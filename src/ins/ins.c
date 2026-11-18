/*
 * ins.c - lifecycle of the loosely-coupled GNSS/SINS module.
 *
 * ins_start() only opens the IMU file; everything else is driven by ins_feed(),
 * which is called from postpos.c once per valid GNSS solution:
 *
 *   init (first valid solution)  ->  seek the IMU file to that epoch, align,
 *                                    build P0/Qt (eskf_initP / eskf_setQt)
 *   every later solution         ->  run_to(t)  : IMU mechanization + time update
 *                                    measurement: position (+ velocity if valid)
 *                                    eskf_meas_update + eskf_feedback
 *                                    one .ins line
 *
 * MSF drives the filter from its own resampled IMU stream (t_gsinskf::_ins_init
 * calls _imudata->interpolate()); here the IMU keeps its native cadence and the
 * residual timing error is absorbed by the ae*dt / ve*dt terms of _set_meas
 * (gins.cpp:951,953), which is what those terms are for.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "rtkins.h"

sins_t   ins_sins;
eskf_t   ins_kf;
insopt_t ins_opt;

#define INS_TMAX 0.1            /* implausible IMU gap -> fall back to nominal ts */

static char ins_msg[256];
static int  ins_on=0;           /* an IMU file was given   */
static int  ins_ready=0;        /* filter initialised      */
static int  ins_err=0;
static FILE *fimu=NULL,*fins=NULL;
static imust_t imust;
static double timu=0.0;         /* time of the last consumed IMU sample */

/* one held-back mechanization group (fetched but not yet propagated) */
static int    ghold=0,gn=0;
static double gt=0.0,gw[INS_IMUN][3],gv[INS_IMUN][3];

/* alignment state, mirrors t_gsinskf::_first_align/_first_pos/_pre_pos/_yaw0 */
static int    first_align=1;
static double first_pos[3],pre_pos[3],align_yaw0;

/* initial velocity when the GNSS file carries no Doppler (sol_t::qv all zero,
 * which is what a position-only .pos gives): the position of an earlier epoch
 * is kept and the two are differenced over INIT_VEL_DT seconds.  MSF always
 * gets a Doppler velocity from t_ggnssdata, so this is a local addition. */
#define INIT_VEL_DT 5.0
static int    pset=0;
static double ppos0[3],ptime0=0.0;

static long nepoch=0,nupdate=0,nreject=0,ndrop=0;

int  ins_enabled(void) { return ins_on&&ins_err==0; }
const char *ins_errmsg(void) { return ins_msg; }

/* ---- IMU group queue ----------------------------------------------------- */
static int fetch_group(double *tg,double w[][3],double v[][3],int *n)
{
    double t[INS_IMUN],wm[INS_IMUN][3],vm[INS_IMUN][3];
    int i,c=0;

    if (ghold) {
        ghold=0;
        *tg=gt; *n=gn;
        for (i=0;i<gn;i++) {
            memcpy(w[i],gw[i],sizeof(double)*3);
            memcpy(v[i],gv[i],sizeof(double)*3);
        }
        return 1;
    }
    for (i=0;i<INS_IMUN;i++) {
        if (imu_next(fimu,&ins_opt,&t[i],wm[i],vm[i])!=1) break;
        c++;
        if (c>=ins_opt.nsub) break;
    }
    if (c<=0) return 0;
    *tg=t[c-1];                     /* the group ends with its last sample */
    for (i=0;i<c;i++) {
        memcpy(w[i],wm[i],sizeof(double)*3);
        memcpy(v[i],vm[i],sizeof(double)*3);
    }
    *n=c;
    return 1;
}
static void hold_group(double tg,const double w[][3],const double v[][3],int n)
{
    int i;
    ghold=1; gt=tg; gn=n;
    for (i=0;i<n;i++) {
        memcpy(gw[i],w[i],sizeof(double)*3);
        memcpy(gv[i],v[i],sizeof(double)*3);
    }
}

/* propagate the mechanization and the filter up to (but not past) t.
 * ins_sins.t and timu must stay equal: every group advances the state by its
 * own tg-timu, so propagating a group that is already behind timu with a
 * substituted nominal ts would leave the state permanently one sample (ts)
 * ahead of the IMU log, which the max-dt gate then rejects at every epoch. */
static void run_to(double t)
{
    double wm[INS_IMUN][3],vm[INS_IMUN][3],tg,ts,dt,phim[3],dvbm[3];
    int n;

    while (fetch_group(&tg,wm,vm,&n)) {
        ts=ins_opt.ts*n;
        dt=tg-timu;
        if (dt<=0.0) continue;             /* stale group, already propagated */
        if (dt>INS_TMAX) dt=ts;            /* implausible gap in the log      */
        if (tg>t+0.5*dt) { hold_group(tg,wm,vm,n); break; }
        imu_comp(wm,vm,n,ins_opt.cps,phim,dvbm,&imust);
        sins_update(&ins_sins,phim,dvbm,dt);
        eskf_time_update(&ins_kf,&ins_sins,dt,ins_opt.infl);
        timu=tg;
    }
}
/* discard the IMU log before the first GNSS epoch, hold the group spanning it.
 * The filter is initialised AT tw, so the group ending at tw itself belongs to
 * the past and is discarded; only a group ending after tw can be propagated. */
static int imu_seek(double t)
{
    double wm[INS_IMUN][3],vm[INS_IMUN][3],tg;
    int n;

    while (fetch_group(&tg,wm,vm,&n)) {
        if (tg>t+0.5*ins_opt.ts) { hold_group(tg,wm,vm,n); return 1; }
    }
    return 0;
}

/* ---- measurement extraction from an RTKLIB solution --------------------- */
/* sol_t::qr and sol_t::qv are float[6] (rtklib.h), never cast them in place */
static void q6d(const float *q,double d[6])
{
    int i;
    for (i=0;i<6;i++) d[i]=q[i];
}
/* diag(C*Q*C') from the packed {vx,vy,vz,vxy,vyz,vzx} RTKLIB covariance */
static void cov6_diag_rot(const double q6[6],const double C[9],double d[3])
{
    static const double I3[9]={1,0,0, 0,1,0, 0,0,1};
    double Q[9];
    int i,j,k;

    Q[0]=q6[0]; Q[4]=q6[1]; Q[8]=q6[2];
    Q[1]=Q[3]=q6[3]; Q[5]=Q[7]=q6[4]; Q[2]=Q[6]=q6[5];
    if (!C) C=I3;
    for (i=0;i<3;i++) {
        d[i]=0.0;
        for (j=0;j<3;j++) for (k=0;k<3;k++) d[i]+=C[i*3+j]*C[i*3+k]*Q[j*3+k];
    }
}
/* absolute ecef position and velocity of a solution */
static void sol_ecef(const sol_t *sol,const double *rb,double pos[3],double vel[3])
{
    int i;

    if (sol->type==0) {
        for (i=0;i<3;i++) { pos[i]=sol->rr[i]; vel[i]=sol->rr[3+i]; }
    } else {
        double blh[3]={0},e[3];
        if (rb) {
            /* enu2ecef only rotates (rtkcmn.c:2188 adds no reference point, and
             * ecef2enu:2174 does subtract it), so the baseline must be put back
             * on top of the reference position here. */
            ecef2pos(rb,blh);
            enu2ecef(blh,sol->rr,e);
            for (i=0;i<3;i++) pos[i]=rb[i]+e[i];
            enu2ecef(blh,sol->rr+3,e);
            for (i=0;i<3;i++) vel[i]=e[i];
        } else {          /* no reference point: keep the numbers as they are */
            for (i=0;i<3;i++) { pos[i]=sol->rr[i]; vel[i]=sol->rr[3+i]; }
        }
    }
}
/* n-frame (east,north,up) velocity of a solution.
 * eth.Cen is the n->e matrix (its columns are the E/N/U ecef unit vectors, the
 * same numbers MSF's t_gbase::Cen returns), so e->n needs the transpose.  For
 * soltype=enu-baseline the velocity is already east-north-up, but in the ENU of
 * the reference position, so rotate it through ecef to the rover ENU (RTKLIB's
 * xyz2enu/ecef2enu/enu2ecef use the identical matrix, rtkcmn.c). */
static void sol_vn(const sol_t *sol,const double *rb,const double lla[3],double vn[3])
{
    double zero[3]={0};
    earth_t eth;

    earth_update(&eth,lla,zero);
    if (sol->type==0) {
        m3_mv(eth.Cne,sol->rr+3,vn);
        return;
    }
    if (rb) {
        double blh[3]={0},ve[3];
        ecef2pos(rb,blh);
        enu2ecef(blh,sol->rr+3,ve);
        m3_mv(eth.Cne,ve,vn);
        return;
    }
    memcpy(vn,sol->rr+3,sizeof(double)*3); /* no reference position to rotate from */
}

/* ---- alignment ----------------------------------------------------------- */
/*
 * INS_ALIGN_AUTO: level frame, yaw so that the body x-axis points along the
 * GNSS velocity.  With n = (east,north,up) and MSF gbase.cpp:16 a2mat, zero
 * pitch and roll give body-x = (cos yaw, sin yaw, 0), hence
 * yaw = atan2(vN,vE) - a mathematical angle from EAST, not a compass heading.
 *
 * INS_ALIGN_VVA / INS_ALIGN_PVA are t_gsinskf::align_vva(gins.cpp:476) and
 * align_pva(gins.cpp:422) transcribed literally, quadrant tables included.
 * DIFF: those tables index the velocity as [0]=north, [1]=east, while the rest
 *       of MSF is east-north-up (gearth.cpp:38, PSINS pp2vn.m).  In the ENU
 *       reading they put body-x 90 deg counter-clockwise of the motion, so use
 *       align = auto unless you are reproducing a MSF run.  Both are kept so the
 *       difference can be measured against the truth file instead of argued.
 *
 * att is returned in degrees, the unit sins_init / sins_align_att expect.
 */
static int init_att(const double vn[3],const double lla[3],double att[3])
{
    double bl[3],dist,e=vn[0],n=vn[1];

    att[0]=att[1]=att[2]=0.0;
    *ins_msg='\0';
    switch (ins_opt.align) {
    case INS_ALIGN_ATT:
        att[0]=ins_opt.att0[0]; att[1]=ins_opt.att0[1]; att[2]=ins_opt.att0[2];
        return 1;
    case INS_ALIGN_VVA:
        if (sqrt(e*e+n*n)<=ins_opt.vtmin) {
            sprintf(ins_msg,"VVA waiting: |vn|=%.2f <= %.2f m/s",sqrt(e*e+n*n),ins_opt.vtmin);
            return 0;
        }
        att[2]=atan2(fabs(e),fabs(n))/D2R;
        if (n>0.0&&e>0.0) att[2]=-att[2];
        if (n<0.0&&e>0.0) att[2]=-(180.0-att[2]);
        if (n<0.0&&e<0.0) att[2]=  (180.0-att[2]);
        return 1;
    case INS_ALIGN_PVA:
        if (first_align) {
            memcpy(first_pos,lla,sizeof(double)*3);
            memcpy(pre_pos ,lla,sizeof(double)*3);
            first_align=0; align_yaw0=0.0;
            sprintf(ins_msg,"PVA waiting for a %.1f m baseline",ins_opt.blmin);
            return 0;
        }
        bl[0]=lla[0]-pre_pos[0]; bl[1]=lla[1]-pre_pos[1]; bl[2]=lla[2]-pre_pos[2];
        memcpy(pre_pos,lla,sizeof(double)*3);
        dist=sqrt(bl[0]*bl[0]+bl[1]*bl[1]);
        if (dist<=ins_opt.blmin/111319.49) {  /* blmin is metres, lla is radians */
            sprintf(ins_msg,"PVA waiting: baseline %.1f m",dist*111319.49);
            return 0;
        }
        att[2]=acos(fabs(bl[1])/dist)/D2R;
        if (bl[1]>0.0&&bl[0]>0.0) att[2]=-att[2];
        if (bl[1]<0.0&&bl[0]>0.0) att[2]=-(180.0-att[2]);
        if (bl[1]<0.0&&bl[0]<0.0) att[2]=  (180.0-att[2]);
        if (align_yaw0==0.0) {
            align_yaw0=att[2]*D2R;
            sprintf(ins_msg,"PVA first baseline, yaw=%.2f deg",att[2]);
            return 0;
        }
        if (fabs(att[2]*D2R-align_yaw0)<10.0*D2R) return 1;
        align_yaw0=0.0; first_align=1;
        sprintf(ins_msg,"PVA failed, keep straight");
        return 0;
    default:
        if (sqrt(e*e+n*n)<=ins_opt.vtmin) {
            if (ins_opt.att0[0]!=0.0||ins_opt.att0[1]!=0.0||ins_opt.att0[2]!=0.0) {
                att[0]=ins_opt.att0[0]; att[1]=ins_opt.att0[1]; att[2]=ins_opt.att0[2];
                return 1;
            }
            sprintf(ins_msg,"alignment waiting: |vn|=%.2f <= %.2f m/s",
                sqrt(e*e+n*n),ins_opt.vtmin);
            return 0;
        }
        att[2]=atan2(n,e)/D2R;
        return 1;
    }
}

/* ---- initialisation ------------------------------------------------------ */
static int init_filter(const sol_t *sol,const double *rb,double tw)
{
    double pos[3],vel[3],lla[3],vn[3],att[3];
    int i;
    
    sol_ecef(sol,rb,pos,vel);
    ecef2pos(pos,lla);
    sol_vn(sol,rb,lla,vn);
    if (!(sol->qv[0]>0.0&&sol->qv[1]>0.0&&sol->qv[2]>0.0)) {   /* no Doppler */
        double dt,bl[3];
        earth_t eth;
        double zero[3]={0};

        if (!pset) {
            pset=1;
            memcpy(ppos0,pos,sizeof(double)*3);
            ptime0=tw;
        }
        dt=tw-ptime0;
        if (dt<INIT_VEL_DT) {
            sprintf(ins_msg,"init velocity waiting: %.1f < %.1f s of positions",
                dt,INIT_VEL_DT);
            return 0;
        }
        for (i=0;i<3;i++) bl[i]=(pos[i]-ppos0[i])/dt;
        earth_update(&eth,lla,zero);
        m3_mv(eth.Cne,bl,vn);
    }
    if (!init_att(vn,lla,att)) return 0;

    /* MSF takes the lever arm off the GNSS position before the filter starts:
     * t_gintegration::_merge_init (gintegration.cpp:626) writes
     *   sins.pos_ecef = gnss_pos - sins.eth.Cen * sins.Cnb * lever
     * because the SINS state is the IMU centre while the GNSS solution is the
     * antenna.  Skipping it leaves a constant Ceb*lever (~0.43 m on campus-01)
     * in every position innovation, and with the cm-level R of a fixed RTK
     * solution the pos-att block [Ceb*lever x] of set_Hk absorbs it as a
     * rotation of up to |Z|/|lever| rad. */
    {
        double Cnb[9],Cen[9],lev_n[3],lev_e[3],zero[3]={0},rad[3];
        earth_t eth;

        for (i=0;i<3;i++) rad[i]=att[i]*D2R;
        earth_update(&eth,lla,zero);
        memcpy(Cen,eth.Cen,sizeof(Cen));   /* n -> e */
        a2m(rad,Cnb);                      /* b -> n */
        m3_mv(Cnb,ins_opt.lever,lev_n);
        m3_mv(Cen,lev_n,lev_e);
        for (i=0;i<3;i++) pos[i]-=lev_e[i];
        ecef2pos(pos,lla);
    }

    if (!imu_seek(tw)) {
        sprintf(ins_msg,"IMU log ends before the first GNSS epoch %.3f",tw);
        ins_err=1;
        return 0;
    }
    sins_init(&ins_sins,&ins_opt,lla,vn,att);
    ins_sins.t=tw;
    timu=tw;

    eskf_initP(&ins_kf,&ins_opt,&ins_sins);
    eskf_setQt(&ins_kf,&ins_opt);

    if (*ins_opt.outfile) {
        if (!(fins=fopen(ins_opt.outfile,"w"))) {
            sprintf(ins_msg,"cannot write %s",ins_opt.outfile);
            ins_err=1;
            return 0;
        }
        fprintf(fins,"# week tow(s) lat(deg) lon(deg) hgt(m) pitch(deg) roll(deg)"
                     " yaw(deg) heading(deg) ve vn vu(m/s) update nr\n");
    }
    ins_ready=1;
    fprintf(stderr,"ins : aligned at %.3f s, lat=%.7f lon=%.7f hgt=%.3f "
        "att=%.3f %.3f %.3f deg, vn=%.3f %.3f %.3f m/s\n",tw,lla[0]*R2D,lla[1]*R2D,lla[2],
        att[0],att[1],att[2],vn[0],vn[1],vn[2]);
    return 1;
}

/* ---- one GNSS epoch ------------------------------------------------------ */
int ins_feed(const sol_t *sol,const double *rb)
{
    double tw,pos[3],vel[3],mpos[3],mvel[3],qp[3],qv[3],C[9],att[3];
    double qp6[6],qv6[6];
    int week=0,has_pos,has_vel,irc,i;

    if (!ins_enabled()||!sol) return 0;
    if (sol->stat==SOLQ_NONE) return 0;

    tw=time2gpst(sol->time,&week);
    nepoch++;

    if (!ins_ready) {
        if (!init_filter(sol,rb,tw)) { ndrop++; return 0; }
    } else run_to(tw);

    sol_ecef(sol,rb,pos,vel);
    memcpy(mpos,pos,sizeof(double)*3);
    memcpy(mvel,vel,sizeof(double)*3);

    q6d(sol->qr,qp6);
    q6d(sol->qv,qv6);
    if (sol->type==0) {
        cov6_diag_rot(qp6,NULL,qp);
        cov6_diag_rot(qv6,NULL,qv);
    } else {
        double lla2[3],zero[3]={0};
        earth_t eth;
        ecef2pos(pos,lla2);
        earth_update(&eth,lla2,zero);
        memcpy(C,eth.Cen,sizeof(C));       /* (e,n,u) -> ecef */
        cov6_diag_rot(qp6,C,qp);
        cov6_diag_rot(qv6,C,qv);
    }
    has_pos=qp[0]>0.0&&qp[1]>0.0&&qp[2]>0.0;
    has_vel=qv[0]>0.0&&qv[1]>0.0&&qv[2]>0.0;
    /* GREAT-MSF takes the measurement noise from its own configuration
     * (t_gsetign pos_noise/vel_noise) rather than from the solution file; the
     * two meas-*-std keys do the same here, 0 = keep the solution covariance */
    if (ins_opt.rpos[0]>0.0||ins_opt.rpos[1]>0.0||ins_opt.rpos[2]>0.0) {
        for (i=0;i<3;i++) qp[i]=ins_opt.rpos[i]*ins_opt.rpos[i];
        has_pos=1;
    }
    if (ins_opt.rvel[0]>0.0||ins_opt.rvel[1]>0.0||ins_opt.rvel[2]>0.0) {
        for (i=0;i<3;i++) qv[i]=ins_opt.rvel[i]*ins_opt.rvel[i];
        has_vel=1;
    }
    switch (ins_opt.meas) {
    case INS_MEAS_POS:     has_vel=0; break;
    case INS_MEAS_VEL:     has_pos=0; break;
    case INS_MEAS_POS_VEL: has_pos=has_vel=1; break;
    default:               break;       /* INS_MEAS_AUTO */
    }
    if (!has_pos&&!has_vel) {
        sprintf(ins_msg,"epoch %.3f unused: no usable solution covariance",tw);
        ndrop++;
        return 0;
    }
    ins_kf.Rpos[0]=qp[0]; ins_kf.Rpos[1]=qp[1]; ins_kf.Rpos[2]=qp[2];
    ins_kf.Rvel[0]=qv[0]; ins_kf.Rvel[1]=qv[1]; ins_kf.Rvel[2]=qv[2];

    irc=eskf_meas_update(&ins_kf,&ins_sins,&ins_opt,mpos,mvel,has_pos,has_vel,tw);
    if (irc>0) {
        eskf_feedback(&ins_kf,&ins_sins,ins_opt.fbsign);
        nupdate++;
    } else if (irc<0) {
        nreject++;
        ndrop++;
        /* rows are [vel(3); pos(3)] when both are present and are e-frame
         * (ECEF) components, so row 2 is the up component of the position */
        sprintf(ins_msg,"epoch %.3f rejected: row %d nr=%.1f > %.1f,"
            " pre-fit |r|=%.3f sqrt(S)=%.3f",tw,ins_kf.gate_row,ins_kf.gate_nr,
            ins_opt.maxnr,fabs(ins_kf.gp_r[ins_kf.gate_row]),
            ins_kf.gp_s[ins_kf.gate_row]);
    } else {
        ndrop++;
        sprintf(ins_msg,"epoch %.3f dropped: |t_ins-t_gnss|=%.3f > max-dt %.3f",
            tw,ins_sins.t-tw,ins_opt.maxdt);
    }
    if (!fins) return 1;

    m2a(ins_sins.Cnb,att);
    /* vn is already (east,north,up): the n-frame of this module IS the local
     * enu frame, so the .ins columns are the state itself, not a rotation.
     * The heading column is the azimuth of the body x-axis, atan2(east,north)
     * of its ENU components: valid for any mounting, whereas 90deg - yaw only
     * is for a roll=0 body (see the frame note in rtkins.h). */
    fprintf(fins,"%5d %10.3f %11.7f %12.8f %9.3f %9.4f %9.4f %9.4f %9.4f"
        " %8.3f %8.3f %8.3f %d %d\n",
        week,tw,ins_sins.pos[0]*R2D,ins_sins.pos[1]*R2D,ins_sins.pos[2],
        att[0]*R2D,att[1]*R2D,att[2]*R2D,
        fmod(atan2(ins_sins.Cnb[0],ins_sins.Cnb[3])*R2D+720.0,360.0),
        ins_sins.vn[0],ins_sins.vn[1],ins_sins.vn[2],irc>0?1:0,ins_kf.nr);
    return 1;
}

/* ---- start / stop -------------------------------------------------------- */
int ins_start(const insopt_t *opt)
{
    ins_on=ins_ready=ins_err=0;
    ghold=0; gn=0;
    first_align=1; align_yaw0=0.0;
    pset=0; ptime0=0.0;
    nepoch=nupdate=nreject=ndrop=0;
    *ins_msg='\0';
    memset(&imust,0,sizeof(imust));
    imust.pf1=imust.pf2=1;              /* one-plus-previous latch */

    ins_opt=*opt;
    if (!*ins_opt.imufile) return 0;    /* no IMU -> module stays completely off */
    if (ins_opt.nsub<1||ins_opt.nsub>INS_IMUN) {
        sprintf(ins_msg,"imu-nsub must be 1..%d",INS_IMUN);
        ins_err=1;
        return -1;
    }
    if (!(fimu=fopen(ins_opt.imufile,"r"))) {
        sprintf(ins_msg,"cannot open IMU file %s",ins_opt.imufile);
        ins_err=1;
        return -1;
    }
    ins_on=1;
    fprintf(stderr,"ins : imu=%s ts=%.4f s inc=%d nsub=%d align=%d\n",
        ins_opt.imufile,ins_opt.ts,ins_opt.inc,ins_opt.nsub,ins_opt.align);
    return 1;
}

void ins_stop(void)
{
    if (fimu) { fclose(fimu); fimu=NULL; }
    if (fins) {
        fprintf(fins,"# %ld solutions, %ld filter updates, %ld rejected, %ld unused\n",
            nepoch,nupdate,nreject,ndrop);
        fclose(fins); fins=NULL;
    }
    if (ins_on) {
        fprintf(stderr,"ins : %ld gnss epochs, %ld updates, %ld outliers, %ld unused\n",
            nepoch,nupdate,nreject,ndrop);
        if (*ins_msg) fprintf(stderr,"ins : last: %s\n",ins_msg);
    }
    ins_on=ins_ready=0;
}
