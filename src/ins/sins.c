/*
 * sins.c - IMU error compensation (coning / sculling) and the three SINS
 *          propagation equations.
 *
 * Coefficients are transcribed from GREAT-MSF LibGREAT/gins/gimu.cpp
 *   phi_cone_crt(27) phi_poly_crt(52) vm_cone_crt(82) vm_poly_crt(104) Update(131)
 * Propagation is a line-by-line port of t_gsins::Update (gins.cpp:107-143),
 * i.e. of PSINS insupdate.m.
 */
#include <math.h>
#include <string.h>
#include "rtkins.h"

static void cross(const double a[3], const double b[3], double c[3])
{
    double t[3];
    t[0]=a[1]*b[2]-a[2]*b[1];
    t[1]=a[2]*b[0]-a[0]*b[2];
    t[2]=a[0]*b[1]-a[1]*b[0];
    c[0]=t[0]; c[1]=t[1]; c[2]=t[2];
}

/* coning compensation, cone mode (gimu.cpp:27) */
static void phi_cone(const double wm[][3], int n, double dphi[3])
{
    double c[3];
    dphi[0]=dphi[1]=dphi[2]=0.0;
    if (n==2) {
        cross(wm[0],wm[1],c); dphi[0]+=2.0/3.0*c[0]; dphi[1]+=2.0/3.0*c[1]; dphi[2]+=2.0/3.0*c[2];
    } else if (n==3) {
        cross(wm[1],wm[2],c); dphi[0]+=27.0/40.0*c[0], dphi[1]+=27.0/40.0*c[1], dphi[2]+=27.0/40.0*c[2];
        cross(wm[0],wm[2],c); dphi[0]+= 9.0/20.0*c[0], dphi[1]+= 9.0/20.0*c[1], dphi[2]+= 9.0/20.0*c[2];
        cross(wm[0],wm[1],c); dphi[0]+=27.0/40.0*c[0], dphi[1]+=27.0/40.0*c[1], dphi[2]+=27.0/40.0*c[2];
    } else if (n==4) {
        cross(wm[2],wm[3],c); dphi[0]+=232.0/315.0*c[0], dphi[1]+=232.0/315.0*c[1], dphi[2]+=232.0/315.0*c[2];
        cross(wm[1],wm[3],c); dphi[0]+= 46.0/105.0*c[0], dphi[1]+= 46.0/105.0*c[1], dphi[2]+= 46.0/105.0*c[2];
        cross(wm[0],wm[3],c); dphi[0]+= 18.0/ 35.0*c[0], dphi[1]+= 18.0/ 35.0*c[1], dphi[2]+= 18.0/ 35.0*c[2];
        cross(wm[1],wm[2],c); dphi[0]+=178.0/315.0*c[0], dphi[1]+=178.0/315.0*c[1], dphi[2]+=178.0/315.0*c[2];
        cross(wm[0],wm[2],c); dphi[0]+= 46.0/105.0*c[0], dphi[1]+= 46.0/105.0*c[1], dphi[2]+= 46.0/105.0*c[2];
        cross(wm[0],wm[1],c); dphi[0]+=232.0/315.0*c[0], dphi[1]+=232.0/315.0*c[1], dphi[2]+=232.0/315.0*c[2];
    } else if (n==5) {
        cross(wm[3],wm[4],c); dphi[0]+=18575.0/24192.0*c[0], dphi[1]+=18575.0/24192.0*c[1], dphi[2]+=18575.0/24192.0*c[2];
        cross(wm[2],wm[4],c); dphi[0]+= 2675.0/ 6048.0*c[0], dphi[1]+= 2675.0/ 6048.0*c[1], dphi[2]+= 2675.0/ 6048.0*c[2];
        cross(wm[1],wm[4],c); dphi[0]+=11225.0/24192.0*c[0], dphi[1]+=11225.0/24192.0*c[1], dphi[2]+=11225.0/24192.0*c[2];
        cross(wm[0],wm[4],c); dphi[0]+=  125.0/  252.0*c[0], dphi[1]+=  125.0/  252.0*c[1], dphi[2]+=  125.0/  252.0*c[2];
        cross(wm[2],wm[3],c); dphi[0]+= 2575.0/ 6048.0*c[0], dphi[1]+= 2575.0/ 6048.0*c[1], dphi[2]+= 2575.0/ 6048.0*c[2];
        cross(wm[1],wm[3],c); dphi[0]+=  425.0/  672.0*c[0], dphi[1]+=  425.0/  672.0*c[1], dphi[2]+=  425.0/  672.0*c[2];
        cross(wm[0],wm[3],c); dphi[0]+=13975.0/24192.0*c[0], dphi[1]+=13975.0/24192.0*c[1], dphi[2]+=13975.0/24192.0*c[2];
        cross(wm[1],wm[2],c); dphi[0]+= 1975.0/ 3024.0*c[0], dphi[1]+= 1975.0/ 3024.0*c[1], dphi[2]+= 1975.0/ 3024.0*c[2];
        cross(wm[0],wm[2],c); dphi[0]+=  325.0/ 1512.0*c[0], dphi[1]+=  325.0/ 1512.0*c[1], dphi[2]+=  325.0/ 1512.0*c[2];
        cross(wm[0],wm[1],c); dphi[0]+=21325.0/24192.0*c[0], dphi[1]+=21325.0/24192.0*c[1], dphi[2]+=21325.0/24192.0*c[2];
    }
}

/* sculling compensation, cone mode (gimu.cpp:82) */
static void vm_cone(const double wm[][3], const double vm[][3], int n, double dvel[3])
{
    double c[3];
    dvel[0]=dvel[1]=dvel[2]=0.0;
    if (n<2||n>4) return;
#define ACC(k) do { dvel[0]+=(k)*c[0]; dvel[1]+=(k)*c[1]; dvel[2]+=(k)*c[2]; } while (0)
    if (n==2) {
        cross(wm[0],vm[1],c); ACC(2.0/3.0);
        cross(vm[0],wm[1],c); ACC(2.0/3.0);
    } else if (n==3) {
        cross(wm[1],vm[2],c); ACC(27.0/40.0);
        cross(vm[1],wm[2],c); ACC(27.0/40.0);
        cross(wm[0],vm[2],c); ACC( 9.0/20.0);
        cross(vm[0],wm[2],c); ACC( 9.0/20.0);
        cross(wm[0],vm[1],c); ACC(27.0/40.0);
        cross(vm[0],wm[1],c); ACC(27.0/40.0);
    } else { /* n==4 */
        cross(wm[2],vm[3],c); ACC(232.0/315.0);
        cross(vm[2],wm[3],c); ACC(232.0/315.0);
        cross(wm[1],vm[3],c); ACC( 46.0/105.0);
        cross(vm[1],wm[3],c); ACC( 46.0/105.0);
        cross(wm[0],vm[3],c); ACC( 18.0/ 35.0);
        cross(vm[0],wm[3],c); ACC( 18.0/ 35.0);
        cross(wm[1],vm[2],c); ACC(178.0/315.0);
        cross(vm[1],wm[2],c); ACC(178.0/315.0);
        cross(wm[0],vm[2],c); ACC( 46.0/105.0);
        cross(vm[0],wm[2],c); ACC( 46.0/105.0);
        cross(wm[0],vm[1],c); ACC(232.0/315.0);
        cross(vm[0],wm[1],c); ACC(232.0/315.0);
    }
#undef ACC
}

/* coning compensation, polynomial mode (gimu.cpp:52).
 * DIFF: gimu.cpp:73 writes wm[0].cross(wm[3]+wm[1].cross(wm[4])) and :119 writes
 * vm[2]+vm[1].cross(wm[2]) - both are malformed (a vector.cross of a vector).
 * The n=5 polynomial term is implemented here as the intended w0xw3 + w1xw4.
 * n=1 uses the one-plus-previous trick, including MSF's asymmetry that phi_poly
 * never refreshes wm_1 (gimu.cpp:59 only latches on the first call), which is
 * kept so that a cps=2 run reproduces MSF bit for bit.
 */
static void phi_poly(const double wm[][3], int n, double dphi[3], imust_t *st)
{
    double c[3];
    dphi[0]=dphi[1]=dphi[2]=0.0;
#define ACC(k) do { dphi[0]+=(k)*c[0]; dphi[1]+=(k)*c[1]; dphi[2]+=(k)*c[2]; } while (0)
    if (n==1) {
        if (st->pf1) { memcpy(st->wm1,wm[0],sizeof(double)*3); st->pf1=0; }
        cross(st->wm1,wm[0],c); ACC(1.0/12.0);
    } else if (n==2) {
        cross(wm[0],wm[1],c); ACC(2.0/3.0);
    } else if (n==3) {
        cross(wm[0],wm[2],c); ACC(33.0/80.0);
        cross(wm[1],wm[2],c); ACC(57.0/80.0);
        cross(wm[1],wm[0],c); ACC(-57.0/80.0);
    } else if (n==4) {
        cross(wm[0],wm[1],c); ACC(736.0/945.0);
        cross(wm[2],wm[3],c); ACC(736.0/945.0);
        cross(wm[0],wm[2],c); ACC(334.0/945.0);
        cross(wm[1],wm[3],c); ACC(334.0/945.0);
        cross(wm[0],wm[3],c); ACC(526.0/945.0);
        cross(wm[1],wm[2],c); ACC(654.0/945.0);
    } else if (n==5) {
        cross(wm[0],wm[1],c); ACC(123425.0/145152.0);
        cross(wm[3],wm[4],c); ACC(123425.0/145152.0);
        cross(wm[0],wm[2],c); ACC( 34875.0/145152.0);
        cross(wm[2],wm[4],c); ACC( 34875.0/145152.0);
        cross(wm[0],wm[3],c); ACC( 90075.0/145152.0);
        cross(wm[1],wm[4],c); ACC( 90075.0/145152.0);
        cross(wm[0],wm[4],c); ACC( 66625.0/145152.0);
        cross(wm[1],wm[2],c); ACC(103950.0/145152.0);
        cross(wm[2],wm[3],c); ACC(103950.0/145152.0);
        cross(wm[1],wm[3],c); ACC( 55400.0/145152.0);
    }
#undef ACC
}

/* sculling compensation, polynomial mode (gimu.cpp:104) */
static void vm_poly(const double wm[][3], const double vm[][3], int n,
                    double dvel[3], imust_t *st)
{
    double c[3];
    dvel[0]=dvel[1]=dvel[2]=0.0;
#define ACC(k) do { dvel[0]+=(k)*c[0]; dvel[1]+=(k)*c[1]; dvel[2]+=(k)*c[2]; } while (0)
    if (n==1) {
        if (st->pf2) { memcpy(st->wm1,wm[0],sizeof(double)*3);
                       memcpy(st->vm1,vm[0],sizeof(double)*3); st->pf2=0; }
        cross(st->wm1,vm[0],c); ACC(1.0/12.0);
        cross(st->vm1,wm[0],c); ACC(1.0/12.0);
        memcpy(st->wm1,wm[0],sizeof(double)*3);
        memcpy(st->vm1,vm[0],sizeof(double)*3);
    } else if (n==2) {
        cross(wm[0],vm[1],c); ACC(2.0/3.0);
        cross(vm[0],wm[1],c); ACC(2.0/3.0);
    } else if (n==3) {
        cross(wm[0],vm[2],c); ACC(33.0/80.0);
        cross(vm[0],wm[2],c); ACC(33.0/80.0);
        cross(wm[0],vm[1],c); ACC(57.0/80.0);
        cross(vm[0],wm[1],c); ACC(57.0/80.0);
        cross(wm[1],vm[2],c); ACC(57.0/80.0);
        cross(vm[1],wm[2],c); ACC(57.0/80.0);
    } else if (n==4) {
        cross(wm[0],vm[1],c); ACC(736.0/945.0);
        cross(wm[2],vm[3],c); ACC(736.0/945.0);
        cross(vm[0],wm[1],c); ACC(736.0/945.0);
        cross(vm[2],wm[3],c); ACC(736.0/945.0);
        cross(wm[0],vm[2],c); ACC(334.0/945.0);
        cross(wm[1],vm[3],c); ACC(334.0/945.0);
        cross(vm[0],wm[2],c); ACC(334.0/945.0);
        cross(vm[1],wm[3],c); ACC(334.0/945.0);
        cross(wm[0],vm[3],c); ACC(526.0/945.0);
        cross(vm[0],wm[3],c); ACC(526.0/945.0);
        cross(wm[1],vm[2],c); ACC(654.0/945.0);
        cross(vm[1],wm[2],c); ACC(654.0/945.0);
    }
#undef ACC
}

/* gimu.cpp:131 */
void imu_comp(const double wm[][3], const double vm[][3], int n, int cps,
              double phim[3], double dvbm[3], imust_t *st)
{
    double dphi[3],scl[3],wmm[3]={0,0,0},vmm[3]={0,0,0},rot[3],c[3];
    int i;

    for (i=0;i<n;i++) {
        wmm[0]+=wm[i][0]; wmm[1]+=wm[i][1]; wmm[2]+=wm[i][2];
        vmm[0]+=vm[i][0]; vmm[1]+=vm[i][1]; vmm[2]+=vm[i][2];
    }
    dphi[0]=dphi[1]=dphi[2]=0.0;
    scl [0]=scl [1]=scl [2]=0.0;
    if (cps==0) {
        phi_cone(wm,n,dphi);
        vm_cone(wm,vm,n,scl);
    } else if (cps==2) {
        phi_poly(wm,n,dphi,st);
        vm_poly(wm,vm,n,scl,st);
    }
    cross(wmm,vmm,c);
    rot[0]=0.5*c[0]; rot[1]=0.5*c[1]; rot[2]=0.5*c[2];
    for (i=0;i<3;i++) {
        phim[i] =wmm[i]+dphi[i];
        dvbm[i] =vmm[i]+rot[i]+scl[i];
    }
}

/* t_gsins::Update (gins.cpp:107-143) - attitude, velocity and position recursion */
void sins_update(sins_t *s, const double phim_in[3], const double dvbm_in[3], double ts)
{
    double nts,nts2,phim[3],dvbm[3],vn12[3],pos12[3],dp[3];
    double wib[3],fb[3],fn[3],an[3],vn1[3],vsum[3],Cbuf[9],qtmp[4],qh[4];
    double qt1[4],qt2[4],rv[3],t[3],tt[3];
    int i;

    nts=fabs(ts);
    nts2=nts/2.0;

    /* calibration: diagonal scale factors, bias removed as increments */
    for (i=0;i<3;i++) {
        phim[i]=s->Kg[i]*phim_in[i]-s->eb[i]*nts;
        dvbm[i]=s->Ka[i]*dvbm_in[i]-s->db[i]*nts;
    }

    /* half-step extrapolation of vn/pos, only to sample the earth at the middle */
    for (i=0;i<3;i++) vn12[i]=s->vn[i]+s->an[i]*nts2;
    earth_v2dp(&s->eth,vn12,nts2,dp);
    for (i=0;i<3;i++) pos12[i]=s->pos[i]+dp[i];
    earth_update(&s->eth,pos12,vn12);

    for (i=0;i<3;i++) { wib[i]=phim[i]/nts; fb[i]=dvbm[i]/nts; }

    /* web = wib - Cbn*wnie, Cbn of the previous epoch (gins.cpp:124) */
    m3_trs(s->Cnb,Cbuf);
    m3_mv(Cbuf,s->eth.wnie,t);
    for (i=0;i<3;i++) s->web[i]=wib[i]-t[i];

    /* wnb = wib - (qh')*wnin, qh = qnb*rv2q(phim/2)  (gins.cpp:126, insupdate.m:28) */
    for (i=0;i<3;i++) rv[i]=phim[i]*0.5;
    rv2q(rv,qtmp);
    q_mul(s->qnb,qtmp,qh);
    q2m(qh,Cbuf);
    m3_tv(Cbuf,s->eth.wnin,tt);
    for (i=0;i<3;i++) s->wnb[i]=wib[i]-tt[i];

    /* fn = Cnb*fb, previous-epoch attitude (gins.cpp:128) */
    q2m(s->qnb,Cbuf);
    m3_mv(Cbuf,fb,fn);

    /* an = rv2q(-wnin*nts/2) * fn + gcc  (gins.cpp:130) */
    for (i=0;i<3;i++) rv[i]=-s->eth.wnin[i]*nts2;
    rv2q(rv,qtmp);
    q_mv(qtmp,fn,an);
    for (i=0;i<3;i++) an[i]+=s->eth.gcc[i];

    /* velocity recursion (gins.cpp:132) */
    for (i=0;i<3;i++) vn1[i]=s->vn[i]+an[i]*nts;

    /* position recursion, trapezoid - this single line is the dead-reckoning
     * (gins.cpp:133, v2dp == PSINS Mpv) */
    for (i=0;i<3;i++) vsum[i]=s->vn[i]+vn1[i];
    earth_v2dp(&s->eth,vsum,nts2,dp);
    for (i=0;i<3;i++) s->pos[i]+=dp[i];
    if (s->pos[1]> PI) s->pos[1]-=2.0*PI;
    if (s->pos[1]<-PI) s->pos[1]+=2.0*PI;
    memcpy(s->vn,vn1,sizeof(double)*3);
    memcpy(s->an,an,sizeof(double)*3);

    /* attitude recursion: qnb = rv2q(-wnin*nts) * qnb * rv2q(phim)  (gins.cpp:134,
     * identical to the header comment of PSINS qupdt2.m) */
    for (i=0;i<3;i++) rv[i]=-s->eth.wnin[i]*nts;
    rv2q(rv,qtmp);
    q_mul(qtmp,s->qnb,qt1);
    rv2q(phim,qt2);
    q_mul(qt1,qt2,s->qnb);
    q_norm(s->qnb);

    s->nts=nts;
    s->t+=nts;
    s->ins_time+=nts;

    /* derived quantities of the new epoch */
    q2m(s->qnb,s->Cnb);
    m2a(s->Cnb,s->att);
    m3_trs(s->Cnb,s->Cbn);
    m3_mv(s->Cbn,s->vn,s->vb);
    memcpy(s->fb,fb,sizeof(double)*3);
    memcpy(s->fn,fn,sizeof(double)*3);
    memcpy(s->wib,wib,sizeof(double)*3);

    earth_update(&s->eth,s->pos,s->vn);
    pos2ecef(s->pos,s->pos_e);
    m3_mul(s->eth.Cen,s->Cnb,s->Ceb);
    m3_trs(s->Ceb,s->Cbe);
    m2q(s->Ceb,s->qeb);
    m3_mv(s->eth.Cen,s->vn,s->ve);
    m3_mv(s->eth.Cen,s->an,s->ae);
}

/* initialise the mechanization state (t_gsins constructor + set_posvel) */
void sins_init(sins_t *s, const insopt_t *opt, const double pos[3], const double vn[3],
               const double att[3])
{
    int i;
    double a[3],C[9];

    memset(s,0,sizeof(sins_t));
    memcpy(s->pos,pos,sizeof(double)*3);
    memcpy(s->vn,vn,sizeof(double)*3);
    s->t=0.0;

    for (i=0;i<3;i++) a[i]=att[i]*D2R;
    a2m(a,C);
    m2q(C,s->qnb);
    q_norm(s->qnb);

    /* biases and scale factors from configuration */
    for (i=0;i<3;i++) {
        s->eb[i]=opt->eb0[i]*(D2R/3600.0);          /* deg/h  -> rad/s */
        s->db[i]=opt->db0[i]*(9.7803267714e-3);     /* mg     -> m/s^2 */
        s->Kg[i]=opt->Kg0[i];
        s->Ka[i]=opt->Ka0[i];
        s->tauG[i]= opt->gbcorr[i]>0.0? -1.0/(opt->gbcorr[i]*3600.0):0.0;
        s->tauA[i]= opt->abcorr[i]>0.0? -1.0/(opt->abcorr[i]*3600.0):0.0;
    }
    /* mounting: Cvb maps body -> vehicle, only used by lever-arm/NHC terms */
    a2m(opt->inm,C);
    m2q(C,s->qvb);
    memcpy(s->Cvb,C,sizeof(double)*9);
    m3_trs(C,s->Cbv);

    earth_update(&s->eth,s->pos,s->vn);
    q2m(s->qnb,s->Cnb);
    m3_trs(s->Cnb,s->Cbn);
    m2a(s->Cnb,s->att);
    m3_mul(s->eth.Cen,s->Cnb,s->Ceb);
    m3_trs(s->Ceb,s->Cbe);
    m2q(s->Ceb,s->qeb);
    m3_mv(s->Cbn,s->vn,s->vb);
    m3_mv(s->eth.Cen,s->vn,s->ve);
    pos2ecef(s->pos,s->pos_e);
    s->nts=opt->ts;
}

void sins_align_att(sins_t *s, const double att[3])
{
    double a[3],C[9];
    int i;
    for (i=0;i<3;i++) a[i]=att[i]*D2R;
    a2m(a,C);
    m2q(C,s->qnb);
    q_norm(s->qnb);
    q2m(s->qnb,s->Cnb);
    m3_trs(s->Cnb,s->Cbn);
    m2a(s->Cnb,s->att);
    m3_mv(s->Cbn,s->vn,s->vb);
    earth_update(&s->eth,s->pos,s->vn);
    m3_mul(s->eth.Cen,s->Cnb,s->Ceb);
    m3_trs(s->Ceb,s->Cbe);
    m2q(s->Ceb,s->qeb);
}
