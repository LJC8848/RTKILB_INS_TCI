/*
 * eskf.c - 15-state error-state Kalman filter, loosely-coupled GNSS/SINS.
 *
 * Transcribed from GREAT-MSF t_gsinskf:
 *   set_Ft(525) set_Hk(550) _set_meas(941) _resize(1002) time_update(895)
 *   _meas_update(907) feedback(688) _ins_init(773) outlier_detect(1039)
 * No matrix inversion: the measurement update is applied row by row as scalar
 * updates, exactly like MSF.
 *
 * Error convention used throughout (and it is the one MSF's own F matrix implies):
 *   phi = attitude error   = estimated - true   (e-frame)
 *   dve = velocity error   = estimated - true   (e-frame)
 *   dpe = position error   = estimated - true   (e-frame)
 *   eb  = gyro  bias error = estimated - true   (b-frame)   [-Ceb in att row]
 *   db  = accel bias error = estimated - true   (b-frame)   [+Ceb in vel row]
 * so every block is removed by subtraction on feedback.
 *
 * DIFF: MSF adds Xk back for the two bias blocks (gins.cpp:694-695) and rotates
 *       qeb by rv2q(+phi) (gquat.cpp:42 "operator-"), i.e. it *adds* for
 *       attitude and bias while it *subtracts* for position and velocity.
 *       Set fb_sign=+1 in the config to reproduce MSF verbatim and compare.
 */
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "rtkins.h"

/* every diagnostic below is gated by the runtime key "debug" (inscfg.c), so a
 * delivery build stays silent without a second compile of the module. */
#define DBGF(...) do { if (ins_opt.dbg) fprintf(stderr,__VA_ARGS__); } while (0)

#define IXA 0
#define IXV 3
#define IXP 6
#define IXB 9
#define IXC 12

#define SQRT_H 60.0                       /* sqrt(3600 s) */
#define DPH    (D2R/3600.0)               /* deg/h  -> rad/s */
#define MG     (9.7803267714e-3)          /* mg     -> m/s^2 */
#define SQR(x) ((x)*(x))                  /* rtklib.h has no sq() */

static void blkadd(double *M,int r,int c,double s,const double *B)   /* 3x3 block += s*B */
{
    int i,j;
    for (i=0;i<3;i++) for (j=0;j<3;j++) M[(r+i)*INS_NX+c+j]+=s*B[i*3+j];
}
static void blkset(double *M,int r,int c,const double *B)
{
    int i,j;
    for (i=0;i<3;i++) for (j=0;j<3;j++) M[(r+i)*INS_NX+c+j]=B[i*3+j];
}
static void eye15(double *M)
{
    int i;
    memset(M,0,sizeof(double)*INS_NX*INS_NX);
    for (i=0;i<INS_NX;i++) M[i*INS_NX+i]=1.0;
}
static void trs15(const double *A,double *B)
{
    int i,j;
    for (i=0;i<INS_NX;i++) for (j=0;j<INS_NX;j++) B[i*INS_NX+j]=A[j*INS_NX+i];
}
static void mv15(const double *A,const double *v,double *w)          /* w = A*v, A 15x15 */
{
    double t[INS_NX];
    int i,k;
    for (i=0;i<INS_NX;i++) {
        t[i]=0.0;
        for (k=0;k<INS_NX;k++) t[i]+=A[i*INS_NX+k]*v[k];
    }
    memcpy(w,t,sizeof(t));
}

/* gins.cpp:833 - psd conversion, then squared (gins.cpp:840) */
void eskf_setQt(eskf_t *kf, const insopt_t *opt)
{
    int i;
    for (i=0;i<3;i++) {
        kf->Qt[IXA+i]=SQR((opt->mispsd[i]*D2R)/SQRT_H);   /* deg/sqrt(h) -> rad^2/s */
        kf->Qt[IXV+i]=SQR(opt->vrwpsd[i]*MG);             /* mg/sqrt(Hz) -> (m/s^2)^2 */
        kf->Qt[IXP+i]=SQR(opt->pospsd[i]/SQRT_H);         /* m/sqrt(h)   -> m^2/s */
        kf->Qt[IXB+i]=SQR((opt->gbpsd[i]*DPH)/SQRT_H);    /* dph/sqrt(h) -> (rad/s)^2/s */
        kf->Qt[IXC+i]=SQR((opt->abpsd[i]*MG)/SQRT_H);     /* mg/sqrt(h)  -> (m/s^2)^2/s */
    }
}

/* _ins_init (gins.cpp:781-811) */
void eskf_initP(eskf_t *kf, const insopt_t *opt, const sins_t *s)
{
    double Patt[9], CenP[9];
    int i;

    memset(kf,0,sizeof(eskf_t));
    for (i=0;i<3;i++) {
        kf->Pk[(IXA+i)*INS_NX+IXA+i]=SQR(opt->misd[i]*D2R);
        kf->Pk[(IXV+i)*INS_NX+IXV+i]=SQR(opt->vstd[i]);
        kf->Pk[(IXP+i)*INS_NX+IXP+i]=SQR(opt->pstd[i]);
        kf->Pk[(IXB+i)*INS_NX+IXB+i]=SQR(opt->gbstd[i]*DPH);
        kf->Pk[(IXC+i)*INS_NX+IXC+i]=SQR(opt->abstd[i]*MG);
    }
    /* misalignment is configured in n, the state lives in e */
    memcpy(Patt,&kf->Pk[IXA*INS_NX+IXA],sizeof(Patt));
    mm(3,3,3,s->eth.Cen,Patt,CenP);
    mm(3,3,3,CenP,s->eth.Cne,Patt);
    blkset(kf->Pk,IXA,IXA,Patt);
    eye15(kf->Phik);
    kf->init=1;
}

/* set_Ft (gins.cpp:525) */
static void set_Ft(eskf_t *kf, const sins_t *s)
{
    double weieSk[9], cenfn[3], fnSk[9];
    int i;

    memset(kf->Ft,0,sizeof(kf->Ft));
    m3_mv(s->eth.Cen,s->fn,cenfn);
    m3_skew(s->eth.weie,weieSk);
    m3_skew(cenfn,fnSk);

    blkadd(kf->Ft,IXA,IXA,-1.0,weieSk);   /* att-att */
    blkadd(kf->Ft,IXA,IXB,-1.0,s->Ceb);   /* att-eb  */
    blkadd(kf->Ft,IXV,IXA, 1.0,fnSk);     /* vel-att */
    blkadd(kf->Ft,IXV,IXV,-2.0,weieSk);   /* vel-vel */
    blkadd(kf->Ft,IXV,IXC, 1.0,s->Ceb);   /* vel-db  */
    for (i=0;i<3;i++) {                   /* pos-vel = I */
        kf->Ft[(IXP+i)*INS_NX+IXV+i]=1.0;
        kf->Ft[(IXB+i)*INS_NX+IXB+i]=s->tauG[i];
        kf->Ft[(IXC+i)*INS_NX+IXC+i]=s->tauA[i];
    }
}

/* time_update (gins.cpp:895) */
void eskf_time_update(eskf_t *kf, sins_t *s, double dt, double infl)
{
    double Fk[INS_NX*INS_NX], FkT[INS_NX*INS_NX], T[INS_NX*INS_NX], X[INS_NX];
    int i,j;

    set_Ft(kf,s);
    for (i=0;i<INS_NX;i++) {
        for (j=0;j<INS_NX;j++) Fk[i*INS_NX+j]=kf->Ft[i*INS_NX+j]*dt;
        Fk[i*INS_NX+i]+=1.0;
    }
    /* Phik accumulates and is never reset - same as MSF (gins.cpp:899) */
    mm(INS_NX,INS_NX,INS_NX,Fk,kf->Phik,T);
    memcpy(kf->Phik,T,sizeof(T));

    mv15(Fk,kf->Xk,X);
    memcpy(kf->Xk,X,sizeof(X));

    mm(INS_NX,INS_NX,INS_NX,Fk,kf->Pk,T);
    trs15(Fk,FkT);
    mm(INS_NX,INS_NX,INS_NX,T,FkT,kf->Pk);
    for (i=0;i<INS_NX;i++) kf->Pk[i*INS_NX+i]+=kf->Qt[i]*dt*infl;
}

/*
 * set_Hk + _set_meas for the loosely-coupled family (gins.cpp:563-578 and
 * 949-969).  Rows are ordered [vel(3); pos(3)] when both are present, exactly
 * as MSF's POS_VEL_MEAS.
 */
static void build_meas(eskf_t *kf, const sins_t *s, const double *lever,
                       const double *mpos, const double *mvel,
                       int has_pos, int has_vel, double dt)
{
    double leverSk[9], vartheta[3], a[9], CebLev[3], tmp[9], sk[9], H[9];
    double eye[9];
    int i,rv,rp,nr;

    m3_eye(eye);
    memset(kf->Hk,0,sizeof(kf->Hk));
    memset(kf->Rk,0,sizeof(kf->Rk));
    kf->Zk[0]=kf->Zk[1]=kf->Zk[2]=kf->Zk[3]=kf->Zk[4]=kf->Zk[5]=0.0;

    /* vartheta = Ceb*[lever x]*wib + [weie x]*Ceb*lever  (gins.cpp:553) */
    m3_skew(lever,leverSk);
    m3_mul(s->Ceb,leverSk,a);
    m3_mv(a,s->wib,vartheta);
    m3_mv(s->Ceb,lever,CebLev);
    m3_skew(s->eth.weie,sk);
    m3_mv(sk,CebLev,tmp);
    for (i=0;i<3;i++) vartheta[i]+=tmp[i];

    rv=0;
    rp=has_vel?3:0;
    nr=(has_vel?3:0)+(has_pos?3:0);
    kf->nr=nr;

    if (has_vel) {
        for (i=0;i<3;i++) {
            kf->Zk[rv+i]=s->ve[i]-mvel[i]-s->ae[i]*dt-vartheta[i];
            kf->Rk[(rv+i)*INS_NR+rv+i]=kf->Rvel[i];
        }
        blkset(kf->Hk,rv,IXV,eye);                     /* vel-vel = I    */
        m3_skew(vartheta,H);
        blkadd(kf->Hk,rv,IXA,-1.0,H);                  /* vel-att        */
        m3_mul(s->Ceb,leverSk,a);                       /* vel-eb = -Ceb[l x] (gins.cpp:566) */
        for (i=0;i<9;i++) kf->Hk[(rv+i/3)*INS_NX+IXB+i%3]=-a[i];
    }
    if (has_pos) {
        for (i=0;i<3;i++) {
            kf->Zk[rp+i]=s->pos_e[i]-mpos[i]-s->ve[i]*dt+CebLev[i];
            kf->Rk[(rp+i)*INS_NR+rp+i]=kf->Rpos[i];
        }
        blkset(kf->Hk,rp,IXP,eye);                    /* pos-pos = I    */
        m3_skew(CebLev,H);
        blkadd(kf->Hk,rp,IXA,1.0,H);                  /* pos-att        */
    }
}

/*
 * _meas_update (gins.cpp:907): row-by-row scalar updates, then MSF's own
 * posterior normalized-innovation gate (gins.cpp:1039), rollback on failure.
 * return:  1 updated, 0 nothing to do, -1 rejected as outlier
 */
int eskf_meas_update(eskf_t *kf, sins_t *s, const insopt_t *opt,
                     const double mpos[3], const double mvel[3],
                     int has_pos, int has_vel, double tmeas)
{
    double Hi[INS_NX], Pxz[INS_NX], Kk[INS_NX], Pkold[INS_NX*INS_NX], Xold[INS_NX];
    double Kall[INS_NR*INS_NX], Ppre[INS_NX], dt,v,Qz,Pzz,p,res;
    int i,j,row,nok=0;

    if (!has_pos&&!has_vel) return 0;
    dt=s->t-tmeas;
    if (fabs(dt)>opt->maxdt) return 0;

    build_meas(kf,s,opt->lever,mpos,mvel,has_pos,has_vel,dt);

    memcpy(Pkold,kf->Pk,sizeof(Pkold));
    memcpy(Xold,kf->Xk,sizeof(Xold));
    if (opt->dbg) for (i=0;i<INS_NX;i++) Ppre[i]=sqrt(kf->Pk[i*INS_NX+i]);

    for (row=0;row<kf->nr;row++) {
        for (j=0;j<INS_NX;j++) Hi[j]=kf->Hk[row*INS_NX+j];
        mv15(kf->Pk,Hi,Pxz);
        Pzz=kf->Rk[row*INS_NR+row];
        for (j=0;j<INS_NX;j++) Pzz+=Hi[j]*Pxz[j];
        if (Pzz<=0.0) continue;
        res=kf->Zk[row];
        for (j=0;j<INS_NX;j++) res-=Hi[j]*kf->Xk[j];
        kf->gp_r[row]=res;                          /* for the reject message */
        kf->gp_s[row]=sqrt(Pzz);
        for (j=0;j<INS_NX;j++) Kk[j]=Pxz[j]/Pzz;
        if (opt->dbg) for (j=0;j<INS_NX;j++) Kall[row*INS_NX+j]=Kk[j];
        for (j=0;j<INS_NX;j++) kf->Xk[j]+=Kk[j]*res;
        for (i=0;i<INS_NX;i++) for (j=0;j<INS_NX;j++)
            kf->Pk[i*INS_NX+j]-=Kk[i]*Pxz[j];
        nok++;
    }
    if (nok<=0) {
        memcpy(kf->Pk,Pkold,sizeof(Pkold));
        memcpy(kf->Xk,Xold,sizeof(Xold));
        return 0;
    }
    /*
     * outlier_detect (gins.cpp:1039) - MSF tests the POSTERIOR: _meas_update()
     * runs the whole sequential update first (gins.cpp:912-924) and only then
     * compares |Zk - Hk*Xk| / sqrt(Hk*Pk*Hk' + Rk) against max_norm, both taken
     * with the updated Xk and Pk.  A failure rolls both back, which is what
     * t_gintegration::_feedback(type,false) does (_global_variance = _gv_sav,
     * Xk = 0, gintegration.cpp:235-237) - the SINS state is left alone because
     * feedback() then subtracts a zero Xk.
     */
    for (row=0;row<kf->nr;row++) {
        for (j=0;j<INS_NX;j++) Hi[j]=kf->Hk[row*INS_NX+j];
        mv15(kf->Pk,Hi,Pxz);
        Qz=kf->Rk[row*INS_NR+row];
        for (j=0;j<INS_NX;j++) Qz+=Hi[j]*Pxz[j];
        v=kf->Zk[row];
        for (j=0;j<INS_NX;j++) v-=Hi[j]*kf->Xk[j];
        if (opt->maxnr>0.0&&Qz>0.0&&fabs(v)/sqrt(Qz)>opt->maxnr) {
            kf->gate_row=row;
            kf->gate_nr=fabs(v)/sqrt(Qz);
            DBGF("gate t=%.3f row %d: |v|=%.4f m sqrt(S)=%.4f nr=%.2f\n",
                tmeas,row,fabs(v),sqrt(Qz),fabs(v)/sqrt(Qz));
            memcpy(kf->Pk,Pkold,sizeof(Pkold));
            memcpy(kf->Xk,Xold,sizeof(Xold));
            return -1;
        }
    }
    /* symmetry */
    for (i=0;i<INS_NX;i++) for (j=i+1;j<INS_NX;j++) {
        p=0.5*(kf->Pk[i*INS_NX+j]+kf->Pk[j*INS_NX+i]);
        kf->Pk[i*INS_NX+j]=kf->Pk[j*INS_NX+i]=p;
    }
    if (opt->dbg) {
        double zp=0.0;
        for (i=0;i<3;i++) zp+=kf->Zk[i]*kf->Zk[i];
        fprintf(stderr,"meas t=%.3f |Z|=%.4f m Z %8.4f %8.4f %8.4f"
            " R %8.4f %8.4f %8.4f\n",tmeas,sqrt(zp),kf->Zk[0],kf->Zk[1],kf->Zk[2],
            sqrt(kf->Rk[0*INS_NR+0]),sqrt(kf->Rk[1*INS_NR+1]),sqrt(kf->Rk[2*INS_NR+2]));
        fprintf(stderr,"     pre-fit sqrt(S) %.4f %.4f %.4f m\n",
            kf->gp_s[0],kf->gp_s[1],kf->gp_s[2]);
        fprintf(stderr,"     std pre  att %7.3f %7.3f %7.3f deg  pos %7.3f m  gb %6.1f dph  ab %6.2f mg\n",
            Ppre[IXA]*R2D,Ppre[IXA+1]*R2D,Ppre[IXA+2]*R2D,Ppre[IXP],
            Ppre[IXB]*3600.0*R2D,Ppre[IXC]/MG);
        fprintf(stderr,"     std post att %7.3f %7.3f %7.3f deg  pos %7.3f m  gb %6.1f dph  ab %6.2f mg\n",
            sqrt(kf->Pk[IXA*INS_NX+IXA])*R2D,sqrt(kf->Pk[(IXA+1)*INS_NX+IXA+1])*R2D,
            sqrt(kf->Pk[(IXA+2)*INS_NX+IXA+2])*R2D,sqrt(kf->Pk[IXP*INS_NX+IXP]),
            sqrt(kf->Pk[IXB*INS_NX+IXB])*3600.0*R2D,sqrt(kf->Pk[IXC*INS_NX+IXC])/MG);
        for (row=0;row<kf->nr;row++) {
            fprintf(stderr,"     K row %d att %8.4f %8.4f %8.4f deg/m  pos %8.4f %8.4f %8.4f"
                "  gb %7.2f %7.2f %7.2f dph/m\n",row,
                Kall[row*INS_NX+IXA]*R2D,Kall[row*INS_NX+IXA+1]*R2D,
                Kall[row*INS_NX+IXA+2]*R2D,Kall[row*INS_NX+IXP],
                Kall[row*INS_NX+IXP+1],Kall[row*INS_NX+IXP+2],
                Kall[row*INS_NX+IXB]*3600.0*R2D,Kall[row*INS_NX+IXB+1]*3600.0*R2D,
                Kall[row*INS_NX+IXB+2]*3600.0*R2D);
        }
        fprintf(stderr,"     Xk att %7.3f %7.3f %7.3f deg  vel %7.3f %7.3f %7.3f m/s"
            "  gb %7.2f %7.2f %7.2f dph  ab %6.2f %6.2f %6.2f mg\n",
            kf->Xk[IXA]*R2D,kf->Xk[IXA+1]*R2D,kf->Xk[IXA+2]*R2D,
            kf->Xk[IXV],kf->Xk[IXV+1],kf->Xk[IXV+2],
            kf->Xk[IXB]*3600.0*R2D,kf->Xk[IXB+1]*3600.0*R2D,kf->Xk[IXB+2]*3600.0*R2D,
            kf->Xk[IXC]/MG,kf->Xk[IXC+1]/MG,kf->Xk[IXC+2]/MG);
    }
    return 1;
}

/* feedback (gins.cpp:688) - sign=+1 reproduces MSF's add-for-att/bias version */
void eskf_feedback(eskf_t *kf, sins_t *s, int sign)
{
    double phi[3], qtmp[4], qn[4], Cne[9], pos[3], vn[3];
    double att0[3]={0}, Xp0[3]={0}, qa[3], Cq[9], qn0[4], qnb0[4];
    int i;
    if (ins_opt.dbg) {
        memcpy(att0,s->att,sizeof(att0));
        memcpy(Xp0,kf->Xk+IXP,sizeof(Xp0));
        /* is the state self-consistent BEFORE anything is touched?  s->att comes
         * from m2a(Cnb), the rebuild below from m2q(Cne)*qeb, so a difference
         * here means qeb/eth and Cnb had already parted company. */
        m2q(s->eth.Cne,qn0);
        q_mul(qn0,s->qeb,qnb0);
        q_norm(qnb0);
        q2m(qnb0,Cq);
        m2a(Cq,qa);
        fprintf(stderr,"fb?  att from Cnb %8.3f %8.3f %8.3f  from Cne*qeb"
            " %8.3f %8.3f %8.3f deg\n",att0[0]*R2D,att0[1]*R2D,att0[2]*R2D,
            qa[0]*R2D,qa[1]*R2D,qa[2]*R2D);
    }

    for (i=0;i<3;i++) phi[i]=(double)sign*kf->Xk[IXA+i];
    rv2q(phi,qtmp);
    q_mul(qtmp,s->qeb,s->qeb);
    q_norm(s->qeb);

    for (i=0;i<3;i++) {
        s->ve[i]    -= kf->Xk[IXV+i];
        s->pos_e[i] -= kf->Xk[IXP+i];
        s->eb[i]    += (double)sign*kf->Xk[IXB+i];
        s->db[i]    += (double)sign*kf->Xk[IXC+i];
    }

    /* rebuild the n-frame state from the corrected e-frame quantities */
    ecef2pos(s->pos_e,pos);
    memcpy(s->pos,pos,sizeof(double)*3);
    m3_trs(s->eth.Cen,Cne);
    m3_mv(Cne,s->ve,vn);
    earth_update(&s->eth,s->pos,vn);
    m3_mv(s->eth.Cne,s->ve,s->vn);

    m2q(s->eth.Cne,qn);
    q_mul(qn,s->qeb,s->qnb);
    q_norm(s->qnb);

    q2m(s->qeb,s->Ceb);
    m3_trs(s->Ceb,s->Cbe);
    q2m(s->qnb,s->Cnb);
    m3_trs(s->Cnb,s->Cbn);
    m2a(s->Cnb,s->att);
    m3_mv(s->Cbn,s->vn,s->vb);

    memset(kf->Xk,0,sizeof(kf->Xk));
    if (ins_opt.dbg) {
        fprintf(stderr,"fb   att %7.3f %7.3f %7.3f -> %7.3f %7.3f %7.3f deg"
            "  dpos ecef %7.3f %7.3f %7.3f m\n",att0[0]*R2D,att0[1]*R2D,att0[2]*R2D,
            s->att[0]*R2D,s->att[1]*R2D,s->att[2]*R2D,-Xp0[0],-Xp0[1],-Xp0[2]);
    }
}
