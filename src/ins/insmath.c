/*
 * insmath.c - small fixed-size linear algebra, quaternion helpers and the
 *             earth-parameter model used by the mechanization.
 *
 * earth_update()/earth_v2dp() are transcriptions of GREAT-MSF
 * LibGREAT/gins/gearth.cpp:32 and :54 (which is PSINS ethupdate.m / Mpv).
 * DIFF: MSF uses f = 1/298.257 (t_gglv::f), RTKLIB's FE_WGS84 = 1/298.257223563
 *       is used here so that the ellipsoid matches rtkcmn.c's pos2ecef/ecef2pos.
 */
#include <math.h>
#include <string.h>
#include "rtkins.h"

#define INS_PI    3.14159265358979
#define INS_DEG   (INS_PI/180.0)
#define INS_A     RE_WGS84                    /* 6378137.0 */
#define INS_F     FE_WGS84                    /* 1/298.257223563 */
#define INS_OMEGA OMGE                        /* 7.2921151467e-5 */
#define INS_G0    9.7803267714                /* normal gravity at equator */

double v_norm(const double v[3])
{
    return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

void m3_eye(double A[9])
{
    int i;
    for (i=0;i<9;i++) A[i]=0.0;
    A[0]=A[4]=A[8]=1.0;
}

/* C = A*B, row-major */
void m3_mul(const double A[9], const double B[9], double C[9])
{
    double t[9];
    int i,j,k;
    for (i=0;i<3;i++) for (j=0;j<3;j++) {
        t[i*3+j]=0.0;
        for (k=0;k<3;k++) t[i*3+j]+=A[i*3+k]*B[k*3+j];
    }
    memcpy(C,t,sizeof(t));
}

void m3_trs(const double A[9], double At[9])
{
    double t[9];
    int i,j;
    for (i=0;i<3;i++) for (j=0;j<3;j++) t[i*3+j]=A[j*3+i];
    memcpy(At,t,sizeof(t));
}

void m3_mv(const double A[9], const double v[3], double w[3])
{
    double t[3];
    int i,k;
    for (i=0;i<3;i++) {
        t[i]=0.0;
        for (k=0;k<3;k++) t[i]+=A[i*3+k]*v[k];
    }
    w[0]=t[0]; w[1]=t[1]; w[2]=t[2];
}

void m3_tv(const double A[9], const double v[3], double w[3])
{
    double t[3];
    int i,k;
    for (i=0;i<3;i++) {
        t[i]=0.0;
        for (k=0;k<3;k++) t[i]+=A[k*3+i]*v[k];
    }
    w[0]=t[0]; w[1]=t[1]; w[2]=t[2];
}

/* askew (gbase.cpp:195) */
void m3_skew(const double v[3], double S[9])
{
    S[0]=0.0;    S[1]=-v[2];  S[2]=v[1];
    S[3]=v[2];   S[4]=0.0;    S[5]=-v[0];
    S[6]=-v[1];  S[7]=v[0];   S[8]=0.0;
}

/* C(m x p) = A(m x n) * B(n x p), row-major */
int mm(int m, int n, int p, const double *A, const double *B, double *C)
{
    double *t;
    int i,j,k;
    if (!(t=(double *)malloc(sizeof(double)*m*p))) return 0;
    for (i=0;i<m;i++) for (j=0;j<p;j++) {
        t[i*p+j]=0.0;
        for (k=0;k<n;k++) t[i*p+j]+=A[i*n+k]*B[k*p+j];
    }
    memcpy(C,t,sizeof(double)*m*p);
    free(t);
    return 1;
}

/* w(m) = A(m x n)' * v(m) */
int mtv(int m, int n, const double *A, const double *v, double *w)
{
    double *t;
    int i,j;
    if (!(t=(double *)malloc(sizeof(double)*n))) return 0;
    for (j=0;j<n;j++) {
        t[j]=0.0;
        for (i=0;i<m;i++) t[j]+=A[i*n+j]*v[i];
    }
    memcpy(w,t,sizeof(double)*n);
    free(t);
    return 1;
}

/* ethupdate.m / gearth.cpp:32 */
void earth_update(earth_t *e, const double pos[3], const double vn[3])
{
    static const double e2 = INS_F*(2.0 - INS_F);
    double sb,cb,tb,sl,cl,sq,sqr,sb2,gn[3],sw[3],S[9];
    int i;

    memcpy(e->pos,pos,sizeof(double)*3);
    sb=sin(pos[0]); cb=cos(pos[0]); tb=sb/cb;
    sl=sin(pos[1]); cl=cos(pos[1]);
    sb2=sb*sb;
    e->Cen[0]=-sl;      e->Cen[1]=-sb*cl; e->Cen[2]= cb*cl;
    e->Cen[3]= cl;      e->Cen[4]=-sb*sl; e->Cen[5]= cb*sl;
    e->Cen[6]= 0.0;     e->Cen[7]= cb;    e->Cen[8]= sb;
    m3_trs(e->Cen,e->Cne);

    sq = 1.0 - e2*sb2;
    sqr= sqrt(sq);
    e->RMh = INS_A*(1.0 - e2)/sq/sqr + pos[2];
    e->RNh = INS_A/sqr + pos[2];

    e->wnie[0]=0.0; e->wnie[1]=INS_OMEGA*cb; e->wnie[2]=INS_OMEGA*sb;
    m3_mv(e->Cen,e->wnie,e->weie);

    e->wnen[0]=-vn[1]/e->RMh;
    e->wnen[1]= vn[0]/e->RNh;
    e->wnen[2]= vn[0]/e->RNh*tb;
    for (i=0;i<3;i++) e->wnin[i]=e->wnie[i]+e->wnen[i];

    e->g0 = INS_G0*(1.0 + 5.27094e-3*sb2 + 2.32718e-5*sb2*sb2) - 3.086e-6*pos[2];
    gn[0]=0.0; gn[1]=0.0; gn[2]=-e->g0;
    /* gcc = gn - (wnie + wnin) x vn   ==  gn - (2*wie + wen) x vn */
    for (i=0;i<3;i++) sw[i]=e->wnie[i]+e->wnin[i];
    m3_skew(sw,S);
    m3_mv(S,vn,e->gcc);
    for (i=0;i<3;i++) e->gcc[i]=gn[i]-e->gcc[i];
}

/* gearth.cpp:54 - velocity to geodetic position rate [dlat/dt; dlon/dt; dh/dt] */
void earth_v2dp(const earth_t *e, const double vn[3], double ts, double dp[3])
{
    double cb=cos(e->pos[0]);
    dp[0]=vn[1]/e->RMh*ts;
    dp[1]=vn[0]/(cb*e->RNh)*ts;
    dp[2]=vn[2]*ts;
}

/* a2mat (gbase.cpp:16), att = {pitch, roll, yaw} */
void a2m(const double att[3], double C[9])
{
    double sp=sin(att[0]),cp=cos(att[0]);
    double sr=sin(att[1]),cr=cos(att[1]);
    double sy=sin(att[2]),cy=cos(att[2]);
    C[0]=cy*cr - sy*sp*sr; C[1]=-sy*cp; C[2]=cy*sr + sy*sp*cr;
    C[3]=sy*cr + cy*sp*sr; C[4]= cy*cp; C[5]=sy*sr - cy*sp*cr;
    C[6]=-cp*sr;           C[7]= sp;    C[8]= cp*cr;
}

/* m2att (gbase.cpp:29) */
void m2a(const double C[9], double att[3])
{
    double t;
    t=C[2*3+1];
    if (t>1.0) t=1.0; else if (t<-1.0) t=-1.0;
    att[0]=asin(t);
    att[1]=atan2(-C[2*3+0],C[2*3+2]);
    att[2]=atan2(-C[0*3+1],C[1*3+1]);
}

/* Hamilton product */
void q_mul(const double a[4], const double b[4], double c[4])
{
    double t[4];
    t[0]=a[0]*b[0]-a[1]*b[1]-a[2]*b[2]-a[3]*b[3];
    t[1]=a[0]*b[1]+a[1]*b[0]+a[2]*b[3]-a[3]*b[2];
    t[2]=a[0]*b[2]+a[2]*b[0]+a[3]*b[1]-a[1]*b[3];
    t[3]=a[0]*b[3]+a[3]*b[0]+a[1]*b[2]-a[2]*b[1];
    c[0]=t[0]; c[1]=t[1]; c[2]=t[2]; c[3]=t[3];
}

void q_conj(const double q[4], double c[4])
{
    c[0]=q[0]; c[1]=-q[1]; c[2]=-q[2]; c[3]=-q[3];
}

void q_norm(double q[4])
{
    double n=sqrt(q[0]*q[0]+q[1]*q[1]+q[2]*q[2]+q[3]*q[3]);
    if (n<=0.0) { q[0]=1.0; q[1]=q[2]=q[3]=0.0; return; }
    q[0]/=n; q[1]/=n; q[2]/=n; q[3]/=n;
}

/* q * v * q' - t_gquat::operator*(Vector3d) (gquat.cpp:77) */
void q_mv(const double q[4], const double v[3], double w[3])
{
    double qm[9],t[3];
    q2m(q,qm);
    m3_mv(qm,v,t);
    w[0]=t[0]; w[1]=t[1]; w[2]=t[2];
}

/* rv2q (gbase.cpp:71), series expansion keeps it exact for tiny rotation vectors */
void rv2q(const double rv[3], double q[4])
{
    double n2=v_norm(rv), c, f, n4;
    n2=n2*n2;
    if (n2 < INS_PI/180.0*INS_PI/180.0) {
        n4=n2*n2;
        c=1.0 - n2*(1.0/8.0) + n4*(1.0/384.0);
        f=0.5 - n2*(1.0/48.0) + n4*(1.0/3840.0);
    } else {
        double h=sqrt(n2)/2.0;
        c=cos(h);
        f=sin(h)/h*0.5;
    }
    q[0]=c; q[1]=f*rv[0]; q[2]=f*rv[1]; q[3]=f*rv[2];
}

/* q2mat (gbase.cpp:140) */
void q2m(const double q[4], double C[9])
{
    double q11=q[0]*q[0],q12=q[0]*q[1],q13=q[0]*q[2],q14=q[0]*q[3],
           q22=q[1]*q[1],q23=q[1]*q[2],q24=q[1]*q[3],
           q33=q[2]*q[2],q34=q[2]*q[3],q44=q[3]*q[3];
    C[0]=q11+q22-q33-q44; C[1]=2*(q23-q14); C[2]=2*(q24+q13);
    C[3]=2*(q23+q14);     C[4]=q11-q22+q33-q44; C[5]=2*(q34-q12);
    C[6]=2*(q24-q13);     C[7]=2*(q34+q12); C[8]=q11-q22-q33+q44;
}

/* m2qua (gbase.cpp:112) */
void m2q(const double C[9], double q[4])
{
    double q0,q1,q2,q3,qq4;
    if (C[0]>=C[4]+C[8]) {
        q1=0.5*sqrt(1.0 + C[0]-C[4]-C[8]); qq4=4.0*q1;
        q0=(C[7]-C[5])/qq4; q2=(C[1]+C[3])/qq4; q3=(C[2]+C[6])/qq4;
    } else if (C[4]>=C[0]+C[8]) {
        q2=0.5*sqrt(1.0 - C[0]+C[4]-C[8]); qq4=4.0*q2;
        q0=(C[2]-C[6])/qq4; q1=(C[1]+C[3])/qq4; q3=(C[5]+C[7])/qq4;
    } else if (C[8]>=C[0]+C[4]) {
        q3=0.5*sqrt(1.0 - C[0]-C[4]+C[8]); qq4=4.0*q3;
        q0=(C[3]-C[1])/qq4; q1=(C[2]+C[6])/qq4; q2=(C[5]+C[7])/qq4;
    } else {
        q0=0.5*sqrt(1.0 + C[0]+C[4]+C[8]); qq4=4.0*q0;
        q1=(C[7]-C[5])/qq4; q2=(C[2]-C[6])/qq4; q3=(C[3]-C[1])/qq4;
    }
    q[0]=q0; q[1]=q1; q[2]=q2; q[3]=q3;
    q_norm(q);
}
