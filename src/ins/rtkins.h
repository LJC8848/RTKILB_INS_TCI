#ifndef RTKINS_H
#define RTKINS_H
/*
 * rtkins.h - SINS mechanization + loosely-coupled (LCI) GNSS/INS EKF, plain C.
 *
 * Formulas are transcribed from GREAT-MSF (C++/Eigen) so results are comparable:
 *   LibGREAT/gins/gimu.cpp   -> imu_comp()      coning / sculling compensation
 *   LibGREAT/gins/gins.cpp   -> sins_update()   t_gsins::Update (mechanization)
 *   LibGREAT/gins/gearth.cpp -> earth_update()  t_gearth::Update / v2dp
 *   LibGREAT/gins/gins.cpp   -> eskf_*()        t_gsinskf set_Ft(525) set_Hk(550)
 *                                             _set_meas(941) time_update(895)
 *                                             _meas_update(907) feedback(688)
 * Deviations from MSF are flagged with "DIFF:" in the code.
 */

#include "rtklib.h"

#ifdef __cplusplus
extern "C" {
#endif

#define INS_NX      15          /* error state dimension */
#define INS_NR      6           /* max loosely-coupled measurement dimension */
#define INS_IMUN    8           /* max IMU samples per mechanization cycle (nsub) */

/*
 * FRAME CONVENTION - read this before changing anything.
 * The n-frame here is the SAME one GREAT-MSF and PSINS use, which is
 * EAST-NORTH-UP, not north-east-down:
 *   PSINS pp2vn.m:  vn = [dLon*clRNh; dLat*RMh; dh]/ts   -> (east, north, up)
 *   PSINS ethupdate.m: wnie = [0; wie*cos(lat); wie*sin(lat)], gn(3) = -g
 *   MSF gearth.cpp:38  the three columns of Cen are the E / N / U unit vectors
 * The attitude follows MSF gbase.cpp:16 a2mat, att = {pitch, roll, yaw}, and at
 * pitch=roll=0 the body x-axis points to (cos(yaw), sin(yaw), 0) in (E,N), i.e.
 * yaw is the mathematical angle counter-clockwise from EAST.  A navigation
 * heading (clockwise from north) of the body x-axis is therefore
 * heading = atan2(Cnb[0],Cnb[3]), which equals 90deg - yaw only for roll=0.
 * A forward-right-down (FRD) IMU is expressible too: with roll=180deg and
 * pitch=0 the a2mat columns are (x,y,z) = (-cos yaw,-sin yaw,0),
 * (-sin yaw,cos yaw,0),(0,0,-1), so heading H maps to yaw = 270deg - H.
 * The .ins heading column is computed from Cnb so it is right for both.
 */

#define INS_ALIGN_AUTO  0       /* level, yaw from the GNSS velocity direction */
#define INS_ALIGN_ATT   1       /* take the initial attitude from the config */
#define INS_ALIGN_PVA   2       /* MSF align_pva: yaw from a GNSS position baseline */
#define INS_ALIGN_VVA   3       /* MSF align_vva: yaw from the GNSS velocity */

/*
 * Body frame the IMU file is written in, i.e. the last three letters of MSF's
 * xml <DataFormat><AxisOrder Type="gaXXX">.  The six matrices below are MSF's
 * gcoders/imufile.cpp:105-133 table pre-multiplied by C(flu<-rfu), because MSF
 * decodes into rfu as its canonical body (Type="garfu" -> identity there) while
 * this module's canonical body is the FLU implied by the ATTITUDE note above.
 * DIFF: MSF leaves the data in the declared frame and its align_pva/align_vva
 *       tables still put body-x along the motion, which is only right for a
 *       forward-axis IMU; here the frame is honoured by the reader instead.
 */
#define INS_FRAME_FLU 0         /* forward-left-up : identity (default)         */
#define INS_FRAME_RFU 1         /* right-forward-up: campus-01-MEMS.txt, garfu  */
#define INS_FRAME_FRD 2         /* forward-right-down: ADIS16465.txt            */
#define INS_FRAME_RBD 3         /* right-behind-down                            */
#define INS_FRAME_LBU 4         /* left-behind-up                               */
#define INS_FRAME_BRU 5         /* behind-right-up                              */

#define INS_MEAS_AUTO    0      /* vel from sol->qv if valid, position always */
#define INS_MEAS_POS     1
#define INS_MEAS_VEL     2
#define INS_MEAS_POS_VEL 3

/* error state layout (identical to t_gsinskf::Xk):
 *   0..2   phi   attitude error, misalignment.  P0 is rotated into the e-frame
 *          by Cen in _ins_init (gins.cpp:810) and fed back onto qeb (gins.cpp:691),
 *          so phi is an e-frame error although it is *configured* in the n-frame.
 *   3..5   dve   velocity error  (e-frame, m/s)
 *   6..8   dpe   position error  (e-frame, m)
 *   9..11  eb    gyro bias error (b-frame, rad/s)
 *   12..14 db    accel bias error(b-frame, m/s^2)
 */
enum {IX_ATT = 0, IX_VEL = 3, IX_POS = 6, IX_EB = 9, IX_DB = 12};

typedef struct {                /* earth parameters (t_gearth) */
    double pos[3];              /* lat, lon (rad), hgt (m) */
    double Cen[9];              /* n -> e */
    double Cne[9];              /* e -> n */
    double RMh, RNh;            /* meridian / prime vertical radius of curvature + h */
    double wnie[3];             /* earth rotation rate in n */
    double wnen[3];             /* n-frame rotation rate wrt earth in n */
    double wnin[3];             /* n-frame rotation rate wrt inertial in n */
    double weie[3];             /* earth rotation rate in e */
    double gcc[3];              /* gravity - (wnie + wnin) x vn  (n-frame) */
    double g0;                  /* normal gravity used for gcc (m/s^2) */
} earth_t;

typedef struct {                /* SINS state (t_gsins) */
    double t;                   /* GPS time of week (s) */
    double nts;                 /* last propagation interval (s) */
    double qnb[4];              /* n <- b quaternion */
    double qeb[4];              /* e <- b quaternion */
    double qvb[4];              /* v <- b quaternion (IMU mounting) */
    double pos[3];              /* lat, lon (rad), hgt (m) */
    double vn[3];               /* n-frame velocity (m/s), n = east, north, up */
    double an[3];               /* n-frame acceleration (m/s^2) */
    double ve[3], ae[3];        /* e-frame velocity / acceleration */
    double pos_e[3];            /* ecef position (m) */
    double Cnb[9], Cbn[9], Ceb[9], Cbe[9], Cen[9];
    double Cvb[9], Cbv[9];      /* body <-> vehicle (IMU mounting) */
    double att[3];              /* pitch, roll, yaw (rad), MSF m2att order */
    double fb[3], fn[3];        /* specific force, b / n frame */
    double wib[3], web[3], wnb[3];
    double vb[3];               /* n-frame velocity expressed in body frame */
    double eb[3];               /* gyro bias (rad/s) */
    double db[3];               /* accel bias (m/s^2) */
    double Kg[3], Ka[3];        /* scale factors (1 = no scale error) */
    double tauG[3], tauA[3];    /* -1/correlation time of the biases (1/s) */
    double ins_time;            /* accumulated pure-INS time (s), t_gsins::pure_ins_time */
    earth_t eth;
} sins_t;

typedef struct {                /* Kalman filter (t_gsinskf) */
    double Xk[INS_NX];
    double Pk[INS_NX * INS_NX];
    double Ft[INS_NX * INS_NX];
    double Phik[INS_NX * INS_NX];
    double Qt[INS_NX];          /* process noise, already squared (variance per second) */
    double Hk[INS_NR * INS_NX];
    double Zk[INS_NR];
    double Rk[INS_NR * INS_NR];
    double Rpos[3], Rvel[3];    /* diagonal measurement variances, e-frame */
    double gp_r[INS_NR];        /* pre-fit residual of each row, last update */
    double gp_s[INS_NR];        /* pre-fit sqrt(H P H' + R) of each row      */
    double gate_nr;             /* normalized innovation that tripped the gate */
    int gate_row;               /* row index of that tripping measurement      */
    int nr;                     /* current measurement dimension */
    int init;                   /* covariance initialised */
} eskf_t;

typedef struct {
    /* imu data */
    char   imufile[1024];       /* ascii IMU file */
    double ts;                  /* imu sample interval (s) */
    int    inc;                 /* 1: increments (dtheta,dvel), 0: rates (omega,f) */
    int    nsub;                /* sub-samples per mechanization cycle */
    int    frame;               /* INS_FRAME_*: body frame the file is written in */
    double inm[3];              /* mounting angle pitch, roll, yaw (deg): vehicle -> body */
    double grsc, arsc;          /* rate-mode scale: deg/s and m/s^2 per data unit */
    /* noise and initial uncertainty, MSF units */
    double misd[3];             /* initial misalignment std (deg) */
    double vstd[3];             /* initial velocity std (m/s) */
    double pstd[3];             /* initial position std (m) */
    double gbstd[3];            /* initial gyro bias std (deg/h) */
    double abstd[3];            /* initial accel bias std (mg) */
    double mispsd[3];           /* misalignment psd (deg/sqrt(h)) */
    double vrwpsd[3];           /* velocity random walk psd (mg/sqrt(Hz)) */
    double pospsd[3];           /* position psd (m/sqrt(h)) */
    double gbpsd[3];            /* gyro bias psd (deg/h/sqrt(h)) */
    double abpsd[3];            /* accel bias psd (mg/sqrt(h)) */
    double gbcorr[3];           /* gyro bias correlation time (h) */
    double abcorr[3];           /* accel bias correlation time (h) */
    double Kg0[3], Ka0[3];      /* initial scale factors */
    double eb0[3];              /* initial gyro bias (deg/h) */
    double db0[3];              /* initial accel bias (mg) */
    /* alignment */
    int    align;               /* INS_ALIGN_* */
    double att0[3];             /* initial attitude pitch, roll, yaw (deg) */
    double blmin;               /* min baseline for pva (m) */
    double vtmin;               /* min speed for vva (m/s) */
    /* loosely coupled */
    double lever[3];            /* GNSS antenna lever arm, body frame (m) */
    double rpos[3];             /* measurement position std (m); all 0 = take sol_t::qr */
    double rvel[3];             /* measurement velocity std (m/s); all 0 = take sol_t::qv */
    int    meas;                /* INS_MEAS_* */
    int    cps;                 /* 0 cone, 1 no compensation, 2 polynomial */
    double maxdt;               /* |t_ins - t_gnss| beyond this the epoch is dropped (s) */
    double maxnr;               /* normalized innovation threshold, 0 = off */
    double infl;                /* process noise inflation */
    int    fbsign;              /* -1 = self-consistent feedback, +1 = MSF gins.cpp:691-695 verbatim */
    int    dbg;                 /* 1 = print the innovation, gain and feedback of every update */
    /* output */
    char   outfile[1024];       /* .ins output */
    double outf;                /* output rate (Hz) */
} insopt_t;

/* ---- public API -------------------------------------------------------- */
EXPORT void ins_defaults(insopt_t *opt);
EXPORT int  ins_readconf(const char *file, insopt_t *opt);
EXPORT int  ins_start(const insopt_t *opt);
EXPORT int  ins_feed(const sol_t *sol, const double *rb);  /* IMU up to sol->time, then LCI update */
EXPORT void ins_stop(void);
EXPORT int  ins_enabled(void);
EXPORT const char *ins_errmsg(void);

/* ---- exposed for the offline harness and for unit checks --------------- */
extern sins_t   ins_sins;
extern eskf_t   ins_kf;
extern insopt_t ins_opt;

/* math (insmath.c) */
void m3_eye(double A[9]);
void m3_mul(const double A[9], const double B[9], double C[9]);
void m3_trs(const double A[9], double At[9]);
void m3_mv(const double A[9], const double v[3], double w[3]);   /* w = A v   */
void m3_tv(const double A[9], const double v[3], double w[3]);   /* w = A' v  */
void m3_skew(const double v[3], double S[9]);                    /* [v x]     */
int  mm(int m, int n, int p, const double *A, const double *B, double *C);
int  mtv(int m, int n, const double *A, const double *v, double *w);
void earth_update(earth_t *e, const double pos[3], const double vn[3]);
void earth_v2dp(const earth_t *e, const double vn[3], double ts, double dp[3]);
void a2m(const double att[3], double C[9]);                      /* att = pitch,roll,yaw */
void m2a(const double C[9], double att[3]);
void q_mul(const double a[4], const double b[4], double c[4]);
void q_conj(const double q[4], double c[4]);
void q_norm(double q[4]);
void q_mv(const double q[4], const double v[3], double w[3]);
void rv2q(const double rv[3], double q[4]);
void q2m(const double q[4], double C[9]);
void m2q(const double C[9], double q[4]);
double v_norm(const double v[3]);

/* mechanization (sins.c) */
typedef struct {                /* one-plus-previous state for the polynomial mode */
    double wm1[3], vm1[3];
    int pf1, pf2;
} imust_t;

void imu_comp(const double wm[][3], const double vm[][3], int n, int cps,
              double phim[3], double dvbm[3], imust_t *st);
void sins_init(sins_t *s, const insopt_t *opt, const double pos[3], const double vn[3],
               const double att[3]);
void sins_update(sins_t *s, const double phim[3], const double dvbm[3], double ts);
void sins_align_att(sins_t *s, const double att[3]);

/* filter (eskf.c) */
void eskf_initP(eskf_t *kf, const insopt_t *opt, const sins_t *s);
void eskf_setQt(eskf_t *kf, const insopt_t *opt);
void eskf_time_update(eskf_t *kf, sins_t *s, double dt, double infl);
int  eskf_meas_update(eskf_t *kf, sins_t *s, const insopt_t *opt,
                      const double meas_pos[3], const double meas_vel[3],
                      int has_pos, int has_vel, double tmeas);
void eskf_feedback(eskf_t *kf, sins_t *s, int sign_att);

/* imu file (imuread.c) */
int  imu_open(FILE **fp, const char *file);
int  imu_next(FILE *fp, const insopt_t *opt, double *t, double w[3], double v[3]);

#ifdef __cplusplus
}
#endif

#endif /* RTKINS_H */
