/*
 * imuread.c - ASCII IMU file reader.
 *
 * Line format (whitespace separated, '#' starts a comment):
 *   t  x1 x2 x3 x4 x5 x6 [x7 [x8]]
 * t is GPS time of week in seconds.  x1..x3 are the gyro channel, x4..x6 the
 * accel channel.  x7/x8 are the OB_GINS odometry channels and are ignored here.
 *
 * The same reader covers the two conventions in use:
 *   inc=1  increments          x1..x3 = dtheta (rad),      x4..x6 = dvel (m/s)
 *          as delivered by ADIS16465.txt / I:\OB_GINS-main\dataset\ob_gins.yaml
 *   inc=0  rates               x1..x3 = omega,             x4..x6 = specific force
 *          scaled to SI by grsc / arsc (units per data unit)
 * imu_next() always returns increments, so the mechanization only ever sees
 * dtheta/dvel, exactly like t_gimu::Update(wm,vm).
 *
 * imu-frame (opt->frame) finally rotates both channels from the body frame the
 * file happens to be written in into this module's FLU, the same job MSF's
 * AxisOrder table does in imufile.cpp:105-133.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "rtkins.h"

#define MAXCOL 9

/* body frame of the file -> our canonical FLU, see the INS_FRAME_* block of
 * rtkins.h.  Row-major, index = INS_FRAME_*. */
static const double frame_rot[6][9] = {
    {1,0,0, 0,1,0, 0,0,1},     /* flu */
    {0,1,0,-1,0,0, 0,0,1},     /* rfu */
    {1,0,0, 0,-1,0, 0,0,-1},   /* frd */
    {0,-1,0,-1,0,0, 0,0,-1},   /* rbd */
    {0,-1,0, 1,0,0, 0,0,1},    /* lbu */
    {-1,0,0, 0,-1,0, 0,0,1}    /* bru */
};

static void frame_mul(int f, const double v[3], double o[3])
{
    const double *C=frame_rot[f];
    int i;
    for (i=0;i<3;i++) o[i]=C[3*i]*v[0]+C[3*i+1]*v[1]+C[3*i+2]*v[2];
}

int imu_open(FILE **fp, const char *file)
{
    if (!file||!*file) return 0;
    if (!(*fp=fopen(file,"r"))) return -1;
    return 1;
}

/*
 * return: 1 data read, 0 end of file, -1 malformed line
 * *t   = seconds of week
 * *dth = delta angle   (rad)
 * *dvl = delta velocity(m/s)
 */
int imu_next(FILE *fp, const insopt_t *opt, double *t, double dth[3], double dvl[3])
{
    char line[1024], *p;
    double d[MAXCOL];
    int i,n=0;

    if (!fp) return 0;

    for (;;) {
        if (!fgets(line,sizeof(line),fp)) return 0;
        for (p=line; *p==' '||*p=='\t'; p++) ;
        if (*p=='#'||*p=='\r'||*p=='\n'||*p=='\0') continue;
        break;
    }
    for (p=line; n<MAXCOL; ) {
        char *q;
        double v=strtod(p,&q);
        if (q==p) break;
        d[n++]=v;
        p=q;
        while (*p==' '||*p=='\t'||*p==',') p++;
    }
    if (n<7) return -1;

    *t=d[0];
    for (i=0;i<3;i++) {
        dth[i]=d[1+i];
        dvl[i]=d[4+i];
    }
    if (!opt->inc) {                  /* rates -> increments over one sample */
        for (i=0;i<3;i++) {
            dth[i]*=opt->grsc*opt->ts;
            dvl[i]*=opt->arsc*opt->ts;
        }
    }
    if (opt->frame!=INS_FRAME_FLU) {  /* file body frame -> FLU */
        double tmp[3];
        frame_mul(opt->frame,dth,tmp);  memcpy(dth,tmp,sizeof(double)*3);
        frame_mul(opt->frame,dvl,tmp);  memcpy(dvl,tmp,sizeof(double)*3);
    }
    return 1;
}
