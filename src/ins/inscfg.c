/*
 * inscfg.c - defaults and "key = value" configuration for the INS module.
 *
 * Defaults are taken from GREAT-MSF t_gsetign::t_gsetign() (gsetign.cpp:159-198)
 * so that a run with no configuration file reproduces MSF's own fallback:
 *   initial_misalignment_std 1 1 10 deg      misalignment_psd  1 deg/sqrt(h)
 *   initial_vel_std          0.1 m/s         vel_psd           0.1 mg/sqrt(Hz)
 *   initial_pos_std          3 m             pos_psd          (MSF leaves it 0)
 *   initial_gyro_std         10 deg/h        gyro_psd          (MSF leaves it 0)
 *   initial_acce_std         5 mg            acce_psd          (MSF leaves it 0)
 *   lever 0 0 0   delay_t 0   TS 1.0
 *
 * DIFF: MSF declares the bias Markov correlation time in t_gsins (_tauG/_tauA,
 *       gins.h:131) but nothing in the whole tree ever assigns it - both stay
 *       zero, so the bias model in MSF is a pure random walk (gins.cpp:544,546).
 *       gyro_corr_time / acce_corr_time are provided here for the Gauss-Markov
 *       form used by PSINS and by I:\OB_GINS-main (corrtime: 1.0 h); the default
 *       is 0, i.e. MSF behaviour.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "rtkins.h"

void ins_defaults(insopt_t *opt)
{
    int i;

    memset(opt,0,sizeof(insopt_t));
    opt->ts      = 0.005;         /* 200 Hz, I:\OB_GINS-main\dataset\ADIS16465.txt */
    opt->inc     = 1;             /* that file carries dtheta/dvel increments     */
    opt->nsub    = 1;
    opt->grsc    = D2R;
    opt->arsc    = 1.0;
    opt->cps     = 0;             /* coning/sculling compensation on              */

    for (i=0;i<3;i++) {
        opt->misd[i]   = i<2?1.0:10.0;
        opt->vstd[i]   = 0.1;
        opt->pstd[i]   = 3.0;
        opt->gbstd[i]  = 10.0;
        opt->abstd[i]  = 5.0;
        opt->mispsd[i] = 1.0;
        opt->vrwpsd[i] = 0.1;
        opt->pospsd[i] = 0.0;
        opt->gbpsd[i]  = 0.0;
        opt->abpsd[i]  = 0.0;
        opt->gbcorr[i]= 0.0;
        opt->abcorr[i]= 0.0;
        opt->Kg0[i]    = 1.0;
        opt->Ka0[i]    = 1.0;
        opt->eb0[i]    = 0.0;
        opt->db0[i]    = 0.0;
        opt->lever[i]  = 0.0;
        opt->rpos[i]   = 0.0;
        opt->rvel[i]   = 0.0;
    }
    opt->align  = INS_ALIGN_AUTO;
    opt->blmin  = 5.0;
    opt->vtmin  = 2.0;
    opt->meas   = INS_MEAS_AUTO;
    opt->maxdt  = 1.0;            /* MSF: |sins.t - tmeas| > delay_t drops the epoch */
    opt->maxnr  = 3.0;            /* MSF outlier_detect(): threshold 3 when unset */
    opt->infl   = 1.0;
    opt->fbsign = -1;           /* see the header comment of eskf.c */
    opt->dbg    = 0;
    opt->outf   = 0.0;            /* 0 = every GNSS epoch + every mechanization epoch */
}

static int iskey(const char *key, const char *name)
{
    return strcmp(key,name)==0;
}

/* keyword values are case insensitive, without pulling in a platform extension */
static int isword(const char *s, const char *w)
{
    int i;
    for (i=0;s[i]&&w[i];i++) {
        char a=s[i],b=w[i];
        if (a>='A'&&a<='Z') a+=(char)('a'-'A');
        if (b>='A'&&b<='Z') b+=(char)('a'-'A');
        if (a!=b) return 0;
    }
    return s[i]=='\0'&&w[i]=='\0';
}
static int getdbl(const char *v, double *d, int n)
{
    int i,c=0;
    const char *p=v;
    char *q;
    for (i=0;i<n&&*p;i++) {
        while (*p==' '||*p=='\t'||*p==','||*p==';') p++;
        if (!*p) break;
        d[c]=strtod(p,&q);
        if (q==p) break;
        c++; p=q;
    }
    return c;
}
static int getint(const char *v,int *i)
{
    double d;
    if (getdbl(v,&d,1)<=0) return 0;
    *i=(int)(d<0?d-0.5:d+0.5);
    return 1;
}
static void getstr(const char *v,char *s,int n)
{
    const char *p=v;
    int i=0;
    while (*p==' '||*p=='\t') p++;
    while (i<n-1&&*p&&*p!='\r'&&*p!='\n') {
        if (*p=='"'||*p=='\'') {p++; continue;}
        s[i++]=*p++;
    }
    while (i>0&&(s[i-1]==' '||s[i-1]=='\t')) i--;
    s[i]='\0';
}
static int getkw(const char *v,const char *w,int *dst)
{
    const char *p=v;
    char s[32];
    int i=0;
    /* first token only, so "meas-type = pos   # comment" parses (getstr keeps
     * everything up to end-of-line because file paths may hold spaces) */
    while (*p==' '||*p=='\t') p++;
    while (i<31&&*p&&*p!=' '&&*p!='\t'&&*p!=','&&*p!=';'&&*p!='#'&&
           *p!='\r'&&*p!='\n') s[i++]=*p++;
    s[i]='\0';
    if (!isword(s,w)) return 0;
    *dst=1;
    return 1;
}
/* one value replicated to three axes, or one value per axis */
static int getstd3(const char *v,double d[3])
{
    int c=getdbl(v,d,3);
    if (c<1) return 0;
    if (c==1) d[1]=d[2]=d[0];
    else if (c==2) d[2]=d[1];
    return c;
}

/*
 * return: 1 ok, 0 unknown key (ignored on purpose so one file can carry the
 * options of other programs), -1 malformed value
 */
static int parse_line(const char *line, insopt_t *opt)
{
    char key[64];
    const char *p=line,*e;
    char *v;
    int n=0;

    while (*p==' '||*p=='\t') p++;
    if (*p=='#'||*p=='/'||*p=='\r'||*p=='\n'||!*p) return 1;
    if (!(e=strchr(p,'='))) return 1;
    for (;p<e&&n<63;p++) if (*p!=' '&&*p!='\t') key[n++]=*p;
    key[n]='\0';
    v=(char *)(e+1);

    /* ---- data ---------------------------------------------------------- */
    if      (iskey(key,"imu-file")||iskey(key,"imufile"))  getstr(v,opt->imufile,sizeof(opt->imufile));
    else if (iskey(key,"imu-ts")||iskey(key,"imudatarate")) {
        double d; if (getdbl(v,&d,1)<1) return -1;
        opt->ts = iskey(key,"imudatarate")?(d>0.0?1.0/d:0.005):d;
    }
    else if (iskey(key,"imu-inc"))    { if (!getint(v,&opt->inc)) return -1; }
    else if (iskey(key,"imu-nsub"))   { if (!getint(v,&opt->nsub)) return -1; }
    else if (iskey(key,"imu-mount"))  { if (getdbl(v,opt->inm,3)<1) return -1; }
    else if (iskey(key,"imu-frame")) {
        static const char *fn[6]={"flu","rfu","frd","rbd","lbu","bru"};
        int f,ok=0;
        double d;
        for (f=0;f<6;f++) {
            if (getkw(v,fn[f],&ok)) { opt->frame=f; break; }
        }
        if (f>=6&&getdbl(v,&d,1)==1&&d>=0&&d<6) opt->frame=(int)d;
        else if (f>=6) { fprintf(stderr,"ins config: unknown imu-frame value: %s",v); return -1; }
    }
    else if (iskey(key,"imu-grsc"))   { if (getdbl(v,&opt->grsc,1)<1) return -1; }
    else if (iskey(key,"imu-arsc"))   { if (getdbl(v,&opt->arsc,1)<1) return -1; }
    else if (iskey(key,"cone-comp"))  { if (!getint(v,&opt->cps)) return -1; }
    /* ---- initial uncertainty (MSF attribute names) ---------------------
     * one value means "the same for all three axes" (getstd3), which is how
     * these keys are written in MSF's own ini files */
    else if (iskey(key,"initial_misalignment_std")) { if (getstd3(v,opt->misd)<1) return -1; }
    else if (iskey(key,"initial_vel_std"))          { if (getstd3(v,opt->vstd)<1) return -1; }
    else if (iskey(key,"initial_pos_std"))          { if (getstd3(v,opt->pstd)<1) return -1; }
    else if (iskey(key,"initial_gyro_std"))         { if (getstd3(v,opt->gbstd)<1) return -1; }
    else if (iskey(key,"initial_acce_std"))         { if (getstd3(v,opt->abstd)<1) return -1; }
    /* ---- process noise psd (MSF attribute names) ----------------------- */
    else if (iskey(key,"misalignment_psd"))         { if (getstd3(v,opt->mispsd)<1) return -1; }
    else if (iskey(key,"vel_psd"))                  { if (getstd3(v,opt->vrwpsd)<1) return -1; }
    else if (iskey(key,"pos_psd"))                  { if (getstd3(v,opt->pospsd)<1) return -1; }
    else if (iskey(key,"gyro_psd"))                 { if (getstd3(v,opt->gbpsd)<1) return -1; }
    else if (iskey(key,"acce_psd"))                 { if (getstd3(v,opt->abpsd)<1) return -1; }
    else if (iskey(key,"gyro_corr_time"))           { if (getstd3(v,opt->gbcorr)<1) return -1; }
    else if (iskey(key,"acce_corr_time"))           { if (getstd3(v,opt->abcorr)<1) return -1; }
    /* ---- fixed IMU errors --------------------------------------------- */
    else if (iskey(key,"initial_gyro_bias"))        { if (getstd3(v,opt->eb0)<1) return -1; }
    else if (iskey(key,"initial_acce_bias"))        { if (getstd3(v,opt->db0)<1) return -1; }
    else if (iskey(key,"initial_gyro_scale"))       { if (getstd3(v,opt->Kg0)<1) return -1; }
    else if (iskey(key,"initial_acce_scale"))       { if (getstd3(v,opt->Ka0)<1) return -1; }
    /* ---- alignment ----------------------------------------------------- */
    else if (iskey(key,"align")) {
        int a=0; double d;
        if      (getkw(v,"auto",&a)) opt->align=INS_ALIGN_AUTO;
        else if (getkw(v,"att",&a))  opt->align=INS_ALIGN_ATT;
        else if (getkw(v,"pva",&a))  opt->align=INS_ALIGN_PVA;
        else if (getkw(v,"vva",&a))  opt->align=INS_ALIGN_VVA;
        else if (getdbl(v,&d,1)>=1)  opt->align=(int)d;
        else return -1;
    }
    else if (iskey(key,"attitude")||iskey(key,"initatt")) {
        if (getdbl(v,opt->att0,3)<3) return -1;
    }
    else if (iskey(key,"min-baseline"))             { if (getdbl(v,&opt->blmin,1)<1) return -1; }
    else if (iskey(key,"min-speed"))                { if (getdbl(v,&opt->vtmin,1)<1) return -1; }
    /* ---- loosely coupled ---------------------------------------------- */
    else if (iskey(key,"lever")||iskey(key,"antlever")) {
        if (getdbl(v,opt->lever,3)<1) return -1;
    }
    else if (iskey(key,"meas-pos-std"))         { if (getstd3(v,opt->rpos)<1) return -1; }
    else if (iskey(key,"meas-vel-std"))         { if (getstd3(v,opt->rvel)<1) return -1; }
    else if (iskey(key,"meas-type")) {
        int m=0; double d;
        if      (getkw(v,"auto",&m))    opt->meas=INS_MEAS_AUTO;
        else if (getkw(v,"pos",&m))     opt->meas=INS_MEAS_POS;
        else if (getkw(v,"vel",&m))     opt->meas=INS_MEAS_VEL;
        else if (getkw(v,"pos-vel",&m)) opt->meas=INS_MEAS_POS_VEL;
        else if (getdbl(v,&d,1)>=1)     opt->meas=(int)d;
        else return -1;
    }
    else if (iskey(key,"max-dt"))                   { if (getdbl(v,&opt->maxdt,1)<1) return -1; }
    else if (iskey(key,"max-nr"))                   { if (getdbl(v,&opt->maxnr,1)<1) return -1; }
    else if (iskey(key,"debug"))                    { if (!getint(v,&opt->dbg)) return -1; }
    else if (iskey(key,"inflation"))                { if (getdbl(v,&opt->infl,1)<1) return -1; }
    else if (iskey(key,"fb-sign"))                  { if (!getint(v,&opt->fbsign)) return -1;
                                                      opt->fbsign=opt->fbsign>=0?1:-1; }
    /* ---- output -------------------------------------------------------- */
    else if (iskey(key,"out-ins"))                  getstr(v,opt->outfile,sizeof(opt->outfile));
    else if (iskey(key,"out-rate"))                 { if (getdbl(v,&opt->outf,1)<1) return -1; }
    return 1;
}

int ins_readconf(const char *file, insopt_t *opt)
{
    char line[1024];
    FILE *fp;
    int i=0;

    if (!file||!*file) return 0;
    if (!(fp=fopen(file,"r"))) return -1;
    while (fgets(line,sizeof(line),fp)) {
        i++;
        if (parse_line(line,opt)<0) {
            fprintf(stderr,"ins config: malformed value at line %d: %s",i,line);
        }
    }
    fclose(fp);
    return 1;
}
