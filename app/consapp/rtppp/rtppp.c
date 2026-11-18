/* ---------------------------------------------------------------------------
 * Name: rtppp.c
 * Function: Windows console shell for real-time PPP over Ntrip streams.
 *           rtkrcv.c cannot be built on Windows (unistd.h, curses vt.c), but
 *           rtksvr.c and stream.c can (rtksvr.c uses CreateThread under WIN32).
 *           This shell does what rtkrcv's startsvr() does, without the POSIX
 *           terminal.
 * Input:  -k conf-file [-c cycle-ms] [-t duration-s] [-x trace-level]
 * Output: solution to conf outstr1 stream + 1 Hz status to stdout
 * -------------------------------------------------------------------------*/
#include "rtklib.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define PROGNM "rtppp"

/* stream options table, copied from app/consapp/rtkrcv/rtkrcv.c -------------*/
#define ISTOPT  "0:off,1:serial,2:file,3:tcpsvr,4:tcpcli,6:ntripcli,7:ftp,8:http"
#define OSTOPT  "0:off,1:serial,2:file,3:tcpsvr,4:tcpcli,5:ntripsvr,9:ntripcas,11:udpcli"
#define FMTOPT  "0:rtcm2,1:rtcm3,2:oem4,4:ubx,5:swift,6:hemis,7:skytraq,8:javad,9:nvs,10:binex,11:rt17,12:sbf,14:unicore,15:rinex,16:sp3,17:clk"
#define NMEOPT  "0:off,1:latlon,2:single"
#define SOLOPT  "0:llh,1:xyz,2:enu,3:nmea,4:stat"

static int   strtype[]={STR_NONE,STR_NONE,STR_NONE,STR_NONE,
                        STR_NONE,STR_NONE,STR_NONE,STR_NONE};
static char  strpath[8][MAXSTRPATH]={"","","","","","","",""};
static int   strfmt[]={STRFMT_RTCM3,STRFMT_RTCM3,STRFMT_RTCM3,SOLF_LLH,SOLF_LLH};
static char  rcvopt[3][256]={"","",""};
static int   navmsgsel=0;
static int   timeout=10000,reconnect=10000,buffsize=32768,fswapmargin=30;

static opt_t rcvopts[]={
    {"inpstr1-type",   3,(void *)&strtype[0],ISTOPT},
    {"inpstr2-type",   3,(void *)&strtype[1],ISTOPT},
    {"inpstr3-type",   3,(void *)&strtype[2],ISTOPT},
    {"inpstr1-path",   2,(void *)strpath [0],""},
    {"inpstr2-path",   2,(void *)strpath [1],""},
    {"inpstr3-path",   2,(void *)strpath [2],""},
    {"inpstr1-format", 3,(void *)&strfmt [0],FMTOPT},
    {"inpstr2-format", 3,(void *)&strfmt [1],FMTOPT},
    {"inpstr3-format", 3,(void *)&strfmt [2],FMTOPT},
    {"inpstr1-rcvopt", 2,(void *)rcvopt[0],""},
    {"inpstr2-rcvopt", 2,(void *)rcvopt[1],""},
    {"inpstr3-rcvopt", 2,(void *)rcvopt[2],""},
    {"outstr1-type",   3,(void *)&strtype[3],OSTOPT},
    {"outstr2-type",   3,(void *)&strtype[4],OSTOPT},
    {"outstr1-path",   2,(void *)strpath [3],""},
    {"outstr2-path",   2,(void *)strpath [4],""},
    {"outstr1-format", 3,(void *)&strfmt [3],SOLOPT},
    {"outstr2-format", 3,(void *)&strfmt [4],SOLOPT},
    {"logstr1-type",   3,(void *)&strtype[5],OSTOPT},
    {"logstr2-type",   3,(void *)&strtype[6],OSTOPT},
    {"logstr3-type",   3,(void *)&strtype[7],OSTOPT},
    {"logstr1-path",   2,(void *)strpath [5],""},
    {"logstr2-path",   2,(void *)strpath [6],""},
    {"logstr3-path",   2,(void *)strpath [7],""},
    {"misc-navmsgsel", 3,(void *)&navmsgsel,"0:all,1:rover,2:base,3:corr"},
    {"",0,NULL,""}
};

static rtksvr_t svr;
static volatile int stop_=0;

/* application defined functions (as in rnx2rtkp.c) --------------------------*/
extern int showmsg(const char *format,...)
{
    va_list arg;
    va_start(arg,format); vfprintf(stderr,format,arg); va_end(arg);
    fprintf(stderr,"\r");
    return 0;
}
extern void settspan(gtime_t ts,gtime_t te) {(void)ts;(void)te;}
extern void settime(gtime_t time) {(void)time;}

#ifdef WIN32
static BOOL WINAPI ctrlc(DWORD sig) {(void)sig; stop_=1; return TRUE;}
#endif

/* load antenna files (rtkrcv.c readant()) -----------------------------------*/
static void readant(prcopt_t *opt, nav_t *nav, pcvs_t *pcvsr,
                    const filopt_t *fopt)
{
    const pcv_t pcv0={0};
    gtime_t time=timeget();
    int i;
    
    opt->pcvr[0]=opt->pcvr[1]=pcv0;
    
    if (*fopt->rcvantp) {
        if (readpcv(fopt->rcvantp,pcvsr)) {
            for (i=0;i<2;i++) {
                pcv_t *pcv;
                if (!*opt->anttype[i]||!strcmp(opt->anttype[i],"*")) continue;
                if (!(pcv=searchpcv(0,opt->anttype[i],time,pcvsr))) {
                    fprintf(stderr,"no antenna %s in %s\n",opt->anttype[i],
                            fopt->rcvantp);
                    continue;
                }
                opt->pcvr[i]=*pcv;
            }
        } else fprintf(stderr,"rcvant open error: %s\n",fopt->rcvantp);
    }
    if (*fopt->satantp) {
        pcvs_t pcvs={0};
        if (readpcv(fopt->satantp,&pcvs)) {
            for (i=0;i<MAXSAT;i++) {
                pcv_t *pcv=searchpcv(i+1,"",time,&pcvs);
                if (!pcv) continue;
                nav->pcvs[i]=*pcv;
            }
            free_pcvs(&pcvs);
            fprintf(stderr,"sat antenna: %s\n",fopt->satantp);
        } else fprintf(stderr,"satant open error: %s\n",fopt->satantp);
    }
}
/* satellites with consistent ssr orbit+clock corrections ---------------------*/
static int count_ssr(const nav_t *nav)
{
    int i,n=0;
    for (i=0;i<MAXSAT;i++) {
        if (nav->ssr[i].t0[0].time&&nav->ssr[i].t0[1].time&&
            nav->ssr[i].iod[0]==nav->ssr[i].iod[1]) n++;
    }
    return n;
}
int main(int argc, char **argv)
{
    prcopt_t prcopt;
    solopt_t sopt[2]={{0},{0}};
    filopt_t filopt={{0}};
    double npos[3]={0,0,0};
    const char *cmds[]={NULL,NULL,NULL};
    const char *cmds_periodic[]={NULL,NULL,NULL};
    const char *ropts[]={rcvopt[0],rcvopt[1],rcvopt[2]};
    char errmsg[2048]="",tstr[64],*conffile=NULL,tracefile[32]="rtppp.trace";
    int stropt[8]={0};
    int i,j,cycle=10,dur=60,tlev=0;
    const char *paths[8]={strpath[0],strpath[1],strpath[2],strpath[3],
                        strpath[4],strpath[5],strpath[6],strpath[7]};
    
    for (i=1;i<argc;i++) {
        if      (!strcmp(argv[i],"-k")&&i+1<argc) conffile=argv[++i];
        else if (!strcmp(argv[i],"-c")&&i+1<argc) cycle=atoi(argv[++i]);
        else if (!strcmp(argv[i],"-t")&&i+1<argc) dur=atoi(argv[++i]);
        else if (!strcmp(argv[i],"-x")&&i+1<argc) tlev=atoi(argv[++i]);
        else if (!strcmp(argv[i],"-f")&&i+1<argc) {
            strcpy(tracefile,argv[++i]);
        } else {
            fprintf(stderr,
                "usage: %s -k conf [-c cycle-ms] [-t sec] [-x trace] [-f file]\n",
                PROGNM);
            return EXIT_FAILURE;
        }
    }
    if (!conffile) {
        fprintf(stderr,"error: no conf file (-k)\n");
        return EXIT_FAILURE;
    }
    /* processing/solution/file options: same table as rnx2rtkp */
    resetsysopts();
    if (!loadopts(conffile,sysopts)) return EXIT_FAILURE;
    getsysopts(&prcopt,&sopt[0],&filopt);
    sopt[1]=sopt[0];
    
    /* stream options: rtkrcv's rcvopts table (sysopts has no inpstr keys) */
    if (!loadopts(conffile,rcvopts)) return EXIT_FAILURE;
    
    fprintf(stderr,"%s RTKLIB %s %s\n",PROGNM,VER_RTKLIB,PATCH_LEVEL);
    fprintf(stderr,"mode=%d sateph=%d ionoopt=%d tropopt=%d navsys=%d\n",
            prcopt.mode,prcopt.sateph,prcopt.ionoopt,prcopt.tropopt,
            prcopt.navsys);
    for (i=0;i<3;i++) {
        fprintf(stderr,"  inpstr%d: type=%d fmt=%d path=%.60s\n",i+1,strtype[i],
                strfmt[i],strpath[i]);
    }
    fprintf(stderr,"  outstr1: type=%d path=%s\n",strtype[3],strpath[3]);
    if (prcopt.mode<PMODE_PPP_KINEMA) {
        fprintf(stderr,"  WARNING: pos1-posmode is not a PPP mode\n");
    }
    if (tlev>0) {
        traceopen(tracefile);
        tracelevel(tlev);
    }
    stropt[0]=timeout; stropt[1]=reconnect; stropt[2]=1000;
    stropt[3]=buffsize; stropt[4]=fswapmargin;
    strsetopt(stropt);
    
    if (!rtksvrinit(&svr)) {
        fprintf(stderr,"rtksvrinit error\n");
        return EXIT_FAILURE;
    }
    /* before rtksvrstart: rtkinit() copies prcopt, update_ssr() fills nav.ssr */
    readant(&prcopt,&svr.nav,&svr.pcvsr,&filopt);
    
    if (!rtksvrstart(&svr,cycle,buffsize,strtype,paths,strfmt,
                     navmsgsel,(const char **)cmds,(const char **)cmds_periodic,
                     ropts,0,0,npos,&prcopt,sopt,NULL,errmsg)) {
        fprintf(stderr,"rtk server start error: %s\n",errmsg);
        rtksvrfree(&svr);
        if (tlev>0) traceclose();
        return EXIT_FAILURE;
    }
#ifdef WIN32
    SetConsoleCtrlHandler(ctrlc,TRUE);
#endif
    printf("server started, run %d s\n",dur);
    printf("%-20s %3s %3s %11s %12s %9s %7s %7s %7s %5s %5s\n","time","Q","ns",
           "lat(deg)","lon(deg)","hgt(m)","sdN(m)","sdE(m)","sdU(m)","eph","ssr");
    
    for (i=0;i<dur&&!stop_;i++) {
        double pos[3];
        sol_t sol;
        int ne,nss,st[3];
        
        sleepms(1000);
        rtksvrlock(&svr);
        sol=svr.rtk.sol;
        ne=svr.nav.n;
        nss=count_ssr(&svr.nav);
        for (j=0;j<3;j++) st[j]=svr.stream[j].state;
        rtksvrunlock(&svr);
        
        ecef2pos(sol.rr,pos);
        time2str(gpst2utc(sol.time),tstr,0); /* utc, to match the .pos file */
        printf("%-20s %3d %3d %11.7f %12.7f %9.3f %7.3f %7.3f %7.3f %5d %5d",
               tstr,sol.stat,sol.ns,pos[0]*R2D,pos[1]*R2D,pos[2],
               sol.qr[0],sol.qr[1],sol.qr[2],ne,nss);
        printf("  str=%d%d%d\n",st[0],st[1],st[2]);
        fflush(stdout);
    }
    printf("\ninput stream message counts (obs,eph,antpos,ssr):\n");
    for (i=0;i<3;i++) {
        static const char *nm[]={"rov ","base","corr"};
        printf("  %s: %d %d %d %d\n",nm[i],svr.nmsg[i][0],svr.nmsg[i][1],
               svr.nmsg[i][4],svr.nmsg[i][7]);
    }
    rtksvrstop(&svr,(const char **)cmds);
    rtksvrfree(&svr);
    if (tlev>0) traceclose();
    printf("stopped\n");
    return EXIT_SUCCESS;
}
