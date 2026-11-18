/* ---------------------------------------------------------------------------
 * Name: rtkrtk.c
 * Function: Windows console shell for real-time RTK (relative positioning over
 *           streams), the RTK twin of rtppp.c and the real-time counterpart of
 *           postrtk.c. It exists for the same reason rtppp.c does: rtkrcv.c
 *           cannot be built on Windows (unistd.h, curses vt.c), but rtksvr.c
 *           and stream.c can (rtksvr.c uses CreateThread under WIN32).
 * Input:  -k conf-file [-c cycle-ms] [-t duration-s] [-x trace-level] [-f file]
 * Output: solution to the conf outstr1 stream + 1 Hz status to stdout
 * Notes:  pos1-posmode = dgps(1) / kinematic(2) / static(3) / static-start(4)
 *         / moving-base(5); inpstr1 = rover, inpstr2 = base, inpstr3 = SSR.
 * -------------------------------------------------------------------------*/
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "rtklib.h"

#define PROGNM "rtkrtk"

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

static const char *modename[]={
    "single","dgps","kinematic","static","static-start","movingbase","fixed",
    "ppp-kine","ppp-static","ppp-fixed"};
static const char *postypname[]={
    "llh","xyz","single","posfile","rinexhead","rtcm"};
static const char *armodename[]={
    "off","continuous","instantaneous","fix-and-hold"};
static const char *ephoptname[]={
    "brdc","precise","sbas","ssr-apc","ssr-com"};
static const char *solqname[]={
    "none","fix","float","sbas","dgps","single","ppp","dr"};
static const char *navselname[]={"all","rover","base","corr"};
static const char *strtpname[]={
    "off","serial","file","tcpsvr","tcpcli","ntripsvr","ntripcli","ftp","http",
    "ntripcas","udpsvr","udpcli","membuf"};
static const char *strfmname[]={
    "rtcm2","rtcm3","oem4","?","ubx","swift","hemis","skytraq","javad","nvs",
    "binex","rt17","sbf","?","unicore","rinex","sp3","clk","sbas","nmea"};
static const char *solfpname[]={"llh","xyz","enu","nmea","stat"};

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
static const char *nm(const char *tab[],int i,int n)
{
    return 0<=i&&i<n?tab[i]:"?";
}
/* one-line stream description -----------------------------------------------*/
static void prnstr(const char *label,int type,int fmt,const char *path,
                   int fmtmax,const char *tab[])
{
    printf("  %-6s: type=%d(%s) fmt=%d(%s) path=%.60s\n",label,type,
           nm(strtpname,type,(int)(sizeof(strtpname)/sizeof(*strtpname))),fmt,
           fmt<0||fmt>=fmtmax?"?":tab[fmt],path);
}
int main(int argc,char **argv)
{
    prcopt_t prcopt;
    solopt_t sopt[2]={{0},{0}};
    filopt_t filopt={{0}};
    double npos[3]={0,0,0},rb0[3]={0,0,0};
    const char *cmds[]={NULL,NULL,NULL};
    const char *cmds_periodic[]={NULL,NULL,NULL};
    const char *ropts[]={rcvopt[0],rcvopt[1],rcvopt[2]};
    char errmsg[2048]="",tstr[64],*conffile=NULL,tracefile[32]="rtkrtk.trace";
    char errbuf[MAXERRMSG];
    int stropt[8]={0};
    int i,j,cycle=10,dur=60,tlev=0,nq[8]={0},nep=0;
    const char *paths[8]={strpath[0],strpath[1],strpath[2],strpath[3],
                        strpath[4],strpath[5],strpath[6],strpath[7]};
    
    for (i=1;i<argc;i++) {
        if      (!strcmp(argv[i],"-k")&&i+1<argc) conffile=argv[++i];
        else if (!strcmp(argv[i],"-c")&&i+1<argc) cycle=atoi(argv[++i]);
        else if (!strcmp(argv[i],"-t")&&i+1<argc) dur=atoi(argv[++i]);
        else if (!strcmp(argv[i],"-x")&&i+1<argc) tlev=atoi(argv[++i]);
        else if (!strcmp(argv[i],"-f")&&i+1<argc) strcpy(tracefile,argv[++i]);
        else {
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
    /* provenance line in the .pos header written by outstr1 */
    sprintf(sopt[0].prog,"%s ver.%s %s",PROGNM,VER_RTKLIB,PATCH_LEVEL);
    sprintf(sopt[1].prog,"%s ver.%s %s",PROGNM,VER_RTKLIB,PATCH_LEVEL);
    
    /* stream options: rtkrcv's rcvopts table (sysopts has no inpstr keys) */
    if (!loadopts(conffile,rcvopts)) return EXIT_FAILURE;
    
    fprintf(stderr,"%s RTKLIB %s %s\n",PROGNM,VER_RTKLIB,PATCH_LEVEL);
    printf("mode=%d(%s) elmask=%.1f ionoopt=%d tropopt=%d navsys=%d nf=%d\n",
           prcopt.mode,nm(modename,prcopt.mode,10),prcopt.elmin*R2D,
           prcopt.ionoopt,prcopt.tropopt,prcopt.navsys,prcopt.nf);
    printf("      sateph=%d(%s) dynamics=%d tidecorr=%d\n",prcopt.sateph,
           nm(ephoptname,prcopt.sateph,5),prcopt.dynamics,prcopt.tidecorr);
    printf("  AR  : mode=%d(%s) glo=%d bds=%d thresar=%.2f minlock=%d minfix=%d\n",
           prcopt.modear,nm(armodename,prcopt.modear,4),prcopt.glomodear,
           prcopt.bdsmodear,prcopt.thresar[0],prcopt.minlock,prcopt.minfix);
    prnstr("rover",strtype[0],strfmt[0],strpath[0],
           (int)(sizeof(strfmname)/sizeof(*strfmname)),strfmname);
    prnstr("base", strtype[1],strfmt[1],strpath[1],
           (int)(sizeof(strfmname)/sizeof(*strfmname)),strfmname);
    prnstr("corr", strtype[2],strfmt[2],strpath[2],
           (int)(sizeof(strfmname)/sizeof(*strfmname)),strfmname);
    prnstr("out1", strtype[3],strfmt[3],strpath[3],
           (int)(sizeof(solfpname)/sizeof(*solfpname)),solfpname);
    prnstr("out2", strtype[4],strfmt[4],strpath[4],
           (int)(sizeof(solfpname)/sizeof(*solfpname)),solfpname);
    printf("  log : %d %d %d (navmsgsel=%d(%s))\n",strtype[5],strtype[6],
           strtype[7],navmsgsel,nm(navselname,navmsgsel,4));
    printf("  base: refpos=%d(%s)",prcopt.refpos,nm(postypname,prcopt.refpos,6));
    if (prcopt.refpos==POSOPT_POS_LLH||prcopt.refpos==POSOPT_POS_XYZ) {
        printf(" pos=%.3f %.3f %.3f (ecef m)",prcopt.rb[0],prcopt.rb[1],
               prcopt.rb[2]);
    }
    printf("\n");
    printf("  ant : sat=%s rcv=%s\n",*filopt.satantp?filopt.satantp:"-",
           *filopt.rcvantp?filopt.rcvantp:"-");
    
    /* relative-mode preflight: warn, never abort, the conf owner decides ----*/
    if (strtype[0]==STR_NONE) {
        printf("ERROR: inpstr1 (rover) is off - there is no observation input, "
               "nothing to process.\n");
        return EXIT_FAILURE;
    }
    /* rcvraw.c input_raw() (rcvraw.c:1431-1443) has no STRFMT_RINEX branch, so
     * a RINEX file stream opens and reads bytes but never yields an epoch. */
    for (i=0;i<3;i++) {
        if (strtype[i]!=STR_NONE&&strfmt[i]==STRFMT_RINEX) {
            printf("ERROR: inpstr%d-format=rinex - rtksvr cannot decode it "
                   "(input_raw() in rcvraw.c has no RINEX branch, only convrnx.c "
                   "reads RINEX). Use the receiver's raw format or rtcm2/rtcm3.\n",
                   i+1);
            return EXIT_FAILURE;
        }
    }
    if (prcopt.mode>=PMODE_PPP_KINEMA) {
        printf("WARNING: pos1-posmode=%s is a PPP mode - use rtppp.\n",
               nm(modename,prcopt.mode,10));
    }
    else if (prcopt.mode<PMODE_DGPS) {
        printf("WARNING: pos1-posmode=%s is not a relative mode; inpstr2 (base) "
               "is ignored.\n",nm(modename,prcopt.mode,10));
    }
    else if (strtype[1]==STR_NONE&&strtype[2]==STR_NONE) {
        printf("ERROR: relative mode %s with both inpstr2 (base) and inpstr3 "
               "(corr) off. rtkpos() (rtkpos.c:2411) hits nr==0, prints \"no "
               "base station observation data for rtk\" and leaves every epoch "
               "Q=0.\n",nm(modename,prcopt.mode,10));
    }
    if (prcopt.mode>=PMODE_DGPS&&prcopt.mode<PMODE_PPP_KINEMA&&
        strtype[1]!=STR_NONE&&strfmt[1]!=STRFMT_RTCM2&&
        strfmt[1]!=STRFMT_RTCM3) {
        printf("WARNING: inpstr2-format=%s, not rtcm2/rtcm3. A base stream is "
               "normally RTCM corrections; a raw receiver format only decodes "
               "if it is the same receiver type as the base.\n",
               nm(strfmname,strfmt[1],(int)(sizeof(strfmname)/sizeof(*strfmname))));
    }
    if (prcopt.sateph==EPHOPT_PREC) {
        printf("ERROR: sateph=precise needs SP3/CLK files; rtksvr has no file "
               "input. Use brdc(0) for RTK, or ssr-apc(3)/ssr-com(4) with "
               "inpstr3.\n");
    }
    if ((prcopt.sateph==EPHOPT_SSRAPC||prcopt.sateph==EPHOPT_SSRCOM)&&
        strtype[2]==STR_NONE) {
        printf("WARNING: sateph=%s but inpstr3 (corr) is off - no SSR "
               "corrections will ever arrive.\n",
               nm(ephoptname,prcopt.sateph,5));
    }
    if (navmsgsel==2&&strtype[1]==STR_NONE) {
        printf("WARNING: misc-navmsgsel=base but inpstr2 is off - no "
               "ephemeris.\n");
    }
    if (navmsgsel==3&&strtype[2]==STR_NONE) {
        printf("WARNING: misc-navmsgsel=corr but inpstr3 is off - no "
               "ephemeris.\n");
    }
    if (prcopt.modear==ARMODE_OFF) {
        printf("NOTE: pos2-armode=off -> float solution only, no Q=1 epochs.\n");
    }
    if (prcopt.mode>=PMODE_DGPS&&prcopt.mode<PMODE_PPP_KINEMA&&
        prcopt.thresar[0]<1.0) {
        printf("ERROR: pos2-arthres=%.2f < 1.0 -> relpos() (rtkpos.c:1840) "
               "skips ambiguity resolution for every epoch, so Q=1 is "
               "impossible. Any value below 1.0 turns AR off; use 3.0.\n",
               prcopt.thresar[0]);
    }
    if ((prcopt.navsys&SYS_GLO)&&prcopt.glomodear==GLO_ARMODE_OFF) {
        printf("NOTE: GLONASS in navsys with pos2-gloarmode=off - GLO satellites "
               "are used for the float solution but never fixed.\n");
    }
    if (prcopt.refpos==POSOPT_FILE||prcopt.refpos==POSOPT_RINEX) {
        printf("ERROR: ant2-postype=%s needs input files, which rtksvr does not "
               "read. Use rtcm(5) for a caster mountpoint that sends RTCM 1006, "
               "or llh/xyz(0|1) with the surveyed position.\n",
               nm(postypname,prcopt.refpos,6));
    }
    if (prcopt.refpos==POSOPT_RTCM) {
        printf("NOTE: base position comes from RTCM 1006 on inpstr2; watch the "
               "antpos counter below - if it stays 0, rtk.rb is never set.\n");
    }
    if (prcopt.refpos==POSOPT_SINGLE&&
        (prcopt.mode==PMODE_STATIC||prcopt.mode==PMODE_STATIC_START)) {
        printf("WARNING: static mode with base position from averaged singles "
               "- the baseline carries that error. Set the surveyed llh/xyz.\n");
    }
    if (strtype[3]==STR_NONE&&strtype[4]==STR_NONE) {
        printf("NOTE: both outstr are off - the solution only goes to this "
               "console.\n");
    }
    if (tlev>0) {
        if (*filopt.trace) strcpy(tracefile,filopt.trace);
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
    printf("%-20s %3s %3s %11s %12s %9s %7s %7s %7s %6s %7s %s\n","time","Q","ns",
           "lat(deg)","lon(deg)","hgt(m)","sdN(m)","sdE(m)","sdU(m)","age(s)",
           "ratio","str");
    
    for (i=0;i<dur&&!stop_;i++) {
        double pos[3],rb[3];
        sol_t sol;
        int st[3],stat;
        
        sleepms(1000);
        rtksvrlock(&svr);
        sol=svr.rtk.sol;
        for (j=0;j<3;j++) {
            st[j]=svr.stream[j].state;
            rb[j]=svr.rtk.rb[j];
        }
        stat=sol.stat;
        rtksvrunlock(&svr);
        
        if (stat>=0&&stat<8) nq[stat]++;
        nep++;
        
        /* the caster sends RTCM 1006 late or not at all - report the change */
        if (norm(rb,3)>RE_WGS84/2.0) {
            double dr[3];
            for (j=0;j<3;j++) dr[j]=rb[j]-rb0[j];
            if (norm(dr,3)>0.001) {
                double l[3];
                ecef2pos(rb,l);
                printf("# base position from stream: %.7f %.7f %.3f\n",l[0]*R2D,
                       l[1]*R2D,l[2]);
                rb0[0]=rb[0]; rb0[1]=rb[1]; rb0[2]=rb[2];
            }
        }
        ecef2pos(sol.rr,pos);
        time2str(gpst2utc(sol.time),tstr,0); /* utc, to match the .pos file */
        printf("%-20s %3d %3d %11.7f %12.7f %9.3f %7.3f %7.3f %7.3f %6.1f %7.1f",
               tstr,sol.stat,sol.ns,pos[0]*R2D,pos[1]*R2D,pos[2],sol.qr[0],
               sol.qr[1],sol.qr[2],sol.age,sol.ratio);
        printf("  str=%d%d%d\n",st[0],st[1],st[2]);
        fflush(stdout);
    }
    printf("\nepochs printed: %d\n",nep);
    for (i=0;i<8;i++) {
        if (nq[i]) printf("  Q=%d %-6s : %d\n",i,solqname[i],nq[i]);
    }
    if (nq[SOLQ_FIX]+nq[SOLQ_FLOAT]>0) {
        printf("  fix rate: %.1f%% of %d fix|float epochs\n",
               100.0*nq[SOLQ_FIX]/(nq[SOLQ_FIX]+nq[SOLQ_FLOAT]),
               nq[SOLQ_FIX]+nq[SOLQ_FLOAT]);
    }
    printf("input stream message counts (obs,eph,glo-eph,antpos,dgps,ssr):\n");
    for (i=0;i<3;i++) {
        static const char *nmn[]={"rov ","base","corr"};
        printf("  %s: %d %d %d %d %d %d\n",nmn[i],svr.nmsg[i][0],svr.nmsg[i][1],
               svr.nmsg[i][6],svr.nmsg[i][4],svr.nmsg[i][5],svr.nmsg[i][7]);
    }
    rtksvrlock(&svr);
    strncpy(errbuf,svr.rtk.errbuf,MAXERRMSG-1);
    errbuf[MAXERRMSG-1]='\0';
    rtksvrunlock(&svr);
    if (*errbuf) printf("last rtk error: %s\n",errbuf);
    
    rtksvrstop(&svr,(const char **)cmds);
    rtksvrfree(&svr);
    if (tlev>0) traceclose();
    printf("stopped\n");
    return EXIT_SUCCESS;
}
