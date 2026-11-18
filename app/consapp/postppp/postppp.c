/* ---------------------------------------------------------------------------
 * Name: postppp.c
 * Function: Windows console shell for POST-PROCESSED PPP, the file-based twin
 *           of rtppp.c (which is the stream-based real-time PPP shell).
 *           rnx2rtkp.c already wraps postpos(), but it takes its input files
 *           from the command line and applies no PPP-specific checks. Here the
 *           files come from the conf file (file-obsfile/navfile/sp3file/
 *           clkfile/erpfile - in this fork those keys are read only by
 *           rnx2rtkp.c:220-224, so no other console app can use them), the
 *           options use the same sysopts table as rnx2rtkp, and the run prints
 *           the option echo + summary that rtppp users get.
 * Input:  -k conf-file [option]... [file]...
 * Output: solution to -o file (stdout if omitted), postpos() progress on
 *         stderr via showmsg(), echo + summary on stdout
 * Notes:  pos1-posmode = ppp-kine(7) / ppp-static(8) / ppp-fixed(9)
 * -------------------------------------------------------------------------*/
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "rtklib.h"

#define PROGNM    "postppp"
#define MAXFILE   64                  /* max number of input files */

static const char *modename[]={
    "single","dgps","kinematic","static","static-start","movingbase","fixed",
    "ppp-kine","ppp-static","ppp-fixed"};
static const char *ephname[]={
    "brdc","precise","brdc+sbas","brdc+ssrapc","brdc+ssrcom"};
static const char *solqname[]={
    "none","fix","float","sbas","dgps","single","ppp","dr"};

/* stream-independent app callbacks required by the rtklib library -----------*/
extern int showmsg(const char *format,...)
{
    va_list arg;
    va_start(arg,format); vfprintf(stderr,format,arg); va_end(arg);
    fprintf(stderr,"\r");
    return 0;
}
extern void settspan(gtime_t ts,gtime_t te) {(void)ts;(void)te;}
extern void settime(gtime_t time) {(void)time;}

static void usage(void)
{
    fprintf(stderr,
        "usage: %s -k conf [-o out.pos] [-ts y/m/d h:m:s] [-te y/m/d h:m:s]\n"
        "          [-ti sec] [-p mode] [-rov ids] [-base ids]\n"
        "          [-x trace-level] [-f trace-file] [file]...\n"
        "  conf   RTKPOST-style option file (same keys as rnx2rtkp -k)\n"
        "  file   overrides the conf file-obsfile/navfile/sp3file/clkfile list\n"
        "         first obs file = rover; PPP needs no second obs file\n",PROGNM);
}
/* split file paths separated by ';' into file list (from rnx2rtkp.c) --------*/
static int split_paths(const char *paths,char (*buf)[MAXSTRPATH],
                       const char **infile,int max)
{
    char tmp[MAXSTRPATH],*p;
    int n=0;
    
    if (*paths=='\0'||max<=0) return 0;
    strcpy(tmp,paths);
    for (p=strtok(tmp,";"); p&&n<max; p=strtok(NULL,";")) {
        if (*p=='\0') continue;
        strcpy(buf[n],p);
        infile[n]=buf[n];
        n++;
    }
    return n;
}
/* case-insensitive extension test, as postpos.c classifies files by suffix --*/
static int ext_is(const char *file,const char *ext)
{
    const char *p=strrchr(file,'.');
    size_t i,n=strlen(ext);
    
    if (!p||strlen(p)!=n) return 0;
    for (i=0;i<n;i++) {
        if (tolower((unsigned char)p[i])!=tolower((unsigned char)ext[i])) {
            return 0;
        }
    }
    return 1;
}
static int count_ext(const char **files,int n,const char *e1,const char *e2)
{
    int i,c=0;
    for (i=0;i<n;i++) {
        if (ext_is(files[i],e1)||(e2&&ext_is(files[i],e2))) c++;
    }
    return c;
}
/* count solution rows by status in the file just written --------------------*/
/* re-reading the output instead of trusting the filter also catches a run that
 * produced zero rows. Column order in this fork (solution.c:1215 outecef,
 * :1259 outpos, :1294 outenu): time(1|2) pos(3) STAT ns std(3) cov(3) age
 * ratio - the status comes BEFORE the std-devs, unlike stock RTKLIB 2.4.3. */
static void sol_summary(const char *file,const solopt_t *sopt)
{
    FILE *fp;
    char line[2048],delim[128],*tok;
    int i,iq,n=0,nq[8]={0};
    
    if (!file||!*file||sopt->posf==SOLF_NMEA||sopt->posf==SOLF_STAT) {
        printf("solutions: written to %s (summary not parsed for this format)\n",
               file&&*file?file:"stdout");
        return;
    }
    iq=(sopt->timef?2:1)+4; /* time + pos(3) + status */
    if (sopt->posf==SOLF_LLH&&sopt->degf) iq+=4; /* dms: lat,lon -> 3 tokens */
    strcpy(delim,*sopt->sep?sopt->sep:" ");
    strcat(delim," \t\r\n");
    if (!(fp=fopen(file,"r"))) {
        printf("solutions: %s not found\n",file);
        return;
    }
    while (fgets(line,sizeof(line),fp)) {
        if (*line=='%') continue;
        n++;
        for (i=0,tok=strtok(line,delim); tok; tok=strtok(NULL,delim)) {
            if (++i==iq) {
                int q=atoi(tok);
                if (q>=0&&q<8) nq[q]++;
                break;
            }
        }
    }
    fclose(fp);
    printf("solutions: %d rows\n",n);
    for (i=0;i<8;i++) {
        if (nq[i]) printf("  Q=%d %-6s : %d\n",i,solqname[i],nq[i]);
    }
}
int main(int argc,char **argv)
{
    prcopt_t prcopt=prcopt_default;
    solopt_t solopt=solopt_default;
    filopt_t filopt={""};
    gtime_t ts={0},te={0};
    double es[6]={2000,1,1,0,0,0},ee[6]={2000,12,31,23,59,59},tint=0.0;
    char filebuf[MAXFILE][MAXSTRPATH],tracefile[64]="postppp.trace";
    const char *infile[MAXFILE],*outfile="",*rov="",*base="",*conffile=NULL;
    int i,n=0,ret,tlev=0;
    
    for (i=1;i<argc;i++) {
        if      (!strcmp(argv[i],"-k")&&i+1<argc) {
            conffile=argv[++i];
            resetsysopts();
            if (!loadopts(conffile,sysopts)) return EXIT_FAILURE;
            getsysopts(&prcopt,&solopt,&filopt);
        }
        else if (!strcmp(argv[i],"-o")&&i+1<argc) outfile=argv[++i];
        else if (!strcmp(argv[i],"-ts")&&i+2<argc) {
            sscanf(argv[++i],"%lf/%lf/%lf",es,es+1,es+2);
            sscanf(argv[++i],"%lf:%lf:%lf",es+3,es+4,es+5);
            ts=epoch2time(es);
        }
        else if (!strcmp(argv[i],"-te")&&i+2<argc) {
            sscanf(argv[++i],"%lf/%lf/%lf",ee,ee+1,ee+2);
            sscanf(argv[++i],"%lf:%lf:%lf",ee+3,ee+4,ee+5);
            te=epoch2time(ee);
        }
        else if (!strcmp(argv[i],"-ti")&&i+1<argc) tint=atof(argv[++i]);
        else if (!strcmp(argv[i],"-p")&&i+1<argc) prcopt.mode=atoi(argv[++i]);
        else if (!strcmp(argv[i],"-rov")&&i+1<argc) rov=argv[++i];
        else if (!strcmp(argv[i],"-base")&&i+1<argc) base=argv[++i];
        else if (!strcmp(argv[i],"-x")&&i+1<argc) tlev=atoi(argv[++i]);
        else if (!strcmp(argv[i],"-f")&&i+1<argc) strcpy(tracefile,argv[++i]);
        else if (!strcmp(argv[i],"-?")||!strcmp(argv[i],"--help")) {
            usage(); return EXIT_SUCCESS;
        }
        else if (*argv[i]=='-') { usage(); return EXIT_FAILURE; }
        else if (n<MAXFILE) infile[n++]=argv[i];
    }
    if (!conffile) {
        fprintf(stderr,"error: no conf file (-k)\n");
        usage();
        return EXIT_FAILURE;
    }
    if (n<=0) { /* input files from the conf, same order as rnx2rtkp.c */
        n+=split_paths(filopt.obs,filebuf+n,infile+n,MAXFILE-n);
        n+=split_paths(filopt.nav,filebuf+n,infile+n,MAXFILE-n);
        n+=split_paths(filopt.sp3,filebuf+n,infile+n,MAXFILE-n);
        n+=split_paths(filopt.clk,filebuf+n,infile+n,MAXFILE-n);
        n+=split_paths(filopt.erp,filebuf+n,infile+n,MAXFILE-n);
    }
    if (n<=0) {
        fprintf(stderr,"error: no input file (conf file-* or command line)\n");
        return EXIT_FAILURE;
    }
    sprintf(solopt.prog,"%s ver.%s %s",PROGNM,VER_RTKLIB,PATCH_LEVEL);
    
    /* echo: what the conf actually selected (rtppp prints the same header) --*/
    fprintf(stderr,"%s RTKLIB %s %s\n",PROGNM,VER_RTKLIB,PATCH_LEVEL);
    printf("mode=%d(%s) sateph=%d(%s) ionoopt=%d tropopt=%d navsys=%d nf=%d soltype=%d\n",
           prcopt.mode,0<=prcopt.mode&&prcopt.mode<=9?modename[prcopt.mode]:"?",
           prcopt.sateph,0<=prcopt.sateph&&prcopt.sateph<=4?ephname[prcopt.sateph]:"?",
           prcopt.ionoopt,prcopt.tropopt,prcopt.navsys,prcopt.nf,prcopt.soltype);
    printf("  obs : %s\n",*filopt.obs?filopt.obs:"(command line)");
    printf("  nav : %s\n",*filopt.nav?filopt.nav:"(command line)");
    printf("  sp3 : %s\n",*filopt.sp3?filopt.sp3:"-");
    printf("  clk : %s\n",*filopt.clk?filopt.clk:"-");
    printf("  ant : sat=%s rcv=%s\n",*filopt.satantp?filopt.satantp:"-",
           *filopt.rcvantp?filopt.rcvantp:"-");
    printf("  out : %s\n",*outfile?outfile:"stdout");
    printf("  %d input file(s)\n",n);
    
    /* PPP-specific preflight: warn, never abort, the conf owner decides -----*/
    if (prcopt.mode<PMODE_PPP_KINEMA) {
        printf("WARNING: pos1-posmode=%s is not a PPP mode (7/8/9). For RTK "
               "use postrtk.\n",modename[prcopt.mode]);
    }
    if (prcopt.sateph==EPHOPT_PREC&&
        count_ext(infile,n,".sp3",".eph")==0) {
        printf("ERROR: sateph=precise but no .sp3/.eph in the input list. "
               "readpreceph() (postpos.c:673) runs readsp3() over the files and "
               "leaves nav.ne=0 -> no precise orbits.\n");
    }
    if (prcopt.sateph==EPHOPT_BRDC) {
        printf("WARNING: sateph=brdc - broadcast orbit/clock PPP is dm-level, "
               "not cm-level. Use sateph=precise with SP3+CLK.\n");
    }
    if (prcopt.sateph==EPHOPT_SSRAPC||prcopt.sateph==EPHOPT_SSRCOM) {
        if (count_ext(infile,n,".rtcm3",NULL)==0) {
            printf("WARNING: sateph=%s needs an SSR/RTCM3 message log in the "
                   "input list; without it no corrections are applied.\n",
                   ephname[prcopt.sateph]);
        }
        if (prcopt.soltype!=SOLTYPE_FORWARD) {
            printf("WARNING: SSR corrections are valid for forward estimation "
                   "only (postpos.c:1379).\n");
        }
    }
    if (prcopt.ionoopt==IONOOPT_OFF) {
        printf("WARNING: ionoopt=off - model_iono() returns 0 for IONOOPT_OFF "
               "(ppp.c:941), ppp_res() then skips every satellite (ppp.c:1007) "
               "and the run yields zero rows.\n");
    }
    if ((prcopt.posopt[0]||prcopt.posopt[1])&&!*filopt.satantp) {
        printf("WARNING: pos1-posopt1/2=on but file-satantfile is empty - PCO "
               "corrections silently become zero.\n");
    }
    if (prcopt.mode==PMODE_PPP_FIXED&&norm(prcopt.ru,3)<=0.0) {
        printf("ERROR: ppp-fixed takes the a-priori straight from prcopt.ru "
               "(ppp.c:536); ant1-postype=llh|xyz + ant1-pos1/2/3 fill it. "
               "They are empty -> the fixed position would be 0,0,0.\n");
    }
    for (i=0;i<n;i++) {
        if (strstr(infile[i],"%r")&&!*rov) {
            printf("ERROR: input file %s has the %%r rover keyword but no -rov "
                   "list was given -> nothing would be processed.\n",infile[i]);
            return EXIT_FAILURE;
        }
    }
    if (tlev>0) {
        if (*filopt.trace) strcpy(tracefile,filopt.trace);
        traceopen(tracefile);
        tracelevel(tlev);
    }
    ret=postpos(ts,te,tint,0.0,&prcopt,&solopt,&filopt,infile,n,outfile,rov,base);
    
    sol_summary(outfile,&solopt);
    if (ret==1) printf("aborted by user\n");
    else if (ret!=0) printf("%s error (ret=%d)\n",PROGNM,ret);
    if (tlev>0) traceclose();
    return ret?EXIT_FAILURE:EXIT_SUCCESS;
}
