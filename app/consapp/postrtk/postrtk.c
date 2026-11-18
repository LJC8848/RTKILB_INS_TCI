/* ---------------------------------------------------------------------------
 * Name: postrtk.c
 * Function: Windows console shell for POST-PROCESSED RTK (relative modes), the
 *           twin of postppp.c and the file-based counterpart of rtppp.c.
 *           Like postppp it drives postpos() through the conf file only -
 *           rnx2rtkp.c falls back to the conf file-* keys only when no file is
 *           given on the command line (rnx2rtkp.c:219-225) and checks nothing
 *           relative-mode specific.
 * Input:  -k conf-file [option]... [file]...
 * Output: solution to -o file (stdout if omitted), postpos() progress on
 *         stderr via showmsg(), echo + summary on stdout
 * Notes:  pos1-posmode = dgps(1) / kinematic(2) / static(3) / static-start(4)
 *         / movingbase(5) / fixed(6); the conf file-obsfile list must be
 *         rover;base (semicolon separated, in that order).
 * -------------------------------------------------------------------------*/
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "rtklib.h"

#define PROGNM    "postrtk"
#define MAXFILE   64                  /* max number of input files */

static const char *modename[]={
    "single","dgps","kinematic","static","static-start","movingbase","fixed",
    "ppp-kine","ppp-static","ppp-fixed"};
static const char *postypname[]={
    "llh","xyz","single","posfile","rinexhead","rtcm"};
static const char *armodename[]={
    "off","continuous","instantaneous","fix-and-hold"};
static const char *solqname[]={
    "none","fix","float","sbas","dgps","single","ppp","dr"};

/* app callbacks required by the rtklib library ------------------------------*/
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
        "  file   overrides the conf file-obsfile/navfile/sp3file/clkfile list;\n"
        "         rover obs first, base obs second, then nav\n",PROGNM);
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
    char filebuf[MAXFILE][MAXSTRPATH],tracefile[64]="postrtk.trace";
    const char *infile[MAXFILE],*outfile="",*rov="",*base="",*conffile=NULL;
    int i,n=0,nobs=-1,ret,tlev=0;
    
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
        nobs=split_paths(filopt.obs,filebuf+n,infile+n,MAXFILE-n);
        n+=nobs;
        n+=split_paths(filopt.nav,filebuf+n,infile+n,MAXFILE-n);
        n+=split_paths(filopt.sp3,filebuf+n,infile+n,MAXFILE-n);
        n+=split_paths(filopt.clk,filebuf+n,infile+n,MAXFILE-n);
        n+=split_paths(filopt.erp,filebuf+n,infile+n,MAXFILE-n);
    }
    if (n<=0) {
        fprintf(stderr,"error: no input file (conf file-obsfile/navfile or "
                "command line)\n");
        return EXIT_FAILURE;
    }
    sprintf(solopt.prog,"%s ver.%s %s",PROGNM,VER_RTKLIB,PATCH_LEVEL);
    
    fprintf(stderr,"%s RTKLIB %s %s\n",PROGNM,VER_RTKLIB,PATCH_LEVEL);
    printf("mode=%d(%s) elmask=%.1f ionoopt=%d tropopt=%d navsys=%d nf=%d soltype=%d\n",
           prcopt.mode,0<=prcopt.mode&&prcopt.mode<=9?modename[prcopt.mode]:"?",
           prcopt.elmin*R2D,prcopt.ionoopt,prcopt.tropopt,prcopt.navsys,
           prcopt.nf,prcopt.soltype);
    printf("  obs : %s\n",nobs>=0?filopt.obs:"(command line)");
    printf("  nav : %s\n",*filopt.nav?filopt.nav:"(command line)");
    printf("  ant : sat=%s rcv=%s\n",*filopt.satantp?filopt.satantp:"-",
           *filopt.rcvantp?filopt.rcvantp:"-");
    printf("  base: postype=%d(%s)",prcopt.refpos,
           0<=prcopt.refpos&&prcopt.refpos<=5?postypname[prcopt.refpos]:"?");
    if (prcopt.refpos==POSOPT_POS_LLH||prcopt.refpos==POSOPT_POS_XYZ) {
        printf(" pos=%.3f %.3f %.3f (ecef m)",prcopt.rb[0],prcopt.rb[1],
               prcopt.rb[2]);
    }
    printf("\n");
    printf("  AR  : mode=%d(%s) glo=%d bds=%d thresar=%.2f\n",
           prcopt.modear,0<=prcopt.modear&&prcopt.modear<=3?
           armodename[prcopt.modear]:"?",prcopt.glomodear,prcopt.bdsmodear,
           prcopt.thresar[0]);
    if (prcopt.baseline[0]>0.0) {
        printf("  baseline constraint: %.3f m +/- %.3f m\n",prcopt.baseline[0],
               prcopt.baseline[1]);
    }
    printf("  out : %s\n",*outfile?outfile:"stdout");
    if (nobs>=0) printf("  %d input file(s), of which %d obs\n",n,nobs);
    else printf("  %d input file(s), from the command line\n",n);
    
    /* relative-mode preflight: warn, never abort, the conf owner decides ----*/
    if (prcopt.mode>=PMODE_PPP_KINEMA) {
        printf("WARNING: pos1-posmode=%s is a PPP mode - use postppp.\n",
               modename[prcopt.mode]);
    }
    else if (prcopt.mode<PMODE_DGPS) {
        printf("WARNING: pos1-posmode=%s is not a relative mode; the base obs "
               "file is then ignored.\n",modename[prcopt.mode]);
    }
    else if (nobs>=0&&nobs<2) {
        printf("ERROR: relative mode %s needs two obs files. postpos.c "
               "readobsnav() numbers receivers in file order - the first obs "
               "file that adds data is rcv=1 (rover), the next rcv=2 (base) - "
               "so file-obsfile must be \"rover;base\".\n",modename[prcopt.mode]);
    }
    if (prcopt.refpos==POSOPT_FILE&&!*filopt.stapos) {
        printf("ERROR: ant2-postype=posfile but file-staposfile is empty - "
               "antpos() (postpos.c:869) stops with \"no position of <name>\".\n");
    }
    if (prcopt.refpos==POSOPT_SINGLE) {
        printf("WARNING: base position is the average of single positions. For "
               "RTK fix, set ant2-postype=llh|xyz with the surveyed position.\n");
    }
    if (prcopt.mode==PMODE_FIXED&&prcopt.rovpos==POSOPT_SINGLE) {
        printf("WARNING: fixed mode with rovpos=single: the rover a-priori is "
               "the averaged single position. Set ant1-pos1/2/3.\n");
    }
    if (prcopt.modear==0) {
        printf("NOTE: pos2-armode=off -> float solution only, no Q=1 rows.\n");
    }
    for (i=0;i<n;i++) {
        if (strstr(infile[i],"%b")&&!*base) {
            printf("ERROR: input file %s has the %%b base keyword but no -base "
                   "list was given -> nothing would be processed.\n",infile[i]);
            return EXIT_FAILURE;
        }
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
