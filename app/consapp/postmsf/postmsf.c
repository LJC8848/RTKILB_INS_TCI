/* ---------------------------------------------------------------------------
 * Name: postmsf.c
 * Function: Windows console shell for post-processed loosely-coupled
 *           GNSS/INS (SINS mechanization + 15-state error-state EKF, src/ins).
 *           It is postrtk.c plus the INS lifecycle: the only RTKLIB change is
 *           one guarded call in postpos.c:procpos(), so this app behaves
 *           exactly like postrtk when no IMU file is given.
 * Input:  -k conf-file [-o out.pos] [option]... [file]...
 *         conf-file carries the RTKPOST keys AND the src/ins keys, so one run
 *         needs one file; -ic ins-conf only splits or overrides them.
 * Output: GNSS solution to -o file (stdout if omitted), INS solution to the
 *         -io file, summary on stdout, postpos() progress on stderr
 * Notes:  the measurement variances come from sol_t::qr/qv, which the
 *         estimators fill regardless of the output options (rtkpos.c:2176-2208,
 *         pntpos.c:451-454,622-627); pos1-posmode may be any postpos mode.
 *         Only the forward analysis feeds the filter, so for a combined
 *         pos4-soltype the .ins file follows the forward pass.
 * -------------------------------------------------------------------------*/
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "rtklib.h"
#include "ins/rtkins.h"

#define PROGNM    "postmsf"
#define MAXFILE   64                  /* max number of input files */

static const char *modename[]={
    "single","dgps","kinematic","static","static-start","movingbase","fixed",
    "ppp-kine","ppp-static","ppp-fixed"};
static const char *alignname[]={"auto","att","pva","vva"};
static const char *measname[]={"auto","pos","vel","pos-vel"};

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
        "usage: %s -k conf [-o out.pos] [-i imu-file]\n"
        "          [-io out.ins] [-ic ins-conf]\n"
        "          [-align auto|att|pva|vva] [-meas auto|pos|vel|pos-vel]\n"
        "          [-fb self|msf] [-ts y/m/d h:m:s] [-te y/m/d h:m:s]\n"
        "            (dates/times also accept - and . as separators)\n"
        "          [-ti sec] [-p mode] [-rov ids] [-base ids]\n"
        "          [-x trace-level] [-f trace-file] [file]...\n"
        "  conf      RTKPOST-style options plus the src/ins keys, one file\n"
        "  ins-conf  separate INS option file, read after conf (\"key = value\"\n"
        "            lines, see src/ins/inscfg.c)\n"
        "  imu-file  overrides imu-file (ASCII: t dtheta(3) dvel(3) or rates)\n"
        "  out.ins   overrides out-ins (\"\" or none = no ins output)\n"
        "  file      overrides the conf file-obsfile/navfile/sp3file/clkfile list;\n"
        "            rover obs first, base obs second, then nav\n"
        "  -ts/-te     both required: postpos() has no data-derived default\n",PROGNM);
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
static int isword(const char *arg,const char *w)
{
    int i;
    if (!arg||!*arg) return 0;
    for (i=0;w[i]!='\0';i++) {
        if (tolower((unsigned char)arg[i])!=tolower((unsigned char)w[i])) return 0;
    }
    return 1;
}
/* -ts/-te: accept y/m/d as well as y-m-d and h:m:s as well as h.m.s.  A pure
 * sscanf("%lf/%lf/%lf") on "2020-10-27" converts only the year and leaves the
 * rest of the default array in place, which gave a 2000-01-01 window, zero
 * solutions and no output file at all - silently.  Return 0 unless all three
 * fields are present so the caller can abort instead. */
static int parsetimef(const char *s,double *v)
{
    char tmp[64];
    int i;
    if (!s||strlen(s)>=sizeof(tmp)) return 0;
    strcpy(tmp,s);
    for (i=0;tmp[i]!='\0';i++) {
        if (tmp[i]=='/'||tmp[i]=='-'||tmp[i]==':'||tmp[i]=='.') tmp[i]=' ';
    }
    return sscanf(tmp,"%lf%lf%lf",v,v+1,v+2)==3;
}
int main(int argc,char **argv)
{
    prcopt_t prcopt=prcopt_default;
    solopt_t solopt=solopt_default;
    filopt_t filopt={""};
    insopt_t insopt;
    gtime_t ts={0},te={0};
    double es[6]={2000,1,1,0,0,0},ee[6]={2000,12,31,23,59,59},tint=0.0;
    char filebuf[MAXFILE][MAXSTRPATH],tracefile[64]="postmsf.trace";
    char insconffile[1024]="",imofile[1024]="",insoutfile[1024]="";
    const char *infile[MAXFILE],*outfile="",*rov="",*base="",*conffile=NULL;
    int i,n=0,nobs=-1,ret,tlev=0,align=-1,meas=-1,fb=0;
    
    ins_defaults(&insopt);
    
    for (i=1;i<argc;i++) {
        if      (!strcmp(argv[i],"-k")&&i+1<argc) {
            conffile=argv[++i];
            resetsysopts();
            if (!loadopts(conffile,sysopts)) {
                fprintf(stderr,"error: cannot open conf %s\n",conffile);
                return EXIT_FAILURE;
            }
            getsysopts(&prcopt,&solopt,&filopt);
        }
        else if (!strcmp(argv[i],"-ic")&&i+1<argc) strcpy(insconffile,argv[++i]);
        else if (!strcmp(argv[i],"-i")&&i+1<argc) strcpy(imofile,argv[++i]);
        else if (!strcmp(argv[i],"-io")&&i+1<argc) strcpy(insoutfile,argv[++i]);
        else if (!strcmp(argv[i],"-align")&&i+1<argc) {
            const char *w=argv[++i];
            if      (isword(w,"auto")) align=INS_ALIGN_AUTO;
            else if (isword(w,"att"))  align=INS_ALIGN_ATT;
            else if (isword(w,"pva"))  align=INS_ALIGN_PVA;
            else if (isword(w,"vva"))  align=INS_ALIGN_VVA;
            else { fprintf(stderr,"error: bad -align %s\n",w); return EXIT_FAILURE; }
        }
        else if (!strcmp(argv[i],"-meas")&&i+1<argc) {
            const char *w=argv[++i];
            if      (isword(w,"pos-vel")) meas=INS_MEAS_POS_VEL;
            else if (isword(w,"pos"))     meas=INS_MEAS_POS;
            else if (isword(w,"vel"))     meas=INS_MEAS_VEL;
            else if (isword(w,"auto"))    meas=INS_MEAS_AUTO;
            else { fprintf(stderr,"error: bad -meas %s\n",w); return EXIT_FAILURE; }
        }
        else if (!strcmp(argv[i],"-fb")&&i+1<argc) {
            const char *w=argv[++i];
            if      (isword(w,"self")) fb=-1;
            else if (isword(w,"msf"))  fb=1;
            else { fprintf(stderr,"error: bad -fb %s (self|msf)\n",w);
                   return EXIT_FAILURE; }
        }
        else if (!strcmp(argv[i],"-o")&&i+1<argc) outfile=argv[++i];
        else if (!strcmp(argv[i],"-ts")&&i+2<argc) {
            const char *sd=argv[++i],*st=argv[++i];
            if (!parsetimef(sd,es)||!parsetimef(st,es+3)) {
                fprintf(stderr,"error: bad -ts %s %s (want y/m/d h:m:s)\n",sd,st);
                return EXIT_FAILURE;
            }
            ts=epoch2time(es);
        }
        else if (!strcmp(argv[i],"-te")&&i+2<argc) {
            const char *sd=argv[++i],*st=argv[++i];
            if (!parsetimef(sd,ee)||!parsetimef(st,ee+3)) {
                fprintf(stderr,"error: bad -te %s %s (want y/m/d h:m:s)\n",sd,st);
                return EXIT_FAILURE;
            }
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
    /* one file carries both option sets: loadopts() ignores keys it does not
     * know (options.c:370) and ins_readconf() ignores the RTKLIB keys
     * (inscfg.c parse_line), so reading the -k conf for the INS options too is
     * safe.  -ic remains available for a separate file and wins on conflict. */
    ins_defaults(&insopt);
    ins_readconf(conffile,&insopt);
    if (*insconffile&&strcmp(insconffile,conffile)) {
        if (ins_readconf(insconffile,&insopt)<0) {
            fprintf(stderr,"error: cannot open ins conf %s\n",insconffile);
            return EXIT_FAILURE;
        }
    }
    if (*imofile) strcpy(insopt.imufile,imofile);
    if (*insoutfile) strcpy(insopt.outfile,insoutfile);
    if (align>=0) insopt.align=align;
    if (meas>=0) insopt.meas=meas;
    if (fb) insopt.fbsign=fb;
    
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
    /* Both ends are required here, and their absence used to be silent:
     * reppaths() returns 0 files when ts.time==0||te.time==0 (rtkcmn.c:3575)
     * and postpos()'s session loop breaks at timediff(tts,te)>0
     * (postpos.c:1427), so a missing -ts or -te - or -te < -ts - gave an empty
     * run, a header-only output file and no message.  This fork's postpos has
     * no tspan_file() to derive the span from the data, unlike stock rnx2rtkp. */
    if (!ts.time||!te.time) {
        fprintf(stderr,"error: -ts and -te must both be given (missing %s)\n",
                ts.time?"-te":"-ts");
        usage();
        return EXIT_FAILURE;
    }
    if (timediff(te,ts)<0.0) {
        fprintf(stderr,"error: -te is before -ts, empty processing window\n");
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
    printf("  out : %s\n",*outfile?outfile:"stdout");
    if (!*insopt.imufile) {
        printf("  ins : OFF (no imu-file key in %s; -ic/-i are required for "
               "loosely coupling)\n",conffile);
    }
    else {
        printf("  ins : imu=%s ts=%.4f s inc=%d nsub=%d\n",insopt.imufile,
               insopt.ts,insopt.inc,insopt.nsub);
        printf("        align=%s att0=%.2f %.2f %.2f lever=%.3f %.3f %.3f (m)\n",
               0<=insopt.align&&insopt.align<=3?alignname[insopt.align]:"?",
               insopt.att0[0],insopt.att0[1],insopt.att0[2],
               insopt.lever[0],insopt.lever[1],insopt.lever[2]);
        printf("        meas=%s max-dt=%.3f max-nr=%.1f fb-sign=%d inflation=%.1f\n",
               0<=insopt.meas&&insopt.meas<=3?measname[insopt.meas]:"?",
               insopt.maxdt,insopt.maxnr,insopt.fbsign,insopt.infl);
        printf("        out=%s\n",*insopt.outfile?insopt.outfile:"(none)");
    }
    if (*insopt.imufile&&prcopt.mode==PMODE_SINGLE) {
        printf("NOTE: single positioning: sol->qv is only filled when Doppler "
               "velocity is estimated (pos1-dop=estvel), so the measurement "
               "type auto may degrade to position-only.\n");
    }
    if (tlev>0) {
        if (*filopt.trace) strcpy(tracefile,filopt.trace);
        traceopen(tracefile);
        tracelevel(tlev);
    }
    if ((ret=ins_start(&insopt))<0) {
        fprintf(stderr,"error: ins: %s\n",ins_errmsg());
        if (tlev>0) traceclose();
        return EXIT_FAILURE;
    }
    ret=postpos(ts,te,tint,0.0,&prcopt,&solopt,&filopt,infile,n,outfile,rov,base);
    
    ins_stop();
    if (ret==1) printf("aborted by user\n");
    else if (ret!=0) printf("%s error (ret=%d)\n",PROGNM,ret);
    if (tlev>0) traceclose();
    return ret?EXIT_FAILURE:EXIT_SUCCESS;
}
