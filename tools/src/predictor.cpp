
/*******************************************************************************

  T-UGO tools, 2006-2018

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

<!-- USE FIND AND REPLACE TO EDIT THIS LINE, SO THAT print_help IS ALSO UPDATED -->
\brief Interpolates constants and predicts tides at given points with given atlases

<!-- A LINK TO main() or print_help() WILL NOT LINK TO THE RIGHT SOURCE ! -->
See the main function for how this works
and the print_help function for how to use this.
*/
/*----------------------------------------------------------------------------*/

#define MAIN_SOURCE

#include "version-macros.def" //for VERSION and REVISION

#include <unistd.h> /* for access */
#include <stdio.h>
#include <glob.h>
#include <time.h> /* for time(), gmtime_r(), strftime() */

#include "poc-netcdf-data.hpp"
#include "tides.h"
#include "functions.h" //safely includes omp.h
#include "matrix.h" // for pos
#include "netcdf-proto.h" /* for poc_get_lon_lat_time() */
#include "poc-time.h"
#include "mgr.h"

#include "tides.def"


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  vector<mgr_t> mgr(vector<hconstant_t> h, spectrum_t spectrum)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  vector<mgr_t> gauges;
  mgr_create(h.size(), gauges, spectrum);
  
  for (int i=0;i<h.size();i++){
    mgr_t mgr;
    for (int k=0;k<spectrum.n;k++) {
      gauges[i].data[k].amp=100*h[i].a[k];
      double G=h[i].G[k];
      if      (G<0.0)   gauges[i].data[k].phi=G+360.0;
      else if (G>360.0) gauges[i].data[k].phi=G-360.0;
      else              gauges[i].data[k].phi=G;
      }
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  vector<mgr_t> mgr_create(hconstant_t *h, double *lon, double *lat, int nlocations, spectrum_t spectrum)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  vector<mgr_t> gauges;
  mgr_create(nlocations, gauges, spectrum, lon, lat, 0);
  
  for (int i=0;i<nlocations;i++){
    sprintf(gauges[i].name, "%s%4.4d","extraction-",i);
    for (int k=0;k<spectrum.n;k++) {
      gauges[i].data[k].amp=100*h[i].a[k];
      double G=h[i].G[k];
      if      (G<0.0)   gauges[i].data[k].phi=G+360.0;
      else if (G>360.0) gauges[i].data[k].phi=G-360.0;
      else              gauges[i].data[k].phi=G;
      }
    }
  
  return(gauges);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void print_help(const char *prog_name)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/** \brief prints help of this programme.

\param prog_name name of the program to be printed by this help function : argv[0]
*/
/*----------------------------------------------------------------------------*/
{
/** The code of the body of this function is :
    \code /**/ // COMPILED CODE BELOW !!!
  printf("\n"
    "NAME AND VERSION\n"
    "  %s version " VERSION " " REVISION "\n", prog_name);
  printf("\n"
    "USE\n"
    "  %s [-p lon_lat_list] -a atlas_convention [-s start -f end] -w wave1 [wave2 ... ] [OPTIONS] \n",prog_name);
  printf("\n"
    "DESCRIPTION\n"
/* USE FIND AND REPLACE TO EDIT THE LINE BELOW, SO THAT THE FILE HEAD IS ALSO UPDATED */
    "  Interpolates constants and predicts tides at given points with given atlases.\n"
    "  If start and end dates are provided, predicts tides.\n"
    "  It activates OpenMP parallelisation only if they are many points per time step.\n"
    "  If the output format is NetCDF, the 'units' attribute of the output variable is set to \"m\".\n"
/* should have linear, spline and GDR */
    "\n"
    "OPTIONS :\n"
    "  -h,--help : show this help and exit\n"
    "  --verbose : force printing constants when in altimetry mode (see option -p)\n"
    "  --nodal=no : do not do nodal corrections\n"
    "  --spring-neap : print spring/neap times of the first point.\n"
    "  --time : followed by the format of the time if the file is ASCII. CNES is replaced by the number of days since 1950/01/01 00:00; ELAPSED, by the number of days since the start date and CALENDAR by the date in yyyy/mm/dd HH:MM:SS.S format. Default: `ELAPSED CALENDAR`\n"
    "  -p : followed by the path of the list of control points. It is either\n"
    "        an ASCII file with the number of control points followed by their coordinates (longitude latitude)\n"
    "        or a pattern of NetCDF files with latitude, longitude and time variables (see also option -o).\n"
    "        If not given, predict on all points of the atlas.\n"
    "  -m : followed by the path of the extracted constants\n"
    "  -a : followed by the atlas file name convention. See below.\n"
    "  -cm : use if atlases are in cm\n"
    "  -s : followed by the start date. See DATE FORMATS below\n"
    "  -f : followed by the end date. See DATE FORMATS below\n"
    "  -i : followed by the date increment concatenated with the unit: s (default), m, h or d. Default increment: 3600s=60m=1h.\n"
    "  -v : followed by the variables names for the amplitude and the phase. Default: Ha Hg\n"
    "  -w : followed by the list of waves to predict for. To optimise speed, sort the waves by type of atlas (sorting by size of atlas should do the same)\n"
    "  -o : followed by the path of the output.\n"
    "        It is\n"
    "          NetCDF if -p is not given or specifies a NetCDF file with a time variable\n"
    "          ASCII otherwise\n"
    "        Default: predictions.dat or predictions.nc or %%s-predictions.nc, with %%s the base name of the patterned input file.\n"
    "        If empty, add prediction to the patterned input file.\n"
    "  --output-var : followed by the name of the NetCDF output variable. Default: prediction\n"
    "\n"
    "TIPS\n"
    "  To get all available atlases :\n"
    "ext=.FES2014b.nc;d=/data/soa/maree/fes2014;f=(`cd $d;ls -S *$ext`);%s -p control.dat -a $d/WAVE$ext -v elevation_a elevation_G -s ... -f ... -i ... -w ${f[@]/$ext}\n",prog_name);
  printf(
    "  In `ipython --pylab=...', plot 1 point with the following command, after replacing yyyy/mm/dd with the start date:\n"
    "t,eta=loadtxt('predictions.dat',usecols=(0,3),unpack=True);plot_date(t+datestr2num('yyyy/mm/dd'),eta,'-')\n"
    );
  print_tide_decode_atlasname_help();
  print_poctime_scan_date_help(0);
  print_OPENMP_help(prog_name);
    /** \endcode */
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void SpringNeap_cycle(astro_angles_t astro_angles, spectrum_t s, hconstant_t *constants, int nodal, double start, double final)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  tidal_wave wave,wM2,wS2;
  double V0,omega;
  double f,Vu;
  double VM2,VS2;
  double time=0;
  
/**--------------------------------------------------------------------------
  compute spring/neap tide calendar */

  int M2id=s.wave_index("M2");
  int S2id=s.wave_index("S2");
  bool abortThis=false;
  if(M2id<0){
    STDERR_BASE_LINE("Can not find M2 in the list of waves:\n");
    abortThis=true;
    }
  if(S2id<0){
    STDERR_BASE_LINE("Can not find S2 in the list of waves:\n");
    abortThis=true;
    }
  if(abortThis){
    STDERR_BASE_LINE("aborting %s\n",__func__);
    return;
    }
  
  wM2=s.wave("M2");
//  wM2.init();
  wS2=s.wave("S2");
  
  double d=wS2.omega-wM2.omega;
//  double r=2*(wS2.omega-wM2.omega)/(wS2.omega+wM2.omega);
  
  double duration=360./24./d;
  printf("spring/neap tides cycle: %lf days\n",duration);
  
  wave=wM2;
  omega=wave.omega*dph2rpd;
  V0=greenwhich_argument(astro_angles,wave);
  if(nodal==1) {
    Vu=nodal_phase(astro_angles,wave);
    f=nodal_factor(astro_angles,wave.formula);
    }
  else {
    Vu=0.;
    f=1.;
    }
  VM2=omega*time/3600.+V0+Vu;
  
  wave=wS2;
  omega=wave.omega*dph2rpd;
  V0=greenwhich_argument(astro_angles,wave);
  if(nodal==1) {
    Vu=nodal_phase(astro_angles,wave);
    f=nodal_factor(astro_angles,wave.formula);
    }
  else {
    Vu=0.;
    f=1.;
    }
  VS2=omega*time/3600.+V0+Vu;
  
/**--------------------------------------------------------------------------
  M2 & S2 phases at start (time=0) */
  VM2-=constants[0].G[M2id]*d2r;
  VS2-=constants[0].G[S2id]*d2r;
  
  double dV=fmod(VM2-VS2,2*M_PI);
  dV-=2*M_PI;
  
/**--------------------------------------------------------------------------
  convert seconds into julian days */
  start/=d2s;
  final/=d2s;
  
  double t=start;
  while(t<final) {
/**--------------------------------------------------------------------------
    in-phase M2 S2 */
    double spring=dV/(d*dph2rpd);
    printf("spring tide : %10.5lf %10.5lf julian day/cnes time\n",spring,spring+start);
    dV+=M_PI;
/**--------------------------------------------------------------------------
    in-opposition M2 S2 */
    double neap=dV/(d*dph2rpd);
    printf("neap tide   : %10.5lf %10.5lf julian day/cnes time\n",neap,neap+start);
    dV+=M_PI;
    t=neap+start;
    }
  
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// function that will be ran when the executable is started
/** See the print_help <a href=#func-members>function</a> for its use.
*/
/*----------------------------------------------------------------------------*/
{
  int pI,wI,n,status,k;        //point index, wave index, argument index, status, variable index
  int pn=-1;                   //numbers of : points, times
  int verbose=0;

  const char *position_file=0;  //path of the lat and lon of the control points
  glob_t globbed;
  globbed.gl_pathv=0;
  vector<int> pnS;
  int fI,fpn;
  
  double *times=0;
  const char *mgrfile=0;        //path of the extracted constants
  hconstant_t *constants=0;     //constants at the control points
  double *buffer;

  FILE *F=0;

  char *atlas_template=NULL;    //atlas path convention
  const char *time_template=NULL;
  date_t start=NADate,final=NADate;
  double startd=NAN,finald=NAN,increment=NAN,t;
  char *startc=0;
  char *unit;                  //unit of time increment
  char *waveNames[100];        //list of wave names
  spectrum_t WaveList;         //list of waves
  int nodal_corrections=1;
  astro_angles_t astro_angles;
  const char *output=NULL;     //path of the output
  bool outputIsNC;
  const char *keyword,*s;      //option and the argument that follows

  const int nvars=2;
  char *varnames[nvars]={strdup("Ha"),strdup("Hg")};
  string outputVarName="prediction";
  int defaultVars=1;
  bool SpringNeap=false;
  double scaleAtlas=1.;
  
  string cmd;
  cmd=fct_echo( argc, argv);

  wI=0;
  //waveNames[wI++]=strdup("Z0");
  waveNames[wI]=NULL;

  n=1;
  while (n < argc) {
    keyword=argv[n];
    if(strcmp(keyword,"--nodal=no")==0) {
      nodal_corrections=0;
      n++;
      }
    else if(strcmp(keyword,"--verbose")==0) {
      verbose=1;
      n++;
      }
    else if(strcmp(keyword,"--spring-neap")==0) {
      SpringNeap=true;
      n++;
      }
    else if(strcmp(keyword,"-cm")==0) {
      scaleAtlas=0.01;
      n++;
      }
    else if(strcmp(keyword,"--time")==0) {
      time_template=argv[n+1];
      n++;
      n++;
      }
    else if(strcmp(keyword,"--output-var")==0) {
      outputVarName=argv[n+1];
      n++;
      n++;
      }
    else switch (keyword[0]) {
      case '-':
        if( strcmp(keyword,"-h")==0 or
            strcmp(keyword,"--help")==0 ){
          print_help(argv[0]);
          return 0;
          }
        
        switch (keyword[1]) {

        case 'p' :
          position_file=argv[n+1];
          n++;
          n++;
          break;

        case 'a' :
          atlas_template=strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 's' :
          s=argv[n+1];
          n++;
          n++;
          status=poctime_scan_date(s,&start,0);
          break;

        case 'f' :
          s=argv[n+1];
          n++;
          n++;
          status=poctime_scan_date(s,0,&final);
          break;

        case 'i' :
          s=argv[n+1];
          n++;
          n++;
          increment=strtod(s,&unit);
          switch(tolower(*unit)){
            case 'd':
              increment*=24.;
            case 'h':
              increment*=60.;
            case 'm':
              increment*=60.;
            case 0:
            case 's':
              printf("Time increment of %g s.\n",increment);
              break;
            default:
              fprintf(stderr,"unit %s not recognised.\n");
              print_help(argv[0]);
              exit(1);
            }
          break;

        case 'v' :
          if(defaultVars)
            defaultVars=0;
          else __FILE_LINE__(stdout,"multiple use of -v : the last one overrides the previous.");
          
          for(k=0;k<nvars;k++){
            free(varnames[k]);
            varnames[k]=strdup(argv[n+1+k]);
            }
          
          n++;
          n+=nvars;
          break;

        case 'w' :
          n++;
          wI=pos((char*)NULL,waveNames,100);
          for(;n<argc;wI++,n++) {
            waveNames[wI]=strdup(argv[n]);
            }
          waveNames[wI]=NULL;
          break;

        case 'o' :
          output=argv[n+1];
          n++;
          n++;
          break;

        case 'm' :
          mgrfile=argv[n+1];
          n++;
          n++;
          break;

        default:
          printf("unknown option %s\n",keyword);
          print_help(argv[0]);
          exit(-1);
        }
        break;

      default:
        printf("unknown option %s\n",keyword);
        print_help(argv[0]);
        exit(-1);
      }
    
    }

  if(isad(start)!=isad(final)){
    fprintf(stderr,"*** Provide start AND end dates for predictions ***\n");
    print_help(argv[0]);
    exit(-1);
    }
/*---------------------------------------------------------------------*//**<h1>
  initialise list of waves </h1>*/
  if(waveNames[0]==NULL){
    fprintf(stderr,"*** Please provide a list of waves. ***\n");
    print_help(argv[0]);
    exit(-1);
    }
#if USE_M1_12
  WaveList.init(initialize_tide(),waveNames,-2);
#else
  WaveList.init(initialize_tide(),waveNames);
#endif
  printf("# number of waves : %d \n",WaveList.n);

/*------------------------------------------------------------------------------
  initialise list of points */
  if(position_file==NULL && (isnad(start) || isnad(final))){
    fprintf(stderr,
      "*** Please provide either : ***\n"
      "*** - a list of points      ***\n"
      "*** - start and end dates.  ***\n");
    print_help(argv[0]);
    exit(-1);
    }
  
  if(atlas_template==NULL){
    fprintf(stderr,"*** Please provide a convention for the name of the atlases. ***\n");
    print_help(argv[0]);
    exit(-1);
    }
  
  poc_var_t tv;/* time variable */
  tv.init("time",NC_DOUBLE,"",(string)"seconds since "+sgetdate(start));
  poc_var_t bv(outputVarName,NC_FLOAT);/* buffer variable */
  
  if(position_file!=0){
/*------------------------------------------------------------------------------
    prediction at given locations */
    double *lat,*lon;
    int strict;
    
    if(strrncasecmp(position_file,".nc")==0){
/*------------------------------------------------------------------------------
      NetCDF */
      
      status=glob(position_file,0,0,&globbed);
      
      if(status!=0 or globbed.gl_pathv==0){
        const char *msg;
        /* see man:/glob */
        switch(status){
        case GLOB_NOSPACE: msg="out of memory";break;
        case GLOB_ABORTED: msg="read error";break;
        case GLOB_NOMATCH: msg="no match";break;
        default: msg="unknown error code";
          }
        TRAP_ERR_EXIT(status,"glob(\"%s\",0,0,) error (%d %s)\n",position_file,status,msg);
        }
      
      pn=0;
      for(fI=0;;fI++){
        position_file=globbed.gl_pathv[fI];
        if(position_file==0)
          break;
        //if(fI>0) printf("%s%s",cuu1,cuu1);
        printf("[%d]",fI+1);
        fpn=poc_get_lon_lat_time(position_file); fflush(stdout);
        pn+=fpn;
        pnS.push_back(fpn);
        }
      
      const int pnSn=pnS.size();
      STDOUT_BASE_LINE("Allocating for %d points (in %d files)... ",pn,pnSn);fflush(stdout);
      exitIfNull(
        lon=new double[pn]
        );
      exitIfNull(
        lat=new double[pn]
        );
      exitIfNull(
        times=new double[pn]
        );
      printf("done.\n");
      
      pn=poc_get_lon_lat_time(globbed.gl_pathv[0],lon,lat,times,&bv.dimensions);
      for(fI=1;fI<pnSn;fI++){
        printf("%s%s",cuu1,cuu1);
        printf("[%d/%d]",fI+1,pnSn);
        position_file=globbed.gl_pathv[fI];
        fpn=poc_get_lon_lat_time(position_file,&lon[pn],&lat[pn],&times[pn]);
        fflush(stdout);
        pn+=fpn;
        if(fpn!=pnS[fI]) NC_TRAP_ERROR(wexit,NC_EDIMSIZE,1,"file %d %s was %d and now is %d\n",fI,position_file,pnS[fI],fpn);
        }
      
      bv.dimensions[0].len=pn;
      bv << poc_att_t("units","m");
      bv << poc_att_t("production","constructed around " __LINE_FILE_PACKAGE_REVISION);
      }
    else{
/*------------------------------------------------------------------------------
      ASCII */
      char *s, line[1024];
      F=fopen(position_file,"r");
      if(F==NULL){
        STDOUT_BASE_LINE("Can not open list of points : %s\n",position_file);
        exit(1);
        }
      
      s=fgets(line,1024,F);
      
      sscanf(line,"%d",&pn);
      printf("# number of points found in %s : %d \n",position_file,pn);
      
      lat=new double[pn];
      lon=new double[pn];
      
      for(pI=0;pI<pn;pI++){
        s=fgets(line,1024,F);
        sscanf(line,"%lf %lf\n",&lon[pI],&lat[pI]);
        printf("point %4d, lat=%10.3lf, lon=%10.3lf\n",pI,lat[pI],lon[pI]);
        }
      
      fclose(F);
      }

/*------------------------------------------------------------------------------
    extract harmonic constants at prediction locations */
    if(times==0)
      status=4|verbose;
    else
      status=verbose;
    strict= 1 or (times!=0);
    constants=tide_atlas2positions(atlas_template,WaveList,varnames[0],varnames[1],lon,lat,pn,NAN,status,strict);
    
    scale_constants(constants,pn,WaveList.n,scaleAtlas);
    
    //if(mgrfile==0) mgrfile=strdup("extraction.mgr");
    if(mgrfile!=0){
      vector<mgr_t> mgr=mgr_create(constants,lon,lat,pn,WaveList);
      mgr_save_ascii(mgrfile, mgr);
      }
    
    delete[]lat;
    delete[]lon;
    
    if(time_template==NULL and F!=0){
      time_template="ELAPSED CALENDAR";
      printf("Default time template : %s\n",time_template);
      }
    
    }
  else{
/*------------------------------------------------------------------------------
    prediction on atlas grid */
    constants=load_atlas(atlas_template,varnames[0],varnames[1],WaveList,&pn,&bv);
    }
  
#if USE_M1_12
#warning USE_M1_12
  wI=WaveList.wave_index("M1");
  
  if(wI>=0){
    int w1,w2;
    w1=WaveList.add(wM1_1);
    w2=WaveList.add(wM1_2);
    STDERR_BASE_LINE("M1_1:%d M1_2:%d\n",w1,w2);
    
    for(pI=0;pI<pn;pI++){
      hconstant_t *cpI=&constants[pI],cp0;
      double aIp,GI;
      
      cp0=*cpI;
      cpI->a=cpI->G=0;
      cpI->init_polar(WaveList.n);
      valcpy(cpI->a,cp0.a,cp0.size);
      valcpy(cpI->G,cp0.G,cp0.size);
      cp0.destroy();
      
      aIp=cpI->a[wI]/wM1.Ap;
      cpI->a[wI]=0.;
      cpI->a[w1]=aIp*wM1_1.Ap;
      cpI->a[w2]=aIp*wM1_2.Ap;
      
      GI=cpI->G[wI];
      cpI->G[w1]=GI;
      cpI->G[w2]=GI;
      }
    }
#endif
  
  outputIsNC= (position_file==0) or times!=0;
  
/*---------------------------------------------------------------------*//**<h1>
  if there are start and end dates, predict with harmonic_prediction()</h1>*/
  if(times==0){
    if(isnad(start) || isnad(final)) TRAP_ERR_EXIT(0,"*** Provide start and end dates for predictions ***\n");
    if(isnan(increment)){
      increment=3600.;
      printf("Default time increment of %g s.\n",increment);
      }
    startd=cnes_time(start,'s');
    finald=cnes_time(final,'s');
    }
  else{
    start=poctime_getdatecnes(times[0   ],'s');
    final=poctime_getdatecnes(times[pn-1],'s');
    }
  init_argument(&astro_angles,start,1);
  startc=sgetdate(start);
  printf("Selected dates from %s to %s.\n",startc,sgetdate(final));fflush(stdout);
  
  buffer=new double[pn];
  
  if(output==NULL){
    if(outputIsNC){
      if(times==0)
        output="predictions.nc";
      else
        output="%s-predictions.nc";
      }
    else
      output="predictions.dat";
    printf("Default output file : %s\n",output);fflush(stdout);
    }
  
  vector<string> token;
  const float mask=1e4f;
  
  if(not outputIsNC){
/*------------------------------------------------------------------------------
    prediction at given locations => ASCII output */
    F=fopen(output,"w");
    if(F==NULL){
      check_error(errno,__LINE__,__FILE__,"error while opening %s",output);
      exit(errno);
      }
    fprintf(F,"#file produced with : "+cmd+"\n#");
    
    token=string_split(time_template, " ");
    for(k=0;k<token.size();k++) {
      if(token[k]=="CNES") {
        fprintf(F,"time(days since 1950-01-01 00:00) ");
        }
      if(token[k]=="ELAPSED") {
        fprintf(F,"time(days since %s) ",startc);
        }
      if(token[k]=="CALENDAR") {
        fprintf(F,"time(human-readable) ");
        }
      }
    
    for(pI=0;pI<pn;pI++){
      fprintf(F," point%d",pI);
      }
    fprintf(F,"\n");
    }
  else if(times==0){
/*------------------------------------------------------------------------------
    gridded NetCDF output */
    poc_global_t global;
    
    poc_dim_t TDim(tv.name,0);
    tv.dimensions << TDim;
    bv.dimensions.insert(0,TDim);
    global << tv;
    
    global << bv;
    status=poc_create(output,global);
    status=poc_def_att(output,poc_att_t("production","constructed around " __LINE_FILE_PACKAGE_REVISION));
    status=poc_def_att(output,poc_att_t("history",cmd));
    }
  else{
/*------------------------------------------------------------------------------
    altimetry */
    
    /* making sure _FillValue will appear as an attribute
        because of python netCDF4 bug... */
    bv << poc_att_t(_FillValue,mask);
    }

  deletep(&startc);
  
/*------------------------------------------------------------------------------
  compute and output tidal prediction */
  t=-increment;
  fI=0;
  fpn=0;
  string outputI;
  
  for(int tI=0;;tI++){
    
    if(times==0){
      t+=increment;
      if(t>finald-startd) break;
      }
    else{
      if(tI>=pn)break;
      t=times[tI]-times[0];
      }
    
    if(F!=0){
      /* ASCII output: time colum(s) */
      for(k=0;k<token.size();k++) {
        if(token[k]=="CNES") {
          fprintf(F,"%.5f ",(t+startd)/d2s);
          }
        if(token[k]=="ELAPSED") {
          fprintf(F,"%.5f ",t/d2s);
          }
        if(token[k]=="CALENDAR") {
          if(startc==0)
            startc=poctime_sdate_cnes(t+startd,'s',' ');
          fprintf(F,"%s ",startc);
          }
        }
      
      deletep(&startc);
      }
    
    /* predictions */
    if(times==0){
      harmonic_prediction(astro_angles,buffer,t,pn,WaveList,constants,nodal_corrections);
      }
    else{
      harmonic_prediction(astro_angles,&buffer[tI],t,1,WaveList,&constants[tI],nodal_corrections);
      }
    
    if(F!=0){
      /* ASCII output: prediction colum(s) */
      for(pI=0;pI<pn;pI++){
        fprintf(F," %12.4g",buffer[pI]);
        }
      fprintf(F,"\n");
      }
    else if(times!=0){
      /* altimetry */
      
      double *bI=&buffer[tI];
      if(isnan(*bI))
        *bI=mask;
      
      if(fpn+1>=pnS[fI]){
        if(strcmp(output,"")!=0){
          string inputName;
          inputName=strrchr0(globbed.gl_pathv[fI],'/');
          inputName.erase(inputName.size()-3);
          
          asprintf(outputI="",output,inputName.c_str());
          }
        else{
          outputI=globbed.gl_pathv[fI];
          }
        
        bv.dimensions[0].len=pnS[fI];
        
        status=access(outputI.c_str(),W_OK);
        
        if(status!=0){
          poc_global_t global;
          global << bv;
          
          printf("[%d/%d:%d]Creating "+outputI+"\n",fI+1,pnS.size(),pnS[fI]);fflush(stdout);
          status=poc_create(outputI,global);
          }
        else{
          printf("[%d/%d:%d]Writing "+outputI+"\n",fI+1,pnS.size(),pnS[fI]);fflush(stdout);
          }
        
        status=poc_update_history(outputI,cmd);
        
        status=poc_put_var(outputI,bv,&buffer[tI-fpn]);
        
        fpn=0;
        fI++;
        }
      else
        fpn++;
      
      }
    else{
      /* gridded NetCDF output */
      for(pI=0;pI<pn;pI++){
        double *bI=&buffer[pI];
        if(isnan(*bI))
          *bI=NC_FILL_DOUBLE;
        }
      
      printf("outputing %s\r",output);fflush(stdout);
      
      status=poc_put_vara(output,tv,tI,&t);
      status=poc_put_vara(output,bv,tI,buffer);
      }
    
    }
  
  if(F!=0)
    fclose(F);
  delete[] buffer;
  deletep(&times);
  if(globbed.gl_pathv!=0)
    globfree(&globbed);
  
/*---------------------------------------------------------------------*//**<h1>
  and SpringNeap_cycle()</h1>*/
  if(SpringNeap)
    SpringNeap_cycle(astro_angles, WaveList, constants, nodal_corrections,startd,finald);
  
  STDOUT_BASE_LINE(" ^^^^^^^^^^^^^ end of computation ^^^^^^^^^^^^^\n");
  exit(0);
}
