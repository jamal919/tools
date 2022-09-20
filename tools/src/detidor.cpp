
/*******************************************************************************

  T-UGO tools, 2006-2014

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

<!-- USE FIND AND REPLACE TO EDIT THIS LINE, SO THAT print_help IS ALSO UPDATED -->
\brief Detide a time serie.

<!-- A LINK TO main() or print_help() WILL NOT LINK TO THE RIGHT SOURCE ! -->
See the main function for how this works
and the print_help function for how to use this.
*/
/*----------------------------------------------------------------------------*/

#define MAIN_SOURCE

#include "version-macros.def" //for VERSION and REVISION


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#include "tools-structures.h"

#include "spectrum.h"
#include "tides.def"
#include "fe.h"
#include "archive.h"
#include "poc-time.h"
#include "mgr-converter.h"
#include "mgr.h"
#include "functions.h"
#include "statistic.h"

extern hconstant_t * atlas_constantsN(const char *atlas_convention,const spectrum_t &s, double *x, double *y, int nb_pos);


//   tseries_t reduce(tseries_t source, date_t start, date_t end) {
//     date_t reference=date_t(1950,1,1,0.);
//     int i, first, last, count;
//
//     double t1 = julian_day(start)-julian_day(reference);
//     double t2 = julian_day(end)  -julian_day(reference);
//
//     count=0;
//     i=0;
//     while(source.t[i]<t1) {
//       if(i==source.n-1) break;
//       i++;
//       }
//     first=i;
//     while(source.t[i]<t2) {
//       if(i==source.n-1) break;
//       i++;
//       count++;
//       }
//     last=i;
//
//     tseries_t *reduced= new tseries_t(count,source.nparam);
//
//     for(i=first;i<last;i++) {
//       reduced->t[i-first]=source.t[i];
//       for (size_t k=0;k<source.nparam;k++) {
//         reduced->x[k][i-first]=source.x[k][i];
//         }
//       }
//     return(*reduced);
//   }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

spectrum_t LPatlas_define()

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// define long period spectrum
/*----------------------------------------------------------------------------*/
{
  spectrum_t atlas;
  int i = 0;
  atlas.nmax = 14;


  atlas.n = 10;
  atlas.waves = new tidal_wave[atlas.n];
  atlas.waves[i++] = wMm;   // FES2004 atlas
  atlas.waves[i++] = wMf;   // FES2004 atlas
  atlas.waves[i++] = wMtm;  // FES2004 atlas
  atlas.waves[i++] = wMSqm; // FES2004 atlas
  atlas.waves[i++] = wMSm;  // spline admittance
  atlas.waves[i++] = wMSf;  // spline admittance
  atlas.waves[i++] = wMStm; // spline admittance
  atlas.waves[i++] = wMqm;  // spline admittance
  atlas.waves[i++] = wSa;   // equilibrium (hard-coded in atlas_constants)
  atlas.waves[i++] = wSsa;  // equilibrium (hard-coded in atlas_constants)

  if (i>atlas.n) TRAP_ERR_EXIT(-1,"Oups, please increase number of waves to be allocated (set  atlas.n)!\n");
  initialize_omega(&atlas);
  
  return(atlas);
}



extern toponyms_t timeseries_formats;


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void print_help(const char *prog_name)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/** \brief prints help of this programme.

\param prog_name name of the program to be printed by this help function : as[0]
*/
/*----------------------------------------------------------------------------*/
{
/** The code of the body of this function is :
    \code /**/ // COMPILED CODE BELOW !!!
  mgr_init_formats();
  printf("\n"
    "NAME AND VERSION\n"
    "  %s version " VERSION " " REVISION "\n", prog_name);
  printf("\n"
    "USE\n"
    "  %s FILE [OPTIONS]\n",prog_name);
  printf("\n"
    "DESCRIPTION\n"
/* USE FIND AND REPLACE TO EDIT THE LINE BELOW, SO THAT THE FILE HEAD IS ALSO UPDATED */
    "  Detide a time serie.\n"
    "  The input file must be in ascii format.\n"
    "\n"
    "OPTIONS\n"
    "Current options:\n"
    "  -h,--help : show this help and exit\n"
    "   -format followed by time serie format within (default: GLOSS, that works also for JASL)\n"
    "     "+get_key_list(timeseries_formats)+"\n"
    "   --residual-format followed by time serie format within (most of) those above (default: GNU)\n"
    "   -context <fomat>  :ctoh to automatically organise the outputs\n"
    "   -n   <ncols>      : number of columns in ASCII GNU format input file (default: 1)\n"
    "   -s : followed by the start date. See DATE FORMATS below\n"
    "   -f : followed by the end date. See DATE FORMATS below\n"
    "   -step <n>y or <n>m: time series spliting per n years or n month (default: all at once)\n"
    "   --mooring=<info>  : input info on mooring, eg: \"lon=<lon> lat=<lat> code=<code> name=<name depth=<depth>\"\n"
    "   -spectrum followed by a predefined spectrum (see SPECTRA below) or a file. Default is COASTAL.\n"
    "   -tides <files>    : convention of the tides atlas         MANDATORY TO ACTIVATE DETIDING\n"
    "   -cm               : use if atlas is in cm\n"
    "   -o  <out_prefix>  : outputs\n"
    "   --no-date         : disable dates in output file names\n"
    "   -w : followed by the list of waves to predict for\n"
    "\n"
    "Backward compatibility option:\n"
    "  --discrete_parsing : use old/obsolete discrete parsing\n"
    "\n"
    "Old syntax options:\n"
    "   -spectrumtype,-spectrumfile,-file are all equivalents of -spectrum\n"
    "   -c   <atlas_type> : specify the atlas file name convention : <prefix>wave<postfix> \n"
    "   -l   <stations_list> : save station list to <stations_list>\n"
    "\n"
    "Not yet documented options:\n"
    "  --time : \n"
    "  -averaging : \n"
    "  -maxratio : \n"
    "  -repetition : \n"
    "  --center : \n"
    "  --debug : \n"
    "  --HF : \n"
    "\n"
    "SPECTRA\n"
    "  Possible values are ESTUARINE, ESTUARINE-HF, COASTAL, COASTAL-HF, SHELF, SHELF-HF, DEEP, DEEP-HF.\n"
    "  Spectra suffixed with `-HF' are without long-period tides.\n"
    "  Use showarg to know what waves are included in a spectrum.\n"
    "\n"
    "EXAMPLES\n"
    "For comparison with altimetry:\n"
    "  %s --no-date -tides /data/soa/maree/fes2014/WAVE.FES2014b.nc -spectrum COASTAL-HF ...\n",prog_name
    );
  print_poctime_scan_date_help(0);
  /** \endcode */
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
#define NB_MAX_INPUTS 300

  char sub_dir[256]="";
  char *stations_list_path=NULL;

  int i,k,nitems,ncol=1;
  int n,status;
  FILE *foutname;
  char *inputs[NB_MAX_INPUTS]={NULL};
  int nb_inputs=0;
  const char *keyword,*s;
  const char *atlas_convention=NULL,*mooring_info=NULL;
  char *time_template=NULL;
  char *spectrum_type=NULL;
  char *onde[10];
  const char *format=FORMAT_NAME_GLOSS,*residual_format=FORMAT_NAME_GNU;
  string time_format, delimiter;
  const char *context="";
  int nonde=0;
  const double *time;
  double **z, mean,mask=1.e+10;
  date_t step=null_date;
  spectrum_t AnalysisList, PredictionList, solved;
  int detiding=0;
  hconstant_t prior;
  double scaleAtlas=1.;
  int nodal=1;
  double x=999999, y=999999, *tides;
  statistic_t stats_input, stats_res;
  date_t start(0,0,0,0.), final(0,0,0,0.);
  date_t start_req(0,0,0,0.), final_req(0,0,0,0.);
  date_t start_file_name(0,0,0,0.), final_file_name(0,0,0,0.);
  date_t one_second(0,0,0,1.);
  int d_date;
  int nb_data;
  string output;
  char  *file_out_prefix=NULL;
  char spectrum_name[256];
  char atlas_name[256]="\0";
  int optimize_atlas_detiding=1;
  bool debug=false, datedNames=true;
  bool discret_parsing=false, HF=false;
  double repetition=0.0;
  double maxratio=0.3;
  char *center=0;
  double averaging_time=0.0;
  mooring_t *mooring;
  mooring_t *moorings;
  vector<mgr_t> *mgr=0;
  
  char mgrfile[1024];
  int nmgr=1;
  char residualsfile[1024];
  char outname[1024];
  string path, basename;
  mgr_data_t **tmp;
  int numf;


  if (argc==1){
    print_help(argv[0]);
    wexit(0);
    };

  fct_echo( argc, argv);

  n=1;
  while (n < argc) {
    keyword=argv[n];
    switch (keyword[0]) {
      case '-':
        if( strcmp(keyword,"-h")==0 or
            strcmp(keyword,"--help")==0 ) {
          print_help(argv[0]);
          wexit(0);
          }
        if(strcmp(keyword,"-spectrum")==0 ||
           /* old */
          strcmp(keyword,"-spectrumtype")==0 || strcmp(keyword,"-file")==0 || strcmp(keyword,"-spectrumfile")==0) {
          spectrum_type= strdup(argv[n+1]);
          n++;
          n++;
          continue;
          }
        if(strcmp(keyword,"--time")==0) {
          time_template=strdup(argv[n+1]);
          n++;
          n++;
          continue;
          }
        if(strcmp(keyword,"-averaging")==0) {
          nitems= sscanf(argv[n+1],"%lf",&averaging_time);
          n++;
          n++;
          continue;
          }
         if(strcmp(keyword,"-maxratio")==0) {
          nitems= sscanf(argv[n+1],"%lf",&maxratio);
          n++;
          n++;
          continue;
          }
       if(strcmp(keyword,"-repetition")==0) {
          nitems= sscanf(argv[n+1],"%lf",&repetition);
          n++;
          n++;
          continue;
          }
        if(strcmp(keyword,"--center")==0) {
          center=strdup(argv[n+1]);
          n++;
          n++;
          continue;
          }
        if(strcmp(keyword,"--discrete_parsing")==0) {
          discret_parsing=true;
          n++;
          continue;
          }
        if(strcmp(keyword,"-format")==0) {
          format= argv[n+1];
          n++;
          n++;
          continue;
          }
        if(strcmp(keyword,"--residual-format")==0) {
          residual_format= argv[n+1];
          n++;
          n++;
          continue;
          }
        if(strcmp(keyword,"-context")==0) {
          context= argv[n+1];
          n++;
          n++;
          continue;
          }
        if(strncmp(keyword,"--mooring=",10)==0){
          mooring_info=strdup(&(argv[n][10]));
          n++;
          continue;
          }
        if(strcmp(keyword,"--debug")==0) {
          debug=true;
          n++;
          continue;
          }
        if(strcmp(keyword,"--no-date")==0) {
          datedNames=false;
          n++;
          continue;
          }
        if(strcmp(keyword,"--HF")==0) {
          HF=true;
          n++;
          continue;
          }
        if(strcmp(keyword,"-tides")==0) {
          atlas_convention = argv[n+1];
          detiding=1;
          n++;
          n++;
          continue;
          }
        if(strcmp(keyword,"-cm")==0) {
          scaleAtlas=0.01;
          n++;
          continue;
          }
        if(strncmp(keyword,"-step",5)==0){
            int nb_step;
            char unit;
            sscanf(argv[n+1],"%d%c",&nb_step,&unit);
            if (unit=='y') step.year = nb_step;
            else if (unit=='m') step.month = nb_step;
            else if (unit=='d') step.day = nb_step;
            else {
              printf("bad syntax %s for option -step\n",argv[n+1]);
              print_help(argv[0]);
              wexit(-1);
              }
            n++;
            n++;
            continue;
            }
        if (strlen(keyword) > 2) {
          printf("Unkown option %s\n", keyword);
          print_help(argv[0]);
          wexit(-1);
          }
        switch (keyword[1]) {

          case 'n' :
            sscanf(argv[n+1],"%d",&ncol);
            if (ncol>2) {
              printf("Sorry, for now cannot consider more than 2 columns data (n=%s)\n",argv[n+1]);
              print_help(argv[0]);
              wexit(-1);
              }
            n++;
            n++;
            break;

          case 's' :
            s=argv[n+1];
            n++;
            n++;
            status=poctime_scan_date(s,&start_req,0);
            break;

          case 'f' :
          case 'e' :
            s=argv[n+1];
            n++;
            n++;
            status=poctime_scan_date(s,0,&final_req);
            break;


            /* OLD way to get tides atlas */
            /* convention */
          case 'c' :
            atlas_convention= strdup(argv[n+1]);
            n++;
            n++;
            detiding=1;
            break;

          case 'w' :
            i=1;
            while(n+i<argc and argv[n+i][0]!='-')  {
              onde[i-1]= strdup(argv[n+i]);
              i++;
              }
            n=n+i;
            nonde=i-1;
            onde[nonde]=0;
            break;
         
         case 'o' :
           file_out_prefix = strdup(argv[n+1]);
           n++;
           n++;
           break;

          case 'l' :
            stations_list_path= strdup(argv[n+1]);
            n++;
            n++;
            break;
          
          default:
            printf("unknown option %s\n",keyword);
            print_help(argv[0]);
            wexit(-1);
          }
        break;

      default:
        if (nb_inputs>=NB_MAX_INPUTS) TRAP_ERR_EXIT(-1,"Sorry too many input files (max %d)\n",NB_MAX_INPUTS);
        inputs[nb_inputs]  = strdup(argv[n]);
        nb_inputs++;
        n++;
        break;
      }
    
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  LOAD SPECTRUMS USED FOR ANALYSIS
      
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if(nonde>0){
    AnalysisList.init(initialize_tide(),onde);
    sprintf(spectrum_name, "custum%02d", AnalysisList.n);
    }
  else{
    AnalysisList.waves=new tidal_wave[100];
    AnalysisList.nmax=100;
    AnalysisList.n = 0;
    
/* -------------------------------------------------------------------------------------
    predifined spectrum */
    if(spectrum_type==0) spectrum_type=strdup("COASTAL");
    
    for (i=0;i<SPECTRUM_DEFAULT_NB;i++) {
      if (strcmp(spectrum_type, spectrum_list_default[i])==0) {
        bool Z0=false;
        AnalysisList = spectrum_init_ref(spectrum_type, Z0);
        /* build spectrum name */
        sprintf(spectrum_name, "%s_%02d", spectrum_type, AnalysisList.n);
        break;
        }
      }
    
/* -------------------------------------------------------------------------------------
  try by file */
    if (i==SPECTRUM_DEFAULT_NB) {
      char *c;
      AnalysisList = spectrum_init_from_file(spectrum_type,1);
/*--------------------------------------------------------------------------------------
      build spectrum name */
      if (c=strrchr(spectrum_type,'/'))
        strcpy(spectrum_name, c+1);
      else
        strcpy(spectrum_name, spectrum_type);
      if (c=strrchr(spectrum_name,'.'))
        *c = '\0';
      }
    }
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  ATLAS DETIDING : 
  
  try to get all mooring positions to read the atlas only once.
  If not possible, the atlas detiding will be done station per station while 
  reading their data
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  hconstant_t *priorN;

  if(optimize_atlas_detiding==1 && detiding==1) {
    double *xN = new double[nb_inputs];
    double *yN = new double[nb_inputs];
    mooring_t mooringp;

/*------------------------------------------------------------------------------
    get all positions of the stations from original file */
    for (numf=0; numf<nb_inputs; numf++) {
//      in=fopen(inputs[numf],"r");
//      if(in==0) {
//        STDERR_BASE_LINE("Probleme a l ouverture du fichier : %s \n",inputs[numf]);
//        TRAP_ERR_EXIT(-1,"exiting\n");
//      }
//      fgets(line,256,in);
//      hawaii_decode_position_p(line, &mooringp);
//      fclose(in);
      const char *header_format=0;
      status=mgr_decode_position(inputs[numf], format, header_format, &mooringp);
      if(status==-1) {
        optimize_atlas_detiding=0;
        STDERR_BASE_LINE("info: failed get all station positions to optimize atlas detiding with format %s. Will be done on the fly", format);
        break;
        }
      xN[numf] = mooringp.lon;
      yN[numf] = mooringp.lat;
      }
    STDERR_BASE_LINE("load model for detiding : %s\n",atlas_convention);
    PredictionList=LPatlas_define();
    priorN=atlas_constantsN(atlas_convention, PredictionList,  xN,  yN, nb_inputs);
    }

/*------------------------------------------------------------------------------
  PREPARE OUTPUTS */
  char postfix[256]="\0";

/* -----------------------------------------------------------------------------
  common file postfix */
  sprintf (postfix, "%s", spectrum_name);
  if(detiding==1) sprintf (postfix, "%s.%s", postfix, atlas_name);
  
/* -----------------------------------------------------------------------------
  build sub-dir */
  if (strncmp(context, "ctoh", 4) == 0) {
    if (file_out_prefix==NULL) sprintf(sub_dir, "../");
    sprintf(sub_dir, "%sanalysis.%s", sub_dir, postfix);
    if (step.year)  sprintf (sub_dir, "%s.%02dy", sub_dir, step.year);
    if (step.month) sprintf (sub_dir, "%s.%02dm", sub_dir, step.month);
    sprintf (sub_dir, "%s/", sub_dir);
    }


/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  loop over moorings (i.e. input files)
      
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  tseries_t **series, residuals;
  int nb_series;

  int nb_data_max=0;
  int ncol_max=0;
  int reconstruct=0;

  moorings = new  mooring_t[nb_inputs];

  for (numf=0; numf<nb_inputs; numf++) {
/* ----------------------------------------------------------------------------
    load time serie */
    printf("#################################################################\n");
    printf("load timeseries file : %s\n",inputs[numf]);
    nb_series = mgr_converter(inputs[numf], format, delimiter, time_format, ncol,
                               NULL, /* header_format  */
                               mooring_info,
                               start_req, final_req, step, reconstruct,0.,
                               NULL, // file_out_prefix, /* file_out */
                               NULL, // "GNU", /* format_out,  */
                               &series);
    if(!nb_series) {
      STDERR_BASE_LINE("could not load properly time serie %s. Ignored.\n", inputs[numf]);
      continue;
      }
 /* -------------------------------------------------------------------------------------
    remove masked data */
    nb_series = mgr_remove_masked_data(series, nb_series);

    x = series[0]->mooring.lon;
    y = series[0]->mooring.lat;
    
/* --------------------------------------------------------------------------------
    LOAD MODEL Would be better not to load it for each mooring  */
    if(detiding==1) {
      if(optimize_atlas_detiding==0) {
        printf("load model for detiding : %s\n",atlas_convention);
        PredictionList=LPatlas_define();
        prior=atlas_constants(atlas_convention, PredictionList,  x,  y);
        }
      else {
        prior = priorN[numf];
        }
      scale_constants(&prior,1,PredictionList.n,scaleAtlas);
      const char *dashes="----------------------------------";
      printf("%10s |%8s |%9s\n","wave","amp(mm)","pha(deg)");
      printf("%.11s+%.9s+%.9s\n",dashes,dashes,dashes);
      for(i=0;i<PredictionList.n;i++)
        printf("%10s |%7.1f  |%8.1f\n",PredictionList.waves[i].name,prior.a[i]*1e3,prior.G[i]);
      }

/* -------------------------------------------------------------------------------------
    prepare  output; local directory if not 'ctoh' context */
    if (strncmp(context, "ctoh", 4) != 0 && file_out_prefix==NULL)  file_out_prefix=strdup("./");
    mgr_build_out_prefix(inputs[numf],file_out_prefix==NULL?"":file_out_prefix, &path, &basename);

/* -------------------------------------------------------------------------------------
    filename (without postfix) */
    output = path + sub_dir + basename;
    if (numf==0) {
      printf ("output=%s\n", output.c_str());
      status=mkdir((path + sub_dir).c_str(),ACCESSPERMS);
      if(status!=0 and errno!=EEXIST) {
        status=errno;
        TRAP_ERR_EXIT(status,"error creating %s with mode %04o (%d %s)",(path + sub_dir).c_str(),ACCESSPERMS,status,strerror(status));
        }
      }
    
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
    allocate data arrays
      
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

    /* loop over series */
    tseries_t *serie;
    int num_s;
    int data_max=0;

/* ----------------------------------------------------------------------
    get data size */
    for (num_s=0; num_s<nb_series; num_s++)
      if (series[num_s]->n > data_max)
        data_max = series[num_s]->n;
    
/* ----------------------------------------------------------------------
    free if data too small */
    if ( nb_data_max != 0 and nb_data_max < data_max ) {
      nb_data_max = 0;
      delete[] tides;
//       delete[] residuals;
      }

/* ----------------------------------------------------------------------
    allocate intermediate data */
    if (nb_data_max == 0) {
      nb_data_max = data_max;
      tides=new double[nb_data_max];
//       residuals=new double[nb_data_max];
      }

/* ----------------------------------------------------------------------
    allocate output arrays */
    if (ncol_max==0) {
/* ----------------------------------------------------------------------
      count nb columns */
      for (num_s=0; num_s<nb_series; num_s++)
        updatemax(&ncol_max,series[num_s]->nparam);
/* ----------------------------------------------------------------------
      alloc data */
      tmp=new mgr_data_t *[ncol_max];
      mgr=new vector<mgr_t> [ncol_max];
      }
      

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
    loop over timeseries
      
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

    for (num_s=0; num_s<nb_series; num_s++) {
      serie = series[num_s];
      mooring = &serie->mooring;
      STDERR_BASE_LINE("assigning moorings[%d](.name=%p)\n",numf,mooring->name);
      moorings[numf] = *mooring;
      nb_data = serie->n;
      ncol = serie->nparam;
      time = serie->t;
      z = serie->x;
      mask = serie->mask;
      if (nb_data<=1) {
        STDOUT_BASE_LINE("No data for %s %04d%02d%02d %04d%02d%02d\n",
                            mooring->name,
                            start_req.year, start_req.month, start_req.day,
                            final_req.year, final_req.month, final_req.day);
        fflush(stdout);
        continue;
        }

/* -------------------------------------------------------------------------------------
      detiding model */
      if(detiding==1) {
        range_t<double> tidalRange;
        printf("model pre-detiding from %s : ",atlas_convention);
        status=harmonic_predictionTS(tides, time, nb_data, PredictionList, prior.a, prior.G, nodal);
        for (k=0;k<nb_data;k++) {
          double *z0k=&z[0][k];
          if(*z0k==serie->mask) continue;
          const double tidesk=tides[k];
          *z0k-=tidesk;
          tidalRange<<tidesk;
          }
        printf("[%.3f;%.3f]m\n",tidalRange.min,tidalRange.max);
        }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
      harmonic analysis for each serie 
      
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

      residuals.allocate(nb_data_max, ncol);

      printf("performs harmonic analysis : %d data\n",nb_data);
      for(k=0;k<ncol;k++) {
        if(!discret_parsing)
          tmp[k]= harmonic_analysis_with_parsing(serie->x[k], serie->mask, time, residuals.x[k], &mean, nb_data, AnalysisList, solved, nodal, stdout);
        else {
          tmp[k]= harmonic_analysis(serie->x[k], residuals.x[k], time, nb_data, AnalysisList, averaging_time, repetition, solved, nodal, maxratio, stdout);
          }
/* -------------------------------------------------------------------------------------
        recomposition */
//         if(detiding==1) {
//           for(w=0;w<PredictionList.n;w++) {
//             l=solved.wave_index(PredictionList.waves[w].name);
//             if(l!=-1) {
//               complex<double> z=polar((double) prior.a[w], prior.G[w]*M_PI/180.)+polar(tmp[k][l].amp,tmp[k][l].phi*M_PI/180.);
//               }
//             }
//           }
        }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
      fill-in mgr data structures
      
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

      start=poctime_getdatecnes(time[0],   'd');
      final=poctime_getdatecnes(time[nb_data-1], 'd');
    
      for(k=0;k<ncol;k++) {
/* -------------------------------------------------------------------------------------
        in case where harmonic analysis failed */
        if(tmp[k]==0) continue;
        for(i=0;i<nmgr;i++) {
          mgr_t gauge;
          gauge.number=i;
          if(mooring->name!=0)
            strcpy(gauge.name,mooring->name);
          else
            strcpy(gauge.name,"no_name");
          gauge.number = mooring->code;
          if(center==0) center=strdup("undocumented");
          strcpy(gauge.origine,center);
          gauge.data=tmp[k];
          gauge.nwave=solved.n;
          gauge.loc.units=strdup("degrees");
          gauge.loc.lon=mooring->lon;
          gauge.loc.lat=mooring->lat;
          gauge.loc.depth=mooring->depth;
          gauge.mindex=1;
          gauge.duree=time[nb_data-1]-time[0];
          sprintf(gauge.debut,"%2.2d/%2.2d/%4d",start.day,start.month,start.year);
          sprintf(gauge.fin,  "%2.2d/%2.2d/%4d",final.day,final.month,final.year);
          mgr[k].push_back(gauge);
          }
        }
      
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
      harmonic constants (*.mgr) archiving 
      
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

      start_file_name = start;
      final_file_name = final;
      char *output_dated=new char[1024];

      /* file names dates start always at begin of the month and of the year */
      if (step.month != 0) {
        start_file_name.day = 1;
        if (start_req.month != 0 and step.month != 1) {
          /* file names dates are modulo of step.month */
          start_file_name.month = start_req.month + floor((start.month - start_req.month)/step.month) * step.month;
          }
        }
      if (step.year != 0) {
        start_file_name.day = 1;
        start_file_name.month = 1;
        if (start_req.year != 0 and step.year != 1) {
          /* file names dates are modulo of step.year */
          start_file_name.year = start_req.year + floor((start.year - start_req.year)/step.year) * step.year;
          }
        }
      
      /* if some seconds, next day started also */
      if (final_file_name.second != 0) {
        final_file_name.second = -1;
        final_file_name.day += 1;
        d_date = cnes_time(final_file_name, 'd');
        calendary(d_date, &final_file_name);
        }
#if DEBUG
      printf("start_real %04d%02d%02d\n", start.year, start.month, start.day);
      printf("start_file %04d%02d%02d\n", start_file_name.year, start_file_name.month, start_file_name.day);
      printf("final_file %04d%02d%02d:%g\n", final_file_name.year, final_file_name.month, final_file_name.day,final_file_name.second);
#endif
      
      if(datedNames){
        if(mooring->name!=0) sprintf(output_dated, "%s_%s", output.c_str(), mooring->name);
        sprintf(output_dated, "%s_%04d-%02d-%02d_%04d-%02d-%02d",
                output.c_str(),
                start_file_name.year, start_file_name.month, start_file_name.day,
                final_file_name.year, final_file_name.month, final_file_name.day);
        sprintf(output_dated, "%s.%s", output_dated, postfix);
        }
      else{
        sprintf(output_dated, "%s", output.c_str());
        }
        
      sprintf(mgrfile, "%s.mgr", output_dated);

/*------------------------------------------------------------------------------
      put harmonic constant in increasing frequency order */
      status=mgr_WaveOrder(mgr[0], mgr[0].size(), (string)"FREQUENCY");
        
      printf("save harmonic constant file : %s\n",mgrfile);
      switch (ncol) {
        case 1:
          status=mgr_save_ascii(mgrfile, mgr[0]);
          break;
          
        case 2:
          status=mgr_save2D(mgrfile, mgr[0].size(), mgr, 10.);
          break;

        default:
          for(int k=0;k<ncol;k++) {
            sprintf(mgrfile, "%s.%2.2d.mgr", output_dated, k);
            status=mgr_save_ascii(mgrfile, mgr[k]);
            }
          break;
          
        }
    
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
      timeseries residuals (*.res) archiving 
      
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

//       printf("save residuals %d timeseries file : %s (original, detided, tides)\n",nb_data,residualsfile);
      tseries_t series(3);
      for(int k=0;k<ncol;k++) {
        switch (ncol) {
          case 1:
            sprintf(residualsfile, "%s", output_dated);
            break;
          
          default:
            sprintf(residualsfile, "%s.%2.2d", output_dated, k);
            break;
          }
     
        series.x[0]=serie->x[k];
        series.x[1]=residuals.x[k];
        series.n=nb_data;
        series.mask=mask;
        series.t=serie->t;
        series.x[2]=new double[series.n];
        series.origin=date_t(1950,1,1,0.);
        series.time_unit='d';
        for(int k=0;k<series.n;k++) {
          if(series.x[0][k]!=series.mask and series.x[1][k]!=series.mask) series.x[2][k]=series.x[0][k]-series.x[1][k];
          else series.x[2][k]=series.mask;
          }
        if(time_template==0) time_template=strdup("");
        printf("save residuals %d timeseries file : %s (original, detided, tides)\n",nb_data,residualsfile);
        mgr_save_timeserie(residualsfile, *mooring, series, 'm', residual_format, false, time_template);
      
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
        variance reduction statistics archiving 
      
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

        switch (ncol) {
          case 1:
            sprintf(outname, "%s.mgr.stat", output_dated, k);
            break;
          
          default:
            sprintf(outname, "%s.%2.2d.mgr.stat", output_dated, k);
            break;
          }
          
        
        printf("save stats file : %s\n",outname);
        foutname=fopen(outname,"w");
        fprintf(foutname,"# mooring: %s\t%lf\t%lf\n",mooring->name, mooring->lon, mooring->lat);
        stats_input = get_statistics(&z[k][0],  mask, nb_data);
        fprintf(foutname,"#nb_samples\t mean\t variance\t deviation\n");
        fprintf(foutname,"%6d\t %8.4lf\t  %8.4lf\t  %8.4lf\n",nb_data, stats_input.mean,stats_input.std*stats_input.std,stats_input.std);
        stats_res= get_statistics(residuals.x[k],  mask, nb_data);
        fprintf(foutname,"%6d\t %8.4lf\t  %8.4lf\t  %8.4lf\n",nb_data, stats_res.mean,stats_res.std*stats_res.std,stats_res.std);
        fclose (foutname);
        delete[] series.x[2];

        }
        
      (*serie).destroy();
      residuals.destroy();
        
      
      fflush(stdout);
      } /* for each time serie in given mooring */

    } /* for each mooring */

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  record moorings list
  
  To get list of stations, use:  scripts/tg_display_stats.py *.stat
      
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if(stations_list_path!=0){
    STDOUT_BASE_LINE("Save mooring list in %s\n", stations_list_path);
    mgr_save_moorings(stations_list_path, numf, moorings);
    }

  STDOUT_BASE_LINE(" ^^^^^^^^^^^^^ end of computation ^^^^^^^^^^^^^\n");
  exit(0);
}
