
/**************************************************************************

  T-UGO tools, 2006-2011

  Unstructured Ocean Grid initiative

***************************************************************************/
/** \file

\author  Florent Lyard      LEGOS/CNRS, Toulouse, France. florent.lyard@legos.obs-mip.fr
\author  Yoann Le Bars      LEGOS, Toulouse, France (PhD)
\author  Laurent Roblou     LEGOS/CNRS, Toulouse, France
\author  Damien Allain      LEGOS/CNRS, Toulouse, France
\author  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
\author  Frédéric Dupont    Université de Laval à Québec, Canada

<!-- USE FIND AND REPLACE TO EDIT THIS LINE, SO THAT print_help IS ALSO UPDATED -->
\brief Predict maximum current.

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
#include <time.h>
#include <unistd.h>

#include <stdint.h> //uint64_t
#include <sys/time.h> //gettimeofday

#include "tools-structures.h"

#include "tides.h"
#include "functions.h" //safely includes omp.h
#include "poc-netcdf-data.hpp"
#include "fe.h"
#include "poc-time.h"
#include "matrix.h"
#include "maths.h"


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
    "  %s -a atlas_convention -v U_amplitude U_phase [V_amplitude V_phase [W_amplitude W_phase]] -w wave1 [wave2 ...]\n",prog_name);
  printf("\n"
    "DESCRIPTION\n"
/* USE FIND AND REPLACE TO EDIT THE LINE BELOW, SO THAT THE FILE HEAD IS ALSO UPDATED */
    "  Predict maximum current.\n"
    "\n"
    "OPTIONS :\n"
    "  -h,--help  Show this help and exit.\n"
    "  -a  followed by the atlas file name convention. See below.\n"
    "  -s  followed by the start date in dd/mm/yyyy format\n"
    "  -f  followed by the end date in dd/mm/yyyy format\n"
    "  -i  followed by the date increment concatenated with the unit: s (default), m, h or d. Default increment: 3600s=60m=1h.\n"
    "  -v  followed by the variables names for the amplitude and the phase.\n"
    "  -w  followed by the list of waves to predict for\n");
  print_tide_decode_atlasname_help();
  print_OPENMP_help(prog_name); /** \endcode */
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
  int i,j,k,l,m;//<wave or time, variable, node, constants and process indexes
  int nn,nv;//number of NetCDF variable names (always even), number of complex variables
  int n,status;//<argument index or number of nodes, NetCDF status
  char *keyword;
  char *waveNames[100];//list of wave names
  char **varnames=NULL;
  
  char *aPC=NULL;//<atlas path convention
  char *aP,*aN,*gN;//<atlas,amplitude,phase names
  poc_var_t var;
  mesh_t mesh;
  
  char *s,*unit;
  date_t start=NADate,final=NADate;
  double startd,finald;
  double increment=NAN;
  
  spectrum_t WaveList;
  astro_angles_t astro_angles;
  
  complex<double> *cbuf=NULL;
  double **predbufs=NULL,**maxbufs=NULL;
  int nproc;
  hconstant_t *constants;
  int nc;//number of constants =n*nv
  
  string cmd=fct_echo( argc, argv);
  
  i=0;
  //waveNames[i++]=strdup("Z0");
  waveNames[i]=NULL;

  n=1;
  while (n < argc) {
    keyword=strdup(argv[n]);
    if(strcmp(keyword,"--help")==0) {
      print_help(argv[0]);
      exit(0);
      }
    switch (keyword[0]) {
      case '-':
        switch (keyword[1]) {

         case 'a' :
          aPC=strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 's' :
            s= strdup(argv[n+1]);
            n++;
            n++;
            status=sscanf(s,"%d/%d/%d",&start.day,&start.month,&start.year);
            if(status==3)
              start.second=0.;
            free(s);
            break;

        case 'f' :
            s= strdup(argv[n+1]);
            n++;
            n++;
            status=sscanf(s,"%d/%d/%d",&final.day,&final.month,&final.year);
            if(status==3)
              final.second=0.;
            free(s);
            break;

        case 'i' :{
          s=strdup(argv[n+1]);
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
              exit(EINVAL);
            }
          free(s);
          break;}

        case 'v' :
          for(j=1;n+j<argc && argv[n+j][0]!='-';j++);
          nn=j-1;
          if(varnames){
            __FILE_LINE__(stdout,"multiple use of -v : the last one overrides the previous.");
            for(j=0;varnames[j]!=NULL;j++){
              free(varnames[j]);
              }
            delete[] varnames;
            }
          varnames=new char*[nn+1];
          for(j=0;j<nn;j++){
            varnames[j]=strdup(argv[n+1+j]);
            }
          varnames[nn]=NULL;
          n++;
          n+=nn;
          break;

        case 'w' :
          n++;
          i=pos((char*)NULL,waveNames,100);
          for(;n<argc;i++,n++) {
            waveNames[i]=strdup(argv[n]);
            }
          waveNames[i]=NULL;
          break;

        case 'h' :
          print_help(argv[0]);
          exit(0);
          break;

        default:
          printf("unknown option %s\n",keyword);
          print_help(argv[0]);
          exit(EINVAL);
        }
        break;

      default:
        printf("unknown option %s\n",keyword);
        print_help(argv[0]);
        exit(EINVAL);
        break;
      }
      free(keyword);
    }

/*-----------------------------------------------------------------------------
  list of waves */
  if(waveNames[0]==NULL){
    fprintf(stderr,"*** Please provide a list of waves. ***\n");
    print_help(argv[0]);
    exit(EINVAL);
    }
  WaveList.init(initialize_tide(),waveNames);
  printf("# number of waves : %d \n",WaveList.n);
  
/*-----------------------------------------------------------------------------
  time boudaries and increment */
  if(isnad(start)){
    start.init(2000);
    keyword=sgetdate(start);
    __ERR_BASE_LINE__("Using default start time of %s\n",keyword);
    delete[]keyword;
    }
  
  if(final<start){
    keyword=sgetdate(start);
    __ERR_BASE_LINE__("Warning : %s after start time!\n",keyword);
    delete[]keyword;
    final=NADate;
    }
  
  if(isnad(final)){
    final.init(2018);
    keyword=sgetdate(final);
    __ERR_BASE_LINE__("Using default final time of %s\n",keyword);
    delete[]keyword;
    }
  
  startd=cnes_time(start,'s');
  finald=cnes_time(final,'s');
  
  if(isnan(increment)){
    increment=1800.;
    __ERR_BASE_LINE__("Using default time increment of %gs\n",increment);
    }
  
  const int nt=(int)floor((finald-startd)/increment)+1;
  __ERR_BASE_LINE__("%d time steps : [%g;%g]\n",nt,startd,finald);
  
/*-----------------------------------------------------------------------------
  list of variables */
  if(varnames==NULL){
    fprintf(stderr,"*** Please provide a list of variables. ***\n");
    print_help(argv[0]);
    exit(EINVAL);
    }
    
  printf("# number of variables : %d \n",nn);
  if(nn&1 || nn>6){
    fprintf(stderr,"*** There must be 2, 4 or 6 variables, %d given. ***\n",nn);
    print_help(argv[0]);
    exit(EINVAL);
    }
  nv=nn/2;
  
  for(i=0;i<WaveList.n;i++) {
    size_t n;
    const char *wave=WaveList.waves[i].name;
    
/*-----------------------------------------------------------------------------
    decode names */
    tide_decode_atlasname(NULL,aPC,wave, 0,&aP);
    __ERR_BASE_LINE__("reading %s :",aP);
    
    for(j=0;j<nv;j++){
      tide_decode_atlasname(NULL,varnames[2*j], 0,  wave,&aN,2);
      tide_decode_atlasname(NULL,varnames[2*j+1], 0,wave,&gN,2);
      
      fprintf(stderr," (%s,%s)",aN,gN);
      
      if(cbuf==NULL){
/*-----------------------------------------------------------------------------
        allocate buffers */
        status=poc_inq_var(aP,aN,&var);
        if(status)__NC_TRAP_ERROR__(wexit,status,1,"poc_inq_var(\"%s\",\"%s\",) error",aP,aN);
        
        status=fe_readgeometry(aP,&mesh);
        
        poc_get_var_length(var,&n);
        fprintf(stderr,"[%d nodes]",n);
        nc=n*nv;
        
        cbuf=new complex<double>[n];
        constants=new hconstant_t[nc];
        for(l=0;l<nc;l++){
          constants[l].init_polar(WaveList.n);
          }
        
        #ifdef OMP_H
        nproc=omp_get_max_threads();//1 set of buffers per process
        #else
        nproc=1;
        #endif
        predbufs=new double*[nproc];
        maxbufs=new double*[nproc];
        for(m=0;m<nproc;m++){
          predbufs[m]=new double[nc];
          maxbufs[m]=new double[n];
          }
        }
      
/*-----------------------------------------------------------------------------
      read atlases */
      status=poc_get_cvara(aP,aN,gN,-1,cbuf,1);
      if(status)__NC_TRAP_ERROR__(wexit,status,1,"poc_get_vara(\"%s\",\"%s\",\"%s\",) error",aP,aN,gN);
      
      for(k=0,l=j;k<n;k++,l+=nv){
        constants[l].a[i]=abs(cbuf[k]);
        constants[l].G[i]=arg(cbuf[k]);
        }
      }
    
    fprintf(stderr," done.\n");
    }
  
  delete[]cbuf;
  
/*-----------------------------------------------------------------------------
  predict */
  init_argument(&astro_angles,start,1);
  
  timeval before;
  gettimeofday(&before);
  
  #pragma omp parallel for schedule(dynamic,1) private(j,k,l,m) if(nproc>1)
  for(i=0;i<nt;i++){
    const double t=startd+increment*i;
    #ifdef OMP_H
    m=omp_get_thread_num();//thread index
    #else
    m=1;
    #endif
    double *predicted=predbufs[m],mod,*maxs=maxbufs[m];
    
    char *ts=poctime_sdate_cnes(t,'s');
    printf("(%d/%d;%05.3fs)%s\r",i,nt,difftime(before),ts);
    fflush(stdout);
    delete[]ts;
    
    harmonic_prediction(astro_angles,predicted,t,nc,WaveList,constants);
    for(k=0,l=0;k<n;k++){
      if(nv>1){
        mod=0.;
        for(j=0;j<nv;j++){
          mod+=square(predicted[l]);
          l++;
          }
        mod=sqrt(mod);
        }
      else{
        mod=predicted[k];
        }
        
/*-----------------------------------------------------------------------------
  statistics */
      if(maxs[k]<mod){
        maxs[k]=mod;
        }
      }
    }
  
  printf("\nmerging...");fflush(stdout);
  gettimeofday(&before);
  
  double *maxs=maxbufs[0];
  
  for(m=1;m<nproc;m++){
    for(k=0;k<n;k++){
      if(maxs[k]<maxbufs[m][k]){
        maxs[k]=maxbufs[m][k];
        }
      }
    delete[]predbufs[m];
    delete[]maxbufs[m];
    }
  
  printf("done in %5.3fs.\n",difftime(before));
  
  poc_global_t global("constructed around " __LINE_FILE_PACKAGE_REVISION);
  var.name="max";
  global<<var;
  global<<poc_att_t("history",cmd);
  string output="max.nc";
  printf("saving to "+output+" :");fflush(stdout);
  if(mesh.vertices){
    printf(" [],");fflush(stdout);
    status=fe_savemesh(output.c_str(),MESH_FILE_FORMAT_NC2D,mesh);
    }
  else{
    status=poc_create(output,global);
    }
  printf(" "+var.name);fflush(stdout);
  status=poc_put_vara(output,var,0,maxs);
  printf(", done.\n");
  
/*-----------------------------------------------------------------------------
  clean-up */
  for(k=0;k<nc;k++){
    constants[k].destroy();
    }
  delete[]constants;
  delete[]predbufs[0];
  delete[]predbufs;
  delete[]maxs;
  delete[]maxbufs;
  
  return 0;
}
