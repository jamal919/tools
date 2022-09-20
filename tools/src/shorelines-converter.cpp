
/*******************************************************************************

  T-UGO tools, 2006-2013

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

\author  Florent Lyard      LEGOS/CNRS, Toulouse, France. florent.lyard@legos.obs-mip.fr
\author  Laurent Roblou     LEGOS/CNRS, Toulouse, France
\author  Damien Allain      LEGOS/CNRS, Toulouse, France
\author  Yves Soufflet      LEGOS/CNRS, Toulouse, France
\author  Clément Mayet      LEGOS, Toulouse, France (PhD)
\author  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
\author  Frédéric Dupont    Université de Laval à Québec, Canada

<!-- USE FIND AND REPLACE TO EDIT THIS LINE, SO THAT print_help IS ALSO UPDATED -->
\brief Convert between different coastline formats.

<!-- A LINK TO main() or print_help() WILL NOT LINK TO THE RIGHT SOURCE ! -->
See the main <a href=#func-members>function</a> for how this works
and the print_help <a href=#func-members>function</a> for how to use this.
*/
/*----------------------------------------------------------------------------*/


#define MAIN_SOURCE

#include "version-macros.def" //for VERSION and REVISION

#include <stdio.h>
#include <string.h>
#include <fstream>
#include <algorithm>
#include <iostream>
#include <limits>
#include <string>
#include <map>

#include "tools-structures.h"

#include "fe.h"
#include "geo.h"
#include "archive.h"
#include "netcdf-proto.h"
#include "poc-time.h"
#include "functions.h"


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int plg_load_raw(const char *filename, plg_t **polygones, int *npolygones, int reorder)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  plg_t *plg_polygones, *tmp;
  int k,i,j,n=0,npoly=0,dum;
  int count;
  FILE *polyfile;
  char line[1024],*s;
  int format;

  polyfile=fopen(filename,"r");
  if(polyfile==NULL) {
    printf("can't open polygone file %s\n",filename);
    return(-1);
    }

/**----------------------------------------------------------------------------
  check format */
  s=fgets(line,1024,polyfile);
  n=strlen(s);
  
  if(strncmp(s,"#NXY",n-1)==0) {
    format=0;
    }
  else if(strncmp(s,"#XY",n-1)==0) {
    format=1;
    }
  else {
    format=1;
    }
    
/**----------------------------------------------------------------------------
  check content */
  npoly=1;
  plg_polygones=new plg_t[npoly];
  plg_polygones[0].npt=0;
  
  do {
    s=fgets(line,1024,polyfile);
    plg_polygones[0].npt++;
    } while(s!=0);

  plg_polygones[0].npt--;
  rewind(polyfile);

/**----------------------------------------------------------------------------
  read file */
  n=0;
  s=fgets(line,1024,polyfile);
  
  plg_polygones[n].t= new double[plg_polygones[n].npt];
  plg_polygones[n].p= new double[plg_polygones[n].npt];
  plg_polygones[n].x= new double[plg_polygones[n].npt];
  plg_polygones[n].y= new double[plg_polygones[n].npt];
  
  for(i=0;i<plg_polygones[n].npt;i++) {
    s=fgets(line,1024,polyfile);
    for(k=0;k<strlen(line);k++) {
      if(s[k]==',')  s[k]=' ';
      if(s[k]=='\n') s[k]=' ';
      }
    switch(format) {
      case 0:
        sscanf(line,"%d %lf %lf",&dum, &(plg_polygones[n].x[i]),&(plg_polygones[n].y[i]));
        break;
      case 1:
        sscanf(line,"%lf %lf",&(plg_polygones[n].x[i]),&(plg_polygones[n].y[i]));
        break;
      }
    plg_polygones[n].t[i]=plg_polygones[n].x[i];
    plg_polygones[n].p[i]=plg_polygones[n].y[i];
    }

  fclose(polyfile);
  if(reorder==0) {
    *polygones=plg_polygones;
    *npolygones=npoly;
    return(-0);
    }
  
/**----------------------------------------------------------------------------
  re-order */
  int *sequence = new int[plg_polygones[n].npt];
  int *passed   = new int[plg_polygones[n].npt];
  double *distance=new double[plg_polygones[n].npt];
  
  for(i=0;i<plg_polygones[n].npt;i++) sequence[i]=i;
  for(i=0;i<plg_polygones[n].npt;i++) passed[i]=0;
  
  double dx,dx2,dy,dy2,dmin,d;
  sequence[0]=0;
  passed[0]=1;
  for(k=0;k<plg_polygones[n].npt-1;k++) {
    i=sequence[k];
    dmin=1.e+20;
    for(j=0;j<plg_polygones[n].npt;j++) {
      if(passed[j]==1) continue;
      dy=plg_polygones[n].y[j]-plg_polygones[n].y[i];
      dy2=dy*dy;
      if(dy2>dmin) continue;
      dx=plg_polygones[n].x[j]-plg_polygones[n].x[i];
      dx2=dx*dx;
      if(dx2>dmin) continue;
      d=dx2+dy2;
      if(d<dmin) {
        dmin=d;
        sequence[k+1]=j;
        distance[k]=sqrt(d);
        }
      }
    passed[sequence[k+1]]=1;
    }

  count=1;
  for(i=0;i<plg_polygones[0].npt-1;i++) {
    d=distance[i];
    if(d>10000.) count++;
    }
    
  tmp=new plg_t[count];
  n=0;
  for(i=0;i<plg_polygones[0].npt-1;i++) {
    d=distance[i];
    if(d>10000.) {
      n++;
      tmp[n].npt=1;
      }
    else {
      tmp[n].npt++;
      }
    }

  int m=0;
  for(n=0;n<count;n++) {
    tmp[n].t= new double[tmp[n].npt];
    tmp[n].p= new double[tmp[n].npt];
    tmp[n].x= new double[tmp[n].npt];
    tmp[n].y= new double[tmp[n].npt];
    for(i=0;i<tmp[n].npt;i++) {
      j=sequence[m];
      tmp[n].t[i]=plg_polygones[0].t[j];
      tmp[n].p[i]=plg_polygones[0].p[j];
      tmp[n].x[i]=plg_polygones[0].x[j];
      tmp[n].y[i]=plg_polygones[0].y[j];
      m++;
      }
    }
  
//  *polygones=plg_polygones;
//  *npolygones=npoly;
  *polygones=tmp;
  *npolygones=count;

  return(-0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void print_help(char *prog_name)

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
    "  %s [OPTIONS]\n",prog_name);
  printf("\n"
    "DESCRIPTION\n"
/* USE FIND AND REPLACE TO EDIT THE LINE BELOW, SO THAT THE FILE HEAD IS ALSO UPDATED */
    "  Convert between different coastline formats.\n"
    "\n"
    "OPTIONS :\n"
    "  -h,--help  Show this help and exit.\n"
    "  -i  followed by the path of the input file.\n"
    "  -o  followed by the path of the output file. Default from the input file and output format if both are given.\n"
    +isFormat_help("  ",
      ". Default from the extension of the input file.",
      ". Default from the extension of the output file.")+
    "  -d  followed by the decimation factor. Default : 1 (no decimation).\n"
    "  -l  followed by the minimum length. Default : 0 (no decimation).\n"
    "  --proj=...  projection parameters. Default : no projection.\n");
  printf("\n"
    "FILE FORMATS :\n");
  plg_print_formats();
  printf("\n"
   "TIP\n"
   "sed -re 's/^[0-9 ]{3}[0-9] .*/nan nan/;t;s/ *[0-9]+ //' <world.scan\n");
  /** \endcode */
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,l;
  int n,status;
  string input,output;
  statistic_t stat;
  string inFormatName,outFormatName;
  int inFormat=PLG_FORMAT_UNKNOWN,outFormat=PLG_FORMAT_UNKNOWN;
  int decimation=1;
  plg_t *polygones;
  int npolygones;
  vector<plg_t> q;
  char *proj4=NULL;
  char *s;
  double MinLength=0.0;
  bool debug=false;
  
  fct_echo( argc, argv);

  n=1;
  while (n < argc) {
    const string keyword=argv[n];
    
    if(strstr(argv[n],"--proj=")!=0){
      proj4=strdup(&argv[n][7]);
      n++;
      continue;
      }
    
    if(isFormat(argv,&n,&inFormatName,&outFormatName)){
      continue;
      }
    
    switch (keyword[0]) {
      case '-':
        switch (keyword[1]) {

        case 'i' :
          input=argv[n+1];
          n++;
          n++;
          break;
        
        case 'o' :
          output=argv[n+1];
          n++;
          n++;
          break;
        
        case 'd' :
          decimation=strtol(argv[n+1],&s,10);
          if(s==argv[n+1]){
            printf("unknown decimation factor %s\n",argv[n+1]);
            print_help(argv[0]);
            exit(-1);
            }
          n++;
          n++;
          break;
        
        case 'l' :
          MinLength=strtol(argv[n+1],&s,10);
          if(s==argv[n+1]){
            printf("unknown decimation factor %s\n",argv[n+1]);
            print_help(argv[0]);
            exit(-1);
            }
          n++;
          n++;
          break;
        
        case 'h' :
          print_help(argv[0]);
          exit(0);

        case '-' :
          if(keyword=="--help"){
            print_help(argv[0]);
            exit(0);
            }
          break;

        default:
          printf("unknown option "+keyword+"\n");
          print_help(argv[0]);
          exit(-1);
        }
        break;

      default:
        printf("unknown option "+keyword+"\n");
        print_help(argv[0]);
        exit(-1);
      }
    }
  
  if(input=="") {
    fprintf(stderr,"*** input missing ***\n");
    print_help(argv[0]);
    exit(-1);
    }
  if(inFormatName!=""){
    inFormat=plg_format_from_name(inFormatName);
    if(inFormat==PLG_FORMAT_UNKNOWN){
      plg_print_formats();
      TRAP_ERR_EXIT(-1,"could not recognise INPUT format "+inFormatName+"\n");
      }
    }
  else{
    inFormat=plg_find_format(input);
    if(inFormat==PLG_FORMAT_UNKNOWN){
      plg_print_formats();
      TRAP_ERR_EXIT(-1,"could not recognise format of "+input+"\n");
      }
    }
  
  if(output=="" && outFormatName=="") {
    fprintf(stderr,"*** output or output format missing ***\n");
    print_help(argv[0]);
    exit(-1);
    }
    
  if(outFormatName!="") {
    outFormat=plg_format_from_name(outFormatName);
    if(outFormat==PLG_FORMAT_UNKNOWN){
      plg_print_formats();
      TRAP_ERR_EXIT(-1,"could not recognise OUTPUT format "+outFormatName+"\n");
      }
    if(output=="") {
      string extension=plg_format_extension(outFormat);
      if(extension==""){
        extension="."+outFormatName;
        }
      output=input+extension;
      }
    }
  else if(output!="" && outFormatName=="") {
    outFormat=plg_find_format(output);
    if(outFormat==PLG_FORMAT_UNKNOWN){
      plg_print_formats();
      TRAP_ERR_EXIT(-1,"could not recognise OUTPUT format of "+output+"\n");
      }
    }
  
  
/*------------------------------------------------------------------------------
  load shorelines data set*/
//   status=plg_load(input, inFormat, q);
//   if(status!=0) NC_TRAP_ERROR(exit,status,1,"plg_load(\""+input+"\",%d,,) error",inFormat);
//   
//   status=plg_decimate(q,5.0,debug);
//   
//   status=plg_decimate(q,5.0,debug);
//   
//   status=plg_decimate(q,5.0,debug);
//   
//   status=plg_decimate(q,5.0,debug);
//   
//   status=plg_save(output.c_str(), outFormat, q);
//   if(status!=NC_NOERR) NC_TRAP_ERROR(wexit,status,1,"plg_save(\""+output+"\",%d,,) error",outFormat);
  
/*------------------------------------------------------------------------------
  load shorelines data set*/
  status=plg_load(input, inFormat, &polygones, &npolygones);
  if(status!=0) NC_TRAP_ERROR(exit,status,1,"plg_load(\""+input+"\",%d,,) error",inFormat);
  
  printf("npolygons=%d\n",npolygones);
  
  if(decimation>1){
    l=0;
    for(n=0;n<npolygones;n++){
      status=plg_decimate(&polygones[n],decimation,debug);
      if(l<n){
        polygones[l]=polygones[l];
        }
      if(polygones[n].npt>2)
        l++;
      else{
        polygones[n].destroy();
        }
      }
    npolygones=l;
    }
  
//   if(proj4!=0) {
//     projPJ ref;
//     double t,p,x,y;
// 
//     ref = init_projection(proj4,true);
//     
//     for(n=0;n<npolygones;n++) {
//       for(k=0;k<polygones[n].npt;k++) {
//         if(polygones[n].x==0) {
//           x=polygones[n].t[k];
//           y=polygones[n].p[k];
//           }
//         else {
//           x=polygones[n].x[k];
//           y=polygones[n].y[k];
//           }
//         projection_to_geo(ref,&p,&t,x,y);
//         polygones[n].t[k]=t;
//         polygones[n].p[k]=p;
//         }
//       }
//     
//     pj_free(ref);
//     }
    
  if(proj4!=0) {
    int nprocs=initialize_OPENMP(-1, 0);
    
    projPJ *refs = init_projection_parallel(proj4,nprocs,true,0);
    
#pragma omp parallel for if(nprocs>1)
    for(n=0;n<npolygones;n++) {
      int id=omp_get_thread_num();
      double x,y,t,p;
      for(int k=0;k<polygones[n].npt;k++) {
        if(polygones[n].x==0) {
          x=polygones[n].t[k];
          y=polygones[n].p[k];
          }
        else {
          x=polygones[n].x[k];
          y=polygones[n].y[k];
          }
        projection_to_geo(refs[id], &p, &t, x, y);
        polygones[n].t[k]=t;
        polygones[n].p[k]=p;
        }
      }
    
    deletep2D(&refs,nprocs,free_threadSafe_projection);
    }
    
  if(MinLength>0.0){
    l=0;
    for(n=0;n<npolygones;n++){
//       double length=polygones[n].length();
      double length=plg_length(polygones[n],0);
      if(length<MinLength){
        polygones[n].destroy();
        }
      }
    }
  
  
  STDERR_BASE_LINE("saving %d polygones to "+output+"\n",npolygones);
  status=plg_save(output.c_str(), outFormat, polygones, npolygones);
  if(status!=NC_NOERR) NC_TRAP_ERROR(wexit,status,1,"plg_save(\""+output+"\",%d,,) error",outFormat);

//  status=plg_load(output, PLG_FORMAT_NETCDF, &polygones, &npolygones);

  STDOUT_BASE_LINE(" ^^^^^^^^^^^^^ end of computation ^^^^^^^^^^^^^\n");
  exit(0);
}
