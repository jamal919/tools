
/**************************************************************************

  T-UGO tools, 2006-2013

  Unstructured Ocean Grid initiative

***************************************************************************/
/** \file

\author  Florent Lyard      LEGOS/CNRS, Toulouse, France. florent.lyard@legos.obs-mip.fr
\author  Yoann Le Bars      LEGOS, Toulouse, France (PhD)
\author  Laurent Roblou     LEGOS/CNRS, Toulouse, France
\author  Damien Allain      LEGOS/CNRS, Toulouse, France
\author  Sara Fleury        LEGOS/CNRS, Toulouse, France
\author  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
\author  Frédéric Dupont    Université de Laval à Québec, Canada

<!-- USE FIND AND REPLACE TO EDIT THIS LINE, SO THAT print_help IS ALSO UPDATED -->
\brief USED BY CTOH

<!-- A LINK TO main() or print_help() WILL NOT LINK TO THE RIGHT SOURCE ! -->
See the main function for how this works
and the print_help function for how to use this.
*/
/*----------------------------------------------------------------------------*/

#define MAIN_SOURCE

#include "version-macros.def" //for VERSION and REVISION

#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "tools-structures.h"

#include "fe.h"
#include "poc-time.h"
#include "archive.h"
#include "rutin.h"     /*  rutin.h contains common utility routines  */
#include "functions.h"


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
    "  %s ...\n",prog_name);
  printf("\n"
    "DESCRIPTION\n"
/* USE FIND AND REPLACE TO EDIT THE LINE BELOW, SO THAT THE FILE HEAD IS ALSO UPDATED */
    "  ... USED BY CTOH.\n"
    "  ...\n"
    "\n"
    "OPTIONS :\n"
    "  -h  Show this help and exit.\n"
    "  -o  \n"
    "  -r  \n"
    "  -f  \n"); /** \endcode */
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  float z,velocity[2],elevation,zero=0,beta[3];
  double dt,cut,mask=0;
  double t0,dT,pulsation;
  double sum[2],a1,p1,a2,p2,zr,zi,d;

  float  **h;
  float  *u,*v;
  float  count,dummy, hours,start,end,w,dum;
  double *hf,*mhf,*mlf,*lf,*vlf,vhf;
  double time;
  double *woce_h,*woce_t,woce_mask,dmask;
  int origine,woce_n, nodes[3];
  int fmt,i,j,k,l,L,ndum,rstatus,lfdepth,mfdepth;
  int minstep,maxstep,l1,l2,n,status,nst;
  int  year,month,day,hour;
  int nnde,nstep,elt,m,resampling=2;
  date_t reference;
  double chk01,chk02;

  size_t size,offset;
  FILE *in,*out,*sout[100];
  FILE *mesh;
  char file1[256],file2[256], wocefile[256], model[256];
  char fout[100][256], outname[256]="\0",rootname[256]="\0";
  pressure_station *sample;
  char *option,*keyword,*s,*input,*datafile,*output, *start_name;
  int filtering=0,fix=0;
  meta_archive_t info,reduced;
  grid_t *grid;
  pressure_station *setptr;

  fct_echo( argc, argv);
  n=1;
  while (n < argc) {
    keyword=strdup(argv[n]);
    switch (keyword[0]) {
      case '-':
        switch (keyword[1]) {
        case 'o' :
          n++;
          output= strdup(argv[n]);
          n++;
          break;

        case 'r' :
          sscanf(argv[n+1],"%d",&resampling);
          n++;
          n++;
          break;

        case 'f' :
          fix=1;
          n++;
          break;

        case 'h':
          print_help(argv[0]);
          exit(0);

        default:
          __OUT_BASE_LINE__("unknown option %s\n",keyword);
          print_help(argv[0]);
          exit(-1);
        }
        break;

      default:
        datafile= strdup(argv[n]);
        n++;
        break;
      }
      free(keyword);
    }

  L=strlen(datafile);
  if(strrchr(datafile,(int) '/') != NULL) {
    l=strlen(strrchr(datafile,(int) '/'));
    strncpy(rootname,datafile,L-l+1);
    }
  else  strcpy(rootname,"./");

  status=clx_archiveinfo(datafile,&info);
  if(status !=0) goto error;

  status=clx_archiveinfo(datafile,&reduced);
  if(status !=0) goto error;

  reduced.sampling=resampling*info.sampling;
  reduced.nframe=info.nframe/resampling;

#if 0 /* sara 21/04/2010 */
/*if (output==NULL) {
    output=malloc(strlen(datafile)+10);
    strcpy(output,rootname);*/
  output=(char *) malloc(1024);
  sprintf(output,"%s.reduced-%d",datafile,resampling);
/*} */
#else

  /* file name without path */
  if ((start_name=strrchr(datafile, '/')) == NULL)
    start_name = datafile;
  else
    start_name += 1;

  output=(char*)malloc(strlen(start_name)+16);
  sprintf(output,"%s.reduced-%d",start_name,resampling);

#endif

  printf("output: %s \n",output);

  reference=info.reference;
  printf("scales : %f %f %f\n",info.scale[0],info.scale[1],info.scale[2]);

  status=fe_list(&info.mesh);
  if(status !=0) goto error;

  printf("mesh setup:\n",status);

  nnde=info.mesh.nvtxs;
  h=smatrix(0,2,0,nnde-1);

  in=fopen(datafile,"r");
  if(in ==NULL) goto error;

  status=clx_archiveread(in,info,1,h,&time);
  if(status !=0) goto error;
  fclose(in);

  n=0;
  for (i=1;i<= (int) info.nframe;i+=resampling) {
    if(n!=0) {
      in=fopen(datafile,"r");
      status=clx_archiveread(in,info,i,h,&time);
      fclose(in);
      if(status !=0) goto error;
      }
    fe_integrale01(info.mesh,h[0],sum);
    chk01=sum[0]/sum[1];
    fe_integrale02(info.mesh,h[0],sum,(double) 65.0);
    chk02=sum[0]/sum[1];
    printf("%s :  %7.4lf %7.4lf\n",sgetnewdate(info.reference,time),chk01,chk02);
//     printf("t= %6.1f hours status= %d %s\n",time/3600.0,status, sgetnewdate(info.reference,time));
    if(fix) {
      for (j=0;j<nnde;j++) {
        h[0][j]=h[0][j]-chk01;
        }
      }
    if(n==0) {
      status=clx_archivewriteheader(output,reduced,h);
      if(status !=0) printf("clx_archivewriteheader done (status=%d) \n",status);
      if(status !=0) goto error;
      }
    else {
      status=archive_update2(output, reduced, h, time);
      if(status !=0) goto error;
      }
    n++;
    }

  free(h);

  __ERR_BASE_LINE__("%s -computation complete ^^^^^^^^^\n",argv[0]);
  exit(0);

 error:
  __ERR_BASE_LINE__("%s -computation aborted ^^^^^^^^^\n",argv[0]);
  exit(-1);
}

