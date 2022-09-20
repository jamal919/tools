/**************************************************************************

  variability

***************************************************************************/

#define MAIN_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include "tools-structures.h"

#include "fe.h"
#include "archive.h"
#include "functions.h"

extern  int quoddy_savec1(char *, int , fcomplex *, char **);


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_statc(mesh_t mesh,char *filelist,char *depthfile, char *rootname)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  float  *rbuffer,count;
  double *rms;
  int status,*elts,fmt,fefmt;
  int i,j,k,l,L,m,n;
  int nitems;
  FILE *in;
  char hfile[1024],output[1024],*comment[3];
  complex<double> *average;
  complex<float> *cbuffer;

  average=new complex<double>[mesh.nvtxs];
  rms    =new double[mesh.nvtxs];

  cbuffer=new complex<float>[mesh.nvtxs];
  rbuffer=new float[mesh.nvtxs];

  for (i=0;i<mesh.nvtxs;i++) {
    average[i]=complex<float>(0.0,0.0);
    rms[i]=0.0;
    }

  if(depthfile!=NULL) {
    status=quoddy_loadr1(depthfile, mesh.nvtxs, rbuffer);
    if(status!=0) goto error;
    }

  in=fopen(filelist,"r");

  count=0;
  for (k=0;k<1000;k++) {
    nitems=fscanf(in, "%s", hfile);
    if(nitems!=1) goto end;
    count++;
    printf("treating %s \n",hfile);
    status=quoddy_loadc1(hfile, mesh.nvtxs, cbuffer);
    if(status!=0) {
      printf("cannot read %s\n",hfile);
      goto error;
      }
    for (i=0;i<mesh.nvtxs;i++) {
      average[i]+=cbuffer[i];
      rms[i]+=norm(cbuffer[i]);
      }
    }

end:

  for (i=0;i<mesh.nvtxs;i++) {
    average[i]/=count;
    rms[i]=sqrt(rms[i]/count-norm(average[i]));
    average[i]/=count;
    cbuffer[i]=average[i];
    rbuffer[i]=rms[i];
    }
  comment[0]=(char *)malloc(256);
  comment[1]=(char *)malloc(256);

  sprintf(comment[0],"average from %s file list",filelist);
  sprintf(comment[1],"harmonic constants (m)");

  sprintf(output,"%s.mean.s2c",rootname);
  status=quoddy_savec1((const char *) output, mesh.nvtxs, cbuffer,comment);

  sprintf(comment[0],"rms from %s file list",filelist);
  sprintf(comment[1],"quadratique rms (m)");

  sprintf(output,"%s.rms.s2r",rootname);
  status=quoddy_saver1(output, mesh.nvtxs, rbuffer,comment);

  delete[] rbuffer;
  delete[] rms;
  delete[] cbuffer;
  delete[] average;

  printf("end of variability ... \n");
  return(0);

error:
  printf("error detected, quit ... \n");
  return(-1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_statr(mesh_t mesh,char *filelist,char *depthfile, char *rootname)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  float  *rbuffer,count;
  double *rms;
  int status,*elts,fmt,fefmt;
  int i,j,k,l,L,m,n;
  int nitems;
  FILE *in;
  char hfile[1024],output[1024],*comment[3];
  double *average;

  average=new double[mesh.nvtxs];
  rms    =new double[mesh.nvtxs];

  rbuffer=new float[mesh.nvtxs];

  for (i=0;i<mesh.nvtxs;i++) {
    average[i]=0.0;
    rms[i]=0.0;
    }

  if(depthfile!=NULL) {
    status=quoddy_loadr1(depthfile, mesh.nvtxs, rbuffer);
    if(status!=0) goto error;
    }

  in=fopen(filelist,"r");

  count=0;
  for (k=0;k<1000;k++) {
    nitems=fscanf(in, "%s", hfile);
    if(nitems!=1) goto end;
    count++;
    printf("treating %s \n",hfile);
    status=quoddy_loadr1(hfile, mesh.nvtxs, rbuffer);
    if(status!=0) {
      printf("cannot read %s\n",hfile);
      goto error;
      }
    for (i=0;i<mesh.nvtxs;i++) {
      average[i]+=rbuffer[i];
      rms[i]+=rbuffer[i]*rbuffer[i];
      }
    }

end:
 
  for (i=0;i<mesh.nvtxs;i++) {
    average[i]/=count;
    rms[i]=sqrt(rms[i]/count-average[i]*average[i]);
    rbuffer[i]=average[i];
    }
  comment[0]=(char *)malloc(256);
  comment[1]=(char *)malloc(256);

  sprintf(comment[0],"average from %s file list",filelist);
  sprintf(comment[1],"harmonic constants (m)");

  sprintf(output,"%s.mean.s2r",rootname);
  status=quoddy_saver1(output, mesh.nvtxs, rbuffer,comment);

  for (i=0;i<mesh.nvtxs;i++) {
    rbuffer[i]=rms[i];
    }

  sprintf(comment[0],"rms from %s file list",filelist);
  sprintf(comment[1],"quadratique rms (m)");

  sprintf(output,"%s.rms.s2r",rootname);
  status=quoddy_saver1(output, mesh.nvtxs, rbuffer,comment);

  for (i=0;i<mesh.nvtxs;i++) {
    rbuffer[i]=100.*rms[i]/MAX(10.,fabs(average[i]));
    }

  sprintf(comment[0],"relative rms from %s file list",filelist);
  sprintf(comment[1],"quadratique rms (m)");

  sprintf(output,"%s.relative-rms.s2r",rootname);
  status=quoddy_saver1(output, mesh.nvtxs, rbuffer,comment);

  delete[] rbuffer;
  delete[] rms;
  delete[] average;

  printf("end of variability ... \n");
  return(0);

error:
  printf("error detected, quit ... \n");
  return(-1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status,fmt,fefmt;
  int i,j,k,l,L,m,n;
  int nitems;
  FILE *in;
  char *keyword;
  char *meshfile=NULL,*depthfile=NULL,*output=NULL,*path=NULL,*filelist=NULL;
  char hfile[1024],ufile[1024],*comment[3];
  mesh_t mesh;

  fct_echo( argc, argv);

  n=1;
  while (n < argc) {
    keyword=strdup(argv[n]);
    switch (keyword[0]) {
      case '-':
        switch (keyword[1]) {

        case 'm' :
          meshfile= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'd' :
          depthfile= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'o' :
          output= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'p' :
          path= strdup(argv[n+1]);
          n++;
          n++;
          break;

        default:
          __OUT_BASE_LINE__("unknown option %s\n",keyword);
          exit(-1);
        }
        break;

      default:
        filelist= strdup(argv[n]);
        n++;
        break;
      }
      free(keyword);
    }

  if(output ==NULL) output=strdup(".");
  if(path ==NULL)   path=strdup(".");

  if(meshfile != NULL) {
    status=fe_readmesh(meshfile,MESH_FILE_FORMAT_TRIGRID,&mesh);
    if(status !=0) goto error;
    }
  else  {
    printf("no mesh file specified; abort...\n");
    goto error;
    }
  status=fe_statr(mesh,filelist,depthfile, output);
//   if(strstr(sfile,".s2r") != NULL)    fefmt=S2R;
//   if(strstr(sfile,".s2c") != NULL)    fefmt=S2C;
//   if(strstr(sfile,".v2r") != NULL)    fefmt=V2R;
//   if(strstr(sfile,".v2c") != NULL)    fefmt=V2C;

  __OUT_BASE_LINE__("end of variability ... \n");
  exit(0);

error:
 __OUT_BASE_LINE__("error detected, quit ... \n");
  exit(-1);
}

