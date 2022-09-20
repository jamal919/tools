#define MAIN_SOURCE

#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
 

#include "tools-structures.h"
     
#include "fe.h"
#include "archive.h"


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  float *rbuf[2],*tmp[4];
  float *rbufx[2],*rbufy[2];
  double t0,dT,pulsation;
  double *serie[500],a1,p1,a2,p2,zr,zi,d;

  float  t,p;
  float  dummy,rmask;
  float  spec;
  int nndes,status,*elts,fmt,fefmt;
  int i,j,k,l,L,m,n;
  int nitems,element;
  FILE *file;
  FILE *out;
  char *keyword,*zone;
  char *meshfile=NULL,*sfile[2]={NULL,NULL},*cfile=NULL,*format=NULL,*nodefile=NULL;
  char file1[256],file2[256];
  char *root,*output=NULL,rootname[1024];
  grid_t grid;
  fcomplex *cbuf[2],cmask;
  fcomplex *cbufx[2],*cbufy[2];
  mesh_t mesh,nodes;
  int count=0;
  char *comment[2],preview[1024];

  spec=999.9;
  comment[0]=(char *)malloc(1024);
  comment[1]=(char *)malloc(1024);

  fct_echo(argc,argv);
 
  n=1;
  while (n < argc)
    {
    keyword=strdup(argv[n]);
    switch (keyword[0])
      {
      case '-':
        switch (keyword[1])
        {
        case 'm' :
          meshfile= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'o' :
          output= strdup(argv[n+1]);
          n++;
          n++;
          break;

        default:
          __OUT_BASE_LINE__("unknown option %s\n",keyword);
          exit(-1);
        }
        break;

      default:
        sfile[count]= strdup(argv[n]);
        printf("input file=%s\n",sfile[count]);
        count++;
        n++;
        break;
      }
      free(keyword);
    }

  if(strstr(sfile[0],".s2r") != NULL)    fefmt=S2R;
  if(strstr(sfile[0],".s2c") != NULL)    fefmt=S2C;
  if(strstr(sfile[0],".v2r") != NULL)    fefmt=V2R;
  if(strstr(sfile[0],".v2c") != NULL)    fefmt=V2C;
    
  if(meshfile != NULL) {
    status=fe_readmesh(meshfile,MESH_FILE_FORMAT_TRIGRID,&mesh);
    if(status !=0) goto error;
    status=fe_list(&mesh);
    if(status !=0) goto error;
    }
  else
    {
    printf("no mesh file specified; abort...\n");
    goto error;
    }
  

  if(strrchr(sfile[0],(int) '/') != NULL) {
    l=strlen(strrchr(sfile[0],(int) '/'));
    strncpy(rootname,strrchr(sfile[0],(int) '/')+1,l);
    }
  else  sprintf(rootname,"%s.operation",sfile[0]);

  if(output==NULL) output=strdup(rootname);
 
  switch (fefmt)
    {
    case S2R:
      rbuf[0]=(float *)malloc(mesh.nvtxs*sizeof(float));
      rbuf[1]=(float *)malloc(nodes.nvtxs*sizeof(float));
      status=quoddy_loadr1(sfile[0], mesh.nvtxs, rbuf[0]);
      if(status!=0) goto error;
      status=quoddy_saver1(output, nodes.nvtxs, rbuf[1],comment);
      if(status!=0) goto error;
      free(rbuf[0]);
      free(rbuf[1]);
      break;
    case V2R:
      sprintf(comment[0],"summation of %s and %s",sfile[0],sfile[1]);
      sprintf(comment[1],"for more details, see files' comment");
      rbufx[0]=(float *)malloc(mesh.nvtxs*sizeof(float));
      rbufx[1]=(float *)malloc(mesh.nvtxs*sizeof(float));
      rbufy[0]=(float *)malloc(mesh.nvtxs*sizeof(float));
      rbufy[1]=(float *)malloc(mesh.nvtxs*sizeof(float));
/*
      status=quoddy_loadr2(sfile[0], mesh.nvtxs, rbufx[0], rbufy[0]);
      status=quoddy_loadr2(sfile[1], mesh.nvtxs, rbufx[1], rbufy[1]);
*/
      status=load_v2r(sfile[0], rbufx[0], rbufy[0],mesh.nvtxs,0);
      status=load_v2r(sfile[1], rbufx[1], rbufy[1],mesh.nvtxs,0);
      if(status!=0) goto error;
      for (n=0;n<mesh.nvtxs;n++) {
        rbufx[1][n]+=rbufx[0][n];
        rbufy[1][n]+=rbufy[0][n];
        }
      status=quoddy_saver2(output, mesh.nvtxs, rbufx[1], rbufy[1],comment);
      free(rbufx[0]);
      free(rbufx[1]);
      free(rbufy[0]);
      free(rbufy[1]);
      break;
    case S2C:
      cbuf[0]=(fcomplex *)malloc(mesh.nvtxs*sizeof(fcomplex));
      tmp[0]=(float *)malloc(mesh.nvtxs*sizeof(float));
      tmp[1]=(float *)malloc(mesh.nvtxs*sizeof(float));
      cbuf[1]=(fcomplex *)malloc(nodes.nvtxs*sizeof(fcomplex));
      status=quoddy_loadc1(sfile[0],mesh.nvtxs , cbuf[0]);
      for (n=0;n<mesh.nvtxs;n++) tmp[0][n]=real(cbuf[0][n]);
      for (n=0;n<mesh.nvtxs;n++) tmp[1][n]=imag(cbuf[0][n]);
      status=quoddy_savec1(output, nodes.nvtxs, cbuf[1],comment);
      free(cbuf[0]);
      free(cbuf[1]);
      free(tmp[0]);
      free(tmp[1]);
      break;
    case V2C:
      cmask=fcomplex(999.9,999.9);
      cbufx[0]=(fcomplex *)malloc(mesh.nvtxs*sizeof(fcomplex));
      cbufx[1]=(fcomplex *)malloc(nodes.nvtxs*sizeof(fcomplex));
      cbufy[0]=(fcomplex *)malloc(mesh.nvtxs*sizeof(fcomplex));
      cbufy[1]=(fcomplex *)malloc(nodes.nvtxs*sizeof(fcomplex));
      tmp[0]=(float *)malloc(mesh.nvtxs*sizeof(float));
      tmp[1]=(float *)malloc(mesh.nvtxs*sizeof(float));
      tmp[2]=(float *)malloc(mesh.nvtxs*sizeof(float));
      tmp[3]=(float *)malloc(mesh.nvtxs*sizeof(float));
      status=quoddy_loadc2(sfile[0],mesh.nvtxs , cbufx[0], cbufy[0]);
      for (n=0;n<mesh.nvtxs;n++) tmp[0][n]=real(cbufx[0][n]);
      for (n=0;n<mesh.nvtxs;n++) tmp[1][n]=imag(cbufx[0][n]);
      for (n=0;n<mesh.nvtxs;n++) tmp[2][n]=real(cbufy[0][n]);
      for (n=0;n<mesh.nvtxs;n++) tmp[3][n]=imag(cbufy[0][n]);
      if(status!=0) goto error;
      status=quoddy_savec2(output, nodes.nvtxs, cbufx[1], cbufy[1],comment);
      free(cbufx[0]);
      free(cbufx[1]);
      free(cbufy[0]);
      free(cbufy[1]);
      free(tmp[0]);
      free(tmp[1]);
      free(tmp[2]);
      free(tmp[3]);
      break;

    default:
      break;
    }


end: printf("end of operation ... \n");
  free(elts);
error:
  __ERR_BASE_LINE__("exiting\n");exit(-1);
}
