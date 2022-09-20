
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
\author  Clément MAYET      LEGOS, Toulouse, France (PhD)
\author  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
\author  Frédéric Dupont    Université de Laval à Québec, Canada

<!-- USE FIND AND REPLACE TO EDIT THIS LINE, SO THAT print_help IS ALSO UPDATED -->
\brief Creates an open boundary file for TUGOm.

<!-- A LINK TO main() or print_help() WILL NOT LINK TO THE RIGHT SOURCE ! -->
See the main function for how this works
and the print_help function for how to use this.
*/
/*----------------------------------------------------------------------------*/

#define MAIN_SOURCE

#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cmath>

#include "tools-structures.h"
#include "constants.h"

#include "functions.h"
#include "map.h"
#include "fe.h"
#include "archive.h"

#include "rutin.h"     /*  rutin.h contains common utility routines  */


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int get_openlimits(char *belfile,char *meshfile, double **x, double **y, int *count)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  mesh_t mesh;
  int *code;
  int m,n,n1,n2,status;
  int nopen=0;

  if(meshfile != NULL) {
    status=fe_readmesh(meshfile,MESH_FILE_FORMAT_TRIGRID,&mesh);
    if(status !=0) goto error;
    status=fe_list(&mesh);
    if(status !=0) goto error;
    }
  else {
    printf("no mesh file specified; abort...\n");
    goto error;
    }

  status= fe_edgetable(&mesh,0,0);

  status= fe_vertex_element_tables(&mesh);

  status=fe_read_boundarycode(belfile,mesh,2);

  code=(int *) malloc(mesh.nvtxs*sizeof(int));

  for (n=0;n<mesh.nvtxs;n++) code[n]=0;

  for (n=0;n<mesh.nedges;n++) {
    if(mesh.edges[n].code==5) {
      n1=mesh.edges[n].extremity[0];
      n2=mesh.edges[n].extremity[1];
      code[n1]=1;
      code[n2]=1;
      nopen++;
      }
    }

  *count=0;
  for (n=0;n<mesh.nvtxs;n++)
    if(code[n]!=0) (*count)++;

  (*x)=(double *) malloc(*count*sizeof(double));
  (*y)=(double *) malloc(*count*sizeof(double));

  m=0;
  for (n=0;n<mesh.nvtxs;n++)
    if(code[n]!=0) {
      (*x)[m]=mesh.vertices[n].lon;
      (*y)[m]=mesh.vertices[n].lat;
      m++;
      }

  free(code);

  return(0);

  error:
  return(-1);
}

static bool expired=expire(20150420,20160420);

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
 
  double *lon,*lat;
  int count,nndes,status,*elts,fmt,fefmt;
  int i,j,k,l,L,m,n;
  int nitems;
  FILE *file;
  FILE *out;
  char *keyword,*zone=NULL,*gridfile=NULL;
  char *meshfile=NULL,*belfile=NULL,*targetfile,*output=NULL,*path=NULL;
  char hfile[1024],ufile[1024],filename[1024];
  char *wave[1024];
  grid_t grid;
  fcomplex *cbuffer,*z,*u,*v;
  fcomplex *cbufferx,*cbuffery;
  float a[3],G[3];
  mesh_t mesh;
  int nwave=0,analysis,ncid;
  spectrum_t spectrum;
  string s,executable,echofile;
  int pos;
  FILE *echo;

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

        case 'b' :
          belfile= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 't' :
          targetfile= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'a' :
          sscanf(argv[n+1],"%d",&analysis);
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
        wave[nwave]= strdup(argv[n]);
        printf("input wave=%s\n",wave[nwave]);
        spectrum.n=nwave+1;
        nwave++;
        n++;
        break;
      }
      free(keyword);
    }

  spectrum.waves=new tidal_wave[spectrum.n];
  for (k=0;k<nwave;k++) {
    strcpy(spectrum.waves[k].name,wave[k]);
    }

  if(path ==NULL) path=strdup(".");

  if(meshfile != NULL) {
    status=fe_readmesh(meshfile,MESH_FILE_FORMAT_TRIGRID,&mesh);
    if(status !=0) goto error;
    status=fe_list(&mesh);
    if(status !=0) goto error;
    }
  else {
    printf("no mesh file specified; abort...\n");
    goto error;
    }

  nndes=mesh.nvtxs;

  status=get_openlimits(belfile,targetfile,&lon,&lat,&count);

  elts=(int *) malloc(count*sizeof(int));
  for(k=0; k<count; k++) {
    elts[k]=fe_whichelement(mesh,lon[k],lat[k]);
    if(elts[k]<0) printf("%lf %lf not in mesh\n",lon[k],lat[k]);
    }

  cbuffer =(fcomplex *) malloc(nndes*sizeof(fcomplex));
  cbufferx=(fcomplex *) malloc(nndes*sizeof(fcomplex));
  cbuffery=(fcomplex *) malloc(nndes*sizeof(fcomplex));

  z=(fcomplex *) malloc(count*sizeof(fcomplex));
  u=(fcomplex *) malloc(count*sizeof(fcomplex));
  v=(fcomplex *) malloc(count*sizeof(fcomplex));

  for (k=0;k<nwave;k++) {
    sprintf(hfile,"%s/%s.ele.%2.2d.s2c",path,wave[k],analysis);
    sprintf(ufile,"%s/%s.uv.%2.2d.v2c",path,wave[k],analysis);
    if(status!=0) goto error;
    printf("treating %s wave, from harmonic analysis #%d \n",wave[k],analysis);
    status=quoddy_loadc1(hfile, nndes, cbuffer);
    if(status!=0) {
      printf("cannot read %s\n",hfile);
      goto error;
      }
    status=quoddy_loadc2(ufile, nndes, cbufferx, cbuffery);
    if(status!=0) {
      printf("cannot read %s\n",ufile);
      goto error;
      }
    for (i=0;i<count;i++) {
      status= fe_intpl_LGP1( mesh, cbuffer ,lon[i],lat[i] ,elts[i], &z[i]);
      status= fe_intpl_LGP1( mesh, cbufferx,lon[i],lat[i] ,elts[i], &u[i]);
      status= fe_intpl_LGP1( mesh, cbuffery,lon[i],lat[i] ,elts[i], &v[i]);
      }
/*---------- Tidal heights and currents at open limits ------------*/
    sprintf(filename,"%s.zuv.%2.2d.obc",wave[k],analysis);
    out=fopen(filename,"w");
    fprintf(out,"%s\n",wave[k]);
    for(i=0; i<count; i++) {
/*---------------------------------------------------------------------
      convert elevation amplitude in cm*/
      a[0]=abs(z[i])*100.0;
      G[0]=360.0-arg(z[i])*r2d;
/*---------------------------------------------------------------------
      convert currents amplitude in cm/s*/
      a[1]=abs(u[i])*100.0;
      G[1]=360.0-arg(u[i])*r2d;
      a[2]=abs(v[i])*100.0;
      G[2]=360.0-arg(v[i])*r2d;
      fprintf(out,"%8f %8f %f %f %f %f %f %f\n",lon[i],lat[i],a[0],G[0],a[1],G[1],a[2],G[2]);
      }
    fclose(out);
    }
  printf("free memory \n");

  free(cbuffer);
  free(cbufferx);
  free(cbuffery);

end: printf("end of %s ... \n",argv[0]);
  free(elts);
  __ERR_BASE_LINE__("exiting\n");exit(0);
error:
 __OUT_BASE_LINE__("error detected, quit ... \n");
  exit(-1);
}
