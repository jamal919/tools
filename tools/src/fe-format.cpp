

/**************************************************************************

  T-UGO tools, 2006-2009

  Unstructured Ocean Grid initiative

Contributors:

  Florent Lyard      LEGOS/CNRS, Toulouse, France
  Yoann Le Bars      LEGOS, Toulouse, France (PhD)
  Laurent Roblou     LEGOS, Toulouse, France
  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
  Frédéric Dupont    Université de Laval à Québec, Canada

E-mail: florent.lyard@legos.obs-mip.fr

***************************************************************************/

#define MAIN_SOURCE

#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "tools-define.h"
#include "tools-structures.h"

#include "functions.h"
#include "fe.h"
#include "bmg.h"
#include "map.h"
#include "ascii.h"
#include "archive.h"
#include "geo.h"
#include "sym-io.h"
#include "netcdf-proto.h"

extern  int surfref_loadr1(char *name, int nvalues, int ncolumns, int pos, int *translation, float *buffer);
extern  int *surfref_reference(char *name, mesh_t mesh, int *ncolumns);
extern  int surfref_reorder(char *name, mesh_t mesh, int ncolumns, int *translation, char *ouput);

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int loadr1 (char* filename,int nndes, int ncolumns, int pos, int *translation, float *buf, float *mask, int fmt )

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{

  int status;
  char *root,output[1024],vname[1024];
  pocgrd_t ncgrid;
  cdfvar_t variable;
  float time;

  switch(fmt) {
    case S2R:
      status=quoddy_loadr1(filename, nndes, buf);
      break;

    case SRF:
      status= surfref_loadr1(filename, nndes, ncolumns, pos, translation, buf);
      break;
    }

  return(status);

}



/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  float *rbuf;
  float  dummy,rmask,*rbuffer,mask;
  float  *rbufferx,*rbuffery;
  float  scale=1.,offset=0.;
  float  spec=1.e+10;
  int nndes,status,*elts,fmt_in,fmt_out,fefmt;
  int *translation,ncolumns,dim;
  int i,j,k,l,L,m,n;
  int nitems;
  FILE *file;
  FILE *out;
  char *meshfile=NULL,*sfile=NULL,*cfile=NULL,*sfmt_in=NULL,*sfmt_out=NULL,*cdl=0;
  char *keyword,*zone=NULL,*gridfile=NULL,*notebook=NULL;
  char file1[256],file2[256];
  char *root,tmp[256],output[1024];
  grid_t grid,cgrid;
  fcomplex *cbuffer,cmask;
  fcomplex *cbufferx,*cbuffery;
  mesh_t mesh;
  cdfvar_t variable;
  pocgrd_t ncgrid;
  string cmd;
  char *comments[2];

  fct_echo( argc, argv);

  comments[0]=(char *)malloc(1024);
  comments[1]=(char *)malloc(1024);
  sprintf(comments[0],"The C->CPP transformation");
  sprintf(comments[1],"The C->CPP transformation");

  spec=999.9;

  n=1;
  while (n < argc) {
    keyword=strdup(argv[n]);
    switch (keyword[0]) {
      case '-':
        switch (keyword[1]) {
        case 't' :
          sfmt_in=strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'f' :
          sfmt_out= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'o' :
          strcpy(output,argv[n+1]);
          n++;
          n++;
          break;

        case 'm' :
          meshfile= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 's' :
          sscanf(argv[n+1],"%f",&spec);
          n++;
          n++;
          break;

        case 'c' :
          cdl= strdup(argv[n+1]);
          n++;
          n++;
          break;

        default:
          __OUT_BASE_LINE__("unknown option %s\n",keyword);
          exit(-1);
        }
        break;

      default:
        if(sfile==NULL) {
          sfile= strdup(argv[n]);
          printf("input file=%s\n",sfile);
          n++;
          }
        else {
          __OUT_BASE_LINE__("unknown option %s\n",keyword);
          exit(-1);
          }
        break;
      }
      free(keyword);
    }

  rmask=spec;
  cmask=fcomplex(spec,spec);

  if(meshfile != NULL) {
    status=fe_readmesh(meshfile,MESH_FILE_FORMAT_TRIGRID,&mesh);
    if(status !=0) goto error;
    status=fe_list(&mesh);
    if(status !=0) goto error;
    status= fe_edgetable(&mesh,0,0);
    if(status !=0) goto error;
    }
  else {
    printf("no mesh file specified; abort...\n");
    goto error;
    }

  nndes=mesh.nvtxs;
  //fe_initaffine(mesh);

  if(sfmt_out==NULL) sfmt_out= strdup("netcdf");

  if(strstr(sfile,".s2r") != NULL)    fefmt=S2R;
  if(strstr(sfile,".s2c") != NULL)    fefmt=S2C;
  if(strstr(sfile,".v2r") != NULL)    fefmt=V2R;
  if(strstr(sfile,".v2c") != NULL)    fefmt=V2C;
  if(strstr(sfile,".dat") != NULL)    fefmt=S2R;

  if(strcmp(sfmt_in,"s2r")==0) {
    fefmt=S2R;
    dim=S2R;
    }

  if(strcmp(sfmt_in,"s2c")==0)   fefmt=S2C;
  if(strcmp(sfmt_in,"v2r")==0)   fefmt=V2R;
  if(strcmp(sfmt_in,"v2c")==0)   fefmt=V2C;


/* *----------------------------------------------------------------------------
  SURFREF specials*/
  if(strcmp(sfmt_in,"srf")==0) {
    fefmt=SRF;
    dim=S2R;
    translation=surfref_reference(sfile, mesh, &ncolumns);
    }
  if(strcmp(sfmt_out,"srf")==0) {
     strcpy(output,sfile);
     root=strstr(output,".dat");
     sprintf(root,"%s",".mesh-array.dat");
     status=surfref_reorder(sfile, mesh, ncolumns, translation, output);
     __ERR_BASE_LINE__("exiting\n");exit(0);
    }

  if(cdl != NULL) {
    cmd=string("ncgen -b -x ")+string(cdl)+string(" -o "+string(output));
    system(cmd.c_str());
    }
  switch (dim) {
    case S2R:
      rbuffer=(float *)malloc(nndes*sizeof(float));
      for (k=0;k<ncolumns;k++) {
        status=loadr1(sfile, nndes, ncolumns, k, translation, rbuffer,&rmask,fefmt);
        if(status!=0) goto error;
        if(strcmp(sfmt_out,"ascii")==0) {
          }
        if(strcmp(sfmt_out,"s2r")==0) {
          sprintf(output,"%s.%2.2d.s2r","ouput",k);
          sprintf(comments[0],"%c POC/Noveltis, extracted from %s",169,sfile);
          sprintf(comments[1],"column %d",k);
          status=quoddy_saver1(output, mesh.nvtxs,rbuffer,comments);
          }
        else if(strcmp(sfmt_out,"netcdf")==0) {
          strcpy(output,sfile);
          root=strstr(output,".s2r");
          sprintf(root,"%s",".nc");
          status= fe_savemeshNC2D(output,mesh, 1);
          variable = poc_variable_UG2D("name_undefined", spec,"units_undefined", scale, offset,"standardname_undefined","N");
          status   = cdf_createvariable(output, &(variable));
          status   = poc_put_UG2D(output, mesh, variable.id, rbuffer);
          variable.destroy();
          }
        }
      free(rbuffer);
      break;
    case V2R:
      rbufferx=(float *)malloc(nndes*sizeof(float));
      rbuffery=(float *)malloc(nndes*sizeof(float));
      status=quoddy_loadr2(sfile, nndes, rbufferx,rbuffery);
      if(strcmp(sfmt_out,"ascii")==0) {
        }
      else if(strcmp(sfmt_out,"netcdf")==0) {
        strcpy(output,sfile);
        root=strstr(output,".s2r");
        sprintf(root,"%s",".nc");
        }
      free(rbufferx);
      free(rbuffery);
      break;
    case S2C:
      cbuffer=(fcomplex *)malloc(nndes*sizeof(fcomplex));
      status=quoddy_loadc1(sfile, nndes, cbuffer);
      if(strcmp(sfmt_out,"ascii")==0) {
        }
      else {
        }
      free(cbuffer);
      break;
    case V2C:
      cmask=fcomplex(999.9,999.9);
      cbufferx=(fcomplex *)malloc(nndes*sizeof(fcomplex));
      cbuffery=(fcomplex *)malloc(nndes*sizeof(fcomplex));
      status=quoddy_loadc2(sfile, nndes, cbufferx, cbuffery);
      if(strcmp(sfmt_out,"ascii")==0) {
        }
      else {
        }
      free(cbufferx);
      free(cbuffery);
      break;

    default:
      break;
    }


end: printf("end of fe-format ... \n");
  free(elts);
error:
  __ERR_BASE_LINE__("exiting\n");exit(-1);
}
