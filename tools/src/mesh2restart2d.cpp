#define MAIN_SOURCE

#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "tools-structures.h"

#include "geo.h"
#include "fe.h"
#include "poc-time.h"
#include "archive.h"

#define LINEAR 0
#define QUADRATIC 1

date_t local_t_reference;
char *local_external_path,*local_external_root;



/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int output_statevector_ascii(double time, double dT, state2D_t *set,int nndes,const char *filename)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,k,index,neq,idum;
  char s[2];
  div_t result;
  char *sdate;
  FILE *out;
  double *ddum,cnestime;
  date_t date;

/*----------------------------------------------------------------------
  Create output files for continuation;
  save model state for further restart, overwrite in 8 turning files
----------------------------------------------------------------------*/

  if ((out = fopen(filename, "w")) == NULL) {
    return(-1);
    }

  date= poctime_getdate(time, local_t_reference, 's', &cnestime);
  sdate=sgetdate(date);
  printf("#Continuation in %s at %s \n",filename,sdate);

  fprintf(out,"%f, %f (model time), h,u,v at %s (t= %f hrs)\n",time,dT,sdate,time/3600.);
  for(i=0;i<nndes;i++) {
    fprintf(out," %lf %lf %lf %f %f", set[i].h,
                    set[i].H, set[i].Hm,
                    set[i].u, set[i].v);
    if(time==0) {
      fprintf(out," %lf %lf %lf %lf %lf\n",0., 0., 0., 0., 0.);
      }
    else
      {
      fprintf(out," %lf %lf %lf %lf %lf",
                    set[i].hmean/time,
                    set[i].umean/time,
                    set[i].vmean/time,
                    set[i].humean/time,
                    set[i].hvmean/time);
      }
    fprintf(out," %d %d\n", set[i].stable,/*suspicious[i]*/0);
    }

  fclose(out);
  free(sdate);
  return(0);

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int nndes,status,*elts,fmt,fefmt;
  int i,j,k,l,L,m,n;
  int option;
  FILE *file;
  FILE *out;  char *keyword,*zone=NULL,*gridfile=NULL;
  char *meshfile=NULL,*depthfile=NULL,*output=NULL,*poly=NULL;
  mesh_t mesh,refined;
  double *depth;
  int *selected;
  char *comment[2],preview[1024];
  float *buffer[2];
  state2D_t *set;
  double time,dT=10.,trestart=0.0,shift=0.0;
 
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

        case 'b' :
          depthfile= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'o' :
          output= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 't' :
          sscanf(argv[n+1],"%lf",&trestart);
          n++;
          n++;
          break;

        case 's' :
          sscanf(argv[n+1],"%lf",&shift);
          n++;
          n++;
          break;

        case '0' :
          fmt=0;
          n++;
          break;

        case '1' :
          fmt=1;
          n++;
          break;

        default:
          __OUT_BASE_LINE__("unknown option %s\n",keyword);
          exit(-1);
        }
        break;

      default:
          local_external_path= strdup(argv[n]);
          n++;
        break;
      }
      free(keyword);
    }

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

 if(local_external_path == NULL) {
   printf("no directory for analysis files specified; abort...\n");
   goto error;
   }

 if(local_external_root == NULL) {
   local_external_root=strdup("analysis");
   printf("%s used as default analysis files rootname\n");
   }

  depth=(double *) malloc(mesh.nvtxs*sizeof(double));

  if(depthfile!=NULL) {
    status=quoddy_loadd1((const char *) depthfile, mesh.nvtxs, depth);
    if(status!=0) goto error;
    }
  else {
    printf("no depth file specified; abort...\n");
    goto error;
    }

  nndes=mesh.nvtxs;
  status= fe_edgetable(&mesh,0,0);

  local_t_reference.day    =01;
  local_t_reference.month  =12;
  local_t_reference.year   =2004;
  local_t_reference.second =0.0;

  set=(state2D_t *) malloc(mesh.nvtxs*sizeof(state2D_t));
  for(n=0;n<mesh.nvtxs;n++) {
    set[n].lon=mesh.vertices[n].lon;
    set[n].lat=mesh.vertices[n].lat;
    set[n].h  =MAX(10.,depth[n]);
    }
  time=600.0;

  time=trestart;

/*   fmt=1; */
  switch (fmt) {
    case 0:
      status= clx_initexternal(time, dT, set,  mesh);
      break;
    case 1:
      status= cdf_initexternal(time, dT, set,  mesh);
      break;
    default:
      __ERR_BASE_LINE__("exiting\n");exit(-1);
    }
/*
  for(n=0;n<mesh.nvtxs;n++) {
    set[n].H  +=set[n].h;
    set[n].Hm +=set[n].h;
    }
*/
  if(fmt==1)
  for(n=0;n<mesh.nvtxs;n++) {
    set[n].h  *=100.;
    set[n].H  *=100.;
    set[n].Hm *=100.;
    set[n].u  *=100.;
    set[n].v  *=100.;
    }

  status=output_statevector_ascii( time+shift,  dT,  set,mesh.nvtxs,"restart.quaker");

  buffer[0]=(float *) malloc(mesh.nvtxs*sizeof(float));
  buffer[1]=(float *) malloc(mesh.nvtxs*sizeof(float));
  comment[0]=(char *) malloc(1024);
  comment[1]=(char *) malloc(1024);

  sprintf(comment[0],"Mog2D simulation, %cPOC/Noveltis, extracted from %s",169,local_external_path);
  sprintf(comment[1],"model elevation (cm)");
  sprintf(preview,"initial-elevation.s2r");

  for(n=0;n<mesh.nvtxs;n++) {
    buffer[0][n]=set[n].H-set[n].h;
    }
  status=quoddy_saver1(preview, mesh.nvtxs, buffer[0],comment);

  sprintf(comment[0],"Mog2D simulation, %cPOC/Noveltis, extracted from %s",169,local_external_path);
  sprintf(comment[1],"model velocity (cm)");
  sprintf(preview,"initial-velocity.v2r");

  for(n=0;n<mesh.nvtxs;n++) {
    buffer[0][n]=set[n].u;
    buffer[1][n]=set[n].v;
    }
  status=quoddy_saver2(preview, mesh.nvtxs, buffer[0], buffer[1],comment);

end: __OUT_BASE_LINE__("end of mesh2mesh ... \n");
  exit(0);
error:
 __OUT_BASE_LINE__("error detected, quit ... \n");
  exit(-1);
}
