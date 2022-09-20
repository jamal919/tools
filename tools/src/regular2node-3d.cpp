
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



#include "tools-structures.h"

#include "fe.h"
#include "map.h"
#include "netcdf-proto.h"

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int nndes,status,*elts,fmt,fefmt;
  int i,j,k,l,L,m,n;
  int nitems;
  int variable=11; /*mean est la 11e variable*/
  char *keyword,*zone;
  char *meshfile=NULL,*inputT=NULL,*inputS=NULL,*nodefile=NULL;
  char *root,output[1024];
  grid_t grid,slice_grid;
  mesh_t mesh,nodes;
  variable_t varinfo;
  int frame;
  float *buf[2],spec[2];
 
  double x,y,t,p;
  float *profile_x[2], *profile_y[2];
  int nloc;
  float *T,*S;

  int nlevel=33;
  float level[33]={0,-10,-20,-30,-50,-75,-100,-125,-150,-200,-250,-300,-400,-500,-600,-700,-800,-900,-1000,-1100,-1200,-1300,-1400,-1500,-1750,-2000,-2500,-3000,-3500,-4000,-4500,-5000,-5500};
  FILE *out;

  fct_echo(argc,argv);
 
  n=1;
  while (n < argc) {
    keyword=strdup(argv[n]);
    switch (keyword[0]) {
      case '-':
        switch (keyword[1]) {
        case 'n' :
          nodefile= strdup(argv[n+1]);
          printf("nodfile=%s\n",nodefile);
          n++;
          n++;
          break;

        default:
          __OUT_BASE_LINE__("unknown option %s\n",keyword);
          exit(-1);
        }
        break;

      default:
        if(inputT==NULL) {
          inputT= strdup(argv[n]);
          printf("input file=%s\n",inputT);
          n++;
          }
        if(inputS==NULL) {
          inputS= strdup(argv[n]);
          printf("input file=%s\n",inputS);
          n++;
          }

        else
          {
          __OUT_BASE_LINE__("unknown option %s\n",keyword);
          exit(-1);
          }
        break;
      }
    free(keyword);
    }

/****************************load FV nodes*******************************/
  status=fe_readnodes_TGL((const char*) nodefile,&nodes);
  if(status !=0) goto error;

/****************************load netcdf grid****************************/
/**----------------------------------------------------------------------------
  obsolete call, will be suppressed */
  status=cdf_loadvargrid_3d (inputT,variable, &grid);/*température*/
  if(status !=0) goto error;
/**----------------------------------------------------------------------------
  obsolete call, will be suppressed */
  status=cdf_loadvargrid_3d (inputS,variable, &grid); /*salinité�*/
  if(status !=0) goto error;
  
/*-----------------------------------------------------------------------
  allocate memory for buffer */
   buf[0]=(float *) malloc(grid.nx*grid.ny*grid.nz*sizeof(float));
   buf[1]=(float *) malloc(grid.nx*grid.ny*grid.nz*sizeof(float));
  
/*-----------------------------------------------------------------------
  load netcdf variable */
  frame=0;
/**----------------------------------------------------------------------------
  obsolete call, will be suppressed */
  status= cdf_loadvar_r1_3d (inputT, variable, frame, grid, buf[0], &spec[0] ,&varinfo);
/**----------------------------------------------------------------------------
  obsolete call, will be suppressed */
  status= cdf_loadvar_r1_3d (inputS, variable, frame, grid, buf[1], &spec[1] ,&varinfo);

  out=fopen("cars_meanTS.dat","w");
  fprintf(out,"Standard level depths\n");
  fprintf(out,"  %d\n",nlevel);
  for (i=0;i<nlevel;i++) fprintf(out,"%8.1f",level[i]*-1);
  fprintf(out,"\n");
  fprintf(out,"Initial temperature and salinity\n");
  fprintf(out,"observed\n");

  T=(float *) malloc(nlevel*sizeof(float));
  S=(float *) malloc(nlevel*sizeof(float));

  for (n=0;n<nodes.nvtxs;n++) {
/*-----------------------------------------------------------------------
    interpolate profiles at node position */
    t=nodes.vertices[n].lon;
    p=nodes.vertices[n].lat;

    status=map_profile02(grid,t,p, buf[0],spec[0],&profile_x[0], &profile_y[0], &nloc);
    status=map_profile02(grid,t,p, buf[1],spec[1],&profile_x[1], &profile_y[1], &nloc);
  
    for (k=0;k<nlevel;k++) {
      status=map_interpolate1D(profile_x[0],profile_y[0],spec[0],nloc,(float) level[k], &T[k],1);
/*-----------------------------------------------------------------------
       save mog3d input file */
      if (T[k]==spec[0]) {
        if(k!=0) T[k]=T[k-1];
        else {
          printf("trouble at %d, t=%f p= %f z=%f, exit\n",n,t,p,level[k]);
/* 	  __ERR_BASE_LINE__("exiting\n");exit(-1); */
          }
        }
      fprintf(out,"%8.1f",T[k]);
      }
    fprintf(out,"\n");
/*--------------------------------------------------------------salinit�---------------*/
   for (k=0;k<nlevel;k++) {
     status=map_interpolate1D(profile_x[1],profile_y[1],spec[1],nloc,(float) level[k], &S[k],1);
/*-----------------------------------------------------------------------
     save mog3d input file*/
      if (S[k]==spec[1]) {
        if(k!=0) S[k]=S[k-1];
        else {
          printf("trouble at %d, t=%f p= %f z=%f, exit\n",n,t,p,level[k]);
/*        __ERR_BASE_LINE__("exiting\n");exit(-1); */
          }
        }
      fprintf(out,"%8.1f",S[k]);
      }
    fprintf(out,"\n");
    }
  fclose(out);

  free(T);
  free(S);
end: printf("end of regular2node-3d ... \n");
  free(&nodes);
error:
  __ERR_BASE_LINE__("exiting\n");exit(-1);

}

