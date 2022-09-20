
/* *@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  topo-merg aims at merging 2 bathymetric database with some controls
  on merging action

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */

#define MAIN_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "tools-structures.h"

#include "rutin.h"
#include "geo.h"
#include "polygones.h"
#include "netcdf-proto.h"
#include "functions.h"
#include "grd.h"
#include "map.def"
#include "map.h"

#include "topo.h"

#define XYZ 0
#define YXZ 1


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_average (grid_t grid, float *buffer, float mask, int i0, int j0, float *z, int range)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *-----------------------------------------------------------------------------

  return average value of buffer over a target cell

------------------------------------------------------------------------------*/
{
  int i,j,k,l,m,n,i_range,j_range;
  int mm;
  int status;
  double count,y;
  float tmp;

  count=0;
  tmp=0;
  y=map_grid_y(grid,i0,j0);
    
/* *-----------------------------------------------------------------------------
  fix range, should be computed from grid and target resolutions*/
  i_range=range*cos(y*M_PI/360.);
  j_range=range;
  
  for(l=MAX(0,j0-j_range);l<MIN(grid.ny-1,j0+j_range);l++) {
    for(k=MAX(0,i0-i_range);k<MIN(grid.nx-1,i0+i_range);k++) {
      mm=l*grid.nx+k;
      if(buffer[mm]!=mask) {
        tmp+=buffer[mm];
        count++;
        }
      }
    }

  if(count!=0) {
    tmp/=count;
    *z=tmp;
    status=0;
    }
  else {
    *z=mask;
    status=-1;
    }
    
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  geo_t projection;

  float  scale=1.0,topo_mask=1.e+10,mask=1.e+10;
  double  x,y,t,p;
  double  dx,dy;

  int i,j,k,l,m,n,signus=0,status;

  char *keyword,*zone;
  int flag,average=0;

  char *rootname=NULL,*output=NULL,*input=NULL,*pathname=NULL,*bathymetry=NULL;
  char *vtopo=0,*vdraft=0;
  char *bathy=NULL,*format=NULL,*hydro=NULL,*poly=NULL;
  char *proj4=NULL;

  grid_t topogrid;
  grid_t grid,spherical;
  short *topo=0,smask=256*127+255;
  float *ftopo=0,*draft=0,ftopomask,*fspherical=0;

  float zmin=0,zmax=0;

  plg_t *polygones=NULL;
  int npolygones=0;

  pocgrd_t ncgrid;
  cdfvar_t variable;
  int dimlgth[4]={10,20,30,40};

  fct_echo( argc, argv);

  n=1;
  while (n < argc) {
    keyword=strdup(argv[n]);
    if(strcmp(keyword,"-inv")==0) {
      scale=-1;
      n++;
      continue;
      }
    if(strcmp(keyword,"-average")==0) {
      scale=-1;
      n++;
      continue;
      }
    if(strcmp(keyword,"-add")==0) {
      hydro= strdup(argv[n+1]);
      signus=+1;
      n++;
      n++;
      continue;
      }
    if(strcmp(keyword,"-substract")==0) {
      hydro= strdup(argv[n+1]);
      signus=-1;
      n++;
      n++;
      continue;
      }
    switch (keyword[0]) {
      case '-':
        switch (keyword[1]) {
          case 'b' :
            bathymetry= strdup(argv[n+1]);
            n++;
            n++;
            break;

          case 'f' :
            format= strdup(argv[n+1]);
            n++;
            n++;
            break;

          case 'p' :
            poly= strdup(argv[n+1]);
            n++;
            n++;
            break;

          case 'o' :
            rootname= strdup(argv[n+1]);
            n++;
            n++;
            break;

          case '-' :
            if(strstr(argv[n],"--proj=")!=0){
              proj4=strdup(argv[n]+7);
              n++;
              }
            break;

          default:
            if(strcmp(keyword,"-zmin")==0) {
              sscanf(argv[n+1],"%f",&zmin);
              n++;
              n++;
              }
            if(strcmp(keyword,"-zmax")==0) {
              sscanf(argv[n+1],"%f",&zmax);
              n++;
              n++;
              }
          }
        break;

      default:
        if(input==NULL) {
          input= strdup(argv[n]);
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

  if(rootname==0) rootname=strdup("out");
  output=new char[1024];

//#define timmerman

#ifdef timmerman
/* *----------------------------------------------------------------------
  load bathymetry grid*/
  printf("#####################################################################\n");
  printf("loading bathymetry in %s\n",input);
  status=map_loadfield(input, "bathy", &topogrid, &ftopo, &ftopomask);

/* *----------------------------------------------------------------------
  load ice bottom grid*/
  printf("#####################################################################\n");
  printf("loading ice bottom position in %s\n",input);
  status=map_loadfield(input, "draft", &topogrid, &draft, &ftopomask);
#elifd albmap
/* *----------------------------------------------------------------------
  load bathymetry grid*/
  printf("#####################################################################\n");
  printf("loading bathymetry in %s\n",input);
  status=map_loadfield(input, "topg", &topogrid, &ftopo, &ftopomask);

/* *----------------------------------------------------------------------
  load ice bottom grid*/
  printf("#####################################################################\n");
  printf("loading ice bottom position in %s\n",input);
  status=map_loadfield(input, "lsrf", &topogrid, &draft, &ftopomask);
  
  ftopomask=-9999.;
#else // IBCAO
/* *----------------------------------------------------------------------
  load bathymetry grid*/
  printf("#####################################################################\n");
  printf("loading bathymetry in %s\n",input);
  status=map_loadfield(input, MAP_FORMAT_CODE_COMODO, "z", 0, 0, 0, &topogrid, &ftopo, &ftopomask);

#endif
  
  if(status!=0) goto error;
  
  printf("#####################################################################\n");
  printf("processing bathymetry\n");
  if(draft!=0) {
/* *----------------------------------------------------------------------
  compute free water column height*/
  for (j=0;j<topogrid.ny;j++) {
    for (i=0;i<topogrid.nx;i++) {
      m=j*topogrid.nx+i;
      if(ftopo[m]!=ftopomask) {
        ftopo[m]-=draft[m];
        }
      }
    }
  }
  
  if(average==1) {
/* *----------------------------------------------------------------------
 smooth column height limits*/
  for (j=0;j<topogrid.ny;j++) {
    for (i=0;i<topogrid.nx;i++) {
      float z;
      m=j*topogrid.nx+i;
      if(ftopo[m]==0.0) {
        status=map_average (topogrid, ftopo, ftopomask, i, j, &z, 4);
        ftopo[m]=z;
        }
      }
    }
    }
    
  if(proj4!=0) {
    projPJ ref;
    char **parms;
    char *s,*tmp;
    int nparms,ntoken;
    double t,p,x,y,x0,y0,t0,p0;
    float z;

    ref = init_projection(proj4,true);
    t=135;
    p=-75;
//     t=135;
//     p=-90;

#ifdef albmap
    t=0;
    p=-90;
    geo_to_projection(ref,p,t,&x0,&y0);
    map_set2Dgrid(&spherical,-180.0,-85.0,+180.0,-60.0,1./60.,1./60.);
#else //IBCAO
    t=0;
    p=90;
    geo_to_projection(ref,p,t,&x0,&y0);
    map_set2Dgrid(&spherical,-180.0,55.0,+180.0,90.0,1./120.,1./120.);
#endif
    x=topogrid.x[0]+x0;
    y=topogrid.y[0]+y0;
    projection_to_geo(ref,&p0,&t0,x,y);

    map_completegridaxis(&spherical,1);
    
    fspherical=new float[spherical.ny*spherical.nx];
    for (j=0;j<spherical.ny;j++) {
      for (i=0;i<spherical.nx;i++) {
//         m=j*spherical.nx+i;
//         t=spherical.x[m];
//         p=spherical.y[m];
        spherical.xy(i,j,t,p);
        geo_to_projection(ref,p,t,&x,&y);
        status=map_interpolation(topogrid,ftopo,ftopomask,x-x0,y-y0,&z);
//        if(z==0.0) z=ftopomask;
        fspherical[m]=z;
        }
      }
    sprintf(output,"%s.regular.nc",rootname);
    status= poc_createfile(output);
    status=poc_sphericalgrid_xy(output,"",spherical,&ncgrid);
    poc_standardvariable_xy(&variable,"bathymetry",ftopomask,"m",1., 0.,"bathymetry","bathymetry","bathymetry",ncgrid);
    status=create_ncvariable(output, &variable);
    status=poc_write_xy(output,  spherical, variable.id,fspherical);
    variable.destroy();

    for (j=0;j<topogrid.ny;j++) {
      for (i=0;i<topogrid.nx;i++) {
        m=j*topogrid.nx+i;
        x=topogrid.x[m]+x0;
        y=topogrid.y[m]+y0;
        projection_to_geo(ref,&p,&t,x,y);
        t=geo_recale(t,t0+180.0,180.0);
        topogrid.x[m]=t;
        topogrid.y[m]=p;
        }
      }
    
    pj_free(ref);
    }


/* *------------------------------------------------------------------------------
  apply change of sign on final topography*/
  for (j=0;j<topogrid.ny;j++) {
    for (i=0;i<topogrid.nx;i++) {
      m=j*topogrid.nx+i;
      if(ftopo[m]!=ftopomask) {
        ftopo[m]*=scale;
        }
      }
    }
 

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  save data, to be replaced with topo_save. Example:
  
  status=topo_save(rootname, output, "netcdf", TOPOgrid, TOPOdata, TOPOmask, debug);
  
  topo_save will hande mirroring...

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

/* *------------------------------------------------------------------------------
  save rectified topo */
  if(format==NULL) format=strdup("grd");


  printf("#####################################################################\n");
  printf("save bathymetry\n");
  if(strcmp(format,"grd")==0) {
    topo=new short[topogrid.nx*topogrid.ny];
    smask=256*127+255;
    for (j=0;j<topogrid.ny;j++) {
      for (i=0;i<topogrid.nx;i++) {
        m=j*topogrid.nx+i;
        n=(topogrid.ny-j-1)*topogrid.nx+i;
        if(ftopo[m]!=ftopomask) {
          topo[n]=(short) floor(scale*ftopo[m]+0.5);
          }
        else {
          topo[n]=smask;
          }
        }
      }
    sprintf(output,"%s.grd",rootname);
    status=grd_save(output,topogrid,topogrid.nx, topo,smask);
    delete[] topo;
    }
  else if(strcmp(format,"netcdf")==0) {
    sprintf(output,"%s.spherical.nc",rootname);
    status= poc_createfile(output);
    status=poc_sphericalgrid_xyzt(output,"","",topogrid,&ncgrid);
    poc_standardvariable_xy(&variable,"bathymetry",ftopomask,"m",1., 0.,"bathymetry","bathymetry","bathymetry",ncgrid);
    status=create_ncvariable(output, &variable);
    status=poc_write_xy(output,  topogrid, variable.id,ftopo);
    variable.destroy();
    }

  delete[] ftopo;
  topogrid.free();

  __OUT_BASE_LINE__("topo_merge sucessfully completed\n");

  exit(0);

 error:
  __OUT_BASE_LINE__("topo_merge aborted\n");
  exit(-1);
}
