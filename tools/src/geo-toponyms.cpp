
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
#include "map.h"

#include "topo.h"

#define XYZ 0
#define YXZ 1


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int geo_select01(grid_t grid, float *topo, float mask, char *poly, signed char *selected, int initialise)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int   k,i,j,m,n,status;
  int   count=-1;
  int   inside;
  edge_t *edges;
  double t,p;
  plg_t *polygones=NULL;
  int npolygones=0;
  int nprocs;
  frame_t frame;

  printf("#################################################################\n");
  printf("select grid nodes in polygon: %s\n",poly);
  
  if(initialise==1) {
    for(j=0;j<grid.ny;j++) {
      for(i=0;i<grid.nx;i++) {
        k=j*grid.nx+i;
        selected[k]=1;
        }
      }
    }

/*-----------------------------------------------------------------------------
  area criterion */

  if(poly!=NULL) status=plg_load_scan((const char *) poly, &polygones, &npolygones);
  else goto end;
  
  if(status!=0) return(-1);

  frame=plg_cartesian_minmax(polygones,npolygones);
 
  nprocs=initialize_OPENMP(-1);

#pragma omp parallel for private(i,j,n,inside,t,p,status) if(nprocs>1)
 
  for (j=0;j<grid.ny;j++) {
    for (i=0;i<grid.nx;i++) {
      n=grid.nx*j+i;
      if(selected[n]==0) continue;
      t=map_grid_x(grid,i,j);
      p=map_grid_y(grid,i,j);
      if(t<frame.xmin) t+=360.0;
      if(t>frame.xmax) t-=360.0;
      if ((t<frame.xmin) || (t>frame.xmax)) {
        selected[n]=0;
        continue;
        }
      if ((p<frame.ymin) || (p>frame.ymax)) {
        selected[n]=0;
        continue;
        }
      inside=plg_TestInterior(t,p,polygones,npolygones);
      if (inside==PLG_POINT_INTERIOR) {
        selected[n]=1;
        }
      else {
        selected[n]=0;
        }
      }
    }

  count=0;
  for (j=0;j<grid.ny;j++) {
    for (i=0;i<grid.nx;i++) {
      n=grid.nx*j+i;
      if (selected[n]==1) {
        count++;
        }
      }
    }

end:
  printf ("number of edges selected %d\n",count);
  return(count);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int set_code(grid_t topogrid, float *topo, float topomask, char *poly, int code, int masked_only)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int   i,j,k,l,m,n,status;
  int initialise;
  float *buffer,mask,z;
  grid_t grid;
  signed char *selected;
  double x,y;

  selected=new signed char[topogrid.nx*topogrid.ny];
  for(j=0;j<topogrid.ny;j++) {
    for(i=0;i<topogrid.nx;i++) {
      k=j*topogrid.nx+i;
      selected[k]=1;
      }
    }

  initialise=1;
  if(poly!=0) {
    status=geo_select01(topogrid, topo, topomask, poly, selected, initialise);
    initialise=0;
    }
  if(status <0) {
    printf("selection failed, abort...\n");
    return(-1);
    }
    
  if(masked_only==1) {
    for(j=0;j<topogrid.ny;j++) {
      for(i=0;i<topogrid.nx;i++) {
        k=j*topogrid.nx+i;
        if(topo[k]!=topomask) selected[k]=0;
        }
      }
    }
    
/* *------------------------------------------------------------------------------
  Interpolate and store topo on symphonie grid*/
//  mask=-31072.0;
  for(j=0;j<topogrid.ny;j++) {
    for(i=0;i<topogrid.nx;i++) {
      k=j*topogrid.nx+i;
      if(selected[k]==0) continue;
      if(topo[k]!=topomask) topo[k]=code;
      }
    }

  delete[] selected;

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  float  scale=1.0,topo_mask=1.e+10,mask=1.e+10;
  double  x,y,t,p;
  double  dx,dy;

  int i,j,k,l,m,n,signus=0,masked_only=0,status;

  size_t size;
  FILE *in,*out;
  char *keyword,*zone;
  int flag;

  char *rootname=NULL,*output=NULL,*input=NULL,*pathname=NULL,*bathymetry=NULL;
  char *bathy=NULL,*format=NULL,*hydro=NULL,*poly=NULL;

  grid_t topogrid;
  grid_t grid,zone_grid;
  mesh_t mesh;
  short *topo,smask=256*127+255;
  float *ftopo,*tmp,*buffer,ftopomask;

  float zmin=0,zmax=0;

  plg_t *polygones=NULL;
  int npolygones=0;

  pocgrd_t ncgrid;
  cdfvar_t variable;
  int dimlgth[4]={10,20,30,40};
  int code;
  
  bool debug=false;

  fct_echo( argc, argv);

  n=1;
  while (n < argc) {
    keyword=strdup(argv[n]);
    if(strcmp(keyword,"-inv")==0) {
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
    if(strcmp(keyword,"--masked-only")==0) {
      masked_only=1;
      n++;
      continue;
      }
    if(strcmp(keyword,"-code")==0) {
      sscanf(argv[n+1],"%d",&code);
      n++;
      n++;
      continue;
      }
    if(strcmp(keyword,"-debug")==0) {
      debug=true;
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

          case 'z' :
            zone= strdup(argv[n+1]);
            n++;
            n++;
            break;

          default:
            printf("unknown option %s\n",keyword);
            break;
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

/*------------------------------------------------------------------------------
  load bathymetry grid*/
  printf("#####################################################################\n");
  printf("load working toponyms file: %s\n",bathymetry);
  status=topo_loadfield((const char*) bathymetry, &topogrid, &ftopo, &ftopomask, debug);
  if(status !=0) {
    printf("cannot load toponyms file=%s\n",bathymetry);
    goto error;
    }

  
  if(input!=NULL) {
    printf("#####################################################################\n");
    printf("import bathymetry depth flag: %s\n", input);
    zmin=0;
    zmax=0;
    status=topo_import(input, 0, topogrid, ftopo, ftopomask, zmin, zmax, (char *) 0, 0, debug);
    if(status !=0) {
      goto error;
      }
    for(m=0;m<topogrid.ny*topogrid.nx;m++) {
      if(ftopo[m]>0) {
        ftopo[m]=ftopomask;
        }
      else {
        ftopo[m]=0;
        }
      }
    }

/*------------------------------------------------------------------------------
  import new depths*/
  if(poly!=NULL) {
    printf("#####################################################################\n");
    printf("set code from polygon: %s\n", poly);
    status= set_code(topogrid, ftopo, ftopomask, poly, code, masked_only);
    if(status !=0) {
      goto error;
      }
    }


/*------------------------------------------------------------------------------
  save rectified topo */
  if(rootname==0) rootname=strdup("out");
  output=new char[1024];

  if(format==NULL) format=strdup("grd");


/*------------------------------------------------------------------------------
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

  printf("#####################################################################\n");
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
    printf("save toponym file: %s\n",output);
    status=grd_save ((const char*) output,topogrid,topogrid.nx, topo,smask);
    delete[] topo;
    }
  else if(strcmp(format,"netcdf")==0) {
    sprintf(output,"%s.spherical.nc",rootname);
    printf("save toponym file: %s\n",output);
    status= poc_createfile(output);
    status=poc_sphericalgrid_xyzt(output,"","",topogrid,&ncgrid);
    poc_standardvariable_xy(&variable,"bathymetry",ftopomask,"m",1., 0.,"bathymetry","bathymetry","bathymetry",ncgrid);
    status=create_ncvariable(output, &variable);
    status=poc_write_xy(output,  topogrid, variable.id,ftopo);
    variable.destroy();
    }

  delete[] ftopo;
  topogrid.free();

  __OUT_BASE_LINE__("geo-toponyms sucessfully completed\n");

  exit(0);

 error:
  __OUT_BASE_LINE__("geo-toponyms aborted\n");
  exit(-1);
}
