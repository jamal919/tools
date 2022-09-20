
/*******************************************************************************

  T-UGO tools, 2006-2019

  Unstructured Ocean Grid initiative

Contributors:

  Florent Lyard      LEGOS/CNRS, Toulouse, France
  Yoann Le Bars      LEGOS, Toulouse, France (PhD)
  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
  Frédéric Dupont    Université de Laval à Québec, Canada

E-mail: florent.lyard@legos.obs-mip.fr

*******************************************************************************/
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  mapxyz-v0 aims to convert complete or parially complete 2D data array in xyz 
  format

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

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
#include "xyz.h"
#include "map.h"

#define XYZ 0
#define YXZ 1

#define ROW    0
#define COLUMN 1

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  geo_t projection;

  float  scale=1.0,topo_mask=1.e+10,mask=1.e+10;
  double  x,y,t,p;
  double  dx,dy;

  int i,j,k,l,m,n,status;

  size_t size;
  FILE *in,*out;
  char *keyword,*zone;
  int flag,inverse=0,ncolumn=0;

  char *rootname=NULL,*output=NULL,*input=NULL,*pathname=NULL,*bathymetry=NULL;
  char *bathy=NULL,*notebook=NULL,*meshbook=NULL,*poly=NULL;
  char *proj4=NULL;
  vector<string> varnames;

  grid_t topogrid, native_grid;
  grid_t grid;
  mesh_t mesh;
  short *topo,smask=256*127+255;
  float *ftopo,*tmp,*buffer;

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
    if(strcmp(keyword,"-mask")==0) {
      scanf("%f",&mask);
      n++;
      n++;
      continue;
      }
    if(strstr(argv[n],"--proj=")!=0){
      proj4=strdup(argv[n]+7);
      n++;
      continue;
      }
    if(strcmp(keyword,"-inverse")==0) {
      inverse=1;
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

        case 'v' :
          varnames.push_back(argv[n+1]);
          n++;
          n++;
          break;

        default:
          __OUT_BASE_LINE__("unknown option %s\n",keyword);
          exit(-1);
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

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  process xyz file

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

/*------------------------------------------------------------------------------
  read xyz file and interpret grid */
  printf("#################################################################\n");
  printf("load grid from file : %s\n",input);
  status=xyz_loadgrid(input, proj4, &native_grid, &ncolumn);
  if(status !=0) {
    __OUT_BASE_LINE__("cannot load grid in bathymetry file=%s %d\n",input,status);
    exit(-1);
    }
  map_printgrid(native_grid);
  
  topogrid.duplicate(native_grid);
 
  if(proj4!=0) {
/*------------------------------------------------------------------------------
    convert cartesian coordinate to spherical ones */
    status=map_projection_backward(topogrid, proj4);
    }
    
/*------------------------------------------------------------------------------
  save rectified topo */
  if(rootname==0) rootname=strdup("out");
  output=new char[1024];
//   sprintf(output,"%s.spherical.nc",rootname);
  sprintf(output,"%s.nc",rootname);

  printf("#################################################################\n");
  printf("save gridded data in netcdf file : %s\n",output);
  status= poc_createfile(output);

//   status=poc_sphericalgrid_xyzt(output,"","",topogrid,&ncgrid);
  status=poc_sphericalgrid_xy(output, "", topogrid, &ncgrid);
  
/*------------------------------------------------------------------------------
  read depth database */
  printf("#################################################################\n");
  printf("load data from file : %s\n",input);
  for(int k=0;k<ncolumn;k++) {
    buffer=new float[native_grid.nx*native_grid.ny];
    status=xyz_loadr1(input, native_grid, k, buffer, &mask);
    if(status !=0) {
      __OUT_BASE_LINE__("cannot load data in bathymetry file=%s\n",input);
      exit(-1);
      }
 
    if(inverse==1) {
      for (j=0;j<native_grid.ny;j++) {
        for (i=0;i<native_grid.nx;i++) {
          m=j*native_grid.nx+i;
          if(buffer[m]!=mask) {
            buffer[m]*=-1.0;
            }
          }
        }
      }
    char *name;
    if(k<varnames.size()) {
      name=strdup(varnames[k].c_str());
      }
    else {
      name=new char[1024];
      sprintf(name, "variable_%d",k);
      }
      
    output=new char[1024];
//     sprintf(output,"%s.spherical.nc",rootname);
    sprintf(output,"%s.nc",rootname);

    poc_standardvariable_xy(&variable,name,mask,"m",1., 0.,name,name,name,ncgrid);
    status=create_ncvariable(output, &variable);
    status=poc_write_xy(output, topogrid, variable.id,buffer);
    variable.destroy();

/*------------------------------------------------------------------------------
    create a degenerated grd file (loss of true grid positions)*/
    frame_t frame=frame_t(topogrid.xmin,topogrid.xmax,topogrid.ymin,topogrid.ymax);
    if(k==0) grid=map_Rgrid(frame, topogrid.dx, topogrid.dy, 2);

    ftopo=new float[grid.nx*grid.ny];
    status=map_export(topogrid,buffer,mask,grid,ftopo,mask,0);

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  save data, to be replaced with topo_save. Example:
  
  status=topo_save(rootname, output, "netcdf", TOPOgrid, TOPOdata, TOPOmask, debug);
  
  topo_save will hande mirroring...

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

    topo=new short[grid.nx*grid.ny];
    for (j=0;j<grid.ny;j++) {
      for (i=0;i<grid.nx;i++) {
        m=j*grid.nx+i;
        n=(grid.ny-j-1)*grid.nx+i;
        if(ftopo[m]!=mask) {
          topo[n]=(short) floor(scale*ftopo[m]+0.5);
          }
        else {
          topo[n]=256*127+255;
          }
        }
      }

    sprintf(output,"%s-%d-short.grd",rootname,k);
    status=grd_save(output,topogrid,topogrid.nx, topo,smask);

    status=grd_mirror_r( grid, grid.nx, ftopo, mask);
    sprintf(output,"%s-%d-float.grd",rootname,k);
    printf("#################################################################\n");
    printf("save gridded data in grd file : %s\n",output);
    status=grd_save(output, topogrid, topogrid.nx, ftopo, mask);

    delete[] topo;
    delete[] ftopo;
    delete[] buffer;
    }

  __OUT_BASE_LINE__("bathymetry sucessfully completed\n");

  exit(0);

error:
  __ERR_BASE_LINE__("exiting\n");exit(-1);

}
