
/*******************************************************************************

  T-UGO tools, 2006-2011

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

\author  Florent Lyard      LEGOS/CNRS, Toulouse, France. florent.lyard@legos.obs-mip.fr
\author  Laurent Roblou     LEGOS/CNRS, Toulouse, France
\author  Damien Allain      LEGOS/CNRS, Toulouse, France
\author  Yves Soufflet      LEGOS/CNRS, Toulouse, France
\author  Clément Mayet      LEGOS, Toulouse, France (PhD)
\author  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
\author  Frédéric Dupont    Université de Laval à Québec, Canada

<!-- USE FIND AND REPLACE TO EDIT THIS LINE, SO THAT print_help IS ALSO UPDATED -->
\brief Calculator.

<!-- A LINK TO main() or print_help() WILL NOT LINK TO THE RIGHT SOURCE ! -->
See the main <a href=#func-members>function</a> for how this works
and the print_help <a href=#func-members>function</a> for how to use this.
*/
/*----------------------------------------------------------------------------*/

#define MAIN_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include "tools-structures.h"

#include "rutin.h"
#include "geo.h"
#include "polygones.h"
#include "functions.h"
#include "grd.h"
#include "netcdf-proto.h"
#include "filter.h"
#include "netcdf-classes.h"

#include "topo.h"

#define XYZ 0
#define YXZ 1

// c **********************************************************
// 	subroutine shell(a,n)
// c
// c shellsort after Sedgewick, Addison-Wesley, 1983
// c at label 30 ge sorts down, le sorts up
// c
// c de mey 1985
// c
// 	dimension a(n)
// 	k= 1
// 10	k= 3*k+1
// 	if (k.le.n) go to 10
// 20	k= k/3
// 	if (k.lt.1) return
// 	do 40 i= k+1, n
// 	  aa= a(i)
// 	  j= i
// 30	  if (a(j-k).ge.aa) go to 40
// 	  a(j)= a(j-k)
// 	  j= j-k
// 	  if (j.gt.k) go to 30 
// 40	  a(j)= aa
// 	go to 20
// 	end


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename TYPE> int topo_rescale(char *bathymetry, grid_t grid, TYPE *topo, TYPE topomask, float maxscale, float scale, char *PolygonsFile)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int   k,l,status;

  float *buffer,*radius,mask,z,d;
  signed char   *selected;

  buffer=new TYPE[grid.nx*grid.ny];
  radius=new float[grid.nx*grid.ny];

  selected=new signed char[grid.nx*grid.ny];
  for(size_t j=0;j<grid.ny;j++) {
    for(size_t i=0;i<grid.nx;i++) {
      size_t m=j*grid.nx+i;
      selected[m]=1;
      radius[m]=0.0;
      }
    }

  if(PolygonsFile!=0) {
    status=topo_PolygonsSelection(grid, PolygonsFile, selected, true, 0);
    if(status!=0) return(-1);
    }

/* *------------------------------------------------------------------------------
  to smooth selected points with prescribed lengthscale*/
  for(size_t j=0;j<grid.ny;j++) {
    for(size_t i=0;i<grid.nx;i++) {
      size_t m=j*grid.nx+i;
      if(selected[m]==0) continue;
      if (topo[m]!=mask) {
        topo[m]*=scale;
        }
      size_t imin=MAX(i-1,0);
      size_t imax=MIN(i+1,grid.nx);
      size_t jmin=MAX(j-1,0);
      size_t jmax=MIN(j+1,grid.ny);
      for(size_t jj=jmin;jj<jmax;jj++) {
        for(size_t ii=imin;ii<imax;ii++) {
          size_t n=jj*grid.nx+ii;
          if(selected[n]==0) radius[m]=maxscale;
          }
        }
      }
    }
  for(size_t j=0;j<grid.ny;j++) {
    for(size_t i=0;i<grid.nx;i++) {
      size_t m=j*grid.nx+i;
      buffer[m]=topo[m];
      }
    }

/* *------------------------------------------------------------------------------
  to smooth surrounding points with decreasing lengthscale*/
  int nRequestedProcs=-1;
  int nprocs=initialize_OPENMP(nRequestedProcs);
//#pragma omp parallel for private(status) if(nprocs>1)
  for(size_t j=0;j<grid.ny;j++) {
    for(size_t i=0;i<grid.nx;i++) {
      size_t m=j*grid.nx+i;
      if (radius[m]==maxscale) {
        size_t imin=MAX(i-maxscale/grid.dx,0);
        size_t imax=MIN(i+maxscale/grid.dx+1,grid.nx);
        size_t jmin=MAX(j-maxscale/grid.dy,0);
        size_t jmax=MIN(j+maxscale/grid.dy+1,grid.ny);
        for(size_t jj=jmin;jj<jmax;jj++) {
          for(size_t ii=imin;ii<imax;ii++) {
            size_t n=jj*grid.nx+ii;
            double dx=(ii-i)*grid.dx;
            double dy=(jj-j)*grid.dy;
            double d=MIN(1.,1.-sqrt(dx*dx+dy*dy)/maxscale);
            d=MAX(d,0.0);
            radius[n]=MAX(radius[n],maxscale*d);
            }
          }
        }
      }
    }

#pragma omp parallel for private(status) if(nprocs>1)
  for(size_t j=0;j<grid.ny;j++) {
    for(size_t i=0;i<grid.nx;i++) {
      size_t m=j*grid.nx+i;
      if (radius[m]!=0) {
        status=loess_filter(grid, topo, mask, radius[m], i, j, buffer[m]);
        }
      }
    }

  for(size_t j=0;j<grid.ny;j++) {
    for(size_t i=0;i<grid.nx;i++) {
      size_t m=j*grid.nx+i;
      topo[m]=buffer[m];
      }
    }

  delete[] buffer;
  delete[] radius;
  delete[] selected;
  
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  vector<plg_t> topo_contours(string rootname, grid_t grid, float*smoothed, vector<float> & isovalues, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  float mask=-10;
  string filename;
  SGfield_t<float> field;
  vector<plg_t> contours, null;
  float *UGarray;
  int verbose=(debug==true);
  float optimal_isovalue=0.2;
  
  frame_t frame;
//   frame=plg_spherical_minmax(polygons);
  
/*------------------------------------------------------------------------------
  optimize contours detection by masking areas where no contours could be expected */
  float *masked=new float[grid.Hsize()];
  for(int m=0;m<grid.Hsize();m++) masked[m]=smoothed[m];

  int border=max(10,10);
  int depth=9;
 
#pragma omp parallel for
  for(int j=border;j<grid.ny-border;j++) {
    for(int i=border;i<grid.nx-border;i++) {
      int k=grid.nx*j+i;
      float value=smoothed[k];
      if( (value!=1) &&  (value!=-1) ) continue;
      int count=0;
/*------------------------------------------------------------------------------
      for a given value, see if neighbour cells show identical values*/
      for(int jj=j-depth;jj<j+depth+1;jj++) {
        for(int ii=i-depth;ii<i+depth+1;ii++) {
          int kk=grid.nx*jj+ii;
          if(kk==k) continue;
          if(smoothed[kk]==mask) continue;
          if(smoothed[kk]!=value) {
            count++;
            goto finish;
            }
          }
        }
    finish:
      if(count==0) {
        masked[k]=mask;
        }
      else {
        masked[k]=value;
        }
      }
    }
  delete[] smoothed;
  smoothed=masked;
  
/*------------------------------------------------------------------------------
  build contours */
  printf("#------------------------------------------------------\n");
  printf("step 3: compute iso_contours\n");
  mesh_t mesh;
    
  for(int k=0;k<isovalues.size();k++) {
    
    plg_destroy_entries(contours);
    contours.clear();
    
    STDOUT_BASE_LINE("treating iso-contour=%lf\n",isovalues[k]);
    
    contours=plg_FindCountours(grid, smoothed, mask, isovalues[k], mesh, UGarray, false, verbose);
    
    if(contours.size()==0) return(contours);
    
    plg_smoothlines(contours);
    
    status=plg_cartesian((projPJ) 0, contours);
    status=plg_recale(contours, frame, PLG_SPHERICAL);
    
    string output;
    std::stringstream tmp;
    tmp.width(4);
    tmp.precision(2);
    tmp.setf( std::ios::fixed, std:: ios::floatfield );
    tmp << fabs(isovalues[k]);
    
    output=rootname+"-contours-" + tmp.str() +".plg";
    status=plg_save(output.c_str(), PLG_FORMAT_SCAN, contours); 
    }

  mesh.destroy();
  if(UGarray!=smoothed) delete[] UGarray;
  delete[] smoothed;
  grid.free();
  
  return(contours);
}

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
  int flag;

  char *rootname=NULL,*output=NULL,*input=NULL,*pathname=NULL,*bathymetry=NULL;
  char *bathy=NULL,*notebook=NULL,*meshbook=NULL,*poly=NULL,*format=NULL;
  char *filename=0;

  grid_t topogrid;
  grid_t grid;
  mesh_t mesh;
  short *topo=0,smask=256*127+255;
  float *ftopo=0,*tmp,*buffer,ftopomask;

  float zmin=-500,zmax=0,radius=0.0;

  plg_t *polygones=NULL;
  int npolygones=0;

  pocgrd_t ncgrid;
  cdfvar_t variable;
  int dimlgth[4]={10,20,30,40};
  
  bool debug=false;
  bool median=false, loess=true;

  fct_echo( argc, argv);

  n=1;
  while (n < argc) {
    keyword=strdup(argv[n]);
    if(strcmp(keyword,"-inv")==0) {
      scale=-1;
      n++;
      continue;
      }
    if(strcmp(keyword,"-scale")==0) {
      sscanf(argv[n+1],"%f",&scale);
      n++;
      n++;
      continue;
      }
    if(strcmp(keyword,"-zmin")==0) {
      sscanf(argv[n+1],"%f",&zmin);
      n++;
      n++;
      continue;
      }
    if(strcmp(keyword,"-zmax")==0) {
      sscanf(argv[n+1],"%f",&zmax);
      n++;
      n++;
      continue;
      }
    if(strcmp(keyword,"--filter=median")==0) {
      median=true;
      loess=false;
      n++;
      continue;
      }
    if(strcmp(keyword,"--filter=loess")==0) {
      median=false;
      loess=true;
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

          case 'r' :
            sscanf(argv[n+1],"%f",&radius);
            n++;
            n++;
            break;


          default:
            if(strcmp(keyword,"-zmin")==0) {
              sscanf(argv[n+1],"%f",&zmin);
              n++;
              n++;
              }
            else if(strcmp(keyword,"-zmax")==0) {
              sscanf(argv[n+1],"%f",&zmax);
              n++;
              n++;
              }
            else {
              __OUT_BASE_LINE__("unknown option %s\n",keyword);
              exit(-1);
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

/* *----------------------------------------------------------------------
  load bathymetry grid*/
  status=grd_loadgrid(bathymetry,&topogrid);
  if(status !=0) {
    printf("cannot load bathymetry file=%s\n",bathymetry);
    goto error;
    }

/* *----------------------------------------------------------------------
  load bathymetry data*/
  ftopo= new float[topogrid.nx*topogrid.ny];
  topogrid.modeH=0;
  topogrid.nz=1;
  status=  grd_loadr1((const char*) bathymetry,topogrid,ftopo,&ftopomask);

  if(radius!=0.0) {
//     zmin=-500;
//    status= topo_smooth_spherical(topogrid, ftopo, ftopomask,  radius, zmin, zmax, poly, input);
    int filter;
    if(median) filter=MEDIAN;
    else filter=LOESS;
    range_t<float> range(zmin,zmax);
   
    bool keep_masked=false;
    
    status=topo_smooth_cartesian(topogrid, ftopo, ftopomask,  radius, range, keep_masked, poly, input, filter);
    if(status!=0) goto error;
    }

  if(scale!=1.0) {
    status= topo_rescale(input, topogrid, ftopo, ftopomask,  radius, scale, poly);
    if(status!=0) goto error;
    }

/*------------------------------------------------------------------------------
  save rectified topo */
  if(rootname==0) rootname=strdup("topo-smooth");
  output=new char[1024];

  topo=new short[topogrid.nx*topogrid.ny];
  for (j=0;j<topogrid.ny;j++) {
    for (i=0;i<topogrid.nx;i++) {
      m=j*topogrid.nx+i;
      n=(topogrid.ny-j-1)*topogrid.nx+i;
      if(ftopo[m]!=mask) {
        topo[n]=(short) floor(ftopo[m]+0.5);
        }
      else {
        topo[n]=256*127+255;
        }
      }
    }
    
  if(format==NULL) format=strdup("grd");
  
  filename=0;
  status=topo_save(rootname, &filename, format, topogrid, ftopo, ftopomask, debug);
  
  status=grd_DuplicateTag(bathymetry, filename);
  
  if(topo!=0) delete[] topo;
  delete[] ftopo;

  topogrid.free();

  __OUT_BASE_LINE__("topo-smooth: bathymetry sucessfully completed\n");

  exit(0);

 error:
  __ERR_BASE_LINE__("topo-smooth error, abort\n");exit(-1);

}
