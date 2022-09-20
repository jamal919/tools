
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
#include "statistic.h"

#include "topo.h"

#define XYZ 0
#define YXZ 1


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int bathymetry_error(grid_t & grid, float *topo, float *delta, float mask, double T, string rootname, int & create)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int i,m,n;
  string varname;
  double g=9.81;
    
  float *Kerror=new float[grid.Hsize()];
  float *Gerror=new float[grid.Hsize()];

  for(m=0;m<grid.Hsize();m++) {
    Kerror[m]=mask;
    Gerror[m]=mask;
    double H=topo[m];
    if(H==mask) continue;
    if(H>-5)     continue;
    double d=delta[m];
    if(d==mask) continue;
    if(H+d>-5)     continue;
    double g=9.81;
    double c=sqrt(-g*H);
    double L=c*T;
    double k=2*M_PI/L;
    
    double c_d=sqrt(-g*(H+d));
    double L_d=c_d*T;
    double k_d=2*M_PI/L_d;
    Kerror[m]=k_d/k-1;
    Gerror[m]=360.0*Kerror[m];
    
    if(isnan(Kerror[m])) {
      printf("trouble \n");
      }
    }
  
  varname="Kerror_"+rootname;
  status=save_SG("bathymetry-error.nc", grid, Kerror, mask, varname.c_str(), "dimensionless","Kerror",create);
  create=0;
  
  varname="Gerror_"+rootname;
  status=save_SG("bathymetry-error.nc", grid, Gerror, mask, varname.c_str(), "dimensionless","Gerror",create);

  delete[] Kerror;
  delete[] Gerror;

  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int bathymetry_error(grid_t & grid, float *topo, float *delta, float mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  double T;
  int create=1;
  
  status=map_completegridaxis(&grid, 1);
  
  printf("#################################################################\n");
  printf("compute bathymetry error\n");

  printf("-----------------------------------------------------------------\n");
  printf("K1 error\n");
  T=24*3600;
  status=bathymetry_error(grid, topo, delta, mask, T, "K1", create);
  printf("-----------------------------------------------------------------\n");
  printf("M2 error\n");
  T=12.5*3600;
  status=bathymetry_error(grid, topo, delta, mask, T, "M2", create);
  printf("-----------------------------------------------------------------\n");
  printf("M4 error\n");
  T=6.25*3600;
  status=bathymetry_error(grid, topo, delta, mask, T, "M4", create);
  printf("-----------------------------------------------------------------\n");
  printf("M6 error\n");
  T=12.5*3600/3.;
  status=bathymetry_error(grid, topo, delta, mask, T, "M6", create);
    
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int SetRegions (const char *filename, grid_t & topogrid, unsigned short *region, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int j,status;
  int nitems,flag;
  FILE *in;
  double value;
  char toponym[64];
  grid_t grid;
  float *code, mask;
    
  pair<map<const char *,int>::iterator,bool> toponym_status;
  toponyms_t *toponyms;
  
  std::map<int, double> wdc;
  pair<map<int, double>::iterator,bool> wdc_status;
  
  toponyms=new toponyms_t;
  status=geo_init_regions(toponyms);
  
  status=grd_loadgrid (filename,&grid);
  if(status!=0) return(-1);
  
  code=new float[grid.nx*grid.ny];
  status=grd_loadr1(filename,  grid, code, &mask);
  if(status!=0) return(-1);
  
  int nprocs=initialize_OPENMP(-1);
#pragma omp parallel for private(status) if(nprocs>1)  
  for(j=0;j<topogrid.ny;j++) {
    for(int i=0;i<topogrid.nx;i++) {
      int n=j*topogrid.nx+i;
      float z;
      double t,p;
      topogrid.xy(i,j,t,p);
      t=map_recale(grid,t);
      status=map_nearestvalue00(grid, grid.nx, code, mask, t, p, &z, 0);
      if(z==mask) {
        if(debug) printf("masked toponyms at t=%lf p=%lf\n",t,p);
        region[n]=(unsigned short) 1100;
        }
      else {
        region[n]=(unsigned short) z;
        }
      }
    }
    
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 int select(grid_t & grid, int *list, unsigned short *region, unsigned short target,int level)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int n,status;
  int count=0;
  unsigned short flag,maxval=1;
  
  while(level>0) {
    maxval*=10;
    level--;
    }
  
//  target=NORTH_ATLANTIC;
  if(target>maxval) return(0);
  
  for(n=0;n<grid.Hsize();n++) {
    flag=region[n]/maxval;
    if(flag==target) {
      list[count]=n;
      count++;
      }
    }
  return(count);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int statistics (const char *filename, grid_t & topogrid, float *delta, float *base, float mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int count, i, j, n, status;
  toponyms_t *toponyms;
  unsigned short *region;
  int *list;
  
      
  for (j=0;j<topogrid.ny;j++) {
    for (i=0;i<topogrid.nx;i++) {
      int m=j*topogrid.nx+i;
      if(base[m]==mask) continue;
      if(base[m]>0.) {
        delta[m]=mask;
        }
      }
    }
  statistic_t s=get_statistics(delta, mask, topogrid.Hsize());
  for (j=0;j<topogrid.ny;j++) {
    for (i=0;i<topogrid.nx;i++) {
      int m=j*topogrid.nx+i;
      if(base[m]==mask) continue;
      if(base[m]>-2000.) {
        delta[m]=mask;
        }
      }
    }
  s=get_statistics(delta, mask, topogrid.Hsize());
  
  region=new unsigned short[topogrid.Hsize()];
  list=new int[topogrid.Hsize()];
  status=SetRegions (filename, topogrid, region, false);
  
  toponyms=new toponyms_t;
  status=geo_init_regions(toponyms);
  int level=2;

  for(toponyms_t::iterator z=toponyms->begin();z!=toponyms->end();z++) {
    count=select(topogrid, list, region, z->second,level);
    if(count==0) continue;
    printf("region %20s   \t \t %d\n",z->first,count);
    s=get_statistics(delta, mask, list, count,1);
//     char regional_output[256];
//     if(output==0) {
//       sprintf(regional_output,"%s",z->first);
//       }
//     else {
//       sprintf(regional_output,"%s.%s",output,z->first);
//       }
    }
  
  delete[] region;
  delete[] list;

  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  float  scale=1.0,topo_mask=1.e+10,mask=1.e+10;
  double  x,y,t,p;
  double  dx,dy;

  int i,j,k,l,m,n,signus=0,masked_only=0,valid_only=0,exclusion=0,status;

  size_t size;
  FILE *in,*out;
  char *keyword,*zone;
  int flag;

  char *rootname=NULL,*output=NULL,*input=NULL,*pathname=NULL,*bathymetry=NULL;
  char *bathy=NULL,*format=NULL,*hydro=NULL,*poly=NULL,*masking=NULL;

  grid_t topogrid;
  grid_t grid,zone_grid;
  mesh_t mesh;
  short *topo,smask=256*127+255;
  float *ftopo=0,*base=0,*tmp,*buffer,ftopomask;
  float value=+INFINITY, shift=+INFINITY;

  int persistence=0;

  float zmin=-INFINITY,zmax=+INFINITY;
  float *fsmoothed;
  
  plg_t *polygones=NULL;
  int npolygones=0;

  pocgrd_t ncgrid;
  cdfvar_t variable;
  bool debug=false;;
  char *filename=0;
  int order=0;

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
    if(strcmp(keyword,"-set")==0) {
      sscanf(argv[n+1],"%f",&value);
      n++;
      n++;
      continue;
      }
     if(strcmp(keyword,"-shift")==0) {
      sscanf(argv[n+1],"%f",&shift);
      n++;
      n++;
      continue;
      }
    if(strcmp(keyword,"-persistence")==0) {
      sscanf(argv[n+1],"%d",&persistence);
      n++;
      n++;
      continue;
      }
    if(strcmp(keyword,"--masked-only")==0) {
      masked_only=1;
      n++;
      continue;
      }
    if(strcmp(keyword,"--exclusion")==0) {
      exclusion=1;
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
//     if(strcmp(keyword,"-fix")==0) {
//       sscanf(argv[n+1],"%f",&zmax);
//       n++;
//       n++;
//       continue;
//       }
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

          case 'm' :
            masking= strdup(argv[n+1]);
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

  if(bathymetry==0){
    printf("base bathymetry not given (-b option)\n");
    goto error;
    }
    
/* *----------------------------------------------------------------------
  load bathymetry grid*/
  printf("#####################################################################\n");
  printf("load working bathymetry: %s\n",bathymetry);
  status=topo_loadfield(bathymetry, &topogrid, &ftopo, &ftopomask, debug);
  if(status !=0) {
    printf("cannot load bathymetry file=%s\n",bathymetry);
    goto error;
    }
  base=new float[topogrid.Hsize()];
  for(int m=0;m<topogrid.Hsize();m++) base[m]=ftopo[m];
  
  if(masking!=NULL) {
    printf("\n#####################################################################\n");
    printf("masking from polygons: %s\n", masking);
    vector<plg_t> polygons;
    status=plg_load(masking,polygons);
    
    status=map_completegridaxis(&topogrid,2);
    float *landmask=set_landmask01(topogrid, polygons);
    for(int m=0;m<topogrid.Hsize();m++) if(landmask[m]==-1) ftopo[m]=ftopomask;
    
    delete[] landmask;
    }


/*------------------------------------------------------------------------------
  convert "zero-hydro related" depths to "mean level related" depths*/
  if(hydro!=NULL) {
    printf("\n#####################################################################\n");
    printf("add/substract: %s\n", hydro);
    int positive_only=0;  /// HERE !!!
    status= topo_operation(hydro, signus, topogrid, ftopo, ftopomask,  zmin,  zmax, poly, masked_only,positive_only,0);
    if(status !=0) {
      goto error;
      }
    }
    
/*------------------------------------------------------------------------------
  set arbitrary value at selected nodes*/
  if(value!=+INFINITY) {
    printf("\n#####################################################################\n");
    printf("set value: %f\n", value);
    short *tag=0, tagvalue=0;
    status=topo_set(topogrid, ftopo, ftopomask, zmin, zmax, poly, masked_only, valid_only, exclusion, value, tag, tagvalue, debug);  
    if(status !=0) {
      goto error;
      }
    }
  
/*------------------------------------------------------------------------------
  set arbitrary value at selected nodes*/
  if(shift!=+INFINITY) {
    printf("\n#####################################################################\n");
    printf("shift value: %f (new=old+shift)\n", shift);
    short *tag=0, tagvalue=0;
    status= topo_shift(topogrid, ftopo, ftopomask, zmin, zmax, poly, masked_only, shift, tag, tagvalue, debug);  
    if(status !=0) {
      goto error;
      }
    }
  
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

  if(persistence!=0) {
    for(int k=0;k<persistence;k++) status=map_persistence(topogrid, ftopo, ftopomask, 0);
    }
    
/* *------------------------------------------------------------------------------
  save rectified topo */
  if(rootname==0) rootname=strdup("topo-operation");
  output=new char[1024];

  if(format==NULL) format=strdup("grd");

  filename=0;
  status=topo_save(rootname, 0, format, topogrid, ftopo, ftopomask, debug);
 
//   status=bathymetry_error(topogrid, base, ftopo, ftopomask);

//   status=statistics ("/home/data/landmask/polygons/geo-toponyms-18.grd", topogrid, ftopo, base, ftopomask);
#if 0  
  {
  grid_t subgrid;
  float *subtopo,*subdelta,submask;
  status=map_remap(topogrid, base, ftopomask,  &subgrid, &subtopo,  &submask, 8, 1);
  status=map_remap(topogrid, ftopo, ftopomask, &subgrid, &subdelta, &submask, 8, 1);
  int nx,ny;
  float *weight,scale=0.2;
  status=loess_filter_init(subgrid, scale, weight, nx, ny);
  fsmoothed=new float[subgrid.Hsize()];
  
  int nRequestedProcs=-1;
  int nprocs=initialize_OPENMP(nRequestedProcs);
  
#pragma omp parallel for private(status) if(nprocs>1)
  for(size_t j=0;j<subgrid.ny;j++) {
    for(size_t i=0;i<subgrid.nx;i++) {
      size_t m=j*subgrid.nx+i;
      status=loess_filter(subgrid, subdelta, submask, weight, nx, ny, i, j, fsmoothed[m]);
      }
    }
  delete[] weight;
  filename=0;
  status=topo_save("smoothed", 0, format, subgrid, fsmoothed, submask, debug);
  status=bathymetry_error(subgrid, subtopo, fsmoothed, submask);
  
  delete[] subdelta;
  delete[] subtopo;
  subgrid.free();
 
  }
#endif
  delete[] base;
  delete[] ftopo;
  topogrid.free();

  __OUT_BASE_LINE__("topo_merge sucessfully completed\n");

  exit(0);

 error:
  __OUT_BASE_LINE__("topo_merge aborted\n");
  exit(-1);
}
