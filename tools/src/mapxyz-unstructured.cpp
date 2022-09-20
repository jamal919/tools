
/*******************************************************************************

  T-UGO tools, 2006-2013

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
\brief Merges structured and random bathymetric data

<!-- A LINK TO main() or print_help() WILL NOT LINK TO THE RIGHT SOURCE ! -->
See the main function for how this works
and the print_help function for how to use this.
*/
/*----------------------------------------------------------------------------*/

#define MAIN_SOURCE

#include "version-macros.def" //for VERSION and REVISION

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "tools-structures.h"

#include "geo.h"
#include "netcdf-proto.h"
#include "grd.h"
#include "map.h"
#include "functions.h"
#include "xyz.h"
#include "topo.h"

#include "map.def"
#include "fe-proto.h"

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int check_duplicate(double *lon, double *lat, double *z, int ndata, char *keep, double resolution)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  resolution in kilometers, unsafe, to be changed

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int count=0;
  double x,y;
    
      
//   printf("#################################################################\n");
//   printf("decimating resolution %lf km\n",resolution);
  for (int m=0;m<ndata;m++) {
    double dmax=resolution;
    if(keep[m]==0) continue;
    x=lon[m];
    y=lat[m];
#pragma omp parallel for
    for (int n=m+1;n<ndata;n++) {
      if(keep[n]==0) continue;
      double xx=lon[n];
      if(fabs(xx-x)>0.001) continue;
      double yy=lat[n];
      if(fabs(yy-y)>0.001) continue;
      double d=geo_haversin_km(x,y,xx,yy);
      if(d>dmax) continue;
      printf("%d %d %lf %lf %lf %lf %lf\n",m,n,x,y,d,z[m],z[n]);
      keep[n]=0;
      }
    count++;
    }
  return(count);
}



/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template<typename T> int add(frame_t frame, double *t, double *p, int ndata, int positive, grid_t SRTMgrid, T *SRTMdata, T SRTMmask, int incr,
                             double * &x_add, double * &y_add, float * &z_add, size_t & n_add, const char *poly)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  interpolate imported DEM on future composite DEM frame (but keeping imported
  DEM native resolution)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int m, status;
  int count;
  char *buffer;
  signed char *selected;
  
  n_add=0;
  if(SRTMdata==0) return(-1);

  grid_t grid=map_Rgrid(frame, SRTMgrid.dx*incr, SRTMgrid.dy*incr, 0);
  
  buffer=new char[grid.Hsize()];
  for (size_t n=0;n<grid.Hsize();n++) buffer[n]=0;
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  flag DTM from polygon set

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  selected=new signed char[grid.Hsize()];
  status=topo_PolygonsSelection(grid, poly, selected, true, 1);

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  check if data already exist in cells, flag buffer in consequence

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  for (m=0;m<ndata;m++) {
    double a,b;
    int i,j;
    status=map_index(grid,t[m],p[m],&i,&j);
    if(status==0) {
      a=(t[m]-grid.xmin-i*grid.dx)/grid.dx;
      b=(p[m]-grid.ymin-j*grid.dy)/grid.dy;
      if(a>0.5) i++;
      if(b>0.5) j++;
      a=(t[m]-grid.xmin-i*grid.dx)/grid.dx;
      b=(p[m]-grid.ymin-j*grid.dy)/grid.dy;
      if(a*a+b*b < 0.25) {
        size_t n=i+grid.nx*j;
        buffer[n]=1;
        }
      }
    }
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  create new sounding nodes

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  count=0;
  for (size_t n=0;n<grid.Hsize();n++) if(buffer[n]==0) count++;

  x_add=new double[count];
  y_add=new double[count];
  z_add=new float [count];
  float *z=new float[grid.Hsize()];
  
  int nRequestedProcs=-1;
  int nprocs=initialize_OPENMP(nRequestedProcs,0);
  
#pragma omp parallel for private(status) if(nprocs>1)
  for(size_t j=0;j<grid.ny;j++) {
    for(size_t i=0;i<grid.nx;i++) {
      size_t n=i+grid.nx*j;
      z[n]=SRTMmask;
/*------------------------------------------------------------------------------
      if data already exist in cells, skip addition */
      if(buffer[n]==1) continue;
      if(selected[n]==1) {
        continue;
        }
      double x=grid.xmin+i*grid.dx;
      double y=grid.ymin+j*grid.dy;
      status=map_interpolation(SRTMgrid, SRTMdata, SRTMmask, x, y, &z[n]);
      }
    }
    
  count=0;
  for(size_t j=0;j<grid.ny;j++) {
    for(size_t i=0;i<grid.nx;i++) {
      size_t n=i+grid.nx*j;
      double x=grid.xmin+i*grid.dx;
      double y=grid.ymin+j*grid.dy;
      if(z[n]!=SRTMmask) {
        switch (positive) {
          case 0:
            if(z[n]<0.0) {
              x_add[count]=x;
              y_add[count]=y;
              z_add[count]=z[n];
              count++;
              }
            break;
          case 1:
//             if(z[n]>20.0) continue;
            if(z[n]>1.0) {
              x_add[count]=x;
              y_add[count]=y;
              z_add[count]=z[n];
              count++;
              }
            break;
          }
        }
      }
    }

  n_add=count;
  delete[] buffer;
  delete[] z;
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int merge(double * &x, double * &y, double * &z, int  &ndata, double *x_add, double *y_add, float *z_add, size_t n_add)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int n, nn, count;
  double *xx,*yy,*zz;
  
  count=ndata+n_add;
    
  xx=new double[count];
  yy=new double[count];
  zz=new double[count];
  
  n=0;
  for(nn=0;nn<ndata;nn++) {
    xx[n]=x[nn];
    yy[n]=y[nn];
    zz[n]=z[nn];
    n++;
    }

  for(nn=0;nn<n_add;nn++) {
    xx[n]=x_add[nn];
    yy[n]=y_add[nn];
    zz[n]=z_add[nn];
    n++;
    }

  delete[] x;
  delete[] y;
  delete[] z;
    
  x=xx;
  y=yy;
  z=zz;
  
  ndata=count;
  
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  statistic_t stats(double t, double p, grid_t & grid, float *data, float mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int m, n, status;
  int i,j;
  int range=10;
  vector<float> buffer;
  statistic_t s;
  
  status=map_index(grid,t,p,&i,&j);
  
  for(size_t jj=max(0,j-range);jj<min(grid.ny-1,j+range);jj++) {
    for(size_t ii=max(0,i-range);ii<min(grid.nx-1,i+range);ii++) {
      size_t m=grid.nx*jj+ii;
      if(data[m]!=mask) buffer.push_back(data[m]);
      }
    }
  
  s=get_statistics(buffer, mask, 0);
  
  return(s);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int xyz_filter (double* & x, double* & y, double* & z, int & ndata, grid_t BASEgrid, float *BASEdata, float BASEmask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int m, n, status;
  char *keep;
  size_t bkp;
  statistic_t s;
  double percent=0.2, lower=1.0-percent, upper=1+percent;
  double range=3.0;
  
  printf("filter from trusted bathymetry: %g%s discrepancy and %g times rms filtering\n",100.*percent,"%",range);
  
  keep=new char[ndata];
  for (m=0;m<ndata;m++) {
    keep[m]=1;
    }
  
  for(n=0;n<ndata;n++) {
    if(keep[n]==0) continue;
    double t=x[n];
    double p=y[n];
    float zz;
    t=map_recale(BASEgrid,t);
    status=map_interpolation(BASEgrid,BASEdata, BASEmask,t,p,&zz);
    if(zz < -50.) {
//       float ratio=fabs(z[n]/zz);
      float ratio=z[n]/zz;
/*------------------------------------------------------------------------------
      data within 20% range with trusted MNT                                  */
      if((ratio-lower)*(ratio-upper)<0) continue;
/*------------------------------------------------------------------------------
      have a closer look based on trusted MNT stats                           */
      s=stats(t, p, BASEgrid, BASEdata, BASEmask);
      if(fabs(z[n]-s.mean)<range*s.std) continue;
      printf("n=%d lon=%f lat=%f : sounding=%f trusted=%f %s\n",n,t,p,z[n],zz,s.print());
      keep[n]=0;
      }
    }
 
  bkp=ndata;
  
//   status=xyz_save("removed.xyz", x, y, z, keep, 0, ndata);
  
  status=xyz_reduce(x, y, z, keep, ndata);
  
  
  if(bkp-ndata!=0) {
    printf("initial count %d, after reduction %d (delta=%d)\n", bkp, ndata, bkp-ndata);
    status=xyz_save("suspicious.xyz", x, y, z, NAN, ndata);
    return(-1);
    }
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int xyz_filter (double * &x, double * &y, double * &z, int  &ndata, const char *filename, bool & sane, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  static grid_t topogrid;
  static float *topo=0, topomask;
  char *bname=0;
  int  format, status;

  if(topo==0) {
    char *s=strstr((char *) filename,".grd");
    if(s==NULL) {
      format=MAP_FORMAT_CODE_COMODO;
      bname=strdup("bathymetry");
      }
    else{
      format=MAP_FORMAT_CODE_GRD;
      }
    status=map_loadfield(filename, format, bname, &topogrid, &topo, &topomask);
    }
  
  status=xyz_filter(x, y, z, ndata, topogrid, topo, topomask);
  
  if(status!=0) sane=false;

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int filter_decimate(const char *filename, string header, char *proj4_options, float scale, float offset, 
                      double resolution, double factor, frame_t frame, const char *polygons, const char *trusted,
                      bool filter, bool paranoiac, bool duplicate, bool cartesian, bool autoscale, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int format, status;
  int m,n,ndata,count,bkp;
  double dmask;
  double *x,*y,*z;
  char *keep, *tmp;
  range_t<double> lon_range,lat_range;
  size_t found;
  string output=filename;
  frame_t plgframe;
  bool sane=true;
  
/*------------------------------------------------------------------------------
  load random data */

  tmp=strdup(filename);
  char *s=strstr(tmp,".shp");
  
  if(s==NULL) {
    format=ASCII;
    }
  else{
    format=SHAPE;
    }
  
  switch (format) {
    case ASCII:
      status=xyz_loadraw (filename, header, proj4_options, x, y, z, &dmask, ndata, debug);
      break;
    case SHAPE:
      status=xyz_load_shp (filename, proj4_options, x, y, z, &dmask, ndata);
      break;
    }
  
  if(status!=0) TRAP_ERR_EXIT(status," error by loading %s\n",filename);
  
  if(autoscale) {
    size_t positive,negative;
    for(size_t m=0;m<ndata;m++) {
      if(z[m]<0.0) negative++;
      else positive++;
      }
    if(verbose==1) printf("%d positive values, %d negative values\n",positive,negative);
    if(positive>negative) scale=-1;
    else scale=+1;
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  select data 
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  bkp=ndata;

  if(polygons!=0) {
    printf("#################################################################\n");
    printf("dataset reduction by polygons, %s)\n",polygons);
    status=xyz_PolygonSelection(polygons, x, y, z, ndata, plgframe);
    printf("initial count %d, after reduction %d (delta=%d)\n", bkp, ndata, bkp-ndata);
    }

  if(scale!=1) {
    for (m=0;m<ndata;m++) {
      if(z[m]!=dmask) z[m]*=scale;
      }
    }
  if(offset!=0.) {
    for (m=0;m<ndata;m++) {
      if(z[m]!=dmask) z[m]+=offset;
      }
    }

  keep=new char[ndata];
  for (m=0;m<ndata;m++) {
//     if(z[m]<=0) keep[m]=1;
//     else keep[m]=0;
    keep[m]=1;
    }
    
  if(frame.initialised()==1) {
    for(n=0;n<ndata;n++) {
      if(frame.inside(x[n],y[n])!=1) keep[n]=0;
      }
    }
  status=xyz_reduce(x, y, z, keep, ndata);
  
  
  if(filter) {
    string bathymetry;
    if(trusted==0) 
      bathymetry="/home/data/topography/smith-sandwell/topo_18.1.grd";
    else 
      bathymetry=(string) trusted;
    if(verbose==1) {
      printf("#################################################################\n");
      printf("dataset filtering, trusted database=%s\n",bathymetry.c_str());
      }
    status=xyz_filter(x,y,z,ndata,bathymetry.c_str(),sane,verbose);
    }
  
  if(not sane and paranoiac) {
    printf("data set %s found unsafe, do not register\n",filename);
    return(-1);
    }
  
  lon_range=range_t<double>(x,ndata);
  lat_range=range_t<double>(y,ndata);

  frame_t dataframe(lon_range,lat_range);
  if(verbose==1) dataframe.print("data range:");
  dataframe.dilatation(0.05);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  decimate data 
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if(resolution!=0) {
    for (m=0;m<ndata;m++) {
      keep[m]=0;
      }
//     int verbose=0;
    count=xyz_decimate(dataframe, x, y, z, keep, ndata, resolution, factor,cartesian, verbose);
    status=xyz_reduce(x, y, z, keep, ndata);
    }
  else {
     for (m=0;m<ndata;m++) {
      keep[m]=1;
      }
    }
    
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  check duplicates 
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if(duplicate) {
    double dmin=0.001;
    printf("#################################################################\n");
    printf("check for duplicates dmin = %lf\n",dmin);
    for (m=0;m<ndata;m++) {
      keep[m]=1;
      }
    status=check_duplicate(x, y, z, ndata, keep, dmin);
    status=xyz_reduce(x, y, z, keep, ndata);
    }
 
  delete[]keep;
  
  if(verbose==1) {
    printf("#################################################################\n");
    printf("dataset reduction, resolution= %lf, scale=%f, offset=%f (offset applies after scaling)\n",resolution,scale,offset);
    }
    
  printf("initial count %d, after reduction %d\n", bkp, ndata);
  
  found=output.find_last_of(".");
  output.insert(found,"-reduced");
  
  found=output.find_last_of("/");
  if(found!=string::npos) output.erase(0,found+1);
  
  status=xyz_save (output.c_str(), x, y, z, dmask, ndata,0);
  
  return(ndata);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 int checks(char *rootname, grid_t topogrid, float * topo, float topomask, bool gradient)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  char *output;
  pocgrd_t ncgrid;
  cdfvar_t variable;
  float *topo_x=0, *topo_y=0;
 
  size_t size=topogrid.nx*topogrid.ny;
  if(gradient) {
    float *topo_x=new float[size];
    float *topo_y=new float[size];
    status= map_gradient(topogrid, topogrid.nx, topo, topomask, GEOCENTRIC, topo_x, topo_y);
    }
  asprintf(&output,"%s.nc",rootname);
  status= poc_createfile(output);

  if(topogrid.modeH==0)  status= map_completegridaxis(&topogrid,1);
  
  status=poc_sphericalgrid_xy(output,"",topogrid,&ncgrid);
  
  poc_standardvariable_xy(&variable,"bathymetry",topomask,"m",1., 0.,"bathymetry","bathymetry","bathymetry",ncgrid);
  status=create_ncvariable(output, &variable);
  status=poc_write_xy(output,  topogrid, variable.id,topo);
  variable.destroy();
  
  if(gradient) {
    poc_standardvariable_xy(&variable,"dhdx",topomask,"dimensionless",1., 0.,"dhdx","dhdx","dhdx",ncgrid);
    status=create_ncvariable(output, &variable);
    status=poc_write_xy(output,  topogrid, variable.id,topo_x);
    variable.destroy();
  
    poc_standardvariable_xy(&variable,"dhdy",topomask,"dimensionless",1., 0.,"dhdy","dhdy","dhdy",ncgrid);
    status=create_ncvariable(output, &variable);
    status=poc_write_xy(output,  topogrid, variable.id,topo_y);
    variable.destroy();
    }
  
  free(output);
  
  deletep(&topo_x);
  deletep(&topo_y);
  
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int xyz_hydro(const char *hydro, int signus, double *x, double *y, double *z, double zmask, int ndata, int persistence)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int   i,j,m,n,status;
  float *buffer,mask,pbma;
  grid_t grid;
  double xx, yy;
  bool debug=false;
//   int ndata;

  status=topo_loadfield(hydro, &grid, &buffer, &mask, debug);
  if(status!=0) return(status);
  
  for (int k=0;k<persistence;k++) status=map_persistence(grid, buffer, mask, 0.0);
      
/*------------------------------------------------------------------------------
  Interpolate minimum low tide level and add to topo*/
  for(j=0;j<ndata;j++) {
    if(z[j]==zmask) {
      continue;
      }
    xx=x[j];
    yy=y[j];
    status=map_interpolation(grid, buffer,mask,xx,yy,&pbma);
    if (pbma!=mask) {
/*------------------------------------------------------------------------------
      assume negative depths and negative lowest astronomical tide*/
      z[j]+=signus*pbma;
      }
    }

  delete[] buffer;

  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int cloud_limits(const string & rootname, double *x, double *y, double *z, int ndata, double ElementMaxSize, double ElementMaxRatio, vector<plg_t> & p, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  size_t n;
  mesh_t mesh, delaunay;
  triangulateio in, out;
  bool meshout=false;
  string output;
  char *keep=0;
  
  keep=new char[ndata];
  for(size_t m=0;m<ndata;m++) {
    keep[m]=0;
    }
    
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  build FE delaunay mesh

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  printf("#################################################################\n");
  printf("generate unstructured mesh\n\n");

  double resolution=1.0, factor=20.;
  bool cartesian=true;
  
  range_t<double> lon_range=range_t<double>(x,ndata);
  range_t<double> lat_range=range_t<double>(y,ndata);
  frame_t frame(lon_range,lat_range);
  
  int verbose=1;
  size_t count;
  
  count=xyz_decimate(frame, x, y, z, keep, ndata, resolution, factor, cartesian, verbose);
//   status=xyz_CheckDuplicate(x, y, z, ndata, keep, resolution, (projPJ) 0);

  status=xyz_reduce(x, y, z, keep, ndata);

  mesh.nvtxs=ndata;
  mesh.vertices=new vertex_t[mesh.nvtxs];
  
  for(size_t n=0;n<ndata;n++) {
//     if(keep[nn]==1) {
      mesh.vertices[n].lon=x[n];
      mesh.vertices[n].lat=y[n];
      mesh.vertices[n].h=z[n];
//       }
    }
    
  status=fe_delaunay_create(rootname.c_str(), mesh, delaunay, meshout, debug);
  mesh.destroy();
  
  ElementMaxSize=1.e+4;
  ElementMaxRatio=10.;
  status=fe_delaunay_improve(rootname.c_str(), delaunay, ElementMaxSize, ElementMaxRatio, debug);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  finalize mesh and get boundary polygons

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  status=fe_edgetable(&delaunay, debug);
  
  int RecycleCodes=0, SetLimits=1;
  
  verbose=0;
  
  status=fe_codetable(delaunay, RecycleCodes, SetLimits, verbose);
  status=fe_limits2poly(delaunay, p, (char *) 0, true);
  
//   if(debug) {
  output=rootname+"delaunay.nei";
  status=fe_savemesh(output.c_str(),MESH_FILE_FORMAT_TRIGRID,delaunay);
//     }
  
  status=plg_save("delaunay.plg", PLG_FORMAT_SCAN, p);

  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  mesh_t *xyz_loadmap_UG(const char *filename, char *rootname, char *bathymetry, vector<string> & elevation,  string PBMA, string distance_to_coast, char *proj4_options, 
                        float scale, float dilatation, float *offset, double resolution, double factor, double grid_resolution,
                        int DecimateNone, int BASEincr,int SRTMincr, double ElementMaxSize, double ElementMaxRatio, 
                        frame_t focus, const char *e_poly, const char *polygons, frame_t mapping, bool CheckDuplicates, bool meshout, bool cartesian, bool isobath, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int positive;
  int k,m,n,ndata,bkp,nn, format;
  size_t count;
  double dmask;
  mesh_t mesh, *finished;
  triangulateio in, out;
  char *keep;
  range_t<double> lon_range,lat_range;
  double *x,*y,*z,d;
  grid_t SRTMgrid;
  float *SRTMdata=0, SRTMmask;
  grid_t DISgrid;
  float *DISdata=0, DISmask;
  grid_t sTOPOgrid;
  grid_t TOPOgrid;
  float *TOPOdata=0, TOPOmask=1.e+10;
  grid_t BASEgrid;
  float *BASEdata=0, BASEmask=1.e+10;
  double *x_add, *y_add;
  float  *z_add;
  size_t n_add;
  char output[1024];
  float d_offset=offset[0], b_offset=offset[1], e_offset=offset[2];
  char *bname=0;
  frame_t plgframe;
  string header="";
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  load gridded bathymetry

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if(bathymetry!=0) {
    char *s=strstr(bathymetry,".grd");
    if(s==NULL) {
      format=MAP_FORMAT_CODE_COMODO;
      bname=strdup("bathymetry");
      }
    else{
      format=MAP_FORMAT_CODE_GRD;
      }
    status=map_loadfield(bathymetry, format, bname, &BASEgrid, &BASEdata, &BASEmask);
    if(status!=0){
      bname=strdup("elevation");
      status=map_loadfield(bathymetry, format, bname, &BASEgrid, &BASEdata, &BASEmask);
      }
    if(status!=0){
       TRAP_ERR_EXIT(-1,"cannot load bathymetry database %s, exiting\n", bathymetry);
      }
    if(b_offset!=0.)
    for (m=0;m<BASEgrid.Hsize();m++) {
      if(BASEdata[m]!=BASEmask) BASEdata[m]+=b_offset;
      }
    }
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  load distance_to_coast

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if(distance_to_coast!="") {
    status=topo_loadfield(distance_to_coast.c_str(), &DISgrid, &DISdata, &DISmask, debug);
    if(status!=0){
       TRAP_ERR_EXIT(-1,"cannot load elevation database, exiting\n");
      }
    }
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  load soundings cloud data 

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  char *tmp=strdup(filename);
  char *s=strstr(tmp,".shp");
  if(s==NULL) {
    format=ASCII;
    }
  else{
    format=SHAPE;
    }
  
  switch (format) {
    case ASCII:
      status=xyz_loadraw (filename, header, proj4_options, x, y, z, &dmask, ndata);
      break;
    case SHAPE:
      status=xyz_load_shp (filename, proj4_options, x, y, z, &dmask, ndata);
      break;
    }
    
  if(status!=0) {
     TRAP_ERR_EXIT(-1,"cannot load souding database, exiting\n");
    }
    
//     sprintf(output, "%s-dump.xyz",rootname);
//     status=xyz_save (output, x, y, z, dmask, ndata);

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  frame selection of depth soundings cloud
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/  

  if(focus.initialised()==1) {
    keep=new char[ndata];
    for (m=0;m<ndata;m++) {
      keep[m]=1;
      }
    bkp=ndata;
    for(n=0;n<ndata;n++) {
      if(focus.inside(x[n],y[n])!=1) keep[n]=0;
      }
    status=xyz_reduce(x, y, z, keep, ndata);
    printf("#################################################################\n");
    printf("dataset reduction by frame : \n");
    printf("initial count %d, after reduction %d\n", bkp, ndata);
    delete[] keep;
    }
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  mapping frame selection of depth soundings cloud
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/  

  if(mapping.initialised()) {
    frame_t frame=mapping;
    frame.dilatation(0.05);
    status=xyz_FrameSelection(frame, x, y, z, ndata);
    }
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  polygon selection of depth soundings cloud
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/  

  if(polygons!=0) {
    bkp=ndata;
    status=xyz_PolygonSelection(polygons, x, y, z, ndata, plgframe);
    printf("#################################################################\n");
    printf("dataset reduction by polygons : %s\n",polygons);
    printf("initial count %d, after reduction %d\n", bkp, ndata);
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  depths correction (sign, vertical reference offset)
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/  

  if(scale!=1) {
    for (m=0;m<ndata;m++) {
      z[m]*=scale;
      }
    }
  if(d_offset!=0.) {
    for (m=0;m<ndata;m++) {
      z[m]+=d_offset;
      }
    }
    
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  convert "zero-hydro related" depths to "mean level related" depths
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/  

//   if(PBMA!="") {
//     int signus=1;
// //     if(negative) signus=+1;
// //     else         signus=-1;
//     status= xyz_hydro(PBMA.c_str(), signus, x, y, z, dmask, ndata, 5);
//     sprintf(output, "%s-meanlevel.xyz",rootname);
//     status=xyz_save (output, x, y, z, dmask, ndata);
//     }
  
  if(debug) {
    sprintf(output, "%s-dump.xyz",rootname);
    status=xyz_save (output, x, y, z, dmask, ndata);
    }
    
//   vector<plg_t> p;
//   status=cloud_limits((const string) rootname, x, y, z, ndata, ElementMaxSize, ElementMaxRatio, p, debug);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  establish completion frame
  
  - if mapping frame not specified, addition frame is a limited dilatation of 
    data frame
  
  - if mapping frame specified, addition frame is a limited dilatation of 
    mapping frame
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/  

/*------------------------------------------------------------------------------
  compute data points geographic range */
  lon_range=range_t<double>(x,ndata);
  lat_range=range_t<double>(y,ndata);

  frame_t dataframe(lon_range,lat_range);
  dataframe.dilatation(dilatation);
    
/*------------------------------------------------------------------------------
  this is an issue when data outside of mapping */
  if(not mapping.initialised()) mapping=dataframe;
  
/*------------------------------------------------------------------------------
  keep cloud extent sligthly larger than mapping one */
//   frame_t addition_frame(lon_range,lat_range);
  frame_t addition_frame=dataframe;
  addition_frame.dilatation(0.05);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  add elevation (dry) data : lidar, SRTM, topography DTM
  
  multi-dataset implementation.
  
  Filtering: positive (z>1m) elevation, optionally distance to coast
  
  distance to coast dataset : http://www.soest.hawaii.edu/wessel/gshhg/

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/  

  double dmax=5.0;
  n_add=0;
  if(elevation.size()!=0) {
    positive=1;
    for(int k=0;k<elevation.size();k++) {
      status=topo_loadfield(elevation[k].c_str(), &SRTMgrid, &SRTMdata, &SRTMmask, debug);
      if(status!=0){
         TRAP_ERR_EXIT(-1,"cannot load elevation database, file=%s, exiting\n", elevation[k].c_str());
        }
      if(SRTMgrid.modeH==1) {
        SRTMgrid.dx=SRTMgrid.x[1]-SRTMgrid.x[0];
        SRTMgrid.dy=SRTMgrid.y[1]-SRTMgrid.y[0];
        }
      if(e_offset!=0.)
      for (m=0;m<SRTMgrid.Hsize();m++) {
        if(SRTMdata[m]!=SRTMmask) SRTMdata[m]+=e_offset;
        }
      printf("#################################################################\n");
      printf("add elevation nodes (positive elevations only) from %s; resolution %9.1lf m \n",elevation[k].c_str(),SRTMincr*SRTMgrid.dx*110000.);
      addition_frame.print("addition frame:");
      int increment=max(1, (int) NINT(grid_resolution/SRTMgrid.dx));
      status= add(addition_frame, x, y, ndata, positive, SRTMgrid, SRTMdata, SRTMmask, increment, x_add, y_add, z_add, n_add, e_poly);
      printf("number of elevation nodes added : %d\n",n_add);
      if(n_add!=0) {
        status=merge(x, y, z, ndata, x_add, y_add, z_add, n_add);
        delete[] x_add;
        delete[] y_add;
        delete[] z_add;
        }
/*------------------------------------------------------------------------------
      filter elevation data following distance to coast indication */
      if(DISdata!=0) {
#pragma omp parallel for
        for (int j=0;j<SRTMgrid.ny;j++) {
          for (int i=0;i<SRTMgrid.nx;i++) {
            double x,y;
            float d;
            SRTMgrid.xy(i,j,x,y);
            status=map_interpolation(DISgrid, DISdata, DISmask, x, y, &d);
            if(d>dmax) SRTMdata[j*SRTMgrid.nx+i]=SRTMmask;
            }
          }
        }
      status= add(addition_frame, x, y, ndata, positive, SRTMgrid, SRTMdata, SRTMmask, SRTMincr, x_add, y_add, z_add, n_add, e_poly);
      printf("number of elevation nodes added : %d\n",n_add);
      if(debug) {
        sprintf(output, "%s-elevation.xyz", rootname);
        status=xyz_save (output, x_add, y_add, z_add, dmask, n_add);
        }
      if(n_add!=0) {
        status=merge(x, y, z, ndata, x_add, y_add, z_add, n_add);
        delete[] x_add;
        delete[] y_add;
        delete[] z_add;
        }
      SRTMgrid.free();
      deletep(&SRTMdata);
      }
    }
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  add bathymetry data
    
  Filtering: negative (z<0m) elevation
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/  

  n_add=0;
  if(bathymetry!=0) {
    printf("#################################################################\n");
    printf("add bathymetry nodes (negative elevations only) from %s; resolution %9.1lf m \n",bathymetry,BASEincr*BASEgrid.dx*110000.);
    addition_frame.print("addition frame:");
    positive=0;
    status= add(addition_frame, x, y, ndata, positive, BASEgrid, BASEdata, BASEmask, BASEincr, x_add, y_add, z_add, n_add, (const char*) 0);
    printf("number of bathymetry nodes added : %d\n",n_add);
    if(debug) {
      sprintf(output, "%s-bathymetry.xyz",rootname);
      status=xyz_save (output, x_add, y_add, z_add, dmask, n_add);
      }
    }
  
  keep=new char[ndata];
  for (m=0;m<ndata;m++) {
    keep[m]=0;
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  decimation : gather data by cells (clustering)
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
//   if(resolution==-1) {
//     resolution=0.025;
//     printf("\ndecimation resolution not set, use automatic value (%lfkm); can be set with -resolution XXX\n\n", resolution);
//     }

  if(DecimateNone==0) {
    printf("#################################################################\n");
    printf("decimate node set\n");
    double r=resolution;
    if(r==0.0) r=grid_resolution/2.0;
    int verbose=1;
    count=xyz_decimate(mapping, x, y, z, keep, ndata, r, factor, cartesian, verbose);
    }
  else {
    for (m=0;m<ndata;m++) {
      keep[m]=1;
      }
    count=ndata;
/*------------------------------------------------------------------------------
    check duplicates */
    if(CheckDuplicates) {
      status=check_duplicate(x, y, z, ndata, keep, 0.001);
      status=xyz_reduce(x, y, z, keep, ndata);
      }
   }
    
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  Remove deep depth data if misfit with bathymetry DTM higher than 20%
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
#if 0
  size_t check;
  check=occurence((char) 0, keep, ndata);
  if(bathymetry!=0) {
    for(nn=0;nn<ndata;nn++) {
      if(keep[nn]==0) continue;
      double t=x[nn];
      double p=y[nn];
      float zz;
      t=map_recale(BASEgrid,t);
      status=map_interpolation(BASEgrid,BASEdata, BASEmask,t,p,&zz);
      if(zz==BASEmask) continue;
      if(zz<-100.) {
        float ratio=fabs(z[nn]/zz);
        if((ratio-0.8)*(ratio-1.2)>0) keep[nn]=0;
        }
      }
    }
  check=occurence((char) 0, keep, ndata);
#endif

  BASEgrid.free();  
  if(BASEdata!=0) delete[] BASEdata;
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  build FE mesh

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  printf("#################################################################\n");
  printf("generate unstructured mesh\n\n");
  
  mesh.nvtxs=count+n_add;
  mesh.vertices=new vertex_t[mesh.nvtxs];
  
  n=0;
  for(nn=0;nn<ndata;nn++) {
    if(keep[nn]==1) {
      mesh.vertices[n].lon=x[nn];
      mesh.vertices[n].lat=y[nn];
      mesh.vertices[n].h=-z[nn];
      n++;
      }
    }
  for(nn=0;nn<n_add;nn++) {
    mesh.vertices[n].lon=x_add[nn];
    mesh.vertices[n].lat=y_add[nn];
    mesh.vertices[n].h =-z_add[nn];
    n++;
    }
  
  delete[] x;
  delete[] y;
  delete[] z;
  
  if(n_add!=0) {
    delete[] x_add;
    delete[] y_add;
    delete[] z_add;
    }
  
  mesh.nvtxs=n;
  
  if(debug) status=fe_savenodes("out.node", NODE_FILE_FORMAT_TRIANGLE, mesh);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  initialize triangle output structure and create triangulation

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  status= fe_triangulateio_init(in);
  status= fe_triangulateio_init(out);
  
  status= fe_node2triangle(mesh, &in);
  mesh.destroy();
  
  triangulate("z", &in, &out, (triangulateio *) NULL);

/*------------------------------------------------------------------------------
  output structure inherit from input structure's pointlist, to handle with care */
  out.pointlist = in.pointlist;

  printf("#################################################################\n");
  printf("convert triangle mesh into sirocco mesh\n");
  finished=new mesh_t;

  status=fe_triangle2mesh(out,finished);
  
/*------------------------------------------------------------------------------
  output structure inherit from input structure's pointlist, to handle with care */
  in.pointlist = out.pointlist;
  fe_FreeTriangleIO(in);
  out.pointlist = 0;
  fe_FreeTriangleIO(out);
  
  if(debug) {
    sprintf(output, "%s-01.nei",rootname);
    status=fe_savemesh(output,MESH_FILE_FORMAT_TRIGRID,*finished);
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  post-processing (mostly for node's cloud stand-alone mesh)

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  filter elements by size (to clean Delaunay mesh), pretty unsafe
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  int filter=(ElementMaxSize>0);
  if(filter==1) {
    printf("#################################################################\n");
    printf("remove oversized element connexion: maxsize=%lf m\n",ElementMaxSize);
    count=0;
    for(n=0;n<finished->nvtxs;n++) {
      for(k=0;k<finished->vertices[n].nngh;k++) {
        m=finished->vertices[n].ngh[k];
        if(m<n) continue;
/*------------------------------------------------------------------------------
        warning : fe_distance now return distance in meters */
        d=fe_distance(*finished, n, m);
        if(d>ElementMaxSize) {
//       if(verbose==1) printf("disconnect vertices %d %d\n",n,m);
          status= fe_disconnectvertices(*finished, n,m);
          k=k-1;
          count++;
          }
        }
      }
    printf("remove oversized edges, discarded connections=%d\n", count);
    status=fe_cleanvertices(finished, false);
    if(count!=0) {
      if(finished->triangles!=0) finished->triangles->destroy();
      if(finished->edges!=0)     finished->edges->destroy();
      finished->triangles=0;
      finished->edges=0;
      }
    if(finished->triangles==0) {
      status=fe_list(finished);
      status=fe_e2n(finished);
      }
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  filter elements by size ratio (to clean Delaunay mesh), pretty unsafe
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  filter=(ElementMaxRatio!=0);
  if(filter==1) {
    printf("#################################################################\n");
    printf("remove oversized element connexion: maxsize=%lf m\n",ElementMaxSize);
    int nboundaries, verbose=0;
//     status=fe_initcodes(*finished, nboundaries, verbose);
    count=0;
    for(n=0;n<finished->nvtxs;n++) {
//       if(finished->vertices[n].code==0) continue;
      double dmin=+1.e+10;
      double dmax=-1.e+10;
      vector<double> lengths;
      if(finished->vertices[n].nngh<5) continue;
      for(k=0;k<finished->vertices[n].nngh;k++) {
        m=finished->vertices[n].ngh[k];
//         if(finished->vertices[m].code==0) continue;
        d=fe_distance(*finished, n, m);
        lengths.push_back(d);
        dmin=min(dmin,d);
        dmax=max(dmax,d);
        }
      for(k=0;k<finished->vertices[n].nngh;k++) {
        m=finished->vertices[n].ngh[k];
        if(m<n) continue;
        d=fe_distance(*finished, n, m);
/*------------------------------------------------------------------------------
        following check should be limited to boundary edges */
        if(d>ElementMaxRatio*dmin) {
          status= fe_disconnectvertices(*finished, n,m);
          k=k-1;
          count++;
          }
        }
      }
    printf("remove unbalanced elements, discarded connections=%d\n", count);
    status=fe_cleanvertices(finished, false);
    if(count!=0) {
      if(finished->triangles!=0) finished->triangles->destroy();
      if(finished->edges!=0)     finished->edges->destroy();
      finished->triangles=0;
      finished->edges=0;
      }
    if(finished->triangles==0) {
      status=fe_list(finished);
      status=fe_e2n(finished);
      }
    }
    
  if(isobath) {

  printf("#################################################################\n");
  printf("construct face tables\n");
  
  finished->nedges=0;
  for(n=0;n<finished->nvtxs;n++) {
    for(k=0;k<finished->vertices[n].nngh;k++) {
      if(n<finished->vertices[n].ngh[k]) finished->nedges++;
      }
    }
   
  finished->faces=new face_t[finished->nedges];
  finished->nedges=0;
  for(n=0;n<finished->nvtxs;n++) {
    for(k=0;k<finished->vertices[n].nngh;k++) {
      if(n<finished->vertices[n].ngh[k]) {
        m=finished->nedges;
        finished->faces[m].extremity[0]=n;
        finished->faces[m].extremity[1]=finished->vertices[n].ngh[k];
        finished->nedges++;
        }
      }
    }
   
  if(debug) {
    sprintf(output, "%s-02.nei",rootname);
    status=fe_savemesh(output,MESH_FILE_FORMAT_TRIGRID,*finished);
    }
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  improve vertex connection to favour isobaths
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  printf("#################################################################\n");
  printf("improve node connections\n");
  count=0;
  for(n=0;n<finished->nedges;n++) {
    int n1,n2,q[4];
    double dh[2],angle;
    n1=finished->faces[n].extremity[0];
    n2=finished->faces[n].extremity[1];
/*------------------------------------------------------------------------------
    form direct order quadrangle*/
    status=fe_quadnodes(*finished, n1, n2, q);
    if(status!=0) continue;
/*------------------------------------------------------------------------------
    check convexity*/
    for(k=0;k<4;k++) {
      angle=fe_angle(*finished,q[k],q[(k+1)%4],q[(k+3)%4]);
      if(angle<M_PI/6.) {
        status=-1;
        }
      }
    if(status!=0) continue;
/*------------------------------------------------------------------------------
    check aspect ratio*/
    k=1;
    angle=fe_angle(*finished,q[k],q[(k+1)%4],q[(k+2)%4]);
    if(angle<M_PI/3.) continue;
    k=3;
    angle=fe_angle(*finished,q[k],q[(k+1)%4],q[(k+2)%4]);
    if(angle<M_PI/3.) continue;
/*------------------------------------------------------------------------------
    check elevation connexion*/
    dh[0]=fabs(finished->vertices[q[0]].h-finished->vertices[q[2]].h);
    dh[1]=fabs(finished->vertices[q[1]].h-finished->vertices[q[3]].h);
    if(dh[0]>dh[1]) {
      status=fe_disconnectvertices(*finished, n1, n2);
      if(status!=0) {
        printf("%s : trouble\n",__func__);
        }
      n1=q[1];
      n2=q[3];
      status=fe_connectvertices(*finished, n1, n2);
      if(status!=0) {
        printf("%s : trouble\n",__func__);
        }
      finished->faces[n].extremity[0]=n1;
      finished->faces[n].extremity[1]=n2;
      count++;
      }
    }

  delete [] finished->faces;
  }
  
  printf("#################################################################\n");
  printf("construct final topography mesh\n");
  finished->nnghm=0;
  for(n=0;n<finished->nvtxs;n++) {
    updatemax(&finished->nnghm,finished->vertices[n].nngh);
    }
  if(finished->edges!=0) finished->edges->destroy();
  
  status=fe_list(finished);

  if(meshout) {
    sprintf(output, "%s.nei",rootname);
    printf("save topography mesh in %s \n",output);
    status=fe_savemesh(output,MESH_FILE_FORMAT_TRIGRID,*finished);
    }
    
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  create bathymetry structured grid 
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  printf("#################################################################\n");
  printf("create bathymetry on structured grid: dx=%lf (%lf m) arc-sec dy=%lf (%lf m) arc-sec\n",
         grid_resolution*3600., grid_resolution*110000., grid_resolution*3600., grid_resolution*110000.);
    
  double dx=grid_resolution;
  double dy=grid_resolution;
  
  TOPOgrid=map_Rgrid(mapping, dx, dy, 0);
  map_printgrid(TOPOgrid);
  
  TOPOdata=new float[TOPOgrid.Hsize()];
  TOPOmask=SRTMmask;
  
  TOPOmask=BASEmask;
  
  int algorithm=0, verbose=1;
  
  status=fe_mapVdepth(*finished, TOPOgrid, TOPOdata, TOPOmask, algorithm, verbose);
  
  status=checks(rootname,TOPOgrid,TOPOdata,TOPOmask,false);
  

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  save data, to be replaced with topo_save. Example:
  
  status=topo_save(rootname, output, "netcdf", TOPOgrid, TOPOdata, TOPOmask, debug);
  
  topo_save will hande mirroring...

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  status=grd_mirror_r( TOPOgrid, TOPOgrid.nx, TOPOdata, TOPOmask);
  sprintf(output, "%s.grd",rootname);
  status=grd_save (output,TOPOgrid, TOPOgrid.nx, TOPOdata, TOPOmask);

  STDOUT_BASE_LINE("Finished with %s\n",__func__);
  
  TOPOgrid.free();
  delete[] TOPOdata;
  
  return(finished);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void print_help(const char *prog_name)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/** \brief prints help of this programme.

\param prog_name name of the program to be printed by this help function : argv[0]
*/
/*----------------------------------------------------------------------------*/
{
/** The code of the body of this function is :
    \code /**/ // COMPILED CODE BELOW !!!
  printf("\n"
    "NAME AND VERSION\n"
    "  %s version " VERSION " " REVISION "\n",prog_name);
  printf("\n"
    "USE\n"
    "  %s [OPTIONS] input\n",prog_name);
  printf("\n"
    "DESCRIPTION\n"
/* USE FIND AND REPLACE TO EDIT THE LINE BELOW, SO THAT THE FILE HEAD IS ALSO UPDATED */
    "  Merges structured and random bathymetric data.\n"
    "  It created a bathymetry from the random input, padded with the base bathymetry and the elevation.\n"
    "  Unless told not to, it will decimate the random data by averaging clusters.\n"
    "\n"
    "OPTIONS\n"
    "  -h,--help  Show this help and exit.\n"
    "  -scale,-offset scale (default:1.) and offset (default:0.) for the random data. Usefull for reference level and depth to altitude conversions\n"
    "  -b_offset,-b_incr offset (default:0.) and point incrementation (default:1) of the base bathymetry\n"
    "  -e_offset,-e_incr offset (default:0.) and point incrementation (default:1) of the elevation\n"
    "  -resolution followed by the resolution of the decimation (default:1/40deg)\n"
    "  -grid followed by the output resolution (default:1/3600deg)\n"
    "  -b  followed by the base bathymetry\n"
    "  -e  followed by the elevations\n"
    "  -o  followed by the root name of the files created. Default : mapxyz-unstructured\n"
    "  -dilatation : followed by\n"
    "  --decimate-only : \n"
    "  --no-decimation : \n"
    "  -maxsize : \n"
    "  -e_poly : \n"
    "  --filter : \n"
    "  --debug : \n"
    "  --mesh : \n"
    "  --check-duplicate : \n"
    "  -p : followed by a selection (inner-wise) polygons set filename\n"
    );
  /** \endcode */
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int n,status;

  char *keyword;

  char *rootname=NULL, *bathymetry=NULL, *trusted=NULL;
  vector<string> elevation;
  vector<string> input;
  string distance_to_coast;
  
  char *polygons=NULL;
  char *e_poly=NULL;
  char *proj4=NULL;
  
  string FrameString="",MappingString="",HeaderString="",PBMA="",list="";

  double resolution,grid_resolution;
  grid_t topogrid,cgrid;
  grid_t grid;
  mesh_t mesh;
  float scale=1.0, dilatation=0.05;
  float b_offset=0.,e_offset=0.,d_offset=0.;
  double factor;

  pocgrd_t ncgrid;
  cdfvar_t variable;
  
  int DecimateOnly=0, DecimateNone=0;
  bool filter=false, paranoiac=false, duplicate=false, debug=false, meshout=false, isobath=true;
  bool cartesian=false, autoscale=false;
  double ElementMaxSize=-1, ElementMaxRatio=0.0;
  
//   frame_t frame, prescribed(-1.e+10,1.e+10,-1.e+10,1.e+10), mapping;
  frame_t frame, prescribed, mapping;

  fct_echo( argc, argv);

//   status=xyz_load_shp("/media/8d0f6fd3-80b8-4696-b618-6f74efb78e3a/data/topography/regional/dalles-SHOM/2013/XYCarreau-14581-20090528-25.shp");

  resolution=-1.0;
  grid_resolution=-1;
  factor=-1;
  
  int BASEincr=1;
  int SRTMincr=1;
  
  int nitems;
  
  n=1;
  while (n < argc) {
    keyword=strdup(argv[n]);
    if(strcmp(keyword,"--help")==0 or strcmp(keyword,"-h")==0){
      print_help(argv[0]);
      exit(0);
      }
    if(strcmp(keyword,"-scale")==0) {
      nitems=sscanf(argv[n+1],"%f",&scale);
      if(nitems!=1) TRAP_ERR_EXIT(-1,"*** missing numerical value, expected after %s ***\n",keyword);
      n++;
      n++;
      continue;
      }
    if(strcmp(keyword,"--autoscale")==0) {
      autoscale=true;
      n++;
      continue;
      }
    if(strcmp(keyword,"-dilatation")==0) {
      nitems=sscanf(argv[n+1],"%f",&dilatation);
      if(nitems!=1) TRAP_ERR_EXIT(-1,"*** missing numerical value, expected after %s ***\n",keyword);
      n++;
      n++;
      continue;
      }
    if(strcmp(keyword,"-maxsize")==0) {
      nitems=sscanf(argv[n+1],"%lf",&ElementMaxSize);
      printf("\n\n\n in addition to maxsize, please use -maxratio 5.0 to retrieve former versions behaviour \n\n\n");
      if(nitems!=1) TRAP_ERR_EXIT(-1,"*** missing numerical value, expected after %s ***\n",keyword);
      n++;
      n++;
      continue;
      }
    if(strcmp(keyword,"-maxratio")==0) {
      nitems=sscanf(argv[n+1],"%lf",&ElementMaxRatio);
      if(nitems!=1) TRAP_ERR_EXIT(-1,"*** missing numerical value, expected after %s ***\n",keyword);
      n++;
      n++;
      continue;
      }
    if(strcmp(keyword,"-offset")==0) {
      nitems=sscanf(argv[n+1],"%f",&d_offset);
      if(nitems!=1) TRAP_ERR_EXIT(-1,"*** missing numerical value, expected after %s ***\n",keyword);
      n++;
      n++;
      continue;
      }
    if(strcmp(keyword,"-b_offset")==0) {
      nitems=sscanf(argv[n+1],"%f",&b_offset);
      if(nitems!=1) TRAP_ERR_EXIT(-1,"*** missing numerical value, expected after %s ***\n",keyword);
      n++;
      n++;
      continue;
      }
    if(strcmp(keyword,"-b_incr")==0) {
      nitems=sscanf(argv[n+1],"%d",&BASEincr);
      if(nitems!=1) TRAP_ERR_EXIT(-1,"*** missing numerical value, expected after %s ***\n",keyword);
      n++;
      n++;
      continue;
      }
     if(strcmp(keyword,"-e_offset")==0) {
      nitems=sscanf(argv[n+1],"%f",&e_offset);
      if(nitems!=1) TRAP_ERR_EXIT(-1,"*** missing numerical value, expected after %s ***\n",keyword);
      n++;
      n++;
      continue;
      }
    if(strcmp(keyword,"-e_incr")==0) {
      nitems=sscanf(argv[n+1],"%d",&SRTMincr);
      if(nitems!=1) TRAP_ERR_EXIT(-1,"*** missing numerical value, expected after %s ***\n",keyword);
      n++;
      n++;
      continue;
      }
    if(strcmp(keyword,"-e_poly")==0) {
      e_poly=strdup(argv[n+1]);
      n++;
      n++;
      continue;
      }
    if( strcmp(keyword,"--isobath=true")==0 ){
      isobath=true;
      n++;
      continue;
      }
    if( strcmp(keyword,"--isobath=false")==0 ){
      isobath=false;
      n++;
      continue;
      }
    if( strcmp(keyword,"--cartesian")==0 ){
      cartesian=true;
      n++;
      continue;
      }
    if( strcmp(keyword,"--pbma")==0 ){
      PBMA=argv[n+1];
      n++;
      n++;
      continue;
      }
    if(strcmp(keyword,"--decimate-only")==0) {
      DecimateOnly=1;
      n++;
      continue;
      }
    if(strcmp(keyword,"--no-decimation")==0) {
      DecimateNone=1;
      n++;
      continue;
      }
    if(strcmp(keyword,"--filter")==0) {
      filter=true;
      n++;
      continue;
      }
    if(strcmp(keyword,"--paranoiac")==0) {
      paranoiac=true;
      n++;
      continue;
      }
   if(strcmp(keyword,"--debug")==0) {
      debug=true;
      n++;
      continue;
      }
   if(strcmp(keyword,"--mesh")==0) {
      meshout=true;
      n++;
      continue;
      }
   if(strcmp(keyword,"--check-duplicate")==0) {
      duplicate=true;
      n++;
      continue;
      }
    if(strcmp(keyword,"-resolution")==0) {
      nitems=sscanf(argv[n+1],"%lf",&(resolution));
      n++;
      n++;
      continue;
      }
    if(strcmp(keyword,"-grid")==0) {
      nitems=sscanf(argv[n+1],"%lf",&(grid_resolution));
      n++;
      n++;
      continue;
      }
    if(strcmp(keyword,"-cluster")==0) {
      nitems=sscanf(argv[n+1],"%lf",&factor);
      n++;
      n++;
      continue;
      }
    if(strncmp("--frame",keyword)==0){
      FrameString=argv[n+1];
      n++;
      n++;
      continue;
      }
    if(strncmp("--mapping",keyword)==0){
      MappingString=argv[n+1];
      n++;
      n++;
      continue;
      }
    if(strncmp("--header",keyword)==0){
      HeaderString=argv[n+1];
      n++;
      n++;
      continue;
      }
    if(strstr(argv[n],"--proj=")!=0){
      proj4=strdup(argv[n]+7);
      n++;
      continue;
      }
    if(strstr(argv[n],"-trusted")!=0){
      trusted=strdup(argv[n+1]);
      n++;
      n++;
      continue;
      }

    switch (keyword[0]) {
      case '-':
        switch (keyword[1]) {
        case 'b' :
          bathymetry=strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'd' :
          distance_to_coast=argv[n+1];
          n++;
          n++;
          break;

        case 'e' :
          elevation.push_back(argv[n+1]);
          n++;
          n++;
          break;

        case 'l' :
          list=argv[n+1];
          n++;
          n++;
          break;

        case 'o' :
          rootname=strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'p' :
          polygons=strdup(argv[n+1]);
          n++;
          n++;
          break;

        default:
          STDOUT_BASE_LINE("unknown option %s\n",keyword);
          exit(-1);
        }
        break;

      default:
        input.push_back(argv[n]);
        n++;
        break;
      }
      free(keyword);
    }

  if(list!="") {
    status=append_filelist(list, input);
    }
  
  if(input.size()==0) {
    STDERR_BASE_LINE("*** input random bathymetric data needed ***\n");
    print_help(argv[0]);
    exit(-1);
    }
  
  if(resolution==-1) {
    resolution=0.025;
    printf("\ndecimation resolution not set, use default value (%lfkm); can be set with -resolution XXX\n\n", resolution);
    }
  
  if(factor==-1) {
    factor=20.0;
    printf("\ndecimation cluster factor not set, use default value (%lf); can be set with -cluster XXX\n\n", factor);
    }
  
  if(grid_resolution==-1) {
    grid_resolution=1./3600.;
    printf("\nstructured grid resolution not set, use default value (%lf arcsec); can be set with -grid XXX\n\n", grid_resolution*3600.);
    }
  
  if(rootname==0) rootname=strdup("mapxyz-unstructured");
  
  if(FrameString!="") {
    status=plg_DecodeFrame(FrameString, prescribed);
    }
  
  if(MappingString!="") {
    status=plg_DecodeFrame(MappingString, mapping);
    if(status!=0) TRAP_ERR_EXIT(-1,"mapping string not understood (may be you forgot \"): %s\n",MappingString.c_str());
    if(debug) mapping.print("mapping frame:");
    }
  
  if(trusted!=0) filter=true;
  
  printf("#################################################################\n");
  printf("\ndecimation resolution=%lf km;", resolution);
  printf(" structured mapping resolution=%lf arc-sec\n\n",grid_resolution*3600);
 
  for(int k=0; k<input.size();k++) {
    printf("treating %s\n", input[k].c_str());
    if(DecimateOnly) {
      int verbose=(input.size()==1);
      status=filter_decimate(input[k].c_str(), HeaderString, proj4, scale, d_offset, resolution, factor, prescribed, 
                             polygons, trusted, filter, paranoiac, duplicate, cartesian, autoscale, verbose, debug);
      }
    else {
      float offset[3];
      offset[0]=d_offset;  // random data offset
      offset[1]=b_offset;  // gridded bathymetry offset
      offset[2]=e_offset;  // gridded elevation offset
      mesh_t *mesh;
      mesh=xyz_loadmap_UG(input[k].c_str(), rootname, bathymetry, elevation, PBMA, distance_to_coast, proj4, 
                          scale, dilatation, offset, resolution, factor, grid_resolution,
                          DecimateNone, BASEincr, SRTMincr, ElementMaxSize, ElementMaxRatio, 
                          prescribed, e_poly, polygons, mapping, duplicate, meshout, cartesian, isobath, debug);
      }
    }
 
// /*------------------------------------------------------------------------------
//   save rectified topo */
//   output=new char[1024];
//   sprintf(output,"%s.spherical.nc",rootname);
// 
//   printf("#################################################################\n");
//   printf("save bathymetry file : %s\n",output);

  STDOUT_BASE_LINE("bathymetry sucessfully completed\n");

  exit(0);

}
