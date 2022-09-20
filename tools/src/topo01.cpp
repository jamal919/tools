
/*******************************************************************************

  T-UGO tools, 2006-2009

  Unstructured Ocean Grid initiative

Contributors:

  Florent Lyard      LEGOS/CNRS, Toulouse, France
  Yoann Le Bars      LEGOS, Toulouse, France (PhD)
  Laurent Roblou     LEGOS/CNRS, Toulouse, France
  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
  Frédéric Dupont    Université de Laval à Québec, Canada

E-mail: florent.lyard@legos.obs-mip.fr

*******************************************************************************/

// #define MAIN_SOURCE
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

#include "poc-netcdf-data.hpp"

#include "bmg.h"
#include "topo.h"

#define XYZ 0
#define YXZ 1
#include "map.def"

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_getgrid2d_0 (const char* filename,cdfgbl_t global, cdfvar_t info, grid_t *grid, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  {
  int status, ncid;
  size_t lengthhp,*dimlgth=NULL;
  int var,vx,vy;
  nc_type xtypep;
  char *text;

  double factor;
  decoded_t decoded;
  cdfvar_t vx_info,vy_info;
  int verbose=0;
  
  if(debug) verbose=1;

/*------------------------------------------------------------------------------
  check variable "associate" attribute */
  var=info.id;
  if(info.ndim <2) goto error;

  grid->circular=0; /*default*/
  status= poc_decode_associates(info, global, &decoded, 1, debug);

  vx=decoded.vx;
  vy=decoded.vy;

//  printf("---->variable x: %d, variable y: %d variable z: %d\n",vx,vy,vz,vt);

  if( (vx==-1) || (vy==-1) ) goto error;

  status=cdf_varinfo(filename,global.variable[vx].name,&vx_info,verbose);
  status=cdf_varinfo(filename,global.variable[vy].name,&vy_info,verbose);

/*------------------------------------------------------------------------------
  open file for reading */
  status=nc_open(filename,0,&ncid);
  if(status != NC_NOERR) goto open_error;

/*------------------------------------------------------------------------------
  assume x dim is last */
  grid->nx=info.dim[info.ndim-1].length;

/*------------------------------------------------------------------------------
  assume y dim is last-1 */
  grid->ny=info.dim[info.ndim-2].length;

  grid->nz=1;

/*------------------------------------------------------------------------------
  get x coordinate data */
  var=vx;
  switch (vx_info.ndim) {
    case 2:
      grid->x=new double [grid->nx*grid->ny];
      status=nc_get_var_double(ncid,var,grid->x);
      if(status != NC_NOERR) goto error;
      grid->modeH=2;
      break;

    case 1:
      grid->x=new double [grid->nx];
      status=nc_get_var_double(ncid,var,grid->x);
      if(status != NC_NOERR) goto error;
      grid->xmin=grid->x[0];
      grid->xmax=grid->x[grid->nx-1];
      grid->modeH=1;
      break;

    default:
      printf("dimension of x=%d ???\n",vx_info.ndim);
      goto error;
      break;
    }

  factor=1.;

  status=nc_inq_att(ncid,var,"units",&xtypep,&lengthhp);
  if(status == NC_NOERR) {
    text=(char *) malloc(lengthhp+1);
    status=nc_get_att_text(ncid,var,"units",text);
    text[lengthhp]='\0';
    if(strcmp(text,"radian_east") ==0) factor=180./M_PI;;
    if(strcmp(text,"degree_east")==0)  factor=1.;
    if(strcmp(text,"degree_east") ==0) factor=1.;
    free(text);
    }

  status=nc_inq_att(ncid,var,"topology",&xtypep,&lengthhp);
  if(status == NC_NOERR) {
    text=(char *) malloc(lengthhp+1);
    status=nc_get_att_text(ncid,var,"topology",text);
    text[lengthhp]='\0';
    if(strcmp(text,"circular")==0) grid->circular=1;
    free(text);
    }

/*------------------------------------------------------------------------------
  get y coordinate data */
  var=vy;
  switch (vy_info.ndim) {
    case 2:
      grid->y=new double [grid->nx*grid->ny];
      status=nc_get_var_double(ncid,var,grid->y);
      grid->ymin=grid->y[0];
      grid->ymax=grid->y[grid->nx*grid->ny-1];
      if(status != NC_NOERR) goto error;
      grid->modeH=2;
      break;

    case 1:
      grid->y=new double [grid->ny];
      status=nc_get_var_double(ncid,var,grid->y);
      if(grid->y[0]<grid->y[grid->ny-1]) {
        grid->ymin=grid->y[0];
        grid->ymax=grid->y[grid->ny-1];
        }
      else {
        grid->ymin=grid->y[grid->ny-1];
        grid->ymax=grid->y[0];
        }
      if(status != NC_NOERR) goto error;
      grid->modeH=1;
      break;

    default:
      printf("dimension of y=%d ???\n",vy_info.ndim);
      goto error;
      break;
    }

  status=map_minmax(grid);
  
  factor=1.;
   
  status=nc_inq_att(ncid,var,"units",&xtypep,&lengthhp);
  if(status == NC_NOERR) {
    text=(char *) malloc(lengthhp+1);
    status=nc_get_att_text(ncid,var,"units",text);
    text[lengthhp]='\0';
    if(strcmp(text,"radian_north")==0) factor=180./M_PI;
    if(strcmp(text,"degree_north")==0) factor=1.;
    if(strcmp(text,"degree_north")==0) factor=1.;
    free(text);
    }

  grid->zmin=0.;
  grid->zmax=0.;
  grid->z=NULL;
  grid->nz=1;
  grid->dz=0.0;

  grid->modeV=-1;


/*------------------------------------------------------------------------------
  default */
  grid->circular=0;
  grid->overlapped=0;
  grid->connex=1;

  grid->connex=mapc_checkconnexity(*grid);
  if(grid->x[grid->nx-1]-grid->x[0]>360.0) {
     grid->circular=1;
     grid->overlapped=1;
     }

  if(grid->xmax-grid->xmin>359.0) {
     grid->circular=1;
     grid->overlapped=1;
     }

  status=nc_close(ncid);
  nc_check_error(status,__LINE__,__FILE__);
  if(status != NC_NOERR) goto close_error;

  status= poc_gettime(filename,  info,  global, &(grid->origine), &(grid->time), &(grid->nt));

  status=0;
  return(status);

  error:
  nc_check_error(status,__LINE__,__FILE__);
  if(dimlgth != NULL) free(dimlgth);
  status=nc_close(ncid);
  nc_check_error(status,__LINE__,__FILE__);
  status=-1;
  printf("loading 2D grid for variable %d failed\n",info.id);
  return(status);

  close_error:
  nc_check_error(status,__LINE__,__FILE__);
  status=3; /* not a file?*/
  return(status);

  open_error:
  nc_check_error(status,__LINE__,__FILE__);
  status=3; /* not a file?*/
  return(status);
  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> int topo_loadfield_grd_template(const char *input,grid_t *grid, T **buffer, T *mask, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  size_t m,count;

  status=grd_loadgrid(input,grid);
  if(status!=0) NC_TRAP_ERROR(wexit,status,1,"cannot load bathymetry file %s",input);
  
  *buffer= new T[grid->Hsize()];
  status=grd_loadr1(input,*grid,*buffer,mask);
  
//   count=occurence(*mask, *buffer, grid->Hsize());

  if(isnan(*mask)) {
    *mask=1.e+10;
/*------------------------------------------------------------------------------
    nan is very inconvenient as mask value*/
    for(m=0;m<grid->Hsize();m++) {
      if(isnan((*buffer)[m])) {
        (*buffer)[m]=*mask;
        }
      }
    }

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int topo_loadfield_grd(const char *input,grid_t *grid, float **buffer, float *mask, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;

  status=topo_loadfield_grd_template( input, grid, buffer, mask, debug);

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int topo_loadfield_grd(const char *input,grid_t *grid, short **buffer, short *mask, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;

  status=topo_loadfield_grd_template( input, grid, buffer, mask, debug);

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int topo_loadfield_cdf(const char *input,grid_t *grid, float **buffer, float *mask, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int    v,status;
  int    verbose,frame=0;
  cdfgbl_t data_info,grid_info;
  cdfvar_t h_info;

  verbose=0;
  status=cdf_globalinfo(input,&data_info,verbose);
  for (v=0;v<data_info.nvarsp;v++) {
    printf("variable %3d: name %s, type %d,ndim %d \n",v,(data_info.variable[v]).name,data_info.variable[v].type,data_info.variable[v].ndim);
    }
  status=cdf_varinfo(input,"bathymetry",&h_info,1);

  status=cdf_globalinfo(input,&grid_info,verbose);
  status= poc_getgrid2d ((const char*) input, grid_info, h_info, grid);

  *buffer=new float[grid->nx*grid->ny];
  status= poc_getvar2d (input, h_info.id, frame, *buffer, mask, h_info);

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int topo_loadfield_cdf(const char *input,const char *varname, grid_t *grid, float **buffer, float *mask, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int    status;
  int    verbose=0,frame=0;
  
#if 0
  cdfgbl_t data_info,grid_info;
  cdfvar_t h_info;

  verbose=0;
  
  if(debug) verbose=1;
  
  status=cdf_globalinfo(input,&data_info,verbose);
  if(status!=0) TRAP_ERR_EXIT(-1,"exiting\n");
  
  if(debug) {
    int v;
    for (v=0;v<data_info.nvarsp;v++) {
      printf("variable %3d: name %15s, type %d,ndim %d \n",v,(data_info.variable[v]).name,data_info.variable[v].type,data_info.variable[v].ndim);
      }
    }
  
  if(debug) verbose=1;
  
  status=cdf_varinfo(input,varname,&h_info,verbose);
  if(status!=0) return(status);

  status=cdf_globalinfo(input,&grid_info,verbose);
  status= poc_getgrid2d_0 (input, grid_info, h_info, grid, debug);
  if(status!=0) return(status);

  *buffer=new float[grid->nx*grid->ny];
  status= poc_getvar2d (input, h_info.id, frame, *buffer, mask, h_info);
#else
  poc_data_t<float> field;
  status=field.init(input,varname,verbose);
  status=field.read_data(input,0,verbose,1);
  status=poc_get_grid(input,field.info,grid,verbose,frame);
  swapValues(buffer,&field.data);
  *mask=field.mask;
#endif

  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int topo_loadfield_cdf(const char *input,const char *varname, grid_t *grid, signed char **buffer, signed char *mask, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int    /*v,*/status;
  int    verbose,frame=0;
  cdfgbl_t data_info,grid_info;
  cdfvar_t h_info;

  verbose=0;
  status=cdf_globalinfo(input,&data_info,verbose);
//   for (v=0;v<data_info.nvarsp;v++) {
//     printf("variable %3d: name %s, type %d,ndim %d \n",v,(data_info.variable[v]).name,data_info.variable[v].type,data_info.variable[v].ndim);
//     }
  status=cdf_varinfo(input,varname,&h_info,1);
  if(status!=0) return(status);

  status=cdf_globalinfo(input,&grid_info,verbose);
  status= poc_getgrid2d_0 (input, grid_info, h_info, grid, debug);
  if(status!=0) return(status);

  *buffer=new signed char[grid->nx*grid->ny];
  status= poc_getvar2d (input, h_info.id, frame, *buffer, mask, h_info);

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int topo_loadfield_cdf(const char *input,const char *varname, grid_t *grid, short **buffer, short *mask, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int    /*v,*/status;
  int    verbose,frame=0;
  cdfgbl_t data_info,grid_info;
  cdfvar_t h_info;

  verbose=0;
  status=cdf_globalinfo(input,&data_info,verbose);
//   for (v=0;v<data_info.nvarsp;v++) {
//     printf("variable %3d: name %s, type %d,ndim %d \n",v,(data_info.variable[v]).name,data_info.variable[v].type,data_info.variable[v].ndim);
//     }
  status=cdf_varinfo(input,varname,&h_info,1);
  if(status!=0) return(status);

  status=cdf_globalinfo(input,&grid_info,verbose);
  status= poc_getgrid2d_0 (input, grid_info, h_info, grid, debug);
  if(status!=0) return(status);

  *buffer=new short[grid->nx*grid->ny];
  status= poc_getvar2d (input, h_info.id, frame, *buffer, mask, h_info);
  if(status!=0) return(status);

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> int topo_loadfield_template(const char *filename, const char *varname, grid_t *grid, T **buffer, T *mask, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  vector<string> variable;

  *buffer=0;

  if(strrncmp(filename,".grd")==0) {
    status=topo_loadfield_grd(filename,grid, buffer, mask, debug);
    }
  else if(strrncmp(filename,".bmg")==0) {
    float time;
    status=bmg_loadgrid (filename,grid);
//    status=map_completegridaxis(grid,2);
    *buffer=new T[grid->nx*grid->ny];
    status= bmg_loadr1(filename,1,1,1,*grid,*buffer,mask,&time);
    }
  else if(strrncmp(filename,".nc")==0) {
    if(varname!=0 and varname[0]!='\0') variable.push_back(varname);
    variable.push_back("bathymetry");
    variable.push_back("z");
    variable.push_back("depth");
    variable.push_back("dst");
    variable.push_back("Band1");
    variable.push_back("LAT");
    for(int k=0; k<variable.size();k++) {
      status=topo_loadfield_cdf(filename,variable[k].c_str(), grid, buffer, mask, debug);
      if(status==0) break;
      }
    }
  else
    status=-1;
  
  if(status==0) grid->circular=map_check_circular(*grid);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int topo_loadfield(const char *filename, grid_t *grid, float **buffer, float *mask, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=topo_loadfield_template(filename, 0, grid, buffer, mask, debug);
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int topo_loadfield(const char *filename, const char *varname, grid_t *grid, float **buffer, float *mask, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=topo_loadfield_template(filename, varname, grid, buffer, mask, debug);
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int topo_loadfield(const char *filename, grid_t *grid, short **buffer, short *mask, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=topo_loadfield_template(filename, 0, grid, buffer, mask, debug);
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int topo_loadfield(const char *filename, const char *varname, grid_t *grid, short **buffer, short *mask, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=topo_loadfield_template(filename, varname, grid, buffer, mask, debug);
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> int topo_savefield_template(const char *filename, grid_t grid, T *buffer, T mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int    status;
  string s;
  pocgrd_t ncgrid;
  cdfvar_t variable;

  if( strrncasecmp(filename,".grd")==0 ) {
//     status=topo_loadfield_grd(filename,grid, buffer, mask);
    }
  if( strrncasecmp(filename,".bmg")==0 ) {
//     float time;
//     status=bmg_loadgrid (filename,grid);
// //    status=map_completegridaxis(grid,2);
//     *buffer=new T[grid->nx*grid->ny];
//     status= bmg_loadr1(filename,1,1,1,*grid,*buffer,mask,&time);
    }
  else if( strrncasecmp(filename,".nc")==0 ) {
    status= poc_createfile(filename);
    status= map_completegridaxis(&grid,2);
    status=poc_sphericalgrid_xy(filename,"",grid,&ncgrid);
    poc_standardvariable_xy(&variable,"bathymetry", mask,"m",1., 0.,"bathymetry","bathymetry","bathymetry",ncgrid);
    status=create_ncvariable(filename, &variable);
    status=poc_write_xy(filename, grid, variable.id, buffer);
    variable.destroy();

    }
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int topo_savefield(const char *filename, grid_t grid, float *buffer, float mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;

  status=topo_savefield_template(filename, grid, buffer, mask);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

frame_t plg_recale_local(plg_t *polygones, int npolygones)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,j;
  frame_t frame;
  double t,t0;

  for(i=0;i<npolygones;i++) {
    if(polygones[i].npt==0) continue;
    t0=polygones[i].t[0];
    for(j=0;j<polygones[i].npt;j++) {
      t=polygones[i].x[j];
      polygones[i].x[j]=geo_recale(t,t0,180.0);
      }
    }
  frame=plg_spherical_minmax(polygones, npolygones);
  return(frame);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

frame_t plg_recale_local(vector<plg_t> & polygons)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,j;
  frame_t frame;
  double t,t0;

  for(i=0;i<polygons.size();i++) {
    if(polygons[i].npt==0) continue;
    t0=polygons[i].t[0];
    for(j=0;j<polygons[i].npt;j++) {
      t=polygons[i].x[j];
      polygons[i].x[j]=geo_recale(t,t0,180.0);
      }
    }
  frame=plg_spherical_minmax(polygons);
  return(frame);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int plg_test_local(double t,double p,plg_t *polygones, frame_t *frame, int npolygones)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  plg_t plg_polygones;
  int tmp,in=0;
  int k;

  for(k=0;k<npolygones;k++) {
    if (t<frame[k].xmin) t+=360.0;
    if (t>frame[k].xmax) t-=360.0;
    if ((t<frame[k].xmin) || (t>frame[k].xmax)) {
      continue;
      }
    tmp=plg_single( polygones[k], t, p,&in);
    if(tmp==PLG_POINT_BOUNDARY) break;
    }
  return(tmp);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  size_t topo_PolygonsSelection(grid_t grid, plg_t *polygones, int npolygones, signed char *selected, bool spherical, int initialise)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int   k,i,j,n;
  size_t   count=-1;
  char *position;
  frame_t frame, *frames;
  double *x,*y;
  bool debug=false;
  int verbose=0;
  size_t *indice;
  
  if(initialise==1) {
    for(size_t m=0;m<grid.Hsize();m++) selected[k]=1;
    }
  
  printf("#################################################################\n");
  printf("select grid nodes in polygons\n");
  
  count=occurence((signed char) 1,selected,grid.Hsize());
  printf ("initial number of nodes selected %d (over %d)\n",count,grid.Hsize());
  
/*-----------------------------------------------------------------------------
  area criterion */
  if(npolygones==0) {
    for(size_t m=0;m<grid.Hsize();m++) selected[k]=0;
    goto end;
    }

  if(spherical)
    frame=plg_recale_local(polygones,npolygones);
  else
    frame=plg_cartesian_minmax(polygones,npolygones);
  
  frames= new frame_t[npolygones];
  for(k=0;k<npolygones;k++) {
    frames[k]=plg_cartesian_minmax(polygones[k]);
    }
 
  x=new double[count];
  y=new double[count];
  indice=new size_t[count];

  n=0;
  for (j=0;j<grid.ny;j++) {
    for (i=0;i<grid.nx;i++) {
      size_t m=grid.nx*j+i;
      if(selected[m]==1) {
        x[n]=map_grid_x(grid,i,j);
        y[n]=map_grid_y(grid,i,j);
        indice[n]=m;
        n++;
        }
      }
    }
  
 position=plg_TestInterior(x, y, count, polygones, npolygones, verbose, debug);
 
 for(n=0;n<count;n++) if(position[n]!=PLG_POINT_INTERIOR) selected[indice[n]]=0;
 
 delete[] x;
 delete[] y;
 delete[] position;
 
//  nprocs=initialize_OPENMP(-1);
//  
// #pragma omp parallel for private(i,j,n,inside,t,p) if(nprocs>1)
//   for (j=0;j<grid.ny;j++) {
//     for (i=0;i<grid.nx;i++) {
//       n=grid.nx*j+i;
//       if(selected[n]==0) continue;
//       t=map_grid_x(grid,i,j);
//       p=map_grid_y(grid,i,j);
// //       if ((t<frame.xmin) || (t>frame.xmax)) {
// //         selected[n]=0;
// //         continue;
// //         }
//       if ((p<frame.ymin) || (p>frame.ymax)) {
//         selected[n]=0;
//         continue;
//         }
// //      inside=plg_TestInterior(t,p,polygones, frames, npolygones);
//       inside=plg_test_local(t,p,polygones, frames, npolygones);
//       if (inside==PLG_POINT_INTERIOR) {
//         selected[n]=1;
//         }
//       else {
//         selected[n]=0;
//         }
//       }
//     }

end:

  count=occurence((signed char) 1,selected,grid.Hsize());
  printf ("resulting number of nodes selected %d (over %d)\n",count,grid.Hsize());
  
  return(count);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  size_t topo_PolygonsSelection(grid_t grid, vector<plg_t> & polygons, signed char *selected, bool spherical, int initialise)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int   k,i,j,n;
  size_t   count=-1;
  char *position;
  frame_t frame, *frames;
  double *x,*y;
  bool debug=false;
  int verbose=0;
  size_t *indice;
  
  if(initialise==1) {
    for(size_t m=0;m<grid.Hsize();m++) selected[k]=1;
    }
  
  printf("#################################################################\n");
  printf("select grid nodes in polygons\n");
  
  count=occurence((signed char) 1,selected,grid.Hsize());
  printf ("initial number of nodes selected %d (over %d)\n",count,grid.Hsize());
  
/*-----------------------------------------------------------------------------
  area criterion */
  if(polygons.size()==0) {
    for(size_t m=0;m<grid.Hsize();m++) selected[k]=0;
    goto end;
    }

  if(spherical)
    frame=plg_recale_local(polygons);
  else
    frame=plg_cartesian_minmax(polygons);
  
  frames= new frame_t[polygons.size()];
  for(k=0;k<polygons.size();k++) {
    frames[k]=plg_cartesian_minmax(polygons[k]);
    }
 
  x=new double[count];
  y=new double[count];
  indice=new size_t[count];

  n=0;
  for (j=0;j<grid.ny;j++) {
    for (i=0;i<grid.nx;i++) {
      size_t m=grid.nx*j+i;
      if(selected[m]==1) {
        x[n]=map_grid_x(grid,i,j);
        y[n]=map_grid_y(grid,i,j);
        indice[n]=m;
        n++;
        }
      }
    }
  
 position=plg_TestInterior(x, y, count, polygons, verbose, debug);
 
 for(n=0;n<count;n++) if(position[n]!=PLG_POINT_INTERIOR) selected[indice[n]]=0;
 
 delete[] x;
 delete[] y;
 delete[] position;
 
end:

  count=occurence((signed char) 1,selected,grid.Hsize());
  printf ("resulting number of nodes selected %d (over %d)\n",count,grid.Hsize());
  
  return(count);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  size_t topo_PolygonsSelection(grid_t grid, const char *PolygonsFile, signed char *selected, bool spherical, int initialise)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  size_t count;
  const char *rootname="anonymous";
  
  vector<plg_t> polygons,q;
  
  if(PolygonsFile==0) return(0);
  
  printf("#----------------------------------------------------------------\n");
  printf("select grid nodes in polygon: %s\n",PolygonsFile);
  
  status=plg_load(PolygonsFile, polygons);
  if(status!=0) return(-1);
  
//   q=plg_dilatation_cartesian(polygons,1.1,true);
//   status=plg_save("dilatation.plg",PLG_FORMAT_UNKNOWN, q);
  
  status=plg_CheckDuplicated(polygons, true, rootname, false);
  
  count=topo_PolygonsSelection(grid, polygons, selected, spherical, initialise);
  
  plg_destroy_entries(polygons);
  
  polygons.clear();
  
  return(count);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int topo_RangeSelection(const grid_t & grid, float *topo, float mask, range_t<float> range, signed char *selected, int initialise)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int   k,i,j,n;
  int   count=-1,nprocs;
  frame_t frame;

  printf("#################################################################\n");
  printf("select grid nodes %f < z < %f\n",range.min,range.max);
  
  count=occurence((signed char) 1,selected,grid.Hsize());
  printf ("initial number of nodes selected %d (over %d)\n",count,grid.Hsize());
  
  if(initialise==1) {
    for(j=0;j<grid.ny;j++) {
      for(i=0;i<grid.nx;i++) {
        k=j*grid.nx+i;
        selected[k]=1;
        }
      }
    }

 nprocs=initialize_OPENMP(-1, 0);
 
#pragma omp parallel for private(i,j,n) if(nprocs>1)
  for (j=0;j<grid.ny;j++) {
    for (i=0;i<grid.nx;i++) {
      n=grid.nx*j+i;
      if(selected[n]==0) continue;
      if ((topo[n]>range.min) and (topo[n]<range.max)) {
        selected[n]=1;
        }
      else {
        selected[n]=0;
        }
      }
    }

  count=occurence((signed char) 1,selected,grid.Hsize());
  printf ("resulting number of nodes selected %d (over %d)\n",count,grid.Hsize());
  
  return(count);
}


// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   int topo_RangeSelection(grid_t grid, float *topo, float mask, range_t<float> range, signed char *selected)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   int result;
//   
//   result=topo_RangeSelection(grid, topo, mask, range, selected, 1);
//   
//   return result;
// }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int topo_operation(const char *bathymetry, int signus, grid_t topogrid, float *topo, float topomask, float zmin, float zmax, char *poly, int masked_only, int positive_only, int persistence)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int    j,status;
  float *buffer,mask;
  grid_t grid;
  bool debug=false;

  
  status=topo_loadfield(bathymetry, &grid, &buffer, &mask, debug);
  if(status!=0) return(status);
//   if(bathymetry!=0) {
//     status=topo_loadfield(bathymetry, &grid, &buffer, &mask);
//     if(status!=0) return(status);
//     }
//   else {
//     }
  
  for (int k=0;k<persistence;k++) status=map_persistence(grid, buffer, mask, 0.0);

/* *------------------------------------------------------------------------------
  Interpolate minimum low tide level and add to topo*/
#define OPEN_MP

#ifdef OPEN_MP
  int nprocs=initialize_OPENMP(-1);
  #pragma omp parallel for private(status) if(nprocs>1)
#endif
  for(j=0;j<topogrid.ny;j++) {
    for(int i=0;i<topogrid.nx;i++) {
      int k=j*topogrid.nx+i;
      float z;
      double x,y;
      if(topo[k]==topomask) {
        continue;
        }
      if(topo[k]>zmax) {
        topo[k]=topomask;
        continue;
        }
      if(topo[k]<zmin) {
        topo[k]=topomask;
        continue;
        }
      topogrid.xy(i,j,x,y);
      x=map_recale(grid,x);
      status=map_interpolation(grid, buffer,mask,x,y,&z);
      if(positive_only==1) {
        if(z>0 and z!=mask) {
          topo[k]=topomask;
          continue;
          }
        }
      if (isnan(z)) {
        topo[k]=topomask;
        continue;
        }
      if (z!=mask) {
/* *------------------------------------------------------------------------------
        assume negative depths and negative lowtide*/
        topo[k]+=signus*z;
        }
      }
    }

  delete[] buffer;

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int topo_set(const grid_t & topogrid, float *topo, float topomask, float zmin, float zmax,
               const char *poly, int masked_only, int valid_only, int exclusion, float value, short *tag, short tagvalue, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int   i,j,k,status;
  int initialise;
  signed char *selected;
  size_t count=0;

//   if(debug) printf("import file: %s, mask=%f\n",bathymetry, mask);

  selected=new signed char[topogrid.Hsize()];
  
  for(size_t m=0;m<topogrid.Hsize();m++) selected[k]=0;
  for(j=0;j<topogrid.ny;j++) {
    for(i=0;i<topogrid.nx;i++) {
      k=j*topogrid.nx+i;
      if(topo[k]!=value) selected[k]=1;
      }
    }

  initialise=0;
  if (zmin!=zmax) {
    range_t<float> range(zmin, zmax);
    status=topo_RangeSelection(topogrid, topo, topomask, range, selected, initialise);
    initialise=0;
    }
  
  if(poly!=0) {
    status=topo_PolygonsSelection(topogrid, poly, selected, initialise);
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
  
  if(exclusion==1) {
    count=occurence((signed char) 1,selected,topogrid.Hsize());
    printf("#exclusion aplies: initial  count=%d\n",count);
    for(size_t m=0;m<topogrid.Hsize();m++) selected[m]=1-selected[m];
    count=occurence((signed char) 1,selected,topogrid.Hsize());
    printf("#exclusion aplies: reversed count=%d\n",count);
    }
  
/*-----------------------------------------------------------------------------
  Interpolate and store new topo */
  int nprocs=initialize_OPENMP(-1);
  
  count=0;
#pragma omp parallel for if(nprocs>1)
  for(int j=0;j<topogrid.ny;j++) {
    for(int i=0;i<topogrid.nx;i++) {
      int k=j*topogrid.nx+i;
      if(selected[k]==0) continue;
      if (topo[k]!=topomask) {
        if(tag!=0) tag[k]=value;
        topo[k]=value;
#pragma omp critical(threadUnsafe)
        {count++;}
        }
      }
    }
  
  printf("#modified soundings: %d\n",count);

  delete[] selected;
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int topo_shift(const grid_t & topogrid, float *topo, float topomask, float zmin, float zmax,
               const char *poly, int masked_only, float value, short *tag, short tagvalue, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int   i,j,k,status;
  int initialise;
  signed char *selected;
  size_t count=0;

//   if(debug) printf("import file: %s, mask=%f\n",bathymetry, mask);

  selected=new signed char[topogrid.Hsize()];
  
  for(size_t m=0;m<topogrid.Hsize();m++) selected[k]=0;
  for(j=0;j<topogrid.ny;j++) {
    for(i=0;i<topogrid.nx;i++) {
      k=j*topogrid.nx+i;
      if(topo[k]!=value) selected[k]=1;
      }
    }

  initialise=0;
  if (zmin!=zmax) {
    range_t<float> range(zmin,zmax);
    status=topo_RangeSelection(topogrid, topo, topomask, range, selected, initialise);
    initialise=0;
    }

  if(poly!=0) {
    status=topo_PolygonsSelection(topogrid, poly, selected, initialise);
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
  
/*-----------------------------------------------------------------------------
  Interpolate and store new topo */
  int nprocs=initialize_OPENMP(-1);
  
#pragma omp parallel for if(nprocs>1)
  for(int j=0;j<topogrid.ny;j++) {
    for(int i=0;i<topogrid.nx;i++) {
      int k=j*topogrid.nx+i;
      if(selected[k]==0) continue;
      if (topo[k]!=topomask) {
        if(tag!=0) tag[k]=value;
        topo[k]+=value;
#pragma omp critical(threadUnsafe)
        {count++;}
        }
      }
    }
  
  printf("#modified soundings: %d\n",count);

  delete[] selected;
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int topo_import(const char *bathymetry, const char *varname, const grid_t & topogrid, float *topo, float topomask, float zmin, float zmax,
                  char *poly, int masked_only, short *tag, short value, const char *input_format, const char *proj4, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int initialise;
  float *buffer,mask;
  grid_t grid;
  signed char *selected;
  size_t count=0;
  frame_t fbase(topogrid),fimport;

  if(input_format==0) {
    status=topo_loadfield(bathymetry, varname, &grid, &buffer, &mask, debug);
    }
  else {
    int fmt=map_get_format(input_format);
    char *xname=0, *yname=0;;
    status=map_loadfield(bathymetry, fmt, varname, xname, yname, proj4, &grid, &buffer, &mask);
    }
  if(status !=0) {
    printf("cannot load bathymetry file=%s\n",bathymetry);
    return(-1);
    }

  map_printgrid(grid);
  
//   fimport=frame_t(grid);
//   fbase.xmax=fbase.x_size()+map_recale(grid,fbase.xmin);
//   fbase.xmin=map_recale(grid,fbase.xmin);
  
  fimport=frame_t(grid);
  double x=map_recale(grid,fbase.x_center());
  
  double width=fbase.x_size();
  
  fbase.xmin=x-width/2.0;
  fbase.xmax=fbase.xmin+width;
  
  double size=fbase.intersection(fimport).size();
  
  if(size==0.0) {
//     if(debug) printf("import file: %s has no intersection, skip\n",bathymetry);
    fimport.print("import file frame: ");
    fbase.print("base file frame: ");
    printf("import file: %s has no intersection, skip\n",bathymetry);
    delete[] buffer;
    grid.free();
    return(0);
    }
  else {
    if(debug) printf("import file: %s, mask=%f\n",bathymetry, mask);
    }

/*-----------------------------------------------------------------------------
  Interpolate and store new topo */
  int nprocs=initialize_OPENMP(-1);
  
  selected=new signed char[topogrid.nx*topogrid.ny];
#pragma omp parallel for if(nprocs>1)
  for(size_t j=0;j<topogrid.ny;j++) {
    for(size_t i=0;i<topogrid.nx;i++) {
      size_t k=j*topogrid.nx+i;
      selected[k]=1;
      }
    }

  initialise=0;
  if(poly!=0) {
    status=topo_PolygonsSelection(topogrid, poly, selected, initialise);
    initialise=0;
    }
  if(status <0) {
    printf("selection failed, abort...\n");
    return(-1);
    }
  
  if (zmin!=zmax) {
    range_t<float> range(zmin, zmax);
    status=topo_RangeSelection(topogrid, topo, topomask, range, selected, initialise);
    initialise=0;
    }

  if(masked_only==1) {
#pragma omp parallel for if(nprocs>1)
    for(size_t j=0;j<topogrid.ny;j++) {
      for(size_t i=0;i<topogrid.nx;i++) {
        size_t k=j*topogrid.nx+i;
        if(topo[k]!=topomask) selected[k]=0;
        }
      }
    }
  
  int verbose=0;
  set_grid_list(&grid,1,verbose);
  
#pragma omp parallel for if(nprocs>1)
  for(size_t j=0;j<topogrid.ny;j++) {
    for(size_t i=0;i<topogrid.nx;i++) {
      size_t k=j*topogrid.nx+i;
      int status=0;
      double x,y;
      float z;
      int64_t prior=-1;
      if(selected[k]==0) continue;
//      x=topogrid.x[k];
//      y=topogrid.y[k];
      x=map_grid_x(topogrid,i,j);
      x=map_recale(grid,x);
      y=map_grid_y(topogrid,i,j);
//       status=map_interpolation(grid, buffer,mask,x,y,&z);
      index_interpolation(grid,x,y,&prior,buffer,mask,&z,verbose);
//      if (z!=mask) {
      if ((status==0)&&(z!=mask)&&(isnan(z)!=1)) {
        if(tag!=0) tag[k]=value;
        topo[k]=z;
#pragma omp critical(threadUnsafe)
        {count++;}
        }
      }
    }
  
  printf("#modified soundings: %d (over %d)\n",count, topogrid.Hsize());

  delete[] buffer;
  delete[] selected;
  grid.free();
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int topo_import(char *bathymetry, const char *varname, const grid_t & topogrid, float *topo, float topomask, float zmin, float zmax, char *poly, int masked_only, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=topo_import(bathymetry, varname, topogrid, topo, topomask, zmin, zmax, poly, masked_only, (short *) 0, (short) 0, (const char *) 0, (const char *) 0, debug);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int topo_import(char *bathymetry, const char *varname, const grid_t & topogrid, float *topo, float topomask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int   i,j,k,status;
  grid_t grid;
  signed char *buffer,mask;
  double x,y;
  float z;
  bool debug=false;

  printf("warning, temporary byte reading...");
  status=topo_loadfield_cdf(bathymetry, varname, &grid, &buffer, &mask, debug);
  if(status!=0) return(status);

/* *------------------------------------------------------------------------------
  Interpolate and store topo on symphonie grid*/
//  mask=-31072.0;
  for(j=0;j<topogrid.ny;j++) {
    for(i=0;i<topogrid.nx;i++) {
      k=j*topogrid.nx+i;
      x=map_grid_x(topogrid,i,j);
      y=map_grid_y(topogrid,i,j);
      status=map_interpolation(grid, buffer,mask,x,y,&z);
//      if (z!=mask) {
        topo[k]=z;
//        }
      }
    }

  delete[] buffer;

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int topo_hydro(char *hydro, int signus, grid_t topogrid, float *topo, float topomask, int persistence)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int   i,j,status;
  float *buffer,mask,z;
  grid_t grid;
  bool debug=false;

  status=topo_loadfield(hydro, &grid, &buffer, &mask, debug);
  if(status!=0) return(status);
  
  for (int k=0;k<persistence;k++) status=map_persistence(grid, buffer, mask, 0.0);

/* *------------------------------------------------------------------------------
  Interpolate minimum low tide level and add to topo*/
  for(j=0;j<topogrid.ny;j++) {
    for(i=0;i<topogrid.nx;i++) {
      int k=j*topogrid.nx+i;
      double x,y;
      if(topo[k]==topomask) {
        continue;
        }
//       if(topo[k]>0) {
//         topo[k]=topomask;
//         continue;
//         }
//       if(topo[k]<-500) {
//         topo[k]=topomask;
//         continue;
//         }
      x=map_grid_x(topogrid,i,j);
      y=map_grid_y(topogrid,i,j);
      status=map_interpolation(grid, buffer,mask,x,y,&z);
//       if(z>0) {
//         topo[k]=topomask;
//         continue;
//         }
      if (z!=mask) {
/* *------------------------------------------------------------------------------
        assume negative depths and negative lowtide*/
        topo[k]+=signus*z;
        }
      }
    }

  delete[] buffer;

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 int topo_checks(const char *rootname, grid_t topogrid, float * topo, float topomask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  char output[1024];
  pocgrd_t ncgrid;
  cdfvar_t variable;
 
  size_t size=topogrid.nx*topogrid.ny;
  float *topo_x=new float[size];
  float *topo_y=new float[size];
  status= map_gradient(topogrid, topogrid.nx, topo, topomask, GEOCENTRIC, topo_x, topo_y);

 
  sprintf(output,"%s.spherical.nc",rootname);
  status= poc_createfile(output);

  if(topogrid.modeH==0)  status= map_completegridaxis(&topogrid,1);
  
  status=poc_sphericalgrid_xy(output,"",topogrid,&ncgrid);
  
  poc_standardvariable_xy(&variable,"bathymetry",topomask,"m",1., 0.,"bathymetry","bathymetry","bathymetry",ncgrid);
  status=create_ncvariable(output, &variable);
  if(status!=0) goto error;
  status=poc_write_xy(output,  topogrid, variable.id,topo);
  variable.destroy();
  
  poc_standardvariable_xy(&variable,"dhdx",topomask,"dimensionless",1., 0.,"dhdx","dhdx","dhdx",ncgrid);
  status=create_ncvariable(output, &variable);
  if(status!=0) goto error;
  status=poc_write_xy(output,  topogrid, variable.id,topo_x);
  variable.destroy();
    
  poc_standardvariable_xy(&variable,"dhdy",topomask,"dimensionless",1., 0.,"dhdy","dhdy","dhdy",ncgrid);
  status=create_ncvariable(output, &variable);
  if(status!=0) goto error;
  status=poc_write_xy(output,  topogrid, variable.id,topo_y);
  variable.destroy();
  
  delete[] topo_x;
  delete[] topo_y;
    
  return(0);
error:
  delete[] topo_x;
  delete[] topo_y;
    
  return(1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int topo_save(const string & output, const char *format, grid_t topogrid, float *ftopo, float ftopomask, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  handle mirroring for grd format

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int i,j,m,n,status;
  short *topo, smask;
  pocgrd_t ncgrid;
  cdfvar_t variable;
  
  if(format==NULL) format=strdup("grd");

//   if(debug) status=topo_checks(rootname.c_str(), topogrid, ftopo, ftopomask);

  if(strcmp(format,"grd")==0) {
    topo=new short[topogrid.nx*topogrid.ny];
    smask=256*127+255;
    for (j=0;j<topogrid.ny;j++) {
      for (i=0;i<topogrid.nx;i++) {
        m=j*topogrid.nx+i;
        n=(topogrid.ny-j-1)*topogrid.nx+i;
        if(ftopo[m]!=ftopomask) {
          topo[n]=(short) floor(ftopo[m]+0.5);
          }
        else {
          topo[n]=smask;
          }
        }
      }
    printf("#################################################################\n");
    printf("create depth file (grd short format) : %s\n",output.c_str());
    status=grd_save(output.c_str(), topogrid, topogrid.nx, topo,smask);
    delete[] topo;
    }
  else if(strcmp(format,"grd-float")==0) {
    printf("#################################################################\n");
    printf("create depth file (grd float format) : %s\n",output.c_str());
    status=grd_mirror_r( topogrid, topogrid.nx, ftopo, ftopomask);
    status=grd_save(output.c_str(),topogrid,topogrid.nx, ftopo, ftopomask);
/*------------------------------------------------------------------------------
    warning 2018-03-22 : now mirror back */
    status=grd_mirror_r( topogrid, topogrid.nx, ftopo, ftopomask);
    }
  else if(strcmp(format,"netcdf")==0) {
    printf("#################################################################\n");
    printf("create depth file (netcdf float format) : %s\n",output.c_str());
    status= poc_createfile(output.c_str());
    if(topogrid.modeH==0) status= map_completegridaxis(&topogrid,1);
    status=poc_sphericalgrid_xy(output.c_str(),"",topogrid,&ncgrid);
    poc_standardvariable_xy(&variable,"bathymetry",ftopomask,"m",1., 0.,"bathymetry","bathymetry","bathymetry",ncgrid);
    status=create_ncvariable(output.c_str(), &variable);
    status=poc_write_xy(output.c_str(),  topogrid, variable.id,ftopo);
    variable.destroy();
    }
      
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int topo_save(string rootname, string & filename, const char *format, grid_t topogrid, float *ftopo, float ftopomask, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  string output;
  
  if(rootname=="") rootname="out";

  if(format==NULL) format=strdup("grd");

  if(debug) status=topo_checks(rootname.c_str(), topogrid, ftopo, ftopomask);

  if(strcmp(format,"grd")==0) {
    output=rootname+".grd";
    }
  else if(strcmp(format,"grd-float")==0) {
    output=rootname+".grd";
    }
  else if(strcmp(format,"netcdf")==0) {
    output=rootname+".nc";
    }
  
  status=topo_save(output, format, topogrid, ftopo, ftopomask, debug);
  
  filename=output;
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int topo_save(const char *rootname, char **output, const char *format, grid_t topogrid, float *ftopo, float ftopomask, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status; 
  string filename;
  
  status=topo_save((string) rootname, filename, format, topogrid, ftopo, ftopomask, debug);
  
  if(output!=0) *output=strdup(filename.c_str());

  return(0);
}
