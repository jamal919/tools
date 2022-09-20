
/*******************************************************************************

  T-UGO tools, 2006-2014

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
\brief Set landmask and interpolate topography for a given grid.

<!-- A LINK TO main() or print_help() WILL NOT LINK TO THE RIGHT SOURCE ! -->
See the main <a href=#func-members>function</a> for how this works
and the print_help <a href=#func-members>function</a> for how to use this.
*/
/*----------------------------------------------------------------------------*/

#define MAIN_SOURCE

#include "config.h"

#include "version-macros.def" //for VERSION and REVISION

#include <stdio.h>
#include <unistd.h> // for getpid
#include <stdlib.h>
#include <string.h>

#include "tools-structures.h"

#include "functions.h"
#include "constants.h"
#include "map.h"
#include "geo.h"
#include "polygones.h"
#include "grd.h"
#include "sym-io.h"
#include "ascii.h"
#include "netcdf-proto.h"
#include "poc-netcdf-data.hpp"
#include "poc-time.h"
#include "datastream.h"

#include "polygons.hpp"
#include "exceptions.hpp"

using namespace Polygons;

extern int map_loadfield_grd(const char *input,grid_t *grid, float **buffer, float *mask);

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int map_loadfield_cdf(const char *input,grid_t *grid, float **buffer, float *mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int    v,status;
  int    verbose,frame=0;
  cdfgbl_t data_info,grid_info;
  cdfvar_t h_info;
  
  *buffer=0;

  verbose=0;
  
  status=cdf_globalinfo(input,&data_info,verbose);
  if(status!=0) return(-1);
  
  if(verbose) {
    for (v=0;v<data_info.nvarsp;v++) {
      printf("variable %3d: name %s, type %d,ndim %d \n",v,(data_info.variable[v]).name,data_info.variable[v].type,data_info.variable[v].ndim);
      }
    }
  status=cdf_varinfo(input,"bathymetry",&h_info,verbose);
  if(status!=0) return(-1);

  status=cdf_globalinfo(input,&grid_info,verbose);
  if(status!=0) return(-1);
  status= poc_getgrid2d (input, grid_info, h_info, grid);
  if(status!=0) return(-1);

  *buffer=new float[grid->nx*grid->ny];
  status= poc_getvar2d (input, h_info.id, frame, *buffer, mask, h_info);
  if(status!=0) return(-1);

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int map_loadgrid(const char *input,const char *variable, grid_t *grid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int    v,status;
  int    verbose;
  cdfgbl_t data_info,grid_info;
  cdfvar_t h_info;

  verbose=0;
  status=cdf_globalinfo(input,&data_info,verbose);
  for (v=0;v<data_info.nvarsp;v++) {
    printf("variable %3d: name %s, type %d,ndim %d \n",v,(data_info.variable[v]).name,data_info.variable[v].type,data_info.variable[v].ndim);
    }
  status=cdf_varinfo(input,variable,&h_info,1);
  if(status!=0) return(-1);

  status=cdf_globalinfo(input,&grid_info,verbose);
  status= poc_getgrid2d (input, grid_info, h_info, grid);
  if(status!=0) return(-1);

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int map_loadfield(const char *filename,grid_t *grid, float **buffer, float *mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;

  if(!strrncmp(filename,".grd")) {
    status=map_loadfield_grd(filename,grid, buffer, mask);
    }
  else if(!strrncmp(filename,".nc")) {
    status=map_loadfield_cdf(filename,grid, buffer, mask);
    if(status!=0) return(-1);
//     for(m=0;m<grid->ny*grid->nx;m++) if((*buffer)[m]!=*mask) (*buffer)[m]*=-1.;
    }
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int topo_interpolate(vector<string> databases, grid_t grid, float * & topo, float mask, int average)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,j,k,n, status;
  grid_t topo_grid;
  float  *topo_base,topo_mask, z;
  double x,y;
  
/* *-----------------------------------------------------------------------------
  read and interpolate topo database */
  topo=new float[grid.Hsize()];
  mask=1.e+35;
  for(j=0;j<grid.ny;j++) {
    for(i=0;i<grid.nx;i++) {
      k=j*grid.nx+i;
      topo[k]=mask;
      }
    }

  for(n=0;n<databases.size();n++) {
/* *------------------------------------------------------------------------------
    Load bottom topography*/
//     printf("#################################################################\n");
    printf("#----------------------------------------------------------------\n");
    printf("load bathymetry file: %s\n",databases[n].c_str());
    status=map_loadfield(databases[n].c_str(),&topo_grid,&topo_base,&topo_mask);
    if(status !=0) goto error;
/* *------------------------------------------------------------------------------
    Interpolate and store topo on symphonie grid*/
    for(j=0;j<grid.ny;j++) {
      for(i=0;i<grid.nx;i++) {
        k=j*grid.nx+i;
        x=grid.x[k];
        y=grid.y[k];
        switch (average) {
          case 0:
            status=map_interpolation(topo_grid, topo_base,topo_mask,x,y,&z);
            break;
          case 1:
            status=map_average(topo_grid, topo_base,topo_mask,grid,i,j,&z);
            break;
          }
        if (z!=topo_mask) {
          topo[k]=z;
          }
//         else {
//           topo[k]=mask;
//           }
        }
      }
    topo_grid.free();
    delete[] topo_base;
    }
  return(0);
  
error:
  printf("file loading error, return with bad status\n");
  return(-1);

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void print_help(char *prog_name)

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
    "  %s version " VERSION " " REVISION "\n", prog_name);
  printf("\n"
    "USE\n"
    "  %s [OPTIONS] notebookFile\n",prog_name);
  printf("\n"
    "DESCRIPTION\n"
/* USE FIND AND REPLACE TO EDIT THE LINE BELOW, SO THAT THE FILE HEAD IS ALSO UPDATED */
    "  Set landmask and interpolate topography for a given grid.\n"
    "\n"
    "OPTIONS :\n"
    "  -h,--help  Show this help and exit.\n"
    "  -inv  switch off altitude to depth conversion\n"
    "  -average/-pointwise  switches between averaging/pointwise interpolations of the bathymetry. Default: pointwise.\n"
    "  -zmin  topography values will be contrained above this minimum. Default:-99999.\n"
    "  -format  fprintf format for the symphonie ASCII file, preferably with leading blank. Default: \" %%6.1f\".\n"
    "  -r  followed by the root name\n"
    "  -b  followed by the path to a bathymetry. Several bathymetries are allowed.\n"
    "  -g  followed by the path to the grid file. Disables -n.\n"
    "  -v  followed by the grid variable name\n"
    "  -n  followed by the path to the notebook file\n"
    "  -p  followed by a path to a land mask line. Disables -s. If not given, compute one from the coast line (see below).\n"
    "  -s  followed by the path to a coast line\n"
    );
  printf("\n"
    "FILE FORMATS :\n");
  plg_print_formats();
  /** \endcode */
}

extern ostringstream assertionCmd;

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  geo_t projection;
  float  mask=1.e+10f;
  int i,j,k,n,status;
  FILE *out;
  char *keyword;
  int flag;
  char *rootname=NULL,*output=NULL;
  vector<string> input;
  char *notebook=NULL,*poly=NULL,*gridfile=NULL,*format=NULL,*varname=NULL;
  vector<string> bathymetry;
  grid_t cgrid,sgrid;

  float  *landmask,*topo;
  cdfgbl_t data_info,grid_info;
  cdfvar_t h_info;
  int niterations=25000;

  double *dx, *dy;
  
  bool do_landmask=false, do_bathymetry=false;
  char *vmask=0, *vtopo=0;
  poc_global_t       meta_file;
  poc_data_t<float>  meta_H_z, meta_H_f;

  string ShorelinesFile, TmpPoly, PointString;

  bool smooth=false;
  bool debug=false;
  
  fct_echo( argc, argv);
  
  if(argc<2) {
    print_help(argv[0]);
    exit(0);
    }
  
  n=1;
  while (n < argc) {
    keyword=argv[n];
    switch (keyword[0]) {
      case '-':
        if(strcmp(keyword,"--debug")==0) {
          debug=true;
          n++;
          continue;
          }
        if(strcmp(keyword,"--smooth")==0) {
          smooth=true;
          n++;
          continue;
          }
        if(strcmp(keyword,"-format")==0) {
          format=strdup(argv[n+1]);
          n++;
          n++;
          continue;
          }
        if(strcmp(keyword,"-iterations")==0) {
          sscanf(argv[n+1],"%d",&niterations);
          n++;
          n++;
          continue;
          }
        if(strcmp(keyword,"-vmask")==0) {
          vmask=strdup(argv[n+1]);
          do_landmask=false;
          n++;
          n++;
          continue;
          }
        if(strcmp(keyword,"-vtopo")==0) {
          vtopo=strdup(argv[n+1]);
          do_bathymetry=false;
          n++;
          n++;
          continue;
          }
        switch (keyword[1]) {
          case 'r' :
            rootname= strdup(argv[n+1]);
            n++;
            n++;
            break;

          case 'b' :
            bathymetry.push_back ((string) argv[n+1]);
            n++;
            n++;
            break;

          case 'g' :
            gridfile= strdup(argv[n+1]);
            n++;
            n++;
            break;

          case 'n' :
            notebook= strdup(argv[n+1]);
            n++;
            n++;
            break;

          case 'p' :
            poly= strdup(argv[n+1]);
            n++;
            n++;
            break;

          case 's' :
            ShorelinesFile=argv[n+1];
            n++;
            n++;
            break;
        
          case 'v' :
            varname= strdup(argv[n+1]);
            n++;
            n++;
            break;

          default:
            STDOUT_BASE_LINE("unknown option %s\n",keyword);
              exit(-1);
          }
        break;

      default:
        if(notebook==NULL) {
	  string tmp=keyword;
          input.push_back(tmp);
          n++;
          }
        else {
          STDOUT_BASE_LINE("unknown option %s\n",keyword);
          exit(-1);
          }
        break;
      }
    }

  output=(char *)malloc(1024);

  if(rootname==NULL) {
    printf("rootname for outputs not specified, use <symphonie>\n");
    rootname= strdup("symphonie");
    }
  else {
    printf("rootname for outputs : %s\n",rootname);
    }

  if(notebook && !gridfile) {
/*-----------------------------------------------------------------------------
    Read notebook data */
    printf("#################################################################\n");
    printf("load notebook file: %s\n",notebook);
    status=load_notebook(notebook, &cgrid, &sgrid, &projection);
    if(status) TRAP_ERR_EXIT(status,"load_notebook(\"%s\",) error\n",notebook);
    printf("notebook successfuly processed\n",notebook);
    }

  if(gridfile) {
    printf("#################################################################\n");
    printf("load grid file: %s\n",gridfile);
    status=poc_getgrid2d(gridfile,varname,&sgrid);
    if(status) TRAP_ERR_EXIT(status,"poc_get*grid*(\"%s\",\"%s\",) error\n",gridfile,varname);
    printf("grid successfuly processed\n");
    double reference_lat=(sgrid.ymax-sgrid.ymin)/2.0;
    double reference_lon=(sgrid.xmax-sgrid.ymin)/2.0;
    projPJ projection=assign_StereoOblique(reference_lat,reference_lon);
    cgrid=map_get_cartesian(projection,sgrid);
    }
   
/*-----------------------------------------------------------------------------
  compute grid resolution (cartesian) */
  float *cresolution;
  
  status=map_cartesian_resolution(cgrid, &dx, &dy);
  cresolution=new float[cgrid.Hsize()];
  for(int m=0;m<cgrid.Hsize();m++) {
    cresolution[m]=sqrt(dx[m]*dy[m]);
    }
  delete[] dx;
  delete[] dy;
  
/*-----------------------------------------------------------------------------
  compute grid resolution (spherical) */
  float *sresolution;

  status=map_spherical_resolution(sgrid, &dx, &dy);
  sresolution=new float[cgrid.Hsize()];
  for(int m=0;m<cgrid.Hsize();m++) {
    sresolution[m]=sqrt(dx[m]*dy[m]);
    }
  delete[] dx;
  delete[] dy;
  
  plg_t GridLimits=plg_GridLimits(sgrid); 
  TmpPoly=(string)rootname +"-grid.plg";
  status=plg_save(TmpPoly.c_str(), PLG_FORMAT_SCAN, &GridLimits, 1);
  
/**------------------------------------------------------------------------------
   */
  {
  poc_var_t zgridVar, fgridVar;
  size_t nv2Dvalues=sgrid.Hsize();
  
  string filename=(string)rootname+"-parameters.nc";
  printf("#################################################################\n");
  printf("create : %s\n",filename.c_str());
  
  poc_data_t<float> meta_m;
  status=meta_m.init(gridfile,"mask_t");
  status=meta_m.read_data(gridfile,0,0,1);
    
  landmask=0;
  if(status==0) {
    int level=meta_m.info.dimensions[0].len-1;
    landmask=aset(sgrid.Hsize(),0.f);
    for(size_t m=0; m< nv2Dvalues; m++) {
      landmask[m]=meta_m.data[level*nv2Dvalues+m];
      }
    }
  else {
    landmask=aset(sgrid.Hsize(),-1.f);
    }
  
  vector<plg_t> MaskGridLimits=plg_GridLimits(sgrid, landmask, debug, 0);
  TmpPoly=(string)rootname +"-maskgrid.plg";
  status=plg_save(TmpPoly.c_str(), PLG_FORMAT_SCAN, MaskGridLimits);
  
  float nu=500.0;
  
/**------------------------------------------------------------------------------
  initialise nc file */
  zgridVar=poc_save_grid(filename.c_str(), sgrid, __FILE__,__LINE__,"t");
  
  poc_var_t var=zgridVar;
  float *buffer=new float[sgrid.Hsize()];
  float *smoothed=new float[sgrid.Hsize()];

//   string polygons="z0.nei",values="z0.dat";
//   vector<float> v =Polygons::initialiseValue<double, float>(polygons, values, sgrid, 1.e-03);
//   for (size_t i = 0; i < v.size(); i++) buffer[i] = v[i];
//   v.clear();
//   
//   for(size_t m=0; m< nv2Dvalues; m++) {
//     if(landmask[m]==0) buffer[m] = mask;
//     }
//   
//   status=map_diffusion_forward(sgrid, buffer, mask, nu, smoothed, niterations);
//   
//   var.init("z0-bottom_t",NC_FLOAT,"z0","m", mask);
//   var.attributes<<poc_att_t("content","YX");
//   status=poc_put_vara(filename.c_str(),var,0,buffer,0);
//   
//   var.init("z0-bottom_diffused_t",NC_FLOAT,"z0","m",mask);
//   var.attributes<<poc_att_t("content","YX");
//   status=poc_put_vara(filename.c_str(),var,0,smoothed,0);
//   
//   polygons="par.nei";
//   values="par.dat"; 
//   v =Polygons::initialiseValue<double, float>(polygons, values, sgrid, 20.0);
//   for (size_t i = 0; i < v.size(); i++) buffer[i] = v[i];
//   v.clear();
//   
//   for(size_t m=0; m< nv2Dvalues; m++) {
//     if(landmask[m]==0) buffer[m] = mask;
//     }
//   
//   status=map_diffusion_forward(sgrid, buffer, mask, nu, smoothed, niterations);
//   
//   var.init("par_t",NC_FLOAT,"par","m",mask);
//   var.attributes<<poc_att_t("content","YX");
//   status=poc_put_vara(filename.c_str(),var,0,buffer,0);
//   
//   var.init("par_diffused_t",NC_FLOAT,"par","m",mask);
//   var.attributes<<poc_att_t("content","YX");
//   status=poc_put_vara(filename.c_str(),var,0,smoothed,0);
  
  for(int k=0;k<input.size();k++) {
    vector<string> tmp=string_split(input[k], " ");
    string polygons=tmp[0],values=tmp[1];
    vector<float> v =Polygons::initialiseValue<double, float>(polygons, values, sgrid, 1.e-03);
    for (size_t i = 0; i < v.size(); i++) buffer[i] = v[i];
    v.clear();
  
    for(size_t m=0; m< nv2Dvalues; m++) {
      if(landmask[m]==0) buffer[m] = mask;
      }
  
    status=map_diffusion_forward(sgrid, buffer, mask, nu, smoothed, niterations);
    
    string varname;
    varname=tmp[2]+"_t";
    var.init(varname.c_str(),NC_FLOAT,varname.c_str(),"m", mask);
    var.attributes<<poc_att_t("content","YX");
    status=poc_put_vara(filename.c_str(),var,0,buffer,0);
    
    varname=tmp[2]+"-diffused_t";
    var.init(varname.c_str(),NC_FLOAT,varname.c_str(),"m", mask);
    var.attributes<<poc_att_t("content","YX");
    status=poc_put_vara(filename.c_str(),var,0,smoothed,0);
    }
    
  delete[] buffer;
  delete[] smoothed;
  }
  
  cgrid.free();
  sgrid.free();

  printf("#################################################################\n");
  STDOUT_BASE_LINE("symtools successfuly completed\n");
  
  return 0;
}
