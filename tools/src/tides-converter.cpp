
/*******************************************************************************

  T-UGO tools, 2006-2014

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

<!-- USE FIND AND REPLACE TO EDIT THIS LINE, SO THAT print_help IS ALSO UPDATED -->
\brief Converts and interpolates structured tidal atlases

<!-- A LINK TO main() or print_help() WILL NOT LINK TO THE RIGHT SOURCE ! -->
See the main <a href=#func-members>function</a> for how this works
and the print_help <a href=#func-members>function</a> for how to use this.
*/
/*----------------------------------------------------------------------------*/

#define MAIN_SOURCE

#include "version-macros.def" //for VERSION and REVISION

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
 

#include "config.h"
#include "tools-structures.h"

#include "rutin.h"    /*  rutin.h contains common utility routines  */
#include "map.h"
#include "poc-netcdf-data.hpp"
#include "bmg.h"
#include "xyz.h"
#include "ascii.h"
#include "topo.h"
#include "tides.h"
#include "statistic.h"
#include "swap.h"
#include "grd.h"

#include "fe.h"

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int netcdf2ascii(mesh_t mesh, vector<string*> varnames, char **wave, int nwave, int iteration, vector<string> input, vector<string>  discretisation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  extern toponyms_t paire_map(void);
  
  int k, status,id;
  char amplitude[256],phase[256];
  discretisation_t *descriptor;
  complex<float> *UGbufC=0, *SGbufC=0;
  float          *UGbufR=0, *SGbufR=0;
  bool is_complex=false;

  poc_var_t input_var[2], av, gv;
  
  toponyms_t Pmap=paire_map();
  
  int nvars=0;


/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  read mesh geometry
  
*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if(mesh.nvtxs==0) {
    printf("#################################################################\n");
    printf("load mesh geometry from %s\n",input[0].c_str());
    status=fe_readgeometry(input[0].c_str(), &mesh);
    if(status) NC_TRAP_ERROR(return,status,1,"fe_readgeometry(\"%s\",) error",input[0].c_str());
    }


/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  identify discretisation
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  for(int v=0;v<varnames.size();v++) {
    id=discretisation_from_name(discretisation[v].c_str(), 0);
    if(id<0) TRAP_ERR_RETURN(status,1,"discretisation_id() error\n");
    printf("#################################################################\n");
    printf("treating %s, %s discretisation\n",varnames[v][0].c_str(), discretisation_name(id));
    status=discretisation_init(&mesh,id,0);
    if(status) NC_TRAP_ERROR(return,status,1,"discretisation_init(,%d,0) error",id);
    
    descriptor=get_descriptor_address(mesh,id);
    if(descriptor->nnodes==0) TRAP_ERR_RETURN(-1,1,"get_descriptor() error\n");
        
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
    setup variable names and variable type (real or complex)
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

    if(varnames[v][1]=="") {
      sprintf(amplitude,"%s",varnames[v][0].c_str());
      UGbufR=new float[descriptor->nnodes];
      is_complex=false;
      nvars=1;
      }
    else {
      sprintf(amplitude,"%s",varnames[v][0].c_str());
      sprintf(phase,    "%s",varnames[v][1].c_str());
      UGbufC=new complex<float>[descriptor->nnodes];
      is_complex=true;
      nvars=2;
      }
    
    for (k=0;k<nwave;k++) {
      printf("treating %s from %s\n",varnames[v][0].c_str(), input[k].c_str());
      
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
      load variable
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

      if(is_complex) {
          status=poc_inq_var(input[k], amplitude, &input_var[0]);
          if(status!=0) TRAP_ERR_EXIT(status,"read error\n");
          status=poc_inq_var(input[k], phase,     &input_var[1]);
          status=poc_get_UG3D(input[k].c_str(), iteration, amplitude, phase, UGbufC);
//           status=quoddy_savec1(const char *name, int nndes, fcomplex *buffer, char **comment);
        }
      else {
        float inputMask;
          status=poc_inq_var(input[k], amplitude, &input_var[0]);
          status=poc_decode_mask(input_var[0], 0, 0, &inputMask);
          status=poc_get_UG3D(input[k].c_str(), iteration, amplitude, UGbufR);
        }
      
      if(status!=0) {
        printf("elevation not done\n");
        continue;
        }
      
      }
     
    }
 
  if(SGbufR!=0) delete[] SGbufR;
  if(UGbufR!=0) delete[] UGbufR;
  if(SGbufC!=0) delete[] SGbufC;
  if(UGbufC!=0) delete[] UGbufC;
  
  descriptor->destroy();
  
  Pmap.clear();
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  grid_t get_zonegrid_1_8(char *zone)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{

  bool zone_initialised;
  grid_t grid=get_zonegrid(zone, &zone_initialised);

//   if(strcmp(zone,"kerguelen")==0)
//   {
//   map_set2Dgrid(&grid,40.0,-75.0,95.0,-33.0);
//   zone_initialised=1;
//   }
 
  grid.dx  =    1./8.;
  grid.dy  =    1./8.;
  grid.nx  = (int)( (grid.xmax-grid.xmin)/grid.dx )+1;
  grid.ny  = (int)( (grid.ymax-grid.ymin)/grid.dy )+1;
  
  return(grid);

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int otis_loadc1(const char* filename, grid_t *grid,fcomplex ***tide, fcomplex *mask, int *nbuffers, spectrum_t *s)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  {
  int i,k,m,n,size,status;
  int nwaves,cdim=4;
  FILE *in;
  range_t<double> range;
  char *names;
  size_t offset, nread;
  float th_lim[2],ph_lim[2];

#define NEED_SWAP  1
  
  in = fopen(filename, "rb");
  offset=4;


/*------------------------------------------------------------------------------
  get binary blocksize */
  nread=fread(&size, sizeof(int), 1, in );
//      read(funit)n,m,nc,th_lim,ph_lim,c_id(1:nc)
  
  nread=fread(&n, sizeof(int), 1, in );
  nread=fread(&m, sizeof(int), 1, in );
  nread=fread(&nwaves, sizeof(int), 1, in );
  nread=fread(&th_lim[0], sizeof(float), 1, in );
  nread=fread(&th_lim[1], sizeof(float), 1, in );
  nread=fread(&ph_lim[0], sizeof(float), 1, in );
  nread=fread(&ph_lim[1], sizeof(float), 1, in );
  
#if NEED_SWAP == 1
  size=lswap(size);
  n=lswap(n);
  m=lswap(m);
  nwaves=lswap(nwaves);
  th_lim[0]=fswap(th_lim[0]);
  th_lim[1]=fswap(th_lim[1]);
  ph_lim[0]=fswap(ph_lim[0]);
  ph_lim[1]=fswap(ph_lim[1]);
#endif

  names=new char[nwaves*cdim];
  nread=fread(names, sizeof(char), cdim*nwaves, in );
  
  fseek(in,offset,SEEK_CUR);
  
  grid->xmin= th_lim[0]*1000;
  grid->ymin=-ph_lim[1]*1000;
  grid->xmax= th_lim[1]*1000;
  grid->ymax=-ph_lim[0]*1000;
  
  grid->nx=m;
  grid->ny=n;
  
  grid->dx=5000.;
  grid->dy=5000.;
  
  double ratio=1.0125;
  
  grid->xmin*=ratio;
  grid->ymin*=ratio;
  grid->xmax*=ratio;
  grid->ymax*=ratio;
  
  grid->dx=(grid->xmax-grid->xmin)/((double) grid->nx-1.);
  grid->dy=(grid->ymax-grid->ymin)/((double) grid->ny-1.);
  
  
  grid->modeH=0;
  
  status=map_completegridaxis(grid,2);
  
  map_minmax(grid);
  
  map_printgrid(*grid);
  
//   grid->dx=(grid->xmax-grid->xmin)/((double) grid->nx-1.);
//   grid->dy=(grid->ymax-grid->ymin)/((double) grid->ny-1.);

//   grid->modeH=2;
  
  *tide=new fcomplex *[nwaves];
  *mask=complex<float>(0.,0.);
  size=grid->nx*grid->ny;
  for(k=0;k<nwaves;k++) {
    (*tide)[k]=new fcomplex[grid->Hsize()];
    nread=fread(&size, sizeof(int), 1, in );
    nread=fread((*tide)[k], sizeof(fcomplex), grid->Hsize(), in);
    fseek(in,offset,SEEK_CUR);
#if NEED_SWAP == 1
    for(m=0;m<grid->Hsize();m++) (*tide)[k][m]=cswap((*tide)[k][m]);
#endif
    status=swap_XY2YX(*grid, (*tide)[k]);
    status=grd_mirror(*grid, grid->Hsize(), (*tide)[k], *mask);
    }
  
  char parameters[1024];
  projPJ pj=assign_StereoNorth(90.0, 90, 0.0, parameters);
  grid->projection=pj;
  grid->proj4options=poc_strdup(parameters);

  fclose(in);
  
  
  *nbuffers=nwaves;
  
  s->n=nwaves;
  s->waves=new tidal_wave[s->n];
  for(k=0;k<nwaves;k++) {
    char tmp[4];
    strncpy(tmp,(const char*) &(names[k*cdim]),4);
    for(i=0;i<cdim;i++) if(tmp[i]==' ') tmp[i]=0;
    strcpy(s->waves[k].name,tmp);
    }
  return(0);
  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int tpxo6_loadc1(const char* filename, grid_t *grid,fcomplex ***tide, fcomplex *mask, int *nbuffers, spectrum_t *s)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,k,m,n,size,status;
  int nwaves,cdim=4;
  FILE *in;
  range_t<double> range;
  char *names;
  size_t offset, nread;
  float th_lim[2],ph_lim[2];

#define NEED_SWAP  1
  
  in = fopen(filename, "rb");
  offset=4;


/*------------------------------------------------------------------------------
  get binary blocksize */
  nread=fread(&size, sizeof(int), 1, in );
//      read(funit)n,m,nc,th_lim,ph_lim,c_id(1:nc)
  
  nread=fread(&m, sizeof(int), 1, in );
  nread=fread(&n, sizeof(int), 1, in );
  nread=fread(&nwaves, sizeof(int), 1, in );
  nread=fread(&ph_lim[0], sizeof(float), 1, in );
  nread=fread(&ph_lim[1], sizeof(float), 1, in );
  nread=fread(&th_lim[0], sizeof(float), 1, in );
  nread=fread(&th_lim[1], sizeof(float), 1, in );
  
#if NEED_SWAP == 1
  size=lswap(size);
  n=lswap(n);
  m=lswap(m);
  nwaves=lswap(nwaves);
  th_lim[0]=fswap(th_lim[0]);
  th_lim[1]=fswap(th_lim[1]);
  ph_lim[0]=fswap(ph_lim[0]);
  ph_lim[1]=fswap(ph_lim[1]);
#endif

  names=new char[nwaves*cdim];
  nread=fread(names, sizeof(char), cdim*nwaves, in );
  
  fseek(in,offset,SEEK_CUR);
  
  grid->xmin= th_lim[0];
  grid->ymin=-ph_lim[1];
  grid->xmax= th_lim[1];
  grid->ymax=-ph_lim[0];
  
  grid->nx=m;
  grid->ny=n;
  
//   grid->dx=5000.;
//   grid->dy=5000.;
//   
//   double ratio=1.0125;
//   
//   grid->xmin*=ratio;
//   grid->ymin*=ratio;
//   grid->xmax*=ratio;
//   grid->ymax*=ratio;
  
  grid->dx=(grid->xmax-grid->xmin)/((double) grid->nx-1.);
  grid->dy=(grid->ymax-grid->ymin)/((double) grid->ny-1.);
  
  
  grid->modeH=0;
  
  status=map_completegridaxis(grid,2);
  
  map_minmax(grid);
  
  map_printgrid(*grid);
  
//   grid->dx=(grid->xmax-grid->xmin)/((double) grid->nx-1.);
//   grid->dy=(grid->ymax-grid->ymin)/((double) grid->ny-1.);

//   grid->modeH=2;
  
  *tide=new fcomplex *[nwaves];
  *mask=complex<float>(0.,0.);
  size=grid->nx*grid->ny;
  for(k=0;k<nwaves;k++) {
    (*tide)[k]=new fcomplex[grid->Hsize()];
    nread=fread(&size, sizeof(int), 1, in );
    nread=fread((*tide)[k], sizeof(fcomplex), grid->Hsize(), in);
    fseek(in,offset,SEEK_CUR);
#if NEED_SWAP == 1
    for(m=0;m<grid->Hsize();m++) (*tide)[k][m]=cswap((*tide)[k][m]);
#endif
//     status=swap_XY2YX(*grid, (*tide)[k]);
//     status=grd_mirror(*grid, grid->Hsize(), (*tide)[k], *mask);
    }
  
//   char parameters[1024];
//   double lon,lat;
//   projPJ pj=assign_StereoNorth(90.0, 90, 0.0, parameters);
//   grid->projection=pj;
//   grid->proj4options=poc_strdup(parameters);

  fclose(in);
  
  
  *nbuffers=nwaves;
  
  s->n=nwaves;
  s->waves=new tidal_wave[s->n];
  for(k=0;k<nwaves;k++) {
    char tmp[4];
    strncpy(tmp,(const char*) &(names[k*cdim]),4);
    for(i=0;i<cdim;i++) if(tmp[i]==' ') tmp[i]=0;
    strcpy(s->waves[k].name,tmp);
    }
  return(0);

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int tpxo9_loadc1(const char* filename,const char* gridfile, grid_t *grid,fcomplex ***tide, fcomplex *mask, int *nbuffers, spectrum_t *s)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,k,m,status;
  int id,ncid,dimid;
  int nwaves,cdim;
  float *realp, *imagp;
  cdfgbl_t global;
  cdfvar_t info;
  decoded_t decoded;
  range_t<double> range;
  char *names;
  
  status=cdf_globalinfo(gridfile,&global,0);
  
  dimid=cdf_identify_dimension(global,"nx");
  grid->nx=global.dimension[dimid].length;
  dimid=cdf_identify_dimension(global,"ny");
  grid->ny=global.dimension[dimid].length;
  
  grid->x=new double[grid->nx*grid->ny];
  grid->y=new double[grid->nx*grid->ny];
  
  status=nc_open(gridfile,0,&ncid);
  nc_check_error(status,__LINE__,__FILE__);
 
  status=cdf_varinfo(global,"lon_z",&info);

  switch(info.ndim) {
    case 1:
      grid->x=new double[grid->nx];
      status=nc_get_var_double(ncid,info.id,grid->x);
      nc_check_error(status,__LINE__,__FILE__);
      break;
    case 2:
      grid->x=new double[grid->nx*grid->ny];
      status=nc_get_var_double(ncid,info.id,grid->x);
      nc_check_error(status,__LINE__,__FILE__);
      status=swap_XY2YX(*grid, grid->x);
      break;
    }
  
  status=cdf_varinfo(global,"lat_z",&info);

  switch(info.ndim) {
    case 1:
      grid->y=new double[grid->ny];
      status=nc_get_var_double(ncid,info.id,grid->y);
      nc_check_error(status,__LINE__,__FILE__);
      grid->modeH=1;
      break;
    case 2:
      grid->y=new double[grid->nx*grid->ny];
      status=nc_get_var_double(ncid,info.id,grid->y);
      nc_check_error(status,__LINE__,__FILE__);
      status=swap_XY2YX(*grid, grid->y);
      grid->modeH=2;
      break;
    }
  
  status=nc_close(ncid);
  
  map_minmax(grid);
  
  grid->dx=(grid->xmax-grid->xmin)/((double) grid->nx-1.);
  grid->dy=(grid->ymax-grid->ymin)/((double) grid->ny-1.);
  
  status=cdf_globalinfo(filename,&global,0);
  
  status=nc_open(filename,0,&ncid);
  nc_check_error(status,__LINE__,__FILE__);
 
  dimid=cdf_identify_dimension(global,"nc");
  
  if(dimid==-1) {
    nwaves=1;
    }
  else {
    nwaves=global.dimension[dimid].length;
    }
  dimid=cdf_identify_dimension(global,"nct");
  cdim=global.dimension[dimid].length;
  
  status=cdf_varinfo(global,"con",&info);
  
  id=cdf_identify(global,"con");
  
  names=new char[nwaves*cdim];
  status=nc_get_var_text(ncid,id,names);
  nc_check_error(status,__LINE__,__FILE__);
  info.destroy();
  
  
  realp=new float[grid->nx*grid->ny*nwaves];
  imagp=new float[grid->nx*grid->ny*nwaves];
   
  status=cdf_varinfo(global,"hRe",&info);
  
  status=nc_get_var_float(ncid,info.id,realp);
  nc_check_error(status,__LINE__,__FILE__);
  info.destroy();
   
  status=cdf_varinfo(global,"hIm",&info);
  
  status=nc_get_var_float(ncid,info.id,imagp);
  nc_check_error(status,__LINE__,__FILE__);
  info.destroy();
  
  *tide=new fcomplex *[nwaves];
  int size=grid->nx*grid->ny;
  for(k=0;k<nwaves;k++) {
    (*tide)[k]=new fcomplex[grid->nx*grid->ny];
    for(m=0;m<grid->nx*grid->ny;m++) {
      (*tide)[k][m]=complex<float>(realp[k*size+m],imagp[k*size+m]);
      }
    status=swap_XY2YX(*grid, (*tide)[k]);
    }
  status = nc_close(ncid);
  nc_check_error(status,__LINE__,__FILE__);

  *mask=complex<float>(0.,0.);
  
  *nbuffers=nwaves;
  
  s->n=nwaves;
  s->waves=new tidal_wave[s->n];
  for(k=0;k<nwaves;k++) {
    char tmp[4];
    strncpy(tmp,(const char*) &(names[k*cdim]),4);
    for(i=0;i<cdim;i++) if(tmp[i]==' ') tmp[i]=0;
    strcpy(s->waves[k].name,tmp);
    }
  
  return(0);
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
    " %s version " VERSION " " REVISION "\n", prog_name);
  printf("\n"
    "USE\n"
    "  %s input [OPTIONS]\n",prog_name);
  printf("\n"
    "DESCRIPTION\n"
/* USE FIND AND REPLACE TO EDIT THE LINE BELOW, SO THAT THE FILE HEAD IS ALSO UPDATED */
    "  Converts and interpolates structured tidal atlases.\n"
    "\n"
    "OPTIONS\n"
    "  -h,--help : show this help and exit\n"
    "  -cm : input is in cm instead of m\n"
    "  -mm : input is in mm instead of m\n"
    "  -bmg : set output format to bmg\n"
    "  -m : followed by mask file\n"
    "  -b : followed by topography file. Usefull to use instead of mask file.\n"
    "  -debug : enable debug mode when reading topography or mask file\n"
    "  -o : followed by output file name\n"
    "  -w : followed by wave name. Targeted wave if multi-waves input\n"
    "  -f : followed by input format : bmg, ascii, netcdf, got, tpxo, otis, dtu (equivalent to ascii), eot, osu or cst\n"
    "    If the extension of the input is .nc, then the default is netcdf,\n"
    "    if the extension of the input is .asc, then the default is ascii,\n"
    "    otherwise it is bmg.\n"
    "  -uv : Read 2 buffers from atlas. Implies ascii input format.\n"
    "  -s : followed by smoothing factor in m. Default is no smoothing.\n"
    "\n"
    "The following options are for netcdf format.\n"
    "\n"
    "  -iv,--input-variables : followed by 2 input variable names for amplitude and phase. Default: Ha Hg\n"
    "  -ig,--input-grid : followed by input grid file. Default is the input file.\n"
    "  -igv,--input-grid-variable : followed by grid variable. Default is amplitude variable.\n"
    "\n"
    "The following options are for modifying the ouput grid. By default, the ouput grid is the same as the input grid.\n"
    "\n"
    );
  print_zone_arg_help();
  printf(
    "  -g : followed by grid file. The default grid is that of the input.\n"
    "  -v : followed by grid variable. For more information see option -g .\n"
    "  -z : followed by grid or zone name or resolution in degrees\n"
    "  -extend : add one column to repeat 0/360 longitude\n"
    "\n"
    "SEE ALSO: gridit-cdf\n"
    ); /** \endcode */
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double x,y;
  float  time=0,scale=1.0;
  float  spec=-9999;
  int i,j,k,n;
  int nbuffers=1,status;
  const char *keyword=NULL,*input=NULL, 
    *input_grid=NULL,*input_grid_var=NULL,
    *zone=NULL,*gridfile=NULL,*gridvar=NULL,*topofile=NULL,*maskfile=NULL;
  frame_t prescribed;
  double dx=NAN,dy=NAN;
  char *output=NULL;
  char *wave=0;
  string format;
  grid_t grid,atlas_grid, topogrid;
  fcomplex **tide=NULL,cmask,z,e;
  spectrum_t spectrum;
  bool remap=1,extend=false,smooth=0, native=true;
  float lscale;
  int bmg=0;
  string Aname="Ha", Gname="Hg";
  
  bool debug=false;
  
  string cmd=fct_echo( argc, argv);

  n=1;
  while (n < argc) {
    keyword=argv[n];
    if(read_zone_arg(keyword,argv[n+1],&prescribed,&dx,&dy)) {
      n++;
      n++;
      native=false;
      continue;
      }
    if( strcmp(keyword,"--help")==0 or
        strcmp(keyword,"-h")==0 ) {
      print_help(argv[0]);
      exit(0);
      }
    if(strcmp(keyword,"-debug")==0) {
      debug=true;
      n++;
      continue;
      }
    if(strcmp(keyword,"-uv")==0) {
      nbuffers=2;
      format="ascii";
      n++;
      continue;
      }
    if(strcmp(keyword,"-cm")==0) {
      scale=0.01;
      n++;
      continue;
      }
    if(strcmp(keyword,"-mm")==0) {
      scale=0.001;
      n++;
      continue;
      }
    if(strcmp(keyword,"-bmg")==0) {
      bmg=1;
      n++;
      continue;
      }
    if( strcmp(keyword,"--input-variables")==0 or
        strcmp(keyword,"-iv")==0 ) {
      Aname=argv[n+1];
      Gname=argv[n+2];
      n+=2;
      n++;
      continue;
      }
    if( strcmp(keyword,"--input-grid")==0 or
        strcmp(keyword,"-ig")==0 ) {
      input_grid=argv[n+1];
      n++;
      n++;
      continue;
      }
    if( strcmp(keyword,"--input-grid-variable")==0 or
        strcmp(keyword,"-igv")==0 ) {
      input_grid_var=argv[n+1];
      n++;
      n++;
      continue;
      }
    if(strcmp(keyword,"-extend")==0) {
      extend=1;
      n++;
      continue;
      }
    switch (keyword[0]) {
      case '-':
        switch (keyword[1]) {

        case 'g' :
          gridfile= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'v' :
          gridvar= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'z' :
          zone= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'm' :
          maskfile= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'o' :
          output= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'w' :
          wave=strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'f' :
          format=argv[n+1];
          n++;
          n++;
          break;

        case 'b' :
          topofile= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 's' :
          sscanf(argv[n+1],"%f",&lscale);
          n++;
          n++;
          smooth=1;
          break;

        default:
          STDOUT_BASE_LINE("unknown option %s\n",keyword);
          print_help(argv[0]);
          wexit(-1);
        }
        break;

      default:
        input= strdup(argv[n]);
        n++;
        break;
      }
    }

  if(input == NULL){
    STDOUT_BASE_LINE("no input given\n");
    print_help(argv[0]);
    wexit(-1);
    }
  
  if(input_grid==0){
    input_grid=input;
    }
  
  if(input_grid_var==0){
    input_grid_var=Aname.c_str();
    }

  if(gridvar!=0){
    if(gridfile==0)
      printf("*** WARNING: GRID VARIABLE %s GIVEN WITHOUT A GRID FILE ! ***\n",gridvar);
    }

//  rmask=spec;
  cmask=fcomplex(spec,spec);

  if(format==""){//get format from default
    const string extension=get_extension(input);
    if(extension==""){
      STDOUT_BASE_LINE("Could not guess format of input file %s\n",input);
      print_help(argv[0]);
      wexit(-1);
      }
    else if(extension=="nc"){
      format="netcdf";
      }
    else if(extension=="asc"){
      format="ascii";
      }
    else{
      format="bmg";
      }
    }
  
  if(format=="dtu")
    format="ascii";
  
  printf("#################################################################\n");
  printf("load input file : %s, format "+format+"\n",input);
  
  if(format=="bmg") {
    tide=new fcomplex *[1];
    tide[0]=0;
    spectrum.n=1;
    status= bmg_loadc1(input,1,1,1,&atlas_grid,&tide[0],&cmask,&time);
    nbuffers=1;
    }
  else if(format=="ascii") {
    spectrum.n=1;
    tide=new fcomplex *[nbuffers];
    if(nbuffers==1)
      status=ascii_loadc2(input,&atlas_grid,&tide[0],0,&cmask);
    else
      status=ascii_loadc2(input,&atlas_grid,&tide[0],&tide[1],&cmask);
    }
  else if(format=="netcdf") {
    tide=new fcomplex *[1];
    spectrum.n=1;
    status=poc_get_grid(input_grid,input_grid_var,&atlas_grid);
    NC_TRAP_ERROR(wexit,status,1,"poc_get_grid(\"%s\",\"%s\",) error",input_grid,input_grid_var);
    tide[0]=new fcomplex[atlas_grid.Hsize()];
    status=poc_get_cvara(input,Aname,Gname,0,tide[0]);
    NC_TRAP_ERROR(wexit,status,1,"poc_get_cvara(\"%s\",\""+Aname+"\",\""+Gname+"\",) error",input);
    cmask=NC_FILL_COMPLEX;
    nbuffers=1;
    }
  else if(format=="got") {
    tide=new fcomplex *[1];
    spectrum.n=1;
//     status= ascii_loadc1_got_complex(input,&atlas_grid,&(tide[0]),&cmask);
    status= ascii_loadc1_got(input,&atlas_grid,&(tide[0]),&cmask);
    nbuffers=1;
    }
  else if(format=="tpxo") {
    status=tpxo_loadc1(input,&atlas_grid,&tide,&cmask,&nbuffers,&spectrum);
//    spectrum.n=1;
    }
  else if(format=="tpxo9") {
    status=tpxo9_loadc1(input, "grid_tpxo9.nc", &atlas_grid,&tide,&cmask,&nbuffers,&spectrum);
//    spectrum.n=1;
    }
  else if(format=="tpxo6") {
    status=tpxo6_loadc1(input,&atlas_grid,&tide,&cmask,&nbuffers,&spectrum);
    }
  else if(format=="otis") {
    status=otis_loadc1(input,&atlas_grid,&tide,&cmask,&nbuffers,&spectrum);
    }
  else if(format=="eot") {
    tide=new fcomplex *[1];
    spectrum.n=1;
    status= eot_loadc1(input,&atlas_grid,&(tide[0]),&cmask);
    nbuffers=1;
    }
  else if(format=="osu") {
    tide=new fcomplex *[1];
    spectrum.n=1;
    status= osu_loadc1(input,&atlas_grid,&(tide[0]),&cmask);
    nbuffers=1;
    }
  else if(format=="cst") {
    int pos;
    if(wave!=0) {
      sscanf(wave,"%d",&pos);
      }
    else {
      pos=-1; //
      }
   status= shomcst_loadc1(input,&atlas_grid,pos,&tide,&cmask,&spectrum);
   nbuffers=spectrum.n;
   }
  else if(format=="boy") {
//   printf("#################################################################\n");
//   printf("load grid from file : %s\n",input);
//     int ncolumn=2;
//     status=xyz_loadgrid(input,0, &atlas_grid,&ncolumn);
//     if(status !=0) {
//       STDOUT_BASE_LINE("cannot load grid from file=%s %d\n",input,status);
//       TRAP_ERR_EXIT(-1,"exiting\n");
//       }
    tide=new fcomplex *[1];
    spectrum.n=1;
    status= boy_loadc1(input,&atlas_grid,&(tide[0]),&cmask);
    nbuffers=1;
    }
  else {
    STDOUT_BASE_LINE("Format "+format+" not recognised.\n");
    print_help(argv[0]);
    wexit(-1);
    }

  if(status != 0) NC_TRAP_ERROR(wexit,status,1,"error while reading %s",input);

//  status=map_completegridaxis(&atlas_grid,2);

  if(extend) {
    int modified=0;
    if(atlas_grid.modeH!=2)
      status=map_completegridaxis(&atlas_grid,2);
    if(nbuffers==1) {
      status=map_extendfield(&atlas_grid, &(tide[0]), &modified);
      }
    else {
      status=map_extendfield(&atlas_grid, tide, nbuffers, &modified);
      }
    }
  
  if(topofile!=NULL) {
    printf("#################################################################\n");
    printf("additional mask prescription from %s\n",topofile);
    float value;
    short rmask,*topo=NULL;
//     status=bmg_loadgrid (topofile,&topogrid);
//     status=map_completegridaxis(&topogrid,2);
//     topo=new float[topogrid.nx*topogrid.ny];
//     status= bmg_loadr1(topofile,1,1,1,topogrid,topo,&rmask,&time);
    status=topo_loadfield(topofile, &topogrid, &topo, &rmask, debug);
    for (j=0;j<atlas_grid.ny;j++) {
      for (i=0;i<atlas_grid.nx;i++) {
        n=i+atlas_grid.nx*j;
        y=atlas_grid.y[n];
        x=map_recale(topogrid,atlas_grid.x[n]);
        status=map_interpolation(topogrid, topo, rmask, x,y,&value);
        if(value>=0){
          for(k=0;k<nbuffers;k++) tide[k][n]=cmask;
          }
        }
      }
    }
  else {
    }

  if(maskfile!=NULL) {
    printf("#################################################################\n");
    printf("additional mask prescription from %s\n",maskfile);
    signed char *buffer,mask;
    float value;
    grid_t maskgrid;
    status=topo_loadfield_cdf((const char*) maskfile, "dst", &maskgrid, &buffer, &mask, debug);
    mask=-2;
    for (j=0;j<atlas_grid.ny;j++) {
      for (i=0;i<atlas_grid.nx;i++) {
        n=i+atlas_grid.nx*j;
        y=atlas_grid.y[n];
        x=map_recale(maskgrid,atlas_grid.x[n]);
        status=map_interpolation(maskgrid, buffer, mask, x,y,&value);
        if(value==0){
          for(k=0;k<nbuffers;k++) tide[k][n]=cmask;
          }
        if(value==-2) {
          for(k=0;k<nbuffers;k++) tide[k][n]=cmask;
          }
        }
      }
    }
  
  if(zone != NULL) {
    grid=get_zonegrid(zone);
    }
  else if(gridfile != NULL){
    if(gridvar)
      status=poc_get_grid(gridfile,gridvar,&grid);
    else
      status=cdf_loadvargrid (gridfile,0,&grid);
    }
  else{
    remap=0;
    grid=atlas_grid;
    }
  
  if(native==false) {
    status=apply_zone_arg(&grid,prescribed,dx,dy);
    if(status!=0){
      grid.modeH=0;
      remap=1;
      }
    }
  
  grid.brief_print();
  
  /* necessary because of buggy tides_savec1 */
  status=map_completegridaxis(&grid,2);
  
  if(remap) {
    printf("#################################################################\n");
    printf("remap tidal solution\n");
#define USE_index_interpolation 1
#if USE_index_interpolation
    set_grid_list(&atlas_grid, 1);
#endif
    for(k=0;k<nbuffers;k++) {
      fcomplex *tmp=new fcomplex[grid.nx*grid.ny];
      for (j=0;j<grid.ny;j++) {
#if USE_index_interpolation
        int64_t accel=-1;
#endif
        for (i=0;i<grid.nx;i++) {
          n=i+grid.nx*j;
          grid.xy(i,j,x,y);
#if USE_index_interpolation
          index_interpolation(atlas_grid,x,y,&accel, tide[k], cmask,&z);
#else
          x=map_recale(atlas_grid,x);
          status=map_interpolation(atlas_grid, tide[k], cmask, x,y,&z);
#endif
          tmp[n]=z;
          }
        }
      
      delete[] tide[k];
      tide[k]=tmp;
      }
    }
  
  if(smooth) {
    printf("#################################################################\n");
    printf("smooth tidal solution : %f km\n",lscale/1000.);
    for(k=0;k<nbuffers;k++) {
      grid_t subgrid;
      complex< float > *smoothed,submask;
      smoothed=new complex<float>[grid.ny*grid.nx];
      
//       status=map_remap(grid, tide[k], cmask, &subgrid, &subtides, &submask, 2, 1);
//       subsmoothed=new complex<float>[subgrid.ny*subgrid.nx];
//       status=map_smooth_latThenLong(subgrid, subtides,    submask, lscale, subsmoothed);
//       status=map_smooth_latThenLong(subgrid, subsmoothed, submask, lscale, subsmoothed);
//       status=map_smooth_latThenLong(subgrid, subsmoothed, submask, lscale, subsmoothed);
//       status=map_smooth_latThenLong(subgrid, subsmoothed, submask, lscale, subsmoothed);
//       status=map_export(subgrid, subsmoothed, submask, grid, smoothed, cmask, 0);
      
      for (j=0;j<grid.ny;j++) {
        for (i=0;i<grid.nx;i++) {
          n=i+grid.nx*j;
          if(tide[k][n]!=cmask) tide[k][n]=abs(tide[k][n]);
          }
        }
      grid.circular=map_check_circular(grid);
      status=map_smooth_latThenLong(grid, tide[k], cmask, lscale, smoothed);
      delete[] tide[k];
      tide[k]=smoothed;
      }
    }
  
  if(output ==NULL) {
    string s=input;
    int npos=s.rfind(".nc");
    string root;
    if(npos==s.size()-3)
      root=s.substr(0,npos)+"-reformatted.nc";
    else{
      npos=s.find_last_of("/.");
      if(npos<0 || s[npos]=='/')
        npos=s.length();
      root=s.substr(0,npos)+".nc";
      }
    output=strdup(root.c_str());
    }

//   spectrum.n=1;
//
//   exitIfNull(spectrum.waves=new tidal_wave[spectrum.n]);
//
//   for (k=0;k<spectrum.n;k++) {
//     strcpy(spectrum.waves[k].name,wave);
//     }
/*
  status=nc_writespectrum(ncid, spectrum);
  if(status !=0) goto error;
*/
{
  
  if(bmg != 0){
    string s="", root="";
    s.assign(input);
    size_t npos=s.rfind(".");
    root=s.substr(0,npos);
    sprintf(output,"%s.bimg",root.c_str());
    char *tmp=strdup(output);
    for (k=0;k<spectrum.n;k++) {
      if(spectrum.n!=1) {
        sprintf(output,"%s.%s",spectrum.waves[k].name,tmp);
        }
      if(scale!=1.0) {
        for(n=0;n<grid.Hsize();n++) if (tide[k][n]!=cmask) tide[k][n]*=scale;
        }
      bmg_savec1(output, 1, 1, 1, grid, tide[k], time, cmask);
      }
    }
  else {
    char *tmp=strdup(output);
    for (k=0;k<spectrum.n;k++) {
      if(spectrum.n!=1) {
        sprintf(output,"%s.%s",spectrum.waves[k].name,tmp);
        }
      if(format=="ascii" and nbuffers==2){
        poc_var_t av,gv;
        printf("#################################################################\n");
        printf("write output file : %s\n",output);
        status=poc_save_grid(output,&av,grid,1,1);
        gv=av.init("Ua",NC_FLOAT,"","m/s",real(cmask));
        gv   .init("Ug",NC_FLOAT,"","deg",real(cmask));
        status=poc_put_cvara(output,av,gv,0,tide[0]);
        gv=av.init("Va",NC_FLOAT,"","m/s",real(cmask));
        gv   .init("Vg",NC_FLOAT,"","deg",real(cmask));
        status=poc_put_cvara(output,av,gv,0,tide[1]);
        }
      else {
        tides_savec1(output, grid,  tide[k], cmask, scale, Aname, Gname);
        }
      status=poc_def_att(output,poc_att_t("production","constructed around " __LINE_FILE_PACKAGE_REVISION));
      status=poc_def_att(output,poc_att_t("history",cmd));
      }
    }
}
  for(k=0;k<nbuffers;k++) delete[] tide[k];
  delete[] tide;
  
  grid.free();
  
  STDOUT_BASE_LINE("end of tides-converter... \n");
  exit(0);
}
