
/**************************************************************************

  T-UGO tools, 2006-2013

  Unstructured Ocean Grid initiative

***************************************************************************/
/** \file

\author  Florent Lyard      LEGOS/CNRS, Toulouse, France. florent.lyard@legos.obs-mip.fr
\author  Laurent Roblou     LEGOS/CNRS, Toulouse, France
\author  Damien Allain      LEGOS/CNRS, Toulouse, France
\author  Yves Soufflet      LEGOS/CNRS, Toulouse, France
\author  Clément Mayet      LEGOS, Toulouse, France (PhD)
\author  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
\author  Frédéric Dupont    Université de Laval à Québec, Canada

<!-- USE FIND AND REPLACE TO EDIT THIS LINE, SO THAT print_help IS ALSO UPDATED -->
\brief Assembles tiles of data

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

#include "poc-netcdf-data.hpp"
#include "grd.h"
#include "map.h"
#include "functions.h"
#include "topo.h"

string cmd;


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void print_topo_assembly_help()

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// print help for topo_assembly() input format
/*----------------------------------------------------------------------------*/
{
/** The code of the body of this function is :
    \code /**/ // COMPILED CODE BELOW !!!
  printf(
    "TILING DESCRIPTION FORMAT\n"
    "  The first line has the X (#columns) and Y (#rows) size of the output, in number of tiles.\n"
    "The following lines are made of :\n"
    "- the X (column index) and Y (row index) indexes of the position of the tile (indices origin is bottom/left),\n"
    "- the file name,\n"
    "- the variable name of the data.\n"
    "\n"
    "  The content of a real life example follows:\n"
    "3 2\n"
    "0 0 n03_e008_3arc_v1.nc Band1\n"
    "1 0 n03_e009_3arc_v1.nc Band1\n"
    "2 0 n03_e010_3arc_v1.nc Band1\n"
    "0 1 n04_e008_3arc_v1.nc Band1\n"
    "1 1 n04_e009_3arc_v1.nc Band1\n"
    "2 1 n04_e010_3arc_v1.nc Band1\n"
    "\n");
  /** \endcode */
}
  

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int topo_assembly(const char *filename,const char *rootname, float mask, bool shrinked, bool grd, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  frame_t frame;
  int status;
  int i,j,k,nitems;
  size_t m,n;
  int nrows, ncols, nplates;
  range_t<double> lon_range,lat_range;
  grid_t SRTMgrid,  grid[128];
//   short *SRTMdata, *data[128], SRTMmask;
  float *SRTMdata, *data[128], SRTMmask;
  grid_t BASEgrid;
  FILE *in;
  vector <size_t> row,col;
  vector <char *> plate,vname;
  char *tmpP,*tmpV;
  string path="";
  int verbose=(debug==true);
  
  in=fopen(filename,"r");
  
  if(in==0) return(-1);
  
  nitems=fscanf(in,"%d %d", &ncols, &nrows);
  if(nitems!=2) return(-1);
  
  k=0;
  do {
    tmpP=new char[1024];
    tmpV=new char[64];
    nitems=fscanf(in,"%d %d %s %s", &i, &j, tmpP, tmpV);
    if(nitems!=4) break;
    row.push_back(j);
    col.push_back(i);
    plate.push_back(tmpP);
    vname.push_back(tmpV);
    } while (nitems==4);
  
  nplates=row.size();
  
  printf("%d tiles found, ncols=%d nrows=%d\n",nplates,ncols,nrows);
  
  range_t<size_t> r=range_t<size_t>(row);
  range_t<size_t> c=range_t<size_t>(col);

  for(k=0;k<nplates;k++) data[k]=0;
  
  
  k=0;
  status=topo_loadfield((const char*) plate[k], vname[k], &grid[k], &data[k], &SRTMmask, debug);
  map_printgrid(grid[k]);
  
/*------------------------------------------------------------------------------
  load gridded elevation */
  for(k=1;k<nplates;k++) {
//     status=topo_loadfield_cdf((const char*) plate[k], vname[k], &grid[k], &data[k], &SRTMmask, debug);
    status=topo_loadfield((const char*) plate[k], vname[k], &grid[k], &data[k], &SRTMmask, debug);
    map_printgrid(grid[k]);
    if(status!=0) {
      printf("loading of %s in %s failed, abort\n",vname[k], plate[k]);
      return(-1);
      }
    }
    
  frame=limits(grid, nplates);
  
  if(grid[0].modeH==0) {
    }
  else if(grid[0].modeH==1) { // can be unaccurate
    grid[0].dx=grid[0].x[1]-grid[0].x[0];
    grid[0].dy=grid[0].y[1]-grid[0].y[0];
    }
  else if(grid[0].modeH==2) { // can be unstable
    grid[0].dx=grid[0].x[1]-grid[0].x[0];
    grid[0].dy=grid[0].y[grid[0].nx]-grid[0].y[0];
    }
  
  SRTMgrid=map_Rgrid(frame, grid[0].dx,grid[0].dy,0);
  status=map_completegridaxis(&SRTMgrid, 1);
  
  SRTMdata=new float[SRTMgrid.Hsize()];
  for(k=0;k<SRTMgrid.Hsize();k++) SRTMdata[k]=SRTMmask;
  
  for(k=0;k<nplates;k++) {
    size_t i0, j0;
    if(shrinked) {
      i0=nint((grid[k].xmin-SRTMgrid.xmin)/SRTMgrid.dx);
      j0=nint((grid[k].ymin-SRTMgrid.ymin)/SRTMgrid.dy);
      }
    else {
      i0=col[k]*(grid[k].nx-1);
      j0=row[k]*(grid[k].ny-1);
      }
    for(j=0;j<grid[k].ny;j++) {
      for(i=0;i<grid[k].nx;i++) {
        n=i+grid[k].nx*j;
/*------------------------------------------------------------------------------
           ok for disjoint tiles */
//         m=i+col[k]*grid[k].nx+SRTMgrid.nx*(j+row[k]*grid[k].ny);
/*------------------------------------------------------------------------------
        ok for overlapping tiles */
//         size_t ii=i+col[k]*(grid[k].nx-1);
//         size_t jj=j+row[k]*(grid[k].ny-1);
        size_t ii=i+i0;
        size_t jj=j+j0;
        m=ii+SRTMgrid.nx*jj;
        SRTMdata[m]=data[k][n];
        }
      }
    }
 
  if(mask!=INFINITY) {
    for(k=0;k<SRTMgrid.Hsize();k++) if(SRTMdata[k]==mask) SRTMdata[k]=SRTMmask;
    }
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  save data with netCDF format
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  path=(string) rootname+".nc";
  poc_var_t gridVar;//<template variable
  poc_var_t dataVar;
  status=poc_save_grid(path,&gridVar,SRTMgrid,1,verbose);
  status=poc_def_att(path,poc_att_t("production","constructed around " __LINE_FILE_PACKAGE_REVISION));
  status=poc_def_att(path,poc_att_t("history",cmd));
  dataVar=gridVar;
//   dataVar.init("Band1",NC_SHORT);
//   dataVar.init("Band1",NC_FLOAT);
  dataVar.init(vname[0],NC_FLOAT);
  dataVar << poc_att_t("content","YX");
  dataVar << poc_att_t("coordinates","latitude longitude");
  dataVar << poc_att_t("_FillValue",SRTMmask);
  status=poc_put_vara(path,dataVar,0,SRTMdata,1);
  

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  save data with historical grd format
  
  to be replaced with topo_save. Example:
  
  status=topo_save(rootname, output, "netcdf", TOPOgrid, TOPOdata, TOPOmask, debug);
  
  topo_save will hande mirroring...

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if(grd) {
    path=(string) rootname+".grd";
    status=grd_mirror_r( SRTMgrid, SRTMgrid.nx, SRTMdata, SRTMmask);
    status=grd_save(path.c_str(), SRTMgrid, SRTMgrid.nx, SRTMdata, SRTMmask);
    }
    
  return(status);
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
    "  %s version " VERSION " " REVISION "\n", prog_name);
  printf("\n"
    "USE\n"
    "  %s [OPTIONS] input\n",prog_name);
  printf("\n"
/* USE FIND AND REPLACE TO EDIT THE LINE BELOW, SO THAT THE FILE HEAD IS ALSO UPDATED */
    "DESCRIPTION\n"
    "  Assembles tiles of data.\n"
    "  This is particulary usefull for tiled data like SRTM.\n"
    "  It created a .nc file and a .grd file from the tiling description contained in the input file.\n"
    "\n"
    "OPTIONS\n"
    "  -h,--help  Show this help and exit.\n"
    "  -o  followed by the root name of the files created. Default : assembly\n"
    "\n");
  print_topo_assembly_help();
  /** \endcode */
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int n,status;

  char *keyword;

  char *rootname=NULL,*input=NULL;
  bool debug=false, grd=false, shrinked=false;
  float mask=+INFINITY;
  
  cmd=fct_echo( argc, argv);

  n=1;
  while (n < argc) {
    keyword=strdup(argv[n]);
    if(strcmp(keyword,"--help")==0 ||
      strcmp(keyword,"-h")==0){
      print_help(argv[0]);
      exit(0);
      }
    if(strcmp(keyword,"-debug")==0) {
      debug=true;
      n++;
      continue;
      }
    if(strcmp(keyword,"-mask")==0) {
      sscanf(argv[n+1],"%f",&mask);
      n++;
      n++;
      continue;
      }
    if(strcmp(keyword,"-grd")==0) {
      grd=true;
      n++;
      continue;
      }
    if(strcmp(keyword,"-shrinked")==0) {
      shrinked=true;
      n++;
      continue;
      }
    switch (keyword[0]) {
      case '-':
        switch (keyword[1]) {
        
        case 'o' :
          rootname= strdup(argv[n+1]);
          n++;
          n++;
          break;
        
        default:
          __OUT_BASE_LINE__("unknown option %s\n",keyword);
          print_help(argv[0]);
          exit(-1);
        }
        break;

      default:
        if(input==NULL) {
          input=strdup(argv[n]);
          n++;
          }
        else {
          __OUT_BASE_LINE__("unknown option %s\n",keyword);
          print_help(argv[0]);
          exit(-1);
          }
        break;
      }
      free(keyword);
    }
  
  if(input==NULL) {
    __ERR_BASE_LINE__("*** Tiling description needed ***\n");
    print_help(argv[0]);
    exit(-1);
    }
  
  if(rootname==NULL) rootname=strdup("assembly");
  
  status=topo_assembly(input, rootname, mask, shrinked, grd, debug);
  
  if(status!=0){
    __ERR_BASE_LINE__("topo-assembly(\"%s\",\"%s\") error\n",input,rootname);
    print_topo_assembly_help();
    }
  
  free(rootname);
  free(input);
  
  exit(status);
}
