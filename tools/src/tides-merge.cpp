
/*******************************************************************************

  T-UGO tools, 2006-2009

  Unstructured Ocean Grid initiative

Contributors:
20
  Florent Lyard      LEGOS/CNRS, Toulouse, France
  Yoann Le Bars      LEGOS, Toulouse, France (PhD)
  Laurent Roblou     LEGOS, Toulouse, France
  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
  Frédéric Dupont    Université de Laval à Québec, Canada

E-mail: florent.lyard@legos.obs-mip.fr

*******************************************************************************/

/* *@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  topo-merg aims at merging 2 bathymetric database with some controls
  on merging action

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */

#define MAIN_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include "tools-structures.h"
#include "version-macros.def" //for VERSION and REVISION

#include "rutin.h"
#include "geo.h"
#include "polygones.h"
#include "poc-netcdf.hpp"
#include "poc-netcdf-data.hpp"
#include "functions.h"
#include "grd.h"
#include "map.h"
#include "filter.h"

#include "tides.h"
#include "topo.h"
#include "map.h"

#define XYZ 0
#define YXZ 1


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
    "DESCRIPTION\n"
/* USE FIND AND REPLACE TO EDIT THE LINE BELOW, SO THAT THE FILE HEAD IS ALSO UPDATED */
    "  merge together MNTs.\n"
    "\n"
    "OPTIONS :\n"
    "  -h,--help : show this help and exit\n"
    "  -b import : merge import MNT in input MNT\n"
    "  -debug : \n"
    "  --masked-only : only import on masked\n"
    "  -persistence : \n"
    "  -tag : \n"
    "  --proj=... : libproj4 projection parameters. Will only have an effect if -iF is specified. Will crash on mode 0 and 1 grids.\n"
    "  -p : followed by selection polygone\n"
    "  -o : followed by output BASE name\n"
    "  -z : \n"
    );
  print_OPENMP_help(prog_name); /** \endcode */
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int merge_template(const grid_t grids[2], complex<double> *buffers[2], complex<double> masks[2], vector<plg_t> polygons, double* & weights)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  merge grid[1]'s buffer into grid[0]'s buffer

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
{
  int status;
  int i,j;
  int initialise=0;
  bool spherical=true;
  size_t m, count;
  double *smoothed=0, dmask=1.e+10;
  signed char *selected;
  complex<double> *tmp, *backup;
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  export buffers[0] on grid[1]

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  tmp=new complex<double> [grids[1].Hsize()];
  status=map_export(grids[0],buffers[0],masks[0],grids[1],tmp,masks[1],0);
    
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  mask polygon-wise inside cells

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  selected=new signed char[grids[1].Hsize()];
  
  for(m=0;m<grids[1].Hsize();m++) selected[m]=1;
    
  status=topo_PolygonsSelection(grids[1], polygons, selected, spherical, initialise);

  for(m=0;m<grids[1].Hsize();m++) {
    if(buffers[1][m]==masks[1] and selected[m]==1) tmp[m]=masks[1];
    }
    
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  compute weights

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if(weights==0) {
    weights=new double[grids[1].Hsize()];
  
    for(m=0;m<grids[1].Hsize();m++) {
      if(selected[m]==0) weights[m]= 0;
      else               weights[m]=+1;
      }
  
    for(int k=0;k<25;k++) status=map_persistence(grids[1], weights, 0.0, 0);
  
    smoothed=new double[grids[1].Hsize()];
    status=map_smooth(LOESS, grids[1], weights, dmask, (float) 250000.0, smoothed);
  
    for(m=0;m<grids[1].Hsize();m++) weights[m]=smoothed[m];
  
    delete[] smoothed;
    }
    
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  merge buffers

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  for(m=0;m<grids[1].Hsize();m++) {
    if(selected[m]==1) {
      continue;
      }
/*------------------------------------------------------------------------------
    if both atlas masked, set mask */
    if(buffers[1][m]==masks[1] and tmp[m]==masks[1]) continue;
    if(buffers[1][m]==masks[1]) {
      buffers[1][m]=tmp[m];
      continue;
      }
    if(tmp[m]==masks[1]) {
      continue;
      }
    if(selected[m]==0) {
      buffers[1][m]=weights[m]*buffers[1][m]+(1.-weights[m])*tmp[m];
      }
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  update buffers[0]

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  backup=new complex<double> [grids[0].Hsize()];
  for(m=0;m<grids[0].Hsize();m++) backup[m]=buffers[0][m];
  
  status=map_export(grids[1],buffers[1],masks[1],grids[0],buffers[0],masks[0],0);
  
  for(m=0;m<grids[0].Hsize();m++) {
    if(buffers[0][m]==masks[0]) buffers[0][m]=backup[m];
    }
    
  delete[] tmp;
  delete[] selected;

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int n,status;

  const char *keyword;

  string rootname, output;
  string PolygonsFile, inputs[2];
  vector<string> waves;
  char *proj4=NULL;
  string varnames[2][2];
  grid_t grids[2];
  complex<double> *buffers[2], masks[2];
  vector<plg_t> polygons;
  double *weights=0;
  poc_var_t av, gv, grid_var;
  
  int persistence=0;
  int tag=-1;
  
  bool debug=false;
  bool tagged=false;
  short *tagbuf=0,tagmask=-1;

  fct_echo( argc, argv);

  n=1;
  while (n < argc) {
    keyword=argv[n];
    if(strcmp(keyword,"--help")==0 or
       strcmp(keyword,"-h")==0) {
      print_help(argv[0]);
      exit(0);
      }
    if(strcmp(keyword,"-debug")==0) {
      debug=true;
      n++;
      continue;
      }
//     if(strcmp(keyword,"--masked-only")==0) {
//       masked_only=1;
//       n++;
//       continue;
//       }
    if(strcmp(keyword,"-persistence")==0) {
      sscanf(argv[n+1],"%d",&persistence);
      n++;
      n++;
      continue;
      }
//     if(strcmp(keyword,"-tag")==0) {
//       sscanf(argv[n+1],"%d",&tag);
//       tagged=true;
//       n++;
//       n++;
//       continue;
//       }
    if(strstr(keyword,"--proj=")!=0){
       proj4=strdup(argv[n]+7);
       n++;
       continue;
       }
    switch (keyword[0]) {
      case '-':
        switch (keyword[1]) {

          case 'p' :
            PolygonsFile=argv[n+1];
            n++;
            n++;
            break;

          case 'o' :
            output= strdup(argv[n+1]);
            n++;
            n++;
            break;

          case 'r' :
            rootname= strdup(argv[n+1]);
            n++;
            n++;
            break;

          default:
            STDOUT_BASE_LINE("unknown option %s\n",keyword);
            exit(-1);
            break;
          }
        break;

      default:
//         if(input==NULL) {
//           input=strdup(argv[n]);
//           n++;
//           }
//         else {
//           STDOUT_BASE_LINE("unknown option %s\n",keyword);
//           exit(-1);
//           }
        break;
      }
    }

//   if(input==0) {
//     }
  
  waves.push_back("M2");
  waves.push_back("S2");
  waves.push_back("N2");
  waves.push_back("K2");
  
  waves.push_back("K1");
  waves.push_back("O1");
  waves.push_back("Q1");
  waves.push_back("P1");

  for(int w=0;w<waves.size();w++) {
  
  string atlas_directory,atlas_convention;
  char *tmp;
  
//   atlas_directory="/home/data/tides/FES2014/structured";
//   atlas_convention="WAVE.FES2014.nc";
//   
//   varnames[0][0]="Ha";
//   varnames[0][1]="Hg";

  atlas_directory=".";
  atlas_convention="WAVE.FES2014b-NWP.nc";
  
  varnames[0][0]="elevation_a";
  varnames[0][1]="elevation_G";
  
  status=tide_decode_atlasname(atlas_directory.c_str(),atlas_convention.c_str(),waves[w].c_str(), 0, &tmp);
  inputs[0]=(string) tmp;
  delete[] tmp;

  atlas_directory="/home/models/sirocco/dev-182";
  atlas_convention="WAVE.generic.nc";
  
  varnames[1][0]="elevation_a";
  varnames[1][1]="elevation_G";

  status=tide_decode_atlasname(atlas_directory.c_str(),atlas_convention.c_str(),waves[w].c_str(), 0, &tmp);
  inputs[1]=(string) tmp;
  delete[] tmp;

/*------------------------------------------------------------------------------
  load inputs */
  for(int k=0;k<2;k++) {
    printf("#####################################################################\n");
    printf("load working file: %s\n",inputs[k].c_str());
    poc_var_t var;
    status=poc_get_grid(inputs[k], varnames[k][0], &(grids[k]), 0);
  
    buffers[k]=new complex<double>[grids[k].Hsize()];
  
    status=poc_get_cvara(inputs[k], varnames[k][0], varnames[k][1],-1, buffers[k]);
    masks[k]=NC_FILL_COMPLEX;
    
//     poc_data_t<double> scalarData;
//     status=scalarData.init(filename,varname);
//     status=scalarData.read_data(filename,i);
//     mask=scalarData.mask;
 
    if(status !=0) {
      printf("cannot load input file=%s\n",inputs[k].c_str());
      goto error;
      }
    }
    
  status=plg_load(PolygonsFile, PLG_FORMAT_UNKNOWN, polygons);

  status=merge_template(grids, buffers, masks, polygons,weights);
    
/*------------------------------------------------------------------------------
  save */
  if(rootname=="") rootname="out";
  
  output=waves[w]+"."+rootname+"-1.nc";
  
  status=poc_save_grid(output,&grid_var,grids[0],1,1);
  av=gv=grid_var;
  av.init(varnames[0][0],NC_FLOAT,comodo_standard_name("",waves[w]),"m",real(masks[0]));
  gv.init(varnames[0][1],NC_FLOAT,comodo_standard_name("",waves[w]),"degrees",real(masks[0]));
  
  status=poc_put_cvara(output,av,gv,0,buffers[0]);
  
  output=waves[w]+"."+rootname+"-2.nc";
  
  status=poc_save_grid(output,&grid_var,grids[1],1,1);
  av=gv=grid_var;
  av.init(varnames[1][0],NC_FLOAT,comodo_standard_name("",waves[w]),"m",real(masks[1]));
  gv.init(varnames[1][1],NC_FLOAT,comodo_standard_name("",waves[w]),"degrees",real(masks[1]));
  
  status=poc_put_cvara(output,av,gv,0,buffers[1]);
  
  av=grid_var;
  av.init("weights",NC_FLOAT,"ponderating weights","none",real(masks[1]));
  
  status=poc_put_vara(output,av,0,weights);
  
  }
    
  STDOUT_BASE_LINE("tides-merge sucessfully completed\n");

  exit(0);

 error:
  STDOUT_BASE_LINE("tides-merge aborted\n");
  exit(-1);
}
