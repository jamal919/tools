
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
#include "functions.h"
#include "grd.h"
#include "map.h"

#include "topo.h"
#include "map.h"

#define XYZ 0
#define YXZ 1

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
 
  add bathymetry replicant 
  
  add grid definition
  
*/

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
    "  -inv : toggle depth/altitude on final topography\n"
    "  -positive : toggle depth/altitude on imported topography\n"
    "  -debug : \n"
    "  -add : followed by a path of a file to convert hydrographic level to mean level by adding\n"
    "  -subtract,-substract : followed by a path of a file to convert hydrographic level to mean level by subtracting\n"
    "  --masked-only : only import on masked\n"
    "  -persistence : \n"
    "  -zmin : followed by minimum imported depth\n"
    "  -zmax : followed by maximum imported depth\n"
    "  -tag : \n"
    "  -iF : followed by anything.\n"
    "  --proj=... : libproj4 projection parameters. Will only have an effect if -iF is specified. Will crash on mode 0 and 1 grids.\n"
    "  -b : followed by working bathymetry\n"
    "  -f : followed by output format: grd, grd-float, nectdf or xyz. Default is grd\n"
    "  -l : followed by list of inputs. Overrides input\n"
    "  -p : followed by selection polygone\n"
    "  -o : followed by output BASE name\n"
    "  -z : \n"
    );
  print_OPENMP_help(prog_name); /** \endcode */
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int get_list(const char *filename, vector <char *> & plates)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int nitems, nplates;
  FILE *in;
  vector <size_t> row,col;
  vector <char *> vname;
  char *tmpP/*,*tmpV*/;
  
  in=fopen(filename,"r");
  if(in==0) return(errno);
  
  do {
    tmpP=new char[1024];
//    tmpV=new char[64];
    nitems=fscanf(in,"%s", tmpP);
    if(nitems!=1) break;
    plates.push_back(tmpP);
//    vname.push_back(tmpV);
    } while (nitems==1);
  
  nplates=plates.size();
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  float  bscale=1.0, scale=1.0;

  int i,j,k,m,n,signus=0,masked_only=0,status;

  const char *keyword;

  char *rootname=NULL,*output=NULL,*input=NULL,*bathymetry=NULL;
  char *format=NULL,*hydro=NULL,*poly=NULL,*list=NULL;
  char *input_format=NULL, *proj4=NULL;

  vector <char *> plates;
  int nplates;
  
  grid_t topogrid;
  grid_t grid,zone_grid;
  mesh_t mesh;
  short *topo,smask=256*127+255;
  float *ftopo,ftopomask;

  float zmin=0,zmax=0;

  pocgrd_t ncgrid;
  cdfvar_t variable;
  
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
    if(strcmp(keyword,"-inv")==0) {
      scale=-1;
      n++;
      continue;
      }
    if(strcmp(keyword,"-positive")==0) {
      bscale=-1;
      n++;
      continue;
      }
    if(strcmp(keyword,"-debug")==0) {
      debug=true;
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
    if(strcmp(keyword,"-substract")==0 or
       strcmp(keyword,"-subtract")==0){
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
    if(strcmp(keyword,"-persistence")==0) {
      sscanf(argv[n+1],"%d",&persistence);
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
    if(strcmp(keyword,"-tag")==0) {
      sscanf(argv[n+1],"%d",&tag);
      tagged=true;
      n++;
      n++;
      continue;
      }
    if(strcmp(keyword,"-iF")==0) {
      input_format=strdup(argv[n+1]);
      n++;
      n++;
      continue;
      }
    if(strstr(keyword,"--proj=")!=0){
       proj4=strdup(argv[n]+7);
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

          case 'l' :
            list= strdup(argv[n+1]);
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

          default:
            STDOUT_BASE_LINE("unknown option %s\n",keyword);
            exit(-1);
            break;
          }
        break;

      default:
        if(input==NULL) {
          input=strdup(argv[n]);
          n++;
          }
        else {
          STDOUT_BASE_LINE("unknown option %s\n",keyword);
          exit(-1);
          }
        break;
      }
    }

  if(input==0) {
    if(list!=0) {
      status=get_list(list, plates);
      if(status !=0) {
        printf("cannot get filenames from list=%s\n",list);
        goto error;
        }
      }
    }
  else{
    plates.push_back(input);
    }
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  load base bathymetry grid
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  printf("#####################################################################\n");
  printf("load working bathymetry: %s\n",bathymetry);
  status=topo_loadfield(bathymetry, &topogrid, &ftopo, &ftopomask, debug);
  if(status !=0) {
    printf("cannot load bathymetry file=%s\n",bathymetry);
    goto error;
    }
  map_printgrid(topogrid);

// if resize resolution
//     int map_remap (grid_t & g_source, float *b_source, float m_source, grid_t *g_target, float **b_target, float *m_target, int incr,int mode)

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  tags handling
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  size_t vlength;
  if(tagged) {
    status=poc_get_var_length( bathymetry, "tag",&vlength);
    if(status==0) {
      tagbuf=aset(topogrid.Hsize(),(short) tagmask);
      poc_get_var( bathymetry, "tag", tagbuf);
      status=grd_mirror_r( topogrid, topogrid.nx, tagbuf, tagmask);
      }
    else tagbuf=aset(topogrid.Hsize(),(short) -1);
    }
  else tagbuf=0;
  
/*------------------------------------------------------------------------------
  apply change of sign on base topography*/
  for (j=0;j<topogrid.ny;j++) {
    for (i=0;i<topogrid.nx;i++) {
      m=j*topogrid.nx+i;
      if(ftopo[m]!=ftopomask) {
        ftopo[m]*=bscale;
        }
      }
    }
 
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  import new depths
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  nplates=plates.size();
  if(nplates!=0) {
    for(k=0;k<nplates;k++) {
      printf("\n#####################################################################\n");
      printf("import extra bathymetry for merging: %s\n", plates[k]);
//      status=map_import(input, topogrid, ftopo, ftopomask,  zmin,  zmax, poly, masked_only);
      range_t<float> base_range;
      range_t<float> inport_range;
      status=topo_import(plates[k], 0, topogrid, ftopo, ftopomask,  zmin,  zmax, poly, masked_only, tagbuf, tag, input_format, proj4, debug);
      if(status !=0) {
        goto error;
        }
      }
    }
    
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  convert "zero-hydro related" depths to "mean level related" depths
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if(hydro!=NULL) {
    status= topo_hydro(hydro, signus, topogrid, ftopo, ftopomask,5);
    if(status !=0) {
      goto error;
      }
    }
  if(persistence!=0) {
    for(int k=0;k<persistence;k++) status=map_persistence(topogrid, ftopo, ftopomask, 0);
    }
    
  if(debug) status=topo_checks(rootname, topogrid, ftopo, ftopomask);

/*------------------------------------------------------------------------------
  apply change of sign on final topography*/
  if(scale!=1.0) {
    for (size_t j=0;j<topogrid.ny;j++) {
      for (size_t i=0;i<topogrid.nx;i++) {
        size_t m=j*topogrid.nx+i;
        if(ftopo[m]!=ftopomask) {
          ftopo[m]*=scale;
          }
        }
      }
    }
 
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  save rectified topo 
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if(rootname==0) rootname=strdup("out.grd");

  if(format==NULL){
    if(strrncasecmp(rootname,".grd")==0){
      format=strdup("grd");
      }
    else if(strrncasecmp(rootname,".nc")==0){
      format=strdup("netcdf");
      }
    else{
      format=strdup("grd");
      }
    }


/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  save data, to be replaced with topo_save. Example:
  
  status=topo_save(rootname, output, "netcdf", TOPOgrid, TOPOdata, TOPOmask, debug);
  
  topo_save will hande mirroring...

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if(strcmp(format,"grd")==0) {
    topo=new short[topogrid.nx*topogrid.ny];
    smask=256*127+255;
    for (size_t j=0;j<topogrid.ny;j++) {
      for (size_t i=0;i<topogrid.nx;i++) {
        size_t m=j*topogrid.nx+i;
        size_t n=(topogrid.ny-j-1)*topogrid.nx+i;
        if(ftopo[m]!=ftopomask) {
          topo[n]=(short) floor(scale*ftopo[m]+0.5);
          }
        else {
          topo[n]=smask;
          }
        }
      }
    if(strrncasecmp(rootname,".grd")==0) asprintf(&output,"%s",rootname);
    else asprintf(&output,"%s.grd",rootname);
    printf("#################################################################\n");
    printf("create depth file (grd short format) : %s\n",output);
    status=grd_save(output,topogrid,topogrid.nx, topo,smask);
    delete[] topo;
    }
  else if(strcmp(format,"grd-float")==0) {
    if(strrncasecmp(rootname,".grd")==0) asprintf(&output,"%s",rootname);
    else asprintf(&output,"%s.grd",rootname);
    printf("#################################################################\n");
    printf("create depth file (grd float format) : %s\n",output);
    status=grd_mirror_r( topogrid, topogrid.nx, ftopo, ftopomask);
    status=grd_save(output,topogrid,topogrid.nx, ftopo, ftopomask);
    }
  else if(strcmp(format,"netcdf")==0) {
    if(topogrid.modeH==0) status= map_completegridaxis(&topogrid,1);
    if(strrncasecmp(rootname,".nc")==0) asprintf(&output,"%s",rootname);
    else asprintf(&output,"%s.nc",rootname);
    printf("#################################################################\n");
    printf("create depth file (netcdf float format) : %s\n",output);
    status= poc_createfile(output);
    status=poc_sphericalgrid_xy(output,"",topogrid,&ncgrid);
    poc_standardvariable_xy(&variable,"bathymetry",ftopomask,"m",1., 0.,"bathymetry","bathymetry","bathymetry",ncgrid);
    status=create_ncvariable(output, &variable);
    status=poc_write_xy(output,  topogrid, variable.id,ftopo);
    variable.destroy();
    }
  else {
    TRAP_ERR_EXIT(-1,"unknown format for output\n");
    }
    
  if(tagged) {
    poc_var_t var;
    status=poc_inq_var( output, "z", &var, 0);
    int id;
    var.init("tag",NC_SHORT,"tag","none",tagmask);
    status=occurence<short>((short) tag, tagbuf, topogrid.Hsize());
    status=poc_def_var( output, var, &id, 0);
    status=grd_mirror_r( topogrid, topogrid.nx, tagbuf, tagmask);
    status=poc_put_var(output, id, tagbuf);
    }

  delete[] ftopo;
  topogrid.free();

  STDOUT_BASE_LINE("topo_merge sucessfully completed\n");

  exit(0);

 error:
  STDOUT_BASE_LINE("topo_merge aborted\n");
  exit(-1);
}
