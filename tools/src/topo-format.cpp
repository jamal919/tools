
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
#include "ascii.h"
#include "map.h"
#include "map.def"

#include "topo.h"

#define XYZ 0
#define YXZ 1



/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  float  scale=1.0,topo_mask=1.e+10,mask=1.e+10;
  double  x,y,t,p;
  double  dx,dy;

  int i,j,k,l,m,n,signus=0,status;

  size_t size;
  FILE *in,*out;
  char *keyword,*zone;
  int flag;

  char *rootname=NULL,*output=NULL,*input=NULL,*pathname=NULL,*bathymetry=NULL;
  char *bathy=NULL,*format=NULL,*hydro=NULL,*poly=NULL;
  char *varname=0;
  char *xname=0,*yname=0;
  char *iF=0;
  char *proj4=NULL;

  grid_t topogrid;
  grid_t grid;
  mesh_t mesh;
  short *topo,smask=256*127+255;
  float *ftopo,*tmp,*buffer, ftopomask;
  bool force_mask=false, flip=false;
  
  float zmin=0,zmax=0;

  plg_t *polygones=NULL;
  int npolygones=0;

  pocgrd_t ncgrid;
  cdfvar_t variable;
  int dimlgth[4]={10,20,30,40};

  fct_echo( argc, argv);

  n=1;
  while (n < argc) {
    keyword=strdup(argv[n]);
    if(strcmp(keyword,"-iF")==0) {
      iF= strdup(argv[n+1]);
      n++;
      n++;
      continue;
      }
    if(strcmp(keyword,"-oF")==0) {
      format= strdup(argv[n+1]);
      n++;
      n++;
      continue;
      }
    if(strcmp(keyword,"-inv")==0) {
      scale=-1;
      n++;
      continue;
      }
    if(strcmp(keyword,"-flip")==0) {
      flip=true;
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

          case 'v' :
            varname= strdup(argv[n+1]);
            n++;
            n++;
            break;

          case 'x' :
            xname= strdup(argv[n+1]);
            n++;
            n++;
            break;

          case 'y' :
            yname= strdup(argv[n+1]);
            n++;
            n++;
            break;

          case '-' :
            if(strstr(argv[n],"--proj=")!=0){
              proj4=strdup(argv[n]+7);
              n++;
              }
            break;

          default:
            if(strcmp(keyword,"-zmin")==0) {
              sscanf(argv[n+1],"%f",&zmin);
              n++;
              n++;
              }
            if(strcmp(keyword,"-scale")==0) {
              sscanf(argv[n+1],"%f",&scale);
              n++;
              n++;
              }
            if(strcmp(keyword,"-zmax")==0) {
              sscanf(argv[n+1],"%f",&zmax);
              n++;
              n++;
              }
            if(strcmp(keyword,"-mask")==0) {
              sscanf(argv[n+1],"%f",&mask);
	          force_mask=true;
              n++;
              n++;
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

/*------------------------------------------------------------------------------
  load bathymetry grid*/
  int fmt;
  
  if(iF!=0) fmt=map_get_format(iF);
  else fmt=-1;
  
/*
  EMODNET
  "PROJCS[\"Mercator\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.000000 , 298.257223563 ]],PRIMEM[\"Greenwich\",0.0],
  UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Mercator\"],PARAMETER[\"False_Easting\",0.000000],PARAMETER[\"False_Northing\",0.000000],PARAMETER[\"Central_Meridian\",0.000000],
  PARAMETER[\"Standard_Parallel_1\",46.000000],UNIT[\"Meter\",1.0]]"*/  
  
//  proj4=strdup("proj=merc ellps=WGS84 lon_0=0E lat_0=49N lat_ts=46N");

/*
 NZ
 -iF netcdf -f netcdf w001001.nc  -o toto -v Band1 "--proj=merc ellps=WGS84 lon_0=100E lat_ts=-41N"
 Projection Mercator (WGS84 Datum) at 41°S true latitude and 100°E central meridian (EPSG 3994).
*/
   status=map_loadfield(input, fmt, varname, xname, yname, proj4, &topogrid, &ftopo, &ftopomask);
//  status=map_loadfield(input, MAP_FORMAT_CODE_GLOBE, (char *) 0, &topogrid, &ftopo, &ftopomask);
  if(status !=0) {
    printf("cannot load bathymetry file=%s\n",input);
    goto error;
    }

  if(force_mask) ftopomask=mask;

/*------------------------------------------------------------------------------
  convert "zero-hydro related" depths to "mean level related" depths*/
  if(hydro!=NULL) {
    status= topo_hydro(hydro, signus, topogrid, ftopo, ftopomask,0);
    }

/*------------------------------------------------------------------------------
  save rectified topo */
  if(rootname==0) rootname=strdup("out");
  output=new char[1024];

  if(format==NULL) format=strdup("grd");
  if(scale!=1.0) {
    for (j=0;j<topogrid.ny;j++) {
      for (i=0;i<topogrid.nx;i++) {
        m=j*topogrid.nx+i;
        if(ftopo[m]!=ftopomask) {
          ftopo[m]=scale*ftopo[m];
          }
        }
      }
    }

  if(flip) status=grd_mirror_r( topogrid, topogrid.nx, ftopo, ftopomask);

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  save data, to be replaced with topo_save. Example:
  
  status=topo_save(rootname, output, "netcdf", TOPOgrid, TOPOdata, TOPOmask, debug);
  
  topo_save will hande mirroring...

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

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
    sprintf(output,"%s.grd",rootname);
    printf("#################################################################\n");
    printf("create depth file (grd short format) : %s\n",output);
    status=grd_save(output,topogrid,topogrid.nx, topo,smask);
    delete[] topo;
    }
  if(strcmp(format,"grd-float")==0) {
    sprintf(output,"%s.grd",rootname);
    printf("#################################################################\n");
    printf("create depth file (grd float format) : %s\n",output);
    status=grd_mirror_r( topogrid, topogrid.nx, ftopo, ftopomask);
    status=grd_save((const char*) output,topogrid,topogrid.nx, ftopo, ftopomask);
    }
  else if(strcmp(format,"netcdf")==0) {
    if(topogrid.modeH==0) status=map_completegridaxis(&topogrid,1); 
    sprintf(output,"%s.nc",rootname);
    printf("#################################################################\n");
    printf("create depth file (netcdf format) : %s\n",output);
    int cmode=NC_CLOBBER | NC_NETCDF4| NC_CLASSIC_MODEL;
    status= poc_createfile(output,cmode);
    status=poc_sphericalgrid_xy(output,"",topogrid,&ncgrid);
    poc_standardvariable_xy(&variable,"bathymetry",ftopomask,"m",1., 0.,"bathymetry","bathymetry","bathymetry",ncgrid);
    status=create_ncvariable(output, &variable);
    status=poc_write_xy(output,  topogrid, variable.id,ftopo);
    variable.destroy();
    }
  else if(strcmp(format,"xyz")==0) {
    sprintf(output,"%s.xyz",rootname);
    printf("#################################################################\n");
    printf("create depth file (xyz format) : %s\n",output);
    status=xyz_save((const char *) output,topogrid, topo,smask);
    }
  else if(strcmp(format,"gis")==0) {
    sprintf(output,"%s.gis",rootname);
    printf("#################################################################\n");
    printf("create depth file (gis format) : %s\n",output);
    status=ascii_save_GIS((const char *) output, topogrid, ftopo, ftopomask, " %7.1f", topogrid.nx);
    }
  delete[] ftopo;

  __OUT_BASE_LINE__("bathymetry sucessfully completed\n");

  exit(0);

 error:
  __ERR_BASE_LINE__("exiting\n");exit(-1);
}
