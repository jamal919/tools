
/***************************************************************************
T-UGO tools, 2006-2009

  Unstructured Ocean Grid initiative
  
E-mail: florent.lyard@legos.obs-mip.fr
***************************************************************************/
/**
\file

\author  Florent Lyard      LEGOS/CNRS, Toulouse, France
\author  Yoann Le Bars      LEGOS, Toulouse, France (PhD)
\author  Laurent Roblou     LEGOS/CNRS, Toulouse, France
\author  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
\author  Frédéric Dupont    Université de Laval à Québec, Canada
\author  Clement MAYET      LEGOS, Toulouse, France (PhD)
\author  Yves Soufflet      LEGOS, Toulouse, France

VERSION :
\date  11/04/2011 : Clement MAYET: add SONEL_HR input format, add a print_help function and some Doxygen comments.
\date  24/04/2011 : Clement MAYET: Bugfix in arguments reading
\date  23/06/2011 : Yves Soufflet: Add input format for profilers

<!-- USE FIND AND REPLACE TO EDIT THIS LINE, SO THAT print_help IS ALSO UPDATED -->
\brief This program is a converter between different sealevel time serie formats

\todo 12/04/2011 : Clement MAYET: re-organize the functions (group the read functions, the write functions, the others... maybe in different files ) ?

***************************************************************************/

#define MAIN_SOURCE

#include "version-macros.def" //for VERSION and REVISION

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <fstream>
#include <algorithm>
#include <iostream>
#include <limits>
#include <string>
#include <map>

#include "tools-structures.h"

#include "tides.h"
#include "fe.h"
#include "archive.h"
#include "netcdf-proto.h"
#include "poc-time.h"
#include "mgr.h"
#include "mgr-converter.h"
#include "functions.h"
#include "map.h"
#include "grd.h"
#include "geo.h"
#include "tides.h"
#include "tides.h"

using namespace std;

extern int mgr_savehawaii(const char *,double,double,double *,double *,int,char *,date_t,int,double);


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
    " NAME AND VERSION\n"
    "  %s version " VERSION " " REVISION "\n", prog_name);
  printf("\n USE");
  printf("\n   %s FILE  [OPTIONS --header=\"lon=lon lat=lat code=code name=name depth=depth\"] \n",prog_name);
  printf("\n DESCRIPTION");
/* USE FIND AND REPLACE TO EDIT THE LINE BELOW, SO THAT THE FILE HEAD IS ALSO UPDATED */
  printf("\n   This program is a converter between different sealevel time serie formats");
  printf("\n   The input file must be in ascii format, and can be of type  SONEL_HR, RMN, PAB (see below for a description of these formats) \n\n");
  printf(
   "\n OPTIONS :\n"
   "   -f        specify the final date (if not specified the whole file is read)\n"
   "   -s        specify the start date (if not specified the whole file is read)\n"
   "   -o        specify the output filename (without extension!)\n"
   "   --header  contains information about the mgr station that can be printed in the output file\n"
    +isFormat_help("   ")+
    "               If no output format is specified, writes a list format with no header\n");
  printf(
   "\n FILE FORMAT : \n"
   "\n    -SONEL_HR       'YYYY-MM-DD HH:MM:SS    float_value_of_slv'   ( ex : 1995-08-10 11:50:00   3.050 ). No header is read (for now ....), \n"
   "                           the file extension is .slv This format is from the 10mn SONEL data (website)\n"
   "\n    -YMD   'YYYY:MM:DD HH:MM    float_value_of_slv'   ( ex : 1995:08:10 11:50   3.050 ). \n"
   "                           No header is read (for now ....), \n"
   "\n    -LIST           'cnes_date    float_value_of_slv'  where cnes_date is a julian date referenced to 1950/01/01 00:00:00  \n"
   "                           ????????? HEADER ???????,   NOT AVAILABLE YET  \n"
   "\n    -RMN            ????????????? DESCRIBE RMN FORMAT ?????????????\n"
   "\n    -PAB            ????????????? DESCRIBE PAB FORMAT ?????????????\n"
   "\n    -PROFILERS       Reads a BINARY file with the following format: 'nbin ntime YYYY MM DD HH MM SS    bin1 bin2 etc.' \n"
   "                         (All values are Float) No header is read (for now ....), \n"
  );
  __ERR_BASE_LINE__("exiting\n");exit(-1); /** \endcode */
/**
* \todo 11/04/2011 : Clement : complete the help function of mgr-converter.cpp, describe RMN and PAB formats
*/
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T> int tide_atlasSG2positions_template(const char *filename, const char *v1, const char *v2, vector<mgr_t> mgr, int nmgr, T *a,T *G, T mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*-----------------------------------------------------------------------------
  extract amplitudes and phases of structured atlas at given locations
-----------------------------------------------------------------------------*/
{
  complex<T> cmask,*z=NULL;
  int status,n;
  int npositions=nmgr;
  
  double *lon=new double[npositions];
  double *lat=new double[npositions];
  for(n = 0; n < npositions; n++) {
    lon[n]=mgr[n].loc.lon;
    lat[n]=mgr[n].loc.lat;
    }
    
  cmask=complex<T>(9999.,9999.);
  z=new complex<T>[npositions];

//   status=tide_atlasSG2positions(filename, v1, v2, lon, lat, npositions, z, cmask,verbose,wave);
  if(status!=NC_NOERR) goto finished;

  for(n = 0; n < npositions; n++) {
    if(z[n]==cmask) {
      a[n] =  mask;
      G[n] =  mask;
      }
    else {
      a[n] =  abs(z[n]);
/*----------------------------------------------------------------------
      convert into phase lag : 0 < G < 360 */
      G[n] = -arg(z[n]) * r2d;
      while(G[n] < 0){
        G[n] += 360;
        }
      while(G[n] > 360.0){
        G[n] -= 360;
        }
      }
    }

finished:
  delete[] z;
  delete[] lon;
  delete[] lat;
  return (status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int mgr_LoadAtlas(vector<mgr_t> mgr, int nmgr, const char **varnames, char **waves, char *atlas_directory, char *atlas_convention, double scale, char *meshfile, 
		   int structured, char *discretisation, int iteration, int level, int strict, double** & a, double** & G, double mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,k,status;
  char *model;
  spectrum_t spectrum;
  
  astro_angles_t astro_angles;
  date_t date;
  
  a=0;
  G=0;

  date.year=1950;
  date.month=1;
  date.day=1;
  date.second=0.0;
  spectrum.init(initialize_tide(&astro_angles,date),waves);

  for (i=0; i<spectrum.n; i++) {
    spectrum.waves[i].init();
    printf ("wave: %10s, pulsation: %12.6f degrees/h \n", spectrum.waves[i].name,spectrum.waves[i].omega);
    }


  exitIfNull(
    a=new double*[spectrum.n]
    );
  
  exitIfNull(
    G=new double*[spectrum.n]
    );

  for (k=0;k<spectrum.n;k++) {
    exitIfNull(
      a[k]=new double[nmgr]
      );
  
    exitIfNull(
      G[k]=new double[nmgr]
      );
    }

/* *-----------------------------------------------------------------------------
  read tidal atlas database and interpolate*/
  atlas_grid_or_mesh_t gm;
  
  for (k=0;k<spectrum.n;k++) {
    status=tide_decode_atlasname(atlas_directory,atlas_convention,waves[k], 0,&model);
//    printf("validate: treating %s wave from %s \n",wave[k],model);
    if(structured==1) {
      printf("validate: treating %s wave from %s \n",waves[k],model);
      status=tide_atlas2mgr(model,varnames[0],varnames[1],mgr,nmgr,a[k],G[k], mask,&gm);
//      status=tide_atlasSG2mgr(model,varnames[0],varnames[1],mgr,nmgr,a[k],G[k], mask);
//      status=tide_atlasSG2mgr(model,varnames[0],varnames[1],mgr,nmgr,a[k],G[k], mask, 0, strict);
      }
    else {
//       if(meshfile==0)
//         status=extract_UGatlas(model,varnames,mgr,nmgr,a[k],G[k],mask,discretisation, iteration, dmax, level, strict, 0);
//       else
//         status=extract_UGatlas_ASCII(model,meshfile,mgr,nmgr,a[k],G[k],mask,discretisation, iteration, dmax, level, strict, 0);
      }
    }
  if(status!=0) goto error;
  
  for (i=0;i<nmgr;i++){
    for (k=0;k<spectrum.n;k++) {
      if(a[k][i]!=mask) a[k][i]*=scale;
      }
    }
  
  return(0);
error:
  return(-1);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int mgr_LoadMultipleAtlas(vector<mgr_t> mgr, int nmgr,char **waves)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,status;
  char *varnames[2], *atlas_directory, *atlas_convention;
  double scale;
  char *meshfile; 
  int structured;
  char *discretisation;
  int iteration, level, strict;
  double **a[10], **G[10];
  double mask=999.9;
  
  iteration=-1;
  level=-1;
  strict=0;
  
  scale=1.0;
  
  meshfile=0;
  discretisation=0;
  
  structured=1;
  varnames[0]=strdup("Ha");
  varnames[1]=strdup("Hg");
  
  atlas_directory=strdup("/home/data/tides/GOT4.8/netcdf");
  atlas_convention=strdup("WAVE.GOT4.8.nc");
  
  k=0;
  status=mgr_LoadAtlas(mgr, nmgr, (const char **) varnames, waves, atlas_directory, atlas_convention, scale, meshfile, 
                       structured, discretisation, iteration, level, strict, a[k], G[k],mask);
  
  atlas_directory=strdup("/home/data/tides/FES2012/atlas-FES2012.2013-01-19");
  atlas_convention=strdup("WAVE.FES2012.elev.nc");
  
  k=1;
  status=mgr_LoadAtlas(mgr, nmgr, (const char **) varnames, waves, atlas_directory, atlas_convention, scale, meshfile, 
                       structured, discretisation, iteration, level, strict, a[k], G[k],mask);
  
  return(status);
  
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int mgr_colocated(vector<mgr_t> mgr, int nmgr, char *keep)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,j,k,l,m,n,count,status;
  double t1,p1,t2,p2,d;
  
  for (m=0;m<nmgr;m++) {
    if(keep[m]==0) continue;
    t1=mgr[m].loc.lon;
    p1=mgr[m].loc.lat;
    for (n=0;n<nmgr;n++) {
      t2=mgr[n].loc.lon;
      p2=mgr[n].loc.lat;
      d=geo_km(t1,p1,t2,p2);
      if(d<5.0) {
        keep[n]=0;
        }
      }
    }
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  vector<int> *mgr_colocated(vector<mgr_t> mgr, int nmgr)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,j,k,l,m,n,status;
  double t1,p1,t2,p2,d;
  vector<int> *clusters=new vector<int>[nmgr];
  bool *eligible=new bool[nmgr];
  int nclusters=-1;
  
  for (m=0;m<nmgr;m++) eligible[m]=true;
  
//  clusters[0].push_back(0);
  
  for (m=0;m<nmgr;m++) {
    if(!eligible[m]) continue;
    nclusters++;
    t1=mgr[m].loc.lon;
    p1=mgr[m].loc.lat;
    clusters[nclusters].push_back(m);
    eligible[m]=0;
    for (n=m+1;n<nmgr;n++) {
      if(!eligible[n]) continue;
      t2=mgr[n].loc.lon;
      p2=mgr[n].loc.lat;
      d=geo_km(t1,p1,t2,p2);
      if(d<5.0) {
        clusters[nclusters].push_back(n);
        eligible[n]=0;
        }
      }
    }
    
  for (k=0;k<nclusters;k++) {
    if(clusters[k].size()==1) continue;
    printf("\n");
    for (l=0;l<clusters[k].size();l++) {
      m=clusters[k][l];
      double t=mgr[m].loc.lon;
      double p=mgr[m].loc.lat;
      printf("%2d %2d %4d %9.3lf %9.3lf",k,l,m,t,p);
      j=mgr[m].wave_index("M2");
      float a,G;
      if(j==-1) {
        a=-1;
        G=-1;
        }
      else {
        a=100*mgr[m].data[j].amp;
        G=mgr[m].data[j].phi;
        }
      printf("    M2 %6.1f %6.1f",a,G);
      j=mgr[m].wave_index("K1");
      if(j==-1) {
        a=-1;
        G=-1;
        }
      else {
        a=100*mgr[m].data[j].amp;
        G=mgr[m].data[j].phi;
        }
      printf("    K1 %6.1f %6.1f",a,G);
      printf(" : %s\n",mgr[m].name);
      }
   }
  
  return(clusters);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int mgr_estuarine(vector<mgr_t> mgr, int nmgr, grid_t topogrid, float*topo, float topomask, char *keep)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,j,k,l,m,n,count,status;
  double t1,p1,t2,p2,d;
  double alpha,radius,max_radius,x,y;
  float positive,negative,average;
  float h;
  int last;
  
  max_radius=0.05;
  
  for (m=0;m<nmgr;m++) {
    if(keep[m]==0) continue;
    t1=mgr[m].loc.lon;
    p1=mgr[m].loc.lat;
    positive=0;
    negative=0;
    average=0;
    for (alpha=0;alpha<2*M_PI;alpha+=0.1) {
      last=0;
      for(radius=max_radius/10.;radius<max_radius;radius+=max_radius/10.) {
        x=t1+radius*cos(alpha);
        y=p1+radius*sin(alpha);
        x=map_recale(topogrid,x);
        status=map_interpolation(topogrid, topo,topomask,x,y,&h);
        if(h!=topomask) {
          if(h>0) {
//            positive+=max_radius-radius;
            positive+=1;
            last=1;
            }
          else  {
//	    negative+=max_radius-radius;
            if(last==1) break;
            negative+=1;
            average+=h;
            last=-1;
            }
          }
        }
      }
//    if(0.5*positive>negative) keep[m]=0;
//    if(average>-20) keep[m]=0;
    if(positive>0) keep[m]=0;
    }
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int consistency(vector<mgr_t> mgr, int nmgr)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,j,k,l,m,n,count,status;
  char *landmask=0,*eligible,*keep;
  
  eligible=new char[nmgr];
  keep    =new char[nmgr];
  for (m=0;m<nmgr;m++) {
    keep[m]=0;
    eligible[m]=0;
    }
    
  for (m=0;m<nmgr;m++) {
    if(mgr[m].loc.depth>-200.0) {
      keep[m]=1;
      eligible[m]=0;
      }
    }
  
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int process_scales(char *bathymetry, const char *scales, vector<mgr_t> & mgr, int nmgr)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,j,k,l,m,n,count,status;
  const char *v1="";
  const char *v2="";
  const char *v3="";
  grid_t l_grid, e_grid, grid;
  fcomplex *glorys, *error, e[2];
  fcomplex cmask;
  float *cardinal, *buffer, *L, rmask, Lmask, z;
  int mode=POLAR, verbose=0;
  double t,p,x,y,xx,yy;
  cdfgbl_t global;
  cdfvar_t info;
  float *lamda1,  *lamda2,  *lamda3;
  char *landmask=0,*eligible,*keep;
  mesh_t mesh;

  grid_t topogrid,Lgrid;
  float *topo=NULL,topomask,h;
  
  char obsfile[1024];
  char *waves[3]={"M2","K1","\0"};
  waves[2]=0;
  
  status=mgr_LoadMultipleAtlas(mgr, nmgr, waves);
  
  vector<int> *clusters=mgr_colocated(mgr, nmgr);
  
/**----------------------------------------------------------------------------
  Load bottom topography */
  if(bathymetry==NULL) bathymetry=strdup("/home/data/topography/gebco/gebco_08.grd");
  printf("#################################################################\n");
  printf("load bathymetry file : %s\n",bathymetry);
  status=grd_loadgrid(bathymetry,&topogrid);
  if(status !=0) {
    __OUT_BASE_LINE__("cannot load bathymetry file=%s\n",bathymetry);
    return(-1);
    }
  topo= new float[topogrid.nx*topogrid.ny];
  if (topo == NULL) {
    __ERR_BASE_LINE__("");perror("topo");
    return(-1);
    }

  status=  grd_loadr1(bathymetry,topogrid,topo,&topomask);

  for (m=0;m<nmgr;m++) {
    if(topo!=0) {
      x=mgr[m].loc.lon;
      y=mgr[m].loc.lat;
      x=map_recale(topogrid,x);
      status=map_interpolation(topogrid, topo,topomask,x,y,&h);
      mgr[m].loc.depth=h;
      }
    }

  L=new float[nmgr];
//  status=map_loadfield("/home/softs/data/tides/FES2004/netcdf/lamda.nc","lamda4_15", & Lgrid, &buffer, &Lmask); 
  status=map_loadfield_cdf("/home/models/FES2013/simulation-tides/lamda.nc", "lamda4_15", Lgrid, buffer, Lmask);
  if(status!=0) return(-1);
  
  
  for(k=0;k<10;k++) status=map_persistence(Lgrid, buffer, Lmask, 0.);

  for (m=0;m<nmgr;m++) {
    x=mgr[m].loc.lon;
    y=mgr[m].loc.lat;
    x=map_recale(Lgrid,x);
    status=map_interpolation(Lgrid, buffer,Lmask,x,y,&z);
    L[m]=z;
    if(z==Lmask) {
      printf("interpolation issue : %lf %lf\n",x,y);
      }
    }
  delete[] buffer;

  printf("#################################################################\n");
  printf("data parsing and editing\n");

  eligible=new char[nmgr];
  keep    =new char[nmgr];
  for (m=0;m<nmgr;m++) {
    keep[m]=0;
    eligible[m]=0;
    }
    
  for (m=0;m<nmgr;m++) {
    if(mgr[m].loc.depth>-200.0) {
      keep[m]=1;
      eligible[m]=0;
      }
    }
  
//   for (m=0;m<nmgr;m++) {
//     if(keep[m]==0) continue;
//     if(L[m]> 50.0) {
//       keep[m]=1;
//       eligible[m]=0;
//       }
//     else {
//       keep[m]=0;
//       eligible[m]=0;
//       }
//     }
  status=mgr_estuarine(mgr, nmgr, topogrid, topo, topomask, keep) ;
  count=0;
  for (m=0;m<nmgr;m++) {
    mgr[m].mindex=keep[m];
    if(keep[m]==1) count++;
    }

  delete[] eligible;
  delete[] keep;
  delete[] topo;
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,j,k,l,nitems,w;
  int n,status;
  char *input=0,*s=0;
  char *output=NULL,*header=0,*wave=0;
  char *bathymetry=NULL, *scales=NULL;
  char *polygons=NULL;
  date_t start,final,reference,start_date;
  statistic_t stat;
  tseries_t *serie=0;
  int nseries;
  mooring_t mooring;
  string input_format="",output_format="";
  double resolution_deep=200.0, resolution_shelf=50.0;
  string rootname,sformat;

  fct_echo( argc, argv);

  n=1;
  while (n < argc) {
    const char *keyword=strdup(argv[n]);
    
    if(isFormat(argv,&n,&input_format,&output_format)){
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

        case 'p' :
          polygons=strdup(argv[n+1]);
          n++;
          n++;
          break;

//         case 's' :
//           s= strdup(argv[n+1]);
//           n++;
//           n++;
//           sscanf(s,"%d/%d/%d",&start.day,&start.month,&start.year);
//           break;
// 
//        case 'f' :
//           s= strdup(argv[n+1]);
//           n++;
//           n++;
//           sscanf(s,"%d/%d/%d",&final.day,&final.month,&final.year);
//           break;

        case 'o' :
          output= strdup(argv[n+1]);
          n++;
          n++;
          break;
          
        case 'w' :
          wave= strdup(argv[n+1]);
          n++;
          n++;
          break;
          
        case 'h' :
          print_help(argv[0]);

        case '-' :
          if(strstr(argv[n],"--header=")!=0){
            header=strdup(argv[n]);
            n++;
            break;
            }
          break;

        default:
          printf("unknown option %s\n",keyword);
          print_help(argv[0]);
          exit(-1);
        }
        break;

      default:
        if(input==NULL) {
          input= strdup(argv[n]);
          n++;
          }
        else {
          printf("redundant option %s\n",keyword);
          print_help(argv[0]);
          exit(-1);
          }
        break;
      }
    
    }
  
  status= init_mgrh_formats();

  printf("#################################################################\n");
  printf("load tide gauge file : %s\n",input);
  vector<mgr_t> mgr;
  int nmgr;
  nmgr=mgr_load((const char *) input, mgr);
  if(nmgr==0) {
    status=-1;
    goto error;
    }
 
  printf("#####################################################################\n");
  printf("process data \n");

  if(polygons!=0) {
    status=mgr_exclude(mgr,polygons,PLG_POINT_INTERIOR);
    if(status==-1) goto error;
    }

  if(bathymetry!=0) {
    status= process_scales(bathymetry, scales, mgr, nmgr);
    if(status==-1) goto error;
    }

  if(output==0) {
    rootname=input;
    size_t pos;
    pos=rootname.rfind(".mgr");
    if (pos!=string::npos) {
      rootname.erase(pos,4);
      }
    pos=rootname.rfind(".nc");
    if (pos!=string::npos) {
      rootname.erase(pos,3);
      }
    }

  if(output==0) {
    if(output_format.compare("ASCII") == 0){
      rootname+=".mgr";
      }
    else if(output_format.compare("LEGOS-ASCII") == 0){
      rootname+=".mgr";
      }
    else if(output_format.compare("NETCDF") == 0){
      rootname+=".nc";
      }
    else if(output_format.compare("LEGOS-NETCDF") == 0){
      rootname+=".nc";
      }
    else if(output_format.compare("LEGOS-OBS") == 0){
      rootname+=".obs";
      }
    output=strdup(rootname.c_str());
    }

  printf("#################################################################\n");
  printf("save tide gauge  file : %s\n",output);
  status=mgr_save((const char *) output, mgr, output_format.c_str(),wave);
  for (size_t m=0;m<nmgr;m++) {
    mgr[m].mindex=(mgr[m].mindex==0);
    }
  status=mgr_save("rejected.mgr", mgr, output_format.c_str(),wave);


  goto end;
error:
end:
  __OUT_BASE_LINE__(" ^^^^^^^^^^^^^ end of computation ^^^^^^^^^^^^^\n");
  exit(status);
}
