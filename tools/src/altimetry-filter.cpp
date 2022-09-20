
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

/**************************************************************************

  altimetry tidal data error

***************************************************************************/

#define MAIN_SOURCE

#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "tools-define.h"
#include "tools-structures.h"

#include "poc-time.h"
#include "list.h"
#include "map.h"
#include "mgr.h"
#include "filter.h"
#include "tides.h"
#include "netcdf-proto.h"
#include "geo.h"
#include "legend.h"
#include "polygones.h"
#include "statistic.h"
#include "filter.h"
#include "functions.h"

#include "legend.def"

#include "grd.h"

#include "matrix.h"

#include "altimetry.h"


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int process(char *bathymetry, const char *glorys_file, const char *error_file, const char *scales_file,
              const char *wave, vector<mgr_t> mgr,
              double resolution_deep, double resolution_shelf, double max_error, double max_ratio, double min_depth, double max_depth, char *poly)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,j,k,l,m,n,w,status;
  const char *v1="";
  const char *v2="";
  const char *v3="";
  grid_t l_grid, e_grid, grid;
  fcomplex *glorys, *error, e[2];
  fcomplex cmask;
  float *scales, *cardinal, length, rmask;
  int mode=POLAR, verbose=0;
  double t,p,x,y,xx,yy;
  cdfgbl_t global;
  cdfvar_t info;
  float *lamda1,  *lamda2,  *lamda3;
  char *landmask=0,*eligible,*keep;
  mesh_t mesh;

  grid_t topogrid;
  float *topo=NULL,topomask,h;
  int *inside;
  char obsfile[1024];
  
  int nmgr=mgr.size();
  
/**----------------------------------------------------------------------------
  horizontal scale, computation */
#if 0
  printf("#################################################################\n");
  printf("load lengthscale characteristic file : %s\n",error_file);
  status=load_atlas(error_file, "a", "G", &l_grid, &error, &cmask, mode, verbose);
  printf("status : %d\n",status);
  if(status!=0) return(-1);
  
  lamda1=lamda2=lamda3=0;
  status=map_LengthScale(l_grid, error, cmask, lamda1, lamda2, lamda3);
  rmask=cmask.real();
  status=save_SG("lamda.nc", l_grid, lamda1, rmask, "lamda1_15", "kilometers", "tidal_length_scale", 1);
  status=save_SG("lamda.nc", l_grid, lamda2, rmask, "lamda2_15", "kilometers", "tidal_length_scale", 0);
  status=save_SG("lamda.nc", l_grid, lamda3, rmask, "lamda3_15", "kilometers", "tidal_length_scale", 0);
#endif

/**----------------------------------------------------------------------------
  observational errors, load */
#if 0
  printf("#################################################################\n");
  printf("load ogcm analysis file : %s\n",glorys_file);
  status=load_atlas(glorys_file, v1, v2, &e_grid, &glorys, &cmask, mode, verbose);
  printf("status : %d\n",status);
  if(status!=0) return(-1);
#endif
  
/**----------------------------------------------------------------------------
  horizontal scale, load */
//  status=map_loadfield_cdf (scales_file, v3, grid, scales, rmask);


    
/**----------------------------------------------------------------------------
  Load bottom topography */
#if 0
  if(bathymetry==NULL) bathymetry=strdup("/home/softs/genesis/data/topography/data/gebco_08.grd");
  printf("#################################################################\n");
  printf("load bathymetry file : %s\n\n",bathymetry);
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
#endif

#if 0
/* *-----------------------------------------------------------------------------
  read tidal atlas database and interpolate*/
  for (k=0;k<spectrum.n;k++) {
    status=tide_decode_atlasname(atlas_directory,atlas_convention,wave[k], 0,&model);
//    printf("validate: treating %s wave from %s \n",wave[k],model);
    if(structured==1) {
      printf("validate: treating %s wave from %s \n",wave[k],model);
      status=extract_SGatlas(model,varnames,mgr,nmgr,a[k],G[k],mask, strict);
      }
    else {
      if(meshfile==0)
        status=extract_UGatlas(model,varnames,mgr,nmgr,a[k],G[k],mask,discretisation, iteration, dmax, level, strict, 0);
      else
        status=extract_UGatlas_ASCII(model,meshfile,mgr,nmgr,a[k],G[k],mask,discretisation, iteration, dmax, level, strict, 0);
      }
    }
  if(status!=0) goto error;
#endif

/**----------------------------------------------------------------------------
  filter */   
  printf("#################################################################\n");
  printf("decimating: deep ocean resolution %lf km, shelf seas resolution %lf km\n",resolution_deep, resolution_shelf);
  
  double *r=new double[nmgr];
  point_t track;
  complex<float> *elevation=new complex<float>[nmgr], *LF=0, mask;
  double  x0=mgr[0].loc.lon;
  double  y0=mgr[0].loc.lat;
  
  
  status=mgr_save("M2.obs", mgr, "LEGOS-OBS", "M2");
//   status=mgr_Constants2legend_01("M2.lgd", "M2", mgr);

  for (m=0;m<nmgr;m++) {
    x=mgr[m].loc.lon;
    y=mgr[m].loc.lat;      
    
    w=mgr[m].wave_index(wave);
    float a=mgr[m].data[w].amp;
    float G=mgr[m].data[w].phi;
    G*=d2r;
    elevation[m]=polar(a,-G);
    r[m]=geo_haversin_km(x0,y0,x,y);
    }
  
  double L=200.;
  status=filter1D_BF(LANCZOS, L, r, nmgr, elevation, mask, LF);
  
  for (m=0;m<nmgr;m++) {
    w=mgr[m].wave_index(wave);
    complex<float> z=elevation[m]-LF[m];
    float a= abs(z);
    float G=-arg(z);
    G/=d2r;
    mgr[m].data[w].amp=a;
    mgr[m].data[w].phi=G;
    }
  status=mgr_save_ascii("HF.mgr", mgr);
  status=mgr_save("M2-HF.obs", mgr, "LEGOS-OBS", "M2");
  
  for (m=0;m<nmgr;m++) {
    w=mgr[m].wave_index(wave);
    complex<float> z=LF[m];
    float a= abs(z);
    float G=-arg(z);
    G/=d2r;
    mgr[m].data[w].amp=a;
    mgr[m].data[w].phi=G;
    }
  status=mgr_save_ascii("LF.mgr", mgr);
  status=mgr_save("M2-LF.obs", mgr, "LEGOS-OBS", "M2");

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double mask=999.9;
  double tau;
  double *residual,*buffer;
  double  zi,zr,mean[10],rms[10],count;
  double  **a=NULL,**G=NULL,a1,p1,a2,p2,d,da,dG;
  int i,j,k,kk,l,m,n,fmt,status,tg;
  FILE *in=NULL,*out=NULL;
  char *report,hfile[1204],tmp[256], datafile[256], lgdfile[1024];
  char *data, *output=NULL, *keyword, *poly=NULL;
  char *mgrfile=NULL,*zone=NULL,*xover_base=NULL,*merge=NULL;
  fcomplex meanc;
  char c;
  spectrum_t spectrum,reference_spectrum;
  char *wave[100],obsfile[1024];
  int found,nwave=0;
  date_t date;
  vector<mgr_t> mgr;
  int i1,i2,m1,m2,n1,n2,nmgr;
  char *atlas_directory=NULL,*atlas_convention=NULL;
  fcomplex cmisfit;
  int *list,indice[4];
  double t1,t2,t_xover,p_xover;
  double dmin;
  int neighbour,nlegends;
  legend_t   *legends;
  legend02_t *legends02;

  char *bathymetry=NULL;
  double x,y,x1[2],y1[2],x2[2],y2[2],z,radius=10.0;
  std::complex<float> tide[2],c1,c2,cmask;
  int maxcycles=1;
  xover_t *xover;
  xoverbase_t full;
  int nxover;
  char *comments;
  char *error_file;
  double resolution_deep=200.0, resolution_shelf=20.0,max_error=0.005, max_ratio=1.e+35, min_depth=0.0, max_depth=-20000.;
  
  fprintf(stderr,"%s  -starting computation *********\n",argv[0]);

  fct_echo( argc, argv);

  n=1;
  while (n < argc) {
    keyword=strdup(argv[n]);
    switch (keyword[0]) {
      case '-':
        if(strcmp(keyword,"-maxcycles")==0) {
          sscanf(argv[n+1],"%d",&maxcycles);
          n++;
          n++;
          continue;
          }
        if(strcmp(keyword,"-min_depth")==0) {
          sscanf(argv[n+1],"%lf",&min_depth);
          n++;
          n++;
          continue;
          }
        if(strcmp(keyword,"-max_depth")==0) {
          sscanf(argv[n+1],"%lf",&max_depth);
          n++;
          n++;
          continue;
          }
        if(strcmp(keyword,"-resolution_deep")==0) {
          sscanf(argv[n+1],"%lf",&resolution_deep);
          n++;
          n++;
          continue;
          }
        if(strcmp(keyword,"-resolution_shelf")==0) {
          sscanf(argv[n+1],"%lf",&resolution_shelf);
          n++;
          n++;
          continue;
          }
        if(strcmp(keyword,"-error")==0) {
          sscanf(argv[n+1],"%lf",&max_error);
          n++;
          n++;
          continue;
          }
       if(strcmp(keyword,"-max_ratio")==0) {
          sscanf(argv[n+1],"%lf",&max_ratio);
          n++;
          n++;
          continue;
          }
        if(strcmp(keyword,"-merge")==0) {
          merge= strdup(argv[n+1]);
          n++;
          n++;
          continue;
          }
        switch (keyword[1]) {
        case 'c' :
          atlas_convention= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'b' :
          bathymetry= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'p' :
          poly= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'o' :
          output= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'g' :
          mgrfile= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'r' :
          sscanf(argv[n+1],"%lf",&radius);
          n++;
          n++;
          break;

        case 'z' :
          zone= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'x' :
          xover_base= strdup(argv[n+1]);
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
        wave[nwave]= strdup(argv[n]);
        printf("input wave=%s\n",wave[nwave]);
        spectrum.n=nwave+1;
        nwave++;
        n++;
        break;
      }
      free(keyword);
    }

  if (nwave==0) {
    printf("no tidal wave specified, abort...\n");
    exit(-1);
    }
  
/* *----------------------------------------------------------------------------
  Initialise tidal library */
  date=date_t(1950,1,1,0.0);
  astro_angles_t astro_angles;
  reference_spectrum = initialize_tide(&astro_angles,date);

  spectrum.nmax=nwave+1;
  spectrum.n=0;
  exitIfNull(spectrum.waves=new tidal_wave[spectrum.nmax]);
  
  for(j = 0; j < nwave; j++) {
    for(k = 0; k < reference_spectrum.n; k++) {
      if(strcmp(reference_spectrum.waves[k].name, wave[j]) == 0) {
        break;
        }
      }
    if(k != reference_spectrum.n) {
      addwave(&spectrum, reference_spectrum.waves[k]);
      }
    }

  for (i=0; i<spectrum.n; i++) {
    spectrum.waves[i].init();
    printf ("wave: %10s, pulsation: %12.6f degrees/h \n", spectrum.waves[i].name,spectrum.waves[i].omega);
    }

  for (i=0; i<spectrum.n; i++) {
    for (j=i+1; j<spectrum.n; j++) {
      tau=deltaOmega2separation(spectrum.waves[j].omega-spectrum.waves[i].omega);
      printf ("wave: %10s %10s, separation: %9.3f days \n", spectrum.waves[i].name,spectrum.waves[j].name,tau);
      }
    }
    
/**----------------------------------------------------------------------------
  Load tide gauge database */
  if(mgrfile==NULL) goto error;
  printf("#################################################################\n");
  printf("import harmonic constants from : %s\n",mgrfile);
  nmgr=mgr_load(mgrfile, mgr);
  printf("#stations : %d\n\n",nmgr);

 
  printf("#####################################################################\n");
  printf("process data \n\n");

//  error_file=strdup("/home/models/FES2012/simulation-tides/spectral-DNP1xLGP2-61-TOPO17x1-WD100-g-variable/M2.vector-difference.nc");
  error_file=strdup("/home/models/FES2012/simulation-tides/spectral-DNP1xLGP2-61-TOPO17LGP1x1WD100-g-varying/M2.vector-difference.nc");
  status= process(bathymetry, "", error_file, "", wave[0], mgr, resolution_deep, resolution_shelf, max_error, max_ratio,  min_depth, max_depth, poly);

  if(merge!=0) {
    vector<mgr_t> additionals;
    int nadditionals;
    printf("#################################################################\n");
    printf("import additional harmonic constants from : %s\n",merge);
    nmgr=mgr_load(merge, mgr);
    printf("#stations : %d\n\n",nmgr);
    status=mgr_exclude(mgr, poly);
    nadditionals=mgr_load("filtered.mgr", additionals);
    status=mgr_concat(mgr,additionals);
    status=mgr_save("concat.mgr", mgr,1.0);
    sprintf(obsfile,"%s.concat.obs",wave[0]);
    status=mgr_save_obs(obsfile, mgr,wave[0]);
    }
  
  __ERR_BASE_LINE__("%s -computation complete ^^^^^^^^^\n",argv[0]);
  exit(0);
error:
  __ERR_BASE_LINE__("exiting\n");exit(-1);
}
