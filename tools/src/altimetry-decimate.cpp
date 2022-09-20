
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

  int data_error(char *poly,char *wave, int kk, char *zone, vector<mgr_t> mgr,int nmgr,xover_t *xover,int nxovers,grid_t topogrid,float *topo,float topomask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int *track,*subtrack,ntracks,ntics;
  int **linked,*nlinked,**list,*card,count,npt,i,j,m,mm,n,nn,k,l,ll,w,status;
  int n1,n2,n3;
  float *alias_error,*coastal_error,*analysis_error,*error;
  double *x,*y,*h,*t,*d,*lon, *lat,*buble,tt;
  float error_min, *rbuffer[10], rmask=data_rmask,a,p;
  std::complex<float> *cbuffer,cmask=-1;
  std::complex<float> *c;
  grid_t grid;
  float **covariance;
  char obsname[1024];

  cdfvar_t variable;
  pocgrd_t ncgrid;
  char output[1024];

  xerror_t xerror;

//  hypermatrix_t Cd;

  buble   = new double [nmgr];
  error   = new float [nmgr];
  alias_error   = new float [nmgr];
  coastal_error = new float [nmgr];
  analysis_error = new float [nmgr];
  subtrack = new int   [nmgr];
  track    = new int   [nmgr];
  x     = new double[nmgr];
  y     = new double[nmgr];


/* *----------------------------------------------------------------------------
  initialisation*/
  for(n = 0; n < nmgr; n++) {
    x[n]=geo_recale(mgr[n].loc.lon,0.,180.);
    y[n]=mgr[n].loc.lat;
    analysis_error[n]=rmask;
    alias_error[n]=rmask;
    coastal_error[n]=0.;
    error[n]=rmask;
    buble[n]=20.;
    }

  xerror.n=nmgr;
  xerror.alias=alias_error;
  xerror.alongtrack=coastal_error;
  xerror.harmonic=analysis_error;

  xerror.cut=new float[nmgr];

  xerror.zraw=new complex<float>[nmgr];
  xerror.zlwf=new complex<float>[nmgr];
  xerror.zhgf=new complex<float>[nmgr];
  xerror.smoothed=new complex<float>[nmgr];

  for(n = 0; n < nmgr; n++) {
    xerror.zraw[n]=data_cmask;
    xerror.zlwf[n]=data_cmask;
    xerror.zhgf[n]=data_cmask;
    }

/* *----------------------------------------------------------------------------
  error estimate: harmonic analysis error bar*/
  status= compute_analysis_error(wave, mgr, nmgr, analysis_error);
  for(n = 0; n < nmgr; n++) {
    if(analysis_error[n]!=rmask) {
      error[n]=analysis_error[n];
      }
    }

/* *----------------------------------------------------------------------------
  error estimate: along-track HF variability*/
  status= compute_coastal_error(wave, kk, mgr, nmgr, xerror ,buble,list,card);
  for(n = 0; n < nmgr; n++) {
    if(coastal_error[n]!=rmask) {
      error[n]=MAX(error[n],coastal_error[n]);
      }
    }

  sprintf(obsname,"%s.obs",wave);
/* *----------------------------------------------------------------------------
  COMAPI : no covariance !!! */
  status=mgr_save_obs4assim(obsname,poly,wave, nmgr ,mgr, covariance, 0, 0);

  sprintf(obsname,"%s.lgd",wave);
  {
  float *all[4];
  all[0]=analysis_error;
  all[1]=alias_error;
  all[2]=coastal_error;
  all[3]=error;
  status=mgr_Constants2Legend_00(obsname,wave, nmgr ,mgr, all, covariance, track, linked, nlinked);
  }

  for(n = 0; n < nmgr; n++) {
    w=mgr[n].wave_index(wave);
    if(w!=-1) {
      a=abs(xerror.zlwf[n]);
      p=arg(xerror.zlwf[n])*r2d;
      mgr[n].data[w].amp=a;
      mgr[n].data[w].phi=p;
      mgr[n].data[w].error=abs(xerror.zhgf[n])/sqrt(xerror.cut[n]/6.0);
      }
    }

  delete[] x;
  delete[] y;
  delete[] buble;
  delete[] error;
  delete[] coastal_error;
  delete[] alias_error;
  delete[] analysis_error;
  delete[] subtrack;
  delete[] track;

  for(k=0;k<4;k++) delete[] rbuffer[k];

  return(0);
  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int process(char *bathymetry, const char *glorys_file, const char *error_file, const char *scales_file,
              const char *wave, vector<mgr_t> & mgr,
              double resolution_deep, double resolution_shelf, double max_error, double max_ratio, double min_depth, double max_depth, char *poly, char *output)

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
//  if(bathymetry==NULL) bathymetry=strdup("/home/softs/genesis/dbg/data/topography/data/GEBCO_08.grd");
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

  for (m=0;m<mgr.size();m++) {
    if(topo!=0) {
      x=mgr[m].loc.lon;
      y=mgr[m].loc.lat;
      x=map_recale(topogrid,x);
      status=map_interpolation(topogrid, topo,topomask,x,y,&h);
      mgr[m].loc.depth=h;
      }
    }

/**----------------------------------------------------------------------------
  cluster cells */
  printf("#################################################################\n");
  printf("create parsing grid\n\n");
  grid=map_Rgrid(frame_t (-180.,180.,-90.,90.), 2., 2., 2);
//  grid=map_Rgrid(frame_t (-180.,180.,-90.,90.), 1., 1.);
//  grid=map_Rgrid(frame_t (-180.,180.,-90.,90.), 1./12, 1./12);
  
  scales=new float[grid.nx*grid.ny];
  status=map_export(l_grid, lamda3, rmask, grid, scales, rmask, 0);
  
  landmask=new char [grid.nx*grid.ny];
  cardinal=new float[grid.nx*grid.ny];
 
  for (j=0;j<grid.ny;j++) {
    for (i=0;i<grid.nx; i++) {
      m=grid.nx*j+i;
      cardinal[m]=0;
      if(scales[m]==rmask) cardinal[m]=rmask;
      }
    }
  for (m=0;m<mgr.size();m++) {
    x=mgr[m].loc.lon;
    y=mgr[m].loc.lat;
    x=map_recale(grid,x);
    status=map_index(grid,x,y,&k,&l);
    if(status!=0) {
      printf("troubles\n");
      }
    n=grid.nx*l+k;
    if(cardinal[n]!=rmask) cardinal[n]++;
    }
  
  status=save_SG("cardinal.nc", grid, cardinal, rmask, "cardinal", "dimensionless", "cardinal", 1);
  status=save_SG("cardinal.nc", grid, scales,   rmask, "scales",   "kilometers",    "tidal_length_scale", 0);
  
#if 0
  landmask=new char[grid.nx*grid.ny];
 
  for (j=0;j<grid.ny;j++) {
    for (i=0;i<grid.nx; i++) {
      x=map_grid_x(grid,i,j);
      y=map_grid_y(grid,i,j);
      x=map_recale(topogrid,x);
      status=map_interpolation(topogrid,topo,topomask,x,y,&h);
      m=grid.ny*j+i;
      if(h>0) landmask[m]=0;
      else    landmask[m]=1;
      }
    }

    status=fe_FGrid2Quadrangle(grid, landmask, mesh,"quadrangle.nc",debug);
#endif
  
#if 0
  vector <size_t> *clusters=new vector <size_t>[nmgr];
  for (m=0;m<nmgr;m++) {
    x=mgr[m].loc.lon;
    y=mgr[m].loc.lat;
    x=map_recale(grid,x);
    status=map_index(grid,x,y,&k,&l);
    n=grid.ny*l+k;
    clusters[m].push_back((size_t) n);
    }
#endif
  
  printf("#################################################################\n");
  printf("data parsing and editing \n\n");

  printf("#################################################################\n");
  printf("error estimate absolute threshold:  %lf mm\n",1000.*max_error);
  printf("error estimate relative threshold:  %lf %\n",100.*max_ratio);
  inside=0;
  if(poly!=0) {
    inside=mgr_select(poly ,mgr);
    if(inside==0) return(-1);
    }
  eligible=new char[mgr.size()];
  keep    =new char[mgr.size()];
  for (m=0;m<mgr.size();m++) {
    keep[m]=0;
    eligible[m]=0;
    w=mgr[m].wave_index(wave);
    if(w==-1) continue;
    if(mgr[m].data[w].error>max_error) {
      eligible[m]=0;
      }
    else {
      eligible[m]=1;
      }
    if(mgr[m].data[w].error/mgr[m].data[w].amp>max_ratio)  eligible[m]=0;
    if(mgr[m].loc.depth>min_depth) eligible[m]=0;
    if(mgr[m].loc.depth<max_depth) eligible[m]=0;
    if(inside!=0) {
      if(inside[m]==PLG_POINT_EXTERIOR) eligible[m]=0;
      }
    }
    
//   for (m=0;m<nmgr;m++) {
//     if(mgr[m].loc.depth>-200.0) {
//       keep[m]=1;
//       eligible[m]=0;
//       }
//     }
  
//   for (m=0;m<nmgr;m++) {
//     if(eligible[m]==1) {
//       keep[m]=1;
//       }
//     }
    
  printf("#################################################################\n");
  printf("decimating: deep ocean resolution %lf km, shelf seas resolution %lf km\n",resolution_deep, resolution_shelf);
  for (m=0;m<mgr.size();m++) {
    vector <size_t> clusters;
    double d,dmin=1.e+10,dmax;
    int    size, target;
    double cx,cy;
    if(keep[m]    ==1) continue;
    if(eligible[m]==0) continue;
    x=mgr[m].loc.lon;
    y=mgr[m].loc.lat;
    if(mgr[m].loc.depth>-500.0) {
      dmax=resolution_shelf;
      }
    else {
      dmax=resolution_deep;
      }
    for (n=0;n<mgr.size();n++) {
      if(eligible[n]==0) continue;
      if(keep[n]    ==1) continue;
      xx=mgr[n].loc.lon;
      yy=mgr[n].loc.lat;
      d=geo_haversin_km(x,y,xx,yy);
      if(d>dmax) continue;
      clusters.push_back((size_t) n);
      }
    cx=cy=0.0;
    for(k=0;k<clusters.size();k++) {
      n=clusters[k];
      xx=mgr[n].loc.lon;
      yy=mgr[n].loc.lat;
      xx=geo_recale(xx,x,180.0);
      cx+=xx;
      cy+=yy;
      }
    size=clusters.size();
    cx/= clusters.size();
    cy/= clusters.size();
    for(k=0;k<clusters.size();k++) {
      n=clusters[k];
      xx=mgr[n].loc.lon;
      yy=mgr[n].loc.lat;
      d=geo_haversin_km(cx,cy,xx,yy);
      if(d<dmin) {
        dmin=d;
        target=n;
        }
      }
    for(k=0;k<clusters.size();k++) {
      n=clusters[k];
      eligible[n]=0;
      }
    keep    [target]=1;
    eligible[target]=0;
    }
  
  
  for (m=0;m<mgr.size();m++) {
    mgr[m].mindex=keep[m];
    }
  status=mgr_save_ascii(output, mgr);

//   for (m=0;m<nmgr;m++) {
//     mgr[m].destroy();
//     }
  int nmgr=mgr_load(output, mgr);
  
  sprintf(obsfile,"%s.obs",wave);
  status=mgr_save_obs(obsfile, mgr,wave);
  
  
  for(n = 0; n < mgr.size(); n++) {
    w=mgr[n].wave_index(wave);
    t=mgr[n].loc.lon;
    p=mgr[n].loc.lat;
    if(w!=-1) {
      t=map_recale(grid,t);
      status=map_interpolation(e_grid, glorys, cmask, t, p, &e[0]);
      status=map_interpolation(l_grid, scales, rmask, t, p, &length);
//       a=abs(xerror.zlwf[n]);
//       p=arg(xerror.zlwf[n])*r2d;
//       mgr[n].data[w].amp=a;
//       mgr[n].data[w].phi=p;
//       mgr[n].data[w].error=abs(xerror.zhgf[n])/sqrt(xerror.cut[n]/6.0);
      }
    }
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

  if(output==0) output=strdup("filtered.mgr");
			      
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
  if(nmgr<1) goto error;

 
  printf("#####################################################################\n");
  printf("process data \n\n");

//  error_file=strdup("/home/models/FES2012/simulation-tides/spectral-DNP1xLGP2-61-TOPO17x1-WD100-g-variable/M2.vector-difference.nc");
  error_file=strdup("/home/models/FES2012/simulation-tides/spectral-DNP1xLGP2-61-TOPO17LGP1x1WD100-g-varying/M2.vector-difference.nc");
  status= process(bathymetry, "", error_file, "", wave[0], mgr, resolution_deep, resolution_shelf, max_error, max_ratio,  min_depth, max_depth, poly, output);
  if(status!=0) goto error;
  
  
  if(merge!=0) {
    vector<mgr_t> additionals;
    int nadditionals;
    printf("#################################################################\n");
    printf("import additional harmonic constants from : %s\n",merge);
    nmgr=mgr_load(merge, mgr);
    printf("#stations : %d\n\n",nmgr);
    status=mgr_exclude(mgr,poly);
    nadditionals=mgr_load("filtered.mgr", additionals);
    status=mgr_concat(mgr,additionals);
    status=mgr_save("concat.mgr", mgr,1.0);
    sprintf(obsfile,"%s.concat.obs",wave[0]);
    status=mgr_save_obs(obsfile, mgr,wave[0]);
    }
  
/* *----------------------------------------------------------------------------
  Warning : code will failed if track number not in station name*/
#if 0
  status=mgr_save("filtered.mgr", nmgr, mgr, (float) 20.0);

  nmgr=mgr_load("filtered.mgr", &mgr);

  list= order(mgr, nmgr);

  if(output==NULL) {
    report=strdup("report.out");
    }
  else {
    report=strdup(output);
    }

  status=mgr_save("smoothed.mgr", nmgr, mgr, (float) 0.10);
#endif

  __ERR_BASE_LINE__("%s -computation complete ^^^^^^^^^\n",argv[0]);
  exit(0);
error:
  __ERR_BASE_LINE__("exiting\n");exit(-1);
}
