
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
#include "grd.h"
#include "filter.h"
#include "functions.h"
#include "tides.h"
#include "netcdf-proto.h"
#include "statistic.h"


extern  int syspack_init(double **in,int **neighbours, int *card, int neq, hypermatrix_t *matrix);
extern  int syspack_inverse(double *in, double *out, int neq, hypermatrix_t matrix);


const char *header="\\documentclass[9pt,a4paper]{article} \n\
\\usepackage{here} \n\
\\usepackage{vmargin} \n\
\\usepackage[counterclockwise]{rotating} \n\
\\usepackage{supertabular} \n\
\\begin{document} \n\
\\newpage \n\
\\thispagestyle{empty} \n\
\\tiny \n";

const char *trailer="\\hline \n\
\\end{supertabular} \n\
\\end{document} \n";

typedef struct {
  char **comments;
  int show_dif,latex;
  } validate_t;
/*-----------------------------------------------------------------------

-----------------------------------------------------------------------*/


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int decode_atlasname(char *atlas_directory,char *atlas_convention,char *wave, char **filename)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=tide_decode_atlasname(atlas_directory,atlas_convention,wave, 0,filename,1);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int read_atlas(char *filename, vector<mgr_t> mgr, int nmgr,double *a,double *G, double mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, j, k, l, m, n, items,node, status;
  int kk,ll,mm;
  int nextk[24]={-1,+1, 0, 0,-1,+1,-1,+1,-2,+2, 0, 0,-2,+2,-2,+2,-1,-1,+1,+1,-2,-2,+2,+2};
  int nextl[24]={ 0, 0,-1,+1,-1,+1,+1,-1, 0, 0,-2,+2,-1,-1,+1,+1,-2,+2,-2,+2,-2,+2,-2,+2};
  float zr, zi;
  float *amp, *pha, *zx, *zy, factor;
  float *buf[2],spec[2];
  fcomplex zz,*z,cmask;
  fcomplex *tide;
  FILE *in;
  float dummy;
  double x, y, dum;
  char units[256];
  char *sunit;
  size_t count[3], start[3], var;
  size_t length;
  int ncid, amp_id, pha_id;
  nc_type xtypep;

  grid_t grid;

  cdfgbl_t global;
  int id;
  variable_t varinfo;

  start[0] = 0;
  start[1] = 0;
  start[2] = 0;

  count[0] = 1;
  count[1] = 1;
  count[2] = 1;

  status= cdf_globalinfo(filename,&global,0);
  if(status !=0) goto error;
  id=cdf_identify(global,"Ha");
  status=cdf_loadvargrid_2d (filename,id, &grid);
  if(status !=0) goto error;
  buf[0] =(float *) malloc(grid.nx*grid.ny*sizeof(float));
  buf[1] =(float *) malloc(grid.nx*grid.ny*sizeof(float));
  tide   =(fcomplex *) malloc(grid.nx*grid.ny*sizeof(fcomplex));
/*-----------------------------------------------------------------------
  load netcdf variable */
  id=cdf_identify(global,"Ha");
  status= cdf_loadvar_r1_2d (filename, id, 0, 0, grid, grid.nx, buf[0], &spec[0] ,&varinfo);
  id=cdf_identify(global,"Hg");
  status= cdf_loadvar_r1_2d (filename, id, 0, 0, grid, grid.nx, buf[1], &spec[1] ,&varinfo);
  cmask=fcomplex(9999.,9999.);
  for (j=0;j<grid.ny;j++)
    for (i=0;i<grid.nx;i++) {
      n=i+grid.nx*j;
      if((buf[0][n]!=spec[0]) && (buf[1][n]!=spec[1])) {
        tide[n]=fcomplex(buf[0][n]*cos(buf[1][n]*d2r),-buf[0][n]*sin(buf[1][n]*d2r));
        }
      else {
        tide[n]=cmask;
        }
      }
  free(buf[0]);
  free(buf[1]);

    for(n = 0; n < nmgr; n++) {
      if(strcmp(mgr[n].loc.units,"degrees")==0) {
        x = mgr[n].loc.lon;
        y = mgr[n].loc.lat;
        }
      else {
        x = mgr[n].loc.lon * r2d;
        y = mgr[n].loc.lat * r2d;
        }
//      y = gLGP1data[0][node].lat * r2d;
      if(x < grid.xmin - grid.dx / 2.)
        x = x + 360.0;
      if(x > grid.xmax + grid.dx / 2.)
        x = x - 360.0;
      status = map_interpolation(grid, tide, cmask, x, y,&zz);
      if(status != 0) {
        status=map_index( grid,  x,  y, &k, &l);
        if(status == 0) {
          mm=0;
          do {
            kk=k+nextk[mm];
            ll=l+nextl[mm];
            mm++;
            status=-1;
            if(kk<0) continue;
            if(ll<0) continue;
            if(kk>grid.nx-1) continue;
            if(ll>grid.ny-1) continue;
            x=map_grid_x(grid,kk,ll);
            y=map_grid_y(grid,kk,ll);
            status = map_interpolation(grid, tide, cmask, x, y,&zz);
            } while((status!=0)&&(mm<24));
          mm=0;
          while((status!=0)&&(mm<24)) {
            kk=k+2*nextk[mm];
            ll=l+2*nextl[mm];
            mm++;
            status=-1;
            if(kk<0) continue;
            if(ll<0) continue;
            if(kk>grid.nx-1) continue;
            if(ll>grid.ny-1) continue;
            x=map_grid_x(grid,kk,ll);
            y=map_grid_y(grid,kk,ll);
            status = map_interpolation(grid, tide, cmask, x, y,&zz);
            }
          }
        }
      if(status != 0) {
        printf("interpolation error: node %d lon=%lf lat=%lf \n", n, x, y);
        a[n] =  mask;
        G[n] =  mask;
        }
      else {
        a[n] =  abs(zz);
/*----------------------------------------------------------------------
      convert into phase lag*/
        G[n] = -arg(zz)* r2d;
        }
      }

  grid.free();
  
  return (status);

error:
  map_printgrid(grid);
  printf("interpolation error: could not read %s\n", filename);
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int mgr_penalty(vector<mgr_t> mgr,int nmgr, double **covariance, int **neighbours, int *card, char *wave, double **a,double **G)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,j,k,l,m,n,status;
  hypermatrix_t Cd;
  double *misfits[2],*x[2],*ratio[2];
  double zi,zr;
  double a1,p1,a2,p2,d,da,dG;
  double penalty,trace,penalty_r,penalty_i;
  statistic_t s[2];
  FILE *out;
  char filename[256];

   for (n=0;n<nmgr;n++) {
//     card[n]=MIN(5,card[n]);
//     card[n]--;
//      printf("%d %lf\n",n,covariance[n][card[n]-1]);
//      covariance[n][card[n]-1]=1.e-06;
     }

  status=syspack_init(covariance,neighbours,card,nmgr, &Cd);
  if(status!=0) goto error;

  misfits[0]=new double[nmgr];
  misfits[1]=new double[nmgr];
  x[0]=new double[nmgr];
  x[1]=new double[nmgr];
  ratio[0]=new double[nmgr];
  ratio[1]=new double[nmgr];

  for (n=0;n<nmgr;n++) {
    a1=mgr[n].data[0].amp;
    p1=mgr[n].data[0].phi;
    a2=a[0][n];
    p2=G[0][n];
    zr=a2*cos(-p2*d2r)-a1*cos(-p1*d2r);
    zi=a2*sin(-p2*d2r)-a1*sin(-p1*d2r);
    misfits[0][n]=zr;
    misfits[1][n]=zi;
    }

  status=syspack_inverse(misfits[0],x[0],nmgr,Cd);
  if(status!=0) goto error;
  status=syspack_inverse(misfits[1],x[1],nmgr,Cd);
  if(status!=0) goto error;

  penalty=0;
  penalty_r=0;
  penalty_i=0;
  for (n=0;n<nmgr;n++) {
    penalty+=x[0][n]*misfits[0][n];
    penalty+=x[1][n]*misfits[1][n];
    penalty_r+=x[0][n]*misfits[0][n];
    penalty_i+=x[1][n]*misfits[1][n];
    ratio[0][n]=x[0][n]*misfits[0][n]+x[1][n]*misfits[1][n];
    }

  penalty/=2.*nmgr;

  trace=0;
  for (n=0;n<nmgr;n++) {
    trace+=misfits[0][n]*misfits[0][n]/covariance[n][0];
    trace+=misfits[1][n]*misfits[1][n]/covariance[n][0];
    ratio[1][n]=(misfits[0][n]*misfits[0][n]+misfits[1][n]*misfits[1][n])/covariance[n][0];
    }

  trace/=2.*nmgr;

  s[0]=get_statistics(ratio[0],nmgr);
  s[1]=get_statistics(ratio[1],nmgr);

  printf("wave %s, penalty=%lf trace=%lf \n",wave,penalty,sqrt(trace));

  sprintf(filename,"penalty.%s",wave);
  out=fopen(filename,"w");
  fprintf(out,"wave %s, penalty=%lf trace=%lf \n",wave,penalty,sqrt(trace));
  fclose(out);

  delete[] misfits[0];
  delete[] misfits[1];
  delete[] x[0];
  delete[] x[1];
  delete[] ratio[0];
  delete[] ratio[1];

  return(0);

error:
  return(-1);

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int mgr_compare(char *output, vector<mgr_t> mgr,int nmgr, int *list, char **wave, int nwave, double **a,double **G, char **comments, validate_t options)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double mask=999.9;
  double zi,zr,mean[10],rms[10],count;
  double a1,p1,a2,p2,d,da,dG;
  int i,j,k,l,m,n,fmt,status;
  FILE *out;
  char *ordering=NULL;
  fcomplex meanc;
  char c;
  fcomplex cmisfit;
  char *tab,*notab,*newline;
  char show_diff=0,latex=0;
  int *keep;

  float hlim=15000.;
  float factor=10.0;

  latex=options.latex;
  show_diff=options.show_dif;

  if(latex) {
    tab  =strdup(" & ");
    notab=strdup("   ");
    newline=strdup("\\\\");
    }
  else {
    tab  =strdup(" ");
    notab=strdup(" ");
    newline=strdup("");
    }

  out=fopen(output,"w");
  if(latex) {
    fprintf(out,"%s",header);
    }

  fprintf(out,"%s",comments[0]);
  if(latex) fprintf(out,"%s",newline);
  fprintf(out,"\n");

  fprintf(out,"%s",comments[1]);
  if(latex) fprintf(out,"%s",newline);
  fprintf(out,"\n");

  if(latex) {
    for(k=0;k<3;k++) fprintf(out,"~ \\\\ \n");
    fprintf(out,"\\begin{supertabular}{|c|c|cc|");
    for (k=0;k<nwave;k++) {
      fprintf(out,"cc");
      if(show_diff) {
        fprintf(out,"ccc");
        }
      fprintf(out,"|");
      }
    fprintf(out,"}\n\\hline\n");
    }

  fprintf(out,"%4s%s%25s%s%15s%s%15s","No",tab,"Station",tab,"longitude",tab,"latitude");
  for (k=0;k<nwave;k++) {
    if(show_diff) {
      fprintf(out,"%s",tab);
      if(latex) fprintf(out,"\\multicolumn{5}{c|} {");
      fprintf(out,"%10s%4s%s",wave[k]," ",notab);
      if(latex) fprintf(out,"}");
      }
    else {
      fprintf(out,"%s",tab);
      if(latex) fprintf(out,"\\multicolumn{2}{c|} {");
      fprintf(out,"%10s%4s%s",wave[k]," ",notab);
      if(latex) fprintf(out,"}");
      }
    }
  fprintf(out,"%s\n",newline);

  fprintf(out,"%s%s%s%59s",tab,tab,tab," ");
  for (k=0;k<nwave;k++) {
    if(show_diff) {
      fprintf(out,"%s%7s%s%7s%s%3s%s%3s%s%6s",tab,"A   ",tab,"G   ",tab,"dA",tab,"dG",tab,"dE");
      }
    else {
      fprintf(out,"%s%7s%s%7s",tab,"A   ",tab,"G   ");
      }
    }
  fprintf(out,"%s\n",newline);
  if(latex) fprintf(out,"\\hline\n");
  fprintf(out,"\n");


  keep=new int[nmgr];

  for (m=0;m<nmgr;m++) {
    string slon,slat;
    n=list[m];
    keep[n]=1;
    slon=dms(mgr[n].loc.lon);
    slat=dms(mgr[n].loc.lat);
    fprintf(out,"%4d%s%32s%s%15s%s%15s",n+1,tab,check_string(mgr[n].name),tab,slon.c_str(),tab,slat.c_str());
    for (k=0;k<nwave;k++) {
      for(i=0;i<mgr[n].nwave;i++) {
        if(strcmp(wave[k],mgr[n].data[i].constituent.name)==0) {
          break;
          }
        }
      if(i==mgr[n].nwave) {
        fprintf(out,"%s---/---%s---/---",tab,tab);
        if(show_diff) {
          fprintf(out,"%s---%s---%s------",tab,tab,tab);
          }
        }
      else {
        a1=100*mgr[n].data[i].amp;
        p1=mgr[n].data[i].phi;
/* *----------------------------------------------------------------------------
        Profondeur limite*/
        if(mgr[n].loc.depth>hlim) {
//          a[k][n]=mask;
          keep[n]=0;
          }
/* *----------------------------------------------------------------------------
        Errrur limite*/
        if(100*mgr[n].data[i].error/a1>0.15) {
//          a[k][n]=mask;
          keep[n]=0;
          };
        if(keep[n]==0) {
          fprintf(out,"%s%3d/---%s%3d/---",
            tab,(int) NINT(a1),tab,(int) NINT(p1));
          if(show_diff) {
            fprintf(out,"%s---%s---%s------",tab,tab,tab);
            }
          }
        else {
          a2=100*a[k][n];
          p2=G[k][n];
          if(p2<  0) p2+=360.0;
          if(p2>360) p2-=360.0;
          zr=a2*cos(-p2*d2r)-a1*cos(-p1*d2r);
          zi=a2*sin(-p2*d2r)-a1*sin(-p1*d2r);
          d=sqrt(zr*zr+zi*zi);
          da=a2-a1;
          dG=p2-p1;
          if(dG>+180.) dG-=360.0;
          if(dG<-180.) dG+=360.0;
          fprintf(out,"%s%3d/%3d%s%3d/%3d%",
            tab,(int) NINT(a1),(int) NINT(a2),tab,(int) NINT(p1),(int) NINT(p2));
          if(show_diff) {
            fprintf(out,"%s%3d%s%3d%s%6.1f",
              tab,(int) NINT(da),tab,(int) NINT(dG),tab,d);
            }
          }
        }
      }
    fprintf(out,"%s\n",newline);
    }

  if(latex) fprintf(out,"\\hline\n\\end{supertabular}\n");

  fprintf(out,"\n");
  if(latex) {
    for(k=0;k<3;k++) fprintf(out,"~ \\\\ \n");
    fprintf(out,"\\begin{supertabular}{|c|cc|cc|ccc|c|} \n");
    fprintf(out,"\\hline\n");
    }
  fprintf(out,"%10s%s","wave",tab);
  if(latex) {
    fprintf(out,"\\multicolumn{2}{c|} {$\\Delta_a$} & \\multicolumn{2}{c|} {$\\Delta_G$} & \\multicolumn{3}{c|} {$e$} & N \\\\ \n");
    }
  else {
//    fprintf(out,"%10s%s%7s%s%7s%s%7s%s%7s%s%6s \\\\ \n","wave",tab,"A   ",tab,"G   ",tab,"E",tab,"N"); LR, comment: 2008/04/08
    fprintf(out,"%10s%s%7s%s%7s%s%7s%s%7s \\\\ \n","wave",tab,"A   ",tab,"G   ",tab,"E",tab,"N"); //LR, change: 2008/04/08
    }
  fprintf(out,"%10s%s%7s%s%7s%s%7s%7s%s%7s%s%7s%s%7s%s%6s \\\\ \n"," ",tab,"mean",tab,"rms",tab,"mean",tab,"rms",tab,"mean",tab,"rms",tab,"rms*",tab,"");
  if(latex) fprintf(out,"\\hline\n");
  for (k=0;k<nwave;k++) {
    for(l=0;l<5;l++) {
      mean[l]=0;
      rms[l]=0;
      }
    count=0;
    for (n=0;n<nmgr;n++) {
      i=mgr[n].wave_index(wave[k]);
      if(i==-1) continue;
      a1=100*mgr[n].data[i].amp;
      p1=mgr[n].data[i].phi;
      if(keep[n]==0) continue;
      a2=100*a[k][n];
      p2=G[k][n];
      if(p2<  0) p2+=360.0;
      if(p2>360) p2-=360.0;
      zr=a2*cos(-p2*d2r)-a1*cos(-p1*d2r);
      zi=a2*sin(-p2*d2r)-a1*sin(-p1*d2r);
      d=sqrt(zr*zr+zi*zi);
      da=a2-a1;
      dG=p2-p1;
      cmisfit=fcomplex(zr,zi);
      count++;
      l=0;
      d=sqrt((zr*zr+zi*zi)/2.);
      mean[l]=mean[l]+d;
      mean[l]=0.0;
      rms[l] =rms[l]+d*d;
      l++;
      d=a1-a2;
      mean[l]=mean[l]+d;
      rms[l] =rms[l]+d*d;
      l++;
      d=p1-p2;
      if(d >  180.) d=d-360;
      if(d < -180.) d=d+360;
      mean[l]=mean[l]+d;
      rms[l] =rms[l]+d*d;
      l++;
      d=zr;
      mean[l]=mean[l]+d;
      rms[l] =rms[l]+d*d;
      l++;
      d=zi;
      mean[l]=mean[l]+d;
      rms[l] =rms[l]+d*d;
      }
    for(l=0;l<5;l++) {
      mean[l]=mean[l]/count;
      rms[l]=sqrt(rms[l]/count - mean[l]*mean[l]);
      }
    fprintf(out,"%10s",wave[k]);
/*----------------------------------------------------------------------
    amplitude mean misfit/rms*/
    fprintf(out,"%s%3d %s %3d",tab,(int) NINT(factor*mean[1]),tab,(int) NINT(factor*rms[1]));
/*----------------------------------------------------------------------
    phase lag mean misfit/rms*/
    fprintf(out,"%s%3d %s %3d",tab,(int) NINT(mean[2]),tab,(int) NINT(rms[2]));
/*----------------------------------------------------------------------
    complex error mean misfit (module) /rms*/
    fprintf(out,"%s%3d %s %3d",tab,(int) NINT(factor*sqrt(mean[3]*mean[3]+mean[4]*mean[4])/sqrt(2.0)),tab,(int) NINT(factor*sqrt(rms[3]*rms[3]+rms[4]*rms[4])/sqrt(2.0)));
/*----------------------------------------------------------------------
    complex error deviation from 0*/
    fprintf(out,"%s%3d%s%3d/%3d",tab,(int) NINT(factor*rms[0]),tab,(int) count,nmgr);
    fprintf(out,"\\\\\n");
    }

  fprintf(out,"\n");
  if(latex) {
    fprintf(out,"%s",trailer);
    }

  fclose(out);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double mask=999.9;
  double tau;
  double *residual,*buffer;

  double  **a,**G;
  int i,j,k,l,m,n,fmt,status;

  FILE *in,*out;
  char *report,*comments[10],tmp[256], datafile[256];
  char *model=NULL,*data, *input,*output=NULL, *keyword, *path=NULL,*mgrfile=NULL;
  char *ordering=NULL;

  spectrum_t spectrum,reference_spectrum;
  char *wave[100];
  int nwave=0;
  date_t date;
  vector<mgr_t> mgr;
  int nmgr;
  char *atlas_directory=NULL,*atlas_convention=NULL;

  char show_diff=0,latex=0;
  int *list;

  char *bathymetry=NULL;
  grid_t topogrid;
  float *topo,topomask,h,hlim=15000.;
  double x,y;
  float factor=10.0;
  double **covariance;
  int **neighbours, *card;

  validate_t options;

  fct_echo( argc, argv);

  n=1;
  while (n < argc) {
    keyword=strdup(argv[n]);
    switch (keyword[0]) {
      case '-':
        switch (keyword[1])
        {

/* *----------------------------------------------------------------------
        bathymetry database*/
        case 'b' :
          bathymetry= strdup(argv[n+1]);
          n++;
          n++;
          break;

/* *----------------------------------------------------------------------
        print data/models differences*/
        case 'd' :
          show_diff=1;
          n++;
          break;

/* *----------------------------------------------------------------------
        output format*/
        case 'l' :
          latex=1;
          n++;
          break;

/* *----------------------------------------------------------------------
        naming convention for tidal atlases*/
        case 'm' :
          atlas_convention= strdup(argv[n+1]);
          n++;
          n++;
          break;

/* *----------------------------------------------------------------------
        path for tidal atlases*/
        case 'p' :
          path= strdup(argv[n+1]);
          n++;
          n++;
          break;

/* *----------------------------------------------------------------------
        output file name*/
        case 'o' :
          output= strdup(argv[n+1]);
          n++;
          n++;
          break;

/* *----------------------------------------------------------------------
        harmonic data input file*/
        case 'c' :
          mgrfile= strdup(argv[n+1]);
          n++;
          n++;
          break;

/* *----------------------------------------------------------------------
        output order: FILE, ALPHA, LAT, LON, IMMERSION*/
        case 'x' :
          ordering= strdup(argv[n+1]);
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
/* *----------------------------------------------------------------------
          tidal wave list*/
          wave[nwave]= strdup(argv[n]);
          spectrum.n=nwave+1;
          nwave++;
          n++;
        break;
      }
      free(keyword);
    }

  if(path!=NULL) atlas_directory=strdup(path);

  options.latex=latex;
  options.show_dif=show_diff;

  date.year=1950;
  date.month=1;
  date.day=1;
  date.second=0.0;
  astro_angles_t astro_angles;
  reference_spectrum = initialize_tide(&astro_angles,date);

  spectrum.nmax=nwave+1;
  spectrum.n=0;
  spectrum.waves=new tidal_wave[spectrum.nmax];

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
//    printf ("wave: %10s, pulsation: %12.6f degrees/h \n", spectrum.waves[i].name,spectrum.waves[i].omega);
    }

  for (i=0; i<spectrum.n; i++) {
    for (j=i+1; j<spectrum.n; j++) {
      tau=deltaOmega2separation(spectrum.waves[j].omega-spectrum.waves[i].omega);
      printf ("wave: %10s %10s, separation: %9.3f days \n",
        spectrum.waves[i].name,spectrum.waves[j].name,tau);
      }
    }

/* *-----------------------------------------------------------------------------
  read tide gauge database */
  strcpy(datafile,mgrfile);
  status= mgr_loadobs(datafile, wave[0], &nmgr, mgr, &covariance, &neighbours, &card);
  if(nmgr==-1) {
    printf("error in reading %s\n",datafile);
    goto error;
    }

  if(ordering==NULL) {
    ordering=strdup("ALPHA");
    }

  list= mgr_order(mgr, nmgr, ordering, 1);

  a=(double **) malloc(nwave*sizeof(double));
  G=(double **) malloc(nwave*sizeof(double));

  for (k=0;k<nwave;k++) {
    a[k]=(double *) malloc(nmgr*sizeof(double));
    G[k]=(double *) malloc(nmgr*sizeof(double));
    }

/* *-----------------------------------------------------------------------------
  read tidal atlas database */
  for (k=0;k<nwave;k++) {
    status=decode_atlasname(atlas_directory,atlas_convention,wave[k],&model);
    printf("validate: treating %s wave from %s \n",wave[k],model);
    status=read_atlas(model,mgr,nmgr,a[k],G[k],mask);
    }

/* *-----------------------------------------------------------------------------
  read topo database */
  if(bathymetry!=NULL) {
    status=grd_loadgrid(bathymetry,&topogrid);
    if(status !=0) {
      __OUT_BASE_LINE__("cannot load bathymetry file=%f\n",bathymetry);
      exit(-1);
      }
    topo= (float *) malloc(topogrid.nx*topogrid.ny*sizeof(float));
    topogrid.modeH=0;
    status=  grd_loadr1(bathymetry,topogrid,topo,&topomask);
    }
  else {
    topo=0;
    }

  for (m=0;m<nmgr;m++) {
    n=list[m];
    if(topo!=0) {
      x=mgr[n].loc.lon;
      y=mgr[n].loc.lat;
      x=map_recale(topogrid,x);
      status=map_interpolation(topogrid, topo,topomask,x,y,&h);
      mgr[n].loc.depth=h;
      }
    }

  if(output==NULL) {
    report=strdup("validate.out");
    }
  else {
    report=strdup(output);
    }

  comments[0]=new char[1024];
  sprintf(comments[0],"model    : %s %s",atlas_directory,atlas_convention);
  comments[1]=new char[1024];
  sprintf(comments[1],"database : %s",clean_string(datafile));

  status=mgr_compare(output,mgr,nmgr,list,wave,nwave,a,G,comments,options);
  status= mgr_penalty(mgr, nmgr, covariance, neighbours, card, wave[0], a, G);

  for (k=0;k<nwave;k++) {
    free(a[k]);
    free(G[k]);
    }
  free(a);
  free(G);

//  __ERR_BASE_LINE__("%s -computation complete ^^^^^^^^^\n",argv[0]);
  exit(0);
error:
  __OUT_BASE_LINE__("validate: error detected, quit ... \n");
  exit(-1);
}
