
/**************************************************************************

  T-UGO tools, 2006-2011

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
\brief Filter spatially model outputs

<!-- A LINK TO main() or print_help() WILL NOT LINK TO THE RIGHT SOURCE ! -->
See the main function for how this works
and the print_help function for how to use this.
You may add more description, but you can also
- document the main function
- or update print_help().

Notes for later :
<h3>Diagram of a biquad filer</h3>
From http://en.wikipedia.org/wiki/Digital_biquad_filter#Direct_Form_2
\dot
digraph biquadFiler {
  rankdir=LR;
  splines=ortho;
  I -> s1 -> B -> s2 -> O;
  I [shape=box];
  B [shape=box];
  O [shape=box];
  {rank=same; {a1 a2}->s1}
  {rank=same; B ->z1->z2}
  {rank=same; {b1 b2}->s2}
  s1 [label="+"];
  s2 [label="+"];
  a1 [label="*a1"];
  a2 [label="*a2"];
  b1 [label="*b1"];
  b2 [label="*b2"];
  z1 [label="1 step delay"];
  z2 [label="1 step delay"];
  z1 -> a1;z1 -> b1;
  z2 -> a2;z2 -> b2;
  }
\enddot
This gives the transfer functions of \f$I\f$, \f$B\f$ and \f$O\f$, with \f$z=e^{j\omega s^{-1}}\f$ the transfer function of a 1 step advance with \f$s\f$ the sampling frequency :
\f{eqnarray*}{
F(B)&=&F(I)+a_1 z^{-1} F(B)+a_2 z^{-2} F(B) \\
F(O)&=&F(B)+b_1 z^{-1} F(B)+b_2 z^{-2} F(B) \\
F(O)&=&F(I)\frac{1+b_1 z^{-1}+b_2 z^{-2}}{1-a_1 z^{-1}-a_2 z^{-2}}
\f}
*/
/*----------------------------------------------------------------------------*/

#define MAIN_SOURCE

#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <vector>
#include <string>

#include "tools-structures.h"
#include "tools-define.h"

#include "tides.h"
#include "functions.h"
#include "netcdf-proto.h"
#include "poc-time.h"
#include "filter.h"
#include "geo.h"
#include "map.h"

#define GLOBAL  0
#define YEARLY  1
#define MONTHLY 2
#define WEEKLY  3
#define DAILY   4
#define HOURLY  5
#define CUSTOM  6

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
//unused function
  int decode_name(date_t actual, const char *varname, const char *name_template, char **filename, int *packing)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int  status, hour;
  char dummy[256], *pointer;
//  date_t actual;
  date_t cdf_reference;
  FILE *out;

//  getmodeldate(t, &actual);

/* *----------------------------------------------------------------------------
  Reconstruct the input file name*/
  (*filename) = new char[strlen(name_template) + 256];
  sprintf((*filename), "%s", name_template);

  out = fopen(*filename, "r");
  status = (out == NULL);

  switch (status) {
    case 0:
/* *----------------------------------------------------------------------------
      file exists, do nothing more*/
      fclose(out);
      *packing=GLOBAL;
      break;

    default:
/*-----------------------------------------------------------------------------
      use format convention information*/
/* *----------------------------------------------------------------------------
      substitute YYYY with current year*/
      pointer = strstr((*filename), "YYYY");
      if(pointer != NULL) {
        sprintf(dummy, "%4d", actual.year);
        strncpy(pointer, dummy, 4);
        *packing=YEARLY;
        }
/* *----------------------------------------------------------------------------
      substitute MM with current month*/
      pointer = strstr((*filename), "MM");
      if(pointer != NULL) {
        sprintf(dummy, "%2.2d", actual.month);
        strncpy(pointer, dummy, 2);
        *packing=MONTHLY;
        }
/* *----------------------------------------------------------------------------
      substitute DD with current day*/
      pointer = strstr((*filename), "DD");
      if(pointer != NULL) {
        sprintf(dummy, "%2.2d", actual.day);
        strncpy(pointer, dummy, 2);
        *packing=DAILY;
        }
/* *----------------------------------------------------------------------------
      substitute HH with current hour*/
      pointer = strstr((*filename), "HH");
      if(pointer != NULL) {
        sprintf(dummy, "%2.2d", NINT(actual.second/3600.));
        strncpy(pointer, dummy, 2);
        *packing=HOURLY;
        }
/* *----------------------------------------------------------------------------
      lookup file fitting ???? with minute and seconds index*/
      pointer = strstr((*filename), "????");
      if(pointer != NULL) {
        int i,j;
        for(j=0;j<60;j++) {
          for(i=0;i<60;i++) {
            sprintf(dummy, "%2.2d%2.2d", j,i);
            strncpy(pointer, dummy, 4);
            status= get_file_size(*filename,0);
/* *----------------------------------------------------------------------------
            status=0 if file found*/
            if(status==0) break;
            }
          if(status==0) break;
          }
        }
/* *----------------------------------------------------------------------------
      substitute VARNAME*/
      pointer = strstr((*filename), "VARNAME");
      if(pointer != NULL) {
        sprintf(dummy, "%s", varname);
        strncpy(pointer, dummy, strlen(varname));
        }
      break;
    }

  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

//unused function
  vector<string> build_filelist(const char **varnames, int nvars, const char *convention, date_t start,date_t final, spectrum_t spectrum)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double t,*time,duration;
  int packing,count;
  int k,n;
  int status;
  int day,month,year,frame;
  char *filename[3],*s;
  grid_t grid;
  cdfgbl_t grid_info[10];
  cdfvar_t info[10],detided[10],variable;
  date_t actual;
  int nvalues,nrhs,nframes,n_spacedims;
  decoded_t decoded;
  vector<string> filelist;

/* *-------------------------------------------------------------------------------------
  record length */
  duration=ellapsed_time(start,  final, 's');

/* *----------------------------------------------------------------------------
  Get variable dimension once and for all (assumes variables have same size !!!) */
  status=decode_name(start, "T", convention, &(filename[0]),&packing);
  status= cdf_globalinfo(filename[0],&grid_info[0],0);
  status= cdf_varinfo  (filename[0],varnames[0],&(info[0]),0);

/* *----------------------------------------------------------------------------
  Get grid once and for all (assumes variables have same discretisation !!!) */
  status= poc_getgrid2d (filename[0], grid_info[0], info[0], &grid);

/* *----------------------------------------------------------------------------
  analyse netcdf informations */
  status= poc_decode_axis(info[0], grid_info[0], &decoded);

/* *----------------------------------------------------------------------------
  retain only space dimensions */
  nvalues=1;
  n_spacedims=0;
  switch (decoded.xlen) {
    case 0:
      break;
    default:
      nvalues*=decoded.xlen;
      n_spacedims++;
      break;
    }
  switch (decoded.ylen) {
    case 0:
      break;
    default:
      nvalues*=decoded.ylen;
      n_spacedims++;
      break;
    }
  switch (decoded.zlen) {
    case 0:
      break;
    default:
      nvalues*=decoded.zlen;
      n_spacedims++;
      break;
    }

  nrhs=nvars;
  nframes=decoded.tlen;

  return(filelist);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int parse_filelist(vector<string> filelist, char **varnames, int nvars, date_t start,date_t final)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double t,*time,duration,tmin,tmax;
  int packing,count;
  int k,n;
  int status;
  int day,month,year,frame,nframes;
  grid_t grid;
  cdfgbl_t global[10];
  cdfvar_t info[10],detided[10],variable;
  date_t actual,origine;
  decoded_t decoded;
  vector<int> nframelist;
 
  tmin=+1.e+10;
  tmax=-1.e+10;
  
  for(vector<string>::const_iterator p = filelist.begin(); p != filelist.end(); p++){
/* *----------------------------------------------------------------------------
    Get variable dimension once and for all (assumes variables have same size !!!) */
    status= cdf_globalinfo (p->c_str(),&global[0],0);
    status= cdf_varinfo  (p->c_str(),varnames[0],&(info[0]),0);
    if(status!=0) return(-1);
/* *----------------------------------------------------------------------------
    analyse netcdf informations */
    status= poc_decode_axis(info[0], global[0], &decoded);
#warning this function does not seam to do anything usefull
    }

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int get_cartesian(grid_t grid, grid_t *cgrid, grid_t *sgrid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int    i,j,k,l,n,status;
  double t,p,x,y;
  projPJ projection;
  double ref_lat[2],ref_lon[2];
  int ptype=0;

  ref_lat[0]=(grid.ymax+grid.ymin)/2;
  ref_lat[0]*=M_PI/180.0;

  projection=assign_projection(ptype,ref_lat, ref_lon);

/* *----------------------------------------------------------------------------
  cgrid is cartesian, regular grid need for filtering */
  cgrid->nx=grid.nx;
  cgrid->ny=grid.ny;
  cgrid->nz=grid.nz;
  cgrid->modeH=0;

  cgrid->x=new double[cgrid->nx];
  cgrid->y=new double[cgrid->ny];

  cgrid->ymin= +1.e+35;
  cgrid->ymax= -1.e+35;
  cgrid->xmin= +1.e+35;
  cgrid->xmax= -1.e+35;

/* *----------------------------------------------------------------------------
  sgrid is spherical counter-part ofgcrid need for interpolation */
  sgrid->nx=grid.nx;
  sgrid->ny=grid.ny;
  sgrid->nz=grid.nz;
  sgrid->modeH=2;

  sgrid->x=new double[sgrid->nx*sgrid->ny];
  sgrid->y=new double[sgrid->nx*sgrid->ny];

  sgrid->ymin= +1.e+35;
  sgrid->ymax= -1.e+35;
  sgrid->xmin= +1.e+35;
  sgrid->xmax= -1.e+35;

  for(j=0;j<cgrid->ny;j++) {
    for(i=0;i<cgrid->nx;i++) {
      k=cgrid->nx*j+i;
      t=grid.x[k];
      p=grid.y[k];
      geo_to_projection(projection,p,t,&x,&y);
      cgrid->ymin=MIN(cgrid->ymin,y);
      cgrid->ymax=MAX(cgrid->ymax,y);
      cgrid->xmin=MIN(cgrid->xmin,x);
      cgrid->xmax=MAX(cgrid->xmax,x);
      }
    }
    
  cgrid->dx=(cgrid->xmax-cgrid->xmin)/(cgrid->nx-1.);
  cgrid->dy=(cgrid->ymax-cgrid->ymin)/(cgrid->ny-1.);

  for(j=0;j<cgrid->ny;j++) {
    for(i=0;i<cgrid->nx;i++) {
      k=cgrid->nx*j+i;
      x=cgrid->xmin+i*cgrid->dx;
      y=cgrid->ymin+j*cgrid->dy;
      cgrid->x[i]=x;
      cgrid->y[j]=y;
      projection_to_geo(projection,&p,&t,x,y);
      sgrid->x[k]=t;
      sgrid->y[k]=p;
      sgrid->ymin=MIN(sgrid->ymin,p);
      sgrid->ymax=MAX(sgrid->ymax,p);
      sgrid->xmin=MIN(sgrid->xmin,t);
      sgrid->xmax=MAX(sgrid->xmax,t);
      }
    }

  for(j=0;j<cgrid->ny;j++) {
    for(i=0;i<cgrid->nx;i++) {
      k=cgrid->nx*j+i;
      }
    }
    
  sgrid->dx=(sgrid->xmax-sgrid->xmin)/(sgrid->nx-1.);
  sgrid->dy=(sgrid->ymax-sgrid->ymin)/(sgrid->ny-1.);
  
  cgrid->projection=projection;
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int filter2D_complex(grid_t grid, grid_t cgrid,grid_t sgrid, float **buffer, float fmask, float *weight, float scale, int nvalues)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,n;
  int status;
  decoded_t decoded[10];

/* *-------------------------------------------------------------------------------------
  convert to real/imag */
  for(n=0;n<nvalues;n++) {
    if((buffer[0][n]!=fmask)&&(buffer[1][n]!=fmask)) {
      buffer[2][n]=buffer[0][n]*cos(buffer[1][n]*d2r);
      buffer[3][n]=buffer[0][n]*sin(buffer[1][n]*d2r);
      }
    else {
      buffer[2][n]=fmask;
      buffer[3][n]=fmask;
      }
    }
    
/* *-------------------------------------------------------------------------------------
  interpolate on evenly spaced grid */
  status=map_export(grid,buffer[2],fmask,sgrid,buffer[4],fmask,0);
  status=map_export(grid,buffer[3],fmask,sgrid,buffer[5],fmask,0);
/* *-------------------------------------------------------------------------------------
  filter */
  buffer[6]=Loess2D_BF(cgrid, buffer[4], fmask, scale, weight, buffer[6]);
  buffer[7]=Loess2D_BF(cgrid, buffer[5], fmask, scale, weight, buffer[7]);
/* *-------------------------------------------------------------------------------------
  back to native grid */
  status=map_export(sgrid,buffer[6],fmask,grid,buffer[8],fmask,0);
  status=map_export(sgrid,buffer[7],fmask,grid,buffer[9],fmask,0);
/* *-------------------------------------------------------------------------------------
  HF */
  for(n=0;n<nvalues;n++) {
    if((buffer[2][n]!=fmask)&&(buffer[8][n]!=fmask)) {
      buffer[6][n]=buffer[2][n]-buffer[8][n];
      buffer[7][n]=buffer[3][n]-buffer[9][n];
      }
    else {
      buffer[6][n]=fmask;
      buffer[7][n]=fmask;
      }
    }

/* *-------------------------------------------------------------------------------------
  LF, convert to amplitude/phase */
  for(n=0;n<nvalues;n++) {
    if((buffer[8][n]!=fmask)&&(buffer[9][n]!=fmask)) {
      buffer[2][n]=sqrt(buffer[8][n]*buffer[8][n]+buffer[9][n]*buffer[9][n]);
      buffer[3][n]=atan2(buffer[9][n],buffer[8][n])*r2d;
      }
    else {
      buffer[2][n]=fmask;
      buffer[3][n]=fmask;
      }
    }
    
/* *-------------------------------------------------------------------------------------
  HF, convert to amplitude/phase */
  for(n=0;n<nvalues;n++) {
    if((buffer[6][n]!=fmask)&&(buffer[7][n]!=fmask)) {
      buffer[4][n]=sqrt(buffer[6][n]*buffer[6][n]+buffer[7][n]*buffer[7][n]);
      buffer[5][n]=atan2(buffer[7][n],buffer[6][n])*r2d;
      }
    else {
      buffer[4][n]=fmask;
      buffer[5][n]=fmask;
      }
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int filter3D_complex(grid_t grid, grid_t cgrid,grid_t sgrid, float **buffer3D, float fmask, float *weight, float scale, int nvalues)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,n;
  int nvalues2D, status;
  decoded_t decoded[10];
  float *buffer2D[10];

  nvalues2D=grid.nx*grid.ny;
  
/* *-------------------------------------------------------------------------------------
  convert to 2D */
  for(n=0;n<grid.nz;n++) {
    for(k=0;k<10;k++) {
      buffer2D[k]=&(buffer3D[k][n*nvalues2D]);
      }
    printf("layer %d\n",n);
    filter2D_complex(grid, cgrid, sgrid, buffer2D, fmask, weight, scale, nvalues2D);
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

///unused function. See filter_01() instead.
  int filter_00(const char **varnames, int nvars, vector<string> filelist, char *gridfile, char *controlfile, date_t start,date_t final)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  float *buffer[10],fmask,*hmean=0;
  double t,*time=0,duration;
  double step;
  int packing,count;
  int k,n;
  int status;
  int day,month,year,frame;
  char *s;
  char test[256];
  string output;
  grid_t grid;
  cdfgbl_t global[10];
  cdfvar_t info[10];
  cdfvar_t filtered_LF[nvars],filtered_HF[nvars],raw[nvars],detided_vtime,vtime;
  decoded_t decoded[10];
  date_t actual,origine;
  bool *mask;
  int nvalues,nrhs,nframes,n_spacedims;
  grid_t cgrid, sgrid;
//float scale_x, float scale_y, float azimuth
  float scale, scale_x, scale_y, azimuth, *weight;

  status= cdf_globalinfo(filelist[0].c_str(),&global[0],0);
  for(k=0;k<nvars;k++) {
/* *----------------------------------------------------------------------------
    get variable dimension once and for all (assumes variables have same size !!!) */
    status= cdf_varinfo  (filelist[0].c_str(),varnames[k],&(info[k]),0);
    if(status==-1) {
      return(0);
      }

/* *----------------------------------------------------------------------------
    analyse netcdf informations */
    status= poc_decode_axis(info[k], global[0], &(decoded[k]));
    status= poc_decode_mask(info[k], &(decoded[k]));
    }
/* *----------------------------------------------------------------------------
  retain only space dimensions */
  nvalues=1;
  n_spacedims=0;
  switch (decoded[0].xlen) {
    case 0:
      break;
    default:
      nvalues*=decoded[0].xlen;
      n_spacedims++;
      break;
    }
  switch (decoded[0].ylen) {
    case 0:
      break;
    default:
      nvalues*=decoded[0].ylen;
      n_spacedims++;
      break;
    }
  switch (decoded[0].zlen) {
    case 0:
      break;
    default:
      nvalues*=decoded[0].zlen;
      n_spacedims++;
      break;
    }

  nrhs=nvars;
  nframes=decoded[0].tlen;
  if(nframes==0) nframes=1;
  
/*-----------------------------------------------------------------------------
  allocate memory for vectors */
  for(k=0;k<5;k++) {
    buffer[k]=new float[nvalues];
    if(buffer[k] == NULL) {
      printf("#memory allocation error for buffer N= %d \n",nvalues);
      check_error(-1, "memory allocation error", __LINE__, __FILE__, 1);;
      }
    }

  if(gridfile==0) gridfile=strdup(filelist[0].c_str());

  switch(n_spacedims) {
    case 2:
      status= poc_getgrid2d (gridfile, global[0], info[0], &grid);
      if(status!=0) {
        printf("cannot load 2D grid from gridfile %s\n",gridfile);
        }
      break;
    case 3:
/* *----------------------------------------------------------------------------
      get grid once and for all (assumes variables have same discretisation !!!) */
      int vdepth;
      status=poc_getgrid3d (gridfile,  global[0], info[0], &grid, &vdepth);
      if(status!=0) {
        printf("cannot load 3D grid from gridfile %s\n",gridfile);
        }
      break;
    default:
      status=-1;
      break;
    }
  status= get_cartesian(grid, &cgrid, &sgrid);

  scale=100000.;
  Loess2D_init(cgrid, scale, &weight);
  
/* *-------------------------------------------------------------------------------------
  construct harmonic matrix and RHS vectors */
  count=0;
  for(vector<string>::const_iterator p = filelist.begin(); p != filelist.end(); p++){
    count++;
/* *----------------------------------------------------------------------------
    get variable dimension once and for all (assumes variables have same size !!!) */
    status= cdf_varinfo  (p->c_str(),varnames[0],&(info[0]),0);
    status= poc_decode_axis(info[0], global[0], &(decoded[0]));
    nframes=decoded[0].tlen;
    if(nframes==0) nframes=1;
    status= poc_gettime(p->c_str(),info[0], global[0], &origine, &time, &nframes);
/* *-------------------------------------------------------------------------------------
    create/open the output files */
    output=string("filtered-");
    output+=*p;
    status= poc_createfile(output.c_str(),global[0]);
/* *-------------------------------------------------------------------------------------
    check filetered variables */
    for(k=0;k<nvars;k++) {
      sprintf(test,"%s",info[k].name);
      status= cdf_varinfo(output.c_str(),test,&(raw[k]),0);
      if(status==-1){
        raw[k]=info[k];
        printf("variable %s does not exist in %s, need to be created...\n",raw[k].name,output.c_str());
        status=create_ncvariable(output.c_str(), &(raw[k]));
        if(status!=NC_NOERR) return(-1);
        printf("done...\n");
        }
      sprintf(test,"%s-%s",info[k].name,"LF");
      status= cdf_varinfo(output.c_str(),test,&(filtered_LF[k]),0);
      if(status==-1){
        filtered_LF[k]=info[k];
        delete[] filtered_LF[k].name;
        filtered_LF[k].name=new char[strlen(info[k].name)+9];
        sprintf(filtered_LF[k].name,"%s-%s",info[k].name,"LF");
        printf("variable %s does not exist in %s, need to be created...\n",filtered_LF[k].name,output.c_str());
        status=create_ncvariable(output.c_str(), &(filtered_LF[k]));
        if(status!=NC_NOERR) return(-1);
        printf("done...\n");
        }
      sprintf(test,"%s-%s",info[k].name,"HF");
      status= cdf_varinfo(output.c_str(),test,&(filtered_HF[k]),0);
      if(status==-1){
        filtered_HF[k]=info[k];
        delete[] filtered_HF[k].name;
        filtered_HF[k].name=new char[strlen(info[k].name)+9];
        sprintf(filtered_HF[k].name,"%s-%s",info[k].name,"HF");
        printf("variable %s does not exist in %s, need to be created...\n",filtered_HF[k].name,output.c_str());
        status=create_ncvariable(output.c_str(), &(filtered_HF[k]));
        if(status!=NC_NOERR) return(-1);
        printf("done...\n");
        }
      }
/* *-------------------------------------------------------------------------------------
    read the output files */
    printf("read %s ...\n",p->c_str());
    for(frame=0;frame<nframes;frame++) {
      for(k=0;k<nvars;k++) {
        size_t *start, *count;
        start=new size_t[info[k].ndim];
        count=new size_t[info[k].ndim];
        switch(n_spacedims) {
          case 2:
            status= poc_getvar2d (p->c_str(), info[k].id, frame,(float *) buffer[0], &fmask ,info[k]);
/* *-------------------------------------------------------------------------------------
            interpolate on evenly spaced grid */
            status=map_export(sgrid,buffer[1],fmask,grid,buffer[0],fmask,0);
/* *-------------------------------------------------------------------------------------
            filter */
            buffer[2]=Loess2D_BF(cgrid, buffer[1], fmask, scale, weight, buffer[2]);
/* *-------------------------------------------------------------------------------------
            back to native grid */
            status=map_export(sgrid,buffer[2],fmask,grid,buffer[1],fmask,0);
            for(n=0;n<nvalues;n++) {
              if((buffer[0][n]!=fmask)&&(buffer[1][n]!=fmask)) {
                buffer[2][n]=buffer[0][n]-buffer[1][n];
                }
              else {
                buffer[2][n]=fmask;
                }
              }
//             start[0]=frame;
//             count[0]=1;
//             start[1]=0;
//             count[1]=decoded[0].ylen;
//             start[2]=0;
//             count[2]=decoded[0].xlen;
            start[0]=0;
            count[0]=decoded[0].ylen;
            start[1]=0;
            count[1]=decoded[0].xlen;
            status=poc_write(output.c_str(), raw[k], start, count, buffer[0]);
            status=poc_write(output.c_str(), filtered_LF[k], start, count, buffer[1]);
            status=poc_write(output.c_str(), filtered_HF[k], start, count, buffer[2]);
            break;
          case 3:
            status= poc_getvar3d (p->c_str(), info[k].id, frame,(float *) buffer[k], &fmask ,info[k]);
            break;
          default:
            status=-1;
            break;
          }
        if(status !=0) check_error(-1, "unable to read file", __LINE__, __FILE__, 1);;
        }
//      t=time[frame]+cnes_time(origine,'s')-cnes_time(start,'s');
      }
    if(time!=0) delete[] time;
    }

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int filter_01(const char **varnames, int nvars, vector<string> filelist, char *gridfile, char *controlfile, date_t start,date_t final)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  float *buffer[10],fmask,*hmean=0;
  double t,*time=0,duration;
  double step;
  int packing,count;
  int k,n;
  int status;
  int day,month,year,frame;
  char *s;
  char test[256];
  string output;
  grid_t grid;
  cdfgbl_t global[10];
  cdfvar_t info[10];
  cdfvar_t filtered_LF[nvars],filtered_HF[nvars],raw[nvars],detided_vtime,vtime;
  decoded_t decoded[10];
  date_t actual,origine;
  bool *mask;
  int nvalues,nrhs,nframes,n_spacedims;
  grid_t cgrid, sgrid;
  float scale, scale_x, scale_y, azimuth, *weight;

  status= cdf_globalinfo(filelist[0].c_str(),&global[0],0);
  for(k=0;k<nvars;k++) {
/* *----------------------------------------------------------------------------
    get variable dimension once and for all (assumes variables have same size !!!) */
    status= cdf_varinfo  (filelist[0].c_str(),varnames[k],&(info[k]),0);
    if(status==-1) {
      return(0);
      }

/* *----------------------------------------------------------------------------
    analyse netcdf informations */
    status= poc_decode_axis(info[k], global[0], &(decoded[k]));
    status= poc_decode_mask(info[k], &(decoded[k]));
    }
/* *----------------------------------------------------------------------------
  retain only space dimensions */
  nvalues=1;
  n_spacedims=0;
  switch (decoded[0].xlen) {
    case 0:
      break;
    default:
      nvalues*=decoded[0].xlen;
      n_spacedims++;
      break;
    }
  switch (decoded[0].ylen) {
    case 0:
      break;
    default:
      nvalues*=decoded[0].ylen;
      n_spacedims++;
      break;
    }
  switch (decoded[0].zlen) {
    case 0:
      break;
    default:
      nvalues*=decoded[0].zlen;
      n_spacedims++;
      break;
    }

  nrhs=nvars;
  nframes=decoded[0].tlen;
  if(nframes==0) nframes=1;
  
/*-----------------------------------------------------------------------------
  allocate memory for vectors */
  for(k=0;k<10;k++) {
    buffer[k]=new float[nvalues];
    if(buffer[k] == NULL) {
      printf("#memory allocation error for buffer N= %d \n",nvalues);
      check_error(-1, "memory allocation error", __LINE__, __FILE__, 1);;
      }
    }

  if(gridfile==0) gridfile=strdup(filelist[0].c_str());

  switch(n_spacedims) {
    case 2:
      status= poc_getgrid2d (gridfile, global[0], info[0], &grid);
      if(status!=0) {
        printf("cannot load 2D grid from gridfile %s\n",gridfile);
        }
      break;
    case 3:
/* *----------------------------------------------------------------------------
      get grid once and for all (assumes variables have same discretisation !!!) */
      int vdepth;
      status=poc_getgrid3d (gridfile,  global[0], info[0], &grid, &vdepth);
      if(status!=0) {
        printf("cannot load 3D grid from gridfile %s\n",gridfile);
        }
      break;
    default:
      status=-1;
      break;
    }
  status= get_cartesian(grid, &cgrid, &sgrid);

  scale=100000.;
  Loess2D_init(cgrid, scale, &weight);
  
/* *-------------------------------------------------------------------------------------
  construct harmonic matrix and RHS vectors */
  count=0;
  for(vector<string>::const_iterator p = filelist.begin(); p != filelist.end(); p++){
    count++;
/* *----------------------------------------------------------------------------
    get variable dimension once and for all (assumes variables have same size !!!) */
    status= cdf_varinfo  (p->c_str(),varnames[0],&(info[0]),0);
    status= poc_decode_axis(info[0], global[0], &(decoded[0]));
    nframes=decoded[0].tlen;
    if(nframes==0) nframes=1;
    status= poc_gettime(p->c_str(),info[0], global[0], &origine, &time, &nframes);
/* *-------------------------------------------------------------------------------------
    create/open the output files */
    output=string("filtered-");
    output+=*p;
    status= poc_createfile(output.c_str(),global[0]);
/* *-------------------------------------------------------------------------------------
    check filetered variables */
    for(k=0;k<nvars;k++) {
      sprintf(test,"%s",info[k].name);
      status= cdf_varinfo(output.c_str(),test,&(raw[k]),0);
      if(status==-1){
        raw[k]=info[k];
        printf("variable %s does not exist in %s, need to be created...\n",raw[k].name,output.c_str());
        status=create_ncvariable(output.c_str(), &(raw[k]));
        if(status!=NC_NOERR) return(-1);
        printf("done...\n");
        }
      sprintf(test,"%s-%s",info[k].name,"LF");
      status= cdf_varinfo(output.c_str(),test,&(filtered_LF[k]),0);
      if(status==-1){
        filtered_LF[k]=info[k];
        delete[] filtered_LF[k].name;
        filtered_LF[k].name=new char[strlen(info[k].name)+9];
        sprintf(filtered_LF[k].name,"%s-%s",info[k].name,"LF");
        printf("variable %s does not exist in %s, need to be created...\n",filtered_LF[k].name,output.c_str());
        status=create_ncvariable(output.c_str(), &(filtered_LF[k]));
        if(status!=NC_NOERR) return(-1);
        printf("done...\n");
        }
      sprintf(test,"%s-%s",info[k].name,"HF");
      status= cdf_varinfo(output.c_str(),test,&(filtered_HF[k]),0);
      if(status==-1){
        filtered_HF[k]=info[k];
        delete[] filtered_HF[k].name;
        filtered_HF[k].name=new char[strlen(info[k].name)+9];
        sprintf(filtered_HF[k].name,"%s-%s",info[k].name,"HF");
        printf("variable %s does not exist in %s, need to be created...\n",filtered_HF[k].name,output.c_str());
        status=create_ncvariable(output.c_str(), &(filtered_HF[k]));
        if(status!=NC_NOERR) return(-1);
        printf("done...\n");
        }
      }
/* *-------------------------------------------------------------------------------------
    read the output files */
    printf("read %s ...\n",p->c_str());
    for(frame=0;frame<nframes;frame++) {
      for(k=0;k<nvars;k+=2) {
        size_t *start, *count;
        start=new size_t[info[k].ndim];
        count=new size_t[info[k].ndim];
        switch(n_spacedims) {
          case 2:
            status= poc_getvar2d (p->c_str(), info[k].id,   frame,(float *) buffer[0], &fmask ,info[k]);
            status= poc_getvar2d (p->c_str(), info[k+1].id, frame,(float *) buffer[1], &fmask ,info[k+1]);
            status= filter2D_complex(grid, cgrid, sgrid,  buffer, fmask, weight, scale, nvalues);
            start[0]=0;
            count[0]=decoded[0].ylen;
            start[1]=0;
            count[1]=decoded[0].xlen;
            status=poc_write(output.c_str(), raw[k],   start, count, buffer[0]);
            status=poc_write(output.c_str(), raw[k+1], start, count, buffer[1]);
            status=poc_write(output.c_str(), filtered_LF[k],   start, count, buffer[2]);
            status=poc_write(output.c_str(), filtered_LF[k+1], start, count, buffer[3]);
            status=poc_write(output.c_str(), filtered_HF[k],   start, count, buffer[4]);
            status=poc_write(output.c_str(), filtered_HF[k+1], start, count, buffer[5]);
            break;
          case 3:
            status= poc_getvar3d (p->c_str(), info[k].id,   frame,(float *) buffer[k],   &fmask ,info[k]);
            status= poc_getvar3d (p->c_str(), info[k+1].id, frame,(float *) buffer[k+1], &fmask ,info[k+1]);
            status= filter3D_complex(grid, cgrid, sgrid,  buffer, fmask, weight, scale, nvalues);
            start[0]=0;
            count[0]=decoded[0].zlen;
            start[1]=0;
            count[1]=decoded[0].ylen;
            start[2]=0;
            count[2]=decoded[0].xlen;
            status=poc_write(output.c_str(), raw[k],   start, count, buffer[0]);
            status=poc_write(output.c_str(), raw[k+1], start, count, buffer[1]);
            status=poc_write(output.c_str(), filtered_LF[k],   start, count, buffer[2]);
            status=poc_write(output.c_str(), filtered_LF[k+1], start, count, buffer[3]);
            status=poc_write(output.c_str(), filtered_HF[k],   start, count, buffer[4]);
            status=poc_write(output.c_str(), filtered_HF[k+1], start, count, buffer[5]);
            break;
          default:
            status=-1;
            break;
          }
        if(status !=0) check_error(-1, "unable to read file", __LINE__, __FILE__, 1);;
        }
      }
    if(time!=0) delete[] time;
    }

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  
  double tau;
  double time,duration;
  int count;
  int fmt,i,j,k,l;
  int n,status;
  char filename[256],*gridfile=NULL,*controlfile=NULL,*keyword,*pathname=NULL,*s;
  char *varlist=NULL,*convention=NULL,*list=0;
  char rootname[256]="\0",tmp[256],output[256],vname[256];
  grid_t grid;
  int nonde;
  date_t start,final,reference,start_date;
  spectrum_t WaveList;
  date_t actual;
  int nvalues,nrhs,nodal_corrections=1,atlas=0;
  harmonic_t harmonic;
  hconstant_t **constants;
  int nvars;
  char **varnames;
  double t;
  vector<string> filelist;

  fct_echo( argc, argv);

  n=1;
  while (n < argc) {
    keyword=strdup(argv[n]);
    if(strcmp(keyword,"--nodal=no")==0) {
      nodal_corrections=0;
      n++;
      }
    else if(strcmp(keyword,"-control")==0) {
      controlfile= strdup(argv[n+1]);  /* directory */
      n++;
      n++;
      continue;
      }
    switch (keyword[0]) {
      case '-':
        switch (keyword[1]) {

        case 'c' :
          convention= strdup(argv[n+1]);  /* directory */
          n++;
          n++;
          break;

        case 'l' :
          list= strdup(argv[n+1]);  /* directory */
          n++;
          n++;
          break;

        case 'g' :
          gridfile= strdup(argv[n+1]);  /* directory */
          n++;
          n++;
          break;

        case 'p' :
          pathname= strdup(argv[n+1]);  /* directory */
          n++;
          n++;
          break;

        case 's' :
          s= strdup(argv[n+1]);
          n++;
          n++;
          sscanf(s,"%d/%d/%d",&start.day,&start.month,&start.year);
          break;

       case 'f' :
          s= strdup(argv[n+1]);
          n++;
          n++;
          sscanf(s,"%d/%d/%d",&final.day,&final.month,&final.year);
          break;

        case 'v' :
          varlist= strdup(argv[n+1]);
          n++;
          n++;
          break;

        default:
          __OUT_BASE_LINE__("unknown option %s\n",keyword);
          exit(-1);
        }
        break;

      default:
        filelist.push_back((string) argv[n]);
        n++;
        break;
      }
      free(keyword);
    }

  if(pathname==NULL) pathname=strdup(".");

  nvars=count_token(varlist);
  varnames=new char*[nvars];
  nvars=get_token(varlist, varnames, nvars);
  printf("#number of variables : %d \n",nvars);
  for(k=0;k<nvars;k++) {
    printf("#variables %d: %s \n",k,varnames[k]);
    }

  if(list!=0) {
    filelist=load_filelist(list);
    }
  status= parse_filelist(filelist, varnames, nvars, start, final);
  if(status==-1) goto error;

  status=filter_01((const char**) varnames,nvars,filelist,gridfile,controlfile, start,final);

  goto end;

error:
  __OUT_BASE_LINE__(" error, abort...\n");
  exit(-1);

end:
  __OUT_BASE_LINE__(" ^^^^^^^^^^^^^ end of computation ^^^^^^^^^^^^^\n");
  exit(0);
}
