
/*******************************************************************************

  T-UGO tools, 2006-2016

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

<!-- USE FIND AND REPLACE TO EDIT THIS LINE, SO THAT print_help IS ALSO UPDATED -->
\brief Produces, from TUGOm outputs, sample files similar to those TUGOm would have done if told to do so.

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
#include <unistd.h>

#include "tools-define.h"
#include "tools-structures.h"

#include "fe.h"
#include "archive.h"
#include "poc-time.h"
#include "poc-netcdf-data.hpp"
#include "map.h"
#include "tides.h"
#include "functions.h"

#include "rutin.h"   /*rutin.h contains common utility routines  */

#include "mgr.h"


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int read_positions(char *input, pressure_station **sample_pos)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*
 \date 2013-01-17 Sara Fleury : modified
 \param Input ASCII file. Various formats authorized: [C] [NAME] LON LAT [NAME]. Lines starting with '#' are ignored. First available line must be the number of stations (integer).
 \param pointer on a  pressure_station table for the output
 \return nb of read stations
*/
{
#define DEBUG 0
  
  int k, nst;
  float lon,lat;
  char line[256], name[64];
  FILE *in;
  char c;
  pressure_station *local;
  enum  FORMAT {LON_LAT, LON_LAT_NAME, NAME_LON_LAT,
                C_LON_LAT, C_LON_LAT_NAME, C_NAME_LON_LAT};
  FORMAT  format;

  in=fopen(input,"r");
  if(in==0)
    TRAP_ERR_RETURN(errno,1,"fopen(\"%s\",\"r\") error (%d %s)\n",input,errno,strerror(errno));

  /* read number of stations */
  do {
    fgets(line,256,in);
#if DEBUG>0
    printf ("%s", line);
#endif
    } while (strlen(line)<1 || line[0] == '#');

  if (sscanf(line, "%d ", &nst) != 1) {
    fclose(in);
    STDOUT_BASE_LINE("error format for file %s (no nb of stations), abort... \n",input);
    }
  printf ("%d stations\n", nst);

  local=new pressure_station[nst];
    
  /* get format with first station */
  do {
    fgets(line,256,in);
    } while (line[0] == '#');

  name[0] = '#'; name[1] = '\0';
  if (sscanf(line, "%c %s %f %f", &c, name, &lon, &lat) == 4)
    format = C_NAME_LON_LAT;
  else if  (sscanf(line, "%c %f %f %s", &c, &lon, &lat, name) == 4) {
    if (name[0] == '#')
      format = C_LON_LAT;
    else
      format = C_LON_LAT_NAME;
    }
  else if  (sscanf(line, "%c %f %f", &c, &lon, &lat) == 3)
    format = C_LON_LAT;
  else if  (sscanf(line, "%s %f %f", name, &lon, &lat) == 3)
    format = NAME_LON_LAT;
  else if  (sscanf(line, "%f %f %s", &lon, &lat, name) == 3) {
    if (name[0] == '#')
      format = LON_LAT;
    else
      format = LON_LAT_NAME;
    }
  else if  (sscanf(line, "%f %f", &lon, &lat) == 2)
    format = LON_LAT;
  else {
    fclose(in);
    STDOUT_BASE_LINE("error station format for file %s, abort... \n",input);
    }

  /* record first station */
  k = 0;
  if (name[0] == '#') sprintf(local[k].name, "%s", k+1);
  else strcpy(local[k].name, name);
  local[k].t = lon;
  local[k].p = lat;
#if DEBUG>0
  printf ("%d %s %f %f\n", k+1, name, lon, lat);
#endif
  /* read each station */
  do {

    /* remove comment or blank lines */
    do {
      fgets(line,256,in);
      } while (strlen(line)<5 || line[0] == '#');

    name[0] = '#'; name[1] = '\0';
    switch(format) {
      case LON_LAT:
        sscanf(line, "%f %f", &lon, &lat);
        break;
      case LON_LAT_NAME:
        sscanf(line, "%f %f %s", &lon, &lat, name);
        break;
      case NAME_LON_LAT:
        sscanf(line, "%s %f %f", name, &lon, &lat);
        break;
      case C_LON_LAT:
        sscanf(line, "%c %f %f", &c, &lon, &lat);
        break;
      case C_LON_LAT_NAME:
        sscanf(line, "%c %f %f %s", &c, &lon, &lat, name);
        break;
      case C_NAME_LON_LAT:
        sscanf(line, "%c %s %f %f", &c, name, &lon, &lat);
        break;
      }

    /* record  station */
    k++;
    if (name[0] == '#') sprintf(local[k].name, "%d", k+1);
    else strcpy(local[k].name, name);
    local[k].t = lon;
    local[k].p = lat;
#if DEBUG>0
    printf ("%d %s %f %f\n", k+1, local[k].name, local[k].t, local[k].p);
#endif
    } while(k<nst-1);

  fclose(in);
  *sample_pos=local;
  
  return(nst);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int write_positions_01(char *output, pressure_station *local, int nst)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/// prepare
{
  int k;
  FILE *out;
//  pressure_station sample;

  out=fopen(output,"w");
  if(out==0)
    TRAP_ERR_RETURN(errno,1,"fopen(\"%s\",\"w\") error (%d %s)\n",output,errno,strerror(errno));

/* *----------------------------------------------------------------------------
  mimic usual input list for extraction softs */
  fprintf(out, "%d ", nst);

  for(k=0; k<nst;k++) {
    fprintf(out," %12lf %12lf %s  %d\n",local[k].t,local[k].p,local[k].name,k+1);
    }
  fclose(out);

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int write_positions_02(char *output, pressure_station *local, int nst)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/// prepare a file needed for site to sample translation
/*
Porto-Corsini      sample.39   12.285  44.492   0.0     0.0     corsini.res.65.mlf.gnu
*/
{
  int k;
  float georef=0.0,ignref=0.0;
  char sample[256];
  FILE *out;

  out=fopen(output,"w");
  if(out==0)
    TRAP_ERR_RETURN(errno,1,"fopen(\"%s\",\"w\") error (%d %s)\n",output,errno,strerror(errno));

/* *----------------------------------------------------------------------------
  mimic translation list for calibration softs */
  fprintf(out, "%d \n", nst);

  for(k=0; k<nst;k++) {
    pressure_station *localk=&local[k];
    fprintf(out,"%-19s ",localk->name);
    snprintf(sample,250uL,"sample.%-20s",localk->name);
    fprintf(out,"%-14s ",sample);
    if(strlen(localk->filename)==0){
      const size_t
        digitsAtStart=strspn(localk->name,"0123456789");
      if(digitsAtStart==3uL)
        sprintf(localk->filename,"h%.3s.gnu",localk->name);
      else
        snprintf(localk->filename,250uL,"%s.gnu",localk->name);
      }
    fprintf(out,"%12lf %12lf %12lf %12lf   %s\n",localk->t,localk->p,georef,ignref,localk->filename);
    }
  fclose(out);

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int truncate (const char *default_dir, pressure_station *sample, int nst,double t)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  float  dum;
  float  count;
  double time;
  int k;

  FILE *in;
  char outname[256]="\0";

  int   nitems;
  long  pos;


/* *---------------------------------------------------------------------
  truncate existing sample file if needed */
  for (k=0;k<nst;k++) {
    count=0;
    sprintf(outname,"%s/sample.%s",default_dir,sample[k].name);
    in = fopen(outname, "r");
    if(in==0) continue;
    pos=ftell(in);
    nitems=fscanf(in,"%lf %f %f %f %f\n", &time,&dum,&dum,&dum,&dum);
    while (!feof(in) && (nitems==5) && (time < t/d2s)) {
      count++;
      pos=ftell(in);
      nitems=fscanf(in,"%lf %f %f %f %f\n", &time,&dum,&dum,&dum,&dum);
      if (time == t/d2s) break;
      }
    fclose(in);
    truncate(outname,pos);
    }
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int resample (const char *default_dir, pressure_station *sample, int nst, int *search_status, vector<float>  *serie[4], vector<double>  tmodel, int nframes, double & tag, double sampling)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  float z,u,v,ib,mask=-9999.9;
  float z_averaged,u_averaged,v_averaged,ib_averaged;

  float  count;
  double t/*,tstore*/;
  double time,average=0.0;
  double tt;

//  int search_status[1000];
  int k;
  int n,status/*,m,mstore*/;

  FILE *out;
  char outname[256]="\0";

//  tstore=t;
//  mstore=m;
  
  n=nframes-1;
  
  printf(" sampling = %lf, from t=%lf to %lf (julian days from 1950-01-01)\n",sampling,tag/(24.0*3600.0),tmodel[n]/(24.0*3600.0));
  fflush(stdout);
  for (k=0;k<nst;k++) {
    t=tag;
//    m=mstore;
    if (search_status[k] !=0) {
      sprintf(outname,"%s/sample.%s.dubious",default_dir,sample[k].name);
      }
    else {
      sprintf(outname,"%s/sample.%s",default_dir,sample[k].name);
      }
    out=fopen(outname,"a");
    while (t<=tmodel[n]) {
      if(average==0.0) {
        status=map_interpolate1D(serie[k][0], tmodel, mask, n+1, t, &z);
        status=map_interpolate1D(serie[k][1], tmodel, mask, n+1, t, &u);
        status=map_interpolate1D(serie[k][2], tmodel, mask, n+1, t, &v);
        status=map_interpolate1D(serie[k][3], tmodel, mask, n+1, t, &ib);
        }
      else {
        z_averaged=0.;
        u_averaged=0.;
        v_averaged=0.;
        ib_averaged=0.;
        count=0.;
        for(tt=t-average/2.;t<=t+average/2.;t+=1./24) {
          status=map_interpolate1D(serie[k][0], tmodel, mask, n+1, t, &z);
          status=map_interpolate1D(serie[k][1], tmodel, mask, n+1, t, &u);
          status=map_interpolate1D(serie[k][2], tmodel, mask, n+1, t, &v);
          status=map_interpolate1D(serie[k][3], tmodel, mask, n+1, t, &ib);
          if(z_averaged != mask) {
            z_averaged+=z;
            u_averaged+=u;
            v_averaged+=v;
            ib_averaged+=ib;
            count++;
            }
          }
        z=z_averaged/count;
        u=u_averaged/count;
        v=v_averaged/count;
        ib=ib_averaged/count;
        }
//       if(z != serie[k][0][m]) {
// //        printf("%12.6f %7.4f %7.4f %7.4f\n",time,z, serie[k][0][m],z-serie[k][0][m]);
//         }
      if(z != mask) {
        time=t/(24.0*3600.0);
        fprintf(out,"%12.6f %7.4f %7.4f %7.4f %7.4f\n",time,z,u,v,ib);
        t+=sampling;
//        m++;
        }
      else {
        time=t/(24.0*3600.0);
//        printf("%12.6f %7.4f %7.4f %7.4f %7.4f status= %d %s\n",time,z,u,v,ib,status,sgetcnesdate(time));
        }
      }
    fclose(out);
    }

  tag=t;
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void fe_nearestvertex_and_beta(mesh_t & mesh,const discretisation_t &descriptor,int nst,pressure_station *sample,int *search_status,int **nodes,double **beta)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,elt;
  
  for (k=0;k<nst;k++) {
    sample[k].code=fe_nearest_vertex(mesh,sample[k].t, sample[k].p, (int *) NULL, 0);
    sample[k].h=mesh.vertices[sample[k].code].h;
    printf("station %25s: x=%7.2f, y=%7.2f, h=%7.2f, node=%6d ", sample[k].name,sample[k].t, sample[k].p, sample[k].h, sample[k].code);
    search_status[k]=fe_beta(mesh, descriptor, sample[k].t, sample[k].p,&elt,nodes[k],beta[k]);
    
    if(search_status[k]!=0) {
      printf("out of domain\n");
      }
    else {
      printf("ok\n");
      }
    
    fflush(stdout);
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int process_binary (const char *pathname, const char *h_root, const char *p_root,  const char *default_dir, pressure_station *sample, date_t start, date_t last, int nst,
                      float scale, int exact, double sampling)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double **beta;
  float  z,u,v,ib,mask=-9999.9,dum;
  float z_averaged,u_averaged,v_averaged,ib_averaged;
  float  *serie[1000][4];

  float  **h,**p;
  float  count;
  double t,t1,t2,tstore;
  double time,average=0.0,ptime,htime;
  double *tmodel,tt;

  int origine, **nodes,search_status[1000];
  int i,j,k;
  int n,status,m,mstore;
  int  year,month;
  int nnde;
  date_t reference;

  FILE *in,*out,*in1,*in2;
  char *hfile,*pfile;
  char outname[256]="\0";

  meta_archive_t hinfo,pinfo;
  discretisation_t descriptor;

  int   nitems;
  long  pos;
  
  if(h_root==NULL) {
    h_root= strdup("analysis");
    printf("use <analysis> as root name for sea state archive files\n");
    }
  if(p_root==NULL) {
    p_root= strdup("forcing");
    printf("use <forcing> as root name for forcing archive files\n");
    }

  hfile=(char *) malloc(strlen(pathname)+256);
  pfile=(char *) malloc(strlen(pathname)+256);
  sprintf(hfile,"%s/%s-%4.4d.%2.2d",pathname,h_root,start.year,start.month);

/*------------------------------------------------------------------------------
  preliminarly fetch information on model mesh */
  status=clx_archiveinfo(hfile,&hinfo);
  if(status!=0) TRAP_ERR_RETURN(status,1,"clx_archiveinfo(\"%s\",) error %d\n",hfile,status);
  printf("compiling %d vertices mesh... ",hinfo.mesh.nvtxs);fflush(stdout);
  status=fe_list(&hinfo.mesh);
  if(status!=0) TRAP_ERR_EXIT(status,"fe_list() error %d\n",status);
  printf("done (%d)\n",status);fflush(stdout);
  
  const int discretisation=LGP1;
  discretisation_init(&hinfo.mesh,discretisation,0);
  descriptor=get_descriptor(hinfo.mesh,discretisation);
  
  nodes=new int*[nst];
  for (k=0;k<nst;k++) {
    nodes[k]=new int[descriptor.nnpe];
    }
  
  beta=new double*[nst];
  for (k=0;k<nst;k++) {
    beta[k]=new double[descriptor.nnpe];
    }
  
/*------------------------------------------------------------------------------
  find FE elements and compute interpolation coefficients*/
  fe_nearestvertex_and_beta(hinfo.mesh,descriptor,nst,sample,search_status,nodes,beta);
  
/*------------------------------------------------------------------------------
  allocate memory for time series */
  if(sampling == 0.) sampling=hinfo.sampling;

  origine=julian_day(1,1,1950);
  t1=(julian_day(start.month,1,start.year)-origine)*24.*3600;
  if(last.month==12)
    t2=(julian_day(1,1,last.year+1)-origine)*24.*3600;
  else
    t2=(julian_day(last.month+1,1,last.year)-origine)*24.*3600;

  count=(float)(t2-t1)/hinfo.sampling;

  tmodel= (double *) malloc( ((int) count) *sizeof(double));
  for (k=0;k<nst;k++)
    for (j=0;j<4;j++) {
      serie[k][j]= (float *) malloc( ((int) count) *sizeof(float));
      if( serie[k][j]==NULL) {
        STDOUT_BASE_LINE("memory allocation failure, abort...\n");
        exit(-1);
        }
      }

  nnde=hinfo.mesh.nvtxs;
  h=smatrix(0,2,0,nnde-1);
  p=smatrix(0,2,0,nnde-1);

/*------------------------------------------------------------------------------
  initialize clock */
  t=hinfo.start+cnes_time(hinfo.reference,'s');

  archive_freeinfo(&hinfo);

/* *---------------------------------------------------------------------
  truncate existing sample file if needed */
  for (k=0;k<nst;k++) {
    count=0;
    sprintf(outname,"%s/sample.%s",default_dir,sample[k].name);
    in = fopen(outname, "r");
    if(in==0) continue;
    pos=ftell(in);
    nitems=fscanf(in,"%lf %f %f %f %f\n", &time,&dum,&dum,&dum,&dum);
    while (!feof(in) && (nitems==5) && (time < t/d2s)) {
      count++;
      pos=ftell(in);
      nitems=fscanf(in,"%lf %f %f %f %f\n", &time,&dum,&dum,&dum,&dum);
      if (time == t/d2s) break;
      }
    fclose(in);
    truncate(outname,pos);
    }

/*------------------------------------------------------------------------------
  start extraction */
  n=-1;
  m=0;
  for (year=start.year;year<=last.year;year++) {
    for (month=1;month<=12;month++) {
      if((year==start.year) && (month < start.month)) continue;
      if((year==last.year)  && (month > last.month))  break;

      sprintf(hfile,"%s/%s-%4.4d.%2.2d",pathname,h_root,year,month);
      sprintf(pfile,"%s/%s-%4.4d.%2.2d",pathname,p_root,year,month);

      status=clx_archiveinfo(hfile,&hinfo);
      if(status!=0) {
        TRAP_ERR_RETURN(status,1,"loading %s failed\n",hfile);
        }
      status=clx_archiveinfo(pfile,&pinfo);
      if(status!=0) {
        TRAP_ERR_RETURN(status,1,"loading %s failed\n",pfile);
        }

/*------------------------------------------------------------------------------
      h and p file consistency check */
      if( hinfo.nframe!=pinfo.nframe ) {
        TRAP_ERR_RETURN(-1,1,"loading frame number in %s (%d) and %s (%d)",hfile,hinfo.nframe,pfile,pinfo.nframe);
        }

      status=fe_list(&hinfo.mesh);
      if(status!=0) {
        TRAP_ERR_RETURN(status,1,"archive mesh error\n");
        }

/*------------------------------------------------------------------------------
      open archive files */
      in1=fopen(hfile,"r");
      in2=fopen(pfile,"r");

      for (i=0;i<(int) hinfo.nframe;i++) {
/*------------------------------------------------------------------------------
        increment the frame number counted from beginning */
        n++;
        status=clx_archiveread(in1,hinfo,i+1,h,&htime);
        status=clx_archiveread(in2,pinfo,i+1,p,&ptime);
        htime=cnes_time(hinfo.reference,'s')+htime;
        ptime=cnes_time(pinfo.reference,'s')+ptime;
/*------------------------------------------------------------------------------
        h and p file consistency check */
        if(ptime != htime) return -1;
/*------------------------------------------------------------------------------
        store time as CNES time */
        tmodel[n]=htime;
        if(status!=0) return -1;
        if(i%48 ==0){
          printf("t= %6.1f hours status= %d %s\n",htime/d2s,status, sgetcnesdate(htime/3600.0));
          fflush(stdout);
          }
        for (k=0;k<nst;k++)
          if((exact==0) || (search_status[k] !=0)) {
            serie[k][0][n] =  h[0][sample[k].code]*scale;
            serie[k][1][n] =  h[1][sample[k].code]*scale;
            serie[k][2][n] =  h[2][sample[k].code]*scale;
            serie[k][3][n] = -p[0][sample[k].code]*scale;
            }
          else {
            serie[k][0][n] = (beta[k][0]*h[0][nodes[k][0]]
                           +  beta[k][1]*h[0][nodes[k][1]]
                           +  beta[k][2]*h[0][nodes[k][2]])*scale;
            serie[k][1][n] = (beta[k][0]*h[1][nodes[k][0]]
                           +  beta[k][1]*h[1][nodes[k][1]]
                           +  beta[k][2]*h[1][nodes[k][2]])*scale;
            serie[k][2][n] = (beta[k][0]*h[2][nodes[k][0]]
                           +  beta[k][1]*h[2][nodes[k][1]]
                           +  beta[k][2]*h[2][nodes[k][2]])*scale;
            serie[k][3][n] = -(beta[k][0]*p[0][nodes[k][0]]
                           +   beta[k][1]*p[0][nodes[k][1]]
                           +   beta[k][2]*p[0][nodes[k][2]])*scale;
            }
        }
      fclose(in1);
      fclose(in2);
      tstore=t;
      mstore=m;
      printf(" sampling = %lf, from t=%lf to %lf \n",sampling,t/(24.0*3600.0),tmodel[n]/(24.0*3600.0));
      fflush(stdout);
      for (k=0;k<nst;k++) {
        t=tstore;
        m=mstore;
        if (search_status[k] !=0) {
          sprintf(outname,"%s/sample.%s.dubious",default_dir,sample[k].name);
          }
        else {
          sprintf(outname,"%s/sample.%s",default_dir,sample[k].name);
          }
        out=fopen(outname,"a");
        while (t<=tmodel[n]) {
          if(average==0.0) {
            status=map_interpolate1D(serie[k][0],tmodel, mask,n+1, t, &z);
            status=map_interpolate1D(serie[k][1],tmodel, mask,n+1, t, &u);
            status=map_interpolate1D(serie[k][2],tmodel, mask,n+1, t, &v);
            status=map_interpolate1D(serie[k][3],tmodel, mask,n+1, t, &ib);
            }
          else {
            z_averaged=0.;
            u_averaged=0.;
            v_averaged=0.;
            ib_averaged=0.;
            count=0.;
            for(tt=t-average/2.;t<=t+average/2.;t+=1./24) {
              status=map_interpolate1D(serie[k][0],tmodel, mask,n+1, t, &z);
              status=map_interpolate1D(serie[k][1],tmodel, mask,n+1, t, &u);
              status=map_interpolate1D(serie[k][2],tmodel, mask,n+1, t, &v);
              status=map_interpolate1D(serie[k][3],tmodel, mask,n+1, t, &ib);
              if(z_averaged != mask) {
                z_averaged+=z;
                u_averaged+=u;
                v_averaged+=v;
                ib_averaged+=ib;
                count++;
                }
              }
            z=z_averaged/count;
            u=u_averaged/count;
            v=v_averaged/count;
            ib=ib_averaged/count;
            }
          if(z != serie[k][0][m]) {
//            printf("%12.6f %7.4f %7.4f %7.4f\n",time,z, serie[k][0][m],z-serie[k][0][m]);
            }
          if(z != mask) {
            time=t/(24.0*3600.0);
            fprintf(out,"%12.6f %7.4f %7.4f %7.4f %7.4f\n",time,z,u,v,ib);
            t+=sampling;
            m++;
            }
          else {
            time=t/(24.0*3600.0);
//            printf("%12.6f %7.4f %7.4f %7.4f %7.4f status= %d %s\n",time,z,u,v,ib,status,sgetcnesdate(time));
            }
          }
        fclose(out);
        }
      archive_freeinfo(&hinfo);
      archive_freeinfo(&pinfo);
/*------------------------------------------------------------------------------
      end of month loop */
      }
/*------------------------------------------------------------------------------
    end of year loop */
    }

  status=clx_archiveinfo(hfile,&hinfo);
  if(status!=0) TRAP_ERR_RETURN(status,1,"clx_archiveinfo(\"%s\",) error %d\n",hfile,status);

  archive_freeinfo(&hinfo);

  free(h);
  free(p);

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int process_netcdf (const char *name_template, const char *default_dir, pressure_station *sample, date_t start, date_t last, int nst,
                      int discretisation, float scale, int exact, double sampling)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int k/*,n*/;
  int month, year, *search_status;
  char *input;
  mesh_t mesh;
  discretisation_t descriptor;
  int **nodes;
  double **beta;
  date_t date;

  vector<float>  **serie;
  vector<double> tmodel;
  poc_global_t global;
  poc_data_t<double> h,u,v,p,time;
  date_t current;
  double t, tag;
  const char *hname="elevation",*uname="ubar",*vname="vbar",*pname="Pa",*tname="time";

  status=decode_name(start, 0, name_template, &input);
  status=fe_readgeometry(input, &mesh);
 
  discretisation_init(&mesh,discretisation,0);
  descriptor=get_descriptor(mesh,discretisation);
  
  nodes=new int*[nst];
  for (k=0;k<nst;k++) {
    nodes[k]=new int[descriptor.nnpe];
    }
  
  beta=new double*[nst];
  for (k=0;k<nst;k++) {
    beta[k]=new double[descriptor.nnpe];
    }
  
  serie=new vector<float> *[nst];
  for (k=0;k<nst;k++) {
    serie[k]=new vector<float> [4];
    }
  
  search_status=new int[nst];
  
  status=poc_inq(input,&global);
      
  if(status!=0) return(-1);
  
  status=time.init(global,tname);
  status=time.read_data(input,0,0,1);
  
/*------------------------------------------------------------------------------
  we need cnes time in seconds */
  t=time.data[0];
  
  tag=t;
  
//  t=cnes_time(hinfo.reference,'s');
  status=truncate (default_dir, sample, nst, t);
  
/*------------------------------------------------------------------------------
  find FE elements and compute interpolation coefficients*/
  fe_nearestvertex_and_beta(mesh,descriptor,nst,sample,search_status,nodes,beta);
  
  for (year=start.year;year<=last.year;year++) {
    for (month=1;month<=12;month++) {
      if((year==start.year) && (month < start.month)) continue;
      if((year==last.year)  && (month > last.month))  break;
      date_t current(year,month,1,0.);
      status=decode_name(current, 0, name_template, &input);
      status=poc_inq(input,&global);
      if(status!=0) break;

      status=h.init(global,hname);
      if(status!=0) NC_TRAP_ERROR(wexit,status,1,"poc_data_t::init(\"%s\",\"%s\") error",input,hname);
      status=u.init(global,uname);
      status=v.init(global,vname);
      
      status=p.init(global,pname);
      
      status=time.init(global,tname);
      
      for (size_t i=0;i<h.nframes;i++) {
        status=h.read_data(input,i,0,1);
        status=u.read_data(input,i,0,1);
        status=v.read_data(input,i,0,1);
        status=p.read_data(input,i,0,1);
        status=time.read_data(input,i,0,1);
        tmodel.push_back(time.data[0]);
        date=poctime_getdatecnes(time.data[0],'s');
        if(date.second==0){
          printf("%4d/%2d/%2d %d\n",date.year,date.month,date.day,i);
          fflush(stdout);
          }
        for (k=0;k<nst;k++) {
          if((exact==0) || (search_status[k] !=0)) {
            size_t node=sample[k].code;
            vector<float> *seriek=serie[k];
            seriek[0].push_back( h.data[node]*scale);
            if(u.data!=0 and v.data!=0){
              seriek[1].push_back( u.data[node]*scale);
              seriek[2].push_back( v.data[node]*scale);
              }
            else{
              seriek[1].push_back(0.);
              seriek[2].push_back(0.);
              }
            if(p.data!=0) seriek[3].push_back(-p.data[node]*scale);
            else seriek[3].push_back(0.);
            }
//           else {
//             vector<float> *seriek=serie[k];
//             seriek[0][n] = (beta[k][0]*h[0][nodes[k][0]]
//                            +  beta[k][1]*h[0][nodes[k][1]]
//                            +  beta[k][2]*h[0][nodes[k][2]])*scale;
//             seriek[1][n] = (beta[k][0]*h[1][nodes[k][0]]
//                            +  beta[k][1]*h[1][nodes[k][1]]
//                            +  beta[k][2]*h[1][nodes[k][2]])*scale;
//             seriek[2][n] = (beta[k][0]*h[2][nodes[k][0]]
//                            +  beta[k][1]*h[2][nodes[k][1]]
//                            +  beta[k][2]*h[2][nodes[k][2]])*scale;
//             seriek[3][n] = -(beta[k][0]*p[0][nodes[k][0]]
//                            +   beta[k][1]*p[0][nodes[k][1]]
//                            +   beta[k][2]*p[0][nodes[k][2]])*scale;
//             }
          }
        }
      status=resample (default_dir, sample, nst, search_status, serie, tmodel, h.nframes, tag, sampling);
      for (k=0;k<nst;k++) for(int l=0;l<4;l++) serie[k][l].clear();
      h.destroy_data();
      u.destroy_data();
      v.destroy_data();
      p.destroy_data();
      tmodel.clear();
      delete[] input;
//       ~h();
//       ~u();
      }
    }
  
  mesh.destroy();
  
  for (k=0;k<nst;k++) {
    delete[] nodes[k];
    }
  delete[] nodes;
  
  for (k=0;k<nst;k++) {
    delete[] beta[k];
    }
  delete[] beta;
  
  for (k=0;k<nst;k++) {
    delete[] serie[k];
    }
  delete[] serie;

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void print_help(const char *prog_name)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/** \brief prints help of this programme.

\param prog_name name of the program to be printed by this help function : as[0]
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
    "  %s [OPTIONS]\n",prog_name);
  printf("\n"
    "DESCRIPTION\n"
/* USE FIND AND REPLACE TO EDIT THE LINE BELOW, SO THAT THE FILE HEAD IS ALSO UPDATED */
    "  Produces, from TUGOm outputs, sample files similar to those TUGOm would have done if told to do so.\n"
    "The files have 5 columns:\n"
    "  - time\n"
    "  - elevation\n"
    "  - zonal velocity\n"
    "  - meridional velocity\n"
    "  - inverse barometer\n"
    "\n"
    "OPTIONS\n"
    "  --help  Show this help and exit.\n"
    "  -i  followed by sation list file, as produced by detidor\n"
    "  -start  followed by start date in mm/yyyy format\n"
    "  -end  followed by end date in mm/yyyy format\n"
    "  -cm  scale to centimeters\n"
    "  -format  followed by anything. If specified, format will be NetCDF; otherwise, binary\n"
    "  -r  followed by directory of archive files (binary format). Default: .\n"
    "  -d  followed by output directory. Default: .\n"
    "  -s  followed by sampling period in seconds. Default: 3600\n"
    "  -h  followed by root name for analysis files (binary format). Default: analysis\n"
    "      When used as last argument, show this help and exit.\n"
    "  -p  followed by root name for forcing files (binary format). Default: forcing\n"
    "  -t  followed by path to translation file. Default: translation\n"
    "  -c  followed by input file name convention (NetCDF format)\n"
    "  -e  With binary format: activate LGP1 interpolation to exact position.\n"
    "      With NetCDF format: analyse random data contained in allocated buffers, i.e. produce rubish.\n"
    "\n"
    "BUGS\n"
    "  See option -e above.\n"
    "  Even when -format is followed by `binary', the format will be NetCDF.\n"
    "\n"
    "EXAMPLES\n"
    "  time ./objects/extract-nf-v1 -i ../h4_list_stat -s 3600 -cm -r ~/tugo/data/archives/global-MR/ -start 01/2001 -end 01/2001\n"
    "  extract-nf-v1 -i /gsa/fles/gloss/hourly-20120625/analysis.ondes65.FES2004.01y/stations_list -s 3600 -cm -r /home/fles/tugo/data/archives/global-MR/ -start 01/2001 -end 02/2001\n"
    );
  /** \endcode */
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double sampling=0.0;
  double scale=1.0;
  int k;
  int n,status,nst;
  int start_year,last_year,start_month,last_month;
  date_t reference, start, last;

  FILE *out;
  char *h_root=NULL,*p_root=NULL,*translation=0,*convention=NULL,*input=NULL,*format=NULL,*pathname=NULL;
  char *default_dir=NULL;
  pressure_station *sample=0;
  const char *keyword,*next_arg;
  int exact=0;
  meta_archive_t hinfo,pinfo;

  if (argc<2) {
    print_help(argv[0]);
    exit(-1);
    }

  fct_echo( argc, argv);

  n=1;
  while (n < argc) {
    
    if(strcmp(argv[n],"--help") ==0) {
      print_help(argv[0]);
      exit(0);
      }

    if(strcmp(argv[n],"-start") ==0) {
      status=sscanf(argv[n+1],"%d/%d",&start_month,&start_year);
      n++;
      n++;
      continue;
      }

    if(strcmp(argv[n],"-end") == 0) {
      status=sscanf(argv[n+1],"%d/%d",&last_month,&last_year);
      n++;
      n++;
      continue;
      }

    if(strcmp(argv[n],"-format") == 0) {
      format=strdup(argv[n+1]);
      n++;
      n++;
      continue;
      }
    
    if(strcmp(argv[n],"-cm")==0) {
      scale=0.01;
      n++;
      continue;
      }

    keyword=argv[n];
    
    switch (keyword[0]) {
      case '-':
        switch (keyword[1]) {
/*------------------------------------------------------------------------------
        root name for default ocean archive */
        case 'c' :
          convention= strdup(argv[n+1]);
          n++;
          n++;
          break;

/*------------------------------------------------------------------------------
        root name for default ocean archive */
        case 'h' :
          next_arg=argv[n+1];
          if(next_arg==0){
            print_help(argv[0]);
            exit(-1);
            }
          h_root= strdup(next_arg);
          n++;
          n++;
          break;

/*------------------------------------------------------------------------------
        root name for default forcing archive */
        case 'p' :
          p_root= strdup(argv[n+1]);
          n++;
          n++;
          break;

/*------------------------------------------------------------------------------
        path to translation file */
        case 't' :
          translation= strdup(argv[n+1]);
          n++;
          n++;
          break;

/*------------------------------------------------------------------------------
        filename for position input file */
        case 'i' :
          input= strdup(argv[n+1]);
          n++;
          n++;
          break;

/*------------------------------------------------------------------------------
        pathname for default output directory */
        case 'd' :
          default_dir= strdup(argv[n+1]);
          n++;
          n++;
          break;

/*------------------------------------------------------------------------------
        pathname for default archive directory */
        case 'r' :
          pathname= strdup(argv[n+1]);
          n++;
          n++;
          break;

/*------------------------------------------------------------------------------
        time sampling of created time series */
        case 's' :
          sscanf(argv[n+1],"%lf",&sampling); /* en secondes */
          n++;
          n++;
          break;

/*------------------------------------------------------------------------------
        exact position for interpolation option */
        case 'e' :
          exact=1;
          n++;
          break;

        default:
          STDOUT_BASE_LINE("unknown option %s\n",keyword);
          exit(-1);
        }
        break;

      default:
        STDOUT_BASE_LINE("unknown option %s\n",keyword);
        exit(-1);
      }
    
    }

  if(translation==0) {
    translation=strdup("translation");
    printf("use <%s> as default translation file path\n",translation);
    }
  
  bool missingOption=false;
  
  if(input==NULL) {
    printf("*** Please specify station list with option -i ***\n");
    missingOption=true;
    }
  
  if(format!=NULL and convention==NULL) {
    printf("*** Please specify convention with option -c ***\n");
    missingOption=true;
    }

  if(missingOption == true) {
    print_help(argv[0]);
    wexit(-1);
    }

  if(default_dir==NULL) {
    default_dir=strdup("./");
    printf("use <./> as default output directory\n");
    }

  if(pathname==NULL and format==NULL) {
    pathname=strdup("./");
    printf("use <./> as default input directory\n");
    }

  if(sampling==0.0) {
    sampling=3600.;
    printf("sampling interval not specified, use 3600s as a default\n");
    }
  
  fflush(stdout);
  
/*------------------------------------------------------------------------------
  read input sampling positions */
  nst =  mgr_read_stations(input,&sample);
//   nst= read_positions(input,&sample);
//   nst= read_hawaii_positions(input,&sample);
//   status= write_positions_01("reformatted.list",sample,nst);
  status= write_positions_02(translation,sample,nst);
  if(status!=0) TRAP_ERR_EXIT(status,"write_positions_02(\"%s\",,) error (%d %s)\n",translation,status,strerror(status));

  start.year =start_year;
  start.month=start_month;
  
  last.year =last_year;
  last.month=last_month;
  
  if(format==NULL) {
    status=process_binary (pathname, h_root, p_root, default_dir, sample, start, last, nst, scale, exact, sampling);
    if(status!=0) TRAP_ERR_EXIT(status,"process_binary() error %d\n");
    }
  else {
    status=process_netcdf (convention, default_dir, sample, start, last, nst, LGP1, scale, exact, sampling);
    if(status!=0) TRAP_ERR_EXIT(status,"process_netcdf() error %d\n");
    }

  out=fopen("sample.info","w");
  for (k=0;k<nst;k++) {
    fprintf(out," %9.3f %9.3f %20s %6d %7.4f   \t%d\n", sample[k].t, sample[k].p, sample[k].name, sample[k].code, sample[k].h, k+1);
    }
  fclose(out);

  TRAP_ERR_EXIT(0,"%s -computation complete ^^^^^^^^^\n",argv[0]);
}

