


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

/* *----------------------------------------------------------------------------


  Extract time series from model's netcdf SG archives


-----------------------------------------------------------------------------*/

#define MAIN_SOURCE

#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "tools-structures.h"

#include "map.h"
#include "fe.h"
#include "archive.h"
#include "rutin.h" /*  rutin.h contains common utility routines  */
#include "poc-time.h"


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double beta[200][3],factor,x,y;
  float  z,u,v,ib,mask[3],dum;
  float  z_averaged,u_averaged,v_averaged,ib_averaged;
  double t0;
  float  *serie[200][4];
  float  **h,**p;
  float  count;
  double t,t1,t2,start,tstore;
  double time,sampling=0.0,average=0.0,ptime,htime;
  double *tmodel,tt;
  int nodes[200][3],search_status[200],zero=0,first=1;
  int i,j,k,l;
  int n,status,nst,m,mstore;
  int  year,month,day,hour;
  int nnde,elt;
  int start_year,last_year,start_month,last_month;
  date_t reference;

  size_t size;
  FILE *in,*out,*in1,*in2;
  char *h_root,*p_root,*pathname,*hfile,*pfile;
  char outname[256]="\0",*rootname=NULL,*default_dir=NULL,directory[256]="\0";
  pressure_station *sample;
  char *option,*keyword,*s,*input,*datafile;
  int filter=0,exact=0;
  meta_archive_t hinfo,pinfo;
  grid_t grid2d;

  int   nitems;
  long  pos;

  char *valid;

  fct_echo(argc,argv);
 
  n=1;
  while (n < argc) {
    if(strcmp(argv[n],"-start") ==0) {
      status=sscanf(argv[n+1],"%d/%d",&start_month,&start_year);
      n++;
      continue;
      }

    if(strcmp(argv[n],"-end") == 0) {
      status=sscanf(argv[n+1],"%d/%d",&last_month,&last_year);
      n++;
      continue;
      }

    keyword=strdup(argv[n]);
    switch (keyword[0]) {
      case '-':
        switch (keyword[1])
        {
/*---------------------------------------------------------------------
        root name for default ocean archive */
        case 'a' :
          sscanf(argv[n+1],"%lf",&average); /* en secondes */
          n++;
          n++;
          break;

/*---------------------------------------------------------------------
        root name for default ocean archive */
        case 'h' :
          h_root= strdup(argv[n+1]);
          n++;
          n++;
          break;

/*---------------------------------------------------------------------
        root name for default forcing archive */
        case 'p' :
          p_root= strdup(argv[n+1]);
          n++;
          n++;
          break;

/*---------------------------------------------------------------------
        filename for position input file */
        case 'i' :
          input= strdup(argv[n+1]);
          n++;
          n++;
          break;

/*---------------------------------------------------------------------
        pathname for default output directory */
        case 'd' :
          default_dir= strdup(argv[n+1]);
          n++;
          n++;
          break;

/*---------------------------------------------------------------------
        pathname for default archive directory */
        case 'r' :
          pathname= strdup(argv[n+1]);
          n++;
          n++;
          break;

/*---------------------------------------------------------------------
        time sampling of created time series */
        case 's' :
          sscanf(argv[n+1],"%lf",&sampling); /* en secondes */
          n++;
          n++;
          break;

/*---------------------------------------------------------------------
        filter option (not used) */
        case 'f' :
          filter=1;
          n++;
          break;

/*---------------------------------------------------------------------
        exact position for interpolation option */
        case 'e' :
          exact=1;
          n++;
          break;


        default:
          __OUT_BASE_LINE__("unknown option %s\n",keyword);
          exit(-1);

        }
        break;

      default:
        datafile= strdup(argv[n]);
        n++;
        break;
      }
      free(keyword);
    }

  if(default_dir==NULL) {
    default_dir=strdup("./");
    printf("use <./> as root name for default output directorys\n");
    }

  if(pathname==NULL) {
    pathname=strdup("./");
    printf("use <./> as root name for default input directorys\n");
    }

  in=fopen(input,"r");
  if(in==NULL) {
    __OUT_BASE_LINE__("cannot open file %s, abort... \n",input);
    exit(-1);
    }
  fscanf(in, "%d ", &nst);

  sample=(pressure_station *) malloc(nst*sizeof(pressure_station));
  valid=(char *) malloc(nst*sizeof(char));

  for(k=0; k<nst;k++) {
    fscanf(in," %lf %lf %s",&sample[k].t,&sample[k].p,sample[k].name);
    valid[k]=1;
    }

  fclose(in);

  if(h_root==NULL) {
    h_root= strdup("AUTO");
    printf("use <AUTO> as root name for regular archive files\n");
    }

  hfile=(char *) malloc(strlen(pathname)+256);
  sprintf(hfile,"%s/%s-%4.4d.%2.2d.huv.nc",pathname,h_root,start_year,start_month);

/*---------------------------------------------------------------------
  preliminarly fetch information on archive*/
  status=archive_info(hfile,2,&hinfo);
  if(status!=0) goto error;

/*---------------------------------------------------------------------
  allocate memory for time series */
  if(sampling == 0.) sampling=hinfo.sampling;

  t1=(julian_day(start_month,1,start_year)-CNES0jd)*24.*3600;
  if(last_month==12)
    t2=(julian_day(1,1,last_year+1)-CNES0jd)*24.*3600;
  else
    t2=(julian_day(last_month+1,1,last_year)-CNES0jd)*24.*3600;

  count=(float)(t2-t1)/hinfo.sampling;

  tmodel= (double *) malloc( ((int) count) *sizeof(double));
  for (k=0;k<nst;k++)
    for (j=0;j<4;j++) {
      serie[k][j]= (float *) malloc( ((int) count) *sizeof(float));
      if( serie[k][j]==NULL) {
        __OUT_BASE_LINE__("memory allocation failure, abort...\n");
        exit(-1);
        }
      }

  nnde=hinfo.grid.nx*hinfo.grid.ny;
  h=smatrix(0,2,0,nnde-1);
  p=smatrix(0,2,0,nnde-1);

/*---------------------------------------------------------------------
  initialize clock */
  t=hinfo.start+cnes_time(hinfo.reference,'s');

  archive_freeinfo(&hinfo);

/*---------------------------------------------------------------------
  truncate existing sample file if needed */
  for (k=0;k<nst;k++) {
    count=0;
    sprintf(outname,"%s/sample.%d",default_dir,k+1);
    if ((in = fopen(outname, "r")) == NULL) continue;
    pos=ftell(in);
    nitems=fscanf(in,"%lf %f %f %f %f\n", &time,&dum,&dum,&dum,&dum);
    while (!feof(in) && (nitems==5) && (time < t/3600.0/24.0)) {
      count++;
      pos=ftell(in);
      nitems=fscanf(in,"%lf %f %f %f %f\n", &time,&dum,&dum,&dum,&dum);
      if (time == t/3600.0/24.0) break;
      }
    fclose(in);
    truncate(outname,pos);
    }

/*---------------------------------------------------------------------
  start extraction */
  n=-1;
  m=0;
  first=1;
  for (year=start_year;year<=last_year;year++) {
    for (month=1;month<=12;month++) {
      if((year==start_year) && (month < start_month)) continue;
      if((year==last_year)  && (month > last_month))  break;

      sprintf(hfile,"%s/%s-%4.4d.%2.2d.huv.nc",pathname,h_root,year,month);

      status=archive_info(hfile,2,&hinfo);
      if(status!=0) goto error;

      grid2d=map_getgrid2d(hinfo.grid);

      for (i=0;i<(int) hinfo.nframe;i++) {
/*---------------------------------------------------------------------
        increment the frame number counted from beginning */
        n++;
        status=archive_read(hfile,TUGO_SG_NETCDF_FORMAT,hinfo,i,h,mask,&htime);
        htime=cnes_time(hinfo.reference,'s')+htime;
/*---------------------------------------------------------------------
        store time as CNES time */
        tmodel[n]=htime;
        if(status!=0) goto error;
        if(i%48 ==0)
          printf("t= %6.1f hours status= %d %s\n",htime/3600.0/24.0,status, sgetcnesdate(htime/3600.0));
/*---------------------------------------------------------------------
        convert to MKS if needed */
        factor=1.;
        for (k=0;k<nst;k++) {
          x=sample[k].t;
          y=sample[k].p;
          status=map_interpolation(grid2d,  h[0], mask[0], x, y, &(serie[k][0][n]));
          if((status !=0 ) && (valid[k] !=0 )) {
            printf("station %d  x=%lf, y=%lf out of domain \n",k,x,y);
            valid[k]=0;
            }
          status=map_interpolation(grid2d,  h[1], mask[1], x, y, &(serie[k][1][n]));
          status=map_interpolation(grid2d,  h[2], mask[2], x, y, &(serie[k][2][n]));
          status=map_interpolation(grid2d,  h[0], mask[0], x, y, &(serie[k][3][n]));
          serie[k][0][n] /=  factor;
          serie[k][1][n] /=  factor;
          serie[k][2][n] /=  factor;
          serie[k][3][n] /= -factor;
          }
        }
      tstore=t;
      mstore=m;
      printf(" sampling = %lf, from t=%lf to %lf \n",sampling,t/(24.0*3600.0),tmodel[n]/(24.0*3600.0));
      for (k=0;k<nst;k++) {
        if (valid[k]==0) continue;
        t=tstore;
        m=mstore;
        sprintf(outname,"%s/sample.%d",default_dir,k+1);
        out=fopen(outname,"a");
        while (t<=tmodel[n]) {
          if(average==0.0) {
            status=map_interpolate1D(serie[k][0],tmodel, mask[0],n+1, t, &z);
            status=map_interpolate1D(serie[k][1],tmodel, mask[1],n+1, t, &u);
            status=map_interpolate1D(serie[k][2],tmodel, mask[2],n+1, t, &v);
            status=map_interpolate1D(serie[k][3],tmodel, mask[0],n+1, t, &ib);
            }
          else {
            z_averaged=0.;
            u_averaged=0.;
            v_averaged=0.;
            ib_averaged=0.;
            count=0.;
            for(tt=t-average/2.;t<=t+average/2.;t+=1./24) {
              status=map_interpolate1D(serie[k][0],tmodel, mask[0],n+1, t, &z);
              status=map_interpolate1D(serie[k][1],tmodel, mask[1],n+1, t, &u);
              status=map_interpolate1D(serie[k][2],tmodel, mask[2],n+1, t, &v);
              status=map_interpolate1D(serie[k][3],tmodel, mask[0],n+1, t, &ib);
              if(z_averaged != mask[0]) {
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
            printf("%12.6f %6.3f %6.3f %6.3f\n",time,z, serie[k][0][m],z-serie[k][0][m]);
            }
          if(z != mask[0]) {
            time=t/(24.0*3600.0);
            fprintf(out,"%12.6f %6.3f %6.3f %6.3f %6.3f\n",time,z,u,v,ib);
            t+=sampling;
            m++;
            }
          else {
            time=t/(3600.0);
/*             printf("%12.6f %6.3f %6.3f %6.3f %6.3f status= %d %s\n",time,z,u,v,ib,status,sgetcnesdate(time)); */
            t+=sampling;
            }
          }
        fclose(out);
        }
      archive_freeinfo(&hinfo);
      archive_freeinfo(&pinfo);
/*---------------------------------------------------------------------
      end of month loop */
      }
/*---------------------------------------------------------------------
    end of year loop */
    }

  free(h);
  free(p);


  __ERR_BASE_LINE__("%s -computation complete ^^^^^^^^^\n",argv[0]);
  exit(0);

 error:
  __ERR_BASE_LINE__("%s -computation aborted ^^^^^^^^^\n",argv[0]);
  exit(-1);
}

