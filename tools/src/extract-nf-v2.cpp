
/*******************************************************************************

  T-UGO tools, 2006-2019

  Unstructured Ocean Grid initiative

*******************************************************************************/

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


  Extract time series from model's netcdf UG archives


xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

#define MAIN_SOURCE

#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "tools-structures.h"

#include "fe.h"
#include "tides.h"
#include "mgr-converter.h"
#include "archive.h"
#include "rutin.h"   /*rutin.h contains common utility routines  */
#include "map.h"
#include "poc-time.h"
#include "functions.h"
#include "filter.h"

#define GLOBAL      0
#define PER_YEAR    1
#define PER_MONTH   2
#define PER_WEEK    3
#define PER_DAY     4
#define PER_HOUR    5
#define PER_MINUTE  6
#define PER_SECOND  7
#define CUSTOM      8


// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   int decode_name(date_t actual, const char *name_template, const char *varname, char **filename, int *packing)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   int  status;
//   char dummy[256], *pointer;
//   FILE *out;
// 
// /*------------------------------------------------------------------------------
//   Reconstruct the input file name*/
//   (*filename) = new char[strlen(name_template) + 256];
//   sprintf((*filename), "%s", name_template);
// 
//   out = fopen(*filename, "r");
//   status = (out == NULL);
// 
//   switch (status) {
//     case 0:
// /*------------------------------------------------------------------------------
//       file exists, do nothing more*/
//       fclose(out);
//       *packing=GLOBAL;
//       break;
// 
//     default:
// /*------------------------------------------------------------------------------
//       use format convention information*/
// /*------------------------------------------------------------------------------
//       substitute YYYY with current year*/
//       pointer = strstr((*filename), "YYYY");
//       if(pointer != NULL) {
//         sprintf(dummy, "%4d", actual.year);
//         strncpy(pointer, dummy, 4);
//         *packing=PER_YEAR;
//         }
// /*------------------------------------------------------------------------------
//       substitute MM with current month*/
//       pointer = strstr((*filename), "MM");
//       if(pointer != NULL) {
//         sprintf(dummy, "%2.2d", actual.month);
//         strncpy(pointer, dummy, 2);
//         *packing=PER_MONTH;
//         }
// /*------------------------------------------------------------------------------
//       substitute DD with current day*/
//       pointer = strstr((*filename), "DD");
//       if(pointer != NULL) {
//         sprintf(dummy, "%2.2d", actual.day);
//         strncpy(pointer, dummy, 2);
//         *packing=PER_DAY;
//         }
// /*------------------------------------------------------------------------------
//       substitute DD with current day*/
//       pointer = strstr((*filename), "NNN");
//       if(pointer != NULL) {
//         int day=day_in_year(actual);
//         sprintf(dummy, "%3.3d", day);
//         strncpy(pointer, dummy, 3);
//         *packing=PER_DAY;
//         }
// /*------------------------------------------------------------------------------
//       substitute HH with current hour*/
//       pointer = strstr((*filename), "HH");
//       if(pointer != NULL) {
//         sprintf(dummy, "%2.2d", (int) floor(actual.second/3600.));
//         strncpy(pointer, dummy, 2);
//         *packing=PER_HOUR;
//         }
// /*------------------------------------------------------------------------------
//       substitute MN with current hour*/
//       pointer = strstr((*filename), "MN");
//       if(pointer != NULL) {
//         double hour=floor(actual.second/3600.);
//         sprintf(dummy, "%2.2d", (int) floor((actual.second-3600.*hour)/60.));
//         strncpy(pointer, dummy, 2);
//         *packing=PER_HOUR;
//         }
// // /*------------------------------------------------------------------------------
// //       lookup file fitting ???? with minute and seconds index*/
// //       pointer = strstr((*filename), "????");
// //       if(pointer != NULL) {
// //         struct stat buf;
// //         int i,j;
// //         for(j=0;j<60;j++) {
// //           for(i=0;i<60;i++) {
// //             sprintf(dummy, "%2.2d%2.2d", j,i);
// //             strncpy(pointer, dummy, 4);
// //             status= get_file_size(*filename, 0);
// // /*------------------------------------------------------------------------------
// //             status=0 if file found*/
// //             if(status==0) break;
// //             }
// //           if(status==0) break;
// //           }
// //         }
// /*------------------------------------------------------------------------------
//       substitute VARNAME*/
//       pointer = strstr((*filename), "VARNAME");
//       if(pointer != NULL) {
//         sprintf(dummy, "%s", varname);
//         strncpy(pointer, dummy, strlen(varname));
//         }
//       break;
// 
//       break;
//     }
// 
//   return (status);
// }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double beta[200][3],factor;
  float  z,u,v,ib,mask=-9999.9,dum;
  float  z_averaged,u_averaged,v_averaged,ib_averaged;
  float  *serie[200][4];

  float  **h,**p;
  float  count;
  double t,t1,t2,tstore;
  double time,sampling=0.0,average=0.0,htime;
  double *tmodel,tt;
  int origine, nodes[200][3],search_status[200],zero=0,first=1;
  int i,j,k;
  int n,status,nst,m,mstore;
  int  year,month;
  int nnde,elt;
  int start_year,last_year,start_month,last_month;
  date_t date;
  int packing;

  FILE *in,*out;
  string h_root, p_root, pathname;
  char *hfile=0, *filename;
  char outname[256]="\0",*default_dir=NULL;
  pressure_station *sample;
  char *keyword,*input,*datafile;
  int filter=0,exact=0;
  meta_archive_t hinfo,pinfo;

  int   nitems;
  long  pos;
  
  status=HarmonicFilter_experiment();
  TRAP_ERR_EXIT(ENOEXEC,"testing\n");

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
        switch (keyword[1]) {
/*------------------------------------------------------------------------------
        root name for default ocean archive */
        case 'a' :
          sscanf(argv[n+1],"%lf",&average); /* en secondes */
          n++;
          n++;
          break;

/*------------------------------------------------------------------------------
        root name for default ocean archive */
        case 'h' :
          h_root=argv[n+1];
          n++;
          n++;
          break;

/*------------------------------------------------------------------------------
        root name for default forcing archive */
        case 'p' :
          p_root=argv[n+1];
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
          pathname=argv[n+1];
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
        filter option (not used) */
        case 'f' :
          filter=1;
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

  if(pathname=="") {
    pathname="./";
    printf("use <./> as root name for default input directorys\n");
    }

  in=fopen(input,"r");
  if(in==NULL) {
    STDOUT_BASE_LINE("cannot open file %s, abort... \n",input);
    exit(-1);
    }
  fscanf(in, "%d ", &nst);

  sample=new pressure_station[nst];
  for(k=0; k<nst;k++) {
    fscanf(in," %lf %lf %s",&sample[k].t,&sample[k].p,sample[k].name);
    }
  fclose(in);

  if(h_root=="") {
    h_root="analysis-YYYY.MM.nc";
    printf("use %s as template name for sea state archive files\n", h_root.c_str());
    }

  size_t position=h_root.find("YYYY");
  if(position==string::npos) printf("\n----> WARNING: code now needs a full filename template (such as analysis.YYYY.MM.nc). You provided %s\n",h_root.c_str());


/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  preliminarly fetch information on model mesh
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  hfile=new char[pathname.size()+h_root.size()+256];
  filename=new char[h_root.size()+256];
  
//   sprintf(hfile,"%s/%s-%4.4d.%2.2d.nc",pathname.c_str(),h_root.c_str(),start_year,start_month);
  
  date.init(start_year,start_month);
  status=decode_name(date, h_root.c_str(), (const char *) 0, &filename, &packing);
  sprintf(hfile,"%s/%s",pathname.c_str(),filename);

  status=archive_info(hfile,1,&hinfo);
  if(status!=0) goto error;
  status=fe_list(&hinfo.mesh);
  if(status!=0) goto error;


/*------------------------------------------------------------------------------
  find FE elements and compute interpolation coefficients*/
  for (k=0;k<nst;k++) {
    sample[k].code=fe_nearest_vertex(hinfo.mesh,sample[k].t, sample[k].p, (int *) NULL, zero);
    printf("station %3d %25s: t=%7.2f, p=%7.2f, node=%7d \n", k+1, sample[k].name,sample[k].t, sample[k].p,sample[k].code);
    search_status[k]=fe_beta(hinfo.mesh, sample[k].t, sample[k].p,&elt,nodes[k],beta[k]);
    }

/*------------------------------------------------------------------------------
  allocate memory for time series */
  if(sampling == 0.) sampling=hinfo.sampling;

  origine=julian_day(1,1,1950);
  t1=(julian_day(start_month,1,start_year)-origine)*24.*3600;
  if(last_month==12)
    t2=(julian_day(1,1,last_year+1)-origine)*24.*3600;
  else
    t2=(julian_day(last_month+1,1,last_year)-origine)*24.*3600;

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

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  truncate existing sample file if needed
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  for (k=0;k<nst;k++) {
    count=0;
    sprintf(outname,"%s/sample.%d",default_dir,k+1);
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

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  start extraction
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  n=-1;
  m=0;
  first=1;
  for (year=start_year;year<=last_year;year++) {
    for (month=1;month<=12;month++) {
      if((year==start_year) && (month < start_month)) continue;
      if((year==last_year)  && (month > last_month))  break;
      
      
//       sprintf(hfile,"%s/%s-%4.4d.%2.2d.nc",pathname.c_str(),h_root.c_str(),year,month);
  
      date.init(year,month);
      status=decode_name(date, h_root.c_str(), (const char *) 0, &filename, &packing);
      sprintf(hfile,"%s/%s",pathname.c_str(),filename);

      status=archive_info(hfile,TUGO_UG_NETCDF_FORMAT,&hinfo);
      if(status!=0) goto error;

      status=fe_list(&hinfo.mesh);
      if(status!=0) goto error;
            

      for (i=0;i<(int) hinfo.nframe;i++) {
/*------------------------------------------------------------------------------
        increment the frame number counted from beginning */
        n++;
        hinfo.flag=FLAG_STATE;
        status=archive_read(hfile,TUGO_UG_NETCDF_FORMAT,hinfo,i,h,&mask,&htime);
        hinfo.flag=FLAG_FORCING;
        status=archive_read(hfile,TUGO_UG_NETCDF_FORMAT,hinfo,i,p,&mask,&htime);
        htime=cnes_time(hinfo.reference,'s')+htime;
/*------------------------------------------------------------------------------
        store time as CNES time */
        tmodel[n]=htime;
        if(status!=0) goto error;
        if(i%48 ==0)
          printf("t= %6.1f hours status= %d %s\n",htime/d2s,status, sgetcnesdate(htime/3600.0));
/*------------------------------------------------------------------------------
        convert to MKS if needed */
        factor=1.;
        for (k=0;k<nst;k++)
          if((exact==0) || (search_status[k] !=0)) {
            serie[k][0][n] =  h[0][sample[k].code]/factor;
            serie[k][1][n] =  h[1][sample[k].code]/factor;
            serie[k][2][n] =  h[2][sample[k].code]/factor;
            serie[k][3][n] = -p[0][sample[k].code]/factor;
            }
          else {
            serie[k][0][n] = (beta[k][0]*h[0][nodes[k][0]] + beta[k][1]*h[0][nodes[k][1]]
                                + beta[k][2]*h[0][nodes[k][2]])/factor;
            serie[k][1][n] = (beta[k][0]*h[1][nodes[k][0]] + beta[k][1]*h[1][nodes[k][1]]
                                + beta[k][2]*h[1][nodes[k][2]])/factor;
            serie[k][2][n] = (beta[k][0]*h[2][nodes[k][0]] + beta[k][1]*h[2][nodes[k][1]]
                                + beta[k][2]*h[2][nodes[k][2]])/factor;
            serie[k][3][n] = -(beta[k][0]*p[0][nodes[k][0]] + beta[k][1]*p[0][nodes[k][1]]
                                + beta[k][2]*p[0][nodes[k][2]])/factor;
            }
        }
      tstore=t;
      mstore=m;
      printf(" sampling = %lf, from t=%lf to %lf \n",sampling,t/(24.0*3600.0),tmodel[n]/(24.0*3600.0));
      for (k=0;k<nst;k++) {
        t=tstore;
        m=mstore;
        sprintf(outname,"%s/sample.%d",default_dir,k+1);
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
            printf("%12.6f %6.3f %6.3f %6.3f\n",time,z, serie[k][0][m],z-serie[k][0][m]);
            }
          if(z != mask) {
            time=t/(24.0*3600.0);
//             fprintf(out,"%12.6f %6.3f %6.3f %6.3f %6.3f\n",time,z,u,v,ib);
            fprintf(out,"%15.9f %6.3f %6.3f %6.3f %6.3f\n",time,z,u,v,ib);
            t+=sampling;
            m++;
            }
          else {
            time=t/(24.0*3600.0);
            printf("%12.6f %6.3f %6.3f %6.3f %6.3f status= %d %s\n",time,z,u,v,ib,status,sgetcnesdate(time));
            }
          }
        fclose(out);
        }
      archive_freeinfo(&hinfo);
/*------------------------------------------------------------------------------
      end of month loop */
      }
/*------------------------------------------------------------------------------
    end of year loop */
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  save extraction informations
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  status=archive_info(hfile,1,&hinfo);
  if(status!=0) goto error;

  out=fopen("sample.info","w");
  for (k=0;k<nst;k++) {
    fprintf(out," %9.3f %9.3f %20s %6d %6.3f\n",sample[k].t,sample[k].p,sample[k].name,sample[k].code,hinfo.mesh.vertices[sample[k].code].h/100.);
    }
  fclose(out);

  archive_freeinfo(&hinfo);

  free(h);
  free(p);
  
//   status=HarmonicFilter_experiment();

  TRAP_ERR_EXIT(0,"%s -computation complete ^^^^^^^^^\n",argv[0]);

 error:
  TRAP_ERR_EXIT(-1,"%s -computation aborted ^^^^^^^^^\n",argv[0]);
}

