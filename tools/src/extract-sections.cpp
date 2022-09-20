


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


  Extract time series from model's netcdf archives


-----------------------------------------------------------------------------*/

#define MAIN_SOURCE

#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "tools-structures.h"
#include "tides.h"

#include "fe.h"
#include "archive.h"
#include "rutin.h"   /*rutin.h contains common utility routines  */
#include "map.h"
#include "poc-time.h"
#include "functions.h"
#include "poc-netcdf-data.hpp"

extern int  sampler_section_init(mesh_t & mesh, char *SectionInputFile);
extern void sampler_section_update(double t, mesh_t mesh, int, state2d_t state, parameter_t data);

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int process_binary (const char *pathname, const char *h_root, const char *p_root,  const char *default_dir, pressure_station *sample, date_t start, date_t last, int nst, 
		      float scale, int exact, double sampling)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double  beta[1000][3];
  float  z,u,v,ib,mask=-9999.9,dum;
  float z_averaged,u_averaged,v_averaged,ib_averaged;
  double t0;
  float  *serie[1000][4];

  float  **h,**p;
  float  count;
  double t,t1,t2,tstore;
  double time,average=0.0,ptime,htime;
  double *tmodel,tt;

  int nodes[1000][3],search_status[1000],zero=0,first=1;
  int i,j,k,l;
  int n,status,m,mstore;
  int  year,month,day,hour;
  int nnde,elt;
  date_t reference;

  FILE *in,*out,*in1,*in2;
  char *hfile,*pfile;
  char outname[256]="\0";

  int filter=0;
  meta_archive_t hinfo,pinfo;

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

/*---------------------------------------------------------------------
  preliminarly fetch information on model mesh */
  status=clx_archiveinfo(hfile,&hinfo);
  if(status!=0) goto error;
  status=fe_list(&hinfo.mesh);
  if(status!=0) goto error;

/*---------------------------------------------------------------------
  find FE elements and compute interpolation coefficients*/
  for (k=0;k<nst;k++) {
    sample[k].code=fe_nearest_vertex(hinfo.mesh,sample[k].t, sample[k].p, (int *) NULL, zero);
    sample[k].h=hinfo.mesh.vertices[sample[k].code].h;
    printf("station %25s: x=%7.2f, y=%7.2f, h=%7.2f, node=%6d", sample[k].name,sample[k].t, sample[k].p, sample[k].h, sample[k].code);
    search_status[k]=fe_beta(hinfo.mesh, sample[k].t, sample[k].p,&elt,nodes[k],beta[k]);
    if(search_status[k]!=0) {
      printf(" out of domain \n");
      }
    else {
      printf("%s ok \n", sample[k].name);
      }
    }

/*---------------------------------------------------------------------
  allocate memory for time series */
  if(sampling == 0.) sampling=hinfo.sampling;

  t1=(julian_day(start.month,1,start.year)-CNES0jd)*24.*3600;
  if(last.month==12)
    t2=(julian_day(1,1,last.year+1)-CNES0jd)*24.*3600;
  else
    t2=(julian_day(last.month+1,1,last.year)-CNES0jd)*24.*3600;

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

  nnde=hinfo.mesh.nvtxs;
  h=smatrix(0,2,0,nnde-1);
  p=smatrix(0,2,0,nnde-1);

/*---------------------------------------------------------------------
  initialize clock */
  t=hinfo.start+cnes_time(hinfo.reference,'s');

  archive_freeinfo(&hinfo);

/* *---------------------------------------------------------------------
  truncate existing sample file if needed */
  for (k=0;k<nst;k++) {
    count=0;
    sprintf(outname,"%s/sample.%s",default_dir,sample[k].name);
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
  for (year=start.year;year<=last.year;year++) {
    for (month=1;month<=12;month++) {
      if((year==start.year) && (month < start.month)) continue;
      if((year==last.year)  && (month > last.month))  break;

      sprintf(hfile,"%s/%s-%4.4d.%2.2d",pathname,h_root,year,month);
      sprintf(pfile,"%s/%s-%4.4d.%2.2d",pathname,p_root,year,month);

      status=clx_archiveinfo(hfile,&hinfo);
      if(status!=0) {
        printf("loading %s failed\n",hfile);
        goto error;
        }
      status=clx_archiveinfo(pfile,&pinfo);
      if(status!=0) {
        printf("loading %s failed\n",pfile);
        goto error;
        }

/*---------------------------------------------------------------------
      h and p file consistency check */
      if( hinfo.nframe!=pinfo.nframe ) {
         printf("inconsistent frame number in %s (%d) and %s (%d)",hfile,hinfo.nframe,pfile,pinfo.nframe);
        goto error;
        }

      status=fe_list(&hinfo.mesh);
      if(status!=0) {
        printf("archive mesh error\n");
        goto error;
        }

/*---------------------------------------------------------------------
      open archive files */
      in1=fopen(hfile,"r");
      if(in1==NULL) goto error;
      in2=fopen(pfile,"r");
      if(in2==NULL) goto error;

      for (i=0;i<(int) hinfo.nframe;i++) {
/*---------------------------------------------------------------------
        increment the frame number counted from beginning */
        n++;
        status=clx_archiveread(in1,hinfo,i+1,h,&htime);
        status=clx_archiveread(in2,pinfo,i+1,p,&ptime);
        htime=cnes_time(hinfo.reference,'s')+htime;
        ptime=cnes_time(pinfo.reference,'s')+ptime;
/*---------------------------------------------------------------------
        h and p file consistency check */
        if(ptime != htime) goto error;
/*---------------------------------------------------------------------
        store time as CNES time */
        tmodel[n]=htime;
        if(status!=0) goto error;
        if(i%48 ==0)
          printf("t= %6.1f hours status= %d %s\n",htime/3600.0/24.0,status, sgetcnesdate(htime/3600.0));
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
/*---------------------------------------------------------------------
      end of month loop */
      }
/*---------------------------------------------------------------------
    end of year loop */
    }

  status=clx_archiveinfo(hfile,&hinfo);
  if(status!=0) goto error;

  archive_freeinfo(&hinfo);

  free(h);
  free(p);

  return(0);

 error:
  return(-1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int process_netcdf (const char *pathname, const char *name_template, const char *default_dir, pressure_station *sample, date_t start, date_t last, int nst, 
		      int discretisation, float scale, int exact, double sampling)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int k,n;
  int month, year, element, *search_status,zero=0;
  char *input;
  mesh_t mesh;
  discretisation_t descriptor;
  int **nodes;
  double **beta;

  vector<float>  **serie;
  vector<double> tmodel;
  poc_global_t global;
  poc_data_t<double> h,u,v,p,time;
  date_t current;
  double t, tag;
  char *hname="elevation",*uname="ubar",*vname="vbar",*pname="Pa",*tname="time";

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
  
/*---------------------------------------------------------------------
  we need cnes time in seconds */
  t=time.data[0];
  
  tag=t;
  
// //  t=cnes_time(hinfo.reference,'s');
//   status=truncate (default_dir, sample, nst, t);
  
/*---------------------------------------------------------------------
  find FE elements and compute interpolation coefficients*/
  for (k=0;k<nst;k++) {
    sample[k].code=fe_nearest_vertex(mesh, sample[k].t, sample[k].p, (int *) NULL, zero);
    sample[k].h=mesh.vertices[sample[k].code].h;
    printf("station %25s: x=%7.2f, y=%7.2f, h=%7.2f, node=%6d ", sample[k].name,sample[k].t, sample[k].p, sample[k].h, sample[k].code);
    search_status[k]=fe_beta(mesh, descriptor, sample[k].t, sample[k].p,&element,nodes[k],beta[k]);
    if(search_status[k]!=0) {
      printf("%s out of domain\n", sample[k].name);
      }
    else {
      printf("%s ok\n", sample[k].name);
      }
    }

  
  for (year=start.year;year<=last.year;year++) {
    for (month=1;month<=12;month++) {
      if((year==start.year) && (month < start.month)) continue;
      if((year==last.year)  && (month > last.month))  break;
      date_t current(year,month,1,0.);
      status=decode_name(current, 0, name_template, &input);
      status=poc_inq(input,&global);
  
      if(status!=0) break;

      status=h.init(global,hname);
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
        for (k=0;k<nst;k++) {
          if((exact==0) || (search_status[k] !=0)) {
            size_t node=sample[k].code;
            serie[k][0].push_back( h.data[node]*scale);
            serie[k][1].push_back( u.data[node]*scale);
            serie[k][2].push_back( v.data[node]*scale);
            if(p.data!=0) serie[k][3].push_back(-p.data[node]*scale);
            else serie[k][3].push_back(0.);
            }
          }
        }
      for (k=0;k<nst;k++) for(int l=0;l<4;l++) serie[k][l].clear();
      h.destroy_data();
      u.destroy_data();
      v.destroy_data();
      p.destroy_data();
      tmodel.clear();
      delete[] input;
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

 error:
  return(-1);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double factor;
  float  z,u,v,ib,mask=-9999.9,dum;
  float  z_averaged,u_averaged,v_averaged,ib_averaged;
  double t0;
  float  *serie[200][4];

  float  **h,**p;
  float  count;
  double t,t1,t2,start,tstore;
  double time,sampling=0.0,average=0.0,ptime,htime;
  double *tmodel;
  int first=1, frame;
  int i,j,k,l;
  int n,status,nst,m,mstore;
  int year,month,day,hour;
  int nnde;
  int start_year,last_year,start_month,last_month;
  date_t reference;

  size_t size;
  char *h_root=NULL,*p_root=NULL,*pathname=0,*hfile,*pfile;
  char outname[256]="\0",*rootname=NULL,*default_dir=NULL,directory[256]="\0";
  char *option,*keyword,*s,*input,*datafile;
  int filter=0,exact=0;
  meta_archive_t hinfo,pinfo;

  int   nitems;
  long  pos;
  
  state2d_t state;
  parameter_t data;
  
  int paire2D=NCP1xLGP1;
  
  discretisation_t z_descriptor, u_descriptor;
  int z_discretisation, u_discretisation;
  
  fct_echo(argc, argv);
  
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

/*---------------------------------------------------------------------
        extraction mode */
//         case 't' :
//           exact=1;
//           n++;
//           break;


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
    
//   if(format==NULL) {
//     status=process_binary (pathname, h_root, p_root, default_dir, sample, start, last, nst, scale, exact, sampling);
//     }
//   else {
//     status=process_netcdf (pathname, convention, default_dir, sample, start, last, nst, LGP1, scale, exact, sampling);
//     }
//   if(status!=0) goto error;


  hfile=new char[strlen(pathname)+256];
  sprintf(hfile,"%s/%s-%4.4d.%2.2d.nc",pathname,h_root,start_year,start_month);

/*---------------------------------------------------------------------
  preliminarly fetch information on model mesh */
  status=archive_info(hfile,1,&hinfo);
  if(status!=0) goto error;
  status=fe_list(&hinfo.mesh);
  if(status!=0) goto error;
  status= init_edge_table(&hinfo.mesh);
  if(status!=0) goto error;

  status= sampler_section_init(hinfo.mesh, input);
  
/*---------------------------------------------------------------------
  allocate memory for time series */
  if(sampling == 0.) sampling=hinfo.sampling;

  t1=(julian_day(start_month,1,start_year)-CNES0jd)*24.*3600;
  if(last_month==12)
    t2=(julian_day(1,1,last_year+1)-CNES0jd)*24.*3600;
  else
    t2=(julian_day(last_month+1,1,last_year)-CNES0jd)*24.*3600;

  count=(float)(t2-t1)/hinfo.sampling;

  tmodel= new double[(int) count];

  status=paire_discretisation_id(paire2D, &z_discretisation,&u_discretisation);
  status= discretisation_init(&hinfo.mesh, z_discretisation, 0);
  status= discretisation_init(&hinfo.mesh, u_discretisation, 0);
  
  z_descriptor=get_descriptor(hinfo.mesh,z_discretisation);
  u_descriptor=get_descriptor(hinfo.mesh,u_discretisation);
  
  nnde=MAX(z_descriptor.nnodes,u_descriptor.nnodes);
  
  h=smatrix(0,2,0,nnde-1);
  p=smatrix(0,2,0,nnde-1);

/*---------------------------------------------------------------------
  initialize clock */
  t=hinfo.start+cnes_time(hinfo.reference,'s');

  archive_freeinfo(&hinfo);

/*---------------------------------------------------------------------
  start extraction */
  frame=-1;
  m=0;
  first=1;
  for (year=start_year;year<=last_year;year++) {
    for (month=1;month<=12;month++) {
      if((year==start_year) && (month < start_month)) continue;
      if((year==last_year)  && (month > last_month))  break;

      sprintf(hfile,"%s/%s-%4.4d.%2.2d.nc",pathname,h_root,year,month);

      status=archive_info(hfile,1,&hinfo);
      if(status!=0) goto error;
      status=fe_list(&hinfo.mesh);
      if(status!=0) goto error;
      status= init_edge_table(&hinfo.mesh);
       if(status!=0) goto error;


      for (i=0;i<(int) hinfo.nframe;i++) {
/*---------------------------------------------------------------------
        increment the frame number counted from beginning */
        frame++;
        status=archive_read(hfile,TUGO_UG_NETCDF_FORMAT,hinfo,i,h,&mask,&htime);
        if(status!=0) {
          goto error;
          }
        htime=cnes_time(hinfo.reference,'s')+htime;
/*---------------------------------------------------------------------
        store time as CNES time */
        tmodel[frame]=htime;
        if(i%48 ==0)
          printf("t= %6.1f hours status= %d %s\n",htime/3600.0/24.0,status, sgetcnesdate(htime/3600.0));
/*---------------------------------------------------------------------
        convert to MKS if needed */
        factor=1.;
        state.colocation=0;
        poc_copy(state.z,h[0],z_descriptor.nnodes);
        poc_copy(state.H,h[0],z_descriptor.nnodes);
        poc_copy(state.u,h[1],u_descriptor.nnodes);
        poc_copy(state.v,h[2],u_descriptor.nnodes);
        for(n=0;n<hinfo.mesh.nvtxs;n++) {
          state.H[n]+=hinfo.mesh.vertices[n].h;
          }
        sampler_section_update(tmodel[frame], hinfo.mesh, paire2D, state, data);
        state.destroy();
        }
      tstore=t;
      mstore=m;
      printf(" sampling = %lf, from t=%lf to %lf \n",sampling,t/(24.0*3600.0),tmodel[n]/(24.0*3600.0));
      archive_freeinfo(&hinfo);
/*---------------------------------------------------------------------
      end of month loop */
      }
/*---------------------------------------------------------------------
    end of year loop */
    }

  status=archive_info(hfile,1,&hinfo);
  if(status!=0) {
    goto error;
    }

  archive_freeinfo(&hinfo);

  free(h);
  free(p);

  __ERR_BASE_LINE__("%s -computation complete ^^^^^^^^^\n",argv[0]);
  exit(0);

 error:
  __ERR_BASE_LINE__("%s -computation aborted ^^^^^^^^^\n",argv[0]);
  exit(-1);
}

