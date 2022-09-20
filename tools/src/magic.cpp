#define MAIN_SOURCE

#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
 
#include "tools-structures.h"
    
#include "rutin.h"     /*  rutin.h contains common utility routines  */
#include "archive.h"
#include "fe.h"
#include "poc-time.h"


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int fe_intpl_xyt(char * datafile, double t, double p, double *time, int n, float *serie[3])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  float  z,z1,z2;
  double beta[3];
  float  **h1,**h2;
  float  count,start,end;
  double time1,time2,sampling=0.0,t0;
  int nodes[3];
  int i,j,k,node;
  int status;
  int nndes,elt;
  date_t reference;

  FILE *in;
  meta_archive_t info;
  
  status=clx_archiveinfo(datafile,&info);
  if(status!=0) {__ERR_BASE_LINE__("exiting\n");exit(-1);}

  if(sampling == 0.) sampling=info.sampling;
  reference=info.reference;

  status=fe_list(&info.mesh);
  if(status!=0) {__ERR_BASE_LINE__("exiting\n");exit(-1);}

  t0=(julian_day(reference.month,reference.day,reference.year)-CNES0jd)*24.*3600.+reference.second;

  start=t0+info.start;
  end  =start+info.sampling*(info.nframe-1);

  status=fe_beta(info.mesh, t, p,&elt,nodes,beta);

  /* temporay, for checking only; to be commented in regular use*/

/*   node=fe_nearest_vertex(info.mesh,t,p, (int *) NULL, 0); */
/*   for (j=0;j<3;j++) if(nodes[j]==node) beta[j]=1.; else beta[j]=0.; */

  /* temporay, for checking only; to be commented in regular use*/

  for (j=0;j<3;j++)
    serie[j]= (float *) malloc(n*sizeof(double));

  nndes=info.mesh.nvtxs;

  h1=smatrix(0,2,0,nndes-1);
  h2=smatrix(0,2,0,nndes-1);
  
  in=fopen(datafile,"r");

  for(k=0;k<n;k++) {
    t=time[k]*24.*3600.;
    i=(int) NINT((t-start)/info.sampling);
    status=clx_archiveread(in,info,i+1,h1,&time1);
    status=clx_archiveread(in,info,i+2,h2,&time2);
    if(status != 0) {
      printf("end of file reached... \n");
      break;
      }
    printf("t= %6.1f hours status= %d %s\n",time1/3600.0,status, sgetnewdate(info.reference,time1));
    printf("t= %6.1f hours status= %d %s\n",time2/3600.0,status, sgetnewdate(info.reference,time2));
    t=t-t0;
    for (j=0;j<3;j++) {
      z1=0;
      z2=0;
      for (i=0;i<3;i++) {
        z1+=h1[j][nodes[i]-1]*beta[i];
        z2+=h2[j][nodes[i]-1]*beta[i];
        }
      serie[j][k] = ((time2-t)*z1+(t-time1)*z2)/info.sampling;
      }
    }

  free_smatrix(h1,0,2,0,nndes-1);
  free_smatrix(h2,0,2,0,nndes-1);
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 
  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double dt;
  double t0;
  float *serie[3];
  
  double  lon,lat;
  float  dummy, hours,start,end;
  double *time,t;
  double time1,time2,sampling=0.0;
  int i,j,k,l,L,count;
  int n,status,nst;
  int year,month,day,hour;
  date_t reference,actual;

  FILE *in;
  char *rootname=NULL;
  char *option,*keyword,*s,*datafile=NULL;
  meta_archive_t info;
  
  fct_echo(argc,argv);
 
  n=1;
  while (n < argc)
    {
    keyword=strdup(argv[n]);
    switch (keyword[0])
      {
      case '-':
        switch (keyword[1])
        {
        case 'f' :
          datafile= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'r' :
          rootname= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 's' :
          sscanf(argv[n+1],"%lf",&sampling);
          n++;
          n++;
          break;

        default:
          __OUT_BASE_LINE__("unknown option %s\n",keyword);
          exit(-1);
        }
        break;

      default:
        sscanf(argv[n],"%f",&lon);
        n++;
        sscanf(argv[n],"%f",&lat);
        n++;
        sscanf(argv[n],"%lf",&t);
        n++;
        break;
      }
      free(keyword);
    }

  actual=poctime_getdatecnes(t*24,'h');

  printf("%s\n", sgetcnesdate(t*24));

  if(datafile==NULL) {
    datafile=(char *)malloc(strlen(rootname)+64);
    sprintf(datafile,"%s-%4.4d.%2.2d",rootname,actual.year,actual.month);
    }

  printf("opening %s\n",datafile);
  status=clx_archiveinfo(datafile,&info);
  if(status!=0) goto error;

  if(sampling == 0.) sampling=info.sampling;

  reference=info.reference;

  status=fe_list(&info.mesh);
  if(status!=0) goto error;

  /* compute time origine t0 in cnes time, seconds ellapsed from 1/1/1950 0h GMT */

/*   t0=(julian_day(reference.month,reference.day,reference.year)-CNES0jd)*24.*3600.+reference.second; */

/*   start=t0+info.start; */
/*   end  =t0+info.start+info.sampling*(info.nframe-1); */
  
/*   count=(end-start)/sampling+1; */
/*   for (j=0;j<count;j++) */
/*     time[j]=t0+info.start+j*sampling; */

  count=1;

  time=(double *)malloc(count*sizeof(double));
  time[0]=t;

  for (j=0;j<3;j++)
    serie[j]= (float *) malloc(count*sizeof(float));

  status=fe_intpl_xyt(datafile,lon,lat,time,count,serie);
  if(status!=0) goto error;
  for (i=0; i< count; i++) {
    printf("%12.6f %6.3f %6.3f %6.3f\n",
                      time[i],serie[0][i]/100.0,serie[1][i]/100.0,serie[2][i]/100.0);
    }
  __ERR_BASE_LINE__("exiting\n");exit(0);

 error:
  __ERR_BASE_LINE__("%s -computation aborted ^^^^^^^^^\n",argv[0]);
  exit(-1);
}

