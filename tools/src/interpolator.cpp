
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

#include "tools-structures.h"

#include "rutin.h"     /*  rutin.h contains common utility routines  */
#include "fe.h"
#include "poc-time.h"
#include "archive.h"

typedef struct {
  double t,p;
  float h;
  int   elt,node[6];
  double beta[6];
  mesh_t *mesh;
  } track_t;


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void fe_intpl_xyt(char *datafile, double t, double p, double *time, int n, float *serie[3])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  float  z,z1,z2;
  double beta[3];
  float  **h1,**h2;
  float  count,start,end;
  double time1,time2,sampling=0.0,t0;
  int nodes[3];
  int i,j,k;
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

  t0=(julian_day(reference.month,reference.day,reference.year)-CNES0jd)*24.*3600.+reference.second;

  start=t0+info.start;
  end  =start+info.sampling*(info.nframe-1);

  status=fe_beta(info.mesh, t, p,&elt,nodes,beta);

  for (j=0;j<3;j++)
    serie[j]= (float *) malloc(n*sizeof(double));

  nndes=info.mesh.nvtxs;

  h1=smatrix(0,2,0,nndes-1);
  h2=smatrix(0,2,0,nndes-1);

  in=fopen(datafile,"r");

  for(k=0;k<n;k++) {
    n++;
    t=time[k];
    i=(int)( (time[k]-start)/info.sampling );
    status=clx_archiveread(in,info,i+1,h1,&time1);
    status=clx_archiveread(in,info,i+2,h2,&time2);
    if(status != 0) {
      printf("end of file reached... \n");
      break;
      }
    printf("t= %6.1f hours status= %d %s\n",time1/3600.0,status, sgetnewdate(info.reference,time1));
    for (j=0;j<3;j++) {
      z1=0;
      z2=0;
      for (i=0;i<3;i++) {
        z1+=h1[j][nodes[i]]*beta[i];
        z2+=h2[j][nodes[i]]*beta[i];
        }
      serie[j][n] = ((time2-t)*z1+(t-time1)*z2)/info.sampling;
      }
    }

  free_smatrix(h1,0,2,0,nndes-1);
  free_smatrix(h2,0,2,0,nndes-1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  float z,z1,z2,beta[3];
  int zero=0;
  double dt,cut,mask=0;
  double t0,dT,pulsation;
  double *serie[200][4];

  float  **h1,**h2;
  float  *u,*v,*buffer;
  float  count,dummy, hours,start,end,w,dum,t;
  double *time;
  double time1,time2,sampling=0.0;
  double sum[2];
  int nodes[3];
  int fmt,i,j,k,l,L,ndum,rstatus,lfdepth,mfdepth;
  int minstep,maxstep,l1,l2,n,status,nst;
  int  year,month,day,hour;
  int nnde,nstep,elt;
  date_t reference;

  size_t size,offset;
  FILE *in,*out,*sout[100];
  FILE *mesh;
  char file1[256],file2[256], wocefile[256], model[256];
  char fout[100][256], outname[256]="\0",*rootname=NULL;
  pressure_station *sample;
  track_t *track;
  char *option,*keyword,*s,meshfile[256]="\0",*input,*datafile;
  int filter=0;
  meta_archive_t info;
  statistic_t stath;

  grid_t *grid;
  pressure_station *setptr;

  fct_echo(argc,argv);
 
  n=1;
  while (n < argc) {
    keyword=strdup(argv[n]);
    switch (keyword[0]) {
      case '-':
        switch (keyword[1]) {
        case 'h' :
          option= strdup(keyword);
          n++;
          break;

        case 'p' :
          option= strdup(keyword);
          n++;
          break;

        case 'm' :
          option= strdup(meshfile);
          n++;
          break;

        case 'i' :
          input= strdup(argv[n+1]);
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

        case 'f' :
          filter=1;
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

  L=strlen(datafile);
  if(strrchr(datafile,(int) '/') != NULL) {
    l=strlen(strrchr(datafile,(int) '/'));
    strncpy(rootname,datafile,L-l+1);
    }
  else  strcpy(rootname,"./");

  printf("rootname: %s \n",rootname);

  in=fopen(input,"r");
  if(in==NULL) {__ERR_BASE_LINE__("exiting\n");exit(-1);}
  fscanf(in, "%s ", meshfile);
  fscanf(in, "%d ", &nst);

  sample=(pressure_station *)malloc(nst*sizeof(pressure_station));
  track=(track_t *)malloc(nst*sizeof(track_t));

  for(k=0; k<nst;k++) {
    fscanf(in," %f %f %s",&sample[k].t,&sample[k].p,sample[k].name);
    printf(" %f %f %s \n",sample[k].t,sample[k].p,sample[k].name);
    }

  fclose(in);

  status=clx_archiveinfo(datafile,&info);
  if(status!=0) {__ERR_BASE_LINE__("exiting\n");exit(-1);}

  if(sampling == 0.) sampling=info.sampling;

  count=info.nframe;
  reference=info.reference;

  status=fe_list(&info.mesh);

/*   status=fe_readmesh(meshfile,MESH_FILE_FORMAT_TRIGRID,&info.mesh); */

/*   status=fe_list(&info.mesh); */

  /* compute time origine t0 in cnes time, seconds ellapsed from 1/1/1950 0h GMT */
  t0=(julian_day(reference.month,reference.day,reference.year)-CNES0jd)*24.*3600.+reference.second;

  start=t0+info.start;
  end  =t0+info.start+info.sampling*(info.nframe-1);

  time=(double *)malloc(info.nframe*sizeof(double));

  for (j=0;j<info.nframe;j++)
    time[j]=t0+info.start+j*sampling;

  count=(end-start)/sampling+1;

  for (k=0;k<nst;k++) {
    sample[k].code=fe_nearest_vertex(info.mesh,sample[k].t, sample[k].p,  NULL, zero);
    status=fe_beta(info.mesh, sample[k].t, sample[k].p,&track[k].elt,track[k].node,track[k].beta);
    printf("station %15s: t=%7.2f, p=%7.2f, node=%6d, elts=%6d \n",
           sample[k].name,sample[k].t, sample[k].p,sample[k].code, track[k].elt);
//     printf(" %6d %6d %6d %8.6f %8.6f %8.6f \n", nodes[0],nodes[1],nodes[2],
//                                                 beta[0],beta[1],beta[2]);
    }

  for (k=0;k<nst;k++)
    for (j=0;j<4;j++)
      serie[k][j]= (double *) malloc( ((int) count) *sizeof(double));

  nnde=info.mesh.nvtxs;

  h1=smatrix(0,2,0,nnde-1);
  h2=smatrix(0,2,0,nnde-1);
  buffer=(float *)malloc(3*nnde*sizeof(float));

  in=fopen(datafile,"r");
  n=-1;
  while(t<end) {
    n++;
    t=info.start+n*sampling;
    i=(int)( n*sampling/info.sampling );
    status=clx_archiveread(in,info,i+1,h1,&time1);
    status=clx_archiveread(in,info,i+2,h2,&time2);
    if(status != 0) {
      printf("end of file reached... \n");
      break;
      }
      printf("t= %6.1f hours status= %d %s\n",time1/3600.0,status, sgetnewdate(info.reference,time1));
    for (k=0;k<nst;k++) {
      for (j=0;j<3;j++) {
        z1=0;
        z2=0;
        for (i=0;i<3;i++) {
          z1+=h1[j][track[k].node[i]]*track[k].beta[i];
          z2+=h2[j][track[k].node[i]]*track[k].beta[i];
          }
        serie[k][j][n] = ((time2-t)*z1+(t-time1)*z2)/info.sampling;
        }
      }
    }

  for (k=0; k<nst; k++) {
    strcpy(outname,rootname);
    strcat(outname,sample[k].name);
    if(strcmp(option,"-p")==0) strcat(outname,".p");
    else strcat(outname,".h");
    if (filter) {
/*       stath=filtering(outname,serie[k][0], &dt,mask, n, dt); */
      }
    else {
      sout[k]=fopen(outname,"w");
      for (i=0; i< n; i++) {
        t=(i*sampling+t0+info.start)/(24.0*3600.0);
        fprintf(sout[k],"%12.6f %6.3f %6.3f %6.3f\n",
                         t,serie[k][0][i]/100.0,serie[k][1][i]/100.0,serie[k][2][i]/100.0);
        }
      fclose(sout[k]);
      }
    }

  free_smatrix(h1,0,2,0,nnde-1);
  free_smatrix(h2,0,2,0,nnde-1);

  __ERR_BASE_LINE__("%s -computation complete ^^^^^^^^^\n",argv[0]);
  exit(0);
}

