#define MAIN_SOURCE

# include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
 
#include "tools-structures.h"
    
#include "poc-time.h"
#include "fe.h"
#include "map.h"
#include "rutin.h"
#include "archive.h"
#include "functions.h"

//define after the main function
extern void load_xyzt(double *lon, double *lat, int *mes, double **sealevel, double **time, FILE *fic_data);
extern void write_list_header(FILE *in,FILE *out,int *nst);


/*---------------------------------------------------------------------------------*/

void get_mean  (meta_archive_t info, char *meanfilename1, char *meanfilename2,float **ms, float **mp, int nst, pressure_station *sample)

/*---------------------------------------------------------------------------------*/
{
  int *search_status,status,i,j,nitems;
  int **nodes,nnde,elt;

  float *MEAN_PA, *MEAN_SL;
  double **beta;

  char line[300];

  FILE *in;



  search_status=(int *)calloc(nst,sizeof(int));
  beta=(double **)malloc(nst*sizeof(double *));
  nodes=(int **)malloc(nst*sizeof(int *));
  for(i=0;i<nst;i++) nodes[i]=(int *)malloc(3*sizeof(int));
  for(i=0;i<nst;i++) beta[i]=(double *)malloc(3*sizeof(double));

  if( (status=fe_list(&info.mesh)) !=0) __TRAP_ERR_EXIT__(status,"can not access mesh infos\n");
  for (i=0;i<nst;i++)  search_status[i]=fe_beta(info.mesh, sample[i].t, sample[i].p,&elt,nodes[i],beta[i]);

  nnde=info.mesh.nvtxs;
  MEAN_PA=(float *)malloc(nnde*sizeof(float));
  MEAN_SL=(float *)malloc(nnde*sizeof(float));
  *mp=(float *) calloc(nst,sizeof(float));
  *ms=(float *) calloc(nst,sizeof(float));

  if( (in=fopen(meanfilename1,"r"))== NULL ) __TRAP_ERR_EXIT__(errno,"can't open mean pressure file %s (%s)\n",meanfilename1,strerror(errno));
  fgets(line,sizeof(line),in);
  fgets(line,sizeof(line),in);
  for(i=0;i<nnde;i++) {
    nitems=fscanf(in,"%d %f",&j,&(MEAN_PA[i]));
    if(nitems != 2) __TRAP_ERR_EXIT__(-1,"can not read mean pressure file");
    }
  fclose(in);
  printf("# Got the mean pressure ... OK\n");

  if( (in=fopen(meanfilename2,"r"))== NULL ) __TRAP_ERR_EXIT__(errno,"can't open mean sea level file %s (%s)\n",meanfilename2,strerror(errno));
  fgets(line,sizeof(line),in);
  fgets(line,sizeof(line),in);
  for(i=0;i<nnde;i++) {
    nitems=fscanf(in,"%d %f",&j,&MEAN_SL[i]);
    if(nitems != 2) __TRAP_ERR_EXIT__(-1,"can not read mean sea level file %s\n",meanfilename2);
    }
  fclose(in);
  printf("# Got the mean sea level ... OK\n");


  for (j=0;j<nst;j++) {
      if(MEAN_PA != NULL)
        for (i=0;i<3;i++)
          (*mp)[j]+=(MEAN_PA[nodes[j][i]]/100)*(beta[j][i]);
      else (*mp)[j]=0.;

      if(MEAN_SL != NULL)
        for (i=0;i<3;i++)
          (*ms)[j]+=(MEAN_SL[nodes[j][i]]/100)*beta[j][i];
      else  (*ms)[j]=0.;
    }

  free(MEAN_PA);
  free(MEAN_SL);

} /* end get_mean */



/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
 
  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double time,*tmodel,t;
  double t0,t1,t2;
  double **data_sealevel,**data_time,*lon,*lat,***corr;

  float **h,**p,ib,z,*mp,*ms;
  double  **beta;
  float  ***serie;

  int nitems,hourly=0,mss=0;
  int i,j,k,count,n,m,status,nst,mask=999999;
  int year,month,day,hour;
  int start_year,last_year,start_month,last_month;
  int nnde,*mes,*cycle;
  int fe_nlistmax=1000;
  int *search_status,zero=0,**nodes,elt,first=1,start;

  int *used, *nodelist, *noderefc, nodecount=0;
  float dum;

  date_t reference,actual,current,current2,start_date;

  FILE *in1,*in2,*out;

  char *rootname_h=NULL,*rootname_p=NULL,*mean_file=NULL,*mean_fileSL=NULL,*default_dir;
  char *option,*keyword,*s,*datafile_h=NULL,*datafile_p=NULL,*input1,*input2,*output;
  char carac,tmp_char[300];
  char filename[300];
  char line[300];
  char *justfilename,*search,dirname[300],rootname[300];

  meta_archive_t info;

  memory_t memory_h,memory_p;

  pressure_station *sample;

  const char *w1_format={"%12.6f %7.4f %7.4f 99.9999 99.9999 99.9999 99.9999 99.9999 %7.4f 99.9999 99.9999 99.9999 99.9999 99.9999 99.9999 99.9999 99.9999 99.9999 99.9999 99.9999 99.9999 99.9999 99.9999 99.9999 99.9999\n"};
  const char *w2_format={"%3d %7.4f %7.4f %12.6f %7.4f %7.4f 99.9999 99.9999 99.9999 99.9999 99.9999 %7.4f 99.9999 99.9999 99.9999 99.9999 99.9999 99.9999 99.9999 99.9999 99.9999 99.9999 99.9999 99.9999 99.9999 99.9999 99.9999 99.9999 %7.4f %7.4f\n"};


/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/* run parameters and inputs */
/*----------------------------------------------------------------------------*/
  
  fct_echo(argc,argv);
 
  n=1;
  while (n < argc)
    {
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

      if(strcmp(argv[n],"-hourly") == 0) {
          hourly=1;
          n++;
          continue;
        }

      if(strcmp(argv[n],"-mss") == 0) {
          mss=1;
          n++;
          continue;
        }

      keyword=strdup(argv[n]);
      switch (keyword[0])
        {
        case '-':
          switch (keyword[1])
            {
            case 'f' :/* data filename */
              s= strdup(argv[n+1]);
              n++;
              n++;
              sscanf(s,"%s",filename);
              break;

            case 'h' :/* root name for 'analysis' file */
              rootname_h= strdup(argv[n+1]);
              n++;
              n++;
              break;

            case 'm' :/* file names for mean sea level and pressure */
              mean_file= strdup(argv[n+1]);
              mean_fileSL= strdup(argv[n+2]);
              n++;
              n++;
              n++;
              break;
          
            case 'p' :/* root name for 'forcing' file */
              rootname_p= strdup(argv[n+1]);
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
          break;
        }
      free(keyword);
    }

  default_dir=strdup("./");
  
  search=strrchr(filename,'.');
  n=strlen(filename)-strlen(search);
  strncpy(rootname,filename,n);
  rootname[n]='\0';
  
  if((in1=fopen(filename,"r"))==NULL) __TRAP_ERR_EXIT__(errno,"can't open data file %s (%s)\n",filename,strerror(errno));
  
  output=(char *)malloc(1024);
  sprintf(output,"%s.corr.list",rootname);
  if( (out=fopen(output,"w")) ==NULL) __TRAP_ERR_EXIT__(errno,"can't open output file %s (%s)\n",output,strerror(errno));
 
  printf("\n\n%s: starting execution ...\n\n",argv[0]);

  if(mss==0)  write_list_header(in1,out,&nst);
  else
    {
      sprintf(line,"wc -l %s  >wc.tmp",filename);
      system(line);
      in2=fopen("wc.tmp","r");
      fscanf(in2,"%d",&nst);
      fclose(in2);
      sprintf(line,"rm wc.tmp");
      system(line);

      do
        {
          fgets(line,200,in1);
          fprintf(out,"%s",line);
        } while (strncmp(line,"#----------------------------------------",20) != 0 );
    }

  /*----------------------------------------------------*/
  /* getting coordinates of all altimetric measurements */
  /*----------------------------------------------------*/

  printf("\n>Getting stations coordinates ...\n\n");

  lon=(double *) calloc(nst,sizeof(double));
  lat=(double *) calloc(nst,sizeof(double));
  mes=(int *) calloc(nst,sizeof(int));
  data_sealevel=(double **) calloc(nst,sizeof(double));
  data_time=(double **) calloc(nst,sizeof(double));

  if(mss==0)
    for (i=0;i<nst;i++) load_xyzt(&(lon[i]),&(lat[i]),&(mes[i]),&(data_sealevel[i]),&(data_time[i]),in1);
  else
    {
      for(i=0;i<nst;i++) mes[i]=1;

      cycle=(int *) calloc(nst,sizeof(int));
      corr=(double ***)malloc(nst*sizeof(double **));
      for (i=0;i<nst;i++) {
          corr[i]=(double **)calloc(mes[i],sizeof(double *));
          for (j=0;j<mes[i];j++) corr[i][j]=(double *)calloc(25,sizeof(double));
        }
      for(i=0;i<nst;i++) data_sealevel[i]=(double *)calloc(mes[i],sizeof(double));
      for(i=0;i<nst;i++) data_time[i]=(double *)calloc(mes[i],sizeof(double));

     for (i=0;i<nst;i++) {
         fgets(line,sizeof(line),in1);
         for (j=0;j<mes[i];j++)
           n=sscanf(line,"%3d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                    &(cycle[i]),
                    &(lon[i]),               /* longitude */
                    &(lat[i]),               /* latitude */
                    &(data_time[i][j]),        /* cnes time in days */
                    &(data_sealevel[i][j]), /* ssh */
                    &(corr[i][j][0]),  /* MOG2D-G model sea level */
                    &(corr[i][j][1]),  /* MOG2D-G model S1 atmospheric tide */
                    &(corr[i][j][2]),  /* MOG2D-G model S2 atmospheric tide */
                    &(corr[i][j][3]),  /* MOG2D-medsea model sea level */
                    &(corr[i][j][4]),  /* MOG2D-medsea model S1 atmospheric tide */
                    &(corr[i][j][5]),  /* MOG2D-medsea model S2 atmospheric tide */
                    &(corr[i][j][6]),  /* inverted barometer */
                    &(corr[i][j][7]),  /* inverted barometer S1 atmospheric tide */
                    &(corr[i][j][8]),  /* inverted barometer S2 atmospheric tide */
                    &(corr[i][j][9]),  /* loading effects */
                    &(corr[i][j][10]), /* solid Earth tide */
                    &(corr[i][j][11]), /* MOG2D-medsea model ocean tide */
                    &(corr[i][j][12]), /* MOG2D-medsea model S1 ocean tide */
                    &(corr[i][j][13]), /* MOG2D-medsea model S2 ocean tide */
                    &(corr[i][j][14]), /* FES99 model ocean tide */
                    &(corr[i][j][15]), /* FES99 model S1 ocean tide */
                    &(corr[i][j][16]), /* FES99 model S2 ocean tide */
                    &(corr[i][j][17]), /* FES2002 model ocean tide */
                    &(corr[i][j][18]), /* FES2002 model S1 ocean tide */
                    &(corr[i][j][19]), /* FES2002 model S2 ocean tide */
                    &(corr[i][j][20]), /* Harmonic analysis ocean tide */
                    &(corr[i][j][21]), /* Harmonic analysis S1 ocean tide */
                    &(corr[i][j][22]), /* Harmonic analysis S2 ocean tide */
                    &(corr[i][j][23]), /* mss */
                    &(corr[i][j][24]));/* geoid */
         if(n!=30) __TRAP_ERR_EXIT__(-1,"error in reading data file");
       }
    }
   fclose(in1);

  sample=(pressure_station *)malloc(nst*sizeof(pressure_station));
  
  for(i=0; i<nst;i++) {
      sample[i].t=lon[i];
      sample[i].p=lat[i];
      sprintf(sample[i].name,"TP-%d",i+1);
    }

  free(lon);
  free(lat);


  /*--------------------------------------------*/
  /* importing mean pressure and mean sea level */
  /*--------------------------------------------*/

  printf("\n>Importing mean levels from global solution ...\nfiles: %s\n       %s\n\n",mean_file,mean_fileSL);

  info.mesh.triangles=NULL;
  info.mesh.vertices=NULL;
  memory_h.activated=0;
  memory_p.activated=0;
  memory_h.file=NULL;
  memory_p.file=NULL;
  
  if(datafile_h==NULL) datafile_h=(char *)malloc(strlen(rootname_h)+64);
  if(datafile_p==NULL) datafile_p=(char *)malloc(strlen(rootname_p)+64);

  archive_freeinfo(&info);

  if(hourly) sprintf(datafile_h,"%s-%4.4d.%2.2d",rootname_h,start_year,start_month);
  else sprintf(datafile_h,"%s-%4.4d.%2.2d.reduced-6",rootname_h,start_year,start_month);

  if ( status=clx_archiveinfo(datafile_h,&info) != 0 ) __TRAP_ERR_EXIT__(-1,"can't access archive infos of %s\n",datafile_h);
 
/*------------------------------------------------------------------------------
  get mean level (sea level and ib) to keep consistent wit usual MSS*/
  get_mean(info,mean_file,mean_fileSL,&ms,&mp,nst,sample);

  archive_freeinfo(&info);
 
  /*------------*/
  /* extraction */
  /*------------*/

  printf("\n>Beginning extraction ...\n\n");

  if(hourly) sprintf(datafile_h,"%s-%4.4d.%2.2d",rootname_h,start_year,start_month);
  else sprintf(datafile_h,"%s-%4.4d.%2.2d.reduced-6",rootname_h,start_year,start_month);

/*------------------------------------------------------------------------------
  open archive file and return mesh information*/
  status=clx_archiveinfo(datafile_h,&info);
  if(status!=0) goto error;
 
/*------------------------------------------------------------------------------
  build the element list from the neighbour list*/
  status=fe_list(&info.mesh);
  if(status!=0) goto error;

  used=(int *)malloc(info.mesh.nvtxs*sizeof(int));
  for (i=0;i<info.mesh.nvtxs;i++) used[i]=0;

/*------------------------------------------------------------------------------
  seek nodes and elements*/
  search_status=(int *)calloc(nst,sizeof(int));
/*------------------------------------------------------------------------------
  Warning, recent change: beta is double !!!*/
  beta =(double **)malloc(nst*sizeof(double *));
  nodes=(int **)malloc(nst*sizeof(int *));
  for (i=0;i<nst;i++) {
    nodes[i]=(int *)    malloc(3*sizeof(int));
    beta[i] =(double *) malloc(3*sizeof(double));
    sample[i].code=fe_nearest_vertex(info.mesh,sample[i].t, sample[i].p, (int *) NULL, zero);
    printf("station %15s: t=%7.2f, p=%7.2f, node=%6d \n", sample[i].name,sample[i].t, sample[i].p,sample[i].code);
    search_status[i]=fe_beta(info.mesh, sample[i].t, sample[i].p,&elt,nodes[i],beta[i]);
    for(j=0;j<3;j++) used[nodes[i][j]]=1;
    }

  t0=(julian_day(info.reference.month,info.reference.day,info.reference.year)-CNES0jd)*24.*3600.+info.reference.second;

  archive_freeinfo(&info);

/*------------------------------------------------------------------------------
  count used mesh nodes*/
  for (i=0;i<info.mesh.nvtxs;i++) if(used[i]==1) nodecount++;
  nodelist=(int *)malloc(nodecount*sizeof(int));
  noderefc=(int *)malloc(info.mesh.nvtxs*sizeof(int));

/*------------------------------------------------------------------------------
  create used nodes list*/
  nodecount=0;
  for (i=0;i<info.mesh.nvtxs;i++)
    if(used[i]==1) {
      nodelist[nodecount]=i;
      noderefc[i]=nodecount;
      nodecount++;
      }
  free(used);

  t1=(julian_day(start_month,1,start_year)-CNES0jd)*24.*3600;
  if(last_month==12) t2=(julian_day(1,1,last_year+1)-CNES0jd)*24.*3600;
  else t2=(julian_day(last_month+1,1,last_year)-CNES0jd)*24.*3600;

  count=(int)( (t2-t1)/info.sampling );

  tmodel= (double *) malloc( ((int) count) *sizeof(double));

  serie=(float ***)malloc(nodecount*sizeof(float **));
  for (i=0;i<nodecount;i++) {
    serie[i]=(float **)malloc(2*sizeof(float *));
    for (j=0;j<2;j++) {
      serie[i][j]= (float *) malloc( ((int) count) *sizeof(float));
      if( serie[i][j]==NULL) {
        printf("i=%d nodecount=%d\n",i,nodecount);
        printf("requested size %d\n",(i*2*count+j)*sizeof(float));
        __TRAP_ERR_EXIT__(-1,"memory allocation failure");
        }
      }
    }

  n=-1;
  first=1;
  for (year=start_year;year<=last_year;year++) {
    for (month=1;month<=12;month++) {
      if((year==start_year) && (month < start_month)) continue;
      if((year==last_year)  && (month > last_month))  break;
      if(hourly) {
        sprintf(datafile_h,"%s-%4.4d.%2.2d",rootname_h,year,month);
        sprintf(datafile_p,"%s-%4.4d.%2.2d",rootname_p,year,month);
         }
       else
         {
         sprintf(datafile_h,"%s-%4.4d.%2.2d.reduced-6",rootname_h,year,month);
         sprintf(datafile_p,"%s-%4.4d.%2.2d.reduced-6",rootname_p,year,month);
         }

       status=clx_archiveinfo(datafile_h,&info);
       if(status!=0) goto error;

       count=info.nframe;
       reference=info.reference;

       status=fe_list(&info.mesh);
       if(status!=0) goto error;

       nnde=info.mesh.nvtxs;
       h=smatrix(0,2,0,nnde-1);
       p=smatrix(0,2,0,nnde-1);
  
       if( (in1=fopen(datafile_h,"r")) == NULL ) __TRAP_ERR_EXIT__(errno,"can't open datafile_h %s (%s)\n",datafile_h,strerror(errno));
       if( (in2=fopen(datafile_p,"r")) == NULL ) __TRAP_ERR_EXIT__(errno,"can't open datafile_p %s (%s)\n",datafile_p,strerror(errno));
       for (i=0;i<(int) info.nframe;i++) {
         n++;
         status=clx_archiveread(in2,info,i+1,p,&tmodel[n]);
         status=clx_archiveread(in1,info,i+1,h,&tmodel[n]);
/* 	 status=clx_archiveread(in2,info,i+1,p,&tmodel[n]); */
         time=tmodel[n];
         if(status!=0) goto error;
         if(fmod(i,50.0) ==0)
           printf("t= %6.1f hours status= %d %s\n",time/3600.0,status, sgetnewdate(info.reference,time));
          
         for (j=0;j<nodecount;j++) {
           serie[j][0][n] = h[0][nodelist[j]];
           serie[j][1][n] = p[0][nodelist[j]];
/*
           serie[j][0][n] = (beta[j][0]*h[0][nodes[j][0]-1] + beta[j][1]*h[0][nodes[j][1]-1]
                                      + beta[j][2]*h[0][nodes[j][2]-1])/100.0;
           serie[j][1][n] = -(beta[j][0]*p[0][nodes[j][0]-1] + beta[j][1]*p[0][nodes[j][1]-1]
                                       + beta[j][2]*p[0][nodes[j][2]-1])/100.0;
*/
            }      	
          }/* end for 0..nframe */
 
        fclose(in1);
        fclose(in2);
/*------------------------------------------------------------------------------
        Warning, bug corrected: memory de-allocation was wrong !!!
                                (wrong place, wrong method)*/
        archive_freeinfo(&info);
        free_smatrix(h,0,2,0,nnde-1);
        free_smatrix(p,0,2,0,nnde-1);
     
        }/* end for month...*/
     }/* end for year ... */
  
/*----------------------------------------------------------------------------*/
/* time interpolation */
/*----------------------------------------------------------------------------*/

/*------------------------------------------------------------------------------
  change CNES time in second into CNES time in days*/
  t1=t1/(24.*3600.);
  t2=t2/(24.*3600.);
  
  for (j=0;j<nst;j++) {
    if(search_status[j] != 0) continue;
    count=0;
    for(m=0;m<mes[j];m++) /* count valid data */
      {
      if ( (data_time[j][m] >= t1) && (data_time[j][m]<= t2) ) {
        t=data_time[j][m]*(24.*3600.)-t0;
/*------------------------------------------------------------------------------
        check validity with first interpolation node only, might be dangerous*/
        k=noderefc[nodes[j][0]];
        status=map_interpolate1D(serie[k][0],tmodel, mask,n+1, t, &z);
        if(z != mask) count++;
        }
      }

    if(mss==0) {
      fprintf(out,"Pt  : %-d\n",j+1);
      fprintf(out,"lon : %-f\n",sample[j].t);
      fprintf(out,"lat : %-f\n",sample[j].p);
      fprintf(out,"Mes : %-d\n",count);
      }

    for(m=0;m<mes[j];m++) /* interpolation */
      {
      if ( (data_time[j][m] >= t1) && (data_time[j][m]<= t2) ) {
        t=data_time[j][m]*(24*3600)-t0;
        z=0;
        ib=0;
        for(i=0;i<3;i++) {
          m=noderefc[nodes[j][i]];
          status=map_interpolate1D(serie[m][0],tmodel, mask,n+1, t, &dum);
          if (dum != mask)  z+=  beta[j][i]*dum/100.;
          else
            {
            z=mask;
            ib=mask;
            break;
            }
          status=map_interpolate1D(serie[m][1],tmodel, mask,n+1, t, &dum);
          if (dum != mask) ib+=-beta[j][i]*dum/100.;
          else
            {
            z=mask;
            ib=mask;
            break;
            }
          }
        if ((z != mask) && (mss==0))
                  fprintf(out,w1_format,
                        data_time[j][m],
                        data_sealevel[j][m],
                        z-ms[j]-mp[j],
                        ib);

         if ((z != mask) && (mss==1))
                fprintf(out,w2_format,
                        cycle[j],
                        sample[j].t,
                        sample[j].p,
                        data_time[j][m],
                        data_sealevel[j][m],
                        z-ms[j]-mp[j],
                        ib,
                        corr[j][m][23],
                        corr[j][m][24]);
          } /* end if t1<t<t2 */
        } /* end mes loop */
          
      if (mss==0) fprintf(out,"#-----------------------------------------------------\n");

    }/* end stations loop */
  
  fclose(out);
  
  free(ms);
  free(mp);

  fprintf(stderr,"%s -computation complete ^^^^^^^^^\n",argv[0]);
  for (j=0;j<nst;j++) {
    free(data_sealevel[j]);
    free(data_time[j]);
    }
  free(data_sealevel);
  free(data_time);
  free(sample);
  if(mss==1) {
      free(corr);
      for (i=0;i<nst;i++) {
          free(corr[i]);
          for (j=0;j<mes[i];j++) free(corr[i][j]);
        }
    }
  free(mes);
  __ERR_BASE_LINE__("exiting\n");exit(0);
  
 error:
  fprintf(stderr,"%s -computation aborted ^^^^^^^^^\n",argv[0]);
  for (j=0;j<nst;j++) {
      free(data_sealevel[j]);
      free(data_sealevel[j]);
      free(mes);
    }
  free(data_sealevel);
  free(data_time);
  free(sample);
  if(mss==1) {
      free(corr);
      for (i=0;i<nst;i++) {
          free(corr[i]);
          for (j=0;j<mes[i];j++) free(corr[i][j]);
        }
    }
  free(mes);
  __ERR_BASE_LINE__("exiting\n");exit(-1);
  
  
}/* end */


/*----------------------------------------------------------------------------*/

void write_list_header(FILE *in,FILE *out,int *nst)
     
/*----------------------------------------------------------------------------*/
{
  int i;
  char line[300];

  for(i=0;i<4;i++) {
      fgets(line,200,in);
      fprintf(out,"%s",line);
    }
  fprintf(out,"# Column 3 : MOG2D-G model sea level (in meters)\n");
  fprintf(out,"# Column 4 : MOG2D-G model S1 atmospheric tide (in meters)\n");
  fprintf(out,"# Column 5 : MOG2D-G model S2 atmospheric tide (in meters)\n");
  fprintf(out,"# Column 6 : MOG2D-medsea model sea level (in meters)\n");
  fprintf(out,"# Column 7 : MOG2D-medsea model S1 atmospheric tide (in meters)\n");
  fprintf(out,"# Column 8 : MOG2D-medsea model S2 atmospheric tide (in meters)\n");
  fprintf(out,"# Column 9 : inverted barometer (in meters)\n");
  fprintf(out,"# Column 10: inverted barometer S1 atmospheric tide (in meters)\n");
  fprintf(out,"# Column 11: inverted barometer S2 atmospheric tide (in meters)\n");
  fprintf(out,"# Column 12: loading effects (in meters)\n");
  fprintf(out,"# Column 13: solid Earth tide (in meters)\n");
  fprintf(out,"# Column 14: MOG2D-medsea model ocean tide (in meters)\n");
  fprintf(out,"# Column 15: MOG2D-medsea model S1 ocean tide (in meters)\n");
  fprintf(out,"# Column 16: MOG2D-medsea model S2 ocean tide (in meters)\n");
  fprintf(out,"# Column 17: FES99 model ocean tide (in meters)\n");
  fprintf(out,"# Column 18: FES99 model S1 ocean tide (in meters)\n");
  fprintf(out,"# Column 19: FES99 model S2 ocean tide (in meters)\n");
  fprintf(out,"# Column 20: FES2002 model ocean tide (in meters)\n");
  fprintf(out,"# Column 21: FES2002 model S1 ocean tide (in meters)\n");
  fprintf(out,"# Column 22: FES2002 model S2 ocean tide (in meters)\n");
  fprintf(out,"# Column 23: harmonic analysis ocean tide (in meters)\n");
  fprintf(out,"# Column 24: harmonic analysis S1 ocean tide (in meters)\n");
  fprintf(out,"# Column 25: harmonic analysis S2 ocean tide (in meters)\n");
  
  for(i=0;i<2;i++) {
      fgets(line,200,in);
      fprintf(out,"%s",line);
    }
 
  fgets(line,200,in);
  sscanf(line,"%d",nst);
  fprintf(out,"%s",line);
  fgets(line,200,in);
  fprintf(out,"%s",line);
  fflush(out);
}



/*--------------------------------------------------------------------------------*/

void load_xyzt(double *lon, double *lat, int *mes, double **sealevel, double **time, FILE *fic_data)

/*--------------------------------------------------------------------------------*/
{
  char a='q',line[200];
  int i, Mes;

  fgets(line,200,fic_data);     /* Pt  : */

  fgets(line,200,fic_data);     /* lon : */
  sscanf(line, "lon : %lf",lon);
  if (*lon<0) *lon+=360;
 
  fgets(line,200,fic_data);     /* lat : */
  sscanf(line, "lat : %lf",lat);
 
  fgets(line,200,fic_data);     /* Mes : */
  sscanf(line, "Mes : %d",mes);
 

  (*sealevel)=(double *)calloc((*mes),sizeof(double));
  (*time)=(double *)calloc((*mes),sizeof(double));

  for (i=0; i<(*mes);i++) {
      fgets(line,200,fic_data);
      sscanf(line,"%lf %lf ",&((*time)[i]),&((*sealevel)[i]));
    }

  for (i=0; i<(*mes);i++)
    if( (*sealevel)[i] == 9.999)
      (*sealevel)[i]=99.9999;

  fgets(line,200,fic_data);    /*  #---------- */

}/*end*/




