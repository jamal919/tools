

/*******************************************************************************

  T-UGO tools, 2006-2009

  Unstructured Ocean Grid initiative

Contributors:

  Florent Lyard      LEGOS/CNRS, Toulouse, France
  Yoann Le Bars      LEGOS, Toulouse, France (PhD)
  Laurent Roblou     LEGOS, Toulouse, France
  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
  Frédéric Dupont    Université de Laval à Québec, Canada

E-mail: florent.lyard@legos.obs-mip.fr

*******************************************************************************/

#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>

#include "tools-structures.h"

#include "mgr.h"
#include "geo.h"
#include "list.h"
#include "polygones.h"
#include "spectrum.h"

#include "legend.h"
#include "legend.def"

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 int mgr_select(vector<mgr_t> mgr, int *list, unsigned short *mgr_region, unsigned short target)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int n;
  int count=0;
  unsigned short flag;
  
//  target=NORTH_ATLANTIC;
  
  for(n=0;n<mgr.size();n++) {
    flag=mgr_region[n]/100;
    if(flag==target) {
      list[count]=n;
      count++;
      }
    }
  return(count);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 int *mgr_select(vector<plg_t> & polygons, vector<mgr_t> & mgr)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,n,status;
  int *inside=NULL,in, count=0;
  plg_t *polygones=NULL;
  int npolygones=0;
  double t,p;
  projPJ proj=0;
  frame_t frame;

  if(polygons.size()==0) return(0);
  
/*------------------------------------------------------------------------------
  checks polygons sanity. Optionally check if polar */
  status=plg_checkAutoSecant(polygons, PLG_SPHERICAL);
  if(status==-1) status=CheckPolar(polygons[0], proj);
 
  inside=new int[mgr.size()];
  frame=plg_cartesian_minmax(polygones, npolygones);

  for (n=0;n<mgr.size();n++) {
    t=mgr[n].loc.lon;
    p=mgr[n].loc.lat;
    in=0;
    for (k=0;k<polygons.size();k++) {
      frame=plg_recale(&(polygons[k]), 1);
      t=degree_recale(mgr[n].loc.lon,0.5*(frame.xmin+frame.xmax));
      
      if(proj==0) {
        inside[n]=plg_single( polygons[k], t, p, &in, PLG_SPHERICAL);
        }
      else {
        double x,y;
        geo_to_projection(proj, p,t, &x, &y);
        inside[n]=plg_single( polygons[k], x, y, &in, PLG_CARTESIAN);
        }
      
      if(inside[n]==PLG_POINT_BOUNDARY) {
        break;
        }
      
      //STDERR_BASE_LINE_FUNC("[%d:%s,%d]: (%g;%g) %d\n",n,mgr[n].name,k,t,p,inside[n]);
      }
    
    if(inside[n]==1) {
      count++;
      }
    
    }
  
  printf("%d station found included over %d\n",count,mgr.size());
  
  return(inside);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 int *mgr_select(const char *poly ,vector<mgr_t> mgr)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int *inside=NULL;
  vector<plg_t> polygons;

  if(poly!=NULL)
    plg_load( poly, PLG_FORMAT_SCAN, polygons);
  else
    return(0);

  if(polygons.size()==0) return(0);

  printf("filter mgr set from %s polygons...\n",poly);
  inside=mgr_select(polygons, mgr);
  
  plg_destroy_entries(polygons);
  
  return(inside);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 int mgr_loadobs(const char *filename, char *wave, int *nmgr ,vector<mgr_t> & mgr, double ***covariance, int ***list, int **card)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int  k,nblines,nitems;
  FILE *in=NULL;
  int i;
  int *n1=NULL,*n2=NULL;
  float a,p,std;
  char mgr_class='X';
  double *cov=NULL;
  char units, tmp[10];
    
  in=fopen(filename ,"r");
  if (in == NULL) TRAP_ERR_EXIT(-1,"cannot open file %s (%d %s)\n", filename, errno, strerror(errno));

  nitems=fscanf(in,"%d %c %s\n",nmgr,&units,tmp);
  
  if(nitems==1) {
    units='m';
    }
  if(nitems==3) {
    strcpy(wave,tmp);
    }

  for(i=0;i<*nmgr;i++) {
    mgr_t gauge;
    exitIfNull(
      gauge.data=(mgr_data_t *)calloc(1,sizeof(mgr_data_t *))
      );
    gauge.nwave=1;
    fscanf(in, "%d %lf %lf %f %f %f %c %s\n",&k,&(gauge.loc.lon),&(gauge.loc.lat),&a,&p,&std,&mgr_class,gauge.name);
    strcpy(gauge.data[0].constituent.name,wave);
    gauge.data[0].amp=a;
    gauge.data[0].phi=p;
    gauge.loc.units=strdup("degrees");
    mgr.push_back(gauge);
    }

  fscanf(in,"%d\n",&nblines);

  if(nblines==0) goto finished;
  
  if(covariance==0)  goto finished;

  n1=new int[nblines];
  n2=new int[nblines];
  cov=new double[nblines];
  for(k=0;k<nblines;k++) {
    fscanf(in, "%d %d %lf\n",&n1[k],&n2[k],&cov[k]);
    }

  *list      =new int*[*nmgr];
  *card      =new int[*nmgr];
  *covariance=new double*[*nmgr];

  for(i=0;i<*nmgr;i++) {
    (*card)[i]=0;
    }

  for(k=0;k<nblines;k++) {
    (*card)[n1[k]-1]++;
    }

  for(i=0;i<*nmgr;i++) {
    (*list)[i]=new int[(*card)[i]];
    (*covariance)[i]=new double[(*card)[i]];
    }

  for(i=0;i<*nmgr;i++) {
    (*card)[i]=0;
    }
  for(k=0;k<nblines;k++) {
    (*list)[n1[k]-1][(*card)[n1[k]-1]]=n2[k]-1;
    (*covariance)[n1[k]-1][(*card)[n1[k]-1]]=cov[k];
    (*card)[n1[k]-1]++;
    }

  for(i=0;i<*nmgr;i++) {
    mgr[i].data[0].error=sqrt((*covariance)[i][0]);
    }
finished:
  fclose(in) ;

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 int mgr_save_obs(const char *filename, vector<mgr_t> & mgr, const char *wave)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int  n;
  FILE *out=NULL;
  int i,j,nobs;
  float a,p,std;
  char mgr_class='X';
  size_t start,end;
  
  out =fopen(filename ,"w");
  if(out==0) TRAP_ERR_EXIT(-1,"file opening issue : %s \n", filename);
    
  nobs=0;
  for(i=0;i<mgr.size();i++) {
    j=mgr[i].wave_index(wave);
    if(j==-1) continue;
    nobs++;
    }

/**----------------------------------------------------------------------------
  temporary */
//  fprintf(out,"%d M\n",nobs);

  n=0;
  for(i=0;i<mgr.size();i++) {
    string s=mgr[i].name,ss;
    j=mgr[i].wave_index(wave);
    if(j==-1) continue;
    a=mgr[i].data[j].amp;
    p=mgr[i].data[j].phi;
//     variance=covariance[i][0];
//     std=sqrt(variance);
    std=mgr[i].data[j].error;
    start=s.rfind("2_t");
    end=s.rfind("_sla");
    if(start!=string::npos) {
      ss=s.substr(start+3, end-start-3);
      }
    else {
      ss=s;
      }
    
/* *----------------------------------------------------------------------------
    fortran: FORMAT(i6,5F12.6,1X,a1,8X,a40,i6,2X,a10)*/
//     fprintf(out, "%6d %11.6f %11.6f %11.6f %11.6f %11.6f %c        %s\n",
//             n,mgr[i].loc.lon,mgr[i].loc.lat,a,p,std,mgr_class, mgr[i].name);
/**----------------------------------------------------------------------------
  temporary */
    fprintf(out, "%6d %11.6f %11.6f %11.6f %11.6f %11.6f %c        %s\n",
            n,mgr[i].loc.lon,mgr[i].loc.lat,a,p,std,mgr_class, ss.c_str());
    n++;
    }

  fclose(out);

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 int mgr_save_obs4assim(const char *filename, char *poly, char *wave, int nmgr ,vector<mgr_t> mgr, float **covariance, int **list, int *card)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int  k,m,n,nblines;
  FILE *out=NULL;
  int i,j,nobs;
  float a,p,variance,std;
  int *inside=NULL,*obsindex=NULL;
  char mgr_class='X';

/* *----------------------------------------------------------------------------
  filter mgr set before saving */
  nobs=0;
  if(poly!=NULL) {
    inside=mgr_select(poly,mgr);
    if(inside==0) return(-1);
    for (n=0;n<nmgr;n++) {
      if(inside[n]==PLG_POINT_INTERIOR) nobs++;
      }
    }
  else {
    inside=new int[nmgr];
    nobs=nmgr;
    for (n=0;n<nmgr;n++) {
      inside[n]=PLG_POINT_INTERIOR;
      }
    }

  obsindex=new int[nmgr];

  nobs=0;
  for(n=0;n<nmgr;n++) {
    obsindex[n]=-1;
    if(inside[n]!=PLG_POINT_INTERIOR) continue;
    j=mgr[n].wave_index(wave);
    if(j==-1) continue;
    obsindex[n]=nobs;
    nobs++;
    }

  out =fopen(filename ,"w");
  if(out==0) TRAP_ERR_EXIT(-1,"file opening issue : %s \n", filename);

  fprintf(out,"%d M\n",nobs);

  for(i=0;i<nmgr;i++) {
    if(inside[i]!=PLG_POINT_INTERIOR) continue;
    j=mgr[i].wave_index(wave);
    if(j==-1) continue;
    a=mgr[i].data[j].amp;
    p=mgr[i].data[j].phi;
    variance=covariance[i][0];
    std=sqrt(variance);
/* *----------------------------------------------------------------------------
    fortran: FORMAT(i6,5F12.6,1X,a1,8X,a40,i6,2X,a10)*/
    fprintf(out, "%6d %11.6f %11.6f %11.6f %11.6f %11.6f %c        %s\n",
            obsindex[i]+1,mgr[i].loc.lon,mgr[i].loc.lat,a,p,std,mgr_class, mgr[i].name);
    }

  if(card==0) {
    nblines=0;
    fprintf(out,"%d\n",nblines);
    goto finished;
    }

  nblines=0;
  for(n=0;n<nmgr;n++) {
    if(inside[n]!=PLG_POINT_INTERIOR) continue;
    for(k=0;k<card[n];k++) {
      m=list[n][k];
      if(inside[m]!=PLG_POINT_INTERIOR) continue;
      nblines++;
      }
    }

  fprintf(out,"%d\n",nblines);
  for(n=0;n<nmgr;n++) {
    if(inside[n]!=PLG_POINT_INTERIOR) continue;
    for(k=0;k<card[n];k++) {
      m=list[n][k];
      if(inside[m]!=PLG_POINT_INTERIOR) continue;
      fprintf(out, "%d %d %e\n",obsindex[n]+1,obsindex[m]+1,covariance[n][k]);
      }
    }

finished:
  fclose(out);
  delete[] inside;
  delete[] obsindex;

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 int mgr_save_obs4assim2(const char *filename, char *wave, int nmgr ,vector<mgr_t> mgr, float **covariance, int **list, int *card)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int  k,nblines;
  FILE *out=NULL;
  int i,j,nbs;
  float a,p,variance,std;
  char mgr_class='A';

  out =fopen(filename ,"w");
  if(out==0) TRAP_ERR_EXIT(-1,"file opening issue : %s \n", filename);

  fprintf(out,"%d M\n",nmgr);

  nbs=0;
  for(i=0;i<nmgr;i++) {
    nbs=nbs+1;
    j=mgr[i].wave_index(wave);
    a=mgr[i].data[j].amp;
    p=mgr[i].data[j].phi;
    variance=covariance[i][0];
    std=sqrt(variance);
    if(std/mgr[i].data[j].amp>0.1) {
      fprintf(out, "%6d %11.6f %11.6f %11.6f %11.6f %11.6f %c        %s\n", i+1,mgr[i].loc.lon,mgr[i].loc.lat,a,p,std,"W", mgr[i].name);
      }
    else {
      fprintf(out, "%6d %11.6f %11.6f %11.6f %11.6f %11.6f %c        %s\n", i+1,mgr[i].loc.lon,mgr[i].loc.lat,a,p,std,mgr_class, mgr[i].name);
      }
    }

  if(card==0) goto finished;

  nblines=0;
  for(i=0;i<nmgr;i++) {
    nblines+=card[i];
    }
/* *----------------------------------------------------------------------------
  neutralized*/
  nblines=0;
  fprintf(out,"%d\n",nblines);
  for(i=0;i<nmgr;i++) {
    for(k=0;k<card[i];k++) {
      fprintf(out, "%d %d %e\n",i+1,list[i][k]+1,covariance[i][k]);
      }
    }

finished:
  fclose(out) ;

  return(0);
}



/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 int mgr_Constants2Legend_00(const char *filename, char *wave, int nmgr ,vector<mgr_t> mgr, float **errors, float **covariance, int *track, int **list, int *card)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int  l,n,nlegends,status;
  int i,j,m,w;
  float variance,std;
  legend_t   *legends=NULL;
  legend02_t *legends02=NULL;
  char *comments=NULL;

  legends=new legend_t;

  l=0;
  legends[l].ID=l;
  legends[l].Type=LGD_GRAPH;
  legends[l].T=0;
  legends[l].P=0;
  legends[l].X=0;
  legends[l].Y=0;

/* *----------------------------------------------------------------------------
  create legend with all valid primary xovers */
  legends02=new legend02_t;
  legends02->np=nmgr;
  legends02->nz=9;
  legends02->pen0=0;
  legends02->pen1=0;
  legends02->points=new point_t[nmgr];

  legends[l].ptr=(char *) legends02;

  nlegends=1;

  for (m=0;m<nmgr;m++) {
    legends02->points[m].z=new float[legends02->nz];
    }

  legends02->ID=0;
  legends02->np=0;
  for(i=0;i<nmgr;i++) {
    w=mgr[i].wave_index(wave);
    if(w==-1) continue;
    legends02->np++;
    n=legends02->np-1;
    legends02->points[n].t=mgr[i].loc.lon;
    legends02->points[n].p=mgr[i].loc.lat;
    variance=covariance[i][0];
    std=sqrt(variance);
    legends02->points[n].N=mgr[i].track;
    j=0;
    legends02->points[n].z[j++]=mgr[i].data[w].amp;
    legends02->points[n].z[j++]=mgr[i].data[w].phi;
    legends02->points[n].z[j++]=errors[0][i];
    legends02->points[n].z[j++]=errors[1][i];
    legends02->points[n].z[j++]=errors[2][i];
    legends02->points[n].z[j++]=std;
    legends02->points[n].z[j++]=100.*std/mgr[i].data[w].amp;
    legends02->points[n].z[j++]=mgr[i].duree;
    legends02->points[n].z[j++]=track[i];
    }

  comments=new char[1024];
  comments[0]=0;

  strcat(comments,"#  0: xover index \n");
  strcat(comments,"#  1: longitude\n");
  strcat(comments,"#  2: latitude\n");
  strcat(comments,"#  3-> 0: amplitude\n");
  strcat(comments,"#  4-> 1: phase lag\n");
  strcat(comments,"#  5-> 2: analysis error\n");
  strcat(comments,"#  6-> 3: alias error\n");
  strcat(comments,"#  7-> 4: coastal error\n");
  strcat(comments,"#  8-> 5: combined error\n");
  strcat(comments,"#  9-> 6: relative combined error (% of amplitude)\n");
  strcat(comments,"# 10-> 7: #cycles\n");
  strcat(comments,"# 11-> 8: track ID\n");

  status=lgd_save(filename, legends, nlegends,NULL,comments);

  delete[] comments;

  for (m=0;m<nmgr;m++) {
    delete[] legends02->points[m].z;
    }
  delete[] legends02->points;

  delete legends02;
  delete legends;

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 int mgr_Constants2legend_01(const char *filename, char *wave, vector<mgr_t> mgr)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int  l,n,nlegends,status;
  int i,j,m,w;
  float /*variance,*/std;
  legend_t   *legends=NULL;
  legend02_t *legends02=NULL;
  char *comments=NULL;
  int nmgr=mgr.size();


  legends=new legend_t;

  l=0;
  legends[l].ID=l;
  legends[l].Type=LGD_GRAPH;
  legends[l].T=0;
  legends[l].P=0;
  legends[l].X=0;
  legends[l].Y=0;

/* *----------------------------------------------------------------------------
  create legend with all valid primary xovers */
  legends02=new legend02_t;
  legends02->np=nmgr;
  legends02->nz=5;
  legends02->pen0=0;
  legends02->pen1=0;
  legends02->points=new point_t[nmgr];

  legends[l].ptr=(char *) legends02;

  nlegends=1;

  for (m=0;m<nmgr;m++) {
    legends02->points[m].z=new float[legends02->nz];
    }

  legends02->ID=0;
  legends02->np=0;
  for(i=0;i<nmgr;i++) {
    legends02->np++;
    n=legends02->np-1;
    legends02->points[n].t=mgr[i].loc.lon;
    legends02->points[n].p=mgr[i].loc.lat;
    w=mgr[i].wave_index(wave);
//     variance=covariance[i][0];
//     std=sqrt(variance);
    legends02->points[n].N=mgr[i].track;
    j=0;
    legends02->points[n].z[j++]=mgr[i].data[w].amp;
    legends02->points[n].z[j++]=mgr[i].data[w].phi;
    legends02->points[n].z[j++]=std;
    legends02->points[n].z[j++]=100.*std/mgr[i].data[w].amp;
//     legends02->points[n].z[j++]=track[i];
    }

  comments=new char[1024];
  comments[0]=0;

  strcat(comments,"#  0: xover index \n");
  strcat(comments,"#  1: longitude\n");
  strcat(comments,"#  2: latitude\n");
  strcat(comments,"#  3: amplitude\n");
  strcat(comments,"#  4: phase lag\n");
  strcat(comments,"#  5: analysis error\n");
  strcat(comments,"#  6: relative analysis error (% of amplitude)\n");
  strcat(comments,"#  7: track ID\n");

  status=lgd_save(filename, legends, nlegends,NULL,comments);

  delete[] comments;

  for (m=0;m<nmgr;m++) {
    delete[] legends02->points[m].z;
    }
  delete[] legends02->points;

  delete legends02;
  delete legends;
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 int mgr_Timeserie2Legend(const char *filename, double lon, double lat, double *time, int nframes, double **values, int nvalues)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int  k,l,n,nlegends,status;
  int i,j;
  legend_t   *legends=NULL;
  legend03_t *legends03=NULL;
  char *comments=NULL;


  legends=new legend_t;

  l=0;
  legends[l].ID=l;
  legends[l].Type=LGD_PROFILE;
  legends[l].T=0;
  legends[l].P=0;
  legends[l].X=0;
  legends[l].Y=0;

/* *----------------------------------------------------------------------------
  create legend with all valid primary xovers */
  legends03=new legend03_t;
  legends03->np=1;
  legends03->ndepths =nframes;
  legends03->ncolumns=nvalues;
  legends03->pen0=0;
  legends03->pen1=0;
  legends03->points=new point3D_t[legends03->np];

  legends[l].ptr=(char *) legends03;

  nlegends=1;

  legends03->points[0].z=new float *[legends03->ndepths];
  legends03->points[0].depths=new float [legends03->ndepths];

  legends03->ID=0;
//  legends03->np=0;
  for(i=0;i<legends03->np;i++) {
    n=i;
    legends03->points[n].t=lon;
    legends03->points[n].p=lat;
    legends03->points[n].N=n;
    for(k=0;k<legends03->ndepths;k++) {
      legends03->points[0].depths[k]=time[k];
      legends03->points[n].z[k]=new float [legends03->ncolumns];
      for(j=0;j<legends03->ncolumns;j++) {
        legends03->points[n].z[k][j]=values[j][k];
        }
      }
    }

  comments=new char[1024];
  comments[0]=0;

  strcat(comments,"#  0: xover index \n");
  strcat(comments,"#  1: longitude\n");
  strcat(comments,"#  2: latitude\n");
  strcat(comments,"#  3: amplitude\n");
  strcat(comments,"#  4: phase lag\n");
  strcat(comments,"#  5: analysis error\n");
  strcat(comments,"#  6: relative analysis error (% of amplitude)\n");
  strcat(comments,"#  7: track ID\n");

  status=lgd_save(filename, legends, nlegends,NULL,comments);

  delete[] comments;

//   for (m=0;m<nmgr;m++) {
//     delete[] legends03->points[m].z;
//     }
//   delete[] legends03->points;

  delete legends03;
  delete legends;
  return(status);
}



/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int mgr_loadLEGOS(const char *filename, vector<mgr_t> & mgr)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  if(!strrncmp(filename,".nc")){
    status=mgr_load_netcdf(filename, mgr);
    }
  else if(!strrncmp(filename,".mgr")){
    status=mgr_load_ascii(filename, mgr);
    }
  else
    status=-1;
  
  return status;
}


// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
// int mgr_load_ascii_obsolete(const char *filename, vector<mgr_t> & mgr_serie)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   FILE *in=NULL;
//   char line1[256],line2[256],line3[256],line4[256],line5[256],line6[256],line7[256],line8[256];
//   char dum1[256],dum2[256];
//   int i,test,count=0,nitems;
//   int track=0,previous=0;
//   char loca[256],carac,nord,est;
//   char *dum=NULL;
//   char *status=NULL;
// 
//   in =fopen(filename ,"r");
//   if(in==0) TRAP_ERR_EXIT(-1,"cannot open file %s\n", filename);
// 
//   List *station_list=NULL;
//   List *waves_list=NULL;
//   mgr_t *mgr_station=NULL;
//   mgr_data_t *data=NULL;
//   station_list=NULL;
// 
// /**----------------------------------------------------------------------------
//   skip first line */
//   status=fgets(line1,sizeof(line1),in);
// 
//   do {
//     exitIfNull(
//       mgr_station=new mgr_t
//       );
//     waves_list=NULL;
//     if(status==NULL)break;
// /**----------------------------------------------------------------------------
//     empty line */
//     status=fgets(line2,sizeof(line2),in);
// /**----------------------------------------------------------------------------
//     station name and number */
//     fgets(line3,sizeof(line3),in);
//     sscanf(line3,"%s %s %d",dum1,dum2,&(mgr_station.number));
// //    printf("%d %s \n",count, line3);track-ref.TP+J1.118.dat
// /* *----------------------------------------------------------------------------
//     fragile patch to supply track index if not present*/
//     if(mgr_station.number < previous) track++;
//     previous=mgr_station.number;
//     mgr_station.track=track;
//     dum=strchr(line3,'\n');
//     dum[0]='\0';
//     dum--;
//     while(dum[0]==' ') {
//       dum[0]='\0';
//       dum--;
//       }
//     dum=strchr(line3,':')+2;
//     while(dum[0]==' ') dum++;
//     strcpy(mgr_station.name,dum);
//     dum=strstr(mgr_station.name,"track-ref.TP+J1.");
//     if(dum!=0) {
//       dum+=strlen("track-ref.TP+J1.");
//       sscanf(dum,"%3d",&(mgr_station.track));
//       }
//     dum=strchr(mgr_station.name,'\n');
// /**----------------------------------------------------------------------------
//     data origin */
//     fgets(line4,sizeof(line4),in);
//     sscanf(line4,"%17s : %50s ",loca,mgr_station.origine);
// //     dum=strchr(line4,':')+2;
// //     while(dum[0]==' ') dum++;
// //     strcpy(mgr_station.origine,dum);
// /**----------------------------------------------------------------------------
//     station location */
//     fgets(line5,sizeof(line5),in);
//     sscanf(line5,"%13s %c %lf%c %lf%c ",loca,&carac,&mgr_station.loc.lat,&nord,&mgr_station.loc.lon,&est);
//     mgr_station.loc.units=strdup("degrees");
//     fgets(line6,sizeof(line6),in);
//     nitems=sscanf(line6,"%17s : Debut: %s     Fin: %s         Duree:  %d",loca,&mgr_station.debut,&mgr_station.fin,&mgr_station.duree);
//     fgets(line7,sizeof(line7),in);
//     test=4;
//     count++;
// /**----------------------------------------------------------------------------
//     harmonic constants */
//     while( ( (test==4)||(test==5)||(test==7) ) && (!feof(in)) ) {
//       exitIfNull(
//         data=(mgr_data_t *) calloc(1,sizeof(mgr_data_t))
//         );
//       fgets(line8,sizeof(line8),in);
//       if(!feof(in)) {
//         test=sscanf(line8,"%d %s %lf %lf %lf",&data->constituent.code,data->constituent.name,&(data->amp),&(data->phi),&(data->error));
// 
// //        test=sscanf(line8,"%d %s %lf %lf %lf %lf %lf",&data->code,data->wave,&(data->amp),&(data->phi),&f1,&f2,&(data->error));
// //        test=sscanf(line8,"%d %s %lf %lf %lf %lf %lf",&data->code,data->wave,&f1,&f2,&(data->amp),&(data->phi),&(data->error));
// 
//         switch(test) {
//           case 4:
// /* *----------------------------------------------------------------------------
//             add wave data to wave list*/
//             waves_list=list_append(waves_list,data);
//             data->error=-1;
//             break;
//           case 5:
// /* *----------------------------------------------------------------------------
//             add wave data to wave list*/
//             waves_list=list_append(waves_list,data);
//             break;
//           case 7:
// /* *----------------------------------------------------------------------------
//             add wave data to wave list*/
//             waves_list=list_append(waves_list,data);
//             break;
//           default:
//             break;
//           }
//         }
//       }
// 
// /* *----------------------------------------------------------------------------
//     create wave list*/
//     exitIfNull(
//       mgr_station.data=(mgr_data_t *)calloc(list_length(waves_list),sizeof(mgr_data_t))
//       );
//     for(i=0;i<list_length(waves_list);i++) {
//       mgr_station.data[i]=*((mgr_data_t *) list_nth_data(waves_list,i));
//       }
//     mgr_station.nwave=list_length(waves_list);
// /* *----------------------------------------------------------------------------
//     free temporary wave list*/
//     waves_list=list_free(waves_list);
// /* *----------------------------------------------------------------------------
//     add station data to station list*/
//     station_list=list_append(station_list,mgr_station);
//     } while (!feof(in));
// 
// /* *----------------------------------------------------------------------------
//   create station list*/
//   exitIfNull(
//     *mgr_serie=(mgr_t **) calloc(list_length(station_list),sizeof(mgr_t *))
//     );
//   for(i=0;i<list_length(station_list);i++) {
//     (*mgr_serie)[i]=(mgr_t *) list_nth_data(station_list,i);
//     }
//   int reout=list_length(station_list);
//   station_list=list_free(station_list);
// 
//   fclose(in) ;
//   return(reout);
// }

/*###############################################################################################*/

  int mgr_load_ascii(const char *filename, vector<mgr_t> & mgr_serie)

/*###############################################################################################*/
{
  FILE *in=NULL;
  char line1[256],line2[256],line3[256],line4[256],line5[256],line6[256],line7[256],line8[256];
  char second_name[256];
  int i,count=0,nitems,lineI;
  int track=0,previous=0;
  char loca[256],carac,nord,est;
  char *dum=NULL;
  char *status=NULL;
  
  in =fopen(filename ,"r");
  if (in == NULL) TRAP_ERR_EXIT(-1,"cannot open file %s (%d %s)\n", filename, errno, strerror(errno));
  lineI=0;

  vector <mgr_t *> station_list;
/**----------------------------------------------------------------------------
  skip first line */
  status=fgets(line1,sizeof(line1),in);lineI++;

  do {
    mgr_t mgr_station;
    vector<mgr_data_t> waves_list;
    if(status==NULL)break;
/**----------------------------------------------------------------------------
    empty line */
    status=fgets(line2,sizeof(line2),in);lineI++;
/**----------------------------------------------------------------------------
    station name and number */
    status=fgets(line3,sizeof(line3),in);lineI++;
//     STDERR_BASE_LINE("%d=fgets(\"%s\",%d,)\n",status!=0,line3,sizeof(line3));
    sscanf(&line3[11],"%d",&mgr_station.number);
//    printf("%d %s \n",count, line3);track-ref.TP+J1.118.dat
/* *----------------------------------------------------------------------------
    fragile patch to supply track index if not present*/
    if(mgr_station.number < previous) track++;
    previous=mgr_station.number;
    mgr_station.track=track;
    dum=strchr(line3,'\n');
    /* if the file does not start with "___________________________________\n\n"
    or the equivalent then it should SegFault here */
    do {
      dum[0]='\0';
      dum--;
      } while(dum[0]==' ');
    dum=strchr(line3,':');
    if(dum==0)
      TRAP_ERR_RETURN(0,1,"%s:%d: no ':' found\n",filename,lineI);
    dum++;
    while(dum[0]==' ') dum++;
    strcpy(mgr_station.name,dum);
    dum=strstr(mgr_station.name,"track-ref.TP+J1.");
    if(dum!=0) {
      dum+=strlen("track-ref.TP+J1.");
      sscanf(dum,"%3d",&(mgr_station.track));
      }
    dum=strchr(mgr_station.name,'\n');
/**----------------------------------------------------------------------------
    data origin */
    status=fgets(line4,sizeof(line4),in);lineI++;
    sscanf(line4,"%17s : %50s ",loca,mgr_station.origine);
//     dum=strchr(line4,':')+2;
//     while(dum[0]==' ') dum++;
//     strcpy(mgr_station.origine,dum);
/**----------------------------------------------------------------------------
    station location */
    status=fgets(line5,sizeof(line5),in);lineI++;
    sscanf(line5,"%13s %c %lf%c %lf%c %lf%c",loca,&carac,&mgr_station.loc.lat,&nord,&mgr_station.loc.lon,&est,&mgr_station.loc.depth,&est);
    
//     double mlon=mgr_station.loc.lon-floor(mgr_station.loc.lon)-1;
//     mgr_station.loc.lon=1+floor(mgr_station.loc.lon)+mlon*100./60.;
//
//     double mlat=mgr_station.loc.lat-floor(mgr_station.loc.lat);
//     mgr_station.loc.lat=floor(mgr_station.loc.lat)+mlat*100./60.;
    
    mgr_station.loc.units=strdup("degrees");
    status=fgets(line6,sizeof(line6),in);lineI++;
    nitems=sscanf(line6,"%17s : Debut: %s     Fin: %s         Duree:  %d",loca,&mgr_station.debut,&mgr_station.fin,&mgr_station.duree);
    status=fgets(line7,sizeof(line7),in);lineI++;
    line7[strlen(line7)-1]=0;
    strcpy(mgr_station.validation, &(line7[20]));
//    nitems=sscanf(line7,"%s : %s",&dum1, &mgr_station.validation);
    count++;
/**----------------------------------------------------------------------------
    harmonic constants */
    do {
      mgr_data_t data;
      data.error=NAN;
      
      status=fgets(line8,sizeof(line8),in);lineI++;
      if(status==0) break;
      
      nitems=sscanf(line8,"%d %s %lf %lf %lf",&data.constituent.code,data.constituent.name,&data.amp,&data.phi,&data.error);
      
//      nitems=sscanf(line8,"%d %s %lf %lf %lf %lf %lf",&data.code,data.wave,&data.amp,&data.phi,&f1,&f2,&data.error);
//      nitems=sscanf(line8,"%d %s %lf %lf %lf %lf %lf",&data.code,data.wave,&f1,&f2,&data.amp,&data.phi,&data.error);
      
      switch(nitems) {
        case 2:
        case 3:
/* *----------------------------------------------------------------------------
          shom funny double wave name*/
          nitems=sscanf(line8,"%d %s %s %lf %lf %lf",&data.constituent.code,data.constituent.name, second_name, &data.amp,&data.phi,&data.error);
          waves_list.push_back(data);
          if(nitems<5)
            TRAP_ERR_RETURN(0,1,"%s:%d: only %d fields recognised\n",filename,lineI,nitems);
          break;
        
        case 4:
/* *----------------------------------------------------------------------------
          no error bar given; add wave data to wave list*/
          waves_list.push_back(data);
          data.error=-1;
          break;
        case 5:
/* *----------------------------------------------------------------------------
          error bar prescribed; add wave data to wave list*/
          waves_list.push_back(data);
          break;
        }
      
      } while ( 4<=nitems and nitems<=6);

/* *----------------------------------------------------------------------------
    create wave list*/
    exitIfNull(
      mgr_station.data=(mgr_data_t *) calloc(waves_list.size(),sizeof(mgr_data_t))
      );
    for(i=0;i<waves_list.size();i++) {
      mgr_station.data[i]=waves_list[i];
      }
    mgr_station.nwave=waves_list.size();
    
/* *----------------------------------------------------------------------------
    add station data to station list*/
    mgr_serie.push_back(mgr_station);
    } while (!feof(in));

  int reout=mgr_serie.size();

  fclose(in) ;
  return(reout);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double sign_from(char c)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  if(c=='-')
    return -1.;
  return 1.;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int mgr_loadSHOM(const char *filename, vector<mgr_t> & mgr_serie)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  FILE *in=NULL;
  const size_t lC=1024;
  char l[lC],*status;
  
  in =fopen(filename ,"r");
  if (in == NULL) TRAP_ERR_EXIT(-1,"cannot open file %s (%d %s)\n", filename, errno, strerror(errno));
  
  do {
    mgr_t mgr_station;
    
/*------------------------------------------------------------------------------
    header line */
    status=fgets(l,lC,in);
    if(status==0) break;
    
    int count,
      lad,lam,las,lod,lom,los,
      tz,tzd,ys,doys,
      missing,method,trust,
      yc,mc,dc,  used_duration;
    char lap,lop,reference_level;
    count=sscanf(l,"  %4d %c%02d%02d%02d   %c%03d%02d%02d  %3d%1d%s  %04d%03d  %d  "
      "%d %d %c%02d%s  %04d%02d%02d  %d",
      &mgr_station.number,  &lap,&lad,&lam,&las,  &lop,&lod,&lom,&los,
      &tz,&tzd,mgr_station.name,  &ys,&doys,  &mgr_station.duree,
      
      &missing,&method, &reference_level,&trust,mgr_station.origine,
      &yc,&mc,&dc, &used_duration);
    
    /* NOTE: THE HEADER LINE MAY NOT BE THAT WELL DOCUMENTED
      There are 24 fields, if only 20 are read, it will be OK. */
    if(count<20) TRAP_ERR_EXIT(-1,"%s:format not recognised, only %d fields parsed\n",__func__,count);
    
    mgr_station.loc.lat=sign_from(lap)*(lad+(lam+las/60.)/60.);
    mgr_station.loc.lon=sign_from(lop)*(lod+(lom+los/60.)/60.);
    mgr_station.loc.depth=NAN;
    
    date_t d(ys,1,doys);
    double dd=cnes_time(d,'s');
    d=poctime_getdatecnes(dd,'s');
    sprintf(mgr_station.debut,"%04d/%02d/%02d",d.year,d.month,d.day);
    
/*------------------------------------------------------------------------------
    harmonic constants */
    vector<mgr_data_t> waves_list;
    
    while(true){
      status=fgets(l,lC,in);
      if(status==0) break;
      
      mgr_data_t data;
      
      /* WARNING: PARSING FORTRAN OUTPUT BELOW... */
      char aw[10],pw[10];
      count=sscanf(l,"%9c%6c%6c",data.constituent.name,&aw,&pw);
      
      vector<string> tokens,keys,values;
      string delimiter=" ";
    
      delimiter=" ";
      tokens=string_split((string) data.constituent.name, delimiter);
      
      if(tokens.size()==0) {
        printf("unable to get wave name from : %s\n", data.constituent.name);
        continue;
        }
      strcpy(data.constituent.name, tokens[0].c_str());
      
      if(count<3) break;
      
      data.amp=atof(aw);
      data.phi=atof(pw);
      
      waves_list.push_back(data);
      }
    
    mgr_station.data=copy(waves_list);
    mgr_station.nwave=waves_list.size();
    
/* *----------------------------------------------------------------------------
    add station data to station list*/
    mgr_serie.push_back(mgr_station);
    
    if(status==0 or count<3) break;
    } while (!feof(in));

  int reout=mgr_serie.size();

  fclose(in) ;
  return(reout);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int mgr_loadDART(const char *filename, vector<mgr_t> & mgr)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  std::ifstream input(filename);
  
  int nitems;
  int nmgr;
  std::string line, sub;
  double factor;
  double frequency;
  char units;
  size_t pos;
  
  int k = 0;
  mgr_t mgr_station;
  mgr_data_t *data;
  
//   pos=line.find ("Headers=", 0);
//   sub = line.substr (pos);
//   nitems=sscanf(sub.c_str(),"%s %d", s1,&nheaders);
  
  cout << endl << "-------- starting conversion --------" << endl << endl;

  do {
    std::getline(input, line );
    } while(line.find("Station =")==string::npos);

  pos=line.rfind ("=");
  sub = line.substr (pos+1);
  nitems=sscanf(sub.c_str(),"%s",mgr_station.name);
  
  nitems=sscanf(filename,"%s",mgr_station.name);
  nitems=sscanf("DART","%s",mgr_station.origine);
  
  nitems=sscanf("01/01/1950","%s",mgr_station.debut);
  nitems=sscanf("01/01/1950","%s",mgr_station.fin);
  
  std::getline(input, line );
  pos=line.rfind ("=");
  sub = line.substr (pos+1);
  nitems=sscanf(sub.c_str(),"%lf%c",&mgr_station.loc.lat,&units);
  std::getline(input, line );
  pos=line.rfind ("=");
  sub = line.substr (pos+1);
  nitems=sscanf(sub.c_str(),"%lf%c",&mgr_station.loc.lon,&units);
  
//   switch (c1) {
//     case 'N':
//      mgr_station.loc.lat=lat+latmin/60.0;
//       break;
//     case 'S':
//       mgr_station.loc.lat=-(lat+latmin/60.0);
//       break;
//     default:
//       TRAP_ERR_EXIT(-1, "format error\n");
//     }
//   switch (c2) {
//     case 'E':
//       mgr_station.loc.lat=lon+lonmin/60.0;
//       break;
//     case 'W':
//       mgr_station.loc.lat=-(lon+lonmin/60.0);
//       break;
//     default:
//       TRAP_ERR_EXIT(-1, "format error\n");
//     }
//  nitems=sscanf(line.c_str(), "%2dd%4lfm%c%3dd%4lfm%c", &lat,&latmin,&c1,&lon,&lonmin,&c2);
  do {
    std::getline(input, line );
    pos=line.find("Number of tidal constituents =");
    } while(line.find("Number of tidal constituents =")==string::npos);
  pos=line.rfind ("=");
  sub = line.substr (pos+1);
  nitems=sscanf(sub.c_str(),"%d",&mgr_station.nwave);
  
  factor=1;
  factor=0.67;
  
  do {
    std::getline(input, line );
    if(line.find("Unit conversion from Z[meters] to Z[psi] Z=Z/")!=string::npos) {
      pos=line.rfind ("/");
      sub = line.substr (pos+1);
      nitems=sscanf(sub.c_str(),"%lf",&factor);
      }
    pos=line.find("Data points for the tidal analysis N");
    if(pos==string::npos) pos=line.find("Data for the tidal analysis N");
    } while(pos==string::npos);
  pos=line.rfind ("=");
  sub = line.substr (pos+1);
  nitems=sscanf(sub.c_str(),"%d",&mgr_station.duree);
  mgr_station.duree/=24;

  do {
    std::getline(input, line );
    } while(line.find("Constituent       Frequency     Amplitude    Greenwich")==string::npos);
  std::getline( input, line );
  
  data=new mgr_data_t[mgr_station.nwave];
  for(k=0;k<mgr_station.nwave;k++) {
    std::getline(input, line );
    nitems=sscanf(line.c_str(), "%s %lf %lf %lf", data[k].constituent.name, &frequency, &(data[k].amp),&(data[k].phi));
    if(strcmp(data[k].constituent.name,"Z0")==0) mgr_station.loc.depth=-data[k].amp;
    data[k].amp*=factor;
    data[k].error=0;
    }
  mgr_station.data=data;
  
  nmgr=1;
  
  mgr.push_back(mgr_station);
  return(nmgr);
}


map<string,int> MgrHarmonic_formats;


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void print_mgrh_formats(void)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  for (map<string,int>::iterator it=MgrHarmonic_formats.begin(); it!=MgrHarmonic_formats.end(); ++it) {
    printf("format name %15s, code=%2d\n",it->first.c_str(),it->second);
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int init_mgrh_formats(void)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
//   int chk;
  
  MgrHarmonic_formats["unknown"]=0;
  
  MgrHarmonic_formats[MGR_FORMAT_NAME_LEGOS_ASCII] =MGR_FORMAT_CODE_LEGOS_ASCII;
  MgrHarmonic_formats[MGR_FORMAT_NAME_LEGOS_NETCDF]=MGR_FORMAT_CODE_LEGOS_NETCDF;
  MgrHarmonic_formats[MGR_FORMAT_NAME_LEGOS_OBS]   =MGR_FORMAT_CODE_LEGOS_OBS;
  MgrHarmonic_formats[MGR_FORMAT_NAME_DART]        =MGR_FORMAT_CODE_DART;
  MgrHarmonic_formats[MGR_FORMAT_NAME_RAY]         =MGR_FORMAT_CODE_RAY;
  MgrHarmonic_formats[MGR_FORMAT_NAME_ATG_DATABASE]=MGR_FORMAT_CODE_ATG_DATABASE;
  MgrHarmonic_formats[MGR_FORMAT_NAME_SHOM]        =MGR_FORMAT_CODE_SHOM;

//   chk=MgrHarmonic_formats[FORMAT_NAME_GLOSS];
//   chk=MgrHarmonic_formats.size();
  
//  print_formats();

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int mgr_load(const char *filename, vector<mgr_t> & mgr, int format, bool reset)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int nmgr,status;
  char wave[10]="";
  
  if(filename==0) {
    printf("please specify input filename\n");
    return(0);
    }
  
  if(reset) mgr.clear();
  
  switch (format) {
    case MGR_FORMAT_CODE_LEGOS_ASCII:
    case MGR_FORMAT_CODE_LEGOS_NETCDF:
      nmgr=mgr_loadLEGOS(filename, mgr);
      return(nmgr);
      break;
      
    case MGR_FORMAT_CODE_DART:
      nmgr=mgr_loadDART(filename, mgr);
      return(nmgr);
      break;
      
    case MGR_FORMAT_CODE_RAY:
      nmgr=mgr_loadRAY(filename, mgr);
      return(nmgr);
      break;
      
    case MGR_FORMAT_CODE_ATG_DATABASE:
      nmgr=mgr_loadATG_DATABASE(filename, mgr);
      return(nmgr);
      break;

    case MGR_FORMAT_CODE_LEGOS_OBS:
      status= mgr_loadobs(filename, wave, &nmgr, mgr, 0, 0, 0);
      return(nmgr);
      break;

    case MGR_FORMAT_CODE_SHOM:
      nmgr=mgr_loadSHOM(filename, mgr);
      return(nmgr);
      break;
      
    default:
      printf("unknown or unspecified input format for %s\n",filename);
      return(0);
      break;
    }
  
  return(-1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int mgr_load(const char *filename, vector<mgr_t> & mgr, bool reset)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int format=MGR_FORMAT_CODE_LEGOS_ASCII;

  status=mgr_load(filename, mgr, format, reset);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int mgr_load(const char *filename, vector<mooring_t> & moorings, vector<spectrum_t> & spectrum, vector<hconstant_t> & constants, int format)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int nmgr;
  vector<mgr_t> mgr;
//   vector<hconstant_t> constants;
//   vector<spectrum_t>  spectrum;
  
  bool Z0=true;
  spectrum_t reference=spectrum_init_ref("ESTUARINE",Z0);
  
  nmgr= mgr_load(filename, mgr, format);
  
  for(int n=0;n<nmgr;n++) {
    mgr_t tmp=mgr[n];
    hconstant_t h=hconstant_t(tmp);
    constants.push_back(h);
    spectrum_t s=spectrum_t(reference,tmp);
    spectrum.push_back(s);
    mooring_t m=mooring_t(tmp);
    moorings.push_back(m);
    }
  
  return(0);
}


int compWaveOmega(const void *a,const void *b){
  const mgr_data_t *a_=(const mgr_data_t *)a;
  const mgr_data_t *b_=(const mgr_data_t *)b;
  return a_->constituent.omega > b_->constituent.omega;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int mgr_save_ascii(const char *filename, const vector<mgr_t> & mgr)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int nmgr=mgr.size();
  FILE *out = NULL;

  out = fopen(filename, "w");
  if(out==0) TRAP_ERR_EXIT(-1,"Unable to open file: %s \n", filename);

  for(size_t i = 0; i < nmgr; i++) {
    if(mgr[i].mindex > 0){
      fprintf(out,"________________________________________________________________________________");
      fprintf(out,"\n \n");
      int number=mgr[i].number;
      if(number < i) number = i;
      fprintf(out," STATION No%5.4d  : %s\n ORIGINE          : %s\n", number, mgr[i].name, mgr[i].origine);
      fprintf(out," LOCALISATION     :%10.3fN %8.3fE %7.2fm Triangle :     0\n", mgr[i].loc.lat, mgr[i].loc.lon, mgr[i].loc.depth);
      fprintf(out," ENREGISTREMENT   : Debut: %s  Fin: %s  Duree:%5d     \n", mgr[i].debut, mgr[i].fin, mgr[i].duree);
      fprintf(out," VALIDATION       : %s\n", mgr[i].validation);
      for (size_t j = 0; j < mgr[i].nwave; j++) mgr[i].data[j].constituent.init();
      qsort(mgr[i].data, mgr[i].nwave, sizeof(*mgr[i].data), compWaveOmega);
      for (size_t j = 0; j < mgr[i].nwave; j++) {
        if(isnan(mgr[i].data[j].error)) {
          mgr[i].data[j].error=1.e+10;
          }
        fprintf(out,"%3d  %-6s        %9.6f   %10.6f   %9.6f\n", mgr[i].data[j].constituent.code, mgr[i].data[j].constituent.name,
                                                                 mgr[i].data[j].amp, mgr[i].data[j].phi, mgr[i].data[j].error);
        }
      }
    }

  fclose(out);

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int mgr_save_BIO(const char *filename, vector<mgr_t> & mgr)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  Dave Greenberg (BIO) format

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int nmgr=mgr.size();
  FILE *out = NULL;
  string wave_file;
  vector<string> waves;
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  get wave list
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  for(size_t i = 0; i < nmgr; i++) {
    if(mgr[i].mindex > 0) {
      for (size_t j = 0; j < mgr[i].nwave; j++) {
        int pos=vpush( (string) mgr[i].data[j].constituent.name, waves);
        }
      }
    }
        
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  save in separate wave files
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  for(int k=0; k<waves.size(); k++) {
    wave_file=waves[k]+"-"+filename;
    out=fopen(wave_file.c_str(), "w");
    printf("%s\n",wave_file.c_str());
    if(out==0) TRAP_ERR_EXIT(-1,"Unable to open file: %s \n", wave_file.c_str());
    for(size_t i = 0; i < nmgr; i++) {
//       string name=mgr[i].name;
      std::stringstream ss;
      ss << i;
      string name="station-"+ss.str();
      if(mgr[i].mindex > 0) {
        int j=mgr[i].wave_index(waves[k].c_str());
        if(j==-1) continue;
        fprintf(out,"%9.3lf %9.3lf %9.3lf %9.1lf %s\n", mgr[i].loc.lon, mgr[i].loc.lat, mgr[i].data[j].amp, mgr[i].data[j].phi, name.c_str());
        }
      }
    fclose(out);
    }

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int mgr_save(const char *filename, vector<mgr_t> & mgr, float seal)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  FILE *out=NULL;
  int i,j,keep;

  out = fopen(filename, "w");
  if(out==0) TRAP_ERR_EXIT(-1,"file opening issue : %s \n", filename);

  for(i=0;i<mgr.size();i++) {
    mgr[i].mindex=1;
    keep=0;
    for (j=0;j<mgr[i].nwave;j++) {
      if(mgr[i].data[j].error > mgr[i].data[j].amp*seal) continue;
      keep=1;
      }
//    if (mgr[i].mindex>0) {
    if (keep==1) {
      for(j=0;j<80;j++) fprintf(out,"_");
      fprintf(out,"\n \n");
      if(mgr[i].number<i)mgr[i].number=i;
      fprintf(out," STATION No%5.4d  : %s\n ORIGINE          : %s\n",mgr[i].number,mgr[i].name,mgr[i].origine);
      fprintf(out," LOCALISATION     :%10.3fN %8.3fE    0.00m Triangle :     0\n",mgr[i].loc.lat,mgr[i].loc.lon);
      fprintf(out," ENREGISTREMENT   : Debut: %s  Fin: %s  Duree: %5d     \n",mgr[i].debut,mgr[i].fin,mgr[i].duree);
      fprintf(out," VALIDATION       : %s\n",mgr[i].validation);
      qsort(mgr[i].data, mgr[i].nwave, sizeof(*mgr[i].data), compWaveOmega);
      for (j=0;j<mgr[i].nwave;j++) {
        if(mgr[i].data[j].error > mgr[i].data[j].amp*seal) continue;
//        if(mgr[i].data[j]->error > seal) continue;
//        fprintf(out,"%3d  %-6s        %9.6f   %10.6f\n",mgr[i].data[j]->code,mgr[i].data[j]->wave,mgr[i].data[j]->amp,mgr[i].data[j]->phi);
        fprintf(out,"%3d  %-6s        %9.6f   %10.6f   %9.6f\n", mgr[i].data[j].constituent.code, mgr[i].data[j].constituent.name,
                                                                 mgr[i].data[j].amp, mgr[i].data[j].phi, mgr[i].data[j].error);
        }
      }
    }

  fclose(out) ;

  return(0);
}


/*###############################################################################################*/

  int mgr_save2D(const char *filename, int nmgr, vector<mgr_t> *mgr, float seal)

/*###############################################################################################*/
{
  FILE *out=NULL;
  int i,j,keep;

  out = fopen(filename, "w");

  if (out == NULL) TRAP_ERR_EXIT(-1,"file opening issue : %s \n", filename);

  for(i=0;i<nmgr;i++) {
    mgr[0][i].mindex=1;
    keep=0;
    for (j=0;j<mgr[0][i].nwave;j++) {
      if(mgr[0][i].data[j].error > mgr[0][i].data[j].amp*seal) continue;
      keep=1;
      }
    if (keep==1) {
      for(j=0;j<80;j++) fprintf(out,"_");
      fprintf(out,"\n \n");
      if(mgr[0][i].number<i)mgr[0][i].number=i;
      fprintf(out," STATION No%5.4d  : %s\n ORIGINE          : %s\n",mgr[0][i].number,mgr[0][i].name,mgr[0][i].origine);
      fprintf(out," LOCALISATION     :%10.3fN %8.3fE    0.00m Triangle :     0\n",mgr[0][i].loc.lat,mgr[0][i].loc.lon);
      fprintf(out," ENREGISTREMENT   : Debut: %s  Fin: %s  Duree: %5d     \n",mgr[0][i].debut,mgr[0][i].fin,mgr[0][i].duree);
      fprintf(out," VALIDATION       : %s\n",mgr[0][i].validation);
      for (j=0;j<mgr[0][i].nwave;j++) {
        if(mgr[0][i].data[j].error > mgr[0][i].data[j].amp*seal) continue;
        fprintf(out,"%3d  %-6s        %9.6f   %10.6f    %9.6f   %10.6f   %9.6f\n", mgr[0][i].data[j].constituent.code, mgr[0][i].data[j].constituent.name,
                                                                 mgr[0][i].data[j].amp, mgr[0][i].data[j].phi,
                                                                 mgr[1][i].data[j].amp, mgr[1][i].data[j].phi, mgr[1][i].data[j].error);
        }
      }
    }

  fclose(out) ;

  return(0);
}

// /*###############################################################################################*/
//
// int mgr_create(int nmgr, vector<mgr_t> & mgr, const spectrum_t& s)
//
// /*###############################################################################################*/
// {
//   int i,k;
//   char *status=NULL;
//
//   if ((*mgr=(mgr_t **) calloc(nmgr,sizeof(mgr_t *)))==NULL) {
//     perror("*mgr");
//     TRAP_ERR_EXIT(-1,"exiting\n");
//     }
//   for(i=0;i<nmgr;i++) {
//     if (((*mgr)[i]=(mgr_t *) calloc(1,sizeof(mgr_t)))==NULL) {
//       perror("(*mgr)[i]");
//       TRAP_ERR_EXIT(-1,"exiting\n");
//       }
//     if (((*mgr)[i]->data=(mgr_data_t *)calloc(s.n,sizeof(mgr_data_t *)))==NULL) {
//       perror("(*mgr)[i]->data");
//       TRAP_ERR_EXIT(-1,"exiting\n");
//       }
//    for(k=0;k<s.n;k++) {
// //      if (((*mgr)[i]->data[k]=(mgr_data_t *)calloc(1,sizeof(mgr_data_t)))==NULL) {
// //        perror("(*mgr)[i]->data");
// //        TRAP_ERR_EXIT(-1,"exiting\n");
// //        }
//       (*mgr)[i]->data[k].code=0;
//       strcpy((*mgr)[i]->data[k].wave,s.waves[k].name);
//       (*mgr)[i]->data[k].amp=0.0;
//       (*mgr)[i]->data[k].phi=0.0;
//       }
//    (*mgr)[i]->nwave=s.n;
//     }
//
//   return(0);
// }
//

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int mgr_create(int nmgr, vector<mgr_t> & mgr, const spectrum_t& s)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  const int& nw = s.n;
  
  for(size_t m = 0; m < nmgr; ++m){
//    gauge = new mgr_t();
    mgr_t gauge;
    gauge.nwave = nw;
    gauge.data = new mgr_data_t[nw];
    gauge.mindex = 1;
    for(size_t k = 0; k < nw; k++) {
      strcpy(gauge.data[k].constituent.name, s.waves[k].name);
      gauge.data[k].constituent = s.waves[k];
      }
//     gauge.loc = XXX ;
//     strcpy(gauge.name, XXX);
    strcpy(gauge.origine,    "IRRELEVANT");
    strcpy(gauge.validation, "IRRELEVANT");
    strcpy(gauge.debut,      "IRRELEVANT");
    strcpy(gauge.fin,        "IRRELEVANT");
    gauge.number = m+1;
//     gauge.track = XXX;
//     gauge.duree = XXX;
    mgr.push_back(gauge);
    }
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int mgr_create(int nmgr, vector<mgr_t> & mgr, const spectrum_t& s, double *lon, double *lat, double *depth)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=mgr_create(nmgr, mgr, s);
  
  for(int n=0;n<nmgr;n++) {
    mgr[n].loc.lon=lon[n];
    mgr[n].loc.lat=lat[n];
    if(depth!=0) mgr[n].loc.depth=depth[n];
    else         mgr[n].loc.depth=0;
    mgr[n].loc.immersion=0;
    mgr[n].loc.units=strdup("degrees");
    }
    
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int mgr_parse(int nmgr, vector<mgr_t> mgr, const spectrum_t& s)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  /// could be useful to parse mgr[] to get constituents within specified spectrum
  /// \todo to be coded
  
  const int& nw = s.n;

  mgr_t **tmp = new mgr_t*[nmgr]; // oversized
  
  for(size_t m = 0; m < nmgr; ++m){
    tmp[m] = new mgr_t[nw];
    tmp[m]->data = new mgr_data_t[nw];

    /// ...
    
  }
  
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int xover_id(xoverbase_t base,xover_t *xover)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int found,n;
  double t,dx,dy,d,dmin;

  dmin=1.e+10;

  for (n=0;n<base.n;n++) {
    t=degree_recale(xover->t,base.xover[n].t);
    dx=t-base.xover[n].t;
    dy=xover->p-base.xover[n].p;
    d=dx*dx+dy*dy;
    if(d<dmin) {
      found=n;
      dmin=d;
      }
    }
/* *----------------------------------------------------------------------------
  verification*/
  d=geo_km(xover->t,xover->p,base.xover[found].t,base.xover[found].p);
  if(d>20.) {
    xover->id=-1;
    xover->track[0]=-1;
    xover->track[1]=-1;
    return(-1);
    }
  else {
    n=found;
    xover->id=base.xover[n].id;
    xover->track[0]=base.xover[n].track[0];
    xover->track[1]=base.xover[n].track[1];
    return(0);
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  xoverbase_t read_xover(char *input)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int n,nitems;
  FILE *in=NULL;
  xoverbase_t base;
  
  in=fopen(input, "r");
  if (in == NULL) TRAP_ERR_EXIT(-1,"cannot open file %s (%d %s)\n", input, errno, strerror(errno));

  fscanf(in,"%d",&n);

  base.xover=new xover_t[n];
  base.n=n;

  for (n=0;n<base.n;n++) {
    nitems=fscanf(in,"%d %d %d %lf %lf",&(base.xover[n].id),
                                        &(base.xover[n].track[0]),&(base.xover[n].track[1]),
                                        &(base.xover[n].t),&(base.xover[n].p));
    }

  fclose(in);

  return(base);
}


// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
//
//   int *order(vector<mgr_t> mgr,int nmgr, char *criterion, int direction)
//
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   int *list;
//   int *position;
//   int m,n;
//   char sign;
//
//   list    = new int [nmgr];
//   position= new int [nmgr];
//
//   for(n = 0; n < nmgr; n++) {
//     position[n]=0;
//     }
//
//   for(n = 0; n < nmgr; n++) {
//     for(m = n+1; m < nmgr; m++) {
//       if(strcasecmp(check_string(mgr[n].name),check_string(mgr[m].name))>0) {
//         position[n]++;
//         }
//       else {
//         position[m]++;
//         }
//       }
//     }
//
//   for(n = 0; n < nmgr; n++) {
//     list[position[n]]=n;
//     }
//
//   return (list);
//   }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int *mgr_order(vector<mgr_t> mgr, int nmgr)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int *list=NULL;
  int *position=NULL;
  int n;

  list    = new int [nmgr];
  position= new int [nmgr];

  for(n = 0; n < nmgr; n++) {
    position[n]=n;
    }

  for(n = 0; n < nmgr; n++) {
    list[position[n]]=n;
    }

  delete [] position;

  return (list);
  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int *mgr_order(vector<mgr_t> mgr,int nmgr,const char *criterion, int direction)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int *list=NULL;
  int *position=NULL;
  int m,n;

  list    = new int [nmgr];
  position= new int [nmgr];
  
  for(n = 0; n < nmgr; n++) {
    position[n]=0;
    }

  if(strcmp(criterion,"ALPHA")==0) {
    char **names=new char*[nmgr];
    for(n = 0; n < nmgr; n++) {
      names[n]=strdup(check_string(mgr[n].name));
      }
    for(n = 0; n < nmgr; n++) {
      for(m = n+1; m < nmgr; m++) {
        if(strcasecmp(names[n],names[m])>0) {
          position[n]++;
          }
        else {
          position[m]++;
          }
        }
      }
    for(n = 0; n < nmgr; n++) {
      free(names[n]);
      }
    delete[] names;
    }
    
  if(strcmp(criterion,"LATITUDE")==0) {
    for(n = 0; n < nmgr; n++) {
      for(m = n+1; m < nmgr; m++) {
        if(mgr[n].loc.lat > mgr[m].loc.lat) {
          position[n]++;
          }
        else {
          position[m]++;
          }
        }
      }
    }
  if(strcmp(criterion,"LONGITUDE")==0) {
    for(n = 0; n < nmgr; n++) {
      for(m = n+1; m < nmgr; m++) {
        if(mgr[n].loc.lon > mgr[m].loc.lon) {
          position[n]++;
          }
        else {
          position[m]++;
          }
        }
      }
    }
  if(strcmp(criterion,"DEPTH")==0) {
    for(n = 0; n < nmgr; n++) {
      for(m = n+1; m < nmgr; m++) {
        if(mgr[n].loc.immersion > mgr[m].loc.immersion) {
          position[n]++;
          }
        else {
          position[m]++;
          }
        }
      }
    }
  else if(strcmp(criterion,"FILE")==0) {
    for(n = 0; n < nmgr; n++) {
      position[n]=n;
      }
    }

  if(direction==1) {
    for(n = 0; n < nmgr; n++) {
      list[position[n]]=n;
      }
    }
  else {
    for(n = 0; n < nmgr; n++) {
      list[nmgr-position[n]-1]=n;
      }
    }

  delete [] position;

  return (list);
}
