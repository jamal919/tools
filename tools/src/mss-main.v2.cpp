#define MAIN_SOURCE

#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "tools-define.h"
#include "tools-structures.h"

#include "matrix.h"
#include "rutin.h"
#include "geo.h"
#include "map.h"
#include "grd.h"
#include "bmg.h"
#include "polygones.h"
#include "tides.h"
#include "filter.h"
#include "mss.v2.h"


#define nstat 12

static bool expired=expire(20130922,20131122);

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
double DEG2RAD = M_PI/180.0;
double RAD2DEG = 180.0/M_PI;

double rTerre = 6378.1363e3;
double omegaTerre = 7.292115e-5;
double j2Terre = 1082.63622070e-6;
double fTerre = 298.257;
double demiGdAxe = 7714.43e3;
double moyenMouv = 2.0*M_PI/6745.72;
double inclin = 66.039*DEG2RAD;
double omegaP = -3.0 / 2.0 * moyenMouv * j2Terre * ( rTerre / demiGdAxe )*( rTerre / demiGdAxe ) * cos( inclin );

  float  *buf=NULL,*error_cls=NULL,z;
  float ftime,mask=1.e+10,dum;
  float  *mss_cls=NULL,*analysis=NULL,*innovation=NULL,*G=NULL,*tmp=NULL;
  float  *mss1=NULL,*mss2=NULL,*mss3=NULL,*mss4=NULL;

  double  x,y;
  double *lon=NULL,*lat=NULL,*error=NULL;

  int i,j,k,m,n;
  int track,status;

  FILE *in;
  char *keyword=NULL, *sat=NULL, *mss_file=NULL,*output=NULL, *alti_file=NULL;
  char sens[3];

  grid_t grid,mssgrid,*nominal_track=NULL;

  serie_t TPdata;
  
  fct_echo(argc,argv);
 
  n=1;
  while (n < argc)
    {
    keyword=strdup(argv[n]);
    if(strcmp(argv[n],"-help") ==0) {
      __OUT_BASE_LINE__("-s satellite -i input_file - t track_number -m mss_bmg_file\n");
      exit(1);
      continue;
      }
    switch (keyword[0]) {
      case '-':
        switch (keyword[1]) {
        case 's' :
          sat= strdup(argv[n+1]);
          printf("the satellite is : %s\n",sat);
         n++;
          n++;
          break;

        case 'i' :
          alti_file= strdup(argv[n+1]);
          printf("the altimetric data are in file is : %s\n",alti_file);
         n++;
          n++;
          break;

        case 't' :
          track= atoi(argv[n+1]);
          if(track%2==0)printf("descending track !!!!! or ascending in the case of gfo : %s\n");
          else printf("ascending track !!!!! or descending in the case of gfo : %s\n");
          n++;
          n++;
          break;


        case 'm' :
          mss_file= strdup(argv[n+1]);
          printf("the mss data are in file is : %s\n",mss_file);
          n++;
          n++;
          break;

        default:
          __OUT_BASE_LINE__("unknown option %s\n",keyword);
          exit(-1);
        }
        break;

      default:
          __OUT_BASE_LINE__("unknown option %s\n",keyword);
          exit(-1);

        break;
      }
    free(keyword);keyword=NULL;
    }

  output=(char *)malloc(1024);


  //check the ascending or descending type
  if(strcmp(sat,"gfo")!=0)
    if(track%2==0) sprintf(sens,"DES"); else sprintf(sens,"ASC");
  else
    if (track%2==0) sprintf(sens,"ASC");
    else sprintf(sens,"DES"); //the GFO case


/*-----------------------------------------------------------------------------
  Read altimeter data */
  TPdata=load_metadata_raw(alti_file);
  printf("readed the altimetric data\n");
  if(TPdata.count==0) {__ERR_BASE_LINE__("exiting\n");exit(-1);}

/*-----------------------------------------------------------------------------
  Select ascending or descending, depends on satellite:
  a GFO ascending track is different than a T/P ascending one!! */
  lon=(double *) calloc(2,sizeof(double));
  lat=(double *) calloc(2,sizeof(double));
  get_extremities(sat,track,TPdata,lon,lat);

/*------------------------------------------------------------------------------
  Create mss grid */

  double interval,longueur;
  double *dtime;

  nominal_track=(grid_t *)calloc(1,sizeof(grid_t));
  interval=5.9;  //interval entre les points de la trace nominale
  
  longueur=geo_km(lon[0],lat[0],lon[1],lat[1]);
  longueur/=interval;

  dtime=(double *)calloc((int)(longueur),sizeof(double));
  nominal_track->x=(double *)calloc((int)(longueur),sizeof(double));
  nominal_track->y=(double *)calloc((int)(longueur),sizeof(double));

  nominal_track->ny=(int)(longueur);
  nominal_track->nx=(int)(longueur);

  for(i=0;i<(int)(longueur);i++)dtime[i]=i;

  interpTrace( sens, lon[0]*DEG2RAD, lat[0]*DEG2RAD, (int)(longueur), dtime,nominal_track->x, nominal_track->y );
  for(i=0;i<(int)(longueur);i++){nominal_track->x[i]*=RAD2DEG;nominal_track->y[i]*=RAD2DEG;}

  
  //plus de grille a gerer donc plus besion de passer dans sioux
  //mssgrid=get_trackgridsioux(lon,lat,TPdata,nominal_track,sens);
  //sprintf(output,"%s.grided_mss.nc",alti_file);
  //status=mss_createfile(output,(size_t) mssgrid.nx,(size_t) mssgrid.ny, mssgrid);
  
  free(lon);lon=NULL;
  free(lat);lat=NULL;
  
  
  /*------------------------------------------------------------------------------
    Load CLS mss */
  //   status=bmg_loadgrid (mss_file,&grid);
  //   grid.mode=0;
  //   status=map_completegridaxis_2(&grid);
  
  //   buf=(float *)malloc(grid.nx*grid.ny*sizeof(float));
  //   if (buf == NULL) {__ERR_BASE_LINE__("exiting\n");exit(8);}
  //   error_cls=(float *)malloc(grid.nx*grid.ny*sizeof(float));
  //   if (error_cls == NULL)exit(11) ;
  
  //   status= bmg_loadr1 (mss_file,1,1,1,grid,buf,&mask,&ftime);
  //   if (status != 0) exit(9) ;
  
  //   status= bmg_loadr1 (mss_file,2,1,1,grid,error_cls,&mask,&ftime);
  //   if (status != 0) exit(10) ;
  
  /*------------------------------------------------------------------------------
    Interpolate and store mss on track grid*/
  
  //     mss_cls=(float *)malloc( mssgrid.nx*mssgrid.ny*sizeof(float));
  //   for(j=0;j<mssgrid.ny;j++)
  //     for(i=0;i<mssgrid.nx;i++)
  //       {
  //       x=mssgrid.x[j*mssgrid.nx+i];
  //       if(x<grid.xmin)x+=360;
  //       y=mssgrid.y[j*mssgrid.nx+i];
  //       status=map_interpolation(grid, buf,mask,x,y,&z);
  //       if (z!=mask) mss_cls[j*mssgrid.nx+i]=z/1000.;
  //       else mss_cls[j*mssgrid.nx+i]=1.0e+35;
  //       }
  
  //   sprintf(output,"%s.grided_mss.nc",alti_file);
  //   status=mss_createfile(output,(size_t) mssgrid.nx,(size_t) mssgrid.ny, mssgrid);
  
  /*-----------------------------------------------------------------------------
    compute the optimal mss*/
  //   error=(double *)calloc(TPdata.count,sizeof(double));
  //   innovation=(float *)calloc(mssgrid.nx*mssgrid.ny,sizeof(float));
  //   analysis=(float *)calloc(mssgrid.nx*mssgrid.ny,sizeof(float));
  //   if( (G =(float *)calloc(mssgrid.nx*mssgrid.ny*TPdata.count,sizeof(float))) == NULL) gmerror("can not allocate G matrix");
  
  //   compute(mssgrid,TPdata,mss_cls,error,innovation,analysis,G);
  
  //   free(error);
  //   error=NULL;
  
  /*-----------------------------------------------------------------------------
    save optimal mss*/
  
  //   for(k=0;k<mssgrid.nx*mssgrid.ny;k++) innovation[k]*=100; /* innovtion in cm */
  //   sprintf(output,"%s.grided_mss.nc",alti_file);
  //   status=mss_createfile(output,(long unsigned int) mssgrid.nx,(long unsigned int) mssgrid.ny, mssgrid);
  //   status=writefile(output, analysis,mss_cls,innovation);
  //   free(innovation);
  //   innovation=NULL;
  
  /*-----------------------------------------------------------------------------*/
  //construction des differentes mss
  
  //   mss1=(float *)malloc(TPdata.count*sizeof(float));
  //   mss2=(float *)malloc(TPdata.count*sizeof(float));
  mss3=(float *)malloc(TPdata.count*sizeof(float));
  //mss4=(float *)malloc(TPdata.count*sizeof(float));
  
  mask=9999.9;
  for(k=0;k<TPdata.count;k++) {
      dum=0.;
      //       for(m=0;m<mssgrid.nx*mssgrid.ny;m++)
      // 	dum+=G[m*TPdata.count+k]*analysis[m];
      //       mss1[k]=dum;
      //       dum=0.;
      //       for(m=0;m<mssgrid.nx*mssgrid.ny;m++)
      // 	dum+=G[m*TPdata.count+k]*mss_cls[m];
      //       mss2[k]=dum;
      mss3[k]=TPdata.data[k].values[24];
      //       x=TPdata.data[k].lon;
      //       if(x<grid.xmin)x+=360.0;
      //       if(x>grid.xmax)x-=360.0;
      //       y=TPdata.data[k].lat;
      //       status=map_interpolation(grid, buf,mask,x,y,&dum);
      //       mss4[k]=dum/1000.0;
    }
  /*-----------------------------------------------------------------------------
    create time series on mean tracks */
  /*   save_meantracks_CLS(TPdata,xover,mss1,polygones,npolygones);  */
  save_meantracks_CTOH(TPdata,alti_file,track,mss3,nominal_track);
  

//  free(analysis);
//   free(mss_cls);
//   free(buf);
//   free(error_cls);
//  free(mss1);
//  free(mss2);
  free(mss3);
  //free(mss4);
  __ERR_BASE_LINE__("exiting\n");exit(0);
  
}
