
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

/* *@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  topo-merg aims at merging 2 bathymetric database with some controls
  on merging action

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */

#define MAIN_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include "tools-structures.h"

#include "rutin.h"
#include "geo.h"
#include "polygones.h"
#include "netcdf-proto.h"
#include "functions.h"
#include "grd.h"
#include "map.h"
#include "xyz.h"
#include "statistic.h"
#include "legend.h"

#include "topo.h"
#include "vector"
#include "iostream"

#define XYZ 0
#define YXZ 1


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int xyz_loadCORA (const char* filename, char *proj4_options, double * &x,double * &y,double * &z, double & mask, int & ndata)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
// CO_DMQCGL01_19900101_PR_XB.nc ; KOELN ATLANTIC                                                   ; 999   ;  10. ;   21. ;  -38.350 ;   49.600 ;   862
  FILE *in;
  int nitems;
  int i,j,k,l,m,n,range;
  int status;
  char line[1024];
  char file[128], chief[128], c[128];
  float dum;
  string s="";
  size_t pos;

  status=0;
  
  ndata=xyz_countdata (filename);

  x=new double[ndata];
  y=new double[ndata];
  z=new double[ndata];

  in=fopen(filename,"rb");
  if (in == NULL) {
    return(-1);
    }
  for(n=0;n<ndata;n++) {
    fgets(line, 1024, in);
    s=(string) line;
    
    pos=s.find(";");
    s.copy(file, pos-1, 0);
    s.erase(0,pos+1);
    
    pos=s.find(";");
    s.copy(chief, pos-1, 0);
    s.erase(0,pos+1);
    
    for(k=0;k<3;k++) {
      pos=s.find(";");
      s.erase(0,pos+1);
      }
    
    pos=s.find(";");
    nitems=sscanf(s.c_str(),"%lf", &x[n]);
    s.erase(0,pos+1);
    
    pos=s.find(";");
    nitems=sscanf(s.c_str(),"%lf", &y[n]);
    s.erase(0,pos+1);
    
    pos=s.find(";");
    nitems=sscanf(s.c_str(),"%lf", &z[n]);
    s.erase(0,pos+1);
    
//     printf("%lf %lf %lf \n",x[n], y[n], z[n]);
//     printf("%s \n",line);
    }

  status=xyz_save("echo.xyz",x,y,z,mask,ndata);
  
  fclose(in);
  exit(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int compare (char *filename, char *proj4_options, grid_t grid, float *topo, float zmask, float scale, float zmin, float zmax)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int option, status, count, format;
  int k,m,n,ndata,ncols=12, i,j;
  double *x,*y,*z,**buf, mask,t,p,mean,variance, tab;
  float zz, remainder;
  statistic_t s;
  range_t<double> r;
  legend_t   *legends=NULL;
  string header="";
 
  char *ss=strstr(filename,".shp");
  
  if(ss==NULL) {
    format=ASCII;
    }
  else{
    format=SHAPE;
    }
  
  switch (format) {
    case ASCII:
      status=xyz_loadraw (filename, header, proj4_options, x, y, z, &mask, ndata);
      break;
    case SHAPE:
      status=xyz_load_shp (filename, proj4_options, x, y, z, &mask, ndata);
      break;
    }

//  status=xyz_loadraw (filename, proj4_options, x, y, z, &mask, ndata);
//  status=xyz_loadCORA (filename, proj4_options, x, y, z, mask, ndata);
    
  r=range_t<double>(x,ndata);
  r.print("x range");
  
  r=range_t<double>(y,ndata);
  r.print("y range");
  
  r=range_t<double>(z,ndata);
  r.print("z range");
  
  printf("#####################################################################\n");
  printf("compare bathymetry MNT with data: \n");
  printf("positive difference means data above MNT floor (quite feasible)\n");
  printf("negative difference means data below MNT floor (MNT is wrong?) \n");
  buf=new double* [ncols];
  for(k=0;k<ncols;k++) {
    buf[k]=new double [ndata];
    }
  
  for(n=0;n<ndata;n++) {
    t=x[n];
    p=y[n];
    t=map_recale(grid,t);
    status=map_interpolation(grid,topo,zmask,t,p,&zz);
    if(zz!=zmask && zz< zmax && zz> zmin) { 
//       buf[0][n]=-z[n]*scale;
      buf[0][n]=z[n]*scale;
      buf[1][n]=zz;
/*------------------------------------------------------------------------------
      data minus MNT : positive mean if MNT deeper than data */
      buf[2][n]=buf[0][n]-buf[1][n];
      if(fabs(buf[1][n])<2.0) {
        buf[3][n]=mask;
        }
      else{
        buf[3][n]=100.*buf[2][n]/fabs(buf[1][n]);
        }
//       if(abs(buf[3][n]) >200.) {
//         printf("%d \n",n);
//         }
      if(fabs(buf[2][n])>0.5) {
        printf("%s : %lf %f\n",__FUNCTION__,zz,z[n]);
        }
      }
    else {
      buf[0][n]=-z[n];
      buf[1][n]=mask;
      buf[2][n]=mask;
      buf[3][n]=mask;
      }
    }
    
  r=poc_minmax(buf[2], ndata, mask);
  s=get_statistics(buf[2], mask, ndata);
 
  r=poc_minmax(buf[3], ndata, mask);
  s=get_statistics(buf[3], mask, ndata);
 
  for(n=0;n<ndata;n++) {
    if(abs(buf[3][n])>25) {
      buf[4][n]=0;
      }
    else {
      buf[4][n]=1;
      }
    remainder=(int) buf[0][n]%100;
    if(remainder==0) {
      buf[5][n]=0;
      }
    else {
      buf[5][n]=1;
      }
/** a revoir */
//     s.std=buf[3][n];
//     if(s.std=0) {
//       buf[6][n]=0;
//       }
//     else {
//       buf[6][n]=1;
//       }
    }
//    get_statistics(
   
  for(n=0;n<ndata;n++) {
    if(buf[4][n]==0) {
      buf[9][n]=mask;
      }
    else {
      buf[9][n]=buf[3][n];
      }
    }
 
  s=get_statistics(buf[9], mask, ndata);
  
  
  int sum=0;
  for(int n=0;n<ndata;n++) {
    if(buf[4][n]==1) sum++;
    }
 
  
 
  vector<int> v(369551);
  //int s=v.size(); //unused and conflicting
  for (int i=0;i<v.size();i++) {
    if (buf[4][n]==1) {
      v[n]=buf[3][n];
      }
    }
  
   
       
   
   s=get_statistics(buf[9], mask, ndata);

   
   
    
  status=xyz_save ("topo-compare.xyz", x, y, buf, mask, ndata, ncols);
  
  legends=lgd_import(x, y, buf, ndata, ncols);
  status=lgd_save("topo-compare.lgd", legends, 1, NULL, NULL);
  
  
    


  for(n=0;n<ndata;n++) {
    if(buf[3][n]<0) {
      buf[10][n]=buf[3][n];
      }
    else {
      buf[10][n]=mask;
      }
    }
    
  
//   for(n=0;n<ndata;n++) {
//     if(buf[10][n]==mask) {
//      for(m=n;m<ndata-1;m++) {
//        for(k=0;k<11;k++){
//          buf[k][m]=buf[k][m+1];
//          }
//        }
//      n--;
//      ndata--;
//      }
//    }

//  int count;
  count=0;
  for(n=0;n<ndata;n++) {
    if(buf[10][n]!=mask) count=count+1;
    }

  int *incidence=new int[count];
  for(n=0;n<count;n++) incidence[n]=-1;
  
  count=0;
  for(n=0;n<ndata;n++) {
    if(buf[10][n]!=mask) {
     incidence[count]=n;
     count++;
     }
   }
  for(n=0;n<count;n++) {
    m=incidence[n];
    x[n]=x[m];
    y[n]=y[m];
    for(k=0;k<11;k++){
      buf[k][n]=buf[k][m];
      }
    }
   
  legends=lgd_import(x, y, buf, count, ncols);
  status=lgd_save("topo-anomalies.lgd", legends, 1, NULL, NULL);

  delete[] x;
  delete[] y;
  delete[] z;
  
  return(status);
  }
  

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  float  scale=1.0,topo_mask=1.e+10,mask=1.e+10;
  double  x,y,t,p;
  double  dx,dy;

  int i,j,k,l,m,n,signus=0,masked_only=0,status;

  size_t size;
  FILE *in,*out;
  char *keyword,*zone;
  int flag;

  char *rootname=NULL,*output=NULL,*input=NULL,*pathname=NULL,*bathymetry=NULL;
  char *bathy=NULL,*format=NULL,*hydro=NULL,*poly=NULL;

  grid_t topogrid;
  grid_t grid,zone_grid;
  mesh_t mesh;
  short *topo,smask=256*127+255;
  float *ftopo,*tmp,*buffer,ftopomask;

  float zmin=0,zmax=0;

  plg_t *polygones=NULL;
  int npolygones=0;

  pocgrd_t ncgrid;
  cdfvar_t variable;
  int dimlgth[4]={10,20,30,40};
  char *proj4=0;
  
  bool debug=false;

  fct_echo( argc, argv) ;

  n=1;
  while (n < argc) {
    keyword=strdup(argv[n]);
    if(strcmp(keyword,"-inv")==0) {
      scale=-1;
      n++;
      continue;
      }
    if(strcmp(keyword,"-zmin")==0) {
      sscanf(argv[n+1],"%f",&zmin);
      n++;
      n++;
      continue;
      }
    if(strcmp(keyword,"-zmax")==0) {
      sscanf(argv[n+1],"%f",&zmax);
      n++;
      n++;
      continue;
      }
    if(strcmp(keyword,"-debug")==0) {
      debug=true;
      n++;
      continue;
      }
    if(strstr(argv[n],"--proj=")!=0){
      proj4=strdup(argv[n]+7);
      n++;
      continue;
      }
    switch (keyword[0]) {
      case '-':
        switch (keyword[1]) {
          case 'b' :
            bathymetry= strdup(argv[n+1]);
            n++;
            n++;
            break;

          case 'f' :
            format= strdup(argv[n+1]);
            n++;
            n++;
            break;

          case 'p' :
            poly= strdup(argv[n+1]);
            n++;
            n++;
            break;

          case 'o' :
            rootname= strdup(argv[n+1]);
            n++;
            n++;
            break;

          case 'z' :
            zone= strdup(argv[n+1]);
            n++;
            n++;
            break;

          default:
            printf("unknown option %s\n",keyword);
            __OUT_BASE_LINE__("unknown option %s\n",keyword);
            exit(-1);
            break;
          }
        break;

      default:
        if(input==NULL) {
          input= strdup(argv[n]);
          n++;
          }
        else {
          __OUT_BASE_LINE__("unknown option %s\n",keyword);
          exit(-1);
          }
        break;
      }
      free(keyword);
    }

    
   
/* *----------------------------------------------------------------------
  load bathymetry grid*/
  printf("#####################################################################\n");
  printf("load working bathymetry: %s\n",bathymetry);
//  status=map_loadfield(bathymetry, (char *) 0, &topogrid, &ftopo, &ftopomask);
  status=topo_loadfield(bathymetry, &topogrid, &ftopo, &ftopomask, debug);
  if(status !=0) {
    printf("cannot load bathymetry file=%s\n",bathymetry);
    goto error;
    }


/* *----------------------------------------------------------------------
  import new depths*/
  if(input!=NULL) {
    printf("\n#####################################################################\n");
    printf("load bathymetry soundings: %s\n", input);
    status= compare (input, proj4, topogrid, ftopo, ftopomask, scale, zmin, zmax);
    if(status !=0) {
      goto error;
      }
    }

// /* *------------------------------------------------------------------------------
//   apply change of sign on final topography*/
//   for (j=0;j<topogrid.ny;j++) {
//     for (i=0;i<topogrid.nx;i++) {
//       m=j*topogrid.nx+i;
//       if(ftopo[m]!=ftopomask) {
//         ftopo[m]*=scale;
//         }
//       }
//     }
//  
//   printf("#####################################################################\n");
//   printf("save bathymetry\n");

  delete[] ftopo;
  topogrid.free();

  __OUT_BASE_LINE__("topo_merge sucessfully completed\n");

  exit(0);

 error:
  __OUT_BASE_LINE__("topo_merge aborted\n");
  exit(-1);
}
