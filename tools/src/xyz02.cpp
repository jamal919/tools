
/**************************************************************************

  T-UGOm hydrodynamic ocean model, 2006

  Unstructured Ocean Grid initiative

Contributors:

  Florent Lyard      LEGOS/CNRS, Toulouse, France
  Laurent Roblou     LEGOS/CNRS, Toulouse, France
  Yoann Le Bars      LEGOS, Toulouse, France
  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
  Frédéric Dupont    Université de Laval à Québec, Canada

***************************************************************************/
#include <stdio.h>
#include <string.h>


#include "tools-structures.h"

#include "rutin.h"
#include "geo.h"
#include "polygones.h"
#include "grd.h"
#include "map.h"
#include "filter.h"
#include "statistic.h"
#include "functions.h"
#include "sym-io.h"
#include "xyz.h"
#include "polygones.h"


#define XYZ 0
#define YXZ 1

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int xyz_skipheader(FILE *in)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int found,k,l,m,n,status;
  char line[1024];

  status=0;

  if (in == NULL) {
    return(-1);
    }

  fgets(line, 1024, in);
  found=(strstr(line, "HEADER")!=0);

  if(found==0) {
    rewind(in);
    return(0);
    }
  else {
     found=0;
     while (found==0) {
       fgets(line, 1024, in);
       found=(strstr(line, "HEADER")!=0);
       if (feof(in)) return(-1);
      }
    }
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int xyz_readgrid_00_from_v1 (const char* filename,grid_t *grid, int mode)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  FILE *in;
  int nitems,ndata,nlon,nlat;
  int i,j,k,l,m,n;
  int downward=0,cartesian=0,complete=0,option,status;
  char line[1024];
  double *x,*y,*z,zz,mask=-9999.0;
  double *xtics,*ytics;
  range_t<double> r;

  status=0;

  in=fopen(filename,"rb");
  if (in == NULL) {
    return(-1);
    }
  status=xyz_skipheader(in);

  fgets(line, 1024, in);

  n=-1;
  while (!feof(in)) {
    fgets(line, 1024, in);
    n++;
    }
  ndata=n;

  x=new double[ndata];
  y=new double[ndata];
  z=new double[ndata];

  rewind (in);

  fgets(line, 1024, in);

  switch (mode) {
    case 0:
      for(n=0;n<ndata;n++) {
        fgets(line, 1024, in);
        nitems=sscanf(line,"%lf %lf %lf",&x[n],&y[n],&zz);
        }
      break;

    case 1:
      for(n=0;n<ndata;n++) {
        fgets(line, 1024, in);
        nitems=sscanf(line,"%lf %lf %lf",&y[n],&x[n],&zz);
        }
      break;
     }

  fclose(in);

  switch (mode) {
    case 0:
      xtics=new double[ndata];
      xtics[0]=x[0];
      nlon=1;
      for(n=0;n<ndata-1;n++) {
        m=n+1;
        while(x[n]==x[m]) {
          if(m==ndata-1) break;
          m++;
          }
        nlon++;
        xtics[nlon-1]=x[n];
        n=m;
        }
      nlon--;
      break;
      ytics=new double[ndata];
      ytics[0]=x[0];
      nlat=1;
      for(n=0;n<ndata-1;n++) {
        k=0;
        while(y[n]!=ytics[k]) {
          k++;
          if(k==nlat) break;
          }
        if(k==nlat) {
          nlat++;
          ytics[nlat-1]=y[n];
          }
        }
//      nlat--;
    case 1:
      ytics=new double[ndata];
      ytics[0]=y[0];
      nlat=1;
      for(n=0;n<ndata-1;n++) {
        m=n+1;
        while(y[n]==y[m]) {
          if(m==ndata-1) break;
          m++;
          }
        nlat++;
        ytics[nlat-1]=y[n];
        n=m;
        }
      nlat--;
      xtics=new double[ndata];
      xtics[0]=x[0];
      nlon=1;
      for(n=0;n<ndata-1;n++) {
        k=0;
        while(x[n]!=xtics[k]) {
          k++;
          if(k==nlon) break;
          }
        if(k==nlon) {
          nlon++;
          xtics[nlon-1]=x[n];
          }
        }
//      nlon--;
      break;
     }

  reorder(xtics,nlon);
  reorder(ytics,nlat);

  grid->xmin=xtics[0];
  grid->xmax=xtics[nlon-1];

  grid->ymin=ytics[0];
  grid->ymax=ytics[nlat-1];

  grid->nx=nlon;
  grid->ny=nlat;
  grid->nz=1;

  grid->dx=(grid->xmax-grid->xmin) / ((double) grid->nx-1.);
  grid->dy=(grid->ymax-grid->ymin) / ((double) grid->ny-1.);

  grid->x=new double[nlon*nlat];
  grid->y=new double[nlon*nlat];
//  grid->z=new double[ndata];

  for (j=0;j<grid->ny;j++) {
    for (i=0;i<grid->nx;i++) {
      m=j*grid->nx+i;
      grid->x[m]=xtics[i];
      grid->y[m]=ytics[j];
//      grid->z[m]=y[n];
      }
    }

  grid->modeH=2;

  delete[] x;
  delete[] y;
  delete[] xtics;
  delete[] ytics;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int xyz_readgrid_01_from_v1 (const char* filename,grid_t *grid, int mode)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  Assume full (no missing values) array structured (sequential array's axes
  ordering) description

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  FILE *in;
  int nitems,ndata,nlon,nlat;
  int i,j,k,l,m,n,range;
  int ii,jj,mm;
  int **list,*count,*tmp,ntest;
  int downward=0,cartesian=0,complete=0,option,status;
  char line[1024];
  double *correlation,*x,*y,*z,zz,mask=-9999.0,xx,yy,mean_x,mean_y;
  range_t<double> r;
  irange_t s;
  int proxies[9][2]={1,0,1,1,0,1,-1,1,-1,0,-1,-1,0,-1,1,-1,0,0};

  status=0;

  in=fopen(filename,"rb");
  if (in == NULL) {
    return(-1);
    }
  status=xyz_skipheader(in);

  fgets(line, 1024, in);

  n=-1;
  while (!feof(in)) {
    fgets(line, 1024, in);
    n++;
    }
  ndata=n;

  x=new double[ndata];
  y=new double[ndata];
  z=new double[ndata];

  rewind (in);

  status=xyz_skipheader(in);

  fgets(line, 1024, in);

  mean_x=0.0;
  mean_y=0.0;

  switch (mode) {
    case 0:
      for(n=0;n<ndata;n++) {
        fgets(line, 1024, in);
        nitems=sscanf(line,"%lf %lf %lf",&x[n],&y[n],&z[n]);
//        mean_x+=x[n];
//        mean_y+=y[n];
        }
      break;

    case 1:
      for(n=0;n<ndata;n++) {
        fgets(line, 1024, in);
        nitems=sscanf(line,"%lf %lf %lf",&y[n],&x[n],&z[n]);
        }
      break;
     }

  fclose(in);

  range=(int) 10*sqrt((double)ndata);

  correlation=get_autocorrelation(y,range*5,range,ndata);
//  correlation=get_autocorrelation(x,range*5,range,ndata);

  nlon=sqrt((double) ndata)/3;
  nlat=sqrt((double) ndata)/3;

  poc_minmax(x,ndata,mask,&grid->xmin,&grid->xmax);
  poc_minmax(y,ndata,mask,&grid->ymin,&grid->ymax);

  grid->nx=nlon;
  grid->ny=nlat;
  grid->nz=1;

  grid->dx=(grid->xmax-grid->xmin) / ((double) grid->nx-1.);
  grid->dy=(grid->ymax-grid->ymin) / ((double) grid->ny-1.);

  grid->xmin-=2*grid->dx;
  grid->xmax+=2*grid->dx;
  grid->ymin-=2*grid->dy;
  grid->ymax+=2*grid->dy;
  grid->nx+=4;
  grid->ny+=4;

  grid->x=new double[grid->Hsize()];
  grid->y=new double[grid->Hsize()];
  grid->z=new double[grid->Hsize()];

  for (j=0;j<grid->ny;j++) {
    for (i=0;i<grid->nx;i++) {
      m=j*grid->nx+i;
      grid->x[m]=grid->xmin+(double) i*grid->dx;
      grid->y[m]=grid->ymin+(double) j*grid->dy;
      grid->z[m]=0;
      }
    }
  grid->modeH=0;

/*------------------------------------------------------------------------------
  count number of points in grid cell*/
  count=new int[grid->Hsize()];
  for (j=0;j<grid->ny;j++) {
    for (i=0;i<grid->nx;i++) {
      m=j*grid->nx+i;
      count[m]=0;
      }
    }

  for(n=0;n<ndata;n++) {
    status=map_index(*grid,x[n],y[n],&i,&j);
    m=j*grid->nx+i;
    grid->z[m]++;
    count[m]++;
    }

/*------------------------------------------------------------------------------
  create list of points in grid cell*/
  list=new int*[grid->Hsize()];
  for (j=0;j<grid->ny;j++) {
    for (i=0;i<grid->nx;i++) {
      m=j*grid->nx+i;
      list[m]=new int[count[m]];
      count[m]=0;
      }
    }

  for(n=0;n<ndata;n++) {
    status=map_index(*grid,x[n],y[n],&i,&j);
    m=j*grid->nx+i;
    list[m][count[m]]=n;
    count[m]++;
    }

  s=poc_minmax(count,grid->Hsize(),-1);
  tmp=new int[9*s.max];

/*------------------------------------------------------------------------------
  detect direct neighbours*/
//   for(n=0;n<ndata;n++) {
//     status=map_index(*grid,x[n],y[n],&i,&j);
//     m=j*grid->nx+i;
//     ntest=0;
//     for(k=0;k<count[m];k++) tmp[ntest+k]=list[m][k];
//     ntest+=count[m];
//     for(k=0;k<8;k++) {
//       ii=i+proxies[k][0];
//       jj=j+proxies[k][1];
//       if(ii<0) continue;
//       if(ii==grid->nx) continue;
//       if(jj<0) continue;
//       if(jj==grid->ny) continue;
//       mm=jj*grid->nx+ii;
//       for(l=0;l<count[mm];l++) tmp[ntest+l]=list[mm][l];
//       ntest+=count[mm];
//       }
//     }

  delete[] x;
  delete[] y;
  delete[] z;

  delete[] count;
  for (m=0;m<grid->ny*grid->nx;m++) delete[] list[m];
  delete[] list;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int xyz_readgrid_02 (const char* filename,grid_t *grid, int mode)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------

  Assume full (no missing values) array structured (sequential array's axes
  ordering) description

------------------------------------------------------------------------------*/
{
  FILE *in;
  int nitems,ndata,nlon,nlat;
  int i,j,k,l,m,n,range;
  int ii,jj,mm;
  int **list,*count,*tmp,ntest;
  int downward=0,cartesian=0,complete=0,option,status;
  char line[1024];
  double *correlation,*x,*y,*z,zz,mask=-9999.0,xx,yy,mean_x,mean_y;
  irange_t s;
  int proxies[9][2]={1,0,1,1,0,1,-1,1,-1,0,-1,-1,0,-1,1,-1,0,0};

  status=0;

  in=fopen(filename,"rb");
  if (in == NULL) {
    return(-1);
    }
  status=xyz_skipheader(in);

  fgets(line, 1024, in);

  n=-1;
  while (!feof(in)) {
    fgets(line, 1024, in);
    n++;
    }
  ndata=n;

  x=new double[ndata];
  y=new double[ndata];
  z=new double[ndata];

  rewind (in);
  status=xyz_skipheader(in);
  fgets(line, 1024, in);

  mean_x=0.0;
  mean_y=0.0;

  switch (mode) {
    case 0:
      for(n=0;n<ndata;n++) {
        fgets(line, 1024, in);
        nitems=sscanf(line,"%lf %lf %lf",&x[n],&y[n],&z[n]);
//        mean_x+=x[n];
//        mean_y+=y[n];
        }
      break;

    case 1:
      for(n=0;n<ndata;n++) {
        fgets(line, 1024, in);
        nitems=sscanf(line,"%lf %lf %lf",&y[n],&x[n],&z[n]);
        }
      break;
     }

  fclose(in);

  range=(int) 10*sqrt((double)ndata);

  nlon=sqrt((double) ndata)/3;
  nlat=sqrt((double) ndata)/3;

  poc_minmax(x,ndata,mask,&grid->xmin,&grid->xmax);
  poc_minmax(y,ndata,mask,&grid->ymin,&grid->ymax);

  grid->nx=nlon;
  grid->ny=nlat;
  grid->nz=1;

  grid->dx=(grid->xmax-grid->xmin) / ((double) grid->nx-1.);
  grid->dy=(grid->ymax-grid->ymin) / ((double) grid->ny-1.);

  grid->xmin-=2*grid->dx;
  grid->xmax+=2*grid->dx;
  grid->ymin-=2*grid->dy;
  grid->ymax+=2*grid->dy;
  grid->nx+=4;
  grid->ny+=4;

  grid->x=new double[grid->Hsize()];
  grid->y=new double[grid->Hsize()];
  grid->z=new double[grid->Hsize()];

  for (j=0;j<grid->ny;j++) {
    for (i=0;i<grid->nx;i++) {
      m=j*grid->nx+i;
      grid->x[m]=grid->xmin+(double) i*grid->dx;
      grid->y[m]=grid->ymin+(double) j*grid->dy;
      grid->z[m]=0;
      }
    }
  grid->modeH=0;

/*------------------------------------------------------------------------------
  count number of points in grid cell*/
  count=new int[grid->Hsize()];
  for (j=0;j<grid->ny;j++) {
    for (i=0;i<grid->nx;i++) {
      m=j*grid->nx+i;
      count[m]=0;
      }
    }

  for(n=0;n<ndata;n++) {
    status=map_index(*grid,x[n],y[n],&i,&j);
    m=j*grid->nx+i;
//    grid->z[m]++;
    grid->z[m]+=z[n];
    count[m]++;
    }

  for (j=0;j<grid->ny;j++) {
    for (i=0;i<grid->nx;i++) {
      m=j*grid->nx+i;
      if(count[m]!=0) {
        grid->z[m]/=count[m];
        }
      else {
        grid->z[m]/9999.9;
        }
      }
    }

/*------------------------------------------------------------------------------
  create list of points in grid cell*/
  list=new int*[grid->Hsize()];
  for (j=0;j<grid->ny;j++) {
    for (i=0;i<grid->nx;i++) {
      m=j*grid->nx+i;
      list[m]=new int[count[m]];
      count[m]=0;
      }
    }

  for(n=0;n<ndata;n++) {
    status=map_index(*grid,x[n],y[n],&i,&j);
    m=j*grid->nx+i;
    list[m][count[m]]=n;
    count[m]++;
    }

  s=poc_minmax(count,grid->Hsize(),-1);
  tmp=new int[9*s.max];

/*------------------------------------------------------------------------------
  detect direct neighbours*/
//   for(n=0;n<ndata;n++) {
//     status=map_index(*grid,x[n],y[n],&i,&j);
//     m=j*grid->nx+i;
//     ntest=0;
//     for(k=0;k<count[m];k++) tmp[ntest+k]=list[m][k];
//     ntest+=count[m];
//     for(k=0;k<8;k++) {
//       ii=i+proxies[k][0];
//       jj=j+proxies[k][1];
//       if(ii<0) continue;
//       if(ii==grid->nx) continue;
//       if(jj<0) continue;
//       if(jj==grid->ny) continue;
//       mm=jj*grid->nx+ii;
//       for(l=0;l<count[mm];l++) tmp[ntest+l]=list[mm][l];
//       ntest+=count[mm];
//       }
//     }

  delete[] x;
  delete[] y;
  delete[] z;

  delete[] count;
  for (m=0;m<grid->ny*grid->nx;m++) delete[] list[m];
  delete[] list;
  
  return(0);

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int xyz_loadmap (const char* filename,grid_t *grid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  FILE *in;
  int downward=0,cartesian=0,mode=0,complete=0,masked=0,option,status;
  char line[1024];

  status=0;

  in=fopen(filename,"rb");
  if (in == NULL) {
    return(-1);
    }
  status=xyz_skipheader(in);
  fgets(line, 1024, in);
  fclose(in);

/*------------------------------------------------------------------------------
  decode header informations*/

  mode     =(strstr(line, "YXZ")!=0);
  complete =(strstr(line, "COMPLETE")!=0);
  cartesian=(strstr(line, "CARTESIAN")!=0);
  downward =(strstr(line, "DOWNWARD")!=0);
  masked   =(strstr(line, "MASK=")!=0);

//   if(masked!=0) {
//     fscanf("%lf",mask);
//     }

  switch (complete) {
    case 0:
      status=xyz_readgrid_02 (filename,grid,mode);
      break;

    case 1:
      status=xyz_readgrid_01_from_v1 (filename,grid,mode);
      break;

    default:
       break;
    }

  return(status);
  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int xyz_read02_r1 (const char* filename,grid_t grid, float *buffer, float mask, int mode)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------

  Assume full (no missing values) array structured (sequential array's axes
  ordering) description

------------------------------------------------------------------------------*/
{
  FILE *in;
  int nitems,ndata,nlon,nlat;
  int i,j,k,l,m,n,range;
  int ii,jj,mm;
  int **list,*count,*tmp,ntest;
  int downward=0,cartesian=0,complete=0,option,status;
  char line[1024];
  double *correlation,*x,*y,*z,zz,xx,yy,mean_x,mean_y;
  range_t<double> r;
  irange_t s;
  int proxies[9][2]={1,0,1,1,0,1,-1,1,-1,0,-1,-1,0,-1,1,-1,0,0};

  status=0;

  in=fopen(filename,"rb");
  if (in == NULL) {
    return(-1);
    }
  status=xyz_skipheader(in);

  fgets(line, 1024, in);

  n=-1;
  while (!feof(in)) {
    fgets(line, 1024, in);
    n++;
    }
  ndata=n;

  x=new double[ndata];
  y=new double[ndata];
  z=new double[ndata];

  rewind (in);

  status=xyz_skipheader(in);

  fgets(line, 1024, in);

  mean_x=0.0;
  mean_y=0.0;

  switch (mode) {
    case 0:
      for(n=0;n<ndata;n++) {
        fgets(line, 1024, in);
        nitems=sscanf(line,"%lf %lf %lf",&x[n],&y[n],&z[n]);
//        mean_x+=x[n];
//        mean_y+=y[n];
        }
      break;

    case 1:
      for(n=0;n<ndata;n++) {
        fgets(line, 1024, in);
        nitems=sscanf(line,"%lf %lf %lf",&y[n],&x[n],&z[n]);
        }
      break;
     }

  fclose(in);

/*------------------------------------------------------------------------------
  count number of points in grid cell*/
  count=new int[grid.nx*grid.ny];
  for (j=0;j<grid.ny;j++) {
    for (i=0;i<grid.nx;i++) {
      m=j*grid.nx+i;
      count[m]=0;
      }
    }

  for(n=0;n<ndata;n++) {
    status=map_index(grid,x[n],y[n],&i,&j);
    m=j*grid.nx+i;
//    grid.z[m]++;
    buffer[m]+=z[n];
    count[m]++;
    }

  for (j=0;j<grid.ny;j++) {
    for (i=0;i<grid.nx;i++) {
      m=j*grid.nx+i;
      if(count[m]!=0) {
        buffer[m]/=count[m];
        }
      else {
        buffer[m]=mask;
        }
      }
    }

  delete[] x;
  delete[] y;
  delete[] z;

  delete[] count;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int xyz_loadr1_from_v1 (const char *filename, grid_t grid, float *buf, float *mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  {
  FILE *in;
  int downward=0,cartesian=0,mode=0,complete=0,option,status;
  char *masked=0;
  char line[1024];

  *mask=-99999.0;

  status=0;

  in=fopen(filename,"rb");
  if (in == NULL) {
    return(-1);
    }
  status=xyz_skipheader(in);
  fgets(line, 1024, in);
  fclose(in);

  masked=strstr(line, "MASK=");

  if(masked!=0) {
    masked+=strlen("MASK=");
    sscanf(masked, "%f",mask);
    }

/*------------------------------------------------------------------------------
  decode header informations*/
  mode=(strstr(line, "YXZ")!=0);

  status=xyz_read02_r1 (filename,grid,buf,*mask,mode);

  return(status);
  }

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int xyz_reduce(double * &x, double * &y, double * &z, char *keep, int  &ndata)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int n, nn, count;
  double *xx,*yy,*zz;
  
  count=0;
  for(nn=0;nn<ndata;nn++) {
    if(keep[nn]==1) count++;
    }
    
  xx=new double[count];
  yy=new double[count];
  zz=new double[count];
  
  n=0;
  for(nn=0;nn<ndata;nn++) {
    if(keep[nn]==1) {
      xx[n]=x[nn];
      yy[n]=y[nn];
      zz[n]=z[nn];
      n++;
      }
    }

  delete[] x;
  delete[] y;
  delete[] z;
    
  x=xx;
  y=yy;
  z=zz;
  
  ndata=count;
  
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int xyz_FrameSelection ( frame_t & frame, double* &x, double* &y, double* &z, int & ndata)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int m,n,status;
  int verbose=0;
  bool debug=false;
  
  char *keep=new char[ndata];
  for (m=0;m<ndata;m++) {
    keep[m]=1;
    }
  if(frame.initialised()==1) {
    for(n=0;n<ndata;n++) {
      if(frame.inside(x[n],y[n])!=1) keep[n]=0;
      }
    }
  status=xyz_reduce(x, y, z, keep, ndata);
  
  delete[] keep;
  
  return(0);
}

  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int xyz_PolygonSelection (vector<plg_t> & p, double* &x, double* &y, double* &z, int & ndata, frame_t & frame)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int m,n,status;
  int verbose=0;
  bool debug=false;
  
  frame=plg_spherical_minmax(p);
  char *keep=new char[ndata];
  for (m=0;m<ndata;m++) {
    keep[m]=1;
    }
  if(frame.initialised()==1) {
    for(n=0;n<ndata;n++) {
      if(frame.inside(x[n],y[n])!=1) keep[n]=0;
      }
    }
  status=xyz_reduce(x, y, z, keep, ndata);
  
  if(ndata==0) {
    printf("no data available inside the polygons frame\n");
    plg_destroy_entries(p);
    p.clear();
    }
  
  char *position=plg_TestInterior(x, y, ndata, p, verbose, debug);
  
  int count=0;
  for(n=0;n<ndata;n++) {
    switch (position[n]) {
      case PLG_POINT_BOUNDARY:
      case PLG_POINT_INTERIOR:
        count++;
        break;
      }
    }
  
  double *xx=new double[count];
  double *yy=new double[count];
  double *zz=new double[count];

  count=0;
  for(n=0;n<ndata;n++) {
    switch (position[n]) {
      case PLG_POINT_BOUNDARY:
      case PLG_POINT_INTERIOR:
        xx[count]=x[n];
        yy[count]=y[n];
        zz[count]=z[n];
        count++;
        break;
      }
    }
  
  delete[] x;
  delete[] y;
  delete[] z;
  
  delete[] keep;
  
  ndata=count;
  x=xx;
  y=yy;
  z=zz;
  
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int xyz_PolygonSelection (const char *polygons, double* &x, double* &y, double* &z, int & ndata, frame_t & frame)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int m,n,status;
  vector<plg_t> p;
  bool debug=false;
  
  status=plg_load(polygons, PLG_FORMAT_SCAN, p);
  if(status!=0) TRAP_ERR_EXIT(status, "file reading error : %s, exit with bad status\n",polygons);
  
  status=xyz_PolygonSelection (p, x, y, z, ndata, frame);
  
  plg_destroy_entries(p);
  p.clear();
  
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int xyz_PolygonSelection (plg_t *polygons, int npolygons, double* &x, double* &y, double* &z, int & ndata, frame_t & frame)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int m,n,status;
  vector<plg_t> p;
  bool debug=false;
  
  p=plg_array2vector(polygons,npolygons);
  
  status=xyz_PolygonSelection (p, x, y, z, ndata, frame);
    
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int xyz_CheckDuplicate_spherical(double *lon, double *lat, double *z, int ndata, char *keep, double resolution)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  resolution in meters

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
{
  int count=0;
  double x,y;
  bool debug=false;
  double dxmax, dymax;
  double dmax=resolution,dtr=M_PI/180.0;
  
  dymax=dmax/110000.;
    
  printf("#################################################################\n");
  printf("decimating resolution %lf m\n",resolution);
  for (int m=0;m<ndata;m++) {
    if(keep[m]==0) continue;
    x=lon[m];
    y=lat[m];
    dxmax=dymax*cos(y*dtr);
#pragma omp parallel for
    for (int n=m+1;n<ndata;n++) {
      if(keep[n]==0) continue;
      double xx=lon[n];
      if(fabs(xx-x)>dxmax) continue;
      double yy=lat[n];
      if(fabs(yy-y)>dymax) continue;
      double d=geo_haversin(x,y,xx,yy);
      if(d>dmax) continue;
      if(debug) printf("%d %d %lf %lf %lf %lf %lf\n",m,n,x,y,d,z[m],z[n]);
      keep[n]=0;
      }
    count++;
    }
  return(count);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int xyz_CheckDuplicate_cartesian(double *lon, double *lat, double *z, int ndata, char *keep, double resolution)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  resolution in meters

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int count=0;
  double x,y;
  bool debug=false;
  double dmax,d2max;
  
  dmax=resolution;
  d2max=dmax*dmax;
      
  printf("#################################################################\n");
  printf("decimating resolution %lf m\n",resolution);
  
  for (int m=0;m<ndata;m++) {
    if(keep[m]==0) continue;
    x=lon[m];
    y=lat[m];
#pragma omp parallel for
    for (int n=m+1;n<ndata;n++) {
      if(keep[n]==0) continue;
      double xx=lon[n];
      const double dx=xx-x;
      if(fabs(dx)>dmax) continue;
      double yy=lat[n];
      const double dy=yy-y;
      if(fabs(dy)>dmax) continue;
      double d=dx*dx+dy*dy;
      if(d>d2max) continue;
      if(debug) printf("%d %d %lf %lf %lf %lf %lf\n",m,n,x,y,d,z[m],z[n]);
      keep[n]=0;
      }
    count++;
    }
  return(count);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int xyz_CheckDuplicate(double *t, double *p, double *z, int ndata, char *keep, double resolution, projPJ projection)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status, count=0;
  range_t<double> lon_range=range_t<double>(t,ndata);
  range_t<double> lat_range=range_t<double>(p,ndata);

  frame_t frame(lon_range,lat_range);
  
  projPJ tmp=0;
  char *params=new char[1024];
//   double *x, *y;

  if(projection==0) {  
    projection=assign_StereoOblique(frame.y_center(), frame.x_center(), params);
    printf("%s : cartesian projection, parameters=%s\n",__func__,params);
//     x=new double[ndata];
//     y=new double[ndata];
//     status=geo_to_projection(params, t, p, x, y, ndata);
    status=geo_to_projection(params, t, p, ndata);
    printf("%s : cartesian projection done\n",__func__);
    }
  
  count=xyz_CheckDuplicate_cartesian(t, p, z, ndata, keep, resolution);
  
  if(tmp!=0) {  
    status=projection_to_geo(params, t, p, ndata);
    }  
   
  
  return(count);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int decimate_spherical(double *lon, double *lat, double *z,char *keep, double resolution, size_t *indices, size_t size, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  
  decimate data subset, substitude original data with average data at resolution
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int k,m,mm,count=0;
  char *eligible;
  double x,y;
  
  size_t ndata=size;
  
  eligible=new char[ndata];
  for (m=0;m<ndata;m++) {
    eligible[m]=1;
    }
      
  for (mm=0;mm<ndata;mm++) {
    vector <size_t> clusters;
    double dmin=1.e+10,dmax=resolution/2.0;
    int    size, target;
    double cx,cy,cz;
    m=indices[mm];
    if(keep[m]     ==1) continue;
    if(eligible[mm]==0) continue;
    x=lon[m];
    y=lat[m];
/*------------------------------------------------------------------------------
    find other data located inside circle centered on m and with radius dmax */
#pragma omp parallel for
    for (int nn=0;nn<ndata;nn++) {
      if(eligible[nn]==0) continue;
      int n=indices[nn];
      if(keep[n]    ==1) continue;
      double xx=lon[n];
      double yy=lat[n];
      double d=geo_haversin_km(x,y,xx,yy);
      if(d>dmax) continue;
#pragma omp critical(clusterIsCritical)
      {
      clusters.push_back((size_t) n);
      }
      eligible[nn]=0;
      }
    if(clusters.size()==1) {
      keep [m]=1;
      count++;
      continue;
      }
      
/*------------------------------------------------------------------------------
    compute mean z value and barycentric position                             */
    cx=cy=cz=0.0;
    for(k=0;k<clusters.size();k++) {
      int n=clusters[k];
      double xx=lon[n];
      double yy=lat[n];
      xx=degree_recale(xx,x);
      cx+=xx;
      cy+=yy;
      cz+=z[n];
      }
    size=clusters.size();
    cx/= clusters.size();
    cy/= clusters.size();
    cz/= clusters.size();
/*------------------------------------------------------------------------------
    then get closer true node position and affect z and position              */
    for(k=0;k<clusters.size();k++) {
      int n=clusters[k];
      double xx=lon[n];
      double yy=lat[n];
      double d=geo_haversin_km(cx,cy,xx,yy);
      if(d<dmin) {
        dmin=d;
        target=n;
        }
      }
    keep [target]=1;
    lon  [target]=cx;
    lat  [target]=cy;
    z    [target]=cz;
    count++;
    }
    
  return(count);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

//   int decimate_cartesian(double *lon, double *lat, double *x, double *y, double *z, char *keep, double resolution, size_t *indices, size_t size)
  int decimate_cartesian(double *x, double *y, double *z, char *keep, double resolution, size_t *indices, size_t size, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  
  decimate data subset, substitude original data with average data at resolution
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int k,m,mm,count=0;
  char *eligible;
  double xm,ym;
  double dmax=resolution/2.0;
  double dmax2=dmax*dmax;
  
  size_t ndata=size;
  
  eligible=new char[ndata];
  for (m=0;m<ndata;m++) {
    eligible[m]=1;
    }

  int nprocs=omp_get_num_procs();
  vector<size_t> *clusters_array=new vector<size_t>[nprocs];
      
  for (mm=0;mm<ndata;mm++) {
    vector <size_t> clusters;
    double dmin=1.e+10;
    int    size, target;
    double cx,cy,cz;
    m=indices[mm];
    if(keep[m]     ==1) continue;
    if(eligible[mm]==0) continue;
    xm=x[m];
    ym=y[m];
/*------------------------------------------------------------------------------
    find other data located inside circle centered on m and with radius dmax */
#pragma omp parallel for
    for (int nn=0;nn<ndata;nn++) {
      if(eligible[nn]==0) continue;
      size_t n=indices[nn];
      if(keep[n]    ==1) continue;
      double dx=x[n]-xm;
      double dy=y[n]-ym;
      double d=dx*dx+dy*dy;
      if(d>dmax2) continue;
      int proc=omp_get_thread_num();
      
//       d=sqrt(d);
//       double d2=1000.*geo_haversin_km(lon[m],lat[m],lon[n],lat[n]);
      
// #pragma omp critical(clusterIsCritical)
//       {
//       clusters.push_back(n);
      clusters_array[proc].push_back(n);
//       }
      eligible[nn]=0;
      }
    
    for(int p=0;p<nprocs;p++) {
      for(int k=0;k<clusters_array[p].size();k++) clusters.push_back(clusters_array[p][k]);
      clusters_array[p].clear();
      }
    if(clusters.size()==1) {
      keep [m]=1;
      count++;
      continue;
      }
      
/*------------------------------------------------------------------------------
    compute mean z value and barycentric position                             */
    cx=cy=cz=0.0;
    for(k=0;k<clusters.size();k++) {
      int n=clusters[k];
//       double dx=x[n]-xm;
//       double dy=y[n]-ym;
      cx+=x[n];
      cy+=y[n];
      cz+=z[n];
      }
    size=clusters.size();
    cx/= clusters.size();
    cy/= clusters.size();
    cz/= clusters.size();
/*------------------------------------------------------------------------------
    then get closer true node position and affect z and position              */
    for(k=0;k<clusters.size();k++) {
      int n=clusters[k];
      double dx=x[n]-xm;
      double dy=y[n]-ym;
      double d=dx*dx+dy*dy;
      if(d<dmin) {
        dmin=d;
        target=n;
        }
      }
    keep [target]=1;
    x    [target]=cx;
    y    [target]=cy;
    z    [target]=cz;
    count++;
    }
    
  delete[] clusters_array;
  delete[] eligible;
  
  return(count);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int xyz_decimate(frame_t frame, double *t, double *p, double *z,char *keep, int ndata, double resolution, double factor, bool cartesian, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  size_t examined, count;
  int *cardinal;
  size_t *cells;
  int step=ndata/100;
  int milestone=step;
  int max_clustersize=400,  max_cardinal=0;
  double ratio;
  
  double resolution_in_degrees=resolution/110;
  
  int nx=MAX(10,frame.x_size()/(factor*resolution_in_degrees));
  int ny=MAX(10,frame.y_size()/(factor*resolution_in_degrees));
  
  if(verbose==1) printf("%s : prior set nx=%d ny=%d\n",__func__,nx,ny);
  ratio=(double) nx/(double) ny;
  
  int nxny=min(nx*ny,max_clustersize*max_clustersize);
  
  nx=sqrt(nxny*ratio);
  ny=nx/ratio;
  
  nx=MAX(10,nx);
  ny=MAX(10,ny);
  
  if(verbose==1) printf("%s : optimised set nx=%d ny=%d\n",__func__,nx,ny);
 
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  decimation working (clustering) grid 
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  grid_t grid=map_Rgrid(frame, (size_t) nx, (size_t) ny, 0);
  
  if(verbose==1) {
    printf("#################################################################\n");
    printf("decimation step, resolution=%lf km; decimation cluster:\n", resolution);
    map_printgrid(grid);
    }
  
  cardinal=new int[grid.Hsize()];
  for(int m=0;m<grid.Hsize();m++) cardinal[m]=0;
  
  cells=new size_t[ndata];

/*------------------------------------------------------------------------------
  count and register data per cluster grid cells */
  for (int m=0;m<ndata;m++) {
    int i,j;
    status=map_index(grid,t[m],p[m],&i,&j);
    if(status==0) {
      size_t n=i+grid.nx*j;
      cells[m]=n;
      cardinal[n]++;
      }
    else {
      cells[m]=-1;
      }
    }
    
  for(size_t n=0;n<grid.Hsize();n++) max_cardinal=max(max_cardinal, cardinal[n]);
  if(verbose==1) printf("max_cardinal=%d\n", max_cardinal);
  
  statistic_t s=get_statistics(cardinal, -1,grid.Hsize(),verbose);
  
  size_t **indices=new size_t*[grid.Hsize()];
  for(size_t n=0;n<grid.Hsize();n++) {
    if(cardinal[n]!=0) indices[n]=new size_t[cardinal[n]];
    cardinal[n]=0;
    }
  for (int m=0;m<ndata;m++) {
    size_t n=cells[m];
    if(n==-1) continue;
    indices[n][cardinal[n]]=m;
    cardinal[n]++;
    }
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  scan data by grid cells 
  
  optimizing questions :
  
    * use array instead of vector for index listing ?                       YES
    * recreate cells propriatery set of data instead of index list ?        ...
    * in that case, deal witl local cell projection ?                       ...
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  
  projPJ projection;
  char *params=new char[1024];
//   double *x, *y;

  if(cartesian) {  
    projection=assign_StereoOblique(frame.y_center(), frame.x_center(), params);
    if(verbose==1) printf("%s : cartesian projection, parameters=%s\n",__func__,params);
//     x=new double[ndata];
//     y=new double[ndata];
//     status=geo_to_projection(params, t, p, x, y, ndata);
    status=geo_to_projection(params, t, p, ndata);
    if(verbose==1) printf("%s : cartesian projection done\n",__func__);
    }
  
  int nprocs=initialize_OPENMP(-1, verbose);
//   size_t *indices=new size_t[ndata];
  count=0;
  examined=0;
  for(int j=0;j<grid.ny;j++) {
    for (int i=0;i<grid.nx;i++) {
      size_t size=0;
      size_t n=i+grid.nx*j;
      if(cardinal[n]==0) continue;
// /*------------------------------------------------------------------------------
//       identify data in current cluster cell n                                 */
// #pragma omp parallel for if(nprocs>1)
//       for (int m=0;m<ndata;m++) {
//         if(cells[m]==n) {
// #pragma omp critical(completedIsCritical)
//           {
//           indices[size]=m;
//           size++;
//           }
//           }
//         }
//       if(size!=0) {
//         printf("i=%d j=%d size=%d\n",i,j,cardinal[n]);
        int done;
        if(cartesian) 
//         done=decimate_cartesian(t, p, x, y, z, keep, 1000.*resolution, indices[n], cardinal[n]);
          done=decimate_cartesian(t, p, z, keep, 1000.*resolution, indices[n], cardinal[n], verbose);
        else
          done=decimate_spherical(t, p, z, keep, resolution, indices[n], cardinal[n], verbose);
        count+=done;
        examined+=cardinal[n];
        if(examined>milestone) {
          if(verbose==1) printf("examined %8d, kept %8d (over %8d)\n",examined,count,ndata);
          milestone+=step;
          }
//         }
      }
    }
    
  if(verbose==1) printf("examined %d, kept %d (over %d)\n",examined,count,ndata);

  if(cartesian) {  
    status=projection_to_geo(params, t, p, ndata);
    }  
  
  delete[] cells;
  delete[] cardinal;
  delete[] indices;
  
  return(count);
}

