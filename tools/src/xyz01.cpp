
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
// #define MAIN_SOURCE
#include <stdio.h>
#include <string.h>

#include <config.h>

#include "tools-structures.h"

#include "rutin.h"
#include "geo.h"
#include "polygones.h"
#include "functions.h"
#include "grd.h"
#include "map.h"
#include "statistic.h"
#include "xyz.h"

#define XYZ 0
#define YXZ 1

#define ROW    0
#define COLUMN 1

// #if HAVE_SHAPEFIL_H || HAVE_LIBSHP
// #include <shapefil.h>
// #endif

#if HAVE_LIBSHP
#if HAVE_SHAPEFIL_H
#include <shapefil.h>
#elif HAVE_LIBSHP_SHAPEFIL_H
#include <libshp/shapefil.h>
#endif
#endif

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void xyz_read_shp(const char *filename, double *x, double *y, double *z, int np)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------*/
/// Read shapes from an SHP file
/**
\param *filename
\param *polygones pointer to the polygones
\param np number of shapes
\returns STUFF
*/
// http://shapelib.maptools.org/shp_api.html
/*----------------------------------------------------------------------*/
{
#if HAVE_SHAPEFIL_H || HAVE_LIBSHP
  SHPHandle shpfile;
  int k,l,oi,p;//object and polygone indexes
  SHPObject *shp;
  
  shpfile=SHPOpen(filename,"rb");
  
  p=0;
  for(oi=0;oi<np;oi++){
    shp=SHPReadObject(shpfile,oi);
    if(shp->nSHPType!=SHPT_POINTZ){
      continue;
      }
    for(l=0;l<shp->nVertices;l++) {
      x[p]=shp->padfX[l];
      y[p]=shp->padfY[l];
      z[p]=shp->padfZ[l];
      p++;
      }
    SHPDestroyObject(shp);
    }
  
  SHPClose(shpfile);
#else
  return;
#endif
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int xyz_inquire_shp(const char *filename, int *np, int *nblocks, int *shapetype)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------*/
/// Get number of shapes and shape types of an SHP file
/**
\param *filename
\param *np pointer to the total number of points.
\param *nblocks pointer to the number of blocks. May be NULL.
\param *shapetype pointer to the type of shapes. May be NULL.
\returns the type of the shapes in the file. See http://shapelib.maptools.org/shp_api.html for more info
*/
/*----------------------------------------------------------------------*/
{
#if HAVE_SHAPEFIL_H || HAVE_LIBSHP
  SHPHandle shpfile;
  int i;
  SHPObject *shp;
  
  shpfile=SHPOpen(filename,"rb");
  __NC_CHKERR_LINE_FILE__(errno,"SHPOpen(\"%s\",,) error",filename);
  
  SHPGetInfo(shpfile,nblocks,shapetype,NULL,NULL);
  for(i=0;i<*nblocks;i++){
    shp=SHPReadObject(shpfile,i);
    switch(shp->nSHPType){
      case SHPT_POINTZ:
        *np+=shp->nVertices;
        break;
      case SHPT_POLYGONZ:
      case SHPT_ARCZ:
        *np+=shp->nParts;
        break;
      default:
        (*nblocks)--;
        break;
      }
//     if(shp->nSHPType!=SHPT_POINTZ){
//       *nblocks--;
//       }
//     else {
//       *np+=shp->nVertices;
//       }
    SHPDestroyObject(shp);
    }
  
  SHPClose(shpfile);
  
  return 0;
#else
  return -1;
#endif
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int xyz_load_shp(const char *filename, char *proj4_options, double * &x, double * &y, double * &z, double *mask, int & ndata)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status, npoints=0, nblocks=0, shapetype;
    
  status=ENOEXEC;

#if HAVE_SHAPEFIL_H || HAVE_LIBSHP
/**----------------------------------------------------------------------------
  one block may contain more than one segement */
  xyz_inquire_shp (filename, &npoints, &nblocks, &shapetype);
  if(shapetype!=SHPT_POINTZ){
    status=NC_EBADTYPE;
    __ERR_BASE_LINE__("Shape %d not supported\n",shapetype);
    return status;
    }

  x=new double[npoints];
  y=new double[npoints];
  z=new double[npoints];

  xyz_read_shp(filename, x, y, z, nblocks);
  
  ndata=npoints;
  status=xyz_save ("echo.xyz", x, y, z, *mask, ndata);

  return status;
#else
  printf("please install libshp-devel....\n");
  return(status);
#endif
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int xyz_countdata(const char* filename)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------


  Assume full (no missing values) array structured (sequential array's axes
  ordering) description

------------------------------------------------------------------------------*/
{
  FILE *in;
  int nitems,ndata,status;
  int i,j,k,l,m,n,range;
  char line[1024];
  char *s;

  status=0;

  in=fopen(filename,"rb");
  if (in == NULL) {
    return(-1);
    }
  status=xyz_skipheader(in);

  s=fgets(line, 1024, in);
  
  if(s==0) return(0);
  
  if(line[0]=='#'|| (strstr(line, "Y")!=0)) n=0;
  else             n=1;
  
  while (!feof(in)) {
    s=fgets(line, 1024, in);
    if(s==0) break;
    n++;
    }
    
  ndata=n;
  
  return(ndata);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int xyz_read(const char* filename,double *x,double *y,double *z, double mask, int & ndata, int mode)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  FILE *in;
  int nitems;
  int i,j,k,l,m,n,range,dum,nlines;
  int status;
  char *line=new char[1024];
  char *s;
  char time[256];
  double e;

  status=0;

  in=fopen(filename,"rb");
  if (in == NULL) {
    return(-1);
    }
  status=xyz_skipheader(in);

//  fgets(line, 1024, in);

  nlines=0;

  switch (mode) {
    case -1:
      for(n=0;n<ndata;n++) {
        nitems=fscanf(in,"%lf %lf %lf",&x[n],&y[n],&z[n]);
        if(nitems!=3) goto error;
        }
      break;

    case 0:
      for(n=0;n<ndata;n++) {
        s=fgets(line, 1024, in);
        nlines++;
        if(s==0) break;
        if(s[0]=='#') {
          n--;
          continue;
          }
        for(k=0;k<strlen(line);k++) if(line[k]==',') line[k]=' ';
        nitems=sscanf(line,"%lf %lf %lf",&x[n],&y[n],&z[n]);
        if(nitems!=3) goto error;
        if(z[n]==mask) {
          n--;
          ndata--;
          }
        }
      break;

    case 1:
      for(n=0;n<ndata;n++) {
        s=fgets(line, 1024, in);
        nlines++;
        if(s==0) break;
        if(s[0]=='#') {
          n--;
          continue;
          }  
        if(strstr(line, "Y")!=0) {
          n--;
          continue;
          }  
        for(k=0;k<strlen(line);k++) if(line[k]==',') line[k]=' ';
        nitems=sscanf(line,"%lf %lf %lf",&y[n],&x[n],&z[n]);
        if(nitems!=3) goto error;
        }
      break;

    case 2:
      for(n=0;n<ndata;n++) {
        s=fgets(line, 1024, in);
        nlines++;
        if(s==0) break;
        if(s[0]=='#') {
          n--;
          continue;
          }
        for(k=0;k<strlen(line);k++) if(line[k]==',')  line[k]=' ';
        for(k=0;k<strlen(line);k++) if(line[k]=='\t') line[k]=' ';
        nitems=sscanf(line,"%s %s %lf %lf %lf %lf",time,time,&x[n],&y[n],&z[n],&e);
        if(nitems!=6) goto error;
        }
      break;

    case 3:
      for(n=0;n<ndata;n++) {
        s=fgets(line, 1024, in);
        nlines++;
        if(s==0) break;
        if(s[0]=='#') {
          n--;
          continue;
          }
        for(k=0;k<strlen(line);k++) if(line[k]==',') line[k]=' ';
        nitems=sscanf(line,"%d %lf %lf %lf",&dum,&y[n],&x[n],&z[n]);
        if(nitems!=4) goto error;
        }
      break;
      
    case 4:
      for(n=0;n<ndata;n++) {
        s=fgets(line, 1024, in);
        nlines++;
        if(s==0) break;
        if(s[0]=='#') {
          n--;
          continue;
          }
        for(k=0;k<strlen(line);k++) if(line[k]==',') line[k]=' ';
        nitems=sscanf(line,"%lf %lf %lf %lf",&x[n],&y[n],&z[n],&e);
        if(nitems!=4) goto error;
        }
      break;
      
    case 5:
      for(n=0;n<ndata;n++) {
        s=fgets(line, 1024, in);
        nlines++;
        if(s==0) break;
        if(s[0]=='#') {
          n--;
          continue;
          }
        for(k=0;k<strlen(line);k++) if(line[k]==',') line[k]=' ';
        nitems=sscanf(line,"%lf %lf",&x[n],&y[n]);
        z[n]=-1;
        if(nitems!=2) goto error;
        }
      break;
      
     }

  fclose(in);
  delete[] line;
  return(0);

error:
  printf("read error at line=%d, exit with status -1\n",nlines);
  fclose(in);
  delete[] line;
  return(-1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int xyz_save_template(const char *filename, double *x, double *y, T *z, char *keep, int flag, int ndata, int header)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  FILE *out;
  size_t n;
  size_t count;
  
  out=fopen(filename,"w");
  
  if(header==1) fprintf(out,"#XYZ\n");
  
  if(flag==0) {
    count=occurence((char) 0,keep,ndata);
    for(n=0;n<ndata;n++) {
      if(keep[n]==0) fprintf(out,"%f %f %g\n",x[n],y[n],z[n]);
      }
    }
  else {
    count=occurence((char) 1,keep,ndata);
    for(n=0;n<ndata;n++) {
      if(keep[n]==1) fprintf(out,"%f %f %g\n",x[n],y[n],z[n]);
      }
    }
  
  fclose(out);
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int xyz_save(const char *filename, double *x, double *y, double *z, char *keep, int flag, int ndata)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=xyz_save_template(filename, x, y, z, keep, flag, ndata, 1);
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int xyz_save_template(const char *filename, double *x, double *y, T *z, T mask, int ndata, int header)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  FILE *out;
  int n;
  
  out=fopen(filename,"w");
  
  if(header==1) fprintf(out,"#XYZ\n");
  for(n=0;n<ndata;n++) {
    fprintf(out,"%f %f %g\n",x[n],y[n],z[n]);
    }
  
  fclose(out);
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int xyz_save(const char *filename, double *x, double *y, double *z, double mask, int ndata)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=xyz_save_template(filename, x, y, z, mask, ndata, 1);
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int xyz_save(const char *filename, double *x, double *y, double *z, double mask, int ndata, int header)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=xyz_save_template(filename, x, y, z, mask, ndata, header);
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int xyz_save(const char *filename, double *x, double *y, short *z, short mask, int ndata)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=xyz_save_template(filename, x, y, z, mask, ndata, 1);
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int xyz_save(const char *filename, double *x, double *y, short *z, short mask, int ndata, int header)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=xyz_save_template(filename, x, y, z, mask, ndata, header);
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int xyz_save(const char *filename, double *x, double *y, float *z, float mask, int ndata)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=xyz_save_template(filename, x, y, z, mask, ndata, 1);
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int xyz_save_template(const char *filename, grid_t grid, T *z, T mask, int header)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  FILE *out;
  int i,j,n;
  int ndata=grid.nx*grid.ny;
  double x,y;
  
  out=fopen(filename,"w");
  
  if(out==0) return(-1);
  
  if(header==1) fprintf(out,"#XYZ\n");

  for(j=0;j<grid.ny;j++) {
    for(i=0;i<grid.nx;i++) {
      n=grid.Hindex( i,  j);
      grid.xy(i,j,x,y);
//      fprintf(out,"%lf %lf %lf\n",x,y, 0.0);
      fprintf(out,"%lf %lf %lf\n",x,y, (double) z[n]);
      }
    }
  
  fclose(out);
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int xyz_save(const char *filename, grid_t grid, double *z, double mask, int header)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=xyz_save_template (filename, grid, z, mask,header);
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int xyz_save(const char *filename, grid_t grid, double *z, double mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=xyz_save_template (filename, grid, z, mask,1);
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int xyz_save(const char *filename, grid_t grid, float *z, float mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=xyz_save_template (filename, grid, z, mask,1);
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int xyz_save(const char *filename, grid_t grid, short *z, short mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=xyz_save_template (filename, grid, z, mask,1);
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int xyz_save(const char *filename, double *x, double *y, double **z, double mask, int ndata, int ncols)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  FILE *out;
  int k,n;
  
  out=fopen(filename,"w");
  
  if(out==0) return(-1);
  
  fprintf(out,"#XYZ\n");
  fprintf(out,"%d %d\n",ndata,ncols);
  for(n=0;n<ndata;n++) {
    fprintf(out,"%lf %lf",x[n],y[n]);
    for(k=0;k<ncols;k++) fprintf(out," %lf",z[k][n]);
    fprintf(out,"\n");
    }
  fclose(out);
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int xyz_loadraw(const char *filename, string header, char *proj4_options, double * &x, double * &y, double * &z, double *mask, int & ndata, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  read XYZ file, and invert projection if projection parameters given
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  FILE *in;
  int downward=0,cartesian=0,mode=0,complete=0,option,status;
  int n,nitems,ntoken;
  char *masked=0;
  char *line,options[16],*s,*tmp;
  double t,p;
  int verbose=1;

  *mask=-99999.0;
  line=new char[1024];
  status=0;

  in=fopen(filename,"r");
  if (in == NULL) {
    fprintf(stderr,"\nxyz_loadraw: opening file %s for reading failed!\n\n",filename);
    return(-1);
    }
  status=xyz_skipheader(in);
  fgets(line, 1024, in);
  fclose(in);

  masked=strstr(line, "MASK=");

  if(masked!=0) {
    masked+=strlen("MASK=");
    sscanf(masked, "%lf",mask);
    }

/*------------------------------------------------------------------------------
  decode header informations*/
  if(header!="") strcpy(line, header.c_str());
  
/*------------------------------------------------------------------------------
  XYZ mode, no checks for mask nor for funny seperator characters*/
  mode=-1;
  
  if     (strstr(line, "TXYZE")!=0) mode=2;
  else if(strstr(line, "NYXZ")!=0)  mode=3;
  else if(strstr(line, "XYZE")!=0)  mode=4;
  else if(strstr(line, "XYZ")!=0)   mode=0;
  else if(strstr(line, "YXZ")!=0)   mode=1;
  else if(strstr(line, "XY")!=0)    mode=5;

  ndata=xyz_countdata (filename);

  x=new double[ndata];
  y=new double[ndata];
  z=new double[ndata];

  
  printf("#################################################################\n");
  printf("load XYZ file : %s, %d data\n",filename,ndata);
  status=xyz_read (filename,x,y,z,*mask,ndata,mode);
  if(status!=0) {
    return(-1);
    delete[] line;
    }
//   printf("load XYZ file : %s, %d data\n",filename,ndata);

  if(proj4_options!=0) {
    status=projection_to_geo(proj4_options, x, y, ndata, verbose);
    }

  if(debug) status=xyz_save ("echo.xyz", x, y, z, *mask, ndata);
  
  delete[] line;
  return(status);
  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int reorder(double *x, int nvalues)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,l,m,n;
  double tmp;

/*------------------------------------------------------------------------------
  re-arrange x array in ascending order*/
  for(k=0;k<nvalues;k++) {
    l=minpos(&(x[k]),nvalues-k)+k;
    if(k!=l) {
      tmp=x[k];
      x[k]=x[l];
      x[l]=tmp;
      }
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int count_columns(char *s)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,l,m,n,count;
  char *token;

  token = strtok(s," ");
  while (token!=0) {
    token = strtok(0," ");
    }

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int xyz_readgrid_00(const char* filename, char *proj4, grid_t *grid, int mode)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  Assumes grid shape is rectangular, and data set is partial (missing values
  not given through a mask flag)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  FILE *in;
  int add,nitems,ndata,nlon,nlat;
  int i,j,k,l,m,n,length;
  int downward=0,cartesian=0,complete=0,option,status;
  char line[1024];
  double *x,*y,*z,zz,mask=-9999.0;
  double *xtics,*ytics,tmp;
  range_t<double> r;

  status=0;

  in=fopen(filename,"rb");
  if (in == NULL) {
    return(-1);
    }

  fgets(line, 1024, in);

  n=-1;
  while (!feof(in)) {
    fgets(line, 1024, in);
    n++;
    }
  ndata=n;

  x=new double[ndata];
  y=new double[ndata];

  rewind (in);

  fgets(line, 1024, in);

  switch (mode) {
    case 0:
      for(n=0;n<ndata;n++) {
        fgets(line,1024,in);
        nitems=sscanf(line,"%lf %lf",&x[n],&y[n]);
        }
      break;

    case 1:
      for(n=0;n<ndata;n++) {
        fgets(line,1024,in);
        nitems=sscanf(line,"%lf %lf",&y[n],&x[n]);
        }
      break;
     }

  fclose(in);
  
  range_t<double> rx=range_t<double>(x,ndata);
  range_t<double> ry=range_t<double>(y,ndata);
  
  rx.print("x range");
  ry.print("y range");

  if(proj4==0) {
    if(x[0]>180.) x[0]-=360.0;
    for(n=0;n<ndata;n++) {
      x[n]=geo_recale(x[n],x[0],180.0);
      }
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  identify all distincts x coordinates

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  xtics=new double[ndata];
  xtics[0]=x[0];
  nlon=1;
  for(n=0;n<ndata-1;n++) {
    m=n+1;
    while(x[n]==x[m]) {
      m++;
      if(m==ndata) break;
      }
    n=m-1;
    add=1;
    for(k=0;k<nlon;k++) {
      if(xtics[k]==x[n]) {
        add=0;
        break;
        }
      }
    if(add==1) {
      nlon++;
      xtics[nlon-1]=x[n];
      }
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  identify all distincts y coordinates

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  ytics=new double[ndata];
  ytics[0]=y[0];
  nlat=1;
  for(n=0;n<ndata-1;n++) {
    m=n+1;
    while(y[n]==y[m]) {
      m++;
      if(m==ndata) break;
      }
    n=m-1;
    add=1;
    for(k=0;k<nlat;k++) {
      if(ytics[k]==y[n]) {
        add=0;
        break;
        }
      }
    if(add==1) {
      nlat++;
      ytics[nlat-1]=y[n];
      }
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  put tics in ascending order

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  reorder(xtics,nlon);
  reorder(ytics,nlat);
  
  double dx=get_timesampling(xtics, NAN, nlon);
  double dy=get_timesampling(ytics, NAN, nlat);

  printf("apparent dx=%f dy=%f\n", dx, dy);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

 set up rectangular grid

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

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

  for (j=0;j<grid->ny;j++) {
    for (i=0;i<grid->nx;i++) {
      m=j*grid->nx+i;
      grid->x[m]=xtics[i];
      grid->y[m]=ytics[j];
      }
    }

  grid->modeH=2;

  delete[] x;
  delete[] y;
  delete[] xtics;
  delete[] ytics;

  return(0);

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int xyz_readgrid_01 (const char* filename,grid_t *grid, int mode, int cartesian)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  Assumes grid shape is rectangular, and data set is complete (may be masked)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  FILE *in;
  int nitems,ndata,nlon,nlat;
  int i,j,k,l,m,n,length,block_length;
  int downward=0,complete=0,arrangement,option,status;
  char line[1024];
  double *x,*y,*z,zz,tmp,mask=-9999.0,dy;
  range_t<double> r;

  status=0;

  in=fopen(filename,"rb");
  if (in == NULL) {
    return(-1);
    }

  fgets(line, 1024, in);

/*------------------------------------------------------------------------------
  rough data count */
  n=-1;
  while (!feof(in)) {
    fgets(line, 1024, in);
    n++;
    }
  ndata=n;

  x=new double[ndata];
  y=new double[ndata];
//  z=new double[ndata];

  rewind (in);

  fgets(line, 1024, in);

  switch (mode) {
    case 0:
      for(n=0;n<ndata;n++) {
/*------------------------------------------------------------------------------
        positions in X Y order */
        fgets(line,1024,in);
        nitems=sscanf(line,"%lf %lf",&x[n],&y[n]);
        }
      break;

    case 1:
/*------------------------------------------------------------------------------
      positions in Y X order */
      for(n=0;n<ndata;n++) {
        fgets(line,1024,in);
        nitems=sscanf(line,"%lf %lf",&y[n],&x[n]);
        }
      break;
     }

  fclose(in);

  if(cartesian==0) {
    if(x[0]>180.) x[0]-=360.0;
//     for(n=0;n<ndata;n++) {
//       x[n]=geo_recale(x[n],x[0]+180.0,180.0);
//       }
    }

  if(x[0]==x[1])
    arrangement=COLUMN;
  else
    arrangement=ROW;

  block_length=0;

/*------------------------------------------------------------------------------

  get array dimensions*/
  switch (arrangement) {
    case COLUMN:
      nlon=1;
      for(n=0;n<ndata-1;n++) {
        m=n+1;
        while(x[n]==x[m]) {
          if(m==ndata-1) break;
          m++;
          }
        nlon++;
        n=m;
        }
      nlon--;
      nlat=ndata/nlon;
      break;

    case ROW:
      nlat=1;
      for(n=0;n<ndata-1;n++) {
        m=n+1;
        while(y[n]==y[m]) {
          if(m==ndata-1) break;
          m++;
          }
        nlat++;
        block_length=m-n;
        n=m;
        }
      nlat--;
      nlon=ndata/nlat;
      break;
     }

  printf("nx=%d ny=%d\n",nlon,nlat);
  poc_minmax(x,ndata,mask,&grid->xmin,&grid->xmax);
  poc_minmax(y,ndata,mask,&grid->ymin,&grid->ymax);

  grid->nx=nlon;
  grid->ny=nlat;
  grid->nz=1;

  grid->dx=(grid->xmax-grid->xmin) / ((double) grid->nx-1.);
  grid->dy=(grid->ymax-grid->ymin) / ((double) grid->ny-1.);

  grid->x=new double[ndata];
  grid->y=new double[ndata];
//  grid->z=new double[ndata];

  for (n=0;n<ndata;n++) {
    dy=y[n+1]-y[n];
    if(dy>0.) {
/*------------------------------------------------------------------------------
      assumes latitude upward (i.e. file starts with lowest latitudes)*/
      downward=0;
      break;
      }
    else {
/*------------------------------------------------------------------------------
      assumes latitude upward (i.e. file starts with highest latitudes)*/
      downward=1;
      break;
      }
    }

  switch (arrangement) {
    case ROW:
      for (j=0;j<grid->ny;j++) {
        for (i=0;i<grid->nx;i++) {
/*------------------------------------------------------------------------------
          assumes latitude upward*/
//          n=i*grid->ny+j;
          if(downward==0)
            n=j*grid->nx+i;
          else
            n=(grid->ny-j-1)*grid->nx+i;
          m=j*grid->nx+i;
          grid->x[m]=x[n];
          grid->y[m]=y[n];
//          grid->z[m]=y[n];
          }
        }
      break;

    case COLUMN:
      for (j=0;j<grid->ny;j++) {
        for (i=0;i<grid->nx;i++) {
          if(downward==0)
/*------------------------------------------------------------------------------
            assumes latitude upward*/
            n=grid->ny*i+j;
          else
/*------------------------------------------------------------------------------
            assumes latitude downward*/
            n=grid->ny-j-1+grid->ny*i;
          m=j*grid->nx+i;
          grid->x[m]=x[n];
          grid->y[m]=y[n];
//          grid->z[m]=y[n];
          }
        }
      break;
     }

  grid->modeH=2;

  delete[] x;
  delete[] y;

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int xyz_loadgrid (const char* filename, char *proj4, grid_t *grid, int *ncol)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  FILE *in;
  int downward=0,cartesian=0,mode=0,complete=0,masked=0,option,status;
  int m;
  double mask,x,y;
  char *multic=0;
  char line[1024];

  status=0;

  in=fopen(filename,"rb");
  if (in == NULL) {
    return(-1);
    }

  fgets(line, 1024, in);

/*------------------------------------------------------------------------------
  decode header informations*/

  mode     =(strstr(line, "YXZ")!=0);
  complete =(strstr(line, "COMPLETE")!=0);
  cartesian=(strstr(line, "CARTESIAN")!=0);
  downward =(strstr(line, "DOWNWARD")!=0);
  masked   =(strstr(line, "MASK=")!=0);
  multic=strstr(line, "NCOL=");

  if(multic!=0) {
    multic+=strlen("NCOL=");
    sscanf(multic, "%d",ncol);
    }
  else {
    *ncol=1;
    } 

  fclose(in);

  switch (complete) {
    case 0:
      status=xyz_readgrid_00 (filename, proj4, grid, mode);
      printf("status=%d\n",status);
      break;

    case 1:
      status=xyz_readgrid_01 (filename, grid, mode, cartesian);
      break;

    default:
       break;
    }
    
//   printf("#################################################################\n");
//   printf("native grid : \n");
//   map_printgrid(*grid);

  if(cartesian==1) {
    projPJ ref;
    const char *parms = "+proj=merc +ellps=WGS84 +lon_0=0E +lat_ts=38N";

    if ( ! (ref = pj_init_plus(parms)) ) {
      __ERR_BASE_LINE__( "Projection initialization failed\n"); exit(1);
      }
    printf("#################################################################\n");
    printf("invert projection with parameters : %s\n",parms);

    for (size_t j=0;j<grid->ny;j++) {
      for (size_t i=0;i<grid->nx;i++) {
        m=j*grid->nx+i;
        x=grid->x[m];
        y=grid->y[m];
        projection_to_geo(ref,&grid->y[m],&grid->x[m],x,y);
        }
      }
    pj_free(ref);
    
    map_minmax(grid);
    
    grid->dx=(grid->xmax-grid->xmin) / ((double) grid->nx-1.);
    grid->dy=(grid->ymax-grid->ymin) / ((double) grid->ny-1.);
    }

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int xyz_read00_r1(FILE *in, grid_t grid, int pos, int ncol, float *buf, float mask, int mode, int cartesian)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int nitems,ndata,nlon,nlat;
  int i,j,k,l,m,n;
  int downward=0,complete=0,option,status;
  char line[1024];
  double xx,yy,tmp;
//  float mask=-9999.0;
  range_t<double> r;
  float *dum,z,w;
  double *array;

  status=0;

  ndata=grid.nx*grid.ny;

  dum=new float[ndata];

  array=new double[ncol+2];

  for (n=0;n<grid.ny*grid.nx;n++) buf[n]=mask;

  for(n=0;n<ndata;n++) {
    for(k=0;k<ncol+2;k++) nitems=fscanf(in,"%lf",&array[k]);
    fgets(line,1024,in);
    dum[n]=array[pos+2];
//    if(nitems==-1) break;
    switch (mode) {
      case 0:
        xx=array[0];
        yy=array[1];
        break;

      case 1:
        xx=array[1];
        yy=array[0];
        break;
      }
    if(cartesian!=1) xx=geo_recale(xx,grid.xmin,180.);
    for (j=0;j<grid.ny;j++) {
      m=j*grid.nx;
      if(yy==grid.y[m]) break;
      }
    if(j==grid.ny) {
      printf("#xyz_read00_r1: %lf %lf\n",xx,yy);
      }
    for (i=0;i<grid.nx;i++) {
      m=i;
      if(xx==grid.x[m]) break;
      }
    if(i==grid.nx) {
      printf("#xyz_read00_r1: %lf %lf\n",xx,yy);
      }
    m=j*grid.nx+i;
    buf[m]=dum[n];
    if(feof(in)) break;
    }

  delete[] dum;
  delete[] array;

  return(0);
  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int xyz_read01_r1(FILE *in, grid_t grid, int pos, int ncol, float *buf,int mode)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int nitems,ndata,nlon,nlat;
  int i,j,k,l,m,n;
  int downward=0,cartesian=0,complete=0,arrangement,option,status;
  char line[1024];
  double *x,yy,tmp,mask=-9999.0;
  range_t<double> r;
  float *dum,z;
  double *array;

  status=0;

  ndata=grid.nx*grid.ny;

  x=new double[ndata];

  array=new double[ncol+2];

  for(n=0;n<ndata;n++) {
    for(k=0;k<ncol+2;k++) nitems=fscanf(in,"%lf",&array[k]);
    buf[n]=array[pos+2];
    switch (mode) {
      case 0:
        x[n]=array[0];
        break;

      case 1:
        x[n]=array[1];
        break;
      }
    }

  if(x[0]==x[1])
    arrangement=COLUMN;
  else
    arrangement=ROW;

  delete[] x;

  switch (arrangement) {
    case COLUMN:
/*------------------------------------------------------------------------------
      assumes latitude upward*/
      dum=new float[grid.ny*grid.nx];
      for (j=0;j<grid.ny;j++) {
        for (i=0;i<grid.nx;i++) {
          n=i*grid.ny+j;
          m=j*grid.nx+i;
          dum[m]=buf[n];
          }
        }
      for (j=0;j<grid.ny;j++) {
        for (i=0;i<grid.nx;i++) {
          m=j*grid.nx+i;
          buf[m]=dum[m];
          }
        }
      delete[] dum;
      break;

    case ROW:
/*------------------------------------------------------------------------------
      assumes latitude downward*/
      for (j=0;j<grid.ny/2;j++) {
        for (i=0;i<grid.nx;i++) {
          n=(grid.ny-j-1)*grid.nx+i;
          m=j*grid.nx+i;
          tmp=buf[m];
          buf[m]=buf[n];
          buf[n]=tmp;
          }
        }
      break;
     }

  delete[] array;

  return(0);
  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int xyz_loadr1_template(const char *filename, grid_t & grid, int pos, T *buf, T *mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  {
  FILE *in;
  int downward=0,cartesian=0,mode=0,complete=0,option,status;
  char *masked=0,*multic=0;
  char line[1024];
  int ncol;

  //*mask=-99999.0;

  status=0;

  in=fopen(filename,"rb");
  if (in == NULL) {
    return(-1);
    }

/*------------------------------------------------------------------------------
  decode header informations*/
  fgets(line, 1024, in);
  masked=strstr(line, "MASK=");

  if(masked!=0) {
    masked+=strlen("MASK=");
    sscanf(masked, "%f",mask);
    }

  multic=strstr(line, "NCOL=");

  if(multic!=0) {
    multic+=strlen("NCOL=");
    sscanf(multic, "%d",&ncol);
    }
  else {
    ncol=1;
    }

  mode=(strstr(line, "YXZ")!=0);
  complete=(strstr(line, "COMPLETE")!=0);
  cartesian=(strstr(line, "CARTESIAN")!=0);
  downward=(strstr(line, "DOWNWARD")!=0);

  switch (complete) {
    case 0:
      cartesian=1;
      status=xyz_read00_r1(in,grid,pos,ncol,buf,*mask,mode,cartesian);
      break;

    case 1:
      status=xyz_read01_r1(in,grid,pos,ncol,buf,mode);
      break;

    default:
       break;
    }

  fclose(in);

  return(status);
  }

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int xyz_loadr1(const char *filename, grid_t & grid, int pos, float *buf, float *mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  {
  int status;
  
  status=xyz_loadr1_template(filename, grid, pos, buf, mask);
  
  return(status);
  }


// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   int xyz_loadr1 (const char *filename, grid_t & grid, int pos, complex<float> *buf, complex<float> *mask)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
//   {
//   int status;
//   
//   status=xyz_loadr1_template(filename, grid, pos, buf, mask);
//   
//   return(status);
//   }


