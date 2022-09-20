
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
#include <stdarg.h>

#include "tools-structures.h"
#include "constants.h"
#include "map.h"
#include "geo.h"
#include "xyz.h"
#include "spectrum.h"

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void ascii_printheader (bmgheader_t header)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{

  printf("#ascii_printheader \n");
  printf(" Remark  : %s \n", header.comment[0]);
  printf(" Remark  : %s \n", header.comment[1]);
  printf(" Remark  : %s \n", header.comment[2]);
  printf(" Remark  : %s \n", header.comment[3]);

  printf(" nx      : %d \n", header.ni);
  printf(" ny      : %d \n", header.nj);
  printf(" nk      : %d \n", header.nk);
  printf(" nt      : %d \n", header.nt);
  printf(" nd      : %d \n", header.nd);

  printf(" xmin    : %f \n", header.xmin);
  printf(" ymin    : %f \n", header.ymin);
  printf(" dx      : %f \n", header.dx);
  printf(" dy      : %f \n", header.dy);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void ascii_skipline(FILE *in)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  char c;
  
  do{
    c = fgetc(in);
    } while((c != '\n') and not feof(in));
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int ascii_readheader (FILE *in,grid_t *grid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int nitems,status=0;

  nitems=fscanf(in,"%lf %lf",&grid->xmin,&grid->xmax);

  nitems=fscanf(in," %lf %lf",&grid->ymin,&grid->ymax);
  nitems=fscanf(in," %lf %lf",&grid->dx,&grid->dy);
  nitems=fscanf(in,"  %d %d",&grid->nx,&grid->ny);
  
  grid->modeH=0;
  
  status= map_completegridaxis(grid,1);

  ascii_skipline(in);

  return(status);
  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int ascii_readheader_GOT (FILE *in, grid_t *grid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int nitems,status=0;

  ascii_skipline(in);
  ascii_skipline(in);
  
  grid->modeH=0;
  
  nitems=fscanf(in,"%d %d",&grid->ny,&grid->nx);
  nitems=fscanf(in,"%lf %lf",&grid->ymin,&grid->ymax);
  nitems=fscanf(in,"%lf %lf",&grid->xmin,&grid->xmax);
  
  grid->dx=(grid->xmax-grid->xmin)/((double) grid->nx-1.);
  grid->dy=(grid->ymax-grid->ymin)/((double) grid->ny-1.);

  return(status);
  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int ascii_readheader_SLIM (FILE *in, grid_t *grid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int nitems,status=0;
  
  nitems=fscanf(in,"%lf %lf",&grid->xmin,&grid->ymin,&grid->zmin);
  nitems=fscanf(in,"%lf %lf",&grid->dx,&grid->dy,&grid->dz);
  nitems=fscanf(in,"%d %d",&grid->nx,&grid->ny,&grid->nz);
  
  grid->xmax=grid->xmin*((double) grid->nx-1.)*grid->dx;
  grid->ymax=grid->ymin*((double) grid->ny-1.)*grid->dy;

  return(status);
  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int ascii_readheader_GIS (FILE *in, grid_t *grid, const char *proj4_options)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int nitems,status=0;
//   char c;
  char key[256];
  double t,p,x,y;

//   ascii_skipline(in);
//   ascii_skipline(in);
  
  nitems=fscanf(in,"%s %d",key,&grid->nx);
  nitems=fscanf(in,"%s %d",key,&grid->ny);
  nitems=fscanf(in,"%s %lf",key,&grid->xmin);
  nitems=fscanf(in,"%s %lf",key,&grid->ymin);
  nitems=fscanf(in,"%s %lf",key,&grid->dx);
//  nitems=fscanf(in,"%s %d",key,&grid->mask);

  grid->dy=grid->dx;

  grid->xmax=grid->xmin+(double) (grid->nx-1.)*grid->dx;
  grid->ymax=grid->ymin+(double) (grid->ny-1.)*grid->dy;
  
  grid->modeH=0;
  
  if(proj4_options!=0) {
    projPJ ref;
    int i,j,n;

    ref = init_projection(proj4_options,true);
//     int nprocs=initialize_OPENMP(-1);
// #pragma omp parallel for private(t,p) if(nprocs>1)

    grid->x=new double[grid->Hsize()];
    grid->y=new double[grid->Hsize()];
    for(j=0;j<grid->ny;j++) {
      y=grid->ymin+j*grid->dy;
      for(i=0;i<grid->nx;i++) {
        n=grid->nx*j+i;
        x=grid->xmin+i*grid->dx;
        projection_to_geo(ref,&p,&t,x,y);
        grid->x[n]=t;
        grid->y[n]=p;
        }
      }
   grid->modeH=2;
   status=map_minmax(grid);
   
   pj_free(ref);
   }

  return(status);
  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int ascii_writeheader_GIS(FILE* in, grid_t & grid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  {
  int status,nitems;
  
  rewind(in);

// NODATA_VALUE -99999.00
  status=-1;
  
  nitems=fprintf(in,"NCOLS %d \n", grid.nx);
  nitems=fprintf(in,"NROWS %d \n", grid.ny);
  nitems=fprintf(in,"XLLCORNER %lf \n", grid.xmin);
  nitems=fprintf(in,"YLLCORNER %lf \n", grid.ymin);
  nitems=fprintf(in,"CELLSIZE %lf \n", grid.dx);

  status=0;

  return(status);
  }

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int ascii_writeheader (FILE *in,bmgheader_t header)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;

  status=0;

  rewind(in);

  fprintf(in," %f %f \n",header.xmin,header.xmin+(header.ni-1)*header.dx);
  fprintf(in," %f %f \n",header.ymin,header.ymin+(header.nj-1)*header.dy);
  fprintf(in," %f %f \n",header.dx,header.dy);
  fprintf(in," %d %d \n",header.ni,header.nj);
/*   fprintf(in," %f \n",header.spec); */

  status=0;

  return(status);
  }

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int ascii_createheader (const char* file,int kmax,int dmax,int tmax,grid_t grid,float *levels,float mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  {
  FILE *out=NULL;
  bmgheader_t header;

  out=fopen(file,"wb");
  if(out==0) {
      TRAP_ERR_EXIT(-1,"file opening issue : %s \n",file);
    }

  strcpy(header.comment[0],"no comment");
  strcpy(header.comment[1],"no comment");
  strcpy(header.comment[2],"no comment");
  strcpy(header.comment[3],"no comment");
  header.ni=grid.nx;
  header.nj=grid.ny;
  header.xmin=grid.xmin;
  header.ymin=grid.ymin;
  header.dx=grid.dx;
  header.dy=grid.dy;
  header.nk=kmax;
  header.nd=dmax;
  header.nt=tmax;
  header.spec=mask;
  header.code=0;
  header.levels=levels;

  ascii_writeheader (out, header);

  fclose(out);
  return(0);
  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void ascii_loadgrid (const char* filename,grid_t *grid, int *status)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  {
  FILE *in=NULL;

  *status=0;

  in=fopen(filename,"rb");
  if (in == NULL) {
    *status=-1;
    return;
    }

 *status=ascii_readheader(in,grid);
  if (*status !=0) {
    *status=-1;
    return;
    }

   fclose(in);

  grid->xmax=grid->xmin+(grid->nx-1)*grid->dx;
  grid->ymax=grid->ymin+(grid->ny-1)*grid->dy;
  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int ascii_getinfo (const char* filename, grid_t *header)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  {
  FILE *in=NULL;
  int status;

  status=0;

  in=fopen(filename,"rb");
  if (in == NULL) {
    status=-1;
    return(status);
    }

  status=ascii_readheader(in,header);
  if (status !=0) {
    status=-1;
    return(status);
    }

  fclose(in);
  status=0;
  return(status);
  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void ascii_checkfile (const char *filename, int *status)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  {
  grid_t grid;

  ascii_loadgrid (filename, &grid, status);
  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int ascii_load (const char* filename, grid_t *grid,float **buf,float *mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  {
  int status=0;
  int i,j,nitems;
  FILE *in=NULL;
  float dum;

  in=fopen(filename,"r");
  if(in==NULL) {
    TRAP_ERR_EXIT(-1,"exiting\n");
    }

  status= ascii_readheader (in,grid);
  
  *buf=new float[grid->nx*grid->ny];
  
  fscanf(in," %f",&dum);
  ascii_skipline(in);

  for (j=0;j<grid->ny;j++) {
    for (i=0;i<grid->nx;i++) {
      int jj=grid->ny-j-1;
      nitems=fscanf(in,"%f",&((*buf)[jj*grid->nx+i]));
      if(nitems!=1) {
        printf("anomaly!!!\n");
        }
      }
    }

  fclose(in);

  *mask=dum;
  return(status);
  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int ascii_load (const char* filename, grid_t *grid,short **buf,short *mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=0;
  int i,j,nitems;
  FILE *in=NULL;
  float dum;

  in=fopen(filename,"r");
  if(in==NULL) {
    TRAP_ERR_EXIT(-1,"exiting\n");
    }

  status= ascii_readheader (in,grid);
  
  *buf=new short[grid->nx*grid->ny];
  
  fscanf(in," %d",&dum);
  ascii_skipline(in);

  for (j=0;j<grid->ny;j++) {
    for (i=0;i<grid->nx;i++) {
      nitems=fscanf(in,"%d",&(buf[j*grid->nx+i]));
      if(nitems!=1) {
        printf("anomaly!!!\n");
        }
      }
    }

  fclose(in);

  *mask=dum;
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> int ascii_load_SLIM_template (const char* filename, grid_t *grid,T **buf,T *mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=0;
  int i,j,nitems;
  FILE *in=NULL;
  float dum;
  char msg[1024]/*,c*/;

  in=fopen(filename,"r");
  if(in==NULL) {
    sprintf(msg, "Cannot open harmonic data file %s \n",filename);
    TRAP_ERR_EXIT(-1,"exiting\n");
    }

  status= ascii_readheader_SLIM (in,grid);
  
  *buf=new T[grid->nx*grid->ny];
  
//   fscanf(in," %d",&dum);
//   ascii_skipline(in);

  for (j=0;j<grid->ny;j++) {
    for (i=0;i<grid->nx;i++) {
      nitems=fscanf(in,"%lf",&dum);
      (*buf)[j*grid->nx+i]=(T) dum;
      if(nitems!=1) {
        printf("anomaly!!!\n");
        }
      }
    }

  fclose(in);

//   *mask=dum;
  *mask=(T)1.e+10;
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int ascii_load_SLIM (const char* filename, grid_t *grid, const char *proj, short **buf, short *mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=0;
  
  status=ascii_load_SLIM_template (filename, grid, buf, mask);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int ascii_load_SLIM (const char* filename, grid_t *grid, const char *proj, float **buf, float *mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=0;
  
  status=ascii_load_SLIM_template (filename, grid, buf, mask);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int ascii_load_SLIM (const char* filename, grid_t *grid, const char *proj, double **buf, double *mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=0;
  
  status=ascii_load_SLIM_template (filename, grid, buf, mask);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> int ascii_load_GIS_template (const char* filename, grid_t *grid, const char *proj, T **buf, T *mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  {
  int status=0;
  int i,j,nitems;
  FILE *in=NULL;
  float dum;
  char msg[1024],key[256];

  in=fopen(filename,"r");
  if(in==NULL) {
    sprintf(msg, "Cannot open harmonic data file %s \n",filename);
    TRAP_ERR_EXIT(-1,"exiting\n");
    }

  status= ascii_readheader_GIS (in, grid, proj);
  
  *buf=new T[grid->nx*grid->ny];
  
  fscanf(in,"%s %f",key, &dum);
  *mask=dum;
  ascii_skipline(in);

  for (j=0;j<grid->ny;j++) {
    for (i=0;i<grid->nx;i++) {
      int jj=grid->ny-j-1;
      nitems=fscanf(in,"%f",&dum);
      (*buf)[jj*grid->nx+i]=dum;
      if(nitems!=1) {
        printf("anomaly!!!\n");
        }
      }
    }

  fclose(in);

  return(status);
  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int ascii_load_GIS (const char* filename, grid_t *grid, const char *proj, float **buf, float *mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  {
  int status=0;
  
  status=ascii_load_GIS_template (filename, grid, proj, buf, mask);
  
  return(status);
  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int ascii_load_GIS(const char* filename, grid_t *grid, const char *proj, short **buf, short *mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=0;
  
  status=ascii_load_GIS_template (filename, grid, proj, buf, mask);
  
  return(status);
  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int ascii_save_r(FILE *out, int t, grid_t grid, float *buf, float mask, double time)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,j,k;

  fprintf(out,"%6.1f %lf (mask,time) \n",mask,time);

  for (j=0;j<grid.ny;j++) {
    for (i=0;i<grid.nx;i+=30) {
      for(k=i;k<min(i+30,grid.nx);k++) fprintf(out," %6.1f",buf[j*grid.nx+k]);
      fprintf(out,"\n");
      }
    }
  return(0);
  }

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int ascii_saver1(const char* filename, grid_t & grid, float *buf, float mask,const char *format, int step)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------

  Comments: ascii_saver1 save an unique time frame,dimension,level

--------------------------------------------------------------------------*/
{
  FILE *out=NULL;
  int i,j,k;
  bmgheader_t header;

  out=fopen(filename,"wb");
  if(out==0) {
      TRAP_ERR_EXIT(-1,"file opening issue : %s \n",filename);
    }

  header.ni=grid.nx;
  header.nj=grid.ny;
  header.xmin=grid.xmin;
  header.ymin=grid.ymin;
  header.dx=grid.dx;
  header.dy=grid.dy;
  header.nk=1;
  header.nd=1;
  header.nt=1;
  header.spec=mask;

  ascii_writeheader (out, header);

  fprintf(out," %f (mask value) \n",mask);

  for (j=0;j<header.nj;j++) {
    for (i=0;i<header.ni;i+=step) {
      for(k=i;k<min(i+step,header.ni);k++) fprintf(out,format,buf[j*header.ni+k]);
      fprintf(out,"\n");
      }
    }
  fclose(out);
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int ascii_save_GIS(const char* filename, grid_t & grid, float *buf, float mask, const char *format, int step)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status, nitems;
  FILE *out=NULL;
  int i,j,k;

  out=fopen(filename,"w");
  if(out==0) {
    TRAP_ERR_EXIT(-1,"file opening issue : %s \n",filename);
    }

  status=ascii_writeheader_GIS(out, grid);

  fprintf(out,"NODATA_VALUE %f\n",mask);

  for (j=0;j<grid.ny;j++) {
    for (i=0;i<grid.nx;i+=step) {
      for(k=i;k<min(i+step,grid.nx);k++) nitems=fprintf(out,format,buf[j*grid.nx+k]);
      fprintf(out,"\n");
      }
    }
  fclose(out);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void polar_buf(int n,const float *a,const float *G,float dum,fcomplex **buffer,fcomplex *mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int m;
  
  *buffer=new fcomplex[n];
  
  *mask=fcomplex(dum,dum);
  
  for(m=0;m<n;m++) {
    if(a[m]!=dum) {
      (*buffer)[m]=polar<float>(a[m],-G[m]*d2r);
      }
    else {
      (*buffer)[m]=*mask;
      }
    }
  
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int ascii_loadc2(const char* filename,grid_t *grid,fcomplex **buf1,fcomplex **buf2,fcomplex *mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,j,k,nitems,status;
  char msg[1024];
  FILE *in=NULL;
  float *a=NULL,*G=NULL,dum;

  in=fopen(filename,"r");
  if(in==NULL) {
    sprintf(msg, "Cannot open harmonic data file %s \n",filename);
    TRAP_ERR_EXIT(-1,"exiting\n");
    }

  status= ascii_readheader (in,grid);
  
  const int n=grid->Hsize();
  
  a=new float[n];
  G=new float[n];

  for(int ibuf=0;ibuf<2;ibuf++){
    
  fscanf(in," %f",&dum,&dum);
  ascii_skipline(in);
  for (j=0;j<grid->ny;j++) {
    for (i=0;i<grid->nx;i+=30) {
      for(k=i;k<min(i+30,grid->nx);k++) {
        nitems=fscanf(in,"%f",&(a[j*grid->nx+k]));
        if(nitems!=1) {
          printf("anomaly!!!\n");
          }
        }
      for(k=i;k<min(i+30,grid->nx);k++) {
        nitems=fscanf(in,"%f",&(G[j*grid->nx+k]));
        if(nitems!=1) {
          printf("anomaly!!!\n");
          }
        }
      }
    }
    
    if(ibuf==0)
      polar_buf(n,a,G,dum,buf1,mask);
    else
      polar_buf(n,a,G,dum,buf2,mask);
    
    if(buf2==0) break;
    }
  
  fclose(in);
  
  delete[] a;
  delete[] G;

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int ascii_loadc1_got(const char* filename,grid_t *grid,fcomplex **buffer,fcomplex *mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  {
  int i,j,k,nitems,status;
  char msg[1024];
  FILE *in=NULL;
  float *a=NULL,*G=NULL,dum;
  int nvpl=11;

  in=fopen(filename,"r");
  if(in==NULL) {
    sprintf(msg, "Cannot open harmonic data file %s \n",filename);
    TRAP_ERR_EXIT(-1,"exiting\n");
    }

  status= ascii_readheader_GOT (in, grid);

  const int n=grid->Hsize();
  
  a=new float[n];
  G=new float[n];

  fscanf(in," %f %f",&dum,&dum);
  ascii_skipline(in);
  ascii_skipline(in);

  for (j=0;j<grid->ny;j++) {
    for (i=0;i<grid->nx;i+=nvpl) {
      for(k=i;k<min(i+nvpl,grid->nx);k++) {
        nitems=fscanf(in,"%f",&(a[j*grid->nx+k]));
        if(nitems!=1) {
          printf("anomaly!!!\n");
          }
        }
      }
    }
  
  for(k=0;k<8;k++) ascii_skipline(in);
  
  for (j=0;j<grid->ny;j++) {
    for (i=0;i<grid->nx;i+=nvpl) {
      for(k=i;k<min(i+nvpl,grid->nx);k++) {
        nitems=fscanf(in,"%f",&(G[j*grid->nx+k]));
        if(nitems!=1) {
          printf("anomaly!!!\n");
          }
        }
      }
    }

  fclose(in);
  
  polar_buf(n,a,G,dum,buffer,mask);
  
  delete[] a;
  delete[] G;

  return(0);
  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int ascii_loadc1_got_complex(const char* filename,grid_t *grid,fcomplex **buffer,fcomplex *mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  {
  int i,j,k,nitems,status;
  char msg[1024];
  FILE *in=NULL;
  float *a=NULL,*G=NULL,*zR=NULL,*zI=NULL,dum;
  int nvpl=11;

  in=fopen(filename,"r");
  if(in==NULL) {
    sprintf(msg, "Cannot open harmonic data file %s \n",filename);
    TRAP_ERR_EXIT(-1,"exiting\n");
    }

  status= ascii_readheader_GOT (in, grid);

  const int n=grid->Hsize();
  
  a=new float[n];
  G=new float[n];

  zR=new float[n];
  zI=new float[n];

  fscanf(in," %f %f",&dum,&dum);
  ascii_skipline(in);
  ascii_skipline(in);

  for (j=0;j<grid->ny;j++) {
    for (i=0;i<grid->nx;i+=nvpl) {
      for(k=i;k<min(i+nvpl,grid->nx);k++) {
        nitems=fscanf(in,"%f",&(zR[j*grid->nx+k]));
        if(nitems!=1) {
          printf("anomaly!!!\n");
          }
        }
      }
    }
  
  for(k=0;k<8;k++) ascii_skipline(in);
  
  for (j=0;j<grid->ny;j++) {
    for (i=0;i<grid->nx;i+=nvpl) {
      for(k=i;k<min(i+nvpl,grid->nx);k++) {
        nitems=fscanf(in,"%f",&(zI[j*grid->nx+k]));
        if(nitems!=1) {
          printf("anomaly!!!\n");
          }
        }
      }
    }

  fclose(in);
  
//   polar_buf(n,a,G,dum,buffer,mask);
  *buffer=new fcomplex[n];
  for (j=0;j<grid->ny;j++) {
    for (i=0;i<grid->nx;i++) {
      size_t m=grid->nx*j+i;
      (*buffer)[m]=fcomplex(zR[m],zI[m]);
      }
    }
  *mask=fcomplex(dum,dum);
  
  delete[] a;
  delete[] G;

  return(0);
  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int ascii_savec1 (const char* filename,grid_t grid,fcomplex *buf,fcomplex mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  {
  FILE *out=NULL;
  int i,j,k;
  bmgheader_t header;
  float *a=NULL,*G=NULL;

  out=fopen(filename,"wb");
  if(out==0) {
      TRAP_ERR_EXIT(-1,"file opening issue : %s \n",filename);
    }

  header.ni=grid.nx;
  header.nj=grid.ny;
  header.xmin=grid.xmin;
  header.ymin=grid.ymin;
  header.dx=grid.dx;
  header.dy=grid.dy;
  header.nk=1;
  header.nd=1;
  header.nt=1;
  header.spec=real(mask);

  ascii_writeheader (out, header);

  fprintf(out," %6.1f %6.1f (mask value) \n",real(mask),imag(mask));

  exitIfNull(
    a=(float *) malloc(header.ni*header.nj*sizeof(float))
    );
  exitIfNull(
    G=(float *) malloc(header.ni*header.nj*sizeof(float))
    );
  for (j=0;j<header.nj;j++) {
    for (i=0;i<header.ni;i++) {
      k=j*header.ni+i;
      if(buf[k] != mask) {
/*
        a[k]=100.*sqrt(buf[k].r*buf[k].r+buf[k].i*buf[k].i);
        G[k]=atan2(-buf[k].i,buf[k].r)*r2d;
*/
        a[k]=abs(buf[k]);
        G[k]=-arg(buf[k])*r2d;
        if (G[k]<0.) G[k]+=360.;
        }
      else {
        a[k]=header.spec;
        G[k]=header.spec;
        }
      }
    }
  
  for (j=0;j<header.nj;j++) {
    for (i=0;i<header.ni;i+=30) {
      for(k=i;k<min(i+30,header.ni);k++) fprintf(out," %6.1f",a[j*header.ni+k]);
      fprintf(out,"\n");
      for(k=i;k<min(i+30,header.ni);k++) fprintf(out," %6.1f",G[j*header.ni+k]);
      fprintf(out,"\n");
      }
    }
  fclose(out);
/*   printf("#ascii_savec1 ... finished\n"); */

  free(a);
  free(G);

  return(0);
  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int ascii_savec2 (const char* filename,grid_t grid,fcomplex *bufx,fcomplex *bufy, fcomplex mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  {
  FILE *out=NULL;
  int i,j,k,step=30;
  bmgheader_t header;
  float levels[1]={0.};
  float *a=NULL,*G=NULL;

  out=fopen(filename,"wb");
  if(out==0) {
      TRAP_ERR_EXIT(-1,"file opening issue : %s \n",filename);
    }

  strcpy(header.comment[0],"no comment");
  strcpy(header.comment[1],"no comment");
  strcpy(header.comment[2],"no comment");
  strcpy(header.comment[3],"no comment");
  header.ni=grid.nx;
  header.nj=grid.ny;
  header.xmin=grid.xmin;
  header.ymin=grid.ymin;
  header.dx=grid.dx;
  header.dy=grid.dy;
  header.nk=1;
  header.nd=1;
  header.nt=1;
  header.spec=real(mask);
  header.levels=levels;

  ascii_writeheader (out, header);

  fprintf(out," %6.1f %6.1f (mask value) \n",real(mask),imag(mask));

  exitIfNull(
    a=(float *) malloc(header.ni*header.nj*sizeof(float))
    );
  exitIfNull(
    G=(float *) malloc(header.ni*header.nj*sizeof(float))
    );

  for (j=0;j<header.nj;j++) {
    for (i=0;i<header.ni;i++) {
      k=j*header.ni+i;
      if(bufx[k] != mask) {
/*
        a[k]=100.*sqrt(bufx[k].r*bufx[k].r+bufx[k].i*bufx[k].i);
        G[k]=atan2(-bufx[k].i,bufx[k].r)*r2d;
*/
        a[k]=abs(bufx[k]);
        G[k]=-arg(bufx[k])*r2d;
        if (G[k]<0.) G[k]+=360.;
        }
        else
        {
        a[k]=header.spec;
        G[k]=header.spec;
        }
      }
    }
  step=30;
  for (j=0;j<header.nj;j++) {
    for (i=0;i<header.ni;i+=30) {
      for(k=i;k<min(i+30,header.ni);k++) fprintf(out," %6.1f",a[j*header.ni+k]);
      fprintf(out,"\n");
      for(k=i;k<min(i+30,header.ni);k++) fprintf(out," %6.1f",G[j*header.ni+k]);
      fprintf(out,"\n");
      }
    }

  fprintf(out," %6.1f %6.1f (mask value) \n",real(mask),imag(mask));

  for (j=0;j<header.nj;j++) {
    for (i=0;i<header.ni;i++) {
      k=j*header.ni+i;
      if(bufy[k] != mask) {
/*
        a[k]=100.*sqrt(bufy[k].r*bufy[k].r+bufy[k].i*bufy[k].i);
        G[k]=atan2(-bufy[k].i,bufy[k].r)*r2d;
        if (G[k]<0.) G[k]+=360.;
*/
        a[k]=abs(bufy[k]);
        G[k]=-arg(bufy[k])*r2d;
        }
        else
        {
        a[k]=header.spec;
        G[k]=header.spec;
        }
      }
    }
  for (j=0;j<header.nj;j++) {
    for (i=0;i<header.ni;i+=30) {
      for(k=i;k<min(i+30,header.ni);k++) fprintf(out," %6.1f",a[j*header.ni+k]);
      fprintf(out,"\n");
      for(k=i;k<min(i+30,header.ni);k++) fprintf(out," %6.1f",G[j*header.ni+k]);
      fprintf(out,"\n");
      }
    }
  fclose(out);
/*   printf("#ascii_savec2 ... finished\n"); */

  free(a);
  free(G);

  return(0);

  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int shomcst_loadgrid (const char* filename,grid_t *grid, int mode)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *-----------------------------------------------------------------------------

  Assumes grid shape is rectangular, and data set is partial (missing values
  not given through a mask flag)

------------------------------------------------------------------------------*/
{
  FILE *in;
  int add,nitems,ndata,nlon,nlat;
  int i,j,k,m,n;
  int status;
  char line[256*256];
  double *x,*y;
  double *xtics,*ytics,dum;
  range_t<double> r;

  status=0;

  in=fopen(filename,"rb");
  if (in == NULL) {
    return(-1);
    }

  fgets(line, 1024, in);
  nitems=sscanf(line,"%lf %lf %lf %lf %d %d",&dum,&dum,&dum,&dum,&(grid->nx),&(grid->ny));

  ndata=grid->nx*grid->ny;

  x=new double[ndata];
  y=new double[ndata];

  for(n=0;n<ndata;n++) {
    fgets(line,1024,in);
    nitems=sscanf(line,"%lf %lf",&y[n],&x[n]);
    fgets(line,256*256,in);
    fgets(line,256*256,in);
    if(feof(in)) {
//      printf("arrgll\n");
//       fclose(in);
//       return(-1);
      break;
      }
    }

  fclose(in);

  ndata=n;

  if(x[0]>180.) x[0]-=360.0;
  for(n=0;n<ndata;n++) {
    x[n]=geo_recale(x[n],x[0],180.0);
    }

  xtics=new double[ndata];
  xtics[0]=x[0];
  nlon=1;
  for(n=0;n<ndata-1;n++) {
    m=n+1;
    while(x[n]==x[m]) {
      if(m==ndata) break;
      m++;
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

  ytics=new double[ndata];
  ytics[0]=y[0];
  nlat=1;
  for(n=0;n<ndata-1;n++) {
    m=n+1;
    while(y[n]==y[m]) {
      if(m==ndata) break;
      m++;
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

  printf("nx =%d ny =%d ndata =%d (full grid =%d)\n",nlon,nlat,ndata,nlon*nlat);

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

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int shomcst_read (FILE *in, grid_t grid, int pos, int ncol, float *buf,float mask,int mode)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  {
  int nitems,ndata;
  int i,j,k,m,n;
  int status;
  char line[1024];
  double xx,yy;
//  float mask=-9999.0;
  range_t<double> r;
  float *dum;
  double *array;

  status=0;

  ndata=grid.nx*grid.ny;

  dum=new float[ndata];

  array=new double[ncol+2];

  for (n=0;n<grid.ny*grid.nx;n++) buf[n]=mask;

  for(n=0;n<ndata;n++) {
    fgets(line,1024,in);
    nitems=sscanf(line,"%lf %lf",&yy,&xx);
    for(k=0;k<ncol;k++) nitems=fscanf(in,"%6lf",&(array[k]));
    dum[n]=array[pos];
    for(k=0;k<ncol;k++) nitems=fscanf(in,"%6lf",&(array[k]));
    fgets(line,1024,in);
    if(feof(in)) {
      printf("ndata =%d\n",n);
      break;
      }
    xx=geo_recale(xx,grid.xmin,180.);
    for (j=0;j<grid.ny;j++) {
      m=j*grid.nx;
      if(yy==grid.y[m]) break;
      }
    if(j==grid.ny) {
      printf("%lf %lf\n",xx,yy);
      }
    for (i=0;i<grid.nx;i++) {
      m=i;
      if(xx==grid.x[m]) break;
      }
    if(i==grid.nx) {
      printf("%lf %lf\n",xx,yy);
      }
    m=j*grid.nx+i;
    buf[m]=dum[n];
    }


  delete[] dum;
  delete[] array;

  return(0);
  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int shomcst_read (FILE *in, grid_t grid, int pos, int ncol, complex<float> *buf, complex<float>  mask,int mode)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  {
  int nitems,ndata;
  int i,j,k,m,n;
  int status;
  char line[1024];
  double xx,yy;
//  float mask=-9999.0;
  range_t<double> r;
  float *amp,*pha;
  double *array;

  status=0;

  ndata=grid.nx*grid.ny;

  amp=new float[ndata];
  pha=new float[ndata];

  array=new double[ncol+2];

  for (n=0;n<grid.ny*grid.nx;n++) buf[n]=mask;

  for(n=0;n<ndata;n++) {
    fgets(line,1024,in);
    nitems=sscanf(line,"%lf %lf",&yy,&xx);
    fgets(line,1024,in);
    for(k=0;k<ncol;k++) {
      char dum[7];
      strncpy(dum, &line[k*6], 6);
      dum[7]=0;
      nitems=sscanf(dum,"%lf",&(array[k]));
      }
    amp[n]=array[pos];
//    for(k=0;k<ncol;k++) nitems=fscanf(in,"%6lf",&(array[k]));
    fgets(line,1024,in);
    for(k=0;k<ncol;k++) {
      char dum[7];
      strncpy(dum, &line[k*6], 6);
      dum[7]=0;
      nitems=sscanf(dum,"%lf",&(array[k]));
      }
    pha[n]=array[pos];
    if(feof(in)) {
      printf("ndata =%d\n",n);
      break;
      }
    xx=geo_recale(xx,grid.xmin,180.);
    for (j=0;j<grid.ny;j++) {
      m=j*grid.nx;
      if(yy==grid.y[m]) break;
      }
    if(j==grid.ny) {
      printf("%lf %lf\n",xx,yy);
      }
    for (i=0;i<grid.nx;i++) {
      m=i;
      if(xx==grid.x[m]) break;
      }
    if(i==grid.nx) {
      printf("%lf %lf\n",xx,yy);
      }
    m=j*grid.nx+i;
//    buf[m]=amp[n]*cos(-pha[n]*M_PI/180.0);
    buf[m]=amp[n]*complex<float>(cos(-pha[n]), sin(-pha[n]));
//    buf[m]=complex<float>(amp[n], pha[n]);
    }


  delete[] amp;
  delete[] pha;
  delete[] array;

  return(0);
  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> int shomcst_read_template (const char *filename, grid_t grid, int pos, T *buf, T *mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  {
  FILE *in;
  int status;
  char line[1024];
  const int ncol=143;

  //*mask=-99999.0;

  status=0;

  in=fopen(filename,"rb");
  if (in == NULL) {
    return(-1);
    }

/* *-----------------------------------------------------------------------------
  decode header informations*/
  fgets(line, 1024, in);
  status=shomcst_read (in,grid,pos,ncol,buf,*mask,1);
  fclose(in);

  return(status);
  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int shomcst_read (const char *filename, grid_t grid, int pos, float *buf, float *mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  {
  int status;
    
  status=shomcst_read_template (filename, grid, pos, buf, mask);
  
  return(status);
  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int shomcst_read (const char *filename, grid_t grid, int pos, complex<float> *buf, complex<float> *mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  {
  int status;
    
  status=shomcst_read_template (filename, grid, pos, buf, mask);
  
  return(status);
  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int shomcst_loadc1 (const char *filename, grid_t *grid, int pos, complex<float> ***buf, complex<float> *mask, spectrum_t *s)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  {
  int w,status;
  int mode=0;
  
  bool Z0=true;
  spectrum_t reference=spectrum_init_ref("ESTUARINE", Z0);
  
  const int ncol=143;
  const char *waves[ncol]={
  "Z0","SA","SSA","MSM","MM","MSF","MF",\
  "2Q1","SIGMA1","Q1","RHO1","O1","MS1","MP1","A19","M1","KHI1","PI1","P1","S1","K1","PSI1","PHI1","THETA1","J1","SO1","OO1","KQ1",\
  "2MN2S2","2NS2","3M2S2","OQ2 MNK2","","MNS2","MNUS2","2MK2","2N2 2NM2","MU2 2MS2","","N2","NU2","OP2 MSK2","","M(SK)2","M2",\
  "M(KS)2","MKS2","LAMBDA2","L2 2MN2","NKM2","T2","S2","R2","K2","MSN2","KJ2 MKN2","2SM2","SKM2",\
  "MQ3","2MK3","M3","SO3","MS3","MK3","A87","SP3","S3","SK3","K3",\
  "2NMS4","2MMUS4","2MNS4","2MNUS4","3MK4","N4","3MS4","MN4","MNU4","2MSK4","M4","2MKS4","SN4","3MN4ML4","NK4","2SMK4","MT4","MS4","MK4","2SNM4","2MSN4","2MKN4","S4","SK4",\
  "3MNK6","3MNS6","3MNUS6","4MK6","2NM6","4MS6","2MN6","2MNU6","3MSK6","M6","3MKS6","MSN6","4MN.2ML6","MNK6","2MT6","2MS6","2MK6","2SN6","3MSN6","3MKN6","2SM6","MSK6",\
  "2MNS8","2(MN)8","3MN8","3MNU8","2MSK8","M8","4MKS8","2MSN8","3ML8","2MNK8","3MS8","3MK8","2SMN8","4MSN8","MSNK8","4MKN8","2(MS)8","2MSK8",\
  "5MNS10","3M2N10","4MN10","M10","3MSN10","4MS10","4MK10","5MSN10","2MSNK10","3M2S10"};
  
  *mask=9999.;
    
  status=shomcst_loadgrid (filename, grid, mode);
  if(status!=0) return(-1);
  
  s->init(ncol);
  
  if(pos==-1) {
    *buf=new complex<float>*[ncol];
    for(int k=0;k<ncol;k++) {
      (*buf)[k]=new complex<float>[grid->nx*grid->ny];
      w=s->add(reference, waves[k],0);
      if(w==-1) continue;
      status=shomcst_read_template (filename, *grid, k, (*buf)[w], mask);
      }
    }
  else {
    *buf=new complex<float>*[1];
    *buf[0]=new complex<float>[grid->nx*grid->ny];
    w=s->add(reference, waves[pos],0);
    if(w!=-1) status=shomcst_read_template (filename, *grid, pos, *buf[0], mask);
    }
  return(status);
  
  }
