#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
 
#include <sys/stat.h> //stat

#include "tools-structures.h"
#include "map.h"

//define de swap function to be homegenous with the Sun Sparc binary format
//-------------------------------------------------------------------------
#if NEED_SWAP == 1
#include "swap.h"
#endif
//-------------------------------------------------------------------------

/*------------------------------------------------------------------------

  Nom         :  fct_filesize
                 return the file size in bytes
--------------------------------------------------------------------------*/

int filesize (char *cmd)
  {
  int status, size;
  struct stat buf;
  /* The structure stat is explained in /usr/include/sys/stat.h */

  /* mode_t   st_mode;      File mode (see mknod(2)) */
  /* ino_t    st_ino;       Inode number */
  /* dev_t    st_dev;       ID of device containing */
  /*                        a directory entry for this file */
  /* dev_t    st_rdev;      ID of device */
  /*                        This entry is defined only for */
  /*                        char special or block special files */
  /* nlink_t  st_nlink;     Number of links */
  /* uid_t    st_uid;       User ID of the file's owner */
  /* gid_t    st_gid;       Group ID of the file's group */
  /* off_t    st_size;      File size in bytes */
  /* time_t   st_atime;     Time of last access */
  /* time_t   st_mtime;     Time of last data modification */
  /* time_t   st_ctime;     Time of last file status change */
  /*                        Times measured in seconds since */
  /*                        00:00:00 UTC, Jan. 1, 1970 */
  /* long     st_blksize;   Preferred I/O block size */
  /* long     st_blocks;    Number of 512 byte blocks allocated*/

  status=stat(cmd, &buf);
  if(status==0) {
    size=buf.st_size;
    }
  else {
    size=0;
    }

  return (size);
  }

/*------------------------------------------------------------------------

  Name        :  bmg_printheader

--------------------------------------------------------------------------*/
void bmg_printheader (bmgheader_t header)
{

  printf("#bmg_printheader \n");
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

/*------------------------------------------------------------------------
  
  Nom         :  bmg_readheader

--------------------------------------------------------------------------*/

int bmg_readheader (FILE *in,bmgheader_t *header)
  {
  int nlevel,i,status, size, s;
  size_t offset, nread;
  float time;

  status=0;

  offset=4;

  nread=fread(&size, sizeof(int), 1, in );
#if NEED_SWAP == 1
  size=lswap(size);
#endif

  if(size !=80) {
    status=-1;
    return(status);
    }
  nread=fread(header->comment[0], 80, 1, in );
  fseek(in,offset,SEEK_CUR);

  fseek(in,offset,SEEK_CUR);
  nread=fread(header->comment[1], 80, 1, in );
  fseek(in,offset,SEEK_CUR);

  fseek(in,offset,SEEK_CUR);
  nread=fread(header->comment[2], 80, 1, in );
  fseek(in,offset,SEEK_CUR);

  fseek(in,offset,SEEK_CUR);
  nread=fread(header->comment[3], 80, 1, in );
  fseek(in,offset,SEEK_CUR);

  nread=fread(&size, sizeof(int), 1, in );
  nread=fread(&header->ni, sizeof(int), 1, in );
  nread=fread(&header->nj, sizeof(int), 1, in );
  nread=fread(&header->nk, sizeof(int), 1, in );
  nread=fread(&header->nt, sizeof(int), 1, in );
  nread=fread(&header->nd, sizeof(int), 1, in );
/*   printf("header.ni=%d, header.nj=%d\n",header->ni,header->nj); */

  nlevel=header->nk;

#if NEED_SWAP == 1
  size=lswap(size);
  nlevel=lswap(nlevel);
#endif

  if(size == 24) {
    nread=fread(&header->code, sizeof(int), 1, in );
    header->size=4*(4+80+4)+(4+24+4)+(4+20+4)+(4+nlevel*4+4);
    }
  else {
    header->size=4*(4+80+4)+(4+20+4)+(4+20+4)+(4+nlevel*4+4);
    }
  fseek(in,offset,SEEK_CUR);

  fseek(in,offset,SEEK_CUR);
  nread=fread(&header->xmin, sizeof(float), 1, in );
  nread=fread(&header->ymin, sizeof(float), 1, in );
  nread=fread(&header->dx,   sizeof(float), 1, in );
  nread=fread(&header->dy,   sizeof(float), 1, in );
  nread=fread(&header->spec, sizeof(float), 1, in );
  fseek(in,offset,SEEK_CUR);

  exitIfNull(
    header->levels=(float * )malloc(nlevel*sizeof(float))
    );

  fseek(in,offset,SEEK_CUR);
  nread=fread(header->levels, sizeof(float), (size_t) nlevel, in );
  fseek(in,offset,SEEK_CUR);

  status=0;

#if NEED_SWAP == 1
/*  header=hearder_swap(header);*/
  header->ni	= lswap(header->ni);
  header->nj	= lswap(header->nj);
  header->nk	= lswap(header->nk);
  header->nt	= lswap(header->nt);
  header->nd	= lswap(header->nd);
  header->code	= lswap(header->code);
  f3swap(&header->xmin);
  f3swap(&header->ymin);
  f3swap(&header->dx);
  f3swap(&header->dy);
  f3swap(&header->spec);
  for(i=0;i<nlevel;i++) f3swap(&header->levels[i]);

#endif
  return(status);
  }

/*------------------------------------------------------------------------

  Nom         :  bmg_writeheader

--------------------------------------------------------------------------*/
int bmg_writeheader (FILE *in,bmgheader_t header)
  {

  int nlevel,i,status;
  size_t offset, nwrite, nitems;
  int    size;
  float time;

  status=0;

  rewind(in);
  offset=4;
  nlevel=header.nk;

#if NEED_SWAP == 1
  header.ni	= lswap(header.ni);
  header.nj	= lswap(header.nj);
  header.nk	= lswap(header.nk);
  header.nt	= lswap(header.nt);
  header.nd	= lswap(header.nd);
  header.code	= lswap(header.code);
  f3swap(&header.xmin);
  f3swap(&header.ymin);
  f3swap(&header.dx);
  f3swap(&header.dy);
  f3swap(&header.spec);
  for(i=0;i<nlevel;i++)
    f3swap(&header.levels[i]);
#endif

  size=80;

#if NEED_SWAP == 1
  size=lswap(size);
#endif

  nwrite=fwrite(&size, sizeof(int), 1, in );
  if(nwrite !=1) {
    status=-1;
    return(status);
    }
  nwrite=fwrite(header.comment[0], 80, 1, in );
  nwrite=fwrite(&size, sizeof(int), 1, in );

  nwrite=fwrite(&size, sizeof(int), 1, in );
  nwrite=fwrite(header.comment[1], 80, 1, in );
  nwrite=fwrite(&size, sizeof(int), 1, in );

  nwrite=fwrite(&size, sizeof(int), 1, in );
  nwrite=fwrite(header.comment[2], 80, 1, in );
  nwrite=fwrite(&size, sizeof(int), 1, in );

  nwrite=fwrite(&size, sizeof(int), 1, in );
  nwrite=fwrite(header.comment[3], 80, 1, in );
  nwrite=fwrite(&size, sizeof(int), 1, in );

  size=6*sizeof(int);
#if NEED_SWAP == 1
  size=lswap(size);
#endif
  nwrite=fwrite(&size, sizeof(int), 1, in );
  nwrite=fwrite(&header.ni, sizeof(int), 1, in );
  nwrite=fwrite(&header.nj, sizeof(int), 1, in );
  nwrite=fwrite(&header.nk, sizeof(int), 1, in );
  nwrite=fwrite(&header.nt, sizeof(int), 1, in );
  nwrite=fwrite(&header.nd, sizeof(int), 1, in );
  nwrite=fwrite(&header.code, sizeof(int), 1, in );
  nwrite=fwrite(&size, sizeof(int), 1, in );
/*   printf("header.ni=%d, header.nj=%d\n",header.ni,header.nj); */

  size=5*sizeof(float);
#if NEED_SWAP == 1
  size=lswap(size);
#endif
  nwrite=fwrite(&size, sizeof(int), 1, in );
  nwrite=fwrite(&header.xmin, sizeof(float), 1, in );
  nwrite=fwrite(&header.ymin, sizeof(float), 1, in );
  nwrite=fwrite(&header.dx,   sizeof(float), 1, in );
  nwrite=fwrite(&header.dy,   sizeof(float), 1, in );
  nwrite=fwrite(&header.spec, sizeof(float), 1, in );
  nwrite=fwrite(&size, sizeof(int), 1, in );

  size=nlevel*sizeof(float);
#if NEED_SWAP == 1
  size=lswap(size);
#endif
  nwrite=fwrite(&size, sizeof(int), 1, in );
  nitems=nlevel;
  nwrite=fwrite(header.levels, sizeof(float), nitems, in );
  nwrite=fwrite(&size, sizeof(int), 1, in );

  status=0;

  return(status);

  }

/*------------------------------------------------------------------------
  
  Nom         :  bmg_createheader

--------------------------------------------------------------------------*/

int bmg_createheader (char* file,int kmax,int dmax,int tmax,grid_t grid,float *levels,float mask)
  {
  FILE *out=NULL;
  float time;
  int i,option,status,d,t,k;
  size_t offset, nwrite, size;
  size_t tsize,dsize,ksize,shift;
  bmgheader_t header;

  if ((out = fopen(file, "wb")) == NULL) {
      __ERR_BASE_LINE__("file opening issue : %s \n", file);
      exit(-1);
    }

  offset=4;

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

  bmg_writeheader (out, header);

  fclose(out);
  return(0);
  }

/*------------------------------------------------------------------------
  
  Nom         :  bmg_read_r

--------------------------------------------------------------------------*/

int bmg_read_r (FILE *in, bmgheader_t header, int k,int d,int t, float *buf, float *time)
  {
  size_t offset=4, nread;
  size_t tsize,dsize,ksize,shift;
  bmgheader_t dum;
  int i,size;
  int status;

/*
  rewind (in);
  bmg_readheader(in,&dum,&status);
*/
 
  ksize=header.ni*header.nj*sizeof(float)+2*offset;
  dsize=header.nk*ksize;
  tsize=header.nd*dsize+sizeof(float)+2*offset;

  shift=header.size;
  fseek(in,shift,SEEK_SET);

  shift=(t-1)*tsize;
  fseek(in,shift,SEEK_CUR);

  fseek(in,offset,SEEK_CUR);
  nread=fread(time, sizeof(float), 1, in );
  if(nread !=1 ) return(-1);
  fseek(in,offset,SEEK_CUR);

  shift=(d-1)*dsize+(k-1)*ksize;
  fseek(in,shift,SEEK_CUR);

  fseek(in,offset,SEEK_CUR);
  nread=fread(buf, sizeof(float), (size_t) header.ni* header.nj, in );
  if(nread != header.ni* header.nj) return(-1);
  fseek(in,offset,SEEK_CUR);
#if NEED_SWAP == 1
  for (i=0;i<header.ni* header.nj;i++) f3swap(&buf[i]);
#endif

  return(0);
  }


/*------------------------------------------------------------------------
  
  Nom         :  bmgf_read_r

--------------------------------------------------------------------------*/

void bmgf_read_r (FILE *in, bmgheader_t *header, int *k,int *d,int *t, float *buf, float *status)
  {
  float time;

  *status=bmg_read_r (in,* header, *k, *d, *t, buf, &time);
  }

/*------------------------------------------------------------------------
  
  Nom         :  bmg_read_d

--------------------------------------------------------------------------*/

int bmg_read_d (FILE *in, bmgheader_t header, int k,int d,int t, double *buf, float *time)
  {
  size_t offset=4, nread;
  size_t tsize,dsize,ksize,shift;
  bmgheader_t dum;
  int i,size;
  int status;

/*
  rewind (in);
  bmg_readheader(in,&dum,&status);
*/
 
  ksize=header.ni*header.nj*sizeof(double)+2*offset;
  dsize=header.nk*ksize;
  tsize=header.nd*dsize+sizeof(float)+2*offset;

  shift=header.size;
  fseek(in,shift,SEEK_SET);

  shift=(t-1)*tsize;
  fseek(in,shift,SEEK_CUR);

  fseek(in,offset,SEEK_CUR);
  nread=fread(time, sizeof(float), 1, in );
  if(nread !=1 ) return(-1);
  fseek(in,offset,SEEK_CUR);

  shift=(d-1)*dsize+(k-1)*ksize;
  fseek(in,shift,SEEK_CUR);

  fseek(in,offset,SEEK_CUR);
  nread=fread(buf, sizeof(double), (size_t) header.ni* header.nj, in );
  if(nread != header.ni* header.nj) return(-1);
  fseek(in,offset,SEEK_CUR);
#if NEED_SWAP == 1
  for (i=0;i<header.ni* header.nj;i++) buf[i]=dswap2(buf[i]);
#endif

  return(0);
  }

/*------------------------------------------------------------------------
  
  Nom         :  bmg_write_r

--------------------------------------------------------------------------*/

int bmg_write_r (FILE *in, bmgheader_t header, int k,int d,int t, float *buf, float time)
  {
  size_t offset=4, nwrite;
  size_t tsize,dsize,ksize,shift;
  int size;
 
  ksize=header.ni*header.nj*sizeof(float)+2*offset;
  dsize=header.nk*ksize;
  tsize=header.nt*dsize+sizeof(float)+2*offset;

  shift=header.size;
  fseek(in,shift,SEEK_SET);

  shift=(t-1)*tsize;
  fseek(in,shift,SEEK_CUR);

  nwrite=fwrite(&offset, sizeof(int), 1, in );
  nwrite=fwrite(&time,  sizeof(float), 1, in );
  nwrite=fwrite(&offset, sizeof(int), 1, in );

  shift=(d-1)*dsize+(k-1)*ksize;
  fseek(in,shift,SEEK_CUR);

  size=header.ni*header.nj*sizeof(float);
  nwrite=fwrite(&size, sizeof(int), 1, in );
  nwrite=fwrite(buf, sizeof(float), (size_t) header.ni* header.nj, in );
  nwrite=fwrite(&size, sizeof(int), 1, in );
  return(0);

  }

/*------------------------------------------------------------------------

  Nom         :  bmg_loadgrid

--------------------------------------------------------------------------*/
int bmg_loadgrid (const char* filename,grid_t *grid)
  {
  FILE *in=NULL;
  int option,status;
  size_t offset, nread, size;
  float time;
  bmgheader_t header;

  status=0;

  if ((in = fopen(filename, "rb")) == NULL) {
      __ERR_BASE_LINE__("file opening issue : %s \n", filename);
      exit(-1);
    }

  status=bmg_readheader(in,&header);
  if (status !=0) {
    status=-1;
    return(status);
    }

  fclose(in);

  grid->xmin=header.xmin;
  grid->ymin=header.ymin;

  grid->dx=header.dx;
  grid->dy=header.dy;

  grid->nx=header.ni;
  grid->ny=header.nj;
  grid->nz=header.nk;

  grid->xmax=grid->xmin+(grid->nx-1)*grid->dx;
  grid->ymax=grid->ymin+(grid->ny-1)*grid->dy;

  return(status);
  }

/*------------------------------------------------------------------------

  Nom         :  bmg_getinfo

--------------------------------------------------------------------------*/

int bmg_getinfo (char* filename, bmgheader_t *header)
  {
  FILE *in=NULL;
  int option,status,s;
  size_t offset, nread, size;
  float time;

  status=0;

  if ((in = fopen(filename, "rb")) == NULL) {
      __ERR_BASE_LINE__("file opening issue : %s \n", filename);
      exit(-1);
    }

  status=bmg_readheader(in,header);
  if (status !=0) {
    status=-1;
    return(status);
    }

  fclose(in);

  s=(filesize(filename)-header->size);
  s=s/header->nt-(sizeof(float) /*time*/+2*4/*rec size*/);
  s=s/header->nd;
  s=s/header->nk-2*4/*rec size*/;
  s=s/(header->ni*header->nj);
  s=s/sizeof(float);
  header->nv=s;

  status=0;

  return(status);

  }

/*------------------------------------------------------------------------
  
  Nom         :  bmg_h2g

--------------------------------------------------------------------------*/

void bmg_h2g (bmgheader_t header, grid_t *grid)
  {
  int status;

  grid->xmin=header.xmin;
  grid->ymin=header.ymin;

  grid->dx=header.dx;
  grid->dy=header.dy;

  grid->nx=header.ni;
  grid->ny=header.nj;

  grid->xmax=grid->xmin+(grid->nx-1)*grid->dx;
  grid->ymax=grid->ymin+(grid->ny-1)*grid->dy;

  status=0;

  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void bmg_checkfile (char *filename, int *status)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  {
  grid_t grid;

  *status=bmg_loadgrid (filename, &grid);
  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int bmg_loadr1 (const char* file,int k,int d,int t,grid_t grid,float *buf,float *mask,float *time)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  {
  FILE *in=NULL;
  int i,option,status;
  size_t offset, nread, size;
  size_t tsize,dsize,ksize;
  long   shift;
  bmgheader_t header;

  in = fopen(file, "rb");
  if (in == NULL) {
    __ERR_BASE_LINE__("file opening issue : %s \n", file);
    exit(-1);
    }

  offset=4;

  status=bmg_readheader(in,&header);

  *mask= header.spec;

  ksize=header.ni*header.nj*sizeof(float)+2*offset;
  dsize=header.nk*ksize;
  tsize=header.nd*dsize+sizeof(float)+2*offset;

/* *------------------------------------------------------------------------
  skip first time frames*/
  shift=(t-1)*tsize;
  fseek(in,shift,SEEK_CUR);

  fseek(in,offset,SEEK_CUR);
  nread=fread(time, sizeof(float), 1, in );
  fseek(in,offset,SEEK_CUR);

/* *------------------------------------------------------------------------
  skip first dimension and layer frames*/
  shift=(d-1)*dsize+(k-1)*ksize;
  fseek(in,shift,SEEK_CUR);

#if NEED_SWAP == 1
/*   *time= f3swap(*time); */
  f3swap(time);
#endif

  fseek(in,offset,SEEK_CUR);
  nread=fread(buf, sizeof(float), (size_t) header.ni* header.nj, in );
  fseek(in,offset,SEEK_CUR);

#if NEED_SWAP == 1
  for (i=0;i<header.ni* header.nj;i++) f3swap(&buf[i]);
/*   for (i=0;i<header.ni* header.nj;i++) buf[i]=f3swap(buf[i]); */
#endif

  fclose(in);

  return(status);
  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int bmg_loadr1 (const char* file, int k,int d,int t,grid_t grid,short *buf,short *mask,float *time)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  {
  int m,status;
  float *rbuf,rmask;
  
  rbuf=new float[grid.nx*grid.ny];
  status=bmg_loadr1 (file, k, d, t, grid, rbuf, &rmask, time);
  for(m=0;m<grid.nx*grid.ny;m++) buf[m]=rbuf[m];
  *mask=rmask;
  delete[] rbuf;
  return(status);
  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int bmgf_loadr1 (char* file, int *k,int *d,int *t,grid_t *grid,float *buf,float *mask,float *time)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  {
  int status;

  status=bmg_loadr1 (file, *k, *d, *t, *grid, buf, mask, time);
  return(status);
  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int bmg_save_r (char* file, int t, grid_t grid, float *buf, float time)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  {
  FILE *out=NULL;
  int i,option,status,d,k,dmax,kmax,tmax, size;
  size_t offset, nwrite, buf_size;
  size_t tsize,dsize,ksize,shift;
  bmgheader_t header;

  if ((out = fopen(file, "rb+wb")) == NULL) {
      __ERR_BASE_LINE__("file opening issue : %s \n", file);
      exit(-1);
    }

  offset=4;

  bmg_readheader (out, &header);
/*   printf("header.ni=%d, header.nj=%d\n",header.ni,header.nj); */
  printf("header.nk=%d, header.nd=%d, header.nt=%d\n",
          header.nk,header.nd,header.nt);

  kmax=header.nk;
  dmax=header.nd;
  tmax=header.nt;

  ksize=header.ni*header.nj*sizeof(float)+2*offset;
  dsize=header.nk*ksize;
  tsize=header.nt*dsize+sizeof(float)+2*offset;

  shift=(t-1)*tsize;
  fseek(out,shift,SEEK_CUR);

#if NEED_SWAP == 1
/*   time= f3swap(time); */
  f3swap(&time);
#endif

  size=sizeof(float);
#if NEED_SWAP == 1
  size=lswap(size);
#endif

  nwrite=fwrite(&size, sizeof(int),   1, out );
  nwrite=fwrite(&time, sizeof(float), 1, out );
  nwrite=fwrite(&size, sizeof(int),   1, out );

  buf_size=header.ni* header.nj*kmax*dmax;

/*   shift=(d-1)*dsize+(k-1)*ksize; */
/*   fseek(in,shift,SEEK_CUR);   */

#if NEED_SWAP == 1
/*   for (i=0;i<buf_size;i++) buf[i]=f3swap(buf[i]); */
  for (i=0;i<buf_size;i++) f3swap(&buf[i]);
#endif

  size=header.ni*header.nj*sizeof(float);
#if NEED_SWAP == 1
  size=lswap(size);
#endif
  for(d=0;d<dmax;d++) {
    for(k=0;k<kmax;k++) {
      nwrite=fwrite(&size, sizeof(int), 1, out );
      nwrite=fwrite(&buf[(k+d*kmax)*header.ni*header.nj],sizeof(float),(size_t)header.ni*header.nj,out);
      nwrite=fwrite(&size, sizeof(int), 1, out );
      }
    }

#if NEED_SWAP == 1
/*   for (i=0;i<buf_size;i++) buf[i]=f3swap(buf[i]); */
  for (i=0;i<buf_size;i++) f3swap(&buf[i]);
#endif

  fclose(out);
  return(0);
  }

/*------------------------------------------------------------------------

  Nom         :  bmgf_save_r

--------------------------------------------------------------------------*/
int bmgf_save_r(char* file, int *t,grid_t *grid,float *buf,float *time)
  {
  int status;
  status=bmg_save_r (file,*t,*grid,buf,*time);
  return(status);
  }

/*------------------------------------------------------------------------

  Nom     : bmg_saver1


  Comments: bmg_saver1 save an unique time frame,dimension,level

--------------------------------------------------------------------------*/

int bmg_saver1 (char* file, int k, int d, int t, grid_t grid, float *buf, float time,float mask)
  {
  FILE *out=NULL;
  int i,option,status,dmax,kmax,tmax;
  size_t offset, nwrite, size, buf_size;
  size_t tsize,dsize,ksize,shift;
  bmgheader_t header;
  float levels[1]={0.};

  out=fopen(file,"r+");

  if(out==NULL) {
    kmax=1;
    dmax=1;
    tmax=t;
    status=bmg_createheader (file, kmax, dmax, tmax, grid, levels, mask);
    out=fopen(file,"r+");
    bmg_readheader(out, &header);
    }
  else {
    bmg_readheader(out, &header);
    header.nt=MAX(t,header.nt);
    status=bmg_writeheader(out, header);
    }

  offset=4;

  kmax=header.nk;
  dmax=header.nd;
  tmax=header.nt;

  ksize=header.ni*header.nj*sizeof(float)+2*offset;
  dsize=header.nk*ksize;
  tsize=header.nd*dsize+sizeof(float)+2*offset;

  shift=(t-1)*tsize;
  fseek(out,shift,SEEK_CUR);

  buf_size=header.ni*header.nj;

  shift=(d-1)*dsize+(k-1)*ksize;
  fseek(out,shift,SEEK_CUR);

  size=4;
#if NEED_SWAP == 1
  size=lswap(size);
#endif
  nwrite=fwrite(&size, sizeof(int), 1, out );
#if NEED_SWAP == 1
  f3swap(&time);
#endif
  nwrite=fwrite(&time,sizeof(float),1,out);
  nwrite=fwrite(&size, sizeof(int), 1, out );

#if NEED_SWAP == 1
/*   for (i=0;i<buf_size;i++) buf[i]=f3swap(buf[i]); */
  for (i=0;i<buf_size;i++) f3swap(&buf[i]);
#endif

  size=header.ni*header.nj*sizeof(float);
#if NEED_SWAP == 1
  size=lswap(size);
#endif
  nwrite=fwrite(&size, sizeof(int), 1, out );
  nwrite=fwrite(buf,sizeof(float),(size_t)header.ni*header.nj,out);
  nwrite=fwrite(&size, sizeof(int), 1, out );

#if NEED_SWAP == 1
/*   for (i=0;i<buf_size;i++) buf[i]=f3swap(buf[i]); */
  for (i=0;i<buf_size;i++) f3swap(&buf[i]);
#endif

  fclose(out);

  return(0);

  }

/*------------------------------------------------------------------------

  Nom         :  bmgf_saver1

--------------------------------------------------------------------------*/

int bmgf_saver1(char* file, int *k,int *d,int *t,grid_t *grid,float *buf,float *mask,float *time)
  {
  int status;
  status=bmg_saver1 (file,*k,*d,*t,*grid,buf,*time,*mask);
  return(status);
  }

/*------------------------------------------------------------------------

  Nom     : bmg_saver2


  Comments: bmg_saver1 save an unique time frame,dimension,level

--------------------------------------------------------------------------*/

int bmg_saver2 (char* file, int k, int d, int t, grid_t grid, float *bufx, float *bufy, float time, float mask)
  {

  FILE *out=NULL;
  int i,option,status,dmax,kmax,tmax;
  size_t offset, nwrite, size, buf_size;
  size_t tsize,dsize,ksize,shift;
  bmgheader_t header;
  float levels[1]={0.};

  out=fopen(file,"r+");

  if(out==NULL) {
    kmax=1;
    dmax=2;
    tmax=t;
    status=bmg_createheader (file, kmax, dmax, tmax, grid, levels, mask);
    out=fopen(file,"r+");
    bmg_readheader(out, &header);
    }
  else
    {
    bmg_readheader(out, &header);
    header.nt=MAX(t,header.nt);
    status=bmg_writeheader(out, header);
    }

  offset=4;

  kmax=header.nk;
  dmax=header.nd;
  tmax=header.nt;

  ksize=header.ni*header.nj*sizeof(float)+2*offset;
  dsize=header.nk*ksize;
  tsize=header.nd*dsize+sizeof(float)+2*offset;

  shift=(t-1)*tsize;
  fseek(out,shift,SEEK_CUR);
  
  buf_size=header.ni*header.nj;

  shift=(d-1)*dsize+(k-1)*ksize;
  fseek(out,shift,SEEK_CUR);

  size=4;
#if NEED_SWAP == 1
  size=lswap(size);
#endif
  nwrite=fwrite(&size, sizeof(int), 1, out );
#if NEED_SWAP == 1
  f3swap(&time);
#endif
  nwrite=fwrite(&time,sizeof(float),1,out);
  nwrite=fwrite(&size, sizeof(int), 1, out );

#if NEED_SWAP == 1
  for (i=0;i<buf_size;i++) f3swap(&bufx[i]);
  for (i=0;i<buf_size;i++) f3swap(&bufy[i]);
#endif

  size=header.ni*header.nj*sizeof(float);
#if NEED_SWAP == 1
  size=lswap(size);
#endif
  nwrite=fwrite(&size, sizeof(int), 1, out );
  nwrite=fwrite(bufx,sizeof(float),(size_t)header.ni*header.nj,out);
  nwrite=fwrite(&size, sizeof(int), 1, out );
  nwrite=fwrite(&size, sizeof(int), 1, out );
  nwrite=fwrite(bufy,sizeof(float),(size_t)header.ni*header.nj,out);
  nwrite=fwrite(&size, sizeof(int), 1, out );

#if NEED_SWAP == 1
  for (i=0;i<buf_size;i++) f3swap(&bufx[i]);
  for (i=0;i<buf_size;i++) f3swap(&bufy[i]);
#endif

  fclose(out);

  return(0);

  }

/*------------------------------------------------------------------------
  
  Nom         :  bmg_loadc1

--------------------------------------------------------------------------*/

int bmg_loadc1(const char* file,int k,int d,int t,grid_t *grid,fcomplex **buffer,fcomplex *mask,float *time)
  {
  
  FILE *in=NULL;
  int option,i,status;
  size_t offset, nread, size;
  size_t tsize,dsize,ksize;
  long shift;
  double *dum=NULL;
  bmgheader_t header;

  if ((in = fopen(file, "rb")) == NULL) {
      __ERR_BASE_LINE__("file opening issue : %s \n", file);
      exit(-1);
    }
  
  offset=4;

  status=bmg_readheader(in,&header);

  grid->xmin=header.xmin;
  grid->ymin=header.ymin;

  grid->dx=header.dx;
  grid->dy=header.dy;

  grid->nx=header.ni;
  grid->ny=header.nj;

  grid->xmax=grid->xmin+(grid->nx-1)*grid->dx;
  grid->ymax=grid->ymin+(grid->ny-1)*grid->dy;

  grid->modeH=0;

  exitIfNull(
    grid->x=(double *) malloc(grid->nx*sizeof(double))
    );
  for (i=0;i<grid->nx;i++) grid->x[i]=grid->xmin+i*grid->dx;

  exitIfNull(
    grid->y=(double *) malloc(grid->ny*sizeof(double))
    );
  for (i=0;i<grid->ny;i++) grid->y[i]=grid->ymin+i*grid->dy;

  ksize=header.ni*header.nj*sizeof(float)+2*offset;
  dsize=header.nk*ksize;
  tsize=header.nd*dsize+sizeof(float)+2*offset;

  shift=(t-1)*tsize;
  fseek(in,shift,SEEK_CUR);

  fseek(in,offset,SEEK_CUR);
  nread=fread(time, sizeof(float), 1, in );
  fseek(in,offset,SEEK_CUR);

  shift=(d-1)*dsize+(k-1)*ksize;
  fseek(in,shift,SEEK_CUR);
/*   fseek(in,offset,SEEK_CUR);   */
/*   nread=fread(time, sizeof(float), 1, in ); */
/*   fseek(in,offset,SEEK_CUR); */

#if NEED_SWAP == 1
/*   *time= f3swap(*time); */
  f3swap(time);
#endif

  if(!*buffer) {
    exitIfNull(
      *buffer=(fcomplex *) malloc(grid->nx*grid->ny*sizeof(fcomplex))
      );
    }
  fseek(in,offset,SEEK_CUR);
  nread=fread(*buffer, sizeof(fcomplex), (size_t) header.ni* header.nj, in );
  fseek(in,offset,SEEK_CUR);

#if NEED_SWAP == 1
  for (i=0;i<header.ni* header.nj;i++)
    (*buffer)[i]=cswap((*buffer)[i]);
#endif

  fclose(in);

/*
  mask->r= header.spec;
  mask->i= header.spec;
*/
  *mask= fcomplex(header.spec,header.spec);

  for (int m=0;m<grid->Hsize();m++) {
    if((*buffer)[m]==*mask) (*buffer)[m]=fcomplex(9999.,9999.);
    }
    
  *mask= fcomplex(9999.,9999.);
  return(0);

  }

/*------------------------------------------------------------------------
  
  Nom         :  bmgf_loadc1

--------------------------------------------------------------------------*/

int bmgf_loadc1(char* file,int *k,int *d,int *t,grid_t *grid,fcomplex *buf,fcomplex *mask,float *time)
  {
  int status;
  grid_t dum;
  status=bmg_loadc1(file,*k,*d,*t,&dum,&buf,mask,time);
  return(status);
  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int bmg_savec1 (char* filename,int k,int t,int d,grid_t grid,fcomplex *buf,float time, fcomplex mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  {
  
  FILE *out=NULL;
  int i,option,status,size;
  size_t offset, nwrite;
  size_t tsize,dsize,ksize,shift;
  bmgheader_t header;
  float levels[1]={0.};

  printf("#bmg_savec1 ... starting\n");

  if ((out = fopen(filename, "wb")) == NULL) {
    __ERR_BASE_LINE__("file opening issue : %s \n", filename);
    exit(-1);
    }

  offset=4;

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
/*
  header.spec=mask.r;
*/
  header.spec=real(mask);
  header.levels=levels;

  bmg_writeheader (out, header);

  ksize=header.ni*header.nj*sizeof(fcomplex)+2*offset;
  dsize=header.nk*ksize;
  tsize=header.nt*dsize+sizeof(float)+2*offset;

  shift=(t-1)*tsize;
  fseek(out,shift,SEEK_CUR);
  
  shift=(d-1)*dsize+(k-1)*ksize;
  fseek(out,shift,SEEK_CUR);

  size=4;
#if NEED_SWAP == 1
  size=lswap(size);
#endif
  nwrite=fwrite(&size, sizeof(int), 1, out );
#if NEED_SWAP == 1
  f3swap(&time);
#endif
  nwrite=fwrite(&time,sizeof(float),1,out);
  nwrite=fwrite(&size, sizeof(int), 1, out );

  size=header.ni*header.nj*sizeof(fcomplex);
#if NEED_SWAP == 1
  size=lswap(size);
  for (i=0;i<header.ni* header.nj;i++) buf[i]=cswap(buf[i]);
#endif

  nwrite=fwrite(&size, sizeof(int), 1, out );
  nwrite=fwrite(buf, sizeof(fcomplex), (size_t) header.ni* header.nj, out );
  nwrite=fwrite(&size, sizeof(int), 1, out );

#if NEED_SWAP == 1
  for (i=0;i<header.ni* header.nj;i++) buf[i]=cswap(buf[i]);
#endif

  fclose(out);
  printf("#bmg_savec1 ... finished\n");
  return(0);

  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int bmg_savec2 (char* filename,int k,int t,int d,grid_t grid,fcomplex *bufx,fcomplex *bufy,float time,fcomplex mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  
  FILE *out=NULL;
  int i,option,status,size;
  size_t offset, nwrite;
  size_t tsize,dsize,ksize,shift;
  bmgheader_t header;
  float levels[1]={0.};

  printf("#bmg_savec2 ... starting\n");

  if ((out = fopen(filename, "wb")) == NULL) {
      __ERR_BASE_LINE__("file opening issue : %s \n", filename);
      exit(-1);
    }
  
  offset=4;

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
  header.nd=2;
  header.nt=1;
/*
  header.spec=mask.r;
*/
  header.spec=real(mask);
  header.levels=levels;

  bmg_writeheader (out, header);

  ksize=header.ni*header.nj*sizeof(fcomplex)+2*offset;
  dsize=header.nk*ksize;
  tsize=header.nt*dsize+sizeof(float)+2*offset;

  shift=(t-1)*tsize;
  fseek(out,shift,SEEK_CUR);

  shift=(d-1)*dsize+(k-1)*ksize;
  fseek(out,shift,SEEK_CUR);

  size=4;
#if NEED_SWAP == 1
  size=lswap(size);
#endif
  nwrite=fwrite(&size, sizeof(int), 1, out );
#if NEED_SWAP == 1
  f3swap(&time);
#endif
  nwrite=fwrite(&time,sizeof(float),1,out);
  nwrite=fwrite(&size, sizeof(int), 1, out );

  size=header.ni*header.nj*sizeof(fcomplex);
#if NEED_SWAP == 1
  size=lswap(size);
  for (i=0;i<header.ni* header.nj;i++) bufx[i]=cswap(bufx[i]);
  for (i=0;i<header.ni* header.nj;i++) bufy[i]=cswap(bufy[i]);
#endif

  nwrite=fwrite(&size, sizeof(int), 1, out );
  nwrite=fwrite(bufx, sizeof(fcomplex), (size_t) header.ni* header.nj, out );
  nwrite=fwrite(&size, sizeof(int), 1, out );

  nwrite=fwrite(&size, sizeof(int), 1, out );
  nwrite=fwrite(bufy, sizeof(fcomplex), (size_t) header.ni* header.nj, out );
  nwrite=fwrite(&size, sizeof(int), 1, out );

#if NEED_SWAP == 1
  for (i=0;i<header.ni* header.nj;i++) bufx[i]=cswap(bufx[i]);
  for (i=0;i<header.ni* header.nj;i++) bufy[i]=cswap(bufy[i]);
#endif

  fclose(out);
  printf("#bmg_savec2 ... finished\n");
  return(0);

  }

