#include "config.h"
   
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
 
#include "tools-structures.h"

#include "fe.def"
#include "fe.h"
#include "map.h"
#include "rutin.h"
#include "geo.h"
#include "poc-time.h"
#include "bmg.h"
#include "archive.h"
#include "swap.h"


/*------------------------------------------------------------------------
  
  Nom         :  fe_filemapr1

--------------------------------------------------------------------------*/

int fe_filemapr1 (char* input,char* output,grid_t grid,int minstep,int maxstep,meta_archive_t info)
  {
  
  FILE *out;
  double time;
  float ftime;
  int i,option,status,d,t,k;
  int dmax=1,tmax,kmax=1, size, shift;
  size_t offset, nwrite,nitems;
  size_t tsize,dsize,ksize;
  bmgheader_t header;
  float *levels,mask=999.9,*buf;
  int *element,nnde;
  date_t date=info.reference;

  status=map_completegridaxis_2(&grid);

  element=fe_scan_elements(info.mesh,grid,0);
  if (element == NULL) return(99);
  //if (element == NULL) return; //modif thierry le 21/10/2005

  buf=(float *)malloc(grid.nx*grid.ny*sizeof(float));
  if (buf == NULL) return(99);
  //if (buf == NULL) return; //modif thierry le 21/10/2005

  out=fopen(output,"wb");
  
  offset=4;

  strcpy(header.comment[0],"mog2d simulation");
  strcpy(header.comment[1],"no comment");
  strcpy(header.comment[2],input);
  strcpy(header.comment[3],"no comment");
  header.ni=grid.nx;
  header.nj=grid.ny;
  header.xmin=grid.xmin;
  header.ymin=grid.ymin;
  header.dx=grid.dx;
  header.dy=grid.dy;
  header.nk=kmax;
  header.nd=dmax;
  header.nt=maxstep-minstep+1;
  header.spec=mask;
  levels=(float *)malloc(kmax*sizeof(float));
  header.levels=levels;

  status=bmg_writeheader (out, header);
  if(status!=0) goto error;

  ksize=header.ni*header.nj*sizeof(float)+2*offset;
  dsize=header.nk*ksize;
  tsize=header.nt*dsize+sizeof(float)+2*offset;

  shift=mjd(date.month,date.day,date.year)-mjd(1,1,1985);

  for(t=minstep;t<=maxstep;t++) {
    status=fe_framemapr1(input, info, grid, mask, element, t, buf, &time);
    if(status!=0) goto error;
    size=sizeof(float);
    ftime=time/(3600.*24.)+shift;
#if NEED_SWAP == 1
/*     time= f3swap(ftime); */
    f3swap(&ftime);
    size=lswap(size);
#endif
    nitems=1;
    nwrite=fwrite(&size, sizeof(int), 1, out );
    if(nwrite!=nitems) goto error;
    nitems=1;
    nwrite=fwrite(&ftime, sizeof(float), 1, out );
    if(nwrite!=nitems) goto error;
    nitems=1;
    nwrite=fwrite(&size, sizeof(int), 1, out );
    if(nwrite!=nitems) goto error;
    size=header.ni*header.nj*sizeof(float);
#if NEED_SWAP == 1
    size=lswap(size);
/*     for (i=0;i<header.ni* header.nj*kmax;i++) buf[i]=f3swap(buf[i]); */
    for (i=0;i<header.ni* header.nj*kmax;i++) f3swap(&buf[i]);
#endif
    for(k=0;k<kmax;k++) {
      nitems=1;
      nwrite=fwrite(&size, sizeof(int), 1, out );
      if(nwrite!=nitems) goto error;
      nitems=(size_t) header.ni* header.nj;
      nwrite=fwrite(&buf[k*header.ni*header.nj], sizeof(float), nitems, out );
      if(nwrite!=nitems) goto error;
      nitems=1;
      nwrite=fwrite(&size, sizeof(int), 1, out );
      if(nwrite!=nitems) goto error;
      }
#if NEED_SWAP == 1
/*  for (i=0;i<header.ni* header.nj*kmax;i++) buf[i]=f3swap(buf[i]); */
    for (i=0;i<header.ni* header.nj*kmax;i++) f3swap(&buf[i]);
#endif
    }

  fclose(out);
  free(element);
  free(buf);
  return(0);

error:
  fclose(out);
  free(element);
  free(buf);
  return(-1);
  }
  
/*------------------------------------------------------------------------
  
  Nom         :  fe_filemapr2

--------------------------------------------------------------------------*/

int fe_filemapr2 (char* input,char* output,grid_t grid,int minstep,int maxstep,meta_archive_t info)
  {
  
  FILE *out;
  double time;
  float ftime;
  int i,option,status,d,t,k;
  int dmax=2,tmax,kmax=1, size, shift;
  size_t offset, nwrite,nitems;
  size_t tsize,dsize,ksize;
  bmgheader_t header;
  float *levels,mask=999.9,*bufx,*bufy;
  int *element, tab[3],ntab;
  date_t date=info.reference;

  status=map_completegridaxis_2(&grid);

  element=fe_scan_elements(info.mesh,grid,0);
  if (element == NULL) return(99);
  //if (element == NULL) return; //modif thierry le 21/10/2005

  bufx=(float *)malloc(grid.nx*grid.ny*sizeof(float));
  if (bufx == NULL) return(99);
  //if (bufx == NULL) return; //modif thierry le 21/10/2005

  bufy=(float *)malloc(grid.nx*grid.ny*sizeof(float));
  if (bufy == NULL) return(99);
  //if (bufy == NULL) return; //modif thierry le 21/10/2005

  out=fopen(output,"wb");
  
  offset=4;

  strcpy(header.comment[0],"mog2d simulation");
  strcpy(header.comment[1],"no comment");
  strcpy(header.comment[2],input);
  strcpy(header.comment[3],"no comment");
  header.ni=grid.nx;
  header.nj=grid.ny;
  header.xmin=grid.xmin;
  header.ymin=grid.ymin;
  header.dx=grid.dx;
  header.dy=grid.dy;
  header.nk=kmax;
  header.nd=dmax;
  header.nt=maxstep-minstep+1;
  header.spec=mask;
  levels=(float *)malloc(kmax*sizeof(float));
  header.levels=levels;

  status=bmg_writeheader (out, header);
  if(status!=0) goto error;

  ksize=header.ni*header.nj*sizeof(float)+2*offset;
  dsize=header.nk*ksize;
  tsize=header.nt*dsize+sizeof(float)+2*offset;

  shift=mjd(date.month,date.day,date.year)-mjd(1,1,1985);

  for(t=minstep;t<=maxstep;t++) {
    status=fe_framemapr2(input, info, grid, mask, element, t, bufx, bufy, &time);
    if(status!=0) goto error;
    size=sizeof(float);
    ftime=time/(3600.*24.)+shift;
#if NEED_SWAP == 1
/*     time= f3swap(ftime); */
    f3swap(&ftime);
    size=lswap(size);
#endif
    nitems=1;
    nwrite=fwrite(&size, sizeof(int), 1, out );
    if(nwrite!=nitems) goto error;
    nitems=1;
    nwrite=fwrite(&ftime, sizeof(float), 1, out );
    if(nwrite!=nitems) goto error;
    nitems=1;
    nwrite=fwrite(&size, sizeof(int), 1, out );
    if(nwrite!=nitems) goto error;
    size=header.ni*header.nj*sizeof(float);
#if NEED_SWAP == 1
    size=lswap(size);
/*     for (i=0;i<header.ni* header.nj*kmax;i++) buf[i]=f3swap(buf[i]); */
    for (i=0;i<header.ni* header.nj*kmax;i++) f3swap(&bufx[i]);
    for (i=0;i<header.ni* header.nj*kmax;i++) f3swap(&bufy[i]);
#endif
    for(k=0;k<kmax;k++) {
      nitems=1;
      nwrite=fwrite(&size, sizeof(int), 1, out );
      if(nwrite!=nitems) goto error;
      nitems=(size_t) header.ni* header.nj;
      nwrite=fwrite(&bufx[k*header.ni*header.nj], sizeof(float), nitems, out );
      if(nwrite!=nitems) goto error;
      nitems=1;
      nwrite=fwrite(&size, sizeof(int), 1, out );
      if(nwrite!=nitems) goto error;
      nitems=1;
      nwrite=fwrite(&size, sizeof(int), 1, out );
      if(nwrite!=nitems) goto error;
      nitems=(size_t) header.ni* header.nj;
      nwrite=fwrite(&bufy[k*header.ni*header.nj], sizeof(float), nitems, out );
      if(nwrite!=nitems) goto error;
      nitems=1;
      nwrite=fwrite(&size, sizeof(int), 1, out );
      if(nwrite!=nitems) goto error;
      }
    }


  fclose(out);
  free(element);
  free(bufx);
  free(bufy);
  return(0);

error:
  fclose(out);
  free(element);
  free(bufx);
  free(bufy);
  return(-1);
  }
  
