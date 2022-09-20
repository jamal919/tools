
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

#include "config.h"

#include <stdio.h>
#include <string.h>

#include "tools-define.h"
#include "tools-structures.h"

#include "functions.h"
#include "fe.h"
#include "bmg.h"
#include "grd.h"
#include "map.h"
#include "ascii.h"
#include "archive.h"
#include "geo.h"
#include "sym-io.h"


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int cross_reference_no_noptimal(mesh_t mesh, double *x, double *y, int *translation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int   m,n;
  double x1,y1,x2,y2,d,dx,dy,*r;

  r=new double[mesh.nvtxs];

/* *----------------------------------------------------------------------------
  return node index in mesh arrangement*/
  for(n=0;n<mesh.nvtxs;n++) {
    translation[n]=-1;
    }

  for(n=0;n<mesh.nvtxs;n++) {
    x1=x[n];
    y1=y[n];
    r[n]=1.e+10;
    for(m=n;m<mesh.nvtxs;m++) {
      x2=mesh.vertices[m].lon;
      dx=(x2-x1)*(x2-x1);
      if(dx>1.e-06) continue;
      y2=mesh.vertices[m].lat;
      dy=(y2-y1)*(y2-y1);
      if(dy>1.e-06) continue;
      d=dx+dy;
      if(d<r[n]) {
        r[n]=d;
        translation[n]=m;
        }
      if(r[n]<1.e-09) break;
      }
    }

  for(n=0;n<mesh.nvtxs;n++) {
    x1=x[n];
    y1=y[n];
    if(r[n]<1.e-09) break;
    for(m=0;m<n;m++) {
      x2=mesh.vertices[m].lon;
      y2=mesh.vertices[m].lat;
      dx=(x2-x1)*(x2-x1);
      if(dx>1.e-06) continue;
      dy=(y2-y1)*(y2-y1);
      if(dy>1.e-06) continue;
      d=dx+dy;
      if(d<r[n]) {
        r[n]=d;
        translation[n]=m;
        }
      if(r[n]<1.e-09) break;
      }
    }

  for(n=0;n<mesh.nvtxs;n++) {
    if(r[n]>1.e-03) goto error;
    }

  for(n=0;n<mesh.nvtxs;n++) {
    if(translation[n]==-1) goto error;
    }

   for(n=0;n<mesh.nvtxs;n++) {
    for(m=n+1;m<mesh.nvtxs;m++) {
      if(translation[n]==translation[m]) {
        printf("research failed...\n");
        goto error;
        }
      }
    }

  delete[] r;
  return (0);

 error:
  return (-1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  frame_t *partition_frame(frame_t frame,int partition)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int  i,j,count,n;
  double dx,dy;
  frame_t *boxes;

  n=partition*partition;

  boxes=new frame_t[n];

  dx=(frame.xmax-frame.xmin)/(double) partition;
  dy=(frame.ymax-frame.ymin)/(double) partition;

  for(j=0;j<partition;j++) {
    for(i=0;i<partition;i++) {
      boxes[count].xmin=frame.xmin+i*dx;
      boxes[count].xmax=frame.xmin+i*dx+dx;
      boxes[count].ymin=frame.ymin+j*dy;
      boxes[count].ymax=frame.ymin+j*dy+dy;
      }
    }

  return(boxes);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int cross_reference(mesh_t mesh, double *x, double *y, int *translation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int   i,j,ii,jj,l,m,mm,n,status;
  int   imin,imax,jmin,jmax;
  double x1,y1,x2,y2,d,dx,dy,ddx,ddy,*r;
  int    **list,*card,nlist,partition=20;
  int    *check;
  frame_t frame,*boxes;

  r=new double[mesh.nvtxs];

  check=new int [mesh.nvtxs];
  for(n=0;n<mesh.nvtxs;n++) {
    check[n]=-1;
    }

/* *----------------------------------------------------------------------------
  return node index in mesh arrangement*/
  for(n=0;n<mesh.nvtxs;n++) {
    translation[n]=-1;
    }

  status=fe_minmax(mesh, frame);
  dx=(frame.xmax-frame.xmin)/(double) partition;
  dy=(frame.ymax-frame.ymin)/(double) partition;

  frame.xmin=frame.xmin-0.1*dx;
  frame.xmax=frame.xmax+0.1*dx;
  frame.ymin=frame.ymin-0.1*dy;
  frame.ymax=frame.ymax+0.1*dy;

  dx=(frame.xmax-frame.xmin)/(double) partition;
  dy=(frame.ymax-frame.ymin)/(double) partition;

//  boxes=partition_frame(frame,partition);

  nlist=partition*partition;

  card=new int[nlist];
  for(l=0;l<nlist;l++) card[l]=0;

  for(n=0;n<mesh.nvtxs;n++) {
    x1=mesh.vertices[n].lon;
    y1=mesh.vertices[n].lat;
    i=(int) floor((x1-frame.xmin)/dx);
    j=(int) floor((y1-frame.ymin)/dy);
    imin=(int) floor((x1-frame.xmin)/dx-0.1);
    imax=(int) floor((x1-frame.xmin)/dx+0.1);
    jmin=(int) floor((y1-frame.ymin)/dy-0.1);
    jmax=(int) floor((y1-frame.ymin)/dy+0.1);
    for (jj=MAX(0,jmin);jj<=MIN(jmax,partition-1);jj++) {
      for (ii=MAX(0,imin);ii<=MIN(imax,partition-1);ii++) {
        l=jj*partition+ii;
        card[l]++;
        }
      }
    }

  list=new int*[nlist];
  for(l=0;l<nlist;l++) list[l]=new int[card[l]];

  for(l=0;l<nlist;l++) card[l]=0;
  for(n=0;n<mesh.nvtxs;n++) {
    x1=mesh.vertices[n].lon;
    y1=mesh.vertices[n].lat;
    i=(int) floor((x1-frame.xmin)/dx);
    j=(int) floor((y1-frame.ymin)/dy);
    imin=(int) floor((x1-frame.xmin)/dx-0.1);
    imax=(int) floor((x1-frame.xmin)/dx+0.1);
    jmin=(int) floor((y1-frame.ymin)/dy-0.1);
    jmax=(int) floor((y1-frame.ymin)/dy+0.1);
    for (jj=MAX(0,jmin);jj<=MIN(jmax,partition-1);jj++) {
      for (ii=MAX(0,imin);ii<=MIN(imax,partition-1);ii++) {
        l=jj*partition+ii;
        list[l][card[l]]=n;
        card[l]++;
        }
      }
    }

  for(n=0;n<mesh.nvtxs;n++) {
    x1=x[n];
    y1=y[n];
    i=(int) floor((x1-frame.xmin)/dx);
    j=(int) floor((y1-frame.ymin)/dy);
    l=j*partition+i;
    r[n]=1.e+10;
    for(mm=0;mm<card[l];mm++) {
      m=list[l][mm];
      x2=mesh.vertices[m].lon;
      y2=mesh.vertices[m].lat;
      ddx=(x2-x1)*(x2-x1);
      if(ddx>1.e-03) continue;
      ddy=(y2-y1)*(y2-y1);
      if(ddy>1.e-03) continue;
      d=ddx+ddy;
      if(d<r[n]) {
        r[n]=d;
        translation[n]=m;
        }
      if(r[n]<1.e-09) break;
      }
    if(check[translation[n]]==-1) {
      check[translation[n]]=n;
      }
    else{
      printf("research failed...\n");
      goto error;
      }
    }

  for(n=0;n<mesh.nvtxs;n++) {
    if(r[n]>1.e-03) goto error;
    }

  for(n=0;n<mesh.nvtxs;n++) {
    if(translation[n]==-1) goto error;
    }

  delete[] card;
  for(l=0;l<nlist;l++) delete[] list[l];
  delete list;

  delete[] r;
  delete[] check;
  return (0);

 error:
  return (-1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int cross_reference_inverse(mesh_t mesh, double *x, double *y, int *translation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int   m,n;
  double x1,y1,x2,y2,d,*r;

  r=new double[mesh.nvtxs];

  for(n=0;n<mesh.nvtxs;n++) {
    translation[n]=-1;
    }

  for(n=0;n<mesh.nvtxs;n++) {
    x1=mesh.vertices[n].lon;
    y1=mesh.vertices[n].lat;
    r[n]=1.e+10;
    for(m=0;m<mesh.nvtxs;m++) {
      x2=x[m];
      y2=y[m];
      d=(x2-x1)*(x2-x1)+(y2-y1)*(y2-y1);
      if(d<r[n]) {
        r[n]=d;
        translation[n]=m;
        }
      if(r[n]<1.e-06) break;
      }
    }

  for(n=0;n<mesh.nvtxs;n++) {
    if(r[n]>1.e-03) goto error;
    }

  for(n=0;n<mesh.nvtxs;n++) {
    if(translation[n]==-1) goto error;
    }

  delete[] r;
  return (0);

 error:
  return (-1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int *surfref_reference(char *name, mesh_t mesh, int *ncolumns)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int   n,k,nitems,status;
  char  line[1000]="";
  string s;
  char c;
  FILE *in;
  double *x,*y;
  int *translation;

  in=fopen(name,"r");
  if(in == NULL) {
    printf("unable to open %s\n",name);
    return (NULL);
    }

  while((strcmp(line,"! -> \n")!=0) && !feof(in)) {
    fgets(line,1024,in);
    }

  nitems=0;
  c=fgetc(in);
  do {
    while (c == ' ') c=fgetc(in);
    if(c=='\n') break;
    nitems++;
    do { c=fgetc(in);   }  while ((c != '\n') && (c != ' ') && !feof(in));
    }  while ((c != '\n') && !feof(in));

//  *ncolumns=13;

/* *----------------------------------------------------------------------------
  do not count lon lat columns*/
  *ncolumns=nitems-2;

  rewind(in);
  strcpy(line,"");
  while((strcmp(line,"! -> \n")!=0) && !feof(in)) {
    fgets(line,1024,in);
    }

  x=new double[mesh.nvtxs];
  y=new double[mesh.nvtxs];
  for(n=0;n<mesh.nvtxs;n++) {
    nitems=fscanf(in,"%lf %lf",&x[n],&y[n]);
    if(nitems ==0) goto error;
    do { c=fgetc(in);   }  while ((c != '\n') && !feof(in));
    if(feof(in)) goto error;
    }

  fclose(in);

  translation=new int[mesh.nvtxs];
  status=cross_reference( mesh, x, y, translation);

  return (translation);

 error:
  fclose(in);
  return (NULL);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int surfref_loadr1(char *name, int nvalues, int ncolumns, int pos, int *translation, float *buffer)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int   n,k,nitems;
  char  line[1000]="";
  string s;
  char c;
  FILE *in;
  double x,y;
  float array[ncolumns];

  in=fopen(name,"r");
  if(in == NULL) {
    printf("unable to open %s\n",name);
    return (-1);
    }

  while((strcmp(line,"! -> \n")!=0) && !feof(in)) {
    fgets(line,1024,in);
    }

  for(n=0;n<nvalues;n++) {
    nitems=fscanf(in,"%lf %lf",&x,&y);
    for(k=0;k<ncolumns;k++) nitems=fscanf(in,"%f",&array[k]);
    if(nitems !=1) goto error;
    buffer[translation[n]]=array[pos];
    }

  fclose(in);

  return (0);

 error:
  fclose(in);
  return (-1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int surfref_reorder(char *name, mesh_t mesh, int ncolumns, int *translation, char *ouput)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int   n,k,nitems;
  char  line[1000]="";
  string s;
  char c;
  FILE *in,*out;
  double x,y;
  float **array;
  int nvalues=mesh.nvtxs;

  in=fopen(name,"r");
  if(in == NULL) {
    printf("unable to open %s\n",name);
    return (-1);
    }

  out=fopen(ouput,"w");
  if(out == NULL) {
    printf("unable to open %s\n",ouput);
    return (-1);
    }

  while((strcmp(line,"! -> \n")!=0) && !feof(in)) {
    fgets(line,1024,in);
    fprintf(out,"%s",line);
    }

  if(feof(in)) {
    printf("missing header (terminating with ! ->) in %s\n",name);
    }

  array=new float *[nvalues];
  for(k=0;k<nvalues;k++) array[k]=new float[ncolumns];

  for(n=0;n<nvalues;n++) {
    nitems=fscanf(in,"%lf %lf",&x,&y);
    for(k=0;k<ncolumns;k++) nitems=fscanf(in,"%f",&array[translation[n]][k]);
    if(nitems !=1) goto error;
    }

  fclose(in);

  for(n=0;n<nvalues;n++) {
    fprintf(out,"%lf %lf",mesh.vertices[n].lon,mesh.vertices[n].lat);
    for(k=0;k<ncolumns;k++) fprintf(out," %f",array[n][k]);
    fprintf(out,"\n");
//    if(nitems !=1) goto error;
    }

  fclose(out);

  for(k=0;k<nvalues;k++) delete[] array[k];
  delete[] array;

  return (0);

 error:
  fclose(out);
  return (-1);
}
