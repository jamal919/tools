
/*******************************************************************************

  T-UGO tools, 2006-2013

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

\author  Florent Lyard      LEGOS/CNRS, Toulouse, France. florent.lyard@legos.obs-mip.fr
\author  Laurent Roblou     LEGOS/CNRS, Toulouse, France
\author  Damien Allain      LEGOS/CNRS, Toulouse, France
\author  Yves Soufflet      LEGOS/CNRS, Toulouse, France
\author  Clément Mayet      LEGOS, Toulouse, France (PhD)
\author  Yoann Le Bars      LEGOS, Toulouse, France (PhD)
\author  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
\author  Frédéric Dupont    Université de Laval à Québec, Canada

\brief mesh geometry initialisation functions + fe_beta interpolation
*/
/*----------------------------------------------------------------------------*/

#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits>

#include "tools-structures.h"
#include "constants.h"

#include "functions.h"

#include "fe.def"
#include "fe.h"

#include "map.h"
#include "geo.h"
#include "rutin.h"
#include "geo.h"
#include "maths.h"


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

inline void push_corners_list(int (*tmp)[3],int *kk,int n1,int n2,int n3)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
//#pragma omp critical
  {
  tmp[*kk][0] = n1;
  tmp[*kk][1] = n2;
  tmp[*kk][2] = n3;
  (*kk)++;
  }
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int fe_spherical_template(mesh_t & mesh, T projection)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=0;
  double *x,*y;
  bool debug=true;
  frame_t frame;
  
  x=new double[mesh.nvtxs];
  y=new double[mesh.nvtxs];
  
  for (int n=0;n<mesh.nvtxs;n++) {
    x[n]=mesh.vertices[n].lon;
    y[n]=mesh.vertices[n].lat;
    }
    
  if(debug) {
    status=fe_minmax(mesh, frame);
    }
    
//   projection_to_geo(projection,x,y,x,y,mesh.nvtxs);
  projection_to_geo(projection,x,y,mesh.nvtxs,0);
  
  for (int n=0;n<mesh.nvtxs;n++) {
    mesh.vertices[n].lon=x[n];
    mesh.vertices[n].lat=y[n];
    }
    
  delete[] x;
  delete[] y;
  
  if(debug) {
    status=fe_minmax(mesh, frame);
    }
    
  mesh.type=SPHERICAL;
  
  return(status);
}


// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   int fe_spherical(mesh_t & mesh, projPJ projection)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   int status=0;
//   
//   status=fe_spherical_template(mesh, projection);
//   
//   return(status);
// }

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_spherical(mesh_t & mesh, const char *projection)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=0;
  
  status=fe_spherical_template(mesh, projection);
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int fe_cartesian_template(mesh_t & mesh, T projection)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=0;
  double *x,*y;
  
  x=new double[mesh.nvtxs];
  y=new double[mesh.nvtxs];
  
  for (int n=0;n<mesh.nvtxs;n++) {
    x[n]=mesh.vertices[n].lon;
    y[n]=mesh.vertices[n].lat;
    }
    
  geo_to_projection(projection,x,y,x,y,mesh.nvtxs);
  
  for (int n=0;n<mesh.nvtxs;n++) {
    mesh.vertices[n].lon=x[n];
    mesh.vertices[n].lat=y[n];
    }
    
  delete[] x;
  delete[] y;
  
  mesh.type=SPHERICAL;
  
  return(status);
}


// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   int fe_cartesian(mesh_t & mesh, projPJ projection)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   int status=0;
//   
//   status=fe_spherical_template(mesh, projection);
//   
//   return(status);
// }

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_cartesian(mesh_t & mesh, const char *projection)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=0;
  
  status=fe_spherical_template(mesh, projection);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_list(mesh_t *mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  Make an element list from nodes' neighbours description.
  If consecutive neighbours of a node are connected, they
  form an element with the node.
    
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int i,j,k,kk,n1,n2,n3,rotation=1;
  int max_nb,nndes;
  vertex_t *vertices=NULL;
  int (*tmp)[3]=NULL;
  int status;
  double angle;

  if(mesh->type==-1) mesh->type=0;

  deletep(&mesh->triangles);
  
  nndes=mesh->nvtxs;
  vertices=mesh->vertices;
  max_nb=mesh->nnghm;
  
#pragma omp parallel for
  for(i=0; i<nndes; i++) {
    status=fe_order(*mesh,i,-1., 0.,rotation);
    }

  for(i=0; i<nndes; i++) {
    if(vertices[i].nngh < 2 ) {
      printf("orphan node %d (less than 2 neighbours, actually %d)\n",i,vertices[i].nngh);
      }
    }

  exitIfNull(
    tmp=new int[2*nndes][3]
    );
  
//   struct timeval b4;
//   gettimeofday(&b4);

  kk = 0;
  for(n1=0;n1<nndes;n1++) {
    vertex_t *vertex=&vertices[n1];
    
    if(vertex->nngh == 2) {
      n2=vertex->ngh[0];
      if(n2 <= n1) continue;
      n3=vertex->ngh[1];
      if(n3 <= n1) continue;
      angle=fe_angle(*mesh, n1,n2,n3);
      if(angle  > 0.0) {
        status=fe_anb(*mesh, n2, n3);
        if(status==-1) {
          printf("not a triangular mesh, or mesh unsafe; abort\n");
          return(-1);
          }
        push_corners_list(tmp,&kk,n1,n2,n3);
        }
      else {
        push_corners_list(tmp,&kk,n1,n2,n3);
        }
      continue;
      }
    
    if(vertex->nngh<=1) continue;
    
    for(k=0; k<vertex->nngh-1; k++) {
      n2=vertex->ngh[k];
      if(n2 <= n1) continue;
      n3=vertex->ngh[k+1];
      if(n3 == -1) break;
      if(n3 <= n1) continue;
      j=fe_anb(*mesh,n3,n2);
      if(j!=-1) {
        push_corners_list(tmp,&kk,n1,n2,n3);
        }
      }
    
    n2=vertex->ngh[0];
    if(n2 <= n1) continue;
    n3=vertex->ngh[k];
    if(n3 <= n1) continue;
    j=fe_anb(*mesh,n3,n2);
    if(j!=-1) {
      angle=fe_angle_cartesian(*mesh, n1,n3,n2);
      if(angle  < 0.0) {
        printf("nodes %d(%g;%g) %d %d do not form a valid triangle (i.e. clock-wise triangle); skip\n",n1,vertex->lon,vertex->lat,n3,n2);
        continue;
        }
      push_corners_list(tmp,&kk,n1,n2,n3);
      }
    
    }
  
//   STDERR_BASE_LINE_FUNC("%gs\n",difftime(b4));
  
  mesh->ntriangles = kk;
  mesh->triangles=new triangle_t[mesh->ntriangles];

  for(i=0; i<mesh->ntriangles; i++) {
    for(j=0; j<3; j++) {
      mesh->triangles[i].vertex[j] = tmp[i][j];
      }
    }

  status=fe_geometry(mesh);
  if(status!=0) return(-1);

  if(mesh->vertices[0].nelmts!=-1) {
    status=fe_vertex_element_tables(mesh);
    }

  delete[] tmp;

  return status;
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_list_new(mesh_t *mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  Make an element list from nodes' neighbours description.
  If consecutive neighbours of a node are connected, they
  form an element with the node.
    
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int i,j,k,kk,n1,n2,n3,rotation=1;
  int max_nb,nndes;
  vertex_t *vertices=NULL;
  int (*tmp)[3]=NULL;
  int status;
  double angle;

  if(mesh->type==-1) mesh->type=0;

  deletep(&mesh->triangles);
  
  nndes=mesh->nvtxs;
  vertices=mesh->vertices;
  max_nb=mesh->nnghm;
  
#pragma omp parallel for
  for(i=0; i<nndes; i++) {
    status=fe_order(*mesh,i,-1., 0.,rotation);
    }

  for(i=0; i<nndes; i++) {
    if(vertices[i].nngh < 2 ) {
      printf("orphan node %d (less than 2 neighbours, actually %d)\n",i,vertices[i].nngh);
      }
    }

//   struct timeval b4;
//   gettimeofday(&b4);

  kk = 0;
  for(n1=0;n1<nndes;n1++) {
    vertex_t *vertex=&vertices[n1];
    
/*------------------------------------------------------------------------------
    special case: 1 neighbour or less */    
    if(vertex->nngh<=1) continue;
    
/*------------------------------------------------------------------------------
    special case: 2 neighbours only */    
    if(vertex->nngh == 2) {
      n2=vertex->ngh[0];
      if(n2 <= n1) continue;
      n3=vertex->ngh[1];
      if(n3 <= n1) continue;
      angle=fe_angle_cartesian(*mesh, n1,n3,n2);
      if(angle  > 0.0) {
        status=fe_anb(*mesh, n2, n3);
        if(status==-1) {
          printf("not a triangular mesh, or mesh unsafe; abort\n");
          return(-1);
          }
        push_corners_list(tmp,&kk,n1,n2,n3);
        }
      else {
        push_corners_list(tmp,&kk,n1,n2,n3);
        }
      continue;
      }
    
/*------------------------------------------------------------------------------
    standard case */    
    for(k=0; k<vertex->nngh-1; k++) {
      n2=vertex->ngh[k];
      if(n2 <= n1) continue;
      n3=vertex->ngh[k+1];
      if(n3 == -1) break;
      if(n3 <= n1) continue;
      j=fe_anb(*mesh,n3,n2);
      if(j!=-1) {
        push_corners_list(tmp,&kk,n1,n2,n3);
        }
      }
    
    n2=vertex->ngh[0];
    if(n2 <= n1) continue;
    n3=vertex->ngh[k];
    if(n3 <= n1) continue;
    j=fe_anb(*mesh,n3,n2);
    if(j!=-1) {
      angle=fe_angle_cartesian(*mesh, n1,n3,n2);
      if(angle  < 0.0) {
        printf("nodes %d(%g;%g) %d %d do not form a valid triangle (i.e. clock-wise triangle); skip\n",n1,vertex->lon,vertex->lat,n3,n2);
        continue;
        }
      push_corners_list(tmp,&kk,n1,n2,n3);
      }
    
    }
  
//   STDERR_BASE_LINE_FUNC("%gs\n",difftime(b4));
  
  mesh->ntriangles = kk;
  mesh->triangles=new triangle_t[mesh->ntriangles];

  for(i=0; i<mesh->ntriangles; i++) {
    for(j=0; j<3; j++) {
      mesh->triangles[i].vertex[j] = tmp[i][j];
      }
    }

  status=fe_geometry(mesh);
  if(status!=0) return(-1);

  if(mesh->vertices[0].nelmts!=-1) {
    status=fe_vertex_element_tables(mesh);
    }

  delete[] tmp;

  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_list(mesh_t *mesh, int ElementType)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  extern int fe_initaffine_spherical(const mesh_t & mesh, quadrangle_t *q, int m);
  int status;
  
  switch(ElementType) {
    case FE_TRIANGLE:
      status=fe_list(mesh);
      break;
    case FE_QUADRANGLE:
      status=quadrangle_list (*mesh, 4, 4);
      for(size_t m=0;m<mesh->nquadrangles;m++) {
        status=fe_initaffine_spherical(*mesh, &(mesh->quadrangles[m]), m);
        }
      break;
    default:
      TRAP_ERR_EXIT(-1, "element type not reckognized\n");
      break;
    }
  return (status);
    
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_geometry(mesh_t *mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------
  Initialise element geometry quantities
----------------------------------------------------------------------*/
{
  int i,m,n;
  int weird=0,severe=0;

  if(mesh->type==-1) mesh->type=0;

  int nprocs __attribute__((unused)) =initialize_OPENMP(-1,0);
#pragma omp parallel for
  for(n=0; n<mesh->ntriangles; n++) {
    int status __attribute__((unused)) =fe_initaffine(mesh,n);
    }

  for(n=0; n<mesh->ntriangles; n++) {
    if(mesh->triangles[n].Area <0.0) {
      printf("element %d has negative area (cw order)\n",n);
      weird++;
      }
    else if(mesh->triangles[n].Area ==0.0) {
      printf("element %d has null area (flat element), vertices :",n);
      for(int j=0; j<3; j++) {
        printf(" %d", mesh->triangles[n].vertex[j]);
        }
      printf("\n");
      severe++;
      }
    }

  if(severe!=0) return(-1);

  for(n = 0; n < mesh->nvtxs; n++) {
    mesh->vertices[n].mw = 0.0;
    }
  for(m = 0; m < mesh->ntriangles; m++) {
    for(i=0;i<3;i++) {
      n = mesh->triangles[m].vertex[i];
      mesh->vertices[n].mw += mesh->triangles[m].Area / 3.;
      }
    }

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_reducebw(mesh_t mesh, int maxpercent)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,m,n,nn,hbw,hbwmin,start,status;
  int nmax;
  mesh_t work;
  int *idx=NULL,*tmp=NULL,*old=NULL;
  int connex=1;

  status=fe_bandwidth(&mesh);
  hbwmin=mesh.hbw;
  
  if(maxpercent==100) {
    nmax=mesh.nvtxs;
    }
  else {
    nmax=mesh.nvtxs*(maxpercent/100.);
    }
  
  idx      =new int[mesh.nvtxs];
  tmp      =new int[mesh.nvtxs];
  old      =new int[mesh.nvtxs];
/*------------------------------------------------------------------------------
  */
  start=-1;

  while(start<nmax) {
next:
    start++;
    if (start==nmax) break;
    hbw=0;
    for(n = 0; n < mesh.nvtxs; n++) tmp[n]=-1;
    work.nvtxs=0;
    tmp[start]=work.nvtxs;
    old[work.nvtxs]=start;
    work.nvtxs++;
//     low=0;
//     first=low;
//     low=work.nvtxs;
    for(nn = 0; nn < work.nvtxs; nn++) {
      n=old[nn];
      for(k=0;k<mesh.vertices[n].nngh;k++) {
        m=mesh.vertices[n].ngh[k];
        if(tmp[m]==-1) {
          int mm=work.nvtxs;
          tmp[m]=mm;
          old[mm]=m;
          work.nvtxs++;
          updatemax(&hbw,mm-nn);
          if(hbw>=hbwmin) goto next;
//           updatemin(&low,nn);
          }
        }
      }
    printf("nvtex %d, start %6d, previous=%6d hbw=%6d\n",mesh.nvtxs,start,hbwmin,hbw);
/* *----------------------------------------------------------------------------
    warning: possible infinite loop if part of mesh not connected with the rest of it*/
//    if(work.nvtxs!=mesh.nvtxs) goto iterate;
    if(work.nvtxs!=mesh.nvtxs) {
      printf("mesh is not connex, exit...\n");
      check_error(-3,"mesh is not connex", __LINE__, __FILE__, 0);
      connex=0;
      }
    hbwmin=hbw;
    for(n = 0; n < mesh.nvtxs; n++) idx[n]=tmp[n];
    if(connex==0) break;
    }

  if(connex==1) work.nvtxs=mesh.nvtxs;
  work.nnghm=mesh.nnghm;
  work.vertices=new vertex_t[work.nvtxs];
  for (n=0; n<work.nvtxs; n++) work.vertices[n].null_value();

  for(n=0;n<work.nvtxs;n++) {
    m=idx[n];
    work.vertices[m].lon =mesh.vertices[n].lon;
    work.vertices[m].lat =mesh.vertices[n].lat;
    work.vertices[m].h   =mesh.vertices[n].h;
    work.vertices[m].code=0;
    work.vertices[m].nngh=mesh.vertices[n].nngh;
    work.vertices[m].ngh=new int[work.vertices[m].nngh];
    for(k=0;k<work.vertices[m].nngh;k++) {
      work.vertices[m].ngh[k]=idx[mesh.vertices[n].ngh[k]];
      }
    }

  if(connex==0) {
    fe_cleanvertices(&work, false);
    }
  status=fe_list(&work);
  if(connex==0) {
    status=fe_e2n (&work);
    }
  if(status!=0) {
    printf("unable to build the element list\n");
    goto error;
    }
  status=fe_bandwidth(&work);
  status=fe_edgetable(&work,0,0);
  status=fe_codetable2(&work,0,1,0);
  status=fe_savemesh("reduced.nei",(int) MESH_FILE_FORMAT_TRIGRID, work);

  work.hbw=hbwmin;

  delete[] idx;
  delete[] tmp;
  delete[] old;

  if(connex==1) return(0);
  else return(-2);

error:
  return(-1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_reducebw_try(mesh_t & mesh, int *ngh, int *nngh, int maxnb, int start, int hbwmin, int & hbw, int *idx)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,m,n,nn;
  int nvtxs;
  int *tmp,*old;
  int next, ptr;
  
  tmp=new int[mesh.nvtxs];
  old=new int[mesh.nvtxs];

/*------------------------------------------------------------------------------
  */
  hbw=0;
  for(n = 0; n < mesh.nvtxs; n++) tmp[n]=-1;
  
  nvtxs=0;
  tmp[start]=nvtxs;
  old[nvtxs]=start;
  nvtxs++;
  
  for(nn = 0; nn < nvtxs; nn++) {
    n=old[nn];
    ptr=maxnb*n;
    for(k=0;k<nngh[n];k++) {
      m=ngh[ptr];
      if(tmp[m]==-1) {
        next=nvtxs;
        if(next-nn>=hbwmin) {
          hbw=next-nn;
          goto abort;
//           ptr++;
//           goto increment;
          }
        updatemax(&hbw,next-nn);
        tmp[m]=next;
        old[next]=m;
        nvtxs++;
        }
      ptr++;
      }
    }
  
  printf("nvtex %d, start %d, previous=%d hbw=%d\n",mesh.nvtxs,start,hbwmin,hbw);
  
/* *----------------------------------------------------------------------------
  warning: possible infinite loop if part of mesh not connected with the rest of it*/
  if(nvtxs!=mesh.nvtxs) {
    printf("mesh is not connex:\n");
    printf("#vertices : %d over %d\n",nvtxs,mesh.nvtxs);
    for(n=0;n<mesh.nvtxs;n++) {
      if(tmp[n]==-1) {
        printf("orphan vertex %d : %lf %lf\n",n,mesh.vertices[n].lon,mesh.vertices[n].lat);
        }
      }
    TRAP_ERR_EXIT(-3,"mesh is not connex\n");
    }
  
  for(n = 0; n < mesh.nvtxs; n++) idx[n]=tmp[n];

abort:
  delete[] tmp;
  delete[] old;

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_reducebw(mesh_t & mesh, int incr, int max_percent, const char *output)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,m,n,hbwmin,start,status;
  mesh_t work;
  int *idx,*ngh,*nngh;
  int maxnb;
  int percent;
  
  status=fe_bandwidth(&mesh);
  hbwmin=mesh.hbw;
  printf("original hbw=%d\n",hbwmin);
    
  maxnb=-1;
  for(n = 0; n < mesh.nvtxs; n++) maxnb=max(mesh.vertices[n].nngh,maxnb);

  idx      =new int[mesh.nvtxs];
  ngh      =new int[mesh.nvtxs*maxnb];
  nngh     =new int[mesh.nvtxs];
  
  for(n = 0; n < mesh.nvtxs; n++) {
    idx[n]=n;
    nngh[n]=mesh.vertices[n].nngh;
    for(k=0;k<mesh.vertices[n].nngh;k++) {
      ngh[maxnb*n+k]=mesh.vertices[n].ngh[k];
      }
    }
/*------------------------------------------------------------------------------
  */
  start=-1;
  while(start<mesh.nvtxs) {
    int hbw;
    start+=incr;
    if (start>=mesh.nvtxs) break;
    percent=(int) floor((100.*start)/mesh.nvtxs+0.5);
    if (percent> max_percent) break;
    status=fe_reducebw_try(mesh, ngh, nngh, maxnb, start, hbwmin, hbw, idx);
    if(hbw<hbwmin) {
      printf("nvtex %d, start %d, previous=%d hbw=%d\n",mesh.nvtxs,start,hbwmin,hbw);
      hbwmin=hbw;
      }
    }

  work.nvtxs=mesh.nvtxs;
  work.nnghm=mesh.nnghm;
  work.vertices=new vertex_t[work.nvtxs];
  for (n=0; n<work.nvtxs; n++) work.vertices[n].null_value();

  for(n=0;n<work.nvtxs;n++) {
    m=idx[n];
    work.vertices[m].lon =mesh.vertices[n].lon;
    work.vertices[m].lat =mesh.vertices[n].lat;
    work.vertices[m].code=0;
    work.vertices[m].nngh=mesh.vertices[n].nngh;
    work.vertices[m].ngh=new int[work.vertices[m].nngh];
    for(k=0;k<work.vertices[m].nngh;k++) {
      work.vertices[m].ngh[k]=idx[mesh.vertices[n].ngh[k]];
      }
    }

  status=fe_list(&work);
  if(status!=0) {
    printf("unable to build the element list\n");
    goto error;
    }
  status=fe_bandwidth(&work);
  status=fe_edgetable(&work,0,0);
  status=fe_codetable2(&work,0,1,0);
  status=fe_savemesh(output,(int) MESH_FILE_FORMAT_TRIGRID, work);

  work.hbw=hbwmin;

  delete[] idx;
  delete[] ngh;
  delete[] nngh;
  
  mesh.destroy();
  mesh=work;

  return(0);

error:
  return(-1);
}



/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_bandwidth(mesh_t *mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,m,n,hbw=0;

/*------------------------------------------------------------------------------
  check for consistency*/
  for(n = 0; n < mesh->nvtxs; n++) {
    for(k=0;k<mesh->vertices[n].nngh;k++) {
      m=mesh->vertices[n].ngh[k];
/*
      if(abs(m-n)>10000) {
        printf("n=%d\n",n);
        }
*/
      hbw=max(hbw,abs(m-n));
      }
    }
  mesh->hbw=hbw;
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_element_bw(mesh_t *mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int  i,ii,n,j,jj,k,elt_nnghmax;
  int  **elist=NULL,*nlist=NULL,nmax;
  int hbw;
  int  ne,nn,maxnb;

  ne=mesh->ntriangles;
  nn=mesh->nvtxs;
  maxnb=mesh->nnghm;

  exitIfNull(
    nlist=(int *) malloc(nn*sizeof(int))
    );

/*------------------------------------------------------------------------------
  Establish element bandwidth
  developped to optimize point detection in external mesh
----------------------------------------------------------------------*/

  nmax=0;
  for(j=0; j<nn; j++) nlist[j] = 0;
  for(j=0; j<ne; j++) {
    for(i=0; i<3; i++) {
      n = mesh->triangles[j].vertex[i];
      nlist[n]++;
      updatemax(&nmax, nlist[n]);
      }
    }
  printf("maximum elements/node = %d \n",nmax);

  elist = (int **) malloc(nn*sizeof(int *));
  for(j=0; j<nn; j++) {
    exitIfNull(
      elist[j]=(int *) malloc(nmax*sizeof(int))
      );
  }
  for(j=0; j<nn; j++) nlist[j] = 0;

  for(j=0; j<ne; j++) {
    for(i=0; i<3; i++) {
      n = mesh->triangles[j].vertex[i];
      elist[n][nlist[n]]=j;
      nlist[n]++;
      }
    }

  hbw=0;

  for(j=0; j<ne; j++) {
    for(i=0; i<3; i++) {
      n = mesh->triangles[j].vertex[i];
      for(ii=0; ii<nlist[n]; ii++) {
        hbw=max(abs(j-elist[n][ii]),hbw);
        }
      }
    }

  printf("element half-bandwidth = %d \n",hbw);
  mesh->hbw=hbw;

  elt_nnghmax=3*maxnb;

  exitIfNull(
    mesh->elt_nngh=(int *) malloc(ne*sizeof(int))
    );
  for(j=0; j<ne; j++) mesh->elt_nngh[j] = 0;
  exitIfNull(
    mesh->elt_nghl=(int **) malloc(ne*sizeof(int *))
    );
  for(j=0; j<ne; j++) {
    exitIfNull(
      mesh->elt_nghl[j]=(int *) malloc(elt_nnghmax*sizeof(int))
      );
  }

  for(j=0; j<ne; j++) {
    for(i=0; i<3; i++) {
      n = mesh->triangles[j].vertex[i];
      for(ii=0; ii<nlist[n]; ii++) {
        jj=elist[n][ii];
        for(k=0;k<mesh->elt_nngh[j];k++)
          if(mesh->elt_nghl[j][k]==jj) goto skip;
        mesh->elt_nngh[j]++;
        if(mesh->elt_nngh[j] > elt_nnghmax) goto error;
        mesh->elt_nghl[j][mesh->elt_nngh[j]-1]=jj;
      skip:
        k=0;
        }
      }
    }

  elt_nnghmax=0;
  for(j=0; j<ne; j++) {
    elt_nnghmax=max(mesh->elt_nngh[j],elt_nnghmax);
    }

  printf("element max neighbours = %d (presumedly less than %d) \n",elt_nnghmax,3*maxnb);

  free(nlist);
  for(j=0; j<nn; j++) free(elist[j]);

  return(0);

 error:
  return(-1);

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_initaffine_spherical(mesh_t *mesh,const int m, bool ignore_CW)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,n1,n2,n3,status,dum;
  double t1,t2,t3,p1,p2,p3;
  double dX, dY, ds;
  double alpha;
  double dlon,dlat;
/*------------------------------------------------------------------------------
  product of local earth radius in meters and of degree to radian conversion */
  const double Rd2r=6.3675e+06*d2r;

/*------------------------------------------------------------------------------
  local*/
  double ef_tref,ef_pref;
  double ef_ctx,ef_cty,ef_cpx,ef_cpy,ef_jacobien;
/*------------------------------------------------------------------------------
  passed through triangle structure*/
  double ef_dxdt,ef_dydt,ef_dxdp,ef_dydp;
  
  triangle_t *triangle=&mesh->triangles[m];// pointer to current triangle

 redo:

  n1=triangle->vertex[0];
  n2=triangle->vertex[1];
  n3=triangle->vertex[2];

  t1=mesh->vertices[n1].lon;
  t2=mesh->vertices[n2].lon;
  t3=mesh->vertices[n3].lon;

  if(mesh->type==0) t2=degree_recale(t2,t1);
  if(mesh->type==0) t3=degree_recale(t3,t1);

  p1=mesh->vertices[n1].lat;
  p2=mesh->vertices[n2].lat;
  p3=mesh->vertices[n3].lat;

  ef_tref=  t1;
  ef_pref=  p1;

  ef_ctx=  t2-ef_tref;
  ef_cty=  t3-ef_tref;
  ef_cpx=  p2-ef_pref;
  ef_cpy=  p3-ef_pref;

  ef_jacobien=ef_ctx*ef_cpy-ef_cty*ef_cpx;

  if (ef_jacobien > 0.) {
    const double r2d_ef_jacobien=r2d/ef_jacobien;
    ef_dxdt= ef_cpy*r2d_ef_jacobien;
    ef_dydt=-ef_cpx*r2d_ef_jacobien;
    ef_dxdp=-ef_cty*r2d_ef_jacobien;
    ef_dydp= ef_ctx*r2d_ef_jacobien;
    status=MESH_STATUS_OK;
    }
  else if(ef_jacobien == 0.) {
    ef_dxdt= 0;
    ef_dydt= 0;
    ef_dxdp= 0;
    ef_dydp= 0;
    status=MESH_STATUS_FLAT_ELEMENT;
    }
  else {
    status=MESH_STATUS_CW_ELEMENT;
    if(ignore_CW) {
      dum=triangle->vertex[2];
      triangle->vertex[2]=triangle->vertex[1];
      triangle->vertex[1]=dum;
      goto redo;
      }
    else {
      printf("clockwise element %d : nodes %d %d %d\n",m,n1,n2,n3);
      }
    }

  triangle->dxdt  =  ef_dxdt;
  triangle->dydt  =  ef_dydt;
  triangle->dxdp  =  ef_dxdp;
  triangle->dydp  =  ef_dydp;

  triangle->dtdx  =  ef_ctx;
  triangle->dpdx  =  ef_cpx;
  triangle->dtdy  =  ef_cty;
  triangle->dpdy  =  ef_cpy;

  triangle->jacobien  =  ef_jacobien;

  triangle->t_base  =  ef_tref;
  triangle->p_base  =  ef_pref;

  triangle->DP[0] =  Rd2r * (t2-t3);
  triangle->DP[1] =  Rd2r * (t3-t1);
  triangle->DP[2] =  Rd2r * (t1-t2);
  triangle->DQ[0] =  Rd2r * (p2-p3);
  triangle->DQ[1] =  Rd2r * (p3-p1);
  triangle->DQ[2] =  Rd2r * (p1-p2);
/*------------------------------------------------------------------------------
  WARNING, PSEUDO area (actual area = pseudo*cos(latitude))
  note: PSEUDO Area=0.5*ef_jacobien*R²*dtr² * /
  triangle->Area = 0.5*R* (t1*triangle->DQ[0]
                                 + t2*triangle->DQ[1]
                                 + t3*triangle->DQ[2])*d2r;
/* *@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
Development notes on pseudo area A :
A= (x0-x1)(y0-y2) - (y0-y1)(x0-x2)
 = (y0-y1)(x2-x0) - (x0-x1)(y2-y0)
 =  DQ[2] *DP[1]  -  DP[2] *DQ[1]
**/
  triangle->Area = 0.5*(triangle->DQ[2]*triangle->DP[1] - triangle->DP[2]*triangle->DQ[1]);

  if(isnan(triangle->Area)) {
    triangle->Area=0.;
    }
  triangle->TrueArea =triangle->Area*cos((p1+p2+p3)*d2r/3.);
  triangle->cosinus  = cos((p1+p2+p3)*M_PI/540.); /* /540 = /180 /3 */

  for(i=0;i<3;i++) {
    int target;
    n1=triangle->vertex[i];
    n2=triangle->vertex[(i+1)%3];
    t1=mesh->vertices[n1].lon;
    t2=mesh->vertices[n2].lon;

    if(mesh->type==0) t2=degree_recale(t2,t1);

    p1=mesh->vertices[n1].lat;
    p2=mesh->vertices[n2].lat;
    alpha = 0.5 * (p1+p2)*d2r;
    dlon=t2-t1;
    dlat=p2-p1;
    dX = Rd2r * cos(alpha) * dlon;
    dY = Rd2r * dlat;
    ds = hypot(dX,dY);
/*------------------------------------------------------------------------------
    outward normal is clockwise rotated */
    target=(i+2)%3;
    triangle->nx[target] = +dY / ds;
    triangle->ny[target] = -dX / ds;
    triangle->l [target] = ds;
    
    ds=geo_haversin(t1,p1,t2,p2);
    triangle->l [target] = ds;
    
    }

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_initaffine_cartesian(mesh_t *mesh,int m, bool ignore_CW)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,n1,n2,n3,status,dum;
  double t1,t2,t3,p1,p2,p3;
  double dX, dY, ds;
  double alpha;
  double dlon,dlat;
/*------------------------------------------------------------------------------
  local earth radius in meters*/
  double R=6.3675e+06;

/*------------------------------------------------------------------------------
  local*/
  double ef_tref,ef_pref;
  double ef_ctx,ef_cty,ef_cpx,ef_cpy,ef_jacobien;
/*------------------------------------------------------------------------------
  passed through triangle structure*/
  double ef_dxdt,ef_dydt,ef_dxdp,ef_dydp;

 redo:

  n1=mesh->triangles[m].vertex[0];
  n2=mesh->triangles[m].vertex[1];
  n3=mesh->triangles[m].vertex[2];

  t1=mesh->vertices[n1].lon;
  t2=mesh->vertices[n2].lon;
  t3=mesh->vertices[n3].lon;

  p1=mesh->vertices[n1].lat;
  p2=mesh->vertices[n2].lat;
  p3=mesh->vertices[n3].lat;

  ef_tref=  t1;
  ef_pref=  p1;

  ef_ctx=  t2-ef_tref;
  ef_cty=  t3-ef_tref;
  ef_cpx=  p2-ef_pref;
  ef_cpy=  p3-ef_pref;

  ef_jacobien=ef_ctx*ef_cpy-ef_cty*ef_cpx;

  if (ef_jacobien > 0.0) {
    ef_dxdt= ef_cpy/ef_jacobien;
    ef_dydt=-ef_cpx/ef_jacobien;
    ef_dxdp=-ef_cty/ef_jacobien;
    ef_dydp= ef_ctx/ef_jacobien;
    status=MESH_STATUS_OK;
    }
  else if(ef_jacobien == 0.0) {
    ef_dxdt= 0;
    ef_dydt= 0;
    ef_dxdp= 0;
    ef_dydp= 0;
    status=MESH_STATUS_FLAT_ELEMENT;
    }
  else {
//     dum=mesh->triangles[m].vertex[2];
//     mesh->triangles[m].vertex[2]=mesh->triangles[m].vertex[1];
//     mesh->triangles[m].vertex[1]=dum;
//     status=MESH_STATUS_CW_ELEMENT;
//     goto redo;
    status=MESH_STATUS_CW_ELEMENT;
    if(ignore_CW) {
      dum=mesh->triangles[m].vertex[2];
      mesh->triangles[m].vertex[2]=mesh->triangles[m].vertex[1];
      mesh->triangles[m].vertex[1]=dum;
      goto redo;
      }
    else {
      printf("clockwise element: %d %d %d\n",n1,n2,n3);
      printf("anti-clockwise element: %d %d %d\n",n1,n2,n3);
      }
    }

  mesh->triangles[m].dxdt  =  ef_dxdt;
  mesh->triangles[m].dydt  =  ef_dydt;
  mesh->triangles[m].dxdp  =  ef_dxdp;
  mesh->triangles[m].dydp  =  ef_dydp;

  mesh->triangles[m].dtdx  =  ef_ctx;
  mesh->triangles[m].dpdx  =  ef_cpx;
  mesh->triangles[m].dtdy  =  ef_cty;
  mesh->triangles[m].dpdy  =  ef_cpy;

  mesh->triangles[m].jacobien  =  ef_jacobien;

  mesh->triangles[m].t_base  =  ef_tref;
  mesh->triangles[m].p_base  =  ef_pref;

  mesh->triangles[m].DP[0] =  (double) R *  (t2-t3);
  mesh->triangles[m].DP[1] =  (double) R *  (t3-t1);
  mesh->triangles[m].DP[2] =  (double) R *  (t1-t2);
  mesh->triangles[m].DQ[0] =  (double) R *  (p2-p3);
  mesh->triangles[m].DQ[1] =  (double) R *  (p3-p1);
  mesh->triangles[m].DQ[2] =  (double) R *  (p1-p2);
/*------------------------------------------------------------------------------
  WARNING, PSEUDO area (actual area = pseudo*cos(latitude))
  note: PSEUDO Area=0.5*ef_jacobien*R²*dtr² */
  mesh->triangles[m].Area = 0.5* (t1*mesh->triangles[m].DQ[0]
                                 + t2*mesh->triangles[m].DQ[1]
                                 + t3*mesh->triangles[m].DQ[2]);

  if(isnan(mesh->triangles[m].Area)==1) {
    mesh->triangles[m].Area=0.;
    }
  mesh->triangles[m].TrueArea=mesh->triangles[m].Area;

  for(i=0;i<3;i++) {
    n1=mesh->triangles[m].vertex[i];
    n2=mesh->triangles[m].vertex[(i+1)%3];
    t1=mesh->vertices[n1].lon;
    t2=mesh->vertices[n2].lon;
    p1=mesh->vertices[n1].lat;
    p2=mesh->vertices[n2].lat;
    alpha = 0.5 * (p1+p2)*d2r;
    dlon=t2-t1;
    dlat=p2-p1;
    dX =  dlon;
    dY =  dlat;
    ds = sqrt(dX * dX + dY * dY);
/*------------------------------------------------------------------------------
   outward normal is clockwise rotated */
    mesh->triangles[m].nx[i] = +dY / ds;
    mesh->triangles[m].ny[i] = -dX / ds;
    mesh->triangles[m].l[i] = ds;
    }
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_initaffine(mesh_t *mesh,int m, bool ignore_CW)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;

  switch(mesh->type) {
    case  0:
      status=fe_initaffine_spherical(mesh, m, ignore_CW);
      break;
    case  1:
      status=fe_initaffine_cartesian(mesh, m, ignore_CW);
      break;
    default:
      STDOUT_BASE_LINE("fe_initaffine: mesh type undefined, set as spherical by default\n");
      mesh->type=0;
      status=fe_initaffine_spherical(mesh, m, ignore_CW);
      break;
    }

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_initaffine(mesh_t *mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int m,status;

  for(m=0;m<mesh->ntriangles; m++) {
    status=fe_initaffine(mesh,m);
    if(status!=0) return(-1);
    }
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int fe_affine_directe(const triangle_t & e,double t,double p,double *ksi,double *eta,int mode)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  calcule les coordonnees ksi,eta dans le triangle droit de reference

     longitude=t(1)+x*[t(2)-t(1)]+y*[t(3)-t(1)]
     latitude =p(1)+x*[p(2)-p(1)]+y*[p(3)-p(1)]

     determinant = [t(2)-t(1)]*[p(3)-p(1)]-[t(3)-t(1)]*[p(2)-p(1)]

     x = ([t-t(1)]*[p(3)-p(1)]-[t(3)-t(1)]*[p-p(1)])/determinant
       = ([t-t(1)]*cpy-cty*[p-p(1)])/determinant
     y = (ctx*[p-p(1)]-[t-t(1)]*cpx)/determinant

     Units : degrees

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  double teta,phi,dt,dp,tref,x,y;

  tref=e.t_base;

  if(mode==0) teta=geo_recale(t,tref,180.0);
  else teta=t;
  
  phi = p;

  dt=(teta-e.t_base)*d2r;
  dp=(phi -e.p_base)*d2r;

  x=e.dxdt*dt+e.dxdp*dp;
  y=e.dydt*dt+e.dydp*dp;

  *ksi= x;
  *eta= y;
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_affine_directe3D(prism_t prism, double t, double p, double h,double *x,double *y, double *z,int mode)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,status;
  double beta[3],bottom,dH;
  triangle_t triangle;

  status=fe_affine_directe(*(prism.triangle), t,p,x,y,mode);
  fe_LGP1base(*x, *y,beta);
  
/*------------------------------------------------------------------------------
  0 : bottom triangle
  1 : top triangle
-----------------------------------------------------------------------*/

  for (k=0;k<3;k++) {
    dH+=beta[k]*(prism.h[1][k]-prism.h[0][k]);
    bottom=beta[k]*prism.h[0][k];
    }
  *z=(h-bottom)/dH;
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_affine_inverse(triangle_t e,double *t, double *p,double ksi,double eta)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
     
  calcule les coordonnees teta,phi dans le triangle reel

     longitude=t(1)+x*[t(2)-t(1)]+y*[t(3)-t(1)]
     latitude =p(1)+x*[p(2)-p(1)]+y*[p(3)-p(1)]

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{

  *t=e.t_base+e.dtdx*ksi+e.dtdy*eta;
  *p=e.p_base+e.dpdx*ksi+e.dpdy*eta;

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_beta(mesh_t & mesh, double t, double p, int hint, int *elt,int *nodes, double *beta)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  find element containing (t,p) point, and return interpolation weights
  
  limited to LGP1 interpolation
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int    i,found,status;
  double tnew,pnew,x,y;
  bool debug=false;

  *elt=-1;
  
  found=fe_whichelement(mesh,t,p,hint);

  if (found ==-1) {
    found= fe_whichnearest(mesh,t,p,&tnew,&pnew);
    if(found<0) goto error;
    if(geo_km(t,p,tnew,pnew) > 20.0) goto error;
    if(debug)  printf ("fe_beta: %d %f %f %f %f\n",found,t,p,tnew,pnew);
    }
  else {
    tnew=t;
    pnew=p;
    }

  status=fe_initaffine(&mesh,found);
  if(status != 0) goto error;

  fe_affine_directe(mesh.triangles[found],tnew,pnew,&x,&y,0);

  if(debug)  printf ("fe_beta: %f %f %f %f %f %f\n",t,p,tnew,pnew,x,y);

  status=fe_LGPbase(mesh, LGP1,  x,  y, beta);
  if(status == 0) {
    for(i=0;i<3;i++) {
      nodes[i]=mesh.triangles[found].vertex[i];
      }
    }
  else {
    printf ("error in fe_beta: %f %f %f %f %f %f %f %f %f\n",t,p,tnew,pnew,x,y,beta[0],beta[1],beta[2]);
    }
  *elt=found;
  return(status);

 error:
  return(-1);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_beta(mesh_t & mesh, double t, double p, int *element, int *nodes, double *beta)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=fe_beta(mesh, t, p,  -1, element, nodes, beta);
  
  return(status);
}
  

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_beta(mesh_t & mesh, const discretisation_t & descriptor, double t, double p, int hint, int *element, int *nodes, double *beta)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  find element containing (t,p) point, and return element and interpolation 
  weights
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int    i,found,discretisation,status;
  double tnew,pnew,x,y;

  found=fe_whichelement(mesh,t,p);

  if (found ==-1) {
    found= fe_whichnearest(mesh,t,p,&tnew,&pnew);
    if(found<0) goto error;
    if(geo_km(t,p,tnew,pnew) > 20.0) goto error;
    }
  else {
    tnew=t;
    pnew=p;
    }

  status=fe_initaffine(&mesh,found);
  if(status != 0) goto error;

  fe_affine_directe(mesh.triangles[found],tnew,pnew,&x,&y,0);

  discretisation=descriptor.type;
  status=fe_LGPbase(mesh, discretisation, x, y, beta);
  if(status == 0) {
    for(i=0;i<descriptor.nnpe;i++) {
//      nodes[i]=mesh.triangles[found].vertex[i];
      nodes[i]=descriptor.NIbE[found][i];
      }
    }
  else {
    printf ("error in ef_beta: %f %f %f %f %f %f %f %f %f\n",t,p,tnew,pnew,x,y,beta[0],beta[1],beta[2]);
    }
  
//   printf ("fe_beta: %lf %lf, %lf %lf : %lf %lf %lf %lf %lf %lf %lf\n",t,p,x,y,beta[0],beta[1],beta[2],beta[3],beta[4],beta[5]);
  
  *element=found;
  return(status);

 error:
  *element=-1;
  return(-1);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_beta(mesh_t & mesh, const discretisation_t & descriptor, double t, double p, int *element, int *nodes, double *beta)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=fe_beta(mesh, descriptor, t, p,  -1, element, nodes, beta);
  
  return(status);
}
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <class T> int fe_beta_inlist_template(const mesh_t & mesh,const T & list,int nlisted,double t,double p,int *elt,int *node,double *beta)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,found,status;
  double tnew,pnew,x,y;

  found=fe_whichelement_inlist(mesh,list,nlisted,t,p);

  if (found ==-1) {
    return(1);
    }
  else {
    tnew=t;
    pnew=p;
    }

  status=fe_test_area(mesh,found);
  if(status)TRAP_ERR_RETURN(status,1,"%s:fe_test_area(,%d)=%d\n",__func__,found,status);

  fe_affine_directe(mesh.triangles[found],tnew,pnew,&x,&y,0);
/*  printf ("fe_beta: %f %f %f %f %f %f\n",t,p,tnew,pnew,x,y);*/
  status=fe_LGPbase(mesh, LGP1,  x,  y, beta);
  if(status == 0) {
    for(i=0;i<3;i++) {
      node[i]=mesh.triangles[found].vertex[i];
      }
    }
  else STDOUT_BASE_LINE("%s:fe_LGPbase() error : %f %f %f %f %f %f %f %f %f\n",__func__,t,p,tnew,pnew,x,y,beta[0],beta[1],beta[2]);
  *elt=found;
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int fe_beta_inlist(const mesh_t & mesh,const int *list,int nlisted,double t,double p,int *elt,int *node,double *beta)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=fe_beta_inlist_template(mesh,list,nlisted,t,p,elt,node,beta);
  
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int fe_beta_inlist(const mesh_t & mesh,const vector<int> & list,double t,double p,int *elt,int *node,double *beta)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=fe_beta_inlist_template(mesh,list,list.size(),t,p,elt,node,beta);
  
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 mesh_t::mesh_t(std::string filename, int format)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{

  vertices = 0;
  triangles = 0;
  quadrangles = 0;
  edges = 0;
  faces = 0;
  limits = 0;
  elt_nghl = 0;
  elt_nngh = 0;
  nvtxs = ntriangles = nquadrangles = nnghm = nedges = nlimits = nlayers = nlevels = 0;
  hbw = level = type = units = circular = -1;
  
  
  int status = 1;

  //read mesh info from file
  if( (status = fe_readmesh( (char*) filename.c_str(), format, this)) != 0){
    *this = null_mesh;
    }
  
  //build ...
  status = 1;
  status = fe_list(this);
  if(status!=0){
    *this = null_mesh;
    };
  
  //build ...
  status = 1;
  int verbose = 0;
  status = fe_edgetable(this, 0, verbose);
  if(status!=0){
    *this = null_mesh;
    };

  
  // build ...
  build_vindex(*this);
  
  //build...
  fe_edge_crosstables01(this);
  
  //build...
  status = 1;
  status = fe_vertex_crosstables01(this);
  if(status!=0){
    *this = null_mesh;
    };
  
  //build...
  fe_vertex_crosstables02(this);
      
  //build...
  fe_element_crosstables(*this);
  
  //build...
  for(size_t m = 0; m < this->ntriangles; m++) {
    status = 1;
    status = fe_initaffine(this, m);
    if(status!=0){
      *this = null_mesh;
      };
    }
  fe_integrale_init();
         
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  mesh_t& mesh_t::initialize(std::string filename, int format)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  // read mesh info from file
  int status = 1;
  if( (status = fe_readmesh( (char*) filename.c_str(), format, this)) != 0){
    return(*this);
    }
  
  // build ...
  status = 1;
  status = fe_list(this);
  if(status!=0){
    return(*this);
    }
  
  // build ...
  status = 1;
  int verbose = 0;
  status = fe_edgetable(this, 0, verbose);
  if(status!=0){
    return(*this);
    }
  
  // build...
  fe_edge_crosstables01(this); /// Use of uninitialised value of size 8 (fe08.cpp:837)
  
  // build...
  status = 1;
  status = fe_vertex_crosstables01(this);
  if(status!=0){
    return(*this);
    }
  
  // build the edges cross-tables for vertices
  fe_vertex_crosstables02(this);
      
  // build...
  fe_element_crosstables(*this);
  
  // build...
  for(size_t m = 0; m < this->ntriangles; m++) {
    status = 1;
    status = fe_initaffine(this, m);
    if(status!=0){
      return(*this);
      }
    }
  fe_integrale_init();
  
  return(*this);
  
}
  
