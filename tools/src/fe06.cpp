

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

#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "tools-structures.h"
#include "constants.h"

#include "fe.def"
#include "fe.h"

#include "map.h"
#include "geo.h"
#include "polygones.h"


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_triangulateio_init(triangulateio & out)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  out.pointlist = (REAL *) NULL;            /* Not needed if -N switch used. */
  out.pointattributelist = (REAL *) NULL;   /* Not needed if -N switch used or number of point attributes is zero: */
  out.pointmarkerlist = (int *) NULL;       /* Not needed if -N or -B switch used. */
  out.trianglelist = (int *) NULL;          /* Not needed if -E switch used. */
  out.triangleattributelist = (REAL *) NULL;/* Not needed if -E switch used or number of triangle attributes is zero: */
  out.neighborlist = (int *) NULL;          /* Needed only if -n switch used. */
  out.segmentlist = (int *) NULL;           /* Needed only if segments are output (-p or -c) and -P not used: */
  out.segmentmarkerlist = (int *) NULL;     /* Needed only if segments are output (-p or -c) and -P and -B not used: */
  out.edgelist = (int *) NULL;              /* Needed only if -e switch used. */
  out.edgemarkerlist = (int *) NULL;        /* Needed if -e used and -B not used. */
  out.holelist      =0;
  out.normlist = 0;
  out.regionlist = 0;

  out.numberofpoints = 0;
  out.numberofpointattributes = 0;
  out.numberofsegments = 0;

  out.numberofholes = 0;
  out.numberofregions  = 0;

  out.numberofedges = 0;
    
  out.numberoftriangleattributes = 0;
  
//   out.trianglelist     = 0;
  out.trianglearealist = 0;
//   out.neighborlist     = 0;
//   out.segmentlist      = 0;
//   out.edgelist         = 0;
//   out.edgemarkerlist   = 0;
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  mesh_t *fe_createnodes(plg_t *polygones, int npolygones, point_t *interiors, int ninteriors)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,j,k,l,m,n,status;
  int count,s, nlimits;
  mesh_t *mesh;
  
  mesh=new mesh_t;

  count=0;
  nlimits=0;
  for(s=0;s<npolygones;s++) {
    if(polygones[s].npt==0) continue;
    count+=polygones[s].npt-1;
    nlimits++;
    }

  count+=ninteriors;

  mesh->nvtxs=count;
  mesh->nlimits=nlimits;

  mesh->vertices=new vertex_t[mesh->nvtxs];
  mesh->limits=new limit_t[mesh->nlimits];

  n=0;
  nlimits=0;
  for(s=0;s<npolygones;s++) {
    if(polygones[s].npt==0) continue;
    mesh->limits[nlimits].nvertex=polygones[s].npt-1;
    mesh->limits[nlimits].vertex=new int[mesh->limits[nlimits].nvertex];
    for(k=0;k<polygones[s].npt-1;k++) {
      mesh->vertices[n].lon=polygones[s].x[k];
      mesh->vertices[n].lat=polygones[s].y[k];
      mesh->vertices[n].nngh=0;
      mesh->vertices[n].code=nlimits+1;
      mesh->limits[nlimits].vertex[k]=n;
      n++;
      }
    nlimits++;
    }

  for(k=0;k<ninteriors;k++) {
    mesh->vertices[n].lon=interiors[k].t;
    mesh->vertices[n].lat=interiors[k].p;
    mesh->vertices[n].code=0;
    mesh->vertices[n].nngh=0;
    n++;
    }

  return(mesh);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_createnodes(vector<plg_t> & polygons, double *x, double *y, double *z, int ninteriors, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,j,k,l,m,n,status;
  int count,s, nlimits;

  count=0;
  nlimits=0;
  for(s=0;s<polygons.size();s++) {
    if(polygons[s].npt==0) continue;
    count+=polygons[s].npt-1;
    nlimits++;
    }

  count+=ninteriors;

  mesh.nvtxs=count;
  mesh.nlimits=nlimits;

  mesh.vertices=new vertex_t[mesh.nvtxs];
  mesh.limits=new limit_t[mesh.nlimits];

  n=0;
  nlimits=0;
  for(s=0;s<polygons.size();s++) {
    if(polygons[s].npt==0) continue;
    mesh.limits[nlimits].nvertex=polygons[s].npt-1;
    mesh.limits[nlimits].vertex=new int[mesh.limits[nlimits].nvertex];
    for(k=0;k<polygons[s].npt-1;k++) {
      mesh.vertices[n].lon=polygons[s].x[k];
      mesh.vertices[n].lat=polygons[s].y[k];
      if(polygons[s].z!=0) mesh.vertices[n].h=polygons[s].z[k];
      mesh.vertices[n].nngh=0;
      mesh.vertices[n].code=nlimits+1;
      mesh.limits[nlimits].vertex[k]=n;
      n++;
      }
    nlimits++;
    }

  for(k=0;k<ninteriors;k++) {
    mesh.vertices[n].lon=x[k];
    mesh.vertices[n].lat=y[k];
    if(z!=0) mesh.vertices[n].h=z[k];
    mesh.vertices[n].code=0;
    mesh.vertices[n].nngh=0;
    n++;
    }

  return(0);

}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_createnodes(plg_t *polygones, int npolygones, double *x, double *y, int ninteriors, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
   
  vector<plg_t> p=plg_array2vector(polygones, npolygones);

  status=fe_createnodes(p, x, y, 0, ninteriors, mesh);
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_createnodes(const vector<plg_t> polygones, point_t *interiors, int ninteriors, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,j,k,l,m,n,status;
  int count,s, nlimits;

  count=0;
  nlimits=0;
  for(s=0;s<polygones.size();s++) {
    if(polygones[s].npt==0) continue;
    count+=polygones[s].npt-1;
    nlimits++;
    }

  count+=ninteriors;

  mesh.nvtxs=count;
  mesh.nlimits=nlimits;

  mesh.vertices=new vertex_t[mesh.nvtxs];
  mesh.limits=new limit_t[mesh.nlimits];

  n=0;
  nlimits=0;
  for(s=0;s<polygones.size();s++) {
    if(polygones[s].npt==0) continue;
    mesh.limits[nlimits].nvertex=polygones[s].npt-1;
    mesh.limits[nlimits].vertex=new int[mesh.limits[nlimits].nvertex];
    for(k=0;k<polygones[s].npt-1;k++) {
      mesh.vertices[n].lon=polygones[s].x[k];
      mesh.vertices[n].lat=polygones[s].y[k];
      mesh.vertices[n].nngh=0;
      mesh.vertices[n].code=nlimits+1;
      mesh.limits[nlimits].vertex[k]=n;
      n++;
      }
    nlimits++;
    }

  for(k=0;k<ninteriors;k++) {
    mesh.vertices[n].lon=interiors[k].t;
    mesh.vertices[n].lat=interiors[k].p;
    mesh.vertices[n].code=0;
    mesh.vertices[n].nngh=0;
    n++;
    }

  return(0);

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_limit_interior (mesh_t mesh,limit_t limit,double lon,double lat,int *in)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/*------------------------------------------------------------------------------

  M=A1 + l A1B1= A2 + l A2AB2

------------------------------------------------------------------------------*/
  int n,size,pos;
  int k,l,ok=0;

  double *x=NULL,*y=NULL,*dx=NULL,*dy=NULL,*dn=NULL,*delta=NULL;
  double epsilon,error,tx,ty,ang,nx,ny,ux,uy,xx,yy,vx,vy;
  double alpha,beta;
/*------------------------------------------------------------------------------
  count the number of intersections of a half-line starting at (lon,lat) with an
  arbitrary direction and the polygone's segments*/

  size=limit.nvertex+1;

  exitIfNull(
    x=new double[size]
    );
  exitIfNull(
    y=new double[size]
    );

  exitIfNull(
    dx=new double[size]
    );
  exitIfNull(
    dy=new double[size]
    );
  exitIfNull(
    dn=new double[size]
    );
  exitIfNull(
    delta=new double[size]
    );

  // error=1e-02; LR: 22/09/2006
  error=1.e-04;
  epsilon=1.e-12;

  for(k=0;k<size-1;k++) {
    n=limit.vertex[k];
    x[k]=mesh.vertices[n].lon;
    y[k]=mesh.vertices[n].lat;
    }
  n=limit.vertex[0];
  x[k]=mesh.vertices[n].lon;
  y[k]=mesh.vertices[n].lat;

  for(k=0;k<size-1;k++) {
    dx[k]=x[k+1]-x[k];
    dy[k]=y[k+1]-y[k];
    dn[k]=sqrt(dx[k]*dx[k]+dy[k]*dy[k]);
    if((dx[k]==0)&&(dy[k]==0)) {
      printf("erreur de polygone\n");
      }
    }

/*------------------------------------------------------------------------------
  seek a convenient arbitrary direction*/

  ang=0.0;
  tx=1.0;
  ty=0.0;
  nx=-ty;
  ny=+tx;

  for(l=0;l<size-1;l++) {
    delta[l]=nx*dx[l]+ny*dy[l];
    if (fabs(delta[l]) < error*dn[l]) {
      ang=ang+0.001;
      tx=cos(ang);
      ty=sin(ang);
      nx=-ty;
      ny=+tx;
      l=-1;
      }
    }
 
/*------------------------------------------------------------------------------
  screen the polygone's segment*/
  for(l=0;l<size-1;l++) {
      ux=lon-x[l];
      uy=lat-y[l];
      alpha=(nx*ux+ny*uy)/delta[l];
/*-------------------------------------------------------------------------------
      (P,t) is the arbitrary testing half line
      AB is one polygone segment
      B' is the ortogonal (along t) projection on (A,n) axis
      P is (lon,lat)
      P' is the ortogonal (along t) projection on (A,n) axis
      I is the intersection point between (P,t) and (A,AB)

      delta[l] is the B and B'  points' coordinate in (A,n) axis
      nx*ux+ny*uy is the P point coordinate in (A,n) axis
      alpha=(nx*ux+ny*uy)/delta[l] is the barycentric coordinate of the
      intersection point P' in (A,AB') axis
      There is an intersection if 0< alpha < 1
------------------------------------------------------------------------------*/
/*
      if (alpha == 0.) {
        pos=plg_point_boundary;
        retun(pos);
        }
*/
      if ((alpha >= 0.0)&&(alpha < 1.0)) {
/*------------------------------------------------------------------------------
        there is an intersection*/
/*------------------------------------------------------------------------------
        x,y is the intersection position*/
        xx=x[l]+alpha*dx[l];
        yy=y[l]+alpha*dy[l];
        vx=xx-lon;
        vy=yy-lat;
/*-------------------------------------------------------------------------------
        beta is the barycentric coordinate of the intersection point in (A,t) axis */
        beta=vx*tx+vy*ty;
        if(fabs(beta)<epsilon) {
          pos=POINT_BOUNDARY;
          goto end;
          }
        if(beta>0.) *in=(!(*in));
        }
     }

  if ((*in)==1) {
    pos=POINT_INTERIOR;
    }
  else {
    pos=POINT_EXTERIOR;
    }

 end:

  delete[] x;
  delete[] y;
  delete[] dx;
  delete[] dy;
  delete[] dn;
  delete[] delta;

  return(pos);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_holestable(mesh_t *mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  FILE *in=NULL;
  vertex_t *set=NULL;
  int i,j,k,l,nndes,position,ninteriors,nexteriors,nsegments;
  int n1,n2,n3,count,interior;
  double dx[2],dy[2],dn[2],alpha;

  set= mesh->vertices;
  nndes=mesh->nvtxs;

  if(set == NULL) {
    return(-1);
    }

  for (l=1; l<mesh->nlimits; l++) {
    if(mesh->limits[l].nvertex==3) {
      mesh->limits[l].x=0.0;
      mesh->limits[l].y=0.0;
      continue;
      }
    for (k=0; k<mesh->limits[l].nvertex-2; k++) {
      n1=mesh->limits[l].vertex[k];
      n2=mesh->limits[l].vertex[k+1];
      n3=mesh->limits[l].vertex[k+2];
      dx[0] = set[n1].lon-set[n2].lon;
      dy[0] = set[n1].lat-set[n2].lat;
      dx[1] = set[n3].lon-set[n2].lon;
      dy[1] = set[n3].lat-set[n2].lat;
      dn[0] = sqrt(dx[0]*dx[0]+dy[0]*dy[0]);
      dn[1] = sqrt(dx[1]*dx[1]+dy[1]*dy[1]);
      dx[0] /= dn[0];
      dy[0] /= dn[0];
      dx[1] /= dn[1];
      dy[1] /= dn[1];
      alpha=asin(dx[0]*dy[1]-dy[0]*dx[1])*r2d;
      mesh->limits[l].x=(set[n1].lon+set[n2].lon+set[n3].lon)/3.0;
      mesh->limits[l].y=(set[n1].lat+set[n2].lat+set[n3].lat)/3.0;
      interior=0;
      position=fe_limit_interior (*mesh,mesh->limits[l],mesh->limits[l].x,mesh->limits[l].y,&interior);
      if(position==POINT_INTERIOR) {
//        printf( "boundary #%d, legth #%d, pts #%d #%d, x=%lf y=%lf alpha=%lf\n", l, mesh->limits[l].nvertex, k,k+1,mesh->limits[l].x,mesh->limits[l].y,alpha);
        break;
        }
      }
    }

  count=0;
  for (l=1; l<mesh->nlimits; l++) {
    if(mesh->limits[l].nvertex>3) {
      count++;
      }
    }

  mesh->nholes=count;

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_node2triangle(mesh_t & mesh, triangulateio *out)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  convert tools node set into triangle node set

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int   k,l,m,i,j,n1,n2,n3,n,status,pass,t;
  int   count=0, mark, total, rotation;
  double t1,t2,p1,p2,d,xx,yy;
  double lon,lat,radius=6300.0,angle;

//  triangleinit(out);

  out->numberofpoints = mesh.nvtxs;
  out->numberofpointattributes = 1;
  exitIfNull(
//    out->pointlist=(REAL *) malloc(out->numberofpoints * 2 * sizeof(REAL))
    out->pointlist=new REAL[out->numberofpoints * 2]
    );

  for (n=0;n<mesh.nvtxs;n++) {
    out->pointlist[2*n]   =mesh.vertices[n].lon;
    out->pointlist[2*n+1] =mesh.vertices[n].lat;
    }

//   if ((out->pointattributelist = (REAL *) malloc(out->numberofpoints *
//                                                  out->numberofpointattributes *
//                                                  sizeof(REAL)))==NULL) {
  if ((out->pointattributelist = new REAL[out->numberofpoints *
                                                 out->numberofpointattributes])==NULL) {
    __ERR_BASE_LINE__("");perror("out->pointattributelist");
    exit(-1);
    }

  for (n=0;n<mesh.nvtxs;n++) {
    out->pointattributelist[n] = mesh.vertices[n].h;
    }

  exitIfNull(
//    out->pointmarkerlist=(int *) malloc(out->numberofpoints * sizeof(int))
    out->pointmarkerlist=new int[out->numberofpoints]
    );
  for (n=0;n<mesh.nvtxs;n++) {
    out->pointmarkerlist[n] = 0;
    }

  count=0;
  for (l=0; l<mesh.nlimits; l++) count+=mesh.limits[l].nvertex;

  out->numberofsegments = count;
  out->segmentmarkerlist=NULL;

  out->numberofholes = 0;
  out->holelist      =0;
  out->numberofregions  = 0;

  out->numberofedges = 0;
  
  out->normlist = 0;
  
  out->numberoftriangleattributes = 0;
  out->triangleattributelist = 0;
  
  out->trianglelist     = 0;
  out->trianglearealist = 0;
  out->neighborlist     = 0;
  out->segmentlist      = 0;
  out->edgelist         = 0;
  out->edgemarkerlist   = 0;
  
  if(count==0) return(0);
  
  out->segmentlist=new int[2*count];

  count=0;
  for (l=0; l<mesh.nlimits; l++) {
    for (k=0; k<mesh.limits[l].nvertex-1; k++) {
      i=mesh.limits[l].vertex[k];
      j=mesh.limits[l].vertex[k+1];
      out->segmentlist[2*count  ]=i;
      out->segmentlist[2*count+1]=j;
      count++;
      }
    i=mesh.limits[l].vertex[k];
    j=mesh.limits[l].vertex[0];
    out->segmentlist[2*count  ]=i;
    out->segmentlist[2*count+1]=j;
    count++;
    }

  status=fe_holestable(&mesh);

  out->numberofholes    = mesh.nholes;
  out->holelist=new double[2*mesh.nholes];
  n=0;
  for (l=1; l<mesh.nlimits; l++) {
    if(mesh.limits[l].nvertex>3) {
      out->holelist[2*n]   =mesh.limits[l].x;
      out->holelist[2*n+1] =mesh.limits[l].y;
      n++;
      }
    }
  out->numberofregions  = 0;

//  out->numberofregions  = 1;
//  out->regionlist = (REAL *) malloc(out->numberofregions * 4 * sizeof(REAL));
//  out->regionlist[0] = 0.5;
//  out->regionlist[1] = 5.0;
//  out->regionlist[2] = 7.0;            /* Regional attribute (for whole mesh). */
//  out->regionlist[3] = 0.1;            /* Area constraint that will not be used. */

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_triangle2mesh(triangulateio & in,mesh_t *mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  convert triangle mesh into tools mesh
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int   k,l,m,i,j,n1,n2,n3,n,status,pass,t;
  int   count=0, mark, total, rotation;
  double t1,t2,p1,p2,d,xx,yy;
  double lon,lat,radius=6300.0,angle;

  mesh->ntriangles=in.numberoftriangles;
  mesh->nvtxs=in.numberofpoints;
  
/*------------------------------------------------------------------------------
  allocate vertex arrays */
  bool *used=new bool[mesh->nvtxs];
  count=0;
  for (i=0; i<mesh->nvtxs; i++) used[i]=false;
  for (m=0;m<mesh->ntriangles;m++) {
    for(k=0;k<3;k++) {
      n=in.trianglelist[3*m+k];
      used[n]=true;
      }
    }
  for (i=0; i<mesh->nvtxs; i++) if(used[i]) count++;
  
  if(mesh->nvtxs-count!=0) printf("#fe_triangle2mesh : found %d unused vertices\n",mesh->nvtxs-count);
  int *ptr=new int[mesh->nvtxs];
  
  count=0;
  for (i=0; i<mesh->nvtxs; i++) {
    if(used[i]) {
      ptr[i]=count;
      count++;
      }
    }
  mesh->nvtxs=count;
  
/*------------------------------------------------------------------------------
  allocate element arrays */
  mesh->triangles=new triangle_t[mesh->ntriangles];
  if(mesh->triangles == NULL) {
    printf("allocation error element set\n");
    return(-1);
    }

  for (m=0;m<mesh->ntriangles;m++) {
    for(k=0;k<3;k++) {
      n=in.trianglelist[3*m+k];
      mesh->triangles[m].vertex[k]=ptr[n];
      }
    }

  delete[] in.trianglelist;
  in.trianglelist=0;
  
/*------------------------------------------------------------------------------
  create neighbour array from element array */
  mesh->vertices= new vertex_t[mesh->nvtxs];
  if(mesh->vertices == NULL) {
    printf("allocation error vertex set\n");
    return(-1);
    }
  for (i=0; i<mesh->nvtxs; i++) {
    mesh->vertices[i].null_value();
    }
    
  for (i=0; i<in.numberofpoints; i++) {
    if(!used[i]) continue;
    n=ptr[i];
    mesh->vertices[n].lon=in.pointlist[2*i];
    mesh->vertices[n].lat=in.pointlist[2*i+1];
    }
  
  delete[] in.pointlist;
  in.pointlist=0;

  status=fe_e2n(mesh);

  if((in.numberofpointattributes==1) && (in.pointattributelist!=0)) {
    for (i=0; i<in.numberofpoints; i++) {
      if(!used[i]) continue;
      n=ptr[i];
      mesh->vertices[n].h=in.pointattributelist[i];
      }
    }

  delete[] used;
  delete[] ptr;
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void fe_FreeTriangleIO(triangulateio & io)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
 
  deletep(&io.pointlist);                                               /* In / out */
  deletep(&io.pointattributelist);                                      /* In / out */
  deletep(&io.pointmarkerlist);                                          /* In / out */

  deletep(&io.trianglelist);                                             /* In / out */
  deletep(&io.triangleattributelist);                                   /* In / out */
  deletep(&io.trianglearealist);                                         /* In only */
  deletep(&io.neighborlist);                                             /* Out only */

  deletep(&io.segmentlist);                                              /* In / out */
  deletep(&io.segmentmarkerlist);                                        /* In / out */

  deletep(&io.holelist);                        /* In / pointer to array copied out */

  deletep(&io.regionlist);                      /* In / pointer to array copied out */

  deletep(&io.edgelist);                                                 /* Out only */
  deletep(&io.edgemarkerlist);            /* Not used with Voronoi diagram; out only */
  deletep(&io.normlist);                /* Used only with Voronoi diagram; out only */
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void fe_InitIO(triangulateio & io)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
 
  io.pointlist=0;                                               /* In / out */
  io.pointattributelist=0;                                      /* In / out */
  io.pointmarkerlist=0;                                          /* In / out */

  io.trianglelist=0;                                             /* In / out */
  io.triangleattributelist=0;                                   /* In / out */
  io.trianglearealist=0;                                         /* In only */
  io.neighborlist=0;                                             /* Out only */

  io.segmentlist=0;                                              /* In / out */
  io.segmentmarkerlist=0;                                        /* In / out */

  io.holelist=0;                        /* In / pointer to array copied out */

  io.regionlist=0;                      /* In / pointer to array copied out */

  io.edgelist=0;                                                 /* Out only */
  io.edgemarkerlist=0;            /* Not used with Voronoi diagram; out only */
  io.normlist=0;                /* Used only with Voronoi diagram; out only */
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  mesh_t *fe_triangulate(mesh_t & nodes, projPJ projection, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int   k,l,m,i,j,n1,n2,n3,n,status,pass;
  int   nelts,nndes;
  int   rotation;
  int   ring[4]={0,1,2,0};
  int   *ncells=NULL,ncellsmax=0,**cells=NULL;
  edge_t *edges=NULL;
  triangle_t *elt=NULL,*ptr=NULL;
  double t,p,x,y;
  mesh_t *finished=NULL;
  triangulateio in,out;
  FILE *chk;

  fe_InitIO(in);
  fe_InitIO(out);
  
#ifdef EXTERN_TRIANGLE
 status=fe_savenodes("refined.poly", NODE_FILE_FORMAT_TRIANGLE, nodes);
 if(status!=0) return(NULL);

 status=fe_savenodes("refined.nod", NODE_FILE_FORMAT_TRIGRID, nodes);
 if(status!=0) return(NULL);
 status=system("node2poly.exe refined.nod refined.poly");
 if(status!=0) return(NULL);

 status=system("triangle.exe -p refined.poly");
 if(status!=0) return(NULL);
#endif

/*------------------------------------------------------------------------------
  transfert node informations to triangle input structure */
  status= fe_node2triangle( nodes, &in);
  
  if((projection!=0) && debug) {
    chk=fopen("holes.xyz","w");
    fprintf(chk,"XYZ\n");
    for (n=0;n<nodes.nlimits;n++) {
      x=nodes.limits[n].x;
      y=nodes.limits[n].y;
      projection_to_geo(projection, &p, &t, x, y);
      fprintf(chk,"%lf %lf %f\n", t, p, 0.0);
      }
    fclose(chk);
    }
//  status=fe_savenodes("tmp.nod", NODE_FILE_FORMAT_TRIGRID, *mesh);

/*------------------------------------------------------------------------------
  initialize triangle output structure */
  out.pointlist = (REAL *) NULL;            /* Not needed if -N switch used. */
  out.pointattributelist = (REAL *) NULL;   /* Not needed if -N switch used or number of point attributes is zero: */
  out.pointmarkerlist = (int *) NULL;       /* Not needed if -N or -B switch used. */
  out.trianglelist = (int *) NULL;          /* Not needed if -E switch used. */
  out.triangleattributelist = (REAL *) NULL;/* Not needed if -E switch used or number of triangle attributes is zero: */
  out.neighborlist = (int *) NULL;          /* Needed only if -n switch used. */
  out.segmentlist = (int *) NULL;           /* Needed only if segments are output (-p or -c) and -P not used: */
  out.segmentmarkerlist = (int *) NULL;     /* Needed only if segments are output (-p or -c) and -P and -B not used: */
  out.edgelist = (int *) NULL;              /* Needed only if -e switch used. */
  out.edgemarkerlist = (int *) NULL;        /* Needed if -e used and -B not used. */

  if(verbose==0) triangulate("pzNQ", &in, &out, (triangulateio *) NULL);
  else triangulate("pzN", &in, &out, (triangulateio *) NULL);

/*------------------------------------------------------------------------------
  output structure inherit from input structure's pointlist */
  out.pointlist = in.pointlist;

  finished=new mesh_t;

#ifdef EXTERN_TRIANGLE
  status=fe_readmesh_TGL ("refined.1.ele","refined.1.node",finished);
#endif

  status=fe_triangle2mesh(out,finished);
  
  in.pointlist = 0;
  fe_FreeTriangleIO(in);
  out.holelist  = 0;
  fe_FreeTriangleIO(out);

  return(finished);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_triangulate(mesh_t & nodes, projPJ projection, mesh_t & finished, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int   k,l,m,i,j,n1,n2,n3,n,status,pass;
  int   nelts,nndes;
  int   rotation;
  int   ring[4]={0,1,2,0};
  int   *ncells=NULL,ncellsmax=0,**cells=NULL;
  edge_t *edges=NULL;
  triangle_t *elt=NULL,*ptr=NULL;
  double t,p,x,y;
  triangulateio in,out;
  FILE *chk;

  fe_InitIO(in);
  fe_InitIO(out);
  
#ifdef EXTERN_TRIANGLE
 status=fe_savenodes("refined.poly", NODE_FILE_FORMAT_TRIANGLE, nodes);
 if(status!=0) return(NULL);

 status=fe_savenodes("refined.nod", NODE_FILE_FORMAT_TRIGRID, nodes);
 if(status!=0) return(NULL);
 status=system("node2poly.exe refined.nod refined.poly");
 if(status!=0) return(NULL);

 status=system("triangle.exe -p refined.poly");
 if(status!=0) return(NULL);
#endif

/*------------------------------------------------------------------------------
  transfert node informations to triangle input structure */
  status= fe_node2triangle( nodes, &in);
  
  if((projection!=0) && debug) {
    chk=fopen("holes.xyz","w");
    fprintf(chk,"XYZ\n");
    for (n=0;n<nodes.nlimits;n++) {
      x=nodes.limits[n].x;
      y=nodes.limits[n].y;
      projection_to_geo(projection, &p, &t, x, y);
      fprintf(chk,"%lf %lf %f\n", t, p, 0.0);
      }
    fclose(chk);
    }
//  status=fe_savenodes("tmp.nod", NODE_FILE_FORMAT_TRIGRID, *mesh);

/*------------------------------------------------------------------------------
  initialize triangle output structure */
  out.pointlist = (REAL *) NULL;            /* Not needed if -N switch used. */
  out.pointattributelist = (REAL *) NULL;   /* Not needed if -N switch used or number of point attributes is zero: */
  out.pointmarkerlist = (int *) NULL;       /* Not needed if -N or -B switch used. */
  out.trianglelist = (int *) NULL;          /* Not needed if -E switch used. */
  out.triangleattributelist = (REAL *) NULL;/* Not needed if -E switch used or number of triangle attributes is zero: */
  out.neighborlist = (int *) NULL;          /* Needed only if -n switch used. */
  out.segmentlist = (int *) NULL;           /* Needed only if segments are output (-p or -c) and -P not used: */
  out.segmentmarkerlist = (int *) NULL;     /* Needed only if segments are output (-p or -c) and -P and -B not used: */
  out.edgelist = (int *) NULL;              /* Needed only if -e switch used. */
  out.edgemarkerlist = (int *) NULL;        /* Needed if -e used and -B not used. */

  if(verbose==0) triangulate("pzNQ", &in, &out, (triangulateio *) NULL);
  else triangulate("pzN", &in, &out, (triangulateio *) NULL);

/*------------------------------------------------------------------------------
  output structure inherit from input structure's pointlist */
  out.pointlist = in.pointlist;
  out.pointattributelist = in.pointattributelist;

#ifdef EXTERN_TRIANGLE
  status=fe_readmesh_TGL ("refined.1.ele","refined.1.node",&finished);
#endif

  status=fe_triangle2mesh(out,&finished);
  
  in.pointlist = 0;
  in.pointattributelist = 0;
  
  fe_FreeTriangleIO(in);
  out.holelist  = 0;
  fe_FreeTriangleIO(out);

  return(0);
}

