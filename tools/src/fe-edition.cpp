/*******************************************************************************

  T-UGO tools, 2006-2015

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

\author  Florent Lyard      LEGOS/CNRS, Toulouse, France. florent.lyard@legos.obs-mip.fr
\author  Laurent Roblou     LEGOS/CNRS, Toulouse, France
\author  Damien Allain      LEGOS/CNRS, Toulouse, France
\author  Yves Soufflet      LEGOS/CNRS, Toulouse, France
\author  Clément Mayet      LEGOS, Toulouse, France (PhD)
\author  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
\author  Frédéric Dupont    Université de Laval à Québec, Canada

\brief finite elements geometry edition functions
*/
/*----------------------------------------------------------------------------*/

#include <stdio.h>

#include "tools-structures.h"

#include "fe.def"
#include "fe.h"

#include "geo.h"


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  mesh_t fe_refine(mesh_t & mesh, int *selected, int limits_locked)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int   n1,n2,n;
  int   nedges,count=0;
  int   nelts,nndes;
  edge_t *edges=NULL;
  triangle_t *elt=NULL;
  double t1,t2,p1,p2;

  mesh_t work;

  elt=mesh.triangles;

  nelts =mesh.ntriangles;
  nndes =mesh.nvtxs;
  nedges=mesh.nedges;
  edges =mesh.edges;

  printf ("number of nodes (original mesh) %d\n",nndes);

  work.type=mesh.type;

  work.vertices=new vertex_t[nedges+nndes];
  for (n=0; n<nedges+nndes; n++) work.vertices[n].null_value();

  for (n=0;n<nndes;n++) {
    work.vertices[n].lon=mesh.vertices[n].lon;
    work.vertices[n].lat=mesh.vertices[n].lat;
    work.vertices[n].h=mesh.vertices[n].h;
    work.vertices[n].code=mesh.vertices[n].code;
    if(work.vertices[n].code==0) {
      work.vertices[n].nngh=0;
      }
    else {
      work.vertices[n].nngh=2;
      exitIfNull(
        work.vertices[n].ngh=new int[work.vertices[n].nngh]
        );
      work.vertices[n].ngh[0]=mesh.vertices[n].ngh[0];
      work.vertices[n].ngh[1]=mesh.vertices[n].ngh[mesh.vertices[n].nngh-1];
//       work.vertices[n].nngh=0;
//       for(k=0;k<mesh.vertices[n].nngh;k++) {
// 	m=mesh.vertices[n].ngh[k];
// 	if(mesh.vertices[m].code!=0) work.vertices[n].nngh++;
//         }
//       work.vertices[n].nngh=0;
//       for(k=0;k<mesh.vertices[n].nngh;k++) {
// 	m=mesh.vertices[n].ngh[k];
// 	if(mesh.vertices[m].code!=0) {
// 	  work.vertices[n].ngh[work.vertices[n].nngh]=m;
// 	  work.vertices[n].nngh++;
//           }
//         }
      }
    }

  count=nndes;
  if(limits_locked==1)  goto interior;

/*------------------------------------------------------------------------------
-
  boundary refinement */
  for (n=0;n<nedges;n++) {
    if(edges[n].nshared==2) continue;
    n1=edges[n].extremity[0];
    n2=edges[n].extremity[1];
    if (mesh.vertices[n1].code!=mesh.vertices[n2].code) continue;
/*     if (mesh.vertices[n1].code!=20) continue; */
    t1=mesh.vertices[n1].lon;
    t2=mesh.vertices[n2].lon;
    p1=mesh.vertices[n1].lat;
    p2=mesh.vertices[n2].lat;
    if(selected[n]) {
/*------------------------------------------------------------------------------
      create a new node at the middle of the edge */
      if(work.type==0) t2=geo_recale(t2,t1,(double) 180.0);
      work.vertices[count].lon=0.5*(t1+t2);
      work.vertices[count].lat=0.5*(p1+p2);
      work.vertices[count].h=0.5*(mesh.vertices[n1].h+mesh.vertices[n2].h);
      work.vertices[count].nngh=2;
      work.vertices[count].code=mesh.vertices[n1].code;
      exitIfNull(
        work.vertices[count].ngh=new int[2]
        );
      work.vertices[count].nngh=2;
/*------------------------------------------------------------------------------
      connect */
      if(work.vertices[n1].ngh[1]==n2) {
//        work.vertices[n2].ngh[0]=count+1;
//        work.vertices[n1].ngh[1]=count+1;
        work.vertices[n2].ngh[0]=count;
        work.vertices[n1].ngh[1]=count;
        work.vertices[count].ngh[0]=n1;
        work.vertices[count].ngh[1]=n2;
        }
      else if(work.vertices[n1].ngh[0]==n2) {
 //       work.vertices[n1].ngh[0]=count+1;
 //       work.vertices[n2].ngh[1]=count+1;
        work.vertices[n1].ngh[0]=count;
        work.vertices[n2].ngh[1]=count;
        work.vertices[count].ngh[0]=n2;
        work.vertices[count].ngh[1]=n1;
        }
      else {
        printf("troubles...\n");
        }

      count++;
      }
    }

/*------------------------------------------------------------------------------
-
  interior refinement */
interior:
  for (n=0;n<nedges;n++) {
    if(edges[n].nshared!=2) continue;
    n1=edges[n].extremity[0];
    n2=edges[n].extremity[1];
    t1=mesh.vertices[n1].lon;
    t2=mesh.vertices[n2].lon;
    p1=mesh.vertices[n1].lat;
    p2=mesh.vertices[n2].lat;
    if(selected[n]) {
/*------------------------------------------------------------------------------
      create a new node at the middle of the edge */
      if(mesh.type==0) t2=geo_recale(t2,t1,(double) 180.0);
      work.vertices[count].lon=0.5*(t1+t2);
      work.vertices[count].lat=0.5*(p1+p2);
      work.vertices[count].h=0.5*(mesh.vertices[n1].h+mesh.vertices[n2].h);
      work.vertices[count].code=0;
      work.vertices[count].nngh=0;
      count++;
      }
    }

// end:
  work.nvtxs=count;
  work.ntriangles=0;

  printf ("number of nodes (refined mesh) %d\n",count);
  return(work);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_reshape_spherical(mesh_t & mesh,int npass, int *targets, int ntargets)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int   k,l,m,n1,n2,n3,n,status,pass,t;
  int   nelts,nndes;
  size_t *touched;;
  double xx,yy,count;
  double lon,lat,lon_backup,lat_backup,angle;
  double *x=NULL,*y=NULL;

#define STEREO

#ifdef STEREO
  projPJ projection;
#else
  geo_t projection;
#endif
  
  nndes=mesh.nvtxs;
  nelts=mesh.ntriangles;

/*------------------------------------------------------------------------------
  work in cartesian coordinates */
  n=targets[0];
  lon=mesh.vertices[n].lon;
  lat=mesh.vertices[n].lat;
  for (n=0;n<nndes;n++) {
//  for (t=1;t<ntargets;t++) {
//    n=targets[t];
//     lon=mesh.vertices[n-1].lon;
//     lat=mesh.vertices[n-1].lat;
    if(mesh.type==SPHERICAL) mesh.vertices[n].lon = geo_recale(mesh.vertices[n].lon, lon, (double) 180.0);
    }

  lon=0;
  lat=0;
  for (t=0;t<ntargets;t++) {
    n=targets[t];
    lon+=mesh.vertices[n].lon;
    lat+=mesh.vertices[n].lat;
    }

  lon/=ntargets;
  lat/=ntargets;

#ifdef STEREO
  projection=assign_StereoOblique(lat,lon);
#else
  projection=geo_mercator_init(lon,lat,radius);
#endif

  exitIfNull(
    x=new double[nndes]
    );
  exitIfNull(
    y=new double[nndes]
    );

  for (n=0;n<nndes;n++) {
#ifdef STEREO
    geo_to_projection(projection,mesh.vertices[n].lat, mesh.vertices[n].lon, &x[n], &y[n]);
#else
    status=geo_mercator_directe(projection,mesh.vertices[n].lon,mesh.vertices[n].lat,&x[n],&y[n]);
#endif
    }

  touched=new size_t[mesh.ntriangles];
  for(m=0; m<nelts; m++) touched[m]=0;
    
  for(pass=0;pass<npass;pass++) {
//    printf("reshape loop: %d over %d\n",pass+1,npass);
    for (t=0;t<ntargets;t++) {
      n=targets[t];
      for(k=0;k<mesh.vertices[n].nelmts;k++) {
        m=mesh.vertices[n].elmts[k];
        touched[m]=1;
        }
      if(mesh.vertices[n].code==0) {
        count=0.;
        xx=0.0;
        yy=0.0;
        lon_backup=mesh.vertices[n].lon;
        lat_backup=mesh.vertices[n].lat;
        for(k=0;k<mesh.vertices[n].nngh;k++) {
          m=mesh.vertices[n].ngh[k];
          xx+=x[m];
          yy+=y[m];
          count++;
          }
        x[n]=xx/count;
        y[n]=yy/count;
#ifdef STEREO
        projection_to_geo(projection,&(mesh.vertices[n].lat),&(mesh.vertices[n].lon),x[n],y[n]);
#else
        status=geo_mercator_inverse(projection,&(mesh.vertices[n].lon),&(mesh.vertices[n].lat),x[n],y[n]);
#endif
/*------------------------------------------------------------------------------
        check elements good shape */
        for(k=0;k<mesh.vertices[n].nngh;k++) {
          n1=n;
          n2=mesh.vertices[n].ngh[k];
          if(k==mesh.vertices[n].nngh-1) l=0;
          else l=k+1;
          n3=mesh.vertices[n].ngh[l];
//          angle=fe_angle(mesh, n1,n2,n3);
          angle=fe_angle_cartesian(mesh, n1,n2,n3, x, y);
/*------------------------------------------------------------------------------
          negative angle means unsafe move for node n; if occures, cancel it*/
          if(angle  < 0.0) {
            printf("reshape alert: node=%d (neighbours %d %d) element=%d angle=%f, move cancelled\n",n1+1,n2+1,n3+1,k,angle);
            mesh.vertices[n].lon=lon_backup;
            mesh.vertices[n].lat=lat_backup;
#ifdef STEREO
            geo_to_projection(projection,mesh.vertices[n].lat, mesh.vertices[n].lon, &x[n], &y[n]);
#else
            status=geo_mercator_directe(projection,mesh.vertices[n].lon,mesh.vertices[n].lat,&x[n],&y[n]);
#endif
            break;
            }
          }
        }
      }
    }

  delete[] x;
  delete[] y;
  
/*------------------------------------------------------------------------------
  final check element good shape */
  for(m=0; m<nelts; m++) {
    if(touched[m]==0) continue;
    status=fe_initaffine(&mesh,m);
    if(mesh.triangles[m].Area <=0.0) {
      printf("element %d has negavtive area (cw order)\n",m);
      }
    }
    
  delete[] touched;

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_reshape_cartesian(mesh_t mesh,int npass)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int   k,l,m,n1,n2,n3,n,status,pass;
  int   nelts,nndes;
  double xx,yy,count;
  double lon,lat,angle;
  double *x=NULL,*y=NULL;

  nndes=mesh.nvtxs;
  nelts=mesh.ntriangles;

/*------------------------------------------------------------------------------
  work in cartesian coordinates */
  x=new double[nndes];
  y=new double[nndes];
  for (n=0;n<nndes;n++) {
    x[n]=mesh.vertices[n].lon;
    y[n]=mesh.vertices[n].lat;
    }

  for(pass=0;pass<npass;pass++) {
//    printf("reshape loop: %d over %d\n",pass+1,npass);
    for (n=0;n<nndes;n++) {
      if(mesh.vertices[n].code==0) {
        count=0.;
        xx=0.0;
        yy=0.0;
        lon=mesh.vertices[n].lon;
        lat=mesh.vertices[n].lat;
/*
        count=1.;
        xx=x[n];
        yy=y[n];
*/
        for(k=0;k<mesh.vertices[n].nngh;k++) {
          m=mesh.vertices[n].ngh[k];
          xx+=x[m];
          yy+=y[m];
          count++;
          }
        x[n]=xx/count;
        y[n]=yy/count;
        mesh.vertices[n].lon=x[n];
        mesh.vertices[n].lat=y[n];
/*------------------------------------------------------------------------------
        check elements good shape */
        for(k=0;k<mesh.vertices[n].nngh;k++) {
          n1=n;
          n2=mesh.vertices[n].ngh[k];
          if(k==mesh.vertices[n].nngh-1) l=0;
          else l=k+1;
          n3=mesh.vertices[n].ngh[l];
          angle=fe_angle_NoRecale(mesh, n1,n2,n3);
/*------------------------------------------------------------------------------
          negative angle means unsafe move for node n; if occures, cancel it*/
          if(angle  < 0.0) {
//            printf("reshape alert: node=%d (neighbours %d %d) element=%d angle=%f, move cancelled\n",n1+1,n2+1,n3+1,k,angle);
            mesh.vertices[n].lon=lon;
            mesh.vertices[n].lat=lat;
            break;
            }
          }
        }
      }
    }

  delete[] x;
  delete[] y;

/*------------------------------------------------------------------------------
  final check element good shape */
  for(n=0; n<nelts; n++) {
    status=fe_initaffine(&mesh,n);
    if(mesh.triangles[n].Area <=0.0) {
      printf("element %d has negavtive area (cw order)\n",n);
      }
   }

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_reshapeall01(mesh_t mesh,int npass)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int   k,l,m,n1,n2,n3,n,status,pass;
  int   nelts,nndes;
  int   rotation;
  double xx,yy,count;
  double lon,lat,radius=6300.0,angle;
  double *x=NULL,*y=NULL;
  geo_t projection;

  nndes=mesh.nvtxs;
  nelts=mesh.ntriangles;

/*------------------------------------------------------------------------------
  work in cartesian coordinates */
  lon=0;
  lat=0;
  for (n=0;n<nndes;n++) {
    lon+=mesh.vertices[n].lon;
    lat+=mesh.vertices[n].lat;
    }

  lon/=nndes;
  lat/=nndes;

  projection=geo_mercator_init(lon,lat,radius);

  exitIfNull(
    x=(double *) malloc(nndes*sizeof(double))
    );
  exitIfNull(
    y=(double *) malloc(nndes*sizeof(double))
    );

  for (n=0;n<nndes;n++) {
    status=geo_mercator_directe(projection,mesh.vertices[n].lon,mesh.vertices[n].lat,&x[n],&y[n]);
    }

  for(pass=0;pass<npass;pass++) {
//    printf("reshape loop: %d over %d\n",pass+1,npass);
    for (n=0;n<nndes;n++) {
      if(mesh.vertices[n].code==0) {
        count=0.;
        xx=0.0;
        yy=0.0;
        lon=mesh.vertices[n].lon;
        lat=mesh.vertices[n].lat;
/*
        count=1.;
        xx=x[n];
        yy=y[n];
*/
        for(k=0;k<mesh.vertices[n].nngh;k++) {
          m=mesh.vertices[n].ngh[k];
          xx+=x[m];
          yy+=y[m];
          count++;
          }
        x[n]=xx/count;
        y[n]=yy/count;
        status=geo_mercator_inverse(projection,&(mesh.vertices[n].lon),&(mesh.vertices[n].lat),x[n],y[n]);
/*------------------------------------------------------------------------------
        check elements good shape */
        for(k=0;k<mesh.vertices[n].nngh;k++) {
          n1=n;
          n2=mesh.vertices[n].ngh[k];
          if(k==mesh.vertices[n].nngh-1) l=0;
          else l=k+1;
          n3=mesh.vertices[n].ngh[l];
          angle=fe_angle(mesh, n1,n2,n3);
/*------------------------------------------------------------------------------
          negative angle means unsafe move for node n; if occures, cancel it*/
          if(angle  < 0.0) {
/*             printf("reshape alert: %d %d %d %d %f\n",n1+1,n2+1,n3+1,k,angle); */
            printf("reshape alert: node=%d (neighbours %d %d) element=%d angle=%f, move cancelled\n",n1+1,n2+1,n3+1,k,angle);
            mesh.vertices[n].lon=lon;
            mesh.vertices[n].lat=lat;
            status=geo_mercator_directe(projection,mesh.vertices[n].lon,mesh.vertices[n].lat,&x[n],&y[n]);
            }
          }
        }
      }
    }
/*
  for (n=0;n<nndes;n++) {
    if(mesh.vertices[n].code==0) {
      status=geo_mercator_inverse(projection,&(mesh.vertices[n].lon),&(mesh.vertices[n].lat),x[n],y[n]);
      }
    }
*/
  free(x);
  free(y);

/*------------------------------------------------------------------------------
  final check element good shape */
  for(n=0; n<nelts; n++) {
    status=fe_initaffine(&mesh,n);
    if(mesh.triangles[n].Area <=0.0) {
      printf("element %d has negative area (cw order)\n",n);
      }
   }

/*------------------------------------------------------------------------------
  need to reorder neighbours */
  for(n=0; n<nndes; n++) {
    rotation=+1;
    status=fe_order(mesh,n,(double) -1., (double) 0.,rotation);
    }

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_reshapeall02_spherical(mesh_t & mesh, int npass)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int   n,status;
  int   nndes;
  double lonmin,lonmax,latmin,latmax,lon_step,lat_step;
  int rotation;
  double lon,lat;
  int *targets=NULL,ntargets;

  nndes=mesh.nvtxs;

/*------------------------------------------------------------------------------
  need to order neighbours */
  for(n=0; n<nndes; n++) {
    rotation=+1;
    status=fe_order(mesh,n,(double) -1., (double) 0.,rotation);
    }

/*------------------------------------------------------------------------------
  seek min/max of nodes positions */

  lonmin=+1.e+10;
  lonmax=-1.e+10;
  latmin=+1.e+10;
  latmax=-1.e+10;

  for (n=0;n<nndes;n++) {
    lon=geo_recale(mesh.vertices[n].lon,0.,180.);
    updatemin(&lonmin,lon);
    updatemin(&latmin,mesh.vertices[n].lat);
    updatemax(&lonmax,lon);
    updatemax(&latmax,mesh.vertices[n].lat);
    }

/*------------------------------------------------------------------------------
  square repartition 10x10 for cartesian projection*/
  lon_step=30.0;
  lat_step=10.0;
  lat_step=30.0;

  lonmin=floor (lonmin/lon_step) *lon_step;
  lonmax=ceil  (lonmax/lon_step) *lon_step;
  latmin=floor (latmin/lat_step) *lat_step;
  latmax=ceil  (latmax/lat_step) *lat_step;

  targets=new int[nndes];

  for (lon=lonmin;lon<lonmax;lon+=lon_step) {
    for (lat=latmin;lat<latmax;lat+=lat_step) {
      ntargets=0;
      for (n=0;n<nndes;n++) {
        mesh.vertices[n].lon = geo_recale(mesh.vertices[n].lon, lon, (double) 180.0);
        if ( (lon<mesh.vertices[n].lon) && (mesh.vertices[n].lon<=lon+lon_step) &&
             (lat<mesh.vertices[n].lat) && (mesh.vertices[n].lat<=lat+lat_step) ) {
          targets[ntargets]=n;
          ntargets++;
          }
        }
      if(ntargets!=0) status=fe_reshape_spherical(mesh, npass, targets, ntargets);
      }
    }
  delete[] targets;

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_reshapeall02_cartesian(mesh_t & mesh,int npass)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int   n,status;
  int   nelts,nndes;
  int rotation;

  nndes=mesh.nvtxs;
  nelts=mesh.ntriangles;

/*------------------------------------------------------------------------------
  need to order neighbours */
  for(n=0; n<nndes; n++) {
    rotation=+1;
    status=fe_order(mesh,n,(double) -1., (double) 0.,rotation);
    }

  status=fe_reshape_cartesian(mesh, npass);

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_reshapeall(mesh_t & mesh,int npass)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;

  switch (mesh.type) {
    case SPHERICAL:
      status=fe_reshapeall02_spherical( mesh, npass);
      break;

    case CARTESIAN:
      status=fe_reshapeall02_cartesian( mesh, npass);
      break;
    }

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_reshape(mesh_t & mesh, int npass, vector<int> vtargets)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int *targets, ntargets;
  
  targets=copy(vtargets);
  ntargets=vtargets.size();
  
  switch (mesh.type) {
    case SPHERICAL:
      status=fe_reshape_spherical(mesh, npass, targets, ntargets);
      break;

    case CARTESIAN:
      status=fe_reshape_cartesian(mesh, npass);  /// HERE !!!
      break;
    }

  if(ntargets!=0) delete[] targets;
  
  return(status);
}

//###########################################################################################################################
//###########################################################################################################################

                                //cleave functions .....

//###########################################################################################################################
//###########################################################################################################################


//this function is static and can not be call out of the following functions -----------------------------
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  static int fe_cleave_heart(mesh_t *mesh, int ndel, int & n1,int & n2, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/*------------------------------------------------------------------------------
     Purpose:  To return indices, nos. of neighbours, indices of neighbours,
               coordinates of neighbours, of the two new vertices replacing
               index NDX in cleave operation.
     delndx  - index of vertex to be replaced
     numnbrs - no. of neighbours of delndx (old vertex)
     nbarray - neighbours of delndx sorted counterclockwise

-----------------------------------------------------------------------------*/
  int   kk,k,k1,k2,kk1,kk2,l,m,m1,m2,n,status;
  int   nndes,rotation=1;
  int   ngh[2][100],nngh[2];
  double lon,lat,radius=6300.0;
  double d,dmin=1.e+06,x,y,xx[2],yy[2];
  float z;
  geo_t projection;
  
  if(debug) {
    m=ndel;
    printf("%s : cleaved vertex %6d code=%2d (nvertices=%d)\n", __func__, m, mesh->vertices[m].code, mesh->nvtxs);
    printf("neighbours: ");
    for(k=0;k<mesh->vertices[m].nngh;k++) printf(" %d",mesh->vertices[m].ngh[k]);
    printf("\n");
    }
    
  nndes=mesh->nvtxs;
  for(k=0;k<100;k++) {
    ngh[0][k]=-1;
    ngh[1][k]=-1;
    }

/*------------------------------------------------------------------------------
  work in local cartesian coordinates */
  lon=0;
  lat=0;
  for(k=0;k<mesh->vertices[ndel].nngh;k++) {
    n=mesh->vertices[ndel].ngh[k];
    mesh->vertices[n].lon=geo_recale(mesh->vertices[n].lon,mesh->vertices[ndel].lon,180.0);
    lon+=mesh->vertices[n].lon;
    lat+=mesh->vertices[n].lat;
    }

  lon/=mesh->vertices[ndel].nngh;
  lat/=mesh->vertices[ndel].nngh;

  projection=geo_mercator_init(lon,lat,radius);

/*------------------------------------------------------------------------------
  in this version, forbid cleave procedure for nodes with 4 or less neighbours */
  if(mesh->vertices[ndel].nngh < 5) return(-1);

/*------------------------------------------------------------------------------
  sort neighbours of old vertex (delndx) into ccw order from east*/
  status=fe_order(*mesh,ndel,(double) -1., (double) 0.,rotation);

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  find opposed pair of neighbours minimum distance apart
   
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  dmin=1.e+10;
  for(k=0;k<mesh->vertices[ndel].nngh;k++) {
    m1=mesh->vertices[ndel].ngh[k];
    m2=fe_oppositeneighbour(*mesh, ndel,m1);
    d=fe_distance(*mesh,m1,m2);
    if(d<dmin) {
      dmin=d;
      kk1=k;
      kk2=fe_anb(*mesh, m2, ndel);
      }
    }
    
  k1=(int) min(kk1,kk2);
  k2=(int) max(kk1,kk2);

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

   note : the two new vertices replacing the old vertex (delndx) both have the
   closest pair of opposed vertices i.e. final values of n1, n2 above, as
   neighbours. i.e.
   
   nb2(1) = nb1(1) = best n1
   nb2(2) = nb1(2) = best n2
   
   allocate neighbours and their coordinates to new vertices
   
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  kk=0;
  for(k=k1;k<=k2;k++) {
    m=mesh->vertices[ndel].ngh[k];
    ngh[0][kk]=m;
    kk++;
    }
  nngh[0]=kk;

  kk=0;
  for(k=k2;k<=k1+mesh->vertices[ndel].nngh;k++) {
    l=k % mesh->vertices[ndel].nngh;
    m=mesh->vertices[ndel].ngh[l];
    ngh[1][kk]=m;
    kk++;
    }
  nngh[1]=kk;

  for(l=0;l<2;l++) {
    xx[l]=0.0;
    yy[l]=0.0;
    for(k=0;k<nngh[l];k++) {
      m=ngh[l][k];
      status=geo_mercator_directe(projection,mesh->vertices[m].lon,mesh->vertices[m].lat,&x,&y);
      xx[l]+=x;
      yy[l]+=y;
      }
    xx[l]/=nngh[l];
    yy[l]/=nngh[l];
    }

  nngh[0]++;
  nngh[1]++;

  if(nngh[0]==4) {
//    printf("fe_disconnectvertex anomaly for vertex %d\n",ndel);
    return(-1);
    }

  if(nngh[1]==4) {
//    printf("fe_disconnectvertex anomaly for vertex %d\n",ndel);
    return(-1);
    }

/**need to give a z value, to be done...*/
  z=9999.999;

  status=fe_ReallocateVertices(mesh, 2);
  
  status=geo_mercator_inverse(projection, &lon,&lat,xx[0],yy[0]);
  n1=fe_setvertex(mesh, lon, lat, z,ngh[0], nngh[0], mesh->nvtxs-2);
  mesh->vertices[n1].code=0;

  status=geo_mercator_inverse(projection, &lon,&lat,xx[1],yy[1]);
  n2=fe_setvertex(mesh, lon, lat, z,ngh[1], nngh[1], mesh->nvtxs-1);
  mesh->vertices[n2].code=0;

  mesh->vertices[n1].ngh[nngh[0]-1]= n2;
  mesh->vertices[n2].ngh[nngh[1]-1]= n1;

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  control future new triangle sanity (positive area) 
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  m=n1;
  for(k=0;k<mesh->vertices[m].nngh;k++) {
    int l=(k+1) % mesh->vertices[m].nngh;
    int m1=mesh->vertices[m].ngh[k];
    int m2=mesh->vertices[m].ngh[l];
    double angle=fe_angle(*mesh, m,m1,m2);
    if(angle <0) {
      printf("negative area, cleave of vertex %d cancelled\n",m);
      status=fe_removevertex(mesh,n1);
      status=fe_removevertex(mesh,n2);
      return(-1);
      }
    }
  m=n2;
  for(k=0;k<mesh->vertices[m].nngh;k++) {
    int l=(k+1) % mesh->vertices[m].nngh;
    int m1=mesh->vertices[m].ngh[k];
    int m2=mesh->vertices[m].ngh[l];
    double angle=fe_angle(*mesh, m,m1,m2);
    if(angle <0) {
      printf("negative area, cleave of vertex %d cancelled\n",m);
      status=fe_removevertex(mesh,n1);
      status=fe_removevertex(mesh,n2);
      return(-1);
      }
    }


  status=fe_insertvertex(mesh,n1);
  status=fe_insertvertex(mesh,n2);

//  status=fe_removevertex(mesh,ndel);
  status=fe_disconnectvertex(*mesh,ndel);
  if(status==-1) {
    printf("fe_disconnectvertex anomaly for vertex %d\n",ndel);
    }
    
  if(debug) {
    m=n1;
    printf("%s : new vertex %6d code=%2d (nvertices=%d)\n", __func__, m, mesh->vertices[m].code, mesh->nvtxs);
    printf("neighbours: ");
    for(k=0;k<mesh->vertices[m].nngh;k++) printf(" %d",mesh->vertices[m].ngh[k]);
    printf("\n");
    m=n2;
    printf("%s : new vertex %6d code=%2d (nvertices=%d)\n", __func__, m, mesh->vertices[m].code, mesh->nvtxs);
    printf("neighbours: ");
    for(k=0;k<mesh->vertices[m].nngh;k++) printf(" %d",mesh->vertices[m].ngh[k]);
    printf("\n");
    }

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_cleave(mesh_t *mesh, int ndel, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
 int n1,n2;
 int status;
 
 status=fe_cleave_heart(mesh,ndel,n1,n2, debug);
 
 return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_cleave(mesh_t *mesh, int ndel,int & n1,int & n2, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
 int status;
 status=fe_cleave_heart(mesh,ndel,n1,n2, debug);
 return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  static int fe_cleavelimit_heart(mesh_t *mesh, int ndel, int *n1,int *n2)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/*------------------------------------------------------------------------------
     Purpose:  To return indices, nos. of neighbours, indices of neighbours,
               coordinates of neighbours, of the two new vertices replacing
               index NDX in cleave operation.
     delndx  - index of vertex to be replaced
     numnbrs - no. of neighbours of delndx (old vertex)
     nbarray - neighbours of delndx sorted counterclockwise

-----------------------------------------------------------------------------*/
  const int nnmax=50;
  int   kk,k,k1,k2,l,m,m1,m2,n,status;
  int   nndes,rotation=1;
  int   count,ngh[2][nnmax],nngh[2];
  double lon,lat,radius=6300.0;
  double dmin=1.e+06,x,y,xx[2],yy[2];
  double alpha=0.05;
  float z;
  geo_t projection;

  nndes=mesh->nvtxs;
  for(k=0;k<nnmax;k++) {
    ngh[0][k]=-1;
    ngh[1][k]=-1;
    }


/*------------------------------------------------------------------------------
  work in local cartesian coordinates */
  lon=0;
  lat=0;
  for(k=0;k<mesh->vertices[ndel].nngh;k++) {
    n=mesh->vertices[ndel].ngh[k];
    lon+=mesh->vertices[n].lon;
    lat+=mesh->vertices[n].lat;
    }

  lon/=mesh->vertices[ndel].nngh;
  lat/=mesh->vertices[ndel].nngh;

  projection=geo_mercator_init(lon,lat,radius);

/*------------------------------------------------------------------------------
  sort neighbours of old vertex (delndx) into ccw order from east*/
  status=fe_order(*mesh,ndel,(double) -1., (double) 0.,rotation);

/*------------------------------------------------------------------------------
  find neighbours of ndel not neighbours for themselves*/
  dmin=1.e+10;
  for(k=0;k<mesh->vertices[ndel].nngh;k++) {
    m1=mesh->vertices[ndel].ngh[k];
    l=(k+1)%mesh->vertices[ndel].nngh;
    m2=mesh->vertices[ndel].ngh[l];
    if(fe_anb(*mesh, m1, m2)==-1) {
      k1=l;
      break;
      }
    }

  for(kk=l;kk<mesh->vertices[ndel].nngh+l;kk++) {
    k=kk%mesh->vertices[ndel].nngh;
    m1=mesh->vertices[ndel].ngh[k];
    l=(k+1)%mesh->vertices[ndel].nngh;
    m2=mesh->vertices[ndel].ngh[l];
    if(fe_anb(*mesh, m1, m2)==-1) {
      k2=kk;
      }
    }

/*------------------------------------------------------------------------------
   note : the two new vertices replacing the old vertex (delndx) both have the
   closest pair of opposed vertices i.e. final values of nt1, nt2 above, as
   neighbours. i.e.
   nb2(1) = nb1(1) = best nt1
   nb2(2) = nb1(2) = best nt2
   allot neighbours and their coordinates to new vertices
*/


  count=0;
  for(k=k1;k<=k2;k++) {
    l=k%mesh->vertices[ndel].nngh;
    m=mesh->vertices[ndel].ngh[l];
    ngh[0][count]=m;
    count++;
    }
  nngh[0]=count;

  count=0;
  for(k=k2+1;k<=k1-1+mesh->vertices[ndel].nngh;k++) {
    l=k%mesh->vertices[ndel].nngh;
    m=mesh->vertices[ndel].ngh[l];
    ngh[1][count]=m;
    count++;
    }
  nngh[1]=count;

  for(l=0;l<2;l++) {
    xx[l]=0.0;
    yy[l]=0.0;
    for(k=0;k<nngh[l];k++) {
      m=ngh[l][k];
      status=geo_mercator_directe(projection,mesh->vertices[m].lon,mesh->vertices[m].lat,&x,&y);
      xx[l]+=x;
      yy[l]+=y;
      }
    xx[l]/=nngh[l];
    yy[l]/=nngh[l];
    }

//   nngh[0]++;
//   nngh[1]++;

  /*need to give a z value, to be done...*/
  z=9999.999; //use a non value z than non assigned value thierry le 26/10/2005

  status=geo_mercator_directe(projection,mesh->vertices[ndel].lon,mesh->vertices[ndel].lat,&x,&y);

  xx[0]=alpha*xx[0]+(1.0-alpha)*x;
  yy[0]=alpha*yy[0]+(1.0-alpha)*y;

  status=geo_mercator_inverse(projection, &lon,&lat,xx[0],yy[0]);
  *n1=fe_addvertex3D(mesh, lon, lat, mesh->vertices[ndel].h, mesh->vertices[ndel].zlevels, ngh[0], nngh[0]);

  xx[1]=alpha*xx[1]+(1.0-alpha)*x;
  yy[1]=alpha*yy[1]+(1.0-alpha)*y;

  status=geo_mercator_inverse(projection, &lon,&lat,xx[1],yy[1]);
  *n2=fe_addvertex3D(mesh, lon, lat, mesh->vertices[ndel].h, mesh->vertices[ndel].zlevels, ngh[1], nngh[1]);

//   mesh->vertices[*n1].ngh[nngh[0]-1]= *n2;
//   mesh->vertices[*n2].ngh[nngh[1]-1]= *n1;

  status=fe_insertvertex(mesh,*n1);
  status=fe_insertvertex(mesh,*n2);

  status=fe_removevertex(mesh,ndel);

  if(*n1>ndel) (*n1)--;
  if(*n2>ndel) (*n2)--;

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_cleavelimit(mesh_t *mesh, int ndel)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Development notes

  Purpose : fix pinched boundary problems

  Check :

  Note:


@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
{
  int n1,n2;
  int status;


  status=fe_savemesh("cleavelimit-00.nei",MESH_FILE_FORMAT_TRIGRID,*mesh);
  status=fe_cleavelimit_heart(mesh,ndel,&n1,&n2);
//  if(status==0) status=fe_disconnectvertices(*mesh, n1,n2);
//  status=-1;
  status=fe_list(mesh);
  status=fe_savemesh("cleavelimit-01.nei",MESH_FILE_FORMAT_TRIGRID,*mesh);

  status=fe_edgetable(mesh,0,0);
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_AddVertexOnEdge(mesh_t *mesh, int edge)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  add vertex on edge n and connect it
  
  required uptodate triangle structure
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int    k,n,status;
  int    n1,n2,nngh;
  double lon, lat;
  float  z;
  size_t m,nn;
  
/*------------------------------------------------------------------------------
  resize vertices array */
  n=fe_ReallocateVertices(mesh, 1);

/*------------------------------------------------------------------------------
  create new vertex */
  n1=mesh->edges[edge].extremity[0];
  n2=mesh->edges[edge].extremity[1];
  
  lon=0.5*(mesh->vertices[n1].lon+mesh->vertices[n2].lon);
  lat=0.5*(mesh->vertices[n1].lat+mesh->vertices[n2].lat);
  z=0.5*(mesh->vertices[n1].h+mesh->vertices[n2].h);
  
  mesh->vertices[n].lon=lon;
  mesh->vertices[n].lat=lat;
  mesh->vertices[n].h=z;
  
/*------------------------------------------------------------------------------
  connect with the 2 side vertices*/
  nngh=2;
  
  mesh->vertices[n].ngh=new int[nngh];
  mesh->vertices[n].nngh=nngh;
  
  mesh->vertices[n].ngh[0]=n1;
  mesh->vertices[n].ngh[1]=n2;

/*------------------------------------------------------------------------------
  connect with the opposite vertices*/
  for(k=0;k<mesh->edges[edge].nshared;k++) {
    m=mesh->edges[edge].shared[k];
    int pos=vpos(edge,mesh->triangles[m].edges,3);
    nn=mesh->triangles[m].vertex[pos];
//    mesh->vertices[n].ngh[k+2]=mesh->triangles[m].vertex[pos];
    status=fe_connectvertices(*mesh, n, nn);
    }

/*------------------------------------------------------------------------------
  substitue the 2 side vertices in neighbours list */
  status=fe_substitute_vertex(*mesh, n1, n2, n);
  
  return(n);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_AddVerticesOnEdge(mesh_t *mesh, int edge, int nvertices)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------

  add vertices on edge n and connect them, to be finalized*/
{
  int    k,n,status;
  int    n1,n2,nngh;
  double lon, lat, L;
  float  z;
  size_t m,nn,start;
  vector<int> nodes;
  
/*------------------------------------------------------------------------------

  resize vertices array */
  start=mesh->nvtxs;
  n=fe_ReallocateVertices(mesh, nvertices);

  n1=mesh->edges[edge].extremity[0];
  n2=mesh->edges[edge].extremity[1];
  
  L=mesh->edges[edge].L;
  
  nodes.push_back(n1);
  for(int k=0;k<nvertices;k++) nodes.push_back(start+k);
  nodes.push_back(n2);
  
  for(int k=0;k<nvertices;k++) {
    double f=(double) k/(double) (nvertices+1);
    lon=mesh->vertices[n1].lon+f*(mesh->vertices[n2].lon-mesh->vertices[n1].lon);
    lat=mesh->vertices[n1].lat+f*(mesh->vertices[n2].lat-mesh->vertices[n1].lat);
    z=  mesh->vertices[n1].h+  f*(mesh->vertices[n2].h  -mesh->vertices[n1].h);
  
    mesh->vertices[n].lon=lon;
    mesh->vertices[n].lat=lat;
    mesh->vertices[n].h  =z;
    
    nngh=2;
  
    mesh->vertices[n].ngh=new int[nngh];
    mesh->vertices[n].nngh=nngh;
  
    mesh->vertices[n].ngh[0]=nodes[k];
    mesh->vertices[n].ngh[1]=nodes[k+1];
    }

  for(k=0;k<mesh->edges[edge].nshared;k++) {
    m=mesh->edges[edge].shared[k];
    int pos=vpos(edge,mesh->triangles[m].edges,3);
    nn=mesh->triangles[m].vertex[pos];
    for(int l=0;l<nvertices;l++) {
      status=fe_connectvertices(*mesh, nodes[k+1], nn);
      }
    }

/**----------------------------------------------------------------------------
  the 2 vertices */
  status=fe_substitute_vertex(*mesh, n1, n2, n);
  
  return(n);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_presplit(mesh_t & mesh, vector<plg_t> & polygons, double threshold, string rootname)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  Purpose : align mesh nodes on polygon's segments before splitting

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int p,n,status;
  mesh_t finished;
  bool *movable, completed;
  double x, y;
  string output;
  size_t count=0;
  
  if(rootname=="") rootname="split-adjusted";
  
  output=(string) rootname+(string) ".nei";
    
  for(n=0;n<mesh.nvtxs;n++) {
    mesh.vertices[n].code=-1;
    }
  start:
  
  movable=new bool[mesh.nvtxs];
  for(n=0;n<mesh.nvtxs;n++) {
    movable[n]=(mesh.vertices[n].code!=0);
    }
  
  for(p=0;p<polygons.size();p++) {
    for(size_t k=0;k<polygons[p].npt-1;k++) {
      line_t a=line_t(polygons[p],k);
      double longtitude_reference=0.5*(a.point[0].t+a.point[1].t);
      for(n=0;n<mesh.nedges;n++) {
        int n1=mesh.edges[n].extremity[0];
        if(!movable[n1]) continue;
        int n2=mesh.edges[n].extremity[1];
        if(!movable[n2]) continue;
        vertex_t v1=mesh.vertices[n1];
        vertex_t v2=mesh.vertices[n2];
        v1.lon=geo_recale(v1.lon, longtitude_reference, 180.0);
        v2.lon=geo_recale(v2.lon, longtitude_reference, 180.0);
        line_t b=line_t(v1.lon,v1.lat,v2.lon,v2.lat);
        int flag=plg_secantpoint(a,b,&x,&y);
        if(flag!=PLG_LINES_SECANT) continue;
        completed=false;
        double L=mesh.edges[n].L;
        L=1000.*geo_km(v1.lon,v1.lat,v2.lon,v2.lat);
        double d1=1000.*geo_km(v1.lon,v1.lat,x,y);
        double d2=1000.*geo_km(v2.lon,v2.lat,x,y);
        if (d1<threshold*L) {
/**----------------------------------------------------------------------------
          move vertex n1 */
          status=fe_displacevertex(mesh, n1, x, y);
          if(status==0) {
            completed=true;
            movable[n1]=false;
//             movable[n2]=false;
            mesh.vertices[n1].code=0;
            }
          }
        if(completed) continue;
        if (d2<threshold*L) {
/**----------------------------------------------------------------------------
          move vertex n2 */
          status=fe_displacevertex(mesh, n2, x, y);
          if(status==0) {
            completed=true;
//             movable[n1]=false;
            movable[n2]=false;
            mesh.vertices[n2].code=0;
            }
          }
        if(completed) continue;
/**----------------------------------------------------------------------------
        add new vertex on edge n */
        size_t node=fe_AddVertexOnEdge(&mesh, n);
/**----------------------------------------------------------------------------
        check neigbours list's consistency (neighbourship symetry)*/
        status=fe_chkvertex_consistency(&mesh,0);
        status=fe_displacevertex(mesh, node, x, y);
        if(status==0) {
//          movable[node]=false;
          mesh.vertices[node].code=0;
          }
        else {
          printf("not suitable\n");
          }
        status=fe_list(&mesh);
        status= fe_edgetable(&mesh,0,0);
        delete[] movable;
        goto start;
        }
      }
    }
  
  int RecycleCodes=0, StopOn_EdgeError=1, StopOn_PinchError=0;
  int SetLimits=1, verbose=0;
  status= fe_codetable(mesh, RecycleCodes, SetLimits, verbose);
  if(status!=0) {
    printf("unable to rebuild the limits table and codes of the original mesh\n");
    }
  
  status=fe_savemesh(output.c_str(),MESH_FILE_FORMAT_TRIGRID,mesh);
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_presplit(mesh_t & mesh, const char *polygonfile, double threshold, string rootname)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Purpose : align mesh nodes on polygon's segments before splitting

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
{
  int status;
  vector<plg_t> polygons;
  
  status=plg_load(polygonfile, PLG_FORMAT_UNKNOWN, polygons);
  if(status!=0) return(-1);
  
  status=fe_presplit(mesh, polygons, threshold, rootname);
  
  plg_destroy_entries(polygons);
  
  polygons.clear();
  
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  mesh_t *fe_split_Q(mesh_t & mesh,int *selected, int *targeted, int *frontier, int size, string rootname, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  selected: edge flag for interior/exterior splitting (1 is interior)
  targeted:
  
  frontier: edge flag, 1 if located along the interior/exterior splitting limit
            setup if non-zero address

            
  size    : number of edges needed to fulfill criterion. if size > 0, then option=0
            else option=1
  
  option 0: select triangle as interior if at least "size" edges are eligible
  
    when size is 3, only triangles strictly included in the selection are taken as
    interior. when size is 1, boarder triangles are taken as interior
  
  option 1: select triangle as exterior if at least "size" edges are not eligible
  
    when size is 3, only triangles strictly excluded from the selection are taken as
    exterior. when size is 1, boarder triangles are taken as exterior
  
    non-symmetric condition
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int   l,m,m1,m2,i,n1,n2,n,status;
  int   count=0;
  int   nndes;
  int   ninternal,nexternal;
  bool  *keep=new bool[mesh.ntriangles];
  int   *used=NULL;
  mesh_t *work=NULL;
  string filename;
  int option;
  
  mesh.set_elements();

  work=new mesh_t[2];

  printf ("number of nodes    (original mesh): %6d\n",mesh.nvtxs);
  printf ("number of elements (original mesh): %6d\n",mesh.nelements);

  ninternal=0;
  nexternal=0;


  if(size>0) option=0;
  else {
    option =1;
    size=-size;
    }

  for (m=0;m<mesh.nelements;m++) {
    keep[m]=false;
    count=0;
    switch (option) {
      case 0:
        for(i=0;i<3;i++) {
          n=mesh.elements[m].edges[i];
          if(selected[n]==1) {
            count++;
            }
          }
        keep[m]=(count>=size);
        break;
      case 1:
        for(i=0;i<3;i++) {
          n=mesh.elements[m].edges[i];
          if(selected[n]==0) {
            count++;
            }
          }
        keep[m]=not (count>=size);
        break;
      }
    switch(keep[m]) {
      case false:
        nexternal++;
        break;
      case true:
        ninternal++;
        break;
      }
    }

  work[0].elements=new element_t [nexternal];
  work[1].elements=new element_t [ninternal];

  work[0].nelements=nexternal;
  work[1].nelements=ninternal;

  work[0].type=SPHERICAL;
  work[1].type=SPHERICAL;
  
/*------------------------------------------------------------------------------
  screen elements*/
  ninternal=0;
  nexternal=0;
  for (m=0;m<mesh.nelements;m++) {
    switch(keep[m]) {
      case false:
        work[0].elements[nexternal].ancestor=m;
        nexternal++;
        break;
      case true:
        work[1].elements[ninternal].ancestor=m;
        ninternal++;
        break;
      }
    }

  used=new int[mesh.nvtxs];
  for(l=0;l<2;l++) {
    nndes=0;
    for (n=0;n<mesh.nvtxs;n++) used[n]=-1;
    for (m=0;m<work[l].nelements;m++) {
      for(i=0;i<3;i++) {
        n=mesh.elements[work[l].elements[m].ancestor].vertex[i];
        if(used[n]==-1) {
          used[n]=nndes;
          nndes++;
          }
        }
      }
    work[l].nvtxs = nndes;
    work[l].vertices=new vertex_t[work[l].nvtxs];
    for (n=0; n<work[l].nvtxs; n++) work[l].vertices[n].null_value();
    for (n=0;n<mesh.nvtxs;n++) {
       if(used[n]!=-1) {
         work[l].vertices[used[n]].lon=mesh.vertices[n].lon;
         work[l].vertices[used[n]].lat=mesh.vertices[n].lat;
         work[l].vertices[used[n]].h  =mesh.vertices[n].h;
         work[l].vertices[used[n]].code=0;
         work[l].vertices[used[n]].ancestor=n;
         }
       }
    for (m=0;m<work[l].nelements;m++) {
      for(i=0;i<3;i++) {
        n=mesh.elements[work[l].elements[m].ancestor].vertex[i];
        if(used[n]!=-1) {
          work[l].elements[m].vertex[i]=used[n];
          }
        else {
          printf("error...\n");
          }
        }
      }
    printf ("number of nodes    (splitted mesh %d): %6d\n",l,work[l].nvtxs);
    printf ("number of elements (splitted mesh %d): %6d\n",l,work[l].nelements);

    status=fe_e2n(&(work[l]));
    if(debug) status=fe_savemesh("tmp-00.nei",MESH_FILE_FORMAT_TRIGRID,work[l]);

    status=fe_geometry(&(work[l]));
    status=fe_edgetable(&(work[l]),0,0);
    status=fe_vertex_crosstables02(&(work[l]));

    int RecycleCodes=0, stopon_EdgeError=1, stopon_PinchError=0;
    status=fe_codetable1(&(work[l]), RecycleCodes, stopon_EdgeError, stopon_PinchError);
    if(status!=0) {
      printf ("code table failed for mesh %d\n",l);
      return(0);
      }
      
    if((l==1) && (frontier!=0)) {
/*------------------------------------------------------------------------------
      no information on edge's ancestor, use elements*/
      for (m1=0;m1<work[l].nelements;m1++) {
        m2=work[l].elements[m1].ancestor;
        for(i=0;i<3;i++) {
          n1=work[l].elements[m1].edges[i];
          n2=mesh.elements[m2].edges[i];
          if((mesh.edges[n2].code==MESH_INTERIOR_EDGE) && (work[l].edges[n1].code==1)) {
/*------------------------------------------------------------------------------
            internal/external limit, do not refine*/
            frontier[n1]=1;
            }
          else {
            frontier[n1]=0;
            }
          }
        }
      }
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  save interior and exterior meshes
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if(rootname=="") rootname="anonymous";
  
  filename=rootname+"-external.nei";
  status=fe_savemesh(filename.c_str(),MESH_FILE_FORMAT_TRIGRID,work[0]);
  
  filename=rootname+"-external.plg";
  status=fe_SaveLimits(filename.c_str(),work[0]);
  
  filename=rootname+"-internal.nei";
  status=fe_savemesh(filename.c_str(),MESH_FILE_FORMAT_TRIGRID,work[1]);
  
  filename=rootname+"-internal.plg";
  status=fe_SaveLimits(filename.c_str(),work[1]);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  create boundary polygon file for further mesh assembly
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  l=1;
  char *flag=new char[work[l].nedges];
  
  count=0;
  for(int n=0;n<work[l].nedges;n++) {
    if(work[l].edges[n].code==MESH_INTERIOR_EDGE) flag[n]='I';
    else flag[n]='T';
    }
    
  for (m1=0;m1<work[l].nelements;m1++) {
    m2=work[l].elements[m1].ancestor;
    for(i=0;i<3;i++) {
      n1=work[l].elements[m1].edges[i];
      n2=mesh.elements[m2].edges[i];
/*------------------------------------------------------------------------------
      02/02/2016: may have side-effects... */
      if((mesh.edges[n2].code==MESH_INTERIOR_EDGE) && (work[l].edges[n1].code>=1)) {
/*------------------------------------------------------------------------------
        internal/external limit, do not refine*/
        flag[n1]='M';
        work[l].edges[n1].code=MESH_FLAGGED_EDGE;
        count++;
        }
      }
    }
    
/*------------------------------------------------------------------------------
  export mesh limits as polygons, passing open/rigid flag */
  vector<plg_t> polygons, splitted, opened;
  status=fe_limits2poly(work[l], polygons, flag, false);
  
/*------------------------------------------------------------------------------
  split polygons to separate open/rigid segments */
  splitted=plg_split(polygons, false);
  
  for(int s=0;s<splitted.size();s++) {
    if(splitted[s].flag[0]=='M') opened.push_back(splitted[s]);
    }
  filename=rootname+"-opened.plg";
  status=plg_save(filename.c_str(), PLG_FORMAT_SCAN, opened);
  
  delete[] used;

  printf("split completed\n");

  return(work);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  mesh_t *fe_split_T(mesh_t & mesh,int *selected, int *targeted, int *frontier, int size, string rootname, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  selected: edge flag for interior/exterior splitting (1 is interior)
  targeted:
  
  frontier: edge flag, 1 if located along the interior/exterior splitting limit
            setup if non-zero address

            
  size    : number of edges needed to fulfill criterion. if size > 0, then option=0
            else option=1
  
  option 0: select triangle as interior if at least "size" edges are eligible
  
    when size is 3, only triangles strictly included in the selection are taken as
    interior. when size is 1, boarder triangles are taken as interior
  
  option 1: select triangle as exterior if at least "size" edges are not eligible
  
    when size is 3, only triangles strictly excluded from the selection are taken as
    exterior. when size is 1, boarder triangles are taken as exterior
  
    non-symmetric condition
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int   l,m,m1,m2,i,n1,n2,n,status;
  int   count=0;
  int   nndes;
  int   ninternal,nexternal;
  bool  *keep=new bool[mesh.ntriangles];
  int   *used=NULL;
  mesh_t *work=NULL;
  string filename;
  int option;

  work=new mesh_t[2];

  printf ("number of nodes    (original mesh): %6d\n",mesh.nvtxs);
  printf ("number of elements (original mesh): %6d\n",mesh.ntriangles);

  ninternal=0;
  nexternal=0;


  if(size>0) option=0;
  else {
    option =1;
    size=-size;
    }

  for (m=0;m<mesh.ntriangles;m++) {
    keep[m]=false;
    count=0;
    switch (option) {
      case 0:
        for(i=0;i<3;i++) {
          n=mesh.triangles[m].edges[i];
          if(selected[n]==1) {
            count++;
            }
          }
        keep[m]=(count>=size);
        break;
      case 1:
        for(i=0;i<3;i++) {
          n=mesh.triangles[m].edges[i];
          if(selected[n]==0) {
            count++;
            }
          }
        keep[m]=not (count>=size);
        break;
      }
    switch(keep[m]) {
      case false:
        nexternal++;
        break;
      case true:
        ninternal++;
        break;
      }
    }

  work[0].triangles=new triangle_t [nexternal];
  work[1].triangles=new triangle_t [ninternal];

  work[0].ntriangles=nexternal;
  work[1].ntriangles=ninternal;

  work[0].type=SPHERICAL;
  work[1].type=SPHERICAL;
  
/*------------------------------------------------------------------------------
  screen elements*/
  ninternal=0;
  nexternal=0;
  for (m=0;m<mesh.ntriangles;m++) {
    switch(keep[m]) {
      case false:
        work[0].triangles[nexternal].ancestor=m;
        nexternal++;
        break;
      case true:
        work[1].triangles[ninternal].ancestor=m;
        ninternal++;
        break;
      }
    }

  used=new int[mesh.nvtxs];
  for(l=0;l<2;l++) {
    nndes=0;
    for (n=0;n<mesh.nvtxs;n++) used[n]=-1;
    for (m=0;m<work[l].ntriangles;m++) {
      for(i=0;i<3;i++) {
        n=mesh.triangles[work[l].triangles[m].ancestor].vertex[i];
        if(used[n]==-1) {
          used[n]=nndes;
          nndes++;
          }
        }
      }
    work[l].nvtxs = nndes;
    work[l].vertices=new vertex_t[work[l].nvtxs];
    for (n=0; n<work[l].nvtxs; n++) work[l].vertices[n].null_value();
    for (n=0;n<mesh.nvtxs;n++) {
       if(used[n]!=-1) {
         work[l].vertices[used[n]].lon=mesh.vertices[n].lon;
         work[l].vertices[used[n]].lat=mesh.vertices[n].lat;
         work[l].vertices[used[n]].h  =mesh.vertices[n].h;
         work[l].vertices[used[n]].code=0;
         work[l].vertices[used[n]].ancestor=n;
         }
       }
    for (m=0;m<work[l].ntriangles;m++) {
      for(i=0;i<3;i++) {
        n=mesh.triangles[work[l].triangles[m].ancestor].vertex[i];
        if(used[n]!=-1) {
          work[l].triangles[m].vertex[i]=used[n];
          }
        else {
          printf("error...\n");
          }
        }
      }
    printf ("number of nodes    (splitted mesh %d): %6d\n",l,work[l].nvtxs);
    printf ("number of elements (splitted mesh %d): %6d\n",l,work[l].ntriangles);

    status=fe_e2n(&(work[l]));
    if(debug) status=fe_savemesh("tmp-00.nei",MESH_FILE_FORMAT_TRIGRID,work[l]);

    status=fe_geometry(&(work[l]));
    status=fe_edgetable(&(work[l]),0,0);
    status=fe_vertex_crosstables02(&(work[l]));

    int RecycleCodes=0, stopon_EdgeError=1, stopon_PinchError=0;
    status=fe_codetable1(&(work[l]), RecycleCodes, stopon_EdgeError, stopon_PinchError);
    if(status!=0) {
      printf ("code table failed for mesh %d\n",l);
      return(0);
      }
      
    if((l==1) && (frontier!=0)) {
/*------------------------------------------------------------------------------
      no information on edge's ancestor, use elements*/
      for (m1=0;m1<work[l].ntriangles;m1++) {
        m2=work[l].triangles[m1].ancestor;
        for(i=0;i<3;i++) {
          n1=work[l].triangles[m1].edges[i];
          n2=mesh.triangles[m2].edges[i];
          if((mesh.edges[n2].code==MESH_INTERIOR_EDGE) && (work[l].edges[n1].code==1)) {
/*------------------------------------------------------------------------------
            internal/external limit, do not refine*/
            frontier[n1]=1;
            }
          else {
            frontier[n1]=0;
            }
          }
        }
      }
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  save interior and exterior meshes
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if(rootname=="") rootname="anonymous";
  
  filename=rootname+"-external.nei";
  status=fe_savemesh(filename.c_str(),MESH_FILE_FORMAT_TRIGRID,work[0]);
  
  filename=rootname+"-external.plg";
  status=fe_SaveLimits(filename.c_str(),work[0]);
  
  filename=rootname+"-internal.nei";
  status=fe_savemesh(filename.c_str(),MESH_FILE_FORMAT_TRIGRID,work[1]);
  
  filename=rootname+"-internal.plg";
  status=fe_SaveLimits(filename.c_str(),work[1]);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  create boundary polygon file for further mesh assembly
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  l=1;
  char *flag=new char[work[l].nedges];
  
  count=0;
  for(int n=0;n<work[l].nedges;n++) {
    if(work[l].edges[n].code==MESH_INTERIOR_EDGE) flag[n]='I';
    else flag[n]='T';
    }
    
  for (m1=0;m1<work[l].ntriangles;m1++) {
    m2=work[l].triangles[m1].ancestor;
    for(i=0;i<3;i++) {
      n1=work[l].triangles[m1].edges[i];
      n2=mesh.triangles[m2].edges[i];
/*------------------------------------------------------------------------------
      02/02/2016: may have side-effects... */
      if((mesh.edges[n2].code==MESH_INTERIOR_EDGE) && (work[l].edges[n1].code>=1)) {
/*------------------------------------------------------------------------------
        internal/external limit, do not refine*/
        flag[n1]='M';
        work[l].edges[n1].code=MESH_FLAGGED_EDGE;
        count++;
        }
      }
    }
    
/*------------------------------------------------------------------------------
  export mesh limits as polygons, passing open/rigid flag */
  vector<plg_t> polygons, splitted, opened;
  status=fe_limits2poly(work[l], polygons, flag, false);
  
/*------------------------------------------------------------------------------
  split polygons to separate open/rigid segments */
  splitted=plg_split(polygons, false);
  
  for(int s=0;s<splitted.size();s++) {
    if(splitted[s].flag[0]=='M') opened.push_back(splitted[s]);
    }
  filename=rootname+"-opened.plg";
  status=plg_save(filename.c_str(), PLG_FORMAT_SCAN, opened);
  
  delete[] used;

  printf("split completed\n");

  return(work);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  mesh_t *fe_split(mesh_t & mesh,int *selected, int *targeted, int *frontier, int size, string rootname, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  wrapper
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int status;
  mesh_t *work;
  
  if(mesh.ntriangles!=0) {
    work=fe_split_T(mesh, selected, targeted, frontier, size, rootname, debug);
    }
  else if(mesh.nquadrangles!=0) {
    work=fe_split_Q(mesh, selected, targeted, frontier, size, rootname, debug);
    }

  return(work);
}
    
    
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  mesh_t *fe_split(mesh_t & mesh, vector<plg_t> & polygons, bool exact, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  wrapper
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int count,k,l,m,n,pos,status;
  int size, nselected;
  int *selected, *targeted;
  mesh_t *splitted;
  
  if(exact) {
    status=fe_presplit(mesh, polygons, 0.5, "");
    if(status!=0) return(0);
    size=+3;
    }
  else {
    size=-3;
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  select edges
 
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  selected=new int[mesh.nedges];
  targeted=0;
  
  for (n=0;n<mesh.nedges;n++) {
    selected[n]=1;
    }
    
  int option =1; 
  nselected=fe_selectedges_01(mesh, polygons, selected, option, false);
  
  plg_destroy_entries(polygons);
  
  polygons.clear();
  
  if(nselected==0) return(0);

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  extract submesh from selection to allow decent cartesian reshape
 
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  printf("#################################################################\n");
  printf("split mesh\n");
  splitted=fe_split(mesh, selected, targeted, 0, size, (string) "anonymous", debug);
  if(splitted==NULL) {
    printf("unable to split the original mesh\n");
    return(0);
    }
    
  deletep(&selected);
  deletep(&targeted);

  printf("#################################################################\n");
  printf("check cut connexity\n");
  status=fe_connex(splitted[0]);
  if(status==-1) {
    printf("splitted mesh not #0 connex\n");
    status=-1;
    }
  status=fe_connex(splitted[1]);
  if(status==-1) {
    printf("splitted mesh not #1 connex\n");
    status=-1;
    }
  
  return(splitted);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_MeshCut(mesh_t & mesh, vector<plg_t> & polygons, bool exact, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
{
  int count,k,l,m,n,pos,status;
  int size, nselected;
  int *selected=0, *targeted=0;
  mesh_t *splitted;
  
  
  if(exact) {
    status=fe_presplit(mesh, polygons, 0.5, "");
    if(status!=0) return(status);
    size=+3;
    }
  else {
    size=-3;
    }

  selected=new int[mesh.nedges];
  targeted=new int[mesh.nedges];
  
  for (n=0;n<mesh.nedges;n++) {
    selected[n]=1;
    }
    
/*------------------------------------------------------------------------------
  select edges to be refined */
  int option =1; 
  nselected=fe_selectedges_01(mesh, polygons, selected, option, false);
  
  plg_destroy_entries(polygons);
  
  polygons.clear();
  
  if(nselected==0) return(0);

/*------------------------------------------------------------------------------
  extract submesh from selection to allow decent cartesian reshape */
  printf("#################################################################\n");
  printf("cut mesh\n");
  splitted=fe_split(mesh, selected, targeted, 0, size, (string) "mesh-reshape", debug);
  if(splitted==NULL) {
    printf("unable to split the original mesh\n");
    return(-1);
    }

  deletep(&selected);
  deletep(&targeted);

  printf("#################################################################\n");
  printf("check cut connexity\n");
  status=fe_connex(splitted[0]);
  if(status==-1) {
    printf("cut mesh not connex\n");
    status=-1;
    }
  
  splitted[1].destroy();
  mesh.destroy();
  
  mesh=splitted[0];

  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_MeshCut(mesh_t & mesh, const char *polygonfile, bool exact, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------
    
------------------------------------------------------------------------------*/
{
  int count,k,l,m,n,pos,status;
  int size, nselected;
  int *selected, *targeted;
  mesh_t *splitted;
  vector<plg_t> polygons;

  status=plg_load(polygonfile, PLG_FORMAT_UNKNOWN, polygons);
  if(status!=0) return(-1);
  
  status=fe_MeshCut(mesh, polygons, exact, debug);
  

  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  mesh_t *fe_node2mesh(mesh_t mesh, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int   n,status;
  int   nelts,nndes;
  double lon,lat;
  double *x=NULL,*y=NULL,xx,yy;
  projPJ projection;
  mesh_t *finished=NULL;

  nndes=mesh.nvtxs;
  nelts=mesh.ntriangles;

/*------------------------------------------------------------------------------
  work in cartesian coordinates */
  lon=mesh.vertices[0].lon;
  lat=mesh.vertices[0].lat;
  for (n=1;n<nndes;n++) {
    if(mesh.type==0) mesh.vertices[n].lon = geo_recale(mesh.vertices[n].lon, lon, (double) 180.0);
    }

  lon=0;
  lat=0;
  for (n=0;n<nndes;n++) {
    lon+=mesh.vertices[n].lon;
    lat+=mesh.vertices[n].lat;
    }
  lon/=nndes;
  lat/=nndes;

//  projection=geo_mercator_init(lon,lat,radius);
  projection=assign_StereoOblique(lat,lon);

  x=new double[nndes];
  y=new double[nndes];

  for (n=0;n<nndes;n++) {
//    status=geo_mercator_directe(projection,mesh.vertices[n].lon,mesh.vertices[n].lat,&x[n],&y[n]);
    geo_to_projection(projection,mesh.vertices[n].lat, mesh.vertices[n].lon, &x[n], &y[n]);
    mesh.vertices[n].lon=x[n];
    mesh.vertices[n].lat=y[n];
    }

// #ifdef EXTERN_TRIANGLE
//  status=fe_savenodes("refined.poly", NODE_FILE_FORMAT_TRIANGLE, mesh);
//  if(status!=0) return(NULL);
//
//  status=fe_savenodes("refined.nod", NODE_FILE_FORMAT_TRIGRID, mesh);
//  if(status!=0) return(NULL);
//  status=system("node2poly.exe refined.nod refined.poly");
//  if(status!=0) return(NULL);
//
//  status=system("triangle.exe -p refined.poly");
//  if(status!=0) return(NULL);
// #endif
//
// /*------------------------------------------------------------------------------

//   transfert node informations to triangle input structure */
//   status=fe_node2triangle( mesh, &in);
//
// /*------------------------------------------------------------------------------

//   initialize triangle output structure */
//   out.pointlist = (REAL *) NULL;            /* Not needed if -N switch used. */
//   out.pointattributelist = (REAL *) NULL;   /* Not needed if -N switch used or number of point attributes is zero: */
//   out.pointmarkerlist = (int *) NULL;       /* Not needed if -N or -B switch used. */
//   out.trianglelist = (int *) NULL;          /* Not needed if -E switch used. */
//   out.triangleattributelist = (REAL *) NULL;/* Not needed if -E switch used or number of triangle attributes is zero: */
//   out.neighborlist = (int *) NULL;          /* Needed only if -n switch used. */
//   out.segmentlist = (int *) NULL;           /* Needed only if segments are output (-p or -c) and -P not used: */
//   out.segmentmarkerlist = (int *) NULL;     /* Needed only if segments are output (-p or -c) and -P and -B not used: */
//   out.edgelist = (int *) NULL;              /* Needed only if -e switch used. */
//   out.edgemarkerlist = (int *) NULL;        /* Needed if -e used and -B not used. */
//
//   triangulate("pzN", &in, &out, (triangulateio *) NULL);
//
// /*------------------------------------------------------------------------------

//   output structure inherit from input structure's pointlist */
//   out.pointlist = in.pointlist;
//
//   finished=new mesh_t;
//
//   status=fe_init(finished);
//
// #ifdef EXTERN_TRIANGLE
//   status=fe_readmesh_TGL ("refined.1.ele","refined.1.node",finished);
// #endif
//   status=fe_triangle2mesh(out,finished);

  finished=fe_triangulate(mesh,0 ,0, debug);
  if(mesh.nvtxs!=nndes) {
    status=fe_savemesh("refined.nei",MESH_FILE_FORMAT_TRIGRID, *finished);
    }
/*------------------------------------------------------------------------------
  back to geocentric coordinates */
  for (n=0;n<nndes;n++) {
    xx=finished->vertices[n].lon;
    yy=finished->vertices[n].lat;
//    status=geo_mercator_inverse(projection,&(finished->vertices[n].lon),&(finished->vertices[n].lat),x[n],y[n]);
    projection_to_geo(projection,&(finished->vertices[n].lat),&(finished->vertices[n].lon), xx, yy);
    }

  lon=finished->vertices[0].lon;
  lat=finished->vertices[0].lat;
  for (n=1;n<nndes;n++) {
    if(mesh.type==0) finished->vertices[n].lon = geo_recale(finished->vertices[n].lon, lon, (double) 180.0);
    }

  finished->type=SPHERICAL;
  status=fe_edgetable(finished,0,0);

  status=fe_codetable1(finished,0,1,0);
  status=fe_savemesh("refined.nei",MESH_FILE_FORMAT_TRIGRID, *finished);

  for (n=0;n<nndes;n++) {
    finished->vertices[n].lon=x[n];
    finished->vertices[n].lat=y[n];
    }
  finished->type=CARTESIAN;
  status=fe_reshapeall(*finished,3);
  
  for (n=0;n<nndes;n++) {
    xx=finished->vertices[n].lon;
    yy=finished->vertices[n].lat;
    projection_to_geo(projection,&(finished->vertices[n].lat),&(finished->vertices[n].lon), xx, yy);
    }
  lon=finished->vertices[0].lon;
  lat=finished->vertices[0].lat;
  for (n=1;n<nndes;n++) {
    if(mesh.type==0) finished->vertices[n].lon = geo_recale(finished->vertices[n].lon, lon, (double) 180.0);
    }

  finished->type=SPHERICAL;
  status=fe_savemesh("refined-reshaped.nei",MESH_FILE_FORMAT_TRIGRID, *finished);

//   for (n=0;n<nndes;n++) {
//     if(finished->vertices[n].nngh>8) status=fe_cleave(finished,n);
//     }
//   status=fe_list(finished);
//   status=fe_edgetable(finished,0);
//   status=fe_codetable1(finished,0);
//   status=fe_reshapeall(*finished,3);
//   status=fe_savemesh("refined-cleaved.nei",MESH_FILE_FORMAT_TRIGRID, *finished);

  delete[] x;
  delete[] y;

  return(finished);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_Triangles2Quadrangle(mesh_t & mesh, int e, quadrangle_t *q)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------

  form a quadrangle from 2 adjacent triangles */
{
  int   k,k0,l,l0,m,m1,m2,n,nodes[2],status=0;
  
/*------------------------------------------------------------------------------
  boundary edge, not suitable */
  if(mesh.edges[e].nshared==0) return(-1);
  
  m1=mesh.edges[e].shared[0];
  m2=mesh.edges[e].shared[1];

  for(k=0;k<2;k++) nodes[k]=-1;
  
  nodes[0]=mesh.edges[e].extremity[0];
  nodes[1]=mesh.edges[e].extremity[1];
  
/*------------------------------------------------------------------------------
  seek extra-edge node (3rd) in element 1 */
  for(k=0;k<3;k++) {
    n=mesh.triangles[m1].vertex[k];
    if((n!=nodes[0]) && (n!=nodes[1])) {
      k0=k;
      break;
      }
    }
    
/*------------------------------------------------------------------------------
  seek extra-edge node (3rd) in element 2 */
  for(k=0;k<3;k++) {
    n=mesh.triangles[m2].vertex[k];
    if((n!=nodes[0]) && (n!=nodes[1])) {
      l0=k;
      break;
      }
    }

/*------------------------------------------------------------------------------
  build quadrangle vertex list */
  k=k0;
  q->vertex[0]=mesh.triangles[m1].vertex[k];
  k=(k+1)%3;
  q->vertex[1]=mesh.triangles[m1].vertex[k];
  q->vertex[2]=mesh.triangles[m2].vertex[l0];
  k=(k+1)%3;
  q->vertex[3]=mesh.triangles[m1].vertex[k];

/*------------------------------------------------------------------------------
  build quadrangle edge list (all but internal T-edge) */

/*------------------------------------------------------------------------------
  edge from k0 to k0+1, i.e opposite to k0+2 */
  k=(k0+2)%3;
  q->edges[0]=mesh.triangles[m1].edges[k];
/*------------------------------------------------------------------------------
  edge from l0+2 to l0, i.e opposite to l0+1 */
  k=(l0+1)%3;
  q->edges[1]=mesh.triangles[m2].edges[k];
/*------------------------------------------------------------------------------
  edge from l0 to l0+1, i.e opposite to l0+2 */
  k=(l0+2)%3;
  q->edges[2]=mesh.triangles[m2].edges[k];
/*------------------------------------------------------------------------------
  edge from k0+2 to k0, i.e opposite to k0+1 */
  k=(k0+1)%3;
  q->edges[3]=mesh.triangles[m1].edges[k];
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_Quadrangles2Triangles(mesh_t & mesh, int e, quadrangle_t q, int *cut)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------

  seek best cut into 2 triangles */
{
  int   k,k0,l,l0,m,m1,m2,n,nodes[2],status;
  double f[4];
  triangle_t *t[2];
  mesh_t work;
  
  work.vertices=mesh.vertices;
  work.nvtxs=2;
  
  work.triangles=new triangle_t[2];
  work.ntriangles=2;
  
  work.type=SPHERICAL;
  
  t[0]=&(work.triangles[0]);
  t[1]=&(work.triangles[1]);
  
/*------------------------------------------------------------------------------
  compute aspect ratio for origninal triangles */
  t[0]->vertex[0]=q.vertex[0];
  t[0]->vertex[1]=q.vertex[1];
  t[0]->vertex[2]=q.vertex[3];
  
  t[1]->vertex[0]=q.vertex[2];
  t[1]->vertex[1]=q.vertex[3];
  t[1]->vertex[2]=q.vertex[1];

  status=fe_initaffine(&work,0, false);
  status=fe_initaffine(&work,1, false);
  
  f[0]=fe_AspectRatio(work,0);
  f[1]=fe_AspectRatio(work,1);
  
/*------------------------------------------------------------------------------
  compute aspect ratio for swapped triangles */
  t[0]->vertex[0]=q.vertex[3];
  t[0]->vertex[1]=q.vertex[0];
  t[0]->vertex[2]=q.vertex[2];
  
  t[1]->vertex[0]=q.vertex[1];
  t[1]->vertex[1]=q.vertex[2];
  t[1]->vertex[2]=q.vertex[0];
  
  status=fe_initaffine(&work,0, false);
  if(status!=0) {
    *cut=0;
    goto terminate;
    }
  status=fe_initaffine(&work,1, false);
  if(status!=0) {
    *cut=0;
    goto terminate;
    }
  
  f[2]=fe_AspectRatio(work,0);
  f[3]=fe_AspectRatio(work,1);
  
/*------------------------------------------------------------------------------
  then compare... */
  if(f[0]+f[1]>f[2]+f[3]) {
    *cut=0;
    }
  else {
    *cut=1;
    }

  terminate:
  
  work.vertices=0;
  work.nvtxs=0;
  
  work.destroy();
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_SubDomainE(mesh_t & mesh, vector<int> & elements)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int m,mm,n,status=0;
  int nelements=elements.size();
  
  for(int k=0;k<nelements;k++) {
    m=elements[k];
    for(int i=0;i<3;i++) {
      n=mesh.triangles[m].vertex[i];
      for(int l=0;l<mesh.vertices[n].nelmts;l++) {
        mm=mesh.vertices[n].elmts[l];
        if(vpos(mm,elements)==-1) elements.push_back(mm);
        }
      }
    }
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_SubDomainE(mesh_t & mesh, vector<int> & elements, int extent)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  for(int k=0;k<extent;k++) status=fe_SubDomainE(mesh, elements);
  
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_SubDomainN(mesh_t & mesh, vector<int> & vertices, int extent)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  vector<int> elements;
  
  for(int k=0;k<vertices.size();k++) {
    int n=vertices[k];
    for(int l=0;l<mesh.vertices[n].nelmts;l++) {
      int mm=mesh.vertices[n].elmts[l];
      if(vpos(mm,elements)==-1) elements.push_back(mm);
      }
    }

  for(int k=0;k<extent;k++) status=fe_SubDomainE(mesh, elements, extent);
  
  for(int k=0;k<elements.size();k++) {
    int m=elements[k];
    for(int i=0;i<3;i++) {
      int n=mesh.triangles[m].vertex[i];
      if(vpos(n,vertices)==-1) vertices.push_back(n);
      }
    }
    
  elements.clear();
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_exchange(mesh_t & mesh, int e, bool reshape, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------

  3-----2     3-----2
  |\    |     |    /|
  | \   |     |   / |
  |  \  |     |  /  |
  |   \ |     | /   |
  |    \|     |/    |
  0-----1     0-----1
  
-------------------------------------------------------------------------------*/
{
  int   k,l,m,m1,m2,i,j,n,n1,n2,status;
  int cut;
  double lat,lon,t1,t2,p1,p2,d;
  projPJ projection;
  quadrangle_t q;
  vector<int> elements, vertices;
  int extent=5,npass=3;
  
  mesh_t work;
  
/*------------------------------------------------------------------------------

  boundary edge, nothing to try */
  if(mesh.edges[e].nshared==1) return(-1);
  
  m1=mesh.edges[e].shared[0];
  m2=mesh.edges[e].shared[1];
  
/*------------------------------------------------------------------------------

  define Oblique Stereographic projection */
//  projection=assign_StereoOblique(lat,lon);

/*------------------------------------------------------------------------------

  form a quadrangle and seek best cut into 2 triangles */
  status=fe_Triangles2Quadrangle(mesh,e,&q);
  if(status!=0) {
    return(-1);
    }
  status=fe_Quadrangles2Triangles(mesh,e,q,&cut);
  if(status!=0) {
    return(-1);
    }
  
//   for(k=0;k<4;k++) n[k]=-1;
  
  switch(cut) {
    case 0:
/*------------------------------------------------------------------------------

      nothing to do, best combination already there*/
      status=-1;
      break;
      
    case 1:
/*------------------------------------------------------------------------------

      exchange edge line*/
      n1=q.vertex[1];
      n2=q.vertex[3];
      status=fe_disconnectvertices(mesh, n1, n2);
      if(status!=0) {
        return(-1);
        }
      n1=q.vertex[0];
      n2=q.vertex[2];
      status=fe_connectvertices(mesh, n1, n2);
      if(status!=0) {
        return(-1);
        }
        
      mesh.edges[e].extremity[0]=n1;
      mesh.edges[e].extremity[1]=n2;
      
      m1=mesh.edges[e].shared[0];
      m2=mesh.edges[e].shared[1];
      
      mesh.triangles[m1].vertex[0]=q.vertex[0];
      mesh.triangles[m1].vertex[1]=q.vertex[1];
      mesh.triangles[m1].vertex[2]=q.vertex[2];
      
      status=fe_initaffine(&mesh,m1);
      
      for(i=0;i<3;i++) {
        n=mesh.triangles[m1].vertex[i];
        status=fe_vertex_Etable(mesh,n);
        }
        
      mesh.triangles[m2].vertex[0]=q.vertex[0];
      mesh.triangles[m2].vertex[1]=q.vertex[2];
      mesh.triangles[m2].vertex[2]=q.vertex[3];
      
      status=fe_initaffine(&mesh,m2);
      
      for(i=0;i<3;i++) {
        n=mesh.triangles[m2].vertex[i];
        status=fe_vertex_Etable(mesh,n);
        }

      n=q.edges[1];
      mesh.triangles[m1].edges[0]=n;      
//       mesh.edges[n].L=mesh.triangles[m1].l[0];
      if(mesh.edges[n].L!=mesh.triangles[m1].l[0]) {
        return(-1);
        }
        

      n=e;
      mesh.triangles[m1].edges[1]=n;
      mesh.edges[n].L=mesh.triangles[m1].l[1];
      
      n=q.edges[0];
      mesh.triangles[m1].edges[2]=n;
//       mesh.edges[n].L=mesh.triangles[m1].l[2];
      if(mesh.edges[n].L!=mesh.triangles[m1].l[2]) {
        return(-1);
        }

      
      n=q.edges[2];
      mesh.triangles[m2].edges[0]=n;
//       mesh.edges[n].L=mesh.triangles[m2].l[0];
      if(mesh.edges[n].L!=mesh.triangles[m2].l[0]) {
        return(-1);
        }
      
      
      n=q.edges[3];
      mesh.triangles[m2].edges[1]=n;
//       mesh.edges[n].L=mesh.triangles[m2].l[1];
      if(mesh.edges[n].L!=mesh.triangles[m2].l[1]) {
        return(-1);
        }
      
      n=e;
      mesh.triangles[m2].edges[2]=n;
//       mesh.edges[n].L=mesh.triangles[m2].l[2];
      if(mesh.edges[n].L!=mesh.triangles[m2].l[2]) {
        return(-1);
        }
      
      if(reshape) {
        for(i=0;i<4;i++) {
          vertices.push_back(q.vertex[i]);
          }
        status=fe_SubDomainN(mesh, vertices, extent);
        status=fe_reshape(mesh, npass, vertices);
        vertices.clear();
        }
      status=0;
      break;
      
    default:
      break;
    }

  return(status);
}
