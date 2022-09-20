
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


#include <stdio.h>
#include <string.h>
 
#include "tools-structures.h"
#include "constants.h"

#include "geo.h"
#include "fe.h"
#include "map.h"
#include "polygones.h"
#include "grd.h"
#include "parallel.h"

#include "zapper.h"     /*  rutin.h contains common utility routines  */


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  mesh_t fe_merge(const mesh_t work[2], double dmax, const char *filename, bool limited, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  merge 2 meshes

\param work array[2] of meshes to merge: for better performance, put the small one as 0th element
\param dmax distance in km
\param *filename Unless 0 (default), path to .nei file where merged mesh will be saved
\param limited whether to limit to the first boundary in first mesh. Default: false
\param debug

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int   k,l,m,i,n1,n2,n,status;
  int   lmin, nndes, nmerged=0, count=0;
  int   *used=NULL;
  double Lmin[2];
   
  if(limited) printf("merge limited to first boundary in first mesh\n");

  mesh_t *merged=new mesh_t;

  merged->type=SPHERICAL;

  merged->ntriangles=work[0].ntriangles+work[1].ntriangles;

  merged->triangles=new triangle_t[merged->ntriangles];
  
  for(l=0;l<2;l++) {
    for(n=0;n<work[l].nvtxs;n++) {
      work[l].vertices[n].lon*=d2r;
      work[l].vertices[n].lat*=d2r;
      }
    }

  for(l=0;l<2;l++) {
    Lmin[l]=1.e+10;
    for(n=0;n<work[l].nedges;n++) {
      Lmin[l]=min(Lmin[l],work[l].edges[n].L);
      if(work[l].edges[n].L<dmax) count++;
      }
    }

  printf("minimum edge length : %lf %lf\n",Lmin[0],Lmin[1]);
  
  if(1000.*dmax>Lmin[0] or 1000.*dmax>Lmin[1]) {
    printf("\n\nWARNING : dmax=%lf larger than %lf or %lf\n\n\n",dmax,Lmin[0],Lmin[1]);
    }
  
  used=new int[work[1].nvtxs];

  nndes=work[0].nvtxs;
  for(n2=0;n2<work[1].nvtxs;n2++) {
    double d,dmin,t1,p1,t2,p2;
    int merging;
    if(work[1].vertices[n2].code==0) {
//    if(work[1].vertices[n2].code!=1) {  ///HERE !!!
      used[n2]=nndes;
      nndes++;
      continue;
      }
    merging=0;
    t2=work[1].vertices[n2].lon;
    p2=work[1].vertices[n2].lat;
    dmin=dmax;
    for(l=0;l<work[0].nlimits;l++) {
      if(limited && l>0) {
        break;
        }
      for(k=0;k<work[0].limits[l].nvertex;k++) {
        n1=work[0].limits[l].vertex[k];
        t1=work[0].vertices[n1].lon;
        p1=work[0].vertices[n1].lat;
        d=geo_haversinR_km(t1,p1,t2,p2);
        if(d<dmin) {
          used[n2]=n1;
          lmin=l;
          dmin=d;
          }
        }
      }
//done:
    merging=(dmin<dmax);
    if(merging) {
//      used[n2]=n1;
      if(debug) printf("codes %d %d\n", work[1].vertices[n2].code, lmin+1);
      nmerged++;
      }
    else {
      used[n2]=nndes;
      nndes++;
      }
    }
  
  status=work[0].nvtxs+work[1].nvtxs-nmerged-nndes;
  printf("merging budget: %d merged nodes of %d in %d totals %d, remains %d (should be 0)\n",nmerged,work[0].nvtxs,work[1].nvtxs,nndes,status);
  if(status!=0) wexit(status);
  
  for(l=0;l<2;l++) {
    for(n=0;n<work[l].nvtxs;n++) {
      work[l].vertices[n].lon/=d2r;
      work[l].vertices[n].lat/=d2r;
      }
    }

  merged->nvtxs=nndes;
  merged->vertices=new vertex_t[merged->nvtxs];
  for (n=0; n<merged->nvtxs; n++) merged->vertices[n].null_value();

  for(n=0;n<work[0].nvtxs;n++) {
    merged->vertices[n].lon =work[0].vertices[n].lon;
    merged->vertices[n].lat =work[0].vertices[n].lat;
    merged->vertices[n].h   =work[0].vertices[n].h;
    merged->vertices[n].code=0;
    }
  for(n=0;n<work[1].nvtxs;n++) {
    merged->vertices[used[n]].lon =work[1].vertices[n].lon;
    merged->vertices[used[n]].lat =work[1].vertices[n].lat;
    merged->vertices[used[n]].h   =work[1].vertices[n].h;
    merged->vertices[used[n]].code=0;
    }
  for(m=0;m<work[0].ntriangles;m++) {
    for(i=0;i<3;i++) {
      merged->triangles[m].vertex[i]=work[0].triangles[m].vertex[i];
      }
    }
  for(m=0;m<work[1].ntriangles;m++) {
    for(i=0;i<3;i++) {
      merged->triangles[m+work[0].ntriangles].vertex[i]=used[work[1].triangles[m].vertex[i]];
      }
    }

  status=fe_e2n(merged);
  
  status=fe_geometry(merged);

//   debug=true;
  if (debug)status=fe_savemesh("merged.nei",MESH_FILE_FORMAT_TRIGRID,*merged);

  status=fe_edgetable(merged,0,0);

/* *---------------------------------------------------------------------------
  to be de-activated ? */
  status=fe_vertex_crosstables02(merged);
  
  int RecycleCodes=0, stopon_EdgeError=1, stopon_PinchError=0;
  status=fe_codetable1(merged, RecycleCodes, stopon_EdgeError, stopon_PinchError);
  if(filename!=0)
    status=fe_savemesh(filename,MESH_FILE_FORMAT_TRIGRID,*merged);
  
  delete[] used;

  printf("merge completed \n");
  return(*merged);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  mesh_t fe_merge(mesh_t & small, mesh_t & big, double dmax)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  mesh_t work[2];
  work[0] = small;
  work[1] = big;
  
  return fe_merge(work, dmax, (const char *) 0, false);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int *fe_read_epartition(const char *filename, const mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int  m,nitems;
  int  *partition=NULL;
  FILE *in=NULL;

  in=fopen(filename, "r");
  if(in==0) {
    TRAP_ERR_EXIT(-1,"cannot open partition file : %s \n",filename);
    }

  partition=new int[mesh.ntriangles];

  for (m=0;m<mesh.ntriangles;m++) {
    nitems=fscanf(in,"%d",&partition[m]);
    }
  fclose(in);

  return(partition);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int *fe_read_npartition(const char *filename, mesh_t mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int  m,nitems;
  int  *partition=NULL;
  FILE *in=NULL;

  in=fopen(filename, "r");
  if(in==0) {
      TRAP_ERR_EXIT(-1,"Probleme a l ouverture du fichier : %s \n",filename);
    }


  partition=new int[mesh.nvtxs];

  for (m=0;m<mesh.nvtxs;m++) {
    nitems=fscanf(in,"%d",&partition[m]);
    }
  fclose(in);

  return(partition);

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  mesh_t *fe_partition_01(mesh_t mesh,int *partition,int npartition)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Partitioned meshes are directly taken from metis element list.

  Not suitable for parallel computing (no meshes overlapping).

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
{
  int   m,i,n,p,status;
  int   nndes=0;
  int   ninternal,nexternal,keep;
  int   *used=NULL;
  mesh_t *work=NULL;
  char output[1024];

  printf ("number of nodes    (original mesh): %6d\n",mesh.nvtxs);
  printf ("number of elements (original mesh): %6d\n",mesh.ntriangles);

  work=new mesh_t[npartition];

  for(p=0;p<npartition;p++) {
    ninternal=0;
    nexternal=0;

    for (m=0;m<mesh.ntriangles;m++) {
      keep=(partition[m]==p);
      switch(keep) {
        case 0:
          nexternal++;
          break;
        case 1:
          ninternal++;
          break;
        }
      }

    work[p].triangles=new triangle_t [ninternal];

    work[p].ntriangles=ninternal;
    work[p].type=SPHERICAL;

    ninternal=0;
    nexternal=0;

/* *-----------------------------------------------------------------------------
    screen elements*/
    for (m=0;m<mesh.ntriangles;m++) {
      keep=(partition[m]==p);
      switch(keep) {
        case 0:
          nexternal++;
          break;
        case 1:
          work[p].triangles[ninternal].ancestor=m;
          ninternal++;
          break;
        }
      }

/* *-----------------------------------------------------------------------------
    count vertices*/
    used=new int[mesh.nvtxs];
    nndes=0;
    for (n=0;n<mesh.nvtxs;n++) used[n]=-1;
    for (m=0;m<work[p].ntriangles;m++) {
      for(i=0;i<3;i++) {
        n=mesh.triangles[work[p].triangles[m].ancestor].vertex[i];
        if(used[n]==-1) {
          used[n]=nndes;
          nndes++;
          }
        }
      }
/* *-----------------------------------------------------------------------------
    build vertices and elements table*/
    work[p].nvtxs = nndes;
    work[p].vertices=new vertex_t[work[p].nvtxs];
    for (n=0; n<work[p].nvtxs; n++) work[p].vertices[n].null_value();
    for (n=0;n<mesh.nvtxs;n++) {
      if(used[n]!=-1) {
        work[p].vertices[used[n]].lon=mesh.vertices[n].lon;
        work[p].vertices[used[n]].lat=mesh.vertices[n].lat;
        work[p].vertices[used[n]].code=0;
        work[p].vertices[used[n]].ancestor=n;
        }
      }
    for (m=0;m<work[p].ntriangles;m++) {
      for(i=0;i<3;i++) {
        n=mesh.triangles[work[p].triangles[m].ancestor].vertex[i];
        if(used[n]!=-1) {
          work[p].triangles[m].vertex[i]=used[n];
          }
        else {
          printf("error...\n");
          }
        }
      }
    printf ("number of nodes    (splitted mesh %d): %6d\n",p,work[p].nvtxs);
    printf ("number of elements (splitted mesh %d): %6d\n",p,work[p].ntriangles);

    status=fe_e2n(&(work[p]));
    status=fe_edgetable(&work[p],0,0);
//     sprintf(output,"%s-raw.%2.2d.nei","internal",p);
//     status=fe_savemesh(output,MESH_FILE_FORMAT_TRIGRID,work[p]);
    status=fe_codetable1(&work[p],0,1,0);
    sprintf(output,"%s-raw.%2.2d.nei","internal",p);
    status=fe_savemesh(output,MESH_FILE_FORMAT_TRIGRID,work[p]);
    delete[] used;
    }

  printf("split completed\n");

  return(work);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  mesh_t *fe_partition_02(mesh_t mesh,int *partition,int npartition)
//unused
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Partitioned meshes are taken from metis element list plus the
  neighbouring elements.

  Partially suitable for parallel computing, but meshes overlap to large.

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
{
  int   k,m,mm,i,n,nn,p,status;
  int   nndes=0;
  int   ninternal,nexternal,keep;
  int   *used=NULL,*flag=NULL,*elist=NULL;
  mesh_t *work=NULL;
  char output[1024];

  printf ("number of nodes    (original mesh): %6d\n",mesh.nvtxs);
  printf ("number of elements (original mesh): %6d\n",mesh.ntriangles);

  used  =new int[mesh.nvtxs];
  flag  =new int[mesh.nvtxs];
  elist =new int[mesh.ntriangles];

  work=new mesh_t[npartition];

  for(p=0;p<npartition;p++) {

    for (n=0;n<mesh.nvtxs;n++) flag[n]=-1;
    for (m=0;m<mesh.ntriangles;m++) elist[m]=-1;

/* *-----------------------------------------------------------------------------
    screen vertices*/
    for (m=0;m<mesh.ntriangles;m++) {
      if(partition[m]==p) {
        elist[m]=0;
        for(i=0;i<3;i++) {
          n=mesh.triangles[m].vertex[i];
          flag[n]=0;
          }
        }
      }
/* *-----------------------------------------------------------------------------
    add vertices and elements neighbours*/
    for (n=0;n<mesh.nvtxs;n++) {
      if(flag[n]==0) {
/* *-----------------------------------------------------------------------------
        */
        for(k=0;k<mesh.vertices[n].nngh;k++) {
          nn=mesh.vertices[n].ngh[k];
          if(flag[nn]==-1) flag[nn]=1;
          }
        for(k=0;k<mesh.vertices[n].nelmts;k++) {
          mm=mesh.vertices[n].elmts[k];
          if(elist[mm]==-1) elist[mm]=1;
          }
        }
      }

    ninternal=0;
    nexternal=0;

    for (m=0;m<mesh.ntriangles;m++) {
      keep=(elist[m]!=-1);
      switch(keep) {
        case 0:
          nexternal++;
          break;
        case 1:
          ninternal++;
          break;
        }
      }

    work[p].triangles=new triangle_t [ninternal];

    work[p].ntriangles=ninternal;

    ninternal=0;
    nexternal=0;

/* *-----------------------------------------------------------------------------
    screen elements*/
    for (m=0;m<mesh.ntriangles;m++) {
      keep=(elist[m]!=-1);
      switch(keep) {
        case 0:
          nexternal++;
          break;
        case 1:
          work[p].triangles[ninternal].ancestor=m;
          ninternal++;
          break;
        }
      }

/* *-----------------------------------------------------------------------------
    count vertices*/
    nndes=0;
    for (n=0;n<mesh.nvtxs;n++) used[n]=-1;
    for (m=0;m<work[p].ntriangles;m++) {
      for(i=0;i<3;i++) {
        n=mesh.triangles[work[p].triangles[m].ancestor].vertex[i];
        if(used[n]==-1) {
          used[n]=nndes;
          nndes++;
          }
        }
      }
/* *-----------------------------------------------------------------------------
    build vertices and elements tables*/
    work[p].nvtxs = nndes;
    work[p].vertices=new vertex_t[work[p].nvtxs];
    for (n=0; n<work[p].nvtxs; n++) work[p].vertices[n].null_value();
    for (n=0;n<mesh.nvtxs;n++) {
      if(used[n]!=-1) {
        work[p].vertices[used[n]].lon=mesh.vertices[n].lon;
        work[p].vertices[used[n]].lat=mesh.vertices[n].lat;
        work[p].vertices[used[n]].code=0;
        }
      }
    for (m=0;m<work[p].ntriangles;m++) {
      for(i=0;i<3;i++) {
        n=mesh.triangles[work[p].triangles[m].ancestor].vertex[i];
        if(used[n]!=-1) {
          work[p].triangles[m].vertex[i]=used[n];
          }
        else {
          printf("error...\n");
          }
        }
      }
    printf ("number of nodes    (splitted mesh %d): %6d\n",p,work[p].nvtxs);
    printf ("number of elements (splitted mesh %d): %6d\n",p,work[p].ntriangles);

    status=fe_e2n(&(work[p]));
    status=fe_edgetable(&work[p],0,0);
    sprintf(output,"%s.%2.2d.nei","internal",p);
    status=fe_savemesh(output,MESH_FILE_FORMAT_TRIGRID,work[p]);
    status=fe_codetable2(&work[p],0,1,0);
    sprintf(output,"%s.%2.2d.nei","internal",p);
    status=fe_savemesh(output,MESH_FILE_FORMAT_TRIGRID,work[p]);
    }

  delete[] used;

  printf("split completed\n");

  return(work);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  mesh_t *fe_partition_03(mesh_t mesh,int *partition,int npartition)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Partitioned meshes are taken from metis element list plus the
  neighbouring elements with overlapping kept as "just as needed".

  Principle:

  1- screen nodes that need to be truly solved in each partiction (starting
     from the last partition down to the first one)
     truly solved: interior node or main mesh's boundary node

  2- extent to node neighbours that need to be solved

  Suitable for parallel computing.

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
{
  int   k,m,mm,i,n,/*nn,*/p,status;
  int   nndes=0;
  int   ninternal,nexternal,keep;
  int   *used=NULL,*flag=NULL,*elist=NULL;
  int   *domain=NULL;
  mesh_t *work=NULL;
  char output[1024];
  connexion_t *table;

  printf ("number of nodes    (original mesh): %6d\n",mesh.nvtxs);
  printf ("number of elements (original mesh): %6d\n",mesh.ntriangles);

  used  =new int[mesh.nvtxs];
  flag  =new int[mesh.nvtxs];
  domain=new int[mesh.nvtxs];
  elist =new int[mesh.ntriangles];

  table=new connexion_t[npartition];
  work=new mesh_t[npartition];

  for(p=npartition-1;p>=0;p--) {
/* *-----------------------------------------------------------------------------
    screen vertices that will be truly solved inside each partition*/
    for (m=0;m<mesh.ntriangles;m++) {
      if(partition[m]==p) {
        for(i=0;i<3;i++) {
          n=mesh.triangles[m].vertex[i];
          domain[n]=p;
          }
        }
      }
    }

  for(p=0;p<npartition;p++) {

    for (n=0;n<mesh.nvtxs;n++) flag[n]=-1;
    for (m=0;m<mesh.ntriangles;m++) elist[m]=-1;

/* *-----------------------------------------------------------------------------
    screen vertices*/
    for (m=0;m<mesh.ntriangles;m++) {
      if(partition[m]==p) {
        elist[m]=0;
        for(i=0;i<3;i++) {
          n=mesh.triangles[m].vertex[i];
          flag[n]=0;
          }
        }
      }
/* *-----------------------------------------------------------------------------
    add vertices and elements neighbours*/
    for (n=0;n<mesh.nvtxs;n++) {
      if(domain[n]==p) {
/* *-----------------------------------------------------------------------------
        */
//         for(k=0;k<mesh.vertices[n].nngh;k++) {
//           nn=mesh.vertices[n].ngh[k];
//           if(flag[nn]==-1) flag[nn]=1;
//           }
        for(k=0;k<mesh.vertices[n].nelmts;k++) {
          mm=mesh.vertices[n].elmts[k];
          if(elist[mm]==-1) elist[mm]=1;
          }
        }
      }

    ninternal=0;
    nexternal=0;

    for (m=0;m<mesh.ntriangles;m++) {
      keep=(elist[m]!=-1);
      switch(keep) {
        case 0:
          nexternal++;
          break;
        case 1:
          ninternal++;
          break;
        }
      }

    work[p].triangles=new triangle_t [ninternal];

    work[p].ntriangles=ninternal;

    ninternal=0;
    nexternal=0;

/* *-----------------------------------------------------------------------------
    screen elements*/
    for (m=0;m<mesh.ntriangles;m++) {
      keep=(elist[m]!=-1);
      switch(keep) {
        case 0:
          nexternal++;
          break;
        case 1:
          work[p].triangles[ninternal].ancestor=m;
          ninternal++;
          break;
        }
      }

/* *-----------------------------------------------------------------------------
    count vertices*/
    nndes=0;
    for (n=0;n<mesh.nvtxs;n++) used[n]=-1;
    for (m=0;m<work[p].ntriangles;m++) {
      for(i=0;i<3;i++) {
        n=mesh.triangles[work[p].triangles[m].ancestor].vertex[i];
        if(used[n]==-1) {
          used[n]=nndes;
          nndes++;
          }
        }
      }
/* *-----------------------------------------------------------------------------
    build vertices and elements tables*/
    work[p].nvtxs = nndes;
    work[p].vertices=new vertex_t[work[p].nvtxs];
    for (n=0; n<work[p].nvtxs; n++) work[p].vertices[n].null_value();
    for (n=0;n<mesh.nvtxs;n++) {
      if(used[n]!=-1) {
        work[p].vertices[used[n]].lon=mesh.vertices[n].lon;
        work[p].vertices[used[n]].lat=mesh.vertices[n].lat;
        work[p].vertices[used[n]].code=0;
        work[p].vertices[used[n]].ancestor=n;
        }
      }
    for (m=0;m<work[p].ntriangles;m++) {
      for(i=0;i<3;i++) {
        n=mesh.triangles[work[p].triangles[m].ancestor].vertex[i];
        if(used[n]!=-1) {
          work[p].triangles[m].vertex[i]=used[n];
          }
        else {
          printf("error...\n");
          }
        }
      }
    printf ("number of nodes    (splitted mesh %d): %6d\n",p,work[p].nvtxs);
    printf ("number of elements (splitted mesh %d): %6d\n",p,work[p].ntriangles);

    status=fe_e2n(&(work[p]));
    status=fe_edgetable(&work[p],0,0);
    sprintf(output,"%s.%2.2d.nei","internal",p);
    status=fe_savemesh(output,MESH_FILE_FORMAT_TRIGRID,work[p]);
    status=fe_codetable2(&work[p],0,1,0);
    sprintf(output,"%s.%2.2d.nei","internal",p);
    status=fe_savemesh(output,MESH_FILE_FORMAT_TRIGRID,work[p]);
    }

  delete[] used;
  delete[] domain;

  printf("split completed\n");

  return(work);
}
