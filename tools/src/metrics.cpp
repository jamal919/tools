
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
#include <time.h>

#include "config.h"

#include "tools-define.h"
#include "tools-structures.h"
#include "tides.h"

#include "fe.h"
#include "geo.h"
#include "zapper.h"
#include "poc-time.h"

#include "legend.def"
#include "legend.h"

section_set_t gSectionSet;
FILE *echo=stdout;

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> T SquareDst(T lon1, T lat1, T lon2, T lat2)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  T SquareDstq;
  T dX,dY;
  double  R = 6.3675e8;

  lon2=geo_recale(lon2,lon1,180.0);

  dX = R*cos(0.5*(lat1+lat2))*(lon1-lon2);
  dY = R*(lat1-lat2);
  SquareDstq = dX*dX + dY*dY;

  return(SquareDstq);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int diag_section(mesh_t mesh, section_set_t SectionSet)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int neighbour,nlegends;
  int j,l,m,n,status;
  legend_t   *legends;
  legend02_t *legends02;
  char *comments;
  char lgdfile[1024];
  edge_t *edgePtr;
  int orientation;

  legends=new legend_t[SectionSet.n];

  for(l=0;l<SectionSet.n;l++) {;
    legends[l].ID=l;
    legends[l].Type=LGD_GRAPH;
    legends[l].T=0;
    legends[l].P=0;
    legends[l].X=0;
    legends[l].Y=0;
    }

  comments=new char[1024];
  comments[0]=0;

  strcat(comments,"#  0: node index \n");
  strcat(comments,"#  1: longitude of xover \n");
  strcat(comments,"#  2: latitude  of xover \n");
  strcat(comments,"#  3->0: normal-x\n");
  strcat(comments,"#  4->1: normal-y\n");

  legends02=new legend02_t[SectionSet.n];

  for(n=0;n<SectionSet.n;n++) {
    legends02[n].ID=n;
    legends02[n].np=SectionSet.sections[n].nsegments;
    legends02[n].nz=2;
    legends02[n].pen0=0;
    legends02[n].pen1=0;
    legends02[n].points=new point_t[legends02[n].np];
    legends[n].ptr=(char *) &(legends02[n]);
    for (m=0;m<legends02[n].np;m++) {
      legends02[n].points[m].z=new float[legends02[n].nz];
      }
    for (l=0;l<SectionSet.sections[n].nsegments;l++) {
      edgePtr     =&(SectionSet.sections[n].segments[l].edge);
//      orientation =SectionSet.sections[n].segments[l].orientation;
      legends02[n].points[l].N=l;
      legends02[n].points[l].t=edgePtr->lon*180./M_PI;
      legends02[n].points[l].p=edgePtr->lat*180./M_PI;
      j=0;
/**----------------------------------------------------------------------------
      right-hande side normal*/
      legends02[n].points[l].z[j++]=  edgePtr->Ty;
      legends02[n].points[l].z[j++]= -edgePtr->Tx;
      }
    }

  sprintf(lgdfile,"diagnostic.lgd");

  nlegends=SectionSet.n;

  status=lgd_save(lgdfile, legends, nlegends,NULL,comments);
  for(n=0;n<SectionSet.n;n++) {
    for (m=0;m<legends02[n].np;m++) {
      delete[] legends02[n].points[m].z;
      }
    delete[] legends02[n].points;
    }

  delete[] legends02;
  delete[] legends;


  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int find_edges(mesh_t mesh, int n1, int n2, edge_t *edge, int *node)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/**@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Find the edge connecting n1 and n2. Algorithm could be more simple

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
  int i, j, k, l, m, n, p, status;
  int j1,j2;
  int  *path,npath;
  double d,dmin,nx,ny;
  char   finished=0;

  if(mesh.vertices[n1].elmts==0) {
    printf("\nfind_edges needs mesh.vertices[b].elmts to be initialized before use  (fe_vertex_crosstables01)\n");
    status=fe_vertex_crosstables01(&mesh);
    }

  for(k=0;k<mesh.vertices[n1].nelmts;k++) {
    m=mesh.vertices[n1].elmts[k];
    l=vpos(m,mesh.vertices[n2].elmts,mesh.vertices[n2].nelmts);
    if(l<0) continue;
/*-----------------------------------------------------------------------------
    */
    j1=vpos(n1,mesh.triangles[m].vertex,3);
    j2=vpos(n2,mesh.triangles[m].vertex,3);
    for(i=0;i<3;i++) {
      p=mesh.triangles[m].edges[i];
//       printf("%d %d %d\n",p,mesh.edges[p].extremity[0],mesh.edges[p].extremity[1]);
//       if(p<0) {
//         printf("\ntroubles...\n");
//         }
//       if(p>mesh.nedges) {
//         printf("\ntroubles...\n");
//         }
      if((mesh.edges[p].extremity[0]==n2)&&(mesh.edges[p].extremity[1]==n1)) {
        break;
        }
      if((mesh.edges[p].extremity[1]==n2)&&(mesh.edges[p].extremity[0]==n1)) {
        break;
        }
      }
    *edge=mesh.edges[p];
    *node=p;
    if((j2-j1+3) % 3 ==1) {
/**-----------------------------------------------------------------------------
      element side goes from n1 to n2*/
      edge->Tx= -mesh.triangles[m].ny[i];
      edge->Ty=  mesh.triangles[m].nx[i];
      }
    else {
/**-----------------------------------------------------------------------------
      element side goes from n2 to n1*/
      edge->Tx=  mesh.triangles[m].ny[i];
      edge->Ty= -mesh.triangles[m].nx[i];
      }
//    printf("%d %d %d %lf %lf\n",n1,n2,p,edge->Tx,edge->Ty);
    goto found;
    }
  return(-1);

found:
  return(0);
  }


/*----------------------------------------------------------------------------*/
/** Finds the linear path connecting a pair
\date documented on 28 Jul 2011
\author Damien Allain

\param mesh
\param connexions indexes of extremities
\param *path array of indexes of nodes. The extremities are respectively the first and last elements of the array.

Old note :
Orientation correct ?
*/
/*----------------------------------------------------------------------------*/
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int find_path(mesh_t mesh, paire_t connexions, int *path)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, j, k, l, next, m, n, n1, n2, status;
  int    npath;
  double d,dmin,nx,ny;
  char   finished=0;

  n1=connexions.value[0];
  n2=connexions.value[1];
  npath=0;
  path[0]=n1;
  npath++;
  do {
    n=path[npath-1];
    dmin=1.e+10;
    for(k=0;k<mesh.vertices[n].nngh;k++) {
      m=mesh.vertices[n].ngh[k];
      if(m==n2) {
        finished=1;
        next=m;
        break;
        }
/**----------------------------------------------------------------------------
      seek closest neighbour*/
      d=geo_distance(mesh.vertices[m].lon,mesh.vertices[m].lat,mesh.vertices[n2].lon,mesh.vertices[n2].lat);
      if(d<dmin) {
        dmin=d;
        next=m;
        }
      }
    path[npath]=next;
    npath++;
    } while(finished==0);

  return(npath);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int section2D_LGP1xLGP1(mesh_t mesh, paire_t *connexions, int nconnexions, state2d_t state, parameter_t data, double *flux)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/**@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Compute transport through a given section

  LGP1 velocities

  Algorithm: seek the shortest path between two given nodes

  NOT bullet-proof, as it can go into dead-ends

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
  int i, j, k, l, next, m, n, n1, n2, status;
  int  *path,npath;
  int   orientation;
  edge_t edge;
  double H,d,dmin,nx,ny;
  char   finished=0;

  path=new int[mesh.nvtxs];

/**----------------------------------------------------------------------------
  find section path*/
  for(l=0;l<nconnexions;l++) {
    n1=connexions[l].value[0];
    n2=connexions[l].value[1];
    npath=find_path(mesh,  connexions[l],path);
    flux[l]=0;
    for(i=0;i<npath-1;i++) {
      n1=path[l];
      n2=path[l+1];
      status=find_edges(mesh, n1, n2, &edge, &orientation);
      if(status!=0) goto error;
/**----------------------------------------------------------------------------
      right-hande side normal*/
      nx= edge.Ty;
      ny=-edge.Tx;
      flux[l]+=0.5*edge.L*(nx*(state.u[n1]+state.u[n2])+ny*(state.v[n1]+state.v[n2]));
      }
    }

  zaparr(path);
  return(0);

error:
  zaparr(path);
  return(-1);

  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int section2D_NCP1xLGP1(mesh_t mesh, paire_t *connexions, int nconnexions, state2d_t state, parameter_t data, double *flux, double *z, double *u, double *v)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/**@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Compute transport through a given section

  LGP1 velocities

  Algorithm: seek the shortest path between two given nodes

  NOT bullet-proof, as it can go into dead-ends

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
  int i, j, k, l, next, m, n, n1, n2, status;
  int  *path,npath;
  int   node;
  edge_t edge;
  double d,dmin,nx,ny,H,F,L,zz;
  char   finished=0;

  path=new int[mesh.nvtxs];

/**----------------------------------------------------------------------------
  find section path*/
  for(l=0;l<nconnexions;l++) {
    n1=connexions[l].value[0];
    n2=connexions[l].value[1];
    npath=find_path(mesh,  connexions[l],path);
    flux[l]=0;
    z[l]=0;
    L=0;
    for(i=0;i<npath-1;i++) {
      n1=path[i];
      n2=path[i+1];
      status=find_edges(mesh, n1, n2, &edge, &node);
      if(status!=0) goto error;
/**----------------------------------------------------------------------------
      right-hand side normal*/
      nx= edge.Ty;
      ny=-edge.Tx;
      H=0.5*(state.H[n1]+state.H[n2]);
      zz=0.5*(state.z[n1]+state.z[n2]);
      F=H*edge.L*(nx*state.u[node]+ny*state.v[node]);
      flux[l]+=F;
      L+=edge.L;
      z[l]+=zz*edge.L;
      }
    z[l]/=L;
    }

  zaparr(path);
  return(0);

error:
  zaparr(path);
  return(-1);

  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int section2D(mesh_t mesh, int paire2D, paire_t *connexions, int nconnexions, state2d_t state, parameter_t data, double *flux, double *z, double *u, double *v)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/**@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Compute transport through a given section

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
  int i, j, k, l, next, m, n, n1, n2, status;

  switch(paire2D) {
    case NCP1xLGP0:
      check_error(-1, "not implemented yet", __LINE__, __FILE__, 1);
      break;

    case NCP1xLGP1:
      status=section2D_NCP1xLGP1(mesh, connexions, nconnexions, state, data, flux, z, u, v);
      break;

    case LGP1xLGP0:
      check_error(-1, "not implemented yet", __LINE__, __FILE__, 1);
      break;

    case LGP1xLGP1:
      status=section2D_LGP1xLGP1(mesh, connexions, nconnexions, state, data, flux);
      break;

    default:
      check_error(-1, "unknown paire", __LINE__, __FILE__, 1);
      break;
    }
  return(0);
  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void sampler_section_update(double t, mesh_t mesh, int paire2D, state2d_t state, parameter_t data)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, n, status;
  char filename[256];
  double time;
  FILE *out;
  section_set_t set=gSectionSet;
  paire_t *connexions;
  int nconnexions;
  double *elevation, *u, *v, *flux, sum;
//  int paire2D=LGP1xLGP1;
  char *output_path=".";

/*-----------------------------------------------------------------------
  Output z,u,v,ib in MKS units at sampling points*/
  if(gSectionSet.n == 0)
    return;

  nconnexions=gSectionSet.n;
  connexions= new paire_t[nconnexions];
  
  elevation= new double[nconnexions];
  u= new double[nconnexions];
  v= new double[nconnexions];
  flux= new double[nconnexions];
  
  for(i = 0; i < gSectionSet.n; i++) {
    connexions[i].value[0]=gSectionSet.sections[i].landmarks[0].node;
    connexions[i].value[1]=gSectionSet.sections[i].landmarks[1].node;
    }

  status=section2D(mesh, paire2D, connexions, nconnexions, state, data, flux, elevation, u, v);

  sum=0;
  for(i = 0; i < gSectionSet.n; i++) {
    sprintf(filename, "%s/section.%d", output_path, i + 1);
    out = fopen(filename, "a+");
//    time = cnestime(t);
    fprintf(out, "%12.6lf %lf %lf", t, flux[i],elevation[i]);
    sum+=flux[i];
    fprintf(out, "\n");
    fclose(out);
    }
    
  sprintf(filename, "%s/section.sum", output_path);
  out = fopen(filename, "a+");
  fprintf(out, "%12.6lf %lf", t, sum);
  fprintf(out, "\n");
  fclose(out);
    
  delete[] connexions;
  delete[] elevation;
  delete[] u;
  delete[] v;
  delete[] flux;

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int sampler_section_complete(mesh_t mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, n, n1, n2, node, status;
  char filename[256];
  double nx,ny;
  FILE *out;
  section_set_t set=gSectionSet;
  paire_t *connexions;
  int nconnexions;
  int npath, *path;
  edge_t edge;
/*-----------------------------------------------------------------------
  */
  if(gSectionSet.n == 0)
    return(0);

  nconnexions=gSectionSet.n;
  connexions= new paire_t[nconnexions];
  path=new int[mesh.nvtxs];

  for(n = 0; n < gSectionSet.n; n++) {
    connexions[n].value[0]=gSectionSet.sections[n].landmarks[0].node;
    connexions[n].value[1]=gSectionSet.sections[n].landmarks[1].node;
/**----------------------------------------------------------------------------
    find section path*/
    npath=find_path(mesh, connexions[n], path);
    gSectionSet.sections[n].segments=new segment_t[npath-1];
    gSectionSet.sections[n].nsegments=npath-1;
    for(i=0;i<npath-1;i++) {
      n1=path[i];
      n2=path[i+1];
      status=find_edges(mesh, n1, n2, &edge, &node);
      if(status!=0) goto error;
/**----------------------------------------------------------------------------
      still room for improvements*/
      gSectionSet.sections[n].segments[i].edge=edge;
      gSectionSet.sections[n].segments[i].node=node;
      }
    }

  delete[] connexions;
  delete[] path;

  status=diag_section( mesh,  gSectionSet);
  return(0);

error:
  return(-1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int sampler_section_init(mesh_t & mesh, char *SectionInputFile)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int    i, j, nitems, k, l, status;
  int    n1,n2;
  long   pos;
  double ds_tst, min_ds, plon, plat, dum;
  FILE   *in, *out;
  double time, model, start;
  char   filename[256],line[1024],msg[1024], sample_format[1024];
  char *output_path=".";

  if(strcmp(SectionInputFile, "NONE") == 0) {
    gSectionSet.n = 0;
    return(-1);
    }

/*----------------------------------------------------------------------
  Read sampling input file*/
  if((in = fopen(SectionInputFile, "r")) == NULL) {
    sprintf(msg,"Cannot open plot file %s",SectionInputFile);
    check_error(-1, msg, __LINE__, __FILE__, 1);
    }

  fgets(line, 1024, in);
  nitems=sscanf(line, "%d %s", &gSectionSet.n,sample_format);
/**----------------------------------------------------------------------
  patch to support old format*/
  if(nitems==1) {
    strcpy(sample_format,"?");
    }

  fprintf(echo, "\n Plotting %d sections\n", gSectionSet.n);

  gSectionSet.sections = new section_t [gSectionSet.n];
  if(gSectionSet.sections == NULL)
    check_error(-1, "allocation error", __LINE__, __FILE__, 1);

  for(i = 0; i < gSectionSet.n; i++) {
    fgets(line, 1024, in);
    nitems=sscanf(line,"%lf %lf %lf %lf", &gSectionSet.sections[i].landmarks[0].lon,
                                          &gSectionSet.sections[i].landmarks[0].lat,
                                          &gSectionSet.sections[i].landmarks[1].lon,
                                          &gSectionSet.sections[i].landmarks[1].lat);
/**----------------------------------------------------------------------
    conversion to radii*/
    for(k=0;k<2;k++) {
      gSectionSet.sections[i].landmarks[k].lon *= d2r;
      gSectionSet.sections[i].landmarks[k].lat *= d2r;
      gSectionSet.sections[i].landmarks[k].node = -1;
      }
    }
  fclose(in);

/*----------------------------------------------------------------------
  Find nearest vertex points for sampling*/
  for(i = 0; i < gSectionSet.n; i++) {
    for(k=0;k<2;k++) {
      min_ds = 1.e30;
      if(gSectionSet.sections[i].landmarks[k].node == -1) {
        plon = gSectionSet.sections[i].landmarks[k].lon;
        plat = gSectionSet.sections[i].landmarks[k].lat;
        for(j = 0; j < mesh.nvtxs; j++) {
/*----------------------------------------------------------------------
          distance in meters*/
          ds_tst = sqrt(SquareDst(plon, plat, mesh.vertices[j].lon*d2r, mesh.vertices[j].lat*d2r));
          if(ds_tst < min_ds) {
            min_ds = ds_tst;
            gSectionSet.sections[i].landmarks[k].node = j;
            }
          }
        }
      }
    }

  status=sampler_section_complete(mesh);

  for(i = 0; i < gSectionSet.n; i++) {
    sprintf(filename, "%s/section.%d", output_path, i + 1);
    if((out = fopen(filename, "w")) == NULL)
      continue;
    fclose(out);
    }

/**----------------------------------------------------------------------
  initialization routine, flushing ok*/
  fflush(echo);
  
  return(0);
}

