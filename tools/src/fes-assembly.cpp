

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

#define MAIN_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <config.h>

#include "tools-define.h"
#include "tools-structures.h"

#include "geo.h"
#include "fe.h"
#include "polygones.h"
#include "functions.h"
#include "mass.h"
#include "matrix.h"
#include "tides.def"
#include "tides.h"
#include "fe-proto.h"

#include "zapper.h"     /*  rutin.h contains common utility routines  */

class ascendance_t {
private :
public :
  int mesh;
  int P1node;
  int type;

  ascendance_t() {
    mesh=0;
    P1node=-1;
    type=-1;
    }
};

  
extern int cefmo_loadc1(char *name, int nndes, fcomplex *buffer);
extern int save_TidalOBCs(spectrum_t s, tidalOBC_t *data, int ndata);

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_gauss7(mesh_t mesh,float *in, double *out)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *-----------------------------------------------------------------------

  Solve : M x = INT()

------------------------------------------------------------------------*/
{
  int j,k,m,n,status;
  double pseudo,C,z;
//  extern int LinearSystem_solve(hypermatrix_t, double *);

  discretisation_t descriptor=mesh.LGP2descriptor;
  double beta[16][6];
  double gauss_x[16], gauss_y[16], gauss_w[16];
  int    gauss_n=7;

  gauss_init(gauss_n,gauss_x,gauss_y,gauss_w);
  for(k=0; k<gauss_n; k++) {
    fe_LGP2base(gauss_x[k],gauss_y[k],beta[k]);
    }

  for(n = 0; n < descriptor.nnodes; n++)
    out[n] = 0;

  if(in==0) return(-1);

  for (m=0;m<mesh.ntriangles;m++) {
    pseudo=mesh.triangles[m].Area;
    C=mesh.triangles[m].cosinus;
    for(k=0; k<gauss_n; k++) {
      for(j=0;j<descriptor.nnpe;j++) {
        n=descriptor.NIbE[m][k];
        z=in[m*7+k]*beta[k][j]*C;
        out[n]+=gauss_w[k]*z;
        }
      }
    }
  status = LinearSystem_solve(descriptor.massmatrix,out);
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int cefmo_loadIPG7(char *name, int nndes, float *buffer)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int   k,n;
  FILE *in;

  in=fopen(name,"r");
  if(in == NULL) {
    printf("unable to open %s\n",name);
    return (-1);
    }

  for(n=0;n<nndes;n++) {
    for(k=0;k<7;k++) {
      fscanf(in,"%f",&(buffer[n*7+k]));
      }
    }

  fclose(in);
  return (0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int P2merge(mesh_t *merged,mesh_t work[2],double dmax)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int   count,k,m,m0,i,j,n;
  int   ii,mm;
  int   cw[3][3] ={{4,3,2},{0,5,4},{2,1,0}};
  int   ccw[3][3]={{2,3,4},{4,5,0},{0,1,2}};
  int   shared,check;
  int   *LGP2index;

  merged->LGP2descriptor.nnodes=merged->nedges+merged->nvtxs;

  merged->LGP2descriptor.NIbE=new int *[merged->ntriangles];
  for(m=0;m<merged->ntriangles;m++) {
    merged->LGP2descriptor.NIbE[m]=new int[6];
    for(i=0;i<6;i++) {
      merged->LGP2descriptor.NIbE[m][i]=-1;
      }
    }

  LGP2index=new int[work[1].LGP2descriptor.nnodes];
  for(n=0;n<work[1].LGP2descriptor.nnodes;n++) {
    LGP2index[n]=-1;
    }

/* *----------------------------------------------------------------------------
  copy full LGP2 descriptor from mesh 0*/
  for(m=0;m<work[0].ntriangles;m++) {
    for(i=0;i<6;i++) {
      n=work[0].LGP2descriptor.NIbE[m][i];
      merged->LGP2descriptor.NIbE[m][i]=work[0].LGP2descriptor.NIbE[m][i];
      }
    }

/* *----------------------------------------------------------------------------
  extend numerotation along shared limits*/
  shared=0;
  check=0;
  for(m=0;m<work[0].ntriangles;m++) {
    for(i=0;i<3;i++) {
      n=work[0].triangles[m].edges[i];
/* *----------------------------------------------------------------------------
      skip if not boundary edge in mesh 0*/
      if(work[0].edges[n].nshared==2) continue;
      n=merged->triangles[m].edges[i];
/* *----------------------------------------------------------------------------
      skip if boundary edge in merged mesh*/
      if(merged->edges[n].nshared==1) continue;
      k=((m==merged->edges[n].shared[0]) ? (1) : (0));
      if(k!=1){
        printf("anomaly...\n");
        }
      mm=merged->edges[n].shared[k];
      ii=merged->edges[n].eindex[k];
      if(merged->edges[n].eindex[0]!=i){
        printf("anomaly...\n");
        }
      if(mm<work[0].ntriangles){
        printf("anomaly...\n");
        }
      shared++;
      for(j=0;j<3;j++) {
        if(merged->LGP2descriptor.NIbE[mm][ccw[ii][j]]==-1) check++;
        merged->LGP2descriptor.NIbE[mm][ccw[ii][j]]=merged->LGP2descriptor.NIbE[m][cw[i][j]];
//        printf("%d ",merged->LGP2descriptor.NIbE[mm][ccw[ii][j]]);
        }
//      printf("\n");
      }
    }

  printf("number of shared edges: %d\n",shared);

/* *----------------------------------------------------------------------------
  identify mesh-1 LGP2 nodes already used*/
  shared=0;
  m0=work[0].ntriangles;
  for(m=0;m<work[1].ntriangles;m++) {
    for(i=0;i<6;i++) {
      if(merged->LGP2descriptor.NIbE[m+m0][i]!=-1) {
/* *----------------------------------------------------------------------------
        LGP2 descriptor already initialized from mesh 0*/
        if(LGP2index[work[1].LGP2descriptor.NIbE[m][i]]==-1) shared++;
        LGP2index[work[1].LGP2descriptor.NIbE[m][i]]=merged->LGP2descriptor.NIbE[m+m0][i];
        }
      }
    }

  printf("number of shared nodes: %d\n",shared);
  count=work[0].LGP2descriptor.nnodes;

/* *----------------------------------------------------------------------------
  copy partial LGP2 descriptor from mesh 1*/
  m0=work[0].ntriangles;
  for(m=0;m<work[1].ntriangles;m++) {
    for(i=0;i<6;i++) {
      if(merged->LGP2descriptor.NIbE[m+m0][i]==-1) {
/* *----------------------------------------------------------------------------
        LGP2 descriptor not initialized from mesh 0*/
        n=LGP2index[work[1].LGP2descriptor.NIbE[m][i]];
        if(n==-1) {
/* *----------------------------------------------------------------------------
          LGP2 node not yet initialized, do it 0*/
          merged->LGP2descriptor.NIbE[m+m0][i]=count;
          LGP2index[work[1].LGP2descriptor.NIbE[m][i]]=count;
          if(count==merged->LGP2descriptor.nnodes){
            printf("LGP2 nodes overflow...\n");
            return(-1);
            }
          count++;
          }
        else {
/* *----------------------------------------------------------------------------
          LGP2 node already exists*/
          merged->LGP2descriptor.NIbE[m+m0][i]=n;
          }
        }
      }
    }

  for(n=0;n<work[1].LGP2descriptor.nnodes;n++) {
    if(LGP2index[n]==-1) {
      printf("P2merge error: LGP2 node %d has no equivalent in merged mesh\n",n);
      return(-1);
      }
    }
  if(count!=merged->LGP2descriptor.nnodes) {
    printf("P2merge error: actual LGP2 nodes count %d differs from predicted one %d\n",count,merged->LGP2descriptor.nnodes);
    return(-1);
    }

  merged->LGP2descriptor.nodes=new node_t[merged->LGP2descriptor.nnodes];
  for(n=0;n<work[0].LGP2descriptor.nnodes;n++) {
    work[0].LGP2descriptor.nodes[n].child=n;
    }
  for(n=0;n<work[1].LGP2descriptor.nnodes;n++) {
    work[1].LGP2descriptor.nodes[n].child=LGP2index[n];
    }

  merged->LGP2descriptor.type=LGP2;
  merged->LGP2descriptor.nnpe=6;
  merged->LGP2descriptor.nelmts=merged->ntriangles;

  printf("P2 merge completed\n");
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int P1merge(mesh_t *merged,mesh_t work[2],double dmax)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int   k,l,m,m0,i,n1,n2,n,status;
  int   nndes;
  int   merging;
  int   *P1index,shared;
  double d,dmin,t1,p1,t2,p2;

  merged->destroy();
  merged->type=0;

  merged->ntriangles=work[0].ntriangles+work[1].ntriangles;

  merged->triangles=new triangle_t[merged->ntriangles];

  P1index=new int[work[1].nvtxs];
  for(n=0;n<work[1].nvtxs;n++) {
    P1index[n]=-1;
    }

  shared=0;
  nndes=work[0].nvtxs;
  for(n2=0;n2<work[1].nvtxs;n2++) {
//    if(work[1].vertices[n2].code==0) {
/* *----------------------------------------------------------------------------
    interior or inner rigid boundary node*/
    if(work[1].vertices[n2].code!=1) {
      P1index[n2]=nndes;
      nndes++;
      continue;
      }
    merging=0;
    t2=work[1].vertices[n2].lon;
    p2=work[1].vertices[n2].lat;
    dmin=dmax;
/* *----------------------------------------------------------------------------
    check nodes in mesh 1 against exterior boundary nodes in mesh 0*/
    for(l=0;l<work[0].nlimits;l++) {
      for(k=0;k<work[0].limits[l].nvertex;k++) {
        n1=work[0].limits[l].vertex[k];
        t1=work[0].vertices[n1].lon;
        p1=work[0].vertices[n1].lat;
        d=geo_distance(t1,p1,t2,p2);
        if(d<dmin) {
          P1index[n2]=n1;
          dmin=d;
          }
        }
      }
//done:
    merging=(dmin<dmax);
    if(merging) {
/* *----------------------------------------------------------------------------
      boundary node is shared by both meshes*/
      shared++;
      }
    else {
/* *----------------------------------------------------------------------------
      boundary node is NOT shared*/
      P1index[n2]=nndes;
      nndes++;
      }
    }

  printf("number of shared vertices: %d\n",shared);
  for(n=0;n<work[0].nvtxs;n++) {
    work[0].vertices[n].child=n;
    }
  for(n=0;n<work[1].nvtxs;n++) {
    work[1].vertices[n].child=P1index[n];
    }

  merged->nvtxs=nndes;
  merged->vertices=new vertex_t[merged->nvtxs];
  for (n=0; n<merged->nvtxs; n++) merged->vertices[n].null_value();

/* *----------------------------------------------------------------------------
  keep all nodes from first mesh*/
  for(n=0;n<work[0].nvtxs;n++) {
    merged->vertices[n].lon =work[0].vertices[n].lon;
    merged->vertices[n].lat =work[0].vertices[n].lat;
    merged->vertices[n].ancestor=n;
    merged->vertices[n].code=0;
    }

  for(n=0;n<work[1].nvtxs;n++) {
    m=P1index[n];
    if(m<work[0].nvtxs) continue;
    merged->vertices[m].lon =work[1].vertices[n].lon;
    merged->vertices[m].lat =work[1].vertices[n].lat;
    merged->vertices[m].ancestor=n;
    merged->vertices[m].code=0;
    }

  for(m=0;m<work[0].ntriangles;m++) {
    for(i=0;i<3;i++) {
      merged->triangles[m].vertex[i]=work[0].triangles[m].vertex[i];
      }
    merged->triangles[m].ancestor=m;
    work[0].triangles[m].child=m;
    }

  m0=work[0].ntriangles;
  for(m=0;m<work[1].ntriangles;m++) {
    for(i=0;i<3;i++) {
      n=P1index[work[1].triangles[m].vertex[i]];
      merged->triangles[m+m0].vertex[i]=n;
      }
    merged->triangles[m+m0].ancestor=m;
    work[1].triangles[m].child=m+m0;
    }

  status= fe_e2n(merged);

  status= init_edge_table(merged);
  status= fe_savemesh("tmp01.nei",MESH_FILE_FORMAT_TRIGRID,*merged);

  status= fe_codetable2(merged,0,1,0);

  for(n=0;n<merged->nedges;n++) {
    merged->edges[n].ancestor=-1;
    }
  for(m=0;m<work[0].ntriangles;m++) {
    for(i=0;i<3;i++) {
      n1=merged->triangles[m].edges[i];
      n2=work[0].triangles[m].edges[i];
      merged->edges[n1].ancestor=n2;
      work[0].edges[n2].child=n1;
      }
    }
  m0=work[0].ntriangles;
  for(m=0;m<work[1].ntriangles;m++) {
    for(i=0;i<3;i++) {
      n1=merged->triangles[m].edges[i];
      n2=work[1].triangles[m].edges[i];
      if(merged->edges[n1].ancestor==-1) {
        merged->edges[n1].ancestor=n2;
        }
      work[1].edges[n2].child=n1;
      }
    }

  status=fe_element_crosstables(*merged);

  status= fe_savemesh("merged.nei",MESH_FILE_FORMAT_TRIGRID,*merged);

  printf("merge completed\n");

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int readmesh_MODULEFP2 (const char *filename, mesh_t *mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  FILE *in;
  int i,m,n,nnodes,nn,np;
  int nitems,n1,n2,n3,n4,n5,n6,count,status;
  char c;
  char *s,dum[32];
  char line[1024]/*,separator[32]*/;
  discretisation_t LGP2_C_numerics;

/*------------------------------------------------------------------------------
  element and node list format*/
  in = fopen(filename, "r");
  if(in==0)
      return(-1);

/*------------------------------------------------------------------------------
  rough estimate of the maximum count of nodes*/
  count=0;
  do {
    fgets(line,1024,in);
    if(feof(in)) break;
    s=strstr(line,"(NTRI ) :");
    } while(s==0);
  sscanf((s+9),"%d", &(mesh->ntriangles));
  mesh->triangles= new triangle_t [mesh->ntriangles];
  if(mesh->triangles == NULL) {
    printf("allocation error vertex set\n");
    return(-1);
    }

  do {
    fgets(line,1024,in);
    if(feof(in)) break;
    s=strstr(line,"(NOE  ) :");
    } while(s==0);
  sscanf((s+9),"%d", &nnodes);

  LGP2_C_numerics.nnodes=nnodes;
  LGP2_C_numerics.type=LGP2_C;
  LGP2_C_numerics.NIbE=new int*[mesh->ntriangles];
  for(m=0;m<mesh->ntriangles;m++) {
    LGP2_C_numerics.NIbE[m]=new int[6];
    }

  do {
    fgets(line,1024,in);
    if(feof(in)) break;
    s=strstr(line,"(NP  ) :");
    } while(s==0);
  sscanf((s+9),"%d", &(mesh->nvtxs));
  mesh->vertices= new vertex_t [mesh->nvtxs];
  if(mesh->vertices == NULL) {
    printf("allocation error vertex set\n");
    return(-1);
    }

  do {
    fgets(line,1024,in);
    if(feof(in)) break;
    s=strstr(line,"| POINT |");
    } while(s==0);
  fgets(line,1024,in);

  for(i=0;i<mesh->nvtxs;i++) {
//     nitems=fscanf(in, "%s %d %s %lf %s %lf %s", &separator, &n1,
//                                                 &separator, &mesh->vertices[i].lon,
//                                                 &separator, &mesh->vertices[i].lat,
//                                                 &separator);
//     if(nitems!=7) {
//       printf("troubles...\n");
//       }
    c=fgetc(in);
    while(c!='|') c=fgetc(in);
    nitems=fscanf(in, "%d", &n1);
    if(nitems!=1) {
      printf("troubles...\n");
      }
    c=fgetc(in);
    while(c!='|') c=fgetc(in);
    nitems=fscanf(in, "%lf", &mesh->vertices[i].lon);
    if(nitems!=1) {
      printf("troubles...\n");
      }
    c=fgetc(in);
    while(c!='|') c=fgetc(in);
    nitems=fscanf(in, "%lf", &mesh->vertices[i].lat);
    if(nitems!=1) {
      printf("troubles...\n");
      }
    c=fgetc(in);
    while(c!='|') c=fgetc(in);
    }

  for(m=0;m<mesh->ntriangles;m++) {
    do {
      fgets(line,1024,in);
      if(feof(in)) break;
      s=strstr(line,"NOMBRE DE  NOEUDS :");
      } while(s==0);
//    sscanf((s+19),"%d", &nn);
    nitems=sscanf((s+19),"%d %d %d %d %d %d %d", &nn,&n1,&n3,&n5,&n2,&n4,&n6);
    if(nitems!=7) {
      nitems=sscanf((s+19),"%d %s %d %d %d %d %d %d", &nn,dum,&n1,&n3,&n5,&n2,&n4,&n6);
      if(nitems!=8) printf("troubles...\n");
      }
    LGP2_C_numerics.NIbE[m][0]=n1-1;
    LGP2_C_numerics.NIbE[m][1]=n2-1;
    LGP2_C_numerics.NIbE[m][2]=n3-1;
    LGP2_C_numerics.NIbE[m][3]=n4-1;
    LGP2_C_numerics.NIbE[m][4]=n5-1;
    LGP2_C_numerics.NIbE[m][5]=n6-1;
    do {
      fgets(line,1024,in);
      if(feof(in)) break;
      s=strstr(line,"NOMBRE DE  POINTS :");
      } while(s==0);
    nitems=sscanf((s+19),"%d %d %d %d", &np,&n1,&n2,&n3);
    if(nitems!=4) {
      nitems=sscanf((s+19),"%d %s %d %d %d", &np,dum,&n1,&n2,&n3);
      if(nitems!=5) printf("troubles...\n");
      }
    mesh->triangles[m].vertex[0]=n1-1;
    mesh->triangles[m].vertex[1]=n2-1;
    mesh->triangles[m].vertex[2]=n3-1;
    }

  LGP2_C_numerics.nodes=new node_t[LGP2_C_numerics.nnodes];
  for(m=0;m<mesh->ntriangles;m++) {
    for(i=0;i<3;i++) {
      n=LGP2_C_numerics.NIbE[m][2*i];
      int nLGP1=mesh->triangles[m].vertex[i];
      LGP2_C_numerics.nodes[n].lon=mesh->vertices[nLGP1].lon;
      LGP2_C_numerics.nodes[n].lat=mesh->vertices[nLGP1].lat;
      }
    for(i=0;i<3;i++) {
      n=LGP2_C_numerics.NIbE[m][2*i+1];
      int n1=LGP2_C_numerics.NIbE[m][2*i];
      int n2=LGP2_C_numerics.NIbE[m][(2*i+2) % 6];
      LGP2_C_numerics.nodes[n].lon=0.5*(LGP2_C_numerics.nodes[n1].lon+LGP2_C_numerics.nodes[n2].lon);
      LGP2_C_numerics.nodes[n].lat=0.5*(LGP2_C_numerics.nodes[n1].lat+LGP2_C_numerics.nodes[n2].lat);
      }
    }
  mesh->LGP2descriptor=LGP2_C_numerics;
  mesh->LGP2descriptor.type=LGP2;
  mesh->LGP2descriptor.nnpe=6;
  mesh->LGP2descriptor.nelmts=mesh->ntriangles;

  status=fe_e2n(mesh);

  fclose(in);

  mesh->type=0;

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void RHS_topo(mesh_t mesh, discretisation_t descriptor, double *rhs, float *topo)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/* *----------------------------------------------------------------------

  Mass matrix for P2 integrale scalar product

----------------------------------------------------------------------*/
  int i, k, m, n, status;
  double C, area;
  
  double gauss_x[16], gauss_y[16], gauss_w[16];
  int    gauss_n=7;
  double **base;
  gauss_init(gauss_n,gauss_x,gauss_y,gauss_w);

  int discretisation=descriptor.type;
  
  for(n = 0; n < descriptor.nnodes; n++)
    rhs[n] = 0;
  
  base=new double*[gauss_n];
  for (k=0;k<gauss_n;k++) {
    base[k]=new double[descriptor.nnpe];
    status=fe_LGPbase(mesh, discretisation, gauss_x[k], gauss_y[k], base[k]);
    }

  for(m = 0; m < mesh.ntriangles; m++) {
    C=mesh.triangles[m].cosinus;
    area = mesh.triangles[m].Area;
    for(i = 0; i < descriptor.nnpe; i++) {
      n = descriptor.NIbE[m][i];
      for (k=0;k<gauss_n;k++) {
        rhs[n] +=topo[m*gauss_n+k]*area*C*gauss_w[k]*base[k][i]*2.0;
        }
      }
    }

  for (k=0;k<gauss_n;k++) {
    delete[] base[k];
    }
  delete[] base;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int load_CEFMOobc(spectrum_t s, tidalOBC_t **data, int *ndata)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,k,nn;
  int nitems;
  FILE *out=NULL;
  char filename[1024];
  float a[3],G[3];
  mesh_t mesh;
  
  for(k=0;k<s.n;k++) {
    sprintf(filename,"m2.cl.p2",s.waves[k].name);
    out=fopen(filename,"r");
    if (out ==NULL) {
      fprintf(stderr,"cannot open obc file : %s \n",filename);
      return(-1);
      }
//    fprintf(out,"%s\n",s.waves[k].name);
    nitems=fscanf(out,"%d",ndata);
    if(*data==0) {
      *data=new tidalOBC_t[*ndata];
      for(i=0; i<*ndata; i++) {
        (*data)[i].Htide.init_polar(s.n);
        }
      }
    for(i=0; i<*ndata; i++) {
/*------------------------------------------------------------------------------
      elevation*/
      nitems=fscanf(out,"%d %f %f\n",&nn,&a[0],&G[0]);
      (*data)[i].node=nn-1;
      (*data)[i].Htide.a[k]=a[0];
      (*data)[i].Htide.G[k]=G[0];
      }
    fclose(out);
    }
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int convert_obc()

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,n,status;
  mesh_t mesh;
  const char *meshfile="tri-45.p2";
  tidalOBC_t *data=0;
  int ndata;
  spectrum_t spectrum;
  discretisation_t descriptor;
  
  status=readmesh_MODULEFP2(meshfile,&mesh);
  
  descriptor=mesh.LGP2descriptor;
  
  spectrum.init(1);
  status=spectrum.add(wM2);
  
  status=load_CEFMOobc(spectrum, &data, &ndata);
  
  for(i=0; i<ndata; i++) {
    n=data[i].node;
    data[i].lon=descriptor.nodes[n].lon;
    data[i].lat=descriptor.nodes[n].lat+45.0;
    }
  
  status=save_TidalOBCs(spectrum, data, ndata);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int j,k,l,n;
  int format,iteration,nmesh,nwave;
  char *keyword;
  char *meshfile[20]={NULL,NULL},*fesfile[20]={NULL,NULL},*output=NULL;
  char *AtlasPath[20];
  char *AtlasConvention=NULL,*AtlasRootPath=NULL,*wave[20];
  mesh_t mesh[20],work[2],tmp,final,merged;
  double dmax=0,cmax=0, mask;
  complex<float> *cbuffer[10],*fes;
  float  *rbuffer[100],*topo;
  double *dLGP2,*dLGP0;
  hypermatrix_t *matrix;

  fct_echo(argc,argv);

  n=1;
  while (n < argc) {
    keyword=strdup(argv[n]);
    switch (keyword[0]) {
      case '-':
        switch (keyword[1]) {
          case 'a' :
            AtlasRootPath = argv[n+1];
            n++;
            n++;
            break;

          case 'c' :
            AtlasConvention = argv[n+1];
            n++;
            n++;
            break;

          case 'd' :
            sscanf(argv[n+1],"%lf",&dmax);
            n++;
            n++;
            break;

          case 'e' :
            sscanf(argv[n+1],"%lf",&cmax);
            n++;
            n++;
            break;

          case 'o' :
            output= strdup(argv[n+1]);
            n++;
            n++;
            break;

          default:
            STDOUT_BASE_LINE("unknown option %s\n",keyword);
            exit(-1);
          }
        break;

      default:
        meshfile[nmesh]= strdup(argv[n]);
        n++;
        nmesh++;
        break;
      }
    free(keyword);
    }
/**----------------------------------------------------------------------------
  for canal obc's, nid de coucou*/
//  convert_obc();

  if(output==0) {
    printf("no filename specified for output (-o [filename]); abort...\n");
    goto error;
    }

  if(AtlasConvention==0) {
    printf("no filename extension specified for FE input (-c [extension]); abort...\n");
    goto error;
    }

  if(AtlasRootPath==0) {
    AtlasRootPath=strdup(".");
    goto error;
    }

  nmesh=12;

  for(k=0;k<nmesh;k++) {
    AtlasPath[k]=new char[strlen(AtlasRootPath)+32];
    }
  k=0;
  sprintf(AtlasPath[k],"%s/%s",AtlasRootPath,"Atlantic");
  k++;
  sprintf(AtlasPath[k],"%s/%s",AtlasRootPath,"Arctic");
  k++;
  sprintf(AtlasPath[k],"%s/%s",AtlasRootPath,"Europe");
  k++;
  sprintf(AtlasPath[k],"%s/%s",AtlasRootPath,"Mediter");
  k++;
  sprintf(AtlasPath[k],"%s/%s",AtlasRootPath,"Labrador");
  k++;
  sprintf(AtlasPath[k],"%s/%s",AtlasRootPath,"Weddel");
  k++;
  sprintf(AtlasPath[k],"%s/%s",AtlasRootPath,"Patagonia");
  k++;
  sprintf(AtlasPath[k],"%s/%s",AtlasRootPath,"Indian");
  k++;
  sprintf(AtlasPath[k],"%s/%s",AtlasRootPath,"Indonesia");
  k++;
  sprintf(AtlasPath[k],"%s/%s",AtlasRootPath,"PacificS");
  k++;
  sprintf(AtlasPath[k],"%s/%s",AtlasRootPath,"Ross");
  k++;
  sprintf(AtlasPath[k],"%s/%s",AtlasRootPath,"PacificN");
  
//   k=0;
//   AtlasPath[k]=strdup("/mnt/pc-titanic/data/tides/FES2002/basins/Atlantic");
//   k++;
//   AtlasPath[k]=strdup("/mnt/pc-titanic/data/tides/FES2002/basins/Arctic");
//   k++;
//   AtlasPath[k]=strdup("/mnt/pc-titanic/data/tides/FES2002/basins/Europe");
//   k++;
//   AtlasPath[k]=strdup("/mnt/pc-titanic/data/tides/FES2002/basins/Mediter");
//   k++;
//   AtlasPath[k]=strdup("/mnt/pc-titanic/data/tides/FES2002/basins/Labrador");
//   k++;
//   AtlasPath[k]=strdup("/mnt/pc-titanic/data/tides/FES2002/basins/Weddel");
//   k++;
//   AtlasPath[k]=strdup("/mnt/pc-titanic/data/tides/FES2002/basins/Patagonia");
//   k++;
//   AtlasPath[k]=strdup("/mnt/pc-titanic/data/tides/FES2002/basins/Indian");
//   k++;
//   AtlasPath[k]=strdup("/mnt/pc-titanic/data/tides/FES2002/basins/Indonesia");
//   k++;
//   AtlasPath[k]=strdup("/mnt/pc-titanic/data/tides/FES2002/basins/PacificS");
//   k++;
//   AtlasPath[k]=strdup("/mnt/pc-titanic/data/tides/FES2002/basins/Ross");
//   k++;
//   AtlasPath[k]=strdup("/mnt/pc-titanic/data/tides/FES2002/basins/PacificN");

  for(k=0;k<nmesh;k++) {
    meshfile[k]=new char[strlen(AtlasPath[k])+8];
    sprintf(meshfile[k],"%s/%s",AtlasPath[k],"tri.p2");
    }

//   status=fe_init(&final);
//   status=fe_init(&tmp);

  format=MESH_FILE_FORMAT_MODULEF_P2;

  for (k=0;k<nmesh;k++) {
    if(meshfile[k] != NULL) {
//      status=fe_readmesh(meshfile[k],MESH_FILE_FORMAT_MODULEF_P2,&mesh[k]);
      status=readmesh_MODULEFP2(meshfile[k],&mesh[k]);
      if(status!=0) {
        printf("unable to read the original mesh in %s\n",meshfile[k]);
        goto error;
        }
      status=fe_geometry(&mesh[k]);
//      status=fe_list(&mesh[k]);
//       if(status!=0) {
//         printf("unable to build the element list from the original mesh\n");
//         goto error;
//         }
      status= fe_edgetable (&mesh[k],0,0);
      status= fe_codetable2(&mesh[k],0,1,0);
      }
    else {
      printf("no %dth mesh filename specified; abort...\n",k);
      goto error;
      }
    }

/**----------------------------------------------------------------------------
    temporary, for testing purposes...*/
//  nmesh=2;
/** temporary, for testing purposes...
------------------------------------------------------------------------------*/
  if(nmesh>1) {
    work[0]=mesh[0];
    work[1]=mesh[1];
    status= P1merge(&merged,work,(double) 2.5);
    status= P2merge(&merged,work, dmax);
    if(status!=0) {
      printf("LGP2 merge failed (mesh %d)\n",1);
      goto error;
      }
    tmp=merged;
    for (k=2;k<nmesh;k++) {
      work[0]=tmp;
      work[1]=mesh[k];
      status= P1merge(&merged,work,(double) 2.5);
      status= P2merge(&merged,work, dmax);
      if(status!=0) {
        printf("LGP2 merge failed (mesh %d)\n",k);
        goto error;
        }
      tmp.destroy();
      tmp=merged;
      }
    }
  else {
    merged=mesh[0];
    }
 
  final=merged;
  status= fe_savemesh("final.nei",MESH_FILE_FORMAT_TRIGRID, final);

/* *----------------------------------------------------------------------------
  LGP2 mass matrix*/
  status=fe_geometry(&final);
  status=fe_element2neighbours (final, final.LGP2descriptor);
  matrix=new hypermatrix_t;
  matrix->ordering=new ordering_t;
  status=dMassMatrix(final, final.LGP2descriptor, matrix);

  k=0;
  wave[k]=strdup("M2");
  k++;
  wave[k]=strdup("S2");
  k++;
  wave[k]=strdup("N2");
  k++;
  wave[k]=strdup("K2");
  k++;
  wave[k]=strdup("2N2");
  k++;
  wave[k]=strdup("K1");
  k++;
  wave[k]=strdup("O1");
  k++;
  wave[k]=strdup("Q1");
  k++;
  wave[k]=strdup("P1");

  nwave=k+1;

  fes=new complex<float>[final.LGP2descriptor.nnodes];

  dLGP2=new double[final.LGP2descriptor.nnodes];
  dLGP0=new double[final.LGP2descriptor.nelmts];

/* *----------------------------------------------------------------------------
  Treat bathymetry */
  topo=new float[final.LGP2descriptor.nelmts*7];
  for (l=0;l<nwave;l++) {
    for(k=0;k<nmesh;k++) {
      fesfile[k]=new char[strlen(AtlasPath[k])+strlen("topo.IPG7")+2];
      sprintf(fesfile[k],"%s/%s",AtlasPath[k],"topo.IPG7");
      }
    for(n=0;n<final.LGP2descriptor.nnodes;n++) {
      fes[n]=complex<float>(1./0.,1./0.);
      }
    output=new char[256];
//    sprintf(output,"%s.FES2004-optimal.nc",wave[l]);
    sprintf(output,"%s.FES2004-prior.nc",wave[l]);

    for (k=0;k<nmesh;k++) {
      rbuffer[k]=new float[mesh[k].LGP2descriptor.nelmts*7];
      status=cefmo_loadIPG7(fesfile[k],mesh[k].LGP2descriptor.nelmts,rbuffer[k]);
      for(n=0;n<mesh[k].ntriangles;n++) {
        for(j=0;j<7;j++) {
          topo[mesh[k].triangles[n].child*7+j]=rbuffer[k][n*7+j];
          }
        }
      delete[] rbuffer[k];
      }
    double *rhs=new double[final.LGP2descriptor.nnodes];
    RHS_topo(final, final.LGP2descriptor, rhs, topo);
    status = LinearSystem_solve(*matrix,rhs);

    for(n=0;n<final.LGP2descriptor.nnodes;n++) {
//       if(fes[n]==complex<float>(1./0.,1./0.)) {
//         printf("trouble...\n");
//         }
      dLGP2[n]=abs(rhs[n]);
      }
    iteration=0;
    status=archiving_UGdummy2D((const char*) output, final, "bathymetry", "m", dLGP2, mask, iteration,(int) LGP2);

    for(k=0;k<nmesh;k++) {
      delete[] fesfile[k];
      }
    delete[] output;
    }


/* *----------------------------------------------------------------------------
  Treat elevations */
//  AtlasConvention=strdup("den.ass.TP_E2_TG.mog2d-1");
  AtlasConvention=strdup("den.b.wdp+0.03");
  for (l=0;l<nwave;l++) {
    for(k=0;k<nmesh;k++) {
      fesfile[k]=new char[strlen(AtlasPath[k])+strlen(wave[l])+strlen(AtlasConvention)+3];
      sprintf(fesfile[k],"%s/%s.%s",AtlasPath[k],wave[l],AtlasConvention);
      }
    for(n=0;n<final.LGP2descriptor.nnodes;n++) {
      fes[n]=complex<float>(1./0.,1./0.);
      }
    output=new char[256];
//    sprintf(output,"%s.FES2004-optimal.nc",wave[l]);
    sprintf(output,"%s.FES2004-prior.nc",wave[l]);

    for (k=0;k<nmesh;k++) {
      cbuffer[k]=new complex<float>[mesh[k].LGP2descriptor.nnodes];
      status=cefmo_loadc1(fesfile[k],mesh[k].LGP2descriptor.nnodes,cbuffer[k]);
      for(n=0;n<mesh[k].LGP2descriptor.nnodes;n++) {
        fes[mesh[k].LGP2descriptor.nodes[n].child]=cbuffer[k][n];
        }
      delete[] cbuffer[k];
      }
    iteration=0;

    for(n=0;n<final.LGP2descriptor.nnodes;n++) {
      if(fes[n]==complex<float>(1./0.,1./0.)) {
        printf("trouble...\n");
        }
      dLGP2[n]=abs(fes[n]);
      }
    status=archiving_UGdummy2D((const char*) output, final, "a_eta_LGP2", "m", dLGP2, mask, iteration,(int) LGP2);
    for(n=0;n<final.LGP2descriptor.nnodes;n++) {
      if(fes[n]==complex<float>(1./0.,1./0.)) {
        printf("trouble...\n");
        }
      dLGP2[n]=360.0-arg(fes[n])*r2d;
      }
    status=archiving_UGdummy2D((const char*) output, final, "G_eta_LGP2", "degrees", dLGP2, mask, iteration,(int) LGP2);
    for(k=0;k<nmesh;k++) {
      delete[] fesfile[k];
      }
    delete[] output;
    }

  STDOUT_BASE_LINE("end of mesh-assembly ... \n");
  exit(0);

error:
  STDOUT_BASE_LINE("error detected, quit ... \n");
  exit(-1);
}
