
/*******************************************************************************

  T-UGOm hydrodynamic ocean model, 2006

  Unstructured Ocean Grid initiative

*******************************************************************************/

#include <stdio.h>
#include <string.h>
#include <time.h>

#include "tools-structures.h"
#include "constants.h"
#include "fe.h"
#include "mass.h"


// static const int gParallel=0;


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void fe_print(mesh_t mesh, int m)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i;
  
  printf("element %d : \n",m);
  for(i=0;i<3;i++) {
    printf("vertex %d = %d\n",i,mesh.triangles[m].vertex[i]);
    }
  for(i=0;i<3;i++) {
    printf("edge %d = %d\n",i,mesh.triangles[m].edges[i]);
    }
  
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename D,typename M> D *get_descriptor_address_template(M &mesh,int discretisation_id)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  D *descriptor=NULL;
  
  switch (discretisation_id) {
    case LGP0:
      descriptor=&mesh.LGP0descriptor;
      break;

    case LGP1:
      descriptor=&mesh.LGP1descriptor;
      break;

    case DGP1:
      descriptor=&mesh.DGP1descriptor;
      break;

    case NCP1:
      descriptor=&mesh.NCP1descriptor;
      break;

    case DNP1:
      descriptor=&mesh.DNP1descriptor;
      break;

    case LGP2:
      descriptor=&mesh.LGP2descriptor;
      break;

    case CQP0:
      descriptor=&mesh.CQP0descriptor;
      break;

    case CQP1:
      descriptor=&mesh.CQP1descriptor;
      break;

    case CQN1:
      descriptor=&mesh.CQN1descriptor;
      break;
  }

  return(descriptor);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  const discretisation_t *get_descriptor_address(const mesh_t &mesh,int discretisation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  const discretisation_t *descriptor;
  
  descriptor=get_descriptor_address_template<const discretisation_t>(mesh,discretisation);

  return(descriptor);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  discretisation_t *get_descriptor_address(mesh_t &mesh,int discretisation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  discretisation_t *descriptor;
  
  descriptor=get_descriptor_address_template<discretisation_t>(mesh,discretisation);

  return(descriptor);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  const discretisation_t & get_descriptor(const mesh_t & mesh, int discretisation_id)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  const discretisation_t *descriptor;
  
  descriptor=get_descriptor_address(mesh,discretisation_id);

  return *descriptor;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int paire_discretisation(const mesh_t& mesh, int paire, discretisation_t *z_descriptor,discretisation_t *u_descriptor)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int z_discretisation, u_discretisation;

  status=paire_discretisation_id(paire, &z_discretisation, &u_discretisation);

  *z_descriptor= get_descriptor(mesh, z_discretisation);
  *u_descriptor= get_descriptor(mesh, u_discretisation);

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  discretisation_t initLGP0(mesh_t mesh, int type)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,k,l,m,n,status;
  int n1,n2;
  int **nghbs;
  int nnghbs_max;
  int acquire;
  discretisation_t local;
  node_t *nodes;
 double t,p;

  local.nnodes=mesh.ntriangles;
  local.nelmts=mesh.ntriangles;

  local.nnpe=1;

  local.NIbE=new int*[mesh.ntriangles];

  nodes=new node_t[local.nnodes];

/* *----------------------------------------------------------------------------
  create LGP0 nodes NIbE list (for each element)*/
  for(m = 0; m < mesh.ntriangles; m++) {
    local.NIbE[m]=new int[local.nnpe];
    local.NIbE[m][0]=m;
    n=m;
    status=fe_position(mesh, mesh.triangles[m], &t, &p,0);
    nodes[n].lon=t;
    nodes[n].lat=p;
//    nodes[n].weight=mesh.triangles[m].Area;
    nodes[n].c  =cos(p*d2r);
    nodes[n].s  =sin(p*d2r);
    nodes[n].h  =-1;
/* *----------------------------------------------------------------------------
    nedded by parallel computing */
//     nodes[n].rank       =mesh.triangles[m].rank;
//     nodes[n].globalized =mesh.triangles[m].origin;
    nodes[n].nedges=3;
    nodes[n].edges=new int[nodes[n].nedges];
    for(i=0;i<nodes[n].nedges;i++) {
      nodes[n].edges[i]=mesh.triangles[m].edges[i];
      }
    nodes[n].nvtces=3;
    nodes[n].vtces=new int[nodes[n].nvtces];
    for(i=0;i<nodes[n].nvtces;i++) {
      nodes[n].vtces[i]=mesh.triangles[m].vertex[i];
      }
    }

/* *----------------------------------------------------------------------------
  create element connectivity (for each LGP0 nodes)*/
  for(n = 0; n < local.nnodes; n++) {
    nodes[n].nelmts=0;
    }
  for(m = 0; m < mesh.ntriangles; m++) {
    for(k=0;k<local.nnpe;k++) {
      n=local.NIbE[m][k];
      nodes[n].nelmts++;
      }
    }
  for(n = 0; n < local.nnodes; n++) {
    nodes[n].elmts=new int[nodes[n].nelmts];
    nodes[n].nelmts=0;
    }
  for(m = 0; m < mesh.ntriangles; m++) {
    for(k=0;k<local.nnpe;k++) {
      n=local.NIbE[m][k];
      nodes[n].elmts[nodes[n].nelmts]=m;
      nodes[n].nelmts++;
      }
    }

/* *----------------------------------------------------------------------------
  create neighbours connectivity (for each LGP0 nodes)*/
  nnghbs_max=mesh.nnghm+1;
  nghbs=new int*[local.nnodes];
  for(n = 0; n < local.nnodes; n++) {
    nodes[n].nnghbs=0;
    nghbs[n]=new int[nnghbs_max];
    for(i=0;i<nnghbs_max;i++) nghbs[n][i]=-1;
    }
  for(m = 0; m < mesh.ntriangles; m++) {
    for(k=0;k<local.nnpe;k++) {
      n1=local.NIbE[m][k];
      for(l=0;l<local.nnpe;l++) {
        if(k==l) continue;
        n2=local.NIbE[m][l];
        acquire=1;
        for(i=0;i<nodes[n1].nnghbs;i++) {
          if(nghbs[n1][i]==n2) {
            acquire=0;
            break;
            }
          }
        if(acquire==1) {
          nghbs[n1][i]=n2;
          nodes[n1].nnghbs++;
          if(nodes[n1].nnghbs>=nnghbs_max) {
            TRAP_ERR_EXIT(-1, "LGP0 discretisation, neigbours overflow");
            }
          }
        }
      }
    }
  for(n = 0; n < local.nnodes; n++) {
    nodes[n].nghbs=new int[nodes[n].nnghbs];
    for(i=0;i<nodes[n].nnghbs;i++) nodes[n].nghbs[i]=nghbs[n][i];
    }

  for(n = 0; n < local.nnodes; n++) {
    delete[] nghbs[n];
    }
  delete[] nghbs;

  local.nodes=nodes;
  local.type=LGP0;
  return(local);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  discretisation_t initLGP1(mesh_t mesh, int type)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,k,l,m,n;
  int n1,n2,nn;
  int **nghbs;
  int nnghbs_max;
  int acquire;
  discretisation_t local;
  node_t *nodes;

  local.nnodes=mesh.nvtxs;
  local.nelmts=mesh.ntriangles;

  local.nnpe=3;

  local.NIbE=new int*[mesh.ntriangles];

  nodes=new node_t[local.nnodes];

/* *----------------------------------------------------------------------------
  create LGP1 nodes NIbE list (for each element)*/
  for(m = 0; m < mesh.ntriangles; m++) {
    local.NIbE[m]=new int[local.nnpe];
    for(k=0;k<local.nnpe;k++) local.NIbE[m][k]=-1;
    for(k=0;k<3;k++) {
      n=mesh.triangles[m].vertex[k];
      local.NIbE[m][k]=n;
//      nodes[n].weight+=mesh.triangles[m].Area/3.0;
      }
    }
  for(n = 0; n < mesh.nvtxs; n++) {
    nn=n;
    nodes[n].lon=mesh.vertices[nn].lon;
    nodes[n].lat=mesh.vertices[nn].lat;
    nodes[n].c  =mesh.vertices[nn].c;
    nodes[n].s  =mesh.vertices[nn].s;
    nodes[n].h  =mesh.vertices[nn].h;
      if(mesh.vertices[nn].nedges>0)
        nodes[n].edges=new int[mesh.vertices[nn].nedges];
      else
        nodes[n].edges=NULL;
      nodes[n].nedges=mesh.vertices[nn].nedges;
      for(i=0;i<nodes[n].nedges;i++) {
        nodes[n].edges[i]=mesh.vertices[nn].edges[i];
        }
      nodes[n].nvtces=mesh.vertices[nn].nngh+1;
      nodes[n].vtces=new int[nodes[n].nvtces];
      nodes[n].vtces[0]=nn;
      for(i=0;i<nodes[n].nvtces-1;i++) {
        nodes[n].vtces[i+1]=mesh.vertices[nn].ngh[i];
        }
    }

/* *----------------------------------------------------------------------------
  create element connectivity (for each LGP1 nodes)*/
  for(n = 0; n < local.nnodes; n++) {
    nodes[n].nelmts=0;
    }
  for(m = 0; m < mesh.ntriangles; m++) {
    for(k=0;k<local.nnpe;k++) {
      n=local.NIbE[m][k];
      nodes[n].nelmts++;
      }
    }
  for(n = 0; n < local.nnodes; n++) {
    nodes[n].elmts=new int[nodes[n].nelmts];
    nodes[n].nelmts=0;
    }
  for(m = 0; m < mesh.ntriangles; m++) {
    for(k=0;k<local.nnpe;k++) {
      n=local.NIbE[m][k];
      nodes[n].elmts[nodes[n].nelmts]=m;
      nodes[n].nelmts++;
      }
    }

/* *----------------------------------------------------------------------------
  create neighbours connectivity (for each LGP1 nodes)*/
  nnghbs_max=mesh.nnghm;
  nghbs=new int*[local.nnodes];
  for(n = 0; n < local.nnodes; n++) {
    nodes[n].nnghbs=0;
    nghbs[n]=new int[nnghbs_max];
    for(i=0;i<nnghbs_max;i++) nghbs[n][i]=-1;
    }
  for(m = 0; m < mesh.ntriangles; m++) {
    for(k=0;k<local.nnpe;k++) {
      n1=local.NIbE[m][k];
      for(l=0;l<local.nnpe;l++) {
        if(k==l) continue;
        n2=local.NIbE[m][l];
        acquire=1;
        for(i=0;i<nodes[n1].nnghbs;i++) {
          if(nghbs[n1][i]==n2) {
            acquire=0;
            break;
            }
          }
        if(acquire==1) {
          nghbs[n1][i]=n2;
          nodes[n1].nnghbs++;
          if(nodes[n1].nnghbs>nnghbs_max) {
            TRAP_ERR_EXIT(-1, "LGP1 discretisation, neigbours overflow");
            }
          }
        }
      }
    }
  for(n = 0; n < local.nnodes; n++) {
    nodes[n].nghbs=new int[nodes[n].nnghbs];
    for(i=0;i<nodes[n].nnghbs;i++) nodes[n].nghbs[i]=nghbs[n][i];
    }

  for(n = 0; n < local.nnodes; n++) {
    delete[] nghbs[n];
    }
  delete[] nghbs;

  local.nodes=nodes;
  local.type=LGP1;
  return(local);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  discretisation_t initDGP1(mesh_t mesh, int type)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,k,l,m,n;
  int n1,n2,nn;
  int **nghbs;
  int nnghbs_max;
  int acquire;
  discretisation_t local;
  node_t *nodes;
  int side[3][2]={{1,2},{2,0},{0,1}};
  int node[3][3]={{2,0,1},{0,1,2},{1,2,0}};

  local.nnodes=3*mesh.ntriangles;
  local.nelmts=mesh.ntriangles;

  local.nnpe=3;

  local.NIbE=new int*[mesh.ntriangles];

  nodes=new node_t[local.nnodes];

/* *----------------------------------------------------------------------------
  create DGP1 nodes NIbE list (for each element)*/
  n=0;
  for(m = 0; m < mesh.ntriangles; m++) {
    local.NIbE[m]=new int[local.nnpe];
    for(k=0;k<local.nnpe;k++) local.NIbE[m][k]=-1;
    for(k=0;k<3;k++) {
      nn=mesh.triangles[m].vertex[k];
      local.NIbE[m][k]=n;
      nodes[n].lon=mesh.vertices[nn].lon;
      nodes[n].lat=mesh.vertices[nn].lat;
//      nodes[n].weight=mesh.triangles[m].Area/3.0;
      nodes[n].c  =mesh.vertices[nn].c;
      nodes[n].s  =mesh.vertices[nn].s;
      nodes[n].h  =mesh.vertices[nn].h;
/* *----------------------------------------------------------------------------
      nedded by parallel computing (to be checked) */
//       nodes[n].rank       =mesh.triangles[m].rank;
//       nodes[n].globalized =mesh.triangles[m].origin*3+k;
      nodes[n].nedges=2;
      nodes[n].edges=new int[nodes[n].nedges];
      for(i=0;i<nodes[n].nedges;i++) {
        nodes[n].edges[i]=mesh.triangles[m].edges[side[k][i]];
        }
      nodes[n].nvtces=3;
      nodes[n].vtces=new int[nodes[n].nvtces];
      for(i=0;i<nodes[n].nvtces;i++) {
        nodes[n].vtces[i]=mesh.triangles[m].vertex[node[k][i]];
        }
      n++;
      }
    }

/* *----------------------------------------------------------------------------
  create element connectivity (for each DGP1 nodes)*/
  for(n = 0; n < local.nnodes; n++) {
    nodes[n].nelmts=0;
    }
  for(m = 0; m < mesh.ntriangles; m++) {
    for(k=0;k<local.nnpe;k++) {
      n=local.NIbE[m][k];
      nodes[n].nelmts++;
      }
    }
  for(n = 0; n < local.nnodes; n++) {
    nodes[n].elmts=new int[nodes[n].nelmts];
    nodes[n].nelmts=0;
    }
  for(m = 0; m < mesh.ntriangles; m++) {
    for(k=0;k<local.nnpe;k++) {
      n=local.NIbE[m][k];
      nodes[n].elmts[nodes[n].nelmts]=m;
      nodes[n].nelmts++;
      }
    }

/* *----------------------------------------------------------------------------
  create neighbours connectivity (for each DGP1 nodes)*/
  nnghbs_max=mesh.nnghm*3;
  nghbs=new int*[local.nnodes];
  for(n = 0; n < local.nnodes; n++) {
    nodes[n].nnghbs=0;
    nghbs[n]=new int[nnghbs_max];
    for(i=0;i<nnghbs_max;i++) nghbs[n][i]=-1;
    }
  for(m = 0; m < mesh.ntriangles; m++) {
    for(k=0;k<local.nnpe;k++) {
      n1=local.NIbE[m][k];
      for(l=0;l<local.nnpe;l++) {
        if(k==l) continue;
        n2=local.NIbE[m][l];
        acquire=1;
        for(i=0;i<nodes[n1].nnghbs;i++) {
          if(nghbs[n1][i]==n2) {
            acquire=0;
            break;
            }
          }
        if(acquire==1) {
          nghbs[n1][i]=n2;
          nodes[n1].nnghbs++;
          if(nodes[n1].nnghbs>=nnghbs_max) {
            TRAP_ERR_EXIT(-1, "DGP1 discretisation, neigbours overflow");
            }
          }
        }
      }
    }
  for(n = 0; n < local.nnodes; n++) {
    nodes[n].nghbs=new int[nodes[n].nnghbs];
    for(i=0;i<nodes[n].nnghbs;i++) nodes[n].nghbs[i]=nghbs[n][i];
    }

  delete[] nghbs;

  local.nodes=nodes;
  local.type=DGP1;
  return(local);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  discretisation_t initDNP1(mesh_t mesh, int type)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,k,l,m,n;
  int n1,n2,nn;
  int **nghbs;
  int nnghbs_max;
  int acquire;
  discretisation_t local;
  node_t *nodes;
  int node[3][2]={{1,2},{2,0},{0,1}};

  local.nnodes=3*mesh.ntriangles;
  local.nelmts=mesh.ntriangles;

  local.nnpe=3;

  local.NIbE=new int*[mesh.ntriangles];

  nodes=new node_t[local.nnodes];

/* *----------------------------------------------------------------------------
  create DGP1 nodes NIbE list (for each element)*/
  n=0;
  for(m = 0; m < mesh.ntriangles; m++) {
    local.NIbE[m]=new int[local.nnpe];
    for(k=0;k<local.nnpe;k++) local.NIbE[m][k]=-1;
    for(k=0;k<3;k++) {
      nn=mesh.triangles[m].edges[k];
      local.NIbE[m][k]=n;
      nodes[n].lon=mesh.edges[nn].lon;
      nodes[n].lat=mesh.edges[nn].lat;
//      nodes[n].weight=mesh.triangles[m].Area/3.0;
      nodes[n].c  =mesh.edges[nn].c;
      nodes[n].s  =mesh.edges[nn].s;
//      nodes[n].h  =mesh.edges[nn].h;
/* *----------------------------------------------------------------------------
      nedded by parallel computing (to be checked) */
//       nodes[n].rank       =mesh.triangles[m].rank;
//       nodes[n].globalized =mesh.triangles[m].origin*3+k;
      nodes[n].nedges=1;
      nodes[n].edges=new int[nodes[n].nedges];
      nodes[n].edges[0]=nn;
      nodes[n].nvtces=2;
      nodes[n].vtces=new int[nodes[n].nvtces];
      for(i=0;i<nodes[n].nvtces;i++) {
        nodes[n].vtces[i]=mesh.triangles[m].vertex[node[k][i]];
        }
      n++;
      }
    }

/* *----------------------------------------------------------------------------
  create element connectivity (for each DGP1 nodes)*/
  for(n = 0; n < local.nnodes; n++) {
    nodes[n].nelmts=0;
    }
  for(m = 0; m < mesh.ntriangles; m++) {
    for(k=0;k<local.nnpe;k++) {
      n=local.NIbE[m][k];
      nodes[n].nelmts++;
      }
    }
  for(n = 0; n < local.nnodes; n++) {
    nodes[n].elmts=new int[nodes[n].nelmts];
    nodes[n].nelmts=0;
    }
  for(m = 0; m < mesh.ntriangles; m++) {
    for(k=0;k<local.nnpe;k++) {
      n=local.NIbE[m][k];
      nodes[n].elmts[nodes[n].nelmts]=m;
      nodes[n].nelmts++;
      }
    }

/* *----------------------------------------------------------------------------
  create neighbours connectivity (for each DGP1 nodes)*/
  nnghbs_max=mesh.nnghm*3;
  nghbs=new int*[local.nnodes];
  for(n = 0; n < local.nnodes; n++) {
    nodes[n].nnghbs=0;
    nghbs[n]=new int[nnghbs_max];
    for(i=0;i<nnghbs_max;i++) nghbs[n][i]=-1;
    }
  for(m = 0; m < mesh.ntriangles; m++) {
    for(k=0;k<local.nnpe;k++) {
      n1=local.NIbE[m][k];
      for(l=0;l<local.nnpe;l++) {
        if(k==l) continue;
        n2=local.NIbE[m][l];
        acquire=1;
        for(i=0;i<nodes[n1].nnghbs;i++) {
          if(nghbs[n1][i]==n2) {
            acquire=0;
            break;
            }
          }
        if(acquire==1) {
          nghbs[n1][i]=n2;
          nodes[n1].nnghbs++;
          if(nodes[n1].nnghbs>=nnghbs_max) {
            TRAP_ERR_EXIT(-1, "DGP1 discretisation, neigbours overflow");
            }
          }
        }
      }
    }
  for(n = 0; n < local.nnodes; n++) {
    nodes[n].nghbs=new int[nodes[n].nnghbs];
    for(i=0;i<nodes[n].nnghbs;i++) nodes[n].nghbs[i]=nghbs[n][i];
    }

  for(size_t n = 0; n < local.nnodes; n++){
    delete [] nghbs[n];
    }
  delete[] nghbs;

  local.nodes=nodes;
  local.type=DNP1;
  return(local);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  discretisation_t NIbE_DNP1(mesh_t mesh, int type)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,m,n;
  int nn;
  discretisation_t local;

  local.nnodes=3*mesh.ntriangles;
  local.nelmts=mesh.ntriangles;

  local.nnpe=3;

  local.NIbE=new int*[mesh.ntriangles];

/* *----------------------------------------------------------------------------
  create DGP1 nodes NIbE list (for each element)*/
  n=0;
  for(m = 0; m < mesh.ntriangles; m++) {
    local.NIbE[m]=new int[local.nnpe];
    for(k=0;k<local.nnpe;k++) local.NIbE[m][k]=-1;
    for(k=0;k<3;k++) {
      nn=mesh.triangles[m].edges[k];
      local.NIbE[m][k]=n;
      n++;
      }
    }

  local.type=DNP1;
  return(local);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  discretisation_t initNCP1(mesh_t mesh, int type)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Development notes:

  Check :

  Note:

  This approach is somehow a waste of memory, as the mesh structure
  contains already all the necessary informations.

  It is used to standardize the multi-discretisation capabilities in
  T-UGOm.

  Mignt be revised later on.

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
{
  int i,k,l,m,n;
  int n1,n2,nn;
  int **nghbs;
  int nnghbs_max;
  int acquire;
  discretisation_t local;
  node_t *nodes;

  local.nnodes=mesh.nedges;
  local.nelmts=mesh.ntriangles;

  local.nnpe=3;

  local.NIbE=new int*[mesh.ntriangles];

  nodes=new node_t[local.nnodes];

/* *----------------------------------------------------------------------------
  create LGP1 nodes NIbE list (for each element)*/
  for(m = 0; m < mesh.ntriangles; m++) {
    local.NIbE[m]=new int[local.nnpe];
    for(k=0;k<local.nnpe;k++) local.NIbE[m][k]=-1;
    for(k=0;k<3;k++) {
      n=mesh.triangles[m].edges[k];
      local.NIbE[m][k]=n;
      }
    }
  for(n = 0; n < local.nnodes; n++) {
    nn=n;
    nodes[n].lon=mesh.edges[nn].lon;
    nodes[n].lat=mesh.edges[nn].lat;
//    nodes[n].weight=mesh.edges[nn].mw;
    nodes[n].c  =mesh.edges[nn].c;
    nodes[n].s  =mesh.edges[nn].s;
/* *----------------------------------------------------------------------------
    nedded by parallel computing */
//     nodes[n].rank       =mesh.edges[nn].rank;
//     nodes[n].globalized =mesh.edges[nn].origin;
//    nodes[n].h  =mesh.edges[nn].h;
//     nodes[n].nvtces=mesh.edges[nn].nvtces;
//     nodes[n].edges=new int[nodes[n].nvtces];
//     for(i=0;i<nodes[n].nvtces;i++) {
//       nodes[n].vtces[i+1]=mesh.edges[nn].vtces[i];
//       }
    nodes[n].nedges=mesh.edges[nn].nngh+1;
    nodes[n].edges=new int[nodes[n].nedges];
    nodes[n].edges[0]=nn;
    for(i=0;i<nodes[n].nvtces-1;i++) {
      nodes[n].edges[i+1]=mesh.edges[nn].ngh[i];
      }
    }

/* *----------------------------------------------------------------------------
  create element connectivity (for each LGP1 nodes)*/
  for(n = 0; n < local.nnodes; n++) {
    nodes[n].nelmts=0;
    }
  for(m = 0; m < mesh.ntriangles; m++) {
    for(k=0;k<local.nnpe;k++) {
      n=local.NIbE[m][k];
      nodes[n].nelmts++;
      }
    }
  for(n = 0; n < local.nnodes; n++) {
    nodes[n].elmts=new int[nodes[n].nelmts];
    nodes[n].nelmts=0;
    }
  for(m = 0; m < mesh.ntriangles; m++) {
    for(k=0;k<local.nnpe;k++) {
      n=local.NIbE[m][k];
      nodes[n].elmts[nodes[n].nelmts]=m;
      nodes[n].nelmts++;
      }
    }

/* *----------------------------------------------------------------------------
  create neighbours connectivity (for each LGP1 nodes)*/
  nnghbs_max=mesh.nnghm;
  nghbs=new int*[local.nnodes];
  for(n = 0; n < local.nnodes; n++) {
    nodes[n].nnghbs=0;
    nghbs[n]=new int[nnghbs_max];
    for(i=0;i<nnghbs_max;i++) nghbs[n][i]=-1;
    }
  for(m = 0; m < mesh.ntriangles; m++) {
    for(k=0;k<local.nnpe;k++) {
      n1=local.NIbE[m][k];
      for(l=0;l<local.nnpe;l++) {
        if(k==l) continue;
        n2=local.NIbE[m][l];
        acquire=1;
        for(i=0;i<nodes[n1].nnghbs;i++) {
          if(nghbs[n1][i]==n2) {
            acquire=0;
            break;
            }
          }
        if(acquire==1) {
          nghbs[n1][i]=n2;
          nodes[n1].nnghbs++;
          if(nodes[n1].nnghbs>nnghbs_max) {
            TRAP_ERR_EXIT(-1, "LGP1 discretisation, neigbours overflow");
            }
          }
        }
      }
    }
  for(n = 0; n < local.nnodes; n++) {
    nodes[n].nghbs=new int[nodes[n].nnghbs];
    for(i=0;i<nodes[n].nnghbs;i++) nodes[n].nghbs[i]=nghbs[n][i];
    }

  for(n = 0; n < local.nnodes; n++) {
    delete[] nghbs[n];
    }
  delete[] nghbs;

  local.nodes=nodes;
  local.type=NCP1;
  return(local);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  discretisation_t initLGP2(mesh_t mesh, int type)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,k,kk,l,m,n;
  int n1,n2,nn;
  int *epassed,*vpassed;
  int **nghbs;
  int nnghbs_max;
  discretisation_t local;
  node_t *nodes;
  
  m=mesh.ntriangles-1;
  if(mesh.triangles[m].vertex==0) {
    TRAP_ERR_EXIT(-1, "vertex tables not initialized ()");
    }
  if(mesh.triangles[m].edges==0 or mesh.nedges==0) {
    TRAP_ERR_EXIT(-1, "edge tables not initialized: call fe_edgetable()\n");
    }
  local.nnodes=mesh.nvtxs+mesh.nedges;
  local.nelmts=mesh.ntriangles;

  local.nnpe=6;

  local.NIbE=new int*[mesh.ntriangles];

  epassed=new int[mesh.nedges];
  vpassed=new int[mesh.nvtxs];

  for(n = 0; n < mesh.nedges; n++) epassed[n]=-1;
  for(n = 0; n < mesh.nvtxs; n++)  vpassed[n]=-1;

  nodes=new node_t[local.nnodes];

/* *----------------------------------------------------------------------------
  create LGP2 nodes NIbE list (for each element)*/
  n=0;
  for(m = 0; m < mesh.ntriangles; m++) {
    local.NIbE[m]=new int[local.nnpe];
    for(k=0;k<local.nnpe;k++) local.NIbE[m][k]=-1;
    for(k=0;k<3;k++) {
      nn=mesh.triangles[m].vertex[k];
      if(vpassed[nn]==-1) {
        local.NIbE[m][2*k]=n;
        nodes[n].lon=mesh.vertices[nn].lon;
        nodes[n].lat=mesh.vertices[nn].lat;
        nodes[n].c  =mesh.vertices[nn].c;
        nodes[n].h  =mesh.vertices[nn].h;
/* *----------------------------------------------------------------------------
        nedded by parallel computing (to be checked) */
//         nodes[n].rank       =mesh.vertices[nn].rank;
// 	nodes[n].globalized =n;
        if(mesh.vertices[nn].nedges>0) {
          nodes[n].edges=new int[mesh.vertices[nn].nedges];
          nodes[n].nedges=mesh.vertices[nn].nedges;
          for(i=0;i<nodes[n].nedges;i++) {
           nodes[n].edges[i]=mesh.vertices[nn].edges[i];
            }
          }
        nodes[n].nvtces=mesh.vertices[nn].nngh+1;
        nodes[n].vtces=new int[nodes[n].nvtces];
        nodes[n].vtces[0]=nn;
        for(i=0;i<nodes[n].nvtces-1;i++) {
          nodes[n].vtces[i+1]=mesh.vertices[nn].ngh[i];
          }
        vpassed[nn]=n;
        n++;
        }
      else {
        local.NIbE[m][2*k]=vpassed[nn];
        }
      nn=mesh.triangles[m].edges[(k+2) % 3];
      kk=(2*k+1) % 6;
      if(epassed[nn]==-1) {
        local.NIbE[m][kk]=n;
//         nodes[n].lon=mesh.edges[nn].lon*r2d;
//         nodes[n].lat=mesh.edges[nn].lat*r2d;
        nodes[n].lon=mesh.edges[nn].lon;
        nodes[n].lat=mesh.edges[nn].lat;
        nodes[n].c  =mesh.edges[nn].c;
        n1=mesh.edges[nn].extremity[0];
        n2=mesh.edges[nn].extremity[1];
        nodes[n].h  =0.5*(mesh.vertices[n1].h+mesh.vertices[n2].h);
/* *----------------------------------------------------------------------------
        nedded by parallel computing (to be checked) */
//         nodes[n].rank       =mesh.edges[nn].rank;
//         nodes[n].globalized =n;

        nodes[n].nedges=mesh.edges[nn].nngh+1;
        nodes[n].edges=new int[nodes[n].nedges];
        nodes[n].edges[0]=nn;
        for(i=0;i<nodes[n].nedges-1;i++) {
          nodes[n].edges[i+1]=mesh.edges[nn].ngh[i];
          }
        if(mesh.edges[nn].nvertex>0) {
          nodes[n].nvtces=mesh.edges[nn].nvertex;
          nodes[n].vtces=new int[nodes[n].nvtces];
          for(i=0;i<nodes[n].nvtces;i++) {
            nodes[n].vtces[i]=mesh.edges[nn].vertex[i];
            }
          }
        epassed[nn]=n;
        n++;
        }
      else {
        local.NIbE[m][kk]=epassed[nn];
        }
      }
    }

  delete[] vpassed;
  delete[] epassed;
  
/* *----------------------------------------------------------------------------
  update globalized numbering */
//   if(gParallel!=0) {
//     for(m = 0; m < mesh.ntriangles; m++) {
//       mm=mesh.triangles[m].origin;
//       for(k=0;k<local.nnpe;k++) {
//         n  = local.NIbE[m][k];
//         nn = mesh.origin->LGP2descriptor.NIbE[mm][k];
//         nodes[n].globalized=nn;
//         }
//       }
/* *----------------------------------------------------------------------------
      do not destroy, needed by output routines*/
//    mesh.origin->LGP2descriptor.destroy();
//    }
    
/* *----------------------------------------------------------------------------
  create element connectivity (for each P2 nodes)*/
  for(n = 0; n < local.nnodes; n++) {
    nodes[n].nelmts=0;
    }
  for(m = 0; m < mesh.ntriangles; m++) {
    for(k=0;k<local.nnpe;k++) {
      n=local.NIbE[m][k];
      nodes[n].nelmts++;
      }
    }
  for(n = 0; n < local.nnodes; n++) {
    nodes[n].elmts=new int[nodes[n].nelmts];
    nodes[n].nelmts=0;
    }
  for(m = 0; m < mesh.ntriangles; m++) {
    for(k=0;k<local.nnpe;k++) {
      n=local.NIbE[m][k];
      nodes[n].elmts[nodes[n].nelmts]=m;
      nodes[n].nelmts++;
      }
    }

/* *----------------------------------------------------------------------------
  create neighbours connectivity (for each LGP2 nodes)*/
  nnghbs_max=mesh.nnghm*3+2;
  nghbs=new int*[local.nnodes];
  for(n = 0; n < local.nnodes; n++) {
    nodes[n].nnghbs=0;
    nghbs[n]=new int[nnghbs_max];
    for(i=0;i<nnghbs_max;i++) nghbs[n][i]=-1;
    }
  
  #pragma omp parallel for private(i,k,l,n1,n2)
  for(m = 0; m < mesh.ntriangles; m++) {
    int *NIbEm=local.NIbE[m];
    for(k=0;k<local.nnpe;k++) {
      n1=NIbEm[k];
      int *nodesn1nnghbs=&nodes[n1].nnghbs;
      for(l=0;l<local.nnpe;l++) {
        if(k==l) continue;
        n2=NIbEm[l];
        int acquire=1;
        for(i=0;i<*nodesn1nnghbs;i++) {
          if(nghbs[n1][i]==n2) {
            acquire=0;
            break;
            }
          }
        if(acquire==1) {
          if(*nodesn1nnghbs>=nnghbs_max) {
            printf("node %d (%d in element %d) : \n",n1,k,m);
            fe_print(mesh, m);
            TRAP_ERR_EXIT(-1, "LGP2 discretisation, neigbours overflow");
            }
          nghbs[n1][i]=n2;
          (*nodesn1nnghbs)++;
          }
        }
      }
    }
  for(n = 0; n < local.nnodes; n++) {
    nodes[n].nghbs=new int[nodes[n].nnghbs];
    for(i=0;i<nodes[n].nnghbs;i++) nodes[n].nghbs[i]=nghbs[n][i];
    }

  for(size_t n = 0; n < local.nnodes; n++){
    delete [] nghbs[n];
    }
  delete[] nghbs;

  local.nodes=nodes;
  local.type=LGP2;
  return(local);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  discretisation_t NIbE_LGP2(mesh_t mesh, int type)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,kk,m,n;
  int nn;
  int *epassed,*vpassed;
  discretisation_t local;

  local.nnodes=mesh.nvtxs+mesh.nedges;
  local.nelmts=mesh.ntriangles;

  local.nnpe=6;

  local.NIbE=new int*[mesh.ntriangles];

  epassed=new int[mesh.nedges];
  vpassed=new int[mesh.nvtxs];

  for(n = 0; n < mesh.nedges; n++) epassed[n]=-1;
  for(n = 0; n < mesh.nvtxs; n++)  vpassed[n]=-1;

/* *----------------------------------------------------------------------------
  create P2 nodes NIbE list (for each element)*/
  n=0;
  for(m = 0; m < mesh.ntriangles; m++) {
    local.NIbE[m]=new int[local.nnpe];
    for(k=0;k<local.nnpe;k++) local.NIbE[m][k]=-1;
    for(k=0;k<3;k++) {
      nn=mesh.triangles[m].vertex[k];
      if(vpassed[nn]==-1) {
        local.NIbE[m][2*k]=n;
        vpassed[nn]=n;
        n++;
        }
      else {
        local.NIbE[m][2*k]=vpassed[nn];
        }
      nn=mesh.triangles[m].edges[(k+2) % 3];
      kk=(2*k+1) % 6;
      if(epassed[nn]==-1) {
        local.NIbE[m][kk]=n;
        epassed[nn]=n;
        n++;
        }
      else {
        local.NIbE[m][kk]=epassed[nn];
        }
      }
    }
  delete[] vpassed;
  delete[] epassed;
  
  local.type=LGP2;
  return(local);
}

// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
//
//   discretisation_t initDGP2_unfinished(mesh_t mesh, int type)
//
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   int i,j,k,kk,l,m,n,status;
//   int mm,n1,n2,nn;
//   int *epassed,*vpassed;
//   int **nghbs;
//   int nnghbs_max;
//   int acquire;
//   discretisation_t local;
//   node_t *nodes;
//
//   local.nnodes=6*mesh.ntriangles;
//   local.nelmts=mesh.ntriangles;
//
//   local.nnpe=6;
//
//   local.NIbE=new int*[mesh.ntriangles];
//
//   epassed=new int[mesh.nedges];
//   vpassed=new int[mesh.nvtxs];
//
//   for(n = 0; n < mesh.nedges; n++) epassed[n]=-1;
//   for(n = 0; n < mesh.nvtxs; n++)  vpassed[n]=-1;
//
//   nodes=new node_t[local.nnodes];
//
// /* *----------------------------------------------------------------------------
//   create P2 nodes NIbE list (for each element)*/
//   n=0;
//   for(m = 0; m < mesh.ntriangles; m++) {
//     local.NIbE[m]=new int[local.nnpe];
//     for(k=0;k<local.nnpe;k++) local.NIbE[m][k]=-1;
//     for(k=0;k<3;k++) {
//       nn=mesh.triangles[m].vertex[k];
//       local.NIbE[m][2*k]=n;
//       nodes[n].lon=mesh.vertices[nn].lon;
//       nodes[n].lat=mesh.vertices[nn].lat;
//       nodes[n].c  =mesh.vertices[nn].c;
//       nodes[n].h  =mesh.vertices[nn].h;
//       n++;
//       }
//     }
//
//
// /* *----------------------------------------------------------------------------
//   create element connectivity (for each P2 nodes)*/
//   for(n = 0; n < local.nnodes; n++) {
//     nodes[n].nelmts=0;
//     }
//   for(m = 0; m < mesh.ntriangles; m++) {
//     for(k=0;k<local.nnpe;k++) {
//       n=local.NIbE[m][k];
//       nodes[n].nelmts++;
//       }
//     }
//   for(n = 0; n < local.nnodes; n++) {
//     nodes[n].elmts=new int[nodes[n].nelmts];
//     nodes[n].nelmts=0;
//     }
//   for(m = 0; m < mesh.ntriangles; m++) {
//     for(k=0;k<local.nnpe;k++) {
//       n=local.NIbE[m][k];
//       nodes[n].elmts[nodes[n].nelmts]=m;
//       nodes[n].nelmts++;
//       }
//     }
//
//   local.nodes=nodes;
//   local.type=DGP2;
//   return(local);
// }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  discretisation_t initDGP2(mesh_t mesh, int type)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,k,l,m,n;
  int n1,n2,nn;
  int *epassed,*vpassed;
  int **nghbs;
  int nnghbs_max;
  int acquire;
  discretisation_t local;
  node_t *nodes;

  local.nnodes=6*mesh.ntriangles;
  local.nelmts=mesh.ntriangles;

  local.nnpe=6;

  local.NIbE=new int*[mesh.ntriangles];

  epassed=new int[mesh.nedges];
  vpassed=new int[mesh.nvtxs];

  for(n = 0; n < mesh.nedges; n++) epassed[n]=-1;
  for(n = 0; n < mesh.nvtxs; n++)  vpassed[n]=-1;

  nodes=new node_t[local.nnodes];

/* *----------------------------------------------------------------------------
  create P2 nodes NIbE list (for each element)*/
  for(m = 0; m < mesh.ntriangles; m++) {
    local.NIbE[m]=new int[local.nnpe];
    for(k=0;k<local.nnpe;k++) local.NIbE[m][k]=-1;
    for(k=0;k<3;k++) {
      nn=mesh.triangles[m].vertex[k];
      n=m*6+2*k;
      local.NIbE[m][2*k]=n;
      nodes[n].lon=mesh.vertices[nn].lon;
      nodes[n].lat=mesh.vertices[nn].lat;
      nodes[n].c  =mesh.vertices[nn].c;
      nodes[n].s  =mesh.vertices[nn].s;
      nodes[n].h  =mesh.vertices[nn].h;
      }
    for(k=0;k<3;k++) {
      nn=mesh.triangles[m].edges[k];
      n=m*6+2*k+1;
      local.NIbE[m][2*k+1]=n;
      nodes[n].lon=mesh.edges[nn].lon;
      nodes[n].lat=mesh.edges[nn].lat;
      nodes[n].c  =mesh.edges[nn].c;
      nodes[n].s  =mesh.edges[nn].s;
//      nodes[n].h  =mesh.edges[nn].h;
      }
    }

  
/* *----------------------------------------------------------------------------
  create element connectivity (for each DGP2 nodes)*/
  for(n = 0; n < local.nnodes; n++) {
    nodes[n].nelmts=0;
    }
  for(m = 0; m < mesh.ntriangles; m++) {
    for(k=0;k<local.nnpe;k++) {
      n=local.NIbE[m][k];
      nodes[n].nelmts++;
      }
    }
  for(n = 0; n < local.nnodes; n++) {
    nodes[n].elmts=new int[nodes[n].nelmts];
    nodes[n].nelmts=0;
    }
  for(m = 0; m < mesh.ntriangles; m++) {
    for(k=0;k<local.nnpe;k++) {
      n=local.NIbE[m][k];
      nodes[n].elmts[nodes[n].nelmts]=m;
      nodes[n].nelmts++;
      }
    }
/* *----------------------------------------------------------------------------
  create neighbours connectivity (for each DGP2 nodes)*/
  nnghbs_max=6;
  nghbs=new int*[local.nnodes];
  for(n = 0; n < local.nnodes; n++) {
    nodes[n].nnghbs=0;
    nghbs[n]=new int[nnghbs_max];
    for(i=0;i<nnghbs_max;i++) nghbs[n][i]=-1;
    }
  for(m = 0; m < mesh.ntriangles; m++) {
    for(k=0;k<local.nnpe;k++) {
      n1=local.NIbE[m][k];
      for(l=0;l<local.nnpe;l++) {
        if(k==l) continue;
        n2=local.NIbE[m][l];
        acquire=1;
        for(i=0;i<nodes[n1].nnghbs;i++) {
          if(nghbs[n1][i]==n2) {
            acquire=0;
            break;
            }
          }
        if(acquire==1) {
          nghbs[n1][i]=n2;
          nodes[n1].nnghbs++;
          if(nodes[n1].nnghbs>=nnghbs_max) {
            TRAP_ERR_EXIT(-1, "DGP2 discretisation, neigbours overflow");
            }
          }
        }
      }
    }
  for(n = 0; n < local.nnodes; n++) {
    nodes[n].nghbs=new int[nodes[n].nnghbs];
    for(i=0;i<nodes[n].nnghbs;i++) nodes[n].nghbs[i]=nghbs[n][i];
    }

  delete[] nghbs;

  local.nodes=nodes;
  local.type=DGP2;
  return(local);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  discretisation_t initCQP0(mesh_t mesh, int type)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,k,l,m,n,status;
  int n1,n2;
  int **nghbs;
  int nnghbs_max;
  int acquire;
  discretisation_t local;
  node_t *nodes;
 double t,p;

  local.nnodes=mesh.nquadrangles;
  local.nelmts=mesh.nquadrangles;

  local.nnpe=1;

  local.NIbE=new int*[mesh.nquadrangles];

  nodes=new node_t[local.nnodes];

/* *----------------------------------------------------------------------------
  create LGP0 nodes NIbE list (for each element)*/
  for(m = 0; m < mesh.nquadrangles; m++) {
    local.NIbE[m]=new int[local.nnpe];
    local.NIbE[m][0]=m;
    n=m;
    status=fe_position(mesh, mesh.quadrangles[m], &t, &p,0);
    nodes[n].lon=t;
    nodes[n].lat=p;
    nodes[n].c  =cos(p*d2r);
    nodes[n].s  =sin(p*d2r);
    nodes[n].h  =-1;
/* *----------------------------------------------------------------------------
    nedded by parallel computing */
//     nodes[n].rank       =mesh.quadrangles[m].rank;
//     nodes[n].globalized =mesh.quadrangles[m].origin;
    nodes[n].nedges=4;
    nodes[n].edges=new int[nodes[n].nedges];
    for(i=0;i<nodes[n].nedges;i++) {
      nodes[n].edges[i]=mesh.quadrangles[m].edges[i];
      }
    nodes[n].nvtces=4;
    nodes[n].vtces=new int[nodes[n].nvtces];
    for(i=0;i<nodes[n].nvtces;i++) {
      nodes[n].vtces[i]=mesh.quadrangles[m].vertex[i];
      }
    }

/* *----------------------------------------------------------------------------
  create element connectivity (for each LGP0 nodes)*/
  for(n = 0; n < local.nnodes; n++) {
    nodes[n].nelmts=0;
    }
  for(m = 0; m < mesh.nquadrangles; m++) {
    for(k=0;k<local.nnpe;k++) {
      n=local.NIbE[m][k];
      nodes[n].nelmts++;
      }
    }
  for(n = 0; n < local.nnodes; n++) {
    nodes[n].elmts=new int[nodes[n].nelmts];
    nodes[n].nelmts=0;
    }
  for(m = 0; m < mesh.nquadrangles; m++) {
    for(k=0;k<local.nnpe;k++) {
      n=local.NIbE[m][k];
      nodes[n].elmts[nodes[n].nelmts]=m;
      nodes[n].nelmts++;
      }
    }

/* *----------------------------------------------------------------------------
  create neighbours connectivity (for each LGP0 nodes)*/
  nnghbs_max=mesh.nnghm+1;
  nghbs=new int*[local.nnodes];
  for(n = 0; n < local.nnodes; n++) {
    nodes[n].nnghbs=0;
    nghbs[n]=new int[nnghbs_max];
    for(i=0;i<nnghbs_max;i++) nghbs[n][i]=-1;
    }
  for(m = 0; m < mesh.nquadrangles; m++) {
    for(k=0;k<local.nnpe;k++) {
      n1=local.NIbE[m][k];
      for(l=0;l<local.nnpe;l++) {
        if(k==l) continue;
        n2=local.NIbE[m][l];
        acquire=1;
        for(i=0;i<nodes[n1].nnghbs;i++) {
          if(nghbs[n1][i]==n2) {
            acquire=0;
            break;
            }
          }
        if(acquire==1) {
          nghbs[n1][i]=n2;
          nodes[n1].nnghbs++;
          if(nodes[n1].nnghbs>=nnghbs_max) {
            TRAP_ERR_EXIT(-1, "LGP0 discretisation, neigbours overflow");
            }
          }
        }
      }
    }
  for(n = 0; n < local.nnodes; n++) {
    nodes[n].nghbs=new int[nodes[n].nnghbs];
    for(i=0;i<nodes[n].nnghbs;i++) nodes[n].nghbs[i]=nghbs[n][i];
    }

  for(n = 0; n < local.nnodes; n++) {
    delete[] nghbs[n];
    }
  delete[] nghbs;

  local.nodes=nodes;
  local.type=CQP0;
  return(local);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  discretisation_t initCQP1(mesh_t mesh, int type)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,k,l,m,n;
  int n1,n2,nn;
  int **nghbs;
  int nnghbs_max;
  int acquire;
  discretisation_t local;
  node_t *nodes;

  local.nnodes=mesh.nvtxs;
  local.nelmts=mesh.nquadrangles;

  local.nnpe=4;

  local.NIbE=new int*[mesh.nquadrangles];

  nodes=new node_t[local.nnodes];

/* *----------------------------------------------------------------------------
  create LGP1 nodes NIbE list (for each element)*/
  for(m = 0; m < mesh.nquadrangles; m++) {
    local.NIbE[m]=new int[local.nnpe];
    for(k=0;k<local.nnpe;k++) local.NIbE[m][k]=-1;
    for(k=0;k<local.nnpe;k++) {
      n=mesh.quadrangles[m].vertex[k];
      local.NIbE[m][k]=n;
      }
    }
  for(n = 0; n < mesh.nvtxs; n++) {
    nn=n;
    nodes[n].lon=mesh.vertices[nn].lon;
    nodes[n].lat=mesh.vertices[nn].lat;
    nodes[n].c  =mesh.vertices[nn].c;
    nodes[n].s  =mesh.vertices[nn].s;
    nodes[n].h  =mesh.vertices[nn].h;
      if(mesh.vertices[nn].nedges>0)
        nodes[n].edges=new int[mesh.vertices[nn].nedges];
      else
        nodes[n].edges=NULL;
      nodes[n].nedges=mesh.vertices[nn].nedges;
      for(i=0;i<nodes[n].nedges;i++) {
        nodes[n].edges[i]=mesh.vertices[nn].edges[i];
        }
      nodes[n].nvtces=mesh.vertices[nn].nngh+1;
      nodes[n].vtces=new int[nodes[n].nvtces];
      nodes[n].vtces[0]=nn;
      for(i=0;i<nodes[n].nvtces-1;i++) {
        nodes[n].vtces[i+1]=mesh.vertices[nn].ngh[i];
        }
    }

/* *----------------------------------------------------------------------------
  create element connectivity (for each LGP1 nodes)*/
  for(n = 0; n < local.nnodes; n++) {
    nodes[n].nelmts=0;
    }
  for(m = 0; m < mesh.nquadrangles; m++) {
    for(k=0;k<local.nnpe;k++) {
      n=local.NIbE[m][k];
      nodes[n].nelmts++;
      }
    }
  for(n = 0; n < local.nnodes; n++) {
    nodes[n].elmts=new int[nodes[n].nelmts];
    nodes[n].nelmts=0;
    }
  for(m = 0; m < mesh.nquadrangles; m++) {
    for(k=0;k<local.nnpe;k++) {
      n=local.NIbE[m][k];
      nodes[n].elmts[nodes[n].nelmts]=m;
      nodes[n].nelmts++;
      }
    }

/* *----------------------------------------------------------------------------
  create neighbours connectivity (for each LGP1 nodes)*/
//  nnghbs_max=mesh.nnghm;
  nnghbs_max=8;
  nghbs=new int*[local.nnodes];
  for(n = 0; n < local.nnodes; n++) {
    nodes[n].nnghbs=0;
    nghbs[n]=new int[nnghbs_max];
    for(i=0;i<nnghbs_max;i++) nghbs[n][i]=-1;
    }
  for(m = 0; m < mesh.nquadrangles; m++) {
    for(k=0;k<local.nnpe;k++) {
      n1=local.NIbE[m][k];
      for(l=0;l<local.nnpe;l++) {
        if(k==l) continue;
        n2=local.NIbE[m][l];
        acquire=1;
        for(i=0;i<nodes[n1].nnghbs;i++) {
          if(nghbs[n1][i]==n2) {
            acquire=0;
            break;
            }
          }
        if(acquire==1) {
          nghbs[n1][i]=n2;
          nodes[n1].nnghbs++;
          if(nodes[n1].nnghbs>nnghbs_max) {
            TRAP_ERR_EXIT(-1, "QP1 discretisation, neigbours overflow");
            }
          }
        }
      }
    }
  for(n = 0; n < local.nnodes; n++) {
    nodes[n].nghbs=new int[nodes[n].nnghbs];
    for(i=0;i<nodes[n].nnghbs;i++) nodes[n].nghbs[i]=nghbs[n][i];
    }

  for(n = 0; n < local.nnodes; n++) {
    delete[] nghbs[n];
    }
  delete[] nghbs;

  local.nodes=nodes;
  local.type=CQP1;
  return(local);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  discretisation_t initCQN1(mesh_t mesh, int type)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,k,l,m,n;
  int n1,n2,nn;
  int **nghbs;
  int nnghbs_max;
  int acquire;
  discretisation_t local;
  node_t *nodes;

  local.nnodes=mesh.nedges;
  local.nelmts=mesh.nquadrangles;

  local.nnpe=4;

  local.NIbE=new int*[mesh.nquadrangles];

  nodes=new node_t[local.nnodes];

/* *----------------------------------------------------------------------------
  create LGP1 nodes NIbE list (for each element)*/
  for(m = 0; m < mesh.nquadrangles; m++) {
    local.NIbE[m]=new int[local.nnpe];
    for(k=0;k<local.nnpe;k++) local.NIbE[m][k]=-1;
    for(k=0;k<local.nnpe;k++) {
      n=mesh.quadrangles[m].edges[k];
      local.NIbE[m][k]=n;
      }
    }
  for(n = 0; n < mesh.nedges; n++) {
    nn=n;
    nodes[n].code =mesh.edges[nn].code;
//     nodes[n].lon  =mesh.edges[nn].lon;
//     nodes[n].lat  =mesh.edges[nn].lat;
    nodes[n].lon=mesh.edges[nn].lon*180./M_PI;
    nodes[n].lat=mesh.edges[nn].lat*180./M_PI;
    nodes[n].c    =mesh.edges[nn].c;
    nodes[n].s    =mesh.edges[nn].s;
//    nodes[n].h  =mesh.edges[nn].h;
//       if(mesh.vertices[nn].nedges>0)
//         nodes[n].edges=new int[mesh.vertices[nn].nedges];
//       else
//         nodes[n].edges=NULL;
//       nodes[n].nedges=mesh.vertices[nn].nedges;
//       for(i=0;i<nodes[n].nedges;i++) {
//         nodes[n].edges[i]=mesh.vertices[nn].edges[i];
//         }
//       nodes[n].nvtces=mesh.vertices[nn].nngh+1;
//       nodes[n].vtces=new int[nodes[n].nvtces];
//       nodes[n].vtces[0]=nn;
//       for(i=0;i<nodes[n].nvtces-1;i++) {
//         nodes[n].vtces[i+1]=mesh.vertices[nn].ngh[i];
//         }
    }

/* *----------------------------------------------------------------------------
  create element connectivity (for each LGP1 nodes)*/
  for(n = 0; n < local.nnodes; n++) {
    nodes[n].nelmts=0;
    }
  for(m = 0; m < mesh.nquadrangles; m++) {
    for(k=0;k<local.nnpe;k++) {
      n=local.NIbE[m][k];
      nodes[n].nelmts++;
      }
    }
  for(n = 0; n < local.nnodes; n++) {
    nodes[n].elmts=new int[nodes[n].nelmts];
    nodes[n].nelmts=0;
    }
  for(m = 0; m < mesh.nquadrangles; m++) {
    for(k=0;k<local.nnpe;k++) {
      n=local.NIbE[m][k];
      nodes[n].elmts[nodes[n].nelmts]=m;
      nodes[n].nelmts++;
      }
    }

/* *----------------------------------------------------------------------------
  create neighbours connectivity (for each LGP1 nodes)*/
//  nnghbs_max=mesh.nnghm;
  nnghbs_max=8;
  nghbs=new int*[local.nnodes];
  for(n = 0; n < local.nnodes; n++) {
    nodes[n].nnghbs=0;
    nghbs[n]=new int[nnghbs_max];
    for(i=0;i<nnghbs_max;i++) nghbs[n][i]=-1;
    }
  for(m = 0; m < mesh.nquadrangles; m++) {
    for(k=0;k<local.nnpe;k++) {
      n1=local.NIbE[m][k];
      for(l=0;l<local.nnpe;l++) {
        if(k==l) continue;
        n2=local.NIbE[m][l];
        acquire=1;
        for(i=0;i<nodes[n1].nnghbs;i++) {
          if(nghbs[n1][i]==n2) {
            acquire=0;
            break;
            }
          }
        if(acquire==1) {
          nghbs[n1][i]=n2;
          nodes[n1].nnghbs++;
          if(nodes[n1].nnghbs>nnghbs_max) {
            TRAP_ERR_EXIT(-1, "QP1 discretisation, neigbours overflow");
            }
          }
        }
      }
    }
  for(n = 0; n < local.nnodes; n++) {
    nodes[n].nghbs=new int[nodes[n].nnghbs];
    for(i=0;i<nodes[n].nnghbs;i++) nodes[n].nghbs[i]=nghbs[n][i];
    }

  for(n = 0; n < local.nnodes; n++) {
    delete[] nghbs[n];
    }
  delete[] nghbs;

  local.nodes=nodes;
  local.type=CQN1;
  return(local);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  discretisation_t init_discretisation(const mesh_t &mesh, int type)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  switch(type){
    case LGP0:
      return initLGP0(mesh,type);
    case LGP1:
      return initLGP1(mesh,type);
    case DGP1:
      return initDGP1(mesh,type);
    case DNP1:
      return initDNP1(mesh,type);
    case NCP1:
      return initNCP1(mesh,type);
    case LGP2:
      return initLGP2(mesh,type);
    case CQP1:
      return initCQP1(mesh,type);
    case CQN1:
      return initCQN1(mesh,type);
    default:
      TRAP_ERR_EXIT(ENOEXEC,"not coded yet for discretisation %d\n",type);
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int discretisation_init(mesh_t *mesh, int discretisation, int compute_massmatrix, int context)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// initialises discretisation and, eventually, mass matrix
/**
\param *mesh
\param discretisation
\param compute_massmatrix If 1(default) : compute mass matrix.
\returns always 0
*/
/*----------------------------------------------------------------------------*/
{
  int status;
  hypermatrix_t *matrix;
  int verbose=0;
//   int context=0;
  discretisation_t *descriptor=get_descriptor_address(*mesh,discretisation);
  if(!descriptor)
    TRAP_ERR_EXIT(ENOEXEC,"could not find descriptor\n");
  const char *disc_name=discretisation_name(discretisation);
  
  switch (discretisation) {
    case LGP0:
/*------------------------------------------------------------------------------
      LGP0 elevation discretisation */
      break;

    case LGP1:
/*------------------------------------------------------------------------------
      Continuous LGP1 discretisation */
      break;

    case DGP1:
/*------------------------------------------------------------------------------
      Discontinuous LGP1 discretisation */
      break;

    case NCP1:
/*------------------------------------------------------------------------------
      Discontinuous NCP1 discretisation */
      break;

    case DNP1:
/*------------------------------------------------------------------------------
      Discontinuous LGP1 discretisation */
      if(gParallel!=0) {
/* *----------------------------------------------------------------------------
        construct unpartitionned LGP2 numbering to access globalized numbering */
        if(mesh->origin==0) {
          TRAP_ERR_EXIT(-1, "unpartitionned mesh not available");
          }
        if(mesh->origin->DNP1descriptor.nnodes==0) {
          mesh->origin->DNP1descriptor=NIbE_DNP1(*(mesh->origin), discretisation);
          }
        }
      break;

    case LGP2:
/*-----------------------------------------------------------------------------
      Continuous LGP2 discretisation */
      if(gParallel!=0) {
/* *----------------------------------------------------------------------------
        construct unpartitionned LGP2 numbering to access globalized numbering */
        if(mesh->origin==0) {
          TRAP_ERR_EXIT(-1, "unpartitionned mesh not available");
          }
        if(mesh->origin->LGP2descriptor.nnodes==0) {
          mesh->origin->LGP2descriptor=NIbE_LGP2(*(mesh->origin), discretisation);
          }
        }
      break;

    case DGP2:
/*-----------------------------------------------------------------------------
      Disontinuous DGP2 discretisation */
//       if(gParallel!=0) {
// /* *----------------------------------------------------------------------------
//         construct unpartitionned DGP2 numbering to access globalized numbering */
//         if(mesh->origin==0) {
//           TRAP_ERR_EXIT(-1, "unpartitionned mesh not available\n");
//           }
//         if(mesh->origin->DGP2descriptor.nnodes==0) {
//           mesh->origin->DGP2descriptor=NIbE_DGP2(*(mesh->origin), discretisation);
//           }
//         }
      break;

    case CQP0:
/*------------------------------------------------------------------------------
      Discontinuous NCP1 discretisation */
      if(mesh->CQP0descriptor.nnodes!=0) {
        if(verbose==1) printf("CQP0 descriptor already initialized, skip instruction\n");
        }
      else {
        mesh->CQP0descriptor=initCQP0(*mesh, discretisation);
//        status=init_distributor(*mesh, &(mesh->CQP0descriptor), context);
        }
      matrix=&(mesh->CQP0descriptor.massmatrix);
      if(matrix->packed!=0) {
        if(verbose==1) printf("CQP0 mass matrix already initialized, skip instruction\n");
        }
      else {
        matrix->ordering=new ordering_t();
        status=dMassMatrix(*mesh, mesh->CQP0descriptor, matrix);
        }
      break;

    case CQP1:
/*------------------------------------------------------------------------------
      Discontinuous NCP1 discretisation */
      break;

    case CQN1:
/*------------------------------------------------------------------------------
      Discontinuous NCP1 discretisation */
      break;

//     case INTRINSIC:
//       printf("intrinsic definition, skip instruction\n");
//       break;

    default:
      TRAP_ERR_EXIT(-1, "illegal discretisation");
      break;
    }
  
  if(descriptor->nnodes!=0) {
    if(verbose==1) printf("%s descriptor already initialized, skip instruction",disc_name);
    }
  else {
    *descriptor=init_discretisation(*mesh, discretisation);
    }
  /// HERE !!! put test on lon lat
  if(compute_massmatrix==1) {
    matrix=&(descriptor->massmatrix);
    if(matrix->packed!=0) {
      if(verbose==1) printf("%s mass matrix already initialized, skip instruction",disc_name);
      }
    else {
      matrix->ordering=new ordering_t();
      status=dMassMatrix(*mesh, *descriptor, matrix);
      }
    }
  
  return(0);
}


