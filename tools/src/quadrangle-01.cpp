
/*******************************************************************************

  T-UGO tools, 2006-2009

  Unstructured Ocean Grid initiative

Contributors:

  Florent Lyard      LEGOS/CNRS, Toulouse, France
  Yoann Le Bars      LEGOS, Toulouse, France (PhD)
  Laurent Roblou     LEGOS, Toulouse, France
  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
  Fr�d�ric Dupont    Universit� de Laval � Qu�bec, Canada

E-mail: florent.lyard@legos.obs-mip.fr

*******************************************************************************/

#include <config.h>

#include <stdio.h>
#include <string.h>

#include "tools-structures.h"

#include "fe.def"
#include "fe.h"

#include "map.h"
#include "geo.h"
#include "swap.h"
#include "maths.h"
#include "matrix.h"
#include "netcdf-proto.h"
#include "poc-netcdf-data.hpp"

#include "polygons.hpp"
#include "exceptions.hpp"

#include "polygones.h"

extern  int fe_savemeshNC2D_new(const char *filename, mesh_t & mesh, int option);
extern  int fe_savemeshNC3D(const char *filename, gmesh_t<quadrangle_t> & mesh, int option);

extern int quadrangle_readneigh (const string &polygonsFileName, mesh_t & final);
extern  int fe_SeekQuadrangles(mesh_t *mesh, Polygons::ElementaryCycles &elementaryCycles,int maxsize);

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int fe_import_Qvertices(grid_t grid, mesh_t *mesh, char *landmask, float *topo, int optimize)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int    i,j,k,m,n;
  int    m1,m2,n1,n2;
  double x,y;
  int *incidence;

  incidence=new int[grid.nx*grid.ny];
  for(m=0;m<grid.nx*grid.ny;m++) incidence[m]=-1;

/* *------------------------------------------------------------------------
  minimize band width by numbering along the shortest side*/
  mesh->nvtxs=0;
  if( (grid.nx<grid.ny) || (optimize==0) ) {
    for(j=0;j<grid.ny;j++) {
      for(i=0;i<grid.nx;i++) {
        m=grid.nx*j+i;
        if(landmask[m]>0) {
          incidence[m]=mesh->nvtxs;
          mesh->nvtxs++;
          }
        }
      }
    }
  else {
    for(i=0;i<grid.nx;i++) {
      for(j=0;j<grid.ny;j++) {
        m=grid.nx*j+i;
        if(landmask[m]>0) {
          incidence[m]=mesh->nvtxs;
          mesh->nvtxs++;
          }
        }
      }
    }

  mesh->vertices= new vertex_t[mesh->nvtxs];
  mesh->nlayers=grid.nz-1;
  mesh->nlevels=grid.nz;

  mesh->nnghm=4;

/* *------------------------------------------------------------------------
  initialize vertices*/
  for(j=0;j<grid.ny;j++) {
    for(i=0;i<grid.nx;i++) {
      m=grid.nx*j+i;
      if(landmask[m]>0.) {
        x= grid.x[m];
        y= grid.y[m];
        n=incidence[m];
        mesh->vertices[n].lon=x;
        mesh->vertices[n].lat=y;
        mesh->vertices[n].nngh=0;
        if(topo!=0) mesh->vertices[n].h=topo[m];
        mesh->vertices[n].code=0;
        mesh->vertices[n].ngh=new int[mesh->nnghm];
        for(k=0;k<mesh->nnghm;k++) mesh->vertices[n].ngh[k]=-1;
//         mesh->vertices[n].sigma  =new double[mesh->nlevels];
//         mesh->vertices[n].zlevels=new double[mesh->nlevels];
//         for (l=0;l<mesh->nlevels;l++) {
//           mesh->vertices[n].zlevels[l]=grid.z[m+l*grid.ny*grid.nx];
//           }
//         for (l=0;l<mesh->nlevels;l++) {
//           mesh->vertices[n].sigma[l]=-grid.z[m+l*grid.ny*grid.nx]/topo[m];
//           }
        }
      }
    }

/* *------------------------------------------------------------------------
  initialize vertices connections*/
  for(j=0;j<grid.ny;j++) {
    for(i=0;i<grid.nx-1;i++) {
      m1=grid.nx*j+i;
      m2=grid.nx*j+i+1;
      if(incidence[m1]==-1) continue;
      if(incidence[m2]==-1) continue;
      n1=incidence[m1];
      n2=incidence[m2];
      k=mesh->vertices[n1].nngh;
      mesh->vertices[n1].ngh[k]=n2;
      k=mesh->vertices[n2].nngh;
      mesh->vertices[n2].ngh[k]=n1;
      mesh->vertices[n1].nngh++;
      mesh->vertices[n2].nngh++;
      }
    }

  for(j=0;j<grid.ny-1;j++) {
    for(i=0;i<grid.nx;i++) {
      m1=grid.nx*j+i;
      m2=grid.nx*j+i+grid.nx;
      if(incidence[m1]==-1) continue;
      if(incidence[m2]==-1) continue;
      n1=incidence[m1];
      n2=incidence[m2];
      k=mesh->vertices[n1].nngh;
      mesh->vertices[n1].ngh[k]=n2;
      k=mesh->vertices[n2].nngh;
      mesh->vertices[n2].ngh[k]=n1;
      mesh->vertices[n1].nngh++;
      mesh->vertices[n2].nngh++;
      }
    }
  
  delete[] incidence;
  
  return(0);
}

// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
// int fe_import_Qvertices(grid_t grid, grid_t zgrid, mesh_t *mesh, char *landmask, float *topo)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   int    i,j,k,l,kk,ll,m,n,nn,I0,J0;
//   int    m1,m2,n1,n2;
//   int    nitems,status;
//   int    nz,type1,type2;
//   FILE   *in,*out;
//   double lon,lat,angle,radius,x,y;
//   double tx,ty,nx,ny,x_shift,y_shift;
//   char   line[1024],nodefile[1024];
//   float  z;
//   int *incidence;
//   int verbose=0;
// 
//   incidence=new int[grid.nx*grid.ny];
//   for(m=0;m<grid.nx*grid.ny;m++) incidence[m]=-1;
// 
// /* *------------------------------------------------------------------------
//   minimize band width by numbering along the shortest side*/
//   mesh->nvtxs=0;
//   if(grid.nx<grid.ny) {
//     for(j=0;j<grid.ny;j++) {
//       for(i=0;i<grid.nx;i++) {
//         m=grid.nx*j+i;
//         if(landmask[m]>0) {
//           incidence[m]=mesh->nvtxs;
//           mesh->nvtxs++;
//           }
//         }
//       }
//     }
//   else {
//     for(i=0;i<grid.nx;i++) {
//       for(j=0;j<grid.ny;j++) {
//         m=grid.nx*j+i;
//         if(landmask[m]>0) {
//           incidence[m]=mesh->nvtxs;
//           mesh->nvtxs++;
//           }
//         }
//       }
//     }
// 
//   mesh->vertices= new vertex_t[mesh->nvtxs];
//   mesh->nlayers=grid.nz-1;
//   mesh->nlevels=grid.nz;
// 
//   mesh->nnghm=4;
// 
// /* *------------------------------------------------------------------------
//   initialize vertices*/
//   for(j=0;j<grid.ny;j++) {
//     for(i=0;i<grid.nx;i++) {
//       m=grid.nx*j+i;
//       if(landmask[m]>0.) {
//         x= grid.x[m];
//         y= grid.y[m];
//         n=incidence[m];
//         mesh->vertices[n].lon=x;
//         mesh->vertices[n].lat=y;
//         mesh->vertices[n].nngh=0;
//         if(topo!=0) mesh->vertices[n].h=topo[m];
//         mesh->vertices[n].code=0;
//         mesh->vertices[n].ngh=new int[mesh->nnghm];
//         for(k=0;k<mesh->nnghm;k++) mesh->vertices[n].ngh[k]=-1;
// //         mesh->vertices[n].sigma  =new double[mesh->nlevels];
// //         mesh->vertices[n].zlevels=new double[mesh->nlevels];
// //         for (l=0;l<mesh->nlevels;l++) {
// //           mesh->vertices[n].zlevels[l]=grid.z[m+l*grid.ny*grid.nx];
// //           }
// //         for (l=0;l<mesh->nlevels;l++) {
// //           mesh->vertices[n].sigma[l]=-grid.z[m+l*grid.ny*grid.nx]/topo[m];
// //           }
//         }
//       }
//     }
// 
// /* *------------------------------------------------------------------------
//   initialize vertices connections*/
//   for(j=0;j<grid.ny;j++) {
//     for(i=0;i<grid.nx-1;i++) {
//       m1=grid.nx*j+i;
//       m2=grid.nx*j+i+1;
//       if(incidence[m1]==-1) continue;
//       if(incidence[m2]==-1) continue;
//       n1=incidence[m1];
//       n2=incidence[m2];
//       k=mesh->vertices[n1].nngh;
//       mesh->vertices[n1].ngh[k]=n2;
//       k=mesh->vertices[n2].nngh;
//       mesh->vertices[n2].ngh[k]=n1;
//       mesh->vertices[n1].nngh++;
//       mesh->vertices[n2].nngh++;
//       }
//     }
// 
//   for(j=0;j<grid.ny-1;j++) {
//     for(i=0;i<grid.nx;i++) {
//       m1=grid.nx*j+i;
//       m2=grid.nx*j+i+grid.nx;
//       if(incidence[m1]==-1) continue;
//       if(incidence[m2]==-1) continue;
//       n1=incidence[m1];
//       n2=incidence[m2];
//       k=mesh->vertices[n1].nngh;
//       mesh->vertices[n1].ngh[k]=n2;
//       k=mesh->vertices[n2].nngh;
//       mesh->vertices[n2].ngh[k]=n1;
//       mesh->vertices[n1].nngh++;
//       mesh->vertices[n2].nngh++;
//       }
//     }
// 
//   return(0);
// }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int quadrangle_list (mesh_t & mesh, int MinSize, int MaxSize) // throw ()

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
using namespace Polygons;

  int k,/*n,*/status;
  size_t add_ntriangles,add_nquadrangles;
  ElementaryCycles elementaryCycles;
  size_t size,count;

//   try {
/* *----------------------------------------------------------------------------
    Read graph from polygonsFileName */
    status=fe_SeekQuadrangles(&mesh, elementaryCycles, MaxSize+1);

    mesh.nquadrangles=0;
    mesh.ntriangles=0;
    mesh.nelements=0;
    
    add_ntriangles=0;
    add_nquadrangles=0;
    
    for (ElementaryCycles::const_iterator i = elementaryCycles.begin(); i != elementaryCycles.end(); i++) {
      const Adjacency cycle = *i;
      size=cycle.size();
      switch(size) {
        case 0:
        case 1:
        case 2:
          printf("Undersized cycle %d\n",size);
          break;
        case 3:
          mesh.ntriangles++;
          mesh.nelements++;
          break;
        case 4:
          mesh.nquadrangles++;
          mesh.nelements++;
          break;
        case 5:
          add_ntriangles++;
          mesh.nelements++;
          add_nquadrangles++;
          mesh.nelements++;
          break;
        default:
          printf("Oversized cycle %d\n",size);
          break;
        }
      }

//     mesh.nvtxs=tmp.nvtxs;
//     mesh.vertices=new vertex_t[mesh.nvtxs];
//     for(n=0;n<mesh.nvtxs;n++) {
//       mesh.vertices[n].lon=tmp.vertices[n].lon;
//       mesh.vertices[n].lat=tmp.vertices[n].lat;
//       }

    if(MinSize<=3  && MaxSize>=3) {
/* *----------------------------------------------------------------------------
      create triangles*/
      mesh.triangles=new triangle_t[mesh.ntriangles+add_ntriangles];

      count=0;
      for (ElementaryCycles::const_iterator i = elementaryCycles.begin(); i != elementaryCycles.end(); i++) {
        const Adjacency cycle = *i;
        size=cycle.size();
        if(size==3) {
          k=0;
          for (Adjacency::const_iterator place = cycle.begin(); place != cycle.end(); place++) {
            mesh.triangles[count].vertex[k]=*place;
            k++;
            }
          count++;
          }
        }
      }

    if(MinSize<=4  && MaxSize>=4) {
/* *----------------------------------------------------------------------------
      create quadrangles*/
      mesh.quadrangles=new quadrangle_t[mesh.nquadrangles+add_nquadrangles];

      count=0;
      for (ElementaryCycles::const_iterator i = elementaryCycles.begin(); i != elementaryCycles.end(); i++) {
        const Adjacency cycle = *i;
        size=cycle.size();
        if(size==4) {
          k=0;
          for (Adjacency::const_iterator place = cycle.begin(); place != cycle.end(); place++) {
            mesh.quadrangles[count].vertex[k]=*place;
            k++;
            }
          count++;
          }
        }
      }
    
    for (ElementaryCycles::iterator i = elementaryCycles.begin(); i != elementaryCycles.end(); i++) {
      (*i).clear();
      }
    
//    ElementaryCycles.clear();
    
//    }

/* *----------------------------------------------------------------------------
  Part of the code executed if a problem occur during file reading*/
//   catch (TugoExceptions::ReadError file) {
//     char errorMessage [256];
//     sprintf(errorMessage, "An error occured when reading the file \"%s\".",file.fileName().c_str());
//     check_error(-2, errorMessage, __LINE__, __FILE__, 1);
//     }
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int fe_convert(grid_t grid, char  *landmask, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int m,count;
//  plg_array_t set;
  mesh_t tmp, final;
  string cleaned;
  float *topo=0, *area=0;
  cdfvar_t variable;
  double mask=-1.;
  char *output;
 
  if(landmask==0) {
    landmask=new char[grid.nx*grid.ny];
    for(m=0;m<grid.nx*grid.ny;m++) landmask[m]=1;
    }
  
  mesh.nvtxs=grid.nx*grid.ny;
  mesh.vertices=new vertex_t[mesh.nvtxs];
  
  status=fe_import_Qvertices(grid, &mesh, landmask, topo, 0);

/* *----------------------------------------------------------------------------
  force neigbours list's consistency (neighbourship symetry)*/
  count=fe_chkvertex_consistency(&tmp,1);
 
/* *----------------------------------------------------------------------------
  remove vertices with no connexions at all*/
  status=fe_cleanvertices(&tmp, false);

/* *----------------------------------------------------------------------------
  save "consistent" mesh*/
  status=fe_savemesh("tmp.nei",MESH_FILE_FORMAT_TRIGRID,tmp);

/* *----------------------------------------------------------------------------
  now process cycles*/
  status= quadrangle_list(mesh);
  for(m=0;m<mesh.nquadrangles;m++) {
    status=fe_initaffine_spherical(mesh, &(mesh.quadrangles[m]), m);
    }
  
  output=strdup("test.nc");
  status=fe_savemeshNC2D_new(output,  mesh, 1);
  
  area=new float[mesh.nquadrangles];
  for(m=0;m<mesh.nquadrangles;m++) {
    area[m]=mesh.quadrangles[m].TrueArea;
    }
  
  variable=poc_variable_UG2D("area", mask, "m²", 1.0, 0.0, "cell_area","NQUADRANGLES");
  status   = cdf_createvariable(output, &(variable));
  status   = poc_put_UG2D(output, mesh, variable.id, area);
  variable.destroy();

  delete[] area;
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int fe_FGrid2Quadrangle(grid_t grid, mesh_t & mesh, const char *echo, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=0;
  int k,i,j,m,n,n1,n2;
  
  mesh.nvtxs=grid.nx*grid.ny;
  mesh.vertices=new vertex_t[mesh.nvtxs];
  
  mesh.nlayers=grid.nz-1;
  mesh.nlevels=grid.nz;

  mesh.nnghm=4;

/* *------------------------------------------------------------------------
  initialize vertices*/
  for(j=0;j<grid.ny;j++) {
    for(i=0;i<grid.nx;i++) {
      m=grid.nx*j+i;
      double x= grid.x[m];
      double y= grid.y[m];
      n=m;
      mesh.vertices[n].lon=x;
      mesh.vertices[n].lat=y;
      mesh.vertices[n].nngh=0;
      mesh.vertices[n].code=0;
      }
    }

  for(j=1;j<grid.ny-1;j++) {
    for(i=1;i<grid.nx-1;i++) {
      n=grid.nx*j+i;
      mesh.vertices[n].nngh=4;
      }
    }
  
  for(j=1;j<grid.ny-1;j++) {
    i=0;
    n=grid.nx*j+i;
    mesh.vertices[n].nngh=3;
    i=grid.nx-1;
    n=grid.nx*j+i;
    mesh.vertices[n].nngh=3;
    }
  
  for(i=1;i<grid.nx-1;i++) {
    j=0;
    n=grid.nx*j+i;
    mesh.vertices[n].nngh=3;
    j=grid.ny-1;
    n=grid.nx*j+i;
    mesh.vertices[n].nngh=3;
    }
  
    
/* *------------------------------------------------------------------------
  initialize vertices connections*/
  for(j=0;j<grid.ny;j++) {
    for(i=0;i<grid.nx-1;i++) {
      n1=grid.nx*j+i;
      n2=grid.nx*j+i+1;
      k=mesh.vertices[n1].nngh;
      mesh.vertices[n1].ngh[k]=n2;
      k=mesh.vertices[n2].nngh;
      mesh.vertices[n2].ngh[k]=n1;
      mesh.vertices[n1].nngh++;
      mesh.vertices[n2].nngh++;
      }
    }

  for(j=0;j<grid.ny-1;j++) {
    for(i=0;i<grid.nx;i++) {
      n1=grid.nx*j+i;
      n2=grid.nx*j+i+grid.nx;
      k=mesh.vertices[n1].nngh;
      mesh.vertices[n1].ngh[k]=n2;
      k=mesh.vertices[n2].nngh;
      mesh.vertices[n2].ngh[k]=n1;
      mesh.vertices[n1].nngh++;
      mesh.vertices[n2].nngh++;
      }
    }
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int fe_FGrid2Quadrangle(grid_t grid, grid_t cgrid, char  *landmask, float *topo, float topomask, mesh_t & mesh, const char *echo, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int m,count;
  int verbose=(debug==true);
  
  mesh_t tmp, final;
  string cleaned;
  float  *area=0;
  cdfvar_t variable;
  char *output;
  
  if(verbose) {
    printf("#################################################################\n");
    printf("convert vorticity structured grid into quadrangle mesh\n");
    }
 
  if(landmask==0) {
    landmask=new char[grid.nx*grid.ny];
    for(m=0;m<grid.nx*grid.ny;m++) landmask[m]=1;
    }
  
  status=fe_import_Qvertices(grid, &mesh, landmask, topo, 0);

/* *----------------------------------------------------------------------------
  force neigbours list's consistency (neighbourship symetry)*/
  count=fe_chkvertex_consistency(&mesh,1);
 
/* *----------------------------------------------------------------------------
  remove vertices with no connexions at all*/
  status=fe_cleanvertices(&mesh, false);

/* *----------------------------------------------------------------------------
  save "consistent" mesh*/
  if(debug) status=fe_savemesh("fe_FGrid2Quadrangle.nei",MESH_FILE_FORMAT_TRIGRID,mesh);

/* *----------------------------------------------------------------------------
  now process cycles*/
  status= quadrangle_list(mesh, 4, 4);
  int nprocs __attribute__((unused)) =initialize_OPENMP(-1,0);
#ifdef OPEN_MP
  #pragma omp parallel for private(status) if(nprocs>1)
#endif
  for(m=0;m<mesh.nquadrangles;m++) {
    status=fe_initaffine_spherical(mesh, &(mesh.quadrangles[m]), m);
    }
  
  status=fe_edgetable_Q(&mesh, verbose);
  if(status!=0) return(status);
  
  if(echo==0) return(status);
  
  output=poc_strdup(echo);
  status=fe_savemeshNC2D_new(output,  mesh, 1);
    
  area=new float[mesh.nquadrangles];
  for(m=0;m<mesh.nquadrangles;m++) {
    area[m]=mesh.quadrangles[m].TrueArea;
    }
  
  variable << poc_variable_UG2D("area", NC_FILL_FLOAT, "m²", 1.f, 0.f, "cell_area","NQUADRANGLES");
  status   = cdf_createvariable(output, &(variable));
  status   = poc_put_UG2D(output, mesh, variable.id, area);
  variable.destroy();
  delete[] area;
  
  float *topography=new float[mesh.nquadrangles];
/*------------------------------------------------------------------------------
  spherical coordinates interpolation will fail in polar regions */

//   for(m=0;m<mesh.nquadrangles;m++) {
//     double t,p;
//     status=fe_position(mesh, mesh.quadrangles[m], &t, &p,0);
//     status=map_interpolation(grid,topo,topomask,t,p,&topography[m]);
//     if(status!=0) {
//       status=map_interpolation(grid,topo,topomask,t,p,&topography[m]);
//       }
//     }

/*------------------------------------------------------------------------------
  thus use cartesian coordinates interpolation */
  double *x=new double[mesh.nquadrangles];
  double *y=new double[mesh.nquadrangles];
  for(int m=0;m<mesh.nquadrangles;m++) {
    status=fe_position(mesh, mesh.quadrangles[m], &x[m], &y[m],0);
    }
  status=geo_to_projection (cgrid.proj4options, x, y, mesh.nquadrangles);
  for(int m=0;m<mesh.nquadrangles;m++) {
    status=map_interpolation(cgrid,topo,topomask, x[m], y[m],&topography[m]);
    if(status!=0) {
      status=map_interpolation(cgrid,topo,topomask, x[m], y[m],&topography[m]);
      }
    }
  delete[] x;
  delete[] y;

  variable << poc_variable_UG2D("topography", topomask, "m", 1.f, 0.f, "topography","NQUADRANGLES");
  status   = cdf_createvariable(output, &(variable));
  status   = poc_put_UG2D(output, mesh, variable.id, topography);
  variable.destroy();
  delete[] topography;
  delete[] output;
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int fe_ZGrid2Quadrangle(const grid_t & zgrid, const grid_t & fgrid, char  *landmask, float *topo, float topomask, mesh_t & mesh, int FixPinched, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  convert structured grid and landmask in quadrangle mesh
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int status;
  int i,j,k,l,m,n,count;

  mesh_t tmp, final;
  string cleaned;
  
  const size_t fs=fgrid.Hsize();
  const size_t zs=zgrid.Hsize();
  
  int *Vincidence=new int[fs];
  int *Qincidence=new int[zs];
  
  int shift=0, remove=0;
  
/*------------------------------------------------------------------------------
  fucking magic, combat structured grid ill-description*/
  if(fs==zs) {
/*------------------------------------------------------------------------------
    eliminate artificial first row and colummn of z-grid (MARS...)*/
    shift=1;
    printf("discard first row and column of tracer grid\n");
    }
  else if(fs<zs) {
/*------------------------------------------------------------------------------
    eliminate artificial first and last row and colummn of z-grid (SYMPHONIE...)*/
    shift=1;
    remove=1;
    printf("discard first and last row and column of tracer grid\n");
    }
 
 
  if(landmask==0) {
    landmask=new char[zs];
    for(m=0;m<zs;m++) landmask[m]=1;
    }
  
  for(n=0;n<fs;n++) {
    Vincidence[n]=-1;
    }
  for(n=0;n<zs;n++) {
    Qincidence[n]=-1;
    }
  
/*------------------------------------------------------------------------------
  preventive anti-pinching changes*/
  if(FixPinched==1)
  for(j=1;j<fgrid.ny-1;j++) {
    for(i=1;i<fgrid.nx-1;i++) {
      m=fgrid.Hindex(i,j);
/*------------------------------------------------------------------------------
      count number of cells */
      int count=0;
      for (l=0;l<2;l++) {
        for (k=0;k<2;k++) {
          n=zgrid.Hindex(i+k,j+l);
          if(n<0 or zs<=n) TRAP_ERR_EXIT(ENOEXEC,"fail safe : %d not in [0;%u[\n",n,zs);
          if (landmask[n]==1) count++;
          }
        }
      if(count==2) {
/*------------------------------------------------------------------------------
        might be pinched, see if diagonal */
        int n1=zgrid.Hindex(i,j);
        int n2=zgrid.Hindex(i+1,j+1);
        if( (landmask[n1]==1) && (landmask[n2]==1) ) {
          if(verbose==1) printf("vorticity node %d (%d,%d) is pinched \n",m,i,j);
/*------------------------------------------------------------------------------
          add the 2 missing cells */
          n1=zgrid.Hindex(i,j+1);
          n2=zgrid.Hindex(i+1,j);
          landmask[n1]=1;
          landmask[n2]=1;
          }
        if( (landmask[n1]==0) && (landmask[n2]==0) ) {
          if(verbose==1) printf("vorticity node %d (%d,%d) is pinched \n",m,i,j);
/*------------------------------------------------------------------------------
          add the 2 missing cells */
          n1=zgrid.Hindex(i,j);
          n2=zgrid.Hindex(i+1,j+1);
          landmask[n1]=1;
          landmask[n2]=1;
          }
        }
      }
    }

  mesh.nnghm=4;
  
/*------------------------------------------------------------------------------
  count quadrangles (valid tracer cells) */
  mesh.nquadrangles=0;
  
  for(j=shift;j<zgrid.ny-remove;j++) {
    for(i=shift;i<zgrid.nx-remove;i++) {
      n=zgrid.Hindex(i,j);
      if (landmask[n]==1) mesh.nquadrangles++;
      }
    }
  
  if(mesh.nquadrangles==0) {
    return(-1);
    }
  
  printf("%s : %d valid quadrangles found over %d\n",__func__,mesh.nquadrangles,zs);
  
  mesh.quadrangles=new quadrangle_t[mesh.nquadrangles];

/*------------------------------------------------------------------------------
  count vertices*/
  mesh.nvtxs=0;

  for(n=0;n<fs;n++) {
    Vincidence[n]=-1;
    }
  
  for(j=shift;j<zgrid.ny-remove;j++) {
    for(i=shift;i<zgrid.nx-remove;i++) {
      m=zgrid.Hindex(i,j);
      if (landmask[m]==0) continue;
/*------------------------------------------------------------------------------
      pressure cell is valid, thus the 4 vorticity nodes arounds */
      for (l=0;l<2;l++) {
        for (k=0;k<2;k++) {
          m=fgrid.Hindex(i+k-shift,j+l-shift);
          if(Vincidence[m]==-1) {
            Vincidence[m]=mesh.nvtxs;
            mesh.nvtxs++;
            }
          }
        }
      }
    }

  if(mesh.nvtxs==0) {
    return(-1);
    }
  
  printf("%s : %d valid vertices found over %d\n",__func__,mesh.nvtxs,fs);
  
/*------------------------------------------------------------------------------
  set up vertices attributes (f-nodes) */
  mesh.vertices=new vertex_t[mesh.nvtxs];
//  for(m=0;m<fs;m++) {
  count=0;
  for(j=0;j<fgrid.ny;j++) {
    for(i=0;i<fgrid.nx;i++) {
      m=fgrid.Hindex(i,j);
      if(Vincidence[m]!=-1) {
        double x= fgrid.x[m];
        double y= fgrid.y[m];
        n=Vincidence[m];
        mesh.vertices[n].lon=x;
        mesh.vertices[n].lat=y;
        mesh.vertices[n].nngh=0;
        if(topo!=0) {
          mesh.vertices[n].h=topo[m];
          if(topo[m]==topomask) {
            mesh.vertices[n].h=0.;
            }
          }
        mesh.vertices[n].code=0;
        mesh.vertices[n].ngh=new int[mesh.nnghm];
        for(k=0;k<mesh.nnghm;k++) mesh.vertices[n].ngh[k]=-1;
        count++;
        }
      }
    }
  
  if(count!=mesh.nvtxs) {
    return(-1);
    }

// /*------------------------------------------------------------------------------
//   debugging verification */
//   status=fe_spherical(mesh, fgrid.proj4options);
//   frame_t frame;
//   status=fe_minmax(mesh, frame);
//   status=map_projection_backward(fgrid, fgrid.proj4options);
  
/*------------------------------------------------------------------------------
  set quadrangle's vertex array*/
  count=0;
  for(j=shift;j<zgrid.ny-remove;j++) {
    for(i=shift;i<zgrid.nx-remove;i++) {
      m=zgrid.Hindex(i,j);
      if (landmask[m]==0) continue;
/*------------------------------------------------------------------------------
      pressure cell is valid, thus is quadrangle */
      n=Vincidence[fgrid.Hindex(i-shift,j-shift)];
      if(n==-1) {
        TRAP_ERR_EXIT(-1, "quadrangle mesh construction failed failed\n");
        }
      mesh.quadrangles[count].vertex[0]=n;
      n=Vincidence[fgrid.Hindex(i-shift+1,j-shift)];
      if(n==-1) {
        TRAP_ERR_EXIT(-1, "quadrangle mesh construction failed failed\n");
        }
      mesh.quadrangles[count].vertex[1]=n;
      n=Vincidence[fgrid.Hindex(i-shift+1,j-shift+1)];
      if(n==-1) {
        TRAP_ERR_EXIT(-1, "quadrangle mesh construction failed failed\n");
        }
      mesh.quadrangles[count].vertex[2]=n;
      n=Vincidence[fgrid.Hindex(i-shift,j-shift+1)];
      if(n==-1) {
        TRAP_ERR_EXIT(-1, "quadrangle mesh construction failed failed\n");
        }
      mesh.quadrangles[count].vertex[3]=n;
/*------------------------------------------------------------------------------
      keep track of C-GRID index */
      mesh.quadrangles[count].origin=m;
      count++;
      }
    }
  
  if(count!=mesh.nquadrangles) {
    return(-1);
    }

/*------------------------------------------------------------------------------
  set up vertices' neighbour list*/
  for(m=0;m<mesh.nquadrangles;m++) {
    for(k=0;k<4;k++) {
      int n1=mesh.quadrangles[m].vertex[k];
      int n2=mesh.quadrangles[m].vertex[(k+1)%4];
      status=fe_connectvertices(mesh, n1, n2);
      if(status!=0) {
        TRAP_ERR_EXIT(-1, "quadrangle mesh construction failed failed\n");
        }
      }
    }

/*------------------------------------------------------------------------------
  force neigbours list's consistency (neighbourship symetry)*/
  count=fe_chkvertex_consistency(&mesh,1);
 
/*------------------------------------------------------------------------------
  remove vertices with no connexions at all*/
  status=fe_cleanvertices(&mesh, false);
  
/*------------------------------------------------------------------------------
  remove isolated (disconnected) mesh parts*/
  status=fe_cleanisolated(mesh, 1);

/*------------------------------------------------------------------------------
  reconstruct elements */
  status= quadrangle_list(mesh,4,4);
  
  for(m=0;m<mesh.nquadrangles;m++) {
    status=fe_initaffine_spherical(mesh, &(mesh.quadrangles[m]), m);
    }
  status=fe_edgetable_Q(&mesh,0);
  
  int nPinched =0;
  for(n=0;n<mesh.nvtxs;n++) {
    if( (mesh.vertices[n].nngh==4) && (mesh.vertices[n].nelmts==2) ) {
      printf("node %d is pinched \n",n);
      printf("position lon=%lf lat=%lf\n",mesh.vertices[n].lon,mesh.vertices[n].lat);
      nPinched++;
      }
    }
  
  int stopon_EdgeError =1;
  int stopon_PinchError=0;
  status=fe_codetable1(&mesh, 0, stopon_EdgeError, stopon_PinchError);
   
/*------------------------------------------------------------------------------
  if grid was not safe, backward incidence is to be reworked*/
  delete[] Vincidence;
  delete[] Qincidence;
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int quadrangle_ImportStructured(const char *gridfile, const char *maskfile, mesh_t & mesh, metagrid_t & input, double tag, string rootname, const char *wprojection, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    build the quadrangle mesh from vorticity nodes.
    
    Normally fgrid is needed for vertices, and zgrid fir landmaks
    
    alternatively, z-landmask may be derived from topo and topo mask

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int status;
  int /*k,*/m,n;
  poc_global_t      meta_gridfile, meta_maskfile, meta_topofile;
  poc_data_t<float> meta_info, meta_H_z, meta_H_f;
  grid_t fgrid,zgrid;
  char  *zlandmask;
  float *ztopo, ztopomask;
  float *ftopo, ftopomask;
  
  resize_t resize;
  
  bool have_fgrid=false, have_zgrid=false;
  bool have_ftopo=false, have_ztopo=false;
  bool have_zmask=false;

  status=poc_inq(gridfile,&meta_gridfile);
  if(status!=0) NC_TRAP_ERROR(return,status,1,"poc_inq(\"file=%s\",) error",gridfile);
  
  if(maskfile!=0) {
    status=poc_inq(maskfile,&meta_maskfile);
    if(status!=0) NC_TRAP_ERROR(return,status,1,"poc_inq(\"file=%s\",) error",maskfile);
    }
  else {
    meta_maskfile=meta_gridfile;
    maskfile=strdup(gridfile);
    }
  
  if(input.topofile!="") {
    status=poc_inq(input.topofile,&meta_topofile);
    if(status!=0) NC_TRAP_ERROR(return,status,1,"poc_inq(\"file=%s\",) error",input.topofile.c_str());
    }
  else {
    input.topofile=gridfile;
    meta_topofile=meta_gridfile;
    }
  
  if(input.f_gnames.vlon!=0) {
    status=meta_info.init(meta_gridfile,input.f_gnames.vlon);
    if(status==0) {
      status=poc_get_grid(gridfile,meta_info.info,&fgrid);
      if(status!=0) return(-1);
      status=map_projection_forward(fgrid, wprojection);
      have_fgrid=true;
      }
    else {
      have_fgrid=false;
      }
    }
  else {
    have_fgrid=false;
    }
  
  if(input.z_gnames.vlon!=0) status=meta_info.init(meta_gridfile,input.z_gnames.vlon);
  else status=-1;
  
  if(status==0) {
    status=poc_get_grid(gridfile,meta_info.info,&zgrid);
    if(status!=0) return(-1);
    status=map_projection_forward(zgrid, wprojection);
// /*------------------------------------------------------------------------------
//     debugging verification */
//     status=map_projection_backward(zgrid, wprojection.c_str());
    have_zgrid=true;
    }
  else {
    have_zgrid=false;
    }
  
  if(!have_zgrid && !have_fgrid) return(-1);
   
/*------------------------------------------------------------------------------
  check if vorticity grid bathymetry is given (it should) or not (then patch it) */
  if(input.f_gnames.vtopo==0) status=-1;
  else
    status=meta_H_f.init(meta_topofile,input.f_gnames.vtopo);
  
  if(status==0) {
    status=meta_H_f.read_data(gridfile,0,0,1);
    if(status!=0) return(-1);
    ftopo=meta_H_f.data;
    ftopomask=meta_H_f.mask;
    ztopomask=ftopomask;
    have_ztopo=false;
    have_ftopo=true;
    }
  else {
    status=meta_H_z.init(meta_topofile,input.z_gnames.vtopo);
    if(status!=0) return(-1);
    if(!have_zgrid) {
      status=poc_get_grid(input.topofile,meta_H_z.info,&zgrid);
      if(status!=0) return(-1);
      have_zgrid=true;
      }
    status=meta_H_z.read_data(input.topofile,0,0,1);
    if(status!=0) return(-1);
    ztopo=meta_H_z.data;
    ztopomask=meta_H_z.mask;
    ftopomask=ztopomask;
    have_ztopo=true;
    have_ftopo=false;
    }
  
  if(!have_zgrid) zgrid=map_f2zgrid(fgrid);

  if(!have_ztopo) {
    ztopo=map_interpolate_v2z(zgrid, fgrid, ftopo, ftopomask);
    }
  
/*------------------------------------------------------------------------------
  check if tracer grid landmask is given (it should) or not (then patch it) */
  if(input.z_gnames.vmask==0) {
    float ftag;
    if(isnan(tag)) ftag=ztopomask;
    else ftag =tag;
    printf("warning: landmask not availbale, detected from bathymetry (land if z=%f)\n",ftag);
    zlandmask=new char [zgrid.Hsize()];
    for(n=0;n<zgrid.Hsize();n++) {
      zlandmask[n]=(ztopo[n]!=ftag);
      }
    }
  else {
    poc_data_t<int8_t> meta_m;
    if(!have_zgrid) {
      status=poc_get_grid(gridfile,meta_info.info,&zgrid);
      if(status!=0) return(-1);
      have_zgrid=true;
      }
    status=meta_m.init(meta_maskfile,input.z_gnames.vmask);
    if(status!=0) return(-1);
    status=meta_m.read_data(maskfile,0,0,1);
//     if(meta_m.nlimited==2) {
//       zlandmask=(char *) meta_m.data;
//       }
//     else {
      zlandmask=new char[zgrid.Hsize()];
      for(int m=0;m<zgrid.Hsize();m++) zlandmask[m]=meta_m.data[m];
//       }
    for(int m=0;m<zgrid.Hsize();m++) if(zlandmask[m]==-1) zlandmask[m]=0;
    have_zmask=true;
    }
  
  int unmasked=occurence('\1',zlandmask,zgrid.Hsize());
     
  if(!have_fgrid) {
/*------------------------------------------------------------------------------
    eleminate masked rows and columns */
    resize=map_CheckZgrid(zgrid, zlandmask);
    status=map_resize(zgrid, resize, zlandmask);
    status=map_resize(zgrid, resize, ztopo);
//     meta_H_z.data=ztopo;
/*------------------------------------------------------------------------------
    vorticity grid not given, apply patch */
    fgrid=map_vgrid(zgrid);
    unmasked=occurence('\1',zlandmask,zgrid.Hsize());
    }
  
    
  if(!have_ftopo) {
    if( (zgrid.nx==fgrid.nx+1) && (zgrid.ny==fgrid.ny+1))
      ftopo=map_interpolate_v2z(fgrid, zgrid, ztopo, ztopomask);
    else
      ftopo=map_extrapolate_z2v(zgrid, fgrid, ztopo, ztopomask);
    }
  if(meta_H_z.data!=ztopo) delete[] ztopo;

/*------------------------------------------------------------------------------
  import f-grid (vorticity) as quadrangle mesh */
  int FixPinched=0;
  unmasked=occurence('\1',zlandmask,zgrid.Hsize());
  
  mesh.type=(fgrid.projection!=0);
  
// /*------------------------------------------------------------------------------
//   debugging verification */
//   status=map_projection_backward(fgrid, wprojection.c_str());
  
  status=fe_ZGrid2Quadrangle(zgrid, fgrid, zlandmask, ftopo, ftopomask, mesh, FixPinched, 0);
  delete[]zlandmask;
  
/*------------------------------------------------------------------------------
  save "consistent" mesh*/
  string echo;
  
  if(rootname=="") rootname="echo";

  status=fe_spherical(mesh, wprojection);
  
  echo=rootname+"-quadrangle.nei";
  status=fe_savemesh(echo.c_str(),MESH_FILE_FORMAT_TRIGRID,mesh);
        
  echo=rootname+"-quadrangle.nc";
  status=fe_savemeshNC2D_new(echo.c_str(),  mesh, 1);
    
  float *area=new float[mesh.nquadrangles];
  float mask=-1;
  for(m=0;m<mesh.nquadrangles;m++) {
    area[m]=mesh.quadrangles[m].TrueArea;
    }
  
  cdfvar_t variable=poc_variable_UG2D("area", mask, "m²", 1.f, 0.f, "cell_area","NQUADRANGLES");
  status   = cdf_createvariable(echo.c_str(), &(variable));
  status   = poc_put_UG2D(echo.c_str(), mesh, variable.id, area);
  variable.destroy();
  delete[] area;
  
  status=map_projection_backward(fgrid, wprojection);
    
  float *topography=new float[mesh.nquadrangles];
  set_grid_list(&fgrid);
  for(m=0;m<mesh.nquadrangles;m++) {
    double t,p;
    int64_t accel;
    status=fe_position(mesh, mesh.quadrangles[m], &t, &p,0);
//    status=map_interpolation(fgrid,topo,topomask,t,p,&topography[m]);
    index_interpolation(fgrid, t, p, &accel, ftopo, ftopomask, &topography[m], 0);
    if(topography[m]==ftopomask) {
      index_interpolation(fgrid, t, p, &accel, ftopo, ftopomask, &topography[m], 0);
      topography[m]=0.0;
      }
    if(isnan(topography[m])!=0) {
      topography[m]=0.0;
      }
    }
  mask=ftopomask;
  variable=poc_variable_UG2D("topography", mask, "m", 1.f, 0.f, "topography","NQUADRANGLES");
  status   = cdf_createvariable(echo.c_str(), &(variable));
  status   = poc_put_UG2D(echo.c_str(), mesh, variable.id, topography);
  variable.destroy();
  delete[] topography;
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int quadrangle_ImportStructured(const char * filename, mesh_t & mesh, metagrid_t input, Cgrid_t & Cgrid, double tag)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Development notes:

  Check :

  Note: 22/12/2012

    build the quadrangle mesh from vorticity nodes.

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
{
  int status;
  int /*k,*/m,n;
  poc_global_t      meta_gridfile, meta_maskfile;
  poc_data_t<float> meta_info, meta_H_z, meta_H_f;
  poc_data_t<int8_t> meta_m;
  grid_t fgrid, zgrid, ugrid, vgrid;
  char  *zlandmask;
  float *ztopo, ztopomask;
  float *ftopo, ftopomask;
  
  resize_t resize;
  
  bool have_fgrid=false, have_zgrid=false;
  bool have_ftopo=false, have_ztopo=false;
  
  status=poc_inq(filename,&meta_gridfile);
  if(status!=0) NC_TRAP_ERROR(return,status,1,"poc_inq(\"file=%s\",) error",filename);
  
  if(input.maskfile!="") {
    status=poc_inq(input.maskfile,&meta_maskfile);
    if(status!=0) NC_TRAP_ERROR(return,status,1,"poc_inq(\"file=%s\",) error",input.maskfile.c_str());
    }
  else {
    meta_maskfile=meta_gridfile;
    input.maskfile=filename;
    }
  
  status=meta_info.init(meta_gridfile,input.f_gnames.vlon);
  if(status==0) {
    status=poc_get_grid(filename,meta_info.info,&fgrid);
    have_fgrid=true;
    }
  else {
    have_fgrid=false;
    }
  
  status=meta_info.init(meta_gridfile,input.z_gnames.vlon);
  if(status==0) {
    status=poc_get_grid(filename,meta_info.info,&zgrid);
    have_zgrid=true;
    }
  else {
    have_zgrid=false;
    }
  
  if(!have_zgrid && !have_fgrid) return(-1);
   
/*------------------------------------------------------------------------------
  check if vorticity grid bathymetry is given (it should) or not (then patch it) */  
  status=meta_H_f.init(meta_gridfile,input.f_gnames.vtopo);
  if(status==0) {
    status=meta_H_f.read_data(filename,0,0,1);
    ftopo=meta_H_f.data;
    ftopomask=meta_H_f.mask;
    ztopomask=ftopomask;
    have_ztopo=false;
    have_ftopo=true;
    }
  else {
    status=meta_H_z.init(meta_gridfile,input.z_gnames.vtopo);
    if(status!=0) return(-1);
    status=meta_H_z.read_data(filename,0,0,1);
    ztopo=meta_H_z.data;
    ztopomask=meta_H_z.mask;
    ftopomask=ztopomask;
    have_ztopo=true;
    have_ftopo=false;
    }
  
  if(!have_zgrid) zgrid=map_f2zgrid(fgrid);

  if(!have_ztopo) {
    ztopo=map_interpolate_v2z(zgrid, fgrid, ftopo, ftopomask);
    }
  
/*------------------------------------------------------------------------------
  check if tracer grid landmask is given (it should) or not (then patch it) */
  if(input.z_gnames.vmask==0) {
    float ftag;
    if(isnan(tag)) ftag=ztopomask;
    else ftag =tag;
    zlandmask=new char [zgrid.Hsize()];
    for(n=0;n<zgrid.Hsize();n++) {
      zlandmask[n]=(ztopo[n]!=ftag);
      }
    }
  else {
    status=meta_m.init(meta_maskfile,input.z_gnames.vmask);
    if(status!=0) return(-1);
    status=meta_m.read_data(input.maskfile,0,0,1);
    if(meta_m.nlimited==2) {
      zlandmask=(char *) meta_m.data;
      }
    else {
      zlandmask=new char[zgrid.Hsize()];
      for(int m=0;m<zgrid.Hsize();m++) zlandmask[m]=meta_m.data[m];
      }
    for(int m=0;m<zgrid.Hsize();m++) if(zlandmask[m]==-1) zlandmask[m]=0;
    }
  
  if(!have_fgrid) {
/*------------------------------------------------------------------------------
    eleminate masked rows and columns */
    resize=map_CheckZgrid(zgrid, zlandmask);
    status=map_resize(zgrid, resize, zlandmask);
    meta_m.data=(signed char *) zlandmask;
    status=map_resize(zgrid, resize, ztopo);
    meta_H_z.data=ztopo;
/*------------------------------------------------------------------------------
    vorticity grid not given, apply patch */
    fgrid=map_vgrid(zgrid);
    }
  if(!have_ftopo) {
    if( (zgrid.nx==fgrid.nx+1) && (zgrid.ny==fgrid.ny+1)) {
      resize.init(zgrid,1,1,1,1);
      status=map_resize_grid(zgrid,resize);
      status=map_resize(zgrid, resize, zlandmask);
      status=map_resize(zgrid, resize, ztopo);
      meta_H_z.data=ztopo;
      }
    ftopo=map_extrapolate_z2v(zgrid, fgrid, ztopo, ztopomask);
    }

// /*------------------------------------------------------------------------------
//   import f-grid (vorticity) as quadrangle mesh */
//   int FixPinched=0;
//   status=fe_ZGrid2Quadrangle(zgrid, fgrid, zlandmask, ftopo, ftopomask, mesh, FixPinched, 0);
    
/*------------------------------------------------------------------------------
  load u grid */
  status=meta_info.init(meta_gridfile,input.u_gnames.vlon);
  if(status==0) {
    status=poc_get_grid(filename,meta_info.info,&ugrid);
    }
  else {
    ugrid=map_f2ugrid(fgrid);
    }
  
/*------------------------------------------------------------------------------
  load v grid */
  status=meta_info.init(meta_gridfile,input.v_gnames.vlon);
  if(status==0) {
    status=poc_get_grid(filename,meta_info.info,&vgrid);
    }
  else {
    vgrid=map_f2vgrid(fgrid);
    }
  
/*------------------------------------------------------------------------------
  import f-grid (vorticity) as quadrangle mesh */
  int FixPinched=1;
  status=fe_ZGrid2Quadrangle(zgrid, fgrid, zlandmask, ftopo, ftopomask, mesh, FixPinched, 0);
  
//   for(n=0;n<mesh.nvtxs;n++) {
//     float z=0;
//     for(k=0;k<mesh.vertices[n].nelmts;k++) {
//       m=mesh.vertices[n].elmts[k];
//       z+=mesh.quadrangles[m].h;
//       }
//     mesh.vertices[n].h=z/mesh.vertices[n].nelmts;
//     }
  
  status=map_minmax(&zgrid);
  
  Cgrid.z_grid=new grid_t;
  Cgrid.z_grid->duplicate(zgrid);
  
  Cgrid.f_grid=new grid_t;
  Cgrid.f_grid->duplicate(fgrid);
  
  Cgrid.u_grid=new grid_t;
  Cgrid.u_grid->duplicate(ugrid);
  
  Cgrid.v_grid=new grid_t;
  Cgrid.v_grid->duplicate(vgrid);
  
 
  Cgrid.z_incidence=aset(zgrid.Hsize(),-1);
//  set_grid_list(&fgrid);
  
#pragma omp parallel for
  for(int n=0;n<mesh.nquadrangles;n++){
    int m;
    int status;
    double lon,lat, distance;
    fe_position(mesh,mesh.quadrangles[n],&lon,&lat);
    status=map_ClosestVertex(zgrid,lon,lat,m,distance, false);
//     if(distance>0.5) {
//       status=map_ClosestVertex(zgrid,lon,lat,m,distance,true);
//       }
    if(status==-1) TRAP_ERR_EXIT(-1, "identification failed\n");
#pragma omp critical(there)
    Cgrid.z_incidence[m]=n;
    }
 
  Cgrid.f_incidence=aset(fgrid.Hsize(),-1);
#pragma omp parallel for
  for(int n=0;n<mesh.nvtxs;n++){
    int m;
    int status;
    double lon,lat, distance;
    fe_position(mesh,mesh.vertices[n],&lon,&lat);
    status=map_ClosestVertex(fgrid,lon,lat,m,distance, false);
    if(status==-1) TRAP_ERR_EXIT(-1, "identification failed\n");
    if(distance>0.1) {
      status=map_ClosestVertex(fgrid,lon,lat,m,distance,true);
      }
    Cgrid.f_incidence[m]=n;
    }

  Cgrid.u_incidence=aset(ugrid.Hsize(),-1);
  Cgrid.v_incidence=aset(vgrid.Hsize(),-1);
#pragma omp parallel for
  for(int n=0;n<mesh.nedges;n++){
    int mu,mv;
    int status;
    double lon,lat, lonu, latu, lonv, latv, distance;
    double dlon=1.e-05;
    double dlat=1.e-05;
    fe_position(mesh,mesh.edges[n],&lon,&lat);
    status=map_ClosestVertex(ugrid,lon,lat,mu,distance, false);
    lonu=ugrid.x[mu];
    latu=ugrid.y[mu];
    if(is_equal(lon,lonu,dlon) && is_equal(lat,latu,dlat)){
      Cgrid.u_incidence[mu]=n;
      continue;
      }
    status=map_ClosestVertex(vgrid,lon,lat,mv,distance, false);
    lonv=vgrid.x[mv];
    latv=vgrid.y[mv];
    if(is_equal(lon,lonv,dlon) && is_equal(lat,latv,dlat)){
      Cgrid.v_incidence[mv]=n;
      }
    else{
      status=map_ClosestVertex(ugrid,lon,lat,mu,distance, true);
      status=map_ClosestVertex(vgrid,lon,lat,mv,distance, true);
      TRAP_ERR_EXIT(8,"programming error: (%g;%g), delta u(%g;%g) v(%g;%g), precision (%g,%g)\n",lon,lat,lon-lonu,lat-latu,lon-lonv,lat-latv);
      }
    }
/*------------------------------------------------------------------------------
  save "consistent" mesh*/
  string echo;
  
  string rootname="tugo";
  
  echo=rootname+"-quadrangle.nei";
  status=fe_savemesh(echo.c_str(),MESH_FILE_FORMAT_TRIGRID,mesh);
        
  echo=rootname+"-quadrangle.nc";
  status=fe_savemeshNC2D_new(echo.c_str(),  mesh, 1);
    
  float *area=new float[mesh.nquadrangles];
  float mask=-1;
  for(m=0;m<mesh.nquadrangles;m++) {
    area[m]=mesh.quadrangles[m].TrueArea;
    }
  
  cdfvar_t variable=poc_variable_UG2D("area", mask, "m²", 1.f, 0.f, "cell_area","NQUADRANGLES");
  status   = cdf_createvariable(echo.c_str(), &(variable));
  status   = poc_put_UG2D(echo.c_str(), mesh, variable.id, area);
  variable.destroy();
  delete[] area;
  
  float *topography=new float[mesh.nquadrangles];
  set_grid_list(&fgrid);
  for(m=0;m<mesh.nquadrangles;m++) {
    double t,p;
    int64_t accel;
    status=fe_position(mesh, mesh.quadrangles[m], &t, &p,0);
//    status=map_interpolation(fgrid,topo,topomask,t,p,&topography[m]);
    index_interpolation(fgrid, t, p, &accel, ftopo, ftopomask, &topography[m], 0);
    if(topography[m]==ftopomask) {
      index_interpolation(fgrid, t, p, &accel, ftopo, ftopomask, &topography[m], 0);
      topography[m]=0.0;
      }
    if(isnan(topography[m])!=0) {
      topography[m]=0.0;
      }
    }
  mask=ftopomask;
  variable=poc_variable_UG2D("topography", mask, "m", 1.f, 0.f, "topography","NQUADRANGLES");
  status   = cdf_createvariable(echo.c_str(), &(variable));
  status   = poc_put_UG2D(echo.c_str(), mesh, variable.id, topography);
  variable.destroy();
  delete[] topography;
  return(0);
}
