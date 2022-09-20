
/***************************************************************************
T-UGO tools, 2006-2009

  Unstructured Ocean Grid initiative
  
E-mail: florent.lyard@legos.obs-mip.fr
***************************************************************************/
/**
\file

*/
/*----------------------------------------------------------------------------*/

#include "version-macros.def" //for VERSION and REVISION

#include <stdio.h>
#include <string.h>
#include <algorithm>

#include <errno.h>

#include "tools-define.h"
#include "tools-structures.h"

#include "fe.h"
#include "poc-time.h"
#include "list.h"
#include "map.h"
#include "map.def"
#include "grd.h"
#include "geo.h"
#include "filter.h"
#include "functions.h"
#include "tides.h"
#include "cefmo.h"
#include "netcdf-proto.h"
#include "statistic.h"
#include "xyz.h"


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int SpArchive_Get3D(const char *filename, int paire, int frame, tide3D_t state)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**-----------------------------------------------------------------------

  Spectral solution archive

 -----------------------------------------------------------------------*/
{
  int i,j,k,l,m,n,status;
  int ncid,varid;
  int zdim,udim;
  int z_discretisation, u_discretisation;
  int wave_id;
  double *buffer[2];
  char varname_a[1024], varname_G[1024];

  status=paire_discretisation_id(paire, &z_discretisation, &u_discretisation);
  if(status!=0) return(-1);
  
  const char *UNAME=discretisation_name(u_discretisation);
  const char *ZNAME=discretisation_name(z_discretisation);

/**----------------------------------------------------------------------------
  elevation */
//   sprintf(varname_a,"a_z3D_%s", ZNAME);
//   sprintf(varname_G,"G_z3D_%s", ZNAME);
//   sprintf(varname_a,"a_z3D");
//   sprintf(varname_G,"G_z3D");
//   status=poc_get_UG4D(filename, frame, (const char *) varname_a, (const char *) varname_G, state.z);
//   if(status!=0) goto error;
  
/**----------------------------------------------------------------------------
  currents */
//   sprintf(varname_a,"a_u_%s", UNAME);
//   sprintf(varname_G,"G_u_%s", UNAME);
  sprintf(varname_a,"a_u3D");
  sprintf(varname_G,"G_u3D");
  status=poc_get_UG4D(filename, frame, (const char *) varname_a, (const char *) varname_G, state.u);
  if(status!=0) goto error;
//   sprintf(varname_a,"a_v_%s", UNAME);
//   sprintf(varname_G,"G_v_%s", UNAME);
  sprintf(varname_a,"a_v3D");
  sprintf(varname_G,"G_v3D");
  status=poc_get_UG4D(filename, frame, (const char *) varname_a, (const char *) varname_G, state.v);
  if(status!=0) goto error;

  return(0);

error:
  return(-1);

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int fe_Zaxis(const char *dataname, mesh_t & mesh, discretisation_t & descriptor, int discretisation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int dim,n,length,fmt,ndim=1;
  int status,frame=NOFRAME;
  int i,j,k,m,m1,m2;
  cdfgbl_t gblinfo;
  int input_id;
  int varid=-1;
  double *z;
  char varname[64];
  
  const char *NAME=discretisation_name(discretisation);

  sprintf(varname,"z-%s",NAME);

  
  switch(discretisation) {
    case LGP0:
//      descriptor=&(mesh.LGP0descriptor);
      varid=cdf_identify(dataname,"z-LGP0");
      break;

    case LGP1:
//       descriptor=&(mesh.LGP1descriptor);
      varid=cdf_identify(dataname,"z-LGP1");
      break;

    case D_LGP1:
//       descriptor=&(mesh.DGP1descriptor);
      varid=cdf_identify(dataname,"z-DGP1");
      break;

    case LGP2:
//       descriptor=&(mesh.LGP2descriptor);
      varid=cdf_identify(dataname,"z-LGP2");
      break;

    case NCP1:
//       descriptor=&(mesh.NCP1descriptor);
      break;

    default:
      return(-1);
      break;
    }
      
  varid=cdf_identify(dataname,varname);
  
  z=new double[descriptor.nnodes*mesh.nlevels];
  if(varid!=-1) {
    if(descriptor.nodes==0) descriptor.nodes=new node_t[descriptor.nnodes];
    status=poc_get_UG3D(dataname, frame, varname, z);
//    status=poc_get_UG4D(dataname, frame, varid, z);
    for(n=0;n<descriptor.nnodes;n++) {
      descriptor.nodes[n].zlevels=new double[mesh.nlevels];
      for(k=0;k<mesh.nlevels;k++) {
        descriptor.nodes[n].zlevels[k]=z[k*descriptor.nnodes+n];
        }
      descriptor.nodes[n].nlevels=mesh.nlevels;
      }
    }
  delete[] z;
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int crozet(const char *filename, const char *polygone, int first, int nframes)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,m,n,status;
  int id,ntargets;
  mesh_t mesh;
  char *targets;
  cefmo_t cefmo;
  float **Vvmax, **Hvmax, *bathy;
  ellipse_t e;
//  int frame, nframes=10, first=24-nframes+1;
//  int frame, nframes=5, first=9-nframes+1;
  int frame;
  statistic_t *s, global;
  double *values;
  int count;
  tide3D_t state;

  int zdim,udim;
  int z_discretisation, u_discretisation;
  
  discretisation_t z_descriptor;
  discretisation_t u_descriptor;

  status=fe_readgeometry(filename, &mesh);
  if(status!=0) return(-1);
  
  cefmo.paire=LGP0xLGP1;

  mesh.nlevels=mesh.nlayers+1;

  status=paire_dimension(mesh,cefmo.paire,&zdim,&udim);
  status=paire_discretisation_id(cefmo.paire, &z_discretisation, &u_discretisation);

  id=discretisation_from_name("LGP1");
  status=fe_readdiscretisation(filename, &mesh, 0, id);
  
  id=discretisation_from_name("LGP0");
  status=fe_readdiscretisation(filename, &mesh, 0, id);
  
  z_descriptor=get_descriptor(mesh,z_discretisation);
  u_descriptor=get_descriptor(mesh,u_discretisation);

  status=fe_Zaxis(filename, mesh, u_descriptor, u_discretisation);
  
  state.allocate(udim,0,mesh.nlayers,0);

  targets=new char[mesh.ntriangles];
  for(m=0;m<mesh.ntriangles;m++) targets[m]=1;
  
  int *selected=new int[mesh.nedges];
  int nselected;
  
  for (n=0;n<mesh.nedges;n++) {
    selected[n]=1;
    }
    
  if(polygone!=0) {
  printf("#################################################################\n");
    printf("select edges from polygons: %s\n",polygone);
    nselected=fe_selectedges_01(mesh, (char *) polygone, selected);
    }
  else {
    nselected=mesh.nedges;
    }
  if(nselected==0) return(-1);
    
  ntargets=0;
  for(m=0;m<mesh.ntriangles;m++) {
    targets[m]=0;
    for(int i=0;i<3;i++) {
      n=mesh.triangles[m].edges[i];
      if(selected[n]==1) {
        targets[m]=1;
        break;
        }
      }
    if(targets[m]==1) ntargets++;
    }
    
  Vvmax=new float *[nframes+1];
  Hvmax=new float *[nframes+1];
  bathy=new float[udim];
  for(k=0;k<nframes+1;k++) {
    Vvmax[k]=new float[udim];
    Hvmax[k]=new float[udim];
    }

/**-----------------------------------------------------------------------------
  compute ellipse major axis for each iteration */
  for(frame=0;frame<nframes;frame++) {
    status=SpArchive_Get3D(filename, cefmo.paire, frame+first,  state);
    for(m=0;m<mesh.ntriangles;m++) {
      if(targets[m]==0) continue;
      float *buf=new float[mesh.nlayers];
      double vvmax=0,hvmax;
      for(k=0;k<mesh.nlayers;k++) {
        e=ellipse_parameter(state.u[m][k],state.v[m][k]);
        buf[k]=100.0*e.a;
        }
      for(k=0;k<mesh.nlayers;k++) {
        double depth=u_descriptor.nodes[m].zlevels[k];
        if(depth<-300.) {
          if(vvmax<buf[k]) {
            vvmax=buf[k];
            hvmax=depth;
            }
          }
        }
      Vvmax[frame][m]=vvmax;
      Hvmax[frame][m]=hvmax;
      }
    }
  
  s=new statistic_t[nframes];
  
/**-----------------------------------------------------------------------------
  compute statistics on ellipse major axis */
  for(frame=0;frame<nframes;frame++) {
//    status=SpArchive_Get3D(filename, cefmo.paire, frame+20,  state);
    double vvmax=0;
    count=0;
    values=new double[ntargets];
    for(m=0;m<mesh.ntriangles;m++) {
      if(targets[m]==0) continue;
      vvmax=MAX(vvmax,Vvmax[frame][m]);
      values[count]=Vvmax[frame][m];
      count++;
      }
    printf("iteration %d :",frame+first);
    s[frame]=get_statistics(values,1.e+10,count);
    delete[] values;
    }

  
/**-----------------------------------------------------------------------------
  compute overall max on ellipse major axis */
  values=new double[nframes*ntargets];
  count=0;
  for(m=0;m<mesh.ntriangles;m++) {
    bathy[m]=u_descriptor.nodes[m].zlevels[0];
    Vvmax[nframes][m]=0.0;
    if(targets[m]==0) continue;
    
    for(frame=0;frame<nframes;frame++) {
      values[count]=Vvmax[frame][m];
      if(Vvmax[nframes][m]<Vvmax[frame][m]) {
        Vvmax[nframes][m]=Vvmax[frame][m];
        Hvmax[nframes][m]=Hvmax[frame][m];
        }
      count++;
      }
    Vvmax[nframes][m]*=2.0;
    }
  global=get_statistics(values,1.e+10,count);
  
  grid_t grid=get_zonegrid("crozet");
  status=map_completegridaxis(&grid,1);
  
  int *elts=fe_scan_elements(mesh,grid,0,1);
  
  float mask=0, *SGbuf=new float[grid.Hsize()];
  
  status=fe_map(mesh, Vvmax[nframes], LGP0, grid, elts, SGbuf, mask);
  
  status=save_SGXY("vmax.nc", grid, SGbuf, mask, "Vvmax", "m/s", "Vvmax", 1);
  status=xyz_save("Vvmax.xyz", grid, SGbuf, mask);
  
  status=fe_map(mesh, Hvmax[nframes], LGP0, grid, elts, SGbuf, mask);
  
  status=save_SGXY("vmax.nc", grid, SGbuf, mask, "Hvmax", "m", "Hvmax", 0);
  status=xyz_save("Hvmax.xyz", grid, SGbuf, mask);
  
  status=fe_map(mesh, bathy, LGP0, grid, elts, SGbuf, mask);
  
  status=save_SGXY("vmax.nc", grid, SGbuf, mask, "bathymetry", "m", "bathymetry", 0);
  status=xyz_save("bathymetry.xyz", grid, SGbuf, mask);
  
  return(0);

error:
  return(-1);

}

