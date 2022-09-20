
/*******************************************************************************

  T-UGO tools, 2006-2017

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

<!-- USE FIND AND REPLACE TO EDIT THIS LINE, SO THAT print_help IS ALSO UPDATED -->
\brief 
*/
/*----------------------------------------------------------------------------*/

#define MAIN_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>

#include "tools-structures.h"

#include "geo.h"
#include "fe.h"
#include "fe-proto.h"
#include "map.h"
#include "polygones.h"
#include "grd.h"
#include "archive.h"
#include "functions.h"
#include "statistic.h"
#include "archive.h"

#include "zapper.h"     /*  rutin.h contains common utility routines  */


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_disassemby(mesh_t mesh, string rootname, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 
  extract base components of mesh
 
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
{
  int status;
  vector<plg_t> polygons;
  string filename;
  
  if(rootname=="") rootname="mesh-disassembly";

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  get mesh limits as polygons
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  status=fe_limits2poly(mesh,polygons,0,debug);
  
  filename=rootname+".plg";
  status=plg_save(filename.c_str(), PLG_FORMAT_SCAN, polygons);
  
  
  filename=rootname+".nod";
  status=fe_savenodes(filename.c_str(), NODE_FILE_FORMAT_TRIGRID, mesh);

  return(status);
  
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_ChkLimits(mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int   l,n,status=0;
  
  for(l=0;l<mesh.nlimits;l++) {
    if(mesh.limits[l].nedges==4) {
      n=mesh.limits[l].vertex[0];
      double L=mesh.limits[l].length/1000.;
      if(L > 75.) {
        printf("limit %4d: 4 nodes island lon=%9.3lf, lat=%9.3lf length=%9.3lf km\n", l, mesh.vertices[n].lon, mesh.vertices[n].lat, L);
        status++;
        }
      }
    }
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int mesh_wavelength_error(mesh_t & mesh, float *resolution, float *topo, double T, string rootname)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int i,m,n;
  string varname;
  double g=9.81;
  
  fe_integrale_init();
  
  status= discretisation_init(&mesh, LGP0);
  status= discretisation_init(&mesh, LGP1);
 
  float *zLGP1=new float[mesh.LGP1descriptor.nnodes], mask=1.e+10;
  float *Kerror_LGP0=new float[mesh.ntriangles];
  float *Gerror_LGP0=new float[mesh.ntriangles];
  float *gAdjustment_LGP0=new float[mesh.ntriangles];

  for(m=0;m<mesh.ntriangles;m++) {
    double H=0;
    for(i=0;i<3;i++) {
      n=mesh.triangles[m].vertex[i];
      H+=topo[n]/3.;
      }
    H=max(10.,H);
    double g=9.81;
    double c=sqrt(g*H);
    double L=c*T;
    double k=2*M_PI/L;
    double dx=resolution[m];
//     if(fabs(k*dx)>1.0) {
//       printf("trouble \n");
//       }
    double r=sqrt((cos(k*dx)+1)/2.);
    double Kdiscrete=2*asin(r*k*dx/2.)/dx;
    Kerror_LGP0[m]=Kdiscrete/k-1;
    Gerror_LGP0[m]=360.0*Kerror_LGP0[m];
    
    double s=r*2.*sin(k*dx/2.)/dx/k;
    gAdjustment_LGP0[m]=g*s*s;
    
    if(isnan(Kerror_LGP0[m])) {
      printf("trouble \n");
      }
    }
  
  status=fe_projection( mesh, Kerror_LGP0, LGP0, zLGP1, LGP1);
  varname="Kerror_LGP1_"+rootname;
  status=archiving_UGdummy2D("mesh-wavelength-error.nc", mesh, varname.c_str(), "dimensionless", zLGP1, mask, LGP1);
  
  status=fe_projection( mesh, Gerror_LGP0, LGP0, zLGP1, LGP1);
  varname="Gerror_LGP1_"+rootname;
  status=archiving_UGdummy2D("mesh-wavelength-error.nc", mesh, varname.c_str(), "degrees", zLGP1, mask, LGP1);
  
  status=fe_projection( mesh, gAdjustment_LGP0, LGP0, zLGP1, LGP1);
  varname="gAdjustment_LGP1_"+rootname;
  status=archiving_UGdummy2D("mesh-wavelength-error.nc", mesh, varname.c_str(), "ms-2", zLGP1, mask, LGP1);
  
  for(m=0;m<mesh.ntriangles;m++) {
    double H=0;
    for(i=0;i<3;i++) {
      n=mesh.triangles[m].vertex[i];
      H+=topo[n]/3.;
      }
    H=max(10.,H);
    double c=sqrt(g*H);
    double L=c*T;
    double k=2*M_PI/L;
    double dx=resolution[m];
    double r=sqrt((cos(k*dx)+4.*cos(k*dx/2)+1)/6.);
    double Kdiscrete=2*asin(r*k*dx/2.)/dx;
    Kerror_LGP0[m]=Kdiscrete/k-1;
    Gerror_LGP0[m]=360.0*Kerror_LGP0[m];
    
    double s=r*dx*k/sin(k*dx/2.)/2;
    gAdjustment_LGP0[m]=g*s*s;
    
    if(isnan(Kerror_LGP0[m])) {
      printf("trouble \n");
      }
    if(s*s>1) {
      printf("trouble \n");
      }
    }
  
  status=fe_projection( mesh, Kerror_LGP0, LGP0, zLGP1, LGP1);
  varname="Kerror_LGP2_"+rootname;
  status=archiving_UGdummy2D("mesh-wavelength-error.nc", mesh, varname.c_str(), "dimensionless", zLGP1, mask, LGP1);
  
  status=fe_projection( mesh, Gerror_LGP0, LGP0, zLGP1, LGP1);
  varname="Gerror_LGP2_"+rootname;
  status=archiving_UGdummy2D("mesh-wavelength-error.nc", mesh, varname.c_str(), "degrees", zLGP1, mask, LGP1);
  
  status=fe_projection( mesh, gAdjustment_LGP0, LGP0, zLGP1, LGP1);
  varname="gAdjustment_LGP2_"+rootname;
  status=archiving_UGdummy2D("mesh-wavelength-error.nc", mesh, varname.c_str(), "ms-2", zLGP1, mask, LGP1);
  
  for(m=0;m<mesh.ntriangles;m++) {
    for(i=0;i<3;i++) {
      n=mesh.triangles[m].vertex[i];
      if(zLGP1[n]/g>1.01) {
        printf("%d %d : %f %f \n",m,n,gAdjustment_LGP0[m]/g,zLGP1[n]/g);
        printf("trouble \n");
        }
      }
    }

  
  delete[] Kerror_LGP0;
  delete[] Gerror_LGP0;
  delete[] gAdjustment_LGP0;
  delete[] zLGP1;

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int mesh_wavelength_error(mesh_t & mesh, string filename)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  float *resolution_LGP0, *topo;
  double T;
  
  printf("#################################################################\n");
  printf("compute mesh discretisation error\n");
//   string filename="/home/models/FES2014-beta-v5/data/bathymetry/etopo-31-meanlevel/topo-LGP1-1.s2r";
  
  topo=new float[mesh.nvtxs];
  
  status=quoddy_loadr1(filename.c_str(), mesh.nvtxs, topo);
  if(status!=0) return(-1);
  
  printf("-----------------------------------------------------------------\n");
  printf("compute mesh equivalent LGP0-resolution\n");
  status=mesh_resolution(mesh, resolution_LGP0);
  if(status!=0) return(-1);

  printf("-----------------------------------------------------------------\n");
  printf("K1 error\n");
  T=24*3600;
  status=mesh_wavelength_error(mesh, resolution_LGP0, topo, T, "K1");
  printf("-----------------------------------------------------------------\n");
  printf("M2 error\n");
  T=12.5*3600;
  status=mesh_wavelength_error(mesh, resolution_LGP0, topo, T, "M2");
  printf("-----------------------------------------------------------------\n");
  printf("M4 error\n");
  T=6.25*3600;
  status=mesh_wavelength_error(mesh, resolution_LGP0, topo, T, "M4");
  printf("-----------------------------------------------------------------\n");
  printf("M6 error\n");
  T=12.5*3600/3.;
  status=mesh_wavelength_error(mesh, resolution_LGP0, topo, T, "M6");
  
  delete[] resolution_LGP0;
  delete[] topo;
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int mesh_bathymetry_error(mesh_t & mesh, float *resolution, float *topo, float *delta, double T, string rootname)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int i,m,n;
  string varname;
  
  fe_integrale_init();
  
  status= discretisation_init(&mesh, LGP0);
  status= discretisation_init(&mesh, LGP1);
 
  float *zLGP1=new float[mesh.LGP1descriptor.nnodes], mask=1.e+10;
  float *Kerror_LGP0=new float[mesh.ntriangles];
  float *Gerror_LGP0=new float[mesh.ntriangles];
  float *gAdjustment_LGP0=new float[mesh.ntriangles];

  for(m=0;m<mesh.ntriangles;m++) {
    double H=0,d=0;
    for(i=0;i<3;i++) {
      n=mesh.triangles[m].vertex[i];
      H+=topo[n]/3.;
      d+=delta[n]/3.;
      }
    H=max(10.,H);
    double g=9.81;
    double c=sqrt(g*H);
    double L=c*T;
    double k=2*M_PI/L;
    
//     double c_d=sqrt(g*(H+d));
//     double L_d=c*T;
    double k_d=2*M_PI/L;
    Kerror_LGP0[m]=k_d/k-1;
    Gerror_LGP0[m]=360.0*Kerror_LGP0[m];
    
    if(isnan(Kerror_LGP0[m])) {
      printf("trouble \n");
      }
    }
  
  status=fe_projection( mesh, Kerror_LGP0, LGP0, zLGP1, LGP1);
  varname="Kerror_LGP1_"+rootname;
  status=archiving_UGdummy2D("mesh_bathymetry-error.nc", mesh, varname.c_str(), "dimensionless", zLGP1, mask, LGP1);
  
  status=fe_projection( mesh, Gerror_LGP0, LGP0, zLGP1, LGP1);
  varname="Gerror_LGP1_"+rootname;
  status=archiving_UGdummy2D("mesh_bathymetry-error.nc", mesh, varname.c_str(), "degrees", zLGP1, mask, LGP1);
    
  delete[] Kerror_LGP0;
  delete[] Gerror_LGP0;
  delete[] gAdjustment_LGP0;
  delete[] zLGP1;

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int mesh_bathymetry_error(mesh_t & mesh, string filename)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  float *resolution_LGP0, *topo, *delta=0;
  double T;
  
  printf("#################################################################\n");
  printf("compute mesh bathymetry error\n");
//   string filename="/home/models/FES2014-beta-v5/data/bathymetry/etopo-37-meanlevel/topo-LGP1-1.s2r";
//   string deltafile="/home/models/FES2014-beta-v5/data/bathymetry/etopo-37-meanlevel/topo-LGP1-1.s2r";
  
  topo=new float[mesh.nvtxs];
  
  status=quoddy_loadr1(filename.c_str(), mesh.nvtxs, topo);
  if(status!=0) return(-1);
  
  printf("-----------------------------------------------------------------\n");
  printf("compute mesh equivalent LGP0-resolution\n");
  status=mesh_resolution(mesh, resolution_LGP0);
  if(status!=0) return(-1);

  printf("-----------------------------------------------------------------\n");
  printf("K1 error\n");
  T=24*3600;
  status=mesh_bathymetry_error(mesh, resolution_LGP0, topo, delta, T, "K1");
  printf("-----------------------------------------------------------------\n");
  printf("M2 error\n");
  T=12.5*3600;
  status=mesh_bathymetry_error(mesh, resolution_LGP0, topo, delta, T, "M2");
  printf("-----------------------------------------------------------------\n");
  printf("M4 error\n");
  T=6.25*3600;
  status=mesh_bathymetry_error(mesh, resolution_LGP0, topo, delta, T, "M4");
  printf("-----------------------------------------------------------------\n");
  printf("M6 error\n");
  T=12.5*3600/3.;
  status=mesh_bathymetry_error(mesh, resolution_LGP0, topo, delta, T, "M6");
  
  delete[] resolution_LGP0;
  delete[] topo;
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int n;
  int channels=0,geometry=0,safety=0,renum=1,reshape=1,nghmax=7;
  int delaunay=0;
  char *keyword,*zone=NULL;
  char *meshfile=NULL,*output=NULL,*poly=NULL;
  string rootname;
  mesh_t mesh,refined,internal,external;
  double dmax=0;
  bool disassembly=false, debug=false;;

  fct_echo( argc, argv);

  n=1;
  while (n < argc) {
    keyword=strdup(argv[n]);
    switch (keyword[0]) {
      case '-':
        if(strcmp(keyword,"--no_renum")==0) {
          renum=0;
          n++;
          continue;
          }
        if(strcmp(keyword,"--no_reshape")==0) {
          reshape=0;
          n++;
          continue;
          }
        if(strcmp(keyword,"--debug")==0) {
          debug=true;
          n++;
          continue;
          }
        if(strcmp(keyword,"--disassembly")==0) {
          disassembly=true;
          n++;
          continue;
          }
        if(strcmp(keyword,"--delaunay")==0) {
          delaunay=1;
          n++;
          continue;
          }
        if(strcmp(keyword,"--nghmax")==0) {
          sscanf(argv[n+1],"%d",&nghmax);
          n++;
          n++;
          continue;
          }
        switch (keyword[1]) {
        case 'm' :
/*-----------------------------------------------------------------------------
          UG mesh filename (input)*/
          meshfile= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'c' :
/*-----------------------------------------------------------------------------
          check channels*/
          channels=1;
          n++;
          break;

        case 'd' :
/*-----------------------------------------------------------------------------
          do not refine if size already smaller than dmax*/
          sscanf(argv[n+1],"%lf",&dmax);
          n++;
          n++;
          break;

        case 'g' :
/*-----------------------------------------------------------------------------
          check geometry*/
          geometry=1;
          n++;
          break;

        case 'p' :
/*-----------------------------------------------------------------------------
          limit action to polygons interior*/
          poly= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'o' :
          output= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'r' :
          rootname=argv[n+1];
          n++;
          n++;
          break;

        case 's' :
          safety=1;
          n++;
          break;

        case 'z' :
/*-----------------------------------------------------------------------------
          limit action to rectuganler zone interior*/
          zone= strdup(argv[n+1]);
          n++;
          n++;
          break;

        default:
          STDOUT_BASE_LINE("unknown option %s\n",keyword);
          exit(-1);
        }
        break;

      default:
          n++;
        break;
      }
      free(keyword);
    }

  if(meshfile == NULL) {
    printf("no mesh file specified; abort...\n");
//     print_help(argv[0]);
    wexit(-1);
    }

  printf("#################################################################\n");
  printf("load mesh and construct related tables: %s\n",meshfile);

  status=fe_readmesh(meshfile,MESH_FILE_FORMAT_TRIGRID,&mesh);
  if(status!=0) {
    printf("unable to read the original mesh in %s\n",meshfile);
//     print_help(argv[0]);
    wexit(-1);
    }

  status=fe_list(&mesh);
  if(status!=0) {
    printf("unable to build the element list from the original mesh\n");
//     print_help(argv[0]);
    wexit(-1);
    }

  printf("#Number of elements   : %d\n",mesh.ntriangles);
  printf("#Number of vertices   : %d\n",mesh.nvtxs);

  status= fe_bandwidth(&mesh);
  printf("#Half-bandwidth       : %d\n",mesh.hbw);

  status= fe_edgetable(&mesh,0,0);
  if(status!=0) {
    printf("unable to build the edge list from the original mesh\n");
//     print_help(argv[0]);
    wexit(-1);
    }
  printf("#Number of edges      : %d\n",mesh.nedges);

  printf("#################################################################\n");
  printf("reconstruct boundary codes and limits table\n");
  status= fe_codetable(mesh,0,1,0);
  if(status!=0) {
    printf("unable to rebuild the limits table and codes of the original mesh\n");
//     print_help(argv[0]);
    wexit(-1);
    }
  printf("#Number of limits     : %d\n",mesh.nlimits);
  
  if(disassembly) {
    status=fe_disassemby(mesh, rootname, debug);
    }

  status=fe_ChkLimits(mesh);

#if 0
  {
  string filename="/home/models/FES2014-beta-v5/data/bathymetry/etopo-31-meanlevel/topo-LGP1-1.s2r";
  status=mesh_wavelength_error(mesh, filename);
  }
#endif

  status=mesh_resolution(mesh);
  
  STDOUT_BASE_LINE("end of mesh-diagnostics ... \n");
  exit(0);
}
