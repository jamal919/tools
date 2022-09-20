
/**************************************************************************

  T-UGO tools, 2006-2011

  Unstructured Ocean Grid initiative

***************************************************************************/
/**
\author  Florent Lyard      LEGOS/CNRS, Toulouse, France
  E-mail: florent.lyard@legos.obs-mip.fr
\author  Yoann Le Bars      LEGOS, Toulouse, France (PhD)
\author  Laurent Roblou     LEGOS/CNRS, Toulouse, France
\author  Damien Allain      LEGOS/CNRS, Toulouse, France
\author  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
\author  Frédéric Dupont    Université de Laval à Québec, Canada
*/

#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <unistd.h>

#include "archive.h"

#include "poc-netcdf.hpp"
#include "solverlib.h"

#include "tides.h"
#include "fe.h"
#include "fe-integration.h"
#include "map.h"
#include "mass.h"

#include "filter.h"
#include "statistic.h"
#include "functions.h"
#include "maths.h"
#include "zapper.h"

#include "matrix.h"

using namespace std;

static int idx = 0;

extern matrix2x2_t  *DragMatrix;
extern zmatrix2x2_t *FricMatrix;
extern void spectral_friction01(mesh_t, parameter_t, int, int);
extern void spectral_friction02(mesh_t, parameter_t, int, int);
extern void spectral_friction03(spectrum_t spectrum,double duration, parameter_t data, int nndes, int target, double **buffer);
extern void spectral_viscosity(spectrum_t ,double ,mesh_t mesh,parameter_t , int , double **, double **);

extern int fe_LGP1_beta_gradient(mesh_t *mesh, int n, int j, double *C, double *D);

/*----------------------------------------------------------------------------*/
///Compute flux divergence
/**
\date redone 19 Jul 2011
\author Damien Allain

\param mesh grid
\param *u longitude component of the speed
\param *v latitude component of the speed
\param u_disc type of discretisation  of the speed
\param *h depth
\param *divergence
\param *Ihu longitude component of integral of fluxes for all triangles
\param *Ihv latitude component of integral of fluxes for all triangles
\param *init Default:true. If true, reset divergence to 0.

Does NCP1xLGP1 and LGP0xLGP1.
\todo other discretisations

It first computes divergence integral over the mesh triangles by integrating fluxes over the triangle's area
then integrates for each nodes using the discretisation gradient.
*/
/*----------------------------------------------------------------------------*/
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T1,typename T2,typename T3> void flux_divergence(mesh_t mesh, T1 *u, T1 *v, int u_disc, T2 *h, T3 *BIhu, T3 *Ihu,T3 *Ihv, int init=1)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int j, l,m,n; //corner or side index; indexes of edges, nodes and triangles
  triangle_t *tp;//pointer to triangle

  double lcB[3],C,D;//cosine of latitude of triangle corners, gradient components
  double NArea;// total area of cornered triangles
  T2 hB[3];//depth of triangle corners
  T1 uB[3],vB[3];//speed components of triangle corners

  //P0 FLUXES
  switch(u_disc){
  case NCP1:
    for(m=0;m<mesh.ntriangles;m++) {// for all triangles
      tp = &mesh.triangles[m];
      for(j=0;j<3;j++) {// for all corners
        n=tp->vertex[j];
        hB[j]=h[n]; // corner depth
        lcB[j]=cos(mesh.vertices[n].lat*d2r); // corner cosine
        }
      for(j=0;j<3;j++) {// for all sides
        l=tp->edges[j];
        uB[j]=u[l]; // side speed
        vB[j]=v[l]; // side speed
        }
      Ihu[m]=fe_integraleLGP1xNCP1_2D(hB,uB);
      Ihv[m]=fe_integraleLGP1xLGP1xNCP1_2D(hB,lcB,vB);
      }
    break;
  case LGP0:
    for(m=0;m<mesh.ntriangles;m++) {// for all triangles
      tp = &mesh.triangles[m];
      for(j=0;j<3;j++) {// for all corners
        n=tp->vertex[j];
        hB[j]=h[n]; // corner depth
        lcB[j]=cos(mesh.vertices[n].lat*d2r); // corner cosine
        }
      Ihu[m]=u[m]*fe_integraleLGP1_2D(hB);
      Ihv[m]=v[m]*fe_integraleLGP1xLGP1_2D(hB,lcB);
      }
    break;
  default:
    __OUT_BASE_LINE__("line %d of %s: %s speed discretisation NOT CODED YET\n",__LINE__,__FILE__,discretisation_name(u_disc));
    exit(1);
    }
  
  if(init)for(n=0;n<mesh.nvtxs;n++){
    BIhu[n]=0.;
    }
  for(n=0;n<mesh.nvtxs;n++){//for all nodes
    NArea=0.;
    for(j = 0; j < mesh.vertices[n].nelmts; j++){//for all cornered triangles
      m=fe_LGP1_beta_gradient(&mesh,n,j,&C,&D);
      BIhu[n]+=Ihu[m]*C+Ihv[m]*D;
      NArea+=mesh.triangles[m].Area;
      }
    BIhu[n]/=NArea;
    }
}


/*----------------------------------------------------------------------------*/
///for use ONLY by callocn()
/**
Do NOT call directly. Use callocn() instead.
\date created 13 May 2011
\author Damien Allain
*/
/*----------------------------------------------------------------------------*/
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> void callocn_(size_t n,T **p,...)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  va_list ap;
  va_start(ap,p);
  do{
    *p=new T[n];
    p=va_arg (ap, T **);
    }while(p);
  va_end (ap);
}
/*----------------------------------------------------------------------------*/
///allocate arrays
/**
\date created 16 May 2011
\author Damien Allain

\param n number of elements to allocate your arrays to
\param "p,..." list of pointers to the arrays you want to allocate

For example :
 \code callocn(n,&a1,&a2)} \endcode
allocates arrays \a a1 and \a a2 to \a n elements.
This macro calls safely callocn_().
*/
/*----------------------------------------------------------------------------*/
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
//*** You _MUST_ have NULL as the last argument to callocn_ ***
#define callocn(n,p,...) callocn_(n,p , ##__VA_ARGS__,NULL)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// info:/cpp-4.5/Pragmas
#pragma GCC poison callocn_

/*----------------------------------------------------------------------------*/
/** Finds the closest mesh node to the given coordinates
\date created 29 Jul 2011
\author Damien Allain

\param mesh
\param lat target latitude in degrees
\param lon target longitude in degrees
\returns the index of the closest node. Negative if error (if mesh.nvtxs<=0).
*/
/*----------------------------------------------------------------------------*/
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int find_closest(mesh_t mesh, double lat, double lon)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int n1=-1,n;//indexes of nodes : last found and current.
  vertex_t *cv;//pointer to current vertex
///Uses sinuses of angles rather than 3D distances.
  double a1=M_2_PI,a;//angles to nodes : last found and current.
///Scans the whole mesh in an optimised way.
  vector3_t M0=math_polar2cartesian(lon*r2d,lat*r2d),M;//coordinates of target point
  
  for(n=0;n<mesh.nvtxs;n++){//for all nodes
    cv=&mesh.vertices[n];
    a=square(M0 & math_polar2cartesian(cv->lon*r2d,cv->lat*r2d));
    if(a1>a){
      n1=n;
      }
    }
  
  return n1;
}


/*----------------------------------------------------------------------------*/
/** Finds the linear path connecting a pair
\date documented on 28 Jul 2011
\author Damien Allain

\param mesh
\param connexions indexes of extremities
\param *path array of indexes of nodes. The extremities are respectively the first and last elements of the array.
\returns the number of nodes in the path

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
/** It iterates by seeking the closest neighbour */
      d=fe_distance(mesh,m,n2);
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


/*----------------------------------------------------------------------------*/
/// Detects bathymetric incoherencies
/**
by measuring the ratio of
the divergence of the flux
over
the time-derivative of the elevation.
Does LGP1xLGP1 (untested), NCP1xLGP1 (untested) and LGP0xLGP1.
\date created 26 Apr 2011
\date last reviewed 28 Jun 2011
\author Damien Allain

\todo other discretisations

\param mesh grid
\param AnalysisList list of tidal frequencies to do the analysis for
\param *output_path path to output directory
\param *path path to directory of netcdf input
*/
/*----------------------------------------------------------------------------*/
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void spectral_depth(mesh_t mesh,spectrum_t AnalysisList, parameter_t data2D, char *output_path, char *path)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  //Indexes
  int l, m, n, w;//of edges, triangles, nodes, waves
  int i, j, k;//of corners, cornered triangles, sides or packed values
  triangle_t *tp;//pointer to triangle
  //Arrays
  const int BL=mesh.nedges;//buffer length : edges are always the most numerous
  double *dB1,*dB2,*dB3,*dB4,mask;//large buffers
  const int n3as=mesh.ntriangles,nEds=mesh.nedges,nNs=mesh.nvtxs;
  double *h,hB[3],lcB[3];//depth and its buffer,cos(lat) buffer
  float *fB;//read buffer of depth
  dcomplex *z,*jwz,*Ijwz,*Ijwz0,*u,*v;//elevation, its time derivative and the integral thereof, its backup, speed componenents
  dcomplex uB[3],vB[3];//buffers of speed componenents
  //results
  double C,D;//Beta gradient components
  dcomplex jw,*ha;
  //
  int status;
  tidal_wave wave;
  /** As it does not seam to work properly on sequential modes, it shows you why. */
  enum{sequential,spectral} mode;
  //netCDF
  char ncfilename[1024];
  int ncid,nvar,iv;
  char **ncvarname,*varNamePart,cB1[128];//,,char buffer
  //discretisation
  int z_discretisation,u_discretisation,*z_OR_u,paire;
  int nu,nz;
  //projection
  hypermatrix_t *LGP1massmatrix;
  hyperzmatrix_t MBIu;//matrix of divergences of speeds for all nodes
  /** It also tries to work out something proper for sequential modes. */
  /** \todo correct effect of non-linearity in the continuation equation */
  struct {dcomplex *z,*u,*v;} tides[AnalysisList.n];
  int kZ0,kM2,kM4;//node index, wave frequency indicator
  /** It also tries to constrain the depth to be a real (or the imaginary part to be 0).
  It works in that the imaginary part of the depth is closer and closer to 0.
  It fails in that the real part of the depth is farther and farther from the truth (or even sense). */
  hyperzmatrix_t ctMBIu,sMBIu,rMBIu;//conjugate transpose, product and regularising form thereof
  int sMBIuNPacked=-1,rMBIuNPacked=-1,*diagonalPacked=new int[4*nNs];//number of non-zero elements of [sr]MBIu.packed, index of values influenced by the constraint factor in rMBIu
  double rMBIuFactor,rMBIuFactors[]={0.,.1,1,10,100,1E3,1E4,1E6,1E7};
  int *clamped,nclamped=-1,npacked=-1;//indexes of nodes at open boundaries, their number, number of non-zero elements of hyper(z)matrix_t::packed
  #define USE_clamped2 0 //does not seam to have any effect
  #if USE_clamped2
  //int clamped2[]={19962,19070,18425,12123,12334,11937,11539,12146},nclamped2=sizeof(clamped2)/sizeof(*clamped2);
  int *clamped2,nclamped2;
  #endif
  /** It is extremely important to set up the right bathymetry at the open boundary! */
  #define ASSERT_clamped_jwz 0
  #define ASSERT_clamped_ha 0
  dcomplex *dcB1,*dcB2,*dcB3,*dcB4;
  bool force_duplicate=false;
  int verbose=0;
  bool debug=false;
  
  //### ALLOCATION ###
  fB=new float[nNs];
  callocn(BL,&dB1,&dB2,&dB3,&dB4);//double
  callocn(2*nNs,&ha);//dcomplex
  callocn(BL,&z,&jwz,&Ijwz,&u,&v,&dcB1,&dcB2,&dcB3,&dcB4,&Ijwz0);//dcomplex
  callocn(nNs,&h);//double
  
  #if USE_clamped2
  clamped2=new int[nNs];
  /// It clamps between hard-coded coordinates! TODO: make flexible.
  nclamped2=find_path(mesh,
    paire_t(
      find_closest(mesh,26.1,51.25),
      find_closest(mesh,27.9,51.75)
      ),
    clamped2
    );
  #endif
  
  fprintf(stderr,"line %d of %s: spectral_depth working.\n",__LINE__,__FILE__);
  
  kZ0=kM2=kM4=-1;
  for(w = 0; w < AnalysisList.n; w++) {
    wave = AnalysisList.waves[w];
    if(strcmp("Z0",wave.name)==0) kZ0=w;
    if(strcmp("M2",wave.name)==0) kM2=w;
    if(strcmp("M4",wave.name)==0) kM4=w;
    tides[w].z=NULL;
    tides[w].u=NULL;
    tides[w].v=NULL;
    }
  
  for(w = 0; w < AnalysisList.n; w++) {
    //### DATA LOADING ###
    wave = AnalysisList.waves[w];
    //Find and open netcdf file
    mode=spectral;
    do{
      //prefer spectral
      sprintf(ncfilename, "%s/%s.spectral.nc" ,path,  wave.name);
      status=nc_open(ncfilename, NC_NOWRITE, &ncid);if(!status)break;
      sprintf(ncfilename, "%s/%s.sequential.nc" ,path,  wave.name);
      status=nc_open(ncfilename, NC_NOWRITE, &ncid);
      mode=sequential;
    }while(0);
    if(status){
      printf("line %d of %s: could not open %s : skipping.\n",__LINE__,__FILE__,ncfilename);
      continue;
      }
    //else
    printf("line %d of %s: Variables of %s are :\n",__LINE__,__FILE__,ncfilename);
    //FIND DISCRETISATION
    z_discretisation=u_discretisation=UNSET;
    varNamePart=NULL;
/**----------------------------------------------------------------------------
    obsolete call, will be suppressed */
    cdf_getvariable(ncid,&nvar,&ncvarname);
    for(iv=0;iv<nvar;iv++){
      if(!strchr("aG",ncvarname[iv][0]))//neither amplitude nor phase
        continue;
      if(!strncmp("_eta_",&ncvarname[iv][1])){
        z_OR_u=&z_discretisation;
        varNamePart=&ncvarname[iv][6];
        }
      else if(!strncmp("_u_",&ncvarname[iv][1]) ||
              !strncmp("_v_",&ncvarname[iv][1])){
        z_OR_u=&u_discretisation;
        varNamePart=&ncvarname[iv][4];
        }
      else // other variables
        continue;
      /*conflicting discretisation*/
#define discFromNameTest(X) if(!strcmp(varNamePart,#X)){ \
        if(*z_OR_u==UNSET)\
          *z_OR_u=X;\
        else if(*z_OR_u!=X)\
          break;\
        continue;\
        }
      discFromNameTest(LGP0)
      discFromNameTest(LGP1)
      discFromNameTest(NCP1)
#undef discFromNameTest
      //discretisation not coded for
      *z_OR_u=UNSET;
      iv=nvar;
      break;
      }
    if(iv<nvar){
      printf("line %d of %s: conflicting discretisation in %s : skipping.\n",__LINE__,__FILE__,ncfilename);
      continue;
      }
    if(z_discretisation==UNSET || u_discretisation==UNSET){
      if(varNamePart)
        printf("line %d of %s: %s discretisation not coded for : skipping %s\n",__LINE__,__FILE__,varNamePart,ncfilename);
      else
        printf("line %d of %s: no known discretisation found in %s : skipping.\n",__LINE__,__FILE__,ncfilename);
      continue;
      }
    for(iv=0;iv<nvar;iv++)
      free(ncvarname[iv]);
    free(ncvarname);
    //CHECK DISCRETISATION and set numbers of elements
    switch(u_discretisation){
    case LGP0:
      nu=n3as;break;
    case LGP1:
      nu=nNs;break;
    case NCP1:
      nu=nEds;break;
    default:
      printf("line %d of %s: %s speed discretisation not coded for : skipping %s\n",__LINE__,__FILE__,discretisation_name(u_discretisation),ncfilename);
      continue;
      }
    switch(z_discretisation){
    case LGP1:
      nz=nNs;break;
    default:
      printf("line %d of %s: %s elevation discretisation not coded for : skipping %s\n",__LINE__,__FILE__,discretisation_name(z_discretisation),ncfilename);
      continue;
      }
    printf("line %d of %s: %sx%s discretisation.\n",__LINE__,__FILE__,discretisation_name(u_discretisation),discretisation_name(z_discretisation));
    //READ NETCDF FILE
    varNamePart=cB1;
    do{
      status=poc_getbuffer(ncid,"bathymetry",fB);if(status)break;
      for(n=0;n<nNs;n++){
        #define FILE_BATHYMETRY 1
        #if FILE_BATHYMETRY // get bathymetry from netcdf
        h[n]=fB[n];
        #else //get bathymetry from mesh
        h[n]=mesh.vertices[n].h;
        #endif
        dB3[n]=fB[n]/mesh.vertices[n].h;
        }
      sprintf(varNamePart,"eta_%s",discretisation_name(z_discretisation));
      status=poc_get_UG3D(ncid,-1,varNamePart,z,dB1,dB2);
      if(status)break;
      sprintf(varNamePart,"u_%s",discretisation_name(u_discretisation));
      status=poc_get_UG3D(ncid,-1,varNamePart,u,dB1,dB2);
      if(status)break;
      sprintf(varNamePart,"v_%s",discretisation_name(u_discretisation));
      status=poc_get_UG3D(ncid,-1,varNamePart,v,dB1,dB2);
      if(status)break;
    }while(0);
    nc_close(ncid);
    if(status){
      fprintf(stderr,"Error while reading %s : %s\n",ncfilename,nc_strerror(status));
      continue;
      }
    printf("line %d of %s: Successfully read %s\n",__LINE__,__FILE__,ncfilename);
    
    sprintf(ncfilename, "%s/%s.hSearch.%2.2d.nc" ,output_path,  wave.name, idx);
    unlink(ncfilename);
    archiving_UGdummy2D(ncfilename, mesh, "h_h0", "m", dB3, mask, LGP1);
    tides[w].z=copy(z,nz);
    tides[w].u=copy(u,nu);
    tides[w].v=copy(v,nu);
    
    //### CALCULATIONS ###
    if(z_discretisation==LGP1 && nclamped<0){
      status=discretisation_init(&mesh, LGP1);//required by mesh.LGP1descriptor.massmatrix
      nclamped=0;
      for(n = 0; n < mesh .nvtxs; n++) {// for all nodes
        if(mesh.vertices[n].code<0){//on open boundary
          nclamped++;
          }
        }
  #if USE_clamped2
      nclamped+=nclamped2;
  #endif
      //point to the same ordering, copy t field, etc...
      //requires a prior call to discretisation_init(&mesh, LGP1) to initialise mesh.LGP1descriptor.massmatrix
      LGP1massmatrix=&mesh.LGP1descriptor.massmatrix;
      npacked=LGP1massmatrix->ordering->pointer[nNs];
      i=0;
      if(nclamped){
        //initialise clamped
        clamped=new int[2*nclamped];
        for(n = 0; n < nNs; n++) {// for all nodes
          if(mesh.vertices[n].code<0){//open boundary
            clamped[i]=n;
            clamped[i+nclamped]=n+nNs;
            i++;
            }
          }
        }
  #if USE_clamped2
      for(l=0;l<nclamped2;l++){
        n=clamped2[l];
        clamped[i]=n;
        clamped[i+nclamped]=n+nNs;
        i++;
        }
  #endif
      }
    if(!MBIu.packed){
      //4 lines inspired from call of init_packed_LGP1xLGP1_0 and context in WE2D_LGP0xLGP1_SpInitialise in cefmo-LGP0xLGP1.cpp
      //See also initLGP1 and dMassMatrix and their calls in discretisation_init in discretisation-initialise.cpp
      //MBIu.ordering=LGP1massmatrix->ordering;//DO NOT SHARE AN ordering BETWEEN AN hypermatrix_t AND AN hyperzmatrix_t !!!
      MBIu.ordering=copy(LGP1massmatrix->ordering);
      MBIu.hbw=mesh.hbw;
      MBIu.packed=new dcomplex[npacked];
      }
    #if 0 //xscan works properly
    for(n=0;n<nu;n++){
      dB1[n]=abs(u[n]);
      dB2[n]=arg(u[n]);
      dB3[n]=abs(v[n]);
      dB4[n]=arg(v[n]);
      }
    archiving_UGdummy2D(ncfilename, mesh, "a_u", dB1,u_discretisation);
    archiving_UGdummy2D(ncfilename, mesh, "G_u", dB2,u_discretisation);
    archiving_UGdummy2D(ncfilename, mesh, "a_v", dB3,u_discretisation);
    archiving_UGdummy2D(ncfilename, mesh, "G_v", dB4,u_discretisation);
    #else //UGLY PATCH !!!
    if(u_discretisation==NCP1){
      for(n=0;n<nNs;n++){
        if(mesh.vertices[n].edges!=0)
          l=mesh.vertices[n].edges[0];
        else
          l=mesh.triangles[mesh.vertices[n].elmts[0]].edges[0];
        dB1[n]=abs(u[l]);
        dB2[n]=arg(u[l]);
        dB3[n]=abs(v[l]);
        dB4[n]=arg(v[l]);
        }
      archiving_UGdummy2D(ncfilename, mesh, "a_u", "m/s",     dB1, mask, LGP1);
      archiving_UGdummy2D(ncfilename, mesh, "G_u", "degrees", dB2, mask, LGP1);
      archiving_UGdummy2D(ncfilename, mesh, "a_v", "m/s",     dB3, mask, LGP1);
      archiving_UGdummy2D(ncfilename, mesh, "G_v", "degrees", dB4, mask, LGP1);
      }
    #endif
    //BATHYMETRY
    jw=dcomplex(0,wave.omega*dph2rps);// degrees/h -> rad/s conversion
    if(z_discretisation==LGP1){
      for(n=0;n<nNs;n++){//for all nodes
        jwz[n]=jw*z[n];
        dB1[n]=abs(jwz[n]);
        dB2[n]=arg(jwz[n]);
        }
      archiving_UGdummy2D(ncfilename, mesh, "a_jwz", "m/s",     dB1, mask, LGP1);
      archiving_UGdummy2D(ncfilename, mesh, "G_jwz", "degrees", dB2, mask, LGP1);
      matrix_operation(*LGP1massmatrix,jwz,Ijwz);
      for(k=0;k<npacked;k++){
        MBIu.packed[k]=0.;
        }
      for(n=0;n<nNs;n++){//for all nodes
        for(j = 0; j < mesh.vertices[n].nelmts; j++){//for all cornered triangles
          m=fe_LGP1_beta_gradient(&mesh,n,j,&C,&D);
          tp = &mesh.triangles[m];
          for(i=0;i<3;i++) {// for all corners and edges
            l=tp->vertex[i];//vertex index
            lcB[i]=cos(mesh.vertices[l].lat*d2r); // corner cosine
            switch(u_discretisation){
            case LGP0:break;
            case NCP1:
              l=tp->edges[i];//edge index
            case LGP1:
              uB[i]=u[l];
              vB[i]=v[l];
              break;
            default:
              printf("line %d of %s: %s speed discretisation BEING CODED : skipping %s\n",__LINE__,__FILE__,discretisation_name(u_discretisation),ncfilename);
              goto EOWave;
              }
            }//EO for all corners and edges
          for(i=0;i<3;i++){//for all corners (also neighbour nodes)
            l=tp->vertex[i];//vertex index
            for(k=MBIu.ordering->pointer[n];k<MBIu.ordering->pointer[n+1];k++) {//look for value index in packed
              if(MBIu.ordering->incidence[k]==l)
                break;
              }
            if(k==MBIu.ordering->pointer[n+1]){//value not found
              printf("line %d of %s: sparse matrix error at %d,%d : skipping %s\n",__LINE__,__FILE__,n,l,ncfilename);
              goto EOWave;
              }
            switch(u_discretisation){
            case LGP0:
              MBIu.packed[k]+=C*u[m]*fe_integraleLGP1_2D(i);
              MBIu.packed[k]+=D*v[m]*fe_integraleLGP1xLGP1_2D(lcB,i);
              break;
            case LGP1:
              MBIu.packed[k]+=C*fe_integraleLGP1xLGP1_2D(uB,i);
              MBIu.packed[k]+=D*fe_integraleLGP1xLGP1xLGP1_2D(lcB,vB,i);
              break;
            case NCP1:
              MBIu.packed[k]+=C*fe_integraleLGP1xNCP1_2D(i,uB);
              MBIu.packed[k]+=D*fe_integraleLGP1xLGP1xNCP1_2D(lcB,i,vB);
              break;
            default:
              printf("line %d of %s: %s speed discretisation BEING CODED : skipping %s\n",__LINE__,__FILE__,discretisation_name(u_discretisation),ncfilename);
              goto EOWave;
              }
            }//EO for all corners (also neighbour nodes)
          }//EO for all triangles
        }//EO for all nodes
      varNamePart=cB1;
      memcpy(Ijwz0,Ijwz,n*sizeof(dcomplex));//backup for later use
      ctMBIu=copy(MBIu);//backup because it will be messed up
      //sparsePPM(MBIu,67);
      for(l=1;l<=2;l++){//for all Modes
        if(l>1){
          MBIu.destroy();
	  MBIu=copy(ctMBIu);//restore original matrix
          memcpy(Ijwz,Ijwz0,n*sizeof(dcomplex));//restore
          }
        printf("line %d of %s: %g%+gi\n",__LINE__,__FILE__,real(Ijwz[9000]),imag(Ijwz[9000]));
        switch(l){
        case 1:
          break;
        case 2://FAKE SPECTRAL and CORRECT FOR NON-LINEARITIES
          printf("line %d of %s: FAKING SPECTRAL for %s\n",__LINE__,__FILE__,ncfilename);
          matrix_operation(MBIu,h,Ijwz);
          for(n=0;n<nNs;n++){//for all nodes
            dB1[n]=real(Ijwz[n]);
            dB2[n]=imag(Ijwz[n]);
            dcB1[n]=jwz[n];//backup
            }
          for(i=0;i<nclamped;i++){
            n=clamped[i];
            dB1[n]=real(jwz[n])/diagonal;
            dB2[n]=imag(jwz[n])/diagonal;
            }
          /*clamped[0]@nclamped LGP1massmatrix->packed[0]@40 dB1[0]@4 dB2[0]@4 */
          status=LinearSystem_initialise(LGP1massmatrix,clamped,nclamped, SOLVER_ID_UMFPACK, force_duplicate, verbose, debug);
          printf("line %d of %s: LinearSystem_initialise said %d for %s\n",__LINE__,__FILE__,status,ncfilename);
          status=LinearSystem_solve(*LGP1massmatrix,dB1);
          printf("line %d of %s: LinearSystem_solve said %d for %s\n",__LINE__,__FILE__,status,ncfilename);
          status=LinearSystem_solve(*LGP1massmatrix,dB2);
          printf("line %d of %s: LinearSystem_solve said %d for %s\n",__LINE__,__FILE__,status,ncfilename);
          #if ASSERT_clamped_jwz
          for(i=0;i<nclamped;i++){
            n=clamped[i];
            printf("line %d of %s: %d %d : %g%+gi%+g%+gi=%g%+gi",__LINE__,__FILE__,i,n,real(dcB1[n]),imag(dcB1[n]),real(-jwz[n]),imag(-jwz[n]),real(dcB1[n]-jwz[n]),imag(dcB1[n]-jwz[n]));
            if(abs(dcB1[n]-jwz[n])/abs(dcB1[n]+jwz[n])>.1){//UGLY PATCH !!!
              printf(" ***CORRECTING***");
              jwz[n]=dcB1[n];
              }
            printf("\n");
            }
          #endif
          for(n=0;n<nNs;n++){//for all nodes
            jwz[n]=dcomplex(dB1[n],dB2[n]);
            dB1[n]=abs(jwz[n]);
            dB2[n]=arg(jwz[n]);
            dcB1[n]-=jwz[n];
            dB3[n]=abs(dcB1[n]);
            dB4[n]=arg(dcB1[n]);
            }
          archiving_UGdummy2D(ncfilename, mesh, "a_jwz_spectral" , "todo", dB1, mask, LGP1);
          archiving_UGdummy2D(ncfilename, mesh, "G_jwz_spectral" , "todo", dB2, mask, LGP1);
          archiving_UGdummy2D(ncfilename, mesh, "a_jwz_nonLinear" , "todo", dB3, mask, LGP1);
          archiving_UGdummy2D(ncfilename, mesh, "G_jwz_nonLinear" , "todo", dB4, mask, LGP1);
          if(w!=kM2)
            continue;//not coded for other frequencies yet
          //Z0xM2
          for(n=0;n<nNs;n++){//for all nodes
            dcB1[n]=0.;
            }
          if(tides[kZ0].z){
            flux_divergence(mesh,u,v,u_discretisation,tides[kZ0].z,dcB2,dcB3,dcB4);
            for(n=0;n<nNs;n++){//for all nodes
              dB1[n]=abs(dcB2[n]);
              dB2[n]=arg(dcB2[n]);
              dcB1[n]+=dcB2[n]*.5;
              }
            archiving_UGdummy2D(ncfilename, mesh, "a_div_uM2_zZ0", "todo" , dB1, mask, LGP1);
            archiving_UGdummy2D(ncfilename, mesh, "G_div_uM2_zZ0", "todo" , dB2, mask, LGP1);
            flux_divergence(mesh,tides[kZ0].u,tides[kZ0].v,u_discretisation,z,dcB2,dcB3,dcB4);
            for(n=0;n<nNs;n++){//for all nodes
              dB1[n]=abs(dcB2[n]);
              dB2[n]=arg(dcB2[n]);
              dcB1[n]+=dcB2[n]*.5;
              }
            archiving_UGdummy2D(ncfilename, mesh, "a_div_uZ0_zM2", "todo" , dB1, mask, LGP1);
            archiving_UGdummy2D(ncfilename, mesh, "G_div_uZ0_zM2", "todo" , dB2, mask, LGP1);
            }
          //M2xM2
          flux_divergence(mesh,u,v,u_discretisation,z,dcB2,dcB3,dcB4);
          for(n=0;n<nNs;n++){//for all nodes
            dB1[n]=abs(dcB2[n]);
            dB2[n]=arg(dcB2[n]);
            //dcB1[n]+=dcB2[n]*.5;
            }
          archiving_UGdummy2D(ncfilename, mesh, "a_div_uM2_zM2", "todo" , dB1, mask, LGP1);
          archiving_UGdummy2D(ncfilename, mesh, "G_div_uM2_zM2", "todo" , dB2, mask, LGP1);
          //M2xM4
          if(tides[kM4].z){
            flux_divergence(mesh,u,v,u_discretisation,tides[kM4].z,dcB2,dcB3,dcB4);
            for(n=0;n<nNs;n++){//for all nodes
              dB1[n]=abs(dcB2[n]);
              dB2[n]=arg(dcB2[n]);
              dcB1[n]+=dcB2[n]*.5;
              }
            archiving_UGdummy2D(ncfilename, mesh, "a_div_uM2_zM4", "todo" , dB1, mask, LGP1);
            archiving_UGdummy2D(ncfilename, mesh, "G_div_uM2_zM4", "todo" , dB2, mask, LGP1);
            flux_divergence(mesh,tides[kM4].u,tides[kM4].v,u_discretisation,z,dcB2,dcB3,dcB4);
            for(n=0;n<nNs;n++){//for all nodes
              dB1[n]=abs(dcB2[n]);
              dB2[n]=arg(dcB2[n]);
              dcB1[n]+=dcB2[n]*.5;
              }
            archiving_UGdummy2D(ncfilename, mesh, "a_div_uM4_zM2", "todo" , dB1, mask, LGP1);
            archiving_UGdummy2D(ncfilename, mesh, "G_div_uM4_zM2", "todo" , dB2, mask, LGP1);
            }
          for(n=0;n<nNs;n++){//for all nodes
            dB3[n]=abs(dcB1[n]);
            dB4[n]=arg(dcB1[n]);
            }
          archiving_UGdummy2D(ncfilename, mesh, "a_jwz_correction", "todo" , dB3, mask, LGP1);
          archiving_UGdummy2D(ncfilename, mesh, "G_jwz_correction", "todo" , dB4, mask, LGP1);
          if(0)for(n=0;n<nNs;n++){//for all nodes
            Ijwz[n];
            }
          //continue;
          break;//EO FAKE SPECTRAL and CORRECT FOR NON-LINEARITIES
          }
        memcpy(ha,Ijwz,n*sizeof(dcomplex));
        for(n=0;n<nNs;n++){//for all nodes
          dB1[n]=abs(Ijwz[n]);
          dB2[n]=arg(Ijwz[n]);
          }
        printf("line %d of %s: %g%+gi\n",__LINE__,__FILE__,real(Ijwz[9000]),imag(Ijwz[9000]));
        sprintf(varNamePart,"a_Ijwz_%d",l);
        archiving_UGdummy2D(ncfilename, mesh, varNamePart, "m" , dB1, mask, LGP1);
        sprintf(varNamePart,"G_Ijwz_%d",l);
        archiving_UGdummy2D(ncfilename, mesh, varNamePart, "degrees" , dB2, mask, LGP1);
        for(i=0;i<nclamped;i++){
          n=clamped[i];
          ha[n]=h[n]/diagonal;
          #if ASSERT_clamped_ha
          dcB1[i]=ha[n];
          #endif
          }
//         status=LinearSystem_initialise(&MBIu,clamped,nclamped, SOLVER_ID_UMFPACK, force_duplicate, verbose, debug);
        status=LinearSystem_initialise(&MBIu,clamped,nclamped, SOLVER_ID_UMFPACK);
        printf("line %d of %s: LinearSystem_initialise said %d for %s\n",__LINE__,__FILE__,status,ncfilename);
        status=LinearSystem_solve(MBIu,ha);
        printf("line %d of %s: LinearSystem_solve said %d for %s\n",__LINE__,__FILE__,status,ncfilename);
        #if ASSERT_clamped_ha
        for(i=0;i<nclamped;i++){
          n=clamped[i];
          printf("line %d of %s: %d %d : %g%+gi%+g%+gi=%g%+gi\n",__LINE__,__FILE__,i,n,real(dcB1[i]),imag(dcB1[i]),real(-ha[n]),imag(-ha[n]),real(dcB1[i]-ha[n]),imag(dcB1[i]-ha[n]));
          }
        #endif
        for(n=0;n<nNs;n++){//for all nodes
          dB1[n]=abs(ha[n]);
          dB2[n]=arg(ha[n]);
          dB3[n]=real(ha[n]);
          dB4[n]=imag(ha[n]);
          }
        sprintf(varNamePart,"a_ha_%d",l);
        archiving_UGdummy2D(ncfilename, mesh, varNamePart, "m", dB1, mask, LGP1);
        sprintf(varNamePart,"G_ha_%d",l);
        archiving_UGdummy2D(ncfilename, mesh, varNamePart, "m" , dB2, mask, LGP1);
        sprintf(varNamePart,"Re_ha_%d",l);
        archiving_UGdummy2D(ncfilename, mesh, varNamePart, "m" , dB3, mask, LGP1);
        sprintf(varNamePart,"Im_ha_%d",l);
        archiving_UGdummy2D(ncfilename, mesh, varNamePart, "m" , dB4, mask, LGP1);
        if(mode==spectral)
          break;
        }//EO for all Modes
      if(mode==spectral)
        break;
      goto EOnc;//SKIPPING CONSTRAINTS
      MBIu.destroy();MBIu=copy(ctMBIu);//restore original matrix
      //initialise rMBIu
    #define USE_Hermitian 1
    #if USE_Hermitian
//       status=LinearSystem_initialise(&MBIu,clamped,nclamped, SOLVER_ID_UMFPACK, force_duplicate, verbose, debug);
      status=LinearSystem_initialise(&MBIu,clamped,nclamped, SOLVER_ID_UMFPACK);
      printf("line %d of %s: LinearSystem_initialise said %d for %s\n",__LINE__,__FILE__,status,ncfilename);
      matrix_transpose(&ctMBIu);
      for(k=0;k<npacked;k++){
        ctMBIu.packed[k]=conj(ctMBIu.packed[k]);
        }
      printf("line %d of %s: ",__LINE__+1,__FILE__);
      matrix_product(ctMBIu,MBIu,&sMBIu);
    #else
      sMBIu.destroy();sMBIu=copy(MBIu);
    #endif
      rMBIu.ordering=new ordering_t;
      rMBIu.ordering->nrows=2*nNs;
      if(sMBIu.ordering->max_cardinal>=0)
        rMBIu.ordering->max_cardinal=sMBIu.ordering->max_cardinal+1;
      rMBIu.ordering->cardinal=new int[rMBIu.neq()];
      rMBIu.ordering->pointer=new int[rMBIu.neq()+1];
      sMBIuNPacked=sMBIu.ordering->pointer[nNs];
      rMBIuNPacked=rMBIu.ordering->pointer[rMBIu.neq()]=2*(sMBIuNPacked+nNs);
      rMBIu.ordering->incidence=new int[rMBIuNPacked];
      rMBIu.packed=new dcomplex[rMBIuNPacked];
      for(n=0;n<nNs;n++){
        l=sMBIu.ordering->cardinal[n];
        rMBIu.ordering->cardinal[n]=l+1;
        j=sMBIu.ordering->pointer[n];
        k=rMBIu.ordering->pointer[n]=j+n;
        memcpy(&rMBIu.packed[k],&sMBIu.packed[j],
          l*sizeof(dcomplex));
        rMBIu.packed[k+l]=0.;
        memcpy(&rMBIu.ordering->incidence[k],&sMBIu.ordering->incidence[j],
          l*sizeof(int));
        i=diagonalPacked[2*n+1]=k+l;
        rMBIu.ordering->incidence[i]=nNs+n;
        for(;k<i && rMBIu.ordering->incidence[k]!=n;k++);
        if(k==i){//value not found
          printf("line %d of %s: sparse matrix error at %d : skipping %s\n",__LINE__,__FILE__,n,ncfilename);
          goto EOWave;
          }
        diagonalPacked[2*n]=k;
        }
      for(n=0;n<nNs;n++){
        l=sMBIu.ordering->cardinal[n];
        rMBIu.ordering->cardinal[n+nNs]=l+1;
        j=sMBIu.ordering->pointer[n];
        k=rMBIu.ordering->pointer[n+nNs]=sMBIuNPacked+nNs+j+n;
        rMBIu.packed[k]=0.;
        rMBIu.ordering->incidence[k]=n;
        diagonalPacked[2*(n+nNs)+1]=k;
        i=k+l+1;
        for(k++;k<i;k++,j++){
          rMBIu.packed[k]=conj(sMBIu.packed[j]);
          rMBIu.ordering->incidence[k]=sMBIu.ordering->incidence[j]+nNs;
          }
        for(k=rMBIu.ordering->pointer[n+nNs]+1;k<i && rMBIu.ordering->incidence[k]!=n+nNs;k++);
        if(k==i){//value not found
          printf("line %d of %s: sparse matrix error at %d : skipping %s\n",__LINE__,__FILE__,n,ncfilename);
          goto EOWave;
          }
        diagonalPacked[2*(n+nNs)]=k;
        }
      rMBIu.ordering->gethbw();
      printf("line %d of %s: made rMBIu for %s\n",__LINE__,__FILE__,ncfilename);
      ctMBIu.destroy();ctMBIu=copy(rMBIu);//backup because it will be messed up
      //EO initialise rMBIu
      if(0)for(l=0;l<sizeof(rMBIuFactors)/sizeof(rMBIuFactor);l++){//constraints
        printf("line %d of %s: %p\n",__LINE__,__FILE__,rMBIu.ordering);
        if(l){//rMBIu is messed up, so
          rMBIu.destroy();rMBIu=copy(ctMBIu);//restore it
          }
    #if USE_Hermitian
        status=matrix_operation(MBIu,Ijwz,ha);
    #else
        memcpy(ha,Ijwz,n*sizeof(dcomplex));
    #endif
        printf("line %d of %s: matrix_operation said %d for %s\n",__LINE__,__FILE__,status,ncfilename);
        for(i=0;i<nclamped;i++){
          n=clamped[i];
          ha[n]=h[n]/diagonal;
          }
        for(n=0;n<nNs;n++){
          ha[nNs+n]=conj(ha[n]);
          }
    #if USE_Hermitian
        rMBIuFactor=pow(rMBIuFactors[l],2);
    #else
        rMBIuFactor=rMBIuFactors[l];
    #endif
        for(i=0;i<2*nNs;i++){
          k=diagonalPacked[2*i];
          rMBIu.packed[k]+=rMBIuFactor;
          k=diagonalPacked[2*i+1];
          rMBIu.packed[k]-=rMBIuFactor;
          }
        #if 0
        if(0 && l==8)
          sparsePPM(rMBIu,67);
        #endif
    #if USE_Hermitian
//         status=LinearSystem_initialise(&rMBIu,clamped,0, SOLVER_ID_UMFPACK, force_duplicate, verbose, debug); 
        status=LinearSystem_initialise(&rMBIu,clamped,0, SOLVER_ID_UMFPACK); 
    #else
        status=LinearSystem_initialise(&rMBIu,clamped,2*nclamped, SOLVER_ID_UMFPACK, force_duplicate, verbose, debug);
    #endif
        printf("line %d of %s: LinearSystem_initialise said %d for %s\n",__LINE__,__FILE__,status,ncfilename);
        status=LinearSystem_solve(rMBIu,ha);
        printf("line %d of %s: LinearSystem_solve said %d for %s\n",__LINE__,__FILE__,status,ncfilename);
        for(n=0;n<nNs;n++){//for all nodes
          dB1[n]=abs(ha[n]);
          dB2[n]=arg(ha[n]);
          dB3[n]=real(ha[n]);
          dB4[n]=imag(ha[n]);
          }
        sprintf(varNamePart,"a_ha_%d_",l);
        archiving_UGdummy2D(ncfilename, mesh, varNamePart, "m", dB1, mask, LGP1);
        sprintf(varNamePart,"G_ha_%d_",l);
        archiving_UGdummy2D(ncfilename, mesh, varNamePart, "m", dB2, mask, LGP1);
        sprintf(varNamePart,"Re_ha_%d_",l);
        archiving_UGdummy2D(ncfilename, mesh, varNamePart, "m", dB3, mask, LGP1);
        sprintf(varNamePart,"Im_ha_%d_",l);
        archiving_UGdummy2D(ncfilename, mesh, varNamePart, "m", dB4, mask, LGP1);
        for(n=0;n<nNs;n++){//for all nodes
          dB1[n]=abs(ha[nNs+n]);
          dB2[n]=arg(ha[nNs+n]);
          dB3[n]=real(ha[nNs+n]);
          dB4[n]=imag(ha[nNs+n]);
          }
        sprintf(varNamePart,"a_ha_%d_x",l);
        archiving_UGdummy2D(ncfilename, mesh, varNamePart, "m", dB1, mask, LGP1);
        sprintf(varNamePart,"G_ha_%d_x",l);
        archiving_UGdummy2D(ncfilename, mesh, varNamePart, "m", dB2, mask, LGP1);
        sprintf(varNamePart,"Re_ha_%d_x",l);
        archiving_UGdummy2D(ncfilename, mesh, varNamePart, "m", dB3, mask, LGP1);
        sprintf(varNamePart,"Im_ha_%d_x",l);
        archiving_UGdummy2D(ncfilename, mesh, varNamePart, "m", dB4, mask, LGP1);
        }//EO constraints
      }//EO z_discretisation==LGP1
EOnc:;
    printf("line %d of %s: wrote to %s\n",__LINE__,__FILE__,ncfilename);
EOWave:;
    rMBIu.destroy();
    sMBIu.destroy();
    ctMBIu.destroy();
    MBIu.destroy();
    }//EO for all waves
  fprintf(stderr,"line %d of %s: spectral_depth worked!\n",__LINE__,__FILE__);
}



/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

#define assertMesh() plotMeshCode(&mesh,__FILE__,__LINE__,58.55,58.7,23.5,23.75)
int plotMeshCode(mesh_t *m,const char *sourcename,int lineNo,double lon_min=-INFINITY,double lon_max=INFINITY,double lat_min=-INFINITY,double lat_max=INFINITY)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**
This function produces a data file of the code values of the mesh.
It is so that I can plot them and work out what fe_read_boundarycode does!

In gnuplot, just do:
  set view map
  splot "2869.dat" every :::1::1 w p ps 0.2 lc palette z
  splot "2869.dat" every :::0::0 w p ps 0.2 lc palette z
*/
{
  FILE *f;//file
  char p[1024];//file path
  int i;//index of nodes or edges
  vertex_t *vp;//pointer to node
  edge_t *ep;//pointer to edge
  #if 0 //mesh range properties are all 0 !
  if(lon_min<m->lon_min)lon_min=m->lon_min;
  if(lat_min<m->lat_min)lat_min=m->lat_min;
  if(lon_max>m->lon_max)lon_max=m->lon_max;
  if(lat_max>m->lat_max)lat_max=m->lat_max;
  #endif
  ///open file
#define plotMeshCode_gplt 0 //produce script. DO NOT! It is so slow!
#if plotMeshCode_gplt
  sprintf(p,"%d.dat",lineNo);
#else
  sprintf(p,"%d.dat",lineNo);
#endif
  printf("%s\n",p);
  f=fopen(p,"w");
  if(!f)return 1;
#if plotMeshCode_gplt
  ///header
  fprintf(f,"%s",
   "unset label\n"
   "l(c,x,y)=sprintf(\"set label '%d' at %g,%g tc lt %d point lt %d\",c,x,y,p,p)\n"
   "set macros;m='eval l'\n"
   );
  ///nodes
  fprintf(f,"\np=1\n");//P1
#endif
  for(i = 0; i < m->nvtxs; i++) {// for all nodes
    vp=&m->vertices[i];
#if plotMeshCode_gplt
    if(vp->lon<lon_min ||
       vp->lon>lon_max ||
       vp->lat<lat_min ||
       vp->lat>lat_max)continue;
    fprintf(f,"@m(%d,%g,%g)\n",vp->code,vp->lon,vp->lat);
#else
    fprintf(f,"%g %g %d\n",vp->lon,vp->lat,vp->code);
#endif
    }
  ///edges
#if plotMeshCode_gplt
  fprintf(f,"\np=7\n");//P2
#else
  fprintf(f,"\n");
#endif
  for(i = 0; i < m->nedges; i++) {// for all edges
    ep=&m->edges[i];
#if plotMeshCode_gplt
    if(ep->lon<lon_min ||
       ep->lon>lon_max ||
       ep->lat<lat_min ||
       ep->lat>lat_max)continue;
    fprintf(f,"@m(%d,%g,%g)\n",ep->code,ep->lon,ep->lat);
#else
    fprintf(f,"%g %g %d\n",ep->lon,ep->lat,ep->code);
#endif
    }
  ///footer
#if plotMeshCode_gplt
  fprintf(f,"\nplot [%g:%g] [%g:%g] sqrt(-1) not\n",lon_min,lon_max,lat_min,lat_max);
#endif
  fclose(f);
  return 0;
}

