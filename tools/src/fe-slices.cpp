
/*******************************************************************************

  T-UGO tools, 2006-2015

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

\brief finite element slices definitions
*/
/*----------------------------------------------------------------------------*/

#include "functions.h" //safely includes omp.h
#include "map.h"
#include "fe.h"


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T> int fe_sliceH_template(const mesh_t & mesh, int discretisation, double *depths, double level, T *buffer, T mask, T *slice)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,status;
  bool ascending;
  double factor;
  
/*------------------------------------------------------------------------------
  compute the horizontal slice array */
  
  const int nproc=omp_get_max_threads();
  double **td;
  T **tv;
  newp2D(&td,nproc,mesh.nlevels);
  newp2D(&tv,nproc,mesh.nlayers);
  
  const discretisation_t *descriptor;
  descriptor=get_descriptor_address(mesh,discretisation);

  status=check_vertical_direction(depths,NAN,descriptor->nnodes,mesh.nlevels,&ascending,&factor);
  if(status)TRAP_ERR_RETURN(status,1,"check_vertical_direction() error\n");
  
  u_int8_t interpMode;
  if(mesh.nlevels>mesh.nlayers)
    interpMode=0x10;
  else
    interpMode=0;
  
  //#pragma omp parallel for
  for(i=0;i<descriptor->nnodes;i++){
    int thI=omp_get_thread_num();
    T *vector=tv[thI];
    double *d=td[thI];
    int k,ld,lv,m;
    
    m=i;
    if(ascending){
      ld=0;
      lv=0;
      }
    else{
      ld=mesh.nlevels-1;
      lv=mesh.nlayers-1;
      }
    
    for(k=0;k<mesh.nlayers;k++){
      d[ld]=factor*depths[m];
      vector[lv]=buffer[m];
      
      if(ascending){
        ld++;
        lv++;
        }
      else{
        ld--;
        lv--;
        }
      
      m+=descriptor->nnodes;
      }
    
    if(ld==0)
      d[ld]=depths[m];
    
    status=map_interpolate1D(vector, d, mask,mesh.nlevels, level, &slice[i], interpMode);
    }
  
  deletep2D(&td,nproc);
  deletep2D(&tv,nproc);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int fe_sliceH(const mesh_t & mesh, int discretisation, double *depths, double level, float *buffer, float mask, float *slice)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=fe_sliceH_template(mesh,discretisation,depths,level,buffer,mask,slice);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int fe_sliceH(const mesh_t & mesh, int discretisation, double *depths, double level, double *buffer, double mask, double *slice)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=fe_sliceH_template(mesh,discretisation,depths,level,buffer,mask,slice);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int fe_sliceH(const mesh_t & mesh, int discretisation, double *depths, double level, double *buffer, double mask, complex<double> *slice)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  const discretisation_t *descriptor;
  descriptor=get_descriptor_address(mesh,discretisation);
  
  const int n=descriptor->nnodes;
  double *d=new double[n];
  
  result=fe_sliceH_template(mesh,discretisation,depths,level,buffer,mask,d);
  
  valcpy(slice,d,n);
  
  delete[]d;
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int fe_sliceH(const mesh_t & mesh, int discretisation, double *depths, double level, complex<double> *buffer, complex<double> mask, complex<double> *slice)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=fe_sliceH_template(mesh,discretisation,depths,level,buffer,mask,slice);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T> int fe_sliceV_template(const mesh_t & mesh, const grid_t &grid, const list_t &list, int discretisation, double *depths, double *x, double *y, int nloc, T *buffer, T mask,grid_t *slice_grid, T **slice, int verbose=1)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,k,m,status;
  bool ascending;
  double factor;

  T *bufferm=NULL,value;
  T *local=NULL;
  double *depthsm=NULL,h;
  
  range_t<double> xR,yR,zR;
  xR.init(x,nloc);
  zR.init(y,nloc);
  
/*------------------------------------------------------------------------------
  define the section grid */

  slice_grid->nx=nloc;
  slice_grid->ny=mesh.nlevels;
  /** \note For easier interpolation later :
  - the mesh_t::x, mesh_t::x and mesh_t::z parameters of \c slice_grid
  - and **slice
  are transposed, therefore the mesh_t::modeH is set to -2. */
  slice_grid->modeH=-2;
  slice_grid->nz=1;
  slice_grid->modeV=-1;

  const int nxy=slice_grid->Hsize();

  exitIfNull(
    slice_grid->x=new double [nxy]
    );
  exitIfNull(
    slice_grid->y=new double [nxy]
    );
  exitIfNull(
    slice_grid->z=new double [nxy]
    );
  
  const int localn=mesh.nlayers*nloc;
  
  exitIfNull(
    local=new T[localn]
    );
  
  const discretisation_t
    *descriptor=get_descriptor_address(mesh,discretisation),
//     *z_descriptor=get_descriptor_address(mesh,LGP1);
    *z_descriptor=descriptor;
  
  status=check_vertical_direction(depths,NAN,z_descriptor->nnodes,mesh.nlevels,&ascending,&factor);
  if(status)TRAP_ERR_RETURN(status,1,"check_vertical_direction() error\n");
  
  if(verbose){
    printf("ascending=%d,factor=%g\n",(int)ascending,factor);
    }
  
  double lon,lat;
  
/*------------------------------------------------------------------------------
  get depths and values along the section */
  int e,last=0;
  
  for(i=0;i<slice_grid->nx;i++){
    lon=x[i];
    lat=y[i];
    e=fe_detect(mesh,mesh.triangles,list,grid,lon,lat,last);
    
    for(k=0;k<mesh.nlevels;k++){
      
      m=z_descriptor->nnodes*k;
      status=fe_interpolate2D(mesh,z_descriptor->type,&depths[m],NAN,lon,lat,e,&h);
      h*=factor;
      yR<<h;
      
      if(ascending)
        m=i*slice_grid->ny+k;
      else
        m=i*slice_grid->ny+mesh.nlevels-k-1;
      
      slice_grid->x[m]=lon;
      slice_grid->z[m]=lat;
      slice_grid->y[m]=h;
      }
    
    for(k=0;k<mesh.nlayers;k++){
      m=descriptor->nnodes*k;
      status=fe_interpolate2D(mesh,discretisation,&buffer[m],mask,lon,lat,e,&value);
      
      if(ascending)
        m=i*mesh.nlayers+k;
      else
        m=i*mesh.nlayers+mesh.nlayers-k-1;
      
      local[m]=value;
      }
    
    }
  
/*---------------------------------------------------------------------*//**<h1>
  fill holes in the depths </h1>*/
  for(i=0;i<slice_grid->nx;i++){
    m=i*slice_grid->ny;
    
    for(k=0;k<mesh.nlevels && isnan(slice_grid->y[m+k]);k++);
    if(k==mesh.nlayers){
      //STDERR_BASE_LINE("fool-proofing for %d\n",i);
      slice_grid->y[m]=INFINITY;
      }
    h=slice_grid->y[m];
    
    for(k=0;k<mesh.nlevels;k++,m++){
      if(isnan(slice_grid->y[m])){
        slice_grid->y[m]=h;
        }
      else{
        h=slice_grid->y[m];
        }
      }
    
    }
  
  slice_grid->xmax=xR.max;
  slice_grid->xmin=xR.min;
  slice_grid->ymax=yR.max;
  slice_grid->ymin=yR.min;
  slice_grid->zmax=zR.max;
  slice_grid->zmin=zR.min;
  
  slice_grid->dx=(slice_grid->xmax-slice_grid->xmin)/(slice_grid->nx-1);
  slice_grid->dy=(slice_grid->ymax-slice_grid->ymin)/(slice_grid->ny-1);
  slice_grid->dz=(slice_grid->zmax-slice_grid->zmin)/(slice_grid->nz-1);
  
  if(verbose){
    slice_grid->print();
    }

  status=0;

  *slice=local;

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int fe_sliceV(const mesh_t & mesh, const grid_t &grid, const list_t &list, int discretisation, double *depths, double *x, double *y, int nloc, float *buffer, float mask,grid_t *slice_grid, float **slice, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=fe_sliceV_template(mesh,grid,list,discretisation,depths,x,y,nloc,buffer,mask,slice_grid,slice,verbose);
  
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int fe_sliceV(const mesh_t & mesh, const grid_t &grid, const list_t &list, int discretisation, double *depths, double *x, double *y, int nloc, double *buffer, double mask,grid_t *slice_grid, double **slice, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=fe_sliceV_template(mesh,grid,list,discretisation,depths,x,y,nloc,buffer,mask,slice_grid,slice,verbose);
  
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int fe_sliceV(const mesh_t & mesh, const grid_t &grid, const list_t &list, int discretisation, double *depths, double *x, double *y, int nloc, double *buffer, double mask,grid_t *slice_grid, complex<double> **slice, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  double *d;
  
  status=fe_sliceV_template(mesh,grid,list,discretisation,depths,x,y,nloc,buffer,mask,slice_grid,&d,verbose);
  
  int n=slice_grid->Hsize();
  
  *slice=new complex<double>[n];
  
  valcpy(*slice,d,n);
  
  delete[]d;
  
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int fe_sliceV(const mesh_t & mesh, const grid_t &grid, const list_t &list, int discretisation, double *depths, double *x, double *y, int nloc, complex<double> *buffer, complex<double> mask,grid_t *slice_grid, complex<double> **slice, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=fe_sliceV_template(mesh,grid,list,discretisation,depths,x,y,nloc,buffer,mask,slice_grid,slice,verbose);
  
  return status;
}
