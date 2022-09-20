
/**************************************************************************

  T-UGO tools, 2006-2011

  Unstructured Ocean Grid initiative

***************************************************************************/
/** \file

\author  Florent Lyard      LEGOS/CNRS, Toulouse, France. florent.lyard@legos.obs-mip.fr
\author  Laurent Roblou     LEGOS/CNRS, Toulouse, France
\author  Damien Allain      LEGOS/CNRS, Toulouse, France
\author  Yves Soufflet      LEGOS/CNRS, Toulouse, France
\author  Clément Mayet      LEGOS, Toulouse, France (PhD)
\author  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
\author  Frédéric Dupont    Université de Laval à Québec, Canada

\brief finite element poc-netcdf definitions
*/
/*----------------------------------------------------------------------------*/

#include <config.h>

#include <stdio.h>
#include <string.h>

#include "tools-structures.h"

#include "fe.def"
#include "fe.h"

#include "map.h"
#include "rutin.h"
#include "geo.h"
#include "poc-time.h"
#include "netcdf-proto.h"
#include "poc-netcdf.def"
#include "functions.h"
#include "parallel.h"
// #include "fe-proto.h"

#define RANK_connectivity 2

// static int gCPU_ID=0,gCPU_MASTER=0;
static char *rootname=".",*output_path=".";


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  mesh_t output_mesh(mesh_t local, int ArchiveExtent)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  mesh_t mesh;
  switch(ArchiveExtent) {
    case SEQUENTIAL_OUTPUT:
      return(local);
      break;
    case GLOBALIZED_OUTPUT:
      return(*(local.origin));
      break;
    case PARTITIONED_OUTPUT:
      return(local);
      break;
    default:
      check_error(-1, "illicit output mode", __LINE__, __FILE__, 1);
    }

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> int archiving_UGdummy2D_template(const char *localname,const mesh_t & local, const char *name, const char *units,const T *buffer, T mask, int frame, int discretisation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int i,k,l,m,n;
  double *ug2D_LGP2,*ug2D_LGP1,*ug2D_LGP0,*ug2D_NCP1,*nvalues;
  FILE *in;
  int count,ncid,variable_id;
  int file_exist;
  cdfvar_t variable;
  T scale,offset;
  char filename[1024];
  mesh_t mesh;
  discretisation_t descriptor,gdescriptor;
  extern discretisation_t DescriptorLocal2Global_generic_template(mesh_t mesh, discretisation_t decriptor, int cpu, int target);

  if(rootname==0) rootname=strdup(output_path);
  
  char *s=strstr( (char *) localname, (char *) rootname);
  if(s==0) {
    sprintf(filename, "%s/%s", rootname, localname);
    }
  else {
//     printf("\n warning, output path found in filename, this is not a good practice\n");
//     printf("output path :%s \n",rootname);
//     printf("filename    :%s \n\n",localname);
    sprintf(filename, "%s", localname);
    }
    
//   switch(gArchiveExtent) {
//     case SEQUENTIAL_OUTPUT:
//     case GLOBALIZED_OUTPUT:
//     case PARTITIONED_OUTPUT:
//       break;
//     default:
//       check_error(-1, "illicit output mode", __LINE__, __FILE__, 1);
//     }
    
//   switch(gArchiveExtent) {
//     case SEQUENTIAL_OUTPUT:
//       sprintf(filename, "%s/%s", rootname, localname);
//       break;
//     case GLOBALIZED_OUTPUT:
//       sprintf(filename, "%s/%s", rootname, localname);
//       break;
//     case PARTITIONED_OUTPUT:
//       sprintf(filename, "%s/%s", rootname, localname);
//       break;
//     default:
//       check_error(-1, "illicit output mode", __LINE__, __FILE__, 1);
//     }

  mesh=local;

/*------------------------------------------------------------------------
  check wether file exist or not*/
  in = fopen(filename, "r");
  if(in == NULL) {
    file_exist=0;
    }
  else {
    file_exist=1;
    fclose(in);
    status=fe_CheckDimensions(filename, mesh);
    if (status!=0) file_exist=0;
    }

  descriptor=get_descriptor(mesh,discretisation);
//  gdescriptor=DescriptorLocal2Global_generic_template(mesh, descriptor, gCPU_ID, 0);
  gdescriptor=descriptor;

  if(file_exist == 0) {
    if(mesh.nlayers==0) mesh.nlayers=1;
    if(mesh.nlevels==0) mesh.nlevels=2;
    if(gCPU_ID==gCPU_MASTER) {
/**----------------------------------------------------------------------------
      get global mesh if needed (globalized output in parallel mode) */
      mesh=output_mesh(local,SEQUENTIAL_OUTPUT);
      if(mesh.nlayers==0) mesh.nlayers=1;
      if(mesh.nlevels==0) mesh.nlevels=2;
      status = fe_savemesh3d(filename, mesh, 1);
      }
    }

  if(gCPU_ID==gCPU_MASTER) {
      mesh=output_mesh(local,SEQUENTIAL_OUTPUT);
  switch (discretisation) {
    case LGP0:
      status =fe_savediscretisation(filename, mesh, gdescriptor);
      break;
    case LGP1:
      status =fe_savediscretisation(filename, mesh, gdescriptor);
      break;
    case DGP1:
      status =fe_savediscretisation(filename, mesh, gdescriptor);
      break;
    case NCP1:
      break;
    case DNP1:
      status =fe_savediscretisation(filename, mesh, gdescriptor);
      break;
    case LGP2:
      status =fe_savediscretisation(filename, mesh, gdescriptor);
      break;
    case CQP0:
      status =fe_savediscretisation(filename, mesh, gdescriptor);
      break;
    case CQN1:
      status =fe_savediscretisation(filename, mesh, gdescriptor);
      break;
    default:
      check_error(-1, "discretisation not implemented in archive procedure", __LINE__, __FILE__, 1);
      break;
      }
  }

/**------------------------------------------------------------------------

  2D section - 2D section - 2D section - 2D section - 2D section - 2D section

-------------------------------------------------------------------------*/
 
  count=frame;
  scale=1.0;
  offset=0.;

  if(gCPU_ID==gCPU_MASTER) {
//   check_error(ENOEXEC, "*** UNITS FORCED TO m, EVEN IF IT IS WRONG TO DO SO ***", __LINE__, __FILE__);
//   #warning ** UNITS FORCED TO m, EVEN IF IT IS WRONG TO DO SO ***
  switch(frame){
    case NOFRAME:
      switch (discretisation) {
        case LGP0:
          variable =poc_variable_UG2D( name, mask, units, scale, offset, name,"M");
          break;

        case LGP1:
          variable =poc_variable_UG2D( name, mask, units, scale, offset, name,"N");
          break;

        case DGP1:
          variable =poc_variable_UG2D( name, mask, units, scale, offset, name,DGP1_standards[NNODES_DIMNAME]);
          break;

        case NCP1:
          variable =poc_variable_UG2D( name, mask, units, scale, offset, name,"E");
          break;

        case DNP1:
          variable =poc_variable_UG2D( name, mask, units, scale, offset, name,DNP1_standards[NNODES_DIMNAME]);
          break;

        case LGP2:
          variable =poc_variable_UG2D( name, mask, units, scale, offset, name,LGP2_standards[NNODES_DIMNAME]);
          break;

        case CQP0:
          variable =poc_variable_UG2D( name, mask, units, scale, offset, name,CQP0_standards[NNODES_DIMNAME]);
          break;

        case CQN1:
          variable =poc_variable_UG2D( name, mask, units, scale, offset, name,CQN1_standards[NNODES_DIMNAME]);
          break;

        default:
          check_error(-1, "discretisation not implemented in archive procedure", __LINE__, __FILE__, 1);
          break;
        }
      break;
    default:
      switch (discretisation) {
    case LGP0:
      variable =poc_variable_UG3D( name, mask, units, scale, offset, name,"T","M");
      break;

    case LGP1:
      variable =poc_variable_UG3D( name, mask, units, scale, offset, name,"T","N");
      break;

    case DGP1:
      variable =poc_variable_UG3D( name, mask, units, scale, offset, name,"T",DGP1_standards[NNODES_DIMNAME]);
      break;

    case NCP1:
      variable =poc_variable_UG3D( name, mask, units, scale, offset, name,"T","E");
      break;

    case DNP1:
      variable =poc_variable_UG3D( name, mask, units, scale, offset, name,"T",DNP1_standards[NNODES_DIMNAME]);
      break;

    case LGP2:
      variable =poc_variable_UG3D( name, mask, units, scale, offset, name,"T",LGP2_standards[NNODES_DIMNAME]);
      break;

    case CQP0:
      variable =poc_variable_UG3D( name, mask, units, scale, offset, name,"T",CQP0_standards[NNODES_DIMNAME]);
      break;

    case CQN1:
      variable =poc_variable_UG3D( name, mask, units, scale, offset, name,"T",CQN1_standards[NNODES_DIMNAME]);
      break;

    default:
      check_error(-1, "discretisation not implemented in archive procedure", __LINE__, __FILE__, 1);
      break;
    }
      break;
    }
  }
  
#ifdef HAVE_MPI
  status = MPI_Barrier(MPI_COMM_WORLD);
#endif
  
   mesh=local;

/**----------------------------------------------------------------------------
  check if variable already exists */
  if(gCPU_ID==gCPU_MASTER) {
  variable.id= cdf_identify((const char *) filename, (const char *) name);
  if(variable.id==-1) {
    status = cdf_createvariable((char *) filename, &(variable));
    }
  variable_id = variable.id;
  variable.destroy();
  }
#ifdef HAVE_MPI
  status = MPI_Barrier(MPI_COMM_WORLD);
#endif
  
  variable_id= cdf_identify(filename, name);
  switch(frame){
    case NOFRAME:
      status = poc_put_UG2D( filename, mesh, variable_id, buffer);
      break;
    default:
      status = poc_put_UG3D((char *) filename, mesh, count, variable_id, buffer);
      break;
    }
terminate:
  return (0);

error:
  printf("archiving_UGarchive failed ..\n");
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int archiving_UGdummy2D(const char *localname,const mesh_t & mesh, const char *name, const char *units,const float *buffer, float mask, int discretisation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
//  float mask=1.e+10;

  status=archiving_UGdummy2D_template(localname, mesh, name, units, buffer, mask, NOFRAME, discretisation);
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int archiving_UGdummy2D(const char *localname,const mesh_t & mesh, const char *name, const char *units,const double *buffer, double mask, int discretisation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;

/**----------------------------------------------------------------------------
  should create file, get var id, then invoke poc_put */
  status=archiving_UGdummy2D_template(localname, mesh, name, units, buffer, mask, NOFRAME, discretisation);
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int archiving_UGdummy2D(const char *localname,const mesh_t & mesh, const char *name, const char *units,const float *buffer, float mask, int frame, int discretisation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
//  float mask=1.e+10;

  status=archiving_UGdummy2D_template(localname, mesh, name, units, buffer, mask, frame, discretisation);
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int archiving_UGdummy2D(const char *localname,const mesh_t & mesh, const char *name, const char *units,const double *buffer, double mask, int frame, int discretisation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
//  double mask=1.e+10;

  status=archiving_UGdummy2D_template(localname, mesh, name, units, buffer, mask, frame, discretisation);
  return (status);
}


// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
// template <typename T> int archiving_UGdummy2D_template(const char *localname, mesh_t & mesh, const char *name1, const char *name2, const char *units, complex <T> *buffer, complex <T> cmask, int frame, int discretisation)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   int status;
//   int n,nnodes;
//   double *tmp, mask;
// 
//   switch (discretisation) {
//     case LGP0:
//       nnodes=mesh.LGP0descriptor.nnodes;
//       break;
// 
//     case LGP1:
//       nnodes=mesh.LGP1descriptor.nnodes;
//       break;
// 
//     case DGP1:
//       nnodes=mesh.DGP1descriptor.nnodes;
//       break;
// 
//     case NCP1:
//       nnodes=mesh.nedges;
//       break;
// 
//     case DNP1:
//       nnodes=mesh.DNP1descriptor.nnodes;
//       break;
// 
//     case LGP2:
//       nnodes=mesh.LGP2descriptor.nnodes;
//       break;
// 
//     default:
//       check_error(-1, "discretisation not implemented in archive procedure", __LINE__, __FILE__, 1);
//       break;
//     }
//     
//   tmp=new double[nnodes];
//   mask=abs(cmask);
//   
//   for (n=0;n<nnodes;n++) {
//     tmp[n]=abs(buffer[n]);
//     }
//   status=archiving_UGdummy2D(localname, mesh, name1, units, tmp, mask, frame, discretisation);
//   
//   for (n=0;n<nnodes;n++) {
//     tmp[n]=-arg(buffer[n])*180./M_PI;
//     }
//   status=archiving_UGdummy2D(localname, mesh, name2, "degrees", tmp, mask, frame, discretisation);
//   
//   delete[] tmp;
//   
//   return (status);
// }

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int archiving_UGdummy2D_template(const char *localname, const mesh_t & mesh, const char *name1, const char *name2, const char *units, complex <T> *buffer, complex <T> cmask, int frame, int discretisation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int n,nnodes;
  T *tmp, mask;

  switch (discretisation) {
    case LGP0:
      nnodes=mesh.LGP0descriptor.nnodes;
      break;

    case LGP1:
      nnodes=mesh.LGP1descriptor.nnodes;
      break;

    case DGP1:
      nnodes=mesh.DGP1descriptor.nnodes;
      break;

    case NCP1:
      nnodes=mesh.nedges;
      break;

    case DNP1:
      nnodes=mesh.DNP1descriptor.nnodes;
      break;

    case LGP2:
      nnodes=mesh.LGP2descriptor.nnodes;
      break;

    case CQP0:
      nnodes=mesh.CQP0descriptor.nnodes;
      break;

    case CQN1:
      nnodes=mesh.CQN1descriptor.nnodes;
      break;

    default:
      check_error(-1, "discretisation not implemented in archive procedure", __LINE__, __FILE__, 1);
      break;
    }
    
  tmp=new T[nnodes];
  mask=(T) abs(cmask);
  
  for (n=0;n<nnodes;n++) {
    tmp[n]=abs(buffer[n]);
    }
  status=archiving_UGdummy2D(localname, mesh, name1, units, tmp, mask, frame, discretisation);
  
  for (n=0;n<nnodes;n++) {
    tmp[n]=-arg(buffer[n])*180./M_PI;
    }
  status=archiving_UGdummy2D(localname, mesh, name2, "degrees", tmp, mask, frame, discretisation);
  
  delete[] tmp;
  
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int archiving_UGdummy2D(const char *localname, const mesh_t & mesh, const char *name1, const char *name2, const char *units, complex <double> *buffer, complex <double> cmask, int frame, int discretisation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=archiving_UGdummy2D_template(localname, mesh, name1, name2, units, buffer, cmask, frame, discretisation);
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int archiving_UGdummy2D(const char *localname, const mesh_t & mesh, const char *name1, const char *name2, const char *units, complex <double> *buffer, complex <double> cmask, int discretisation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=archiving_UGdummy2D_template(localname, mesh, name1, name2, units, buffer, cmask, NOFRAME, discretisation);
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int archiving_UGdummy2D(const char *localname, const mesh_t & mesh, const char *name1, const char *name2, const char *units, complex <float> *buffer, complex <float> cmask, int discretisation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=archiving_UGdummy2D_template(localname, mesh, name1, name2, units, buffer, cmask, NOFRAME, discretisation);
  return (status);
}

// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
//
//   int archiving_UGdummy2D(const char *localname, mesh_t & mesh, const char *name1, const char *name2, const char *units, complex <double> *buffer, int discretisation)
//
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   int status;
//   complex <double> cmask=1.e+10;
//
//   status=archiving_UGdummy2D(localname, mesh, name1, name2, units, buffer, cmask, (int) 0, discretisation);
//   return (status);
// }




/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int read_UG2Datlas(char *filename, mesh_t & mesh,fcomplex * &tide,fcomplex & mask,char *discretisation,int iteration,int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, j, k, l, n,node,frame, status;
  float spec[2];
  fcomplex zz,*z=NULL,cmask;
  char units[256],amplitude[256],phase[256];
  char *sunit=NULL;
  int nnodes;
  discretisation_t descriptor;

  cdfgbl_t global;
  int id;
  variable_t varinfo;

  printf("#################################################################\n");
  printf("load harmonic file : %s\n",filename);
  status= cdf_globalinfo(filename,&global,0);
  if(status !=0) {
    printf("cannot open %s\n",filename);
    goto error;
    }

  id=discretisation_from_name(discretisation);
  status=fe_readgeometry(filename, &mesh);
  status=fe_readdiscretisation(filename, &mesh, 0, id);
  
  descriptor=get_descriptor(mesh,id);
  
  exitIfNull(
    tide=new complex<float>[descriptor.nnodes]
    );
/*-----------------------------------------------------------------------
  load netcdf variable */
  frame=iteration;
  sprintf(amplitude,"a_eta_%s",discretisation);
  sprintf(phase,"G_eta_%s",discretisation);
  
  status=poc_get_UG3D(filename, frame, amplitude, phase, tide);
  if(status !=0) {
    printf("cannot load frame %d in %s\n",frame,filename);
    goto error;
    }
  cmask=fcomplex(9999.,9999.);

  return (status);

error:
  status=-1;
  return (status);
}



/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

//  int archiving_UGdummy3D(const char *localname, mesh_t & mesh,const char *name, const char *units, double **buffer, double mask, int frame, int h_discretisation, int v_discretisation, int l_discretisation)
  int archiving_UGdummyTZH(const char *localname, mesh_t & mesh,const char *name, const char *units, double **buffer, double mask, int frame, int h_discretisation, int v_discretisation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int i,k,l,m,n;
  double *nvalues;
  FILE *in;
//  double mask=1.e+10;
  int count,ncid,variable_id;
  int file_exist;
  cdfvar_t variable;
  double scale,offset;
  char filename[1024];

  if(rootname==0) rootname=strdup(output_path);
  char *s=strstr( (char *) localname, (char *) rootname);
  if(s==0) {
    sprintf(filename, "%s/%s", rootname, localname);
    }
  else {
//     printf("\n warning, output path found in filename, this is not a good practice\n");
//     printf("output path :%s \n",rootname);
//     printf("filename    :%s \n\n",localname);
    sprintf(filename, "%s", localname);
    }
/*------------------------------------------------------------------------
  see if file exist*/
  in = fopen(filename, "r");
  if(in == NULL) {
    file_exist=0;
    }
  else {
    file_exist=1;
    }

  if(file_exist == 0) {
    if(mesh.nlayers==0) mesh.nlayers=1;
    if(mesh.nlevels==0) mesh.nlevels=2;
    status = fe_savemesh3d(filename, mesh, 1);
    }
/**----------------------------------------------------------------------------
  temporary, for plotting in xscan (to be fixed in xscan) */
  status =fe_savediscretisation(filename, mesh, mesh.LGP1descriptor);

  switch (h_discretisation) {
    case LGP0:
      status =fe_savediscretisation(filename, mesh, mesh.LGP0descriptor);
      break;
    case LGP1:
//      status =fe_savediscretisation(filename, mesh, mesh.LGP1descriptor);
      break;
    case DGP1:
      status =fe_savediscretisation(filename, mesh, mesh.DGP1descriptor);
      break;
    case NCP1:
      break;
    case DNP1:
      status =fe_savediscretisation(filename, mesh, mesh.DNP1descriptor);
      break;
    case LGP2:
      status =fe_savediscretisation(filename, mesh, mesh.LGP2descriptor);
      break;
      }

/**------------------------------------------------------------------------

  2D section - 2D section - 2D section - 2D section - 2D section - 2D section

-------------------------------------------------------------------------*/

  count=frame;

  scale=1.0;
  offset=0.;


  if(frame==-1) {
    variable =poc_variable_UG4D( name, mask, units, scale, offset, name, "T",h_discretisation, v_discretisation);
    }
  else if(frame==NOFRAME) {
    switch (v_discretisation) {
      case LAYERS:
        variable =poc_variable_UG3D( name, mask, units, scale, offset, name, "K",h_discretisation);
        break;
      case LEVELS:
        variable =poc_variable_UG3D( name, mask, units, scale, offset, name, "L",h_discretisation);
        break;
        }
    }
  else {
    variable =poc_variable_UG4D( name, mask, units, scale, offset, name, "T", h_discretisation, v_discretisation);
    }
    
//   switch (l_discretisation) {
//     case LGP0:
//       variable.att[4].data=(char *) poc_strdup("time z-LGP0 lat lon");
//       break;
//     case LGP1:
//       variable.att[4].data=(char *) poc_strdup("time z-LGP1 lat lon");
//       break;
//     case DGP1:
//       variable.att[4].data=(char *) poc_strdup("time z-DGP1 lat lon");
//       break;
//     case NCP1:
//       variable.att[4].data=(char *) poc_strdup("time z-NCP1 lat lon");
//       break;
//     case DNP1:
//       variable.att[4].data=(char *) poc_strdup("time z-DNP1 lat lon");
//       break;
//     case LGP2:
//       variable.att[4].data=(char *) poc_strdup("time z-LGP2 lat lon");
//       break;
//       }

  switch (v_discretisation) {
    case LAYERS:
      variable.att[4].data=(char *) poc_strdup("time z-LGP0 lat lon");
      break;
    case LEVELS:
      variable.att[4].data=(char *) poc_strdup("time z-LGP1 lat lon");
      break;
      }
      
      
/**----------------------------------------------------------------------------
  check if variable already exists */
  variable.id= cdf_identify(filename, name);
  if(variable.id==-1) {
    status = cdf_createvariable((char *) filename, &(variable));
    }
  variable_id = variable.id;

  variable.destroy();
  
  if(frame==-1) {
    status = poc_put_UG4D((char *) filename, mesh, count, variable_id, buffer);
    }
  else if(frame==NOFRAME) {
    status = poc_put_UG3D((const char *) filename, mesh, variable_id, buffer);
    }
  else {
    status = poc_put_UG4D((char *) filename, mesh, count, variable_id, buffer);
    }

terminate:
  return (0);

error:
  printf("archiving_UGarchive failed ..\n");
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int archiving_UGdummy3D(const char *localname, mesh_t & mesh, const char *name, const char *units, double **buffer, double mask, int frame, int h_discretisation, int v_discretisation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int k,n,nnodes,nlevels;
  double **tmp;
//  int frame=NOFRAME;
  
  discretisation_t descriptor=get_descriptor(mesh,h_discretisation);
  nnodes=descriptor.nnodes;
  
  switch (v_discretisation) {
    case LAYERS:
      nlevels=mesh.nlayers;
      break;
    case LEVELS:
      nlevels=mesh.nlevels;
      break;
    default:
      check_error(-1, "discretisation not implemented in archive procedure", __LINE__, __FILE__, 1);
      break;
    }

/**----------------------------------------------------------------------------
  re-organize buffer order to save layer slices */
  tmp=new double*[nlevels];

  for(k=0;k<nlevels;k++) {
    tmp[k]=new double[nnodes];
    for (n=0;n<nnodes;n++) {
      tmp[k][n]=buffer[n][k];
      }
    }
    
  status=archiving_UGdummyTZH(localname, mesh, name, units, tmp, mask, frame, h_discretisation, v_discretisation);
  
  for(k=0;k<nlevels;k++) {
    if (tmp[k] != NULL) {
      delete[] tmp[k];
      tmp[k]=NULL;
      }
    }
  if (tmp != NULL){
    delete[] tmp;
    tmp=NULL;
    }
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int archiving_UGdummy3D(const char *localname, mesh_t & mesh, const char *name, const char *units, double **buffer, double mask, int h_discretisation, int v_discretisation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int frame=NOFRAME;
  status=archiving_UGdummy3D(localname, mesh, name, units, buffer, mask, frame, h_discretisation, v_discretisation);
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int archiving_UGdummy3D(const char *localname, mesh_t & mesh, const char *name1, const char *name2, const char *units, complex <double> **buffer, complex <double> cmask,  int frame, int h_discretisation, int v_discretisation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int k,n,nnodes,nlevels;
  double **tmp, mask=1.e+10;

  switch (h_discretisation) {
    case LGP2:
      nnodes=mesh.LGP2descriptor.nnodes;
      break;

    case LGP1:
      nnodes=mesh.LGP1descriptor.nnodes;
      break;

    case DGP1:
      nnodes=mesh.DGP1descriptor.nnodes;
      break;

    case LGP0:
      nnodes=mesh.LGP0descriptor.nnodes;
      break;

    case NCP1:
      nnodes=mesh.nedges;
      break;

    default:
      check_error(-1, "discretisation not implemented in archive procedure", __LINE__, __FILE__, 1);
      break;
    }
  switch (v_discretisation) {
    case LAYERS:
      nlevels=mesh.nlayers;
      break;
    case LEVELS:
      nlevels=mesh.nlevels;
      break;
    default:
      check_error(-1, "discretisation not implemented in archive procedure", __LINE__, __FILE__, 1);
      break;
    }

/**----------------------------------------------------------------------------
  re-arrange by layers */
  tmp=new double*[nlevels];
  for(k=0;k<nlevels;k++) {
    tmp[k]=new double[nnodes];
    for (n=0;n<nnodes;n++) {
      tmp[k][n]=abs(buffer[n][k]);
      }
    }
    
  status=archiving_UGdummyTZH(localname, mesh, name1, units, tmp, mask, frame, h_discretisation, v_discretisation);

  for(k=0;k<nlevels;k++) {
//    tmp[n]=new double[nnodes];
    for (n=0;n<nnodes;n++) {
      tmp[k][n]=-arg(buffer[n][k])*180./M_PI;
      }
    }
    
  status=archiving_UGdummyTZH(localname, mesh, name2, "degrees", tmp, mask, frame, h_discretisation, v_discretisation);

  for(k=0;k<nlevels;k++) {
    delete[] tmp[k];
    }
  delete[] tmp;
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int archiving_UGdummy3D(const char *localname, mesh_t & mesh, const char *name1, const char *name2, const char *units, complex <double> **buffer, complex <double> mask, int h_discretisation, int v_discretisation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int frame=NOFRAME;
//  complex <double> cmask=1.e+10;
  
  status= archiving_UGdummy3D(localname, mesh, name1, name2, units, buffer, mask, frame, h_discretisation, v_discretisation);
  
  return (status);
  
}


