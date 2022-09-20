
/**************************************************************************

  T-UGO tools, 2006-2011

  Unstructured Ocean Grid initiative

***************************************************************************/
/** \file

\author  Florent Lyard      LEGOS/CNRS, Toulouse, France. florent.lyard@legos.obs-mip.fr
\author  Yoann Le Bars      LEGOS, Toulouse, France (PhD)
\author  Laurent Roblou     LEGOS/CNRS, Toulouse, France
\author  Damien Allain      LEGOS/CNRS, Toulouse, France
\author  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
\author  Frédéric Dupont    Université de Laval à Québec, Canada

VERSION :
\date Aug 11, 2011 : Yves Soufflet : First draft

<!-- USE FIND AND REPLACE TO EDIT THIS LINE, SO THAT print_help IS ALSO UPDATED -->
\brief Extract a vertical profile from a 3D Unstructured NETCDF and print it into an ascii file

\todo EVERYTHING

*/
/*----------------------------------------------------------------------------*/

#define MAIN_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "version-macros.def" //for VERSION and REVISION

#include "tools-define.h"
#include "tools-structures.h"

#include "poc-netcdf.def"

#include "fe.h"
#include "map.h"
#include "geo.h"
#include "sym-io.h"
#include "polygones.h"
#include "grd.h"
#include "netcdf-proto.h"
#include "poc-time.h"

#include <netcdf.h>
#include "fe.def"
#include "rutin.h"
#include "functions.h"

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void print_help(const char *prog_name)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/** \brief prints help of this programme.

\param prog_name name of the program to be printed by this help function : as[0]
*/
/*----------------------------------------------------------------------------*/
{
/** The code of the body of this function is :
    \code /**/ // COMPILED CODE BELOW !!!
  printf("\n"
    " NAME AND VERSION\n"
    "  %s version " VERSION " " REVISION "\n", prog_name);
  printf("\n USE");
  printf("\n   %s -f FILE -p\"lat lon\" \n",prog_name);
  printf("\n DESCRIPTION");
/* USE FIND AND REPLACE TO EDIT THE LINE BELOW, SO THAT THE FILE HEAD IS ALSO UPDATED */
  printf("\n   This program extract a vertical profile from a unstructured NETCDF file");
    printf(
   "\n OPTIONS :\n"
   "   -f       specify the file name of the netcdf\n"
   "   -p       specify the geographical point where the profile is extracted\n"
   "   -v       list of variables name to be extracted\n"
   "   -o       output file name\n"
   "   -h       print the help\n");

  __ERR_BASE_LINE__("exiting\n");exit(-1); /** \endcode */
}



/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  geo_t projection;
  int n,status,ele;
  char *keyword;
  char *output=NULL,*input=NULL,*lon, *lat,*mooring_name,**varname,**dimname;
  nc_type *vartype;
 // typedef int size_t;
  size_t *dimlgth;
  //int *dimlgth;
  mesh_t mesh;
  int nlayers,varNum,verbose=0;
  int ncid;
  cdfgbl_t g_info;
  cdfvar_t vU_info;
 double factor, *variable,*var2,*var3,*var4,*var5;
 decoded_t decoded;
 int *ids,id,dimv,ldim,tdim,pdim;
 int ndimsp,nvarsp,ngattsp,unlimdimidp,dim,att,*nattsp,*dimids,ndim;
 FILE *out;
 double **var_array,**var_array2,**var_array3,**var_array4,**var_array5;
 
  fct_echo(argc,argv);
 
 n=1;
 while (n < argc) {
    keyword=strdup(argv[n]);
    switch (keyword[0]) {
      case '-':
        switch (keyword[1])
        {
        case 'f' :
          input= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'p' :
          lon= strdup(argv[n+1]);
          lat=strdup(argv[n+2]);
          n++;
          n++;
          n++;
          break;

        case 'o' :
          output= strdup(argv[n+1]);
          n++;
          n++;
          break;
        case 'm' :
          mooring_name= strdup(argv[n+1]);
          n++;
          n++;
          break;
            
        case 'h' :
           print_help(argv[0]);
           break;

        default:
          printf("unknown option %s\n",keyword);
          print_help(argv[0]);
          exit(-1);
        }
        break;

      default:
         printf("unknown option %s\n",keyword);
         print_help(argv[0]);
         exit(-1);
         break;
      }
    free(keyword);
    }



/*-----------------------------------------------------------------------------
 check input line */

  if((lon!=NULL)&&(lat!=NULL)){
    
    printf("Extraction of profile at position: %s %s \n",lon,lat);
    
    }

  if(input!=NULL) {
    printf("reading %s \n",input);
    status=fe_readmesh3d(input,&mesh, 0);
    
    }
  status=cdf_globalinfo((const char*) input,&g_info,verbose);
 
  status=nc_open(input,0, &ncid);

  nvarsp=g_info.nvarsp;
  varname = new char *[nvarsp];
  nattsp  = new int[nvarsp];
  vartype = new nc_type[nvarsp];
  //dimids  = new int [nvarsp];
  //ndim    = new int[nvarsp];
  dimname = new char *[ndimsp];
  dimlgth = new size_t[ndimsp];


  /*------------------------------------------------------------------------
   *   inquire dimensions name and length*/
 /* for(dim = 0; dim < ndimsp; dim++) {
    dimname[dim] = new char[1024];
    //status = nc_inq_dim(ncid, dim, dimname[dim], &dimlgth[dim]);
    if(status != NC_NOERR) printf("ERROR WHILE INQUIRING NC FILE");
    //if (verbose) printf("dimension %d,name1 %s, length %d \n",dim,dimname[dim],dimlgth[dim]);
   }
*/
  for(varNum=1;varNum<g_info.nvarsp;varNum++) {
        /*----------------------------------------------------------------------------
         *     get unstructured netcdf variable details*/
    status=cdf_varinfo((const char*) input, varNum, &vU_info);
    status=poc_decode_axis(vU_info,g_info,&decoded);
    status=poc_decode_mask(vU_info,&decoded);
    status=poc_decode_names(vU_info,&decoded);
    status=poc_decode_units(vU_info,&decoded,&factor);
  
    tdim=1;//dimlgth[ids[0]];
    ldim=mesh.nlevels;//dimlgth[ids[1]];
    nlayers=mesh.nlayers;
    pdim=mesh.ntriangles;//[ids[2]];
    
    printf("treating var %s \n",vU_info.name);
   if (strcmp(vU_info.name,"a_u3D" )==0){
      varname[varNum]= new char[1024];
      status = nc_inq_varndims(ncid, varNum, &ndim);
      if(status != NC_NOERR) printf("ERROR WHILE INQUIRING NC FILE");
      dimids = new int[vU_info.ndim];
      ids = new int[vU_info.ndim];
      status = nc_inq_var(ncid, varNum, varname[varNum], &vartype[varNum], &ndim,
      dimids, &nattsp[varNum]);
    variable= new double[tdim*ldim*pdim];
    status=nc_get_var_double(ncid,varNum,variable);
   }
   if (strcmp(vU_info.name,"G_u3D" )==0){
      varname[varNum]= new char[1024];
      status = nc_inq_varndims(ncid, varNum, &ndim);
      if(status != NC_NOERR) printf("ERROR WHILE INQUIRING NC FILE");
      dimids = new int[vU_info.ndim];
      ids = new int[vU_info.ndim];
      status = nc_inq_var(ncid, varNum, varname[varNum], &vartype[varNum], &ndim,
      dimids, &nattsp[varNum]);
    var2= new double[tdim*ldim*pdim];
    status=nc_get_var_double(ncid,varNum,var2);
   }
   if (strcmp(vU_info.name,"a_v3D" )==0){
      varname[varNum]= new char[1024];
      status = nc_inq_varndims(ncid, varNum, &ndim);
      if(status != NC_NOERR) printf("ERROR WHILE INQUIRING NC FILE");
      dimids = new int[vU_info.ndim];
      ids = new int[vU_info.ndim];
      status = nc_inq_var(ncid, varNum, varname[varNum], &vartype[varNum], &ndim,
      dimids, &nattsp[varNum]);
    var3= new double[tdim*ldim*pdim];
    status=nc_get_var_double(ncid,varNum,var3);
   }
   if (strcmp(vU_info.name,"G_v3D" )==0){
      varname[varNum]= new char[1024];
      status = nc_inq_varndims(ncid, varNum, &ndim);
      if(status != NC_NOERR) printf("ERROR WHILE INQUIRING NC FILE");
      dimids = new int[vU_info.ndim];
      ids = new int[vU_info.ndim];
      status = nc_inq_var(ncid, varNum, varname[varNum], &vartype[varNum], &ndim,
      dimids, &nattsp[varNum]);
    var4= new double[tdim*ldim*pdim];
    status=nc_get_var_double(ncid,varNum,var4);
   }
   if (strcmp(vU_info.name,"Kv3D" )==0){
      varname[varNum]= new char[1024];
      status = nc_inq_varndims(ncid, varNum, &ndim);
      if(status != NC_NOERR) printf("ERROR WHILE INQUIRING NC FILE");
      dimids = new int[vU_info.ndim];
      ids = new int[vU_info.ndim];
      status = nc_inq_var(ncid, varNum, varname[varNum], &vartype[varNum], &ndim,
      dimids, &nattsp[varNum]);
    var5= new double[tdim*ldim*pdim];
    status=nc_get_var_double(ncid,varNum,var5);
   }
   }


  var_array=new double* [nlayers];
  var_array2=new double* [nlayers];
  var_array3=new double* [nlayers];
  var_array4=new double* [nlayers];
  var_array5=new double* [nlayers];
  for (int i=0;i<tdim;i++){
  for (int k=0;k<nlayers;k++){
          var_array[k]=new double[pdim];
          var_array2[k]=new double[pdim];
          var_array3[k]=new double[pdim];
          var_array4[k]=new double[pdim];
          var_array5[k]=new double[pdim];
          for (int j=0;j<pdim;j++){
          var_array[k][j]=variable[i*pdim*nlayers+k*pdim+j];
          var_array2[k][j]=var2[i*pdim*nlayers+k*pdim+j];
          var_array3[k][j]=var3[i*pdim*nlayers+k*pdim+j];
          var_array4[k][j]=var4[i*pdim*nlayers+k*pdim+j];
          var_array5[k][j]=var5[k*pdim+j];
          }
  }
  }
  ele=fe_whichelement(mesh,strtod(lon,NULL),strtod(lat,NULL));
  printf("vertical profile will be made over element number %d \n",ele);
  
  //retrieve level's depth
  int n1,n2,n3;
  n1=mesh.triangles[ele].vertex[0];
  n2=mesh.triangles[ele].vertex[1];
  n3=mesh.triangles[ele].vertex[2];
  
  double levels_depth[nlayers];
  
  
  for (int k=0;k<ldim-1;k++){
  levels_depth[k]=((mesh.vertices[n1].zlevels[k]+mesh.vertices[n2].zlevels[k]+mesh.vertices[n3].zlevels[k])
          +(mesh.vertices[n1].zlevels[k+1]+mesh.vertices[n2].zlevels[k+1]+mesh.vertices[n3].zlevels[k+1]))/6.0;
  }
  
  out=fopen(output,"w");
  fprintf(out,"# %s %s %s\n",mooring_name,lon,lat);
  fprintf(out,"# levels depth a_u3D G_u3D a_v3D G_v3D Kv3D\n");

  for (int k=0;k<nlayers;k++){
    printf("level %d %lf %lf %lf %lf %lf %lf \n",k,levels_depth[k],var_array[k][ele],var_array2[k][ele],var_array3[k][ele],var_array4[k][ele],var_array5[k][ele]);
    fprintf(out,"%d %lf %lf %lf %lf %lf %lf\n",k,levels_depth[k],var_array[k][ele],var_array2[k][ele],var_array3[k][ele],var_array4[k][ele],var_array5[k][ele]);
  }
  /*
  char* Varnames[]={"a_u3D","G_u3D","a_v3D","G_v3D","Kv3D"};
  
  for (int var=0;var<5;var++){
            printf("treating variable: %s \n",Varnames[var]);
  }
*/
 status = nc_close(ncid);

}
