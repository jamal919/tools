/**************************************************************************

  T-UGO tools, 2019

  Unstructured Ocean Grid initiative

Contributors:

  Simon Barbot       LEGOS/CNRS-CNES-CLS-UPS, Toulouse, France
  Florent Lyard      LEGOS/CNRS, Toulouse, France
  Damien Allain      LEGOS/CNRS, Toulouse, France

E-mail: florent.lyard@legos.obs-mip.fr

***************************************************************************/

#define MAIN_SOURCE 1

#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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
#include "poc-netcdf-data.hpp"

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*
interpolation of a slice grid (like COMODO) from listing profile observation
Options :
  -i    :input file name (.txt) [longitude, latitude, depth, variable]
  -m    :bathymetry file of the wanted domain (.nc) - regular grid
  -o    :output file name (.nc)
*/
{
  float *lon,*lat,*depth,*buffer,*bufferLon,*bufferLat,*bufferD,*bufferV,*var;
  float *lonprofiles,*priori;
  float mask=1.e+10,shift=-13.1,dlon,sum;
// CAUTION ! shift value must be set for each different file
// TODO : find the shift used in mesh-academics for the bathymetry
//   float  *buffer[2];

  int status=0;
  int fmt,id,count,delta,idrows,idnonzero,card,left,right;
  int i,j,k,l,n;
  int nobs,nlevels,nprofiles,nrows,nparams,nmaskobs,nmaskparams;
  int *tmp,*idprofile;
  int boundary;
  
  bool test_lon;

  FILE *in,*out;
  char text[256];
  char varname[256];
  char *keyword,*input,*output,*bathyfile;
//   char boundary;
  string standards[2]={"depth","bathymetry"};
  mesh_t mesh;
  poc_global_t Dglobal;
  string name="";
  poc_data_t <float> vD;
  grid_t grid;
  size_t *sortlist;
//   size_t *sortlist;
//   cdfvar_t info;
//   cdfgbl_t data_info,grid_info;
//   date_t origine,reference;

  float *X,*Y;
//   float **A1,**A2;
  hypermatrix_t A,A1,A2;
  
#define DEBUG 1

#define NO_BOUNDARY 0
#define LEFT_BOUNDARY 1
#define RIGHT_BOUNDARY 2

  fct_echo(argc,argv);
 
  n=1;
  while (n < argc) {
    keyword=strdup(argv[n]);
    switch (keyword[0]) {
      case '-':
        switch (keyword[1]) {
        case 'i' :
          input= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'b' :
          bathyfile= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'o' :
          output= strdup(argv[n+1]);
          n++;
          n++;
          break;

        default:
          __OUT_BASE_LINE__("unknown option %s\n",keyword);
          exit(-1);
        }
        break;

      default:
        __OUT_BASE_LINE__("not enough options\n");
        exit(-1);
      }
      free(keyword);
    }

  printf("Loading the observations -- ");

  in=fopen(input,"r");
  if(in==NULL) {__ERR_BASE_LINE__("exiting\n");exit(-1);}
  fscanf(in, "%d", &nobs);
  fgets(varname,256,in);
  depth = new float[nobs];
  buffer = new float[nobs];
  bufferLon = new float[nobs];
  bufferLat = new float[nobs];
  bufferD = new float[nobs];
  bufferV = new float[nobs];
  nprofiles = 0;
  tmp = new int[nobs];
  for(n=0; n<nobs;n++) {
    fscanf(in,"%f\t%f\t%f\t%f", &bufferLon[n], &bufferLat[n], &bufferD[n], &bufferV[n]);
    fgets(text,256,in);
    bufferLon[n]-=shift;
    if(bufferD[n]==0) {
      tmp[nprofiles]=n;
      nprofiles++;
      }
    else if(bufferD[n]>0){
      bufferD[n]=-bufferD[n];
      }
    }
  fclose(in);
  idprofile = new int[nprofiles];
  for(n=0; n<nprofiles;n++) {
    idprofile[n]=tmp[n];
    }
  
  delete [] tmp;
  
  printf("Done\n");
  printf("Reshape the profiles to negative ascending -- ");
  
  lonprofiles = new float[nprofiles]; 

  for(n=0; n<nprofiles;n++) {
    if (n<nprofiles-1) nlevels=idprofile[n+1]-idprofile[n];
    else nlevels=nobs-idprofile[n];
    
    lonprofiles[n]=bufferLon[idprofile[n]];
    sortlist=sort(&bufferD[idprofile[n]],nlevels);
    for(k=0;k<nlevels;k++){
      depth[idprofile[n]+k]=bufferD[idprofile[n]+sortlist[k]];
      buffer[idprofile[n]+k]=bufferV[idprofile[n]+sortlist[k]];
      }
    delete [] sortlist;
  }
  
  delete [] bufferD;
  delete [] bufferV;
  
  printf("Done\n");
  printf("Loading the bathymetry -- ");
  
  if(bathyfile != NULL) {
    status=poc_inq(bathyfile,&Dglobal);
    }
  else {
    printf("no mesh file specified; abort...\n");
    goto error;
    }
// try standard names
  count=0;
  id=-1;
  while(id==-1) {
    name=standards[count];
    id=Dglobal.variables.find(name.c_str());
    count++;
    if(count==2) break;
    }
  if(id==-1)  check_error(-1, "depth initialisation failed", __LINE__, __FILE__, 1);

  status=vD.init(Dglobal,name.c_str());
  if(status != 0) goto error;

// load grid
  status=poc_get_grid(bathyfile,vD.info,&grid,0,0);
  if(status != 0) goto error;
  
// load bathymetry variable
  status=vD.read_scaled_data(bathyfile);
  if(status != 0) goto error;
  
  printf("Done\n");
  printf("Build vertical grid -- ");
  
#if DEBUG
  grid.nz = 5;
  grid.z = new double[grid.nz];
  grid.z[0]=-3500;
  grid.z[1]=-1000;
  grid.z[2]=-100;
  grid.z[3]=-50;
  grid.z[4]=0;
#else  
  grid.nz = 104;
  grid.z = new double[grid.nz];
  grid.z[0]=round(vD.data[0]/100)*100;
  for(k=1; k<grid.nz;k++) {
    if(grid.z[k-1] < -1000) delta = 100;
    else if(grid.z[k-1] >= -1000 && grid.z[k-1] < -500) delta = 50;
    else if(grid.z[k-1] >= -500 && grid.z[k-1] < -100) delta = 10;
    else if(grid.z[k-1] >= -100 && grid.z[k-1] < -10) delta = 5;
    else if(grid.z[k-1] >= -10 && grid.z[k-1] < 0) delta = 1;
    
    grid.z[k]=grid.z[k-1]+delta;
    }
#endif

// Number of grid point in the mask
  nmaskparams=0;
  for(j=0;j<grid.nx;j++) {
    for(k=0;k<grid.nz;k++) {
      if(grid.z[k]<vD.data[j]) nmaskparams++;
      }
    }
#if DEBUG
  float *lonmaskparams;
  int *idmaskparams;
  lonmaskparams = new float[nmaskparams];
  idmaskparams = new int[nmaskparams];
  count=0;
  for(j=0;j<grid.nx;j++) {
    for(k=0;k<grid.nz;k++) {
      if(grid.z[k]<vD.data[j]) {
        lonmaskparams[count]=grid.x[j];
        idmaskparams[count]=j;
        count++;
        break;
        }
      }
    }
  printf("\nAll longitude inside the mask (%i) :\n",count);
  for(n=0;n<count;n++) printf("%f\t",grid.x[idmaskparams[n]]);
  printf("\n");
#endif

  printf("Done\n");
  printf("Interpolate observation profiles on grid levels -- ");
  
  var = new float[nprofiles*grid.nz];
  lon = new float[nprofiles*grid.nz];
  lat = new float[nprofiles*grid.nz];
  nmaskobs = 0;
  sortlist=sort(lonprofiles,nprofiles);
  for(n=0; n<nprofiles;n++) {
    if (n<nprofiles-1) nlevels=idprofile[n+1]-idprofile[n];
    else nlevels=nobs-idprofile[n];
        
    for(k=0; k<grid.nz;k++){
      status=map_interpolate1D(&buffer[idprofile[k]], &depth[idprofile[k]], mask, nlevels, grid.z[k], &var[sortlist[n]*grid.nz+k], 0);
      lon[sortlist[n]*grid.nz+k]=bufferLon[idprofile[n]];
      lat[sortlist[n]*grid.nz+k]=bufferLat[idprofile[n]];
      if(var[sortlist[n]*grid.nz+k]==mask) nmaskobs++;
      }
    }
  
  quick_sort(lonprofiles,nprofiles);
  delete [] sortlist;
  
  printf("Done\n");
  printf("Building priori solution -- ");
  
  priori = new float[grid.nz];
//  method : horizontal mean
  for(k=0;k<grid.nz;k++){
    sum=0;
    for(i=0;i<nprofiles;i++)sum=+var[i*k];
    priori[k]=sum/nprofiles;
    }
  
  printf("Done\n");
  printf("Building working matrices -- ");
// ----------------------------------------------------
// ---------------------TABLE SHAPE--------------------
// ----------------------------------------------------
//   nobs = nprofiles*grid.nz;
//   nparams=grid.nz*grid.nx;
//   X = new float[nparams];
//   Y = new float[nobs+nparams];
//   A1 = new float*[nobs];
//   for(i=0;i<nobs;i++) {
//     A1[i] = new float[nparams];
//     aset(A1[i],nparams,0.f);
//     }
//   A2 = new float*[nparams];
//   for(i=0;i<nparams;i++) {
//     A2[i] = new float[nparams];
//     aset(A2[i],nparams,0.f);
//     }
// //   Y = new float[nobs+2*nparams];
// //   A = new float*[nobs+2*nparams];
// //   for(i=0;i<nobs+2*nparams;i++) A[i] = new float[nparams];
//   
//   index=0;
//   for(i=0;i<nobs;i+=grid.nz){
//     for(j=index;j<grid.nx-1;j+=grid.nz){
//       test_lon=(abs(lon[i]-grid.x[j])<=abs(lon[i]-grid.x[j+1]));
//       if(test_lon){
//         for(n=0;n<grid.nz;n++) {A1[i+n][j+n]=var[i+n];}
//         index+=grid.nz;
//         break;
//         }
//       else{
//         index+=grid.nz;
//         }
//       }
//     } 
// 
//   for(i=0;i<nparams;i++){
//     if(i==0){
//       }
//     else if(i==nparams-1){
//       }
//     else{
//       }
//     }

// ----------------------------------------------------
// -------------------MATRIX SHAPE---------------------
// ----------------------------------------------------
// Interpolation matrix
  nrows = nprofiles*grid.nz-nmaskobs;
  nparams=grid.nz*grid.nx-nmaskparams;
  card = 2;
  X = new float[nparams];
  Y = new float[nrows+nparams];
  A1.ordering=new ordering_t;
  A1.ordering->allocate(nrows);
  
// method : nearest grid point = obs point
//   for(i=0;i<nrows;i++) A1.ordering->cardinal[i]=1;
//   A1.ordering->allocate_finalize();
//   A1.allocate();
//   index=0;
//   for(i=0;i<nrows;i+=grid.nz){
//     for(j=index;j<grid.nx-1;j+=grid.nz){
//       test_lon=(abs(lon[i]-grid.x[j])<=abs(lon[i]-grid.x[j+1]));
//       if(test_lon){
//         for(n=0;n<grid.nz;n++) {
//           if(var[i+n]==mask) break;
//           A1.ordering->pointer[i+n]=j+n;
//           A1.ordering->incidence[i+n]=j+n;
//           A1.packed[i+n]=var[i+n];
//           }
//         index+=grid.nz;
//         break;
//         }
//       else{
//         index+=grid.nz;
//         }
//       }
//     } 

// method : interpolation, weight on the two closest nodes
  for(i=0;i<nrows;i++) A1.ordering->cardinal[i]=card;
  A1.ordering->allocate_finalize();
  A1.allocate();
  idnonzero=0;
  idrows=0;
  count=0;
  j=0;
  for(i=0;i<nprofiles;i++){
    while(j<grid.nx-1){
      dlon = abs(lon[i*grid.nz]-grid.x[j]);
      if(dlon<grid.dx){
        for(k=0;k<grid.nz;k++) {
          if(var[i*grid.nz+k]==mask) {
            if(grid.z[k]>vD.data[j]) count++;
            continue;
            }
          A1.ordering->pointer[idrows]=count;
          A1.ordering->incidence[idnonzero]=count;
          A1.packed[idnonzero]=var[i*grid.nz+k]*(1-(dlon/grid.dx));
          idnonzero++;
          idrows++;
          count++;
          }
        j++;
        for(k=0;k<grid.nz;k++) {
          if(var[i*grid.nz+k]==mask){
            if(grid.z[k]>vD.data[j]) count++;
            continue;
            }
          A1.ordering->incidence[idnonzero]=count;
          A1.packed[idnonzero]=var[i*grid.nz+k]*(dlon/grid.dx);
          idnonzero++;
          count++;
          }
        j++;
//         if(i==nprofiles-1) continue;
        break;
        }
      else{
        for(k=0;k<grid.nz;k++) {
          if(grid.z[k]<vD.data[j]) {
            continue;
            }
          else count++;
          }
        }
      j++;
#if DEBUG
      if(i*grid.nz > nprofiles*grid.nz){
        printf("\nobs index value too big !\n");
        goto error;
        }
      if(idnonzero > nrows*card){
        printf("\nnon-zero index value too big !\n");
        goto error;
        }
      if(count > nparams) {
        printf("\ncolumn index value too big !\n");
        goto error;
        }
#endif
      }
    } 
  
  status = matrix_check_Zero(A1);
  
//  Horizontal gradient matrix
  nparams=grid.nz*grid.nx-nmaskparams;
  nrows=nparams;
  card = 2;
  A2.ordering=new ordering_t;
  A2.ordering->allocate(nrows);  
//   for(i=0;i<nrows;i++) A2.ordering->cardinal[i]=card;
// looking for null gradient nodes
  for(i=0;i<grid.nx;i++){
    for(k=0;k<grid.nz;k++) {
      if(i==0){
        if(grid.z[k]<vD.data[i+1]) A2.ordering->cardinal[i*k]=0;
        else A2.ordering->cardinal[i*k]=card;
        idnonzero++;
        }
      else if(i==grid.nx-1){
        if(grid.z[k]<vD.data[i-1]) A2.ordering->cardinal[i]=0;
        else A2.ordering->cardinal[i]=card;
        idnonzero++;
        }
      else{
        if(grid.z[k]<vD.data[i-1] && grid.z[k]<vD.data[i+1]) A2.ordering->cardinal[i]=0;
        else A2.ordering->cardinal[i]=card;
        idnonzero++;
        }
      }
    }  
  
  A2.ordering->allocate_finalize();
  A2.allocate();
  idnonzero=0;
  idrows=0;
  

  for(i=0;i<grid.nx;i++){
    for(k=0;k<grid.nz;k++) {
      if(grid.z[k]<vD.data[i]) continue;
      left = 0;
      right = 0;
      
      if(i==0){
        if(grid.z[k]<vD.data[i+1]){
          idrows++; 
          continue;
          }
        else boundary = LEFT_BOUNDARY;
//         else boundary = 'l';
        for(n=0;n<grid.nz;n++){
          if(grid.z[n]<vD.data[i]){
            n++;
            break;
            }
          }
        left = n;
        }
      else if(i==grid.nx-1){
        if(grid.z[k]<vD.data[i-1]){
          idrows++; 
          continue;
          }
        else boundary = RIGHT_BOUNDARY;
//         else boundary = 'r';
        for(n=0;n<grid.nz;n++){
          if(grid.z[n]<vD.data[i-1]){
            n++;
            break;
            }
          }
        right = n;
        }
      else{
        if(grid.z[k]<vD.data[i-1] && grid.z[k]<vD.data[i+1]){
          idrows++; 
          continue;
          }
        
        for(n=0;n<grid.nz;n++){
          if(grid.z[n]<vD.data[i]){
            n++;
            break;
            }
          }
        left = n;
        for(n=0;n<grid.nz;n++) {
          if(grid.z[n]<vD.data[i-1]) {
            n++;
            break;
            }
          }
        right = n; 
        
        if(grid.z[k]<vD.data[i-1]){
          boundary = LEFT_BOUNDARY;
//           boundary = 'l';
          }
        else if(grid.z[k]<vD.data[i+1]){
          boundary = RIGHT_BOUNDARY;
//           boundary = 'r';
          }
        else{
          boundary = NO_BOUNDARY;
//           boundary = 'n';
          }
        }
      switch(boundary){
        case NO_BOUNDARY :
//         case 'n' :
          A2.ordering->pointer[idrows]=idnonzero-left;
          A2.ordering->incidence[idnonzero]=idnonzero-left;
          A2.packed[idnonzero]= -1/(2*grid.dx);
          idnonzero++;
          A2.ordering->incidence[idnonzero]=idnonzero+right;
          A2.packed[idnonzero]= 1/(2*grid.dx);
          idnonzero++;
          idrows++;
          break;
        case LEFT_BOUNDARY :
//         case 'l' :
          A2.ordering->pointer[idrows]=idnonzero;
          A2.ordering->incidence[idnonzero]=idnonzero;
          A2.packed[idnonzero]= -1/grid.dx;
          idnonzero++;
          A2.ordering->incidence[idnonzero]=idnonzero+right;
          A2.packed[idnonzero]= 1/grid.dx;
          idnonzero++;
          idrows++;
          break;
        case RIGHT_BOUNDARY :
//         case 'r' :
          A2.ordering->pointer[idrows]=idnonzero-left;
          A2.ordering->incidence[idnonzero]=idnonzero-left;
          A2.packed[idnonzero]= -1/grid.dx;
          idnonzero++;
          A2.ordering->incidence[idnonzero]=idnonzero;
          A2.packed[idnonzero]= 1/grid.dx;
          idnonzero++;
          idrows++;
          break;
        default :
          printf("\nboundary not define for the node i=%i ; k=%i\n",i,k);
          goto error;
        }
      if(idnonzero > nrows*card){
        printf("\nnon-zero index value too big at i=%i ; k=%i\n",i,k);
        goto error;
        }
      }
    }
  
  status = matrix_check_Zero(A2);


// Create matrix A from A1 and A2 :
  nparams=grid.nz*grid.nx-nmaskparams;
  nrows=nparams+nprofiles*grid.nz-nmaskobs;
  A.ordering=new ordering_t;
  A.ordering->allocate(nrows);  
  for(i=0;i<nrows;i++){
      A.ordering->cardinal[i]=A1.ordering->cardinal[i];
    }
  for(i=0;i<nparams;i++){
      A.ordering->cardinal[nrows+i]=A2.ordering->cardinal[i];
    }
  
  A.ordering->allocate_finalize();
  A.allocate();
  idnonzero=0;
  idrows=0;

  status=ordering_Vconcat(A1.ordering,A2.ordering,A.ordering);
  
  printf("Done\n");
  printf("Inversion -- ");
  
//produit matriciel                   status=matrix_product(A, M, AM, 1);
//transposée                          status=matrix_CSR2CSC(A1)
//Resolution d'un systeme matriciel   LinearSystem_solve(hypermatrix_t & M, double *B)
//Concatenate matrix                  ?  
  printf("Done\n");
  printf("Saving...\n");
  
  printf("-----------------------------\n");
  printf("-----------------------------\n");
  printf("        Exiting\n");
  
  return(status);

 error:
  return(-1);
}
