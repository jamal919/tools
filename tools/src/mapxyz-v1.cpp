
/**************************************************************************

  T-UGO tools, 2006-2009

  Unstructured Ocean Grid initiative

Contributors:

  Florent Lyard      LEGOS/CNRS, Toulouse, France
  Yoann Le Bars      LEGOS, Toulouse, France (PhD)
  Laurent Roblou     LEGOS, Toulouse, France
  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
  Frédéric Dupont    Université de Laval à Québec, Canada

E-mail: florent.lyard@legos.obs-mip.fr

***************************************************************************/

/* *@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  rectify aims to re-format ifremer bathymetric database.

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */

#define MAIN_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include "tools-structures.h"
#include "constants.h"

#include "rutin.h"
#include "geo.h"
#include "polygones.h"
#include "netcdf-proto.h"
#include "grd.h"
#include "map.h"
#include "filter.h"
#include "statistic.h"
#include "functions.h"
#include "sym-io.h"
#include "xyz.h"

#define XYZ  0
#define YXZ  1
#define NYXZ 3


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_fill(const grid_t grid, float *buf, float mask, const float scale, float *out, int niteration)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  Smoothing by diffusion operator:
  --------------------------------
  
    dC/dt = nu (d²C/dx²+d²C/dy²)
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int status=0;
  int j,k, l;
  double zz;
  const double R=6.3e+06,dx=R*grid.dx*d2r;
  const double nu=scale;
  double *state[3],*swap,dmask=mask;
  const double L=2*M_PI*sqrt(nu);
  const double CFL=1./(4.*nu/dx/dx);
//  const float tau=180.;
  const double tau=CFL/1.13;
  const double coef=nu*tau/(dx*dx);
  
  for(k=0;k<2;k++) state[k]=new double [grid.nx*grid.ny];
  
  printf("\ntypical lengths=%f m smoothing lengths=%f km\n",L,L*sqrt(tau*niteration)/1000.0);
  printf("CFL=%f s time step=%f s\n\n",CFL,tau);

  int nRequestedProcs=-1;
  int nprocs=initialize_OPENMP(nRequestedProcs);
  
#pragma omp parallel for if(nprocs>1)
  for(j=0;j<grid.ny;j++) {
    int i;
    for(i=0;i<grid.nx;i++) {
      const int m=j*grid.nx+i;
      state[0][m]=(double) buf[m];
      state[1][m]=(double) buf[m];
      if(out!=buf) out[m]=buf[m];
      }
    }
  for(int iteration=0;iteration<niteration;iteration++) {
    printf("iteration %6d\n",iteration);
#pragma omp parallel for private(zz) if(nprocs>1)
    for(j=1;j<grid.ny-1;j++) {
      int i;
      for(i=1;i<grid.nx-1;i++) {
        const int m=j*grid.nx+i;
        int done=0;
        state[1][m]=state[0][m];
        if((state[0][m]!=mask)&&(out[m]==mask)) {
          done=0;
          const int m1=m-1;
          const int m2=m+1;
          const int n1=m-grid.nx;
          const int n2=m+grid.nx;
          const int frame=0;
          const double z1= state[frame][m1];
          if(z1==dmask) continue;
          const double z2= state[frame][m2];
          if(z2==dmask) continue;
          const double z3= state[frame][n1];
          if(z3==dmask) continue;
          const double z4= state[frame][n2];
          if(z4==dmask) continue;
          zz=z1+z2+z3+z4-4.* state[frame][m];
          state[1][m]=state[0][m]+coef*zz;
//           state[1][m]=state[0][m];
//           done=1;
// //           const int mm1=m-grid.nx-1;
// //           const int mm2=m+grid.nx+1;
// //           const int nn1=m-grid.nx+1;
// //           const int nn2=m+grid.nx-1;
// //           const double zz1= state[frame][mm1];
// //           if(zz1==dmask) continue;
// //           const double zz2= state[frame][mm2];
// //           if(zz2==dmask) continue;
// //           const double zz3= state[frame][nn1];
// //           if(zz3==dmask) continue;
// //           const double zz4= state[frame][nn2];
// //           if(zz4==dmask) continue;
// //           zz=zz1+zz2+zz3+zz4-4.* state[frame][m];
// //           state[1][m]+=coef*zz/2.0;
          }

//        if(done==0) {
        if((state[0][m]!=mask)&&(out[m]!=mask)) {
          state[1][m]=state[0][m];
          }
        if(state[0][m]==mask) {
          double w=0.0;
          zz=0;
          const int m1=m-1;
          const int m2=m+1;
          const int n1=m-grid.nx;
          const int n2=m+grid.nx;
          const int frame=0;
          const double z1= state[frame][m1];
          if(z1!=dmask) {
            zz+=z1;
            w++;
            }
          const double z2= state[frame][m2];
          if(z2!=dmask) {
            zz+=z2;
            w++;
            }
          const double z3= state[frame][n1];
          if(z3!=dmask) {
            zz+=z3;
            w++;
            }
          const double z4= state[frame][n2];
          if(z4!=dmask) {
            zz+=z4;
            w++;
            }
          const int mm1=m-grid.nx-1;
          const int mm2=m+grid.nx+1;
          const int nn1=m-grid.nx+1;
          const int nn2=m+grid.nx-1;
          const double zz1= state[frame][mm1];
          if(zz1!=dmask) {
            zz+=zz1;
            w++;
            }
          const double zz2= state[frame][mm2];
          if(zz2!=dmask) {
            zz+=zz2;
            w++;
            }
          const double zz3= state[frame][nn1];
          if(zz3!=dmask) {
            zz+=zz3;
            w++;
            }
          const double zz4= state[frame][nn2];
          if(zz4!=dmask) {
            zz+=zz4;
            w++;
            }
          if(state[frame][m]!=dmask) {
            zz+=state[frame][m];
            w++;
            }
          if(w!=0) state[1][m]=zz/w;
          else state[1][m]=mask;
          }
        }
      }
    swap=state[0];
    state[0]=state[1];
    state[1]=swap;
    }
    
#pragma omp parallel for if(nprocs>1)
  for(j=0;j<grid.ny;j++) {
    int i;
    for(i=0;i<grid.nx;i++) {
      const int m=j*grid.nx+i;
      out[m]=state[0][m];
      }
    }
    
  for(k=0;k<2;k++) delete[] state[k];
  
  return(status);
}



/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int xyz_map (double *x, double *y, double *z, int ndata, grid_t grid, float *buffer, float mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  Assume full (no missing values) array structured (sequential array's axes
  ordering) description

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int i,j,k,l,m,n,range;
  int status;
  int *count;

  status=0;

/*------------------------------------------------------------------------------
  count number of points in grid cell*/
  count=new int[grid.nx*grid.ny];
  for (j=0;j<grid.ny;j++) {
    for (i=0;i<grid.nx;i++) {
      m=j*grid.nx+i;
      count[m]=0;
      buffer[m]=0;
      }
    }

  for(n=0;n<ndata;n++) {
    status=map_index(grid,x[n],y[n],&i,&j);
    if(status==0) {
      m=j*grid.nx+i;
      buffer[m]+=z[n];
      count[m]++;
      }
    }

  for (j=0;j<grid.ny;j++) {
    for (i=0;i<grid.nx;i++) {
      m=j*grid.nx+i;
      if(abs(buffer[m])>1.e+6) {
        buffer[m]/=count[m];
        }
      if(count[m]!=0) {
        buffer[m]/=count[m];
        }
      else {
        buffer[m]=mask;
        }
      }
    }

  delete[] count;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int xyz_loadmap (const char *filename, grid_t grid, char *proj4_options, float *buf, float *mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  FILE *in;
  int downward=0,cartesian=0,mode=0,complete=0,option,status;
  int n,ndata,nitems,ntoken;
  char *masked=0;
  char line[1024],options[16],*s,*tmp;
  double *x,*y,*z,t,p;
  double dmask;

  *mask=-99999.0;

  string header="";
  status=xyz_loadraw (filename, header, proj4_options, x, y, z, &dmask, ndata);
  
  *mask=dmask;
  
  status=xyz_map (x, y, z, ndata, grid, buf, *mask);

  delete[] x;
  delete[] y;
  delete[] z;
  return(status);
  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  geo_t projection;
  float  topo_mask=1.e+10,mask=1.e+10;
  double  x,y,t,p;
  double  dx,dy;

  int i,j,k,l,m,n,status;

  size_t size;
  FILE *in,*out;
  char *keyword;
  int flag;

  char *rootname=NULL,*output=NULL,*input=NULL,*pathname=NULL,*bathymetry=NULL;
  char *bathy=NULL,*notebook=NULL,*meshbook=NULL,*poly=NULL,*zone=NULL;
  char *proj4=NULL;

  grid_t topogrid,cgrid;
  grid_t grid;
  mesh_t mesh;
  short *topo,smask=256*127+255;;
  float *ftopo,*tmp,*buffer,scale=1.0;

  plg_t *polygones=NULL;
  int npolygones=0;

  pocgrd_t ncgrid;
  cdfvar_t variable;
  int dimlgth[4]={10,20,30,40};

  fct_echo( argc, argv);

  n=1;
  while (n < argc) {
    keyword=strdup(argv[n]);
    if(strcmp(keyword,"-inv")==0) {
      scale=-1;
      n++;
      continue;
      }
    switch (keyword[0]) {
      case '-':
        switch (keyword[1]) {
        case 'b' :
          bathymetry= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'n' :
          notebook= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'p' :
          poly= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'o' :
          rootname= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'z' :
          zone= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case '-' :
          if(strstr(argv[n],"--proj=")!=0){
            proj4=strdup(argv[n]+7);
            n++;
            }
          break;

        default:
          __OUT_BASE_LINE__("unknown option %s\n",keyword);
          exit(-1);
        }
        break;

      default:
        if(input==NULL) {
          input= strdup(argv[n]);
          n++;
          }
        else {
          __OUT_BASE_LINE__("unknown option %s\n",keyword);
          exit(-1);
          }
        break;
      }
      free(keyword);
    }

 mask=9999.9;

 if(notebook != NULL) {
    printf("#################################################################\n");
    printf("load notebook file : %s\n",notebook);
    status=load_notebook(notebook, &cgrid, &topogrid, &projection);
    printf("%s (notebook file) processed\n",notebook);
    buffer=new float[topogrid.nx*topogrid.ny];
    status= xyz_loadmap (input, topogrid, proj4, buffer, &mask);
    }
  else if(zone != NULL) {
/*------------------------------------------------------------------------------
    Apply zone definition by name*/
    topogrid=get_zonegrid(zone);
    status=map_completegridaxis(&topogrid,0);
    buffer=new float[topogrid.nx*topogrid.ny];
    status= xyz_loadmap (input, topogrid, proj4, buffer, &mask);
    }
  else {
/*------------------------------------------------------------------------------
    read xyz file and interprete grid */
    status=xyz_loadmap(input,&topogrid);
    if(status !=0) {
      __OUT_BASE_LINE__("cannot load grid in bathymetry file=%s\n",input);
      exit(-1);
      }
    buffer=new float[topogrid.nx*topogrid.ny];
    for (j=0;j<topogrid.ny;j++) {
      for (i=0;i<topogrid.nx;i++) {
        m=j*topogrid.nx+i;
        buffer[m]=(float) topogrid.z[m];
        }
      }
    }

/*------------------------------------------------------------------------------
  check */
  for(m=0;m<topogrid.ny*topogrid.nx;m++) {
    if(isnan(buffer[m])) {
      printf("Nan value out of mask in topography\n");
      }
    }
  for (j=0;j<topogrid.ny;j++) {
    for (i=0;i<topogrid.nx;i++) {
      m=j*topogrid.nx+i;
      if(buffer[m]!=mask) {
        buffer[m]*=scale;
        }
      }
    }

/*------------------------------------------------------------------------------
  save rectified topo */
  if(rootname==0) rootname=strdup("test");
  output=new char[1024];
  sprintf(output,"%s.spherical.nc",rootname);

  printf("#################################################################\n");
  printf("save bathymetry file : %s\n",output);
  status= poc_createfile(output);
  
  
  if(topogrid.modeH==2) {
    int imid=topogrid.nx/2;
    int jmid=topogrid.ny/2;
    dx=topogrid.x[topogrid.nx*jmid+imid+1]-topogrid.x[topogrid.nx*jmid+imid-1];
    dy=topogrid.y[topogrid.nx*jmid+imid+1]-topogrid.y[topogrid.nx*jmid+imid-1];
    topogrid.dx=sqrt(dx*dx+dy*dy);
    dx=topogrid.x[topogrid.nx*(jmid+1)+imid]-topogrid.x[topogrid.nx*(jmid-1)+imid];
    dy=topogrid.y[topogrid.nx*(jmid+1)+imid]-topogrid.y[topogrid.nx*(jmid-1)+imid];
    topogrid.dy=sqrt(dx*dx+dy*dy);
    }
    
  status=poc_sphericalgrid_xy(output,"",topogrid,&ncgrid);
  poc_standardvariable_xy(&variable,"bathymetry",mask,"m",1., 0.,"bathymetry","bathymetry","bathymetry",ncgrid);
  status=create_ncvariable(output, &variable);
  status=poc_write_xy(output,  topogrid, variable.id,buffer);
  variable.destroy();

  status=map_fill(topogrid,buffer, mask, 1.0, buffer,20);

//  status=poc_sphericalgrid_xyzt(output,"","",topogrid,&ncgrid);
  poc_standardvariable_xy(&variable,"bathymetry-no-gap",mask,"m",1., 0.,"bathymetry","bathymetry","bathymetry",ncgrid);
  status=create_ncvariable(output, &variable);
  status=poc_write_xy(output,  topogrid, variable.id,buffer);
  variable.destroy();

  __OUT_BASE_LINE__("bathymetry sucessfully completed\n");

  exit(0);

 error:
  __ERR_BASE_LINE__("exiting\n");exit(-1);
  
}
