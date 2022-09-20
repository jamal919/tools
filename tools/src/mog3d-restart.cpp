
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

#define MAIN_SOURCE

#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "tools-structures.h"
     
#include "fe.h"
#include "map.h"
#include "netcdf-proto.h"

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int load_state(const char *datafile, const char *gridfile, grid_t *grid, const char *varname, float **buffer, float *mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  float fmask;
  double t,*time,duration;
  date_t actual,origine;
  int k,n,i;
  int frame,nvalues,n_spacedims,status;
  cdfgbl_t global;
  cdfatt_t aB;//attribute buffer
  cdfvar_t vtime, info;
  decoded_t decoded;
  int nlayers;

  status= cdf_globalinfo(datafile,&global,0);
  status= cdf_varinfo  (datafile,varname,&info,0);
  
  status= poc_decode_axis(info, global, &decoded);
  status= poc_decode_mask(info, &(decoded));
  
  nvalues=1;
  n_spacedims=0;
  switch (decoded.xlen) {
    case 0:
      break;
    default:
      nvalues*=decoded.xlen;
      n_spacedims++;
      break;
    }
  switch (decoded.ylen) {
    case 0:
      break;
    default:
      nvalues*=decoded.ylen;
      n_spacedims++;
      break;
    }
  switch (decoded.zlen) {
    case 0:
      break;
    default:
      nvalues*=decoded.zlen;
      n_spacedims++;
      break;
    }

  *buffer=new float[nvalues];
  
  frame=0;
  
  switch(n_spacedims) {
    case 2:
      status= poc_getvar2d (datafile, info.id, frame,(float *) *buffer, &fmask ,info);
      break;
    case 3:
      status= poc_getvar3d (datafile, info.id, frame,(float *) *buffer, &fmask ,info);
      break;
    default:
      status=-1;
      break;
    }

  switch(n_spacedims) {
    case 2:
      status= poc_getgrid2d (gridfile, global, info, grid);
      break;
    case 3:
/* *----------------------------------------------------------------------------
      get grid once and for all (assumes variables have same discretisation !!!) */
      int vdepth;
      status=poc_getgrid3d (gridfile,  global, info, grid, &vdepth);
      break;
    default:
      status=-1;
      break;
    }
  
  *mask=fmask;
  return(status);
}

// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
//
// int fe_e2edges (mesh_t *mesh)
//
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// /*------------------------------------------------------------------------
//   build the neigh list from the element list*/
// {
//   float dummy;
//   int i,j,k,m,n,maxn,ndum;
//   int nitems,n1,n2,n3,count;
//   char c;
//
//   for(n=0;n<mesh->nedges;n++) mesh->edges[n].nshared=0;
//
//   for(m=0;m<mesh->ntriangles;m++)
//     for(i=0;i<3;i++) {
//       n=mesh->triangles[m].edges[i];
//       mesh->edges[n].nshared+=1;
//       n1=mesh->triangles[m].vertex[(i+1)%3];
//       n2=mesh->triangles[m].vertex[(i+2)%3];
//       mesh->edges[n].extremity[0]=n1;
//       mesh->edges[n].extremity[1]=n2;
//       mesh->edges[n].lon=0.5*(mesh->vertices[n1].lon+mesh->vertices[n2].lon);
//       mesh->edges[n].lat=0.5*(mesh->vertices[n1].lat+mesh->vertices[n2].lat);
//       }
//
//   for(m=0;m<mesh->nquadrangles;m++)
//     for(i=0;i<4;i++) {
//       n=mesh->quadrangles[m].edges[i];
//       mesh->edges[n].nshared+=1;
//       n1=mesh->quadrangles[m].vertex[(i+1)%3];
//       n2=mesh->quadrangles[m].vertex[(i+2)%3];
//       mesh->edges[n].extremity[0]=n1;
//       mesh->edges[n].extremity[1]=n2;
//       mesh->edges[n].lon=0.5*(mesh->vertices[n1].lon+mesh->vertices[n2].lon);
//       mesh->edges[n].lat=0.5*(mesh->vertices[n1].lat+mesh->vertices[n2].lat);
//       }
//
//   return(0);
// }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int fe_readlevels(const char *dataname, mesh_t mesh, int discretisation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,n;
  int status,frame=0;
  cdfvar_t info;
  cdfgbl_t gblinfo;
  int input_id;
  int varid=-1;
  discretisation_t *descriptor;
  float *z=0,mask;
//  extern int fe_read3d(int ncid, int varid, int frame, float *buffer, float *mask);

  switch(discretisation) {
    case LGP0:
      descriptor=&(mesh.LGP0descriptor);
      varid=cdf_identify(dataname,"z-LGP0");
      break;

    case LGP1:
      descriptor=&(mesh.LGP1descriptor);
      varid=cdf_identify(dataname,"z-LGP1");
      break;

    case D_LGP1:
      descriptor=&(mesh.DGP1descriptor);
      varid=cdf_identify(dataname,"z-DGP1");
      break;

    case LGP2:
      descriptor=&(mesh.LGP2descriptor);
      varid=cdf_identify(dataname,"z-LGP2");
      break;

    case NCP1:
      descriptor=&(mesh.NCP1descriptor);
      break;

    default:
      return(-1);
      break;
    }
    
  if(varid!=-1) {
    if(descriptor->nodes==0) descriptor->nodes=new node_t[descriptor->nnodes];
    status=fe_read3d(dataname,varid,frame,&z,&mask);
//    status=poc_get_UG3D(dataname,varid,frame,&z,&mask);
    for(n=0;n<descriptor->nnodes;n++) {
      descriptor->nodes[n].zlevels=new double[mesh.nlevels];
      for(k=0;k<mesh.nlevels;k++) {
        descriptor->nodes[n].zlevels[k]=z[k*descriptor->nnodes+n];
        }
      }
    }

  return(0);
}

// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
//
// int fe_readgeometry(const char *filename,mesh_t *mesh)
//
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   int n,k,option,status, ncid,verbose=1;
//   double time;
//   int var,varid,dim;
//   double *x=0,*y=0,*z=0,*h=0,H;
//   float *sigma=0,*elevation=0,mask;
//   int *v;
//   variable_t details;
//   cdfgbl_t fileinfo;
//
//
//   status=nc_open(filename,0,&ncid);
//   if(status != NC_NOERR) goto error;
//
// /* *------------------------------------------------------------------------
//   inquire file informations */
//   status= cdf_globalinfo(ncid, &fileinfo, 0);
//
//   status=fe_init(mesh);
//
// /* *----------------------------------------------------------------------------
//   T-UGOM*/
//
//   dim=cdf_identify_dimension(fileinfo,"N");
//   if(dim!=-1) mesh->nvtxs=fileinfo.dimension[dim].length;
//
//   dim=cdf_identify_dimension(fileinfo,"M");
//   if(dim!=-1) mesh->ntriangles=fileinfo.dimension[dim].length;
//
//   dim=cdf_identify_dimension(fileinfo,"NQUADRANGLES");
//   if(dim!=-1) mesh->nquadrangles=fileinfo.dimension[dim].length;
//
//   dim=cdf_identify_dimension(fileinfo,"E");
//   if(dim!=-1) mesh->nedges=fileinfo.dimension[dim].length;
//
//   dim=cdf_identify_dimension(fileinfo,"K");
//   if(dim!=-1) mesh->nlayers=fileinfo.dimension[dim].length;
//
//   dim=cdf_identify_dimension(fileinfo,"L");
//   if(dim!=-1) mesh->nlevels=fileinfo.dimension[dim].length;
//
// /* *----------------------------------------------------------------------------
//   WW3*/
//   dim=cdf_identify_dimension(fileinfo,"node");
//   if(dim!=-1) mesh->nvtxs=fileinfo.dimension[dim].length;
//
//   dim=cdf_identify_dimension(fileinfo,"element");
//   if(dim!=-1) mesh->ntriangles=fileinfo.dimension[dim].length;
//
// /* *----------------------------------------------------------------------------
//   FVCOM*/
//   dim=cdf_identify_dimension(fileinfo,"node");
//   if(dim!=-1) mesh->nvtxs=fileinfo.dimension[dim].length;
//
//   dim=cdf_identify_dimension(fileinfo,"nele");
//   if(dim!=-1) mesh->ntriangles=fileinfo.dimension[dim].length;
//
//   dim=cdf_identify_dimension(fileinfo,"E");
//   if(dim!=-1) mesh->nedges=fileinfo.dimension[dim].length;
//
//   dim=cdf_identify_dimension(fileinfo,"siglay");
//   if(dim!=-1) mesh->nlayers=fileinfo.dimension[dim].length;
//
//   dim=cdf_identify_dimension(fileinfo,"siglev");
//   if(dim!=-1) mesh->nlevels=fileinfo.dimension[dim].length;
//
//   mesh->vertices = new vertex_t[mesh->nvtxs];
//   if(mesh->nedges!=0)       mesh->edges       = new edge_t[mesh->nedges];
//   if(mesh->ntriangles!=0)   mesh->triangles   = new triangle_t[mesh->ntriangles];
//   if(mesh->nquadrangles!=0) mesh->quadrangles = new quadrangle_t[mesh->nquadrangles];
//
//   x=new double[mesh->nvtxs];
//   y=new double[mesh->nvtxs];
//   h=new double[mesh->nvtxs];
//   z=0;
//   if(mesh->nlevels!=0)  z=new double[mesh->nvtxs*mesh->nlevels];
//
// /* *----------------------------------------------------------------------------
//   T-UGOM geometry*/
//   status=nc_inq_varid (ncid,"lon",&varid);
//   if(status!=0) {
// /* *----------------------------------------------------------------------------
//     WW3*/
//     status=nc_inq_varid (ncid,"longitude",&varid);
//     }
//   status=nc_get_var_double(ncid,varid,x);
//
//   status=nc_inq_varid (ncid,"lat",&varid);
//   if(status!=0) {
// /* *----------------------------------------------------------------------------
//     WW3*/
//     status=nc_inq_varid (ncid,"latitude",&varid);
//     }
//   status=nc_get_var_double(ncid,varid,y);
//
//   mesh->type=SPHERICAL;
//
// /* *----------------------------------------------------------------------------
//   T-UGOM bathymetry*/
//   status=nc_inq_varid     (ncid,"bathymetry",&varid);
//   if(status==NC_NOERR) {
//     status=nc_get_var_double(ncid,varid,h);
//     }
//
// /* *----------------------------------------------------------------------------
//   FVCOM bathymetry*/
//   status=nc_inq_varid     (ncid,"h",&varid);
//   if(status==NC_NOERR) {
//     status=nc_get_var_double(ncid,varid,h);
//     }
//
//   for (n=0;n<mesh->nvtxs;n++) {
//     mesh->vertices[n].lon=x[n];
//     mesh->vertices[n].lat=y[n];
//     }
//
//   if(x!=0) delete[] x;
//   if(y!=0) delete[] y;
//   if(z!=0) delete[] z;
//
//   free(h);
//
// /* *----------------------------------------------------------------------------
//   T-UGOM, former format*/
//   status=nc_inq_varid     (ncid,"element",&varid);
//   if(status==NC_NOERR) {
//     v=new int[mesh->ntriangles*3];
//     status=nc_get_var_int   (ncid,varid,v);
//     for (n=0;n<mesh->ntriangles;n++) {
//       for(k=0;k<3;k++) {
//         mesh->triangles[n].vertex[k]=v[n*3+k];
//         }
//       }
//     delete[] v;
//     }
//
// /* *----------------------------------------------------------------------------
//   T-UGOM, new format*/
//   status=nc_inq_varid     (ncid,"triangles",&varid);
//   if(status==NC_NOERR) {
//     v=new int[mesh->ntriangles*3];
//     status=nc_get_var_int   (ncid,varid,v);
//     for (n=0;n<mesh->ntriangles;n++) {
//       for(k=0;k<3;k++) {
//         mesh->triangles[n].vertex[k]=v[n*3+k];
//         }
//       }
//     delete[] v;
//     }
//
//   status=nc_inq_varid     (ncid,"quadrangles",&varid);
//   if(status==NC_NOERR) {
//     v=new int[mesh->nquadrangles*4];
//     status=nc_get_var_int   (ncid,varid,v);
//     for (n=0;n<mesh->ntriangles;n++) {
//       for(k=0;k<4;k++) {
//         mesh->quadrangles[n].vertex[k]=v[n*4+k];
//         }
//       }
//     delete[] v;
//     }
//
// /* *----------------------------------------------------------------------------
//   WW3*/
//   status=nc_inq_varid     (ncid,"tri",&varid);
//   if(status==NC_NOERR) {
//     v=new int[mesh->ntriangles*3];
//     status=nc_get_var_int   (ncid,varid,v);
//     for (n=0;n<mesh->ntriangles;n++) {
//       for(k=0;k<3;k++) {
//         mesh->triangles[n].vertex[k]=v[k+3*n]-1;
//         }
//       }
//     delete[] v;
//     }
//
// /* *----------------------------------------------------------------------------
//   FVCOM*/
//   status=nc_inq_varid     (ncid,"nv",&varid);
//   if(status==NC_NOERR) {
//     v=new int[mesh->ntriangles*3];
//     status=nc_get_var_int   (ncid,varid,v);
//     for (n=0;n<mesh->ntriangles;n++) {
//       for(k=0;k<3;k++) {
//         mesh->triangles[n].vertex[k]=v[n+mesh->ntriangles*k]-1;
//         }
//       }
//     delete[] v;
//     }
//
//   status = nc_inq_varid(ncid, "edges", &varid);
//   if(status==NC_NOERR) {
//     v=new int[mesh->ntriangles*3+mesh->nquadrangles*4];
//     status = nc_get_var_int(ncid, varid, v);
//     for(n = 0; n < mesh->ntriangles; n++) {
//       for(k = 0; k < 3; k++) {
//         mesh->triangles[n].edges[k] = v[n * 3 + k];
//         }
//       }
//     for(n = 0; n < mesh->nquadrangles; n++) {
//       for(k = 0; k < 4; k++) {
//         mesh->quadrangles[n].edges[k] = v[mesh->ntriangles*3 + n * 4 + k];
//         }
//       }
//     delete[] v;
//     }
//
// /* *----------------------------------------------------------------------------
//   build vertex neighbour list from element list*/
//   status= fe_e2n (mesh);
//
//   if(mesh->nedges==0) {
// /* *----------------------------------------------------------------------------
//     edge table need to be built from scratch*/
// //    status=build_edgetable(mesh);
//     status=init_edge_table(mesh);
//     }
//   else {
// /* *----------------------------------------------------------------------------
//     use file's edge table*/
//     status= fe_e2edges (mesh);
//     }
//   for(n=0; n<mesh->ntriangles; n++) {
//     status=fe_initaffine(mesh,n);
//     if(mesh->triangles[n].Area <=0.0) {
//       printf("element %d has negative area (cw order)\n",n);
//       }
//     }
//
//   status = nc_close(ncid);
//   return(status);
//
// error:
//   return(status);
// }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int initialize_LGP1(mesh_t mesh, grid_t Tgrid, grid_t Sgrid, float *buf[2], float spec[2])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int nndes,status,fmt,verbose=1;
  int i,j,k,l,m,n;
  int nitems;
  int ncid,varid;
  char *discretisation=NULL;
  char *root,output[1024];
  grid_t slice_grid;
  variable_t varinfo;
  int frame;
  float *T,  *S;
 
  double x,y,t,p;
  float *profile_x[2]={NULL,NULL}, *profile_y[2]={NULL,NULL};
  int nloc;

  int nlevels=mesh.nlevels,vdepth;
  float level[100];

  cdfgbl_t gT_info,gS_info;
  cdfvar_t vT_info,vS_info;
  cdfvar_t Tvar,Svar;
  
  discretisation_t descriptor=mesh.LGP1descriptor;

//   T=new float[descriptor.nnodes*descriptor.nlevels];
//   S=new float[descriptor.nnodes*descriptor.nlevels];
  T=new float[descriptor.nnodes*nlevels];
  S=new float[descriptor.nnodes*nlevels];
  
  profile_x[0]=new float[Tgrid.nz];
  profile_x[1]=new float[Tgrid.nz];
  profile_y[0]=new float[Tgrid.nz];
  profile_y[1]=new float[Tgrid.nz];
  
  for (n=0;n<descriptor.nnodes;n++) {
/*-----------------------------------------------------------------------
    interpolate profiles at node position */
//     mesh.vertices[n].sigma=(double *) malloc(nlevels*sizeof(double));
//     for (k=0;k<nlevels;k++) {
//       mesh.vertices[n].sigma[k]=level[k];
//       }
 
/*-----------------------------------------------------------------------
    create temperature profile at database levels */
    t=descriptor.nodes[n].lon;
    p=descriptor.nodes[n].lat;
    status=map_profile02(Tgrid,t,p, buf[0],spec[0],&(profile_x[0]), &(profile_y[0]), &nloc);
    if(status!=0) {
/*       printf("node %d: t=%lf p=%lf isolated, try neighbours\n",n,t,p); */
/*-----------------------------------------------------------------------
      create profile at database levels from nearest neighbour*/
      for(i=0;i< mesh.vertices[n].nngh;i++) {
        m=mesh.vertices[n].ngh[i];
        t=descriptor.nodes[m].lon;
        p=descriptor.nodes[m].lat;
        status=map_profile02(Tgrid,t,p, buf[0],spec[0],&(profile_x[0]), &(profile_y[0]), &nloc);
        if(status==0) break;
        }
      }
    if(status!=0) {
      printf("node %d: t=%lf p=%lf isolated, no profile available\n",n,t,p);
      }

/*-----------------------------------------------------------------------
    create salinity profile at database levels */
    t=descriptor.nodes[n].lon;
    p=descriptor.nodes[n].lat;
    status=map_profile02(Sgrid,t,p, buf[1],spec[1],&profile_x[1], &profile_y[1], &nloc);
    if(status!=0) {
/*-----------------------------------------------------------------------
      create profile at database levels from nearest neighbour*/
      for(i=0;i< mesh.vertices[n].nngh;i++) {
        m=mesh.vertices[n].ngh[i];
        t=descriptor.nodes[m].lon;
        p=descriptor.nodes[m].lat;
        status=map_profile02(Tgrid,t,p, buf[0],spec[0],&profile_x[0], &profile_y[0], &nloc);
        if(status==0) break;
        }
      }
  
/*-----------------------------------------------------------------------
    create profile at mesh levels */
    for (k=0;k<nlevels;k++) {
      status=map_interpolate1D(profile_x[0],profile_y[0],spec[0],nloc,(float) descriptor.nodes[n].zlevels[k], &T[k*descriptor.nnodes+n],1);
      status=map_interpolate1D(profile_x[1],profile_y[1],spec[1],nloc,(float) descriptor.nodes[n].zlevels[k], &S[k*descriptor.nnodes+n],1);
      }
//     for (k=0;k<nlevels;k++) {
//       T[k*descriptor.nnodes+n]=profile_x[0][k];
//       S[k*descriptor.nnodes+n]=profile_x[1][k];
//       }
/*-----------------------------------------------------------------------
    get bathymetry from climatology */
//     for (k=0;k<nlevels;k++) {
//       if (T[k*mesh.nvtxs+n]!=spec[0]) {
// 	descriptor.nodes[n].h=level[k];
// 	break;
//         }
//       }
    }

//   for (n=0;n<descriptor.nnodes;n++) {
//      for (k=0;k<nlevels;k++) {
//        T[k*descriptor.nnodes+n]=1;
//        S[k*descriptor.nnodes+n]=1;
//        }
//     }
  mesh.nlayers=nlevels;
  sprintf(output,"%s","initialisation-TS.nc");
  status= fe_savemeshNC3D(output, mesh, 1);

  status= fe_savediscretisation(output, mesh, descriptor);

//   Tvar=poc_floatvariable_nz("T", spec[0],"C",1.,0.,"water temperature");
//   Svar=poc_floatvariable_nz("S", spec[1],"psu",1.,0.,"water salinity");
//
//   status=create_ncvariable(output, &Tvar);
//   status=create_ncvariable(output, &Svar);
//
//   status= poc_putfloat_nz(output, mesh, Tvar.id, T);
//   status= poc_putfloat_nz(output, mesh, Svar.id, S);
  float scale=1., offset=0;
  Tvar = poc_variable_UG4D("T", spec[0], "C", scale, offset,"potential_sea_water_temperature","T","K","N");
  status   = cdf_createvariable(output, &(Tvar));
  Svar = poc_variable_UG4D("S",  spec[1], "PSU", scale, offset,"sea_water_sanility","T","K","N");
  status   = cdf_createvariable(output, &(Svar));

  int count=0;
  status   = poc_put_UG4D(output, mesh, count, Tvar.id, T);
  Tvar.destroy();

  status   = poc_put_UG4D(output, mesh, count, Svar.id, S);
  Svar.destroy();

  delete[] T;
  delete[] S;

error:
  __ERR_BASE_LINE__("exiting\n");exit(-1);

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int nndes,status,*elts,fmt,fefmt,verbose=1;
  int i,j,k,l,L,m,n;
  int nitems;
  int ncid,varid;
  char *keyword,*zone;
  char *meshfile=NULL,*inputT=NULL,*inputS=NULL,*discretisation=NULL,*gridfile=NULL;
  char *root,output[1024];
  grid_t Tgrid,Sgrid,slice_grid;
  mesh_t mesh,nodes;
  variable_t varinfo;
  int frame;
  float *buf[2],spec[2];
 
  double x,y,t,p;
  float *profile_x[2]={NULL,NULL}, *profile_y[2]={NULL,NULL};
  int nloc;
  float *T,*S;

  int nlevel=33,vdepth;
  float level[100];

  cdfgbl_t gT_info,gS_info;
  cdfvar_t vT_info,vS_info;
  cdfvar_t Tvar,Svar;

  fct_echo(argc,argv);
 
  n=1;
  while (n < argc) {
    keyword=strdup(argv[n]);
    switch (keyword[0]) {
      case '-':
        switch (keyword[1]) {
        case 'd' :
          discretisation= strdup(argv[n+1]);
          printf("discretisation=%s\n",discretisation);
          n++;
          n++;
          break;

        case 'm' :
          meshfile= strdup(argv[n+1]);
          printf("meshfile=%s\n",meshfile);
          n++;
          n++;
          break;

        case 'g' :
          gridfile= strdup(argv[n+1]);  /* directory */
          n++;
          n++;
          break;

        default:
          __OUT_BASE_LINE__("unknown option %s\n",keyword);
          exit(-1);
        }
        break;

      default:
        if(inputT==NULL) {
          inputT= strdup(argv[n]);
          printf("input file=%s\n",inputT);
          n++;
          }
        if(inputS==NULL) {
          inputS= strdup(argv[n]);
          printf("input file=%s\n",inputS);
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

  if(inputT==NULL) {
//    inputT= strdup("/home/ocean/softs/genesis/sparc3/data/climatology/clim.med.temp.nc");
    inputT= strdup("/data1/models/sirocco/dev-49/data-v1/model/20090110_000036_domconcat.nc");
    }
  if(inputS==NULL) {
//    inputS= strdup("/home/ocean/softs/genesis/sparc3/data/climatology/clim.med.psal.nc");
    inputS= strdup("/data1/models/sirocco/dev-49/data-v1/model/20090110_000036_domconcat.nc");
    }
    
  if(gridfile==NULL) {
    gridfile= strdup("/data1/models/sirocco/dev-49/data-v1/grille_domconcat.nc");
    }
    
  printf("temperature input file=%s\n",inputT);
  printf("temperature input file=%s\n",inputS);

/* *-----------------------------------------------------------------------
  load mesh */
//   status=fe_init(&mesh);
//   if(status !=0) goto error;
//
//   status=fe_readmesh(meshfile,MESH_FILE_FORMAT_TRIGRID,&mesh);
//   if(status !=0) goto error;
//
//   status=fe_list(&mesh);
//   if(status !=0) goto error;
  status=fe_readgeometry(meshfile,&mesh);
  status=fe_geometry(&mesh);
    
  fe_integrale_init();

  fmt=0;
//  status=fe_readdiscretisation(meshfile,&mesh,fmt);
  
 status=discretisation_init(&mesh, LGP1);
 status=fe_readlevels(meshfile, mesh, LGP1);
 
/* *-----------------------------------------------------------------------
  load climatology */
  status=load_state(inputT, gridfile, &Tgrid, "tem", &(buf[0]), &spec[0]);
  if(status != 0) goto error;
  status=load_state(inputS, gridfile, &Sgrid, "sal", &(buf[1]), &spec[1]);
  if(status != 0) goto error;

/*-----------------------------------------------------------------------
  Z level climatology ?*/
  status= map_chklevels(Tgrid);
  status= map_chklevels(Sgrid);

//   nlevel=Tgrid.nz;
//   for (k=0;k<nlevel;k++) {
//     level[k]=-Tgrid.z[k*Tgrid.nx*Tgrid.ny];
//     }

  status=map_diffusion01(Tgrid, buf[0],spec[0]);
  status=map_diffusion01(Sgrid, buf[1],spec[1]);

  status=map_diffusion01(Tgrid, buf[0],spec[0]);
  status=map_diffusion01(Sgrid, buf[1],spec[1]);
  
  status=initialize_LGP1( mesh, Tgrid, Sgrid, buf, spec);
  
end: printf("end of  ... \n");

error:
  __ERR_BASE_LINE__("exiting\n");exit(-1);

}

