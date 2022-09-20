
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

#include "tools-define.h"
#include "tools-structures.h"

#include "functions.h"
#include "fe.h"
#include "bmg.h"
#include "grd.h"
#include "map.h"
#include "ascii.h"
#include "archive.h"
#include "geo.h"
#include "sym-io.h"
#include "netcdf-proto.h"


extern  int surfref_loadr1(char *name, int nvalues, int ncolumns, int pos, int *translation, float *buffer);
extern  int *surfref_reference(char *name, mesh_t mesh, int *ncolumns);

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  float *rbuf;
  float *rbufx,*rbufy;
  int *translation, ncolumns;
  float  dummy,rmask,*rbuffer;
  float  *rbufferx,*rbuffery;
  float  spec,scale_factor=1.0;
  int status,*elts,fmt,fefmt;
  int i,j,k,l,L,m,n;
  int nitems;
  FILE *file;
  FILE *out;
  char *meshfile=NULL,*sfile=NULL,*cfile=NULL,*format=NULL,*cdl=0;
  char *keyword,*zone=NULL,*gridfile=NULL,*notebook=NULL;
  char file1[256],file2[256],vname[256];
  char *root,tmp[256],output[1024];
  grid_t cgrid,grid;
  fcomplex *cbuffer,*cbuf,cmask;
  fcomplex *cbufferx,*cbuffery,*cbufx,*cbufy;
  mesh_t mesh;
  geo_t projection;
  cdfvar_t variable;
  pocgrd_t ncgrid;
  string cmd;
  char *comments[2];

  comments[0]=new char[1024];
  comments[1]=new char[1024];

  fct_echo( argc, argv);

  spec=999.9;

  n=1;
  while (n < argc) {
    keyword=strdup(argv[n]);
    if(strcmp(keyword,"-scale")==0) {
      sscanf(argv[n+1],"%f",&scale_factor);
      n++;
      n++;
      continue;
      }
    switch (keyword[0]) {
      case '-':
        switch (keyword[1]) {
        case 'f' :
          format= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'g' :
          gridfile= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'n' :
          notebook= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'o' :
          strcpy(output,argv[n+1]);
          n++;
          n++;
          break;

        case 'm' :
          meshfile= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 's' :
          sscanf(argv[n+1],"%f",&spec);
          n++;
          n++;
          break;

        case 'c' :
          cdl= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'z' :
          zone= strdup(argv[n+1]);
          n++;
          n++;
          break;

        default:
          __OUT_BASE_LINE__("unknown option %s\n",keyword);
          exit(-1);
        }
        break;

      default:
        if(sfile==NULL) {
          sfile= strdup(argv[n]);
          printf("input file=%s\n",sfile);
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

  rmask=spec;
  cmask=fcomplex(spec,spec);

  if(format==NULL) format= strdup("bmg");

  fefmt=S2R;

  if(strstr(sfile,".s2r") != NULL)    fefmt=S2R;
  if(strstr(sfile,".s2c") != NULL)    fefmt=S2C;
  if(strstr(sfile,".v2r") != NULL)    fefmt=V2R;
  if(strstr(sfile,".v2c") != NULL)    fefmt=V2C;

  if(zone != NULL) {
    grid=get_zonegrid(zone);
    status=map_completegridaxis_2(&grid);
    }

 if(gridfile != NULL) {
   status=cdf_loadvargrid (gridfile,0,&grid);
//    cmd=string("ncdum -h ")+string(gridfile)+string(" > gridit.cdl");
//    system(cmd.c_str());
   }

 if(cdl != NULL) {
   cmd=string("ncgen -b -x ")+string(cdl)+string(" -o "+string(output));
   system(cmd.c_str());
   }

 if(notebook != NULL) {
    status=load_notebook(notebook, &cgrid, &grid, &projection);
    printf("%s (notebook file) processed\n",notebook);
    }


 if(meshfile != NULL) {
   status=fe_readmesh(meshfile,MESH_FILE_FORMAT_TRIGRID,&mesh);
   if(status !=0) goto error;
   status=fe_list(&mesh);
   if(status !=0) goto error;
   }
 else {
   printf("no mesh file specified; abort...\n");
   goto error;
   }

  //fe_initaffine(&mesh);

  elts=fe_scan_elements(mesh,grid,0);
  if(elts==NULL) goto error;

  translation=surfref_reference(sfile, mesh, &ncolumns);

  switch (fefmt) {
    case S2R:
      rbuffer=(float *)malloc(mesh.nvtxs*sizeof(float));
      rbuf=(float *)malloc(grid.nx*grid.ny*sizeof(float));
      for(k=0;k<ncolumns;k++) {
        status=surfref_loadr1(sfile, mesh.nvtxs, ncolumns, k, translation, rbuffer);
        if(status!=0) goto error;
        if(strrchr(sfile,'/')==0)
          sprintf(output,"%s.nc",sfile);
        else
          sprintf(output,"%s.nc",(char *) (strrchr(sfile,'/')+1));
        root=strstr(output,".dat");
        sprintf(root,"-%2.2d.s2r",k);
        sprintf(comments[0],"%c POC/Noveltis, extracted from %s",169,sfile);
        sprintf(comments[1],"column %d",k);
        status=quoddy_saver1(output, mesh.nvtxs, rbuffer,comments);
        status=fe_map(mesh,rbuffer,grid,elts,rbuf,rmask);
        if(status!=0) goto error;
        for(n=0;n<grid.nx*grid.ny;n++) {
          if(rbuf[n]!=rmask) {
            rbuf[n]*=scale_factor;
            }
          }
        if(strcmp(format,"ascii")==0) {
          strcpy(output,sfile);
          root=strstr(output,".s2r");
          sprintf(root,"%s",".asc");
          status=ascii_saver1(output,grid, rbuf, rmask,"%5.0f", grid.nx);
          }
        else if(strcmp(format,"bmg")==0) {
          strcpy(output,sfile);
          root=strstr(output,".s2r");
          sprintf(root,"%s",".bmg");
          status=bmg_saver1(output, 1, 1, 1, grid, rbuf, 0.,rmask);
          }
        else if(strcmp(format,"grd")==0) {

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  save data, to be replaced with topo_save. Example:
  
  status=topo_save(rootname, output, "netcdf", TOPOgrid, TOPOdata, TOPOmask, debug);
  
  topo_save will hande mirroring...

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

          strcpy(output,sfile);
          root=strstr(output,".s2r");
          sprintf(root,"%s",".grd");
          status=grd_mirror_r( grid, grid.nx, rbuf, rmask);
          status=grd_save(output, grid, grid.nx, rbuf, rmask);
          }
        else if(strcmp(format,"netcdf")==0) {
          if(strrchr(sfile,'/')==0)
            sprintf(output,"%s.nc",sfile);
          else
            sprintf(output,"%s.nc",(char *) (strrchr(sfile,'/')+1));
          root=strstr(output,".dat");
          sprintf(root,"%s",".nc");
          if(k==0) {
            status= poc_createfile(output);
            status=poc_sphericalgrid_xy(output,"",grid,&ncgrid);
//             for(n=0;n<grid.nx*grid.ny;n++) rbuf[n]=elts[n];
//             poc_standardvariable_xy(&variable,"element-index",-1.,"no_units",1., 0.,"element-index","element-index","element-index",ncgrid);
//             status=create_ncvariable(output, &variable);
//             status=poc_write_xy(output, grid, variable.id, rbuf);
//             status=free_ncvariable(&variable);
            }
          sprintf(vname,"variable%2.2d",k);
          poc_standardvariable_xy(&variable,vname,rmask,"no_units",1., 0.,"no_name","no_name","no_name",ncgrid);
          status=create_ncvariable(output, &variable);
          status=poc_write_xy(output, grid, variable.id, rbuf);
          variable.destroy();
          }
        }
      break;
    case V2R:
      rbufferx=(float *)malloc(mesh.nvtxs*sizeof(float));
      rbuffery=(float *)malloc(mesh.nvtxs*sizeof(float));
      rbufx=(float *)malloc(grid.nx*grid.ny*sizeof(float));
      rbufy=(float *)malloc(grid.nx*grid.ny*sizeof(float));
      status=quoddy_loadr2(sfile, mesh.nvtxs, rbufferx,rbuffery);
      status=fe_map(mesh,rbufferx,grid,elts,rbufx,rmask);
      status=fe_map(mesh,rbuffery,grid,elts,rbufy,rmask);
      if(strcmp(format,"ascii")==0) {
        strcpy(output,sfile);
        root=strstr(output,".v2r");
        sprintf(root,"%s",".asc");
//        status=ascii_saver2(output, grid, rbufx, rbufy,rmask);
        }
      else if(strcmp(format,"bmg")==0) {
        strcpy(output,sfile);
        root=strstr(output,".v2r");
        sprintf(root,"%s",".bmg");
        status=bmg_saver2(output, 1, 1, 1, grid, rbufx,rbufy,999.99, rmask);
        }
      else if(strcmp(format,"netcdf")==0) {
        strcpy(output,sfile);
        root=strstr(output,".s2r");
        sprintf(root,"%s",".nc");
        status= poc_createfile(output);
        status=poc_sphericalgrid_xy(output,"",grid,&ncgrid);
        poc_standardvariable_xy(&variable,"variable01",rmask,"no_units",1., 0.,"no_name","no_name","no_name",ncgrid);
        status=create_ncvariable(output, &variable);
        status=poc_write_xy(output, grid, variable.id, rbufx);
        variable.destroy();
        poc_standardvariable_xy(&variable,"variable02",rmask,"no_units",1., 0.,"no_name","no_name","no_name",ncgrid);
        status=create_ncvariable(output, &variable);
        status=poc_write_xy(output, grid, variable.id, rbufy);
        variable.destroy();
        }
      free(rbufferx);
      free(rbuffery);
      free(rbufx);
      free(rbufy);
      break;
    case S2C:
      cbuffer=(fcomplex *)malloc(mesh.nvtxs*sizeof(fcomplex));
      cbuf=(fcomplex *)malloc(grid.nx*grid.ny*sizeof(fcomplex));
      status=quoddy_loadc1(sfile, mesh.nvtxs, cbuffer);
      status=fe_map(mesh,cbuffer,grid,elts,cbuf,cmask);
      if(strcmp(format,"ascii")==0) {
        strcpy(output,sfile);
        root=strstr(output,".s2c");
        sprintf(root,"%s",".asc");
        status=ascii_savec1(output, grid, cbuf, cmask);
        }
      else
        {
        strcpy(output,sfile);
        root=strstr(output,".s2c");
        sprintf(root,"%s",".bmg");
        status=bmg_savec1(output, 1, 1, 1, grid, cbuf, 0.,cmask);
        }
      free(cbuffer);
      free(cbuf);
      break;
    case V2C:
      cmask=fcomplex(999.9,999.9);
      cbufferx=(fcomplex *)malloc(mesh.nvtxs*sizeof(fcomplex));
      cbuffery=(fcomplex *)malloc(mesh.nvtxs*sizeof(fcomplex));
      cbufx=(fcomplex *)malloc(grid.nx*grid.ny*sizeof(fcomplex));
      cbufy=(fcomplex *)malloc(grid.nx*grid.ny*sizeof(fcomplex));
      status=quoddy_loadc2(sfile, mesh.nvtxs, cbufferx, cbuffery);
      status=fe_map(mesh,cbufferx,grid,elts,cbufx,cmask);
      status=fe_map(mesh,cbuffery,grid,elts,cbufy,cmask);
      if(strcmp(format,"ascii")==0) {
        strcpy(output,sfile);
        root=strstr(output,".v2c");
        sprintf(root,"%s",".asc");
        status=ascii_savec2(output, grid, cbufx, cbufy, cmask);
        }
      else {
        strcpy(output,sfile);
        root=strstr(output,".v2c");
        sprintf(root,"%s",".bmg");
        status=bmg_savec2(output, 1, 1, 1, grid, cbufx, cbufy, 0.,cmask);
        }
      free(cbufferx);
      free(cbuffery);
      free(cbufx);
      free(cbufy);
      break;

    default:
      break;
    }


end: printf("end of gridit ... \n");
  free(elts);
error:
  __ERR_BASE_LINE__("exiting\n");exit(-1);
}
