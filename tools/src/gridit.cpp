
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

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_DecodeGrid(const string input, grid_t & frame)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
//  [tmin:tmax:dt; pmin:pmax:dp]
  int status;
  string s=input;
  string substring,longitude, latitude;
  size_t pointer;
  int count,nitems;
  
  pointer = s.find('[');
  if(pointer==string::npos) return (-1);
  pointer = s.find(';');
  if(pointer==string::npos) return (-1);
  pointer = s.find(']');
  if(pointer==string::npos) return (-1);
  
  pointer = s.find('[');
  s[pointer]=' ';
  pointer = s.find(']');
  s[pointer]=' ';
  
  pointer = s.find(';');
  s[pointer]=' ';
  longitude =s.substr(0,pointer);
  latitude  =s.substr(pointer+1);
  
  s=longitude;
  pointer = s.find(':');
  if(pointer==string::npos) return (-1);
  
  count=0;
  while(pointer!=string::npos) {
    s[pointer]=' ';
    pointer = s.find(':');
    count++;
    }
  switch(count) {
    case 1:
      nitems=sscanf(s.c_str(),"%lf %lf", &frame.xmin, &frame.xmax);
      if(nitems!=2) return (-1);
      break;
    case 2:
      nitems=sscanf(s.c_str(),"%lf %lf %lf", &frame.xmin, &frame.dx, &frame.xmax);
      if(nitems!=3) return (-1);
      break;
    default:
     return (-1);
    }
    
  s=latitude;
  pointer = s.find(':');
  if(pointer==string::npos) return (-1);
  
  count=0;
  while(pointer!=string::npos) {
    s[pointer]=' ';
    pointer = s.find(':');
    count++;
    }
  switch(count) {
    case 1:
      nitems=sscanf(s.c_str(),"%lf %lf", &frame.ymin, &frame.ymax);
      if(nitems!=2) return (-1);
      break;
    case 2:
      nitems=sscanf(s.c_str(),"%lf %lf %lf", &frame.ymin, &frame.dy, &frame.ymax);
      if(nitems!=3) return (-1);
      break;
    default:
     return (-1);
    }
     
  status=map_set2Dgrid(&frame, frame.xmin, frame.ymin, frame.xmax, frame.ymax, frame.dx, frame.dy);
     
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  float *rbuf=NULL;
  float *rbufx=NULL,*rbufy=NULL;
  float  dummy,rmask,*rbuffer=NULL;
  float  *rbufferx=NULL,*rbuffery=NULL;
  float  spec,scale_factor=1.0;
  int status,*elts=NULL,fmt,fefmt;
  int i,j,k,l,L,m,n;
  int nitems;
  FILE *file=NULL;
  FILE *out=NULL;
  char *meshfile=NULL,*sfile=NULL,*cfile=NULL,*format=NULL,*cdl=0;
  char *keyword=NULL,*zone=NULL,*gridfile=NULL,*notebook=NULL;
  char file1[256],file2[256];
  char *root=NULL,tmp[256],output[1024]="";
  grid_t cgrid,grid;
  fcomplex *cbuffer=NULL,*cbuf=NULL,cmask;
  fcomplex *cbufferx=NULL,*cbuffery=NULL,*cbufx=NULL,*cbufy=NULL;
  mesh_t mesh;
  geo_t projection;
  cdfvar_t variable;
  pocgrd_t ncgrid;
  string cmd;
  int persistence=0;
  string MappingString="";

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
    if(strcmp(keyword,"--persistence")==0) {
      sscanf(argv[n+1],"%d",&persistence);
      n++;
      n++;
      continue;
      }
    if(strncmp("--mapping",keyword)==0){
      MappingString= argv[n+1];
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

 if(sfile == NULL) {
   printf("no FE file specified; abort...\n");
   goto error;
   }
   
  rmask=spec;
  cmask=fcomplex(spec,spec);

  if(format==NULL) format= strdup("netcdf");

  if(strstr(sfile,".s2r") != NULL)    fefmt=S2R;
  if(strstr(sfile,".s2c") != NULL)    fefmt=S2C;
  if(strstr(sfile,".v2r") != NULL)    fefmt=V2R;
  if(strstr(sfile,".v2c") != NULL)    fefmt=V2C;

  if(zone != NULL) {
    bool identified;
    grid=get_zonegrid(zone, &identified);
    if(not identified) goto error;
    status=map_completegridaxis_2(&grid);
    if(status !=0) goto error;
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
    
  if(MappingString!="") {
    status=map_DecodeGrid(MappingString, grid);
    if(status!=0) TRAP_ERR_EXIT(-1,"mapping string not understood (may be you forgot \"): %s\n",MappingString.c_str());
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

  elts=fe_scan_elements(mesh,grid,0);
  if(elts==NULL) goto error;

  switch (fefmt) {
    case S2R:
      exitIfNull(
        rbuffer=(float *)malloc(mesh.nvtxs*sizeof(float))
        );
      exitIfNull(
        rbuf=(float *)malloc(grid.nx*grid.ny*sizeof(float))
        );
      status=quoddy_loadr1(sfile, mesh.nvtxs, rbuffer);
      if(status!=0) goto error;
      status=fe_map(mesh,rbuffer,grid,elts,rbuf,rmask);
      if(status!=0) goto error;
      for(n=0;n<grid.nx*grid.ny;n++) {
        if(rbuf[n]!=rmask) {
          rbuf[n]*=scale_factor;
          }
        }
      if(persistence!=0) {
        for(int k=0;k<persistence;k++) status=map_persistence(grid, rbuf, rmask, 0);
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

        if (strlen(output)==0) {
          strcpy(output,sfile);
          }
        else {
          sprintf(output,"%s-%s",output,sfile);
          }
        root=strstr(output,".s2r");
        sprintf(root,"%s",".grd");
        status=grd_mirror_r( grid, grid.nx, rbuf, rmask);
        status=grd_save(output, grid, grid.nx, rbuf, rmask);
        }
      else if(strcmp(format,"netcdf")==0) {
        if (strlen(output)==0) {
          strcpy(output,sfile);
          }
        else {
          sprintf(output,"%s-%s",output,sfile);
          }
        root=strstr(output,".s2r");
        sprintf(root,"%s",".nc");
        status= poc_createfile(output);
        status=poc_sphericalgrid_xy(output,"",grid,&ncgrid);
        poc_standardvariable_xy(&variable,"variable01",rmask,"no_units",1., 0.,"no_name","no_name","no_name",ncgrid);
        status=create_ncvariable(output, &variable);
        status=poc_write_xy(output, grid, variable.id, rbuf);
        variable.destroy();
        for(n=0;n<grid.nx*grid.ny;n++) rbuf[n]=elts[n];
        poc_standardvariable_xy(&variable,"element-index",-1.,"no_units",1., 0.,"element-index","element-index","element-index",ncgrid);
        status=create_ncvariable(output, &variable);
        status=poc_write_xy(output, grid, variable.id, rbuf);
        variable.destroy();
        }
      break;
    case V2R:
      exitIfNull(
        rbufferx=(float *)malloc(mesh.nvtxs*sizeof(float))
        );
      exitIfNull(
        rbuffery=(float *)malloc(mesh.nvtxs*sizeof(float))
        );
      exitIfNull(
        rbufx=(float *)malloc(grid.nx*grid.ny*sizeof(float))
        );
      exitIfNull(
        rbufy=(float *)malloc(grid.nx*grid.ny*sizeof(float))
        );
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
        if (strlen(output)==0) {
          strcpy(output,sfile);
          }
        else {
          sprintf(output,"%s-%s",output,sfile);
          }
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
      exitIfNull(
        cbuffer=(fcomplex *)malloc(mesh.nvtxs*sizeof(fcomplex))
        );
      exitIfNull(
        cbuf=(fcomplex *)malloc(grid.nx*grid.ny*sizeof(fcomplex))
        );
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
      exitIfNull(
        cbufferx=(fcomplex *)malloc(mesh.nvtxs*sizeof(fcomplex))
        );
      exitIfNull(
        cbuffery=(fcomplex *)malloc(mesh.nvtxs*sizeof(fcomplex))
        );
      exitIfNull(
        cbufx=(fcomplex *)malloc(grid.nx*grid.ny*sizeof(fcomplex))
        );
      exitIfNull(
        cbufy=(fcomplex *)malloc(grid.nx*grid.ny*sizeof(fcomplex))
        );
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
  delete [] elts;
  exit(0);

error:
  __ERR_BASE_LINE__("error detected, exiting\n");
  exit(-1);
}
