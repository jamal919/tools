
/*******************************************************************************

  T-UGO tools, 2006-2009

  Unstructured Ocean Grid initiative

Contributors:

  Florent Lyard      LEGOS/CNRS, Toulouse, France
  Yoann Le Bars      LEGOS, Toulouse, France (PhD)
  Laurent Roblou     LEGOS, Toulouse, France
  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
  Frédéric Dupont    Université de Laval à Québec, Canada

E-mail: florent.lyard@legos.obs-mip.fr

*******************************************************************************/

#define MAIN_SOURCE

#include <config.h>

#include <stdio.h>
#include <string.h>
 
#include "tools-structures.h"
#include "version-macros.def" //for VERSION and REVISION

#include "functions.h"
#include "map.h"
#include "xyz.h"
#include "netcdf-proto.h"
#include "geo.h"
#include "sym-io.h"
#include "polygones.h"
#include "grd.h"
#include "topo.h"
#include "fe-proto.h"
#include "fe.def"

extern   grid_t get_topogrid(grid_t topogrid);
extern   grid_t get_topogrid(grid_t topogrid);

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void print_help(const char *prog_name)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/** \brief prints help of this programme.

\param prog_name name of the program to be printed by this help function : argv[0]
*/
/*----------------------------------------------------------------------------*/
{
/** The code of the body of this function is :
    \code /**/ // COMPILED CODE BELOW !!!
  printf("\n"
    "NAME AND VERSION\n"
    "  %s version " VERSION " " REVISION "\n", prog_name);
  printf("\n"
    "USE\n"
    "  %s -o output [OPTIONS]\n",prog_name);
  printf("\n"
    "DESCRIPTION\n"
/* USE FIND AND REPLACE TO EDIT THE LINE BELOW, SO THAT THE FILE HEAD IS ALSO UPDATED */
    "  create a new MNT.\n"
    "\n"
    "OPTIONS\n"
    "  Grid specification options are:\n"
    "\n"
    );
  print_zone_arg_help();
  printf(
    "\n"
    "  -z [zone]     : zone name for grid specification   \n"
    "\n"
    "  -n [notebook] : notebook file for grid specification \n"
    "\n"
    "  -g : followed by grid file path\n"
    "\n"
    "  Other options are:\n"
    "\n"
    "  --help,-h : show this help and exit\n"
    "  -mask [value] : specify mask value\n"
    "  -b [import]   : merge import MNT in new MNT\n"
    "  -v : folloowed by MNT variable name\n"
    "  -o [output]   : specify output filename\n"
    "  -f : followed by output format: grd, grd-float, nectdf or xyz. Default is netcdf if output file ends with .nc, otherwise is grd\n"
    "  -inv : toggle depth/altitude\n"
    "  -zmin : followed by minimum imported depth\n"
    "  -zmax : followed by maximum imported depth\n"
    "  -debug : \n"
    "  -s : followed by a mask value. It may be a bug that this option has NO effect.\n"
    "  -m : followed by mesh file path\n"
    "  --mapping : followed by frame and resolution, typically \"[-10.0:.0083333333:0.0;42.5:.0083333333:50]\" or equivalently \"[-10.0:30 arcsec:0.0;42.5:30 arcsec:50]\"\n"
    );
//  print_OPENMP_help(prog_name); /** \endcode */
}

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
  double factor;

/*------------------------------------------------------------------------------
  check if basic pattern [ ; ] present */
  pointer = s.find('[');
  if(pointer==string::npos) {
    printf("could not find [ in %s, mapping string decode abort\n",input.c_str());
    return (-1);
    }
  pointer = s.find(';');
  if(pointer==string::npos) return (-1);
  pointer = s.find(']');
  if(pointer==string::npos) {
    printf("could not find ] in %s, mapping string decode abort\n",input.c_str());
    return (-1);
    }
  
  pointer = s.find("[:");
  if(pointer!=string::npos) {
    s.replace(pointer, 2, "[nan:");
    }
  pointer = s.find("::");
  if(pointer!=string::npos) {
    s.replace(pointer, 2, ":nan:");
    }
  pointer = s.find(":;");
  if(pointer!=string::npos) {
    s.replace(pointer, 2, ":nan;");
    }
  pointer = s.find(";:");
  if(pointer!=string::npos) {
    s.replace(pointer, 2, ";nan:");
    }
  pointer = s.find("::");
  if(pointer!=string::npos) {
    s.replace(pointer, 2, ":nan:");
    }
  pointer = s.find(":]");
  if(pointer!=string::npos) {
    s.replace(pointer, 2, ":nan]");
    }
  
 /*------------------------------------------------------------------------------
  remove brackets */
  pointer = s.find('[');
  s[pointer]=' ';
  pointer = s.find(']');
  s[pointer]=' ';
    
 /*------------------------------------------------------------------------------
  remove semi-column, split longitude/latitude setting */
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
    
  factor=1.0;
  
  pointer = s.find("arcsec"); 
  if(pointer!=string::npos) {
    for(int k=0;k<strlen("arcsec");k++) s[pointer+k]=' ';
    factor=1./3600.;
    }
    
  pointer = s.find("mn"); 
  if(pointer!=string::npos) {
    for(int k=0;k<strlen("mn");k++) s[pointer+k]=' ';
    factor=1./60.;
    }
    
  switch(count) {
    case 1:
      nitems=sscanf(s.c_str(),"%lf %lf", &frame.xmin, &frame.xmax);
      if(nitems!=2) return (-1);
      break;
    case 2:
      nitems=sscanf(s.c_str(),"%lf %lf %lf", &frame.xmin, &frame.dx, &frame.xmax);
      if(nitems!=3) return (-1);
      frame.dx*=factor;
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
    
  factor=1.0;
  pointer = s.find("arcsec");
  
  if(pointer!=string::npos) {
    for(int k=0;k<strlen("arcsec");k++) s[pointer+k]=' ';
    factor=1./3600.;
    }
    
  switch(count) {
    case 1:
      nitems=sscanf(s.c_str(),"%lf %lf", &frame.ymin, &frame.ymax);
      if(nitems!=2) return (-1);
      break;
    case 2:
      nitems=sscanf(s.c_str(),"%lf %lf %lf", &frame.ymin, &frame.dy, &frame.ymax);
      if(nitems!=3) return (-1);
      frame.dy*=factor;
     break;
    default:
     return (-1);
    }
     
//   status=map_set2Dgrid(&frame, frame.xmin, frame.ymin, frame.xmax, frame.ymax, frame.dx, frame.dy);
  status=0;
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  float  spec=-9999;
  int status;
  int i,j,n,count;
  size_t m;
  const char *keyword,*notebook=NULL;
  char *poly=NULL, *depthfile=NULL,*output=NULL,*path=NULL, *gridfile=NULL, *zone=NULL;
  char *meshfile=NULL,*format=NULL,*varname=NULL;
  grid_t topogrid;
  grid_t cartesian_topogrid;
  string cmd;
  cdfvar_t variable;
  pocgrd_t z_ncgrid, u_ncgrid, v_ncgrid;
  geo_t projection;
  frame_t frame,prescribed(1.e+10,1.e+10,1.e+10,1.e+10);
  pocgrd_t ncgrid;
  double dx=1.e+10,dy=1.e+10,delta;
  mesh_t mesh;
  int *elts=0;
  string MappingString;
  
/* *-----------------------------------------------------------------------------------------------------
                     15"       20"       30"       1'     2'     5'     6'     10'     15'     30'    */
  double preset[10]={1./60./4.,1./60./3.,1./60./2.,1./60.,2./60.,5./60.,6./60.,10./60.,15./60.,30./60.};
  short *topo,smask=256*127+255;
  float *ftopo,ftopomask,scale=1.0;

  float zmin=0,zmax=0;

  bool debug=false;

  fct_echo( argc, argv);

  ftopomask=spec;

  n=1;
  while (n < argc) {
    keyword=argv[n];
    if(strcmp(keyword,"--help")==0 or
       strcmp(keyword,"-h")==0 ) {
      print_help(argv[0]);
      exit(0);
      }
    if(read_zone_arg(keyword,argv[n+1],&prescribed,&dx,&dy)) {
      n++;
      n++;
      continue;
      }
    if(strcmp(keyword,"-inv")==0) {
      scale=-1;
      n++;
      continue;
      }
    if(strcmp(keyword,"-zmin")==0) {
      sscanf(argv[n+1],"%f",&zmin);
      n++;
      n++;
      continue;
      }
    if(strcmp(keyword,"-zmax")==0) {
      sscanf(argv[n+1],"%f",&zmax);
      n++;
      n++;
      continue;
      }
    if(strcmp(keyword,"-mask")==0) {
      sscanf(argv[n+1],"%f",&ftopomask);
      n++;
      n++;
      continue;
      }
    if(strcmp(keyword,"-debug")==0) {
      debug=true;
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
        case 'g' :
          gridfile= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'z' :
          zone= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'f' :
          format= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'n' :
          notebook= argv[n+1];
          n++;
          n++;
          break;

        case 'b' :
          depthfile= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 's' :
          sscanf(argv[n+1],"%f",&spec);
          n++;
          n++;
          break;

        case 'o' :
          output= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'm' :
          meshfile= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'v' :
          varname= strdup(argv[n+1]);
          n++;
          n++;
          break;

        default:
          STDOUT_BASE_LINE("unknown option %s\n",keyword);
          exit(-1);
        }
        break;

      default:
        STDOUT_BASE_LINE("unknown option %s\n",keyword);
        exit(-1);
        break;
      }
    }

  if(path ==NULL) path=strdup(".");

/* *-----------------------------------------------------------------------------
  Regular grid prescription: priority sequence*/
  if(notebook!=0) {
/* *-----------------------------------------------------------------------------
    Read notebook data */
    status=load_notebook(notebook, &cartesian_topogrid, &topogrid, &projection);
    printf("%s (notebook file) processed\n",notebook);
    }
  else if(gridfile != NULL) {
/* *-----------------------------------------------------------------------------
    Read grid file */
    short *sbuffer,smask;
    status=topo_loadfield(gridfile, &topogrid, &sbuffer, &smask, debug);
    map_printgrid(topogrid);
    delete[] sbuffer;
    }
  else if(zone != NULL) {
/* *-----------------------------------------------------------------------------
    Apply zone definition by name*/
    topogrid=get_zonegrid(zone);
    status=map_completegridaxis(&topogrid,0);
    }
/*----------------------------------------------------------------------------
  Read FE mesh */
  else if(meshfile != NULL) {
    status=fe_readmesh(meshfile,MESH_FILE_FORMAT_TRIGRID,&mesh);
    if(status !=0) TRAP_ERR_EXIT(status,"fe_readmesh(\"%s\") error\n",meshfile);
    status=fe_list(&mesh);
    if(status !=0) TRAP_ERR_EXIT(status,"fe_list((\"%s\")) error\n",meshfile);
    status=fe_minmax(mesh, frame);
    topogrid.xmin=frame.xmin;
    topogrid.ymin=frame.ymin;
    topogrid.xmax=frame.xmax;
    topogrid.ymax=frame.ymax;
// /* *-----------------------------------------------------------------------------
//     default resolution (preset)*/
//     count=0;
//     topogrid.dx=preset[0];
//     topogrid.dy=preset[0];
//     topogrid.nx=floor((topogrid.xmax-topogrid.xmin)/topogrid.dx+0.5)+1;
//     topogrid.ny=floor((topogrid.ymax-topogrid.ymin)/topogrid.dy+0.5)+1;
//     while (sqrt(topogrid.nx*topogrid.ny)>1200) {
//       count++;
//       topogrid.dy=preset[count];
//       topogrid.ny=floor((topogrid.ymax-topogrid.ymin)/topogrid.dy+0.5)+1;
//       topogrid.dx=preset[count];
//       topogrid.nx=floor((topogrid.xmax-topogrid.xmin)/topogrid.dx+0.5)+1;
//       }
// 
//     delta=max(topogrid.dx,topogrid.dy);
//     topogrid.dx=delta;
//     topogrid.dy=delta;

/* *-----------------------------------------------------------------------------
    apply prescribed resolution*/
    if(dx!=1.e+10) topogrid.dx=dx;
    if(dy!=1.e+10) topogrid.dy=dy;
    else if(dx!=1.e+10) topogrid.dy=dx;

    topogrid.xmin=floor(topogrid.xmin/topogrid.dx)*topogrid.dx;
    topogrid.xmax=floor(topogrid.xmax/topogrid.dx+1)*topogrid.dx;
    topogrid.ymin=floor(topogrid.ymin/topogrid.dy)*topogrid.dy;
    topogrid.ymax=floor(topogrid.ymax/topogrid.dy+1)*topogrid.dy;

    topogrid.nx=floor((topogrid.xmax-topogrid.xmin)/topogrid.dx+0.5)+1;
    topogrid.ny=floor((topogrid.ymax-topogrid.ymin)/topogrid.dy+0.5)+1;
    topogrid.modeH=0;
    status=map_completegridaxis(&topogrid, 1);
    }
  else if(MappingString!="") {
    status=map_DecodeGrid(MappingString, topogrid);
    if(status!=0) TRAP_ERR_EXIT(-1,"mapping string not understood (may be you forgot \"): %s\n",MappingString.c_str());
    status=map_set2Dgrid(&topogrid, topogrid.xmin, topogrid.ymin, topogrid.xmax, topogrid.ymax, topogrid.dx, topogrid.dy);
    status=map_completegridaxis(&topogrid, 1);
    }
  else {
/* *-----------------------------------------------------------------------------
    Apply frame definition*/

/* *-----------------------------------------------------------------------------
    default frame (mesh limits)*/
    topogrid.xmin=frame.xmin;
    topogrid.ymin=frame.ymin;
    topogrid.xmax=frame.xmax;
    topogrid.ymax=frame.ymax;

/* *-----------------------------------------------------------------------------
    apply prescribed limits*/
    if(prescribed.xmin!=1.e+10) topogrid.xmin=prescribed.xmin;
    if(prescribed.xmax!=1.e+10) topogrid.xmax=prescribed.xmax;
    if(prescribed.ymin!=1.e+10) topogrid.ymin=prescribed.ymin;
    if(prescribed.ymax!=1.e+10) topogrid.ymax=prescribed.ymax;

/* *-----------------------------------------------------------------------------
    default resolution (preset)*/
    count=0;
    topogrid.dx=preset[0];
    topogrid.dy=preset[0];
    topogrid.nx=floor((topogrid.xmax-topogrid.xmin)/topogrid.dx+0.5)+1;
    topogrid.ny=floor((topogrid.ymax-topogrid.ymin)/topogrid.dy+0.5)+1;
    while (sqrt(topogrid.nx*topogrid.ny)>1200) {
      count++;
      topogrid.dy=preset[count];
      topogrid.ny=floor((topogrid.ymax-topogrid.ymin)/topogrid.dy+0.5)+1;
      topogrid.dx=preset[count];
      topogrid.nx=floor((topogrid.xmax-topogrid.xmin)/topogrid.dx+0.5)+1;
      }

    delta=max(topogrid.dx,topogrid.dy);
    topogrid.dx=delta;
    topogrid.dy=delta;

/* *-----------------------------------------------------------------------------
    apply prescribed resolution*/
    if(dx!=1.e+10) topogrid.dx=dx;
    if(dy!=1.e+10) topogrid.dy=dy;
    else if(dx!=1.e+10) topogrid.dy=dx;

    topogrid.xmin=floor(topogrid.xmin/topogrid.dx)*topogrid.dx;
    topogrid.xmax=floor(topogrid.xmax/topogrid.dx+1)*topogrid.dx;
    topogrid.ymin=floor(topogrid.ymin/topogrid.dy)*topogrid.dy;
    topogrid.ymax=floor(topogrid.ymax/topogrid.dy+1)*topogrid.dy;
    
    topogrid.modeH=0;

    topogrid.nx=floor((topogrid.xmax-topogrid.xmin)/topogrid.dx+0.5)+1;
    topogrid.ny=floor((topogrid.ymax-topogrid.ymin)/topogrid.dy+0.5)+1;
//     status=map_completegridaxis_2(&topogrid);
    status=map_completegridaxis(&topogrid,1);
    }
  
  if (status!=0) TRAP_ERR_EXIT(-1,"map_completegridaxis() error\n");

/* *----------------------------------------------------------------------
  initialize bathymetry data*/
  ftopo=aset(topogrid.Hsize(),ftopomask);
  topogrid.nz=1;
   
  printf("#################################################################\n");
  printf("create depth grid : \n");
  map_printgrid(topogrid);

  if(depthfile!=0) {
    printf("#################################################################\n");
    printf("import depth from : %s\n",depthfile);
    status= topo_import(depthfile, varname, topogrid, ftopo, ftopomask,  zmin,  zmax, poly, 0, debug);
    if (status!=0) TRAP_ERR_EXIT(-1,"topo_import() error %d\n",status);
    }
  
/*----------------------------------------------------------------------------
  Read FE mesh */
  if(meshfile != NULL) {
    printf("#################################################################\n");
    printf("import depth from : %s\n",meshfile);
    elts=fe_scan_elements(mesh,topogrid,0);
    if(elts ==0) TRAP_ERR_EXIT(-1,"fe_scan_elements() error\n");
    float *depths=new float[mesh.nvtxs];
    for(n=0;n<mesh.nvtxs;n++) depths[n]=-mesh.vertices[n].h;
    status=fe_map(mesh,depths,topogrid,elts,ftopo,ftopomask);
    if(status!=0) TRAP_ERR_EXIT(-1,"fe_map() error\n");
    delete depths;
    }
  
//     for (j=0;j<topogrid.ny;j++) {
//       for (i=0;i<topogrid.nx;i++) {
//         m=j*topogrid.nx+i;
//         ftopo[m]=ftopomask;
//         }
//       }

  if(scale!=1.) {
    for (j=0;j<topogrid.ny;j++) {
      for (i=0;i<topogrid.nx;i++) {
        m=j*topogrid.nx+i;
        if(ftopo[m]!=ftopomask) {
          ftopo[m]*=scale;
          }
        }
      }
    }

  if(format==NULL){
    if(output==0) {
      format=strdup("grd");
      printf("format not specified, use default grd\n");
      } 
    else if(strrncasecmp(output,".grd")==0) {
      format=strdup("grd");
      printf("format not specified, follows extension, use grd\n");
      }
    else if(strrncasecmp(output,".nc")==0) {
      format=strdup("netcdf");
      printf("format not specified, follows extension, use netcdf\n");
      }
    else {
      format=strdup("grd");
      printf("format not specified, use default grd\n");
      }
    }
    
  if (output==0) {
    asprintf(&output,"%s.grd","topo-create");
    }


/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  save data, to be replaced with topo_save. Example:
  
  status=topo_save(rootname, output, "netcdf", TOPOgrid, TOPOdata, TOPOmask, debug);
  
  topo_save will hande mirroring...

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if(strcmp(format,"grd")==0) {
    topo=new short[topogrid.Hsize()];
    smask=256*127+255;
    for (j=0;j<topogrid.ny;j++) {
      for (i=0;i<topogrid.nx;i++) {
        m=(size_t)j*(size_t)topogrid.nx+i;
        size_t n=(size_t) (topogrid.ny-j-1)*(size_t) topogrid.nx+i;
        if(ftopo[m]!=ftopomask) {
          topo[n]=(short) floor(ftopo[m]+0.5);
          }
        else {
          topo[n]=smask;
          }
        }
      }
    printf("#################################################################\n");
    printf("create depth file (grd short format) : %s\n",output);
    status=grd_save((const char*) output,topogrid,topogrid.nx, topo,smask);
    delete[] topo;
    }
  else if(strcmp(format,"grd-float")==0) {
    printf("#################################################################\n");
    printf("create depth file (grd float format) : %s\n",output);
    status=grd_mirror_r( topogrid, topogrid.nx, ftopo, ftopomask);
    status=grd_save((const char*) output,topogrid,topogrid.nx, ftopo, ftopomask);
    }
  else if(strcmp(format,"netcdf")==0) {
    printf("#################################################################\n");
    printf("create depth file (netcdf format) : %s\n",output);
    status= poc_createfile(output);
//     status=poc_sphericalgrid_xyzt(output,"","",topogrid,&ncgrid);
    status=poc_sphericalgrid_xy(output,"",topogrid,&ncgrid);
    poc_standardvariable_xy(&variable,"bathymetry",ftopomask,"m",1., 0.,"bathymetry","bathymetry","bathymetry",ncgrid);
    status=create_ncvariable(output, &variable);
    status=poc_write_xy(output,  topogrid, variable.id,ftopo);
    variable.destroy();
    }
  else if(strcmp(format,"xyz")==0) {
    if (output==0) {
      asprintf(&output,"%s.nc","topo-create");
      }
    printf("#################################################################\n");
    printf("create depth file (xyz format) : %s\n",output);
    status=xyz_save((const char *) output,topogrid, ftopo,ftopomask);
    }

  delete[] ftopo;

  TRAP_ERR_EXIT(0,"exiting\n");
}
