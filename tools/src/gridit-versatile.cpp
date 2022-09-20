
/*******************************************************************************

  T-UGO tools, 2006-2013

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

<!-- USE FIND AND REPLACE TO EDIT THIS LINE, SO THAT print_help IS ALSO UPDATED -->
\brief Converts unstructured to structured tidal atlases

<!-- A LINK TO main() or print_help() WILL NOT LINK TO THE RIGHT SOURCE ! -->
See the main function for how this works
and the print_help function for how to use this.
*/
/*----------------------------------------------------------------------------*/

#define MAIN_SOURCE

#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "version-macros.def" //for VERSION and REVISION

#include "tools-structures.h"

#include "functions.h"
#include "map.h"
#include "xyz.h"
#include "fe.h"
#include "poc-netcdf-data.hpp"
#include "archive.h"
#include "geo.h"
#include "sym-io.h"
#include "tides.h"
#include "cefmo.h"

extern   grid_t get_zgrid(grid_t zgrid);
extern   grid_t get_zgrid(grid_t zgrid);

string cmd;

class foptions_t {
private:
public:
  bool create, append, overwrite;
  foptions_t() {
    create=append=overwrite=false;
    }
};

// class decode_t {
// private:
// public:
//   bool create, append, overwrite;
//   foptions_t() {
//     create=append=overwrite=false;
//     }
// };

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

/*------------------------------------------------------------------------------
  check if basic pattern [ ; ] present */
  pointer = s.find('[');
  if(pointer==string::npos) return (-1);
  pointer = s.find(';');
  if(pointer==string::npos) return (-1);
  pointer = s.find(']');
  if(pointer==string::npos) return (-1);
  
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
  
//   status=map_set2Dgrid(&frame, frame.xmin, frame.ymin, frame.xmax, frame.ymax, frame.dx, frame.dy);
  status=0;
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int map_CheckMasked(const grid_t & grid, T* buffer, T spec, signed char *mask, size_t & missing, size_t & extra)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  bool debug=false;
  
  missing=extra=0;
  
/*------------------------------------------------------------------------------
  no pre-defined mask, return 0 missing */
  if(mask==0) return(0);
  
/*------------------------------------------------------------------------------
  check values against pre-defined mask */
  for(size_t m=0;m<grid.Hsize();m++) {
    if(buffer[m]==spec and mask[m]==1) {
      if(debug) printf("missing value at %d\n",m);
      missing++;
      }
    else if(buffer[m]!=spec and mask[m]==0) {
      if(debug) printf("extra value at %d\n",m);
      buffer[m]=spec;
      extra++;
      }
    }
 
  printf("#missing : %d\n",missing);
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int gridit_NETCDF_vector1D(mesh_t mesh, vector<string*> varnames, vector<string*> outnames, const grid_t & grid, signed char *mask, char **wave, int nwave, int iteration,
                  vector<string> input, vector<string> output, const char *format, poc_var_t grid_var, foptions_t foptions, vector<string>  discretisation, int persistence)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  extern toponyms_t paire_map(void);
  struct timeval before;
  
  int k, status,id;
  char amplitude[256],phase[256];
  tide2D_t state;
  atlas2D_t atlas;
//   mesh_t mesh;
  discretisation_t *descriptor;
  int *elts=0;
  complex<float> *UGbufC=0, *SGbufC=0, cmask=1.e+10;
  float          *UGbufR=0, *SGbufR=0, rmask=1.e+10;
  bool is_complex=false;
  size_t missing, missing_bkp=-1,extra;

  poc_var_t input_var[2], av, gv;
  
  toponyms_t Pmap=paire_map();
  
  int nvars=0;


/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  read mesh geometry
  
*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if(mesh.nvtxs==0) {
    printf("#################################################################\n");
    printf("load mesh geometry from %s\n",input[0].c_str());
    status=fe_readgeometry(input[0].c_str(), &mesh);
    if(status) NC_TRAP_ERROR(return,status,1,"fe_readgeometry(\"%s\",) error",input[0].c_str());
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  detect elements corresponding to grid nodes
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  printf("scanning elements (may take ages): ");
  fflush(stdout);
  gettimeofday(&before);
  
  elts=fe_scan_elements(mesh,grid,0,1);
  if(elts==NULL) TRAP_ERR_RETURN(-1,1,"fe_scan_elements() error\n");
    
  printf(" Took %gs.\n",__func__,difftime(before));


/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  identify discretisation
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
//   id=discretisation_from_name(discretisation, 0);
//   if(id<0) TRAP_ERR_RETURN(status,1,"discretisation_id() error\n");
    
  for(int v=0;v<varnames.size();v++) {
    id=discretisation_from_name(discretisation[v].c_str(), 0);
    if(id<0) TRAP_ERR_RETURN(status,1,"discretisation_id() error\n");
    printf("#################################################################\n");
    printf("treating %s, %s discretisation\n",varnames[v][0].c_str(), discretisation_name(id));
    status=discretisation_init(&mesh,id,0);
    if(status) NC_TRAP_ERROR(return,status,1,"discretisation_init(,%d,0) error",id);
    
    descriptor=get_descriptor_address(mesh,id);
    if(descriptor->nnodes==0) TRAP_ERR_RETURN(-1,1,"get_descriptor() error\n");
      
    size_t nxy=grid.Hsize();
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
    setup variable names and variable type (real or complex)
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

    if(varnames[v][1]=="") {
      sprintf(amplitude,"%s",varnames[v][0].c_str());
      SGbufR=new float[nxy];
      UGbufR=new float[descriptor->nnodes];
      is_complex=false;
      nvars=1;
      }
    else {
      sprintf(amplitude,"%s",varnames[v][0].c_str());
      sprintf(phase,    "%s",varnames[v][1].c_str());
      SGbufC=new complex<float>[nxy];
      UGbufC=new complex<float>[descriptor->nnodes];
      is_complex=true;
      nvars=2;
      }
    
    for (k=0;k<nwave;k++) {
//       status=tide_decode_atlasname(directory, convention, wave[k], 0, &input);
      printf("treating %s from %s\n",varnames[v][0].c_str(), input[k].c_str());
      
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
      load variable
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

      if(is_complex) {
        if(strcmp(format, "ascii")==0) {
          status=quoddy_loadc1(input[k].c_str(),descriptor->nnodes,UGbufC);
          }
        else {
          status=poc_inq_var(input[k], amplitude, &input_var[0]);
          if(status!=0) TRAP_ERR_EXIT(status,"read error\n");
          status=poc_inq_var(input[k], phase,     &input_var[1]);
          status=poc_get_UG3D(input[k].c_str(), iteration, amplitude, phase, UGbufC);
          }
        }
      else {
        float inputMask;
        if(strcmp(format, "ascii")==0) {
          status=quoddy_loadr1(input[k].c_str(),descriptor->nnodes,UGbufR);
          }
        else {
          status=poc_inq_var(input[k], amplitude, &input_var[0]);
          status=poc_decode_mask(input_var[0], 0, 0, &inputMask);
          status=poc_get_UG3D(input[k].c_str(), iteration, amplitude, UGbufR);
          }
        }
      
      if(status!=0) {
        printf("elevation not done\n");
        continue;
        }
      
      gettimeofday(&before);
      
      if(is_complex) {
        status=fe_map(mesh, UGbufC, id, grid, elts, SGbufC, cmask);
        if(persistence!=0) for(int loop=0; loop<persistence; loop++) status=map_persistence(grid, SGbufC, cmask, 0);
        status=map_CheckMasked(grid, SGbufC, cmask, mask, missing, extra);
        
        while(missing!=0) {
          missing_bkp=missing;
          for(int loop=0; loop<1; loop++) status=map_persistence(grid, SGbufC, cmask, 0);
          status=map_CheckMasked(grid, SGbufC, cmask, mask, missing, extra);
          if(missing==missing_bkp) break;
          }
        }
      else {
        status=fe_map(mesh, UGbufR, id, grid, elts, SGbufR, rmask);
        if(persistence!=0) for(int loop=0; loop < persistence; loop++) status=map_persistence(grid, SGbufR, rmask, 0);
        status=map_CheckMasked(grid, SGbufR, rmask, mask, missing, extra);

        while(missing!=0) {
          missing_bkp=missing;
          for(int loop=0; loop<1; loop++) status=map_persistence(grid, SGbufR, rmask, 0);
          status=map_CheckMasked(grid, SGbufR, rmask, mask, missing, extra);
          if(missing==missing_bkp) break;
          }
        }
     
      STDERR_BASE_LINE_FUNC("fe_map() took %gs\n",difftime(before));

      av=grid_var;
      gv=grid_var;
      if(is_complex){
        if(outnames.size()<=v) outnames=varnames;
        if(outnames[v][0]=="") outnames[v][0]=varnames[v][0];
        if(outnames[v][1]=="") outnames[v][1]=varnames[v][1];
        av.init(outnames[v][0], NC_FLOAT, comodo_standard_name("",wave[k]), "m",rmask);
        gv.init(outnames[v][1], NC_FLOAT, comodo_standard_name("",wave[k]), "degrees",rmask);
        poc_att_t *att;
        att=input_var[0].attributes.findP("units");
        if(att!=0){
          av << *att;
          }
        att=input_var[0].attributes.findP("coordinates");
        if(att!=0){
          av << *att;
          }
        att=input_var[0].attributes.findP("standard_name");
        if(att!=0){
          av << *att;
          }
        att=input_var[1].attributes.findP("units");
        if(att!=0){
          gv << *att;
          }
        att=input_var[1].attributes.findP("coordinates");
        if(att!=0){
          av << *att;
          }
        att=input_var[1].attributes.findP("standard_name");
        if(att!=0){
          gv << *att;
          }
        poc_put_cvara(output[k],av,gv,0,SGbufC);
        }
      else{
        if(outnames.size()<=v) outnames=varnames;
        if(outnames[v][0]=="") outnames[v][0]=varnames[v][0];
        av.init(outnames[v][0],NC_FLOAT,"","",rmask);
        poc_att_t *att;
        att=input_var[0].attributes.findP("units");
        if(att!=0){
          av << *att;
          }
        att=input_var[0].attributes.findP("coordinates");
        if(att!=0){
          av << *att;
          }
        att=input_var[0].attributes.findP("standard_name");
        if(att!=0){
          av << *att;
          }
        av << poc_att_t("production","many attributes not coded yet at " __LINE_FILE_PACKAGE_REVISION);
        poc_put_vara(output[k],av,0,SGbufR);
        }
      }
    }
 
  if(SGbufR!=0) delete[] SGbufR;
  if(UGbufR!=0) delete[] UGbufR;
  if(SGbufC!=0) delete[] SGbufC;
  if(UGbufC!=0) delete[] UGbufC;
  
  descriptor->destroy();
    
  Pmap.clear();
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int process_GridOptions(const string gridfile, const char *vlon, const char *vmask, grid_t & grid, signed char* & mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int verbose=0, status;
  poc_data_t<int8_t> t_mask,u_mask,v_mask;
  
  status=poc_get_grid(gridfile, vlon, &grid, verbose, -1);
  if(status !=0) TRAP_ERR_EXIT(status,"parse_GridOptions error\n");;
  
  if(vmask!=0) {
    status=t_mask.init(gridfile,vmask,1);
    status=t_mask.read_data(gridfile,-1,1);
    if(t_mask.nlimited==3) {
      printf("actually 3D mask, extract 2D\n");
      int ndim=t_mask.info.dimensions.size();
      int k;
      k=0;
      size_t lower=occurence((int8_t) 0, &t_mask.data[k*grid.Hsize()], grid.Hsize());
      k=t_mask.info.dimensions[ndim-3].len-1;
      size_t upper=occurence((int8_t) 0, &t_mask.data[k*grid.Hsize()], grid.Hsize());
      mask=new signed char[grid.Hsize()];
      if(lower<upper) {
        k=0;
        printf("upper layer is for k=%d (%5.1f\%)\n",k,100.*lower/grid.Hsize());
        for(size_t m=0; m<grid.Hsize(); m++) mask[m]=t_mask.data[k*grid.Hsize()+m];
        }
      else {
        k=t_mask.info.dimensions[ndim-3].len-1;
        printf("upper layer is for k=%d (%5.1f\%)\n",k,100.*upper/grid.Hsize());
        for(size_t m=0; m<grid.Hsize(); m++) mask[m]=t_mask.data[k*grid.Hsize()+m];
        }
      }
    }

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int duplicate_grid(const string gridfile, const char *vlon, const char *vlat, const string output, poc_var_t & grid_var)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  poc_data_t<double> x,y;
  
  status=x.init(gridfile,vlon,1);
  status=x.read_data(gridfile,-1,1);
  status=x.write_data(output,-2,1);

  status=y.init(gridfile,vlat,1);
  status=y.read_data(gridfile,-1,1);
  status=y.write_data(output,-2,1);
  
  grid_var=x.info;

  return(0);
}


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
    " %s version " VERSION " " REVISION "\n", prog_name);
  printf("\n"
    "USE\n"
    "  %s input [OPTIONS] wave1 [wave2 ...]\n",prog_name);
  printf("\n"
    "DESCRIPTION\n"
/* USE FIND AND REPLACE TO EDIT THE LINE BELOW, SO THAT THE FILE HEAD IS ALSO UPDATED */
    "  Converts unstructured to structured tidal atlases.\n"
    "\n"
    "OPTIONS :\n"
    "  -h,--help  show this help and exit\n"
    );
  print_zone_arg_help();
  printf(
    "  -discretisation  followed by discretisation pair of netcdf input atlas\n"
    "  -c  followed by netcdf input atlas file name convention. See below.\n"
    "  -p  followed by directory containing the atlas files. Default: .\n"
    "  -v  followed by an empty-terminated list of up to 6 variable names for the amplitudes and phases of e, u and v. Defaults computed from discretisation as TUGO would.\n"
    "  -a  followed by frame index. Default: -1 (last).\n"
    "  -g  followed by grid file AND by an empty-terminated list of up to 3 gridded variable names for e, u and v\n"
    "  -n  followed by notebook file\n"
    "  -m  followed by mesh file\n"
    "  -z  followed by zone name\n"
    "  -o  followed by some output root name. Default: generic.\n"
    "  -f  followed by input format. Default: ascii. If anything else: netcdf, in which case discretisation will be compulsory.\n"
    "\n"
    "EXAMPLES :\n"
    "  %s -z global -f nc -c WAVE.spectral.nc -o FES2012 M2 -discretisation DNP1xLGP2\n"
    "  %s -dx 1m -dy 1m -f nc -c WAVE.spectral.energy.nc -v WaveDragRow '' -o energy-RG M2 -discretisation LGP1xLGP0\n"
    "\n"
    "SEE ALSO: tides-converter\n"
    , prog_name, prog_name);
  
  print_tide_decode_atlasname_help();
  /** \endcode */
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int k,n;

  char *keyword=NULL,*notebook=NULL;
  char *meshfile=NULL, *depthfile=NULL,*rootname=NULL,*path=NULL, *gridfile=NULL, *griddedvar[3]={NULL,NULL,NULL}, *zone=NULL;
  char *convention=0,*format=0;
  string MappingString;
  char *wave[1024];
  grid_t grid;
  grid_t cartesian_zgrid;
  mesh_t mesh;
  int nwave=0,analysis=-1;
  geo_t projection;
  frame_t frame,prescribed;
  double dx=NAN,dy=NAN;
  int persistence=0;
  string grid_options="";
  vector<string *> varnames;
  vector<string *> outnames;
  signed char *mask=0;
  foptions_t foptions;
  vector<string> input,output,discretisation;
  poc_var_t grid_var;

  cmd=fct_echo( argc, argv);
  
  foptions.create=true;

  n=1;
  while (n < argc) {
    keyword=strdup(argv[n]);
    if(read_zone_arg(keyword,argv[n+1],&prescribed,&dx,&dy)) {
      n++;
      n++;
      continue;
      }
    if(strcmp(keyword,"--append")==0) {
      foptions.create=false;
      foptions.append=true;
      n++;
      continue;
      }
    if(strcmp(keyword,"--rename")==0) {
      string *tmp=new string[2];
      for(k=0;k<2;k++){
        if(argv[n+1+k][0]==0){
          k++;
          break;
          }
        tmp[k]=argv[n+1+k];
        }
      outnames.push_back(tmp);
      n++;
      n+=k;
      continue;
      }
    if(strcmp(keyword,"--persistence")==0) {
      sscanf(argv[n+1],"%d",&persistence);
      n++;
      n++;
      continue;
      }
    if(strcmp(keyword,"--help")==0) {
      print_help(argv[0]);
      exit(0);
      }
    if(strstr(argv[n],"--grid")!=0){
      grid_options=strstr(argv[n],"--grid")+7;
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
//         if(strcmp(keyword,"-discretisation")==0) {
//           discretisation= strdup(argv[n+1]);
//           n++;
//           n++;
//           break;
//           }
        switch (keyword[1]) {
/*------------------------------------------------------------------------------
        naming convention for tidal atlases*/
        case 'c' :
          convention= strdup(argv[n+1]);
          n++;
          n++;
          break;
          
        case 'g' :
          n++;
          gridfile= strdup(argv[n]);
          n++;
/* ----------------------------------------------------------------------------
          empty-terminated list of up to 3 variable names */
          for(k=0;k<3;k++){
            if(argv[n+k][0]==0){
              k++;
              break;
              }
            griddedvar[k]= strdup(argv[n+k]);
            }
          n+=k;
          break;

        case 'h' :
          print_help(argv[0]);
          exit(0);

        case 'z' :
          zone= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'n' :
          notebook= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'm' :
          meshfile= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'd' :
          depthfile= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'a' :
          sscanf(argv[n+1],"%d",&analysis);
          n++;
          n++;
          break;

        case 'o' :
          rootname= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'f' :
          format= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'p' :
          path= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'v' :
/* ----------------------------------------------------------------------------
          empty-terminated list of 1 (real) up to 2 (complex) variable names */
          {
          string *tmp=new string[2];
          for(k=0;k<2;k++){
            if(argv[n+1+k][0]==0){
              k++;
              break;
              }
            tmp[k]=argv[n+1+k];
            }
/* ----------------------------------------------------------------------------
          discretisation */
          discretisation.push_back(argv[n+1+k]);
          varnames.push_back(tmp);
          n++;
          n++;
          n+=k;
          }
          break;

        default:
          STDOUT_BASE_LINE("unknown option %s\n",keyword);
          print_help(argv[0]);
          exit(-1);
        }
        break;

      default:
          wave[nwave]= strdup(argv[n]);
          printf("input wave=%s\n",wave[nwave]);
          nwave++;
          n++;
        break;
      }
      free(keyword);
    }

//   if(format!=0 && strcmp(format,"ascii")!=0 && discretisation==0) {
//     fprintf(stderr,"*** discretisation missing even though input format is NetCDF ***\n");
//     print_help(argv[0]);
//     exit(-1);
//     }

  if(nwave==0) {
    fprintf(stderr,"*** empty tidal wave list ***\n");
    print_help(argv[0]);
    exit(-1);
    }

  if(rootname ==NULL) rootname=strdup("generic");
  if(path ==NULL) path=strdup(".");
  if(format==0) format=strdup("ascii");
  
/*----------------------------------------------------------------------------
  Read FE mesh */
  if(meshfile != NULL) {
    printf("#################################################################\n");
    printf("load mesh from ascii file : %s\n",meshfile);
    status=fe_readmesh(meshfile,MESH_FILE_FORMAT_TRIGRID,&mesh);
    if(status) TRAP_ERR_EXIT(status,"fe_readmesh() error\n");
    status=fe_list(&mesh);
    if(status) TRAP_ERR_EXIT(status,"fe_list() error\n");
    }
  else {
    if(strcmp(format,"ascii")==0) {
      fprintf(stderr,"*** no mesh file specified ***\n");
      print_help(argv[0]);
      exit(-1);
      }
    }

  for (k=0;k<nwave;k++) {
    char *tmp;
    status=tide_decode_atlasname(path, convention, wave[k], 0, &tmp);
    input.push_back(tmp);
    }
  for (k=0;k<nwave;k++) {
    char *tmp=new char[1024];
    sprintf(tmp,"%s.%s.nc",wave[k],rootname);
    output.push_back(tmp);
    }


/*----------------------------------------------------------------------------
  Regular grid prescription: priority sequence*/
  if(grid_options!="") {
    printf("#################################################################\n");
    printf("load grid from meta option: %s\n",grid_options.c_str());
    metagrid_t meta;
    poc_data_t<int8_t> t_mask,u_mask,v_mask;
    
    status= parse_GridOptions(grid_options, meta);
    
    if(meta.target!="")
      status=process_GridOptions(meta.gridfile, meta.target.c_str(), 0, grid, mask);
    else if(meta.u_gnames.vlon!=0)
      status=process_GridOptions(meta.gridfile, meta.u_gnames.vlon, meta.u_gnames.vmask, grid, mask);
    else if(meta.v_gnames.vlon!=0)
      status=process_GridOptions(meta.gridfile, meta.v_gnames.vlon, meta.v_gnames.vmask, grid, mask);
    else if(meta.z_gnames.vlon!=0)
      status=process_GridOptions(meta.gridfile, meta.z_gnames.vlon, meta.z_gnames.vmask, grid, mask);
    else if(meta.f_gnames.vlon!=0)
      status=process_GridOptions(meta.gridfile, meta.f_gnames.vlon, meta.f_gnames.vmask, grid, mask);
    else
      STDOUT_BASE_LINE("unable to load notebook file: %s\n",notebook);
    }
  else if(notebook!=0) {
/*----------------------------------------------------------------------------
    Read notebook data */
    printf("#################################################################\n");
    printf("load grid from notebook file : %s\n",notebook);
    status=load_notebook(notebook, &cartesian_zgrid, &grid, &projection);
    printf("%s (notebook file) processed\n",notebook);
    if(status!=0) {
      STDOUT_BASE_LINE("unable to load notebook file: %s\n",notebook);
      exit(-1);
      }
    }
  else if(gridfile != NULL) {
/*----------------------------------------------------------------------------
    Read grid file */
    printf("#################################################################\n");
    printf("load regular grid from %s for",gridfile);
//     for(k=0;k<3;k++){
//       if(griddedvar[k]==0) break;
//       printf(" %s",griddedvar[k]);fflush(stdout);
//       status=poc_get_grid(gridfile,griddedvar[k],&grid[k],1);
//       if(status) NC_TRAP_ERROR(wexit,status,1,"poc_get_grid(\"%s\",\"%s\",) error",gridfile,griddedvar[k]);
//       if(!grid[k].nx || !grid[k].ny) TRAP_ERR_EXIT(status,"empty grid\n");
//       switch(k){
//         case 0:printf(",");break;
//         case 1:printf(" and");break;
//         case 2:printf("\n");break;}
//       }
    }
  else if(zone != NULL) {
/*----------------------------------------------------------------------------
    Apply zone definition by name*/
    printf("#################################################################\n");
    printf("define grid from zone definition : %s\n",zone);
    grid=get_zonegrid(zone);
/**----------------------------------------------------------------------------
    allow for lighter grid storage */
//    status=map_completegridaxis(&grid,2);
    status=map_completegridaxis(&grid,1);
    }
  else if(MappingString!="") {
    status=map_DecodeGrid(MappingString, grid);
    if(status!=0) TRAP_ERR_EXIT(-1,"mapping string not understood (may be you forgot \"): %s\n",MappingString.c_str());
    status=fe_minmax(mesh, frame);
    if(isnan(grid.xmin)) grid.xmin=frame.xmin;
    if(isnan(grid.ymin)) grid.ymin=frame.ymin;
    if(isnan(grid.xmax)) grid.xmax=frame.xmax;
    if(isnan(grid.ymax)) grid.ymax=frame.ymax;
    status=map_set2Dgrid(&grid, grid.xmin, grid.ymin, grid.xmax, grid.ymax, grid.dx, grid.dy);
    status=map_completegridaxis(&grid, 1);
    }
  else {
    if(!mesh.vertices){
//       status=tide_decode_atlasname(path,convention,wave[0], 0,&meshfile);
//       printf("#################################################################\n");
      printf("load mesh geometry from %s\n",input[0].c_str());
      status=fe_readgeometry(input[0].c_str(), &mesh);
      if(status) NC_TRAP_ERROR(return,status,1,"fe_readgeometry(\"%s\",) error",meshfile);
      }
/*----------------------------------------------------------------------------
    Apply frame definition*/
    printf("#################################################################\n");
    printf("define grid from arguments: ");fflush(stdout);
    
    status=fe_minmax(mesh, frame);
    if(status)TRAP_ERR_EXIT(status,"fe_minmax() error\n");
    mesh.destroy();
    
/*----------------------------------------------------------------------------
    default frame (mesh limits)*/
    grid.xmin=frame.xmin;
    grid.ymin=frame.ymin;
    grid.xmax=frame.xmax;
    grid.ymax=frame.ymax;
    
/*----------------------------------------------------------------------------
    apply prescribed limits and resolution*/
    apply_zone_arg(&grid,prescribed,dx,dy);
    grid.modeH=0;
    
    grid.brief_print();
    
    status=map_completegridaxis_2(&grid);
    }

/**----------------------------------------------------------------------------
  save structured grid(s) */
  printf("#################################################################\n");
  for (k=0;k<nwave;k++) {
    printf("creating file %s\n", output[k].c_str());
    string loc;
    int create=(foptions.create);
    string production;
    int compression=0;
    int verbose=0;
/*------------------------------------------------------------------------------
    check wether file exist or not*/
    FILE *in;
    bool file_exist;
    in = fopen(output[k].c_str(), "r");
    if(in == NULL) {
      file_exist=false;
      }
    else {
      file_exist=true;
      fclose(in);
      }
    
    production="gridit-versatile";
    
    if(!file_exist or foptions.overwrite) {
      status=poc_create(output[k], production, verbose, compression);
      }
    
    if(grid_options!="") {
//     printf("#################################################################\n");
//     printf("load grid from meta option: %s\n",grid_options.c_str());
      metagrid_t meta;
      poc_data_t<int8_t> t_mask,u_mask,v_mask;
    
      status= parse_GridOptions(grid_options, meta);
    
//       if(meta.target!="")
//         status=process_GridOptions(meta.gridfile, meta.target.c_str(), 0, grid, mask);
//       else if(meta.u_gnames.vlon!=0)
      if(meta.u_gnames.vlon!=0)
        status=duplicate_grid(meta.gridfile, meta.u_gnames.vlon, meta.u_gnames.vlat, output[k], grid_var);
      else if(meta.v_gnames.vlon!=0)
        status=duplicate_grid(meta.gridfile, meta.v_gnames.vlon, meta.v_gnames.vlat, output[k], grid_var);
      else if(meta.z_gnames.vlon!=0)
        status=duplicate_grid(meta.gridfile, meta.z_gnames.vlon, meta.z_gnames.vlat, output[k], grid_var);
      else if(meta.f_gnames.vlon!=0)
        status=duplicate_grid(meta.gridfile, meta.f_gnames.vlon, meta.f_gnames.vlat, output[k], grid_var);
      else
        STDOUT_BASE_LINE("unable to load notebook file: %s\n",notebook);
      
      if(meta.target!="") {
        status=poc_inq_var(meta.gridfile, meta.target.c_str(), &grid_var);
        }
      }
    else {
      status=poc_save_grid(output[k],&grid_var,grid,create,1,loc);
      status=poc_def_att(output[k],poc_att_t("production","constructed around " __LINE_FILE_PACKAGE_REVISION));
      status=poc_def_att(output[k],poc_att_t("history",cmd));
      }
    }

  if(strcmp(format,"ascii")==0) {
    printf("#################################################################\n");
    printf("treat ascii files\n");
    
//     status=gridit_ASCII(path, meshfile, depthfile, grid, wave, nwave, atlas_convention, analysis, rootname);
    status=gridit_NETCDF_vector1D(mesh, varnames, outnames, grid, mask, wave, nwave, analysis, input, output, format, grid_var, foptions, discretisation, persistence);
    if(status) TRAP_ERR_EXIT(status,"gridit_ASCII() error\n");
    }
  else {
    printf("#################################################################\n");
    printf("treat netcdf files\n");
    status=gridit_NETCDF_vector1D(mesh, varnames, outnames, grid, mask, wave, nwave, analysis, input, output, format, grid_var, foptions, discretisation, persistence);
    if(status) TRAP_ERR_EXIT(status,"gridit_NETCDF_vector1D() error\n");
    }
  
  return 0;
}
