
/*******************************************************************************
T-UGO tools, 2006-2009

  Unstructured Ocean Grid initiative
  
E-mail: florent.lyard@legos.obs-mip.fr
*******************************************************************************/
/**
\file

\author  Florent Lyard      LEGOS/CNRS, Toulouse, France
\author  Yoann Le Bars      LEGOS, Toulouse, France (PhD)
\author  Laurent Roblou     LEGOS/CNRS, Toulouse, France
\author  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
\author  Frédéric Dupont    Université de Laval à Québec, Canada
\author  Yves Soufflet      LEGOS, Toulouse, France

<!-- USE FIND AND REPLACE TO EDIT THIS LINE, SO THAT print_help IS ALSO UPDATED -->
\brief 
*/
/*----------------------------------------------------------------------------*/

#define MAIN_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <errno.h>

#include <config.h>

#include "tools-define.h"
#include "tools-structures.h"

#include "version-macros.def" //for VERSION and REVISION

#include "fe.h"
#include "poc-time.h"
#include "list.h"
#include "map.h"
#include "mgr.h"
#include "grd.h"
#include "geo.h"
#include "filter.h"
#include "functions.h"
#include "tides.h"
#include "poc-netcdf-data.hpp"
#include "archive.h"

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int tide_decode_atlasname(const char *atlas_directory, const char *atlas_convention,
                            const char *wave, const char *member, char **filename, int mode)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///computes a name from a convention
/**
\param *atlas_directory directory convention : see atlas_convention below for the replacements
\param *atlas_convention file name convention : "WAVE" will be replaced with the name of the wave (which is all upper-case, but see mode below)
\param *wave name of the wave
\param **filename pointer to computed file name. Will be allocated with \b new
\param mode Default: 3.
mode&1 : if the file exists without any replacements made, do not do any replacement and return 0.
mode&2 : replace "wave" with lower-cased name of the wave and "Wave" with all-but-first-letter-lower-cased name of the wave
\return 0 if file is found or -1

\sa print_tide_decode_atlasname_help() decode_atlasname()
*/
/*----------------------------------------------------------------------------*/
{
  int i;
  char dummy[256], *pointer=NULL, *p2=NULL;
  FILE *out=NULL;

  if(atlas_convention==0) {
    STDERR_BASE_LINE("atlas naming template not given, quit \n");
    return(-1);
    }
/*------------------------------------------------------------------------------
  build the atlas file name*/
  *filename=0;
  if(atlas_directory==0) {
    (*filename) = new char[strlen(atlas_convention)+10];
    sprintf((*filename), "%s", atlas_convention);
    }
  else {
    (*filename) = new char[strlen(atlas_directory) + 1 + strlen(atlas_convention)+10];
    sprintf((*filename), "%s/%s", atlas_directory, atlas_convention);
    }

  if(mode&1){
/*------------------------------------------------------------------------------
    check existence name*/
    out = fopen(*filename, "r");
    if(out!=NULL){
/*------------------------------------------------------------------------------
      file exists, do nothing more*/
      fclose(out);
      return 0;
      }
    }

  ///\note Only the first occurrence of WAVE, wave or Wave, in this order, is replaced!
  for(i=0;i<3;i++){
/*------------------------------------------------------------------------------
    use format information*/
    switch(i){
    case 0:
      pointer = strstr((*filename), "WAVE");
      break;
    case 1:
      pointer = strstr((*filename), "wave");
      break;
    case 2:
      pointer = strstr((*filename), "Wave");
      break;
      }
    if(pointer != NULL) {
      p2 =strdup(pointer);
      sprintf(dummy, "%s", wave);
      int k,kMax=strlen(dummy);
      switch(i){
      case 0:
        k=kMax;
        break;
      case 1:
        k=0;
        break;
      case 2:
        k=1;
        break;
        }
      for(;k<kMax;k++) {
        dummy[k]=tolower(dummy[k]);
        }
      strcpy(pointer, dummy);
      strcat(pointer, p2+4);
      free(p2);
      break;
      }
    if(!(mode&2))
      break;
    }
   
  if(member!=0) {
    pointer = strstr((*filename), "MEMBER");
    if(pointer != NULL) {
      p2 =strdup(pointer);
      sprintf(dummy, "%s", member);
      strcpy(pointer, dummy);
      strcat(pointer, p2+6);
      }
    }

  /* check existence */
  if ((out = fopen(*filename, "r"))) {
    fclose(out);
    return 0;
  }
  return -1;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int extract(vector<string> waves, const char* path, const char* input_template, const char* output_template, int discretisation, const char *polygons, int nmembers, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k, m, n, status;
  int verbose=0;
  mesh_t mesh, *splitted;
  int *selected, nselected, *targeted=0;
  char *filename;
  discretisation_t descriptor, child;
  complex<double> *global, *buffer, cmask;
  char member[4];
  

  sprintf(member,"%3.3d",0);
  tide_decode_atlasname(path, input_template, waves[0].c_str(), member, &filename);
  
/*------------------------------------------------------------------------------
  load original mesh */
  status=fe_readmesh3d(filename, &mesh,verbose); 
  fe_init_from_elements(&mesh, verbose);

  status=fe_geometry(&mesh);
  
  status=discretisation_init(&mesh, discretisation, 0);
  
  descriptor=get_descriptor(mesh, discretisation);
  
  if(polygons!=0) {
  printf("#################################################################\n");
    printf("select edges from polygons: %s\n",polygons);
    selected=new int[mesh.nedges];
    for(n=0;n<mesh.nedges;n++) selected[n]=1;
    nselected=fe_selectedges_01(mesh, (char *) polygons, selected);
    }
  
  splitted=fe_split(mesh, selected, targeted, 0, 1, (string) "extraction", debug);

  status=discretisation_init(&(splitted[1]), discretisation, 0);
  
  child=get_descriptor(splitted[1],discretisation);
  
  for (m=0;m<splitted[1].ntriangles;m++) {
    int mm=splitted[1].triangles[m].ancestor;
    for(k=0;k<child.nnpe; k++) {
      n=child.NIbE[m][k];
      int nn=descriptor.NIbE[mm][k];
      child.nodes[n].ancestro=nn;
      }
    }
  
  global=new complex<double> [descriptor.nnodes];
  buffer=new complex<double> [child.nnodes];
  
  for(k=0;k<waves.size();k++) {
    for(m=0;m<nmembers;m++) {
      char *output;
      sprintf(member,"%3.3d",m);
      tide_decode_atlasname(path, input_template, waves[k].c_str(), member, &filename);
/**----------------------------------------------------------------------------
      to be turned more flexible later on*/
      char varname_a[1024], varname_G[1024];
      sprintf(varname_a,"a_eta_%s", discretisation_name(discretisation));
      sprintf(varname_G,"G_eta_%s", discretisation_name(discretisation));
      int frame=0;
      status=poc_get_UG3D(filename, frame, (const char *) varname_a, (const char *) varname_G, global);
//     if(status!=0) goto error;
      for (n=0;n<child.nnodes;n++) {
        int nn=child.nodes[n].ancestro;
        buffer[n]=global[nn];
        }
      tide_decode_atlasname(path, output_template, waves[0].c_str(), member, &output);
      status=archiving_UGdummy2D(output, splitted[1], (const char *) varname_a, (const char *) varname_G, "m/s", buffer, cmask, frame, discretisation);
      delete[] filename;
      delete[] output;
      } 
    }
    
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------
  function that will be ran when the executable is started
  See the print_help <a href=#func-members>function</a> for its use.
----------------------------------------------------------------------*/
{
  int k,n,status;
  char *output=NULL, *keyword=NULL, *mgrfile=NULL;
  char *ordering=NULL,*regions=NULL,*polygons=NULL;
  vector<string> waves;
  int nwave=0;
  char *atlas_directory=NULL,*input_template=NULL,*output_template=NULL,*meshfile=NULL;
  string FrameString;
  int nmembers=0;
  vector<plg_t> limits;
  
  char *format=NULL;

  const int nvars=2;

  char *varnames[nvars]={0,0},*discretisation;
  
  bool debug=true;
  int verbose=0;

  fct_echo( argc, argv);

  n=1;
  while (n < argc) {
    keyword=strdup(argv[n]);
    switch (keyword[0]) {
      case '-':
        if(strcmp(keyword,"-discretisation")==0) {
          discretisation=strdup(argv[n+1]);
          n++;
          n++;
          break;
          }
        if(strcmp(keyword,"--regions")==0) {
          regions= strdup(argv[n+1]);
          n++;
          n++;
          break;
          }
        if(strcmp(keyword,"--polygons")==0) {
          polygons= strdup(argv[n+1]);
          n++;
          n++;
          break;
          }
        if(strncmp("--frame",keyword)==0){
          FrameString= argv[n+1];
          n++;
          n++;
          break;
          }
        switch (keyword[1]) {

/* *----------------------------------------------------------------------
        naming convention for tidal atlases*/
        case 'c' :
          input_template= strdup(argv[n+1]);
          n++;
          n++;
          break;

/* *----------------------------------------------------------------------
        path for tidal atlases*/
        case 'p' :
          atlas_directory= strdup(argv[n+1]);
          n++;
          n++;
          break;

/* *----------------------------------------------------------------------
        output file name*/
        case 'o' :
          output_template= strdup(argv[n+1]);
          n++;
          n++;
          break;

/* *----------------------------------------------------------------------
        number of members*/
        case 'n' :
          sscanf(argv[n+1],"%d",&nmembers);
          n++;
          n++;
          break;

        case 'v' :
          if(varnames[0]){
//             if(defaultVars)
//               defaultVars=0;
//             else{
//               __FILE_LINE__(stdout,"multiple use of -v : the last one overrides the previous.");
//               for(k=0;k<nvars;k++){
//                 free(varnames[k]);
//                 }
//               }
            }
          for(k=0;k<nvars;k++){
            varnames[k]=strdup(argv[n+1+k]);
            }
          n++;
          n+=nvars;
          break;

//         case 'h' :
//           print_help(argv[0]);
//           exit(0);
        default:
//           printf("unknown option %s\n",keyword);
//           print_help(argv[0]);
          exit(-1);
        }
        break;

      default:
/* *----------------------------------------------------------------------
          tidal wave list*/
          waves.push_back(argv[n]);
//          printf("input wave=%s\n",wave[nwave]);
          //spectrum.n=nwave+1;
          nwave++;
          n++;
        break;
      }
    free(keyword);
    }

  if(atlas_directory==NULL)  atlas_directory=strdup(".");
  if(input_template==NULL) {
    printf("*** Please specify atlas filename template with -c ***\n");
//     print_help(argv[0]);
    exit(-1);
    }
  if(output_template==NULL) {
    printf("*** Please specify output filename template with -o ***\n");
//     print_help(argv[0]);
    exit(-1);
    }
  if(strcmp(input_template, output_template)==0) {
    printf("*** Please specify different input/output filename template ***\n");
//     print_help(argv[0]);
    exit(-1);
    }
  
  if(FrameString!="") {
    frame_t frame;
    status=plg_DecodeFrame(FrameString, frame);
    limits.push_back(plg_t(frame, PLG_INIT_SEPARATE));
    status=plg_cartesian((projPJ) 0, limits);
    }

  int discretisation_id=discretisation_from_name(discretisation, verbose);
//   const char *polygons="selection.plg";
  
  status=extract(waves, atlas_directory, input_template, output_template, discretisation_id, polygons, nmembers, debug);
  
  TRAP_ERR_EXIT(0,"%s -computation complete ^^^^^^^^^\n",argv[0]);
}
