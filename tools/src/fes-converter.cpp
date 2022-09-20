


/**************************************************************************

  T-UGO tools, 2006-2013

  Unstructured Ocean Grid initiative

***************************************************************************/
/** \file

\author  Florent Lyard      LEGOS/CNRS, Toulouse, France. florent.lyard@legos.obs-mip.fr
\author  Yoann Le Bars      LEGOS, Toulouse, France (PhD)
\author  Damien Allain      LEGOS/CNRS, Toulouse, France
\author  Sara Fleury        LEGOS/CNRS, Toulouse, France
\author  Clément MAYET      LEGOS, Toulouse, France (PhD)
\author  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
\author  Frédéric Dupont    Université de Laval à Québec, Canada

<!-- USE FIND AND REPLACE TO EDIT THIS LINE, SO THAT print_help IS ALSO UPDATED -->
\brief Creates an open boundary file for TUGOm.

<!-- A LINK TO main() or print_help() WILL NOT LINK TO THE RIGHT SOURCE ! -->
See the main function for how this works
and the print_help function for how to use this.
*/
/*----------------------------------------------------------------------------*/

#define MAIN_SOURCE

#include "version-macros.def" //for VERSION and REVISION

#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cmath>

#include "tools-structures.h"

#include "functions.h"
#include "map.h"
#include "fe.h"
#include "tides.h"
#include "archive.h"
#include "poc-netcdf.hpp"



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
    "  %s input [OPTIONS]\n",prog_name);
  printf("\n"
    "DESCRIPTION\n"
/* USE FIND AND REPLACE TO EDIT THE LINE BELOW, SO THAT THE FILE HEAD IS ALSO UPDATED */
    "  Converts and interpolates structured tidal atlases.\n"
    "\n"
    "OPTIONS\n"
    "  -h,--help  show this help and exit\n"
    "  -cm  input is in cm instead of m\n"
    "  -mm  input is in mm instead of m\n"
    "  -bmg  set output format to bmg\n"
    "  -m   followed by mask file\n"
    "  -b   followed by topography file. Usefull to use instead of mask file.\n"
    "  -debug  enable debug mode when reading topography or mask file\n"
    "  -o   followed by output file name\n"
    "  -w   followed by wave name. Targeted wave if multi-waves input\n"
    "  -f   followed by input format : bmg, ascii, netcdf, got, tpxo, otis, dtu (equivalent to ascii), eot, osu or cst\n"
    "       If the extension of the input is .nc, then the default is netcdf,\n"
    "       if the extension of the input is .asc, then the default is ascii,\n"
    "       otherwise it is bmg.\n"
    "  -uv  Read 2 buffers from atlas. Implies ascii input format.\n"
    "  -s   followed by smoothing factor in m. Default is no smoothing.\n"
    "\n"
    "The following options are for modifying the ouput grid. By default, the ouput grid is the same as the input grid.\n"
    "\n"
    );
  print_zone_arg_help();
  printf(
    "  -g   followed by grid file. The default grid is that of the input.\n"
    "  -v   followed by grid variable. For more information see option -g .\n"
    "  -z   followed by grid or zone name or resolution in degrees\n"
    "  -extend add one column to repeat 0/360 longitude\n"
    "\n"
    "SEE ALSO: gridit-cdf\n"
    ); /** \endcode */
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int read_UG2D(const char *filename, const string & varname, mesh_t & mesh, complex<float>* & tide, complex<float> &mask, const char *discretisation,int iteration, bool init, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, j, k, l, n,node,frame, status;
  float spec[2];
  fcomplex zz,*z=NULL,cmask;
  char units[256];
  string amplitude,phase;
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
  if(init) {
    status=fe_readgeometry(filename, &mesh);
    status=fe_readdiscretisation(filename, &mesh, 0, id);
    }
    
  descriptor=get_descriptor(mesh,id);
  
  if(tide==0) tide=new complex<float>[descriptor.nnodes];
  
/*-----------------------------------------------------------------------
  load netcdf variable */
  frame=iteration;
  
  amplitude="a_"+varname+"_"+discretisation;
  phase="G_"+varname+"_"+discretisation;
  
  status=poc_get_UG3D(filename, frame, amplitude.c_str(), phase.c_str(), tide);
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

  int process(string assimilated, string & rootname, const char *target, string out)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  mesh_t mesh;
  float *buffer, mask, scale, offset;
  string filename, output;
  complex<float> *LGP2buffer=0,*DNP1buffer=0,*LGP1buffer=0,cmask;
  int verbose=0;
  int iteration=-1;
  string discretisation;
  vector<string> sequential;
  vector<string> spectral;
  vector<string> optimised;
  string vname[2],varname;
  string in;
  bool init;
  int id;
  string extension="."+rootname+".nc";

/*------------------------------------------------------------------------------
  LGP1xLGP1 sequential, no atmospheric forcing (7) */
  sequential.push_back("Mf");
  sequential.push_back("Mm");
  sequential.push_back("MSf");
  sequential.push_back("MSqm");
  sequential.push_back("Mtm");
  sequential.push_back("Sa");
  sequential.push_back("Ssa");
  
/*------------------------------------------------------------------------------
  LGP1xLGP1 sequential, with atmospheric forcing (6) */
  sequential.push_back("J1");
  sequential.push_back("M3");
  sequential.push_back("M8");
  sequential.push_back("MKS2");
  sequential.push_back("N4");
  sequential.push_back("S1");
  sequential.push_back("S4");
//   sequential.push_back("R2");
//   sequential.push_back("T2");
  
/*------------------------------------------------------------------------------
  DNP1xLGP2 spectral, hydrodynamic (5) */
  spectral.push_back("M6");
  spectral.push_back("MN4");
  spectral.push_back("MS4");
  spectral.push_back("R2");
  spectral.push_back("T2");
  
/*------------------------------------------------------------------------------
  DNP1xLGP2 spectral, optimised (15) */
  optimised.push_back("2N2");
  optimised.push_back("L2");
  optimised.push_back("Mu2");
  optimised.push_back("P1");
  optimised.push_back("E2");
  optimised.push_back("La2");
  optimised.push_back("N2");
  optimised.push_back("Q1");
  optimised.push_back("K1");
  optimised.push_back("M2");
  optimised.push_back("Nu2");
  optimised.push_back("K2");
  optimised.push_back("M4");
  optimised.push_back("O1");
  optimised.push_back("S2");
  
  in="/home/data/tides/FES2014/distribution_fes2014a/HYDRO/unstructured/";
//   out="/home/data/tides/FES2014/distribution_fes2014b/synthesis/unstructured/";
  
  filename=in+"Mf.sequential.nc";
  
  iteration=-1;
  
  status=fe_readgeometry(filename.c_str(), &mesh);
  
//   id=discretisation_from_name(discretisation.c_str());
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  generate needed discretisations
  
  does assume archived discretisation and generated ones will be identical...

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  int compute_massmatrix=1, context=SEQUENTIAL_COMPUTING;
  
  status=discretisation_init(&mesh, LGP2, compute_massmatrix, context);
  LGP2buffer=new complex<float>[mesh.LGP2descriptor.nnodes];
  
  status=discretisation_init(&mesh, DNP1, compute_massmatrix, context);
  DNP1buffer=new complex<float>[mesh.DNP1descriptor.nnodes];
  
  status=discretisation_init(&mesh, LGP1, 0, context);
  LGP1buffer=new complex<float>[mesh.LGP1descriptor.nnodes];
  
  cmask=fcomplex(9999.,9999.);
  
/*------------------------------------------------------------------------------
  LGP1xLGP1 sequential, hydrodynamic */
//   status=fe_readdiscretisation(filename.c_str(), &mesh, 0, LGP1);
  
  in="/home/data/tides/FES2014/distribution_fes2014a/HYDRO/unstructured/";
  for(int k=0;k<sequential.size();k++) {
    filename=in+sequential[k]+".sequential.nc";
    printf("load harmonic file : %s\n",filename.c_str());
    
    discretisation="LGP1";
    
    varname="eta";
    
    vname[0]="a_"+varname+"_"+discretisation;
    vname[1]="G_"+varname+"_"+discretisation;
  
    status=poc_get_UG3D(filename.c_str(), iteration, vname[0].c_str(), vname[1].c_str(), LGP1buffer);
    if(status !=0) {
      printf("cannot load frame %d in %s\n",iteration,filename.c_str());
      }
    status=fe_projection( mesh, LGP1buffer, LGP1, LGP2buffer, LGP2);
//   variable = poc_variable_UG2D("name_undefined", mask,"units_undefined", scale, offset,"standardname_undefined","N");
//   status   = cdf_createvariable(output, &(variable));
//   status   = poc_put_UG2D(output, mesh, variable.id, buffer);
    output=out+sequential[k]+extension;
    remove(output.c_str());
    
    int option=1;
    bool edges=false, levels=false, code=false, bathymetry=false;
    status=fe_savemesh3d(output.c_str(), mesh, option, edges, levels, code, bathymetry);
    
    status=archiving_UGdummy2D(output.c_str(), mesh, "a_eta_LGP2", "G_eta_LGP2", "m",   LGP2buffer, cmask, LGP2);
    
    init=false;
    
    varname="u";
    
    vname[0]="a_"+varname+"_"+discretisation;
    vname[1]="G_"+varname+"_"+discretisation;
  
    status=poc_get_UG3D(filename.c_str(), iteration, vname[0].c_str(), vname[1].c_str(), LGP1buffer);
    status=fe_projection( mesh, LGP1buffer, LGP1, DNP1buffer, DNP1);
    status=archiving_UGdummy2D(output.c_str(), mesh, "a_u_DNP1", "G_u_DNP1", "m",   DNP1buffer, cmask, DNP1);
    
    varname="v";
    
    vname[0]="a_"+varname+"_"+discretisation;
    vname[1]="G_"+varname+"_"+discretisation;
  
    status=poc_get_UG3D(filename.c_str(), iteration, vname[0].c_str(), vname[1].c_str(), LGP1buffer);
    status=fe_projection( mesh, LGP1buffer, LGP1, DNP1buffer, DNP1);
    status=archiving_UGdummy2D(output.c_str(), mesh, "a_v_DNP1", "G_v_DNP1", "m",   DNP1buffer, cmask, DNP1);
    }
    
/*------------------------------------------------------------------------------
  DNP1xLGP2 spectral, hydrodynamic */
  in="/home/data/tides/FES2014/distribution_fes2014a/HYDRO/unstructured/";
  for(int k=0;k<spectral.size();k++) {
    filename=in+spectral[k]+".spectral.nc";
    printf("load harmonic file : %s\n",filename.c_str());
    
    discretisation="LGP2";
    
    varname="eta";
    
    vname[0]="a_"+varname+"_"+discretisation;
    vname[1]="G_"+varname+"_"+discretisation;
    
    status=poc_get_UG3D(filename.c_str(), iteration, vname[0].c_str(), vname[1].c_str(), LGP2buffer);
    output=out+spectral[k]+extension;
    remove(output.c_str());
    status=archiving_UGdummy2D(output.c_str(), mesh, "a_eta_LGP2", "G_eta_LGP2", "m",   LGP2buffer, cmask, LGP2);
        
    discretisation="DNP1";
    
    varname="u";
    
    vname[0]="a_"+varname+"_"+discretisation;
    vname[1]="G_"+varname+"_"+discretisation;
    
    status=poc_get_UG3D(filename.c_str(), iteration, vname[0].c_str(), vname[1].c_str(), DNP1buffer);
    status=archiving_UGdummy2D(output.c_str(), mesh, "a_u_DNP1", "G_u_DNP1", "m",   DNP1buffer, cmask, DNP1);

    varname="v";
    
    vname[0]="a_"+varname+"_"+discretisation;
    vname[1]="G_"+varname+"_"+discretisation;

    status=poc_get_UG3D(filename.c_str(), iteration, vname[0].c_str(), vname[1].c_str(), DNP1buffer);
    status=archiving_UGdummy2D(output.c_str(), mesh, "a_v_DNP1", "G_v_DNP1", "m",   DNP1buffer, cmask, DNP1);
    }
    
/*------------------------------------------------------------------------------
  DNP1xLGP2 spectral, optimised */
//   in="/home/data/tides/FES2014/distribution_fes2014b/ASSIMILE+COURANTS/unstructured/";
  in=assimilated;
  for(int k=0;k<optimised.size();k++) {
    filename=in+optimised[k]+".EnOI.nc";
    printf("load harmonic file : %s\n",filename.c_str());
    
    discretisation="LGP2";
    
    varname="eta";
    
    vname[0]="analysis_Ha";
    vname[1]="analysis_Hg";
  
    status=poc_get_UG3D(filename.c_str(), iteration, vname[0].c_str(), vname[1].c_str(), LGP2buffer);
    output=out+optimised[k]+extension;
    remove(output.c_str());
    status=archiving_UGdummy2D(output.c_str(), mesh, "a_eta_LGP2", "G_eta_LGP2", "m",   LGP2buffer, cmask, LGP2);
        
    discretisation="DNP1";
    
    varname="u";
    
    vname[0]="analysis_Ua";
    vname[1]="analysis_Ug";
    
    status=poc_get_UG3D(filename.c_str(), iteration, vname[0].c_str(), vname[1].c_str(), DNP1buffer);
    status=archiving_UGdummy2D(output.c_str(), mesh, "a_u_DNP1", "G_u_DNP1", "m",   DNP1buffer, cmask, DNP1);

    varname="v";
    
    vname[0]="analysis_Va";
    vname[1]="analysis_Vg";

    status=poc_get_UG3D(filename.c_str(), iteration, vname[0].c_str(), vname[1].c_str(), DNP1buffer);
    status=archiving_UGdummy2D(output.c_str(), mesh, "a_v_DNP1", "G_v_DNP1", "m",   DNP1buffer, cmask, DNP1);
    }
    
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double x,y;
  float  time=0,scale=1.0;
  float  spec=-9999;
  int i,j,k,n;
  int nbuffers=1,status;
  char *keyword=NULL,*gridfile=NULL,*gridvar=NULL,*maskfile=NULL;
  double dx=NAN,dy=NAN;
  char *output=NULL,*input=NULL;
  char *wave=0;
  string format;
  spectrum_t spectrum;
  string Aname, Gname;
  string rootname;
  string assimilated="/home/data/tides/FES2014/distribution_fes2014b/ASSIMILE+COURANTS/unstructured/";
  string out="/home/data/tides/FES2014/distribution_fes2014b/synthesis/unstructured/";

  bool debug=false;
  rootname="FES2014-unified";
  
  string cmd=fct_echo( argc, argv);

  n=1;
  while (n < argc) {
    keyword=argv[n];
    if(strcmp(keyword,"--help")==0) {
      print_help(argv[0]);
      exit(0);
      }
    if(strcmp(keyword,"-debug")==0) {
      debug=true;
      n++;
      continue;
      }
    if(strcmp(keyword,"--assimilated")==0) {
      assimilated=argv[n+1];
      n++;
      n++;
      continue;
      }
    if(strcmp(keyword,"--output")==0) {
      out=argv[n+1];
      n++;
      n++;
      continue;
      }
    switch (keyword[0]) {
      case '-':
        switch (keyword[1]) {
        case 'h' :
          print_help(argv[0]);
          exit(0);

        case 'g' :
          gridfile= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'v' :
          gridvar= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'm' :
          maskfile= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'o' :
          output= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'w' :
          wave=strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'f' :
          format=argv[n+1];
          n++;
          n++;
          break;

        case 'r' :
          rootname=argv[n+1];
          n++;
          n++;
          break;


        default:
          STDOUT_BASE_LINE("unknown option %s\n",keyword);
          print_help(argv[0]);
          exit(-1);
        }
        break;

      default:
        input= strdup(argv[n]);
        n++;
        break;
      }
    }
  
  fe_integrale_init();

  status=process(assimilated, rootname, "LGP2", out);
    
}
