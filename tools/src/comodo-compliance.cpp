
/*******************************************************************************

  T-UGO tools, 2006-2011

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

\author  Florent Lyard      LEGOS/CNRS, Toulouse, France. florent.lyard@legos.obs-mip.fr
\author  Yoann Le Bars      LEGOS, Toulouse, France (PhD)
\author  Laurent Roblou     LEGOS/CNRS, Toulouse, France
\author  Damien Allain      LEGOS/CNRS, Toulouse, France
\author  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
\author  Frédéric Dupont    Université de Laval à Québec, Canada

<!-- USE FIND AND REPLACE TO EDIT THIS LINE, SO THAT print_help IS ALSO UPDATED -->
\brief Test compliance of NetCDF files to comodo standard.

<!-- A LINK TO main() or print_help() WILL NOT LINK TO THE RIGHT SOURCE ! -->
See the main function for how this works
and the print_help function for how to use this.
*/
/*----------------------------------------------------------------------------*/

#define MAIN_SOURCE
#include "version-macros.def" //for VERSION and REVISION

#include <stdio.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include "tools-structures.h"

#include "functions.h"
#include "poc-netcdf-data.hpp"
#include "poc-time.h"
#include "map.h"


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
    "  %s [OPTIONS] file1 var11 [var12 ... ] [ file2 var21 [var22 ... ] ... ]\n",prog_name);
  printf("\n"
    "DESCRIPTION\n"
/* USE FIND AND REPLACE TO EDIT THE LINE BELOW, SO THAT THE FILE HEAD IS ALSO UPDATED */
    "  Test compliance of NetCDF files to comodo standards : " COMODO_URL "\n"
    "  It takes a file followed by variable names and tries and get the grid for each variable.\n"
    "  If the file is a .grd file, the variable names are prepended with z .\n"
    "\n"
    "OPTIONS :\n"
    "  -h,--help  show this help and exit\n"
    "  -v  verbose mode\n"
    "\n"
    "TIPS\n"
    "  Use NCOs to correct your files.\n"
    "  Change attributes with : ncatted -a att,var,mode,type,new_value [-a ...] [-O] file [optional_new_file]\n"
    "Replacing :\n"
    "- att and var as the attribute name and variable regular expression (`global' for global attributes, blank for all variables)\n"
    "- mode with one of a,c,d,m,o for append, create, delete, modify (no effect if attribute does not exist), overwrite (default)\n"
    "- type with one of f,d,l,s,c,b for float, double, short, char, byte.\n"
    "Option -O is for overwriting the new file.\n"
    " e.g. : ncatted -a \"coordinates,ssh_.,o,char,lat lon\" -a \"content,ssh_.,o,char,YX\" -O M2-ssh-atlas.nc M2-ssh-atlas-.nc\n"
    );
  /** \endcode */
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int aC, char *as[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// function that will be ran when the executable is started
/** See the print_help <a href=#func-members>function</a> for its use.
*/
/*----------------------------------------------------------------------------*/
{
  int i,status;
  char *a;//argument
  vector<string> files,*varList=NULL;
  vector<vector<string> > varLists;
  string report;
  int verbose=0;
  
  fct_echo( aC, as);

  i=1;
/*---------------------------------------------------------------------*//**<h1>
  Parse the argument list </h1>*/
  while (i < aC) {
    a=as[i];
    if(strcmp(a,"--help")==0) {
      print_help(as[0]);
      return 0;
      }
    switch (a[0]) {
      case '-':
        switch (a[1]) {

        case 'h' :
          print_help(as[0]);
          return 0;

        case 'v' :
          verbose++;
          i++;
          break;

        default:
          printf("unknown option %s\n",a);
          print_help(as[0]);
          return -1;
        }
        break;

      default:
/*---------------------------------------------------------------------------*/
        if(strrncasecmp(a,".nc")==0){
          /* Assume it is a file and add it to the file list */
          files.push_back(a);
          varLists.push_back(vector<string>());
          varList=&varLists.back();
          }
/*---------------------------------------------------------------------------*/
        else if(strrncasecmp(a,".grd")==0){
          /* Assume it is a file and add it to the file list */
          files.push_back(a);
          varLists.push_back(vector<string>());
          varList=&varLists.back();
          varList->push_back("z");
          }
/*---------------------------------------------------------------------------*/
        else{
          if(varList==NULL){/// <h3>If it is the first argument</h3>
            /// abort, calling print_help().
            __ERR_BASE_LINE__("*** You must have a file name before the variable name ***\n");
            print_help(as[0]);
            return -1;
            }
          /// <h3>Otherwise, add it to the current variable list</h3>
          varList->push_back(a);
          }
        i++;
        break;
      }
    }
  
  if(files.size()<1){
    print_help(as[0]);
    return -1;
    }
  
  for(i=0;i<files.size();i++){/// <h1>For all the files</h1>
    ///check compliance with comodo_compliance()
    status=comodo_compliance(files[i],varLists[i],&report,verbose);
    printf(report);
    }
  
  return 0;
}
