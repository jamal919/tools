
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
\brief Show build options for programmes using tools library

<!-- A LINK TO main() or print_help() WILL NOT LINK TO THE RIGHT SOURCE ! -->
See the main function for how this works
and the print_help function for how to use this.
*/
/*----------------------------------------------------------------------------*/

#define MAIN_SOURCE

#include "functions.h"
#include "poc-assertions.h"

#include "tools-define.h"
#include "version-macros.def"

#include "flags.def"

#include <stdio.h>
#include <errno.h>
#include <string.h>

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
    "  %s OPTION\n",prog_name);
  printf("\n"
    "DESCRIPTION\n"
/* USE FIND AND REPLACE TO EDIT THE LINE BELOW, SO THAT THE FILE HEAD IS ALSO UPDATED */
    "  Show build options for programmes using tools library.\n"
    "\n"
    "OPTION\n"
    "  -h,--help  show this help and exit\n"
    "  -v         show version. TIP: compare with m4_version_compare() (See info:/autoconf/Number processing Macros) in your configure.ac scripts.\n"
    "  -i         show include options\n"
    "  -l         show link options\n"
    "  -d         If uninstalled, give path to libtools.a . This is usefull for dependencies checks in Makefile's. Gives nothing if installed.\n"
    "  -dT        If uninstalled, give path to lib4tugo.a . See note above.\n"
    "\n"
    "WARNING\n"
    "  It will detect whether it is installed. Unless so, it will add the tools sources or build directory respectively to the include or the link options and it will also print relevant warnings.\n"
    "Please note that if you:\n"
    "- are using a hashing shell like bash or csh and\n"
    "- are calling this programme without an absolute path to test what it does and\n"
    "- also have the tools build directory in your PATH\n"
    "then you need to ask your shell to redo its hash table of executables with the rehash command (for bash or csh) right after a make install.\n"
    "For more information, see your shell's documentation.\n"
    ); /** \endcode */
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// function that will be ran when the executable is started
/** See the print_help <a href=#func-members>function</a> for its use.
*/
/*----------------------------------------------------------------------------*/
{
/*---------------------------------------------------------------------*//**<h1>
  checks whether the programme is installed </h1>*/
  string execPath=canonicalize_file_name((string)"/proc/self/exe");
  int afterLastSep=execPath.rfind("/")+1;
  string name=execPath.substr(afterLastSep);
  string installPath=canonicalize_file_name(TOOLS_BINDIR "/"+name);
  string buildPath=canonicalize_file_name(TOOLS_BUILDDIR "/"+name);
//   STDERR_BASE_LINE(execPath+"\n");
//   STDERR_BASE_LINE(installPath+"\n");
//   STDERR_BASE_LINE(buildPath+"\n");
  int isInstalled=execPath==installPath;
  if(!isInstalled){
    if(execPath==buildPath)
      STDERR_BASE_LINE("warning: "+name+" is not installed!\n");
    else
      STDERR_BASE_LINE("WARNING: "+name+" APPEARS TO HAVE BEEN INSTALLED MANUALLY:\n"
        "DOING AS IF IT WAS NOT INSTALLED!!!\n");
    }
  
/*---------------------------------------------------------------------*//**<h1>
  then </h1>*/
  if(argc!=2){
    if(argc<2) STDOUT_BASE_LINE("*** Missing the argument ***\n");
    else STDOUT_BASE_LINE("*** Only one argument allowed ***\n");
    print_help(argv[0]);
    return 1;
    }
  
  switch (argv[1][0]) {
    
  case '-':
    if(strcmp(argv[1],"-h")==0 || strcmp(argv[1],"--help")==0 ) {
      print_help(argv[0]);
      return 0;
      }
    
    switch (argv[1][1]) {
      
    case 'v' :
      printf(VERSION "\n");
      return 0;
      
    case 'l' :
      if(isInstalled){
        printf("-ltools"  TOOLS_LDFLAGS "\n");
        }
      else{
        STDERR_BASE_LINE("WARNING : GIVING TOOLS BUILD DIRECTORIES AMONG LINK FLAGS!\n");
        printf(TOOLS_BUILDDIR "/libtools.a " TOOLS_LDFLAGS "\n");
        }
      return 0;
      
    case 'd' :
      if(isInstalled){
        }
      else{
        if(argv[1][2]=='\0')
          printf(TOOLS_BUILDDIR "/libtools.a\n");
        else
          printf(TOOLS_BUILDDIR "/lib4tugo.a\n");
        }
      return 0;
      
    case 'i' :
      if(isInstalled){
        printf(TOOLS_INCFLAGS "\n");
        }
      else{
        STDERR_BASE_LINE("WARNING : GIVING TOOLS SOURCES DIRECTORY AMONG INCLUDE FLAGS!\n");
        printf(TOOLS_INCFLAGS "\n");
        }
      return 0;
      
      }
    
    }
  
  TRAP_ERR_EXIT(1,"unknown option %s. Run with --help for help or update your version of the tools!\n",argv[1]);
}
