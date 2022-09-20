
/*******************************************************************************

  T-UGO tools, 2006-2017

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

<!-- USE FIND AND REPLACE TO EDIT THIS LINE, SO THAT print_help IS ALSO UPDATED -->
\brief 
*/
/*----------------------------------------------------------------------------*/

#define MAIN_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <sstream>


#include "tools-structures.h"

#include "geo.h"
#include "fe.h"
#include "map.h"
#include "polygones.h"
#include "grd.h"
#include "archive.h"
#include "functions.h"

#include "zapper.h"     /*  rutin.h contains common utility routines  */


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int n;
  char *keyword;
  char *meshfile=NULL,*belfile=NULL,*output=NULL;
  mesh_t mesh;
  ostringstream report;
  int incr=1;
  int max_percent=10;

  fct_echo( argc, argv);

  n=1;
  while (n < argc) {
    keyword=strdup(argv[n]);
    switch (keyword[0]) {
      case '-':
        switch (keyword[1]) {
        case 'm' :
/*-----------------------------------------------------------------------------
          UG mesh filename (input)*/
          meshfile= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'b' :
/*-----------------------------------------------------------------------------
          optional bel file filename (output)*/
          belfile= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'o' :
/*-----------------------------------------------------------------------------
          UG mesh filename (output)*/
          output= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'i' :
          sscanf(argv[n+1],"%d",&incr);
          n++;
          n++;
          break;

        case 'p' :
          sscanf(argv[n+1],"%d",&max_percent);
          n++;
          n++;
          break;

        default:
          STDOUT_BASE_LINE("unknown option %s\n",keyword);
          exit(-1);
        }
        break;

      default:
          n++;
        break;
      }
      free(keyword);
    }

  if(meshfile == NULL) {
    printf("no mesh file specified; abort...\n");
//     print_help(argv[0]);
    wexit(-1);
    }

  printf("#################################################################\n");
  printf("load mesh %s\n",meshfile);
  status=fe_readmesh(meshfile,MESH_FILE_FORMAT_TRIGRID,&mesh);
  if(status!=0) {
    printf("unable to read the original mesh in %s\n",meshfile);
//     print_help(argv[0]);
    wexit(-1);
    }

  if(output==0) output=strdup("mesh-renum.reduced.nei");

  report << "#Number of nodes      : " << mesh.nvtxs <<endl;

  printf("#################################################################\n");
  printf("reduce bandwidth with node increment=%d (-i option), max percentage=%d (-p option) \n",incr,max_percent);
  status= fe_reducebw(mesh,incr,max_percent,output);

  if(belfile!=0) {
    printf("#################################################################\n");
    printf("rebuild boundary code file\n");
/* *----------------------------------------------------------------------------
    rebuild codes from edge description*/
    int option=0;
    status=fe_readmesh(output,MESH_FILE_FORMAT_TRIGRID,&mesh);
    status=fe_list(&mesh);
    status= fe_edgetable(&mesh,0,0);
    status=fe_codetable1(&mesh, option,1,0);
    status= fe_write_boundarycode(belfile, mesh, 1);
    }

  STDOUT_BASE_LINE("end of mesh-renum ... \n");
  exit(0);
}
