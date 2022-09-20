
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

/* *@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  topo-merg aims at merging 2 bathymetric database with some controls
  on merging action

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */

#define MAIN_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include "tools-structures.h"
#include "version-macros.def" //for VERSION and REVISION

#include "rutin.h"
#include "geo.h"
#include "polygones.h"
#include "netcdf-proto.h"
#include "functions.h"
#include "grd.h"
#include "map.h"
#include "fe.def"
#include "fe.h"

#include "topo.h"

#define XYZ 0
#define YXZ 1


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
    "  %s -o output input [OPTIONS]\n",prog_name);
  printf("\n"
    "DESCRIPTION\n"
/* USE FIND AND REPLACE TO EDIT THE LINE BELOW, SO THAT THE FILE HEAD IS ALSO UPDATED */
    "  merge together MNTs.\n"
    "\n"
    "OPTIONS :\n"
    "  -b import : merge import MNT in input MNT\n");
  print_OPENMP_help(prog_name); /** \endcode */
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int get_list(const char *filename, vector <char *> & plates)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int nitems, nplates;
  FILE *in;
  vector <size_t> row,col;
  vector <char *> vname;
  char *tmpP,*tmpV;
  
  in=fopen(filename,"r");
  
  if(in==0) return(-1);
    
  do {
    tmpP=new char[1024];
//    tmpV=new char[64];
    nitems=fscanf(in,"%s", tmpP);
    if(nitems!=1) break;
    plates.push_back(tmpP);
//    vname.push_back(tmpV);
    } while (nitems==1);
  
  nplates=plates.size();
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  float  bscale=1.0, scale=1.0,topo_mask=1.e+10,mask=1.e+10;
  double  x,y,t,p;
  double  dx,dy;

  int i,j,k,l,m,n,signus=0,masked_only=0,status;

  size_t size;
  FILE *in,*out;
  char *keyword,*zone;
  int flag;

  char *rootname=NULL,*output=NULL,*input=NULL,*pathname=NULL,*bathymetry=NULL;
  char *bathy=NULL,*format=NULL,*hydro=NULL,*poly=NULL,*list=NULL;
  
  vector <char *> plates;
  int nplates;
  
  grid_t topogrid;
  grid_t grid,zone_grid;

//   short *topo,smask=256*127+255;
  float *ftopo,*base,*buffer,ftopomask;

  float zmin=0,zmax=0;

  plg_t *polygones=NULL;
  int npolygones=0;

  bool debug=false;
  vector<plg_t> r;
  mesh_t mesh;

  fct_echo( argc, argv);

  n=1;
  while (n < argc) {
    keyword=strdup(argv[n]);
    if(strcmp(keyword,"--help")==0) {
      print_help(argv[0]);
      exit(0);
      }
    if(strcmp(keyword,"-inv")==0) {
      scale=-1;
      n++;
      continue;
      }
    if(strcmp(keyword,"-positive")==0) {
      bscale=-1;
      n++;
      continue;
      }
    if(strcmp(keyword,"-debug")==0) {
      debug=true;
      n++;
      continue;
      }
    if(strcmp(keyword,"-add")==0) {
      hydro= strdup(argv[n+1]);
      signus=+1;
      n++;
      n++;
      continue;
      }
    if(strcmp(keyword,"-substract")==0) {
      hydro= strdup(argv[n+1]);
      signus=-1;
      n++;
      n++;
      continue;
      }
    if(strcmp(keyword,"--masked-only")==0) {
      masked_only=1;
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
//     if(strcmp(keyword,"-fix")==0) {
//       sscanf(argv[n+1],"%f",&zmax);
//       n++;
//       n++;
//       continue;
//       }
    switch (keyword[0]) {
      case '-':
        switch (keyword[1]) {
          case 'b' :
            bathymetry= strdup(argv[n+1]);
            n++;
            n++;
            break;

          case 'f' :
            format= strdup(argv[n+1]);
            n++;
            n++;
            break;

          case 'l' :
            list= strdup(argv[n+1]);
            n++;
            n++;
            break;

          case 'p' :
            poly= strdup(argv[n+1]);
            n++;
            n++;
            break;

          case 'o' :
            rootname= strdup(argv[n+1]);
            n++;
            n++;
            break;

          case 'z' :
            zone= strdup(argv[n+1]);
            n++;
            n++;
            break;

          default:
            printf("unknown option %s\n",keyword);
            break;
          }
        break;

      default:
        if(input==NULL) {
          input=strdup(argv[n]);
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

  switch ((unsigned long) input) {
    case 0 :
      if(list!=0) {
        status=get_list(list, plates);
        if(status !=0) {
          printf("cannot get filenames from list=%s\n",list);
          goto error;
          }
        }
      break;
    default :
      plates.push_back(input);
      break;
    }
    
//  status=fe_readmesh(poly,MESH_FILE_FORMAT_TRIGRID,&mesh);

  r=plg_readneigh ((string) poly, true);
  
  status=plg_save("cycles.plg",PLG_FORMAT_SCAN,r);
  
/* *----------------------------------------------------------------------
  load bathymetry grid*/
  printf("#####################################################################\n");
  printf("load working bathymetry: %s\n",bathymetry);
  status=topo_loadfield(bathymetry, &topogrid, &base, &ftopomask, debug);
  if(status !=0) {
    printf("cannot load bathymetry file=%s\n",bathymetry);
    goto error;
    }

/* *------------------------------------------------------------------------------
  apply change of sign on base topography*/
  for (j=0;j<topogrid.ny;j++) {
    for (i=0;i<topogrid.nx;i++) {
      m=j*topogrid.nx+i;
      if(base[m]!=ftopomask) {
        base[m]*=bscale;
        }
      }
    }
  
  ftopo=new float[topogrid.Hsize()];
  
  for(k=17;k<r.size();k++) {
    char filename[1024];
    for (j=0;j<topogrid.ny;j++) {
      for (i=0;i<topogrid.nx;i++) {
        m=j*topogrid.nx+i;
        ftopo[m]=base[m];
        }
      }
/* *----------------------------------------------------------------------
    import new depths*/
    nplates=plates.size();
    if(nplates!=0) {
      for(int p=0;p<nplates;p++) {
        printf("\n#####################################################################\n");
        printf("import extra bathymetry for merging: %s\n", plates[p]);
        status=plg_save("tmp.plg",PLG_FORMAT_SCAN,&(r[k]),1);
        status=topo_import(plates[p], 0, topogrid, ftopo, ftopomask,  zmin,  zmax, "tmp.plg", masked_only, debug);
        if(status !=0) {
          goto error;
          }
        }
      }

/* *----------------------------------------------------------------------
    convert "zero-hydro related" depths to "mean level related" depths*/
    if(hydro!=NULL) {
      status= topo_hydro(hydro, signus, topogrid, ftopo, ftopomask,5);
      if(status !=0) {
        goto error;
        }
      }
    sprintf(filename,"%s-%2.2d",rootname,k);
//     char *dummy=0;
    status=topo_save(filename, 0, format, topogrid, ftopo, ftopomask, debug);
    }
  
  delete[] ftopo;
  delete[] base;
  topogrid.free();

  __OUT_BASE_LINE__("topo_merge sucessfully completed\n");

  exit(0);

 error:
  __OUT_BASE_LINE__("topo_merge aborted\n");
  exit(-1);
}
