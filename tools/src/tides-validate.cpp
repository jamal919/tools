
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
\brief This program creates a table of comparison of waves between model and observations. USED BY CTOH
**/
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

/*
const char *header="\\documentclass[9pt,a4paper]{article} \n\
\\usepackage{here} \n\
\\usepackage{vmargin} \n\
\\usepackage[counterclockwise]{rotating} \n\
\\usepackage{supertabular} \n\
\\begin{document} \n\
\\newpage \n\
\\thispagestyle{empty} \n\
\\tiny \n";
*/
const char *header=
  "\\documentclass[10pt,a4paper]{article}\n"
  "\\setlength{\\textheight}{257mm}\\setlength{\\textwidth}{170mm}\n"
  "\\setlength{\\voffset}{-5.4mm}\\setlength{\\hoffset}{-5.4mm} %The TeX top and left margins are 1in=25.4mm\n"
  "\\setlength{\\topmargin}{0mm}\\setlength{\\headheight}{0mm}\\setlength{\\headsep}{0mm}\n"
  "\\setlength{\\tabcolsep}{3pt} %The LaTeX default value is 6pt\n"
  "\\setlength{\\evensidemargin}{0mm}\n"
  "\\setlength{\\oddsidemargin}{\\evensidemargin}\n"
  "\\pagestyle{empty}\n"
  "\\usepackage{longtable}\n"
  "\\usepackage{textcomp}\n"
  "\\usepackage{amsmath}\n"
  "\\begin{document}\n"
  "\\tiny\\centering\n";
/*
const char *tableTrailer="\\hline \n\
\\end{supertabular} \n";
*/
const char *tableTrailer="\\end{longtable}\n";

/*------------------------------------------------------------------------------

  It reads the sample file from the MOG2D model and filter/detide the time
  series

  input:  sample.# where # is the index of the sampling point
          extract.input where the list position/ name is kept
          (presumedly in /export/bartok4/medsea/extract.input)

  output: XXX.h and XXX.p where XXX is the sampling point name

-----------------------------------------------------------------------*/


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int extract_SGatlas(char *filename, const char **varnames, vector<mgr_t> mgr, int nmgr,double *a,double *G, double mask, int strict, int verbose)
//similar to tide_atlasSG2positions
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
///extract amplitudes and phases of structured atlas at tide gauges locations
{
  int k, l, n, nerrors, status;
  int kk,ll,mm;
  int64_t accel;
  const int nnext=24;
  int nextk[nnext]={-1,+1, 0, 0,-1,+1,-1,+1,-2,+2, 0, 0,-2,+2,-2,+2,-1,-1,+1,+1,-2,-2,+2,+2};
  int nextl[nnext]={ 0, 0,-1,+1,-1,+1,+1,-1, 0, 0,-2,+2,-1,-1,+1,+1,-2,+2,-2,+2,-2,+2,-2,+2};
  fcomplex zz,cmask;
  fcomplex *tide=NULL;
  double x, y;
  grid_t grid;

  status=poc_get_grid(filename, varnames[0], &grid);
  if(status)NC_TRAP_ERROR(return,status,1,"poc_get_grid(\"%s\",\"%s\",) error",filename,varnames[0]);

/*------------------------------------------------------------------------------
  load netcdf variable */
  exitIfNull(tide   =new fcomplex[grid.nx*grid.ny]);
  cmask=NC_FILL_FLOAT;
  status=poc_get_cvara(filename,varnames[0],varnames[1],-1,tide);
  if(status)NC_TRAP_ERROR(return,status,1,"poc_get_cvara(\"%s\",\"%s\",\"%s\",) error",filename,varnames[0],varnames[1]);
  
  set_grid_list(&grid);

/*------------------------------------------------------------------------------
  interpolate */
  nerrors=0;
  for(n = 0; n < nmgr; n++) {
    if(strcmp(mgr[n].loc.units,"degrees")==0) {
      x = mgr[n].loc.lon;
      y = mgr[n].loc.lat;
      }
    else {
      x = mgr[n].loc.lon * r2d;
      y = mgr[n].loc.lat * r2d;
      }
//      y = gLGP1data[0][node].lat * r2d;
    if(x < grid.xmin - grid.dx / 2.)
      x = x + 360.0;
    if(x > grid.xmax + grid.dx / 2.)
      x = x - 360.0;
//     status = map_interpolation(grid, tide, cmask, x, y,&zz);
    index_interpolation(grid, x, y, &accel, tide, cmask, &zz, 0);
    if(zz==cmask) status=-1;
    else status=0;
    if (strict==0)
    if((status != 0) || (zz==cmask) ) {
      status=map_index( grid,  x,  y, &k, &l);
      ///\todo 2011-10-26 Damien Allain : ask Florent to explain this difference with tide_atlasSG2positions
      if(status == 0) {
        mm=0;
        do {
          kk=k+nextk[mm];
          ll=l+nextl[mm];
          mm++;
          status=-1;
          if(kk<0) continue;
          if(ll<0) continue;
          if(kk>grid.nx-1) continue;
          if(ll>grid.ny-1) continue;
          x=map_grid_x(grid,kk,ll);
          y=map_grid_y(grid,kk,ll);
//           status = map_interpolation(grid, tide, cmask, x, y,&zz);
          index_interpolation(grid, x, y, &accel, tide, cmask, &zz, 0);
          if(zz==cmask) status=-1;
          else status=0;
          } while((status!=0)&&(mm<24));
        mm=0;
        while((status!=0)&&(mm<24)) {
          kk=k+2*nextk[mm];
          ll=l+2*nextl[mm];
          mm++;
          status=-1;
          if(kk<0) continue;
          if(ll<0) continue;
          if(kk>grid.nx-1) continue;
          if(ll>grid.ny-1) continue;
          x=map_grid_x(grid,kk,ll);
          y=map_grid_y(grid,kk,ll);
//           status = map_interpolation(grid, tide, cmask, x, y,&zz);
          index_interpolation(grid, x, y, &accel, tide, cmask, &zz, 0);
          if(zz==cmask) status=-1;
          else status=0;
          }
        }
      }
    if(status != 0) {
      if(verbose) printf("interpolation error: node %d lon=%lf lat=%lf \n", n, x, y);
      a[n] =  mask;
      G[n] =  mask;
      nerrors++;
      }
    else {
      a[n] =  abs(zz);
/*------------------------------------------------------------------------------
      convert into phase lag*/
      G[n] = -arg(zz)* r2d;
      }
    }
  printf("%d points = %d interpolation errors + %d ok\n", nmgr, nerrors, nmgr-nerrors);
  
  if(nerrors and verbose==1)
    map_printgrid(grid);

  grid.free();
  delete[]tide;

  return (0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int *mgr_detect(mesh_t & mesh, vector<mgr_t> mgr, int nmgr)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int n;
  double *x, *y;
  int *elements;
  
  x=new double[nmgr];
  y=new double[nmgr];
  
  for(n = 0; n < nmgr; n++) {
    if(strcmp(mgr[n].loc.units,"degrees")==0) {
      x[n] = mgr[n].loc.lon;
      y[n] = mgr[n].loc.lat;
      }
    else {
      x[n] = mgr[n].loc.lon * r2d;
      y[n] = mgr[n].loc.lat * r2d;
      }
    }
  elements=fe_detect( mesh, x, y, nmgr);
  
  delete[] x;
  delete[] y;
  
  return(elements);

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int extract_UGatlas(const char *filename,const char **varnames, vector<mgr_t> mgr, int nmgr,double *a,double *G, double mask,const char *disc_name, int iteration, double dmax, int level, int strict, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/**
\param dmax maximum distance in m
*/
/*----------------------------------------------------------------------------*/
{
  int n,frame, status;
  fcomplex zz,cmask;
  fcomplex *tide=NULL;
  double x, y;
  char amplitude[256],phase[256];
  mesh_t mesh;
  discretisation_t descriptor;
  int *elements;
  frame_t mesh_frame;

  cdfgbl_t global;
  int id;
  variable_t varinfo;
  
  frame=iteration;

  if(verbose==1) printf("load harmonic file : %s (variables %s %s, iteration %d, level %d)\n",filename, varnames[0],varnames[1], frame, level);
  status= cdf_globalinfo(filename,&global,0);
  if(status !=0) {
    printf("cannot open %s\n",filename);
    goto error;
    }

  id=discretisation_from_name(disc_name);
  if(id ==-1) {
    goto error;
    }
    
  status=fe_readgeometry(filename, &mesh);
  if(status !=0) {
    goto error;
    }
    
  status=fe_readdiscretisation(filename, &mesh, 0, id);
  
  if(status!=0) {
    if(verbose==1) printf("connectivity not found in file %s for discretisation %s\n",filename,disc_name);
    status=discretisation_init(&mesh, id, 0);
    }
  
  descriptor=get_descriptor(mesh,id);
  
  exitIfNull(
    tide=new complex<float>[descriptor.nnodes]
    );
/*------------------------------------------------------------------------------
  load netcdf variable */
  frame=iteration;
  
  if(varnames[0]==0) {
    sprintf(amplitude,"a_eta_%s",disc_name);
    }
  else {
    sprintf(amplitude,"%s",varnames[0]);
    }
  if(varnames[1]==0) {
    sprintf(phase,"G_eta_%s",disc_name);
    }
  else {
    sprintf(phase,"%s",varnames[1]);
    }
  
  if(level==-1) {
    status=poc_get_UG3D(filename, frame, amplitude, phase, tide);
    }
  else  {
    status=poc_get_UG3D(filename, frame, level, amplitude, phase, tide);
    }
    
  if(status !=0) {
    printf("cannot load frame %d in %s\n",frame,filename);
    goto error;
    }
  cmask=fcomplex(9999.,9999.);

//  nprocs=initialize_OPENMP(-1, 0);

//#pragma omp parallel for private(x,y,zz,status) if(nprocs>1)

  elements=mgr_detect(mesh, mgr, nmgr);
        
  status=fe_minmax(mesh,  mesh_frame);
  
  for(n = 0; n < nmgr; n++) {
    if(strcmp(mgr[n].loc.units,"degrees")==0) {
      x = mgr[n].loc.lon;
      y = mgr[n].loc.lat;
      }
    else {
      x = mgr[n].loc.lon * r2d;
      y = mgr[n].loc.lat * r2d;
      }
    int m=elements[n];
    a[n] =  mask;
    G[n] =  mask;
    if(m == -1) {
      if (strict==0) {
        bool inside=mesh_frame.inside(x,y);
        if(not inside) continue;
/*------------------------------------------------------------------------------
        try to use nearest node in mesh */
        double d=fe_extrapolate(descriptor,tide,x,y,dmax,&zz,verbose);
        
        if(d<dmax) {
          a[n] =  abs(zz);
/*------------------------------------------------------------------------------
          convert into phase lag*/
          G[n] = -arg(zz)* r2d;
          }
        else {
          if(verbose==1) {
            printf("interpolation issue: station %4d %s, lon=%9.3lf lat=%9.3lf out of mesh\n", n, mgr[n].name,  x, y);
            }
          a[n] =  mask;
          G[n] =  mask;
          }
        }
      }
    else {
      status =fe_interpolate2D(mesh, id, tide, cmask, x, y, m, &zz);
      if(status != 0) {
        printf("interpolation error: node %d lon=%lf lat=%lf \n", n, x, y);
        a[n] =  mask;
        G[n] =  mask;
        }
      else {
        a[n] =  abs(zz);
/*------------------------------------------------------------------------------
        convert into phase lag*/
        G[n] = -arg(zz)* r2d;
        }
      }
    }

  delete[] tide;
  return (status);

error:
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int extract_UGatlas_ASCII(const char *filename, const char *meshfile, vector<mgr_t> mgr, int nmgr,double *a,double *G, double mask,const char *disc_name, int iteration, double dmax, int level, int strict, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/**
\param dmax maximum distance in m
*/
/*----------------------------------------------------------------------------*/
{
  int n,frame, status;
  fcomplex zz,cmask;
  fcomplex *tide=NULL;
  double x, y;
  mesh_t mesh;
  discretisation_t descriptor;
  int *elements;

  cdfgbl_t global;
  int id;
  variable_t varinfo;
  
  frame=iteration;

/* *-----------------------------------------------------------------------------
  Read FE mesh */
  if(meshfile != NULL) {
    status=fe_readmesh(meshfile,MESH_FILE_FORMAT_TRIGRID,&mesh);
    if(status !=0) goto error;
    status=fe_list(&mesh);
    if(status !=0) goto error;
    status=fe_edgetable(&mesh,0,0);
    if(status !=0) goto error;
    }
  else {
    printf("no mesh file specified; abort...\n");
    status=ENOENT;
    goto error;
    }
    
  id=discretisation_from_name(disc_name);

  status=discretisation_init(&mesh, id, 0);
  descriptor=get_descriptor(mesh, id);
  
  exitIfNull(
    tide=new complex<float>[descriptor.nnodes]
    );
/*------------------------------------------------------------------------------
  load netcdf variable */
  frame=iteration;
  
  status=quoddy_loadc1(filename, descriptor.nnodes, tide);
    
  if(status !=0) {
    printf("cannot load file %s\n",filename);
    goto error;
    }
  cmask=fcomplex(9999.,9999.);

  elements=mgr_detect(mesh, mgr, nmgr);
  
  for(n = 0; n < nmgr; n++) {
    if(strcmp(mgr[n].loc.units,"degrees")==0) {
      x = mgr[n].loc.lon;
      y = mgr[n].loc.lat;
      }
    else {
      x = mgr[n].loc.lon * r2d;
      y = mgr[n].loc.lat * r2d;
      }
    int m=elements[n];
    a[n] =  mask;
    G[n] =  mask;
    if(m == -1) {
      if (strict==0) {
/*------------------------------------------------------------------------------
        try to use nearest node in mesh */
        double d=fe_extrapolate(descriptor,tide,x,y,dmax*1e3,&zz,verbose);
        
        if(d<dmax) {
          a[n] =  abs(zz);
/*------------------------------------------------------------------------------
          convert into phase lag*/
          G[n] = -arg(zz)* r2d;
          }
        else {
          if(verbose==1) {
            printf("interpolation issue: station %4d %s, lon=%9.3lf lat=%9.3lf out of mesh\n", n, mgr[n].name,  x, y);
            }
          a[n] =  mask;
          G[n] =  mask;
          }
        }
      }
    else {
      status =fe_interpolate2D(mesh, id, tide, cmask, x, y, m, &zz);
      if(status != 0) {
        printf("interpolation error: node %d lon=%lf lat=%lf \n", n, x, y);
        a[n] =  mask;
        G[n] =  mask;
        }
      else {
        a[n] =  abs(zz);
/*------------------------------------------------------------------------------
        convert into phase lag*/
        G[n] = -arg(zz)* r2d;
        }
      }
    }

  delete[] tide;
  return (status);

error:
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void print_help(char *prog_name)

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
    " NAME AND VERSION\n"
    "  %s version " VERSION " " REVISION "\n", prog_name);
  printf("\n"
    " USE\n"
    "   %s -g file.mgr -a WAVE-file.nc wave1 [ wave2 ...]\n",prog_name);
  printf("\n"
    " DESCRIPTION\n"
/* USE FIND AND REPLACE TO EDIT THE LINE BELOW, SO THAT THE FILE HEAD IS ALSO UPDATED */
    "   This program creates a table of comparison of waves between model and observations. USED BY CTOH.\n"
    "\n"
    "OPTIONS :\n"
    "  -h,--help : show this help and exit\n"
    "  -b : bathymetry database\n"
    "  -d : print data/models differences \n"
    "  -f : followed by analysis file format\n"
    "  -l : output format in LaTeX\n"
    "  -a : naming convention for tidal atlases, see below\n"
    "  -m : ASCII mesh file. (I do not know how to use this...)\n"
    "  -p : path for tidal atlases\n"
    "  -o : output root name. Default output: validate.out or validate.tex\n"
    "  -g : harmonic data input file\n"
    "  -v : followed by 2 variables' name conventions. Default for structured atlases: Ha Hg\n"
    "  -i : iteration number if any (-1 for last one)\n"
    "  -x : output order: FILE, ALPHA, LAT, LON, IMMERSION\n"
    "  -unstructured [discretisation] : if the netcdf file is unstructured: please provide the discretisation type (usually LGP1 or LGP2)\n"
    "  -strict : only compare with data strictly included in grid (no extrapolation)\n"
    "  --frame [lonmin:lonmax;latmin:latmax] : limits validation to stations included in frame\n"
    "  --polygons : followed by path of polygons containing stations validation is limited to\n"
    "  -level : \n"
    "  --regions : \n"
    "  -cm : analysis file unit is cm\n"
   );
  
  print_tide_decode_atlasname_help();
  /** \endcode */

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double rounded(double factor, double value)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double tmp=NINT(10.*factor*value)/10.;
  return(tmp);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void mgr_print(spectrum_t spectrum, const char *rootname, vector<mgr_t> mgr, int nmgr, int *list, grid_t topogrid, float *topo, float topomask,
                  double  **a,double  **G,const char *atlas_directory,const char *atlas_convention, char latex, char show_diff, bool do_legend, double mean[10],double rms[10])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double mask=999.9;
  complex<double> z;
//   double  mean[10],rms[10],count;
  double  count;
  double  a1,p1,a2,p2,d,da,dG;
  int i,k,l,m,n,status;
  FILE *out=NULL,*stream[2];
  int s;
  char report[1024];
  fcomplex meanc;
  date_t date;

  const char *tab=NULL,*notab=NULL,*newline=NULL,*percent=NULL;

  float h;
  double x,y;
  const float factor=10.f,hlim=15e3f;
  
  if(rootname==NULL) {
    rootname="validate";
    }
  if(latex)
    sprintf(report,"%s.tex",rootname);
  else
    sprintf(report,"%s.out",rootname);
  
//   STDOUT_BASE_LINE("writing ");
//   if(latex)printf("LaTeX");else printf("ascii");
//   printf(" report to %s\n",report);

  if(latex) {
    tab  =" & ";
    notab="   ";
    newline="\\\\";
    percent="\\%";
    }
  else {
    tab  =" ";
    notab=" ";
    newline="";
    percent="%";
    }
  
  out=fopen(report,"w");
  if(out==0) TRAP_ERR_EXIT(errno,"Failed to open %s for writing (%d %s)\n", report, errno, strerror(errno));
//   printf("#################################################################\n");
//   printf("output report to %s\n",report);

  stream[0]=stdout;
  stream[1]=out;

  if(latex) {
    fprintf(out,"%s",header);
    }
  fprintf(out,"model    : %s/%s%s\n",atlas_directory,atlas_convention,newline);
  #if 0
  ///\todo 2011-10-19 Damien Allain : show bathymetry database ?
  fprintf(out,"database : %s%s\n",clean_string(datafile),newline);
  #endif
  fprintf(out,"\n");

  if(latex) {
    fprintf(out,"\\begin{longtable}{|c|c|cc|");
    for (k=0;k<spectrum.n;k++) {
      fprintf(out,"cc");
      if(show_diff) {
        fprintf(out,"ccc");
        }
      fprintf(out,"|");
      }
    fprintf(out,"}\n\\hline\n");
    }
//  fprintf(out,"%s\n",newline);

  const char *NoStationLonLatFormat="%4s%s%32s%s%14s%s%14s";
  fprintf(out,NoStationLonLatFormat,"No",tab,"Station",tab,"longitude ",tab,"latitude  ");
  for (k=0;k<spectrum.n;k++) {
    fprintf(out,"%s",tab);
    if(latex) fprintf(out,"\\multicolumn{%d}{c|} {",show_diff?5:2);
    fprintf(out,"%9s%5s%*s",spectrum.waves[k].name,"(o/m)",show_diff?16:1,notab);
    if(latex) fprintf(out,"}");
    }
  fprintf(out,"%s\n",newline);

  fprintf(out,NoStationLonLatFormat,"",tab,"",tab,"",tab,"");
  for (k=0;k<spectrum.n;k++) {
    fprintf(out,"%s%7s%s%7s",tab,"A(cm) ",tab,"G(deg)");
    if(show_diff) {
      fprintf(out,"%s%3s%s%3s%s%6s",tab,"dA",tab,"dG",tab,"dE");
      }
    }
  fprintf(out,"%s\n",newline);
  if(latex)
    fprintf(out,"\\hline\\endhead\n"
      "\\hline\\endfoot\n\n");

  ///For all tides gauges
  for (m=0;m<nmgr;m++) {
    n=list[m];
    ///If it is in the domain
    if(a[0][n]==mask)continue;
    ///It prints :
    char *Number;
    asprintf(&Number,"%d",n+1);
    string slon,slat;
    slon=dms(mgr[n].loc.lon);
//     mgr[n].loc.lat=-0.33;
    slat=dms(mgr[n].loc.lat);
    if(latex){
      replace(&slon,"°","\\textdegree");
      replace(&slat,"°","\\textdegree");
      }
    ///- its number, name and position
    fprintf(out,NoStationLonLatFormat,Number,tab,check_string(mgr[n].name),tab,slon.c_str(),tab,slat.c_str());
    free(Number);
    ///\todo 2011-10-18 Damien Allain : finish doxygenation
    if(topo!=0) {
      x=mgr[n].loc.lon;
      y=mgr[n].loc.lat;
      x=map_recale(topogrid,x);
      status=map_interpolation(topogrid, topo,topomask,x,y,&h);
      }
    for (k=0;k<spectrum.n;k++) {
      for(i=0;i<mgr[n].nwave;i++) {
        //getting the indices of wave in mgr
         if(strcmp(spectrum.waves[k].name,mgr[n].data[i].constituent.name)==0) {
          break;
          }
        }
      if(i==mgr[n].nwave) {
        fprintf(out,"%s---/---%s---/---",tab,tab);
        if(show_diff) {
          fprintf(out,"%s---%s---%s------",tab,tab,tab);
          }
        }
      else {
        a1=100*mgr[n].data[i].amp;
        p1=mgr[n].data[i].phi;
        if(topo != (size_t) 0) {
/* *----------------------------------------------------------------------------

          Profondeur limite

-----------------------------------------------------------------------------*/
          if(h>hlim) {
            a[k][n]=mask;
            }
          };
        if(a[k][n]==mask) {
          fprintf(out,"%s%3.0f/---%s%3.0f/---",tab,a1,tab,p1);
          if(show_diff) {
            fprintf(out,"%s---%s---%s------",tab,tab,tab);
            }
          }
        else {
          a2=100*a[k][n];
          p2=G[k][n];
          if(p2<  0) p2+=360.0;
          if(p2>360) p2-=360.0;
          fprintf(out,"%s%3.0f/%3.0f%s%3.0f/%3.0f",tab,a1,a2,tab,p1,p2);
          if(show_diff) {
            z=polar(a2,-p2*d2r)-polar(a1,-p1*d2r);
            d=abs(z);
            da=a2-a1;
            dG=p2-p1;
            if(dG>+180.) dG-=360.0;
            if(dG<-180.) dG+=360.0;
            fprintf(out,"%s%6.2g%s%3.0f%s%6.1f", tab,da,tab,dG,tab,d);
            }
          }
        }
      }
    fprintf(out,"%s\n",newline);
    }
  if(latex) fprintf(out,"\\end{longtable}\n");

  fprintf(out,"\n");
//  fprintf(out,"%s\n",newline);
  
  if(do_legend)
  for(s=0;s<2;s++){
    if(latex) {
      fprintf(stream[s],"\\begin{longtable}{|c|cc|cc|ccc|c|c|} \n");
      fprintf(stream[s],"\\hline\n");
      }
    fprintf(stream[s],"%10s%2s","wave",tab);
    if(latex) {
      fprintf(stream[s],"\\multicolumn{2}{c|} {$\\Delta a$(mm)} & \\multicolumn{2}{c|} {$\\Delta G$(deg)} & \\multicolumn{3}{c|} {$e$(mm)} & e/a & N \\\\ \n");
      }
    else {
      fprintf(stream[s],"%5s%s%5s%2s","A",tab,"(mm)",tab);
      fprintf(stream[s],"%5s%s%5s%2s","G",tab,"(deg)",tab);
      fprintf(stream[s],"%5s%s%5s%s%5s%2s%5s%2s%s%s\n","E",tab,"(mm)",tab,"",tab,"E/A",tab,"N (model/obs)",newline);
      }
    fprintf(stream[s],"%10s","");
    for(k=0;k<3;k++)
      fprintf(stream[s],"%2s%5s%s%5s",tab,"mean",tab,"rms");
    fprintf(stream[s],"%s%5s%s%s%s\n",tab,"rms*",tab,tab,newline);
    if(latex)
      fprintf(stream[s],"\\hline\\endhead\n"
        "\\hline\\endfoot\n\n");
    }
  
/*------------------------------------------------------------------------------
  print statistics */
  string equations;
  
  for (k=0;k<spectrum.n;k++) {
    for(l=0;l<10;l++) {
      mean[l]=0;
      rms[l]=0;
      }
    count=0;
    for (m=0;m<nmgr;m++) {
      n=list[m];
      for(i=0;i<mgr[n].nwave;i++) {
        if(strcmp(spectrum.waves[k].name,mgr[n].data[i].constituent.name)==0) {
          break;
          }
        }
      if(i==mgr[n].nwave) continue;
      a1=100*mgr[n].data[i].amp;
      p1=mgr[n].data[i].phi;
      if(a[k][n]==mask) continue;
      a2=100*a[k][n];
      p2=degree_recale(G[k][n],0.);
      
      equations="With:\n"
        "\\begin{gather*}\n";
      
      z=polar(a2,-p2*d2r)-polar(a1,-p1*d2r);
      //"\\\sigma = \\sqrt{N{-1} \\sum_{i=1}^N (x_i - \\\mu)^2}, {\\rm \\ \\ where\\ \\ } \\\mu = N{-1} \\sum_{i=1}^N x_i"
      asprintf(equations,"z=z_m-z_o %s\n",newline);
      
/*------------------------------------------------------------------------------
      complex error deviation from 0 (RMS)*/
      l=0;
      d=abs(z)*M_SQRT1_2;
      rms[l]+=d*d;
      //asprintf(equations,"d_%d=\\sqrt{0.5}\\left|z\\right| %s\n",l,newline);/* WRONG because mean[0] is kept at 0. */
      
/*------------------------------------------------------------------------------
      amplitude mean misfit/rms*/
      l++;
      d=a1-a2;
      mean[l]+=d;
      rms[l]+=d*d;
      asprintf(equations,"d_%d=\\left|z_m\\right|-\\left|z_o\\right| %s\n",l,newline);
      
/*------------------------------------------------------------------------------
      phase lag mean misfit/rms*/
      l++;
      d=degree_recale(p1-p2,0.);
      mean[l]+=d;
      rms[l]+=d*d;
      asprintf(equations,"d_%d=\\arg\\left(z_m\\right)-\\arg\\left(z_o\\right) %s\n",l,newline);
      
/*------------------------------------------------------------------------------
      complex error mean misfit (module) / standard deviation*/
      l++;
      d=real(z);
      mean[l]+=d;
      rms[l]+=d*d;
      asprintf(equations,"d_%d=\\Re\\left(z\\right) %s\n",l,newline);
      
      l++;
      d=imag(z);
      mean[l]+=d;
      rms[l]+=d*d;
      asprintf(equations,"d_%d=\\Im\\left(z\\right) %s\n",l,newline);
      
/*------------------------------------------------------------------------------
      complex error deviation from 0 (RMS), percent of amplitude*/
      l++;
      d=abs(z)*M_SQRT1_2/(0.5*(a1+a2));
      mean[l]+=d;
      asprintf(equations,"d_%d=\\frac{\\sqrt{0.5}\\left|z\\right|}{0.5\\left(\\left|z_m\\right|+\\left|z_o\\right|\\right)} %s\n",l,newline);
      
      count++;
      }
    
    for(l=0;l<6;l++) {
      mean[l]/=count;
      rms[l]=sqrt(rms[l]/count - square(mean[l]) );
      }
    asprintf(equations,"\\mu_n=N^{-1}\\sum_N d_n %s\n",newline);
    asprintf(equations,"\\sigma_n=\\sqrt{N^{-1}\\sum \\left(d_n-\\mu_n\\right)^2} %s\n",newline);
    equations+="\\end{gather*}\n";
    
    static const char *meanRMSFormat="%2s%5.1f%s%5.1f";
    for(s=0;s<2;s++) {
      fprintf(stream[s],"%10s",spectrum.waves[k].name);
/*------------------------------------------------------------------------------
      amplitude mean misfit/rms*/
      fprintf(stream[s],meanRMSFormat,tab,rounded(factor,mean[1]),tab,rounded(factor,rms[1]));
      if(s==0)asprintf(equations,"The values in the $\\Delta a$ column are $\\mu_1$ and $\\sigma_1$.%s\n",newline);
/*------------------------------------------------------------------------------
      phase lag mean misfit/rms*/
      fprintf(stream[s],meanRMSFormat,tab,rounded(1.,mean[2]),tab,rounded(1.,rms[2]));
      if(s==0)asprintf(equations,"The values in the $\\Delta G$ column are $\\mu_2$ and $\\sigma_2$.%s\n",newline);
/*------------------------------------------------------------------------------
      complex error mean misfit (module) / standard deviation*/
      fprintf(stream[s],meanRMSFormat,tab,rounded(factor,hypot(mean[3],mean[4])/M_SQRT2),tab,rounded(factor,hypot(rms[3],rms[4])/M_SQRT2));
      if(s==0)asprintf(equations,"The values in the $e$ column are $c_1=\\sqrt{0.5\\left(\\mu_3^2+\\mu_4^2\\right)}$, $c_2=\\sqrt{0.5\\left(\\sigma_3^2+\\sigma_4^2\\right)}$");
/*------------------------------------------------------------------------------
      complex error deviation from 0 (RMS)*/
      fprintf(stream[s],"%s%5.1f%2s",tab,rounded(factor,rms[0]),tab);
      //if(s==0)asprintf(equations,"\\sigma_0,");
      if(s==0)asprintf(equations," and $c_3=\\sqrt{N^{-1}\\sum 0.5\\left|z\\right|^2}$. Note that $c_1^2+c_2^2=c_3^2$.%s\n",newline);
/*------------------------------------------------------------------------------
      complex error deviation from 0 (RMS), percent of amplitude*/
      fprintf(stream[s],"%5.1f%s%2s",rounded(100.,mean[5]),percent,tab);
      if(s==0)asprintf(equations,"The value in the e/a column is $\\mu_5$.\n");
/*------------------------------------------------------------------------------
      nb valid comparisons / nb stations*/
//       fprintf(stream[s],"%3d/%3d%s\n",(int) count,nmgr, newline);
      
      fprintf(stream[s],"%3d/%3d%2s",(int) count,nmgr,tab);
      fprintf(stream[s],"real / imaginary : %5.1f%s %5.1f%s\n",rounded(factor,mean[3]),tab,rounded(factor,mean[4]), newline);

      }
    
    }

  fprintf(out,"\n");
  if(latex) {
    fprintf(out,"%s",tableTrailer);
    fprintf(out,"%s%s%s","\n",equations.c_str(),"\n"
      "\\end{document}\n");
    }
  
  fclose(out);
  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void mgr_print2(spectrum_t spectrum,const char *rootname, vector<mgr_t> mgr, int nmgr, int *list, grid_t topogrid, float *topo, float topomask,
                  double  **a,double  **G,const char *atlas_directory,const char *atlas_convention, char latex, char show_diff)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  print data/model misfit (complex and polar)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  double  mask=999.9;
  double  zi,zr;
  double  a1,p1,a2,p2,d,da,dG;
  int i,k,m,n;
  FILE *out=NULL;
  char report[1024];

  fcomplex cmisfit;

  sprintf(report,"%s.diags",rootname);
  
  for (k=0;k<spectrum.n;k++) {
    sprintf(report,"%s.%s.diags",rootname,spectrum.waves[k].name);
    out=fopen(report,"w");
    for (m=0;m<nmgr;m++) {
      n=list[m];
      if(a[0][n]==mask)continue;
      i=mgr[n].wave_index(spectrum.waves[k].name);
      if(i==-1) continue;
      a1=100*mgr[n].data[i].amp;
      p1=mgr[n].data[i].phi;
      a2=100*a[k][n];
      p2=G[k][n];
      if(p2<  0) p2+=360.0;
      if(p2>360) p2-=360.0;
      zr=a2*cos(-p2*d2r)-a1*cos(-p1*d2r);
      zi=a2*sin(-p2*d2r)-a1*sin(-p1*d2r);
      d=sqrt(zr*zr+zi*zi);
      da=a2-a1;
      dG=p2-p1;
      if(dG>+180.) dG-=360.0;
      if(dG<-180.) dG+=360.0;
      fprintf(out,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", mgr[n].loc.lon,mgr[n].loc.lat, zr, zi, da, dG, a1, a2, p1, p2 );
      }
    fclose(out);
    }

  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int mgr_SetRegions (const char *filename, vector<mgr_t> mgr, int nmgr, unsigned short *region)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int n,status;
  double t,p;
  grid_t grid;
  float *code, mask, z;
    
  pair<map<const char *,int>::iterator,bool> toponym_status;
  toponyms_t *toponyms;
  
  std::map<int, double> wdc;
  pair<map<int, double>::iterator,bool> wdc_status;
  
  toponyms=new toponyms_t;
  status=geo_init_regions(toponyms);
  
  status=grd_loadgrid (filename,&grid);
  if(status!=0) return(-1);
  
  code=new float[grid.nx*grid.ny];
  status=grd_loadr1(filename,  grid, code, &mask);
  if(status!=0) return(-1);
    
  for(n=0;n<nmgr;n++) {
    t = mgr[n].loc.lon;
    p = mgr[n].loc.lat;
    t=map_recale(grid,t);
    status=map_nearestvalue00(grid, grid.nx,code,mask,t,p,&z);
    if(z==mask) {
      printf("masked toponyms at t=%lf p=%lf\n",t,p);
      region[n]=(unsigned short) 1100;
      }
    else {
      region[n]=(unsigned short) z;
      }
    }
    
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int mgr_validate(const char **varnames, char **wave,const char *output,const char *atlas_directory,const char *atlas_convention, double scale,const char *meshfile,
		   const char *mgrfile, int format,const char *ordering, const char *bathymetry,const char *regionsfile, vector<plg_t> & limits, char latex, 
		   char show_diff, int structured,const char *discretisation, int iteration, int level, int strict, string tag, string year, bool silent)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double mask=999.9;
  int count,id;
  double  **a=NULL,**G=NULL;
  int i,k,status;
  const double dmax=10e3;/*< maximum distance in m */
  double  mean[10],rms[10];
  double  summary_mean[20],summary_rms[20];
  FILE *dump;
  spectrum_t spectrum;
  date_t date;
  vector<mgr_t> mgr;
  string filename;

  int nmgr;
  char *model;

  int *list=NULL;
  
  unsigned short *region;
  toponyms_t *toponyms;
  
  grid_t topogrid;
  float *topo=NULL,topomask;

  date.year=1950;
  date.month=1;
  date.day=1;
  date.second=0.0;
  astro_angles_t astro_angles;
  spectrum.init(initialize_tide(&astro_angles,date),wave);

  for (i=0; i<spectrum.n; i++) {
    spectrum.waves[i].init();
    if(not silent) printf ("wave: %10s, pulsation: %12.6f degrees/h \n", spectrum.waves[i].name,spectrum.waves[i].omega);
    }

//   for (i=0; i<spectrum.n; i++) {
//     for (j=i+1; j<spectrum.n; j++) {
//       tau=deltaOmega2separation(spectrum.waves[j].omega-spectrum.waves[i].omega);
//       printf ("wave: %10s %10s, separation: %9.3f days \n",
//         spectrum.waves[i].name,spectrum.waves[j].name,tau);
//       }
//     }

/* *-----------------------------------------------------------------------------
  read tide gauge database */
  printf("#################################################################\n");
  printf("load observation file : %s\n",mgrfile);
  nmgr=mgr_load(mgrfile, mgr, format);

  if(limits.size()!=0) {
    status=mgr_exclude(mgr,limits,PLG_POINT_EXTERIOR);
    nmgr=mgr.size();
    }

/** to nest in a separate routine -------------------------------------------*/
  vector<mgr_t> mgr_model;
  mgr_create(nmgr, mgr_model, spectrum);
  
  for(size_t m = 0; m < nmgr; ++m){
    mgr_model[m].loc = mgr[m].loc;
    strcpy(mgr_model[m].name, mgr[m].name);
    strcpy( mgr_model[m].validation, mgr[m].validation);
    mgr_model[m].number = mgr[m].number;
    mgr_model[m].track = mgr[m].track;
    mgr_model[m].duree = mgr[m].duree;
    }
    
/** to nest in a separate routine -------------------------------------------*/
  if(nmgr==-1) {
    printf("error in reading %s\n",mgrfile);
    goto error;
    }

  if(ordering==NULL) {
    ordering=strdup("ALPHA");
    }

  list= mgr_order(mgr, nmgr,ordering,1);

  exitIfNull(
    a=new double*[spectrum.n]
    );
  
  exitIfNull(
    G=new double*[spectrum.n]
    );

  for (k=0;k<spectrum.n;k++) {
    exitIfNull(
      a[k]=new double[nmgr]
      );
  
    exitIfNull(
      G[k]=new double[nmgr]
      );
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  read tidal atlas database and interpolate

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  char *decodedNames[2];
  printf("#################################################################\n");
  printf("load solution files\n");
  for (k=0;k<spectrum.n;k++) {
    status=tide_decode_atlasname(atlas_directory,atlas_convention,wave[k], 0,&model);
    decodedNames[0]=decode_atlasname(varnames[0],wave[k],0);
    decodedNames[1]=decode_atlasname(varnames[1],wave[k],0);
    if(not silent) printf("treating %s wave from %s(%s,%s)  ",wave[k],model,varnames[0],varnames[1]);
    if(structured==1) {
      status=extract_SGatlas(model,(const char**)decodedNames,mgr,nmgr,a[k],G[k],mask, strict,0);
      }
    else {
      if(not silent) printf("\n");
      if(meshfile==0)
        status=extract_UGatlas(model,(const char**)decodedNames,mgr,nmgr,a[k],G[k],mask,discretisation, iteration, dmax, level, strict, (silent==false));
      else
        status=extract_UGatlas_ASCII(model,meshfile,mgr,nmgr,a[k],G[k],mask,discretisation, iteration, dmax, level, strict, 0);
      }
    
    delete[] model;
    delete[] decodedNames[0];
    delete[] decodedNames[1];
    }
  if(status!=0) goto error;
  
  for (i=0;i<nmgr;i++){
    for (k=0;k<spectrum.n;k++) {
      if(a[k][i]!=mask) a[k][i]*=scale;
      }
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  read topo database

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if(bathymetry!=NULL) {
    status=grd_loadgrid(bathymetry,&topogrid);
    if(status !=0) {
      STDOUT_BASE_LINE("cannot load bathymetry file=%s\n",bathymetry);
      exit(-1);
      }
    exitIfNull(
      topo=(float *) malloc(topogrid.nx*topogrid.ny*sizeof(float))
      );
    topogrid.modeH=0;
    status=  grd_loadr1(bathymetry,topogrid,topo,&topomask);
    }
  else {
    topo=0;
    }
  
  printf("\n");
  mgr_print(spectrum, output, mgr, nmgr, list, topogrid, topo, topomask, a,  G, atlas_directory, atlas_convention, latex, show_diff,true, mean, rms);
  id=0;
  summary_mean[id]=mean[0];
  summary_rms[id]=rms[0];

/** to nest in a separate routine -------------------------------------------*/
  for (i=0;i<nmgr;i++){
    for (k=0;k<spectrum.n;k++) {
      mgr_model[i].data[k].amp=100*a[k][i];
      if (G[k][i]<0.0) mgr_model[i].data[k].phi=G[k][i]+360.0;
      else if (G[k][i]>360.0) mgr_model[i].data[k].phi=G[k][i]-360.0;
      else mgr_model[i].data[k].phi=G[k][i];
      }
    }
/** to nest in a separate routine -------------------------------------------*/
    

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  process regional statistics
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  region=new unsigned short[nmgr];
  
  if(regionsfile==NULL) goto end;
  printf("#################################################################\n");
  printf("load partition file : %s\n",regionsfile);
  status=mgr_SetRegions (regionsfile, mgr, nmgr, region);
  if(status !=0) {
    printf("cannot load region file=%s\n",regionsfile);
    goto error;
    }
  
  toponyms=new toponyms_t;
  status=geo_init_regions(toponyms);
  
  for(toponyms_t::iterator z=toponyms->begin();z!=toponyms->end();z++) {
    count=mgr_select(mgr, list, region, z->second);
    if(count==0) continue;
    printf("region %20s   \t \t",z->first);
    char regional_output[256];
    if(output==0) {
      sprintf(regional_output,"%s",z->first);
      }
    else {
      sprintf(regional_output,"%s.%s",output,z->first);
      }
    mgr_print(spectrum,  regional_output, mgr, count, list, topogrid, topo, topomask, a,  G, atlas_directory, atlas_convention, latex, show_diff, false, mean, rms);
    mgr_print2(spectrum, regional_output, mgr, count, list, topogrid, topo, topomask, a,  G, atlas_directory, atlas_convention, latex, show_diff);
    id++;
    summary_mean[id]=mean[0];
    summary_rms[id]=rms[0];
    }

    filename=tag+"-"+year+".sum";
    dump=fopen(filename.c_str(), "w");
    fprintf(dump, "%s %s ", tag.c_str(), year.c_str());
    for(int k=0;k<id+1;k++) fprintf(dump, "%lf ", summary_rms[k]);
    fprintf(dump,"\n");
    fclose(dump);
    
end:

/** to nest in a separate routine -------------------------------------------*/
   mgr_save_ascii("model.mgr", mgr_model);
/** to nest in a separate routine -------------------------------------------*/


  for (k=0;k<spectrum.n;k++) {
    delete[] a[k];
    delete[] G[k];
    }
  delete[] a;
  delete[] G;
  delete [] list;

  free(spectrum.waves);
  
  for(i=0;i<nmgr;i++) {
//     for(j=0;j<spectrum.n;j++) {
//       free(mgr[i]->data[j]);
//       }
    free(mgr[i].data);
    free(mgr[i].loc.units);
    }
  mgr.clear();
  
  return(0);
error:
  return(-1);
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
  const char *output=NULL, *keyword=NULL, *mgrfile=NULL;
  const char *ordering=NULL,*regions=NULL;
  char *wave[101];
  int nwave=0;
  date_t date;
  char *atlas_directory=NULL,*atlas_convention=NULL,*meshfile=NULL;
  string FrameString,plg_path,tag="X",year="1950";
  char show_diff=0,latex=0;
  double scale=1.0;
  bool silent=false;
  
  vector<plg_t> limits;
  
  char *bathymetry=NULL, *format=NULL;
  int structured=1;

  const int nvars=2;

  char *varnames[nvars]={0,0},*discretisation;
  int defaultVars=1,iteration=-1,level=-1;
  int strict=0;

  fct_echo( argc, argv);

//  fprintf(stderr,"%s  -starting computation *********\n",argv[0]);
  n=1;
  while (n < argc) {
    keyword=argv[n];
    switch (keyword[0]) {
      case '-':
        if(strcmp(keyword,"--silent")==0) {
          silent=true;
          n++;
          break;
          }
        if(strcmp(keyword,"-unstructured")==0) {
          discretisation= strdup(argv[n+1]);
          structured=0;
          n++;
          n++;
          break;
          }
        if(strcmp(keyword,"-strict")==0) {
          strict=1;
          n++;
          break;
          }
        if(strcmp(keyword,"-level")==0) {
          sscanf(argv[n+1],"%d",&level);
          n++;
          n++;
          break;
          }
        if(strcmp(keyword,"--regions")==0) {
          regions= argv[n+1];
          n++;
          n++;
          break;
          }
        if(strcmp(keyword,"--tag")==0) {
          tag= argv[n+1];
          n++;
          n++;
          break;
          }
        if(strcmp(keyword,"--year")==0) {
          year= argv[n+1];
          n++;
          n++;
          break;
          }
        if(strcmp(keyword,"-cm")==0) {
          scale=0.01;
          n++;
          break;
          }
        if(strncmp("--frame",keyword)==0){
          FrameString= argv[n+1];
          n++;
          n++;
          break;
          }
        if(strncmp("--polygons",keyword)==0){
          plg_path= argv[n+1];
          n++;
          n++;
          break;
          }
        switch (keyword[1]) {

/* *----------------------------------------------------------------------
        bathymetry database*/
        case 'b' :
          bathymetry= strdup(argv[n+1]);
          n++;
          n++;
          break;

/* *----------------------------------------------------------------------
        print data/models differences*/
        case 'd' :
          show_diff=1;
          n++;
          break;

/* *----------------------------------------------------------------------
        harmonic database format*/
        case 'f' :
          format= strdup(argv[n+1]);
          n++;
          n++;
          break;

/* *----------------------------------------------------------------------
        iteration index*/
        case 'i' :
          sscanf(argv[n+1],"%d",&iteration);
          n++;
          n++;
          break;

/* *----------------------------------------------------------------------
        output format*/
        case 'l' :
          latex=1;
          n++;
          break;

/* *----------------------------------------------------------------------
        naming convention for tidal atlases*/
        case 'm' :
          meshfile= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'a' :
          atlas_convention= strdup(argv[n+1]);
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
          output= argv[n+1];
          n++;
          n++;
          break;

/* *----------------------------------------------------------------------
        harmonic data input file*/
        case 'g' :
          mgrfile= argv[n+1];
          n++;
          n++;
          break;

/* *----------------------------------------------------------------------
        output order: FILE, ALPHA, LAT, LON, IMMERSION*/
        case 'x' :
          ordering= argv[n+1];
          n++;
          n++;
          break;

        case 'v' :
          if(varnames[0]){
            if(defaultVars)
              defaultVars=0;
            else{
              __FILE_LINE__(stdout,"multiple use of -v : the last one overrides the previous.");
              for(k=0;k<nvars;k++){
                free(varnames[k]);
                }
              }
            }
          for(k=0;k<nvars;k++){
            varnames[k]=strdup(argv[n+1+k]);
            }
          n++;
          n+=nvars;
          break;

        case 'h' :
          print_help(argv[0]);
          exit(0);
        default:
          printf("unknown option %s\n",keyword);
          print_help(argv[0]);
          exit(-1);
        }
        break;

      default:
/* *----------------------------------------------------------------------
          tidal wave list*/
          wave[nwave]= strdup(argv[n]);
//          printf("input wave=%s\n",wave[nwave]);
          //spectrum.n=nwave+1;
          nwave++;
          n++;
        break;
      }
    
    }

  if(atlas_directory==NULL)  atlas_directory=strdup(".");
  if(atlas_convention==NULL) {
    printf("*** Please specify atlas filename convention with -a ***\n");
    print_help(argv[0]);
    exit(-1);
    }
  wave[nwave]=NULL;
  
  status=init_mgrh_formats();
  if(format==0) format=strdup("LEGOS-ASCII");
  
  if(FrameString!="" and plg_path!="") {
    printf("*** Please specify EITHER --frame OR --polygons ***\n");
    print_help(argv[0]);
    exit(-1);
    }
  
  if(FrameString!="") {
    frame_t frame;
    status=plg_DecodeFrame(FrameString, frame);
    limits.push_back(plg_t(frame, PLG_INIT_SEPARATE));
    status=plg_cartesian((projPJ) 0, limits);
    }
  
  if(plg_path!="") {
    status=plg_load(plg_path,limits);
    }

  if(structured==1) {
    if(varnames[0]==0) varnames[0]=strdup("Ha");
    if(varnames[1]==0) varnames[1]=strdup("Hg");
    }
  else if(varnames[0]==0 or varnames[1]==0){
    printf("*** Please specify amplitude and phase variable names with -v ***\n");
    print_help(argv[0]);
    exit(-1);
    }
  
  int informat=MgrHarmonic_formats[format];
  status=mgr_validate((const char**)varnames, wave, output, atlas_directory, atlas_convention, scale, meshfile, mgrfile, informat, ordering, bathymetry, regions, limits,
                      latex, show_diff, structured, discretisation, iteration, level, strict, tag, year, silent);

  printf("\n%s - computation completed ^^^^^^^^^\n",argv[0]);
}
