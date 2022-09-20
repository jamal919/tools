

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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "version-macros.def" //for VERSION and REVISION

#include "tools-structures.h"

#include "geo.h"
#include "fe.h"
#include "map.h"
#include "polygones.h"
#include "grd.h"
#include "functions.h"

#include "zapper.h"     /*  rutin.h contains common utility routines  */


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  mesh_t *fe_split_G(mesh_t & mesh,int *selected, int *targeted, int *frontier, int size, string rootname, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  selected: edge flag for interior/exterior splitting (1 is interior)
  targeted:
  
  frontier: edge flag, 1 if located along the interior/exterior splitting limit
            setup if non-zero address

            
  size    : number of edges needed to fulfill criterion. if size > 0, then option=0
            else option=1
  
  option 0: select triangle as interior if at least "size" edges are eligible
  
    when size is 3, only triangles strictly included in the selection are taken as
    interior. when size is 1, boarder triangles are taken as interior
  
  option 1: select triangle as exterior if at least "size" edges are not eligible
  
    when size is 3, only triangles strictly excluded from the selection are taken as
    exterior. when size is 1, boarder triangles are taken as exterior
  
    non-symmetric condition
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int   l,m,m1,m2,i,n1,n2,n,status;
  int   count=0,element_size;
  int   nndes;
  int   ninternal,nexternal;
  bool  *keep;
  int   *used=NULL;
  mesh_t *work=NULL;
  string filename;
  int option;
  
  mesh.set_elements();
  
/*------------------------------------------------------------------------------
  to revise in case of composite mesh (triangles+quadrangles) */
  element_size=mesh.elements[0].size;

  work=new mesh_t[2];

  printf ("number of nodes    (original mesh): %6d\n",mesh.nvtxs);
  printf ("number of elements (original mesh): %6d\n",mesh.nelements);

  ninternal=0;
  nexternal=0;

  keep=new bool[mesh.nelements];
  
  if(size>0) option=0;
  else {
    option =1;
    size=-size;
    }

  for (m=0;m<mesh.nelements;m++) {
    keep[m]=false;
    count=0;
    switch (option) {
      case 0:
        for(i=0;i<element_size;i++) {
          n=mesh.elements[m].edges[i];
          if(selected[n]==1) {
            count++;
            }
          }
        keep[m]=(count>=size);
        break;
      case 1:
        for(i=0;i<element_size;i++) {
          n=mesh.elements[m].edges[i];
          if(selected[n]==0) {
            count++;
            }
          }
        keep[m]=not (count>=size);
        break;
      }
    switch(keep[m]) {
      case false:
        nexternal++;
        break;
      case true:
        ninternal++;
        break;
      }
    }

  switch(element_size) {
    case 3:
      work[0].triangles=new triangle_t [nexternal];
      work[1].triangles=new triangle_t [ninternal];
      work[0].ntriangles=nexternal;
      work[1].ntriangles=ninternal;
      work[0].type=SPHERICAL;
      work[1].type=SPHERICAL;
      break;
    case 4:
      work[0].quadrangles=new quadrangle_t [nexternal];
      work[1].quadrangles=new quadrangle_t [ninternal];
      work[0].nquadrangles=nexternal;
      work[1].nquadrangles=ninternal;
      work[0].type=SPHERICAL;
      work[1].type=SPHERICAL;
      break;
    }
    
//   work[0].elements=new element_t [nexternal];
//   work[1].elements=new element_t [ninternal];
// 
//   work[0].nelements=nexternal;
//   work[1].nelements=ninternal;
// 
//   work[0].type=SPHERICAL;
//   work[1].type=SPHERICAL;
    
  work[0].set_elements();
  work[1].set_elements();

/*------------------------------------------------------------------------------
  screen elements*/
  ninternal=0;
  nexternal=0;
  for (m=0;m<mesh.nelements;m++) {
    switch(keep[m]) {
      case false:
        work[0].elements[nexternal].ancestor=m;
        nexternal++;
        break;
      case true:
        work[1].elements[ninternal].ancestor=m;
        ninternal++;
        break;
      }
    }

  used=new int[mesh.nvtxs];
  for(l=0;l<2;l++) {
    nndes=0;
    for (n=0;n<mesh.nvtxs;n++) used[n]=-1;
    for (m=0;m<work[l].nelements;m++) {
      for(i=0;i<element_size;i++) {
        n=mesh.elements[work[l].elements[m].ancestor].vertex[i];
        if(used[n]==-1) {
          used[n]=nndes;
          nndes++;
          }
        }
      }
    work[l].nvtxs = nndes;
    work[l].vertices=new vertex_t[work[l].nvtxs];
    for (n=0; n<work[l].nvtxs; n++) work[l].vertices[n].null_value();
    for (n=0;n<mesh.nvtxs;n++) {
       if(used[n]!=-1) {
         work[l].vertices[used[n]].lon=mesh.vertices[n].lon;
         work[l].vertices[used[n]].lat=mesh.vertices[n].lat;
         work[l].vertices[used[n]].h  =mesh.vertices[n].h;
         work[l].vertices[used[n]].code=0;
         work[l].vertices[used[n]].ancestor=n;
         }
       }
    for (m=0;m<work[l].nelements;m++) {
      for(i=0;i<element_size;i++) {
        n=mesh.elements[work[l].elements[m].ancestor].vertex[i];
        if(used[n]!=-1) {
          work[l].elements[m].vertex[i]=used[n];
          }
        else {
          printf("error...\n");
          }
        }
      }
    printf ("number of nodes    (splitted mesh %d): %6d\n",l,work[l].nvtxs);
    printf ("number of elements (splitted mesh %d): %6d\n",l,work[l].nelements);

    status=fe_e2n(&(work[l]));
    if(debug) status=fe_savemesh("tmp-00.nei",MESH_FILE_FORMAT_TRIGRID,work[l]);

    status=fe_geometry(&(work[l]));
    status=fe_edgetable(&(work[l]),0,0);
    status=fe_vertex_crosstables02(&(work[l]));

    int RecycleCodes=0, stopon_EdgeError=1, stopon_PinchError=0;
    status=fe_codetable1(&(work[l]), RecycleCodes, stopon_EdgeError, stopon_PinchError);
    if(status!=0) {
      printf ("code table failed for mesh %d\n",l);
      return(0);
      }
      
    if((l==1) && (frontier!=0)) {
/*------------------------------------------------------------------------------
      no information on edge's ancestor, use elements*/
      for (m1=0;m1<work[l].nelements;m1++) {
        m2=work[l].elements[m1].ancestor;
        for(i=0;i<3;i++) {
          n1=work[l].elements[m1].edges[i];
          n2=mesh.elements[m2].edges[i];
          if((mesh.edges[n2].code==MESH_INTERIOR_EDGE) && (work[l].edges[n1].code==1)) {
/*------------------------------------------------------------------------------
            internal/external limit, do not refine*/
            frontier[n1]=1;
            }
          else {
            frontier[n1]=0;
            }
          }
        }
      }
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  save interior and exterior meshes
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if(rootname=="") rootname="anonymous";
  
  filename=rootname+"-external.nei";
  status=fe_savemesh(filename.c_str(),MESH_FILE_FORMAT_TRIGRID,work[0]);
  
  filename=rootname+"-external.plg";
  status=fe_SaveLimits(filename.c_str(),work[0]);
  
  filename=rootname+"-internal.nei";
  status=fe_savemesh(filename.c_str(),MESH_FILE_FORMAT_TRIGRID,work[1]);
  
  filename=rootname+"-internal.plg";
  status=fe_SaveLimits(filename.c_str(),work[1]);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  create boundary polygon file for further mesh assembly
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  l=1;
  char *flag=new char[work[l].nedges];
  
  count=0;
  for(int n=0;n<work[l].nedges;n++) {
    if(work[l].edges[n].code==MESH_INTERIOR_EDGE) flag[n]='I';
    else flag[n]='T';
    }
    
  for (m1=0;m1<work[l].nelements;m1++) {
    m2=work[l].elements[m1].ancestor;
    for(i=0;i<element_size;i++) {
      n1=work[l].elements[m1].edges[i];
      n2=mesh.elements[m2].edges[i];
/*------------------------------------------------------------------------------
      02/02/2016: may have side-effects... */
      if((mesh.edges[n2].code==MESH_INTERIOR_EDGE) && (work[l].edges[n1].code>=1)) {
/*------------------------------------------------------------------------------
        internal/external limit, do not refine*/
        flag[n1]='M';
        work[l].edges[n1].code=MESH_FLAGGED_EDGE;
        count++;
        }
      }
    }
    
/*------------------------------------------------------------------------------
  export mesh limits as polygons, passing open/rigid flag */
  vector<plg_t> polygons, splitted, opened;
  status=fe_limits2poly(work[l], polygons, flag, false);
  
/*------------------------------------------------------------------------------
  split polygons to separate open/rigid segments */
  splitted=plg_split(polygons, false);
  
  for(int s=0;s<splitted.size();s++) {
    if(splitted[s].flag[0]=='M') opened.push_back(splitted[s]);
    }
  filename=rootname+"-opened.plg";
  status=plg_save(filename.c_str(), PLG_FORMAT_SCAN, opened);
  
  delete[] used;

  printf("split completed\n");

  return(work);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int delaunay(mesh_t & mesh, const char *rootname)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int n,status;
  mesh_t finished;
  triangulateio in, out;
  char output[1024];
  
  for(n=0;n<mesh.nvtxs;n++) {
    if(mesh.vertices[n].code==1) mesh.vertices[n].code=0;
    }
  status= fe_triangulateio_init(in);
  status= fe_triangulateio_init(out);
  
  status= fe_node2triangle(mesh, &in);
  
  triangulate("z", &in, &out, (triangulateio *) NULL);

/* *-----------------------------------------------------------------------------
  output structure inherit from input structure's pointlist */
  out.pointlist = in.pointlist;

  status=fe_triangle2mesh(out,&finished);
  
  sprintf(output, "%s-delaunay.nei",rootname);
  status=fe_savemesh(output,MESH_FILE_FORMAT_TRIGRID,finished);
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_create_circle(plg_array_t *plg, double radius, double t0, double p0)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int n,status;
  int i,j,k,utm_zone,ntoken;
  double *x,*y,x0,y0;
  double L,H;
  char *proj4_options;
  const char *hemisphere;
  projPJ proj;
  char **parms;
  
  L=2.*M_PI*radius;

  plg->n=1;
  plg->p=new plg_t[plg->n];

  plg->p[0].npt=100;
  
  x=new double[plg->p[0].npt];
  y=new double[plg->p[0].npt];

  plg->p[0].t=new double[plg->p[0].npt];
  plg->p[0].p=new double[plg->p[0].npt];
  
  if(p0<0) hemisphere="+south";
  else hemisphere="";
  
  utm_zone=(t0+180)/6+1;
  
  asprintf(&proj4_options,"+proj=utm +zone=%d +ellps=WGS84 +datum=WGS84 +units=m +no_defs %s",utm_zone,hemisphere);

  proj = pj_init_plus(proj4_options);
  if (!proj) __TRAP_ERR_EXIT__(1,"Projection initialization failed\n");
  
  geo_to_projection(proj, p0, t0, &x0, &y0);
  pj_free(proj);
  
  for(n=0;n<plg->p[0].npt;n++) {
    double alpha=n*2*M_PI/(plg->p[0].npt-1.);
    x[n]=x0+radius*cos(alpha);
    y[n]=y0+radius*sin(alpha);
    }
    
  status=projection_to_geo (proj4_options, x, y, plg->p[0].npt);
  
  free(proj4_options);
  
  for(n=0;n<plg->p[0].npt;n++) {
    plg->p[0].t[n]=x[n];
    plg->p[0].p[n]=y[n];
    }

  status=plg_save("circle.plg",PLG_FORMAT_SCAN,plg->p,1);
  
  return(status);
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
    "NAME AND VERSION\n"
    "  %s version " VERSION " " REVISION "\n", prog_name);
  printf("\n"
    "USE\n"
    "  %s file1 [ file2 ... ] [OPTIONS]\n",prog_name);
  printf("\n"
    "DESCRIPTION\n"
/* USE FIND AND REPLACE TO EDIT THE LINE BELOW, SO THAT THE FILE HEAD IS ALSO UPDATED */
    "  refine mesh.\n"
    "\n"
    "OPTIONS :\n"
    "  -h,--help  Show this help and exit.\n"
    "  -m  followed by the path of the mesh input file.\n"
    "  -p  followed by the path of the polygon input file.\n"
    "  -o  followed by the path of the output file. Default from the input file and output format if both are given.\n"
    "  --exact  exactly follow polygons limits\n");
  printf("\n"
    "FILE FORMATS :\n");
  plg_print_formats();
  /** \endcode */
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int nndes,status,*elts,fmt,fefmt;
  int i,j,k,l,L,m,n;
  int nselected;
  int option,channels=0;
  FILE *file;
  FILE *out;

  char *meshfile=NULL,*depthfile=NULL,*output=NULL,*poly=NULL;
  mesh_t mesh,refined,internal,external,*splitted,*final;
  double dmax=0,cmax=0;
  int *selected,*targeted;
  bool exact=false;
  string keyword="",zone="",rootname,filename;
  plg_array_t plg;
  int RecycleCodes=0, StopOn_EdgeError=1, StopOn_PinchError=0, SetLimits=1;
  int element=FE_TRIANGLE;
  int verbose=0;
  bool debug=false, polar=false;

  fct_echo(argc,argv);

  n=1;
  while (n < argc) {
    keyword=argv[n];
    switch (keyword[0]) {
      case '-':
        switch (keyword[1]) {
        case 'z' :
          zone=argv[n+1];
          n++;
          n++;
          break;

        case 'm' :
          meshfile= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'o' :
          output= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'p' :
          poly= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'r' :
          rootname=argv[n+1];
          n++;
          n++;
          break;

        case '-' :
          if(keyword=="--help"){
            print_help(argv[0]);
            exit(0);
            }
          else if(keyword=="--quadrangle"){
            element=FE_QUADRANGLE;
            n++;
            }
          else if(strncmp("--exact",keyword)==0){
            exact=true;
            n++;
            }
          else if(strncmp("--polar",keyword)==0){
            polar=true;
            n++;
            }
          else if(strncmp("--debug",keyword)==0){
            debug=true;
            n++;
            }
          else {
            printf("unknown option "+keyword+"\n");
            print_help(argv[0]);
            exit(-1);
            }
          break;

        default:
          printf("unknown option "+keyword+"\n");
          print_help(argv[0]);
          exit(-1);
          break;
        }
        break;
      default:
          __OUT_BASE_LINE__("unknown option %s\n",keyword.c_str());
          exit(-1);
      }
    }
  
//   status= plg_create_circle (&plg, 90000., 52, -46.5);
//   status= plg_create_ellipse(&plg, 60000., 90000., 52, -46.5);

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  load mesh
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  printf("#################################################################\n");
  printf("load mesh ant initialize related tables\n");
  if(meshfile != NULL) {
    status=fe_readmesh(meshfile,MESH_FILE_FORMAT_TRIGRID,&mesh);
    if(status!=0) {
      printf("unable to read the original mesh in %s\n",meshfile);
      goto error;
      }
    status=fe_list(&mesh, element);
    if(status!=0) {
      printf("unable to build the element list from the original mesh\n");
      goto error;
      }
    }
 else {
   printf("no mesh file specified; abort...\n");
   goto error;
   }

  status= fe_edgetable(&mesh,0,0);
  if(status!=0) {
    printf("unable to build the edge list from the original mesh\n");
    goto error;
    }

//   RecycleCodes=0; StopOn_EdgeError=1; StopOn_PinchError=0;
//   status= fe_codetable2(&mesh, RecycleCodes, StopOn_EdgeError, StopOn_PinchError);

  RecycleCodes=0; SetLimits=1; verbose=1;
  status=fe_codetable(mesh, RecycleCodes, SetLimits, verbose);
  if(status!=0) {
    printf("unable to rebuild the limits table and codes of the original mesh\n");
    goto error;
    }
    
  {
  string LimitsFile=meshfile;
  LimitsFile=replace(LimitsFile,(string) ".nei",(string) ".plg",-1);
  
  vector<plg_t> limits, selection;
  status=fe_limits2poly(mesh, limits, (char *) 0, true);  
  status=plg_save(LimitsFile.c_str(), PLG_FORMAT_SCAN, limits);
  
  status=plg_load(poly, PLG_FORMAT_SCAN, selection);
  status=plg_CutAtIntersections(limits, selection, PLG_SPHERICAL);    
  status=plg_save("dump-01.plg", PLG_FORMAT_SCAN, limits);
  status=plg_save("dump-02.plg", PLG_FORMAT_SCAN, selection);
  }
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  optionally, prepare mesh nodes to fit exactly polygons
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if(exact) {
    status=fe_presplit(mesh, poly, 0.5, "split-adjusted");
    if(status!=0) goto error;
    }
    
  selected=new int[mesh.nedges];
  for (n=0;n<mesh.nedges;n++) {
    selected[n]=1;
    }
  
  if(zone!="") {
    frame_t frame;
    point2D_t point;
    vector<plg_t> limits;
    status=plg_SetZone(zone.c_str(), frame, rootname, point);
    if(status!=0) goto error;
    limits.push_back(plg_t(frame, PLG_INIT_SEPARATE));
    status=plg_cartesian((projPJ) 0, limits);
    status=plg_save("mesh-split.plg", PLG_FORMAT_SCAN, limits);
    poly=strdup("mesh-split.plg");
//    extract=1;
    }
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  identify edges in polygon
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  printf("#################################################################\n");
  printf("select edges from polygons\n");
  
  nselected=fe_selectedges_01(mesh, poly, selected, polar);

  if(nselected<=0) goto error;
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  create mesh partition
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  if(rootname=="") rootname="mesh-split";
    
  printf("#################################################################\n");
  printf("split mesh\n");
//   splitted=fe_split(mesh, selected, selected, 0, 3, rootname, debug);
  splitted=fe_split_G(mesh, selected, selected, 0, 3, rootname, debug);
  
  if(exact) {
    delaunay(splitted[0],"external");
    delaunay(splitted[1],"internal");
    }
    
end: __OUT_BASE_LINE__("end of mesh-split ... \n");
  exit(0);

error:
  __OUT_BASE_LINE__("error detected, quit ... \n");
  exit(-1);
}
