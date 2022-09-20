

/**************************************************************************

  T-UGO tools, 2006-2009

  Unstructured Ocean Grid initiative

Contributors:

  Florent Lyard      LEGOS/CNRS, Toulouse, France
  Yoann Le Bars      LEGOS, Toulouse, France (PhD)
  Laurent Roblou     LEGOS, Toulouse, France
  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
  Frédéric Dupont    Université de Laval à Québec, Canada
  vector<int> nodes;

E-mail: florent.lyard@legos.obs-mip.fr

***************************************************************************/

#define MAIN_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include "tools-structures.h"

#include "geo.h"
#include "fe.h"
#include "fe-proto.h"
#include "map.h"
#include "polygones.h"
#include "grd.h"
#include "functions.h"
#include "statistic.h"

#include "zapper.h"     /*  rutin.h contains common utility routines  */


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int remove_duplicates(vector<int> & v)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,n,nvalues,status;
  int *used;
  range_t<int> r=poc_minmax(v,-1);
  
//   printf("initial size=%d\n",v.size());
  
  nvalues=r.max+1;
  
  used=new int[nvalues];
  
  for(n=0;n<nvalues;n++)     used[n]=-1;

//   for(k=0;k<v.size();k++) {
  for(std::vector<int>::iterator s=v.begin(); s<v.end();s++) {
    int m=*s;
    if(used[m]==-1)  {
      used[m]=1;
      }
    else {
      v.erase(s);
      }
    }

  delete[] used;
//   for(std::vector<int>::iterator s=v.begin(); s<v.end();s++) {
//     int m=*s;
//     for(std::vector<int>::iterator r=s+1;r<v.end();r++) {
//       int n=*r;
//       if(m==n) v.erase(r);
//       }
//     }
//   printf("final size=%d\n",v.size());
  
  return(status);
  
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_submesh(mesh_t & mesh, vector<int> elements, const char *filename, mesh_t *work)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int m, n, status;
  int s, smax;
  
//  mesh_t *work=new mesh_t;
  
  vector<int> nodes;
  vector <vector<int> > partition;

  int *kept=new int[mesh.nvtxs];
  for(n=0;n<mesh.nvtxs;n++) kept[n]=-1;

  for(s=0;s<elements.size();s++) {
    m=elements[s];
    for(int k=0;k<3;k++) {
      n=mesh.triangles[m].vertex[k];
      if(kept[n]==-1) {
        kept[n]=nodes.size();
        nodes.push_back(n);
        }
      }
    }
  
  work->nvtxs=nodes.size();
  work->vertices=new vertex_t[work->nvtxs];
  
  work->ntriangles=elements.size();
  work->triangles=new triangle_t[work->ntriangles];
  
  work->type=mesh.type;
 
  for(s=0;s<elements.size();s++) {
    m=elements[s];
    for(int k=0;k<3;k++) {
      n=mesh.triangles[m].vertex[k];
      work->triangles[s].vertex[k]=kept[n];
      }
    work->triangles[s].ancestor=m;
    }
  
  for(s=0;s<nodes.size();s++) {
    n=nodes[s];
    work->vertices[s].lon=geo_recale(mesh.vertices[n].lon, 0.0, 180.);
    work->vertices[s].lat=mesh.vertices[n].lat;
    work->vertices[s].code=0;
    work->vertices[s].ancestor=n;
    }
  
  status=fe_e2n(work);

  status=fe_edgetable(work,0,0);
  
//   status=fe_vertex_element_tables(work);
//   
//   status=fe_GetConnex(*work, partition);
  
  if(filename!=0) status=fe_savemesh(filename,MESH_FILE_FORMAT_TRIGRID,*work);
  
  delete[] kept;
  
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int vertex_copy(vertex_t & target, const vertex_t & vertex)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int m, n, status;
  
  target.destroy();
  
  target.lon   = vertex.lon;
  target.lat   = vertex.lat;
  target.h     = vertex.h;
  target.code  = vertex.code;
  target.nngh  = vertex.nngh;
  
  target.ngh=new int[target.nngh];
  
  for(int k=0;k<vertex.nngh; k++) {
    target.ngh[k]=vertex.ngh[k];
    }

  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int edge_copy(edge_t & target, const edge_t & edge)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int m, n, status;
  
  target.destroy();
  
  target.lon   = edge.lon;
  target.lat   = edge.lat;
  target.extremity[0]  = edge.extremity[0];
  target.extremity[1]  = edge.extremity[1];

  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int triangle_copy(triangle_t & target, const triangle_t & triangle, bool reverse)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int m, n, status;
  
  target.destroy();
  
  if(reverse) {
    for(int k=0;k<3; k++) {
      target.vertex[k]=triangle.vertex[2-k];
      }
    }
  else {    
    for(int k=0;k<3; k++) {
      target.vertex[k]=triangle.vertex[k];
      }
    }
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_flip(const mesh_t & mesh, mesh_t & flip)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int m, n, status;
  double offset=0;
  
  flip.nvtxs=mesh.nvtxs;
  flip.ntriangles=mesh.ntriangles;

  flip.vertices  = new vertex_t[flip.nvtxs];
  flip.triangles = new triangle_t[flip.ntriangles];
 
  for (n=0;n<mesh.nvtxs;n++) {
    vertex_copy(flip.vertices[n], mesh.vertices[n]);
    flip.vertices[n].lat=-(flip.vertices[n].lat+offset);
    }
  
  for (n=0;n<mesh.ntriangles;n++) {
    triangle_copy(flip.triangles[n], mesh.triangles[n], true);
    }
  
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_translate(const mesh_t & mesh, mesh_t & flip, double offset)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int m, n, status;
  
  flip.nvtxs=mesh.nvtxs;
  flip.ntriangles=mesh.ntriangles;

  flip.vertices  = new vertex_t[flip.nvtxs];
  flip.triangles = new triangle_t[flip.ntriangles];
 
  for (n=0;n<mesh.nvtxs;n++) {
    vertex_copy(flip.vertices[n], mesh.vertices[n]);
    flip.vertices[n].lat=flip.vertices[n].lat+offset;
    }
  
  for (n=0;n<mesh.ntriangles;n++) {
    triangle_copy(flip.triangles[n], mesh.triangles[n], false);
    }
  
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_symetric(const mesh_t & mesh, mesh_t & symetric)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int m, n, status;
  double dmax=0.01;
  mesh_t flip, translate;
  double offset=0;
  int RecycleCodes=0, stopon_EdgeError=1, stopon_PinchError=0;
  
  offset=7.99999982e-02;
  
  status=fe_translate(mesh, translate, offset);
  status=fe_edgetable(&translate,0,0);
  status=fe_codetable1(&translate, RecycleCodes, stopon_EdgeError, stopon_PinchError);
  
  status=fe_flip(translate, flip);
  
  status=fe_edgetable(&flip,0,0);
  status=fe_codetable1(&flip, RecycleCodes, stopon_EdgeError, stopon_PinchError);
    
  printf("#################################################################\n");
  printf("merge meshes, maximum merging distance=%lf km\n",dmax);
  symetric=fe_merge(translate, flip, dmax);

  status=fe_edgetable(&symetric,0,0);
  status=fe_codetable1(&symetric, RecycleCodes, stopon_EdgeError, stopon_PinchError);
  
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_relocate(mesh_t & mesh, double offset)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int m, n, status;
  frame_t frame;

//   double scale=1./cos(offset*M_PI/180.);
  
  double scale=110.0/100.0;
  offset=acos(1./scale)*180./M_PI;
  
  for (n=0;n<mesh.nvtxs;n++) {
    mesh.vertices[n].lon*=scale;
    mesh.vertices[n].lat+=offset;
    }
  
  status=fe_minmax(mesh, frame);
  
//   frame.print();
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int nndes,status,*elts,fmt,fefmt;
  int i,j,k,l,L,m,n;
  int option,channels=0;
  FILE *file;
  FILE *out;
  char *keyword,*zone=NULL,*gridfile=NULL;
  char *meshfile=NULL,*depthfile=NULL,*output=NULL,*poly=NULL,*belfile=NULL;
  mesh_t mesh,final;
  double dmax=0,cmax=0;
  double offset;
  int *selected=0, *targeted=0, *frontier=0, nselected;
  bool debug;

  int RecycleCodes, stopon_EdgeError, stopon_PinchError;
  
  fct_echo(argc,argv);

  n=1;
  while (n < argc) {
    keyword=strdup(argv[n]);
    switch (keyword[0]) {
      case '-':
        switch (keyword[1]) {
        case 'z' :
          zone= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'm' :
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

        case 'd' :
          sscanf(argv[n+1],"%lf",&dmax);
          n++;
          n++;
          break;

        case 'e' :
          sscanf(argv[n+1],"%lf",&cmax);
          n++;
          n++;
          break;

        case 'c' :
          channels=1;
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

        case '-' :
//           if(keyword=="--help"){
//             print_help(argv[0]);
//             exit(0);
//             }
          if(strncmp("--debug",keyword)==0){
            debug=true;
            n++;
            }
          else {

            exit(-1);
            }
          break;

           
        default:
          __OUT_BASE_LINE__("unknown option %s\n",keyword);
          exit(-1);
        }
        break;

      default:
        n++;
        break;
      }
    free(keyword);
    }

  if(output==0) output=strdup("mesh-editor.nei");

  if(meshfile != NULL) {
    status=fe_readmesh(meshfile,MESH_FILE_FORMAT_TRIGRID,&mesh);
    if(status!=0) {
      printf("unable to read the original mesh in %s\n",meshfile);
      goto error;
      }
    status=fe_list(&mesh);
    if(status!=0) {
      printf("unable to build the element list from the original mesh\n");
      goto error;
      }
    }
  else {
    printf("no mesh file specified; abort...\n");
    goto error;
    }

  nndes=mesh.nvtxs;
  status=fe_edgetable(&mesh,0,0);
  if(status!=0) {
    printf("unable to build the edge list from the original mesh\n");
    goto error;
    }

  option=3;
     
  switch (option) {
    case 0:

      mesh.destroy();
      break;
    }

//   printf("#################################################################\n");
//   printf("reshape final mesh\n");
//   status= fe_reshapeall(final,3);

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  COMODO mesh construction : symetric mesh 

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

//   status=fe_symetric(mesh, final);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  COMODO mesh construction : relocate mesh

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  offset=46.95;
  status=fe_relocate(mesh, offset);
  final=mesh;
  
  printf("#################################################################\n");
  printf("save reshaped mesh : %s\n", output);
  status= fe_savemesh(output,MESH_FILE_FORMAT_TRIGRID, final);

end: __OUT_BASE_LINE__("end of refine ... \n");
  exit(0);

error:
  __OUT_BASE_LINE__("error detected, quit ... \n");
  exit(-1);
}
