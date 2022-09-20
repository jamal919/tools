

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

/* *----------------------------------------------------------------------------

  mesh-doctor aims at improving a finite element mesh (nei format)


-----------------------------------------------------------------------------*/

#define MAIN_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
 

#include "tools-structures.h"
#include "constants.h"

#include "geo.h"
#include "fe.h"
#include "map.h"
#include "polygones.h"
#include "grd.h"
#include "archive.h"
#include "functions.h"
#include "statistic.h"
#include "matrix.h"

#include "zapper.h"     /*  rutin.h contains common utility routines  */

#define S2R 0
#define S2C 1
#define V2R 3
#define V2C 4

static float *flags;


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int clean_angles(mesh_t & mesh, double threshold, bool reshape, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *----------------------------------------------------------------------------
  Detect small angles (vertex-wise)
-----------------------------------------------------------------------------*/
{
  int   e,e1,e2,i,m,n,n1,n2,n3,status;
  int   count=0;
  char *passed;
  double l1,l2,angle[3],r;

  passed=new char[mesh.ntriangles];
  
  for (m=0;m<mesh.ntriangles;m++) passed[m]=0;
   
/*-----------------------------------------------------------------------------
  screen large angles*/
  for (m=0;m<mesh.ntriangles;m++) {
    if(passed[m]==1) continue;
    for (i=0;i<3;i++) {
      n1=mesh.triangles[m].vertex[i];
      n2=mesh.triangles[m].vertex[(i+2)%3];
      n3=mesh.triangles[m].vertex[(i+1)%3];
      angle[i]=fe_angle(mesh, n1, n2, n3);
      if(fabs(angle[i])>threshold*d2r) {
        e1=mesh.triangles[m].edges[(i+1)%3];
        e2=mesh.triangles[m].edges[(i+2)%3];
        l1=mesh.edges[e1].L;
        l2=mesh.edges[e2].L;
        r=abs(l2-l1)/(l2+l1);
/*-----------------------------------------------------------------------------
        recut triangles*/
        if(r<0.2) {
          e=mesh.triangles[m].edges[i];
	  if(mesh.edges[e].nshared!=2) continue;
//           if(passed[mesh.edges[e].shared[0]]==1) continue;
//           if(passed[mesh.edges[e].shared[1]]==1) continue;
          status=fe_exchange(mesh, e,reshape, debug);
          if(status==0) {
            count++;
//             passed[mesh.edges[e].shared[0]]=1;
//             passed[mesh.edges[e].shared[1]]=1;
            break;
            }
          }
        }
      }
finish:
    passed[m]=1;
    }
  printf ("number of elements checked %d\n",count);
  
  delete[] passed;
  return(count);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int clean_delaunay01(mesh_t & mesh, double maxsize, bool exterior_only, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  Detect small angles (vertex-wise)
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int   k,m,n,n1,n2,n3,status;
  int   count=0,*passed;
  double d;

/**-----------------------------------------------------------------------------
  pre-filter elements*/
//  resolution=0.1;
  count=0;
  for(n=0;n<mesh.nvtxs;n++) {
    for(k=0;k<mesh.vertices[n].nngh;k++) {
      m=mesh.vertices[n].ngh[k];
      d=fe_distance(mesh, n, m);
/*------------------------------------------------------------------------------
      warning : fe_distance now return distance in meters */
      if(d>maxsize) {
        int e=fe_isedge(mesh, n, m);
        if(mesh.edges[e].code==MESH_INTERIOR_EDGE and exterior_only) continue;
        status=fe_disconnectvertices(mesh, n,m);
        k=k-1;
        count++;
        }
      }
    }

  printf("remove oversized connections (L > %fm) : %d\n", maxsize, count);
  status=fe_cleanvertices(&mesh,debug);
  status=fe_list(&mesh);
  status=fe_e2n(&mesh);
  status=fe_edgetable(&mesh,0,0);
  int RecycleCodes=0, stopon_EdgeError=1, stopon_PinchError=0;
  status=fe_codetable2(&mesh, RecycleCodes, stopon_EdgeError, stopon_PinchError);
  if(debug)
    status=fe_savemesh("out-01.nei",MESH_FILE_FORMAT_TRIGRID,mesh);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int clean_delaunay02(mesh_t & mesh, bool exterior_only, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  Detect small angles (vertex-wise)
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int   k,l,m,n,n1,n2,n3,status;
  int   count=0,*passed;
  double *esize;
  
  if(mesh.nedges==0) return(-1);
    
  esize=new double[mesh.nedges];
  
  for(n=0;n<mesh.nedges;n++) esize[n]=mesh.edges[n].L;
  statistic_t s=get_statistics(esize, (double) -1, mesh.nedges, 0);

/**-----------------------------------------------------------------------------
  pre-filter elements*/
  count=0;
  for(n=0;n<mesh.nvtxs;n++) {
    vector<double> lengths;
    if (mesh.vertices[n].nngh<4) continue;
    double d;
    for(k=0;k<mesh.vertices[n].nngh;k++) {
      m=mesh.vertices[n].ngh[k];
/*------------------------------------------------------------------------------
      warning : fe_distance now return distance in meters */
      d=fe_distance(mesh, n, m);
      lengths.push_back(d);
      }
    statistic_t s=get_statistics(lengths, (double) -1, 0);
    for(l=0;l<mesh.vertices[n].nngh;l++) {
      if (mesh.vertices[n].nngh<4) break;
      m=mesh.vertices[n].ngh[l];
      if(m<n) continue;
      if (mesh.vertices[m].nngh<4) continue;
      int e=fe_isedge(mesh, n, m);
      if(mesh.edges[e].code==MESH_INTERIOR_EDGE and exterior_only) continue;
      d=fe_distance(mesh, n, m);
      if(d > 5.0*s.min) {
        status=fe_disconnectvertices(mesh, n, m);
        l=l-1;
        count++;
        }
      }
    lengths.clear();
    }

  printf("remove oversized connections (L > %f x Lmin) : %d\n", 5.0, count);
  
  if(count==0) {
    return(0);
    }
    
  status=fe_cleanvertices(&mesh, false);
  status=fe_list(&mesh);
  status=fe_e2n(&mesh);
  status=fe_edgetable(&mesh,0,0); 
  int RecycleCodes=0, stopon_EdgeError=1, stopon_PinchError=0;
  status=fe_codetable2(&mesh, RecycleCodes, stopon_EdgeError, stopon_PinchError);
  if(debug)
    status=fe_savemesh("out-02.nei",MESH_FILE_FORMAT_TRIGRID,mesh);
  
  return(count);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int clean_delaunay03(mesh_t & mesh, bool exterior_only, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  Detect small edges ()
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int   k,l,m,n,n1,n2,n3,status;
  int   count=0,*passed;
  double *esize;
  
  if(mesh.nedges==0) return(-1);
    
  esize=new double[mesh.nedges];
  
  for(n=0;n<mesh.nedges;n++) esize[n]=mesh.edges[n].L;
  statistic_t s=get_statistics(esize, (double) -1, mesh.nedges, 0);

/**-----------------------------------------------------------------------------
  */
  count=0;
  for(n=0;n<mesh.nvtxs;n++) {
    redo:
    if (mesh.vertices[n].code==0) continue;
    if (mesh.vertices[n].nngh==0) continue;
    if (mesh.vertices[n].nngh<5) continue;
    double *d=new double[mesh.vertices[n].nngh];
    for(k=0;k<mesh.vertices[n].nngh;k++) {
      m=mesh.vertices[n].ngh[k];
      d[k]=fe_distance(mesh, n, m);
      }
    statistic_t s=get_statistics(d, (double) -1, mesh.vertices[n].nngh, 0);
    vector<double> lengths;
    for(k=0;k<mesh.vertices[n].nngh;k++) {
      if(d[k]!=s.max) lengths.push_back(d[k]);
      }
    s=get_statistics(lengths, (double) -1, 0);
    double median=median_value(lengths);
    for(l=0;l<mesh.vertices[n].nngh;l++) {
      if (mesh.vertices[n].nngh<4) break;
      if(d[l] > 2.0*median) {
        m=mesh.vertices[n].ngh[l];
        if (mesh.vertices[m].nngh<4) break;
        int e=fe_isedge(mesh, n, m);
        if(mesh.edges[e].code==MESH_INTERIOR_EDGE and exterior_only) continue;
        status=fe_disconnectvertices(mesh, n,m);
        count++;
        goto redo;
        }
      }
    delete[] d;
    }

  printf("remove oversized connections (L > mean + %f x rms) : %d\n", 2.5, count);
  if(count==0) {
    return(0);
    }
    
  status=fe_cleanvertices(&mesh, false);
  status=fe_list(&mesh);
  status=fe_e2n(&mesh);
  status=fe_edgetable(&mesh,0,0);
  if(debug)
    status=fe_savemesh("out-03.nei",MESH_FILE_FORMAT_TRIGRID,mesh);
  
  return(count);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_analyseBCflags(mesh_t & mesh, int target)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  Retrieve mesh boundary structure from mesh.edges[].code and mesh.limits[]
  
  Useful to reconstruct bel file after mesh editing
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/ 
{
  int k,kk,l,m,n;
  int start,upto;
  size_t count;
  int n1, n2;
  
  m=0;
  count=0;
  for(l=0;l<mesh.nlimits;l++) {
    vector< range_t<int> > partition;
    int chk=0;
    printf("analysing mesh limit %d\n",l);
    int nedges=mesh.limits[l].nedges;
    vector<int> codes;
    for(k=0;k<mesh.limits[l].nedges;k++) {
      n=mesh.limits[l].edges[k];
      codes.push_back(mesh.edges[n].code);
      }
/*------------------------------------------------------------------------------
    identify different types of flag in the current limit */
    vector<int> flags;
    for(k=0;k<nedges;k++) {
      int flag=codes[k];
      if(vpos(flag,flags)==-1) {
        flags.push_back(flag);
        }
      }
    printf("#number of computational flags=%d\n",flags.size());
    if(flags.size()==1) {
      range_t<int> r(0,nedges);
      continue;
      }
/*------------------------------------------------------------------------------
    identify differents types of segments in the current limit */
    k=0;
    bool completed=false;
    while(!completed) {
      start=(k+1) % nedges;
      int flag=codes[k];
      while(codes[start]==flag) {
        start=(start+1) % nedges;
        if(start==0) break;
        }
      upto=start;
      flag=codes[start];
      count=0;
      while(codes[upto]==flag) {
        upto=(upto+1) % nedges;
        count++;
        if(partition.size()!=0) {
          if(upto==partition[0].min) {
            completed=true;
            break;
            }
          }
        }
      upto--;
      range_t<int> r(start,upto);
      partition.push_back(r);
      k=(upto) % nedges;
      printf("computational flag=%d, start=%d upto=%d (%d edges)\n",flag,start,upto, count);
      int m1=mesh.limits[l].edges[r.min];
      int m2=mesh.limits[l].edges[r.max];
      n1=mesh.edges[m1].extremity[0];
      n2=mesh.edges[m2].extremity[1];
      printf("%lf %lf %lf %lf    %1d\n", mesh.vertices[n1].lon, mesh.vertices[n1].lat, mesh.vertices[n2].lon, mesh.vertices[n2].lat, flag);
      chk+=count;
      }
    printf("done\n");
    }

//   printf("%d edges received computational flag=%d\n",count,flag);
  
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int redep(const mesh_t & input, mesh_t & ouput)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  grid_t grid;
  list_t list;
  double *buffer=new double[input.nvtxs];
  
  for (int n=0;n<input.nvtxs;n++) buffer[n]=input.vertices[n].h;
    
  fe_Allocate_and_CreateList(input ,&grid, &list);
  
#pragma omp parallel for
  for (int n=0;n<ouput.nvtxs;n++) {
    int status;
    double t,p,z;
    int k,l,element;
    t=ouput.vertices[n].lon;
    p=ouput.vertices[n].lat;
    t=map_recale(grid,t);
    k=(int)( floor((t-grid.xmin)/grid.dx) );
    l=(int)( floor((p-grid.ymin)/grid.dy) );
    element=fe_whichelement_inlist(input,list.elements[k][l], (double) t,(double) p);
    if(element>0) {
      status= fe_intpl_LGP1(input, buffer, t,p,element,&z);
      }
    else {
      int m;
      m=fe_nearest_vertex (input,t,p,0,0);
      z=buffer[m];
      }
    ouput.vertices[n].h=z;
    }

  list.destroy();
  grid.free();
  
  delete[] buffer;
  
  return(0);
  
}
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int nndes,status,*elts,fmt,fefmt;
  int i,j,k,l,L,m,n;
  int option,channels=0,geometry=0,safety=0,renum=1,reshape=1,nghmax=7;
  int delaunay=0;
  FILE *file;
  FILE *out;
  char *keyword,*zone=NULL,*gridfile=NULL;
  char *meshfile=NULL,*belfile=NULL,*depthfile=NULL,*output=NULL,*poly=NULL;
  mesh_t mesh,refined,internal,external,*splitted,*final;
  double dmax=0;
  int *selected=0,*targeted=NULL;
  ostringstream report;
  int histogram[100];
  char *comment[2];
  int io_status;
  int InteriorsOnly=1;
  bool debug=false, aspect=false;
  float maxsize=5000, MinSharpAngle=10.0;
  int nitems;
  vector<int> watch;

  fct_echo( argc, argv);

  n=1;
  while (n < argc) {
    keyword=strdup(argv[n]);
    switch (keyword[0]) {
      case '-':
        if(strcmp(keyword,"--no_renum")==0) {
          renum=0;
          n++;
          continue;
          }
        if(strcmp(keyword,"--no_reshape")==0) {
          reshape=0;
          n++;
          continue;
          }
        if(strcmp(keyword,"--no-reshape")==0) {
          reshape=0;
          n++;
          continue;
          }
        if(strcmp(keyword,"--delaunay")==0) {
          delaunay=1;
          n++;
          continue;
          }
        if(strcmp(keyword,"--debug")==0) {
          debug=true;;
          n++;
          continue;
          }
        if(strcmp(keyword,"--aspect")==0) {
          aspect=true;;
          n++;
          continue;
          }
        if(strcmp(keyword,"--nghmax")==0) {
          nitems=sscanf(argv[n+1],"%d",&nghmax);
          if(nitems!=1) TRAP_ERR_EXIT(-1,"*** missing numerical value, expected after %s ***\n",keyword);
           n++;
          n++;
          continue;
          }
        if(strcmp(keyword,"--maxsize")==0) {
          nitems=sscanf(argv[n+1],"%f",&maxsize);
          if(nitems!=1) TRAP_ERR_EXIT(-1,"*** missing numerical value, expected after %s ***\n",keyword);
          n++;
          n++;
          continue;
          }
        if(strcmp(keyword,"--minsharp")==0) {
          nitems=sscanf(argv[n+1],"%f",&MinSharpAngle);
          if(nitems!=1) TRAP_ERR_EXIT(-1,"*** missing numerical value, expected after %s ***\n",keyword);
          n++;
          n++;
          continue;
          }
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
          original bel file, optional*/
          belfile= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'c' :
/*-----------------------------------------------------------------------------
          check channels*/
          channels=1;
          n++;
          break;

        case 'd' :
/*-----------------------------------------------------------------------------
          do not refine if size already smaller than dmax*/
          sscanf(argv[n+1],"%lf",&dmax);
          n++;
          n++;
          break;

        case 'g' :
/*-----------------------------------------------------------------------------
          check geometry*/
          geometry=1;
          n++;
          break;

        case 'p' :
/*-----------------------------------------------------------------------------
          limit action to polygons interior*/
          poly= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'o' :
          output= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 's' :
          safety=1;
          n++;
          break;

        case 'z' :
/*-----------------------------------------------------------------------------
          limit action to rectuganler zone interior*/
          zone= strdup(argv[n+1]);
          n++;
          n++;
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

  if(meshfile == NULL) {
    printf("*** Please specify mesh with -m ***\n");
    exit(-1);
    }

  printf("#################################################################\n");
  printf("load mesh and construct related tables: %s\n",meshfile);

  status=fe_readmesh(meshfile,MESH_FILE_FORMAT_TRIGRID,&mesh);
  if(status!=0) {
    printf("unable to read the original mesh in %s\n",meshfile);
    goto error;
    }
  status=fe_cleanvertices(&mesh, debug);
  
  report << "#Number of nodes      : " << mesh.nvtxs <<endl;

  selected=new int[mesh.nvtxs];

  status=fe_chkvertex_00( mesh, selected, histogram);
  if(status!=0) {
    printf("unsafe mesh description, abort\n");
    goto error;
    }
  status=fe_chkvertex_01( mesh, selected, histogram);

  delete[] selected;
  selected=0;

  status=fe_list(&mesh);
  if(status!=0) {
    printf("unable to build the element list from the original mesh\n");
    goto error;
    }

  report << "#Number of elements   : " << mesh.ntriangles <<endl;

  status=fe_bandwidth(&mesh);
  report << "#Half-bandwidth       : " << mesh.hbw <<endl;
//  cout << report.str() << flush;

  status=fe_edgetable(&mesh,0,0);
  if(status!=0) {
    printf("unable to build the edge list from the original mesh\n");
    goto error;
    }
  report << "#Number of edges      : " << mesh.nedges <<endl;
//  cout << report.str() << flush;

  if(safety==1) {
    if(poly!=NULL) {
      targeted=new int[mesh.nedges];
      status=fe_selectedges_01( mesh, poly, targeted);
      }
    status=fe_Echk_crossing ( mesh, targeted);
    if(targeted!=NULL )delete[] targeted;
    }

  
  printf("#################################################################\n");
  printf("reconstruct boundary codes and limits table\n");
  status=fe_codetable2(&mesh,0,1,0,1);
  if(status!=0) {
    printf("unable to rebuild the limits table and codes of the original mesh\n");
    goto error;
    }
  report << "#Number of limits     : " << mesh.nlimits <<endl;
  cout << report.str() << flush;
    
  if(delaunay==1) {
    printf("#################################################################\n");
    printf("Delaunay mesh improvements, affects boundary edges only (thus islands must be initialised by editing the mesh)\n");
    
    while(true) {
      int count=0;
      debug=true;
      count=clean_angles(mesh,(double) 120., false, debug);
//       status=fe_cleanvertices(&mesh, debug);
//       status=fe_list(&mesh);
//       status=fe_edgetable(&mesh,0,0);
      if (count==0) break;
      }
    final=&mesh;
    goto save;
     
//     status=fe_Mesh2DElasticityT(mesh, debug);
//     final=&mesh;
//     goto save;
    
/*------------------------------------------------------------------------------
    detect and remove oversized boundary edges */
    bool exterior_only=false;
    status=clean_delaunay01(mesh, maxsize, exterior_only, debug);
    
    exterior_only=true;
    selected=new int[mesh.nedges];
    while(true) {
      status=fe_codetable2(&mesh,0,1,1.0);
/*------------------------------------------------------------------------------
      detect and remove flat boundary elements */
      status=fe_Echk_angles(mesh, selected, 1.0, 0, exterior_only, debug);
      int count=0;
      for(n=0;n<mesh.nedges;n++) {
        if(selected[n]==0) continue;
        int n1=mesh.edges[n].extremity[0];
        int n2=mesh.edges[n].extremity[1];
        status=fe_disconnectvertices(mesh, n1, n2);
        count++;
        }
      if (count==0) break;
      status=fe_cleanvertices(&mesh, debug);
      status=fe_list(&mesh);
      status=fe_edgetable(&mesh,0,0);
      }
    delete[] selected;
    selected=0;
  
    status=fe_codetable2(&mesh,0,1,0,1);
    while(true) {
      int count=0;
      count=clean_delaunay02(mesh, exterior_only, debug);
      if (count==0) break;
      }
      
    selected=new int[mesh.nedges];
    for(k=0;k<100;k++) {
      status=fe_codetable2(&mesh,0,1,1.0);
/*------------------------------------------------------------------------------
      detect and remove sharp boundary elements */
      status=fe_Echk_angles(mesh, selected, MinSharpAngle, 1, exterior_only, debug);
      int count=0;
      for(n=0;n<mesh.nedges;n++) {
        if(selected[n]==0) continue;
        int n1=mesh.edges[n].extremity[0];
        int n2=mesh.edges[n].extremity[1];
        status=fe_disconnectvertices(mesh, n1, n2);
        count++;
        }
      if (count==0) break;
      status=fe_cleanvertices(&mesh, debug);
      status=fe_list(&mesh);
      status=fe_edgetable(&mesh,0,0);
      }
    delete[] selected;
    selected=0;
//     status=fe_codetable2(&mesh,0,1,0,1);
//     while(true) {
//       count=clean_delaunay02(mesh, debug);
//       if (count==0) break;
//       }
//     status=fe_codetable2(&mesh,0,1,1.0);
    while(true) {
      int count=0;
      count=clean_delaunay03(mesh, exterior_only, debug);
      if (count==0) break;
      status=fe_codetable2(&mesh,0,1,1.0);
      }
    status=fe_Tchk_angles( mesh, selected, (double) 5.0, debug);
    status=fe_codetable2(&mesh,0,1,1.0);
    final=&mesh;
    goto save;
    }
  
  if(belfile!=0) {
    status=fe_read_boundarycode(belfile, mesh, 0);
    status=fe_analyseBCflags(mesh, 0);
    }

  selected=new int[mesh.nedges];

  printf("#################################################################\n");
  printf("check for channel edges and bridging edges\n");
  status=fe_selectedges_02( mesh, selected, true);
  status=fe_selectedges_03( mesh, selected, true);
//  status=fe_selectedges_04( mesh, selected);

  if(geometry==1) {
    status=fe_Tchk_angles( mesh, selected, (double) 15.0, debug);
    status=fe_fix_RiverMouth(&mesh);
    status=fe_clean03(&mesh, (double) 0.1);
    status=fe_Tchk_AspectRatio( mesh, selected, (double) 0.75, NULL);
    }

  if(channels==1) {
    status  =fe_selectedges_02(mesh, selected, true);
    status  =fe_selectedges_03(mesh, selected, false);
    refined =fe_refine(mesh,selected,0);
    status  =fe_codetable2(&refined,1,1,0);
    final=fe_node2mesh(refined, debug);
    status=fe_savemesh("refined.nei",MESH_FILE_FORMAT_TRIGRID,*final);
    }
  else {
    final=&mesh;
    }

  printf("#################################################################\n");
  printf("cleave over-connected nodes\n");
  InteriorsOnly=1;
  
//   watch.push_back(6659);
//   watch.push_back(6660);
  
  status=fe_CureOverconnections(final, nghmax, InteriorsOnly, watch, debug);

  if(debug or status!=0) io_status=fe_savemesh("mesh-doctor.00.nei",MESH_FILE_FORMAT_TRIGRID,*final);
  if(status!=0) goto error;

  printf("#################################################################\n");
  printf("remove diamonds (interior nodes with 4 neighbours)\n");
  status=fe_supress_diamond(final, 4, 2, debug);

  if(debug) io_status=fe_savemesh("mesh-doctor.01.nei",MESH_FILE_FORMAT_TRIGRID,*final);

  if(aspect) {
    printf("#################################################################\n");
    printf("improve aspect ratio\n");
    status= clean_angles(*final,(double) 120., true, debug);
    if(debug) io_status=fe_savemesh("mesh-doctor.02.nei",MESH_FILE_FORMAT_TRIGRID,*final);
    if(status!=0) goto error;
    }
    
  flags=new float[final->nvtxs];
  for(n=0;n<final->nvtxs;n++) {
    flags[n]=0;
    }
  status=fe_Tchk_AspectRatio(*final, selected, (double) 0.75, flags);
  comment[0]= strdup("...");
  comment[1]= strdup("...");
  io_status=quoddy_saver1("flags.s2r", final->nvtxs,flags,comment);

  delete[] flags;

  if(renum==1) {
    printf("#################################################################\n");
    printf("optimize numbering\n");
    int  increment=MAX(10, 0.001*final->nvtxs);
    status=fe_reducebw(*final, increment, 100, (const char *) 0);
    }

  if(reshape==1) {
    printf("#################################################################\n");
    printf("reshape mesh\n");
    status=fe_reshapeall(*final,3);
    }

save:

  status=redep(mesh, *final);

  if(output!=0) {
    status=fe_savemesh(output,MESH_FILE_FORMAT_TRIGRID, *final);
    }
  else status=fe_savemesh("mesh-doctor.nei",MESH_FILE_FORMAT_TRIGRID,*final);
  
end: __OUT_BASE_LINE__("end of mesh-doctor ... \n");
  exit(0);

error:
  __OUT_BASE_LINE__("error detected, quit ... \n");
  exit(-1);
}
