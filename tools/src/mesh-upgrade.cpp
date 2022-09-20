

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


#include "tools-structures.h"

#include "geo.h"
#include "fe.h"
#include "map.h"
#include "polygones.h"
#include "grd.h"
#include "functions.h"

#include "zapper.h"     /*  rutin.h contains common utility routines  */

#include "statistic.h"


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  mesh_t refine(mesh_t mesh,int *selected)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int   k,l,m,i,j,j1,j2,n1,n2,n,status;
  int   nedges,chk=0,count=0,tmp=0;
  int   nelts,nndes;
  int   interior=0,boundary=0,weird=0;
  int   ring[4]={0,1,2,0};
  int   *ncells=NULL,ncellsmax=0,**cells=NULL;
  edge_t *edges=NULL;
  triangle_t *elt=NULL,*ptr=NULL;
  double t1,t2,p1,p2,d;

  mesh_t work;

  elt=mesh.triangles;

  nelts =mesh.ntriangles;
  nndes =mesh.nvtxs;
  nedges=mesh.nedges;
  edges =mesh.edges;

  printf ("number of nodes (original mesh) %d\n",nndes);

  work.type=mesh.type;

  work.vertices=new vertex_t[nedges+nndes];
  for (n=0; n<nedges+nndes; n++) work.vertices[n].null_value();

  for (n=0;n<nndes;n++) {
    work.vertices[n].lon=mesh.vertices[n].lon;
    work.vertices[n].lat=mesh.vertices[n].lat;
    work.vertices[n].h=mesh.vertices[n].h;
    work.vertices[n].code=mesh.vertices[n].code;
    if(work.vertices[n].code==0) {
      work.vertices[n].nngh=0;
      }
    else {
      work.vertices[n].nngh=2;
      exitIfNull(
        work.vertices[n].ngh=new int[work.vertices[n].nngh]
        );
      work.vertices[n].ngh[0]=mesh.vertices[n].ngh[0];
      work.vertices[n].ngh[1]=mesh.vertices[n].ngh[mesh.vertices[n].nngh-1];
//       work.vertices[n].nngh=0;
//       for(k=0;k<mesh.vertices[n].nngh;k++) {
// 	m=mesh.vertices[n].ngh[k];
// 	if(mesh.vertices[m].code!=0) work.vertices[n].nngh++;
//         }
//       work.vertices[n].nngh=0;
//       for(k=0;k<mesh.vertices[n].nngh;k++) {
// 	m=mesh.vertices[n].ngh[k];
// 	if(mesh.vertices[m].code!=0) {
// 	  work.vertices[n].ngh[work.vertices[n].nngh]=m;
// 	  work.vertices[n].nngh++;
//           }
//         }
      }
    }

  count=nndes;
  goto interior;

/* *-----------------------------------------------------------------------------
  boundary refinement */
  for (n=0;n<nedges;n++) {
    if(edges[n].nshared==2) continue;
    n1=edges[n].extremity[0];
    n2=edges[n].extremity[1];
    if (mesh.vertices[n1].code!=mesh.vertices[n2].code) continue;
/*     if (mesh.vertices[n1].code!=20) continue; */
    t1=mesh.vertices[n1].lon;
    t2=mesh.vertices[n2].lon;
    p1=mesh.vertices[n1].lat;
    p2=mesh.vertices[n2].lat;
    if(selected[n]) {
/*-----------------------------------------------------------------------------
      create a new node at the middle of the edge */
      if(work.type==0) t2=geo_recale(t2,t1,(double) 180.0);
      work.vertices[count].lon=0.5*(t1+t2);
      work.vertices[count].lat=0.5*(p1+p2);
      work.vertices[count].h=0.5*(mesh.vertices[n1].h+mesh.vertices[n2].h);
      work.vertices[count].nngh=2;
      work.vertices[count].code=mesh.vertices[n1].code;
      exitIfNull(
        work.vertices[count].ngh=new int[2]
        );
      work.vertices[count].nngh=2;
/*-----------------------------------------------------------------------------
      connect */
      if(work.vertices[n1].ngh[1]==n2) {
//        work.vertices[n2].ngh[0]=count+1;
//        work.vertices[n1].ngh[1]=count+1;
        work.vertices[n2].ngh[0]=count;
        work.vertices[n1].ngh[1]=count;
        work.vertices[count].ngh[0]=n1;
        work.vertices[count].ngh[1]=n2;
        }
      else if(work.vertices[n1].ngh[0]==n2) {
 //       work.vertices[n1].ngh[0]=count+1;
 //       work.vertices[n2].ngh[1]=count+1;
        work.vertices[n1].ngh[0]=count;
        work.vertices[n2].ngh[1]=count;
        work.vertices[count].ngh[0]=n2;
        work.vertices[count].ngh[1]=n1;
        }
      else {
        printf("troubles...\n");
        }
      count++;
      }
    }

/* *-----------------------------------------------------------------------------
  interior refinement */
interior:
  for (n=0;n<nedges;n++) {
    if(edges[n].nshared!=2) continue;
    n1=edges[n].extremity[0];
    n2=edges[n].extremity[1];
    t1=mesh.vertices[n1].lon;
    t2=mesh.vertices[n2].lon;
    p1=mesh.vertices[n1].lat;
    p2=mesh.vertices[n2].lat;
    if(selected[n]) {
/*-----------------------------------------------------------------------------
      create a new node at the middle of the edge */
      if(mesh.type==0) t2=geo_recale(t2,t1,(double) 180.0);
      work.vertices[count].lon=0.5*(t1+t2);
      work.vertices[count].lat=0.5*(p1+p2);
      work.vertices[count].h=0.5*(mesh.vertices[n1].h+mesh.vertices[n2].h);
      work.vertices[count].code=0;
      work.vertices[count].nngh=0;
      count++;
      }
    }

 end:
  work.nvtxs=count;
  work.ntriangles=0;
  work.nlimits=mesh.nlimits;
  work.limits=new limit_t[work.nlimits];
  for (l=0;l<work.nlimits;l++) {
    work.limits[l].nvertex=mesh.limits[l].nvertex;
    work.limits[l].vertex=new int[work.limits[l].nvertex];
    for (k=0;k<work.limits[l].nvertex;k++) {
      work.limits[l].vertex[k]=mesh.limits[l].vertex[k];
      }
    }
  printf ("number of nodes (refined mesh) %d\n",count);
  return(work);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_checkpartition(mesh_t mesh,int *partition, int npartition)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *----------------------------------------------------------------------------
  check mesh connexity */
{
  int i,k,l,p,q,m,n,n1,n2,nndes,status;
  int *used,previous;
  
  used=new int[mesh.nvtxs];
  nndes=0;
  for (n=0;n<mesh.nvtxs;n++) used[n]=-1;
  
  used[0]=1;
  nndes=1;
  previous=0;
  while(nndes>previous) {
    previous=nndes;
    for (n=0;n<mesh.nvtxs;n++) {
      if(used[n]==-1) continue;
      for(i=0;i<mesh.vertices[n].nngh;i++) {
        m=mesh.vertices[n].ngh[i];
        if(used[m]==-1) {
          used[m]=nndes;
          nndes++;
          }
        }
      }
    }
    
  if( nndes!=mesh.nvtxs) {
    printf("partition not connex\n");
    for (m=0;m<mesh.ntriangles;m++) {
      for(i=0;i<3;i++) {
        n=mesh.triangles[m].vertex[i];
        if(used[n]==-1) partition[mesh.triangles[m].ancestor]=npartition;
        break;
        }
      }
    delete[] used;
    return(-1);
    }
  
  delete[] used;
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  mesh_t fe_global(mesh_t mesh, int global, int npartition)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,l,p,q,n,n1,n2,status;
  mesh_t refurbished,test,refined;
  mesh_t *splitted,submesh[2];
  int *selected,*targeted, nselected;
  double d,t1,t2,p1,p2;
  char *cmdline,*partition_file;
  int  *partition;
  int **connected,*sequence;
  bool debug=false;

  cmdline=new char[64];
  partition_file=new char[256];
  
  status=fe_savemesh("mesh.metis",MESH_FILE_FORMAT_METIS,mesh);

  sprintf(cmdline,"%s %s %d","/home/softs/metis-4.0/partnmesh","mesh.metis",npartition);
  status=system(cmdline);

  sprintf(partition_file,"mesh.metis.epart.%d",npartition);
  
  partition=fe_read_epartition(partition_file,mesh);
redo:
  splitted=fe_partition_01(mesh,partition,npartition);

  for(p=0;p<npartition;p++) {
    printf("checking partition %d\n",p);
    status=fe_checkpartition(splitted[p],partition,npartition);
    if (status!=0) {
      for(q=0;q<npartition;q++) splitted[q].destroy();
      npartition++;
      goto redo;
      }
    }

  sequence=new int[npartition];
  connected=new int*[npartition];
  for(p=0;p<npartition;p++) {
    connected[p]=new int[npartition];
    }

  for(p=0;p<npartition;p++) {
    for(q=p+1;q<npartition;q++) {
      connected[p][q]=0;
      connected[q][p]=0;
      for(k=0;k<splitted[p].limits[0].nvertex;k++) {
        n1=splitted[p].limits[0].vertex[k];
        for(l=0;l<splitted[q].limits[0].nvertex;l++) {
          n2=splitted[q].limits[0].vertex[l];
          if(splitted[p].vertices[n1].ancestor==splitted[q].vertices[n2].ancestor) {
            connected[p][q]=1;
            connected[q][p]=1;
            printf("partition %d is connected with %d\n",p,q);
            break;
            }
          }
        if(connected[p][q]==1) break;
        }
      }
    }

  for(p=0;p<npartition;p++) {
//    nselected=fe_selectedges_04( mesh, criteria.maxsize, selected);
    selected=new int[splitted[p].nedges];
    for (n=0;n<splitted[p].nedges;n++) {
      selected[n]=1;
      }
    nselected=fe_selectedges_04(splitted[p], 100., selected);
    if(nselected==0) continue;
    refined=refine(splitted[p], selected);
//    status =fe_codetable1(&refined,1);
    splitted[p].destroy();
    splitted[p]=*(fe_node2mesh(refined, debug));
//    status=fe_savemesh("temporary.nei",MESH_FILE_FORMAT_TRIGRID,*final);
    delete[] selected;
    }


  sequence[0]=0;
  n=0;
  while(n<npartition-1) {
    p=sequence[n];
    for(q=0;q<npartition;q++) {
      if(p==q) continue;
      if(connected[p][q]==1) {
        sequence[n+1]=q;
        for(k=0;k<npartition;k++) connected[q][k]=MAX(connected[p][k],connected[q][k]);
        for(k=0;k<npartition;k++) connected[k][p]=0;
        n++;
        break;
        }
      }
    }
  
//   partition=fe_read_epartition("mesh.metis.epart.4",mesh);
//   splitted=fe_partition_03(mesh,partition,npartition);
  
//   selected=new int[mesh.nedges];
//   targeted=new int[mesh.nedges];
//
//   for (n=0;n<mesh.nedges;n++) {
//     selected[n]=0;
//     }
//   for (n=0;n<mesh.nedges;n++) {
//     n1=mesh.edges[n].extremity[0];
//     n2=mesh.edges[n].extremity[1];
//     t1=mesh.vertices[n1].lon;
//     t2=mesh.vertices[n2].lon;
//     p1=mesh.vertices[n1].lat;
//     p2=mesh.vertices[n2].lat;
//     if(p1>60.0) {
//       selected[n]=1;
//       }
//     if(p2>60.0) {
//       selected[n]=1;
//       }
//     }
//
//   splitted=fe_split(mesh,selected,targeted);
//   if(splitted==NULL) {
//     printf("unable to split the original mesh\n");
//     goto error;
//     }

  refurbished.type=SPHERICAL;
  submesh[0]=splitted[0];
  for(p=0;p<npartition-1;p++) {
    submesh[1]=splitted[sequence[p+1]];
    printf("merging with partition %d\n",sequence[p+1]);
    refurbished=fe_merge(submesh,(double) 0.1, (const char*) 0);
    submesh[0].destroy();
    submesh[1].destroy();
    submesh[0]=refurbished;
    }
    
  return(test);
error:
  return(test);

}



/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int nndes,status,*elts,fmt,fefmt;
  int i,j,k,l,L,m,n;
  int option,channels=0,resample=0;
  FILE *file;
  FILE *out;
  char *keyword,*zone=NULL,*gridfile=NULL;
  char *meshfile=NULL,*depthfile=NULL,*output=NULL,*poly=NULL,*critere=NULL;
  char *meshsize=NULL;
  mesh_t mesh,refined,internal,external,*splitted,final;
  double dmax=0,cmax=0,maxrate=0.75;
  int *selected,*targeted,nselected;
  criteria_t criteria;
  plg_t *polygones;
  int npolygones;
  grid_t grid;
  float *density,mask;
  bool debug=false, exact=false;
  int size=1;

  fct_echo(argc,argv);

  n=1;
  while (n < argc) {
    keyword=strdup(argv[n]);
/* *----------------------------------------------------------------------
    boundary re-sampling*/
    if(strcmp(keyword,"--resample")==0) {
      resample=1;
      n++;
      continue;
      }
    if(strcmp(keyword,"--debug")==0) {
      debug=true;
      n++;
      continue;
      }
    if(strcmp(keyword,"--exact")==0) {
      exact=true;
      n++;
      continue;
      }
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

        case 'd' :
          sscanf(argv[n+1],"%lf",&dmax);
          n++;
          n++;
          break;

        case 'r' :
          sscanf(argv[n+1],"%lf",&maxrate);
          n++;
          n++;
          break;

        case 'e' :
          sscanf(argv[n+1],"%lf",&cmax);
          n++;
          n++;
          break;

        case 'c' :
          critere= strdup(argv[n+1]);
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

/* *----------------------------------------------------------------------
  load and initialize criteria structure*/
  printf("#################################################################\n");
  printf("load meshing criteria : %s\n",critere);
  status=fe_load_criteria(critere, &criteria);
  if(status!=0) {
    printf("unable to load meshing criteria, abort...\n");
    goto error;
    }
  
/* *----------------------------------------------------------------------
  boundary re-sampling*/
  criteria.resample_obsolete=resample;

/* *----------------------------------------------------------------------
  element size increasing/decreasing max rate*/
//  criteria.maxrate=maxrate;

/**-----------------------------------------------------------------------------
  echo criteria file */
  status= fe_save_criteria("mesh-upgrade-echo.crt", criteria);
  
/* *----------------------------------------------------------------------
  load mesh*/
  printf("#################################################################\n");
  printf("load mesh ant initialize related tables\n");
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
  status= fe_edgetable(&mesh,0,0);
  if(status!=0) {
    printf("unable to build the edge list from the original mesh\n");
    goto error;
    }
    
  status=fe_vertex_crosstables02(&mesh);
  
  status= fe_codetable2(&mesh,0,1,0);
  if(status!=0) {
    printf("unable to rebuild the limits table and codes of the original mesh\n");
    goto error;
    }
    
  if(exact) {
    status=fe_presplit(mesh, poly, 0.5, "");
    if(status!=0) goto error;
    size=3;
    }
  else {
    size=1;
    }

  selected=new int[mesh.nedges];
  targeted=new int[mesh.nedges];

/* *----------------------------------------------------------------------
  refine global mesh */
//  refined=fe_global(mesh, 1, 30);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  select edges to be refined 
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  for (n=0;n<mesh.nedges;n++) {
    selected[n]=1;
    }

if(poly!=0) {
  printf("#################################################################\n");
    printf("select edges from polygons: %s\n",poly);
    nselected=fe_selectedges_01(mesh, poly, selected);
    }
  else {
    nselected=mesh.nedges;
    }
  if(nselected==0) goto end;


/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  extract submesh from selection to allow decent cartesian reshape 
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  printf("#################################################################\n");
  printf("split mesh\n");
  splitted=fe_split(mesh, selected, targeted, 0, size, (string) "mesh-upgrade", debug);
  if(splitted==NULL) {
    printf("unable to split the original mesh\n");
    goto error;
    }
    
  printf("#################################################################\n");
  printf("check split connexity\n");
  status=fe_connex(splitted[1]);
  if(status==-1) {
    printf("internal mesh not connex\n");
    goto error;
    }
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  transform mesh limits into polygons
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  status=fe_limits2poly(splitted[1], &polygones, &npolygones,true);
 
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  mesh size mapping
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  criteria.resample_obsolete=( (criteria.resample_openlimits==1) && (criteria.resample_rigidlimits==1) );
   
  if(meshsize==0) {
    printf("#################################################################\n");
    printf("compute mesh density\n");
    meshsize=strdup("mesh-size.nc");
    status=fe_ComputeMeshsize(meshsize, criteria, &grid, &density, &mask, polygones, npolygones,1,false);
    if(status!=0) {
      goto error;
      }
    }
  else {
    printf("#################################################################\n");
    printf("load mesh density from %s\n",meshsize);
    status=fe_ReloadMeshsize(meshsize, &grid, &density, &mask, polygones, npolygones);
    if(status!=0) {
      goto error;
      }
    status= defproj(&grid, polygones, npolygones);
    }

//  debug=true;  
  if(criteria.resample_rigidlimits) {
    printf("#################################################################\n");
    printf("rigid boundary adjustment\n");
/* *---------------------------------------------------------------------------
    re-sample polygons by following meshing criteria */
    for(int s=0; s<  npolygones; s++) {
      if(polygones[s].npt==0) {
        printf("empty polygon %d\n",s);
        continue;
        }
      if(polygones[s].flag!=0) {
        plg_t q=plg_resample_rigid_obsolete(polygones[s], grid.projection, criteria.shelf_minsize, 1);
        polygones[s].duplicate(q);
        }
      }
    status=plg_spherical(grid.projection, polygones, npolygones);
    status=plg_checkAutoSecant(polygones, npolygones, PLG_SPHERICAL);
    status=plg_save("mesh-generator-rigid-adjusted.plg", PLG_FORMAT_SCAN,  polygones, npolygones);
    status=plg_load_scan("mesh-generator-rigid-adjusted.plg", &polygones, &npolygones);
    grid.free();
    delete[] density;
    density=0;
    criteria.resample_obsolete=0;
    criteria.resample_rigidlimits=0;
    status=("criteria-adjusted.nc", criteria, &grid, &density, &mask, polygones, npolygones,1, false);
/* *---------------------------------------------------------------------------
    next step necessary because of projection/spherical issue */
    for(int k=0;k<20;k++) status=map_persistence(grid, density, mask, 0.);
    }
 
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  mesh generation
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  printf("#################################################################\n");
  printf("create mesh nodes\n");
  refined=fe_nodit(criteria, grid, density, mask, polygones, npolygones, debug);
    
  printf("#################################################################\n");
  printf("merge meshes\n");
  splitted[1]=refined;
  final=fe_merge(splitted,(double) 0.01, (char *) 0);
  
  if(debug) 
    status= fe_savemesh("mesh-upgrade-no-reshape.nei",MESH_FILE_FORMAT_TRIGRID, final);
  
  printf("#################################################################\n");
  printf("reshape refined mesh\n");
  status= fe_reshapeall(final,3);

  if(output!=0) {
    status= fe_savemesh(output,MESH_FILE_FORMAT_TRIGRID, final);
    }
  else
    status= fe_savemesh("mesh-upgrade.nei",MESH_FILE_FORMAT_TRIGRID, final);
 
  printf("#################################################################\n");
  printf("summary:\n");
  printf("before upgrade: %d vertices %d elements\n",mesh.nvtxs,mesh.ntriangles);
  printf("after  upgrade: %d vertices %d elements\n",final.nvtxs,final.ntriangles);
  
end: __OUT_BASE_LINE__("end of mesh-upgrade ... \n");
  exit(0);

error:
  __OUT_BASE_LINE__("error detected, quit ... \n");
  exit(-1);
}
