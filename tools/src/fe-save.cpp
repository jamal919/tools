
#include <config.h>

#include <stdio.h>
#include <string.h>
#include <fstream>

#include "tools-structures.h"
#include "constants.h"

#include "fe.def"
#include "fe.h"

#include "map.h"
#include "netcdf-proto.h"
#include "poc-time.h"
#include "rutin.h"
#include "geo.h"
#include "swap.h"

  /**
   * \brief Convert from little or big endian to the other.
   * \param x Data to be converted.
   * \return The data in the new endian style.
   */
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> inline T littleBigEndian (T x) {

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
    /// \brief The new value.
    T tmp;
    
    /// \brief To interpret the value as an array of chars.
    unsigned char *toConvert = reinterpret_cast<unsigned char *>(&x);
    
    /// \brief The new value as an array of chars.
    unsigned char *converted = reinterpret_cast<unsigned char *>(&tmp);

/* *---------------------------------------------------------------------------
    Unix is big-endian, Linux is little-endian */

    for (size_t i = 0; i < sizeof(T); ++i)
      converted[i] = toConvert[sizeof(T) - i - 1];

    return tmp;
  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  inline unsigned readRecordSize (ifstream &file, /*const fs::path &path,*/ bool lbe) {

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  unsigned recordSize;

  file.read(reinterpret_cast<char*>(&recordSize), sizeof(unsigned));
  
//  if (!file) throw MeshTools::SerafinError (path.string());
  if (lbe) recordSize = littleBigEndian(recordSize);
  return recordSize;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  inline bool matchSize (ifstream &file, /*const fs::path &path,*/ unsigned recordSize, bool lbe)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  return (recordSize == readRecordSize(file, /*path,*/ lbe));
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> void write_telemac_BinaryBlock (ofstream &file, /*const fs::path &path,*/ bool lbe, T *vector, int length)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  const unsigned recordSize = (lbe)? littleBigEndian(static_cast<unsigned>(length * sizeof(T))):
                                                     static_cast<unsigned>(length * sizeof(T));
 
  file.write(reinterpret_cast<const char*>(&recordSize), sizeof(unsigned));
//  if (!file) throw SerafinWriteError (path.string());

  for (size_t i = 0; i < length; ++i) {

      const T n = (lbe)? littleBigEndian(static_cast<T>(vector[i])): static_cast<T>(vector[i]);

      file.write(reinterpret_cast<const char*>(&n), sizeof(T));
      
//      if (!file) throw SerafinWriteError (path.string());
      }

  file.write(reinterpret_cast<const char*>(&recordSize), sizeof(unsigned));
//  if (!file) throw SerafinWriteError (path.string());
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_savenodes_TGD(const char *filename,mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  TRIGRID format
  
  vulnerable to pinched boundaries ?
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  FILE *in=NULL;
  vertex_t *set=NULL;
  int i,k,l,ninteriors,nexteriors;

  in = fopen(filename, "w");
  if(in==0)
      return(-1);

  set= mesh.vertices;

  if(set == NULL) {
    return(-1);
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  count interior and exterior nodes
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  ninteriors=0;
  for (i=0; i<mesh.nvtxs; i++) {
    if(set[i].code==0) ninteriors++;
    }

  nexteriors=0;
  for (i=0; i<mesh.nvtxs; i++) {
    if(set[i].code>0) nexteriors++;
    }

  fprintf(in, "%d \n", ninteriors+nexteriors);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  first save exterior nodes by following mesh limits structure
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  fprintf(in, "%d \n", mesh.nlimits);

  for (l=0; l<mesh.nlimits; l++) {
    fprintf(in, "%d \n", mesh.limits[l].nvertex);
    for (k=0; k<mesh.limits[l].nvertex; k++) {
      i=mesh.limits[l].vertex[k];
      fprintf(in,"%lf %lf %f\n",set[i].lon,set[i].lat,set[i].h);
      }
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  then save interior nodes
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  fprintf(in, "%d \n", ninteriors);
  for (i=0; i<mesh.nvtxs; i++) {
    if(set[i].code==0) fprintf(in,"%lf %lf %f\n",set[i].lon,set[i].lat,set[i].h);
    }
  fclose(in);

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_savepoly(const char *filename,mesh_t mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  TRIANGLE format
  
  vulnerable to pinched boundaries ?
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  FILE *in=NULL;
  vertex_t *set=NULL;
  int i,j,k,l,nndes,nexteriors;
  int n1,n2,n3,count,interior;
  double dx[2],dy[2],dn[2],alpha;

  in = fopen(filename, "w");
  if(in==0)
      return(-1);

  set= mesh.vertices;
  nndes=mesh.nvtxs;

  if(set == NULL) {
    return(-1);
    }

  nexteriors=0;
  for (i=0; i<nndes; i++) {
    if(set[i].code!=0) nexteriors++;
    }

  fprintf(in, "%d 2 1 0\n", mesh.nvtxs);
  for (i=0; i<nndes; i++) {
    fprintf(in,"%d %lf %lf %f\n",i,set[i].lon,set[i].lat,set[i].h);
    }

  fprintf(in, "%d 0\n",nexteriors);

  count=0;
  for (l=0; l<mesh.nlimits; l++) {
    for (k=0; k<mesh.limits[l].nvertex-1; k++) {
      i=mesh.limits[l].vertex[k];
      j=mesh.limits[l].vertex[k+1];
      fprintf(in,"%d %d %d \n",count,i,j);
      count++;
      }
    i=mesh.limits[l].vertex[k];
    j=mesh.limits[l].vertex[0];
    fprintf(in,"%d %d %d \n",count,i,j);
    count++;
    }

//  fprintf(in, "%d \n", mesh.nlimits-1);

  for (l=1; l<mesh.nlimits; l++) {
    if(mesh.limits[l].nvertex==3) {
      mesh.limits[l].x=0.0;
      mesh.limits[l].y=0.0;
      continue;
      }
    for (k=0; k<mesh.limits[l].nvertex-1; k++) {
      n1=mesh.limits[l].vertex[k];
      n2=mesh.limits[l].vertex[k+1];
      n3=mesh.limits[l].vertex[k+2];
      dx[0] = set[n1].lon-set[n2].lon;
      dy[0] = set[n1].lat-set[n2].lat;
      dx[1] = set[n3].lon-set[n2].lon;
      dy[1] = set[n3].lat-set[n2].lat;
      dn[0] = sqrt(dx[0]*dx[0]+dy[0]*dy[0]);
      dn[1] = sqrt(dx[1]*dx[1]+dy[1]*dy[1]);
      dx[0] /= dn[0];
      dy[0] /= dn[0];
      dx[1] /= dn[1];
      dy[1] /= dn[1];
      alpha=asin(dx[0]*dy[1]-dy[0]*dx[1])*r2d;
      mesh.limits[l].x=(set[n1].lon+set[n2].lon+set[n3].lon)/3.0;
      mesh.limits[l].y=(set[n1].lat+set[n2].lat+set[n3].lat)/3.0;
      interior=0;
      if(fe_limit_interior (mesh,mesh.limits[l],mesh.limits[l].x,mesh.limits[l].y,&interior)==POINT_INTERIOR) {
        printf( "boundary #%d, legth #%d, pts #%d #%d, x=%lf y=%lf alpha=%lf\n", l, mesh.limits[l].nvertex, k,k+1,mesh.limits[l].x,mesh.limits[l].y,alpha);
        break;
        }
      }
    }

  count=0;
  for (l=1; l<mesh.nlimits; l++) {
    if(mesh.limits[l].nvertex>3) {
      count++;
      }
    }

  fprintf(in, "%d \n", count);
  for (l=1; l<mesh.nlimits; l++) {
    if(mesh.limits[l].nvertex>3) {
      fprintf(in, "%d %lf %lf\n", l,mesh.limits[l].x,mesh.limits[l].y);
      }
    }

  fclose(in);
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_savenodes(const char *filename, int fmt, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;

  switch(fmt) {
    case NODE_FILE_FORMAT_TRIGRID:
      status=fe_savenodes_TGD(filename,mesh);
      break;

    case NODE_FILE_FORMAT_TRIANGLE:
      status=fe_savepoly(filename,mesh);
      break;

    case NODE_FILE_FORMAT_QUODDY:
      status=-1;
      break;

//    case NODE_FILE_FORMAT_GOM:
      status=-1;
      break;

//    case NODE_FILE_FORMAT_NC3D:
      status=-1;
      break;

//    case NODE_FILE_FORMAT_MODULEF_P1:
      status=-1;
      break;

//    case NODE_FILE_FORMAT_MODULEF_P2:
      status=-1;
      break;

    case NODE_FILE_FORMAT_GMSH:
      status=-1;
      break;

    default:
      status=-1;
      break;

    }
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_savemesh_NGH(const char *filename,const mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  NEIGHBOUR format (TRIGRID)
  
  vulnerable to mono-triangle islands ?
 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  FILE *in=NULL;
  vertex_t *set=NULL;
  int i,j,nndes,maxn,status;
  frame_t frame;

  nndes=mesh.nvtxs;
  maxn=mesh.nnghm;
  set=mesh.vertices;
  
  maxn=0;
  for(int n=0;n<mesh.nvtxs;n++) {
    maxn=max(maxn, mesh.vertices[n].nngh);
    }

  if(nndes==0) return(-1);
  if(set == NULL) {
    printf("allocation error vertex set\n");
    return(-1);
    }

  if(maxn<=0){
    for (i=0; i<nndes; i++) {
      if(maxn<set[i].nngh)
        maxn=set[i].nngh;
      }
    }

  status=fe_minmax( mesh, frame);

/**----------------------------------------------------------------------------
  MPI unsafe */
  in = fopen(filename, "w");
  if(in==0)
      return(-1);

  fprintf(in, "%d \n",  nndes);
  fprintf(in, "%d \n",  maxn);
  fprintf(in, "%f %f %f %f \n", frame.xmin,frame.xmax,frame.ymin,frame.ymax);

  for (i=0; i<nndes; i++) {
    fprintf(in, "%d %12.8lf %12.8lf %d %f", i+1,set[i].lon,set[i].lat, set[i].code,set[i].h);
    for(j=0; j<set[i].nngh; j++) fprintf(in, " %d",set[i].ngh[j]+1);
    for(   ; j<maxn       ; j++) fprintf(in, " %d",0);
    fprintf(in, "\n");
    }
  fclose(in);

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_savemesh_METIS(const char *filename,mesh_t mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  METIS format
   
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  FILE *out=NULL;
  int i,k,nitems;

/*------------------------------------------------------------------------------
  neighbour list format*/

  out = fopen(filename, "w");
  if(out==0)
    return(-1);

  fprintf(out, "%d 1\n",  mesh.ntriangles);
  for (i=0; i<mesh.ntriangles; i++) {
    for (k=0; k<3; k++) {
      nitems=fprintf(out, " %d", mesh.triangles[i].vertex[k]+1);
      }
    fprintf(out, "\n");
    }
  fclose(out);

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_savemesh_SCHISM(const char *filename, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  FILE *in;
  int i,j,k,elt_count=0;
  int nitems;
  int type,ntags,tag1,tag2;
  vertex_t *set;

  in = fopen(filename, "w");
  if(in==0) {
    return(-1);
    }

  set= mesh.vertices;
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  mesh format : version, ascii, size of double
 
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  nitems=fprintf(in, "%s\n", "$MeshFormat");
  nitems=fprintf(in, "%s %d %d\n", "2",0,sizeof(double));
  nitems=fprintf(in, "%s\n", "$EndMeshFormat");

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  vertices
 
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  nitems=fprintf(in, "%s\n", "$Nodes");
  nitems=fprintf(in, "%12d\n", mesh.nvtxs);
  for (i=0; i<mesh.nvtxs; i++) {
    fprintf(in, "%10d %15.10lf %15.10lf %10.2f", i+1,set[i].lon,set[i].lat,set[i].h);
    fprintf(in, "\n");
    }
  nitems=fprintf(in, "%s\n", "$EndNodes");

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  elements
  
  pre-count the total number of "elements" : 
  
  triangle elements (type 2) + single point elements (type 15)
 
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  for (i=0; i<mesh.nlimits; i++) {
    for (j=0; j<mesh.limits[i].nvertex; j++){
      elt_count++;
      }
    }
  elt_count+=mesh.ntriangles;
  
  nitems=fprintf(in, "%s\n", "$Elements");
  nitems=fprintf(in, "%d\n ", elt_count);

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  single point elements (type 15)
 
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  elt_count=0;
  type=15;
  ntags=2;
  int buf_vtx_number=-1;
  
  for (i=0; i<mesh.nlimits; i++) {
    for (j=0; j<mesh.limits[i].nvertex; j++){
      elt_count++;
      nitems=fprintf(in,"%10d %10d %10d", elt_count, type, ntags);
      tag1=0;
      buf_vtx_number=mesh.limits[i].vertex[j];
      if (mesh.vertices[buf_vtx_number].code==5) {
        tag1=1;
        }
      tag2=mesh.limits[i].code -1 ;
      nitems=fprintf(in,"%10d %10d %10d\n", tag1, tag2, buf_vtx_number+1);
      }
    }
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  triangles (type 2)
 
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  type=2;
  ntags=3;
  for (i=0; i<mesh.ntriangles; i++) {
    elt_count++;
    nitems=fprintf(in, "%10d %10d %10d", elt_count,type,ntags);
    nitems=fprintf(in, "%10d %10d %10d", 0,i+1,5);
    for (k=0; k<3; k++) {
      nitems=fprintf(in, " %10d", mesh.triangles[i].vertex[k]+1);
      }
    fprintf(in, "\n");
    }
  nitems=fprintf(in, "%s\n", "$EndElements");

  fclose(in);

  return(0);
}
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_savemesh_GMSH_WW (const char *filename, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  GMSH WaveWatch (WW3) format
  
  see http://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format-version-2-_0028Legacy_0029
  
  $MeshFormat
  version-number file-type data-size
  $EndMeshFormat
  
  $PhysicalNames
  number-of-names
  physical-dimension physical-tag "physical-name"
  …
  $EndPhysicalNames
  
  $Nodes
  number-of-nodes
  node-number x-coord y-coord z-coord
  …
  $EndNodes
  
  $Elements
  number-of-elements
  elm-number elm-type number-of-tags < tag > … node-number-list
  …
  $EndElements
  
  $Periodic
  number-of-periodic-entities
  dimension entity-tag master-entity-tag
  number-of-nodes
  node-number master-node-number
  …
  $EndPeriodic
  
  $NodeData
  number-of-string-tags
  < "string-tag" >
  …
  number-of-real-tags
  < real-tag >
  …
  number-of-integer-tags
  < integer-tag >
  …
  node-number value …
  …
  $EndNodeData
  
  $ElementData
  number-of-string-tags
  < "string-tag" >
  …
  number-of-real-tags
  < real-tag >
  …
  number-of-integer-tags
  < integer-tag >
  …
  elm-number value …
  …
  $EndElementData
  
  $ElementNodeData
  number-of-string-tags
  < "string-tag" >
  …
  number-of-real-tags
  < real-tag >
  …
  number-of-integer-tags
  < integer-tag >
  …
  elm-number number-of-nodes-per-element value …
  …
  $EndElementNodeData
  
  $InterpolationScheme
  "name"
  number-of-element-topologies
  elm-topology
  number-of-interpolation-matrices
  num-rows num-columns value …
  …
  $EndInterpolationScheme
   
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  FILE *in;
  int i,j,k,elt_count=0;
  int nitems;
  int type,ntags,tag1,tag2;
  vertex_t *set;

  in = fopen(filename, "w");
  if(in==0) {
    return(-1);
    }

  set= mesh.vertices;
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  mesh format : version, ascii, size of double
 
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  nitems=fprintf(in, "%s\n", "$MeshFormat");
  nitems=fprintf(in, "%s %d %d\n", "2",0,sizeof(double));
  nitems=fprintf(in, "%s\n", "$EndMeshFormat");

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  vertices
 
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  nitems=fprintf(in, "%s\n", "$Nodes");
  nitems=fprintf(in, "%12d\n", mesh.nvtxs);
  for (i=0; i<mesh.nvtxs; i++) {
    fprintf(in, "%10d %15.10lf %15.10lf %10.2f", i+1,set[i].lon,set[i].lat,set[i].h);
    fprintf(in, "\n");
    }
  nitems=fprintf(in, "%s\n", "$EndNodes");

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  elements
  
  pre-count the total number of "elements" : 
  
  triangle elements (type 2) + single point elements (type 15)
 
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  for (i=0; i<mesh.nlimits; i++) {
    for (j=0; j<mesh.limits[i].nvertex; j++){
      elt_count++;
      }
    }
  elt_count+=mesh.ntriangles;
  
  nitems=fprintf(in, "%s\n", "$Elements");
  nitems=fprintf(in, "%d\n ", elt_count);

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  single point elements (type 15)
 
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  elt_count=0;
  type=15;
  ntags=2;
  int buf_vtx_number=-1;
  
  for (i=0; i<mesh.nlimits; i++) {
    for (j=0; j<mesh.limits[i].nvertex; j++){
      elt_count++;
      nitems=fprintf(in,"%10d %10d %10d", elt_count, type, ntags);
      tag1=0;
      buf_vtx_number=mesh.limits[i].vertex[j];
      if (mesh.vertices[buf_vtx_number].code==5) {
        tag1=1;
        }
      tag2=mesh.limits[i].code -1 ;
      nitems=fprintf(in,"%10d %10d %10d\n", tag1, tag2, buf_vtx_number+1);
      }
    }
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  triangles (type 2)
 
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  type=2;
  ntags=3;
  for (i=0; i<mesh.ntriangles; i++) {
    elt_count++;
    nitems=fprintf(in, "%10d %10d %10d", elt_count,type,ntags);
    nitems=fprintf(in, "%10d %10d %10d", 0,i+1,5);
    for (k=0; k<3; k++) {
      nitems=fprintf(in, " %10d", mesh.triangles[i].vertex[k]+1);
      }
    fprintf(in, "\n");
    }
  nitems=fprintf(in, "%s\n", "$EndElements");

  fclose(in);

  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_savemesh_GMSH (const char *filename, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  GMSH format
   
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  FILE *in=NULL;
  int i,k;
  int nitems;
  int type,ntags,tag;
  vertex_t *set=NULL;

  in = fopen(filename, "w");
  if(in==0) {
    return(-1);
    }

  set= mesh.vertices;

  nitems=fprintf(in, "%s\n", "$MeshFormat");
  nitems=fprintf(in, "%s %d %d\n", "2.0",0,sizeof(double));
  nitems=fprintf(in, "%s\n", "$EndMeshFormat");

  nitems=fprintf(in, "%s\n", "$Nodes");
  nitems=fprintf(in, "%d\n", mesh.nvtxs);
  for (i=0; i<mesh.nvtxs; i++) {
    fprintf(in, "%d %lf %lf %f", i+1,set[i].lon,set[i].lat,set[i].h);
    fprintf(in, "\n");
    }
  nitems=fprintf(in, "%s\n", "$EndNodes");

  nitems=fprintf(in, "%s\n", "$Elements");
  nitems=fprintf(in, "%d\n ", mesh.ntriangles);

  type=2;
  ntags=0;
  for (i=0; i<mesh.ntriangles; i++) {
    nitems=fprintf(in, "%d %d %d", i+1,type,ntags);
//     if(nitems != 3) {
//       printf("error in fe_readmesh\n");
//       return(-1);
//       }
    for (k=0; k<ntags; k++) {
      nitems=fprintf(in, "%d", tag);
      }
//    fprintf(in, "\n");
    for (k=0; k<3; k++) {
      nitems=fprintf(in, " %d", mesh.triangles[i].vertex[k]+1);
      }
    fprintf(in, "\n");
    }
  nitems=fprintf(in, "%s\n", "$EndElements");

  fclose(in);

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_savemesh_TELEMAC_ASCII (const char *rootname, mesh_t mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  TELEMAC format
   
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  FILE *in1=NULL,*in2=NULL;
  int i,m;
//   int dum;
  int fmt_ele,fmt_nod;
  int nitems,n1,n2,n3;
  char *elements=NULL,*nodes=NULL;
//   char  line[1000];

  elements=new char[strlen(rootname)+4];
  nodes=new char[strlen(rootname)+4];

  sprintf(elements,"%s.ele",rootname);
  sprintf(nodes,"%s.nod",rootname);

/*------------------------------------------------------------------------------
  element and node list format*/
  in1 = fopen(elements, "w");
  if(in1==0)
      return(-1);
  
//   fgets(line,999,in1);
//   nitems=sscanf(line, "%d %d %d %d", &dum, &n1, &n2, &n3);
//   switch(nitems) {
//     case 3:
//       fmt_ele=0;
//       break;
//     case 4:
//       fmt_ele=1;
//       break;
//     default:
//       return(-1);
//       break;
//     }
//   rewind (in1);
//
  in2 = fopen(nodes, "w");
  if(in2==0)
      return(-1);
//
//   fgets(line,999,in2);
//   if(strcmp(line,"NXY")==0) {
//     fmt_nod=0;
//     }
//   else if(strncmp(line,"XYZ",3)==0) {
//     fmt_nod=1;
//     }
//   else {
//     fmt_nod=0;
//     rewind (in2);
//     }
//

  for(m=0;m<mesh.ntriangles;m++) {
    n1=mesh.triangles[m].vertex[0]+1;
    n2=mesh.triangles[m].vertex[1]+1;
    n3=mesh.triangles[m].vertex[2]+1;
    switch(fmt_ele) {
      case 0:
        nitems=fprintf(in1, "%d %d %d", n1, n2, n3);
        if(nitems !=3) break;
        break;
      case 1:
        nitems=fprintf(in1, "%d %d %d %d", m+1, n1, n2, n3);
        if(nitems !=4) break;
        break;
      }
    }

  for (i=0; i<mesh.nvtxs; i++) {
    switch(fmt_nod) {
      case 0:
        nitems=fprintf(in2, "%d %lf %lf", i+1, mesh.vertices[i].lon,mesh.vertices[i].lat);
        if(nitems !=3) break;
        break;
      case 1:
        nitems=fprintf(in2, "%lf %lf %f", mesh.vertices[i].lon,mesh.vertices[i].lat,mesh.vertices[i].h);
        if(nitems !=3) break;
        break;
      }
    }

  fclose(in1);
  fclose(in2);

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_savemesh_TELEMAC_BINARY (const char *filename, mesh_t & mesh, const char *proj4, bool lbe)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  TELEMAC format
   
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  FILE *in=NULL;
  int i,k,status;
  char *title=NULL,*name=NULL;
  int nval;
  int *ibuf=NULL;
  float *rbuf=NULL;
  int *nbv=NULL;
  int cartesian=1;
  
  double *xx,*yy;
  
//  bool lbe=true;
  ofstream file (filename, std::ios_base::binary);
    
  in = fopen(filename, "r");
  if(in==0) {
    return(-1);
    }
/* *----------------------------------------------------------------------------
  Title */
  title=new char[80];
  strcpy(title,"TELEMAC TEST");
  write_telemac_BinaryBlock (file, lbe, title, 80);
  delete[] title;

/* *----------------------------------------------------------------------------
  Number of variables */
  nval=2;
  nbv=new int[nval];
  nbv[0]=1;
  nbv[1]=0;
  write_telemac_BinaryBlock (file, lbe, nbv, nval);

/* *----------------------------------------------------------------------------
  Variables names */
  for(i=0;i<nbv[0];i++) {
    name=new char[33];
//    for(int l=0; l<32; l++) name[l]=' ';
    strcpy(name,"FOND            m               ");
    write_telemac_BinaryBlock (file, lbe, name, 32);
    delete [] name;
    }

/* *----------------------------------------------------------------------------
  IPARAM */
  nval=10;
  ibuf=new int[nval];
  ibuf[0]=1;
  for(i=1;i<nval;i++) ibuf[i]=0;
  
  write_telemac_BinaryBlock (file, lbe, ibuf, nval);
  delete[] ibuf;

/* *----------------------------------------------------------------------------
  CARDINAL */
  nval=4;
  ibuf=new int[nval];
  ibuf[0]=mesh.ntriangles;
  ibuf[1]=mesh.nvtxs;
  ibuf[2]=3;
  ibuf[3]=1;
  write_telemac_BinaryBlock (file, lbe, ibuf, nval);
  delete[] ibuf;

/* *----------------------------------------------------------------------------
  CONNECTIVITY */
  nval=3*mesh.ntriangles;
  ibuf=new int[nval];
  for (i=0; i<mesh.ntriangles; i++) {
    for (k=0; k<3; k++) {
      ibuf[i*3+k]=mesh.triangles[i].vertex[k]+1;
      }
    }
  write_telemac_BinaryBlock (file, lbe, ibuf, nval);
  delete[] ibuf;

/* *----------------------------------------------------------------------------
  BOUNDARY CODES */
  nval=mesh.nvtxs;
  ibuf=new int[nval];
/**----------------------------------------------------------------------------
  temporary */
//   for (i=0; i<mesh.nvtxs; i++) {
//     ibuf[i]=0;
//     }
//
//   int count=1;
//   for(int l=0;l<mesh.nlimits;l++) {
//     for(k=0;k<mesh.limits[l].nvertex;k++) {
//       ibuf[mesh.limits[l].vertex[k]]=count;
//       count++;
//       }
//     }
    
/**----------------------------------------------------------------------------
  temporary */
  for (i=0; i<mesh.nvtxs; i++) {
    ibuf[i]=mesh.vertices[i].code;
    }

  write_telemac_BinaryBlock (file, lbe, ibuf, nval);
  delete[] ibuf;

  cartesian=1;
  xx=new double[mesh.nvtxs];
  yy=new double[mesh.nvtxs];
    
/* *----------------------------------------------------------------------------
  Coordinates and projection */
  if(cartesian==1) {
//     projPJ ref=NULL;
/* *----------------------------------------------------------------------------
    Transverse Mercator */
//     const char *parms = "proj=utm +zone=32 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs";
    
//     int utm_zone;
//     char proj4_options[128], hemisphere[32];
//     if(p0<0) sprintf(hemisphere,"+south");
//     else sprintf(hemisphere,"");
//     utm_zone=(t0+180)/6+1;
//     sprintf(proj4_options,"proj=utm +zone=%d +ellps=WGS84 +datum=WGS84 +units=m +no_defs %s",utm_zone,hemisphere);
    
//     if ( ! (ref = pj_init_plus(parms)) ) {
//       TRAP_ERR_EXIT(1, "Projection initialization failed\n");
//       }
    for(i=0;i<mesh.nvtxs;i++) {
      xx[i]=mesh.vertices[i].lon;
      yy[i]=mesh.vertices[i].lat;
      }
    
    if(proj4!=0) {
      status=geo_to_projection (proj4, xx, yy, mesh.nvtxs);
      }
    else {
//       double *tt=new double[mesh.nvtxs];
//       double *pp=new double[mesh.nvtxs];
//       CONV_LAMBERT_DEGDEC_TAB(mesh.nvtxs,xx, yy, tt, pp );
//       for(i=0;i<mesh.nvtxs;i++) {
//         mesh.vertices[i].lon=tt[i];
//         mesh.vertices[i].lat=pp[i];
//         }
//       delete[] tt;
//       delete[] pp;
      }
    }
  else {
    for(i=0;i<mesh.nvtxs;i++) {
      xx[i]=mesh.vertices[i].lon;
      yy[i]=mesh.vertices[i].lat;
      }
    }

/* *----------------------------------------------------------------------------
  LONGITUDES */
  nval=mesh.nvtxs;
  rbuf=new float[nval];
  for (i=0; i<mesh.nvtxs; i++) {
    rbuf[i]=xx[i];
    }
  write_telemac_BinaryBlock (file, lbe, rbuf, nval);
  delete[] rbuf;

/* *----------------------------------------------------------------------------
  LATITUDES */
  nval=mesh.nvtxs;
  rbuf=new float[nval];
  for (i=0; i<mesh.nvtxs; i++) {
    rbuf[i]=yy[i];
    }
  write_telemac_BinaryBlock (file, lbe, rbuf, nval);
  delete[] rbuf;

  delete[] xx;
  delete[] yy;
  
/* *----------------------------------------------------------------------------
  TEMPS */
  nval=1;
  rbuf=new float[nval];
  for (i=0; i<nval; i++) {
    rbuf[i]=0.;
    }
  write_telemac_BinaryBlock (file, lbe, rbuf, nval);
  delete[] rbuf;

/* *----------------------------------------------------------------------------
  DEPTHS */
  nval=mesh.nvtxs;
  rbuf=new float[nval];
  for (i=0; i<mesh.nvtxs; i++) {
    rbuf[i]=mesh.vertices[i].h;
    }
  write_telemac_BinaryBlock (file, lbe, rbuf, nval);
  delete[] rbuf;
  
  file.close();
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_find_format(const string &filename)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  if(strrncasecmp(filename,".nei")==0)
    return MESH_FILE_FORMAT_TRIGRID;
  
  if(strrncasecmp(filename,".nc")==0)
    return MESH_FILE_FORMAT_NC2D;
  
  return MESH_FILE_FORMAT_UNKNOWN;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_savemesh(const char *filename, int fmt, mesh_t & mesh, const char *proj4)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  char *rootname=NULL;
  
  if(fmt==MESH_FILE_FORMAT_UNKNOWN)
    fmt=fe_find_format(filename);

  switch(fmt) {
    case MESH_FILE_FORMAT_TRIGRID:
/*------------------------------------------------------------------------------
      neighbour list format*/
      status=fe_savemesh_NGH(filename,mesh);
      break;

    case MESH_FILE_FORMAT_TRIANGLE:
      status=-1;
      break;

    case MESH_FILE_FORMAT_QUODDY:
      status=-1;
      break;

    case MESH_FILE_FORMAT_GOM:
      status=-1;
      break;

    case MESH_FILE_FORMAT_NC2D:
      status=fe_savemeshNC2D(filename,mesh,1);
      break;

    case MESH_FILE_FORMAT_NC3D:
      status=fe_savemeshNC3D(filename,mesh,1);
      break;

    case MESH_FILE_FORMAT_MODULEF_P1:
      status=-1;
      break;

    case MESH_FILE_FORMAT_MODULEF_P2:
      status=-1;
      break;

    case MESH_FILE_FORMAT_GMSH:
      status=fe_savemesh_GMSH(filename,mesh);
      break;

    case MESH_FILE_FORMAT_GMSH_WW:
      status=fe_savemesh_GMSH_WW(filename,mesh);
      break;
      
    case MESH_FILE_FORMAT_SCHISM:
      status=fe_savemesh_SCHISM(filename,mesh);
      break;
      
    case MESH_FILE_FORMAT_METIS:
      status=fe_savemesh_METIS(filename,mesh);
      break;

    case MESH_FILE_FORMAT_TELEMAC_ASCII:
      rootname=strdup(filename);
      rootname[strlen(filename)-4]=0;
      status=fe_savemesh_TELEMAC_ASCII(rootname,mesh);
      break;

    case MESH_FILE_FORMAT_TELEMAC_BINARY:
      status=fe_savemesh_TELEMAC_BINARY(filename,mesh,proj4, (bool) false);
      break;

    case MESH_FILE_FORMAT_TELEMAC_SWAPPED:
      status=fe_savemesh_TELEMAC_BINARY(filename,mesh,proj4, (bool) true);
      break;

    default:
      printf("%s: undefined mesh format %d for %s\n",__func__,fmt,filename);
      status=-1;
      break;

    }
  return(status);
}
