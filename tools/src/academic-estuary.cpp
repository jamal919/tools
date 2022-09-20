

/*******************************************************************************

  T-UGO tools, 2006-2009

  Unstructured Ocean Grid initiative

Contributors:

  Florent Lyard      LEGOS/CNRS, Toulouse, France
  Yoann Le Bars      LEGOS, Toulouse, France (PhD)
  Laurent Roblou     LEGOS, Toulouse, France
  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
  Frédéric Dupont    Université de Laval à Québec, Canada

E-mail: florent.lyard@legos.obs-mip.fr

*******************************************************************************/

/*******************************************************************************

  2D mesh generator (academic, simple geometry meshes)

  F. Lyard, 2009, CNRS/LEGOS, Toulouse, France

  e-mail: florent.lyard@legos.obs-mip.fr

*******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <list>

using namespace std;

#include "tools-structures.h"

#include "functions.h"
#include "geo.h"
#include "fe.h"
#include "archive.h"
#include "map.h"
#include "polygones.h"
#include "grd.h"
#include "list.h"
#include "maths.h"
#include "sym-io.h"
#include "map.def"
#include "mass.h"
#include "topo.h"
#include "filter.h"
#include "statistic.h"

#include "zapper.h"     /*  rutin.h contains common utility routines  */

extern  void fe_LGP0_to_LGP1(double *in, double *out, mesh_t mesh);
#include "statistic.h"

#include "academic.h"

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_spherical(mesh_t & mesh, projPJ projection)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=0;
  for (int n=0;n<mesh.nvtxs;n++) {
    double x,y;
    x=mesh.vertices[n].lon;
    y=mesh.vertices[n].lat;
    projection_to_geo(projection,&(mesh.vertices[n].lat),&(mesh.vertices[n].lon),x,y);
    }
  
  mesh.type=SPHERICAL;
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_rotate_cartesian(plg_t & p, point2D_t center, double angle)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int status=0;
  double x,y;
  
  matrix2x2_t rotation=math_rotation2D_init(angle);
  
  for(int i=0;i<p.npt;i++) {
    vector2D_t u(p.x[i]-center.x, p.y[i]-center.y), v;
    v=math_rotation2D(rotation, u);
    p.x[i]=v.x+center.x;
    p.y[i]=v.y+center.y;
    }
    
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_rotate_cartesian(vector<plg_t> & polygons, point2D_t center, double angle)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int status=0;
  double x,y;
  
  matrix2x2_t rotation=math_rotation2D_init(angle);
  
  for(int i=0;i<polygons.size();i++) {
    status=plg_rotate_cartesian(polygons[i], center, angle);
    }
    
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_translate_cartesian(plg_t & p, vector2D_t & u)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int status=0;
  double x,y;
    
  for(int i=0;i<p.npt;i++) {
    p.x[i]+=u.x;
    p.y[i]+=u.y;
    }
    
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_translate_cartesian(vector<plg_t> & polygons, vector2D_t & u)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int status=0;
  double x,y;
  
  
  for(int i=0;i<polygons.size();i++) {
    status=plg_translate_cartesian(polygons[i], u);
    }
    
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int process_01(river_t & backbone, plg_t & p, double L0, double H0, projPJ projection)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  process backbone central polygon 
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int status=0;
  size_t size,count;
  double x0,y0;
  
  backbone.parse();
  size=backbone.n+1;
  
  p.init(size, PLG_INIT_SEPARATE);
 
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  curved shape given by:
  
    k=2.*M_PI*2.0/L0   ->    y=H0*(1.0-cos(k*l))
    
  for continuity allowance, rather use dy definition:
      
    dy/dl=H0*k*sin(k*l)
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
    
  p.x[0]=0.0;
  p.y[0]=0.0;
  
  p.z=new double[p.npt];
  
  double G=M_PI/4.0;
  double lamda,l=0;
  
  count=1;
  for(int i=0;i<backbone.parts.size();i++) {
    double k=backbone.parts[i].wavenumber;
    if(not isnan(backbone.parts[i].phase)) G=backbone.parts[i].phase;
    double atom=backbone.parts[i].resolution;
    for(int j=0; j<backbone.parts[i].n;j++) {
      double dx,dy,y;
      int n=count;
      l=j*atom;
      dx=atom;
      p.x[n]=p.x[n-1]+dx;
      dy=H0*k*sin(k*l+G)*atom;
      p.y[n]=p.y[n-1]+dy;
      p.z[n]=backbone.parts[i].interpolate(j, backbone.parts[i].datum);
      count++;
      }
    G=k*l+G;
    }
  
  for(int step=0;step<2;step++) {
    for(int i=0;i<backbone.parts.size()-1;i++) {
      for(int k=0;k<5;k++) {
        int n;
        n=backbone.parts[i].n+k;
        p.x[n]=0.5*(p.x[n-1]+p.x[n+1]);
        p.y[n]=0.5*(p.y[n-1]+p.y[n+1]);
        p.z[n]=0.5*(p.z[n-1]+p.z[n+1]);
        n=backbone.parts[i].n-k;
        p.x[n]=0.5*(p.x[n-1]+p.x[n+1]);
        p.y[n]=0.5*(p.y[n-1]+p.y[n+1]);
        p.z[n]=0.5*(p.z[n-1]+p.z[n+1]);
        }
      }     
   }     

  status=plg_spherical(projection, &p, 1);
  
  status=plg_save("backbone.plg", PLG_FORMAT_SCAN, &p, 1);
     
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int process_LinearWidth(river_t & backbone, plg_t & p, double *s)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  process backbone width : linearly varying and symetrical width
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int status=0;
  size_t size,count;
  double x0,y0;
  
  backbone.parse();
  
  double lamda,l=0;
  
  count=0;
  for(int i=0;i<backbone.parts.size();i++) {
    double atom=backbone.parts[i].resolution;
    for(int j=0; j<backbone.parts[i].n;j++) {
      double dx,dy,y;
      int n=count;
      l=j*atom;
/*------------------------------------------------------------------------------
      geoid elevation above ellipsoïd  */
      p.z[n]=backbone.parts[i].interpolate(j, backbone.parts[i].datum);
/*------------------------------------------------------------------------------
      river width  */
      s[n]=backbone.parts[i].interpolate(j, backbone.parts[i].width);
      count++;
      }
    }
         
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int process_CosineWidth(river_t & backbone, plg_t & p, double *s)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  process backbone width : linearly varying and symetrical width
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int status=0;
  size_t size,count;
  double x0,y0;
  
  backbone.parse();
  
  double lamda,l=0;
  
  count=0;
  for(int i=0;i<backbone.parts.size();i++) {
    double atom=backbone.parts[i].resolution;
    for(int j=0; j<backbone.parts[i].n;j++) {
      double dx,dy,y;
      int n=count;
      l=j*atom;
/*------------------------------------------------------------------------------
      geoid elevation above ellipsoïd  */
      p.z[n]=backbone.parts[i].interpolate(j, backbone.parts[i].datum);
/*------------------------------------------------------------------------------
      river width  */
      s[n]=backbone.parts[i].interpolate_cosine(j, backbone.parts[i].width);
      count++;
      }
    }
         
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int process_depths(river_t & backbone, plg_t & p, double* & x, double* & y, double* & z, vector<plg_t> & limits, size_t & count)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  process backbone depths
  
  memory allocation to be checked
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int i, status=0;
  size_t size, n;
  double x0,y0;
  plg_t extremities;
  bool equalize=true,debug=false;
  
  backbone.parse();
 
  double lamda,l=0;
  
  for(int k=0;k<4;k++) limits[k].z=new double[limits[k].npt];

/*------------------------------------------------------------------------------
  first count the number of interior nodes */
  count=0;
  
  n=0;
  for(int i=0;i<backbone.parts.size();i++) {
    int start=0,finish=backbone.parts[i].n;
    if(i==0) start++;
    if(i==backbone.parts.size()-1) finish--;
    double atom=backbone.parts[i].resolution;
    for(int j=start; j<finish;j++) {
      extremities=plg_t(limits[0].x[n],limits[0].y[n],limits[1].x[n],limits[1].y[n]);
      plg_t q=plg_resample(extremities,atom,equalize,debug);
      extremities.destroy();
      count+=q.npt-2;
      q.destroy();
      n++;
      }
    }

/*------------------------------------------------------------------------------
  allocate arrays */
  x=new double[count];
  y=new double[count];
  z=new double[count];
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  process nodes depths 
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  count=0;
  n=1;
  for(i=0;i<backbone.parts.size();i++) {
    int start=0,finish=backbone.parts[i].n;
    if(i==0) start++;
    if(i==backbone.parts.size()-1) finish--;
    double atom=backbone.parts[i].resolution;
    double L,lamda=0.0;
    for(int j=start; j<finish;j++) {
/*------------------------------------------------------------------------------
      create interior nodes */
      extremities=plg_t(limits[0].x[n],limits[0].y[n],limits[1].x[n],limits[1].y[n]);
      plg_t q=plg_resample(extremities,atom,equalize,debug);
      extremities.destroy();
      L=q.cartesian_size();
//       backbone.parts[i].shape.parse(L);
/*------------------------------------------------------------------------------
      interpolate first boundary node depth */
      lamda=0.0;
      status=backbone.parts[i].compute_shape(j, L, lamda*L, limits[0].z[n]);
      limits[0].z[n]+=p.z[n];
      lamda=q.cartesian_size(0)/L;
      for(int k=1;k<q.npt-1;k++) {
        x[count]=q.x[k];
        y[count]=q.y[k];
/*------------------------------------------------------------------------------
        interpolate interior node depth */
        status=backbone.parts[i].compute_shape(j, L, lamda*L, z[count]);
        if(status!=0) {
          printf("%s : interpolation failed\n", __func__);
          }
/*------------------------------------------------------------------------------
        add geoid elevation above ellipsoïd, fragile... */
        z[count]+=p.z[n];
        lamda+=q.cartesian_size(k)/L;
        count++;
        }
 /*------------------------------------------------------------------------------
      interpolate second boundary node depth */
      lamda=1.0;
      status=backbone.parts[i].compute_shape(j, L, lamda*L, limits[1].z[n]);
      if(status!=0) {
        printf("%s : interpolation failed\n", __func__);
        status=backbone.parts[i].compute_shape(j, L, lamda*L, limits[1].z[n]);
        }
      limits[1].z[n]+=p.z[n];
      q.destroy();
      n++;
      }
    }
  
  {
  i=0;
  double atom=backbone.parts[i].resolution;
  plg_t & q=limits[2];
  lamda=0.0;
  n=0;
  int j=0;
  double L=q.cartesian_size();
  for(int k=0; k<q.npt; k++) {
/*------------------------------------------------------------------------------
    interpolate interior node depth */
    status=backbone.parts[i].compute_shape(j, L, lamda*L, q.z[k]);
    if(status!=0) {
      printf("%s : interpolation failed\n", __func__);
      }
/*------------------------------------------------------------------------------
    add geoid elevation above ellipsoïd, fragile... */
    q.z[k]+=p.z[n];
    lamda+=q.cartesian_size(k)/L;
    } 
  }
  
  {
  i=backbone.parts.size()-1;
  double atom=backbone.parts[i].resolution;
  plg_t & q=limits[3];
  lamda=0.0;
  n=p.npt-1;
  int j=backbone.parts[i].n-1;
  double L=q.cartesian_size();
  for(int k=0;k<q.npt;k++) {
/*------------------------------------------------------------------------------
    interpolate interior node depth */
    status=backbone.parts[i].compute_shape(j, L, lamda*L, q.z[k]);
    if(status!=0) {
      printf("%s : interpolation failed\n", __func__);
      }
/*------------------------------------------------------------------------------
    add geoid elevation above ellipsoïd, fragile... */
    q.z[k]+=p.z[n];
    lamda+=q.cartesian_size(k)/L;
    } 
  }
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int process_estuary_shorelines(vector<plg_t> & p, vector<plg_t> & q, vector<plg_t> & shorelines, projPJ projection)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  get overall external polygon
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int status=0;
  vector<plg_t> islands, external, targeted, selection;
  bool debug=false;
  string rootname="debug", filename;
  
  for(int k=0;k<p.size();k++){
    vector<plg_t> pp;
    plg_t tmpp;
    tmpp.duplicate(p[k]);
    pp.push_back(tmpp);
    if(debug) status=plg_save("debug-selection-interim-0.plg", PLG_FORMAT_SCAN, pp);
    status=plg_CutAtIntersections(pp, q, PLG_SPHERICAL);
/*------------------------------------------------------------------------------
    if nothing was cut, polygon is inside frame, nothing else to do*/
    if(pp.size()==1 && plg_isclosed(pp[0])==1) {
      islands.push_back(pp[0]);
      continue;
      }
    if(debug) status=plg_save("debug-selection-interim-1.plg", PLG_FORMAT_SCAN, pp);
    
    for(int kk=0;kk<pp.size();kk++) {
      if(pp[kk].npt>1)
        external.push_back(pp[kk]);
      else
        pp[kk].destroy();
      }
    targeted.push_back(p[k]);
    }
    
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  assembly external shorelines and frame in closed polygons (elementary cycles)
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  status=plg_merge(external, q);
  
  status=plg_cartesian(projection, external);
  
  status=plg_save("combo.plg", PLG_FORMAT_SCAN, external);
  
  if(debug) {
    filename=rootname+"-external.plg";
    status=plg_save(filename.c_str(), PLG_FORMAT_SCAN, external);
    }
  
  vector<plg_t> r;
  if(external.size()!=1) {
    
    mesh_t mesh;
    bool strict=false, do_edges=true, debug=false;
    status=plg_export2nei(external, mesh, strict, do_edges, debug);

//     int RecycleCodes=0, SetLimits=1, verbose=0;
//     status=fe_codetable(mesh, RecycleCodes, SetLimits, verbose);

//     int StopOnPinched=0;
//     status=fe_codetable_obsolete(&mesh, RecycleCodes, StopOnPinched, verbose);

    bool check_islands=false;
    int verbose=1;
    r=Nei2Polygons(mesh, check_islands, verbose, debug);
    
    if(debug) {
      filename=rootname+"-elementary-cycles.plg";
      status=plg_save(filename.c_str(), PLG_FORMAT_SCAN, r);
      }
    mesh.destroy();
    }
  else {
    r=plg_duplicate(q);
    }
  
  for(int s=0;s<r.size();s++){
    double area=r[s].area();
    printf("polygon %d: area=%lf\n",s, area);
    if(area<0.0) {
      plg_t *tmp=new plg_t;
      tmp->duplicate(r[s]);
      shorelines.push_back(*tmp);
      }
    }
    
  if(debug) status=plg_save("cycle.plg", PLG_FORMAT_SCAN, r);
  
  status=plg_save("shorelines.plg", PLG_FORMAT_SCAN, shorelines);
         
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int process_mesh(const string & rootname, mesh_t & mesh, river_t & backbone, plg_t & p, vector<plg_t> & polygons, projPJ projection, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  
  process mesh
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int status;
  double *x, *y, *z;
  size_t count;
  mesh_t work;
  string filename;
  
  status=process_depths(backbone, p, x, y, z, polygons, count);

  bool repair=true;
  int verbose=0;
  
  status=plg_CheckClosure(polygons, repair, verbose, debug);
  
  status=fe_createnodes(polygons, x, y, z, count, work);
  delete[] x;
  delete[] y;
  delete[] z;
  
  status=fe_triangulate(work, projection, mesh, 0, debug);
  
  status=fe_spherical(work, projection);
  
  filename=rootname+"-depth.nod";
  status=fe_savenodes(filename.c_str(), NODE_FILE_FORMAT_TRIGRID, work);
  
  work.destroy();
  
  status=fe_spherical(mesh, projection);
  
  filename=rootname+"-depth.nei";
  status=fe_savemesh(filename.c_str(),MESH_FILE_FORMAT_TRIGRID, mesh);
         
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int process_MouthBathymetry(river_t & backbone, plg_t & p, projPJ projection, mesh_t & mesh, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  river mouth bathymetry
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int i, status=0;
  double c, atom=20.0;
  double L0,H0;
  double *s;
  atom_river_t atom_river;
  plg_t limits[4];
  vector<plg_t> polygons;
  
  backbone.destroy();
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  estuary mouth horizontal shape settings
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
   
/*------------------------------------------------------------------------------
  river mouth backbone                                                        */
  double L_mouth=50000.0;
  c=0.0;
  atom_river.set(0, L_mouth, atom, c, M_PI/4.0);
  atom_river.set3D(1.5e+04, 5.e+02, 0.0, 0.0);
  backbone.parts.push_back(atom_river);
  
/*------------------------------------------------------------------------------
  process information to create backbone polygon p                            */
  status=process_01(backbone, p, L0, H0, projection);
  
  limits[0].init(p.npt-1, PLG_INIT_SEPARATE);
  limits[1].init(p.npt-1, PLG_INIT_SEPARATE);
  
  s=new double[p.npt];
  
  status=process_CosineWidth(backbone, p, s);
  
  for(int n=0; n<p.npt-1;n++) {
    vector2D_t u(p.x[n],p.y[n],p.x[n+1],p.y[n+1]);
    vector2D_t v=u.normal();
    double ss=s[n];
    limits[0].x[n]=0.5*(p.x[n]+p.x[n+1])+v.x*ss;
    limits[0].y[n]=0.5*(p.y[n]+p.y[n+1])+v.y*ss;
    ss=500.;
    limits[1].x[n]=0.5*(p.x[n]+p.x[n+1])-v.x*ss;
    limits[1].y[n]=0.5*(p.y[n]+p.y[n+1])-v.y*ss;
    }
    
/*------------------------------------------------------------------------------
  create lateral segments to close polygon, order matters                     */
  plg_t extremities;
  int m;
 
  m=0;
  extremities=plg_t(limits[0].x[m],limits[0].y[m],limits[1].x[m],limits[1].y[m]);
  limits[2]=plg_resample(extremities,atom,true,debug);
  extremities.destroy();
  
  m=p.npt-2;
  extremities=plg_t(limits[0].x[m],limits[0].y[m],limits[1].x[m],limits[1].y[m]);
  limits[3]=plg_resample(extremities,atom,true,debug);
  extremities.destroy();

//   polygons.push_back(p);
  polygons.push_back(limits[0]);
  polygons.push_back(limits[1]);
  polygons.push_back(limits[2]);
  polygons.push_back(limits[3]);
  
  point2D_t center(0.0,0.0);
  double angle=90.0;
  vector2D_t u(+5000.0, -25000.0);
  
  status=plg_rotate_cartesian(polygons, center, angle);

  status=plg_translate_cartesian(polygons, u);

  status=plg_spherical(projection, polygons);
  
  status=plg_save("topo-mouth.plg", PLG_FORMAT_SCAN, polygons);

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  estuary mouth topography settings
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
   
  shape_t      mouth_shape;
  atom_shape_t section;
  
/*------------------------------------------------------------------------------
                                          length, slope, depth (m) depth (m)  */
  section.set_linear(TOPO_FIXED,          2500.0,   NAN,    -50.0,  -15.0);
  mouth_shape.transects.push_back(section);
  section.set_linear(TOPO_RESIZABLE,      1000.0,   NAN,    -15.0,  -15.0);
  mouth_shape.transects.push_back(section);
  
  backbone.parts[0].shapes[0]=mouth_shape;
  backbone.parts[0].shapes[1]=mouth_shape;
  
  status=process_mesh((string) "mouth", mesh, backbone, p, polygons, projection, debug);
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_create_Estuary_Geometry(vector<plg_t> & shorelines, mesh_t & estuary_mesh, mesh_t & coastal_mesh, mesh_t & mouth_mesh, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  
  create typical estuarine shape meshes
    
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/  
{
  int status=0;
  vector<plg_t> polygons;
  plg_t p, limits[4];
  double ref_lat=45.,ref_lon=-1.0;
  double L0=1.e+05, H0=1.e+04;
  double atom=10.0;
  char *pj_parameters=new char[512];
  projPJ projection;
  river_t backbone;
  atom_river_t atom_river;
  
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  estuarine configuration parameters
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  double L_coast= 1.e+04;
  double L_lower= 5.e+04;
  double L_upper= 5.e+04;
  double L_river= 2.e+04;
  
  atom=50.0;
  
  double L_estuary=L_lower+L_upper+L_river;

  size_t size=L_estuary/atom;

  int n1=L_lower/atom;
  int n2=(L_lower+L_upper)/atom;
  int n3=(L_lower+L_upper+L_river)/atom;

  projection=assign_StereoOblique(ref_lat, ref_lon, pj_parameters);
  
  double c;
  double lamda;
  
  status=process_MouthBathymetry(backbone, p, projection, mouth_mesh, debug);
  backbone.destroy();
  
//   c=0;
//   atom_river.set(0, L_coast, atom, c, 0.0);
//   atom_river.set3D(2.e+04, 2.e+04, 0.0, 0.0);
//   backbone.parts.push_back(atom_river);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  estuary/river shoreline backbone line
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
/*------------------------------------------------------------------------------
  lower estuary, weak curvature                                               */
  c=2.*M_PI*2.0/L0/2.0;
  atom_river.set(0, L_lower, atom, c, M_PI/4.0);
  atom_river.set3D(1.e+04, 5.e+02, 0.0, 0.0);
  backbone.parts.push_back(atom_river);
  
/*------------------------------------------------------------------------------
  upper estuary, larger curvature (x2)                                        */
  c=2.*M_PI*2.0/L0;
  atom_river.set(0, L_upper, atom, c, NAN);
  atom_river.set3D(5.e+02, 1.5e+02, 0.0, 1.0);
  backbone.parts.push_back(atom_river);
  
/*------------------------------------------------------------------------------
  river, larger curvature (x2)                                                */
  c=2.*M_PI*2.0/L0;
  atom_river.set(0, L_river, atom, c, NAN);
  atom_river.set3D(1.5e+02, 1.e+02, 1.0, 4.0);
  backbone.parts.push_back(atom_river);
  
/*------------------------------------------------------------------------------
  process information to create backbone polygon p                            */
  status=process_01(backbone, p, L0, H0, projection);
  
  polygons.push_back(p);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  shoreline side line
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  limits[0].init(p.npt-1, PLG_INIT_SEPARATE);
  limits[1].init(p.npt-1, PLG_INIT_SEPARATE);
  
  double *s=new double[p.npt];
  
/*------------------------------------------------------------------------------
  process information to create river width array                             */
  status=process_LinearWidth(backbone, p, s);
  
  for(int n=0; n<p.npt-1;n++) {
    vector2D_t u(p.x[n],p.y[n],p.x[n+1],p.y[n+1]);
    vector2D_t v=u.normal();
    double ss=s[n];
    limits[0].x[n]=0.5*(p.x[n]+p.x[n+1])+v.x*ss;
    limits[0].y[n]=0.5*(p.y[n]+p.y[n+1])+v.y*ss;
    limits[1].x[n]=0.5*(p.x[n]+p.x[n+1])-v.x*ss;
    limits[1].y[n]=0.5*(p.y[n]+p.y[n+1])-v.y*ss;
    }
    
  polygons.push_back(limits[0]);
  polygons.push_back(limits[1]);
  
/*------------------------------------------------------------------------------
  create lateral segments to close polygon, order matters                     */
  plg_t extremities;
  int m;
 
  m=0;
  extremities=plg_t(limits[0].x[m],limits[0].y[m],limits[1].x[m],limits[1].y[m]);
  limits[2]=plg_resample(extremities,atom,true,debug);
  extremities.destroy();
  
  m=p.npt-2;
  extremities=plg_t(limits[0].x[m],limits[0].y[m],limits[1].x[m],limits[1].y[m]);
  limits[3]=plg_resample(extremities,atom,true,debug);
  extremities.destroy();

  polygons.push_back(limits[2]);
  polygons.push_back(limits[3]);
  
  status=plg_spherical(projection, polygons);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  open coast
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/  

  plg_t coast[5];
  
/*------------------------------------------------------------------------------
  should be set through a frame prescription                                  */
  double xmin= -50000.0, ymin= -30000.0;
  double xmax=  +5000.0, ymax= +30000.0;

  frame_t continental(xmin, xmax, ymin, ymax);

//   extremities=plg_t(limits[0].x[0],limits[0].y[0], 0.0, ymax);
//   coast[0]=plg_resample(extremities,atom,true,debug);
//   extremities.destroy();
//   
//   extremities=plg_t(0.0, ymax, xmin, ymax);
//   coast[1]=plg_resample(extremities,atom,true,debug);
//   extremities.destroy();
//   
//   extremities=plg_t(xmin, ymax, xmin, ymin);
//   coast[2]=plg_resample(extremities,atom,true,debug);
//   extremities.destroy();
//   
//   extremities=plg_t(xmin, ymin, 0.0, ymin);
//   coast[3]=plg_resample(extremities,atom,true,debug);
//   extremities.destroy();
//   
//   extremities=plg_t(0.0, ymin,limits[1].x[0],limits[1].y[0]);
//   coast[4]=plg_resample(extremities,atom,true,debug);
//   extremities.destroy();
//   
//   status=plg_concat(coast[0], coast[1], CARTESIAN);
//   status=plg_concat(coast[0], coast[2], CARTESIAN);
//   status=plg_concat(coast[0], coast[3], CARTESIAN);
//   status=plg_concat(coast[0], coast[4], CARTESIAN);

  coast[0]=plg_t(continental, PLG_INIT_SEPARATE);
  status=plg_spherical(projection, &coast[0], 1);
  
  int equalize=1;
  coast[1]=plg_subdivide(coast[0], atom, equalize, PLG_CARTESIAN, debug);
  
  status=plg_spherical(projection, &coast[1], 1);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  finalize shorelines
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
 
  status=plg_save("material.plg", PLG_FORMAT_SCAN, polygons);
  
  polygons.clear();

  plg_t estuary_shorelines, tmp;
  vector<plg_t> coastal, final;
  
  estuary_shorelines.duplicate(limits[0]);
  
  tmp.duplicate(limits[2]);
  status=plg_concat(estuary_shorelines, tmp, CARTESIAN);
  
  tmp.duplicate(limits[1]);
  status=plg_concat(estuary_shorelines, tmp, CARTESIAN);
  
  tmp.duplicate(limits[3]);
  status=plg_concat(estuary_shorelines, tmp, CARTESIAN);
  
  polygons.push_back(estuary_shorelines);
  status=plg_spherical(projection, polygons);
    
  tmp.destroy();
  
  tmp.duplicate(coast[1]);
  coastal.push_back(tmp);
  
  status=process_estuary_shorelines(polygons, coastal, shorelines, projection);

  status=plg_spherical(projection, polygons);
  
  status=plg_save("estuary.plg", PLG_FORMAT_SCAN, polygons);

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  estuarine topography settings
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
 
  shape_t      estuary_shape, river_shape, coastal_shape, lower_shape;
  atom_shape_t section;
  
/*------------------------------------------------------------------------------
                                          length, slope, depth (m) depth (m)  */
  section.set_linear(TOPO_FIXED,           500.0,   NAN,     10.0,  -10.0);
  lower_shape.transects.push_back(section);
  section.set_parabolic(TOPO_RESIZABLE,   1000.0, -10.0,    -15.0,  -10.0);
  lower_shape.transects.push_back(section);
  section.set_linear(TOPO_FIXED,           500.0,   NAN,    -10.0,   10.0);
  lower_shape.transects.push_back(section);

  section.set_linear(TOPO_FIXED,           100.0,   NAN,     10.0,  -10.0);
  estuary_shape.transects.push_back(section);
  section.set_parabolic(TOPO_RESIZABLE,   1000.0, -10.0,    -15.0,  -10.0);
  estuary_shape.transects.push_back(section);
  section.set_linear(TOPO_FIXED,           100.0,   NAN,    -10.0,   10.0);
  estuary_shape.transects.push_back(section);

  section.set_linear(TOPO_FIXED,            20.0,   NAN,     10.0,  -10.0);
  river_shape.transects.push_back(section);
  section.set_linear(TOPO_RESIZABLE,      1000.0,   NAN,    -10.0,  -10.0);
  river_shape.transects.push_back(section);
  section.set_linear(TOPO_FIXED,            20.0,   NAN,    -10.0,   10.0);
  river_shape.transects.push_back(section);
 
  section.set_linear(TOPO_FIXED,           500.0,   NAN,     10.0,  -10.0);
  coastal_shape.transects.push_back(section);
  section.set_linear(TOPO_FIXED,          5000.0,   NAN,    -10.0,  -50.0);
  coastal_shape.transects.push_back(section);
  section.set_linear(TOPO_RESIZABLE,     50000.0,   NAN,    -50.0,  -100.0);
  coastal_shape.transects.push_back(section);
 
  backbone.parts[0].shapes[0]=lower_shape;
  backbone.parts[0].shapes[1]=estuary_shape;
  
  backbone.parts[1].shapes[0]=estuary_shape;
  backbone.parts[1].shapes[1]=river_shape;
  
  backbone.parts[2].shapes[0]=river_shape;
  backbone.parts[2].shapes[1]=river_shape;
 
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  estuary depth mesh
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
 
//   double *x, *y, *z;
//   size_t count;
//   mesh_t work;
  
  plg_destroy_entries(polygons);
//   polygons.clear();
  
  for(int k=0;k<4;k++) polygons.push_back(limits[k]);
  status=plg_spherical(projection, polygons);
  
  status=process_mesh((string) "estuary", estuary_mesh, backbone, p, polygons, projection, debug);
//   status=process_depths(backbone, p, x, y, z, polygons, count);
// 
//   bool repair=true;
//   int verbose=0;
//   
//   status=plg_CheckClosure(polygons, repair, verbose, debug);
//   
//   status=fe_createnodes(polygons, x, y, z, count, work);
//   delete[] x;
//   delete[] y;
//   delete[] z;
//   
//   status=fe_triangulate(work, projection, estuary_mesh, 0, debug);
//   
//   status=fe_spherical(work, projection);
//   status=fe_savenodes("estuary-depth.nod", NODE_FILE_FORMAT_TRIGRID, work);
//   
//   work.destroy();
//   
//   status=fe_spherical(estuary_mesh, projection);
//   status=fe_savemesh("estuary-depth.nei",MESH_FILE_FORMAT_TRIGRID, estuary_mesh);
    
  backbone.destroy();
  
  plg_destroy_entries(polygons);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  coastal depth mesh
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  frame_t frame=continental;
   
  atom=100.0;
   
  extremities=plg_t(frame.xmax, frame.ymin, frame.xmax, frame.ymax);
  limits[0]=plg_resample(extremities,atom,true,debug);
  extremities.destroy();
  
  extremities=plg_t(frame.xmin, frame.ymin, frame.xmin, frame.ymax);
  limits[1]=plg_resample(extremities,atom,true,debug);
  extremities.destroy();
    
/*------------------------------------------------------------------------------
  create lateral segments to close polygon, order matters                     */
 
  m=0;
  extremities=plg_t(limits[0].x[m],limits[0].y[m],limits[1].x[m],limits[1].y[m]);
  limits[2]=plg_resample(extremities,atom,true,debug);
  extremities.destroy();
  
  m=limits[0].npt-1;
  extremities=plg_t(limits[0].x[m],limits[0].y[m],limits[1].x[m],limits[1].y[m]);
  limits[3]=plg_resample(extremities,atom,true,debug);
  extremities.destroy();

  status=plg_add_entries(polygons, limits, 4);
  
  status=plg_spherical(projection, polygons);
   
/*------------------------------------------------------------------------------
  "fake" coastal backbone                                                     */
  c=0.0;
  atom_river.set(0, limits[0].cartesian_size(), atom, c, 0.0);
  atom_river.set3D(1.e+04, 5.e+02, 0.0, 0.0);
  backbone.parts.push_back(atom_river);
  
  backbone.parts[0].shapes[0]=coastal_shape;
  backbone.parts[0].shapes[1]=coastal_shape;
  
  status=process_mesh((string) "coastal", coastal_mesh, backbone, p, polygons, projection, debug);
//   status=process_depths(backbone, p, x, y, z, polygons, count);
// 
//   status=plg_CheckClosure(polygons, repair, verbose, debug);
//   
//   status=fe_createnodes(polygons, x, y, z, count, work);
//   delete[] x;
//   delete[] y;
//   delete[] z;
//   
//   status=fe_triangulate(work, projection, coastal_mesh, 0, debug);
//   
//   status=fe_spherical(work, projection);
//   status=fe_savenodes("coastal-depth.nod", NODE_FILE_FORMAT_TRIGRID, work);
//   
//   work.destroy();
//   
//   status=fe_spherical(coastal_mesh, projection);
//   status=fe_savemesh("coastal-depth.nei",MESH_FILE_FORMAT_TRIGRID, coastal_mesh);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  create unstructured mesh "esthetic" open limit
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  plg_t opened, fermeture;
  
  double x0, y0;
  double radius=30.e+03;
  range_t<double> angles=range_t<double>(M_PI/2.,3*M_PI/2);
  int npt=100;
  
  x0=0.5*(limits[0].x[0]+limits[0].x[limits[0].npt-1]);
  y0=0.5*(limits[0].y[0]+limits[0].y[limits[0].npt-1]);
  
  double a=40.e+03,b=25.e+03;
  status=plg_cartesian_ellipse(x0, y0, a, b, angles, npt, opened);
  
  fermeture.init(4,PLG_INIT_SEPARATE);
  
  fermeture.x[0]=opened.x[0];
  fermeture.y[0]=opened.y[0];
  
  fermeture.x[1]=opened.x[0]+2500.;
  fermeture.y[1]=opened.y[0];
  
  fermeture.x[2]=opened.x[opened.npt-1]+2500.;
  fermeture.y[2]=opened.y[opened.npt-1];
  
  fermeture.x[3]=opened.x[opened.npt-1];
  fermeture.y[3]=opened.y[opened.npt-1];
  
  status=plg_spherical(projection, &opened, 1);
  status=plg_spherical(projection, &fermeture, 1);
  
  status=plg_save("opened.plg", PLG_FORMAT_SCAN, &opened, 1);
  
  status=plg_concat(opened, fermeture);
  
  status=plg_save("coastal-mesh.plg", PLG_FORMAT_SCAN, &opened, 1);
   
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_create_Estuary(bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  
  create typical estuarine shape meshes
    
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/  
{
  int status=0;
  grid_t TOPOgrid;
  float *TOPOdata[4], TOPOmask=1.e+10;
  mesh_t coastal_mesh, estuary_mesh, mouth_mesh;
  vector<plg_t> shorelines;
  frame_t frame;
  double grid_resolution=0.0001;
  string rootname="academic", output, filename;
  
//   bool debug=false;

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  create estuary bathymetry
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  status=fe_create_Estuary_Geometry(shorelines, estuary_mesh, coastal_mesh, mouth_mesh, debug);

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  create bathymetry structured grid 
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  printf("#################################################################\n");
  printf("create bathymetry on structured grid: dx=%lf (%lf m) arc-sec dy=%lf (%lf m) arc-sec\n",
         grid_resolution*3600., grid_resolution*110000., grid_resolution*3600., grid_resolution*110000.);
    
  double dx=grid_resolution;
  double dy=grid_resolution;
  
  frame=plg_spherical_minmax(shorelines);
  
  TOPOgrid=map_Rgrid(frame, dx, dy, 0);
  map_printgrid(TOPOgrid);
  
  for(int k=0; k<4; k++) TOPOdata[k]=new float[TOPOgrid.Hsize()];
  
  int algorithm=0, verbose=0;  
 
/*------------------------------------------------------------------------------
  estuary depths structured file */
  status=fe_initaffine(&estuary_mesh);
  status=fe_mapVdepth(estuary_mesh, TOPOgrid, TOPOdata[0], TOPOmask, algorithm, verbose);
  
  output=rootname+".estuary.grd";
  status=topo_save(output.c_str(), "grd-float", TOPOgrid, TOPOdata[0], TOPOmask, debug);
  
/*------------------------------------------------------------------------------
  coastal depths structured file */
  status=fe_initaffine(&coastal_mesh);
  status=fe_mapVdepth(coastal_mesh, TOPOgrid, TOPOdata[1], TOPOmask, algorithm, verbose);
  
//   status=checks(rootname,TOPOgrid,TOPOdata,TOPOmask,false);
  
  output=rootname+".coastal.grd";
  status=topo_save(output.c_str(), "grd-float", TOPOgrid, TOPOdata[1], TOPOmask, debug);
   
/*------------------------------------------------------------------------------
  estuary mouth depths structured file */
  status=fe_initaffine(&mouth_mesh);
  status=fe_mapVdepth(mouth_mesh, TOPOgrid, TOPOdata[2], TOPOmask, algorithm, verbose);
  
  output=rootname+".mouth.grd";
  status=topo_save(output.c_str(), "grd-float", TOPOgrid, TOPOdata[2], TOPOmask, debug);
  
/*------------------------------------------------------------------------------
  topography merging */
  for(size_t m=0;m<TOPOgrid.Hsize(); m++) {
    TOPOdata[3][m]=TOPOmask;
/*------------------------------------------------------------------------------
    keep shallowest of mouth and coastal in modified coastal */
    if(TOPOdata[2][m]!=TOPOmask and TOPOdata[1][m]!=TOPOmask) TOPOdata[1][m]=min(TOPOdata[1][m], TOPOdata[2][m]);
/*------------------------------------------------------------------------------
    keep deepest of estuary and modified coastal */
    if(TOPOdata[0][m]!=TOPOmask and TOPOdata[1][m]!=TOPOmask) TOPOdata[3][m]=max(TOPOdata[0][m], TOPOdata[1][m]);
    else if(TOPOdata[0][m]==TOPOmask) TOPOdata[3][m]=TOPOdata[1][m];
    else if(TOPOdata[1][m]==TOPOmask) TOPOdata[3][m]=TOPOdata[0][m];
/*------------------------------------------------------------------------------
    change sign */
    if(TOPOdata[3][m]!=TOPOmask)  TOPOdata[3][m]=-TOPOdata[3][m]; 
/*------------------------------------------------------------------------------
    copy and set dry altitude (20 m)*/
//     if(TOPOdata[3][m]==TOPOmask) TOPOdata[3][m]=20.0;
    }
  
  rootname="academic-merged";
  
//   status=grd_mirror_r( TOPOgrid, TOPOgrid.nx, TOPOdata[3], TOPOmask);
  status=topo_save(rootname, output, "netcdf", TOPOgrid, TOPOdata[3], TOPOmask, debug);
    
/*------------------------------------------------------------------------------
  smooth merging */
#if 0
  frame_t smoothing_frame;
  status=fe_minmax(mouth_mesh, smoothing_frame);
  
  vector<plg_t> polygons;
//   plg_t p(smoothing_frame);
//   polygons.push_back(p);

  float radius, zmin, zmax;
  int filter=LOESS;
  char *trusted=0;
  
  zmin=-100.0;
  zmax=   0.0;
  
  range_t<float> range(zmin,zmax);
  
  radius=2000.0/110000.0;
  
  bool keep_masked=true;
  status=plg_load("smoothing.plg", PLG_FORMAT_UNKNOWN, polygons);
  
  status=topo_smooth_cartesian(TOPOgrid, TOPOdata[0], TOPOmask, radius, range, keep_masked, polygons, trusted, filter);

  rootname="academic-smoothed";
  
  for(size_t m=0;m<TOPOgrid.Hsize(); m++) {
    if(TOPOdata[0][m]!=TOPOmask) TOPOdata[3][m]=TOPOdata[0][m];
    else TOPOdata[3][m]=20.0;
    }
#else

  vector<plg_t> p,q, intersection;
  int RecycleCodes=0, SetLimits=1;
  
  verbose=0;
  
  status=fe_edgetable(&estuary_mesh, 0, 0);
  status=fe_codetable(estuary_mesh, RecycleCodes, SetLimits, verbose);
  status=fe_limits2poly(estuary_mesh, p, (char *) 0, true);
  
  status=fe_edgetable(&mouth_mesh, 0, 0);
  status=fe_codetable(mouth_mesh, RecycleCodes, SetLimits, verbose);
  status=fe_limits2poly(mouth_mesh,   q, (char *) 0, true);
  
  int target=0;
  intersection=plg_extract(p, q, rootname, PLG_SPHERICAL, target, debug);
  status=plg_save("intersection.plg", PLG_FORMAT_SCAN, intersection);
  
/*------------------------------------------------------------------------------
  identify open pieces */
  vector<plg_t> splitted, opened;
  
  splitted=plg_split(intersection, debug);
  
  for(int s=0;s<splitted.size();s++) {
    if(splitted[s].flag[0]!='T') opened.push_back(splitted[s]);
    }

  status=plg_save("limit.plg", PLG_FORMAT_SCAN, opened);
  
  frame_t f=plg_spherical_minmax(intersection);
  for(size_t m=0;m<TOPOgrid.Hsize(); m++) {
    double x,y;
//     TOPOdata[3][m]=TOPOmask;
    if(TOPOdata[0][m]==TOPOmask) continue;
    TOPOgrid.xy(m,x,y);
    if(not f.inside(x,y)) continue;
    int inside=0;
    inside=plg_single(intersection[0],x,y,&inside,PLG_SPHERICAL,0);
    if(inside!=PLG_POINT_INTERIOR) continue;
/*------------------------------------------------------------------------------
    copy and set dry altitude (20 m)*/
    int n1,n2;
    vector2D_t u;
    double d, scale=0.01;
    status=plg_NearestSegment(opened[0], x, y, n1, n2, d, u);
    double r=exp(-d/scale)/exp(0.0);
    float depth;
    depth=r*TOPOdata[0][m]+(1-r)*TOPOdata[1][m];
    if(TOPOdata[0][m]!=TOPOmask and TOPOdata[1][m]!=TOPOmask) TOPOdata[3][m]=-max(depth, TOPOdata[1][m]);
    }

  for(size_t m=0;m<TOPOgrid.Hsize(); m++) {
    if(TOPOdata[3][m]==TOPOmask) TOPOdata[3][m]=10.0;
    }
    
#endif

  rootname="academic-smoothed";
  status=topo_save(rootname, output, "netcdf", TOPOgrid, TOPOdata[3], TOPOmask, debug);
    
  STDOUT_BASE_LINE("Finished with %s\n",__func__);
  
  TOPOgrid.free();
  delete[] TOPOdata[0];
  delete[] TOPOdata[1];
  
  return(status);
}
