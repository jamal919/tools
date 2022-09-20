

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

/**-************************************************************************

  2D mesh generator

  F. Lyard, 2009, CNRS/LEGOS, Toulouse, France

  e-mail: florent.lyard@legos.obs-mip.fr

***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string>


#include "tools-structures.h"

#include "constants.h"

#include "functions.h"
#include "geo.h"
#include "fe.h"
#include "archive.h"
#include "map.h"
#include "polygones.h"
#include "grd.h"
#include "list.h"
#include "maths.h"
#include "map.def"
#include "mass.h"

#include "zapper.h"     /*  rutin.h contains common utility routines  */

using namespace std;

// #define MSH_MAXSIZE_OPEN     1
// #define MSH_MAXSIZE_SHELF    2
// #define MSH_TIDAL_WAVELENGTH 3
// #define MSH_TOPO_SLOPE       4
// #define MSH_SURFWAVE_LENGTH  5
// #define MSH_SURFWAVE_DZ      6
// #define MSH_SURFWAVE_CFL     7

extern  void fe_LGP0_to_LGP1(double *in, double *out, mesh_t mesh);
#include "statistic.h"


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  float plg_resolution(plg_t & target, int k, grid_t & grid, float *density, float mask, criteria_t & criteria,  int equalize)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  float *z,size=1.e+10;
  equalize=1;
  double km2m=1000.0;
  
  plg_t p=plg_t(target.x[k],target.y[k],target.x[k+1],target.y[k+1]);
  double radius=criteria.cellsize*km2m;
    
/*------------------------------------------------------------------------------
  sample interval at base grid granularity */
  plg_t q=plg_resample(p, radius, equalize);
  
  z=new float[q.npt];
  for(size_t l=0;l<q.npt;l++) {
    status=map_interpolation(grid, density,mask,q.x[l],q.y[l],&z[l]);
    if(z[l]==mask) {
      double t,p;
      projection_to_geo(grid.projection, &p,&t,q.x[l],q.y[l]);
      printf("%s: masked value at t=%lf p=%lf (x=%lf y=%lf)\n",__func__,t,p,q.x[l],q.y[l]);
      continue;
      }
    }
  for(size_t l=0;l<q.npt;l++) {
    if(z[l]==mask) {
//       printf("%lf %lf\n",q.t[l],q.p[l]);
      continue;
      }
    size=MIN(size,z[l]);
    }

  p.destroy();
  q.destroy();
  
  delete[] z;
  
  return(size);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  plg_t plg_sample(plg_t target, grid_t & grid, float *density, float mask, criteria_t & criteria, int equalize, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  re-sample M flag sections of polygon with tuned sampling radius 

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int status;
  plg_t *tmp,out;
  int concat[2];
  double km2m=1000.0;
  double radius=criteria.cellsize*km2m;
  
  equalize=1;
  
  tmp=new plg_t[target.npt-1];
  
  concat[0]=0;
  
  for(size_t k=0;k<target.npt-1;k++) {
    if(target.flag[k]=='M') {
/*------------------------------------------------------------------------------
      next step necessary to insure new boundary to lon/lat aligned */
      double check=target.cartesian_size(k);
      int n=7*check/radius+2;
      plg_t p=plg_t(target.t[k],target.p[k],target.t[k+1],target.p[k+1],n,PLG_SPHERICAL);
      
      status=plg_cartesian(grid.projection, &p, 1);
      status=plg_spherical(grid.projection, &p, 1);
    
/*------------------------------------------------------------------------------
      sample interval at base grid granularity, working polygon */
      plg_t q=plg_resample(p, radius, equalize);
      
      if(q.npt<2) {
        printf("sampling with radius=%lf m failed: %lf %lf %lf %lf \n",radius,target.t[k],target.p[k],target.t[k+1],target.p[k+1]);
        return(out);
        }
      status=plg_spherical(grid.projection, &q, 1);
    
/*------------------------------------------------------------------------------
      then decimate over-densified points */
      for(size_t l=1;l<q.npt-1;l++) {
        double L1=q.cartesian_size(l-1);
        double L2=q.cartesian_size(l);
        double L=MAX(L1,L2);
        double C1=plg_resolution(q, l-1, grid, density, mask, criteria,  equalize);
        double C2=plg_resolution(q, l,   grid, density, mask, criteria,  equalize);
        double C=MIN(C1,C2);
        if(L<0.95*C) {
          status= plg_delete_point(&q,l);
          l--;
          }
        }
/*------------------------------------------------------------------------------
      equalize sampling */
      for(int pass=0;pass<10;pass++) 
        for(size_t l=1;l<q.npt-1;l++) {
          double L1=q.cartesian_size(l-1);
          double L2=q.cartesian_size(l);
          double L=MAX(L1,L2);
          double C1=plg_resolution(q, l-1, grid, density, mask, criteria,  equalize);
          double C2=plg_resolution(q, l,   grid, density, mask, criteria,  equalize);
          double C=C1+C2, r=1.-C1/C;
          q.x[l]=r*q.x[l-1]+(1.-r)*q.x[l+1];
          q.y[l]=r*q.y[l-1]+(1.-r)*q.y[l+1];
          }
      status=plg_spherical(grid.projection, &q, 1);
      q.SetFlag('M');
      
      tmp[k].duplicate(q);
      if(k!=0) {
        concat[1]=k;
/*------------------------------------------------------------------------------
        aggregate new pieces in first polygon*/
        status=plg_concat(concat, tmp, k+1, PLG_CARTESIAN);
        if(status!=0) {
          printf("concatenation issue : %d\n",k);
          }
        }
      p.destroy();
      q.destroy();
      }
    else if(target.flag[k]=='T') {
      plg_t p=plg_t(target.x[k],target.y[k],target.x[k+1],target.y[k+1]);
      p.SetFlag('T');
      tmp[k].duplicate(p);
      if(k!=0) {
        concat[1]=k;
        status=plg_concat (concat, tmp, k+1, PLG_CARTESIAN);
        if(status!=0) {
          printf("concatenation issue : %d\n",k);
          }
        }
      p.destroy();      
      }
    else if(target.flag[k]=='X') {
      plg_t p=plg_t(target.x[k],target.y[k],target.x[k+1],target.y[k+1]);
      p.SetFlag('M');
      tmp[k].duplicate(p);
      if(k!=0) {
        concat[1]=k;
        status=plg_concat (concat, tmp, k+1, PLG_CARTESIAN);
        if(status!=0) {
          printf("concatenation issue : %d\n",k);
          }
        }
      p.destroy();      
      }
    }
  
  out.duplicate(tmp[0]);
  
  plg_deletep(&tmp,target.npt-1);
  
  return(out);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  plg_t plg_resample_rigid_obsolete(plg_t target, projPJ projection, double radius, int equalize, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  re-sample T flag sections of polygon with fixed sampling radius 

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int status;
  plg_t *tmp,out;
  int concat[2];
  double km2m=1000.0;
  size_t npassed=0, newlines=0, nlines=target.npt-1;
  bool debug=false;
  
  equalize=1;
  
  tmp=new plg_t[nlines];
  
  concat[0]=0;
  
  radius=radius*km2m;

/*-------------------------------------------------------------------------------
  reorganize polygon with possibly multiple codes in multiple polygons of similar code */
  size_t k=0;
  while (npassed!=nlines) {
    if(target.flag[k]=='M') {
      npassed++;
      size_t l=(k+1) % (nlines);
      while (target.flag[l]=='M') {
        l=(l+1) % (nlines);
        npassed++;
        if (npassed==nlines) break;
        }
/*------------------------------------------------------------------------------
      duplicate open sections */
      tmp[newlines].duplicate(target,k,l);
      
      if(newlines!=0) {
        concat[1]=newlines;
        status=plg_concat (concat, tmp, newlines+1, PLG_CARTESIAN);
        if(status!=0) {
          printf("concatenation issue : %d\n",newlines);
          }
        }
      k=l;
      newlines++;
      }
    else {
      npassed++;
      size_t l=(k+1) % (nlines);
      while (target.flag[l]=='T') {
        l=(l+1) % (nlines);
        npassed++;
        if (npassed==nlines) break;
        }

/*------------------------------------------------------------------------------
      extract rigid section*/
      plg_t p;
      p.duplicate(target,k,l);
      status=plg_spherical(projection, &p, 1);
      if(debug) status=plg_save ("mesh-tmp.plg", PLG_FORMAT_SCAN, &p, 1);
    
/*------------------------------------------------------------------------------
      sample interval at base grid granularity, working polygon */
      double L=plg_length(p, PLG_CARTESIAN);
//       if(verbose==1) {
//         printf("re-sampling from vertex %d to vertex %d, \n",k,l);
//         printf("initial length/number of vertices : %lf km/%5d (average resolution=%lf km)\n",p.npt,L/1000.,L/(p.npt-1)/1000.);
//         }
/*------------------------------------------------------------------------------
      re-sample if needed rigid sections */
      plg_t q;
      if(L<radius) {
        equalize=0;
        q.npt=0;
        }
      else if(L<6*radius) {
        equalize=1;
        q=plg_resample(p, L/10.0, equalize);
        }
      else {
        equalize=1;
        q=plg_resample(p, radius, equalize);
        }
      status=plg_spherical(projection, &q, 1);
      double LL=plg_length(q, PLG_CARTESIAN);
// /*------------------------------------------------------------------------------
//       fractal nature of shorelines, check final size as it can be much smaller than original one */
//       if(LL<2*radius) {
//         q.destroy();
//         }
      if(q.npt>3) {
        q.SetFlag('T');
//         double LL=plg_length(q, PLG_CARTESIAN);
// //        L=plg_length(q, PLG_SPHERICAL);
        if(verbose==1) {
          printf("re-sampling from vertex %d to vertex %d, \n",k,l);
          printf("initial length/number of vertices : %lf km/%5d (average resolution=%lf km)\n",p.npt,L/1000.,L/(p.npt-1)/1000.);
          printf("final   length/number of vertices : %lf km/%5d (average resolution=%lf km)\n",q.npt,LL/1000.,LL/(q.npt-1)/1000.);
          }
        tmp[newlines].duplicate(q);
        if(newlines!=0) {
          concat[1]=newlines;
          status=plg_concat(concat, tmp, newlines+1, PLG_CARTESIAN);
          if(status!=0) {
            printf("concatenation issue : %d\n",k);
            }
          }
        }
      p.destroy();
      q.destroy();
      k=l;
      newlines++;
      }
    }

  out.duplicate(tmp[0]);
  
  plg_deletep(&tmp,target.npt-1);
  
  return(out);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  plg_t plg_resample_rigid(plg_t target, projPJ projection, vector<double> radius, vector<double> length, char mask, int equalize, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 
  re-sample non mask flag sections of polygon with uniform sampling radius
  
  now expect distance in meters
  
  in earlier version (<2017), mask was implicitely set as 'M' 

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
{
  int status;
  plg_t *tmp,out;
  int concat[2];
  size_t npassed=0, nlines=target.npt-1;
  size_t additional_plg=0;
  bool debug=false;
  char flag;
  
  if(nlines==1) {
    if(verbose==1) printf("sampling limited to one segment...\n");
    }
  
  equalize=1;
  
  tmp=new plg_t[nlines];
  
  concat[0]=0;
  
/*-------------------------------------------------------------------------------
  reorganize polygon with possibly multiple flags in multiple polygons of similar flag */
  size_t k=0;
  while (npassed!=nlines) {
    if(target.flag[k]==mask) {
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      no sampling
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
      npassed++;
      size_t l=(k+1) % (nlines);
      while (target.flag[l]==mask) {
        l=(l+1) % (nlines);
        npassed++;
        if (npassed==nlines) break;
        }
/*------------------------------------------------------------------------------
      duplicate open sections */
      tmp[additional_plg].duplicate(target,k,l);
      
      if(additional_plg!=0) {
        concat[1]=additional_plg;
        status=plg_concat (concat, tmp, additional_plg+1, PLG_CARTESIAN);
        if(status!=0) {
          printf("concatenation issue : %d\n",additional_plg);
          }
        }
      k=l;
      additional_plg++;
      }
    else if(target.flag[k]=='M') {
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      "open limits" sampling : straight line sub-division at radius resolution
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
      npassed++;
      flag=target.flag[k];
      size_t l=(k+1) % (nlines);
      while ((target.flag[l]==flag) && (npassed<nlines)){
        l=(l+1) % (nlines);
        npassed++;
        if (npassed==nlines) break;
        }
/*------------------------------------------------------------------------------
      extract open section*/
      plg_t p;
      p.duplicate(target,k,l);
      status=plg_spherical(projection, &p, 1);
      plg_t q;
      q=plg_subdivide(p, radius[0], equalize, PLG_CARTESIAN, debug);
      status=plg_spherical(projection, &q, 1);
      q.SetFlag(flag);
      tmp[additional_plg].duplicate(q);
      if(additional_plg!=0) {
        concat[1]=additional_plg;
        status=plg_concat(concat, tmp, additional_plg+1, PLG_CARTESIAN);
        if(status!=0) {
          printf("concatenation issue : %d\n",k);
          }
        }
      p.destroy();
      q.destroy();
      additional_plg++;
      k=l;
      }
    else if(target.flag[k]=='T') {
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      "shorelines" sampling : resampling following existing segments
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
      npassed++;
      flag=target.flag[k];
      size_t l=(k+1) % (nlines);
      while ((target.flag[l]==flag) && (npassed<nlines)){
        l=(l+1) % (nlines);
        npassed++;
        if (npassed==nlines) break;
        }

/*------------------------------------------------------------------------------
      extract rigid section*/
      plg_t p;
      p.duplicate(target,k,l);
      status=plg_spherical(projection, &p, 1);
      if(debug) status=plg_save ("mesh-tmp.plg", PLG_FORMAT_SCAN, &p, 1);
    
/*------------------------------------------------------------------------------
      sample interval at base grid granularity, working polygon */
      double L=plg_length(p, PLG_CARTESIAN);
/*------------------------------------------------------------------------------
      re-sample if needed rigid sections */
      plg_t q;
      int r=-1;
      for(int s=0;s<length.size();s++) {
        if(L>length[s]) {
          r=s;
          break;
          }
        }

      if(r==-1) {
        equalize=0;
        q.npt=0;
        }
      else {
        equalize=1;
        q=plg_resample(p, radius[r], equalize);
        }
      status=plg_spherical(projection, &q, 1);
      double LL=plg_length(q, PLG_CARTESIAN);
/*------------------------------------------------------------------------------
      fractal nature of shorelines, check final size as it can be much smaller than original one */
      if(LL<2*radius[0]) {
        q.destroy();
        }
      if(q.npt>3) {
        if(verbose==1) {
          printf("re-sampling from vertex %d to vertex %d, \n",k,l);
          printf("initial length/number of vertices : %lf km/%5d (average resolution=%lf km)\n",p.npt,L/1000.,L/(p.npt-1)/1000.);
          printf("final   length/number of vertices : %lf km/%5d (average resolution=%lf km)\n",q.npt,LL/1000.,LL/(q.npt-1)/1000.);
          }
        }
      else {
        q.destroy();
        q.duplicate2(p);
        }
      q.SetFlag(flag);
      tmp[additional_plg].duplicate(q);
      if(additional_plg!=0) {
        concat[1]=additional_plg;
        status=plg_concat(concat, tmp, additional_plg+1, PLG_CARTESIAN);
        if(status!=0) {
          printf("concatenation issue : %d\n",k);
          }
        }
      additional_plg++;
      p.destroy();
      q.destroy();
      k=l;
      }
    }

  out.duplicate(tmp[0]);
  
/*------------------------------------------------------------------------------
  projection got involved, accuracy issue may appear at extremities; fix that */
  plg_copy_point(&out, 0, target, 0);
  plg_copy_point(&out, out.npt-1, target, target.npt-1);
  
  plg_deletep(&tmp,nlines);
  
  return(out);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  plg_t plg_resample_rigid_standard(plg_t target, projPJ projection, double radius, char mask, int equalize, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  plg_t out;
  vector<double> Vradius, Vlength;
  
  Vradius.push_back(radius);
  Vlength.push_back(radius*6.0);
  
  Vradius.push_back(radius/2.0);
  Vlength.push_back(radius*3.0);
  
  Vradius.push_back(radius/4.0);
  Vlength.push_back(radius*1.0);
  
  out=plg_resample_rigid(target, projection, Vradius, Vlength, mask, equalize, verbose);
  
  return(out);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  plg_t plg_resample_rigid_variable(const plg_t & target, projPJ projection,const SGfield_t<float> & resolution, double radius, char mask, int equalize, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------
  
  sample target polygons with variable sampling radius; radius transitions might 
  be improved in further developments
  
  expect distance in meters

-------------------------------------------------------------------------------*/
{
  int status;
  plg_t out;
  vector<double> Vradius, Vlength;
  float r,s;
  int k,l;
  double t,p;
  vector<plg_t>  subtargets;
  vector<float> subradius;
  extern int field_interpolation(const SGfield_t<float> & field, double x, double y, float *z);
  
  k=0;
  t=0.5*(target.t[k]+target.t[k+1]);
  p=0.5*(target.p[k]+target.p[k+1]);
  status=field_interpolation(resolution, t, p, &r);
  if(r==resolution.mask) r=radius;
         
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  split/aggregate sub-segments with similar sampling prescription 
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  for(l=1; l<target.npt-1;l++) {
    t=0.5*(target.t[l]+target.t[l+1]);
    p=0.5*(target.p[l]+target.p[l+1]);
    status=field_interpolation(resolution, t, p, &s);
    if(s==resolution.mask) s=radius;
    double L=plg_length(target, k, l, PLG_CARTESIAN);
    double rmax=max(r,s);
//     if( (fabs(r-s)/r > 0.1) and (L>5*r) ) {
    if( fabs(r-s)/(r+s)/2. > 0.1 and L>5*r ) {
/*------------------------------------------------------------------------------
      significant resolution change */
      plg_t tmp;
      tmp.duplicate(target,k,l);
      subtargets.push_back(tmp);
      subradius.push_back(r);
      k=l;
      r=s;
      }
    }
    
  plg_t tmp;
  tmp.duplicate(target,k,l);
  subtargets.push_back(tmp);
  subradius.push_back(r);
 
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  sample sub-segments 
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  for(k=0;k<subtargets.size();k++){
    plg_t tmp;
    double rr=subradius[k];
    Vradius.push_back(rr);
    Vlength.push_back(rr*6.0);
    
    Vradius.push_back(rr/1.2);
    Vlength.push_back(rr*5.0);
  
    Vradius.push_back(rr/1.5);
    Vlength.push_back(rr*4.0);
  
    Vradius.push_back(rr/2.0);
    Vlength.push_back(rr*3.0);
  
    Vradius.push_back(rr/4.0);
    Vlength.push_back(rr*1.0);
  
    tmp=plg_resample_rigid(subtargets[k], projection, Vradius, Vlength, mask, equalize, verbose);
    status=plg_concat(out,tmp, PLG_CARTESIAN);
    Vradius.clear();
    Vlength.clear();
    }

//   for(k=0;k<subtargets.size();k++){
//     plg_t tmp;
//     double rr=subradius[k];
//     Vradius.push_back(rr);
//     Vlength.push_back(rr*2.0);
//     
//     Vradius.push_back(rr/1.2);
//     Vlength.push_back(rr*2.0);
//   
//     Vradius.push_back(rr/1.5);
//     Vlength.push_back(rr*2.0);
//   
//     Vradius.push_back(rr/2.0);
//     Vlength.push_back(rr*2.0);
//   
//     Vradius.push_back(rr/4.0);
//     Vlength.push_back(rr*2.0);
//   
//     tmp=plg_resample_rigid(subtargets[k], projection, Vradius, Vlength, mask, equalize, verbose);
//     status=plg_concat(out,tmp, PLG_CARTESIAN);
//     Vradius.clear();
//     Vlength.clear();
//     }
    
  return(out);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  plg_t plg_check(plg_t target, grid_t & grid, float *density, float mask, criteria_t & criteria,  int equalize)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  plg_t *tmp,out;
//  point2D_t point;
  plg_point_t point;
  int concat[2];
  double km2m=1000.0;
  
  equalize=1;
  
  tmp=new plg_t[target.npt-1];
  
  concat[0]=0;
  
  for(size_t k=0;k<target.npt-1;k++) {
    
    plg_t p=plg_t(target.x[k],target.y[k],target.x[k+1],target.y[k+1]);
    
    if(target.flag[k]=='M') {
    
      double check=p.cartesian_size(k);
//       double L1=q.cartesian_size(l-1);
//       double L2=q.cartesian_size(l);
//       double L=MAX(L1,L2);
    
      double radius=criteria.cellsize*km2m;
      double C;
      
/*------------------------------------------------------------------------------
      sample interval at base grid garnularity */
      plg_t q=plg_resample(p, radius, equalize);
    
      for(size_t l=0;l<q.npt-1;l++) {
        double C1=plg_resolution(q, l, grid, density, mask, criteria,  equalize);
        double C=MIN(C1,C);
        }
      p.destroy();
      q.destroy();
      }
    else {
      p.destroy();      
      }
    }
  
  plg_deletep(&tmp,target.npt-1);
  
  return(out);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  plg_t plg_sample_rigid(const char *filename, plg_t target, criteria_t & criteria,  int equalize)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  plg_t *tmp,out;
  int concat[2];
  double km2m=1000.0, radius=1000.0;
  plg_t *shorelines, *polygones;
  int format=PLG_FORMAT_UNKNOWN, nshorelines, npolygones;
//  int equalize;
  double t1, p1, t2, p2;
  
  npolygones=2;
  polygones=new plg_t[npolygones];
  
  equalize=1;
  
  status=plg_load(filename, format, &shorelines, &nshorelines);
  
  t1 = -4.3488908;
  p1 = 47.845745;
  t2 = -4.1802444;
  p2 = 48.687477;
  
  polygones[0]=plg_sample_shorelines(t1, p1, t2, p2, shorelines, nshorelines, radius, equalize);
  
  t1 = -5.1347866;
  p1 = 48.454182;
  t2 = -5.1347866;
  p2 = 48.454182;
  
  polygones[1]=plg_sample_shorelines(t1, p1, t2, p2, shorelines, nshorelines, radius, equalize);
  
  status=plg_save ("test.plg", PLG_FORMAT_SCAN,  polygones, npolygones);
  
  return(out);
  
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  vector<int> plg_select(frame_t frame, plg_t *shorelines, int nshorelines)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  vector<int> found;
  
  for (size_t s=0;s<nshorelines;s++) {
    bool keep=true;
    for(size_t k=0;k<shorelines[s].npt;k++) {
      double t=geo_recale(shorelines[s].t[k],frame.xmin+180,180.);
      double a=t-frame.xmin;
      double b=t-frame.xmax;
      if(a*b>0.0) {
        keep=false;
        break;
        }
      a=shorelines[s].p[k]-frame.ymin;
      b=shorelines[s].p[k]-frame.ymax;
      if(a*b>0.0) {
        keep=false;
        break;
        }
      }
    if(keep) found.push_back(s);
    }
  
  return(found);
  
}
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  plg_t plg_sample_islands(const char *filename, plg_t target, criteria_t & criteria,  int equalize)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int count, status;
  plg_t *tmp,out;
  double km2m=1000.0, radius=500.0;
  plg_t *shorelines, *polygones;
  int format=PLG_FORMAT_UNKNOWN, nshorelines, npolygones;
//  int equalize;
  double t1, p1, t2, p2;
  
  frame_t frame;
  
  frame.xmin=-5.5;
  frame.xmax=-3.5;
  
  frame.ymin= 47.5;
  frame.ymax= 49.0;
  
  equalize=1;
  
  status=plg_load(filename, format, &shorelines, &nshorelines);
  
  vector<int> list=plg_select(frame, shorelines, nshorelines);
  
  for(int s=0;s<list.size();s++) {
    size_t n=list[s];
    if(plg_isclosed(shorelines[n])==0) {
      list.erase(list.begin()+s);
      s--;
      }
    }
    
  npolygones=list.size();
  polygones=new plg_t[npolygones];
    
  for(size_t s=0;s<npolygones;s++) {
    size_t n=list[s];
    t1 = shorelines[n].t[0];
    p1 = shorelines[n].p[0];
    t2 = t1;
    p2 = p1;
    polygones[s]=plg_sample_shorelines(t1, p1, t2, p2, &(shorelines[n]), 1, radius, equalize);
    }
  
  status=plg_save ("test.plg", PLG_FORMAT_SCAN,  polygones, npolygones);
  
  return(out);
  
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_decimate(plg_t *polygone, int decimation, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  aimed to decimate polygon with given increment (decimation parameter)
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  if(decimation<2) return(-1);
  
  int i,j0;
  double j;
  double *p,*t,*x,*y;
  
  ///Ensure you have first and last point
  const int n=max(polygone->npt/decimation,2);
  const double d=(double)(polygone->npt-1)/(n-1);
  
  p=new double[n];
  t=new double[n];
  x=new double[n];
  y=new double[n];
  
  for(i=0,j=0.;i<n;i++,j+=d){
    j0=j;
    p[i]=polygone->p[j0];
    t[i]=polygone->t[j0];
    x[i]=polygone->x[j0];
    y[i]=polygone->y[j0];
    }
  
  delete[]polygone->p;
  delete[]polygone->t;
  delete[]polygone->x;
  delete[]polygone->y;
  
  polygone->npt=n;
  polygone->p=p;
  polygone->t=t;
  polygone->x=x;
  polygone->y=y;
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_decimate(plg_t & polygon, double threshold, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  aimed to removed quasi-aligned consecutive points (with respect to a given
  threshold)
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int status;
  double epsilon=threshold;
  statistic_t s;
  int i,j;
  double *p,*t,*x,*y,*length,mask=-1;
  plg_t tmp;
  int verbose=0;
  
  if(threshold < 0.0) return(0);
  if(polygon.npt < 3) return(0);
  
  length=new double[polygon.npt-1];
  for(j=0;j<polygon.npt-1;j++) {
    length[j]=plg_length(polygon, j, j+1, PLG_SPHERICAL)*1000.0;
    }
  
  if(debug) verbose=1;
  
  if(debug) printf("polygon segments statistics (m):\n");
  s=get_statistics(length, mask, polygon.npt-1, verbose);
    
/*------------------------------------------------------------------------------
  initialise tmp as a copy of polygon */ 
  tmp.duplicate(polygon);
    
  j=1;
  i=j+1;
  while(i<polygon.npt) {
/*------------------------------------------------------------------------------
    checking angle of lines j-1 and j i.e. [j-1 -> j] and [j -> j+1] */ 
    bool discard, flat, last;
    
    double angle=tmp.angle(j);
    
    if(fabs(angle*r2d)<=epsilon)
      discard=true;
    else 
      discard=false;
    
    if(discard) {
/*------------------------------------------------------------------------------
      if aligned lines, point j discarded: */ 
      tmp.t[j]=tmp.t[j+1];
      tmp.p[j]=tmp.p[j+1];
      tmp.x[j]=tmp.x[j+1];
      tmp.y[j]=tmp.y[j+1];
      }
    else {
/*------------------------------------------------------------------------------
      else, point j kept */
      j++;
      }
    if(i==polygon.npt-1) {
      tmp.npt=j+1; 
      break;
      }
/*------------------------------------------------------------------------------
    initialise j+1 from i */ 
    i++;
    tmp.t[j+1]=polygon.t[i];
    tmp.p[j+1]=polygon.p[i];
    tmp.x[j+1]=polygon.x[i];
    tmp.y[j+1]=polygon.y[i];
    if(i==polygon.npt-1) {
      tmp.npt=j+2; 
      break;
      }
    }
    
//   tmp.npt=j+2;    
  
  polygon.destroy();
  
  polygon.duplicate(tmp);
  
  tmp.destroy();
  
  for(j=0;j<polygon.npt-1;j++) {
    length[j]=plg_length(polygon, j, j+1, PLG_SPHERICAL)*1000.0;
    }
  s=get_statistics(length, mask, polygon.npt-1, verbose);
  
  delete[] length;
  
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_decimate(vector<plg_t> & polygons, double threshold, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int bkp, status;
  
  for(int p=0;p<polygons.size();p++) {
    bkp=polygons[p].npt;
    status=plg_decimate(polygons[p],threshold,debug);
    printf("polygon %3d size: original=%5d decimated=%5d\n", p, bkp, polygons[p].npt);
    }
    
  return(0);
}
