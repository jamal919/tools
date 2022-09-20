
/*******************************************************************************

  T-UGO tools, 2006-2011

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

\author  Florent Lyard      LEGOS/CNRS, Toulouse, France. florent.lyard@legos.obs-mip.fr
\author  Laurent Roblou     LEGOS/CNRS, Toulouse, France
\author  Damien Allain      LEGOS/CNRS, Toulouse, France
\author  Yves Soufflet      LEGOS/CNRS, Toulouse, France
\author  Clément Mayet      LEGOS, Toulouse, France (PhD)
\author  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
\author  Frédéric Dupont    Université de Laval à Québec, Canada

\brief poc_data_t functions
*/
/*----------------------------------------------------------------------------*/


#include <math.h>

#include "poc-data-operators.hpp"


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void isMasked(poc_deque_t<poc_cdata_t*> &vars,poc_cdata_t * output,const poc_cdata_t & input)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  POC_FORMULA_FUNC_HELP("whether the value is masked\n");
  
  output->info.dimensions=input.info.dimensions;
  output->init(input.axes);
  
  for(int i=0;i<input.length;i++){
    output->data[i]= input.data[i]==input.mask;
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void zeroIfNan(poc_deque_t<poc_cdata_t*> &vars,poc_cdata_t * output,const poc_cdata_t & input)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  POC_FORMULA_FUNC_HELP("copy input, setting masked and nan values to 0.\n");
  
  output->info.dimensions=input.info.dimensions;
  output->init(input.axes);
  
  for(int i=0;i<input.length;i++){
    const complex<double> *ii=&input.data[i];
    complex<double> *oi=&output->data[i];
    if(isnan(*ii) or *ii==input.mask)
      *oi=0.;
    else
      *oi=*ii;
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void nearest(poc_deque_t<poc_cdata_t*> &vars,poc_cdata_t * output,const poc_cdata_t & input)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// wrapper for distance_to_nearest_unmasked()
/*----------------------------------------------------------------------------*/
{
  POC_FORMULA_FUNC_HELP( "value of nearest unmasked node that has a value (not NAN), assuming the grid is orthonormal.\n"
    "  "+funcName+"(z[,doDistance])\n"
    "   doDistance : if different from 0, the default, return distance in pixels to nearest unmasked node, assuming the grid is orthonormal\n"
    );
  
  const poc_list_t<poc_dim_t> &dimensions=input.info.dimensions;
  const int dl=dimensions.size();
  if(dl!=2) TRAP_ERR_RETURN(,1,"%s:%d dimensions when 2 are expected\n",__func__,dl);
  
  grid_t grid;
  grid.ny=dimensions[0].len,
  grid.nx=dimensions[1].len;
  
  double *distance;
  distance=new double[input.length];
  
  output->info.dimensions=dimensions;
  output->init(input.axes);
  output->mask=input.mask;
  
/*------------------------------------------------------------------------------
  optional argument */
  bool doDistance=false;
  
  string name;
  const poc_att_t *att;
  att=input.info.attributes.findP("");
  if(att!=0)
    name=att->as_string();
  
  if(att!=0){
    char *c=&name[name.size()-1];
    
    const poc_cdata_t
      *var=poc_formula_check_var(vars,name);
    const complex<double>
      *var0=&var->data[0];
    
    if(not isnan(*var0) and *var0!=0.)
      doDistance=true;
    
    if(*c>'a'){
      STDERR_BASE_LINE_FUNC("NOTE:all %d last arguments are ignored\n",*c-'a');
      }
    
    }
  
/*----------------------------------------------------------------------------*/
  if(doDistance){
    distance_to_nearest_unmasked(grid,input.data,input.mask,distance);
    valcpy(output->data,distance,output->length);
    }
  else{
    distance_to_nearest_unmasked(grid,input.data,input.mask,distance,output->data);
    }
  
  delete[]distance;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void sum(poc_deque_t<poc_cdata_t*> &vars,poc_cdata_t * output,const poc_cdata_t & input)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  POC_FORMULA_FUNC_HELP("sum of all the unmasked elements\n");
  
  int i;
  complex<double> s=0;
  
  for(i=0;i<input.length;i++){
    complex<double> *idi=&input.data[i];
    if(*idi==input.mask or isnan(*idi)) continue;
    s+=*idi;
    }
  
  *output=s;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void average(poc_deque_t<poc_cdata_t*> &vars,poc_cdata_t * output,const poc_cdata_t & input)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  POC_FORMULA_FUNC_HELP("average of all the unmasked elements\n");
  
  int i,c=0;
  complex<double> s=0;
  
  for(i=0;i<input.length;i++){
    complex<double> *idi=&input.data[i];
    if(*idi==input.mask or isnan(*idi)) continue;
    s+=*idi;
    c++;
    }
  
  *output=s*(1./c);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  complex<double> abs_toComplex(const complex<double> &z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  complex<double> a;
  
  a=abs(z);
  
  return a;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  complex<double> arg_toComplex(const complex<double> &z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  complex<double> a;
  
  a=arg(z);
  
  return a;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  complex<double> sign_toComplex(const complex<double> &z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  const complex<double>
    a(sign(real(z)),sign(imag(z)));
  
  return a;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  complex<double> real_toComplex(const complex<double> &z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  complex<double> a;
  
  a=real(z);
  
  return a;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  complex<double> imag_toComplex(const complex<double> &z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  complex<double> a;
  
  a=imag(z);
  
  return a;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  complex<double> floor_toComplex(const complex<double> &z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  complex<double> a(floor(real(z)),floor(imag(z)));
  
  return a;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  complex<double> round_toComplex(const complex<double> &z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  complex<double> a(round(real(z)),round(imag(z)));
  
  return a;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  complex<double> ceil_toComplex(const complex<double> &z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  complex<double> a(ceil(real(z)),ceil(imag(z)));
  
  return a;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  complex<double> isnan_toComplex(const complex<double> &z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  complex<double> a=isnan(z);
  
  return a;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  complex<double> nanIf0(const complex<double> &z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  complex<double> a;
  
  if(z==0.)
    a=NAN;
  else
    a=z;
  
  return a;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void dimmax(poc_deque_t<poc_cdata_t*> &vars,poc_cdata_t * output,const poc_cdata_t & input)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  POC_FORMULA_FUNC_HELP("maximum along the " POC_DIM_INDEX_VAR_NAME "-th dimension\n");
  
  const poc_cdata_t
    *dimIndexVarP=poc_formula_check_var(vars,POC_DIM_INDEX_VAR_NAME);
  double dimIndexD0=NAN;
  if(dimIndexVarP!=0)
    dimIndexD0=real(dimIndexVarP->data[0]);
  int i,j,dimIndex=-1;
  if(isfinite(dimIndexD0))
    dimIndex=dimIndexD0;
  STDERR_BASE_LINE_FUNC("%d\n",dimIndex);
  
  poc_list_t<poc_dim_t> &oDs=output->info.dimensions;
  const poc_list_t<poc_dim_t> &iDs=input.info.dimensions;
  
  if(dimIndex>=0)
    for(i=0;i<iDs.size();i++){
      if(i==dimIndex) continue;
      oDs<<iDs[i];
      }
  
  output->init();
  aset(output->data,output->length,(complex<double>)-INFINITY);
  
  for(i=0;i<input.length;i++){
    const double ie=real(input.data[i]);/*<input element*/
    
    j=output->translate_index(iDs,i);
    complex<double> *oe=&output->data[j];/*<output element*/
    
    if(ie==input.mask) continue;
    
    if(real(*oe)<ie)
      *oe=ie;
    
    }
  
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void grad_xORy(poc_deque_t<poc_cdata_t*> &vars,poc_cdata_t * output,const poc_cdata_t & input,char axis)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,j,m;
  int m0,m1;
  
  /* check whether dummy run */
  if(input.length<=1) return;
  
  poc_cdata_t
    &lon=poc_formula_get_var(vars,POC_LON_VAR_NAME),
    &lat=poc_formula_get_var(vars,POC_LAT_VAR_NAME);
  
  const poc_list_t<poc_dim_t> &iDs=input.info.dimensions;
  output->info.dimensions=iDs;
  output->init(input.axes);
  
  const int nd=iDs.size();
  int nx,ny;
  bool skipDistance=false;
  
  switch(nd){
  case 1:
    ny=1,
    nx=iDs[0].len;
    skipDistance=true;
    break;
  case 2:
    ny=iDs[0].len,
    nx=iDs[1].len;
    if( ny!=lon.info.dimensions[0].len or
        nx!=lon.info.dimensions[1].len )
      skipDistance=true;
    break;
  default:
    TRAP_ERR_EXIT(ENOEXEC,"%s not coded for %d dimensions\n",__func__,nd);
    }

  axis=tolower(axis);
  
  for(j=0;j<ny;j++){
    
    switch(axis){
    case 'x':
      m=j*nx;
      m0=m-1;
      m1=m+1;
      break;
    case 'y':
      m=j*nx;
      m0=m-nx;
      m1=m+nx;
      break;
    default:
      TRAP_ERR_EXIT(ENOEXEC,"coding error, %s should only be called with axis as x or y, not %c\n",__func__,axis);
      }
    
    for(i=0;i<nx;i++,m++,m0++,m1++){
      complex<double> *oe=&output->data[m];
      
      if( ( axis=='x' and ( i<=0 or i>=nx-1 ) ) or
          ( axis=='y' and ( j<=0 or j>=ny-1 ) ) ) {
        *oe=output->mask;
        continue;
        }
      
      complex<double>
        *ie0=&input.data[m0],
        *ie1=&input.data[m1];
      
      if(*ie0==input.mask || *ie1==input.mask){
        *oe=NC_FILL_COMPLEX;
        continue;
        }
      
      *oe=*ie1-*ie0;
      
      if(skipDistance)
        continue;
      
      const double
        lo0=real(lon.data[m0]),
        lo1=real(lon.data[m1]),
        la0=real(lat.data[m0]),
        la1=real(lat.data[m1]);
      
      double
        a=geo_distance(lo0,la0,lo1,la1);
      
      if( ( axis=='x' and lo1<lo0 ) or
          ( axis=='y' and la1<la0 ) )
        a*=-1.;
      
      *oe/=a;
      }
    }
  
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void gradx(poc_deque_t<poc_cdata_t*> &vars,poc_cdata_t * output,const poc_cdata_t & input)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  POC_FORMULA_FUNC_HELP("X-axis derivative, in meters\n");
  
  grad_xORy(vars,output,input,'x');
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void grady(poc_deque_t<poc_cdata_t*> &vars,poc_cdata_t * output,const poc_cdata_t & input)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  POC_FORMULA_FUNC_HELP("Y-axis derivative, in meters\n");
  
  grad_xORy(vars,output,input,'y');
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int cmp_complex_double(const void *a,const void *b)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  const complex<double>
    *za=(const complex<double> *)a,
    *zb=(const complex<double> *)b;
  
  return cmp(real(*za),real(*zb));
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void ellipseM(poc_deque_t<poc_cdata_t*> &vars,poc_cdata_t * output,const poc_cdata_t & input)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* See also ellipse_Madmp() and ellipse_Madmp3D() */
{
  POC_FORMULA_FUNC_HELP("norm of vector or maximum current of ellipse\n"
    "  "+funcName+"(u[,v[,w]])\n"
    "If u and v are real numbers, this is equivalent to:\n"
    "  (u^2+v^2)^.5\n"
    );
  
/*------------------------------------------------------------------------------
  arguments, including optional ones */
  int i;
  poc_deque_t< complex<double>* > args;
  poc_deque_t< complex<double> > masks;
  args.push_back(input.data);
  masks.push_back(input.mask);
  
  string name;
  const poc_att_t *att;
  att=input.info.attributes.findP("");
  if(att!=0)
    name=att->as_string();
  const string name0=name;
  
  if(att!=0){
    name=name0;
    char *c=&name[name.size()-1];
    ostringstream argLengths;
    bool mismatch=false;
    
    argLengths<<"1:"<<input.length;
    
    for(i=2;i<=3 && *c>='a';i++,(*c)--){
      poc_cdata_t *var;
      var=poc_formula_check_var(vars,name);
      
      args.push_back(var->data);
      
      if(var->length!=input.length)
        mismatch=true;
      argLengths<<","<<i<<":"<<var->length;
      }
    
    if(mismatch)
      TRAP_ERR_RETURN(,1,"%s:argument length mismatch: "+argLengths.str()+"\n",__func__);
    }
  
/*----------------------------------------------------------------------------*/
  output->info.dimensions=input.info.dimensions;
  output->init(input.axes);
  
  int m;
  
  for(m=0;m<input.length;m++){
    complex<double> sumS=0,z,z2;
    double sumAS=0,aSumS;
    
    for(i=0;i<args.size();i++){
      z=args[i][m];
      if(z==masks[i])
        break;
      z2=square(z);
      sumS+=z2;
      sumAS+=abs(z2);
      }
    
    if(i<args.size()){
      output->data[m]=output->mask;
      continue;
      }
    
    aSumS=abs(sumS);
    output->data[m]=sqrt((sumAS+aSumS)*.5);
    }
  
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void sorted(poc_deque_t<poc_cdata_t*> &vars,poc_cdata_t * output,const poc_cdata_t & input)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  POC_FORMULA_FUNC_HELP("sorted vector, ignoring nan and masked values\n"
    "  "+funcName+"(x)\n"
    "If x is complex, return an empty array.\n"
    "The sorted array contains none of the nan or masked values.\n"
    );
  
  int run,i,unmasked;
  
  if(input.isComplex())
    TRAP_ERR_RETURN(,1,"%s:complex input\n",__func__);
  
  for(run=0;;run++){
    
    unmasked=0;
    
    /* run==0: count unmasked values */
    /* run==1: fill with unmasked values */
    for(i=0;i<input.length;i++){
      complex<double> *z=&input.data[i];
      
      if(isnan(*z) or *z==input.mask)
        continue;
      
      if(run>=1)
        output->data[unmasked]=*z;
      
      unmasked++;
      }
    
    if(run>=1)
      break;
    
    output->info.dimensions<<poc_dim_t(input.info.name+"_sorted",unmasked);
    output->init();
    }
  
  qsort(output->data,output->length,sizeof(complex<double>),cmp_complex_double);
}
