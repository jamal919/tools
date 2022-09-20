
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

\brief poc_data_t operators
*/
/*----------------------------------------------------------------------------*/


#include "poc-data-operators.hpp"


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  complex<double> complex_operator_lt(const complex<double> &x1,const complex<double> &x2)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  complex<double> result;
  
  result=real(x1)<real(x2);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  poc_cdata_t operator<(const poc_cdata_t & d1,const poc_cdata_t & d2)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  poc_cdata_t d0;
  
  implement_operator(&d0,d1,d2,"<",complex_operator_lt);
  
  return d0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  complex<double> complex_operator_le(const complex<double> &x1,const complex<double> &x2)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  complex<double> result;
  
  result=real(x1)<=real(x2);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  poc_cdata_t operator<=(const poc_cdata_t & d1,const poc_cdata_t & d2)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  poc_cdata_t d0;
  
  implement_operator(&d0,d1,d2,"<=",complex_operator_le);
  
  return d0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  complex<double> complex_operator_gt(const complex<double> &x1,const complex<double> &x2)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  complex<double> result;
  
  result=real(x1)>real(x2);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  poc_cdata_t operator>(const poc_cdata_t & d1,const poc_cdata_t & d2)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  poc_cdata_t d0;
  
  implement_operator(&d0,d1,d2,">",complex_operator_gt);
  
  return d0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  complex<double> complex_operator_ge(const complex<double> &x1,const complex<double> &x2)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  complex<double> result;
  
  result=real(x1)>=real(x2);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  poc_cdata_t operator>=(const poc_cdata_t & d1,const poc_cdata_t & d2)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  poc_cdata_t d0;
  
  implement_operator(&d0,d1,d2,">=",complex_operator_ge);
  
  return d0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  complex<double> complex_operator_eq(const complex<double> &x1,const complex<double> &x2)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  complex<double> result;
  
  result=real(x1)==real(x2);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  poc_cdata_t operator==(const poc_cdata_t & d1,const poc_cdata_t & d2)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  poc_cdata_t d0;
  
  implement_operator(&d0,d1,d2,"==",complex_operator_eq);
  
  return d0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  complex<double> complex_operator_ne(const complex<double> &x1,const complex<double> &x2)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  complex<double> result;
  
  result=real(x1)!=real(x2);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  poc_cdata_t operator!=(const poc_cdata_t & d1,const poc_cdata_t & d2)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  poc_cdata_t d0;
  
  implement_operator(&d0,d1,d2,"!=",complex_operator_ne);
  
  return d0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  poc_cdata_t operator+(const poc_cdata_t & d1,const poc_cdata_t & d2)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  poc_cdata_t d0;
  
  implement_operator(&d0,d1,d2,"+",operator_add<complex<double>,complex<double>,complex<double> >);
  
  return d0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  poc_cdata_t operator-(const poc_cdata_t & d1,const poc_cdata_t & d2)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  poc_cdata_t d0;
  
  implement_operator(&d0,d1,d2,"-",operator_subtract<complex<double>,complex<double>,complex<double> >);
  
  return d0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  poc_cdata_t operator-(const poc_cdata_t & d1)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  poc_cdata_t d0;
  
  implement_unary(&d0,d1,"-",operator_minus<complex<double> >);
  
  return d0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  poc_cdata_t operator*(const poc_cdata_t & d1,const poc_cdata_t & d2)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  poc_cdata_t d0;
  
  implement_operator(&d0,d1,d2,"*",operator_mult<complex<double>,complex<double>,complex<double> >);
  
  return d0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  poc_cdata_t operator%(const poc_cdata_t & d1,const poc_cdata_t & d2)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  poc_cdata_t d0;
  
  implement_operator(&d0,d1,d2,"%",operator_scalar<complex<double>,complex<double>,complex<double> >);
  
  return d0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  poc_cdata_t operator/(const poc_cdata_t & d1,const poc_cdata_t & d2)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  poc_cdata_t d0;
  
  implement_operator(&d0,d1,d2,"/",operator_divide<complex<double>,complex<double>,complex<double> >);
  
  return d0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  poc_cdata_t operator^(const poc_cdata_t & d1,const poc_cdata_t & d2)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  poc_cdata_t d0;
  
  implement_operator(&d0,d1,d2,"^",operator_pow<complex<double>,complex<double>,complex<double> >);
  
  return d0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void get2dimLen(const poc_cdata_t & d,int *l0,int *l1)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// row concatenation
/*----------------------------------------------------------------------------*/
{
  const poc_list_t<poc_dim_t> &dimensions=d.info.dimensions;
  const int dl=dimensions.size();
  
  switch(dl){
  case 0:
  case 1:
    *l0=d.length;
    *l1=1;
    break;
  case 2:
    *l0=dimensions[0].len;
    *l1=dimensions[1].len;
    break;
  default:
    TRAP_ERR_EXIT(ENOEXEC,"%s should not have been called with %d dimensions\n",__FUNCTION__,dl);
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  poc_cdata_t operator&(const poc_cdata_t & d1,const poc_cdata_t & d2)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// column concatenation
/*----------------------------------------------------------------------------*/
{
  poc_cdata_t d0;
  
  const int
    d1l=d1.info.dimensions.size(),
    d2l=d2.info.dimensions.size();
  
  if(2<d1l || 2<d2l) TRAP_ERR_RETURN(d0,1,"%s:at least 1 wrong number of dimensions: %d or %d\n",__FUNCTION__,d1l,d2l);
  
  int d0l0,d0l1,d1l0,d1l1,d2l0,d2l1;
  get2dimLen(d1,&d1l0,&d1l1);
  get2dimLen(d2,&d2l0,&d2l1);
  
  if(d1l1==1 && d2l1==1){
    swapValues(&d1l0,&d1l1);
    swapValues(&d2l0,&d2l1);
    }
  
//   if(d1l0!=d2l0 && d1.length!=1 && d2.length!=1) TRAP_ERR_RETURN(d0,1,"%s:first dimension mismatch: %dx%d and %dx%d\n",__FUNCTION__,d1l0,d1l1,d2l0,d2l1);
  
  d0l0=max(d1l0,d2l0);
  
  int i,j,k,l,*l1,*l2,zero=0;
  
  if(d1.length==1){
    l1=&zero;
    d1l0=d0l0;
    }
  else
    l1=&l;
  
  if(d2.length==1){
    l2=&zero;
    d2l0=d0l0;
    }
  else
    l2=&l;
  
  const string
    formula=getNameOrFormula(d1.info)+"&"+getNameOrFormula(d2.info);
  d0.info<<poc_att_t("formula",formula);
  
  d0l1=d1l1+d2l1;
  
  if(d0l0==1)
    d0.info<<poc_dim_t("n",d0l1);
  else
    d0.info<<poc_dim_t("nj",d0l0)<<poc_dim_t("ni",d0l1);
  
  d0.init();
  d0.info.type=d1.info.type;
  
  for(j=0;j<d0l0;j++){
    k=j*d0l1;
    
    if(j<d1l0){
      l=j*d1l1;
      for(i=0;i<d1l1;i++,k++,l++)
        d0.data[k]=d1.data[*l1];
      }
    else
      for(i=0;i<d1l1;i++,k++)
        d0.data[k]=NAN;
    
    if(j<d2l0){
      l=j*d2l1;
      for(i=0;i<d2l1;i++,k++,l++)
        d0.data[k]=d2.data[*l2];
      }
    else
      for(i=0;i<d2l1;i++,k++)
        d0.data[k]=NAN;
    }
  
  return d0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  poc_cdata_t operator|(const poc_cdata_t & d1,const poc_cdata_t & d2)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// row concatenation
/*----------------------------------------------------------------------------*/
{
  poc_cdata_t d0;
  
  const int
    d1l=d1.info.dimensions.size(),
    d2l=d2.info.dimensions.size();
  
  if(2<d1l || 2<d2l) TRAP_ERR_RETURN(d0,1,"%s:at least 1 wrong number of dimensions: %d or %d\n",__FUNCTION__,d1l,d2l);
  
  int d0l0,d0l1,d1l0,d1l1,d2l0,d2l1;
  get2dimLen(d1,&d1l0,&d1l1);
  get2dimLen(d2,&d2l0,&d2l1);
  
  if(d1l1==1)
    swapValues(&d1l0,&d1l1);
  
  if(d2l1==1)
    swapValues(&d2l0,&d2l1);
  
//   if(d1l1!=d2l1 && d1.length!=1 && d2.length!=1) TRAP_ERR_RETURN(d0,1,"%s:second dimension mismatch: %dx%d and %dx%d\n",__FUNCTION__,d1l0,d1l1,d2l0,d2l1);
  
  d0l1=max(d1l1,d2l1);
  
  int i,j,k,l=0,*k1,*k2,zero=0;
  
  if(d1.length==1){
    k1=&zero;
    d1l1=d0l1;
    }
  else
    k1=&k;
  
  if(d2.length==1){
    k2=&zero;
    d2l1=d0l1;
    }
  else
    k2=&k;
  
  const string
    formula=getNameOrFormula(d1.info)+"|"+getNameOrFormula(d2.info);
  d0.info<<poc_att_t("formula",formula);
  
  d0l0=d1l0+d2l0;
  d0.info<<poc_dim_t("nj",d0l0)<<poc_dim_t("ni",d0l1);
  d0.init();
  
  for(j=0,k=0;j<d1l0;j++){
    for(i=0;i<d1l1;i++,k++,l++)
      d0.data[l]=d1.data[*k1];
    for(;i<d0l1;i++,l++)
      d0.data[l]=NAN;
    }
  
  for(j=0,k=0;j<d2l0;j++){
    for(i=0;i<d2l1;i++,k++,l++)
      d0.data[l]=d2.data[*k2];
    for(;i<d0l1;i++,l++)
      d0.data[l]=NAN;
    }
  
  return d0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  poc_cdata_t operator&&(const poc_cdata_t & d1,const poc_cdata_t & d2)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  poc_cdata_t d0;
  
  const int
    d1l=d1.info.dimensions.size(),
    d2l=d2.info.dimensions.size();
  
  if(d1l!=d2l) TRAP_ERR_RETURN(d0,1,"%s:dimension mismatch: %d!=%d\n",__FUNCTION__,d1l,d2l);
  if(d1.length!=d2.length)
    TRAP_ERR_RETURN(d0,1,"%s:size mismatch: %d!=%d\n",__FUNCTION__,d1.length,d2.length);
  
  const string
    formula=getNameOrFormula(d1.info)+"&&"+getNameOrFormula(d2.info);
  d0.info<<poc_att_t("formula",formula);
  
  d0.info.dimensions=d1.info.dimensions;
  
  d0.init();
  
  for(int i=0;i<d0.length;i++){
    const complex<double>
      *d1i=&d1.data[i],
      *d2i=&d2.data[i];
    complex<double> *d0i=&d0.data[i];
    
    if(*d1i==d1.mask)
      *d0i=NC_FILL_COMPLEX;
    else{
      if(*d2i==d2.mask)
        *d0i=NC_FILL_COMPLEX;
      else
        *d0i=*d2i;
      }
    }
  
  return d0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  poc_cdata_t operator||(const poc_cdata_t & d1,const poc_cdata_t & d2)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,zero=0,*i2=&i;
  poc_cdata_t d0;
  
  const int
    d1l=d1.info.dimensions.size(),
    d2l=d2.info.dimensions.size();
  
  if(d2l!=0){
    if(d1l!=d2l) TRAP_ERR_RETURN(d0,1,"%s:dimension mismatch: %d!=%d\n",__FUNCTION__,d1l,d2l);
    if(d1.length!=d2.length)
      TRAP_ERR_RETURN(d0,1,"%s:size mismatch: %d!=%d\n",__FUNCTION__,d1.length,d2.length);
    }
  else{
    i2=&zero;
    }
  
  const string
    formula=getNameOrFormula(d1.info)+"||"+getNameOrFormula(d2.info);
  d0.info<<poc_att_t("formula",formula);
  
  d0.info.dimensions=d1.info.dimensions;
  d0.init();
  
  for(i=0;i<d0.length;i++){
    const complex<double>
      *d1i=&d1.data[i],
      *d2i=&d2.data[*i2];
    complex<double> *d0i=&d0.data[i];
    
    if(*d1i==d1.mask){
      if(*d2i==d2.mask)
        *d0i=NC_FILL_COMPLEX;
      else
        *d0i=*d2i;
      }
    else
      *d0i=*d1i;
    }
  
  return d0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  poc_cdata_t semicolon_operator(const poc_cdata_t & d1,const poc_cdata_t & d2)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// reverse index
/*----------------------------------------------------------------------------*/
{
  poc_cdata_t d0;
  const int d2l=d2.info.dimensions.size();
  
  const string
    formula=getNameOrFormula(d1.info)+':'+getNameOrFormula(d2.info);
  
  if(d2l!=1) TRAP_ERR_RETURN(d0,1,"%s("+formula+"):wrong number of dimensions: %d!=1\n",__FUNCTION__,d2l);
  
  d0.info<<poc_att_t("formula",formula);
  d0.info.dimensions=d1.info.dimensions;
  d0.init();
  
  int m,k0=-1;
  
  double *indexes,*data;
  indexes=new double[d2.length];
  data=new double[d2.length];
  const double mask=real(d2.mask);
  
  for(m=0;m<d2.length;m++){
    indexes[m]=m;
    data[m]=real(d2.data[m]);
    }
  
  for(m=0;m<d1.length;m++)
    map_interpolate1D(indexes,data,mask,d2.length,real(d1.data[m]),&d0.data[m],1,&k0);
  
  delete[]indexes;
  delete[]data;
  
  return d0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void poc_data_index(poc_cdata_t * output,const poc_cdata_t & array,const poc_cdata_t & index)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  /* When using poc_cdata_t::operator[const poc_cdata_t &], the result was
  destroyed before being used for the assignement, causing crashes.
  So this function has been defined instead. */
  
  const poc_list_t<poc_dim_t> &dims=array.info.dimensions;
  const int ndims=dims.size();
  
  int nid,njd;
  switch(ndims){
  case 0:
    *output=array;
    return;
  case 1:
    nid=dims[0].len;
    njd=1;
    break;
  case 2:
    njd=dims[0].len;
    nid=dims[1].len;
    break;
  default:
    TRAP_ERR_EXIT(ENOEXEC,"%s not coded yet for %d dimensions",__FUNCTION__,ndims);
    }
  
  output->info.attributes.clear();
  const string
    formula=getNameOrFormula(array.info)+'['+getNameOrFormula(index.info)+']';
  output->info<<poc_att_t("formula",formula);
  output->info.dimensions=index.info.dimensions;
  output->init();
  
  int mi;          /* index of index.data and output->data */
  int id,jd,md;       /* indexes of array.data */
  int di,dj;
  double sumw,w,wj,wi;     /* weights */
  complex<double> sum;
  
  for(mi=0;mi<index.length;mi++){
    const complex<double> *idmi=&index.data[mi];
    complex<double> *odmi=&output->data[mi];
    const double
      idmir=real(*idmi),
      idmii=imag(*idmi);
    
    if(isnan(*idmi) || idmir<0 || nid-1<idmir || idmii<0 || njd-1<idmii){
      *odmi=NAN;
      continue;
      }
    
    sum=0;
    sumw=0;
    jd=idmii;
    for(dj=0;dj<=1;dj++,jd++){
      if(jd>=njd) continue;
      
      id=idmir;
      md=jd*nid+idmir;
      wj=1-fabs(idmii-jd);
      
      for(di=0;di<=1;di++,id++,md++){
        if(id>=nid) continue;
        
        const complex<double> *admd=&array.data[md];
        if(*admd==array.mask) continue;
        
        wi=1-fabs(idmir-id);
        w=wi*wj;
        
        if(w==0.) continue;
        
        sum+=*admd*w;
        sumw+=w;
        }
      }
    
    if(sumw<.5){
      *odmi=NAN;
      continue;
      }
    
    *odmi=sum/sumw;
    }
  
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void transpose(poc_deque_t<poc_cdata_t*> &vars,poc_cdata_t * output,const poc_cdata_t & input)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  POC_FORMULA_FUNC_HELP("transpose a matrix\n");
  
  const poc_list_t<poc_dim_t> &dimensions=input.info.dimensions;
  const int dl=dimensions.size();
  
  if(dl==1){
    *output=input;
    return;
    }
  
  if(dl!=2) TRAP_ERR_RETURN(,1,"%s:%d dimensions when 1 or 2 are expected\n",__FUNCTION__,dl);
  
  int dl0,dl1;
  get2dimLen(input,&dl0,&dl1);
  
  output->info.dimensions<<dimensions[1]<<dimensions[0];
  output->init();
  
  int i,j,k=0,l;
  
  for(j=0;j<dl0;j++){
    for(i=0;i<dl1;i++,k++){
      l=i*dl0+j;
      output->data[l]=input.data[k];
      }
    }
}
