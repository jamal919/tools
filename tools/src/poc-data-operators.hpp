
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


#if POC_DATA_OPERATORS_HPP == 0
#define POC_DATA_OPERATORS_HPP 1

#include "maths.h" //operator %
#include "constants.h" //MeanEarthRadius
#include "poc-netcdf-data.hpp"
#include "polygones.h"


extern string getNameOrFormula(const poc_var_t & var);


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T0,typename T1>
  void implement_unary(poc_data_t<T0> *d0,
    const poc_data_t<T1> & d1,const string &name,
    T0 (*op)(const T1 &)
    )

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// unary operator
/*----------------------------------------------------------------------------*/
{
  d0->init(d1.info,d1.axes);
  d0->info.name="";
  
  for(int i=0;i<d1.length;i++){
    T0 *d0i=&d0->data[i];
    T1 *d1i=&d1.data[i];
    
    if(*d1i==d1.mask){
      *d0i=d0->mask;
      continue;
      }
    
    *d0i=op(*d1i);
    }
  
  string d1formula=getNameOrFormula(d1.info);
  const bool nameIsVar=isLetter(name[0]);
  if(nameIsVar && d1formula[0]!='(')
    d1formula='('+d1formula+')';
  
  const string
    formula=name+d1formula;
  
  d0->info<<poc_att_t("formula",formula);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T0,typename T1,typename T2,typename To>
  void implement_operator(poc_data_t<T0> *d0,
    const poc_data_t<T1> & d1,const T2 & d2,const string &name,
    To (*op)(const T1 &,const T2 &)
    )

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// binary operator
/*----------------------------------------------------------------------------*/
{
  d0->init(d1.info);
  d0->info.name="";
  
  for(int i=0;i<d1.length;i++){
    double *d0i=&d0->data[i],*d1i=&d1.data[i];
    if(*d1i==d1.mask){
      *d0i=d0->mask;
      continue;
      }
    *d0i=op(*d1i,d2);
    }
  
  ostringstream formula;
  formula<<getNameOrFormula(d1.info);
  formula<<name;
  formula<<d2;
  
  d0->info<<poc_att_t("formula",formula.str());
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T0,typename T1,typename T2,typename To>
  void implement_operator(poc_data_t<T0> *d0,
    const poc_data_t<T1> & d1,const poc_data_t<T2> & d2,const string &name,
    To (*op)(const T1 &,const T2 &)
    )

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
{
  int i,zero=0,translated,*i1=&i,*i2=&i;
  if(d1.length>1 && d2.length>1 && d1.length!=d2.length){
    d0->info.dimensions=d1.info.dimensions;
    d0->info.dimensions<<d2.info.dimensions;
    
    STDERR_BASE_LINE(""+getNameOrFormula(d1.info)+name+getNameOrFormula(d2.info)+":\n");
    poc_print(d0->info,cerr);
    d0->init();
    
    if(d1.length<d0->length) i1=&translated;
    if(d2.length<d0->length) i2=&translated;
    }
  else if(d1.length>d2.length){
    d0->init(d1.info,d1.axes);
    }
  else{
    d0->init(d2.info,d2.axes);
    }
  if(d1.length==1) i1=&zero;
  if(d2.length==1) i2=&zero;
  
  d0->info.name="";
  
  for(i=0;i<d0->length;i++){
    T0 *d0i=&d0->data[i];
    
    if(i1==&translated)
      translated=d1.translate_index(d0->info.dimensions,i);
    T1 *d1i=&d1.data[*i1];
    
    if(i2==&translated)
      translated=d2.translate_index(d0->info.dimensions,i);
    T2 *d2i=&d2.data[*i2];
    
    if(*d1i==d1.mask || *d2i==d2.mask){
      *d0i=d0->mask;
      continue;
      }
    
    *d0i=op(*d1i,*d2i);
    }
  
  const string
    formula=getNameOrFormula(d1.info)+name+getNameOrFormula(d2.info);
  d0->info<<poc_att_t("formula",formula);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename To,typename T1,typename T2>
  To operator_gt(const T1 &x1,const T2 &x2)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  To result;
  
  result=x1>x2;
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T1,typename T2>
  poc_data_t<int8_t> operator>(const poc_data_t<T1> & d1,const T2 & d2)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  poc_data_t<int8_t> d0;
  
  implement_operator(&d0,d1,d2,">",operator_gt<int8_t,T1,T2>);
  
  return d0;
}


extern poc_cdata_t operator< (const poc_cdata_t & d1,const poc_cdata_t & d2);
extern poc_cdata_t operator<=(const poc_cdata_t & d1,const poc_cdata_t & d2);
extern poc_cdata_t operator> (const poc_cdata_t & d1,const poc_cdata_t & d2);
extern poc_cdata_t operator>=(const poc_cdata_t & d1,const poc_cdata_t & d2);
extern poc_cdata_t operator==(const poc_cdata_t & d1,const poc_cdata_t & d2);
extern poc_cdata_t operator!=(const poc_cdata_t & d1,const poc_cdata_t & d2);


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename To,typename T1,typename T2>
  To operator_add(const T1 &x1,const T2 &x2)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  To result;
  
  result=x1+x2;
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T>
  poc_data_t<double> operator+(const poc_data_t<T> & d1,const double & d2)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  poc_data_t<double> d0;
  
  implement_operator(&d0,d1,d2,"+",operator_add<double,T,double>);
  
  return d0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T>
  poc_data_t<double> operator+(const poc_data_t<double> & d1,const T & d2)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  poc_data_t<double> d0;
  
  implement_operator(&d0,d1,d2,"+",operator_add<double,double,T>);
  
  return d0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T>
  poc_data_t<double> operator+(const poc_data_t<T> & d1,const poc_data_t<double> & d2)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  poc_data_t<double> d0;
  
  implement_operator(&d0,d1,d2,"+",operator_add<double,T,double>);
  
  return d0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T>
  poc_data_t<double> operator+(const poc_data_t<double> & d1,const poc_data_t<T> & d2)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  poc_data_t<double> d0;
  
  implement_operator(&d0,d1,d2,"+",operator_add<double,double,T>);
  
  return d0;
}


extern poc_cdata_t operator+(const poc_cdata_t & d1,const poc_cdata_t & d2);


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename To,typename T1,typename T2>
  To operator_subtract(const T1 &x1,const T2 &x2)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  To result;
  
  result=x1-x2;
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T1,typename T2>
  poc_data_t<double> operator-(const poc_data_t<T1> & d1,const T2 & d2)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  poc_data_t<double> d0;
  
  implement_operator(&d0,d1,d2,"-",operator_subtract<double,T1,T2>);
  
  return d0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T1,typename T2>
  poc_data_t<double> operator-(const poc_data_t<T1> & d1,const poc_data_t<T2> & d2)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  poc_data_t<double> d0;
  
  implement_operator(&d0,d1,d2,"-",operator_subtract<double,T1,T2>);
  
  return d0;
}


extern poc_cdata_t operator-(const poc_cdata_t & d1,const poc_cdata_t & d2);


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T> T operator_minus(const T & x)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  T result;
  
  result=-x;
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T> poc_data_t<T> operator-(const poc_data_t<T> & d1)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  poc_data_t<T> d0;
  
  implement_unary(&d0,d1,"-",operator_minus<T>);
  
  return d0;
}


extern poc_cdata_t operator-(const poc_cdata_t & d1);


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename To,typename T1,typename T2>
  To operator_mult(const T1 &x1,const T2 &x2)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  To result;
  
  result=x1*x2;
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T1,typename T2>
  poc_data_t<double> operator*(const poc_data_t<T1> & d1,const T2 & d2)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  poc_data_t<double> d0;
  
  implement_operator(&d0,d1,d2,"*",operator_mult<double,T1,T2>);
  
  return d0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T1,typename T2>
  poc_data_t<double> operator*(const poc_data_t<T1> & d1,const poc_data_t<T2> & d2)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  poc_data_t<double> d0;
  
  implement_operator(&d0,d1,d2,"*",operator_mult<double,T1,T2>);
  
  return d0;
}


extern poc_cdata_t operator*(const poc_cdata_t & d1,const poc_cdata_t & d2);


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename To,typename T1,typename T2>
  To operator_scalar(const T1 &x1,const T2 &x2)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  To result;
  
  result=x1%x2;
  
  return result;
}


extern poc_cdata_t operator%(const poc_cdata_t & d1,const poc_cdata_t & d2);


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename To,typename T1,typename T2>
  To operator_divide(const T1 &x1,const T2 &x2)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  To result;
  
  result=x1/x2;
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T1,typename T2>
  poc_data_t<double> operator/(const poc_data_t<T1> & d1,const T2 & d2)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  poc_data_t<double> d0;
  
  implement_operator(&d0,d1,d2,"/",operator_divide<double,T1,T2>);
  
  return d0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T1,typename T2>
  poc_data_t<double> operator/(const poc_data_t<T1> & d1,const poc_data_t<T2> & d2)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  poc_data_t<double> d0;
  
  implement_operator(&d0,d1,d2,"/",operator_divide<double,T1,T2>);
  
  return d0;
}


extern poc_cdata_t operator/(const poc_cdata_t & d1,const poc_cdata_t & d2);


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename To,typename T1,typename T2>
  To operator_pow(const T1 &x1,const T2 &x2)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  To result;
  
  result=pow(x1,x2);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T1,typename T2>
  poc_data_t<double> operator^(const poc_data_t<T1> & d1,const T2 & d2)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  poc_data_t<double> d0;
  
  implement_operator(&d0,d1,d2,"^",operator_pow<double,T1,T2>);
  
  return d0;
}


extern poc_cdata_t operator^(const poc_cdata_t & d1,const poc_cdata_t & d2);
extern poc_cdata_t operator&(const poc_cdata_t & d1,const poc_cdata_t & d2);
extern poc_cdata_t operator|(const poc_cdata_t & d1,const poc_cdata_t & d2);
extern poc_cdata_t operator&&(const poc_cdata_t & d1,const poc_cdata_t & d2);
extern poc_cdata_t operator||(const poc_cdata_t & d1,const poc_cdata_t & d2);
extern poc_cdata_t semicolon_operator(const poc_cdata_t & d1,const poc_cdata_t & d2);
extern void poc_data_index(poc_cdata_t * output,const poc_cdata_t & input,const poc_cdata_t & index);
extern poc_cdata_t operator<<(const poc_cdata_t & d1,const poc_cdata_t & d2);

extern void var2plg(const poc_cdata_t & var,plg_t *plg);

extern void poc_save_grid(poc_deque_t<poc_cdata_t*> &vars,const grid_t& grid);
extern int get_poc_cdata_dims(const poc_cdata_t & var,int *nx,int *ny);
extern int poc_cdata2grid(grid_t * grid,const poc_cdata_t & var);
extern int poc_save_vars(const string & path,const poc_deque_t<poc_cdata_t*> & vars,int verbose=0,bool saveAsGrid=false,bool overwrite=true);
extern int poc_load_vars(const string & path,poc_deque_t<poc_cdata_t*> *vars,int verbose=0,string *formula=0);

#define POC_DIM_INDEX_VAR_NAME "dimIndex"
#define POC_LAT_VAR_NAME "lat"
#define POC_LON_VAR_NAME "lon"

#define HASH_LINE "################################################################################\n"
#define POC_FORMULA_FUNC_HELP(s, args... ) if(output==0){string funcName=input.info.name;printf(HASH_LINE "# "+input.info.name+"(): " s HASH_LINE, ##args );return;}


extern poc_cdata_t *poc_formula_check_var(const poc_deque_t<poc_cdata_t*> &vars,const string & name, int *index=0);
extern poc_cdata_t & poc_formula_get_var(poc_deque_t<poc_cdata_t*> &vars,const string & name);
extern complex<double> * poc_formula_init_var(poc_deque_t<poc_cdata_t*> &vars,const string & name);

extern poc_cdata_t comma_operator(poc_deque_t<poc_cdata_t*> &vars,const poc_cdata_t & d1,const poc_cdata_t & d2);
extern void poc_formula_keep_consts(poc_deque_t<poc_cdata_t*> &vars);
extern void poc_formula_implement_func(poc_deque_t<poc_cdata_t*> &vars,poc_cdata_t * output,const string & name,const poc_cdata_t & input);

extern int poc_formula_parse(poc_deque_t<poc_cdata_t*> &vars,const char *str);

#endif
