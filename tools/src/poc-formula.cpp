
/*******************************************************************************

  T-UGO tools, 2006-2019

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

\brief formula function definitions
*/
/*----------------------------------------------------------------------------*/


#include <unistd.h> /* for close() */
#include <fcntl.h> /* for O_* */

#include "poc-data-operators.hpp"


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void open(poc_deque_t<poc_cdata_t*> &vars,poc_cdata_t * output,const poc_cdata_t & input)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  POC_FORMULA_FUNC_HELP("open a file in append mode and return a descriptor\n");
  
  int i;
  
  if(input.info.type!=NC_CHAR)
    TRAP_ERR_RETURN(,"%s:argument must be a string\n",__func__);
  
  char *s=new char[input.length+1];
  
  for(i=0;i<input.length;i++)
    s[i]=real(input.data[i]);
  s[i]='\0';
  
  int f;
  f=open(s,O_APPEND|O_CREAT|O_WRONLY,0666);
  delete[]s;
  
  *output=f;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void close(poc_deque_t<poc_cdata_t*> &vars,poc_cdata_t * output,const poc_cdata_t & input)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  POC_FORMULA_FUNC_HELP("close a file descriptor\n");
  
  if(input.info.type==NC_CHAR)
    TRAP_ERR_RETURN(,"%s:file must be a descriptor\n",__func__);
  
  int f;
  f=real(input.data[0]);
  
  close(f);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void write(poc_deque_t<poc_cdata_t*> &vars,poc_cdata_t * output,const poc_cdata_t & input)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  POC_FORMULA_FUNC_HELP("write\n"
    "  count="+funcName+"(file,buffer)\n"
    );
  
  if(input.info.type==NC_CHAR)
    TRAP_ERR_RETURN(,"%s:file must be a descriptor\n",__func__);
  
  const complex<double> *buffer=0;
  size_t c0=0,c1;
  
/*------------------------------------------------------------------------------
  buffer argument */
  string name;
  const poc_att_t *att;
  att=input.info.attributes.findP("");
  if(att!=0)
    name=att->as_string();
  
  if(att!=0){
    char *c=&name[name.size()-1];
    
    const poc_cdata_t
      *var=poc_formula_check_var(vars,name);
    
    buffer=var->data;
    c0=var->length*sizeof(complex<double>);
    
    if(*c>'a'){
      STDERR_BASE_LINE_FUNC("NOTE:all %d last arguments are ignored\n",*c-'a');
      }
    
    }
  
  if(buffer==0) TRAP_ERR_RETURN(,"%s:you must give a buffer argument\n",__func__);
  
/*----------------------------------------------------------------------------*/
  
  int f;
  f=real(input.data[0]);
  
  c1=write(f,buffer,c0);
  
  *output=c1;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  string getNameOrFormula(const poc_var_t & var)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
{
  if(var.name.length())
    return var.name;
  
  const poc_att_t *att=var.attributes.findP("formula");
  if(!att)
    return "()";
  else
    return "("+att->as_string()+")";
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void print(poc_deque_t<poc_cdata_t*> &vars,poc_cdata_t * output,const poc_cdata_t & input)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  POC_FORMULA_FUNC_HELP("print the content of an input as an eventually reusable formula and return that input\n");
  
  int i;
  
  printf(""+getNameOrFormula(input.info));
  
  const poc_list_t<poc_dim_t> &dimensions=input.info.dimensions;
  
  putchar('[');
  for(i=0;i<dimensions.size();i++){
    const poc_dim_t *dim=&dimensions[i];
    if(i) putchar(',');
    printf(""+dim->name+"=%d",dim->len);
    }
  putchar(']');
  
  if(input.info.type==NC_CHAR){
    char *s=new char[input.length];
    
    for(i=0;i<input.length;i++)
      s[i]=real(input.data[i]);
    
    printf("=\"%s\"",s);
    delete[]s;
    }
  else
   for(i=0;i<input.length;i++){
    printf("%c",i?'&':'=');
    
    complex<double> *z=&input.data[i];
    const bool
      hasReal=real(*z)!=0,
      hasImag=imag(*z)!=0;
    
    if(hasReal || not hasImag)
      printf("%.5g",real(*z));
    
    if(hasImag){
      if(hasReal)
        printf("%+.5gj",imag(*z));
      else
        printf("%.5gj",imag(*z));
      }
    }
  
  printf(";\n");
  
  *output=input;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void printInfo(poc_deque_t<poc_cdata_t*> &vars,poc_cdata_t * output,const poc_cdata_t & input)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  POC_FORMULA_FUNC_HELP("print the informations of an input and return that input\n");
  
  poc_print(input.info);
  
  *output=input;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void renameDim(poc_deque_t<poc_cdata_t*> &vars,poc_cdata_t * output,const poc_cdata_t & input)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  POC_FORMULA_FUNC_HELP("return input with dimension renamed\n"
    "  input="+funcName+"(input,newName)\n"
    );
  
  if(input.info.type==NC_CHAR)
    TRAP_ERR_RETURN(,"%s:file must be a descriptor\n",__func__);
  
  string newName;
  
/*------------------------------------------------------------------------------
  buffer argument */
  string name;
  const poc_att_t *att;
  att=input.info.attributes.findP("");
  if(att!=0)
    name=att->as_string();
  
  if(att!=0){
    char *c=&name[name.size()-1];
    
    const poc_cdata_t
      *var=poc_formula_check_var(vars,name);
    
    newName=var->toStr();
    
    if(*c>'a'){
      STDERR_BASE_LINE_FUNC("NOTE:all %d last arguments are ignored\n",*c-'a');
      }
    
    }
  
  if(newName=="") TRAP_ERR_RETURN(,"%s:you must give a non-empty newName argument\n",__func__);
  
/*----------------------------------------------------------------------------*/
  
  *output=input;
  
  output->info.dimensions[0].name=newName;
  output->info.attributes.erase("");
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void dims(poc_deque_t<poc_cdata_t*> &vars,poc_cdata_t * output,const poc_cdata_t & input)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  POC_FORMULA_FUNC_HELP("dimension lengths\n");
  
  int i;
  
  const poc_list_t<poc_dim_t> &iDs=input.info.dimensions;
  const int
    n=iDs.size();
  
  output->info << poc_dim_t("dims",n);
  output->init();
  
  for(i=0;i<n;i++){
    output->data[i]=iDs[i].len;
    }
  
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  poc_cdata_t comma_operator(poc_deque_t<poc_cdata_t*> &vars,const poc_cdata_t & d1,const poc_cdata_t & d2)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// reverse index
/*----------------------------------------------------------------------------*/
{
  poc_cdata_t d0;
  
  const string
    formula=getNameOrFormula(d1.info)+','+getNameOrFormula(d2.info);
  
  static uint64_t uses=0;
  
  d0=d1;
  d0.info.name="";
  d0.info << poc_att_t("formula",formula);
  
  string name;
  
  const poc_att_t *att;
  /* In a comma-separated list of arguments,
  the last comma will be parsed first. */
  att=d2.info.attributes.findP("");
  
  if(att==0){
    struct timeval t;
    gettimeofday(&t);
    
    char *s;
    
    #pragma omp critical(comma_operator_uses)
    {
    asprintf(&s,"%lx,%c",uses,'a'-1);
    if(uses>=0xfffffffffffffff0uL)
      TRAP_ERR_EXIT(ENOEXEC,"line below not coded for more than %lu uses.\n",uses);
    uses++;
    }
  
    name=s;
    free(s);
    }
  else{
    name=att->as_string();
    }
  
  char *c=&name[name.size()-1];
  
  if(*c>='z') TRAP_ERR_EXIT(ENOEXEC,"not coded for more than %d arguments\n",*c-'a'+2);
  (*c)++;
  
#define TRACK_COMMA_STACK 0
#if TRACK_COMMA_STACK
  STDERR_BASE_LINE_FUNC("creating "+name+"\n");
#endif
  poc_cdata_t *stacked=&poc_formula_get_var(vars,name);
  *stacked=d2;
  stacked->info.name=name;
  
  d0.info << poc_att_t("",name);
  
  return d0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void poc_formula_make_const(poc_deque_t<poc_cdata_t*> &vars,const string & name,complex<double> z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  poc_cdata_t *var=&poc_formula_get_var(vars,name);
  
  if(var->length==1 && var->data[0]==z) return;
  
  *var=z;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void poc_formula_keep_consts(poc_deque_t<poc_cdata_t*> &vars)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  poc_formula_make_const(vars,"e", M_E);
  poc_formula_make_const(vars,"nan",NAN);
  poc_formula_make_const(vars,"pi",M_PI);
  
  poc_formula_make_const(vars,"d2r",d2r);
  poc_formula_make_const(vars,"r2d",r2d);
}


map<string,complex<double> (*)(const complex<double> &)> unaryFunctions;
map<string,void (*)(poc_deque_t<poc_cdata_t*> &vars,poc_cdata_t*,const poc_cdata_t&)> fieldFunctions;

extern complex<double> abs_toComplex(const complex<double> &z);
extern complex<double> arg_toComplex(const complex<double> &z);
extern complex<double> sign_toComplex(const complex<double> &z);
extern complex<double> real_toComplex(const complex<double> &z);
extern complex<double> imag_toComplex(const complex<double> &z);
extern complex<double> floor_toComplex(const complex<double> &z);
extern complex<double> round_toComplex(const complex<double> &z);
extern complex<double> ceil_toComplex(const complex<double> &z);
extern complex<double> isnan_toComplex(const complex<double> &z);
extern complex<double> nanIf0(const complex<double> &z);

extern void isMasked(poc_deque_t<poc_cdata_t*> &vars,poc_cdata_t * output,const poc_cdata_t & input);
extern void zeroIfNan(poc_deque_t<poc_cdata_t*> &vars,poc_cdata_t * output,const poc_cdata_t & input);
extern void nearest(poc_deque_t<poc_cdata_t*> &vars,poc_cdata_t * output,const poc_cdata_t & input);
extern void transpose(poc_deque_t<poc_cdata_t*> &vars,poc_cdata_t * output,const poc_cdata_t & input);
extern void sum(poc_deque_t<poc_cdata_t*> &vars,poc_cdata_t * output,const poc_cdata_t & input);
extern void average(poc_deque_t<poc_cdata_t*> &vars,poc_cdata_t * output,const poc_cdata_t & input);
extern void dimmax(poc_deque_t<poc_cdata_t*> &vars,poc_cdata_t * output,const poc_cdata_t & input);
extern void gradx(poc_deque_t<poc_cdata_t*> &vars,poc_cdata_t * output,const poc_cdata_t & input);
extern void grady(poc_deque_t<poc_cdata_t*> &vars,poc_cdata_t * output,const poc_cdata_t & input);
extern void ellipseM(poc_deque_t<poc_cdata_t*> &vars,poc_cdata_t * output,const poc_cdata_t & input);
extern void sorted(poc_deque_t<poc_cdata_t*> &vars,poc_cdata_t * output,const poc_cdata_t & input);

extern void plg2grd(poc_deque_t<poc_cdata_t*> &vars,poc_cdata_t * output,const poc_cdata_t & input);
extern void resample(poc_deque_t<poc_cdata_t*> &vars,poc_cdata_t * output,const poc_cdata_t & input);
extern void orthogonalise(poc_deque_t<poc_cdata_t*> &vars,poc_cdata_t * output,const poc_cdata_t & input);
extern void crange(poc_deque_t<poc_cdata_t*> &vars,poc_cdata_t * output,const poc_cdata_t & input);
extern void fmerc(poc_deque_t<poc_cdata_t*> &vars,poc_cdata_t * output,const poc_cdata_t & input);
extern void imerc(poc_deque_t<poc_cdata_t*> &vars,poc_cdata_t * output,const poc_cdata_t & input);
extern void geokm(poc_deque_t<poc_cdata_t*> &vars,poc_cdata_t * output,const poc_cdata_t & input);
extern void fliph(poc_deque_t<poc_cdata_t*> &vars,poc_cdata_t * output,const poc_cdata_t & input);
extern void flipv(poc_deque_t<poc_cdata_t*> &vars,poc_cdata_t * output,const poc_cdata_t & input);


#if __cplusplus < 201103L
/* From http://en.wikipedia.org/wiki/Inverse_hyperbolic_function */
complex<double> asinh(const complex<double> &z){return log(z+sqrt(z*z+1.));}
#endif


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void print_poc_formula_help()

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// show help for all functions
/*----------------------------------------------------------------------------*/
{
  map<string,complex<double> (*)(const complex<double> &)>::iterator unaryI;
  map<string,void (*)(poc_deque_t<poc_cdata_t*> &vars,poc_cdata_t*,const poc_cdata_t&)>::iterator fieldI;
  string unaryList,fieldList;
  
  for(unaryI=unaryFunctions.begin();unaryI!=unaryFunctions.end();unaryI++)
    unaryList+=" "+unaryI->first;
  
  for(fieldI=fieldFunctions.begin();fieldI!=fieldFunctions.end();fieldI++)
    fieldList+=fieldI->first+"();";
  
  printf(HASH_LINE
    "# Calling an unknown function (with or without argument) prints THIS help.\n"
    "#\n"
    "# The most intuitive functions are:"+unaryList+"\n"
    "# Calling a most intuitive function without any argument ALSO prints THIS help.\n"
    "#\n"
    "# The not so intuitive functions are: "+fieldList+"\n"
    "# Calling a not so intuitive function without any argument prints ITS help.\n"
    "#\n"
    "# The most intuitive operators are: + - * /\n"
    "# The comparison operators are: < <= > >= == !=\n"
    "# %% is the scalar product operator.\n"
    "# ^ is the exponent operator.\n"
    "# && and || are mask operators.\n"
    "# << is the grid deformation operator.\n"
    "# = is the assignement operator.\n"
    "# ? is the eventual assignement operator.\n"
    "# [] is the index operator. If the index is not an interger, it interpolates!\n"
    "# : is the reverse index operator. It also interpolates, e.g. a=2&3&4;a[2.5:a] gives 2.5 .\n"
    "#\n"
    "# The following variables are constants: e nan pi d2r=pi/180 r2d=180/pi\n"
    "#\n"
    "# The column and row concatenation operators, respectively & and |, will, in the case of a dimension mismatch:\n"
    "#  - repeat scalars\n"
    "#  - or padd with nan.\n"
    HASH_LINE);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void poc_formula_implement_func(poc_deque_t<poc_cdata_t*> &vars,poc_cdata_t * output,const string & name,const poc_cdata_t & input)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/*------------------------------------------------------------------------------
  populate function maps */
  static bool functionMapsInitialised=false;
  
  if(!functionMapsInitialised){
    unaryFunctions["cos"]=cos;
    unaryFunctions["sin"]=sin;
    unaryFunctions["log"]=log;
    unaryFunctions["sinh"]=sinh;
    unaryFunctions["asinh"]=asinh;
    unaryFunctions["abs"]=abs_toComplex;
    unaryFunctions["arg"]=arg_toComplex;
    unaryFunctions["sign"]=sign_toComplex;
    unaryFunctions["real"]=real_toComplex;
    unaryFunctions["imag"]=imag_toComplex;
    unaryFunctions["conj"]=conj;
    unaryFunctions["floor"]=floor_toComplex;
    unaryFunctions["round"]=round_toComplex;
    unaryFunctions["ceil"]=ceil_toComplex;
    unaryFunctions["isnan"]=isnan_toComplex;
    unaryFunctions["nanIf0"]=nanIf0;
    
    fieldFunctions["isMasked"]=isMasked;
    fieldFunctions["zeroIfNan"]=zeroIfNan;
    fieldFunctions["nearest"]=nearest;
    fieldFunctions["T"]=transpose;
    fieldFunctions["sum"]=sum;
    fieldFunctions["average"]=average;

    fieldFunctions["open"]=open;
    fieldFunctions["close"]=close;
    fieldFunctions["write"]=write;
    
    fieldFunctions["print"]=print;
    fieldFunctions["printInfo"]=printInfo;
    fieldFunctions["renameDim"]=renameDim;
    fieldFunctions["dims"]=dims;
    fieldFunctions["dimmax"]=dimmax;

    fieldFunctions["gradx"]=gradx;
    fieldFunctions["grady"]=grady;
    fieldFunctions["ellipseM"]=ellipseM;
    fieldFunctions["sorted"]=sorted;
    
    /* map edition functions */
    fieldFunctions["plg2grd"]=plg2grd;
    fieldFunctions["resample"]=resample;
    fieldFunctions["ortho"]=orthogonalise;
    fieldFunctions["crange"]=crange;
    fieldFunctions["fmerc"]=fmerc;
    fieldFunctions["imerc"]=imerc;
    fieldFunctions["geokm"]=geokm;
    fieldFunctions["fliph"]=fliph;
    fieldFunctions["flipv"]=flipv;
    
    functionMapsInitialised=true;
    }
  
/*------------------------------------------------------------------------------
  search function maps */
  if(unaryFunctions.count(name)>0){
    if(output==0){
      print_poc_formula_help();
      return;
      }
    implement_unary(output,input,name,unaryFunctions[name]);
    }
  else if(fieldFunctions.count(name)>0){
    
    if(output!=0){
      string formula;
      formula=name+"("+getNameOrFormula(input.info)+")";
      output->info << poc_att_t("formula",formula);
      }
    
    fieldFunctions[name](vars,output,input);
    }
  else{
    STDERR_BASE_LINE("function name not recognised : "+name+"\n");
    print_poc_formula_help();
    }
  
  const poc_att_t *att;
  att=input.info.attributes.findP("");
  
  if(att!=0){
    string name=att->as_string();
    char *c=&name[name.size()-1];
    
    for(;*c>='a';(*c)--){
      int i;
      poc_cdata_t *var;
#if TRACK_COMMA_STACK
      STDERR_BASE_LINE_FUNC("deleting "+name+"\n");
#endif
      var=poc_formula_check_var(vars,name,&i);
      delete var;
      vars.erase(i);
      }
    
    }
  
}
