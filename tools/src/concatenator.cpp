
/**************************************************************************

  T-UGO tools, 2006-2011

  Unstructured Ocean Grid initiative

***************************************************************************/
/** \file

\author  Florent Lyard      LEGOS/CNRS, Toulouse, France. florent.lyard@legos.obs-mip.fr
\author  Yoann Le Bars      LEGOS, Toulouse, France (PhD)
\author  Laurent Roblou     LEGOS/CNRS, Toulouse, France
\author  Damien Allain      LEGOS/CNRS, Toulouse, France
\author  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
\author  Frédéric Dupont    Université de Laval à Québec, Canada

<!-- USE FIND AND REPLACE TO EDIT THIS LINE, SO THAT print_help IS ALSO UPDATED -->
\brief Concatenate model grid and output files.

<!-- A LINK TO main() or print_help() WILL NOT LINK TO THE RIGHT SOURCE ! -->
See the main function for how this works
and the print_help function for how to use this.
*/
/*----------------------------------------------------------------------------*/

#define MAIN_SOURCE

#include "version-macros.def" //for VERSION and REVISION

#include <stdio.h>
#include <unistd.h>

#include "poc-netcdf-data.hpp"


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<class V>  poc_dim_t * findConcatenationDim_template(V & info,const string & name)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// find time dimension
/**
\param *info variable
\param *index set to the index of the dimension
\returns a pointer to the dimension or NULL if none found
*/
/*----------------------------------------------------------------------------*/
{
  poc_dim_t *dim;
  
  if(name==""){
    dim=findTimeDim(info);
    }
  else{
    dim=info.dimensions.findP(name);
    }
  
  return dim;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  poc_dim_t * findConcatenationDim(poc_var_t & info,const string & name)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  poc_dim_t *dim;
  
  dim=findConcatenationDim_template(info,name);
  
  return dim;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  poc_dim_t * findConcatenationDim(poc_global_t & glob,const string & name)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  poc_dim_t *dim;
  
  dim=findConcatenationDim_template(glob,name);
  
  return dim;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void print_help(const char *prog_name)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/** \brief prints help of this programme.

\param prog_name name of the program to be printed by this help function : argv[0]
*/
/*----------------------------------------------------------------------------*/
{
/** The code of the body of this function is :
    \code /**/ // COMPILED CODE BELOW !!!
  printf("\n"
    "NAME AND VERSION\n"
    "  %s version " VERSION " " REVISION "\n", prog_name);
  printf("\n"
    "USE\n"
    "  %s [OPTION] input1.nc [ input2.nc ... ] output\n",prog_name);
  printf("\n"
    "DESCRIPTION\n"
/* USE FIND AND REPLACE TO EDIT THE LINE BELOW, SO THAT THE FILE HEAD IS ALSO UPDATED */
    "  Concatenate model grid and output files.\n"
    "  Will exit if the output file already exists or if one or more of the input can not be read.\n"
    "  With less than 2 arguments, print this help.\n"
    "\n"
    "OPTION\n"
    "  -a : followed by dimension. Default is whatever is recognised as a time dimension.\n"
    "     See also option -t.\n"
    "  -v : followed by comma-separated list of names of variables to concatenate. Default is everything.\n"
    "  -t : followed by comma-separated list of names of variables to which a time dimension should be added.\n"
    "     The name of the dimension is then what is specified with option -a. Default : time.\n"
    "     BUG : if used, the count of time frames will be the count of input files!\n"
    ); /** \endcode */
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// function that will be ran when the executable is started
/** See the print_help <a href=#func-members>function</a> for its use.
*/
/*----------------------------------------------------------------------------*/
{
  const int lastArgI=argc-1;
  int firstArgI=1;
  string dimName;
  char *outPath=argv[lastArgI];
  const char *keyword;
  poc_list_t<poc_name_id_t> varNameList,timedList;
  
  /* file, variable and attribute indexes */
  int i,j,k,status;
  
  /* command line */
  string cmd;
  
  bool do_help=false;
  
  cmd=fct_echo( argc, argv);
  
/*------------------------------------------------------------------------------
  check arguments */
  
  if(argc<=2){
    print_help(argv[0]);
    return 0;
    }
  
  if(access(outPath,F_OK)==0){
    fprintf(stderr,"*** %s already exists ***\n",outPath);
    do_help=true;
    }
  
  firstArgI=1;
  while (firstArgI < argc) {
    keyword=argv[firstArgI];
    
    if(strcmp(keyword,"-a")==0){
      dimName=argv[firstArgI+1];
      firstArgI+=2;
      continue;
      }
    
    if(strcmp(keyword,"-v")==0){
      istringstream parsedList(argv[firstArgI+1]);
      string token;
      
      while(parsedList.good()){
        getline(parsedList,token,',');
        
        if(token=="")
          break;
        
        varNameList<<poc_name_id_t(token);
        }
      
      firstArgI+=2;
      continue;
      }
    
    if(strcmp(keyword,"-t")==0){
      istringstream parsedList(argv[firstArgI+1]);
      string token;
      
      while(parsedList.good()){
        getline(parsedList,token,',');
        
        if(token=="")
          break;
        
        timedList<<poc_name_id_t(token);
        }
      
      firstArgI+=2;
      continue;
      }
    
    break;
    }
  
  for(i=firstArgI;i<lastArgI;i++){
    char *inPath=argv[i];
    if(access(inPath,R_OK)==0) continue;
    fprintf(stderr,"*** can not read %s : %s ***\n",inPath,strerror(errno));
    do_help=true;
    }
  
  if(do_help){
    print_help(argv[0]);
    return -1;
    }
  
  if(timedList.size()>0 and dimName==""){
    dimName="time";
    }
  
/*------------------------------------------------------------------------------
  header */
  
  poc_global_t gin,gout("constructed around " __LINE_FILE_PACKAGE_REVISION);
  bool withTimeDim;
  int maxLen=-1;
  
  /* frame indexes */
  int ofi=0,ifi;
  
  for(i=firstArgI;i<lastArgI;i++){
    char *inPath=argv[i];
    status=poc_inq(inPath,&gin);
    NC_CHKERR_BASE_LINE(status,"poc_inq(\"%s\",) error",inPath);
    
    const size_t nVar=gin.variables.size();
    
    for(j=0;j<nVar;j++){
      poc_var_t *var=&gin.variables[j];
      
      if(varNameList.size()>0 and varNameList.findP(var->name)==0)
        continue;
      
      size_t len;
      status=poc_get_var_length(*var,&len);
      updatemax(&maxLen,(int)len);
      
      if(timedList.size()>0){
        ifi=1;
        }
      else if(dimName==""){
        withTimeDim=unlimitTimeDim(var);
        }
      else{
        poc_dim_t *vdim;
        vdim=var->dimensions.findP(dimName);
        
        if(vdim==0){
          withTimeDim=false;
          }
        else{
          withTimeDim=true;
          ifi=vdim->len;
          }
        
        }
      gout.variables<<*var;
      }
    
    ofi+=ifi;
    
    for(j=0;j<gin.attributes.size();j++){
      poc_att_t *att=&gin.attributes[j];
      gout.attributes<<*att;
      }
    
    }
  
  if(timedList.size()>0){
    poc_dim_t timeDim(dimName,ofi,true);
    for(j=0;j<gout.variables.size();j++){
      poc_var_t *var=&gout.variables[j];
      if(timedList.findP(var->name)==0)
        continue;
      var->dimensions.insert(0,timeDim);
      }
    }
  else if(dimName!=""){
    for(j=0;j<gout.variables.size();j++){
      poc_var_t *var=&gout.variables[j];
      poc_dim_t *vdim;
      vdim=var->dimensions.findP(dimName);
      if(vdim!=0)
        vdim->init(ofi);
      }
    }
  
  gout<<poc_att_t("history",cmd);
  printf("creating %s ...",outPath);fflush(stdout);
  status=poc_create(outPath,gout,1);
  if(status!=NC_NOERR) NC_TRAP_ERROR(return,status,1,"poc_create(\"%s\",,) error");
  printf(" done (max length=%d).\n",maxLen);
  
/*------------------------------------------------------------------------------
  data */
  
  ofi=0;
  double *data=new double[maxLen];
  int inNc,outNc;
  
  status=nc_open(outPath,NC_WRITE,&outNc);
  
  for(i=firstArgI;i<lastArgI;i++){
    char *inPath=argv[i];
    
    printf("processing %s",inPath);fflush(stdout);
    
    status=nc_open(inPath,NC_NOWRITE,&inNc);
    if(status!=NC_NOERR) NC_TRAP_ERROR(wexit,status,1,"nc_open(\"%s\",NC_NOWRITE,) error",inPath);
    status=poc_inq(inNc,&gin);
    if(status!=NC_NOERR) NC_TRAP_ERROR(wexit,status,1,"poc_inq((\"%s\"),) error",inPath);
    
    /* number of input frames */
    int nif;
    poc_dim_t *gdim,*vdim;
    bool incrementFrames;
    
    gdim=findConcatenationDim(gin,dimName);
    incrementFrames=(gdim!=0 or timedList.size()>0);
    
    if(incrementFrames){
      if(gdim!=0)
        nif=gdim->len;
      else
        nif=1;
      printf(" (%d+%d frames)",ofi,nif);fflush(stdout);
      }
    else
      nif=1;
    
    for(j=0;j<gin.variables.size();j++){
      poc_var_t *var=&gin.variables[j];
      poc_var_t *outVar=gout.variables.findP(var->name);
      int run;
      
      if(outVar==0)
        continue;
      
      printf("%c "+var->name,j==0?':':',');fflush(stdout);
      
      if(dimName!=""){
        for(run=0;run<2;run++){
          poc_var_t *thisVar=(run==0)?var:outVar;
          vdim=findConcatenationDim(*thisVar,dimName);
          
          for(k=0;k<thisVar->dimensions.size();k++){
            poc_dim_t *dim=&thisVar->dimensions[k];
            dim->isunlimited=(dim==vdim);
            /* make sure vdim is now the only time dim */
            if(not dim->isunlimited)
              dim->name="";
            }
          }
        }
      
      for(ifi=0;ifi<max(1,nif);ifi++){
        status=poc_get_vara(inNc ,*var   ,ifi    ,data);
        status=poc_put_vara(outNc,*outVar,ifi+ofi,data);
        if(vdim==0)break;
        }
      
      }
    
    status=nc_close(inNc);
    //status=nc_sync(outNc);
    
    if(incrementFrames)
      ofi+=nif;
    
    printf(", done.\n");
    }
  
  status=nc_close(outNc);
  
  return 0;
}
