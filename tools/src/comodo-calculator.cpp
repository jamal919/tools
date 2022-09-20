
/*******************************************************************************

  T-UGO tools, 2006-2011

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

<!-- USE FIND AND REPLACE TO EDIT THIS LINE, SO THAT print_help IS ALSO UPDATED -->
\brief Calculator.

<!-- A LINK TO main() or print_help() WILL NOT LINK TO THE RIGHT SOURCE ! -->
See the main <a href=#func-members>function</a> for how this works
and the print_help <a href=#func-members>function</a> for how to use this.
*/
/*----------------------------------------------------------------------------*/


#define MAIN_SOURCE

#include "version-macros.def" //for VERSION and REVISION

#include <stdio.h>

#include "functions.h" //fct_echo
#include "poc-data-operators.hpp"
#include "map.h"


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void print_help(const char *prog_name)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/** \brief prints help of this programme.

\param prog_name name of the program to be printed by this help function : as[0]
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
    "  %s -i file1.nc [var1]=ncVarName1 [ncPhaseVarName1] [ [var2]=ncVarName2 [ncPhaseVarName2] ...] [-i file2.nc ...] [OPTIONS]\n",prog_name);
  printf("\n"
    "DESCRIPTION\n"
/* USE FIND AND REPLACE TO EDIT THE LINE BELOW, SO THAT THE FILE HEAD IS ALSO UPDATED */
    "  Calculator.\n"
    );
  printf("\n"
    "OPTIONS :\n"
    "  -h,--help  Show this help and exit.\n"
    "  -i  followed by path to netcdf input file, then variables\n"
    "  -o  followed by path to netcdf output file\n"
    "  -f  followed by formula. If it starts with `/', `./' or `../', it is the file path of the formula\n"
    "  -d  do not standardise dimensions\n"
    );
  print_OPENMP_help(prog_name); /** \endcode */
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int narg, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/** \brief main function

Document this to show how this works.
See the print_help <a href=#func-members>function</a> for how to use this.
*/
/*----------------------------------------------------------------------------*/
{
  int i,status=0;//general index, status

  char *oP=NULL;//output path
  const char *a,*iP=NULL;//input path
  char *formula=NULL;
  
  grid_t grid;
  poc_var_t gridVar;
  bool keepDims=false;
  poc_deque_t<poc_cdata_t*> vars;

  struct timeval before;
  string cmd=fct_echo(narg,argv);

  if(narg<2){
    print_help(argv[0]);
    exit(0);
    }
  
  for(i=1;i<narg;i++){
    a=argv[i];
    if(a[0]=='-'){
      if(strcmp(a,"--help")==0 ||
         strcmp(a,"-h")==0){
        print_help(argv[0]);
        exit(0);
        }
      //else
      if(a[2]!='\0'){
        printf("unknown option %s\n",a);
        print_help(argv[0]);
        exit(-1);
        }
      //else
      switch(a[1]){
      
/*------------------------------------------------------------------------------
      input */
      case 'i':
        i++;
        iP=argv[i];
        printf("From %s, reading:\n",iP);
        
        for(i++;i<narg;i++){
          const string va=argv[i];
          if(va[0]=='-') break;
          const int pos=va.find_first_of('=');
          if(pos<0) break;
          
          const string ncVarName=va.substr(pos+1);
          string ncPhaseVarName;
          if(i+1<narg){
            const char *pa=argv[i+1];
            if(pa[0]!='-' && strchr(pa,'=')==0){
              printf("phase: %s\n",pa);
              ncPhaseVarName=pa;
              i++;
              }
            }
          
          poc_cdata_t *varp=new poc_cdata_t;
          printf("("+ncVarName+","+ncPhaseVarName+")\n");
//           status=varp->init(iP,ncVarName,ncPhaseVarName);
          status=varp->init(iP,ncVarName,ncPhaseVarName,"");
          NC_TRAP_ERROR(wexit,status,1,"poc_cdata_t::init(\"%s\",\""+ncVarName+"\",\""+ncPhaseVarName+"\",\"\") error",iP);
          poc_print(varp->info,cout,0);
          printf("mask=%g\n",real(varp->mask));
          fflush(stdout);
          
          if(!vars.size()){
            /* get grid from 1st variable */
            status=poc_get_grid(iP,varp->info,&grid);
            if(status!=0)
              grid.print(cout,__FILE__ ":" TOSTRING(__LINE__) ":" "  ");
            }
          
          if(pos>0)
            varp->info.name=va.substr(0,pos);
          else
            /* default to NetCDF variable name */
            varp->info.name=ncVarName;
          
          printf(" -> "+varp->info.name+"\n");
          
          vars.push_back(varp);
          }
        
        i--;
        break;
      
/*------------------------------------------------------------------------------
      output */
      case 'o':
        i++;
        oP=strdup(argv[i]);
        break;
      
/*------------------------------------------------------------------------------
      formula */
      case 'f':
        i++;
        formula=argv[i];
        break;
      
/*------------------------------------------------------------------------------
      formula */
      case 'd':
        keepDims=true;
        break;
      
/*----------------------------------------------------------------------------*/
      default:
        printf("unknown option %s\n",a);
        print_help(argv[0]);
        exit(-1);
        }
      continue;
      }
    
    printf("unknown option %s\n",a);
    print_help(argv[0]);
    exit(-1);
    }
  
  if(vars.size()==0 || oP==0){
    if(vars.size()==0)
      printf("*** input variable missing ***\n");
    if(oP==0)
      printf("*** output file missing ***\n");
    }
  
  if(formula==0){
    printf("*** formula missing ***\n");
    print_help(argv[0]);
    exit(-1);
    }
  
  if( strncmp(  "/",formula)==0 or
      strncmp( "./",formula)==0 or
      strncmp("../",formula)==0 ){
    
    const char *fpath=formula;
    size_t fsize;
    
    status=get_file_size(fpath,&fsize);
    formula=new char[fsize+1];
    formula[fsize]='\0';
    
    printf(HASH_LINE "Reading formula from %s :\n",fpath);
    
    FILE *f;
    f=fopen(fpath,"r");
    
    status=fread(formula,1,fsize,f);
    printf("%s\n" HASH_LINE,formula);
    
    fclose(f);
    }
  
  gettimeofday(&before);
  
/*---------------------------------------------------------------------*//**<h1>
  grid output </h1>*/
  if(oP!=0){
    printf("Preparing %s:",oP);fflush(stdout);
    
    int verbose=0;
    if(grid.Hsize()>0){
      printf(" grid,");fflush(stdout);
      status=poc_save_grid(oP,&gridVar,grid,1,verbose);
      }
    else{
      status=poc_create(oP,poc_global_t(),verbose);
      }
    
    printf(" attributes,");fflush(stdout);
    status=poc_def_att(oP,poc_att_t("production","constructed around " __LINE_FILE_PACKAGE_REVISION));
    status=poc_def_att(oP,poc_att_t("history",cmd));
    
    printf(" done in %gs.\n",difftime(&before));
    }
  
/*---------------------------------------------------------------------*//**<h1>
  calculation </h1>*/
  printf("Calculations...\n");fflush(stdout);
  /*
  h.info=gridVar;
  h.info.init("h",NC_DOUBLE,"depth","m");
  
  g.init(h.info);
  g.info.init("g",NC_DOUBLE,"gravity_at_mean_sea_level","ms-2");
  map_completegridaxis_2(&grid);
  for(i=0;i<g.length;i++){
    g.data[i]=gravitational_constant_01deg(grid.y[i]);
    }
  
  vars.push_back(&h);
  vars.push_back(&g);
  */
  
  poc_formula_parse(vars,formula);
  
  printf("done in %gs.\n",difftime(&before));fflush(stdout);
  
/*---------------------------------------------------------------------*//**<h1>
  grid output </h1>*/
  if(oP==0) return status;
  
  printf("Saving to %s:\n",oP);fflush(stdout);
  for(i=0;i<vars.size();i++){
    poc_cdata_t *varp=vars[i];
    int length=varp->length;
    
    const string varname=varp->info.name;
    cout << varp->info.name << "(";
    poc_print(",",varp->info.dimensions,cout);
    cout << "):" << length << " values\n";
    
    if(length==grid.Hsize()){
      poc_att_t att;
      varp->info.attributes.find("formula",&att);
      
      if(keepDims){
        varp->info.attributes=gridVar.attributes;
        }
      else{
        varp->info=gridVar;
        varp->info.init(varname,NC_DOUBLE);
        }
      
      if(att.name!="") varp->info << att;
      printf("mask=%g\n",varp->mask);
      if(varp->mask!=NC_FILL_COMPLEX)
        varp->info << poc_att_t(_FillValue,real(varp->mask));
      }
    
    //poc_print(varp->info);
    fflush(stdout);
    }
  
  status=poc_save_vars(oP,vars,0,0,0);
  
  printf("done in %gs.\n",difftime(&before));
  
  return status;
}
