
/*******************************************************************************

  T-UGO tools, 2006-2019

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

<!-- USE FIND AND REPLACE TO EDIT THIS LINE, SO THAT print_help IS ALSO UPDATED -->
\brief Convert a SHOM grid file to NetCDF or gridit-cdf atlas files to SHOM's .a files.

<!-- A LINK TO main() or print_help() WILL NOT LINK TO THE RIGHT SOURCE ! -->
See the main <a href=#func-members>function</a> for how this works
and the print_help <a href=#func-members>function</a> for how to use this.
*/
/*----------------------------------------------------------------------------*/


#define MAIN_SOURCE

#include "version-macros.def" //for VERSION and REVISION

#include <stdio.h>
#include <fstream>

#include "tools-structures.h"

#include "functions.h" //fct_echo
#include "swap.h"
#include "poc-netcdf-data.hpp"


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_fgets(char *s, int n, FILE* f,int verbose=0)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  if(fgets(s,n,f)==0)
    return EOF;
  
  int count;
  
  count=max(0,(int)strlen(s)-1);
  s[count]='\0';
  
//   if(verbose) STDERR_BASE_LINE_FUNC("read \"%s\"\n",s);
  if(verbose) printf("%s, read: \"%s\"\n",__func__,s);
 
  return 0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  size_t shom_length(size_t ni,size_t nj,int verbose=0)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  size_t nij,n;
  
  nij=ni*nj;
  n=ceil(nij/4096.)*4096;
  
  if(verbose) STDOUT_BASE_LINE("%u*%u=%u<=%u (*4=%u)\n",nj,ni,nij,n,n*4);
  
  return n;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void swapBuffer(float *buffer,size_t n)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// MSB to LSB conversion
/**
\param[in,out] buffer array[n]
\param n
*/
/*----------------------------------------------------------------------------*/
{
  for(int i=0;i<n;i++){
    float *bi=&buffer[i];
    *bi=swap(*bi);
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int SHOM_data2nc(string inPath)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  inPath.erase(inPath.end()-1);
  
  const string
    ascPath=inPath+"b",
    binPath=inPath+"a",
    ncPath=inPath+"nc";
  STDOUT_BASE_LINE("ascii : "+ascPath+"\n");
  STDOUT_BASE_LINE("binary: "+binPath+"\n");
  STDOUT_BASE_LINE("output: "+ncPath+"\n");
  
  FILE *ascFile,*binFile;
  
  ascFile=fopen(ascPath.c_str(),"r");
  if(ascFile==0) TRAP_ERR_EXIT(errno,"Could not open "+ascPath+": %s\n",strerror(errno));
  
  binFile=fopen(binPath.c_str(),"r");
  if(binFile==0) TRAP_ERR_EXIT(errno,"Could not open "+binPath+": %s\n",strerror(errno));
  
  const size_t l=100;
  char *s, *name;
  
  s   =new char[l+1];
  name=new char[l];
  
  int status=0,mode=-9,count;
  size_t n,ni=0,nj=0;
  
/*------------------------------------------------------------------------------
  find mode and sizes */
  while(mode==-9 or ni==0 or nj==0){
    
    if(status) TRAP_ERR_EXIT(status,"Reached EOF in "+ascPath+"\n");
    status=poc_fgets(s,l,ascFile,1);
    
    count=0;
    
/*------------------------------------------------------------------------------
    field file */
    if(mode==-9) {
      count=sscanf(s,"%u x %u",&ni,&nj);
      }
    if(count==2) {
      mode=1;
      size_t pos=vpos('=',s,strlen(s));
      if(pos!=-1) {
        char *tmp=strdup(&s[pos+1]);
        pos=vpos(',',tmp,strlen(tmp));
        if(pos!=-1) tmp[pos]=' ';
        count=sscanf(tmp,"%u %u",&ni,&nj);
        }
      break;
      }
    
/*------------------------------------------------------------------------------
    grid file */
    count=sscanf(s,"%d '%s",&n,name);
    if(count!=2) continue;
    mode=2;
    switch(name[0]){
      case 'i':  ni=n; break;
      case 'j':  nj=n; break;
      }
    }
  
/*------------------------------------------------------------------------------
  compute record size */
  n=shom_length(ni,nj,1);
  STDOUT_BASE_LINE("mode=%d\n",mode);
  
  poc_data_t<float> d;
  double mi,ma;
  
#if 0
  printf("create variable\n");
  d.info << poc_dim_t("nj",nj);
  d.info << poc_dim_t("ni",ni);
  d.init("T");
  printf("reallocate : %u %d\n",d.data,n);
  deletep(&d.data);
  printf("reallocate : %u %d\n",d.data,n);
  d.data=new float[n];
#endif
  
  switch(mode){
  
/*------------------------------------------------------------------------------
  bathymetry */
  case 1:{
#if 0
    grid_t grid;
    grid.nx=ni;
    grid.ny=nj;
    
    /* NOTE: scanning strings because otherwise
      `E' (standing for East) is scanned as part of the number representation */
    char xmin[l],xmax[l],ymin[l],ymax[l];
    
    while(true){
      status=poc_fgets(s,l,ascFile);
      if(status!=0) TRAP_ERR_EXIT(-2,"format not recognised\n");
      
      count=sscanf(s,"lon: %[-+0-9.]E to %[-+0-9.]E. lat: %[-+0-9.]N to %[-+0-9.]N", xmin,xmax,ymin,ymax);
      if(count!=4)continue;
      
      sscanf(xmin,"%lg",&grid.xmin);
      sscanf(xmax,"%lg",&grid.xmax);
      sscanf(ymin,"%lg",&grid.ymin);
      sscanf(ymax,"%lg",&grid.ymax);
      grid.modeH=0;
      map_set_dxdy(&grid);
      STDERR_BASE_LINE("%d:",count);grid.brief_print(cerr);
      break;
      }
    
    status=map_completegridaxis(&grid,1);
    
    STDOUT_BASE_LINE("output: "+ncP+"\n");
    status=poc_save_grid(ncP,&d.info,grid,1);
    if(status!=0) NC_TRAP_ERROR(wexit,status,1,"poc_save_grid(\""+ncP+"\",,,...) error");
    status=poc_def_att(ncP,poc_att_t("production","constructed around " __LINE_FILE_PACKAGE_REVISION));
    if(status!=0) NC_CHKERR_BASE_LINE(status, "poc_def_att(\""+ncP+"\",poc_att_t(\"production\",...),) error");
#else
    const string loc="p";
    
    d.info << poc_dim_t("nj"+loc,nj);
    d.info << poc_dim_t("ni"+loc,ni);
    d.info << poc_att_t("coordinates",loc+"lat "+loc+"lon");
    d.info << poc_att_t("associates",loc+"lat "+loc+"lon");
    d.info << poc_att_t("axis","YX");
    
    status=poc_create(ncPath,"constructed around " __LINE_FILE_PACKAGE_REVISION);
    if(status!=0) NC_CHKERR_BASE_LINE(status,"poc_create(\""+ncPath+"\",\"...\") error");
#endif
    
    while(true){
      status=poc_fgets(s,l,ascFile);
      if(status!=0) TRAP_ERR_EXIT(-2,"format not recognised\n");
      
      count=sscanf(s," min,max %s = %lg %lg",name,&mi,&ma);
      if(count==3) break;
      }
    
    printf("%s[%g;%g]. ",name,mi,ma);fflush(stdout);
    d.info.init(name,NC_FLOAT,"depth_at_"+loc+"_location","m",mi);
    d.init();
    if(status!=0) NC_CHKERR_BASE_LINE(status, "poc_def_var(\""+ncPath+"\",) error");
    
    printf("Reading "+binPath+"... ");fflush(stdout);
    count=fread(d.data,sizeof(float),d.length,binFile);
    if(count!=d.length) TRAP_ERR_EXIT(errno,"fread() error: %s\n",strerror(errno));
    
    printf("swapping endianness... ");fflush(stdout);
    swapBuffer(d.data,d.length);
    
    printf(" writing to "+ncPath+"... ");fflush(stdout);
    status=d.write_data(ncPath,0,1);
    if(status!=0) NC_TRAP_ERROR(wexit,status,1,"poc_data_t<>::write_data(\""+ncPath+"\") error");
    
    printf(" done!\n");
    }break;
  
/*------------------------------------------------------------------------------
  fields */
  case 2:{
    status=poc_create(ncPath,"constructed around " __LINE_FILE_PACKAGE_REVISION);
    if(status!=0) NC_CHKERR_BASE_LINE(status,"poc_create(\""+ncPath+"\",\"...\") error");
    
    printf("create variable\n");
    d.info << poc_dim_t("nj",nj);
    d.info << poc_dim_t("ni",ni);
    d.init("T");
    printf("reallocate : %u %d\n",d.data,n);
    deletep(&d.data);
    printf("reallocate : %u %d\n",d.data,n);
    d.data=new float[n];
    
    printf("start scanning binary file\n");
    while(true){
      status=poc_fgets(s,l,ascFile);
      if(status!=0) break;
      
      count=sscanf(s,"%[^:]:  min,max = %lg %lg",name,&mi,&ma);
      if(count!=3) continue;
      
      printf("%s[%g;%g]. ",name,mi,ma);fflush(stdout);
      
      const string loc(name,1);
      
      if(strcasecmp(&name[1],"lon")==0){
        d.info.init(name,NC_FLOAT,"longitude_at_"+loc+"_location","degrees");
        }
      else if(strcasecmp(&name[1],"lat")==0){
        d.info.init(name,NC_FLOAT,"latitude_at_"+loc+"_location","degrees");
        }
      else{
        d.info.init(name,NC_FLOAT);
        }
      
      d.info.dimensions.clear();
      d.info << poc_dim_t("nj"+loc,nj);
      d.info << poc_dim_t("ni"+loc,ni);
      d.info << poc_att_t("coordinates",loc+"lat "+loc+"lon");
      d.info << poc_att_t("associates",loc+"lat "+loc+"lon");
      d.info << poc_att_t("axis","YX");
      
      printf("Reading "+binPath+"... ");fflush(stdout);
      count=fread(d.data,sizeof(float),n,binFile);
      if(count!=n) TRAP_ERR_EXIT(errno,"fread() error: %s\n",strerror(errno));
      
      printf("swapping... ");fflush(stdout);
      swapBuffer(d.data,d.length);
      
      printf(" writing to "+ncPath+"... ");fflush(stdout);
      status=d.write_data(ncPath,0,1);
      if(status!=0) NC_TRAP_ERROR(wexit,status,1,"poc_data_t<>::write_data(\""+ncPath+"\") error");
      
      printf(" done.\n");
      }
    }break;
  
  default:
    TRAP_ERR_EXIT(ENOEXEC,"programming error\n");
    }
  
  fclose(ascFile);
  fclose(binFile);
  
  return 0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void SHOM_nc2a(char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int j;
  
  for(j=0;argv[j]!=0;j++){
    const string
      ncPath=argv[j];
    if(strrncasecmp(ncPath,".nc")!=0) continue;
    
    int i,n,d,status;
    
    n=ncPath.length();
    string binPath(ncPath,0,n-2),ascPath=binPath+"b";
    binPath+="a";
    FILE *binFile,*ascFile;
    
    binFile=fopen(binPath.c_str(),"w");
    if(binFile==0) TRAP_ERR_EXIT(errno,"Could not open "+binPath+" (%d %s)\n",errno,strerror(errno));
    
    ascFile=fopen(ascPath.c_str(),"w");
    if(ascFile==0) STDERR_BASE_LINE("Could not open "+ascPath+" (%d %s)\n",errno,strerror(errno));
    
    poc_global_t glob;
    status=poc_inq(ncPath,&glob);
    
    float *buffer=0;
    const poc_var_t *var;
    
    /* find a variable on a T point */
    printf("Scanning "+ncPath+" for a variable \"_at_t_location\" : ");fflush(stdout);
    int verbose=1;
    i=findAnyVariableThasIs(glob,2,"_at_t_location",verbose);
    if(i<0) TRAP_ERR_EXIT(NC_ENOTVAR,"no variable \"_at_t_location\" found.\n");
    
    var=&glob.variables[i];
    printf(""+var->name+" .\n");fflush(stdout);
    
    const poc_list_t<poc_dim_t> &dims=var->dimensions;
    
    const size_t tsize=dims[0].len*dims[1].len;
    
    n=shom_length(dims[0].len,dims[1].len);
    buffer=new float[n];
    
    printf("Converting "+ncPath+" to "+binPath);
    if(ascFile!=0){
      printf(" and "+ascPath);
      status=fprintf(ascFile,"%5u    'idm   ' = longitudinal array size\n",dims[0].len);
      status=fprintf(ascFile,"%5u    'jdm   ' = latitudinal array size\n" ,dims[1].len);
      }
    printf(" :\n");
    
    for(i=0;i<glob.variables.size();i++){
      var=&glob.variables[i];
      
      if(var->dimensions.size()!=2)
        continue;
      
      for(d=0;d<2;d++)
        if(var->dimensions[d]!=dims[d])
          break;
      
      if(d<2)
        continue;
      
      printf(""+var->name);fflush(stdout);
      
      status=poc_get_vara(ncPath,*var,0,buffer);
      if(status!=0) NC_TRAP_ERROR(wexit,status,1,"poc_get_var(\""+ncPath+"\",(\""+var->name+"\"),) error");
      
      if(ascFile!=0){
        float mask;
        status=poc_decode_mask(*var,0,0,&mask);
        range_t<float> bR;
        bR.init(buffer,mask,tsize);
        status=fprintf(ascFile,""+var->name+":  min,max = %12g %12g\n",bR.min,bR.max);
        }
      
      swapBuffer(buffer,n);
      
      status=fwrite(buffer,sizeof(float),n,binFile);
      if(status<n) TRAP_ERR_EXIT(errno,"fwrite(,,,("+binPath+")) returned %d<%d: %s\n",status,n,strerror(errno));
      
      printf(", ");fflush(stdout);
      }
    
    delete[]buffer;
    
    if(ascFile!=0)
      status=fclose(ascFile);
    
    status=fclose(binFile);
    
    printf("done.\n");
    }
  
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
    "  %s [OPTIONS] file\n",prog_name);
  printf(
    "  %s file1.nc [ file2.nc ... ]\n",prog_name);
  printf("\n"
    "DESCRIPTION\n"
/* USE FIND AND REPLACE TO EDIT THE LINE BELOW, SO THAT THE FILE HEAD IS ALSO UPDATED */
    "  Convert a SHOM grid file to NetCDF or gridit-cdf NetCDF atlas files to SHOM's .a files.\n"
    "The file(s) is(are) either:\n"
    "  - one SHOM ascii .grd file\n"
    "  - one of the SHOM pair of binary .a and ASCII .b files\n"
    "  - one or several gridit-cdf NetCDF file(s).\n"
    "\n"
    "OPTIONS\n"
    "  -h,--help  Show this help and exit.\n"
    "  -3  create NetCDF 3 file, that should be huge with the dummy variables, instead of a ~20x compressed NetCDF 4 file\n"
    "  --small3  create NetCDF 3 file without the dummy variables\n"
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
  int i,n,status;//general index, argument index, status
  
  vector<poc_data_t<double>* > vars;
  
  string cmd=fct_echo(argc,argv);

  const char *inPath=NULL,*keyword=NULL;
  int compressed=NC_NETCDF4|NC_CLASSIC_MODEL,small=0;
  
  if(argc<2 ||
    strcmp(argv[1],"--help")==0 ||
    strcmp(argv[1],"-h")==0){
    print_help(argv[0]);
    wexit(0);
    }
  
  n=1;
  while (n < argc) {
    keyword=argv[n];
    
    if(!strcmp(keyword,"--small3")){
      small=1;
      compressed=0;
      n++;
      continue;
      }
    
    switch (keyword[0]) {
      case '-':
        switch (keyword[1]) {
        case '3' :
          compressed=0;
          n++;
          break;
        
        default:
          STDERR_BASE_LINE("unknown option %s\n",keyword);
          print_help(argv[0]);
          wexit(-1);
          }
        
        break;
        
      default:
        if(n==1 and strrncasecmp(keyword,".nc")==0){
          SHOM_nc2a(&argv[n]);
          TRAP_ERR_EXIT(0,"exiting\n");
          }
        if(inPath){
          if(strcmp(inPath,keyword)){
            STDERR_BASE_LINE("ONLY ONE INPUT FILE ACCEPTED: CHOOSE BETWEEN %s AND %s\n",inPath,keyword);
            print_help(argv[0]);
            wexit(-1);
            }
          
          STDERR_BASE_LINE("WARNING : same input file given twice\n");
          }
        else
          inPath=keyword;
        
        n++;
        break;
      }
    }
  
  if(strrncasecmp(inPath,".a")==0 or strrncasecmp(inPath,".b")==0){
    status=SHOM_data2nc(inPath);
    return status;
    }
  
  STDOUT_BASE_LINE("reading %s\n",inPath);
  
  ifstream inFile(inPath);
  string line,word;
  
  poc_data_t<double> *data;
  poc_global_t global;
  double value;
  
/*------------------------------------------------------------------------------
  LINE */
  while(inFile.good()){
    getline(inFile,line);
    if(strncmp("LINE ",line)) break;
    STDOUT_BASE_LINE("got "+line+"\n");
    istringstream lineStream(line);
    lineStream>>word;
    
    poc_dim_t dim;
    lineStream>>dim.name;
    lineStream>>dim.len;
    global.dimensions << dim;
    
    lineStream>>word;/* WE or NS */
    
    string axis,location;
    lineStream>>word;
    if(word=="LONGITUDE")
      axis="X";
    else if(word=="LATITUDE")
      axis="Y";
    else
      TRAP_ERR_EXIT(ENOEXEC,"not coded yet for "+word+" axis type\n");
    location=dim.name.substr(dim.name.length()-1);
    
    data=new poc_data_t<double>();
    vars.push_back(data);
    data->info.init(dim.name,NC_DOUBLE,word+"_at_"+location+"_location");
    data->info << dim;
    data->info << poc_att_t("axis",axis);
    poc_print(data->info);
    
    lineStream>>word;
    if(word[word.length()-1]==':') word.resize(word.length()-1);
    if(word=="GIVEN_BELOW") continue;
    if(word!="START,DELTA") TRAP_ERR_EXIT(ENOEXEC,"data format "+word+" not recognised\n");
    
    data->init("");
    
    double increment;
    lineStream>>value;
    if(axis=="X")
      value=degree_recale(value,0.);
    lineStream>>increment;
    STDOUT_BASE_LINE("from %g ",value);
    
    for(i=0;i<dim.len;i++){
      data->data[i]=value;
      value+=increment;
      }
    
    printf("to %g step %g\n",value,increment);
    }
  
/*------------------------------------------------------------------------------
  GRID */
  while(!strncmp("GRID ",line)){
    if(small)
      STDOUT_BASE_LINE("skipping "+line+"\n");
    else{
      STDOUT_BASE_LINE("got "+line+"\n");
      istringstream lineStream(line);
      lineStream>>word;
      
      poc_var_t dummy;
      lineStream>>word;
      dummy.init(word,NC_BYTE,"dummy");
      
      while(lineStream.good()){
        lineStream>>word;
        
        lineStream>>word;
        if(word=="UNKNOWN") break;
        
        poc_dim_t *dim=global.dimensions.findP(word);
        if(!dim) NC_TRAP_ERROR(wexit,NC_EBADDIM,1,word+" dimension not found");
        
        dummy << *dim;
        }
      
      global.variables << dummy;
      poc_print(dummy);
      }
    
    getline(inFile,line);
    }
  
/*------------------------------------------------------------------------------
  COORDS */
  while(inFile.good()){
    STDOUT_BASE_LINE("got "+line+"\n");
    istringstream lineStream(line);
    
    lineStream>>word;
    if(word!="COORDS") TRAP_ERR_EXIT(ENOEXEC,"data format "+word+" not recognised\n");
    
    lineStream>>word;
    for(i=0;i<vars.size() && vars[i]->info.name!=word;i++);
    
    data=vars[i];
    data->init("");
    STDOUT_BASE_LINE("reading %d values for "+data->info.name+"\n",data->length);
    
    for(i=0;i<data->length;i++){
      inFile>>data->data[i];
      inFile>>value;
      }
    
    do{
      getline(inFile,line);
      }while(line=="" && inFile.good());
    }
  
/*------------------------------------------------------------------------------
  output */
  string oPath=inPath;
  oPath.replace(oPath.find(".grd"),string::npos,".nc");
  STDOUT_BASE_LINE("writing ");
  printf(compressed?"compressed NetCDF4":"NetCDF3");
  printf(" "+oPath+": ");fflush(stdout);
  
  poc_create(oPath,global,0,compressed);
  status=poc_def_att(oPath,poc_att_t("production","constructed around " __LINE_FILE_PACKAGE_REVISION));
  status=poc_def_att(oPath,poc_att_t("history",cmd));
  
  for(i=0;i<vars.size();i++){
    data=vars[i];
    if(i) printf(", ");
    printf(data->info.name);fflush(stdout);
    vars[i]->write_data(oPath);
    }
  
  if(vars.size()) printf(", ");
  printf("done.\n");
  
  return 0;
}
