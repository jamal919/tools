
/*******************************************************************************

  T-UGO tools, 2006-2014

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

\brief functions definitions for COMODO (COmmunauté de MODélisation Océanique) grid parsing
*/
/*----------------------------------------------------------------------------*/

#if TUGO

#include "tugo-prototypes.h"

#include "poc-netcdf-data.hpp"
int poc_formula_parse(poc_deque_t<poc_cdata_t*> &vars,const char *str){
  TRAP_ERR_EXIT(ENOEXEC,"%s not defined in TUGO\n",__func__); }
poc_cdata_t *poc_formula_check_var(const poc_deque_t<poc_cdata_t*> &vars,const string & name, int *index=0){
  TRAP_ERR_EXIT(ENOEXEC,"%s not defined in TUGO\n",__func__); }
#else

#include "version-macros.def" //for VERSION and REVISION

#include "fe-proto.h"

#include "poc-data-operators.hpp"

#endif

#include "poc-grib.h"

#include "grd.h" //for grd_loadgrid
#include "statistic.h" //for poc_minmax


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_parse_coordinates_attribute(const poc_var_t &var,const poc_global_t & global,const poc_var_t *vs[],int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// parse coordinates (or associate) attribute
/** mainly for, e.g., C-grid files with
  P, U, V and F variables that all have the same stupid dimensions.
As the badness of these files should also extend
to the absence of standard_name attribute in the coordinate variable,
the template for the command to salvage them is :\code
ncatted -a standard_name,<badlonvarname>,o,char,longitude[_at_X_location] -a ... -a coordinates,var,o,char,"badlonvarname ..." -a ... [-O] file [optional_new_file]
\endcode

\returns NC_NOERR if success or the NC_ENOTATT error code if error.
*/
/*----------------------------------------------------------------------------*/
{
  const poc_att_t *att=var.attributes.findP("coordinates");
  
  if(att==0) {
    if(verbose>1) STDERR_BASE_LINE("trying associate\n");
    att=var.attributes.findP("associate");
    }
  
  if(att==0){
    NC_TRAP_ERROR(return,NC_ENOTATT,verbose,"No coordinates attribute for "+var.name);
    }
  
  const string coordList=att->as_string();
  if(verbose>1) STDERR_BASE_LINE("using "+var.name+":"+att->name+" = \""+coordList+"\"\n");
  istringstream coordStream(coordList);
  
/* NOTE: the latest CF Conventions documents DO NOT give ANY indication
about the order of coordinate variable names in the coordinates attribute,
which makes sense because of time series with time-dependent coordinates */
  string coordName;
  do{
    coordStream>>coordName;
    bool isLon,isLat,isDepth;
    int check;
    
    if(verbose>1) STDERR_BASE_LINE(" eof=%d good=%d "+coordName+"\n",coordStream.eof(),coordStream.good());
    
    /* time variable is not in grid files */
    if(isT(coordName.c_str()))
      /* so make sure you do not check for its existence */
      continue;
    
    const poc_var_t *coordVar=global.variables.findP(coordName);
    if(coordVar==0) NC_TRAP_ERROR(return,NC_ENOTVAR,verbose,"error with "+coordName+" in "+var.name+":"+att->name+" attribute");
    
    check=0;
    
    if(check==0) {
      isLon=nameS_start_or_end_with(*coordVar,"longitude");
      isLat=nameS_start_or_end_with(*coordVar,"latitude");
      isDepth=nameS_start_or_end_with(*coordVar,"depth");
      check=(int)isLon+isLat+isDepth;
      }
    
    if(check==0) {
      isLon=nameS_start_or_end_with(*coordVar,"COLUMNS");
      isLat=nameS_start_or_end_with(*coordVar,"LINES");
      isDepth=nameS_start_or_end_with(*coordVar,"LAYERS");
      isDepth=isDepth or nameS_start_or_end_with(*coordVar,"z-");
      check=(int)isLon+isLat+isDepth;
      }
    
    if(check!=1){
      if(verbose>=0)
        STDERR_BASE_LINE_FUNC(
          "COORDINATES ERROR WITH "+var.name+":"+att->name+" ATTRIBUTE : "+coordName+"'S NAMES NOT UNDERSTANDABLE: %d COORDINATE TYPE RECOGNISED.\n"
          "FIX WITH EITHER OF:\n"
          "  ncatted -a standard_name,"+coordName+",o,c,\"longitude\" <file>.nc\n"
          "  ncatted -a standard_name,"+coordName+",o,c,\"latitude\" <file>.nc\n"
          "  ncatted -a standard_name,"+coordName+",o,c,\"depth\" <file>.nc\n"
          ,check);
      continue;
      }
    
    int i;
    if(isLon)
      i=0;
    else if(isLat)
      i=1;
    else if(isDepth)
      i=2;
    else TRAP_ERR_EXIT(ENOEXEC,"programming error\n");
    
    vs[i]=coordVar;
    
    }while(coordName!="" && not coordStream.eof());
  
  return 0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_parse_grid_data(const string &path,const poc_var_t &var,poc_grid_data_t *gdata, int verbose, poc_global_t *global=0)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///Load grid data for a variable from a NetCDF file
/**
\param var variable
\param *gdata grid data
\param path NetCDF file path
\param verbose
\param *global if specified, pointer to use when calling poc_inq(string,poc_global_t*,int)
\returns NC_NOERR if success or the NetCDF error code if error.

The \c *gdata will have:
  - for meshes, \c gdata->tri and \c gdata->quad set to,
    respectively, \c "triangles" and \c "quadrangles" variables
  - for structured grids, \c gdata->xd and \c gdata->yd :
    + always with one dimension
    + optionally with a name
    + with data that may be fake
  - for reduced Gaussian grids, see poc_grib_get_grid_data()
*/
/*----------------------------------------------------------------------------*/
{
  int i,status;//general index,NetCDF status
  poc_global_t global_;
  const poc_dim_t *dim;
  const int ndim=var.dimensions.size();
  int nlimited=-1;
  const int vcount=gdata->vc;//<maximum number of coordinate variables
  const poc_var_t *vs[vcount],
    *&xv=vs[0],*&yv=vs[1],*&zv=vs[2],
    *&xd=vs[3],*&yd=vs[4],*&zd=vs[5],
    *&tri=vs[3],*&quad=vs[4];
  int v;
  const int sNC=6;
  int nx=-1,ny=-1,nz=-1;//numbers related to each grid dimension
  int modeH,modeV;

  aset(vs,vcount,(const poc_var_t*)0);

/*------------------------------------------------------------------------------
  poc_inq() global file informations */
  
  if(!global){
    global=&global_;
    }
  status=poc_inq(path,global,verbose);
  if(status!=NC_NOERR)NC_TRAP_ERROR(return,status,verbose,"poc_inq(\""+path+"\",) error");
  
/*------------------------------------------------------------------------------
  first look at a valid "coordinates" attribute*/
  status=poc_parse_coordinates_attribute(var,*global,vs,verbose);
  
  if(status!=0 and verbose>=0) {
    STDERR_BASE_LINE_FUNC(
      "No valid coordinates attribute found. If actually missing, it can be easily set by using:\n"
      "  ncatted -O -a coordinates,"+var.name+",o,c,\"<last axis variable name> ... <1st axis variable name>\" "+path+"\n"
      );
    }
  
/*---------------------------------------------------------------------*//**<h1>
  check number of dimensions </h1>*/
  
  if(ndim<1 or ndim>4) NC_TRAP_ERROR(return,NC_EBADDIM,verbose,"%d dimensions when only 1 to 4 are allowed",ndim);
  
  /** \NOTE As this uses poc_get_var_length(), which calls isT(),
  this will be influenced by ::timeNames . */
  poc_get_var_length(var,NULL,&nlimited);
  if(nlimited<1 or nlimited>3) NC_TRAP_ERROR(return,NC_EBADDIM,verbose,"%d space dimensions when only 1 to 3 are allowed.",nlimited);
  
/*---------------------------------------------------------------------*//**<h1>
  Find the coordinates variables </h1>*/
  {
/*---------------------------------------------------------------------*//**<h2>
  try with findVariableThasIs() </h2>
  The prefix for findVariableThasIs() is taken from \c vns :
  \code /**/ // COMPILED CODE BELOW !!!
  const char vns[][3][16]={{"lon","longitude","x"},{"lat","latitude","y"},{"dep","depth","z"}};
  /** \endcode
  For the last dimension (the depth), if findVariableThasIs() fails,
  try and find, with findStandardName(), a variable whose standard name is in \c sNs :
  \code /**/ // COMPILED CODE BELOW !!!
  const string sNs[sNC]={"ocean_levels","ocean_sigma_coordinate", "ocean_s_coordinate","ocean_sigma_z_coordinate","ocean_double_sigma_coordinate",
    "levels_depths" //<patch
    };
  /** \endcode
  see the latest CF Conventions documents http://cfconventions.org/documents.html
  */
/*---------------------------------------------------------------------*//**<h3>
  patched for unstructured grids </h3>*/
  if(nlimited<=2){
    
    for(v=0;v<2;v++){
      ///use find1DVariableThasIs() instead
      i=find1DVariableThasIs(*global,vns[v][0],verbose);
      if(i<-1){
        if(verbose>0) STDERR_BASE_LINE("trying tighter search criterion...\n");
        i=find1DVariableThasIs(*global,vns[v][1],verbose);
        if(verbose && i>=0)STDERR_BASE_LINE("\""+global->variables[i].name+"\" is the only one that complies with the tighter search criterion.\n");
        }
      if(i<0) break;
      vs[v]=&global->variables[i];
      }
    
    bool hasOnly1Dim=(nlimited==1);
    if(v<2){
      if(hasOnly1Dim)
        NC_TRAP_ERROR(return,NC_ENOTVAR,verbose,"[%d]neither %s nor %s found with find1DVariableThasIs(,,): %d",v,vns[v][0],vns[v][1],i);
      else
        goto SG;
      }
    else if(!hasOnly1Dim){
      nx=xv->dimensions.size();
      ny=yv->dimensions.size();
      if(verbose>0)STDERR_BASE_LINE("making sure this is unstructured: %d %d... ",nx,ny);
      if(nx==1 && ny==1 && xv->dimensions[0].name==yv->dimensions[0].name){
        if(verbose>0)fprintf(stderr,"Maybe?\n");
        }
      else{
        if(verbose>0)fprintf(stderr,"NO!***\n");
        goto SG;
        }
      }
    
    i=find2DVariableThasIs(*global,3,verbose);
    if(i<-1){
      if(verbose>0) STDERR_BASE_LINE("trying patch...\n");
      i=global->variables.find("element");
      if(i<0)
        i=global->variables.find("triangles");
      if(i<0)
        i=global->variables.find("tri");
      if(verbose>0 && i>=0)STDERR_BASE_LINE("\""+global->variables[i].name+"\" is the only one that complies with the patch.\n");
      }
    if(i<0){
      if(hasOnly1Dim){
        if(verbose>0)STDERR_BASE_LINE("no triangle variable found, trying quadrangles\n");
        }
      else{
        if(verbose)STDERR_BASE_LINE("no triangle variable found. Trying SG.\n");
        /* cancel */
        vs[0]=vs[1]=0;
        goto SG;
        }
      }
    if(i>=0)
      tri=&global->variables[i];
    
quadrangles:
    if(i<0){
      i=find2DVariableThasIs(*global,4,verbose);
      if(i<-1){
        if(verbose) STDERR_BASE_LINE("trying patch...\n");
        i=global->variables.find("quadrangles");
        if(verbose && i>=0)STDERR_BASE_LINE("\""+global->variables[i].name+"\" is the only one that complies with the patch.\n");
        }
      if(i<0){
        if(verbose>0)STDERR_BASE_LINE("no quadrangle variable found (%d). Trying SG.\n",i);
        /* cancel */
        vs[0]=vs[1]=0;
        goto SG;
        }
      quad=&global->variables[i];
      }
    
    if(nlimited==2){
      /* 3D */
      poc_list_t<poc_dim_t> zvDimensions=var.dimensions;
      
      for(i=0;i<zvDimensions.size();i++){
        if(isT(zvDimensions[i])){
          zvDimensions.erase(i);
          break;
          }
        }
      fakeDimVar(&gdata->zd,zvDimensions[0],"Z",verbose>6);
      
      v=2;
      int attempt;
      for(attempt=1;/**/;attempt++,zvDimensions[0].len++){
        i=findVariableThasIs(*global, zvDimensions, 4, vns[v], verbose, 1);
        if(i>=0 || attempt>=2) break;
        }
      if(i<0)STDERR_BASE_LINE("[%d]even after %d attempts, neither %s nor %s nor %s found with findVariableThasIs(): %d\n",v,attempt,vns[v][0],vns[v][1],vns[v][2],i);
      else
        zv=&global->variables[i];
      }
    
    goto init_grid_info;
    }
  
SG:
/*---------------------------------------------------------------------*//**<h3>
  for all space dimensions </h3>*/
  for(v=0;v<nlimited;v++) {
    
    /* check whether poc_parse_coordinates_attribute() has done this part of the job */
    if(vs[v]!=0) continue;
    
    i=findVariableThasIs(*global, var.dimensions, v<2?2:4, vns[v], verbose);
    if(v==2 && i>=0){
      poc_get_var_length(global->variables[i],0,&nz);
      if(nz!=1 && nz!=3){
        if(verbose>0)STDERR_BASE_LINE("%d dimensions\n",nz);
        i=-1;
        }
      }
    if(v==2 && i<0){
      if(verbose>0)STDERR_BASE_LINE("trying ocean_*_coordinate...\n");
      i=findStandardName(*global,sNs,sNC,var.dimensions,verbose);
      if(i<0) {
        NC_TRAP_ERROR(return,NC_ENOTVAR,verbose,"[%d]not found with findVariableThasIs(): %d",v,i);
        }
      }
    if(i<0) {
      NC_TRAP_ERROR(return,NC_ENOTVAR,verbose,"[%d]neither %s nor %s nor %s found with findVariableThasIs(): %d",v,vns[v][0],vns[v][1],vns[v][2],i);
      }
    
    vs[v]=&global->variables[i];
    }
  
  }/* end of sNs and co */
  
/*---------------------------------------------------------------------*//**<h1>
  Check the match between the number of dimensions of the coordinates variables </h1>*/
  
  /* use poc_get_var_length because EVERY coordinates can have 1 time frame... */
  poc_get_var_length(*xv,0,&nx);
  poc_get_var_length(*yv,0,&ny);
  
  if(nx!=ny) NC_TRAP_ERROR(return,NC_EDIMSIZE,verbose,"number of dimensions mismatch "+xv->name+":%d "+zv->name+":%d dimensions",nx,ny);
  
  if(nx>3) NC_TRAP_ERROR(return,NC_EMAXDIMS,verbose,"horizontal coordinates may have at most 2 dimensions, not %d like "+xv->name,nx);
  if(zv){
    poc_get_var_length(*zv,0,&nz);
    if(nz!=1 && nz!=3)NC_TRAP_ERROR(return,NC_EMAXDIMS,verbose,"vertical coordinates may have either 1 or 3 dimensions, not %d like "+zv->name,nz);
    }
  else{
    nz=-1;
    }
  /// Set grid_t::modeH and grid_t::modeV to the number of dimensions of the relevant coordinates variables
  modeH=nx;
  modeV=nz;
  
  if(modeH==2){/// <h2>If grid_t::modeH==2</h2>
    /// <h3>Check also dimension sizes between the surface coordinates</h3>
    for(i=0;i<nx;i++){
      /** If one coordinate size differs from the other coordinate,
      return NC_EDIMSIZE */
      if(xv->dimensions[i].len!=yv->dimensions[i].len)
        NC_TRAP_ERROR(return,NC_EDIMSIZE,verbose,"dimension %d size mismatch x:%d y:%d",i,xv->dimensions[i].len,yv->dimensions[i].len);
      }
    }
  
/*---------------------------------------------------------------------*//**<h1>
  find the dimension variables, if any </h1>*/
  switch(modeH){
  case 1:
    for(i=0;i<vcount/2;i++){
      if(vs[i]==NULL)continue;
      dim=&vs[i]->dimensions[0];//only one dimension
      vs[i+vcount/2]=global->variables.findP(dim->name);//dimension variable
      if(vs[i+vcount/2]==NULL){
        (*gdata)[i+vcount/2].info=poc_var_t() << *dim;
        }
      }
    break;
  
  case 2:
/*---------------------------------------------------------------------*//**<h2>
    If grid_t::modeH == 2 </h2>
    Assume the coordinates have the same dimensions.
    Calls findAxisDim(), patched to recognise "axis" or "content" attributes.
    If the dimension variable is not found, call fakeDimVar().
    */
    
    /* X dimension */
    dim=findAxisDim(*xv,*global,'X',&var,-1);
    if(dim==0){
      dim=&xv->dimensions[1];
      if(verbose>0)STDERR_BASE_LINE_FUNC("Failed to find X dimension. Assuming it's "+dim->name+"\n");
      }
    
    xd=global->variables.findP(dim->name);/* dimension variable */
    if(xd && xd->dimensions.size()!=1) {STDERR_BASE_LINE("\""+xd->name+"\" has %d dimensions. PATCHING\n",xd->dimensions.size()); xd=NULL;}
    if(xd==NULL){
      if(verbose>6)NC_CHKERR_LINE_FILE(NC_ENOTVAR,"Failed to find dimension variable "+dim->name);
      fakeDimVar(&gdata->xd,*dim,"X",verbose>6);
      }
    
    /* Y dimension */
    dim=findAxisDim(*xv,*global,'Y',&var,-1);
    if(dim==0){
      dim=&xv->dimensions[0];
      if(verbose>0)STDERR_BASE_LINE_FUNC("Failed to find Y dimension. Assuming it's "+dim->name+"\n");
      }
    
    yd=global->variables.findP(dim->name);/* dimension variable */
    if(yd && yd->dimensions.size()!=1) {STDERR_BASE_LINE("\""+yd->name+"\" has %d dimensions. PATCHING\n",yd->dimensions.size()); yd=NULL;}
    if(yd==NULL){
      if(verbose>6)NC_CHKERR_LINE_FILE(NC_ENOTVAR,"Failed to find dimension variable "+dim->name);
      fakeDimVar(&gdata->yd,*dim,"Y",verbose>6);
      }
    
    /* Z dimension */
    if(modeV==1){
      if(zv==NULL)break;//2D data
      dim=&zv->dimensions[0];//only one dimension
      }
    else if(modeV==3){
      dim=findAxisDim(*zv,*global,'Z',&var,-1);
      if(dim==0){
        int j;
        for(i=0;i<3;i++){
          dim=&zv->dimensions[i];
          j=xv->dimensions.find(dim->name);
          if(j<0)break;
          }
        if(i==3)NC_TRAP_ERROR(return,NC_EBADDIM,verbose,"Failed to find Z dimension.");
        if(verbose>0)STDERR_BASE_LINE_FUNC("Failed to find Z dimension. Assuming it's "+dim->name+"\n");
        }
      }
    else break;
    
    zd=global->variables.findP(dim->name);/* dimension variable */
    if(zd && zd->dimensions.size()!=1) {STDERR_BASE_LINE("\""+zd->name+"\" has %d dimensions. PATCHING\n",zd->dimensions.size()); zd=NULL;}
    if(zd==NULL){
      if(verbose>6)NC_CHKERR_LINE_FILE(NC_ENOTVAR,"Failed to find dimension variable "+dim->name);
      fakeDimVar(&gdata->zd,*dim,"Z",verbose>6);
      }
    
    break;
  
  default:
/*----------------------------------------------------------------------------*/
    TRAP_ERR_EXIT(ENOEXEC,"%s programming error: modeH=%d\n",__func__,modeH);
    }
  
init_grid_info:
/*------------------------------------------------------------------------------
  conclusion */
  for(i=0;i<vcount;i++){
    if(vs[i]==NULL) continue;///if there is data to read
    (*gdata)[i].info=*vs[i];
    }
  
  return 0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_get_grid_data(const string &path,const poc_var_t &var,poc_grid_data_t *gdata, int verbose,poc_global_t *global)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///Load grid data for a variable from a NetCDF file
/**
\param var variable
\param *gdata grid data
\param path NetCDF file path
\param verbose
\param *global if specified, pointer to use when calling poc_inq(string,poc_global_t*,int)
\returns NC_NOERR if success or the NetCDF error code if error.
*/
/*----------------------------------------------------------------------------*/
{
  int status,i;
  
  if(is_grib(path)){
#ifdef HAVE_LIBGRIB_API
    status=poc_grib_get_grid_data(path,var,gdata,verbose);
    return status;
#else
    TRAP_ERR_RETURN(ENOEXEC,1,"Please compile with GRIB_API\n");
#endif
    }
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  parse grid informations

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  status=poc_parse_grid_data(path,var,gdata,verbose,global);
  if(status) NC_TRAP_ERROR(return,status,verbose,"poc_parse_grid_data(\""+path+"\",(\""+var.name+"\"),,,) error");
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  load grid-related buffers

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  const int frame=0;
  
  for(i=0;i<gdata->vc;i++){
    poc_data_t<double> *vdata=&(*gdata)[i];
    if(vdata->info.name=="" or vdata->data!=0)
      continue;
    
    vdata->init();
    status=vdata->read_data(path,frame,0,1);
    if(status!=NC_NOERR)
      NC_TRAP_ERROR(return,status,verbose,"poc_data_t<>::read_scaled_data(\""+path+"\") error with "+vdata->info.name);
    
    /* PATCH for partially masked grids */
    
    if(vdata==&gdata->zv)
      /* skip patch for z variable */
      continue;
    
    for(int j=0;j<vdata->length;j++){
      if(vdata->data[j]==vdata->mask)
        vdata->data[j]=NAN;
      }
    vdata->mask=NAN;
    }
  
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_get_axes(const string &path,poc_var_t *var,int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// set axes attribute
/**
\param path NetCDF file path
\param *var variable
\param verbose
\return output of poc_parse_grid_data()
*/
/*----------------------------------------------------------------------------*/
{
  int status,i,c;
  poc_grid_data_t gd;
  string axes;
  
  poc_list_t<poc_dim_t> &dimensions=var->dimensions;
  
  var->axes="";
  
  for(i=0;i<dimensions.size();i++){
    if(isT(dimensions[i]))
      var->axes+="T";
    else
      var->axes+=" ";
    }
  
  if(path=="")
    return 0;
  
  if(is_grib(path))
    return 0;
  
  status=poc_parse_grid_data(path,*var,&gd,verbose);
  if(status!=NC_NOERR)NC_TRAP_ERROR(return,status,verbose,"poc_parse_grid_data(\""+path+"\",) error with "+var->name);
  
  axes=var->axes;
  
  for(i=3;i<gd.vc;i++){
    poc_data_t<double> *gdi=&gd[i];
    
    if(gdi->info.dimensions.size()>2)
      NC_TRAP_ERROR(return,NC_EDIMSIZE,verbose,"%s: "+gdi->info.name+" has %d dimensions, i.e. "+var->name+" is not structured in an unknown way",__func__,gdi->info.dimensions.size());
    
    c=dimensions.find(gdi->info.name);
    
    if(gdi->info.dimensions.size()==2){
      for(c=0;c<var->dimensions.size() and isT(var->dimensions[c]);c++);
      axes[c]='0'+gdi->info.dimensions[1].len;
      continue;
      }
    
    axes[c]=gd.axes()[i];
    }
  
  var->axes=axes;
  
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int findWord(const string &s,const string &match,int pos=0,const char *notWord=NULL)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  if(!notWord)
    notWord="([+-*/, ])";
  do{
    pos=s.find(match,pos);
    if(pos>=s.length())break;
    //STDERR_BASE_LINE("%c\n",s[pos+match.length()]);
    if(!*notWord)break;
    if(strchr(notWord,s[pos+match.length()]) &&
       (pos<=0 || strchr(notWord,s[pos-1])) )
      break;
    pos++;
    }while(1);
  
  return pos;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void replaceWord(string *s,const string &match,const string &rep,const char *notWord=NULL)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  for(int i=0;(i=findWord(*s,match,i,notWord))<s->length();i+=rep.length()){
    s->replace(i,match.length(),rep);
    //STDERR_BASE_LINE(*s+"\"\n");
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void transfer_data(double **buffer,poc_data_t<double> *src)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  *buffer=src->data;
  src->data=NULL;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void copy_data(double **buffer,const poc_data_t<double> *src)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  poc_copy(*buffer,src->data,src->length);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<class G,typename D> int poc_grid_data_to_grid_template(const string &path,const poc_global_t &global,const poc_var_t &var,G *gd,grid_t *grid,int verbose,int frame,void (*d2d)(double **,D *))

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///Convert grid data to grid
/**
\param path NetCDF file path
\param global
\param var
\param gd
\param *grid
\param verbose
\param frame
 - if <0(default), by-pass level computation
 - if >=0, compute levels at \c frame
\returns NC_NOERR if success or the NetCDF error code if error.
*/
/*----------------------------------------------------------------------------*/
{
  int status=NC_NOERR;//NetCDF status
  
  const poc_dim_t
    *xdim=&gd->xd.info.dimensions[0],
    *ydim=&gd->yd.info.dimensions[0];
  
  /* see poc_parse_grid_data() and poc_grib_get_grid_data() */
  const bool
    isReducedGaussian=
      xdim->len==ydim->len and
      xdim->name==ydim->name and
      gd->xd.data!=0 and
      gd->yd.data==0;
  
/*------------------------------------------------------------------------------
  set defaults */
  
  grid->free();
  
  grid->circular=0;
  grid->overlapped=0;
  
/*------------------------------------------------------------------------------
  transfer data */
  
  /* use poc_get_var_length because EVERY coordinates can have 1 time frame... */
  if(isReducedGaussian){
    const poc_dim_t *dim=&var.dimensions.back();
    if(dim->name=="reduced_gg")
      grid->modeH=MODEH_REDUCED_GG;
    else
      TRAP_ERR_EXIT(ENOEXEC,"not coded yet for gridType = "+dim->name+"\n");
    }
  else{
    poc_get_var_length(gd->xv.info,0,&grid->modeH);
    }
  poc_get_var_length(gd->zv.info,0,&grid->modeV);
  
  if(isReducedGaussian){
    grid->ny=xdim->len;
    poc_copy(grid->reduced_nx,gd->xd.data,grid->ny);
    grid->nx=0;
    int sum=0,i;
    grid->reduced_xi=new int[grid->ny];
    
    double *oldy=gd->yv.data;
    gd->yv.data=new double[grid->ny];
    gd->yv.length=grid->ny;
    
    for(i=0;i<grid->ny;i++){
      const int *ri=&grid->reduced_nx[i];
      gd->yv.data[i]=oldy[sum];
      grid->reduced_xi[i]=sum;
      sum+=*ri;
      updatemax(&grid->nx,*ri);
      }
    if(sum!=gd->xv.length)
      TRAP_ERR_EXIT(ENOEXEC,"programming error : %u!=%d\n",sum,gd->xv.length);
    
    delete[]oldy;
    }
  else{
    grid->nx=xdim->len;
    grid->ny=ydim->len;
    }
  if(grid->modeV>0){
    grid->nz=gd->zd.info.dimensions[0].len;
    }
  else{
    grid->modeV=-1;
    }
  
  if(gd->xv.data){
    d2d(&grid->x,&gd->xv);
    
    ///\todo if units is radians, convert to degrees
    map_recale_all_x(grid);
    poc_minmax(grid->x,gd->xv.length,gd->xv.mask,&grid->xmin,&grid->xmax);
    grid->dx=(grid->xmax-grid->xmin)/(grid->nx-1);
    }
  if(gd->yv.data){
    d2d(&grid->y,&gd->yv);
    
    ///\todo if units is radians, convert to degrees
    poc_minmax(grid->y,gd->yv.length,gd->yv.mask,&grid->ymin,&grid->ymax);
    grid->dy=(grid->ymax-grid->ymin)/(grid->ny-1);
    //TODO:all as above but for map_recale_all_x()
    }
  
  if(grid->modeH==1 &&
      var.dimensions.find(gd->xd.info.name)<var.dimensions.find(gd->yd.info.name)
    ){
    STDERR_BASE_LINE("transposing grid\n");
    map_completegridaxis_2(grid,1);
    }
  
  if(gd->zv.data){
    const poc_list_t<poc_att_t> *attributes=&gd->zv.info.attributes;
    const poc_att_t *att=attributes->findP("standard_name");
    string sN="";//standard_name
    if(att){
      sN=att->as_string();
      }
    if(frame>=0 and is_Dimensionless_Vertical_Coordinates(sN)){
      int i;
      string formula;
      poc_deque_t<poc_cdata_t*> vars;
      poc_cdata_t *varp;
      
      if(verbose)STDERR_BASE_LINE("computing sigma levels formula\n");
      if(att=attributes->findP("formula")){
        formula=att->as_string();
        }
      else if(att=attributes->findP("formula_terms")){
        istringstream terms(att->as_charp());
        string term,varname;
        /// see #Dimensionless_Vertical_Coordinates_URL
        /* NOTE: as the operators will fill indexes the way they will be found,
        make sure to re-order some operands accordingly.
        See also bug note after formula parse. */
        if(sN=="ocean_s_coordinate"){
          formula="C(k) = (1-b)*sinh(a*s(k))/sinh(a) + b*[tanh(a*(s(k)+0.5))/(2*tanh(0.5*a)) - 0.5]\n"
            "z(n,k,j,i) = eta(n,j,i)*(1+s(k)) + depth_c*s(k) + C(k)*(depth(j,i)-depth_c)";
          }
        else if(sN=="ocean_s_coordinate_g1"){
          formula="S(k,j,i) = depth_c * s(k) + C(k) * (depth(j,i) - depth_c)\n"
            "z(n,k,j,i) = S(k,j,i) + eta(n,j,i) * (1 + S(k,j,i) / depth(j,i))";
          }
        else if(sN=="ocean_s_coordinate_g2"){
          formula="S(k,j,i) = (depth_c * s(k) + C(k) * depth(j,i)) / (depth_c + depth(j,i))\n"
            "z(n,k,j,i) = eta(n,j,i) + S(k,j,i) * (eta(n,j,i) + depth(j,i))";
          }
        else TRAP_ERR_EXIT(ENOEXEC,"no support yet for \""+sN+"\". For an update, please see:\n" Dimensionless_Vertical_Coordinates_URL "\n");
        if(verbose>1)STDERR_BASE_LINE("parsing formula_terms=\"%s\":\n",att->as_charp());
        while(terms.good()){
          terms>>term>>varname;
          term.erase(term.end()-1);
          if(term=="eta"){
            if(verbose>0)STDERR_BASE_LINE("setting "+term+"=0\n");
            poc_formula_parse(vars,"eta=0");
            continue;
            }
          if(verbose>0)STDERR_BASE_LINE("loading, from \""+path+"\", \""+varname+"\" as \""+term+"\"\n");
          varp=new poc_cdata_t;
          status=varp->init(path,varname,"","T",indexes_t() << frame);
          varp->info.name=term;
          vars.push_back(varp);
          if(term=="depth" and varp->length>grid->Hsize()){
            ostringstream indexing;
            indexing << "depth=depth[crange(0&" << (grid->nx-1) << "+" << (grid->ny-1) << "j&" << grid->nx << "+" << grid->ny << "j)]";
            if(verbose>0)STDERR_BASE_LINE("PATCHING WITH "+indexing.str()+"\n");
            poc_formula_parse(vars,indexing.str().c_str());
            }
          }
        }
      else NC_TRAP_ERROR(return,NC_ENOTATT,verbose>=0,"No level formula attribute found");
      
      //clean up formula
      replaceWord(&formula,"i","");
      replaceWord(&formula,"j","");
      replaceWord(&formula,"k","");
      replaceWord(&formula,"n","");
      replaceWord(&formula,",","");
      replaceWord(&formula,"[","(","");
      replaceWord(&formula,"]",")","");
      replaceWord(&formula,"()","","");
      if(verbose>0)STDERR_BASE_LINE("parsing formula "+formula+"\n");
      
      poc_formula_parse(vars,formula.c_str());
      
      varp=poc_formula_check_var(vars,"z");
      /// \bug should transpose indexes
      if(verbose>0)STDERR_BASE_LINE("z.length=%d\n",varp->length);
      
      grid->modeV=3;
      
      grid->z=new double[varp->length];
      
      for(i=0;i<varp->length;i++){
        grid->z[i]=real(varp->data[i]);
        }
      
      gd->zv.length=varp->length;
      gd->zv.mask=real(varp->mask);
      
      for(i=0;i<vars.size();i++){
        delete vars[i];
        }
      }
    else{
      d2d(&grid->z,&gd->zv);
      }
    
    poc_minmax(grid->z,gd->zv.length,gd->zv.mask,&grid->zmin,&grid->zmax);
    grid->dz=(grid->zmax-grid->zmin)/(grid->nz-1);
    
    grid->zmask=gd->zv.mask;
    }
  else{
    grid->nz=1;
    }
  
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_grid_data_to_grid(const string &path,const poc_global_t &global,const poc_var_t &var,const poc_grid_data_t &gd,grid_t *grid,int verbose,int frame)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_grid_data_to_grid_template(path,global,var,&gd,grid,verbose,frame,copy_data);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_grid_data_to_grid(const string &path,const poc_global_t &global,const poc_var_t &var,poc_grid_data_t *gd,grid_t *grid,int verbose,int frame)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_grid_data_to_grid_template(path,global,var,gd,grid,verbose,frame,transfer_data);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_grid_data_to_mesh(const string &path,const poc_global_t &global,const poc_var_t &var,const poc_grid_data_t &gd,mesh_t *mesh,int verbose,int frame)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///Load structured grid for a variable
/**
\param path ignored NetCDF file path
\param global ignored
\param var ignored
\param gd
\param *mesh
\param verbose
\param frame ignored
\returns always 0
*/
/*----------------------------------------------------------------------------*/
{
  int i;
  
  struct timeval before;
  gettimeofday(&before);
  
  const size_t zDimSize=gd.zd.info.dimensions.size();
  const size_t triDimSize=gd.tri.info.dimensions.size();
  const size_t quadDimSize=gd.quad.info.dimensions.size();
  
  mesh->nvtxs=gd.xv.info.dimensions[0].len;
  if(triDimSize>0)
    mesh->ntriangles=gd.tri.info.dimensions[0].len;
  if(quadDimSize>0)
    mesh->nquadrangles=gd.quad.info.dimensions[0].len;
  if(zDimSize>0){
    mesh->nlayers=gd.zd.info.dimensions[0].len;
    mesh->nlevels=gd.zv.info.dimensions[0].len;
    }
  
  if(verbose>0) STDERR_BASE_LINE("%gs:mesh->vertices[%d]\n",difftime(&before),mesh->nvtxs);
  mesh->vertices=new vertex_t[mesh->nvtxs];
  
  for(i=0;i<mesh->nvtxs;i++){
    mesh->vertices[i].lon=gd.xv.data[i];
    mesh->vertices[i].lat=gd.yv.data[i];
    }
  
  if(gd.tri.data!=0){
    i=gd.tri.data[0];
    if(i>0)
      i=amin(gd.tri.data,mesh->ntriangles*3);
    if(verbose>0) STDERR_BASE_LINE("%gs:fe_set_triangles()\n",difftime(&before),mesh->nvtxs);
    fe_set_triangles(mesh,gd.tri.data,i);
    }
  
  if(gd.quad.data!=0){
    i=gd.quad.data[0];
    if(i>0)
      i=amin(gd.quad.data,mesh->nquadrangles*4);
    if(verbose>0) STDERR_BASE_LINE("%gs:fe_set_quadrangles()\n",difftime(&before),mesh->nvtxs);
    fe_set_quadrangles(mesh,gd.quad.data,i);
    }
  
  mesh->type=0;//spherical
  
  if(verbose>0) STDERR_BASE_LINE("%gs:fe_init_from_elements()\n",difftime(&before),mesh->nvtxs);
  fe_init_from_elements(mesh,verbose);
  
  if(verbose>0) STDERR_BASE_LINE("%gs.\n",difftime(before));
  return 0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_grid_data_to_mesh(const string &path,const poc_global_t &global,const poc_var_t &var,poc_grid_data_t *gd,mesh_t *mesh,int verbose,int frame,double **depths)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  const size_t zDimSize=gd->zd.info.dimensions.size();
  
  if(zDimSize>0){
    if(depths!=0){
      *depths=gd->zv.data;
      gd->zv.data=0;
      }
    }
  
  status=poc_grid_data_to_mesh(path,global,var,*gd,mesh,verbose,frame);
  
  return status;
}
