
/*******************************************************************************

  T-UGO tools, 2006-2011

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

\brief functions definitions for COMODO (COmmunauté de MODélisation Océanique) grid input/output
*/
/*----------------------------------------------------------------------------*/

#if TUGO

#include "tugo-prototypes.h"

#else

#include "version-macros.def" //for VERSION and REVISION

#include "fe-classes.h"
#include "fe-proto.h"
#include "constants.h"

#endif

#include "grd.h" //for grd_loadgrid
#include "poc-netcdf-data.hpp"


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_get_lon_lat_time(const char *path,double *lon,double *lat,double *times,poc_list_t<poc_dim_t> *dimensions)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// read longitude, latitude and times
/**
\param *path to the file
\param *lon can be \c NULL, in which case only return
\param *lat can be \c NULL, in which case only return
\param *times can be \c NULL, in which case only return
\param *dimensions can be \c NULL
\return the number of points or crash if error.
/*----------------------------------------------------------------------------*/
{
  poc_global_t global;
  int status,ncid,vid;
  int pI,pn=-1;
  const poc_var_t *var;
  const size_t *len;
  
  printf("In %s :\n",path);
  
  status=nc_open(path,NC_NOWRITE,&ncid);
  NC_TRAP_ERROR(wexit,status,1,"nc_open(\"%s\",NC_NOWRITE,) error",path);
  
  status=poc_inq(ncid,&global);
  NC_TRAP_ERROR(wexit,status,1,"poc_inq((\"%s\"),) error",path);
  
/*------------------------------------------------------------------------------
  time */
  vid=find1DVariableThasIs(global,"time",0);
  
  if(vid<=-2)
    vid=find1DVariableThasIs(global,"time_",0);
  
  if(vid<0)
    NC_TRAP_ERROR(wexit,NC_ENOTVAR,1,"No time variable found in %s",path);
  
  size_t upn;
  date_t start;
  double startd,*rtime=0,**timep;
  const bool
    doRead= (lon!=0 or lat!=0 or times!=0);
  
  if(times!=0)
    timep=&rtime;
  else
    timep=0;
  
  status=poc_gettime(ncid,vid,&start,timep,&upn);
  
  if(status!=0 or rtime==0)
    NC_TRAP_ERROR(wexit,status,1,"poc_gettime((\"%s\"),%d,,%p=&%p,) error",path,vid,timep,rtime);
  
  pn=upn;
  STDOUT_BASE_LINE("Found "+global.variables[vid].name+"[%d]. ",pn);
  
  if(not doRead)
    goto cleanUp;
  
  if(times!=0){
    valcpy(times,rtime,pn);
    deletep(&rtime);
    
    startd=cnes_time(start,'s');
    
    for(pI=0;pI<pn;pI++)
      times[pI]+=startd;
    }
  
/*------------------------------------------------------------------------------
  lat */
  vid=find1DVariableThasIs(global,"lat",1);
  if(vid<0)
    TRAP_ERR_EXIT(status,"Could not find latitude variable.\n");
  var=&global.variables[vid];
  len=&var->dimensions[0].len;
  printf("Found "+var->name+"[%u]: ",*len);fflush(stdout);
  if(pn>0 and *len!=pn)
    NC_TRAP_ERROR(wexit,NC_EDIMSIZE,1,"%d!=%u\n",pn,*len);
  printf("reading, ");fflush(stdout);
  status=poc_get_var(ncid,vid,lat);
  NC_TRAP_ERROR(wexit,status,1,"poc_get_var((\"%s\"),(\""+var->name+"\"),) error",path);
  status=poc_scale_data(lat,pn,*var);
  if(status!=NC_ENOTATT) NC_TRAP_ERROR(wexit,status,1,"poc_scale_data(,%d,(\""+var->name+"\")) error",pn);
  printf("done. ");
  
/*------------------------------------------------------------------------------
  lon */
  vid=find1DVariableThasIs(global,"lon",1);
  if(vid<0)
    TRAP_ERR_EXIT(status,"Could not find longitude variable.\n");
  var=&global.variables[vid];
  len=&var->dimensions[0].len;
  printf("Found "+var->name+"[%u]: ",*len);fflush(stdout);
  if(pn>0 and *len!=pn)
    NC_TRAP_ERROR(wexit,NC_EDIMSIZE,1,"%d!=%u\n",pn,*len);
  printf("reading, ");fflush(stdout);
  status=poc_get_var(ncid,vid,lon);
  NC_TRAP_ERROR(wexit,status,1,"poc_get_var((\"%s\"),(\""+var->name+"\"),) error",path);
  status=poc_scale_data(lon,pn,*var);
  if(status!=NC_ENOTATT) NC_TRAP_ERROR(wexit,status,1,"poc_scale_data(,%d,(\""+var->name+"\")) error",pn);
  printf("done.");
  
cleanUp:
  
  printf("\n");fflush(stdout);
  
  status=nc_close(ncid);
  NC_CHKERR_BASE_LINE(status,"nc_close((\"%s\")) error",path);
  
  if(dimensions!=0)
    *dimensions=var->dimensions;
  
  return pn;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_get_grd_PATCH(const string &path,const poc_var_t &var,grid_t *grid,int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  if(verbose>0)STDERR_BASE_LINE("grd_loadgrid PATCH\n");
  status=grd_loadgrid(path.c_str(),grid);
  if(verbose>=0)NC_CHKERR_BASE_LINE(status,"grd_loadgrid(\""+path+"\",) error");
  
  if(verbose>0)STDERR_BASE_LINE("PATCHING grid\n");
  map_completegridaxis(grid,1);
  for(int j=0;j<grid->ny/2;j++)
    swapValues(&grid->y[j],&grid->y[grid->ny-1-j]);
  
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_get_grid(const string &path,const poc_var_t &var,grid_t *grid,int verbose,int frame)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///Load structured grid for a variable
/**
\param var
\param *grid
\param path NetCDF file path
\param verbose
\param frame
 - if <0(default), by-pass level computation
 - if >=0, compute levels at \c frame
\returns NC_NOERR if success or the NetCDF error code if error.
*/
/*----------------------------------------------------------------------------*/
{
  int status,nxd,nyd;//NetCDF status
  poc_grid_data_t gdata;
  poc_global_t global;
  
  if(var.dimensions.size()==1 && var.dimensions[0].name=="xysize"){
    status=poc_get_grd_PATCH(path,var,grid,verbose);
    return status;
    }
  
  status=poc_get_grid_data(path,var,&gdata,verbose,&global);
  if(status) NC_TRAP_ERROR(return,status,verbose,"poc_get_grid_data() error");
  
  nxd=gdata.xd.info.dimensions.size();
  nyd=gdata.yd.info.dimensions.size();
  
  if(nxd!=1 || nyd!=1)
    NC_TRAP_ERROR(return,NC_EBADDIM,verbose,"wrong number of dimensions for structured grid: %d and %d",nxd,nyd);
  
  status=poc_grid_data_to_grid(path,global,var,&gdata,grid,verbose,frame);
  
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_get_discretisation(const poc_var_t &var,const mesh_t *mesh,int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int discretisation=UNSET;
  const char *discname=strrchr(var.name.c_str(),'_');
  
  if(discname!=NULL){
    discretisation=discretisation_from_name(&discname[1],verbose);
    if(discretisation==UNSET)
      discname=NULL;
    else
      if(verbose>0) STDERR_BASE_LINE("%s discretisation.\n",&discname[1]);
    }
  
  if(discname==NULL){
    string discdim=var.dimensions.back().name;
    if(replace(&discdim,"_NNODES","")>0){
      discname=discdim.c_str();
      discretisation=discretisation_from_name(discname);
      if(verbose>0) STDERR_BASE_LINE("%s discretisation.\n",discname);
      }
    }
  
  if(discname==NULL && mesh!=0){
    int length=var.dimensions.back().len;
    verbose++;
    
    if(verbose>0) STDERR_BASE_LINE("guessing discretisation from %d and ... ",length);
    if(length==mesh->nvtxs){
      if(verbose>0) fprintf(stderr,"%d=LGP1\n",mesh->nvtxs);
      discretisation=LGP1;
      }
    else if(length==mesh->ntriangles){
      if(verbose>0) fprintf(stderr,"%d=LGP0\n",mesh->ntriangles);
      discretisation=LGP0;
      }
    else if(length==mesh->nedges){
      if(verbose>0) fprintf(stderr,"%d=NCP1\n",mesh->nedges);
      discretisation=NCP1;
      }
    else if(length==mesh->nedges+mesh->nvtxs){
      if(verbose>0) fprintf(stderr,"%d=LGP2\n",mesh->nedges+mesh->nvtxs);
      discretisation=LGP2;
      }
    else{
      if(verbose>0) fprintf(stderr,"could not guess dimension!\n");
      }
    }
  
  return discretisation;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int test_grid_data(const poc_grid_data_t &gdata,bool *isGrid,bool *isMesh,int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  *isGrid=false;
  *isMesh=false;
  
  const size_t xDimSize=gdata.xd.info.dimensions.size();
  const size_t yDimSize=gdata.yd.info.dimensions.size();
  const size_t triDimSize=gdata.tri.info.dimensions.size();
  const size_t quadDimSize=gdata.quad.info.dimensions.size();
  
  if(xDimSize==1 and yDimSize==1){
    *isGrid=true;
    }
  else if(triDimSize==2 or quadDimSize==2){
    *isMesh=true;
    }
  else
    NC_TRAP_ERROR(return,NC_EVARSIZE,verbose,"error: %d %d",xDimSize,yDimSize);
  
  return 0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_get_grid_or_mesh(const string &path,const poc_var_t &var,grid_t *grid,mesh_t *mesh,int verbose,int frame,int *discretisation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///Load structured or unstructured grid for a variable
/**
\param path NetCDF file path
\param var
\param *grid will be freed if not loaded
\param *mesh will be freed if not loaded
\param verbose
\param frame
 - if <0(default), by-pass level computation
 - if >=0, compute levels at \c frame
\returns NC_NOERR if success or the NetCDF error code if error.

In the case of an unstructured 3D grid, the depths may be put in \c grid->z !
*/
/*----------------------------------------------------------------------------*/
{
  int status;//NetCDF status
  poc_grid_data_t gdata;
  poc_global_t global;
  
  grid->free();
  mesh->destroy();
  
  if(var.dimensions.size()==1 && var.dimensions[0].name=="xysize"){
    status=poc_get_grd_PATCH(path,var,grid,verbose);
    return status;
    }
  
  struct timeval before;
  gettimeofday(&before);
  
  if(verbose>0)STDERR_BASE_LINE("poc_get_grid_data(\""+path+"\",(\""+var.name+"\"),...)\n");
  
  status=poc_get_grid_data(path,var,&gdata,verbose,&global);
  if(status)NC_TRAP_ERROR(return,status,verbose,"poc_get_grid_data(\""+path+"\",...) error");
  
  bool isGrid,isMesh;
  
  status=test_grid_data(gdata,&isGrid,&isMesh,verbose);
  NC_TRAP_ERROR(return,status,verbose,"test_grid_data() error");
  
  if(isGrid){
    if(verbose>0)
      STDERR_BASE_LINE("%gs:poc_grid_data_to_grid(,,,[\""+gdata.xv.info.name+"\",\""+gdata.yv.info.name+"\"],...)\n",difftime(before));
    status=poc_grid_data_to_grid(path,global,var,&gdata,grid,verbose,frame);
    }
  else if(isMesh){
    if(verbose>0)
      STDERR_BASE_LINE("%gs:poc_grid_data_to_mesh ...\n",difftime(before));
    status=poc_grid_data_to_mesh(path,global,var,&gdata,mesh,verbose,frame,&grid->z);
    
    if(discretisation!=0){
      *discretisation=poc_get_discretisation(var,mesh,verbose);
      if(*discretisation==UNSET)
        status=NC_EDIMSIZE;
      }
    
    }
  else
    TRAP_ERR_EXIT(ENOEXEC,"programming error\n");
  
  if(verbose>0)STDERR_BASE_LINE("%gs.\n",difftime(before));
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_get_grid(const string &path,const string &varname,grid_t *grid, int verbose,int frame)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------
  wrapper for poc_get_grid(poc_var_t,grid_t*,string,int)
------------------------------------------------------------------------------*/
{
  int status;
  poc_var_t var;
  
  status=poc_inq_var(path,varname,&var,verbose);
  if(status!=NC_NOERR)
    NC_TRAP_ERROR(return,status,verbose,"poc_inq_var(\""+path+"\",\""+varname+"\",,) error");
  
  status= poc_get_grid(path,var,grid,verbose,frame);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_get_mesh(const string &path,const poc_var_t &var,mesh_t *mesh,int verbose,int frame,int *discretisation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///Load unstructured grid for a variable
/**
\param path NetCDF file path
\param var
\param *mesh will be freed if not loaded
\param verbose
\param frame
 - if <0(default), by-pass level computation
 - if >=0, compute levels at \c frame
\param *discretisation
\returns NC_NOERR if success or the NetCDF error code if error.
------------------------------------------------------------------------------*/
{
  int status;//NetCDF status
  poc_grid_data_t gdata;
  poc_global_t global;
  
  mesh->destroy();
  
  struct timeval before;
  gettimeofday(&before);
  
  if(verbose>0) STDERR_BASE_LINE("poc_get_grid_data(\""+path+"\",...)\n");
  
  status=poc_get_grid_data(path,var,&gdata,verbose,&global);
  if(status!=0) NC_TRAP_ERROR(return,status,verbose,"poc_get_grid_data(\""+path+"\",...) error");
  
  if(gdata.tri.info.dimensions.size()!=2) NC_TRAP_ERROR(return,NC_EVARSIZE,1,"error: %d %d",gdata.xd.info.dimensions.size(),gdata.yd.info.dimensions.size());
  
  if(verbose>0) STDERR_BASE_LINE("%gs:poc_grid_data_to_mesh ...\n",difftime(before));
  status=poc_grid_data_to_mesh(path,global,var,&gdata,mesh,verbose,frame);
  
  if(discretisation!=0){
    *discretisation=poc_get_discretisation(var,mesh,verbose);
    if(*discretisation==UNSET)
      status=NC_EDIMSIZE;
    }
  
  if(verbose>0) STDERR_BASE_LINE("%gs.\n",difftime(before));
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int comodo_compliance(const string &path,const vector<string> &vars,string *report,int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///Test compliance of NetCDF file to comodo standard
/**
\param path
\param vars list of variables
\param *report will contain the details
\param verbose
\returns NC_NOERR if compliant or the <i>last</i> NetCDF error code

\sa #COMODO_URL
*/
/*----------------------------------------------------------------------------*/
{
  string base=strrchr0(__FILE__,'/');
  #define BASE_LINE_STR base+":" TOSTRING(__LINE__) ":"
  #define NC_CHKERR(status) " ("+nc_strerror(status)+")"
  int i,j,status;//general indexes,NetCDF status
  poc_global_t global;
  const string *name;
  grid_t grid;
  poc_grid_data_t gdata;
  poc_list_t<poc_name_id_t> frame_names;
  
  *report="Checking compliance of "+path+"\n";
  
/*---------------------------------------------------------------------*//**<h1>
  Call poc_inq(string,poc_global_t*,int) </h1>*/
  status=poc_inq(path,&global,verbose);
  if(status!=NC_NOERR){
    *report+=BASE_LINE_STR "poc_inq() error" NC_CHKERR(status) "\n";
    return status;
    }
  
  for(i=0;i<global.dimensions.size();i++){/// <h1>For all dimensions</h1>
    name=&global.dimensions[i].name;
    /// <h2>Check if it exists</h2>
    j=global.variables.find(*name);
    if(j<0){
      *report+=BASE_LINE_STR+"\""+*name+"\" dimension variable not found\n";
      status=NC_ENOTVAR;
      }
    }
  
  if(vars.size()==0){/// <h1>If there is no variables in the list</h1>
    /// report
    *report+=BASE_LINE_STR+"no variables given to check. The variables are:";
    for(i=0;i<global.variables.size();i++){
      *report+=" "+global.variables[i].name;
      }
    *report+="\n";
    }
  
  for(i=0;i<vars.size();i++){/// <h1>For all variables</h1>
    name=&vars[i];
    /// <h2>Check if it exists</h2>
    *report+="Checking "+*name+" within "+path+"\n";
    j=global.variables.find(*name);
    if(j<0){
      *report+=BASE_LINE_STR+" variable not found\n";
      status=NC_ENOTVAR;
      continue;
      }
    
    /// <h2>Call poc_getgrid()</h2>
    status=poc_get_grid(path,global.variables[j],&grid,verbose,0);
    if(status==NC_NOERR){
      ostringstream gridString;
      grid.print(gridString,"  ");
      *report+=gridString.str();
      continue;
      }
    *report+=BASE_LINE_STR " poc_get_grid() error" NC_CHKERR(status) "\n";
    }
  
  status=poc_get_frame_names(path,&frame_names);
  if(status==NC_NOERR){
    const int n=frame_names.size();
    
    *report+="frames:";
    
    for(i=0;i<n;i++){
      asprintf(*report,"[%d]"+frame_names[i].name,i);
      }
    
    *report+="\n";
    }
  
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void poc_grid_dim(const grid_t& grid,poc_dim_t *nx,poc_dim_t *ny,const string & nameSuf)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  nx->init("nx"+nameSuf,grid.nx);
  ny->init("ny"+nameSuf,grid.ny);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_save_grid(const string& path, poc_var_t* gridVar,const grid_t& grid, int overwrite, int verbose,const string & loc,float xdimoffset,float ydimoffset,float zdimoffset,const poc_list_t<poc_att_t> *attributes)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///Save a grid
/**
\param grid
\param path a NetCDF file path
\param *gridVar set to the template variable that has the right dimensions and attributes if not NULL
\param overwrite If non-0, overwrite the file with poc_create()
\param verbose
\param loc location string
\param xdimoffset offset for X dimension variable
\param ydimoffset offset for Y dimension variable
\param zdimoffset offset for Z dimension variable
\returns NC_NOERR if success or the NetCDF error code if error.
*/
/*----------------------------------------------------------------------------*/
{
  string nameSuf,stdNameSuf;
  if(loc!=""){
    nameSuf="_"+loc;
    stdNameSuf="_at_"+loc+"_location";
    }
  
  int i,status;                 //<general index,NetCDF status
  
  poc_dim_t nx,ny,nz;           //<dimensions
  poc_data_t<float> xd,yd,zd;   //<dimension variables
  poc_var_t xv,yv,zv;           //<variables
  const poc_att_t degrees("units","degrees"),xaxis("axis","X"),yaxis("axis","Y"),zaxis("axis","Z");
  
  bool oneDim=grid.ny==1 && grid.modeH==2 && grid.modeV==0;
  
  if(oneDim){
    nx.init("n"+nameSuf,grid.nx);
    }
  else
    poc_grid_dim(grid,&nx,&ny,nameSuf);
  
  if(!oneDim){
    xd.info=poc_var_t(nx.name,NC_FLOAT);
    xd.info.attributes << xaxis;
    xd.info.dimensions << nx;
    xd.init();
    for(i=0;i<nx.len;i++){
      xd.data[i]=i+xdimoffset;
      }
    }
  xv=poc_var_t("longitude"+nameSuf,NC_DOUBLE);
  if(attributes!=0)
    xv.attributes=*attributes;
  xv.attributes << poc_att_t("standard_name","longitude"+stdNameSuf);
  xv.attributes << degrees << xaxis;
  
  if(!oneDim){
    yd.info=poc_var_t(ny.name,NC_FLOAT);
    yd.info.attributes << yaxis;
    yd.info.dimensions << ny;
    yd.init();
    for(i=0;i<ny.len;i++){
      yd.data[i]=i+ydimoffset;
      }
    }
  yv=poc_var_t("latitude"+nameSuf,NC_DOUBLE);
  if(attributes!=0)
    yv.attributes=*attributes;
  yv.attributes << poc_att_t("standard_name","latitude"+stdNameSuf);
  yv.attributes << degrees << yaxis;
  
  if(gridVar){
    *gridVar=poc_var_t();
    gridVar->dimensions << ny << nx;
    /* See the NOTE in poc_parse_coordinates_attribute() documentation
    about how to fill the coordinates attribute... */
    *gridVar << poc_att_t("coordinates",yv.name+" "+xv.name);
    }
  
  if(grid.modeV>0){
    nz.init("nz"+nameSuf,grid.nz,0);
    zd.info=poc_var_t(nz.name,NC_FLOAT);
    if(attributes!=0)
      zd.info.attributes=*attributes;
    zd.info.attributes << zaxis;
    zd.info.dimensions << nz;
    zd.init();
    for(i=0;i<nz.len;i++){
      zd.data[i]=i+zdimoffset;
      }
    zv.init("depth"+nameSuf,NC_DOUBLE,"depth"+stdNameSuf,"m");
    zv << zaxis;
    }
  
  switch(grid.modeH){
  case 1:
    xv.dimensions << nx;
    yv.dimensions << ny;
    break;
  case 2:
    if(oneDim){
      xv.dimensions << nx;
      yv.dimensions << nx;
      }
    else{
      xv.dimensions << ny << nx;
      xv.attributes << poc_att_t("content","YX");
      xv.attributes << poc_att_t("coordinates",yv.name+" "+xv.name);
      yv.dimensions << ny << nx;
      yv.attributes << poc_att_t("content","YX");;
      yv.attributes << poc_att_t("coordinates",yv.name+" "+xv.name);
      }
    break;
  default:
    TRAP_ERR_EXIT(ENOEXEC,"%s not coded for grid_t::mode=%d, use map_completegridaxis() beforehand ?\n",__func__,grid.modeH);
    }
  
  if(overwrite){
    poc_global_t global("constructed around " __LINE_FILE_PACKAGE_REVISION);
    if(!oneDim)
      global.variables << xd.info << yd.info;
    global.variables << xv << yv;
    if(grid.modeV>0){
      global.variables << zd.info << zv;
      }
    status=poc_create(path,global,verbose);
    }
  
  if(!oneDim){
    status=xd.write_data(path,0,verbose,1);//descaling
    status=yd.write_data(path,0,verbose,1);//descaling
    }
  status=poc_put_vara(path,xv,0,grid.x,verbose,1);//descaling
  status=poc_put_vara(path,yv,0,grid.y,verbose,1);//descaling
  if(grid.modeV>0){
    zv.dimensions << nz << ny << nx;
    status=poc_put_vara(path,zv,0,grid.z,verbose,1);//descaling
    status=zd.write_data(path,0,verbose,1);//descaling
    }
  
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  string getProduction(const string & def,const char *srcFile,int srcLine)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  string production;
  
  if(srcFile==0)
    production="constructed around " __LINE_FILE_PACKAGE_REVISION;
  else if(srcLine<=0)
    production="constructed in "+(string)srcFile+" of " PACKAGE_STRING " " REVISION;
  else{
    asprintf(production,"constructed around line %d of %s of " PACKAGE_STRING " " REVISION,srcLine,srcFile);
    }
  
  return production;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  poc_var_t poc_save_grid(const string& path,const grid_t& grid, const char *srcFile, int srcLine, const string & loc, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///save a grid with 0 offsets
/**
*/
/*----------------------------------------------------------------------------*/
{
  int status;
  poc_var_t gridVar;
  
  int create=(srcFile!=0);
  status=poc_save_grid(path,&gridVar, grid, create, max(verbose-1,0), loc);
  
  if(verbose && status) NC_CHKERR_LINE_FILE(status,"poc_save_grid(\""+path+"\",) error");
  
  const string
    production=getProduction(__LINE_FILE_PACKAGE_REVISION,srcFile,srcLine);
  status=poc_def_att(path,poc_att_t("production",production));
  
  return gridVar;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_save_mesh(const string & path,const mesh_t & mesh,const char *srcFile,int srcLine,int verbose,int discretisation,bool withTime,poc_var_t * var)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,j,m,status;
  
  const poc_dim_t
    timed("time",NC_UNLIMITED),
    vd("N",mesh.nvtxs),
    td("M",mesh.ntriangles),
    three("P",3);
  
  var->attributes.clear();
  var->dimensions.clear();
  
  if(withTime)
    var->dimensions << timed;
  
  switch(discretisation){
  case LGP0: var->dimensions << td; break;
  case LGP1: var->dimensions << vd; break;
  default:
    discretisation_t descriptor;
    descriptor=init_discretisation(mesh, discretisation);
    
    poc_dim_t discdim;
    discdim.name=discretisation_name(discretisation)+string("_NNODES");
    discdim.len=descriptor.nnodes;
    
    descriptor.destroy();
    
    var->dimensions << discdim;
    }
  
  poc_var_t tv;
  status=poc_inq_var(path,"triangles",&tv,-1);
  if(status==0)
    return status;
  
/*------------------------------------------------------------------------------
  overwrite */
  poc_data_t<int> triangles;
  triangles.info.init("triangles",NC_INT,"indexes_of_corners_of_elements")
    .dimensions << td << three;
  
  poc_data_t<double> lon,lat;
  lon.info.init("lon",NC_DOUBLE,"longitude","degrees")
    .dimensions << vd;
  lat.info.init("lat",NC_DOUBLE,"latitude","degrees")
    .dimensions << vd;
  
  const string
    production=getProduction(__LINE_FILE_PACKAGE_REVISION,srcFile,srcLine);

  poc_global_t global(production);
  
  if(withTime){
    poc_dim_t td("time",NC_UNLIMITED);
    poc_var_t timev;
    timev.init(td,NC_DOUBLE);
    timev << poc_att_t("units","seconds since 1950-01-01 00:00:00");
    global.variables << timev;
    }
  
  global.variables << lon.info << lat.info << triangles.info;
  status=poc_create(path,global,verbose);
  
/*------------------------------------------------------------------------------
  data */
  lon.init();
  lat.init();
  for(i=0;i<mesh.nvtxs;i++){
    vertex_t *vertex=&mesh.vertices[i];
    lon.data[i]=vertex->lon;
    lat.data[i]=vertex->lat;
    }
  status=lon.write_data(path);
  status=lat.write_data(path);
  
  triangles.init();
  j=0;
  for(m=0;m<mesh.ntriangles;m++){
    triangle_t *triangle=&mesh.triangles[m];
    for(i=0;i<3;i++,j++)
      triangles.data[j]=triangle->vertex[i];
    }
  status=triangles.write_data(path);
  
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T> int poc_put_mesh_vara_template(const string & path,const mesh_t & mesh,int discretisation,poc_var_t *var,int frame,const T *buffer,int verbose=1)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  if(verbose){STDOUT_BASE_LINE("saving to "+path+"("+var->name+"[%d])... ",frame);fflush(stdout);}
  
  const bool
    withTime=not updatemax(&frame,0);
  
  status=poc_save_mesh(path,mesh,__FILE__,__LINE__,verbose,discretisation,withTime,var);
  if(status!=NC_NOERR) NC_TRAP_ERROR(return,status,verbose,"poc_save_mesh(\""+path+"\",...) error");
  
  status=poc_def_var(path,*var);
  if(status!=NC_NOERR) NC_TRAP_ERROR(return,status,verbose,"poc_def_var(\""+path+"\",(\""+var->name+"\"),...) error");
  
  status=poc_put_vara(path,*var,frame,buffer,verbose);
  if(status!=NC_NOERR) NC_TRAP_ERROR(return,status,verbose,"poc_put_vara(\""+path+"\",(\""+var->name+"\"),...) error");
  
  if(verbose) {printf("%d done.\n",status);fflush(stdout);}
  return status;
}
#if 0
extern int poc_put_mesh_vara(const string & path,const mesh_t & mesh,int discretisation,poc_var_t *var,int frame,const T *buffer,int verbose=1);


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_put_mesh_vara(const string & path,const mesh_t & mesh,int discretisation,poc_var_t *var,int frame,const T *buffer,int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_put_mesh_vara_template(path,mesh,discretisation,var,frame,buffer,verbose);
  return status;
}
#endif


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_put_mesh_vara(const string & path,const mesh_t & mesh,int discretisation,poc_var_t *var,int frame,const double *buffer,int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_put_mesh_vara_template(path,mesh,discretisation,var,frame,buffer,verbose);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_put_mesh_vara(const string & path,const mesh_t & mesh,int discretisation,poc_var_t *var,int frame,const float *buffer,int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_put_mesh_vara_template(path,mesh,discretisation,var,frame,buffer,verbose);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_put_mesh_vara(const string & path,const mesh_t & mesh,int discretisation,poc_var_t *var,int frame,const int *buffer,int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_put_mesh_vara_template(path,mesh,discretisation,var,frame,buffer,verbose);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T> int poc_put_mesh_vara_template(const string & path,const mesh_t & mesh,int discretisation,const string & name,int frame,const T *buffer,nc_type type,const string & units="",int verbose=1)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  poc_var_t var;
  
  var.init(name,type,comodo_standard_name(),units);
  
  status=poc_put_mesh_vara(path,mesh,discretisation,&var,frame,buffer,verbose);
  
  return status;
}
#if 0
extern int poc_put_mesh_vara(const string & path,const mesh_t & mesh,int discretisation,const string & name,int frame,const T *buffer,const string & units="",int verbose=1);


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_put_mesh_vara(const string & path,const mesh_t & mesh,int discretisation,const string & name,int frame,const T *buffer,const string & units,int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_put_mesh_vara_template(path,mesh,discretisation,name,frame,buffer,NC_T,units,verbose);
  return status;
}
#endif


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_put_mesh_vara(const string & path,const mesh_t & mesh,int discretisation,const string & name,int frame,const double *buffer,const string & units,int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_put_mesh_vara_template(path,mesh,discretisation,name,frame,buffer,NC_DOUBLE,units,verbose);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_put_mesh_vara(const string & path,const mesh_t & mesh,int discretisation,const string & name,int frame,const float *buffer,const string & units,int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_put_mesh_vara_template(path,mesh,discretisation,name,frame,buffer,NC_FLOAT,units,verbose);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_put_mesh_vara(const string & path,const mesh_t & mesh,int discretisation,const string & name,int frame,const int *buffer,const string & units,int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=poc_put_mesh_vara_template(path,mesh,discretisation,name,frame,buffer,NC_INT,units,verbose);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  string comodo_standard_name(const string & varName,const string waveName)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// compute standard name from common variable name
/**
\param varName common variable name, e.g. "H","U,"V","Hg","Ha","Ug","V_a",...
\return standard name
*/
/*----------------------------------------------------------------------------*/
{
  if(varName=="")
    return "_"+waveName;
  
  string standard_name;
  const int vnl=varName.length();
  const bool
    hasPartName=vnl>1,
    hasWaveName=waveName.length()>0;
  
  switch( toupper(varName[0]) ){
  case 'H':
  case 'E':
  case 'P':
    standard_name="sea_surface_height_above_geoid";
    break;
  case 'U':
    standard_name="sea_water_x_velocity";
    break;
  case 'V':
    standard_name="sea_water_y_velocity";
    break;
  default:
    return "";
    }
  
  if(hasPartName)
    switch( toupper(varName[vnl-1]) ){
    case 'A':
      standard_name+="_amplitude";
      break;
    case 'G':
      standard_name+="_phaselag";
      break;
    default:
      return "";
      }
  
  if(hasWaveName)
    standard_name+="_at_"+waveName+"_frequency";
  
  if(hasPartName || hasWaveName)
    standard_name+="_due_to_non_equilibrium_ocean_tide";
  
  return standard_name;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  string long_name_from_varname(const string & varName)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// compute long name for Ferret from common variable name
/**
\param varName common variable name, e.g. "H","U,"V","Hg","Ha","Ug","V_a",...
\return long name
*/
/*----------------------------------------------------------------------------*/
{
  if(varName=="")
    return "_";
  
  string long_name;
  const int vnl=varName.length();
  const bool hasPartName=vnl>1;
  
  switch( toupper(varName[0]) ){
  case 'H':
  case 'E':
  case 'P':
    long_name="sea_surface_height";
    break;
  case 'U':
    long_name="x_velocity";
    break;
  case 'V':
    long_name="y_velocity";
    break;
  default:
    return "";
    }
  
  if(hasPartName)
    switch( toupper(varName[vnl-1]) ){
    case 'A':
      long_name+="_amplitude";
      break;
    case 'G':
      long_name+="_phase";
      break;
    default:
      return "";
      }
  
  return long_name;
}
