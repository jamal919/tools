
/*******************************************************************************

  T-UGO tools, 2006-2011

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

\brief poc_data_t operators
*/
/*----------------------------------------------------------------------------*/


#include <unistd.h> //unlink

#include "poc-data-operators.hpp"
#include "map.h"


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void poc_save_grid(poc_deque_t<poc_cdata_t*> &vars,const grid_t& grid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  poc_dim_t nx,ny;           //<dimensions
  poc_var_t xv,yv;           //<variables
  const poc_att_t
    degrees("units","degrees"),xaxis("axis","X"),yaxis("axis","Y");
  
  poc_grid_dim(grid,&nx,&ny);
  
  xv.dimensions<<ny<<nx;
  xv.attributes<<poc_att_t("standard_name","longitude");
  xv.attributes<<degrees<<xaxis;
  
  yv.dimensions<<ny<<nx;
  yv.attributes<<poc_att_t("standard_name","latitude");
  yv.attributes<<degrees<<yaxis;
  
  const int nxy=grid.Hsize();
  
  poc_cdata_t *lon,*lat;
  lon=&poc_formula_get_var(vars,POC_LON_VAR_NAME);
  lat=&poc_formula_get_var(vars,POC_LAT_VAR_NAME);
  
  lon->info.dimensions=xv.dimensions;
  lon->info.attributes=xv.attributes;
  lon->init();
  if(grid.modeH!=2)
    TRAP_ERR_EXIT(ENOEXEC,"not coded yet for grid_t::modeH=%d\n",grid.modeH);
  valcpy(lon->data,grid.x,nxy);
  
  lat->info.dimensions=yv.dimensions;
  lat->info.attributes=yv.attributes;
  lat->init();
  valcpy(lat->data,grid.y,nxy);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int get_poc_cdata_dims(const poc_cdata_t & var,int *nx,int *ny)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///
/**
\param var
\return 0 on success and -1 if \c var has more than 2 dimensions
*/
/*----------------------------------------------------------------------------*/
{
  const poc_list_t<poc_dim_t> &dimensions=var.info.dimensions;
  const int n=dimensions.size();
  
  switch(n){
  case 0:
    *nx=1;
    *ny=1;
    break;
  case 1:
    *nx=dimensions[0].len;
    *ny=1;
    break;
  case 2:
    *nx=dimensions[1].len;
    *ny=dimensions[0].len;
    break;
  default:
    TRAP_ERR_RETURN(-1,1,"%s called with a %d-dimension variable when only 0, 1 or 2-dimension variables are allowed.\n",__func__,n);
    }
  
  return 0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_cdata2grid(grid_t * grid,const poc_cdata_t & var)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///
/**
\param *grid
\param var
\return 0 on success and -1 if \c var has more than 2 dimensions
*/
/*----------------------------------------------------------------------------*/
{
  int status=0,nx,ny,k;
  
  status=get_poc_cdata_dims(var,&nx,&ny);
  if(status!=0) return status;
  
  grid->free();
  map_allocate_x_y(grid,2,nx,ny);
  
  for(k=0;k<var.length;k++){
    const complex<double> *d=&var.data[k];
    grid->x[k]=real(*d);
    if(isnan(real(*d)))
      grid->y[k]=NAN;
    else
      grid->y[k]=imag(*d);
    }
  
  return 0;
}


#define VARIABLE_LIST_NAME "variable_list"


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_save_vars(const string & path,const poc_deque_t<poc_cdata_t*> & vars,int verbose,bool saveAsGrid,bool overwrite)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status,length,i,j;
  bool isComplex;
  vector<string> savedVars;
  
  if(overwrite){
    poc_global_t global("constructed around " __LINE_FILE_PACKAGE_REVISION);
    status=unlink(path.c_str());
    status=poc_create(path,global,verbose,0);
    }
  
/*------------------------------------------------------------------------------
  save data */
  for(i=0;i<vars.size();i++){
    const poc_cdata_t *varp=vars[i];
    const poc_var_t *info=&varp->info;
    range_t<double> rR;
    double f;/* fractional part */
    
    length=varp->length;
    isComplex=varp->isComplex(&rR,&f);
    
    if(isComplex){
      /* is complex */
      savedVars.push_back(info->name);
      
      if(saveAsGrid){
        grid_t grid;
        poc_cdata2grid(&grid,*varp);
        poc_save_grid(path,0,grid,0,verbose,info->name,0.,0.,0.,&info->attributes);
        grid.free();
        continue;
        }
      
      poc_var_t ainfo=*info;
      ainfo.name+=AMPSUF;
      poc_var_t ginfo=*info;
      ginfo.name+=PHASUF;
      STDERR_BASE_LINE(""+info->name+"[%d] is a complex\n",varp->length);
      status=poc_put_cvara(path,ainfo,ginfo,0,varp->data,verbose);
      }
    else{
      /* is real */
      double *rdata=new double[length];
      bool allnan=true;
      for(j=0;j<length;j++){
        rdata[j]=real(varp->data[j]);
        if(allnan && !isnan(rdata[j]))
          allnan=false;
        }
      if(allnan)
        STDERR_BASE_LINE("skipping "+info->name+" as it is made of nan\n");
      else{
        if(f==0.){
          /* is integer */
          poc_var_t iinfo=*info;
          rR<<1.;
          const double rRa=max(-rR.min,rR.max);
          const size_t
            size=ceil((log2(rRa)+1.)/8.),
            log2size=ceil(log2(size));
          STDERR_BASE_LINE(""+info->name+"[%d] is a %u byte (%u) integer in [%g;%g](%g)\n",varp->length,size,log2size,rR.min,rR.max,rRa);
          switch(log2size){
          case 0:
            iinfo.type=NC_BYTE;
            break;
          case 1:
            iinfo.type=NC_SHORT;
            break;
          case 2:
            iinfo.type=NC_INT;
            break;
          case 3:
            iinfo.type=NC_INT64;
            break;
          default:
            TRAP_ERR_EXIT(ENOEXEC,"not coded for %g bytes yet",exp2(log2size));
            }
          poc_att_t *att=iinfo.attributes.findP(_FillValue);
          if(att!=0)
            iinfo.type=att->type;
          status=poc_put_vara(path,iinfo,-2,rdata,verbose,1);
          }
        else if( info->dimensions.size()==2 and (info->name=="lon" or info->name=="lat") ){
          /* is 2D coordinates */
          const int
            ni=info->dimensions[1].len,
            nj=info->dimensions[0].len;
          int i,j,m,mj,m0;
          bool compressible=true;
          const bool isLon= info->name=="lon";
          
          m=ni;
          for(j=1;j<nj and compressible;j++){
            mj=m;
            
            for(i=0;i<ni;i++,m++){
              
              if(isLon)
                m0=i;
              else
                m0=mj;
              
              if(rdata[m0]==rdata[m]) continue;
              
              compressible=false;
              break;
              }
            }
          
          poc_var_t iinfo=*info;
          size_t length=varp->length;
          
          if(compressible){
            STDERR_BASE_LINE("COMODO-COMPRESSING "+info->name+"\n");
            if(isLon){
              iinfo.dimensions.erase(0);
              length=ni;
              }
            else{
              iinfo.dimensions.erase(1);
              length=nj;
              
              m=ni;
              for(j=1;j<nj;j++,m+=ni){
                rdata[j]=rdata[m];
                }
              }
            }
          
          STDERR_BASE_LINE(""+info->name+"[%u] is a floating point coordinate\n",length);
          status=poc_put_vara(path,iinfo,-2,rdata,verbose,1);
          }
        else{
          /* is floating point */
          STDERR_BASE_LINE(""+info->name+"[%d] is a floating point\n",varp->length);
          status=poc_put_vara(path,*info,-2,rdata,verbose,1);
          }
        savedVars.push_back(info->name);
        }
      delete[]rdata;
      }
    
    }
  
/*------------------------------------------------------------------------------
  save variable name list */
  string varNameList;
  
  for(i=0;i<savedVars.size();i++){
    const string &name=savedVars[i];
    j=name.find(' ');
    if(j==string::npos)
      varNameList+=" "+name;
    }
  
  if(varNameList.length()>1)
    poc_def_att(path,poc_att_t(VARIABLE_LIST_NAME,varNameList.substr(1)));
  
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_load_vars(const string & path,poc_deque_t<poc_cdata_t*> *vars,int verbose,string *formula)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status,m;
  
  poc_global_t global;
  status=poc_inq(path,&global,verbose);
  if(status!=0) NC_TRAP_ERROR(return,status,1,"poc_inq(\""+path+"\",) error");
  
  poc_att_t *varNameListAtt;
  varNameListAtt=global.attributes.findP(VARIABLE_LIST_NAME);
  if(varNameListAtt==0) NC_TRAP_ERROR(return,NC_ENOTATT,1,"no \"" VARIABLE_LIST_NAME "\" attribute in "+path);
  if(formula!=0){
    poc_att_t *formulaAtt;
    formulaAtt=global.attributes.findP("formula");
    if(formulaAtt!=0)
      *formula=formulaAtt->as_string();
    }
  
  istringstream varNameStream(varNameListAtt->as_charp());
  
  while(varNameStream.good() && !varNameStream.eof()){
    string varName;
    varNameStream>>varName;
    
    poc_var_t *fileVar0,*fileVar1=0;
    poc_cdata_t *varp=&poc_formula_get_var(*vars,varName);
    
    do{
      fileVar0=global.variables.findP(varName);
      if(fileVar0!=0){
        if(verbose) STDERR_BASE_LINE_FUNC("loading "+varName+" from "+path+"("+fileVar0->name+") ");
        status=varp->init(path,fileVar0->name);
        if(verbose) fprintf(stderr,"(%d)\n",status);
        break;
        }
      
      fileVar0=global.variables.findP(varName+AMPSUF);
      fileVar1=global.variables.findP(varName+PHASUF);
      if(fileVar0!=0 && fileVar1!=0){
        if(verbose) STDERR_BASE_LINE_FUNC("loading "+varName+" from "+path+"("+fileVar0->name+","+fileVar1->name+") ");
        status=varp->init(path,fileVar0->name,fileVar1->name);
        varp->info.name=varName;
        if(verbose) fprintf(stderr,"(%d)\n",status);
        break;
        }
      
      fileVar0=global.variables.findP("longitude_"+varName);
      fileVar1=global.variables.findP("latitude_" +varName);
      if(fileVar0!=0 && fileVar1!=0){
        if(verbose) STDERR_BASE_LINE_FUNC("loading "+varName+" as grid from "+path+"("+fileVar0->name+","+fileVar1->name+")\n");
        
        poc_data_t<double> lon,lat;
        status=lon.init(*fileVar0,"");
        status=lat.init(*fileVar1,"");
        if(lon.info.dimensions.size()!=lat.info.dimensions.size()) TRAP_ERR_EXIT(ENOEXEC,"not coded yet for this case\n");
        
        poc_list_t<poc_dim_t> &dimensions=varp->info.dimensions;
        dimensions=lon.info.dimensions;
        varp->init();
        
        status=lon.read_data(path);
        status=lat.read_data(path);
        for(m=0;m<varp->length;m++){
          complex<double> *z=&varp->data[m];
          *z=complex<double>(lon.data[m],lat.data[m]);
          }
        
        fileVar0->attributes.erase("axis");
        
        break;
        }
      
      /* just in case only the latitude has been found */
      fileVar0=0;
      
      }while(0);
    
    if(fileVar0!=0){
      fileVar0->attributes.erase("standard_name");
      varp->info.attributes=fileVar0->attributes;
      continue;
      }
    
    status=NC_ENOTVAR;
    NC_CHKERR_BASE_LINE(status,"no file variable for formula variable \""+varName+"\"");
    }
  
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  poc_cdata_t *poc_formula_check_var(const poc_deque_t<poc_cdata_t*> &vars,const string & name, int *index)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  poc_cdata_t *varp;
  
  for(int i=0;i<vars.size();i++){
    varp=vars[i];
    if(varp->info.name==name){
      if(index!=0)
        *index=i;
      return varp;
      }
    }
  
  if(index!=0)
    *index=-1;
  
  return 0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  poc_cdata_t & poc_formula_get_var(poc_deque_t<poc_cdata_t*> &vars,const string & name)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  poc_cdata_t *varp;
  
  /* search existing */
  varp=poc_formula_check_var(vars,name);
  if(varp!=0)
    return *varp;
  
  /* create new */
  varp=new poc_cdata_t;
  varp->info.init(name,NC_DOUBLE);
  
  vars.push_back(varp);
  return *varp;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  complex<double> * poc_formula_init_var(poc_deque_t<poc_cdata_t*> &vars,const string & name)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  poc_cdata_t &var=poc_formula_get_var(vars,name);
  
  if(var.length==0)
    var=NAN;
  
  return &var.data[0];
}
