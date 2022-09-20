
/*******************************************************************************

  T-UGO tools, 2006-2014

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

\brief poc_data_t and related classes
*/
/*----------------------------------------------------------------------------*/
/** \page tips Coding tips

\section oneLineRead read a buffer in one line of code
\code
  int status,frame;
  string path,varname,ampname,phaname;
  grid_t grid;
  poc_data_t<double> scalarData;
  
  status=poc_get_grid(path,varname,grid);
  status=scalarData.read_data(path,varname,frame);
  for(i=0;i<scalarData.length;i++){
    ... scalarData.data[i] ...
    }
  
  complex<double> *complexBuf;
  complexBuf=new complex<...
  status=poc_get_vara(path,ampname,phaname,frame,complexBuf);
\endcode

\section oneLineSave save a buffer in one line of code
\code
  int status;
  string outputPath;
  double *buffer;
  poc_var_t var;
  grid_t grid;
  var=poc_save_grid(outputPath,grid,__FILE__,__LINE__);
  ...
  status=poc_put_var(outputPath,var.init("name1",NC_DOUBLE,"standard_name1","units1",mask1),buffer,1);
  ...
  status=poc_put_var(outputPath,var.init("name2",NC_DOUBLE,"standard_name2","units2",mask2),buffer,1);
\endcode
*/
/*----------------------------------------------------------------------------*/

#if POC_NETCDF_DATA_HPP == 0
#define POC_NETCDF_DATA_HPP 1

#include "constants.h"

#include "poc-netcdf-io.h"
#include "archive.h" /* for quoddy_load*() */


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T> class poc_data_t

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///
/**
\note You may NOT put a poc_data_t directly or indirectly in a vector<>,
as resizing the vector<> will call poc_data_t::~poc_data_t() !
USE poc_data_t* INSTEAD !!!
*/
/*----------------------------------------------------------------------------*/
  {
public:
  poc_var_t info;
  T *data,mask;
  double scale,offset;
  size_t length;     ///<length of data
  int nlimited;   ///<number of limited dimensions
  int nframes;    ///<number of frames
  string axes;    ///< "T":time frame,"XYZ":time series,"":all the data

  poc_data_t()//it is not necessary to initialise info as its default constructor will be called
    {
    data=0;
    scale=offset=NAN;
    mask=NAN;
    length=0;
    axes="T";
    }
  
  int decode_mask(){
    int status;
    status=poc_decode_mask(info,&scale,&offset,&mask);
    return status;
    }

protected:
  
  int init_info(const poc_global_t &global, const string &var,const string & axes0){
    int status=global.variables.find(var,&info);
    if(status<0)return NC_ENOTVAR;
    return init_info(axes0);
    }
  
  int init_info(const poc_var_t &info0,const string & axes0){
    info=info0;
    return init_info(axes0);
    }
  
  int init_info(const string & axes0){
    int status,i;
    status=decode_mask();
    
//     printf("%s #1\n",__func__);
    
    status=poc_get_var_length(info,&length,&nlimited,&nframes);
    axes=axes0;
    
//     printf("%s #2\n",__func__);
    
    if(axes=="T")
      return status;

//     printf("%s #3\n",__func__);
    
    length*=nframes;
    nframes=1;
    nlimited++;
    
    for(i=0;i<axes.size();i++){
      const char c=axes[i];
      const size_t ic=info.axes.find(c);
      if(ic==string::npos) continue;
      
      const int l=info.dimensions[ic].len;
      length/=l;
      nframes*=l;
      nlimited--;
      }
    
    return 0;
    }
  
  void init_data(){
    destroy_data();
    if(length>0)
      data=new T[length];
    }

public:
  
  int init(const poc_global_t &global, const string &var,const string & axes0="T"){
    int status=init_info(global,var,axes0);
    if(status!=0) return status;
    init_data();
    return status;
    }
  
  int init(const poc_var_t &info0,const string & axes0="T"){
    int status=init_info(info0,axes0);
    init_data();
    return status;
    }
  
  int init(int mode0) __attribute__(( warning("badly deprecated") )) {
    switch(mode0){
    case 1:
      axes="T";
      break;
    case 0:
      axes="";
      break;
    default:
      TRAP_ERR_EXIT(ENOEXEC,"programming error:%d\n",mode0);
      }
    return init(axes);
    }
  
  int init(const string & axes0="T"){
    int status=init_info(axes0);
    init_data();
    return status;
    }
  
  int init(const string &filename, const string &var,int verbose=0,const string & axes0="T"){
    int status;
    status=poc_inq_var(filename,var,&info,verbose);
    if(status!=0) return status;
    status=init_info(axes0);
    if(status!=0) return status;
    init_data();
    return status;
    }
  
  int scale_data(int verbose=0){
    int status;
    
    status=poc_scale_data(data,length,scale,offset,&mask,verbose);
    
    return status;
    }
  
  int read_data(const string & filename,const indexes_t & indexes= indexes_t(), int verbose=0, int do_scale_data=0){
    int status;
    
    status=poc_get_vara(filename,info,data,axes,indexes,verbose);
    
    if(do_scale_data)
      scale_data(verbose);
    
    return status;
    }
  
  int read_data(int ncid,const indexes_t & indexes= indexes_t(), int verbose=0, int do_scale_data=0){
    int status;
    
    status=poc_get_vara(ncid,info,data,axes,indexes);
    
    if(do_scale_data)
      scale_data(verbose);
    
    return status;
    }
  
  int read_data(const string & filename,int index, int verbose=0, int do_scale_data=0){
    int status;
    status=read_data(filename,indexes_t()<<index,verbose,do_scale_data);
    return status;
    }
  
  int read_data(int ncid,int index, int verbose=0, int do_scale_data=0){
    int status;
    status=read_data(ncid,indexes_t()<<index,verbose,do_scale_data);
    return status;
    }
  
  template<typename F> int read_scaled_data(F filenameORncid, int frame=0,int verbose=0){
    int status;
    timeval before;
    
    if(verbose){
      gettimeofday(&before);
      }
    
    status=read_data(filenameORncid, frame, verbose);
    if(status!=0) return status;
    
    status=decode_mask();
    
    if(verbose)STDERR_BASE_LINE("%s:%gs\n",__func__,difftime(before));
    
    poc_scale_data(data,length,scale,offset,&mask,verbose);
    
    return 0;
    }
  
  int write_data(const string &filename,int frame=0,int verbose=0,int unscale=0) const{
    int status;
    status=poc_put_vara(filename,info,frame,data,verbose,unscale);
    return status;
    }
  
  void transfer(poc_data_t<T> * src){
    if(src->data!=0)
      init_info(src->info,src->axes);
    data=src->data;
    src->data=0;
    }
  
  void destroy_data(){
    deletep(&data);
    }
  
  ~poc_data_t(){
    this->destroy_data();
    }
  
  void init(const poc_data_t<T> & src){
    /** You will need to call poc_data_t::info.init() if you want to change the name */
    init(src.info,src.axes);
    if(src.data)
      valcpy(data,src.data,length);
    else if(length>0)
      *data=NAN;
    
    /* this is because read_data may change the mask if the scale option is set */
    mask=src.mask;
    }
  
  int init(const T & x){
    if(info.name==""){
      ostringstream oss;
      oss<<x;
      info.name=oss.str();
      }
    info=poc_var_t(info.name,NC_NAT);
    info<<poc_att_t("formula","=");
    int status=init();
    data[0]=x;
    return status;
    }
  
  poc_data_t<T> & operator=(const poc_data_t<T> & src){
    init(src);
    return *this;
    }
  
  poc_data_t<T> & operator=(const T & x){
    init(x);
    return *this;
    }
  
  int translate_index(const poc_list_t<poc_dim_t> & dimensions,int index) const{
    int translated=0,i,j,l;
    const int n1=dimensions.size();
    const int n2=info.dimensions.size();
    vector<int> indexes(n1);
    
    for(i=n1-1;i>=0;i--){
      const poc_dim_t *dim=&dimensions[i];
      
      if(axes=="T"){
        if( isT(*dim) )
          continue;
        }
      else TRAP_ERR_EXIT(ENOEXEC,"not coded yet for .axes="+axes+"\n");
      
      l=dim->len;
      indexes[i]=index % l;
      index/=l;
      }
    
    for(i=0;i<n2;i++){
      const poc_dim_t *dim=&info.dimensions[i];
      
      if(axes=="T"){
        if( isT(*dim) )
          continue;
        }
      else TRAP_ERR_EXIT(ENOEXEC,"not coded yet for .axes="+axes+"\n");
      
      if(i>0)
        translated*=dim->len;
      
      j=dimensions.find(dim->name);
      if(j<0)continue;
      
      translated+=indexes[j];
      }
    
    if(translated<0 or length<=translated)
      TRAP_ERR_EXIT(ENOEXEC,"%d not in [0;%d[\n",translated,length);
    
    return translated;
    }
  };


extern void fakeDimVar(poc_data_t<double>*dv,const poc_dim_t&dim,const char*axis,int verbose=0);


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  class poc_grid_data_t

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  {
public:
  static const int vc=6;
  static const char *axes(){return "XYZXYZ";}
  
protected:
  poc_data_t<double> vs[vc];

public:
  poc_data_t<double> &xv,&yv,&zv;  // variables
  poc_data_t<double> &xd,&yd,&zd;  // dimension variables
  poc_data_t<double> &tri,&quad;   // unstructured grid variables
  
  poc_grid_data_t():
    //initialising references
    xv(vs[0]),yv(vs[1]),zv(vs[2]),
    xd(vs[3]),yd(vs[4]),zd(vs[5]),
    tri(vs[3]),quad(vs[4])
    {}
  
  void destroy(){
    int i;
    for(i=0;i<vc;i++){
      poc_data_t<double> *vsi=&vs[i];
      vsi->destroy_data();
      vsi->info.dimensions.clear();
      }
    }
  
  void transfer(poc_grid_data_t * src){
    int i;
    
    for(i=0;i<vc;i++){
      vs[i].transfer(&src->vs[i]);
      }
    }
  
  poc_data_t<double>& operator[](int i){
    return vs[i];
    }
  
  const poc_data_t<double>& operator[](int i) const {
    return vs[i];
    }
  
  bool operator==(const poc_grid_data_t & src) const{
    int i,m;
    
    for(i=0;i<vc;i++){
      const poc_data_t<double> *vi=&vs[i];
      const poc_data_t<double> *svi=&src.vs[i];
      
      const int length=vi->length;
      
      if(length!=svi->length)
        return false;
      
      const double *vid=vi->data;
      const double *svid=svi->data;
      
      for(m=0;m<length;m++){
        if(vid[m]!=svid[m])
          return false;
        }
      }
    
    return true;
    }
  
  bool operator!=(const poc_grid_data_t & src) const{
    bool equal;
    equal=(src==*this);
    return not equal;
    }
  
  };


#define AMPSUF "_a"
#define PHASUF "_G"


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

class poc_cdata_t : public poc_data_t<complex<double> >

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  {
private:
  
  void init_type(){
    info.type=NC_DOUBLE;
    }
  
public:
  
  poc_cdata_t(){
    init_type();
    }
  
  int init(const string &axes0="T"){
    return poc_data_t<complex<double> >::init(axes0);
    }
  
  poc_cdata_t & operator=(const complex<double> & x){
    if(info.name==""){
      const bool
        hasReal=real(x)!=0.,
        hasImag=imag(x)!=0.;
      
      ostringstream oss;
      if(hasReal or not hasImag)
        oss<<real(x);
      if(hasReal and hasImag)
        oss<<showpos;
      if(hasImag)
        oss<<imag(x)<<'j';
      info.name=oss.str();
      }
    poc_data_t<complex<double> >::init(x);
    init_type();
    return *this;
    }
  
  poc_cdata_t & operator=(const string & s){
    info.attributes.clear();
    info.dimensions.clear();
    poc_dim_t dim(s.size());
    info.dimensions<<dim;
    init("");
    info.type=NC_CHAR;
    valcpy(data,s,dim.len);
    }
  
  string toStr() const{
    if(info.type!=NC_CHAR)
      return "";
    
    const size_t n=info.dimensions[0].len;
    
    string s;
    s.resize(n);
    
    for(int i=0;i<n;i++)
      s[i]=real(data[i]);
    
    return s;
    }
  
  int init(const string &filename, const string &aname, const string &gname="",const string & axes0="T",const indexes_t & indexes= indexes_t()<<0){
    int status;
    
    status=poc_inq_var(filename,aname,&info);
    if(status!=0) return status;
    
    status=init_info(axes0);
    if(status!=0) return status;
    init_data();
    
    if(gname==""){
      double *rdata,rmask=NC_FILL_DOUBLE;
      rdata=new double[length];
      
      status=poc_get_vara(filename,info,rdata,axes,indexes);
      
      rmask=real(mask);
      poc_scale_data(rdata,length,scale,offset,&rmask);
      mask=rmask;
      valcpy(data,rdata,length);
      delete[]rdata;
      }
    else{
      if(axes!="T")
        TRAP_ERR_EXIT(ENOEXEC,"poc_get_cvara(,,,string,indexes_t) not coded yet\n");
      
      const int frame=indexes[0];
      
      status=poc_get_cvara(filename,aname,gname,frame,data);
      
      mask=NC_FILL_COMPLEX;
      }
    
    return status;
    }
  
  bool isComplex(range_t<double> *range=0,double *maxFrac=0) const{
    int j;
    
    if(range!=0)
      range->init();
    
    if(maxFrac!=0)
      *maxFrac=0.;
    
    for(j=0;j<length;j++){
      complex<double> *z=&data[j];
      double
        re=abs(real(*z)),
        im=abs(imag(*z));
      
      if(range!=0)
        (*range)<<re;
      
      if(maxFrac!=0){
        const double
          a=fabs(re);
        const u_int64_t
          i=a;
        updatemax(maxFrac,a-i);
        }
      
      if(im>1e-3*re)/* is complex */
        return true;
      }
    
    return false;
    }
  
  };


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

class atlas_grid_or_mesh_t

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  {
public:
  poc_grid_data_t data;
  grid_t grid;
  mesh_t mesh;
  int discretisation;
  
  void destroy(){
    data.destroy();
    grid.free();
    mesh.destroy();
    }
  };

extern int poc_inq_dimvars(const string &filename,const poc_list_t<poc_dim_t> &dims,poc_data_t<double> *dimvars,int verbose=0);
extern int findVariableThasIs(const poc_global_t &global,const poc_list_t<poc_dim_t> &list,int nDimMax,const char (*prefixes)[16],int verbose=0,int mode=0);
extern int findAnyVariableThasIs(const poc_global_t &global,int ndim,const char *suffix,int verbose=0);
extern const poc_dim_t *findAxisDim(const poc_var_t &cvar,const poc_global_t &global,const char axis,const poc_var_t *var=NULL,int verbose=0);

extern int poc_get_lon_lat_time(const char *path,double *lon=0,double *lat=0,double *times=0,poc_list_t<poc_dim_t> *dimensions=0);

extern int poc_grid_data_to_grid(const string &path,const poc_global_t &global,const poc_var_t &var,const poc_grid_data_t &gd,grid_t *grid,int verbose,int frame);
extern int poc_grid_data_to_grid(const string &path,const poc_global_t &global,const poc_var_t &var,poc_grid_data_t *gd,grid_t *grid,int verbose,int frame);
extern int poc_grid_data_to_mesh(const string &path,const poc_global_t &global,const poc_var_t &var,const poc_grid_data_t &gd,mesh_t *mesh,int verbose,int frame);
extern int poc_grid_data_to_mesh(const string &path,const poc_global_t &global,const poc_var_t &var,poc_grid_data_t *gd,mesh_t *mesh,int verbose,int frame,double **depths=0);

extern int poc_get_grid_data(const string &path,const poc_var_t &var,poc_grid_data_t *gd,int verbose=0,poc_global_t *global=NULL);
extern int poc_get_axes(const string &path,poc_var_t *var,int verbose=0);

extern int poc_get_discretisation(const poc_var_t &var,const mesh_t *mesh=0,int verbose=0);
extern int test_grid_data(const poc_grid_data_t &gdata,bool *isGrid,bool *isMesh,int verbose=0);
extern int poc_get_grid_or_mesh(const string &path,const poc_var_t &var,grid_t *grid,mesh_t *mesh,int verbose,int frame,int *discretisation=0);
extern int poc_get_mesh(const string &path,const poc_var_t &var,mesh_t *mesh,int verbose,int frame,int *discretisation=0);

extern int poc_get_grid(const string &path,const poc_var_t &var,grid_t *grid,int verbose=0,int frame=-1);
extern int poc_get_grid(const string &path,const string &varname,grid_t *grid,int verbose=0,int frame=-1);

extern void poc_grid_dim(const grid_t& grid, poc_dim_t *nx, poc_dim_t *ny, const string & nameSuf="");
extern int poc_save_grid(const string& path, poc_var_t *gridVar, const grid_t& grid, int overwrite, int verbose, const string & loc="", float xdimoffset=0.f, float ydimoffset=0.f, float zdimoffset=0.f, const poc_list_t<poc_att_t> *attributes=0);
extern poc_var_t poc_save_grid(const string& path,const grid_t& grid,const char *srcFile=0,int srcLine=-1,const string & loc="",int verbose=1);

extern int poc_save_mesh(const string& path,const mesh_t & mesh,const char *srcFile=0,int srcLine=-1,int verbose=0,int discretisation=-1,bool withTime=false,poc_var_t *var=0);

extern int poc_put_mesh_vara(const string & path,const mesh_t & mesh,int discretisation,poc_var_t *var,int frame,const double *buffer,int verbose=1);
extern int poc_put_mesh_vara(const string & path,const mesh_t & mesh,int discretisation,poc_var_t *var,int frame,const float *buffer,int verbose=1);
extern int poc_put_mesh_vara(const string & path,const mesh_t & mesh,int discretisation,poc_var_t *var,int frame,const int *buffer,int verbose=1);

extern int poc_put_mesh_vara(const string & path,const mesh_t & mesh,int discretisation,const string & name,int frame,const double *buffer,const string & units="",int verbose=1);
extern int poc_put_mesh_vara(const string & path,const mesh_t & mesh,int discretisation,const string & name,int frame,const float *buffer,const string & units="",int verbose=1);
extern int poc_put_mesh_vara(const string & path,const mesh_t & mesh,int discretisation,const string & name,int frame,const int *buffer,const string & units="",int verbose=1);


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T> int poc_put_mesh_cvara(const char *path,const mesh_t & mesh,int discretisation,const char *aname,const char *gname,int frame,const complex<T> *buffer,const string & units="",int verbose=1)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// tidal loading with an unstructured grid
/*----------------------------------------------------------------------------*/
{
  int status;
  poc_var_t var,tv;
  
  if(verbose){STDOUT_BASE_LINE("saving to %s(%s,%s)... ",path,aname,gname);fflush(stdout);}
  
  const bool
    withTime=not updatemax(&frame,0);
  
  status=poc_save_mesh(path,mesh,__FILE__,__LINE__,verbose,discretisation,withTime,&var);
  
  var.init(aname,NC_DOUBLE,comodo_standard_name(),units);
  status=poc_def_var(path,var);
  if(status!=NC_NOERR) NC_TRAP_ERROR(return,status,verbose,"poc_def_var(\"%s\",(\""+var.name+"\"),...) error",path);
  
  var.init(gname,NC_DOUBLE,comodo_standard_name(),"degrees");
  status=poc_def_var(path,var);
  if(status!=NC_NOERR) NC_TRAP_ERROR(return,status,verbose,"poc_def_var(\"%s\",(\""+var.name+"\"),...) error",path);
  
  status=poc_put_cvara(path,aname,gname,frame,buffer,verbose);
  if(status!=NC_NOERR) NC_TRAP_ERROR(return,status,verbose,"poc_put_cvara(\"%s\",\"%s\",\"%s\",...) error",path,aname,gname);
  
  if(verbose) printf("%d done.\n",status);
  return status;
}


#define COMODO_URL "http://pycomodo.forge.imag.fr/norm.html"
extern int comodo_compliance(const string &path,const vector<string> &vars,string *report, int verbose=0);

#endif
