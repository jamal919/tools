#if !defined(NETCDF_CLASSES)
#define NETCDF_CLASSES

#include <cfloat>               //for DBL_MAX
#include "functions.h"
#include "poc-netcdf-assertions.h"

#include "map-classes.h"

/*------------------------------------------------------------------------
cdf dimension type*/
class cdfdim_t { // cdfdim_t has moved to class
  private:
  public:
    int     id;
    size_t  length;
    bool isunlimited;
    char    *name;

  // Constructors
  private:
    void init(){
      id=NC_EBADDIM;/* Invalid dimension id or name */
      length = 0;
      isunlimited=false;
      name = NULL;
      }
    void init(const cdfdim_t &src){
      id=src.id;
      length=src.length;
      isunlimited=src.isunlimited;
      name = new char[strlen(src.name)+1];
      strcpy(name, src.name);
      }
  public:
    cdfdim_t(){
      init();
      }
    cdfdim_t(const char* nname,const size_t nlength){
      id=NC_EBADDIM;/* Invalid dimension id or name */
      length = nlength;
      isunlimited=(nlength==NC_UNLIMITED);
      name = poc_strdup(nname);
      }
    cdfdim_t(const cdfdim_t &src) {
      init(src);
      }
    cdfdim_t &operator = (const cdfdim_t &src) {
      destroy();
      init(src);
      }
    friend int operator == (const cdfdim_t &a,const cdfdim_t &b) {
      return (a.isunlimited || b.isunlimited || a.length==b.length) && a.name!=NULL && b.name!=NULL && !strcmp(a.name,b.name);
      }
    void destroy(){
      if(name!=NULL){delete[] name;name=NULL;}
      }

} ;

/*------------------------------------------------------------------------
cdf attribute type*/
class cdfatt_t
{
  private:

  public:
    int     id;
    nc_type type;
    size_t  length;
    char    *name;
    char    *data;

  // Constructors
  cdfatt_t(){
    id=NC_ENOTATT; /* Attribute not found */
    type = NC_NAT; /* NAT = 'Not A Type' (c.f. NaN) */
    length=0;
    name=NULL;
    data=NULL;
    }
  private:
  void init_name_and_id(const char *n){
    id=NC_ENOTATT; /* Attribute not found */
    name =poc_strdup(n);
    }
  public:
  void init(const char *n,const char *d){
    init_name_and_id(n);
    type = NC_CHAR; /* NAT = 'Not A Type' (c.f. NaN) */
    length=strlen(d);
    data =poc_strdup(d);
    }
  cdfatt_t(const char *n,const char *d){
    init(n,d);
    }
  void init(const char *n,const size_t dn,const double *d){
    init_name_and_id(n);
    type = NC_DOUBLE;
    length=dn;
    size_t size=length*sizeof(double);
    data = new char[size];
    memcpy(data,d,size);
    }
  void init(const char *n,const size_t dn,const float *d){
    init_name_and_id(n);
    type = NC_FLOAT;
    length=dn;
    size_t size=length*sizeof(float);
    data = new char[size];
    memcpy(data,d,size);
    }
  template<typename T> void init(const char *n,const T d){
    init(n,1,&d);
    }
  template<typename T> cdfatt_t(const char *n,const T d){
    init(n,d);
    }

  cdfatt_t &operator = (const cdfatt_t &src) {
    id=src.id;
    type=src.type;
    length=src.length;

    size_t size = 0;
    switch (type) {
      case NC_CHAR:
        size=length+1;
        break;

      case NC_BYTE:
        size=length;
        break;

      case NC_SHORT:
        size=length*sizeof(short);
        break;

      case NC_INT:
        size=length*sizeof(int);
        break;

      case NC_FLOAT:
        size=length*sizeof(float);
        break;

      case NC_DOUBLE:
        size=length*sizeof(double);
        break;
      }
    if(name != 0) delete [] name;
    name = new char[strlen(src.name)+1];
    strcpy(name, src.name);
    if(data!=0) {
      delete[] data;
      }
    data =new char[size];
    for(size_t n=0;n<length;n++) data[n]=src.data[n];
    if(type==NC_CHAR) data[length]=0;
    return(*this);
    }

    void reset(){
      id = -1;
      type = NC_NAT;
      length = 0;
      destroy();
      }// end reset()

    void destroy(){
      deletep(&data);
      deletep(&name);
      }
    ~cdfatt_t(){//called by functions that take the class as argument (EVEN const-QUALIFIED!!!)
      //destroy();
      }
};

class cdfvar_t
{
  private:

  public:
    int     id;
    nc_type type;
    int     ndim;
    cdfdim_t *dim;
    int     natt;
    cdfatt_t *att;
    char    *name;

    // Constructors
  cdfvar_t(){
    init();
    }
  void init(){
    id=NC_ENOTVAR;
    type = NC_NAT;
    ndim=0;
    dim=NULL;
    natt=0;
    att =NULL;
    name=NULL;
    }
  void init(const char *n,const nc_type t){
    id=NC_ENOTVAR;
    ndim=0;
    dim=NULL;
    natt=0;
    att =NULL;
    type = t;
    name=poc_strdup(n);
    }
  cdfvar_t(const char *n,const nc_type t){
    init(n,t);
    }

  void initatt(int nnatt=-1){
  ///Initialises natt and (re-)allocates att accordingly
  /**
  \param nnatt new value for natt. Default: natt itself.
  */
    if(att!=NULL) {
      for(size_t k = 0; k < natt; k++){
        att[k].destroy();
        }
      delete[] att;
      }
    if(nnatt>=0)
      natt=nnatt;
    if(natt==0){
      att=NULL;
      return;
      }
    if((att=new cdfatt_t[natt])==NULL){
      __ERR_BASE_LINE__("\n%s:%d: null pointer allocation. Exiting with an error code.\n",__FILE__,__LINE__);
      exit(-1);
      }
    }

  void add_att(cdfatt_t &att_){
    if(att==NULL)
      natt=0;
    cdfatt_t *new_att=new cdfatt_t[natt+1];
    if(att!=NULL){
      for(int i=0;i<natt;i++)
        new_att[i]=att[i];
      delete[]att;
      }
    att=new_att;
    att[natt]=att_;
    att[natt].id=natt;
    natt++;
    }
  void add_att(const char *n,char *d){
    add_att(n,(const char *)d);
    }
  void add_att(const char *n,const char *d){
    cdfatt_t att_(n,d);
    add_att(att_);
    }
  template<typename T> void add_att(const char *n,const T d){
    cdfatt_t att_(n,d);
    add_att(att_);
    }

  void initdim(int nndim=-1){
  ///Initialises ndim and (re-)allocates dim accordingly
  /**
  \param nndim new value for ndim. Default: ndim itself.
  */
    if(dim!=NULL) {
      for(size_t k = 0; k < ndim; k++){
        dim[k].destroy();
        }
      delete[] dim;
      }
    if(nndim>=0)
      ndim=nndim;
    if(ndim==0){
      dim=NULL;
      return;
      }
    if((dim=new cdfdim_t[ndim])==NULL){
      __ERR_BASE_LINE__("\n%s:%d: null pointer allocation. Exiting with an error code.\n",__FILE__,__LINE__);
      exit(-1);
      }
    }

  size_t size(){
    size_t s=1;
    if(dim!=NULL) {
      for(size_t k = 0; k < ndim; k++){
        s*=dim[k].length;
        }
      }
    return(s);
    }

  void init(const cdfvar_t &src) {
    id=src.id;
    type=src.type;
    
    ///This reinitialises dim and its elements
    initdim(src.ndim);
    for(size_t n=0;n<ndim;n++){
      dim[n]=src.dim[n];
      }
    
    initatt(src.natt);
    for(size_t n=0;n<natt;n++) att[n]=src.att[n];
    
    if(src.name!=NULL) {
      if(name!=NULL)delete[]name;
      name=poc_strdup(src.name);
      }
    }

  cdfvar_t &operator = (const cdfvar_t &src) {
    init(src);
    return(*this);
    }

  cdfvar_t & operator << (cdfvar_t src) {
    swapValues(&id,&src.id);
    swapValues(&type,&src.type);
    swapValues(&ndim,&src.ndim);
    swapValues(&dim,&src.dim);
    swapValues(&natt,&src.natt);
    swapValues(&att,&src.att);
    swapValues(&name,&src.name);
    return(*this);
    }

  friend int operator == (const cdfvar_t &a,const cdfvar_t &b){
    int status=a.ndim==b.ndim && a.name!=NULL && b.name!=NULL && !strcmp(a.name,b.name);
    for(int i=0;status && i<a.ndim;i++){
      status=a.dim[i]==b.dim[i];
      }
    return status;
    }

  friend int operator != (const cdfvar_t &a,const cdfvar_t &b){
    return !(a==b);
    }

  int findattr(const char *name) {
    int n;
    for(n=0;n<this->natt;n++) {
      if(strcmp(this->att[n].name,name)==0){
        return n;
        }
      }
    return -1;
    }

  char *getattr(const char *name) {
    int n=findattr(name);
    if(0<=n && n<=natt)
      return this->att[n].data;
    return NULL;
    }

  void destroy_dim() {
    if(dim!=NULL){
      for(size_t k = 0; k < ndim; k++){
        dim[k].destroy();
        }
      delete[]dim;dim=NULL;
      }
    }

  void destroy() {
    deletep(&name);
    destroy_dim();
    if(att!=NULL){
      for(size_t k = 0; k < natt; k++){
        att[k].destroy();
        }
      delete[]att;
      att=NULL;
      }
    }

//   ~cdfvar_t() {//called by functions that take the class as non-reference argument (EVEN const-QUALIFIED!!!)
//     //destroy();
//     }

};

/*------------------------------------------------------------------------
cdf global type*/
class cdfgbl_t{
private:
public:
  int     nvarsp,ndimsp,ngattsp;
  cdfdim_t *dimension;
  cdfvar_t *variable;
  cdfatt_t *attribute;
  int  unlimdimid;
  char *production;
  // Constructors
  cdfgbl_t(){
    nvarsp = ndimsp = ngattsp = 0;
    dimension  = NULL;
    variable   = NULL;
    attribute  = NULL;
    unlimdimid = NC_EBADDIM;
    }

  void initattribute(int nngattsp=-1){
  ///Initialises ngattsp and (re-)allocates attribute accordingly
  /**
  \param nngattsp new value for ngattsp. Default: ngattsp itself.
  */
    if(nngattsp>0)
      ngattsp=nngattsp;
    if(attribute!=NULL) {
      delete[] attribute;
      }
    if((attribute=new cdfatt_t[ngattsp])==NULL){
      __ERR_BASE_LINE__("\n%s:%d: null pointer allocation. Exiting with an error code.\n",__FILE__,__LINE__);
      exit(-1);
      }
    }

  cdfgbl_t(const cdfvar_t &info){
    add_var(info);
    }
  int add_var(const cdfvar_t &info){
    ///\returns the id of the variable or -1 if error.
    int i;//variable index
    for(i=0;i<nvarsp;i++){
      if(!strcmp(info.name,variable[i].name)){
        if(info!=variable[i])
          return -1;
        return i;//variable already exist
        }
      }
    ///It adds the variable
    cdfvar_t *new_variable=new cdfvar_t[nvarsp+1];
    for(i=0;i<nvarsp;i++){
      new_variable[i]=variable[i];
      }
    new_variable[i].id=i;
    new_variable[i]=info;
    ///and the dimensions.
    for(int j=0;j<info.ndim;j++)
      new_variable[i].dim[j].id=add_dim(info.dim[j]);
    delete[]variable;
    variable=new_variable;
    nvarsp++;
    return i;
    }
  int add_dim(const cdfdim_t &dim){
    ///\returns the id of the dimension
    int i;
    for(i=0;i<ndimsp;i++){
      if(dimension[i]==dim)return i;
      }
    cdfdim_t *new_dimension=new cdfdim_t[ndimsp+1];
    for(i=0;i<ndimsp;i++){
      new_dimension[i]=dimension[i];
      }
    new_dimension[i]=dim;
    new_dimension[i].id=i;
    if(dim.isunlimited)unlimdimid=i;
    delete[]dimension;
    dimension=new_dimension;
    ndimsp++;
    return i;
    }

//    ~cdfgbl_t()
    void destroy() {

    if(variable != 0){
      for(size_t k = 0; k < nvarsp; k++){
        this->variable[k].destroy();
        }
      delete [] variable; variable = NULL;
      }

    if(attribute != 0){
      for(size_t k = 0; k < ngattsp; k++){
        this->attribute[k].destroy();
        }
      delete [] attribute; attribute = NULL;
      }

    if(dimension != 0){
      for(size_t k = 0; k < ndimsp; k++){
        this->dimension[k].destroy();
        }
      delete [] dimension; dimension = NULL;
      }

    /// uncoment here when possible

//     this->dimension.destroy();
//     this->variable.destroy();
//     this->attribute.destroy();
    nvarsp = ndimsp = ngattsp = unlimdimid = -1;
    }
};


extern void delete_cdfvar_p(cdfvar_t **var);


/*------------------------------------------------------------------------
poc grid type*/

class pocgrd_t{
  private:
  public:
  int      id;
  char     *name;
  cdfvar_t *lon,*lat,*x,*y,*z,*t;
  char     *filename;

    // Constructors
  pocgrd_t(){
    id = -1;
    name = 0; filename = 0;
    lon = 0; lat = 0;
    x = 0; y = 0; z = 0; t = 0;
    }
  
  void destroy(){
    id = -1;
    deletep(&name);
    delete_cdfvar_p(&lon);
    delete_cdfvar_p(&lat);
    delete_cdfvar_p(&x);
    delete_cdfvar_p(&y);
    delete_cdfvar_p(&z);
    delete_cdfvar_p(&t);
    deletep(&filename);
    }

    // Destructors
/*    ~pocgrd_t(){
      destroy();
    }*/
};

/*------------------------------------------------------------------------
decoding type*/
class decoded_t{
  private:
  public:
  cdfvar_t  info;
  char      *axis;
  char      *associate;
  int       xdim,ydim,zdim,tdim,fdim;
  size_t    xlen,ylen,zlen,tlen,flen;
  int       vx,vy,vz,vt,vf;
  char      *units;
  grid_t    *grid;
  double    *time;
  float     *buffer,offset,scale;
  float spec;///< mask value, from "missing_value" or "_FillValue" attribute. See poc_decode_mask()
//  char      *buffer,*offset,*scale,*spec;
  float valid_min,valid_max;;
  char *standard_name, *long_name, *short_name;
  char  *production;

  decoded_t() {
    axis=0;
    associate=0;
    xdim=ydim=zdim=tdim=fdim=-1;
    xlen=ylen=zlen=tlen=flen=0;
    vx=vy=vz=vt=vf=-1;
    units=0;
    grid=0;
    time=0;
    buffer=0;
    standard_name=0; long_name=0; short_name=0;
    production=0;
    }

  void destroy() {
    info.destroy();
    deletep(&axis);
    deletep(&associate);
    deletep(&units);
    deletep(&time);
    deletep(&buffer);
    deletep(&standard_name);
    deletep(&long_name);
    deletep(&short_name);
    xdim=ydim=zdim=tdim=fdim=-1;
    xlen=ylen=zlen=tlen=flen=-1;
    vx=vy=vz=vt=vf=-1;
    deletep(&production);
    }
};


class UGdecoded_t{

  private:
  public:
  cdfvar_t  info;
  char      *axis;
  char      *associate;
  int       hdim,vdim,tdim,fdim;
  int       hdiscretisation,vdiscretisation;
  size_t    hlen,vlen,tlen,flen;
  int       vx,vy,vz,vt,vmask;
  char      *units;
  grid_t    *grid;
  double    *time;
  float     *buffer,offset,scale,spec;
  
  UGdecoded_t () {
    hdim=vdim=tdim=fdim=-1;
    hlen=vlen=tlen=flen=-1;
  }
} ;

/*------------------------------------------------------------------------
variable type*/
// typedef struct {
//   int     id;
//   int     ndim;
//   int     *dim,*axis;
//   char    **dimname,**axisname;
//   int     natt;
//   char    **attname;
//   double  time,lon,lat;
//   int     layer,frame;
//   char    *name,*standard_name,*long_name;
//   char    *units;
//   char    *origin;
//   grid_t  *grid;
//   float   *buffer;
// } variable_t;


class variable_t{ // variable has moved to class

  private:
  public:
    int     id;
    int     ndim;
    int     *dim,*axis;
    char    **dimname,**axisname;
    int     natt;
    char    **attname;
    double  time,lon,lat;
    int     layer,frame;
    char    *name,*standard_name,*long_name;
    char    *units;
    char    *origin,*production;
    grid_t  *grid;
    float   *buffer;

  variable_t(){
    id = ndim = natt = -1;
    time = lon = lat = DBL_MAX;
    dim = 0;
    axis = 0;
    dimname = 0;
    axisname = 0;
    attname = 0;
    name = 0;
    standard_name = 0;
    long_name = 0;
    units = 0;
    origin = production = 0;
    grid = 0;
    buffer = 0;
  }

  void reset(){ // free memory blocks
    id = ndim = natt = -1;
    time = lon = lat = DBL_MAX;
    // free *dim
  if(dim != 0){
    free(dim);
    dim = 0;
    }
  //free *axis
  if(axis != 0){
    free(axis);
    axis = 0;
    }
  // free **dimname
  for(int i = 0; i< ndim; i++){// ndim not sure
    if(dimname[i] != 0){
      free(dimname[i]);
      }
    }
  if(dimname != 0){
    free(dimname);
    dimname = 0;
    }
  // free **axisname
  for(int i = 0; i< ndim; i++){ // ndim not sure
    if(axisname[i] != 0){
      free(axisname[i]);
      }
    }
  if(axisname != 0){
    free(axisname);
    axisname = 0;
    }
  // free **attname
  for(int i = 0; i< natt; i++){ // natt not sure
    if(attname[i] != 0){
      free(attname[i]);
      }
    }
  if(attname != 0){
    free(attname);
    attname = 0;
    }
  // free *name
  if(name != 0){
    free(name);
    name = 0;
    }
  // free *standard_name
  if(standard_name != 0){
    free(standard_name);
    standard_name = 0;
    }
  // free *long_name
  if(long_name != 0){
    free(long_name);
    long_name = 0;
    }
  // free *units
  if(units != 0){
    free(units);
    units = 0;
    }
  // free *origin
  if(origin != 0){
    free(origin);
    origin = 0;
    }
  // free *grid (grid_t *)
/*  for(int i = 0; i< ???; i++){
  if(){
  grid[i].free();
}
}*/
  // *buffer
  if(origin != 0){
    free(origin);
    origin= 0;
    }
  }
};

#endif
