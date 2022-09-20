
/*******************************************************************************

  T-UGO tools, 2006-2017

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

\brief declarations of usefull lists
*/
/*----------------------------------------------------------------------------*/

#if POC_LIST_HPP == 0
#define POC_LIST_HPP 1

#include <dirent.h> //for DT_UNKNOWN and scandir

#include <vector>
#include <string>

using namespace std;


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<class T> class poc_deque_t : public vector<T> {

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// adds wrappers to erase() and insert() to vector
/**
list allows constant time insert and erase operations anywhere within the sequence
but deque and vector have an operator[]
and deque messes-up debugger
so using vector.

\note for printing as arrays in KDevelop, replace the 2 XXXX by the object name in the code below
\code
XXXX.v()[0]@XXXX.size()
\endcode
*/
/*----------------------------------------------------------------------------*/
public:
  
  /* THIS MUST BE CALLED SOMEWHERE TO BE ENABLED IN DEBUGGER */
  const T * v () const __attribute__((noinline)) {
    return &(*this)[0];
    }
  
  poc_deque_t(){
    /* SEE NOTE ABOVE */
    v();
    }
  
  poc_deque_t(const vector<T> & src):
  vector<T>(src){
    }
  
  poc_deque_t(const size_t n):
  vector<T>(n){
    }
  
  poc_deque_t(const T &src){
    (*this)=poc_deque_t<T>();
    this->set(src);
    }
  
  void clear(int n=0) {
    vector<T>::clear();
    if(n>0)
      vector<T>::reserve(n);
    }
  
  poc_deque_t<T>& erase(const int i){
    if(i>=0){
      vector<T>::erase(this->begin()+i);
      }
    return *this;
    }
  
  poc_deque_t<T>& insert(const int i, const T &x){
    if(i>=0){
      vector<T>::insert(this->begin()+i,x);
      }
    return *this;
    }
  
  poc_deque_t<T>& operator<<(const T &x){
    this->push_back(x);
    return *this;
    }
  
  int findI(const T *x) const {
    for(int i=0;i<this->size();i++){
      if(x==&(*this)[i])
        return i;
      }
    return -1;
    }
  
  int find(const T &x) const {
    for(int i=0;i<this->size();i++){
      if(x==(*this)[i])
        return i;
      }
    return -1;
    }
  
  };


typedef poc_deque_t<int> indexes_t;


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

class poc_name_id_t{

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
public:
  string name;
  int    id;

protected:
  
template<class T> inline void init_name_and_id(const T &src){
    name=src.name;
    id=src.id;
    }
  
  /* NOT CALLING init_name_and_id() `init()` FOR DESCENDANTS */
  void init_name_and_id(const string & name0,int id0=-1){
    name=name0;
    id=id0;
    }
  
public:
  
  poc_name_id_t(const string & name0="",int id0=-1){
    init_name_and_id(name0,id0);
    }
  };


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<class T> class poc_list_t : public poc_deque_t<T> {

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
public:
  
  poc_list_t(){
    }
  
  poc_list_t(const T &src){
    (*this)=poc_list_t<T>();
    this->set(src);
    }
  
  int find(const char *searched, T *dest=NULL) const{
    if(searched==NULL)
      return -1;
    string s=searched;
    return find(s,dest);
    }
  
  int find(const string &searched, T *dest=NULL) const{
    int i;
    for(i=0;i<this->size();i++){
      if((*this)[i].name==searched){
        if(dest!=NULL)
          *dest=(*this)[i];
        return i;
        }
      }
    return -1;
    }
  
  const T* findP(const string &searched) const{
    int i=find(searched);
    if(i<0)return NULL;
    return &(*this)[i];
    }
  
  T* findP(const string &searched){
    int i=find(searched);
    if(i<0)return NULL;
    return &(*this)[i];
    }
  
  const T* findP(const char *searched) const{
    if(searched==NULL)
      return NULL;
    string s=searched;
    return findP(s);
    }
  
  T* findP(const char *searched){
    if(searched==NULL)
      return NULL;
    string s=searched;
    return findP(s);
    }
  
  const T* findP(const char *s1,const char *s2) const{///\todo variable number of arguments
    const T *r;
    r=findP(s1);
    if(r)return r;
    r=findP(s2);
    return r;
    }
  
  int set(const T &src){
    int i=find(src.name);
    if(i<0){
      i=this->size();
      this->push_back(src);
      }
    else{
      (*this)[i]=src;
      }
    (*this)[i].id=i;
    return i;
    }
  
  poc_list_t<T>& erase(const int i){
    poc_deque_t<T>::erase(i);
    return *this;
    }
  
  poc_list_t<T>& erase(const string &searched){
    poc_deque_t<T>::erase(find(searched));
    return *this;
    }
  
  poc_list_t<T>& operator<<(const T &src){
    this->set(src);
    return *this;
    }
  
  poc_list_t<T>& operator<<(const poc_list_t<T> &src){
    int i;
    for(i=0;i<src.size();i++)
      this->set(src[i]);
    return *this;
    }
  };


class poc_entry_list_t;

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

class poc_entry_t : public poc_name_id_t {

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
public:
  unsigned char type;
  poc_entry_list_t *list;

private:
  void nullPointers(){
    type=DT_UNKNOWN;
    list=0;
    }

public:
  poc_entry_t(){
    nullPointers();
    }

  poc_entry_t(const struct dirent *entry){
    nullPointers();
    type=entry->d_type;
    name=entry->d_name;
    }

  ~poc_entry_t();

};


class poc_entry_list_t : public poc_list_t<poc_entry_t> {
  };

extern int poc_list_dir(const string path,poc_entry_list_t *list,int onlyVisible=1,int verbose=0);

#endif
