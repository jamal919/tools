
/*******************************************************************************

  T-UGO tools, 2006-2014

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

\author  Florent Lyard      LEGOS/CNRS, Toulouse, France. florent.lyard@legos.obs-mip.fr
\author  Laurent Roblou     LEGOS/CNRS, Toulouse, France
\author  Damien Allain      LEGOS/CNRS, Toulouse, France
\author  Yves Soufflet      LEGOS/CNRS, Toulouse, France
\author  Clément Mayet      LEGOS, Toulouse, France (PhD)
\author  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
\author  Frédéric Dupont    Université de Laval à Québec, Canada

\brief definition of usefull string functions missing from standard libraries
*/
/*----------------------------------------------------------------------------*/

#include <stdio.h> //for FILE
#include <stdarg.h> //variadic arrays
#include <sys/stat.h> /* for stat() */

#include <netcdf.h> //for NC_NOERR
#include <string.h>
#include <stdlib.h> //for exit,malloc

#include <fstream>

#include "poc-assertions.h"

#include "functions.h"


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  string get_extension(const string & path)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// get file extension
/**
\param path
\return the string after the last '.' or "" if there is no '.'
\sa strrncasecmp() and strncmp(), that can also be used instead.
*/
/*----------------------------------------------------------------------------*/
{
  const char *dot=strrchr(path.c_str(),'.');
  
  if(dot==0)
    return "";
  
  return &dot[1];
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  string replace_extension(const string & path,const string & extension)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  replace extension of path (if any) by the new extension

  path : original path
  extension : new extension
  return the new path

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  size_t pos;
  string newPath;
  
  pos=path.find_last_of('.');
  
  if(pos==string::npos){
    /* no extension : just append the new one */
    newPath=path+".";
    }
  else{
    pos++;
    newPath=path.substr(0,pos);
    }
  
  newPath+=extension;
  
  return newPath;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int get_file_size(const string & path,size_t *size)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// get file size in bytes
/**
\param path
\param *size can be \c NULL to only test the existence of the file
\return 0 on success and \c errno value on error
*/
/*----------------------------------------------------------------------------*/
{
  int status;
  struct stat fileStat;
  
  status=stat(path.c_str(),&fileStat);
  if(status!=0)
    status=errno;
  else if(size!=0)
    *size=fileStat.st_size;
  
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 bool isLetter(char c)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// rather than isalpha, because it is locale-dependent
/* see also strchr(const char*,int) */
/*----------------------------------------------------------------------------*/
{
  bool result;
  
  result = c=='_'||('A'<=c && c<='Z')||('a'<=c && c<='z');
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  const char *strrchr0(const char *s, int c)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// gets first character, after the last of a character if found
/*----------------------------------------------------------------------------*/
{
  const char *r=strrchr(s,c);
  if(r==NULL)
    return s;
  return r+1;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int replace(string *s,const string & match,const string & replacement,int count)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// does what you think it does
/*----------------------------------------------------------------------------*/
{
  int i,pos=0;
  if(count<0)count=s->length();
  
  for(i=0;i<count;i++){
    pos=s->find(match,pos);
    if(pos>=s->length()) break;
    s->replace(pos,match.length(),replacement);
    pos+=replacement.length();
    }
  
  return i;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  string replace(string s,const string & match,const string & replacement,int count)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// does what you think it does
/*----------------------------------------------------------------------------*/
{
  replace(&s,match,replacement,count);
  return s;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int strncmp(const string &prefix, const string &s)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// strncmp overload
/*----------------------------------------------------------------------------*/
{
  int status;
  status=(prefix!=s.substr(0,prefix.length()));
  
  return(status);
  
//   /** With this declaration : \code strncmp(prefix,s) \endcode */
//   return /** is equivalent to
//     \code /**/ // COMPILED CODE BELOW !!!
//     prefix!=s.substr(0,prefix.length()) /** \endcode */
//     ;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int strrncmp(const string &s1,const string &s2)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// test suffix
/*----------------------------------------------------------------------------*/
{
  const size_t l1=s1.length();
  const size_t l2=s2.length();
  if(l2>=l1)
    return s1.compare(s2.substr(l2-l1));
  else
    return -s2.compare(s1.substr(l1-l2));;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int strncasecmp(const string & prefix, const string & str)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// strncasecmp overload
/*----------------------------------------------------------------------------*/
{
//   const size_t ls=str.length();
  const size_t lpref=prefix.length();
  int result;
  
  result=strncasecmp(prefix.c_str(),str.c_str(),lpref);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int strrncasecmp(const string & str, const string & suf)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  test case insensitive suffix
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  const size_t lstr=str.length();
  const size_t lsuf=suf.length();
  int result=-1;
  
  if(lstr>=lsuf)
    result=strncasecmp(&str.c_str()[lstr-lsuf],suf.c_str(),lsuf);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void tolower(string *str)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  for(int i=0;i<str->length();i++){
    char *c=&(*str)[i];
    *c=tolower(*c);
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  char *poc_strdup(const char *source)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  size_t l;
  char *s=NULL;

  if(source==0) return(0);
  
  l=strlen(source);
  s=new char[l+1];
  strcpy(s,source);
  s[l]='\0';
  
  return s;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  char *poc_fortran_strdup(const char *source)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  size_t l;
  char *s=NULL;

  if(source==NULL) return s;
  
  l=strchr(source,' ')-source;
  s=new char[l+1];
  strncpy(s,source,l);
  s[l]='\0';
  
  return s;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  string poc_fortran_string(const char *source)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  size_t l;
  string s;

  if(source==NULL) return s;
  
  l=strchr(source,' ')-source;
  s.assign(source,l);
  
  return s;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  char *check_string(const char* word)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/** converts a string in a certain way
\param *word the word to convert
\returns the converted string
*/
/*----------------------------------------------------------------------------*/
{
  char *string=NULL,*p=NULL;

  string= strdup(word);

/*------------------------------------------------------------------------------
  replaces all '_' by ' ' */
  do {
    p=strchr(string,'_');
    if(p!=NULL) {
      *p=' ';
      }
    } while (p!=NULL);
  
/*------------------------------------------------------------------------------
  then capitalises every first letter of every word */
  p=string;
  if((*p>96)&&(*p<123)) {
    *p-=32;
    }
  
  do {
    p=strchr(p,' ');
    if(p!=NULL) {
      p+=1;
      if((*p>96)&&(*p<123)) {
        *p-=32;
        p+=1;
        }
      }
    } while (p!=NULL);
  
  return (string);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

char *clean_string(const char *word)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  char *string=NULL,*p=NULL;

  string= strdup(word);

  do {
    p=strchr(string,'_');
    if(p!=NULL) {
      *p='-';
      }
    }
  while (p!=NULL);

  return (string);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int decode_format(const string & format, string *input, string *output)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/** \brief decodes the file format names given by the user in arguments list (format="..."") and returns two strings corresponding to the input and output format names, as defined in the main function.

@param format the format string as read from char *argv[]
@param input  a pointer to the input file format name string
@param output a pointer to the output file format name string
@return zero
*/
/*----------------------------------------------------------------------------*/
{
  int nitems;
  char buf[32];
  string s,dum;
  size_t pos;

  s=format;

//rootname.erase(pos,4);

  if(s.c_str()==0) return(-1);

  pos=s.find("--format=");
  s.erase(pos,9);

  pos=s.find("input=");
  if(pos!=string::npos) {
    dum=s.substr(pos+6,string::npos);
    nitems=sscanf(dum.c_str(),"%s",buf);
    *input=string(buf);
    }

  pos=s.find("output=");
  if(pos!=string::npos) {
    dum=s.substr(pos+7,string::npos);
    nitems=sscanf(dum.c_str(),"%s",buf);
    *output=string(buf);
    }

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  string isFormat_help(const string & indent,const string & inputDefault,const string & outputDefault)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  return
    indent+"--format=... deprecated\n"+
    indent+"-iF  followed by input format"+inputDefault+"\n"+
    indent+"-oF  followed by output format"+outputDefault+"\n";
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  bool isFormat(char **argv, int *n, string *input, string *output)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/** \brief decodes the file format names given by the user in arguments list (format="..."") and returns 1 or 2 strings corresponding to the input and output format names, as defined in the main function.

\param **argv
\param *n argument index
\param *input input format, modified only if necessary
\param *output output format, modified only if necessary
\return zero
*/
/*----------------------------------------------------------------------------*/
{
  const string argvn=argv[*n];
  
  if( strncmp("--format=",argvn)==0 ){
    decode_format(argvn,input,output);
    (*n)++;
    
    fprintf(stderr,"WARNING: YOU SHOULD REPLACE\n"
      "  "+argvn+"\n"
      "BY:\n"
      "  -iF "+*input+" -oF "+*output+"\n");
    
    return true;
    }
  
  if( strncmp("--format_in=",argvn)==0 ){
    *input=&argvn[12];
    (*n)++;
    
    fprintf(stderr,"WARNING: YOU SHOULD REPLACE\n"
      "  "+argvn+"\n"
      "BY:\n"
      "  -iF "+*input+"\n");
    
    return true;
    }
  
  if( strncmp("--format_out=",argvn)==0 ){
    *output=&argvn[12];
    (*n)++;
    
    fprintf(stderr,"WARNING: YOU SHOULD REPLACE\n"
      "  "+argvn+"\n"
      "BY:\n"
      "  -oF "+*input+"\n");
    
    return true;
    }
  
  if( argvn=="-iF" ){
    (*n)++;
    *input=argv[*n];
    (*n)++;
    return true;
    }
  
  if( argvn=="-oF" ){
    (*n)++;
    *output=argv[*n];
    (*n)++;
    return true;
    }
  
  return false;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  vector<string> string_split(const  string  & theString, const  string  & theDelimiter)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
//    UASSERT( theDelimiter.size(), >, 0); // My own ASSERT macro.
  vector<string> theStringVector;
  size_t  start = 0, end = 0;

  while ( end != string::npos) {
    end = theString.find( theDelimiter, start);

    // If at end, use length=maxLength.  Else use length=end-start.
    theStringVector.push_back( theString.substr( start, (end == string::npos) ? string::npos : end - start));

    // If at end, use start=maxSize.  Else use start=end+delimiter.
    start = (   ( end > (string::npos - theDelimiter.size()) ) ?  string::npos  :  end + theDelimiter.size());
    }
  return(theStringVector);
}

// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   vector<string> string_split2(const string & line, const  string  & separator)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   int nitems;
//   size_t pointer, next;
//   vector<string> array;
//   string substring;
//   string tmp;
//   
//   substring=line;
//   
// /*------------------------------------------------------------------------------
//   remove initial occurence of delimiter */
//   pointer = substring.find(separator);
//   while(pointer==0) {
//     tmp=substring.substr(pointer+separator.length());
//     substring=tmp;
//     pointer=substring.find(separator);
//     }
//     
//     
//   while(true) {
//     pointer = substring.find(separator);
// /*------------------------------------------------------------------------------
//     remove repetitive occurence of delimiter */
//     while(pointer==0) {
//       tmp=substring.substr(pointer+separator.length());
//       substring=tmp;
//       pointer=substring.find(separator);
//       }
// /*------------------------------------------------------------------------------
//     at this point, substring first character is valid (i.e. NOT a delimiter) */
//     next=substring.find(separator);
//     if(next == string::npos) {
//       string *token=new string;
//       *token=substring.substr(0);
//       array.push_back(*token);
//       }
//     else {
//       string *token=new string;
//       *token=substring.substr(0, next);
//       array.push_back(*token);
//       }
//     if(next == string::npos)  break;
//     tmp=substring.substr(next);
//     substring=tmp;
//     }
// 
//   return (array);
// }
// 

// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   vector<string> string_split(const  string  & theString, const  vector<string>  & theDelimiters)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
// //    UASSERT( theDelimiter.size(), >, 0); // My own ASSERT macro.
//   vector<string> theStringVector;
//   size_t  start = 0, end = 0;
// 
//   while ( end != string::npos) {
//     end = theString.find( theDelimiter, start);
// 
//     // If at end, use length=maxLength.  Else use length=end-start.
//     theStringVector.push_back( theString.substr( start, (end == string::npos) ? string::npos : end - start));
// 
//     // If at end, use start=maxSize.  Else use start=end+delimiter.
//     start = (   ( end > (string::npos - theDelimiter.size()) ) ?  string::npos  :  end + theDelimiter.size());
//     }
//   return(theStringVector);
// }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int count_token(const char *s, const char *delim)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/** Counts the space-separated words in a string.
\date reviewed 21 Jul 2011
\author  Damien Allain

\param *s list of words
\returns the number of words
*/
/*----------------------------------------------------------------------------*/
{
  int count=0;
  char *token;
  char *dummy,*bkp;

  if(s==0) return(0);

  if(strlen(s)==0) return(0);

  dummy=strdup(s);
  bkp=dummy;

  // man:/strtok(3) says ``the tokens returned by strtok() are always nonempty strings''
  token = strtok(dummy,delim);
  while (token!=0) {
    count++;
    token = strtok(0,delim);
    }

  free(bkp);
  return(count);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int get_token(const char *s, char **params, int nparams)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/** Separates words from a string into an array of strings.
\date reviewed 21 Jul 2011
\author  Damien Allain

\param *s list of words
\param **params array of strings to fill in
\param nparams maximum number of words
\returns the number of words
*/
/*----------------------------------------------------------------------------*/
{
  int count;
  char *token;
  char *dummy,*bkp;

  if(s==NULL || strlen(s)==0)
    return(0);

  dummy=strdup(s);
  bkp=dummy;

  token = strtok(dummy," ");
  for(count=0;count<nparams;count++){
    if(token==NULL)break;
    params[count]=strdup(token);
    token = strtok(0," ");
    }

  free(bkp);
  return(count);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int get_token(const char *s, char ***params)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/** Separates words from a string into an array of strings.
\date reviewed 21 Jul 2011
\author  Damien Allain

\param *s list of words
\param ***params array of strings to fill in
\returns the number of words
*/
/*----------------------------------------------------------------------------*/
{
  int n,i;
  char *token;
  char *dummy,*bkp;

  n=count_token(s);
  if(n==0){
    return(0);
    }
  *params=new char*[n];

  dummy=strdup(s);
  bkp=dummy;

  token = strtok(dummy," ");
  for(i=0;i<n;i++){
    if(token==NULL)break;
    (*params)[i]=poc_strdup(token);
    token = strtok(0," ");
    }
  if(i<n)
    (*params)[i]=NULL;

  free(bkp);
  return n;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  char *get_nth_token(const char *s, int parami)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/** Get one particular word from a string.
\param *s string
\param parami word index
\returns the word or NULL if error
*/
/*----------------------------------------------------------------------------*/
{
  if(parami<0)
    return NULL;
  char **params;
  int n=get_token(s,&params);
  char *param;
  if(parami>=n){
    param=NULL;
    parami=-1;
    }
  else
    param=params[parami];
  
  for(int i=0;i<n;i++){
    if(i==parami)
      continue;
    //else
    delete[]params[i];
    }
  delete[]params;
  
  return param;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  vector<double> get_tokens(string & line, const char* separator)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int nitems;
  size_t pointer;
  vector<double> array;
  double q;
  string substring;
  
  pointer = line.find(separator);
  while(pointer ==0) {
    substring=line.substr(pointer+strlen(separator));
    line.assign(substring);
    pointer = line.find(separator);
    }
    
  nitems = 1;
  
  while(nitems == 1) {
    nitems  = sscanf(line.c_str(), "%lf", &q);
    if(nitems!=0) array.push_back(q);
    pointer = line.find(separator);
    while(pointer ==0) {
      substring=line.substr(pointer+strlen(separator));
      line.assign(substring);
      pointer = line.find(separator);
      }
    if(pointer != string::npos) {
      substring=line.substr(pointer+strlen(separator));
      line.assign(substring);
      }
    else break;
    }

  return (array);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  regex_t poc_regcomp_(const char *regex, int cflags,const char *file,int line)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// wrapper for regcomp(). See man:/regcomp
/**
Do no call directly. Call #poc_regcomp instead.
\param file replace this with __FILE__
\param line replace this with __LINE__
\return the pattern buffer storage area \c preg
*/
/*----------------------------------------------------------------------------*/
{
  regex_t preg;
  int status;
  
  status=regcomp(&preg,regex,cflags);
  
  if(status!=0){
    const int l=1024;
    char msg[l];
    
    regerror(status,&preg,msg,l);
    
    fprintf(stderr,"%s:%d:regex coding error (%d %s)\n",file,line,status,msg);
    wexit(ENOEXEC);
    }
  
  return preg;
}

