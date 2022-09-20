
/*******************************************************************************

  T-UGO tools, 2006-2013

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

\brief definitions of usefull lists
*/
/*----------------------------------------------------------------------------*/


#include <errno.h>
#include <string.h> // for strerror !
#include <stdlib.h> // for free

#include "poc-assertions.h"
#include "poc-array-stats.hpp"
#include "poc-list.hpp"


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int isfinite(const range_t<double> &r)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=isfinite(r.min) and isfinite(r.max);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

poc_entry_t::~poc_entry_t()

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
if(list!=0){
  list->~poc_entry_list_t();
  list=0;
  }
}


static int one(const struct dirent *unused){return 1;}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

static int is_visible(const struct dirent *entry)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  const bool visible=entry->d_name[0]!='.';
  return visible;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_list_dir(const string path,poc_entry_list_t *list,int onlyVisible,int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  /* see info:/libc/Accessing%20Directories */
  struct dirent **entries;
  int (*selector)(const struct dirent *entry);
  int n,i;
  
  if(onlyVisible)
    selector=is_visible;
  else
    selector=one;
  
  n=scandir(path.c_str(),&entries,selector,alphasort);
  if(n<0) TRAP_ERR_RETURN(errno,verbose,"scandir(\""+path+"\",) error (%d %s)\n",errno,strerror(errno));
  
  list->clear();
  
  for(i=0;i<n;i++){
    struct dirent *entry=entries[i];
    
    (*list)<<poc_entry_t(entry);
    }
  
  free(entries);
  return 0;
}
