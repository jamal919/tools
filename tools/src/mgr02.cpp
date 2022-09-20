

/*******************************************************************************

  T-UGO tools, 2006-2009

  Unstructured Ocean Grid initiative

Contributors:

  Florent Lyard      LEGOS/CNRS, Toulouse, France
  Yoann Le Bars      LEGOS, Toulouse, France (PhD)
  Laurent Roblou     LEGOS/CNRS, Toulouse, France
  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
  Frédéric Dupont    Université de Laval à Québec, Canada

E-mail: florent.lyard@legos.obs-mip.fr

*******************************************************************************/

#include <config.h>

#include <list>
#include <vector>
#include <iomanip>
#include <fstream>

#include "tools-structures.h"
#include "netcdf-proto.h"
#include "mgr.h"
#include "tides.h"
#include "tides.def"
#include "list.h"
#include "zapper.h"
#include "polygones.h"


inline bool wcompare(tidal_wave first, tidal_wave second){
  return(first.omega < second.omega);
  }

inline bool compare_nocase (tidal_wave first, tidal_wave second){
  return -strcasecmp(first.name,second.name);
  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int mgr_compact(vector<mgr_t> & mgr, int& nmgr)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/* *----------------------------------------------------------------------------
  check empty mgr_t[] */
  
  for(size_t m = 0; m < mgr.size(); ++m){
    if(mgr[m].nwave == 0){
      mgr.erase(mgr.begin()+m);
      m--;
      }
    }
  
  nmgr=mgr.size();
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int *order_ALPHA(mgr_t **mgr,int nmgr, char *criterion, int direction)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int *list = 0;
  int *position = 0;
  char *s1 = 0 ,* s2 = 0;

  list    = new int [nmgr];
  position= new int [nmgr];

  for(size_t n = 0; n < nmgr; n++) {
    position[n]=0;
    }

  for(size_t n = 0; n < nmgr; n++) {
    s1 = check_string(mgr[n]->name);
    for(size_t m = n+1; m < nmgr; m++) {
      s2 = check_string(mgr[m]->name);
      if(strcasecmp(s1, s2) > 0) {
        position[n]++;
        }
      else {
        position[m]++;
        }
      free(s2);
      }
    free(s1);
    }

  for(size_t n = 0; n < nmgr; n++) {
    list[position[n]]=n;
    }

  return (list);
  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int mgr_WaveOrder(vector<mgr_t> mgr_serie, int nmgr, string orderCrit)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  for(size_t m = 0; m < nmgr; m++){
    
    tidal_wave w = wNUL;
    list<tidal_wave> current(0);
    
    for(size_t j = 0; j != mgr_serie[m].nwave; j++){
      w = wave_from_name(mgr_serie[m].data[j].constituent.name);
      if(strcmp(w.name, "NUL") != 0){
        w.init();
        current.push_back(w);
        }
      }
    if(current.size() != mgr_serie[m].nwave){ // internal check
#if defined(DEBUG)
      printf("ERROR in mgr_order(): %d != %s", current.size(), mgr_serie[m].nwave);
#endif
      return(1);
      }

    if(orderCrit=="INVERSE"){
      current.reverse();
      }

    if(orderCrit=="FREQUENCY"){
      current.sort(wcompare);
      }

    if(orderCrit=="ALPHANUM"){
/*      current.sort(compare_nocase);*/
      return(1);
      }

    std::vector<size_t> perm(0);
    for(list<tidal_wave>::iterator it = current.begin(); it != current.end(); it++){
      int i = mgr_serie[m].wave_index(it->name);
      if(i != -1) {
        perm.push_back(i);
        }
      else{
        return(1);
        }
      }

    if(perm.size() != mgr_serie[m].nwave){ // internal check
#if defined(DEBUG)
      printf("ERROR in mgr_order(): %d != %s", perm.size(), mgr_serie[m].nwave);
#endif
      return(1);
      }

    mgr_data_t *bck = new mgr_data_t[mgr_serie[m].nwave];
    for(size_t i = 0; i != mgr_serie[m].nwave; i++){
      bck[i] = mgr_serie[m].data[i];
      }
    for(size_t i = 0; i != mgr_serie[m].nwave; i++){
      mgr_serie[m].data[i] = bck[perm[i]];
      }

    current.clear();
    delete [] bck;
    perm.clear();
    }

  return(0);

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int mgr_reduce(vector<mgr_t> & mgr, char *keep)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int m;
  vector<mgr_t> tmp;
  
  int count=0;
  for (m=0;m<mgr.size();m++) {
    if(keep[m]!=0) {
      count++;
      }
    }
  
  count=0;
  for (m=0;m<mgr.size();m++) {
    if(keep[m]!=0) {
      tmp.push_back(mgr[m]);
      count++;
      }
    else {
      mgr[m].clear();
      }
    }
  
  mgr=tmp;
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int mgr_concat(vector<mgr_t> & mgr, vector<mgr_t> additionals)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int m;
  
  for (m=0;m<additionals.size();m++) {
    mgr.push_back(additionals[m]);
    }
  
  printf("mgr_concat: %d mgr's\n",mgr.size());
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int mgr_exclude(vector<mgr_t> & mgr, vector<plg_t> & polygons, int position)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  excludes gauges with respect to their position relative to a polygon
  
  polygons whose elements will be plg_recale()'d
  
  position if PLG_POINT_EXTERIOR, exclude those outside
            if PLG_POINT_INTERIOR exclude those inside

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int m,status;
  int *inside=0;
  char *keep=0;
  
  if(polygons.size()==0) return(0);
  
  inside=mgr_select(polygons ,mgr);
  
  keep = new char[mgr.size()];
  
  for (m=0;m<mgr.size();m++) {
    char *keepm=&keep[m];
    
    if(inside!=0 and inside[m]==position)
      *keepm=0;
    else
      *keepm=1;
    
    }
  
  status=mgr_reduce(mgr, keep);
  
  delete[] keep;
  delete[] inside;
  
  printf("mgr_exclude: %d mgr's left\n",mgr.size());
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int mgr_exclude(vector<mgr_t> & mgr, const char *poly, int position)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  vector<plg_t> polygons;
  int *inside;
  
  inside=0;
  if(poly!=0) {
    plg_load( poly, PLG_FORMAT_SCAN, polygons);
    }
  
  status=mgr_exclude(mgr,polygons,position);
  
  plg_destroy_entries(polygons);
  
  return status;
}


