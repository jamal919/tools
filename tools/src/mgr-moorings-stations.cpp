
/*******************************************************************************
T-UGO tools, 2006-2013

  Unstructured Ocean Grid initiative
  
E-mail: florent.lyard@legos.obs-mip.fr
*******************************************************************************/
/**
\file

\author  Florent Lyard      LEGOS/CNRS, Toulouse, France
\author  Yoann Le Bars      LEGOS, Toulouse, France (PhD)
\author  Laurent Roblou     LEGOS/CNRS, Toulouse, France
\author  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
\author  Frédéric Dupont    Université de Laval à Québec, Canada
\author  Clement MAYET      LEGOS, Toulouse, France (PhD)
\author  Yves Soufflet      LEGOS, Toulouse, France
\author  Sara Fleury        LEGOS, Toulouse, France

*******************************************************************************/

#include "config.h"

#include <stdio.h>
#include <string.h>

#include "tools-structures.h"

int mgr_check_morring_exists(mooring_t *mooring_new, int nb_moorings_old, mooring_t *moorings_old);

static int isfloat(char *str)
{
  char c;
  for (c=0;c<strlen(str);c++)
    if (!isdigit(str[c]) && str[c] != '-' && str[c] != '.')
      return 0;
  return 1;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int mgr_read_moorings(const char *file_in, mooring_t **sample_pos)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*
 \date 2013-01-17 Sara Fleury : created
 \param Input ASCII file. Various formats authorized: [C] [NAME] LON LAT [NAME].
 Lines starting with '#' are ignored. First available line must be the number of stations (integer).
 \param pointer on a  mooring_t table for the output
 \return nb of read stations
*/
#define DEBUG 0
{
  int k, nst=0;
  float lon,lat;
  char *line, name[64], num_name[64];
  FILE *in;
  char c;
  char c_str[64]="";
  mooring_t *local;
  enum  FORMAT {LON_LAT=0,
                LON_LAT_NAME=1,
                NAME_LON_LAT=2,
                C_LON_LAT=3,
                C_LON_LAT_NAME=4,
                C_NAME_LON_LAT=5};
  FORMAT  format;

  line =new char[256];
  
  in=fopen(file_in,"r");
  if(in==NULL) {
    STDOUT_BASE_LINE_FUNC("fopen(\"%s\",\"r\") error (%d %s)\n",file_in,errno,strerror(errno));
    return(0);
    }

/*-------------------------------------------------------------------------------
  read number of stations */
  do {
    fgets(line,256,in);
#ifdef DEBUG
    STDOUT_BASE_LINE(" %s", line);
#endif
    } while (strlen(line)<1 || line[0] == '#');

  if (sscanf(line, "%d ", &nst) != 1) {
    fclose(in);
    STDOUT_BASE_LINE("error format for file %s (no nb of stations), abort... \n",file_in);
    return(0);
    }
  STDOUT_BASE_LINE("%d stations\n", nst);
  
  local=new mooring_t[nst];
  
/*-------------------------------------------------------------------------------
  detect format from first station */
  do {
    fgets(line,256,in);
    } while (line[0] == '#');

  name[0] = '#'; name[1] = '\0';
  if (sscanf(line, "%s %s %f %f", &c_str, name, &lon, &lat) == 4 && strlen(c_str) == 1
      && !isfloat(c_str) && !isfloat(name))
    format = C_NAME_LON_LAT;
  else if  (sscanf(line, "%s %f %f %s", &c_str, &lon, &lat, name) == 4 && strlen(c_str) == 1
            && !isfloat(c_str) && !isfloat(name)) {
    if (name[0] == '#')
      format = C_LON_LAT;
    else
      format = C_LON_LAT_NAME;
    }
  else if  (sscanf(line, "%s %f %f", &c_str, &lon, &lat) == 3 && strlen(c_str) == 1 && !isfloat(c_str))
    format = C_LON_LAT;
  else if  (sscanf(line, "%s %f %f", name, &lon, &lat) == 3 && !isfloat(name))
    format = NAME_LON_LAT;
  else if  (sscanf(line, "%f %f %s", &lon, &lat, name) == 3  && !isfloat(name)) {
    if (name[0] == '#') {
      format = LON_LAT;
      name[0] = '#'; name[1] = '\0';
      }
    else
      format = LON_LAT_NAME;
    }
  else if  (sscanf(line, "%f %f", &lon, &lat) == 2){
    format = LON_LAT;
    name[0] = '#'; name[1] = '\0';
    }
  else {
    fclose(in);
    STDOUT_BASE_LINE("error station format for file %s, abort... \n",file_in);
    }
  if ( strlen(c_str) == 1) c=c_str[0];
#if DEBUG>0
  STDOUT_BASE_LINE(" detected format %d\n", format);
#endif

/*-------------------------------------------------------------------------------
  record first station */
  k = 0;
  if (name[0] == '#') {
    sprintf(num_name,"%d", k+1); local[k].name = strdup(num_name);
    }
  else
    local[k].name = strdup(name);
  local[k].lon = lon;
  local[k].lat = lat;
#if DEBUG>0
  STDOUT_BASE_LINE(" %d %s %f %f\n", k+1, name, lon, lat);
#endif
  
/*-------------------------------------------------------------------------------
  read each station */
  do {
    if(k==nst-1) break;
/*-------------------------------------------------------------------------------
    remove comment or blank lines */
    do {
      strcpy(line,"");
      char *s=fgets(line,256,in);
      if(s==0) break;
      } while (strlen(line)<5 || line[0] == '#');
    
    if(strlen(line)==0) {
      fclose(in);
      STDOUT_BASE_LINE("error in reading file %s, abort... \n",file_in);
      return(-1);
      }
    
    name[0] = '#'; name[1] = '\0';
    switch(format) {
    case LON_LAT:
      sscanf(line, "%f %f", &lon, &lat);
      break;
    case LON_LAT_NAME:
      sscanf(line, "%f %f %s", &lon, &lat, name);
      break;
    case NAME_LON_LAT:
      sscanf(line, "%s %f %f", name, &lon, &lat);
      break;
    case C_LON_LAT:
      sscanf(line, "%c %f %f", &c, &lon, &lat);
      break;
    case C_LON_LAT_NAME:
      sscanf(line, "%c %f %f %s", &c, &lon, &lat, name);
      break;
    case C_NAME_LON_LAT:
      sscanf(line, "%c %s %f %f", &c, name, &lon, &lat);
      break;
    }

/*-------------------------------------------------------------------------------
    record  station */
    k++;
    if (name[0] == '#') {
      sprintf(num_name,"%d", k+1); local[k].name = strdup(num_name);
      }
    else local[k].name = strdup(name);
    local[k].lon = lon;
    local[k].lat = lat;
#if DEBUG>0
    STDOUT_BASE_LINE(" %d %s %f %f\n", k+1, local[k].name, local[k].lon, local[k].lat);
#endif
  } while(k<nst-1);

  fclose(in);
  *sample_pos=local;
  return(nst);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int mgr_save_moorings(const char *file_out, int nb_moorings_new, mooring_t *moorings_new)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*
 \date 2013-01-17 Sara Fleury : created
 \param Output ASCII file NAME LON LAT. Lines starting with '#' are ignored. First available line must be the number of stations (integer).
 \param pointer on a  mooring_t table
 \return nb of read stations
*/

/* Save positions only ! */
{
  int nb_moorings_old=0;
  int nb_moorings_new_to_keep=0;
  int *num_moorings_new;
  mooring_t *moorings_old=0;
  int add;
  FILE *out;
  int nb_tot,n;
  char name[64];
#if 0
  char file_save[1024];
#endif

  num_moorings_new = new int[nb_moorings_new];
  nb_moorings_new_to_keep = nb_moorings_new;
  for (n=0;n<nb_moorings_new;n++)
    num_moorings_new[n] = n;

#if 0
  /* check if file exists */
  out=fopen(file_out,"r");
  if(out!=0) {
    fclose(out);

    /* save old file */
    sprintf(file_save, "%s.save", file_out);
    STDOUT_BASE_LINE("File already exists: %s. Saved in %s\n",  file_out, file_save);
    nb_moorings_old = mgr_read_moorings(file_out,  &moorings_old);
    mgr_save_moorings(file_save, nb_moorings_old, moorings_old);

    /* keep only unknown new moorings */
    nb_moorings_new_to_keep=0;
    for (n=0;n<nb_moorings_new;n++) {
      if (!mgr_check_morring_exists(&moorings_new[n], nb_moorings_old, moorings_old)) {
        num_moorings_new[nb_moorings_new_to_keep] = n;
        nb_moorings_new_to_keep++;
      }
    }
  }
#endif
  nb_tot = nb_moorings_old + nb_moorings_new_to_keep;

  out=fopen(file_out,"w");
  fprintf(out, "# lon lat [CODE]NAME\n");
  fprintf(out, "%d XYN\n", nb_tot);

  /* loop over old moorings to keep */
  for (n=0; n<nb_moorings_old; n++) {
    if (moorings_old[n].code == -1)
      sprintf(name, "%s", moorings_old[n].name);
    else
      sprintf(name, "%03d%s", moorings_old[n].code, moorings_old[n].name);
      fprintf(out, "%10.6f  %10.6f   %s\n", moorings_old[n].lon, moorings_old[n].lat, name);
  }  /* for old */


  /* loop over new moorings to add */
  for (add=0;add<nb_moorings_new_to_keep;add++) {
    n =  num_moorings_new[add];
    const mooring_t *moorings_newn=&moorings_new[n];
    if(moorings_newn->name==0){
      STDERR_BASE_LINE_FUNC("skipping moorings_new[%d].name==%p\n",n,moorings_newn->name);
      continue;
      }
    if (moorings_newn->code == -1)
      sprintf(name, "%s", moorings_newn->name);
    else
      sprintf(name, "%03d%s", moorings_newn->code, moorings_newn->name);
      fprintf(out, "%10.6f  %10.6f   %s\n", moorings_newn->lon, moorings_newn->lat, name);
  }  /* for new */

  fclose(out);

  return 1;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int mgr_check_morring_exists(mooring_t *mooring_new, int nb_moorings_old, mooring_t *moorings_old)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k;

  /* check if already in file list */
  for (k=0;k<nb_moorings_old;k++) {

    /* already in list */
    if (strcmp(moorings_old[k].name, mooring_new->name)==0) {
      
      /* pos identical ? */
      if (moorings_old[k].lon == mooring_new->lon && moorings_old[k].lat == mooring_new->lat)
        break;
      
      /* not identical !!! */
      STDOUT_BASE_LINE("2 positions for same mooring %s: old %f %f -> new %f %f !!\n",  mooring_new->name,
                        moorings_old[k].lon, moorings_old[k].lat,
                        mooring_new->lon, mooring_new->lat);
    }
  } /* for old */
  
  /* new mooring */
  if (k == nb_moorings_old) {
    printf("ADD %f  %f %s\n", mooring_new->lon, mooring_new->lat, mooring_new->name);
    return 0;
  }
  /* known mooring */
  return 1;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int mgr_stations_2_moorings(int nst, pressure_station *stations, mooring_t **moorings)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*
 \date 2013-03-17 Sara Fleury : created
 \comment position converter
 \param nb of stations
 \param pointer on a  pressure_station table
 \return moorings table
*/
{
  int k;

  /* alloc */
  *moorings = new mooring_t[nst];

  /* convert stations to moorings */
  for (k=0; k<nst; k++) {
    (*moorings)[k].name = strdup(stations[k].name);
    (*moorings)[k].lon = stations[k].t;
    (*moorings)[k].lat = stations[k].p;
  }

  return nst;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int mgr_moorings_2_stations(int nst, mooring_t *moorings, pressure_station **stations)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*
 \date 2013-03-17 Sara Fleury : created
 \comment position converter
 \param nb of stations
 \param pointer on a  mooring_t table
 \return pressure_station table
*/
{
  int k;

  /* alloc */
  *stations=new pressure_station[nst];

  /* convert moorings to stations */
  for (k=0; k<nst; k++) {
    const mooring_t *mooringk=&moorings[k];
    pressure_station *stationk=&(*stations)[k];
    
    strcpy(stationk->name, mooringk->name);
    stationk->t = mooringk->lon;
    stationk->p = mooringk->lat;
    }
  
  return nst;
}



/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int mgr_read_stations(const char *filename, pressure_station **sample_pos)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/*------------------------------------------------------------------------------
 \date 2013-01-17 Sara Fleury : created
 \comment interface to mgr_read_moorings
 \param Input ASCII file. Various formats authorized: [C] [NAME] LON LAT [NAME].
 Lines starting with '#' are ignored. First available line must be the number of stations (integer).
 \param pointer on a  mooring_t table for the output
 \return nb of read stations
------------------------------------------------------------------------------*/
  int nst=0;
  pressure_station *stations=NULL;

  mooring_t *moorings=0;

/*------------------------------------------------------------------------------
  read moorings */
  nst = mgr_read_moorings(filename, &moorings);
  
  if(nst==0) {
    STDOUT_BASE_LINE("opening %s failed...\n",filename);
    exit(0);
    }

/*------------------------------------------------------------------------------
  convert to stations */
  nst = mgr_moorings_2_stations(nst, moorings, &stations);

/*------------------------------------------------------------------------------
  remove moorings */
  if (moorings !=0) delete [] moorings;

  *sample_pos = stations;
  return (nst);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int mgr_save_stations(char *file_out, int nb_moorings_new, pressure_station *stations)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* Save positions only ! */
/*
 \date 2013-01-17 Sara Fleury : created
 \param Output ASCII file NAME LON LAT. Lines starting with '#' are ignored.
 First available line must be the number of stations (integer).
 \param pointer on a  mooring_t table
 \return nb of read stations
*/

{
  mooring_t *moorings=NULL;
  int ret;

  nb_moorings_new = mgr_stations_2_moorings(nb_moorings_new, stations, &moorings);

  ret = mgr_save_moorings(file_out, nb_moorings_new, moorings);

  if (moorings !=0) delete [] moorings;

  return ret;
}
