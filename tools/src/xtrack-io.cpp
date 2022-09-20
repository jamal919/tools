
/*******************************************************************************

  T-UGO tools, 2006-2017

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

\brief X-TRACK I/O function definitions
*/
/*----------------------------------------------------------------------------*/

#include <config.h>

#include <iostream>
#include <vector>

#include "tools-structures.h"

#include "tides.h"
#include "fe.h"
#include "archive.h"
#include "poc-time.h"
#include "poc-netcdf-data.hpp"
#include "mgr.h"
#include "functions.h"
#include "statistic.h"
#include "xtrack-io.h"

#ifndef VERBOSE
#define VERBOSE 1
#endif


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int formatage_CTOH (int nb_traces, track_data *data, const char *dir_ah, char *convention)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,j,k,cpt;
  char *filename=NULL;
  FILE *f_data;
  double *Hmer=NULL;
  double *Temps=NULL;
  double *DAC=NULL;
  double *Tide=NULL;
  double temp=99.9999;
  int num_trace;

  char mesg[70]= "Problème d'allocation mémoire - Programme formatage";
  if ((filename=(char *) malloc(100*sizeof(char))) == NULL) perror(mesg);

  printf("Formation des fichiers d'entree de l'outil de l analyse harmonique \n\n");

  for (k=0;k<nb_traces;k++) {
    num_trace=data[k].num_trace;
    printf("Num trace : %d\n", num_trace);
    if (num_trace<100 && num_trace>=10) {
      sprintf(filename,"%s/%s.0%d.dat",dir_ah,convention,num_trace);
      }
    else if (num_trace <10) {
      sprintf(filename,"%s/%s.00%d.dat",dir_ah,convention,num_trace);
      }
    else {
      sprintf(filename,"%s/%s.%d.dat",dir_ah,convention,num_trace);
      }

    printf("Nom de fichier : %s \n", filename);

    f_data=fopen(filename,"w");
    if(f_data==0) TRAP_ERR_RETURN(errno,1,"Could not (re)create .dat file for track %d (%d %s)\n",num_trace,errno,strerror(errno));

    fprintf(f_data,"#-- HEADER -------------------------------------------\n");
    fprintf(f_data,"# Column 1 : date in days referred to\n");
    fprintf(f_data,"#            CNES date (01-JAN-1950 00:00:00.0)\n");
    fprintf(f_data,"# Column 2 : sea surface height anomaly (in meters)\n");
    fprintf(f_data,"# Column 3 : MOG2D-G model sea level (in meters)\n");
    fprintf(f_data,"# Column 4 : MOG2D-G model S1 atmospheric tide (in meters)\n");
    fprintf(f_data,"# Column 5 : MOG2D-G model S2 atmospheric tide (in meters)\n");
    fprintf(f_data,"# Column 6 : MOG2D regional model sea level (in meters)\n");
    fprintf(f_data,"# Column 7 : MOG2D regional model S1 atmospheric tide (in meters)\n");
    fprintf(f_data,"# Column 8 : MOG2D regional model S2 atmospheric tide (in meters)\n");
    fprintf(f_data,"# Column 9 : inverted barometer (in meters)\n");
    fprintf(f_data,"# Column 10: inverted barometer S1 atmospheric tide (in meters)\n");
    fprintf(f_data,"# Column 11: inverted barometer S2 atmospheric tide (in meters)\n");
    fprintf(f_data,"# Column 12: loading effects (in meters)\n");
    fprintf(f_data,"# Column 13: solid Earth tide (in meters)\n");
    fprintf(f_data,"# Column 14: MOG2D regional model ocean tide (in meters)\n");
    fprintf(f_data,"# Column 15: MOG2D regional model S1 ocean tide (in meters)\n");
    fprintf(f_data,"# Column 16: MOG2D regional model S2 ocean tide (in meters)\n");
    fprintf(f_data,"# Column 17: GOT00 model ocean tide (in meters)\n");
    fprintf(f_data,"# Column 18: GOT00 model S1 ocean tide (in meters)\n");
    fprintf(f_data,"# Column 19: GOT00 model S2 ocean tide (in meters)\n");
    fprintf(f_data,"# Column 20: FES2004 model ocean tide (in meters)\n");
    fprintf(f_data,"# Column 21: FES2004 model S1 ocean tide (in meters)\n");
    fprintf(f_data,"# Column 22: FES2004 model S2 ocean tide (in meters)\n");
    fprintf(f_data,"# Column 23: harmonic analysis ocean tide (in meters)\n");
    fprintf(f_data,"# Column 24: harmonic analysis S1 ocean tide (in meters)\n");
    fprintf(f_data,"# Column 25: harmonic analysis S2 ocean tide (in meters)\n");
    fprintf(f_data,"#-- HEADER END ---------------------------------------\n");
    fprintf(f_data,"Number of points :\n");
    fprintf(f_data,"%d\n",data[k].nbpoints);
    fprintf(f_data,"#-----------------------------------------------------\n");

    for (i=0;i<data[k].nbpoints;i++) {
      cpt=0;
      /*Suppression des données manquantes, ie valeur par défaut*/
      for (j=0;j<data[k].nbcycles;j++) {
        /*exclusion des données manquantes pour faire l'analyse harmonique : la valeur à écarter est 99.99*/
        if (data[k].time[i][j]<99 || data[k].time[i][j]>100) {
          cpt++;
          if ((Hmer=(double*)realloc(Hmer,sizeof(double)*(cpt))) == NULL) perror(mesg);
          if ((Temps=(double*)realloc(Temps,sizeof(double)*(cpt))) == NULL) perror(mesg);
          if ((Tide=(double*)realloc(Tide,sizeof(double)*(cpt))) == NULL) perror(mesg);
          if ((DAC=(double*)realloc(DAC,sizeof(double)*(cpt))) == NULL) perror(mesg);
          Tide[cpt-1]=(double)data[k].tide[i][j];
          DAC[cpt-1]=(double)data[k].mog2d[i][j];
          if(Tide[cpt-1]>=99) {
            Hmer[cpt-1]=99.9999;
            }
          else {
            Hmer[cpt-1]=(double)data[k].sla[i][j]+(double)data[k].tide[i][j]+(double)data[k].mog2d[i][j];
            }
          Temps[cpt-1]=data[k].time[i][j];
          }
        }
      /*Ecriture dans les fichiers .dat*/
      fprintf(f_data,"Pt  : %d\n",i+1);
      fprintf(f_data,"lon : %f\n",data[k].lon[i]);
      fprintf(f_data,"lat : %f\n",data[k].lat[i]);
      fprintf(f_data,"Mes : %d\n",cpt);
      for (j=0;j<cpt;j++) {
        fprintf(f_data,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
                Temps[j],Hmer[j],DAC[j],temp,temp,temp,temp,temp,temp,temp,temp,0.0,0.0,temp,temp,temp,temp,temp,temp,Tide[j],temp,temp,temp,temp,temp);
        }

      fprintf(f_data,"#-----------------------------------------------------\n");
      }
    }

  free(filename);
  free(Hmer);
  free(Temps);
  free(Tide);
  free(DAC);

  return (0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int formatage_DUACS (int nb_traces, track_data *data, char *dir_ah)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,j,k,cpt, cpt_points;
  char *filename=NULL;
  FILE *f_data;
  double *Hmer=NULL;
  double *Temps=NULL;
  double *DAC=NULL;
  double *Tide=NULL;
  double temp=0;
  int num_trace;

  char mesg[70]= "Problème d'allocation mémoire - Programme formatage";
  if ((filename=(char *) malloc(100*sizeof(char))) == NULL) perror(mesg);

  printf("Formation des fichiers d'entree de l'outil de l analyse harmonique \n\n");

  for (k=0;k<nb_traces;k++) {
    num_trace=data[k].num_trace;
    cpt_points=0;
    printf("Num trace : %d\n", num_trace);
    if (num_trace<100 && num_trace>=10) {
      sprintf(filename,"%s/track-ref.TP+J1.0%d.dat",dir_ah,num_trace);
    }
    else if (num_trace <10) {
      sprintf(filename,"%s/track-ref.TP+J1.00%d.dat",dir_ah,num_trace);
    }
    else {
      sprintf(filename,"%s/track-ref.TP+J1.%d.dat",dir_ah,num_trace);
    }

    printf("Nom de fichier : %s \n", filename);

    f_data=fopen(filename,"w");
    if(f_data==0) TRAP_ERR_RETURN(errno,1,"Could not (re)create .dat file for track %d (%d %s)\n",num_trace,errno,strerror(errno));

    fprintf(f_data,"#-- HEADER -------------------------------------------\n");
    fprintf(f_data,"# Column 1 : date in days referred to\n");
    fprintf(f_data,"#            CNES date (01-JAN-1950 00:00:00.0)\n");
    fprintf(f_data,"# Column 2 : sea surface height anomaly (in meters)\n");
    fprintf(f_data,"# Column 3 : MOG2D-G model sea level (in meters)\n");
    fprintf(f_data,"# Column 4 : MOG2D-G model S1 atmospheric tide (in meters)\n");
    fprintf(f_data,"# Column 5 : MOG2D-G model S2 atmospheric tide (in meters)\n");
    fprintf(f_data,"# Column 6 : MOG2D regional model sea level (in meters)\n");
    fprintf(f_data,"# Column 7 : MOG2D regional model S1 atmospheric tide (in meters)\n");
    fprintf(f_data,"# Column 8 : MOG2D regional model S2 atmospheric tide (in meters)\n");
    fprintf(f_data,"# Column 9 : inverted barometer (in meters)\n");
    fprintf(f_data,"# Column 10: inverted barometer S1 atmospheric tide (in meters)\n");
    fprintf(f_data,"# Column 11: inverted barometer S2 atmospheric tide (in meters)\n");
    fprintf(f_data,"# Column 12: loading effects (in meters)\n");
    fprintf(f_data,"# Column 13: solid Earth tide (in meters)\n");
    fprintf(f_data,"# Column 14: MOG2D regional model ocean tide (in meters)\n");
    fprintf(f_data,"# Column 15: MOG2D regional model S1 ocean tide (in meters)\n");
    fprintf(f_data,"# Column 16: MOG2D regional model S2 ocean tide (in meters)\n");
    fprintf(f_data,"# Column 17: GOT00 model ocean tide (in meters)\n");
    fprintf(f_data,"# Column 18: GOT00 model S1 ocean tide (in meters)\n");
    fprintf(f_data,"# Column 19: GOT00 model S2 ocean tide (in meters)\n");
    fprintf(f_data,"# Column 20: FES2004 model ocean tide (in meters)\n");
    fprintf(f_data,"# Column 21: FES2004 model S1 ocean tide (in meters)\n");
    fprintf(f_data,"# Column 22: FES2004 model S2 ocean tide (in meters)\n");
    fprintf(f_data,"# Column 23: harmonic analysis ocean tide (in meters)\n");
    fprintf(f_data,"# Column 24: harmonic analysis S1 ocean tide (in meters)\n");
    fprintf(f_data,"# Column 25: harmonic analysis S2 ocean tide (in meters)\n");
    fprintf(f_data,"#-- HEADER END ---------------------------------------\n");
    fprintf(f_data,"Number of points :\n");
    fprintf(f_data,"%d\n",data[k].nbpoints);
    fprintf(f_data,"#-----------------------------------------------------\n");


    for (i=0;i<data[k].nbpoints;i++) {
      cpt=0;
      /*Suppression des données manquantes, ie valeur par défaut*/
      for (j=0;j<data[k].nbcycles;j++) {
        /*exclusion des données manquantes pour faire l'analyse harmonique : la valeur à écarter est 99.99*/
        if (data[k].time[i][j]<99 || data[k].time[i][j]>100) {
          cpt++;
          if ((Hmer=(double*)realloc(Hmer,sizeof(double)*(cpt))) == NULL) perror(mesg);
          if ((Temps=(double*)realloc(Temps,sizeof(double)*(cpt))) == NULL) perror(mesg);
          if ((Tide=(double*)realloc(Tide,sizeof(double)*(cpt))) == NULL) perror(mesg);
          if ((DAC=(double*)realloc(DAC,sizeof(double)*(cpt))) == NULL) perror(mesg);
          Tide[cpt-1]=0;
          DAC[cpt-1]=0;
          Hmer[cpt-1]=(double)data[k].sla[i][j];
          Temps[cpt-1]=data[k].time[i][j];
        }
      }
      /*Ecriture dans les fichiers .dat*/
      // fprintf(f_data,"Pt  : %d\n",i+1);
      fprintf(f_data,"Pt  : %d\n",cpt_points++);
      fprintf(f_data,"lon : %f\n",data[k].lon[i]);
      fprintf(f_data,"lat : %f\n",data[k].lat[i]);
      fprintf(f_data,"Mes : %d\n",cpt);
      for (j=0;j<cpt;j++) {
        fprintf(f_data,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
                Temps[j],Hmer[j],DAC[j],temp,temp,temp,temp,temp,temp,temp,temp,temp,temp,temp,temp,temp,temp,temp,temp,Tide[j],temp,temp,temp,temp,temp);
      }
      fprintf(f_data,"#-----------------------------------------------------\n");
    }

  }

  free(filename);
  free(Hmer);
  free(Temps);
  free(Tide);
  free(DAC);

  return (0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  string XTRACK_get_fileFormat(const string & input)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  string format = "";
  size_t pos = string::npos;

  // choose format
  pos = input.rfind(".dat");
  if( pos != string::npos){
    format = "ASCII";
    }

  pos = input.rfind(".list");
  if( pos != string::npos){
    format = "ASCII";
    }

  pos = input.rfind(".nc");
  if( pos != string::npos){
    format = "NETCDF";
    }
  
  return(format);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int XTRACK_get_asciiFormat(const string & input)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  char line[200];
  FILE *f = fopen(input.c_str(), "r");
  
  // skip header
  do  fgets(line, sizeof(line), f);
  while (strncmp(line,"#-- HEADER END ---------------------------------------",20) != 0 );
  fgets(line, sizeof(line), f); /* Number of... */
  fgets(line, sizeof(line), f); /* value */
  fgets(line, sizeof(line), f); /*#--------------------- */

  // read first record
  fgets(line, sizeof(line), f); /* Pt : */
  fgets(line, sizeof(line), f); /* lon : */
  fgets(line, sizeof(line), f); /* lat : */
  fgets(line, sizeof(line), f); /* Mes : */

  int ncol = 0;
  fgets(line, sizeof(line), f); /* dataset */
  double a = 1e-35, b = 1e-35, c = 1e-35;
  int nitems = sscanf(line, " %lf %lf %lf", &a, &b, &c);

  // assume more than 3 values means full line is written (31 values)
  // TODO do it better later
  ncol = (nitems == 2) ? nitems : X_TRACK_RAW_NUMVAL;
  
  fclose(f);
  
  return(ncol);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int XTRACK_ref_get_nRecords_ascii(const string & input)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  FILE *fic_data = 0;

  fic_data = fopen(input.c_str(),"r");
  if(fic_data!=0) {

    int nRecords;
    char line[300];

    do  fgets(line,sizeof(line),fic_data);
    while (strncmp(line,"#-- HEADER END ---------------------------------------",20) != 0 );

    fgets(line,sizeof(line),fic_data);
    fgets(line,sizeof(line),fic_data);
    sscanf(line, "%d", &nRecords);

    fclose(fic_data);
    return(nRecords);
    }
  else {
    cout << "ERROR : can not open : "<< input << endl;
    return(0);
    }

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int XTRACK_ref_get_nRecords_netcdf(const string & input)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status = -1;
  int ncid = -1;
  status = nc_open(input.c_str(),NC_NOWRITE,&ncid);
  if(status!=0) {
    fprintf(stderr, "%s\n", nc_strerror(status));
    return(0);
    }

  size_t nRecords = 0;
  int nRecords_id = -1;
  string nRecords_dimName = "nbpoints";

  status = nc_inq_dimid(ncid,nRecords_dimName.c_str(),&nRecords_id);
  if(status!=0) {
    fprintf(stderr, "%s\n", nc_strerror(status));
    return(0);
   }

  status = nc_inq_dimlen(ncid,nRecords_id,&nRecords);
  if(status!=0) {
    fprintf(stderr, "%s\n", nc_strerror(status));
    return(0);
    }

  status = nc_close(ncid);
  if(status!=0) {
    fprintf(stderr, "%s\n", nc_strerror(status));
    return(0);
    }
  else{
    return(nRecords);
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int XTRACK_ref_get_nRecords(const string & input)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  string format = XTRACK_get_fileFormat(input);
  
  if(format.empty()){
    return(1);
    }
  
  // get nRecords
  if(format=="ASCII"){
    format.clear();
    return(XTRACK_ref_get_nRecords_ascii(input));
    }

  if(format=="NETCDF"){
    format.clear();
    return(XTRACK_ref_get_nRecords_netcdf(input));
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

std::vector<int> XTRACK_ref_get_indexRecords_ascii(const string & input)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  FILE *fic_data = 0;
  std::vector<int> index(0);
  int nRecords = 0;
  char line[300], dummy[10];
  int nitems = 0, count = 0;

  fic_data = fopen(input.c_str(),"r");
  if(fic_data!=0) {
    do  fgets(line,sizeof(line),fic_data);
    while (strncmp(line,"#-- HEADER END ---------------------------------------",20) != 0 );

    fgets(line,sizeof(line),fic_data);
    fgets(line,sizeof(line),fic_data);
    sscanf(line, "%d", &nRecords);
    }
  else {
    cout << "ERROR : can not open : "<< input << endl;
    return(index);
    }

  fgets(line,sizeof(line),fic_data);
  while(!feof(fic_data)){
    if( strncmp(line,"Pt ",3) == 0 ){
      if( (nitems = sscanf(line,"%s : %d",dummy,&count)) != 1){
        index.push_back(count);
        }
      }
    fgets(line,sizeof(line),fic_data);
    }

  // internal check
  if(index.size() == nRecords){
    fclose(fic_data);
    return(index);
    }
  else {
    index.clear();
    fclose(fic_data);
    return(index);
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  std::vector<int> XTRACK_ref_get_indexRecords_netcdf(const string & input)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status = -1;
  int ncid = -1;
  
  status = nc_open(input.c_str(),NC_NOWRITE,&ncid);
  NC_TRAP_ERROR(wexit,status,1, " file=%s, nc_open() failed", input.c_str());
 
  // get records number
  int nRecords_id = -1;
  string s = "nbpoints";
  status = nc_inq_dimid(ncid,s.c_str(),&nRecords_id);
  NC_TRAP_ERROR(wexit,status,1, " file=%s, nRecords_id, nc_inq_dimid() failed", input.c_str());

  size_t nRecords = 0;
  status = nc_inq_dimlen(ncid,nRecords_id,&nRecords);
  NC_TRAP_ERROR(wexit,status,1," file=%s, nRecords_id, nc_inq_dimlen() failed", input.c_str());
  
  /// LOOSE INDICES FOR NETCDF/ASCII COMPATIBILITY REASONS
    vector<int> index(nRecords);
    for(size_t i = 0; i < nRecords; ++i){
      index[i] = i;
      }
  
    s.clear();

    return(index);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  std::vector<int> XTRACK_ref_get_indexRecords(const string & input)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  string format = XTRACK_get_fileFormat(input);
  vector<int> indices;
  
  if(format.empty()){
    return(std::vector<int>(0));
    }

  if(format=="ASCII"){
    format.clear();
    indices=XTRACK_ref_get_indexRecords_ascii(input);
    return(indices);
    }
  if(format=="NETCDF"){
    format.clear();
    indices=XTRACK_ref_get_indexRecords_netcdf(input);
    return(indices);
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  vector<serie_t> ctoh_load_metadata(const string & input, int& status)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  string format = XTRACK_get_fileFormat(input);
  
  vector<serie_t> series;
  
  // get nRecords
  if(format=="ASCII"){
    series=ctoh_load_metadata_ASCII(input.c_str(), status);
    return(series);
    }
  else if(format=="NETCDF"){
    series=ctoh_load_metadata_NetCDF(input);
    return(series);
    }
  else  TRAP_ERR_EXIT(-1,"unknown format\n");
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  serie_t XTRACK_ref_get_record_ascii(const char *input, int point, int& status)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int    i;
  int    npoints,nitems,nmes,count;
  double lon,lat;
  char   line[1024],dummy[1024];
  FILE   *in=NULL;

  in = fopen(input,"r");
  if(in==0) {
    cout << "ERROR : can not open : "<< input << endl;
    }

  do  fgets(line, sizeof(line), in);
  while (strncmp(line,"#-- HEADER END ---------------------------------------",20) != 0 );
  fgets(line, sizeof(line), in); /* Number of... */
  if( (nitems = fscanf(in,"%d\n",&npoints)) !=1) {
    cout << "ERROR : can not read : "<< line << endl;
    }

 read:
  for(i = 0; i < 1; i++) fgets(line,sizeof(line),in);

  for(i = 0; i < 1; i++) fgets(line,sizeof(line),in);
  if( (nitems = sscanf(strstr(line," : "),"%s %d",dummy,&count)) != 2){
    cout << "ERROR : can not read : "<< line << endl;
    }
  for(i = 0; i < 1; i++) fgets(line,sizeof(line),in);
  if( (nitems = sscanf(strstr(line," : "),"%s %lf",dummy,&lon)) != 2){
    cout << "ERROR : can not read : "<< line << endl;
    }
  for(i = 0; i < 1; i++) fgets(line,sizeof(line),in);
  if( (nitems = sscanf(strstr(line," : "),"%s %lf",dummy,&lat)) != 2){
    cout << "ERROR : can not read : "<< line << endl;
    }
  for(i = 0; i < 1; i++) fgets(line,sizeof(line),in);
  if( (nitems = sscanf(strstr(line," : "),"%s %d",dummy,&nmes)) != 2){
    cout << "ERROR : can not read : "<< line << endl;
    }

  serie_t metadata;
  if(count < point) {
    /* *----------------------------------------------------------------------
       skip current point*/
    for(i = 0; i < nmes;  i++) fgets(line,sizeof(line),in);
    goto read;
    }
  else {
    if(count == point) {
      //printf("file %s : %d %d %lf %lf\n",input,count,nmes,lon,lat); //DBG
      data_t *buffer = new data_t[nmes];
      int j=0;
      int k = 0;
      while(k < nmes) {
        buffer[k].lon = lon;
        buffer[k].lat = lat;
        nitems = fscanf(in,"%lf",&(buffer[k].time));
        if(nitems != 1) {
          printf("********** %d %d\n",count,k);
          break;
          }

        for(j = 0; j < 24; j++) {
          nitems = fscanf(in,"%f",&(buffer[k].values[j]));
          if(nitems != 1) {
            printf("************** %d %d %d\n",count,k,j);
            break;
            }
          }
        k++;
        }
#if VERBOSE > 0
      printf("%s starts %s, finishes %s \n",input,sgetcnesdate(buffer[0].time*24.),sgetcnesdate(buffer[nmes-1].time*24.));
      fflush(stdin);
#endif
      metadata.count = nmes;
      metadata.mask = 99.9999f;
      metadata.data = new data_t[nmes];
      for(k = 0; k < nmes; k++) metadata.data[k] = buffer[k];
      delete [] buffer;

      fgets(line,sizeof(line),in);
      fclose(in);
      status=0;
      return(metadata);
      }
    else {
      if(count>point) {  // nominal location #point is not in the file!
        status=1;
        fclose(in);
        metadata.count=0;
        metadata.data=0;
        return(metadata);
        }
      }
    }

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  vector<serie_t> ctoh_load_metadata_ASCII(const char *input, int& status)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
// dealing with 2 types of ASCII format
  serie_t record;
  vector<serie_t> track;
  int n = XTRACK_ref_get_nRecords_ascii(input);
  
  int ncol = XTRACK_get_asciiFormat(input);
  switch(ncol){
    
    case 2: {
      // LIST format (2 columns)
      // init some variables
      int nmesMax = 2000; // TODO compute max number of data in time series (use XTRACK_ref_get_header_ascii() )
      double *ttime = new double[nmesMax];
      double *sealevel = new double[nmesMax];
      int scale_time = 1; // there is no time unit conversion when constant is 1
      float mask = 99.9999f;
      double lon, lat;
      FILE *file = fopen(input, "r");
      
      // skip file header (no need to returned value here)
      skipHeader(file);
    
      // read record #i
      for(size_t i = 0; i < n; ++i){

        int nmes = aktarus_load_serie(ttime, sealevel, &lon, &lat, file, scale_time);
            
        if(nmes < 0){ // reading failed, close, free and exit
          status = 1;
          fclose(file);
          delete [] ttime;
          delete [] sealevel;
          return(track);
          }
        else { // complete class instance
          record.data = new data_t[nmes];
          record.count = nmes;
          record.mask = mask;
          for(size_t k = 0; k < nmes; ++k){
            record.data[k].lon = lon;
            record.data[k].lat = lat;
            record.data[k].time = ttime[k];
            record.data[k].values[X_TRACK_SERIE_T_SSHA] = sealevel[k];
            for(size_t l = 1; l < X_TRACK_SERIE_T_NUMVAL; ++l){
              record.data[k].values[l] = 0;
              }
            }
          track.push_back(record);   // push back record in list
          delete [] record.data;       // clean local memory
          record.data = 0;
          }
        }
      break;
      }
    
    case X_TRACK_RAW_NUMVAL: {
      // X-TRACK format (31 columns)
      std::vector<int> idx = XTRACK_ref_get_indexRecords(input);    // get record indices
      for(size_t i = 0; i < n; ++i){
        record = XTRACK_ref_get_record_ascii(input, idx[i], status); // read record #i
        track.push_back(record);   // push back record in list
        delete [] record.data;     // clean local memory
        record.data = 0;
        }
      break;
      }
    
    default: {
      status = 1;
      return(track);
      }
    
    }
  
    if(track.size() != n){
    status = 1;
    return(track);
    }
  status = 0;
  return(track);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  serie_t XTRACK_ref_get_record_netcdf(const char *input, int point, int& status)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  vector<serie_t>track = ctoh_load_metadata_NetCDF(input);
 
  serie_t record;
  if(track.size() == 0){
    status = 1;
    record.count = 0;
    record.data  = 0;
    }
  else {
    record = track[point];
    status = 0;
    }

  track.clear();
  return(record);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  serie_t XTRACK_ref_get_record(const string & input, int point, int& status)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  string format = XTRACK_get_fileFormat(input);
  
  // get nRecords
  if(format=="ASCII"){
    format.clear();
    return(XTRACK_ref_get_record_ascii(input.c_str(), point, status));
    }
  else if(format=="NETCDF"){
    format.clear();
    return(XTRACK_ref_get_record_netcdf(input.c_str(), point, status));
    }
  else  TRAP_ERR_EXIT(-1,"unknown format\n");
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int XTRACK_ref_get_header_ascii(const string & input, int point, double &x, double &y, int &n)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int    i;
  int    npoints,nitems,nmes,count;
  double lon,lat;
  char   line[1024],dummy[1024];
  FILE   *in=NULL;
  int status = 1;

  in = fopen(input.c_str(),"r");
  if(in==0) {
    cout << "ERROR : can not open : "<< input << endl;
    return(status);
    }

  // read number of points in file
  do  fgets(line, sizeof(line), in);
  while (strncmp(line,"#-- HEADER END ---------------------------------------",20) != 0 );
  fgets(line,sizeof(line),in); /* Number of points : */
  if( (nitems = fscanf(in, "%d\n", &npoints)) !=1) {
    cout << "ERROR : can not read : "<< line <<endl;
    return(status);
    }

  // read block data
 read:

  fgets(line,sizeof(line),in);

  fgets(line,sizeof(line),in);
  if( (nitems = sscanf(strstr(line," : "),"%s %d",dummy,&count)) != 2){
    cout << "ERROR : can not read : "<< line << endl;
    return(status);
    }
  fgets(line,sizeof(line),in);
  if( (nitems = sscanf(strstr(line," : "),"%s %lf",dummy,&lon)) != 2){
    cout << "ERROR : can not read : "<< line << endl;
    return(status);
    }
  fgets(line,sizeof(line),in);
  if( (nitems = sscanf(strstr(line," : "),"%s %lf",dummy,&lat)) != 2){
    cout << "ERROR : can not read : "<< line << endl;
    return(status);
    }
  fgets(line,sizeof(line),in);
  if( (nitems = sscanf(strstr(line," : "),"%s %d",dummy,&nmes)) != 2){
    cout << "ERROR : can not read : "<< line << endl;
    return(status);
    }

  if(count < point) {
    /* *----------------------------------------------------------------------
       skip current point*/
    for(i = 0; i < nmes;  i++) fgets(line,sizeof(line),in);
    goto read;
    }
  else {
    if(count == point) {
      x = lon;
      y = lat;
      n = nmes;
      fgets(line,sizeof(line),in);
      fclose(in);
      return(0);
      }
    else {
      if(count > point) {  // nominal location #point is not in the file!
        status = 1;
        fclose(in);
        return(status);
        }
      }
    }

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int XTRACK_ref_get_header_netcdf(const string & input, int point, double &x, double &y, int &n)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  vector<serie_t> track = ctoh_load_metadata_NetCDF(input);
  x = track[point].data[0].lon;
  y = track[point].data[0].lat;
  n = track[point].count;
  
  track.clear();

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void XTRACK_ref_get_header(const string & input, int point, double &x, double &y, int &n, int &status)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  string format = XTRACK_get_fileFormat(input);

  // get header info
  if(format=="ASCII"){
    status = XTRACK_ref_get_header_ascii(input,point,x,y,n);
    }
  if(format=="NETCDF"){
    status = XTRACK_ref_get_header_netcdf(input,point,x,y,n);
    }
  format.clear();
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int XTRACK_ref_get_AllHeader_ascii(const string & input, std::vector<int> id, double *x, double *y, int *n)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status = 1; // error status
  char line[1024],dummy[1024];
  int nitems, npoints;
  
  FILE   *in=NULL;
  in = fopen(input.c_str(),"r");
  if(in==0) {
    cout << "ERROR : can not open : "<< input << endl;
    return(status);
    }

  // read number of points in file
  do  fgets(line, sizeof(line), in);
  while (strncmp(line,"#-- HEADER END ---------------------------------------",20) != 0 );
  
  fgets(line,sizeof(line),in); /* Number of points : */
  if( (nitems = fscanf(in, "%d\n", &npoints)) !=1) {
    cout << "ERROR : can not read : "<< line <<endl;
    return(status);
    }
  if(npoints != id.size()){
    return(status);
    }
  
  // read record headers
  for(size_t record = 0; record < npoints; record++){
    fgets(line,sizeof(line),in); //#-----------------------------------------------------
    fgets(line,sizeof(line),in); // Pt:
    
    fgets(line,sizeof(line),in); // lon :
    if( (nitems = sscanf(strstr(line," : "), "%s %lf", dummy, &(x[record]))) != 2){
      cout << "ERROR : can not read : "<< line << endl;
      return(status);
      }
    fgets(line,sizeof(line),in); // lat :
    if( (nitems = sscanf(strstr(line," : "), "%s %lf", dummy, &(y[record]))) != 2){
      cout << "ERROR : can not read : "<< line << endl;
      return(status);
      }
    fgets(line,sizeof(line),in); // nmes :
    if( (nitems = sscanf(strstr(line," : "), "%s %d", dummy, &(n[record]))) != 2){
      cout << "ERROR : can not read : "<< line << endl;
      return(status);
      }
    
    for(size_t m = 0; m < n[record]; ++m){
      fgets(line,sizeof(line),in);
      }

    }
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void XTRACK_ref_get_AllHeaders(const string & input, std::vector<int> id, double *x, double *y, int *nrecords, int &status)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  string format = XTRACK_get_fileFormat(input);

  // get header info
  if(format=="ASCII"){
     status = XTRACK_ref_get_AllHeader_ascii (input, id, x, y, nrecords);
    }
  if(format=="NETCDF"){
    vector<serie_t> track = ctoh_load_metadata_NetCDF(input);
    for(size_t n = 0; n < id.size(); n++) {
      x[n] = track[id[n]].data[0].lon;
      y[n] = track[id[n]].data[0].lat;
      nrecords[n] = track[id[n]].count;
      }
    track.clear();
    }

  format.clear();
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int XTRACK_get_meteo(serie_t *currentRecord)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i = 0, j = 0;

/*------------------------------------------------------------------------------
  check if regional DAC is available for all points */
  for(i = 0, j = 0; i <currentRecord->count; i++){
    if(!is_equal(currentRecord->data[j].values[4],currentRecord->mask,float(1e-3))) j++;
    }

/*------------------------------------------------------------------------------
  if regional DAC available for all points, use it */
  if(j == currentRecord->count){
    printf("use regional solutions for dealiasing atmospheric effects\n");
    return(X_TRACK_SERIE_T_DAC_REGIONAL);
    }
  else {
    printf("use global solutions for dealiasing atmospheric effects\n");
    return(X_TRACK_SERIE_T_DAC_GLOBAL);
    }

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int XTRACK_get_tide(serie_t *currentRecord)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i = 0, j = 0;

  for(i = 0, j = 0; i < currentRecord->count; i++){
    if(!is_equal(currentRecord->data[i].values[X_TRACK_SERIE_T_TIDE_REGIONAL],currentRecord->mask,float(1e-3))) j++;
    }
  
  if(j == currentRecord->count){
    printf("use regional solutions for dealiasing tides effects\n");
    return(X_TRACK_SERIE_T_TIDE_REGIONAL);
    }
  else {
    printf("use FES2004 solutions for dealiasing tides effects\n");
    return(X_TRACK_SERIE_T_TIDE_FES2004);
    }

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void XTRACK_ref_build_SLA(serie_t *currentRecord, int use_tide, int use_ib)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  float epsilon = float(1e-3);
  
  // automatic detection for dac correction
  float scale = 1.0;
  int DAC = X_TRACK_SERIE_T_DAC_GLOBAL;
  
  if(use_ib) {
//     printf ("warning-warning-warning-warning-warning-warning-warning-warning\n");
//     printf ("scale=%f , should be one except for Persian Gulf\n", scale);
    DAC = XTRACK_get_meteo(currentRecord);
    }

  // automatic detection for tide correction*/
  int TIDE = X_TRACK_SERIE_T_TIDE_FES2004;
  if(use_tide){
    TIDE = XTRACK_get_tide(currentRecord);
    }

  for(int i = 0; i < currentRecord->count; i++){
/*------------------------------------------------------------------------------
    SLA */
    if(!is_equal(currentRecord->data[i].values[X_TRACK_SERIE_T_SSHA],currentRecord->mask, epsilon)) {
      currentRecord->data[i].values[X_TRACK_SERIE_T_RESIDUAL] = currentRecord->data[i].values[X_TRACK_SERIE_T_SSHA];
      }
    else{
      currentRecord->data[i].values[X_TRACK_SERIE_T_RESIDUAL] = currentRecord->mask;
      continue;
      }

/*------------------------------------------------------------------------------
    meteo effects */
    if(use_ib){
      if(i==0) printf("eliminate DAC correction from ssha\n");
      if(!is_equal(currentRecord->data[i].values[DAC],currentRecord->mask, epsilon)) {
        currentRecord->data[i].values[X_TRACK_SERIE_T_RESIDUAL] -= currentRecord->data[i].values[DAC] * scale;
        }
      else {
        currentRecord->data[i].values[X_TRACK_SERIE_T_RESIDUAL] = currentRecord->mask;
        continue;
        }
      }

/*------------------------------------------------------------------------------
    LSA */
    if(!is_equal(currentRecord->data[i].values[X_TRACK_SERIE_T_DAC_LSA],currentRecord->mask, epsilon)) {
      currentRecord->data[i].values[X_TRACK_SERIE_T_RESIDUAL] -= currentRecord->data[i].values[X_TRACK_SERIE_T_DAC_LSA];
      }
    else {
      currentRecord->data[i].values[X_TRACK_SERIE_T_RESIDUAL] = currentRecord->mask;
      continue;
      }

/*------------------------------------------------------------------------------
    Solid Earth tide */
    if(!is_equal(currentRecord->data[i].values[X_TRACK_SERIE_T_DAC_SET],currentRecord->mask, epsilon)) {
      currentRecord->data[i].values[X_TRACK_SERIE_T_RESIDUAL] -= currentRecord->data[i].values[X_TRACK_SERIE_T_DAC_SET];
      }
    else {
      currentRecord->data[i].values[X_TRACK_SERIE_T_RESIDUAL] = currentRecord->mask;
      continue;
      }

/*------------------------------------------------------------------------------
    Ocean tide */
    if(use_tide){
      if(i==0) printf("eliminate OCEAN TIDE correction from ssha\n");
      if(!is_equal(currentRecord->data[i].values[TIDE],currentRecord->mask, epsilon)) {
        currentRecord->data[i].values[X_TRACK_SERIE_T_RESIDUAL] -= currentRecord->data[i].values[TIDE];
        }
      else {
        currentRecord->data[i].values[X_TRACK_SERIE_T_RESIDUAL] = currentRecord->mask;
        continue;
        }
      }
    //if(fabs(currentRecord->data[i].values[X_TRACK_SERIE_T_RESIDUAL9])>1.5){
    //  currentRecord->data[i].values[X_TRACK_SERIE_T_RESIDUAL]=currentRecord->mask;
    // }
    }

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void XTRACK_ref_write_header(FILE *out,int n)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{

  fprintf(out,"#-- HEADER -------------------------------------------\n");
  fprintf(out,"# Column 1 : date in days referred to\n");
  fprintf(out,"#            CNES date (01-JAN-1950 00:00:00.0)\n");
  fprintf(out,"# Column 2 : sea surface height anomaly (in meters)\n");
  fprintf(out,"# Column 3 : MOG2D-G model sea level (in meters)\n");
  fprintf(out,"# Column 4 : MOG2D-G model S1 atmospheric tide (in meters)\n");
  fprintf(out,"# Column 5 : MOG2D-G model S2 atmospheric tide (in meters)\n");
  fprintf(out,"# Column 6 : MOG2D regional model sea level (in meters)\n");
  fprintf(out,"# Column 7 : MOG2D regional model S1 atmospheric tide (in meters)\n");
  fprintf(out,"# Column 8 : MOG2D regional model S2 atmospheric tide (in meters)\n");
  fprintf(out,"# Column 9 : inverted barometer (in meters)\n");
  fprintf(out,"# Column 10: inverted barometer S1 atmospheric tide (in meters)\n");
  fprintf(out,"# Column 11: inverted barometer S2 atmospheric tide (in meters)\n");
  fprintf(out,"# Column 12: loading effects (in meters)\n");
  fprintf(out,"# Column 13: solid Earth tide (in meters)\n");
  fprintf(out,"# Column 14: MOG2D regional model ocean tide (in meters)\n");
  fprintf(out,"# Column 15: MOG2D regional model S1 ocean tide (in meters)\n");
  fprintf(out,"# Column 16: MOG2D regional model S2 ocean tide (in meters)\n");
  fprintf(out,"# Column 17: GOT4.7 model ocean tide (in meters)\n");
  fprintf(out,"# Column 18: GOT4.7 model S1 ocean tide (in meters)\n");
  fprintf(out,"# Column 19: GOT4.7 model S2 ocean tide (in meters)\n");
  fprintf(out,"# Column 20: FES2004 model ocean tide (in meters)\n");
  fprintf(out,"# Column 21: FES2004 model S1 ocean tide (in meters)\n");
  fprintf(out,"# Column 22: FES2004 model S2 ocean tide (in meters)\n");
  fprintf(out,"# Column 23: harmonic analysis ocean tide (in meters)\n");
  fprintf(out,"# Column 24: harmonic analysis S1 ocean tide (in meters)\n");
  fprintf(out,"# Column 25: harmonic analysis S2 ocean tide (in meters)\n");
  fprintf(out,"#-- HEADER END ---------------------------------------\n");
  fprintf(out,"Number of points :\n");
  fprintf(out,"%d\n",n);
  fprintf(out,"#-----------------------------------------------------\n");

  fflush(out);

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void XTRACK_WriteAsciiRecords(FILE *stream, serie_t *record, int *index, size_t size, const char *rootname)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  FILE *gnu;
  char gnufile[256];
  
  for(size_t count = 0; count < size ; count++){
    fprintf(stream,"Pt  : %d\n",index[count]);
    fprintf(stream,"lon : %lf\n",record[count].data[0].lon);
    fprintf(stream,"lat : %lf\n",record[count].data[0].lat);
    fprintf(stream,"Mes : %d\n",record[count].count);

    for(int i = 0; i < record[count].count; i++){
      fprintf(stream,"%12.6f",record[count].data[i].time);
      for(int l = 0; l < X_TRACK_SERIE_T_NUMVAL; l++){
        fprintf(stream," %7.4f",record[count].data[i].values[l]);
        }
      fprintf(stream,"\n");
      }
    fprintf(stream,"#-----------------------------------------------------\n");
    
    if(rootname==0) continue;
    
    sprintf(gnufile,"%s.%3.3d.gnu",rootname,index[count]);
    gnu=fopen(gnufile,"w");
    for(int i = 0; i < record[count].count; i++){
//       if(is_equal(record[count].data[i].values[X_TRACK_SERIE_T_RESIDUAL],record[count].mask,float(1e-3)) ) continue;
      if(is_equal(record[count].data[i].values[X_TRACK_SERIE_T_RESIDUAL],record[count].mask,float(1e-3)) ){
        continue;
        fprintf(gnu,"\n");
        }
      fprintf(gnu,"%12.6f",record[count].data[i].time);
      for(int l = 0; l < X_TRACK_SERIE_T_NUMVAL; l++){
        fprintf(gnu," %7.4f",record[count].data[i].values[l]);
        }
      fprintf(gnu,"\n");
      }
    fclose(gnu);
    }

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int XTRACK_CheckRecords(serie_t & record, int index, const char *rootname)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
//   FILE *gnu;
//   char gnufile[256];
  vector<double> time;
  vector<int>  clusters;
  vector<int>  *nodes;
  bool agregate;
  
//   sprintf(gnufile,"%s.%3.3d.gnu",rootname,index);
  
  nodes=new vector<int> [record.count];
  
//   gnu=fopen(gnufile,"w");
  for(int i = 0; i < record.count; i++){
    if(is_equal(record.data[i].values[X_TRACK_SERIE_T_SSHA],record.mask,float(1e-3)) ) continue;
    double t=record.data[i].time;
    agregate=false;
    if(time.size()==0) {
      nodes[0].push_back(i);
      time.push_back(t);
      continue;
      }
    for(int l = 0; l < time.size(); l++){
      if(fabs(t-time[l])<1.e-03) {
        nodes[l].push_back(i);
        agregate=true;
        break;
        }
      }
    if(!agregate) {
      time.push_back(t);
      int l = time.size()-1;
      nodes[l].push_back(i);
      }
//     fprintf(gnu,"\n");
    }
  
//   fclose(gnu);
  for(int l = 0; l < time.size(); l++){
    if(nodes[l].size()==1) continue;
    int n=nodes[l][0];
    for(int k = 1; k < nodes[l].size(); k++){
      int m=nodes[l][k];
      for(int i = 0; i < X_TRACK_SERIE_T_NUMVAL; i++) {
        record.data[n].values[i]+=record.data[m].values[i];
        record.data[m].values[i]=record.mask;
        }
      }
    for(int i = 0; i < X_TRACK_SERIE_T_NUMVAL; i++) record.data[n].values[i]/=nodes[l].size();
    }
  
  return(0);

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 int ctoh_get_variable(int nc_id, int var_id, int size, int **buffer, int *scale, int *offset, int *mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int n,status;
  
  status = nc_get_att_int(nc_id, var_id, "scale_factor", scale);
  if(status != NC_NOERR){
    *scale = 1;
    }
  
  status = nc_get_att_int(nc_id, var_id, "offset", offset);
  if(status != NC_NOERR){
    *offset = 0;
    }

  status = nc_get_att_int(nc_id, var_id, "_FillValue", mask);
  if(status != NC_NOERR){
//    NC_TRAP_ERROR(wexit,status,1, "nc_get_att_float() failed");
    *mask=-1;
    }
  
  *buffer = new int[size];
  status = nc_get_var_int(nc_id, var_id, *buffer);
  NC_TRAP_ERROR(wexit,status,1, "nc_get_var_float() failed");
  for(n=0;n<size;n++) {
    if((*buffer)[n]!=*mask) (*buffer)[n]=(*buffer)[n] * (*scale) + *offset;
    }
  status=0;
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 int ctoh_get_variable(int nc_id, int var_id, int size, float **buffer, float *scale, float *offset, float *mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int n,status;
  
  status = nc_get_att_float(nc_id, var_id, "scale_factor", scale);
  if(status != NC_NOERR){
    *scale = 1;
    }
  
  status = nc_get_att_float(nc_id, var_id, "offset", offset);
  if(status != NC_NOERR){
    *offset = 0;
    }

  status = nc_get_att_float(nc_id, var_id, "_FillValue", mask);
  NC_CHKERR_BASE_LINE(status,"nc_get_att_float(,,_FillValue,) failed");
  if(status!=0)*mask=NC_FILL_FLOAT;
  
  *buffer = new float[size];
  status = nc_get_var_float(nc_id, var_id, *buffer);
  NC_TRAP_ERROR(wexit,status,1, "nc_get_var_float() failed");
  for(n=0;n<size;n++) {
    if((*buffer)[n]!=*mask) (*buffer)[n]=(*buffer)[n] * (*scale) + *offset;
    }
  status=0;
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 int ctoh_get_variable(int nc_id, int var_id, int size, double **buffer, double *scale, double *offset, double *mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int n,status;
  
  status = nc_get_att_double(nc_id, var_id, "scale_factor", scale);
  if(status != NC_NOERR){
    *scale = 1;
    }
  
  status = nc_get_att_double(nc_id, var_id, "offset", offset);
  if(status != NC_NOERR){
    *offset = 0;
    }

  status = nc_get_att_double(nc_id, var_id, "_FillValue", mask);
  NC_TRAP_ERROR(wexit,status,1, "nc_get_att_float() failed");
  
  *buffer = new double[size];
  status = nc_get_var_double(nc_id, var_id, *buffer);
  NC_TRAP_ERROR(wexit,status,1, "nc_get_var_float() failed");
  for(n=0;n<size;n++) {
    if((*buffer)[n]!=*mask) (*buffer)[n]=(*buffer)[n] * (*scale) + *offset;
    }
  status=0;
  return(status);
}


#define ctoh_round(a) (floor(a * 10000. +0.5) / 10000.)
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  vector<serie_t> ctoh_load_metadata_NetCDF(const string & filename)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status = NC_NOERR + 1;

  int nc_id;
  status = nc_open(filename.c_str(), NC_NOWRITE, &nc_id);
  if(status!=0) NC_TRAP_ERROR(wexit,status,1,"nc_open(\""+filename+"\",NC_NOWRITE,) error");

/* *----------------------------------------------------------------------------
  get dimensions */
  int dim_cycle_id;
  status = nc_inq_dimid(nc_id, "nbcycles", &dim_cycle_id);
  NC_TRAP_ERROR(wexit,status,1, "nc_inq_dimid() failed");
  size_t i_CYCLES_dim;
  status = nc_inq_dimlen(nc_id, dim_cycle_id, &i_CYCLES_dim);
  NC_TRAP_ERROR(wexit,status,1, "nc_inq_dimlen() failed");

  int dim_points_id;
  status = nc_inq_dimid(nc_id, "nbpoints", &dim_points_id);
  NC_TRAP_ERROR(wexit,status,1, "nc_inq_dimid() failed");
  size_t i_POINTS_dim;
  status = nc_inq_dimlen(nc_id, dim_points_id, &i_POINTS_dim);
  NC_TRAP_ERROR(wexit,status,1, "nc_inq_dimlen() failed");

/* *----------------------------------------------------------------------------
  Laurent, il faudra un jour que tu m'explique pourquoi tu utilises des types qui
  te compliquent la mise au point !!!
  Tes codes sont deja assez confus comme ca... */
  vector<serie_t> track(i_POINTS_dim);
  for(size_t i = 0; i < i_POINTS_dim; ++i){
    track[i].data = new data_t[i_CYCLES_dim];
    track[i].count = i_CYCLES_dim;
    track[i].mask = 9999.f;
    }
  
  float epsilon = 1.e-3f;
  float scale;
  float offset;
  float fillValue;

  double dscale = 1;
  double doffset = 0;
  double dfillValue = 0;
  int iscale, ioffset, ifillValue;
  size_t size;
  
  int *cycles;
  float *lon,*lat,*tide=0,*sla=0,*dac=0;
  double *time;
  
/**----------------------------------------------------------------------------
  get longitudes */
  int lon_id;
  status = nc_inq_varid(nc_id, "lon", &lon_id);
  NC_TRAP_ERROR(wexit,status,1, "nc_inq_varid() failed");
  size=i_POINTS_dim;
  status=ctoh_get_variable(nc_id,lon_id,size,&lon,&scale,&offset,&fillValue);
  NC_TRAP_ERROR(wexit,status,1, "nc_get_var_float() failed");

  for(size_t i = 0; i < track.size(); ++i){
    bool valid=!is_equal(lon[i], fillValue, epsilon);
    for(size_t j = 0; j < track[i].count; ++j){
      if(valid ) {
        track[i].data[j].lon = lon[i];
        }
      else{
        track[i].data[j].lon = track[i].mask;
        }
      }
    }
  delete [] lon;


/**----------------------------------------------------------------------------
  get latitudes */
  int lat_id;
  status = nc_inq_varid(nc_id,"lat",&lat_id);
  NC_TRAP_ERROR(wexit,status,1, "nc_inq_varid() failed");
  size=i_POINTS_dim;
  status=ctoh_get_variable(nc_id,lat_id,size,&lat,&scale,&offset,&fillValue);
  NC_TRAP_ERROR(wexit,status,1, "nc_get_var_float() failed");
  for(size_t i = 0; i < track.size(); ++i){
    bool valid=!is_equal(lat[i], fillValue, epsilon);
    for(size_t j = 0; j < track[i].count; ++j){
      if(valid ) {
        track[i].data[j].lat = lat[i];
        }
      else{
        track[i].data[j].lat = track[i].mask;
        }
      }
    }
  delete [] lat;

 
/**----------------------------------------------------------------------------
  get time */
  int time_id;
  status = nc_inq_varid(nc_id,"time",&time_id);
  NC_TRAP_ERROR(wexit,status,1, "nc_inq_varid() failed");
  size=i_POINTS_dim * i_CYCLES_dim;
  status=ctoh_get_variable(nc_id,time_id,size,&time,&dscale,&doffset,&dfillValue);
  NC_TRAP_ERROR(wexit,status,1, "nc_get_var_float() failed");
  for(size_t i = 0; i < track.size(); ++i){
    for(size_t j = 0; j < track[i].count; ++j){
      size_t m=i * i_CYCLES_dim + j;
      bool valid=!is_equal(time[m], dfillValue, (double) epsilon);
      if(valid){
        track[i].data[j].time = time[m];
        }
      else{
        track[i].data[j].time = track[i].mask;
        }
      }
    }
  delete [] time;

/**----------------------------------------------------------------------------
  get cycles */
  int cycle_id;
  status = nc_inq_varid(nc_id,"cycle",&cycle_id);
  if( status != NC_NOERR){
    status = nc_inq_varid(nc_id,"cycles",&cycle_id);
    }
  NC_CHKERR_BASE_LINE(status,"nc_inq_varid((\""+filename+"\"),[\"cycle\"|\"cycles\"],) error");
  
  if( status != NC_NOERR){
    for(size_t i = 0; i < track.size(); ++i){
      for(size_t j = 0; j < track[i].count; ++j){
        track[i].data[j].cycle = j;
        }
      }
    }
  else {
    size=i_CYCLES_dim;
    status=ctoh_get_variable(nc_id,cycle_id,size,&cycles,&iscale,&ioffset,&ifillValue);
    NC_TRAP_ERROR(wexit,status,1, "ctoh_get_variable() failed");
  
    for(size_t i = 0; i < track.size(); ++i){
      for(size_t j = 0; j < track[i].count; ++j){
        size_t m=j;
        bool valid=!is_equal(cycles[m], ifillValue, (int) 0);
        if(valid){
          track[i].data[j].cycle = cycles[m];
          }
        else{
          track[i].data[j].cycle = track[i].mask;
          }
        }
      }
    delete [] cycles;
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  get tides
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  int tide_id;
  status = nc_inq_varid(nc_id,"tide",&tide_id);
  if( status != NC_NOERR){
    status = nc_inq_varid(nc_id,"tidefilt",&tide_id);
    }
  NC_CHKERR_BASE_LINE(status,"nc_inq_varid((\""+filename+"\"),[\"tide\"|\"tidefilt\"],) error");
  
  if(status==0){
    size=i_POINTS_dim * i_CYCLES_dim;
    status=ctoh_get_variable(nc_id,tide_id,size,&tide,&scale,&offset,&fillValue);
    NC_TRAP_ERROR(wexit,status,1, "ctoh_get_variable() failed");
    }
  
  for(size_t i = 0; i < track.size(); ++i){
    for(size_t j = 0; j < track[i].count; ++j){
      if(tide==0){
        track[i].data[j].values[X_TRACK_SERIE_T_TIDE_REGIONAL] = 0.f;
        track[i].data[j].values[X_TRACK_SERIE_T_TIDE_GOT4_7]   = 0.f;
        track[i].data[j].values[X_TRACK_SERIE_T_TIDE_FES2004]  = 0.f;
        continue;
        }
      size_t m=i * i_CYCLES_dim + j;
      bool valid=!is_equal(tide[m], fillValue, epsilon);
      if(valid ) {
        if(fabs(tide[m])>10.) {
          printf("trouble\n");
          }
        track[i].data[j].values[X_TRACK_SERIE_T_TIDE_REGIONAL] = ctoh_round(tide[m]);
        track[i].data[j].values[X_TRACK_SERIE_T_TIDE_GOT4_7]   = ctoh_round(tide[m]);
        track[i].data[j].values[X_TRACK_SERIE_T_TIDE_FES2004]  = ctoh_round(tide[m]);
        }
      else{
        track[i].data[j].values[X_TRACK_SERIE_T_TIDE_REGIONAL] = track[i].mask;
        track[i].data[j].values[X_TRACK_SERIE_T_TIDE_GOT4_7]   = track[i].mask;
        track[i].data[j].values[X_TRACK_SERIE_T_TIDE_FES2004]  = track[i].mask;
        }
      }
    }
  deletep(&tide);

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  get DAC
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  int dac_id;
  status = nc_inq_varid(nc_id,"mog2d",&dac_id);
  if( status != NC_NOERR) status = nc_inq_varid(nc_id,"dac",&dac_id);
  
  if( status != NC_NOERR){
    for(size_t i = 0; i < track.size(); ++i){
      for(size_t j = 0; j < track[i].count; ++j){
//         track[i].data[j].values[X_TRACK_SERIE_T_DAC_GLOBAL]   = track[i].mask;
//         track[i].data[j].values[X_TRACK_SERIE_T_DAC_REGIONAL] = track[i].mask;
        track[i].data[j].values[X_TRACK_SERIE_T_DAC_GLOBAL]   = 0;
        track[i].data[j].values[X_TRACK_SERIE_T_DAC_REGIONAL] = 0;
        }
      }
    }
  else {
    size=i_POINTS_dim * i_CYCLES_dim;
    status=ctoh_get_variable(nc_id,dac_id,size,&dac,&scale,&offset,&fillValue);
    NC_TRAP_ERROR(wexit,status,1, "ctoh_get_variable() failed");
    
    for(size_t i = 0; i < track.size(); ++i){
      for(size_t j = 0; j < track[i].count; ++j){
        size_t m=i * i_CYCLES_dim + j;
        bool valid=!is_equal(dac[m], fillValue, epsilon);
        if(valid ) {
          track[i].data[j].values[X_TRACK_SERIE_T_DAC_GLOBAL]   = ctoh_round(dac[m]);
          track[i].data[j].values[X_TRACK_SERIE_T_DAC_REGIONAL] = ctoh_round(dac[m]);
          }
        else{
          track[i].data[j].values[X_TRACK_SERIE_T_DAC_GLOBAL]   = track[i].mask;
          track[i].data[j].values[X_TRACK_SERIE_T_DAC_REGIONAL] = track[i].mask;
          }
        }
      }
    delete [] dac;
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  get SLA
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  int sla_id;
  status = nc_inq_varid(nc_id,"sla",&sla_id);
  if( status != NC_NOERR){
    status = nc_inq_varid(nc_id,"slafilt",&sla_id);
    }
  NC_TRAP_ERROR(wexit,status,1, "nc_inq_varid() failed on sla variable\n");
  size=i_POINTS_dim * i_CYCLES_dim;
  status=ctoh_get_variable(nc_id,sla_id,size,&sla,&scale,&offset,&fillValue);
  NC_TRAP_ERROR(wexit,status,1, "ctoh_get_variable() failed");
  for(size_t i = 0; i < track.size(); ++i){
    for(size_t j = 0; j < track[i].count; ++j) {
      size_t m=i * i_CYCLES_dim + j;
      bool valid=!is_equal(sla[m], fillValue, epsilon);
      if(valid ) {
        track[i].data[j].values[X_TRACK_SERIE_T_SSHA] = ctoh_round(sla[m]); // caution tides and DAC already applied
        }
      else{
        track[i].data[j].values[X_TRACK_SERIE_T_SSHA] = track[i].mask;
        }
      track[i].data[j].values[X_TRACK_SERIE_T_RESIDUAL] = track[i].data[j].values[X_TRACK_SERIE_T_SSHA];
      }
    }
  delete [] sla;

  status = nc_close(nc_id);
  NC_TRAP_ERROR(wexit,status,1, "nc_close() failed");

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  accumulate data with valid values
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  for(size_t r = 0; r < track.size(); ++r) {
    serie_t *trackr=&track[r];
    
    size_t n = trackr->count;
    data_t* tmp = new data_t[n]; //over alloc
    size_t count = 0, time_count=0,ssha_count=0,tide_count=0;
    for(size_t i = 0; i < n; i++){
      float mask  = trackr->mask;
      double time = trackr->data[i].time;
      float ssha  = trackr->data[i].values[X_TRACK_SERIE_T_SSHA];
      float tide  = trackr->data[i].values[X_TRACK_SERIE_T_TIDE_GOT4_7];
      bool time_is_valid = isfinite(time) and not is_equal(time, (double)mask, 1e-3);
      bool ssha_is_valid = not is_equal(ssha, mask, 1e-3f);
      bool tide_is_valid = not is_equal(tide, mask, 1e-3f);
      if(ssha_is_valid != tide_is_valid) {
//        printf("time serie anomaly record %d, point %d: ssh flag=%d tide flag=%d\n", r, i, ssha_is_valid, tide_is_valid);
        }
      if(time_is_valid) time_count++;
      if(ssha_is_valid) ssha_count++;
      if(tide_is_valid) tide_count++;
      if(time_is_valid && ssha_is_valid && tide_is_valid) {
        tmp[count] = trackr->data[i];
        tmp[count].cycle = i;
        count++;
        }
      }
    trackr->count = count;
    if(count == 0){
      STDOUT_BASE_LINE("track[%u]: %u time, %u SSHA and %u tide out of %u\n",r,time_count,ssha_count,tide_count,n);
      delete[] tmp;
      continue;
      }
    delete [] trackr->data;
    trackr->data = new data_t[count];
    for(size_t i = 0; i < count; ++i){
      trackr->data[i] = tmp[i];
      }
    delete[] tmp;
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  set back tides and dac, already applied in SLA
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  for(size_t r = 0; r < track.size(); ++r) {
    if(track[r].count==0) continue;
//    printf("track %d: treating %d valid data\n",r,track[r].count);
    for(size_t i = 0; i < track[r].count; ++i){
      float ssha =track[r].data[i].values[X_TRACK_SERIE_T_SSHA];
      float dac  =track[r].data[i].values[X_TRACK_SERIE_T_DAC_GLOBAL];
      float tide =track[r].data[i].values[X_TRACK_SERIE_T_TIDE_GOT4_7];
      if(fabs(ssha)>10.) {
        printf("trouble on ssha\n");
        }
      if(fabs(tide)>10.) {
        printf("trouble on tide\n");
        }
      if(fabs(dac)>10.) {
        printf("trouble on dac\n");
        }
      ssha+=dac+tide;
      track[r].data[i].values[X_TRACK_SERIE_T_SSHA] =ssha;
      }
   }
  return(track);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 int ctoh_put_variable(const char *filename, const char *name, vector<poc_dim_t> dim, float *buffer, float scale, float offset, float mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int id,ncid,status;
  poc_var_t var(name, NC_FLOAT);
  
  status = nc_open(filename, NC_WRITE, &ncid);
  NC_TRAP_ERROR(wexit,status,1,"nc_open() failed");
  
  status=nc_redef(ncid);
  NC_TRAP_ERROR(wexit,status,1,"nc_open() failed");
  
  var.attributes << poc_att_t("_FillValue",mask);
  var.attributes << poc_att_t("scale",scale);
  var.attributes << poc_att_t("offset",offset);
    
  for(int k=0;k<dim.size();k++) {
    var.dimensions << dim[k];
    }
  
  status=poc_def_var(ncid,var,&id);
  
  status=nc_enddef(ncid);
    
  status=nc_put_var_float(ncid, id, buffer);
  
  status=nc_close(ncid);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 int ctoh_put_variable(const char *filename, const char *name, vector<poc_dim_t> dim, double *buffer, double scale, double offset, double mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int id,ncid,status;
  poc_var_t var(name, NC_DOUBLE);
  
  status = nc_open(filename, NC_WRITE, &ncid);
  NC_TRAP_ERROR(wexit,status,1,"nc_open() failed");
  
  status=nc_redef(ncid);
  NC_TRAP_ERROR(wexit,status,1,"nc_open() failed");
  
  var.attributes << poc_att_t("_FillValue",mask);
  var.attributes << poc_att_t("scale",scale);
  var.attributes << poc_att_t("offset",offset);
  
  var.name=name;
  
  for(int k=0;k<dim.size();k++) {
    var.dimensions << dim[k];
    }
  
  status=poc_def_var(ncid, var, &id);
  
  status=nc_enddef(ncid);
  
  status=nc_put_var_double(ncid, id, buffer);
    
  status=nc_close(ncid);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int ctoh_save_metadata_NetCDF(const string & filename, vector<serie_t> & track)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status = NC_NOERR + 1;
  int ncid;
 
  int dimid;

  poc_dim_t ncycles, npoints;
  
  vector<poc_dim_t> dim1D, dim2D;

  float scale = 1;
  float offset = 0;
  float fillValue=track[0].mask;

  double dscale = 1;
  double doffset = 0;
  double dfillValue = 0;
  
  
  status=nc_create(filename.c_str(), NC_CLOBBER, &ncid);
  NC_TRAP_ERROR(wexit,status,1,"nc_open() failed");
  
  ncycles.len  = 0;
  for(size_t k=0;k<track.size();k++) {
    updatemax(&ncycles.len,track[k].count);
    }
  ncycles.name = "nbcycles";
  
  status=poc_def_dim(ncid,(const poc_dim_t) ncycles, &dimid);
  
  npoints.len  = track.size();
  npoints.name = "nbpoints";
  
  status=poc_def_dim(ncid,(const poc_dim_t) npoints, &dimid);
  
  status=nc_close(ncid);
  
  dim1D.push_back(npoints);
  
  dim2D.push_back(ncycles);
  dim2D.push_back(npoints);
 
  float  *buffer = new float[ncycles.len];
  double *time   = new double[ncycles.len*npoints.len];
  
/**----------------------------------------------------------------------------
  put longitudes */
  for(size_t i = 0; i < track.size(); ++i){
    buffer[i]=track[i].data[0].lon;
    }
  status=ctoh_put_variable(filename.c_str(),"lon",dim1D,buffer,scale,offset,fillValue);
  NC_TRAP_ERROR(wexit,status,1, "nc_get_var_float() failed");

/**----------------------------------------------------------------------------
  put latitudes */
  for(size_t i = 0; i < track.size(); ++i){
    buffer[i]=track[i].data[0].lat;
    }
  status=ctoh_put_variable(filename.c_str(),"lat",dim1D,buffer,scale,offset,fillValue);
  NC_TRAP_ERROR(wexit,status,1, "nc_get_var_float() failed");
  
  delete [] buffer;

/**----------------------------------------------------------------------------
  put time */
  for(size_t i = 0; i < track.size(); ++i){
    for(size_t j = 0; j < track[i].count; ++j){
      size_t m=i * ncycles.len + j;
      time[m]=track[i].data[j].time;
      }
    }
  status=ctoh_put_variable(filename.c_str(),"time",dim2D,time,dscale,doffset,dfillValue);
  NC_TRAP_ERROR(wexit,status,1, "nc_get_var_float() failed");
  delete [] time;

/**----------------------------------------------------------------------------
  put sla */
  buffer  = new float[ncycles.len*npoints.len];;
//   for(size_t i = 0; i < track.size(); ++i){
//     for(size_t j = 0; j < track[i].count; ++j){
//       size_t m=i * ncycles.len + j;
//       buffer[m]=track[i].data[j].values[X_TRACK_SERIE_T_SSHA];
//       }
//     }
//   status=ctoh_put_variable(filename.c_str(),"sla_0",dim2D,buffer,scale,offset,fillValue);
//   if( status != 0){
//     NC_TRAP_ERROR(wexit,status,1, "nc_get_var_float() failed");
//     }
  
/**----------------------------------------------------------------------------
  put dac */
  for(size_t i = 0; i < track.size(); ++i){
    for(size_t j = 0; j < track[i].count; ++j){
      size_t m=i * ncycles.len + j;
      buffer[m]=track[i].data[j].values[X_TRACK_SERIE_T_DAC_GLOBAL];
      }
    }
  status=ctoh_put_variable(filename.c_str(),"mog2d",dim2D,buffer,scale,offset,fillValue);
  NC_TRAP_ERROR(wexit,status,1, "nc_get_var_float() failed");
  
/**----------------------------------------------------------------------------
  put tidal prediction (from analysis)*/
  for(size_t i = 0; i < track.size(); ++i){
    for(size_t j = 0; j < track[i].count; ++j){
      size_t m=i * ncycles.len + j;
      buffer[m]=track[i].data[j].values[X_TRACK_SERIE_T_TIDE_ANALYSIS];
      }
    }
  status=ctoh_put_variable(filename.c_str(),"tide",dim2D,buffer,scale,offset,fillValue);
  NC_TRAP_ERROR(wexit,status,1, "nc_get_var_float() failed");
  
/**----------------------------------------------------------------------------
  put residuals (from analysis)*/
  buffer  = new float[ncycles.len*npoints.len];;
  for(size_t i = 0; i < track.size(); ++i){
    for(size_t j = 0; j < track[i].count; ++j){
      size_t m=i * ncycles.len + j;
      buffer[m]=track[i].data[j].values[X_TRACK_SERIE_T_RESIDUAL];
      }
    }
  status=ctoh_put_variable(filename.c_str(),"sla",dim2D,buffer,scale,offset,fillValue);
  NC_TRAP_ERROR(wexit,status,1, "nc_get_var_float() failed");
  
  delete [] buffer;

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 int poc_register_variable(const char *filename, poc_var_t var)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int id,ncid,status;
  
  status = nc_open(filename, NC_WRITE, &ncid);
  NC_TRAP_ERROR(wexit,status,1,"nc_open() failed");
  
  status=nc_redef(ncid);
  NC_TRAP_ERROR(wexit,status,1,"nc_open() failed");
    
  status=poc_def_var(ncid,var,&id);
  
  status=nc_enddef(ncid);
  
  status=nc_close(ncid);
  
  return(status);
}



/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int ctoh_save_metadata_NetCDF(const string & filename,const string & source, vector<serie_t> & track)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status = NC_NOERR + 1;
//   int ncid;
 
//   int dimid;

  poc_dim_t ncycles, npoints;
  
  poc_list_t<poc_dim_t> dim1D, dim2D;
  
  poc_var_t vlon, vlat, vtime;
  poc_var_t vmssh, vcycle, vpoint;;
  poc_var_t vsla, vdac, vtide;
  poc_global_t info;
  
  status=poc_inq(source,&info);
 
  
  status=poc_inq_var(source,"lon",  &vlon);
  status=poc_inq_var(source,"lat",  &vlat);
  status=poc_inq_var(source,"time", &vtime);
  status=poc_inq_var(source,"mssh", &vmssh);
  status=poc_inq_var(source,"cycle", &vcycle);
  status=poc_inq_var(source,"point", &vpoint);

  status=poc_create(filename.c_str(), info, 1);
  NC_TRAP_ERROR(wexit,status,1,"nc_open() failed");
  
  ncycles=vtime.dimensions[1];
  npoints=vtime.dimensions[0];
  
//   status=nc_create(filename.c_str(), NC_CLOBBER, &ncid);
//   status=poc_def_dim(ncid,(const poc_dim_t) ncycles, &dimid);
//   
//   status=poc_def_dim(ncid,(const poc_dim_t) npoints, &dimid);
//   
//   status=nc_close(ncid);
  
  dim1D.push_back(npoints);
  
  dim2D.push_back(ncycles);
  dim2D.push_back(npoints);
 
/**----------------------------------------------------------------------------
  put longitudes */
  poc_data_t<float> lon;
  
  lon.init(vlon,"");
  status=lon.read_data(source);
  
  status=lon.write_data(filename);

/**----------------------------------------------------------------------------
  put latitudes */
  poc_data_t<float> lat;
  
  lat.init(vlat,"");
  status=lat.read_data(source);
  
  status=lat.write_data(filename);

/**----------------------------------------------------------------------------
  put time */
  poc_data_t<float> time;

  time.init(vtime,"");
  status=time.read_data(source);

  status=time.write_data(filename);

/**----------------------------------------------------------------------------
  put mssh */
  poc_data_t<float> mssh;
  
  mssh.init(vmssh,"");
  status=mssh.read_data(source);
  
  if(status==0)
    status=mssh.write_data(filename);

/**----------------------------------------------------------------------------
  put cycle */
  poc_data_t<float> cycle;
  
  cycle.init(vcycle,"");
  status=cycle.read_data(source);
  
  if(status==0)
    status=cycle.write_data(filename);

/**----------------------------------------------------------------------------
  put point */
  poc_data_t<float> point;
  
  point.init(vpoint,"");
  status=point.read_data(source);
  
  if(status==0)
    status=point.write_data(filename);

/**----------------------------------------------------------------------------
  put sla */
  poc_data_t<float> sla;
  
  status=poc_inq_var(source, "sla",  &vsla);
  
  if(status==0){
    sla.init(vsla,"");
    
    for(size_t i = 0; i < npoints.len; ++i){
      for(size_t j = 0; j < ncycles.len; ++j){
        size_t m=i * ncycles.len + j;
        sla.data[m]=sla.mask;
        }
      for(size_t j = 0; j < track[i].count; ++j){
        size_t m=i * ncycles.len + track[i].data[j].cycle;
        sla.data[m]=track[i].data[j].values[X_TRACK_SERIE_T_RESIDUAL];
        }
      }
    status=sla.write_data(filename);
    }

/**----------------------------------------------------------------------------
  put std */
  poc_data_t<float> std;
  std.info.init("std_sla", NC_FLOAT, "", "", sla.mask);
  std.info.dimensions=dim1D;
  std.info.attributes
    <<poc_att_t("scale",sla.scale)
    <<poc_att_t("offset",sla.offset);
  std.init("");
  
  for(size_t i = 0; i < npoints.len; ++i) {
    float *buffer=new float[ncycles.len];
    statistic_t s;
    for(size_t j = 0; j < ncycles.len; ++j) {
      buffer[j]=std.mask;
      }
    for(size_t j = 0; j < track[i].count; ++j) {
      size_t m=track[i].data[j].cycle;
      buffer[m]=track[i].data[j].values[X_TRACK_SERIE_T_RESIDUAL];
      }
    s=get_statistics(buffer, std.mask, ncycles.len);
    std.data[i]=s.std;
    delete[] buffer;
    }
  status=std.write_data(filename);
  
/**----------------------------------------------------------------------------
  put dac */
  poc_data_t<float> dac;
  
  status=poc_inq_var(source, "dac",  &vdac);
  
  if(status==0){
    dac.init(vdac,"");
    
    for(size_t i = 0; i < npoints.len; ++i) {
      for(size_t j = 0; j < ncycles.len; ++j) {
        size_t m=i * ncycles.len + j;
        dac.data[m]=dac.mask;
        }
      for(size_t j = 0; j < track[i].count; ++j){
        size_t m=i * ncycles.len + track[i].data[j].cycle;
        dac.data[m]=track[i].data[j].values[X_TRACK_SERIE_T_DAC_GLOBAL];
        }
      }
    
    status=dac.write_data(filename);
    }
  
/**----------------------------------------------------------------------------
  put tidal prediction (from analysis)*/
  poc_data_t<float> tide;
  
  status=poc_inq_var(source, "tide",  &vtide);
  
  if(status==0){
    vtide.attributes << poc_att_t("long_name","harmonic analysis");
    
    tide.init(vtide,"");
    
    for(size_t i = 0; i < npoints.len; ++i) {
      for(size_t j = 0; j < ncycles.len; ++j) {
        size_t m=i * ncycles.len + j;
        tide.data[m]=tide.mask;
        }
      for(size_t j = 0; j < track[i].count; ++j) {
        size_t m=i * ncycles.len + track[i].data[j].cycle;
        tide.data[m]=track[i].data[j].values[X_TRACK_SERIE_T_TIDE_ANALYSIS];
        }
      }
    
    status=tide.write_data(filename);
    }
  
/**----------------------------------------------------------------------------
  put residuals (from analysis)*/

  return(status);
}
