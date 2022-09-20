/**************************************************************************

  X-TRACK tools, 2006-2009

Contributors:

  Florent Lyard      LEGOS/CNRS, Toulouse, France
  Yoann Le Bars      LEGOS, Toulouse, France (PhD)
  Laurent Roblou     LEGOS, Toulouse, France
  Anne Vernier       NOVELTIS, Ramonville-Saint-Agne, France

E-mail: florent.lyard@legos.obs-mip.fr

***************************************************************************/
#include <config.h>

#include <iostream>
#include <vector>

#include "tools-structures.h"

#include "tides.h"
#include "fe.h"
#include "archive.h"
#include "poc-time.h"
#include "mgr.h"
#include "functions.h"
#include "xtrack-io.h"

#ifndef VERBOSE
#define VERBOSE 1
#endif


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void free_liste_ctoh(track_data *liste, int nb_traces, int nb_ondes)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,k;

  for (k=0;k<nb_traces;k++) { //boucle sur le nombre de traces contenues dans le tableau
    free(liste[k].lon);
    free(liste[k].lat);
    free(liste[k].mssh);
    free(liste[k].cycle);
    free(liste[k].point);
    free(liste[k].lat_ah);
    free(liste[k].lon_ah);
    for(i=0;i<liste[k].nbpoints;i++) { //boucle sur le nombre de points
      free(liste[k].time[i]);
    }
    free(liste[k].time);
    for(i=0;i<liste[k].nbpoints;i++) {//boucle sur le nombre de points
      free(liste[k].sla[i]);
    }
    free(liste[k].sla);

    for(i=0;i<liste[k].nbpoints;i++) {//boucle sur le nombre de points
      free(liste[k].tide[i]);
    }
    free(liste[k].tide);

    for(i=0;i<liste[k].nbpoints;i++) {//boucle sur le nombre de points
      free(liste[k].mog2d[i]);
    }
    free(liste[k].mog2d);

    for(i=0;i<nb_ondes;i++) {//boucle sur le nombre d'ondes analys�es
      free(liste[k].amplitude[i]);
    }
    free(liste[k].amplitude);

    for(i=0;i<nb_ondes;i++) { //boucle sur le nombre d'ondes analys�es
      free(liste[k].phase[i]);
    }
    free(liste[k].phase);
  }

}



/*Routine pour lancer l'analyse harmonique*/
/* int analyse_harmonique (track_data *data, int *nb_ondes, int nb_traces, int choix_zone, int choix_type, char *dir_tide,char *dir_ah, char *dir_ch) */

char  **analyse_harmonique (track_data *data,int nb_traces, int choix_zone, int choix_type, char *dir_tide,char *dir_ah, char *dir_ch)


{
  int cpt=0,k=0;
  int a,i,j;
  double b;
  char c[50],f[50];
  char line[256];
  char e[2];
  float amp,phas;
  double lat, lon;
  char c_lat[1], c_lon[1];
  int nb_ondes,status;

  char **wave_name=NULL; //tableau contenant les noms des ondes � �tudier
  char mesg[70]= "Problème d'allocation mémoire - Programme Analyse harmonique";

  char *commande=NULL, *filename=NULL, *file_ondes=NULL;
  char *file_ini=NULL, *file_input=NULL, *file_output=NULL;
  FILE *f_data=NULL,*f_dir=NULL, *f_harmo=NULL;
  char *date_debut=NULL, *date_fin=NULL, *dir_name=NULL;


  if ((commande=(char *) malloc(250*sizeof(char))) == NULL) perror(mesg);
  if ((filename=(char *) malloc(50*sizeof(char))) == NULL) perror(mesg);
  if ((file_ondes=(char *) malloc(50*sizeof(char))) == NULL) perror(mesg);
  if ((file_ini=(char *) malloc(50*sizeof(char))) == NULL) perror(mesg);
  if ((file_input=(char *) malloc(50*sizeof(char))) == NULL) perror(mesg);
  if ((file_output=(char *) malloc(50*sizeof(char))) == NULL) perror(mesg);
  if ((date_debut=(char *) malloc(20*sizeof(char))) == NULL) perror(mesg);
  if ((date_fin=(char *) malloc(20*sizeof(char))) == NULL) perror(mesg);
  if ((dir_name=(char *) malloc(100*sizeof(char))) == NULL) perror(mesg);

  /*Bornes temporelles de l'analyse harmonique*/
  sprintf(date_debut,"-start 01/1992");
  sprintf(date_fin,"-end 12/2008");


  /*Recherche du fichier contenant les ondes � analyser : ondes.XX.dut*/
  sprintf(filename,"nom_fichier.txt");
  sprintf(commande, "ls *.dut > %s",filename);

  printf("Recherche du fichier contenant le nom des ondes : %s \n", commande);
  system(commande);
  printf("Nom du fichier : %s \n", filename);

  if((f_data=fopen(filename,"r")) == NULL) {
    fprintf(stderr,"Probleme a l ouverture du fichier des traces %s\n", filename);
  }

  fscanf(f_data,"%s",file_ondes);

  fclose(f_data);


  /*Recherche du fichier de configuration ini.hut.V...*/
  sprintf(commande, "ls *ini* > %s",filename);

  printf("Recherche du fichier contenant le nom des ondes : %s \n", commande);
  system(commande);
  printf("Nom du fichier : %s \n", filename);

  if((f_data=fopen(filename,"r")) == NULL) {
    fprintf(stderr,"Probleme a l ouverture du fichier des traces %s\n", filename);
  }

  fscanf(f_data,"%s",file_ini);

  fclose(f_data);

  /*Ouverture du fichier contenant les ondes et lecture du nombre d'ondes � analyser et de leur nom*/
  if((f_data=fopen(file_ondes,"r")) == NULL) {
    fprintf(stderr,"Probleme a l ouverture du fichier des traces %s\n", file_ondes);
  }

  fscanf(f_data,"%d",&nb_ondes);
  printf("Nombre d'ondes a analyser : %d\n", nb_ondes);

  exitIfNull(
    wave_name=(char **) malloc((nb_ondes)*sizeof(char *))
    );

  for(i=0;i<(nb_ondes);i++) {
    if ((wave_name[i]=(char *) malloc(6*sizeof(char))) == NULL) perror(mesg);
  }

  fgets(line,sizeof(line),f_data);

  for (i=0;i<(nb_ondes);i++) {
    fscanf(f_data,"%d %3s %c %lf %s",&a,wave_name[i],e,&b,c);
    printf("Nom de l'onde Num %d : %s \n\n", (i+1),wave_name[i]);
  }

  fclose(f_data);

  /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  /*Execution de l'analyse harmonique*/
  /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  /*Boucle sur le nombre de traces*/
  for (i=0;i<nb_traces;i++) {
    /*Allocation des tableaux de la structure pour amplitude, phase et erreur*/
    printf("Nombre de points de la trace %d : %d \n",data[i].num_trace,data[i].nbpoints);

    if((data[i].lat_ah = (float*)malloc(sizeof(float))) == NULL) perror(mesg);
    if((data[i].lon_ah = (float*)malloc(sizeof(float))) == NULL) perror(mesg);

    /*la 2eme dimension est allou�e plus loin au fur et a mesure de la lecture du fichier de sortie de l'ah*/
    if ((data[i].amplitude=(double **) malloc((nb_ondes)*sizeof(double))) == NULL) perror(mesg);
    if ((data[i].phase=(double **) malloc((nb_ondes)*sizeof(double))) == NULL) perror(mesg);

    for(j=0;j<(nb_ondes);j++) {
      if ((data[i].amplitude[j]=(double*)malloc(sizeof(double))) == NULL) perror(mesg);
    }

    for(j=0;j<(nb_ondes);j++) {
      if ((data[i].phase[j]=(double*)malloc(sizeof(double))) == NULL) perror(mesg);
    }

    /*Commande permettant de lancer l'outil d'analyse harmonique*/

    /*Nom du fichier d'entr�e pour la trace*/
    if (data[i].num_trace<100 && data[i].num_trace>=10) {
      sprintf(file_input,"%s/track-ref.TP+J1.0%d.dat",dir_ah,data[i].num_trace);
    }
    else if (data[i].num_trace <10) {
      sprintf(file_input,"%s/track-ref.TP+J1.00%d.dat",dir_ah,data[i].num_trace);
    }
    else {
      sprintf(file_input,"%s/track-ref.TP+J1.%d.dat",dir_ah,data[i].num_trace);
    }

    printf("Fichier ini : %s \n", file_ini);
    printf("Fichier input : %s \n", file_input);
    printf("repertoire tide : %s \n", dir_tide);

    /*Execution de l analyse harmonique*/
    sprintf(commande,"./analyse.1.2/bin/analyse.1.2 -i %s -f %s -d %s %s %s -alti",file_ini,file_input,dir_tide,date_debut,date_fin);
    printf("Execution de la commande : %s \n\n", commande);
    system(commande);

    /*Le fichier de sortie *.mgr est deplace vers repertoire commun de stockage pour toutes les traces dir_ch*/
    if((f_dir=fopen("directory.txt","r")) == NULL) {
      fprintf(stderr,"Probleme a l ouverture du fichier de repertoire %d \n",data[i].num_trace);
      return(NULL);
    }

    fscanf(f_dir,"%s", dir_name);
    printf("Repertoire de stockage des données AH de la trace %d  :  %s \n",data[i].num_trace, dir_name);
    fclose(f_dir);


    if (data[i].num_trace<100 && data[i].num_trace>=10) {
      sprintf(commande, "mv ./%s/track-ref.TP+J1.0%d.mgr %s/track-ref.TP+J1.0%d.mgr", dir_name,data[i].num_trace,dir_ch,data[i].num_trace);
      sprintf(file_output,"%s/track-ref.TP+J1.0%d.mgr",dir_ch,data[i].num_trace);
    }
    else if (data[i].num_trace <10) {
      sprintf(commande, "mv ./%s/track-ref.TP+J1.00%d.mgr %s/track-ref.TP+J1.00%d.mgr", dir_name,data[i].num_trace,dir_ch,data[i].num_trace);
      sprintf(file_output,"%s/track-ref.TP+J1.00%d.mgr",dir_ch,data[i].num_trace);
    }
    else {
      sprintf(commande, "mv ./%s/track-ref.TP+J1.%d.mgr %s/track-ref.TP+J1.%d.mgr", dir_name,data[i].num_trace,dir_ch,data[i].num_trace);
      sprintf(file_output,"%s/track-ref.TP+J1.%d.mgr",dir_ch,data[i].num_trace);
    }

    printf("Execution de la commande : %s \n\n", commande);
    status=system(commande);

    if (status != 0) { // pour certaines traces, le fichier mgr genere est directement dans ./src/ : traces 134 141 et 236 (???)
      if (data[i].num_trace<100 && data[i].num_trace>=10) {
        sprintf(commande, "mv track-ref.TP+J1.0%d.mgr %s/track-ref.TP+J1.0%d.mgr",data[i].num_trace,dir_ch,data[i].num_trace);
      }
      else if (data[i].num_trace <10) {
        sprintf(commande, "mv track-ref.TP+J1.00%d.mgr %s/track-ref.TP+J1.00%d.mgr",data[i].num_trace,dir_ch,data[i].num_trace);
      }
      else {
        sprintf(commande, "mv track-ref.TP+J1.%d.mgr %s/track-ref.TP+J1.%d.mgr",data[i].num_trace,dir_ch,data[i].num_trace);
      }
      printf("Execution de la commande 2 : %s\n",commande);
      system(commande);
    }
    /*Ouverture du fichier de sortie de l'AH*/
    if((f_harmo=fopen(file_output,"r")) == NULL) {
      fprintf(stderr,"Probleme a l ouverture du fichier de sortie harmo trace %d \n", data[i].num_trace);
      return(NULL);
    }

    data[i].nbpoints_ah=0;

    /*Lecture du fichier de sortie pour mettre les valeurs d'amplitude/phase dans la structure*/
    do {
      /*Entete*/
      fgets(line,sizeof(line),f_harmo);//_______
      fgets(line,sizeof(line),f_harmo); //espace
      fgets(line,sizeof(line),f_harmo); //station
      fgets(line,sizeof(line),f_harmo); //origine

      fgets(line,sizeof(line),f_harmo); //localisation
      sscanf(line,"%12s : %lf%1s  %lf%1s %20s",f,&lat,c_lat,&lon,c_lon,c);
      // fscanf(f_harmo, "%12s : %f%1s  %f%1s %20s",f,&lat,c_lat,&lon,c_lon,c);
      // printf("%d  %f %s   %f %s   \n",k,lat,c_lat,lon,c_lon);

      /*Allocation des vecteurs de la structure contenant les lon,lat des points ou l'ah est appliqu�e*/
      if((data[i].lat_ah = (float*)realloc(data[i].lat_ah,sizeof(float)*((data[i].nbpoints_ah)+1))) == NULL) perror(mesg);
      if((data[i].lon_ah = (float*)realloc(data[i].lon_ah,sizeof(float)*((data[i].nbpoints_ah)+1))) == NULL) perror(mesg);

      /*Allocation des vecteurs de la structure contenant les amplitudes et phases*/
      for(j=0;j<(nb_ondes);j++) {
          if ((data[i].amplitude[j]=(double*)realloc(data[i].amplitude[j],((data[i].nbpoints_ah)+1)*sizeof(double))) == NULL) perror(mesg);
        }

      for(j=0;j<(nb_ondes);j++) {
          if ((data[i].phase[j]=(double*)realloc(data[i].phase[j],((data[i].nbpoints_ah)+1)*sizeof(double))) == NULL) perror(mesg);
        }

      data[i].lat_ah[data[i].nbpoints_ah]=(float)lat;
      data[i].lon_ah[data[i].nbpoints_ah]=(float)lon;
                 

      fgets(line,sizeof(line),f_harmo); //enregistrement
      fgets(line,sizeof(line),f_harmo); //validation
                 
      /*Constantes harmoniques des ondes etudiees*/
      for(j=0;j<(nb_ondes);j++) {
          fgets(line,sizeof(line),f_harmo);
          sscanf(line, " %d  %s %f %f",&a,c,&amp,&phas);
          data[i].amplitude[j][data[i].nbpoints_ah]=((double)amp);
          data[i].phase[j][data[i].nbpoints_ah]=((double)phas);
          //printf("%d %d  %s  %lf %lf\n",j,k,c,data[i].amplitude[j][data[i].nbpoints_ah],data[i].phase[j][data[i].nbpoints_ah]);
        }

      if (k=!(feof(f_harmo))) //test si on est � la fin du fichier
        {
          data[i].nbpoints_ah++;
        }
        
    } while(k);

    printf("Nombre de points de l'AH pour la trace %d : %d\n", data[i].num_trace, data[i].nbpoints_ah);

    fclose(f_harmo);
  }

  free(commande);
  free(filename);
  free (file_ondes);
  free(file_ini);
  free(file_input);
  free(file_output);
  free(date_debut);
  free(date_fin);
  free(dir_name);

  return wave_name;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void XTRACK_WriteAsciiRecords_obsolete(FILE *stream, serie_t *record, int *index, size_t size)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
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
    }
  fflush(stream);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void XTRACK_WriteAsciiRecords_obsolete(FILE *stream, serie_t *record, int *index, size_t size, const char *rootname)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  FILE *gnu;
  char gnufile[256];
  
  for(size_t count = 0; count < size ; count++){
    sprintf(gnufile,"%s.%3.3d.gnu",rootname,index[count]);
    gnu=fopen(gnufile,"w");
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
    
    for(int i = 0; i < record[count].count; i++){
      if(is_equal(record[count].data[i].values[X_TRACK_SERIE_T_RESIDUAL],record[count].mask,float(1e-3)) ) continue;
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

vector<serie_t> ctoh_load_metadata_ref_obsolete(string filename)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status = NC_NOERR + 1;

  int nc_id;
  if( (status = nc_open((char*) filename.c_str(), NC_NOWRITE, &nc_id)) != NC_NOERR){
    check_error(status,"nc_open() failed", __LINE__,__FILE__, 1);
    }


/* *----------------------------------------------------------------------------
  get dimensions */
  int dim_cycle_id;
  if( (status = nc_inq_dimid(nc_id, "nbcycles", &dim_cycle_id)) != NC_NOERR){
    check_error(status, "nc_inq_dimid() failed", __LINE__,__FILE__, 1);
    }
  size_t i_CYCLES_dim;
  if( (status = nc_inq_dimlen(nc_id, dim_cycle_id, &i_CYCLES_dim)) != NC_NOERR){
    check_error(status, "nc_inq_dimlen() failed", __LINE__,__FILE__, 1);
    }

  int dim_points_id;
  if( (status = nc_inq_dimid(nc_id, "nbpoints", &dim_points_id)) != NC_NOERR){
    check_error(status, "nc_inq_dimid() failed", __LINE__,__FILE__, 1);
    }
  size_t i_POINTS_dim;
  if( (status = nc_inq_dimlen(nc_id, dim_points_id, &i_POINTS_dim)) != NC_NOERR){
    check_error(status, "nc_inq_dimlen() failed", __LINE__,__FILE__, 1);
    }

/* *----------------------------------------------------------------------------
  Laurent, il faudra un jour que tu m'explique pourquoi tu utilises des types qui
  te compliquent la mise au point !!!
  Tes codes sont deja assez confus comme ca... */
  vector<serie_t> track(i_POINTS_dim);
  for(size_t i = 0; i < i_POINTS_dim; ++i){
    track[i].data = (data_t *) new data_t[i_CYCLES_dim];
    track[i].count = i_CYCLES_dim;
    track[i].mask = static_cast<float>(9999);
    }
    
  float epsilon = static_cast<float>(1.e-3);

/* *----------------------------------------------------------------------------
  get longitudes */
  int lon_id;
  if( (status = nc_inq_varid(nc_id, "lon", &lon_id)) != NC_NOERR){
    check_error(status, "nc_inq_varid() failed", __LINE__,__FILE__, 1);
    }
  float scale;
  status = NC_NOERR + 13;
  if( (status = nc_get_att_float(nc_id, lon_id, "scale_factor", &scale)) != NC_NOERR){
    scale = 1.;
    }
  float offset;
  status = NC_NOERR + 13;
  if( (status = nc_get_att_float(nc_id, lon_id, "offset", &offset)) != NC_NOERR){
    offset = 0.;
    }
  float fillValue;
  status = NC_NOERR + 13;
  if( (status = nc_get_att_float(nc_id, lon_id, "_FillValue", &fillValue)) != NC_NOERR){
    check_error(status, "nc_get_att_float() failed", __LINE__,__FILE__, 1);
    }
  float *lon = new float[i_POINTS_dim];
  if( (status = nc_get_var_float(nc_id,lon_id,lon)) != NC_NOERR){
    check_error(status, "nc_get_var_float() failed", __LINE__,__FILE__, 1);
    }

  for(size_t i = 0; i < track.size(); ++i){
    for(size_t j = 0; j < track[i].count; ++j){
      if(!is_equal(lon[i], fillValue, epsilon) ){
        track[i].data[j].lon = lon[i] * scale + offset;
        }
      else{
        track[i].data[j].lon = track[i].mask;
        }
      }
    }
  delete [] lon;


/* *----------------------------------------------------------------------------
  get latitudes */
  int lat_id;
  if( (status = nc_inq_varid(nc_id,"lat",&lat_id)) != NC_NOERR){
    check_error(status, "nc_inq_varid() failed", __LINE__,__FILE__, 1);
    }
  status = NC_NOERR + 13;
  if( (status = nc_get_att_float(nc_id, lat_id, "scale_factor", &scale)) != NC_NOERR){
    scale = 1;
    }
  status = NC_NOERR + 13;
  if( (status = nc_get_att_float(nc_id, lat_id, "offset", &offset)) != NC_NOERR){
    offset = 0;
    }
  status = NC_NOERR + 13;
  if( (status = nc_get_att_float(nc_id, lat_id, "_FillValue", &fillValue)) != NC_NOERR){
    check_error(status, "nc_get_att_float() failed", __LINE__,__FILE__, 1);
    }
  float *lat = new float[i_POINTS_dim];
  if( (status = nc_get_var_float(nc_id,lat_id,lat)) != NC_NOERR){
    check_error(status, "nc_get_var_float() failed", __LINE__,__FILE__, 1);
    }
  for(size_t i = 0; i < track.size(); ++i){
    for(size_t j = 0; j < track[i].count; ++j){
      if(!is_equal(lat[i], fillValue, epsilon) ){
        track[i].data[j].lat = lat[i] * scale + offset;
        }
      else{
        track[i].data[j].lat = track[i].mask;
        }
      }
    }
  delete [] lat;

 
/* *----------------------------------------------------------------------------
  get times */
  int time_id;
  if( (status = nc_inq_varid(nc_id,"time",&time_id)) != NC_NOERR){
    check_error(status, "nc_inq_varid() failed", __LINE__,__FILE__, 1);
    }
    
  double dscale = 1;
  status = NC_NOERR + 13;
  if( (status = nc_get_att_double(nc_id, time_id, "scale_factor", &dscale)) != NC_NOERR){
    dscale = 1;
    }
  double doffset = 0;
  status = NC_NOERR + 13;
  if( (status = nc_get_att_double(nc_id, time_id, "offset", &doffset)) != NC_NOERR){
    doffset = 0;
    }
  double dfillValue = 0;
  status = NC_NOERR + 13;
  if( (status = nc_get_att_double(nc_id, time_id, "_FillValue", &dfillValue)) != NC_NOERR){
    check_error(status, "nc_get_att_double() failed", __LINE__,__FILE__, 1);
    }
  double *time = new double[i_POINTS_dim * i_CYCLES_dim];
  if( (status = nc_get_var_double(nc_id,time_id,time)) != NC_NOERR){
    check_error(status, "nc_get_var_double() failed", __LINE__,__FILE__, 1);
    }

  for(size_t i = 0; i < track.size(); ++i){
    for(size_t j = 0; j < track[i].count; ++j){
      if(!is_equal(time[i * i_CYCLES_dim + j], dfillValue, (double) epsilon) ){
        track[i].data[j].time = time[i * i_CYCLES_dim + j] * dscale + doffset;
        }
      else{
        track[i].data[j].time = track[i].mask;
        }
      }
    }
  delete [] time;

/* *----------------------------------------------------------------------------
  get tides */
  int tide_id;
  if( (status = nc_inq_varid(nc_id,"tide",&tide_id)) != NC_NOERR){
    check_error(status, "nc_inq_varid() failed", __LINE__,__FILE__, 1);
    }
  status = NC_NOERR + 13;
  if( (status = nc_get_att_float(nc_id, tide_id, "scale_factor", &scale)) != NC_NOERR){
    scale = 1;
    }
  status = NC_NOERR + 13;
  if( (status = nc_get_att_float(nc_id, tide_id, "offset", &offset)) != NC_NOERR){
    offset = 0;
    }
  status = NC_NOERR + 13;
  if( (status = nc_get_att_float(nc_id, tide_id, "_FillValue", &fillValue)) != NC_NOERR){
    check_error(status, "nc_get_att_float() failed", __LINE__,__FILE__, 1);
    }
  float *tide = new float[i_POINTS_dim * i_CYCLES_dim];
  if( (status = nc_get_var_float(nc_id, tide_id, tide)) != NC_NOERR){
    check_error(status, "nc_get_var_float() failed",__LINE__,__FILE__, 1);
    }

  for(size_t i = 0; i < track.size(); ++i){
    for(size_t j = 0; j < track[i].count; ++j){
      if(!is_equal(tide[i * i_CYCLES_dim + j], fillValue, epsilon) ){
        tide[i * i_CYCLES_dim + j] = tide[i * i_CYCLES_dim + j] * scale + offset;
        track[i].data[j].values[X_TRACK_SERIE_T_TIDE_REGIONAL] = floor(tide[i * i_CYCLES_dim + j] * 10000. +0.5) / 10000.;
        track[i].data[j].values[X_TRACK_SERIE_T_TIDE_GOT4_7]   = floor(tide[i * i_CYCLES_dim + j] * 10000. +0.5) / 10000.;
        track[i].data[j].values[X_TRACK_SERIE_T_TIDE_FES2004]  = floor(tide[i * i_CYCLES_dim + j] * 10000. +0.5) / 10000.;
        }
      else{
        track[i].data[j].values[X_TRACK_SERIE_T_TIDE_REGIONAL] = track[i].mask;
        track[i].data[j].values[X_TRACK_SERIE_T_TIDE_GOT4_7]   = track[i].mask;
        track[i].data[j].values[X_TRACK_SERIE_T_TIDE_FES2004]  = track[i].mask;
        }
      }
    }
  delete [] tide;

/* *----------------------------------------------------------------------------
  get DAC */
  int dac_id;
  
  status = nc_inq_varid(nc_id,"mog2d",&dac_id);
  if( status != NC_NOERR){
//    check_error(status, "nc_inq_varid() failed",__LINE__,__FILE__, 1);
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
  status = NC_NOERR + 13;
  if( (status = nc_get_att_float(nc_id, dac_id, "scale_factor", &scale)) != NC_NOERR){
    scale = 1;
    }
  status = NC_NOERR + 13;
  if( (status = nc_get_att_float(nc_id, dac_id, "offset", &offset)) != NC_NOERR){
    offset = 0;
    }
  status = NC_NOERR + 13;
  if( (status = nc_get_att_float(nc_id, dac_id, "_FillValue", &fillValue)) != NC_NOERR){
    check_error(status, "nc_get_att_float() failed", __LINE__,__FILE__, 1);
    }
  float *dac = new float[i_POINTS_dim * i_CYCLES_dim];
  status = nc_get_var_float(nc_id,dac_id,dac);
  if( status != NC_NOERR){
    check_error(status, "nc_get_var_float() failed",__LINE__,__FILE__, 1);
    }
  for(size_t i = 0; i < track.size(); ++i){
    for(size_t j = 0; j < track[i].count; ++j){
      if(!is_equal(dac[i * i_CYCLES_dim + j], fillValue, epsilon) ){
        dac[i * i_CYCLES_dim + j] = dac[i * i_CYCLES_dim + j] * scale + offset;
        track[i].data[j].values[X_TRACK_SERIE_T_DAC_GLOBAL]   = floor(dac[i * i_CYCLES_dim + j] * 10000. +0.5) / 10000.;
        track[i].data[j].values[X_TRACK_SERIE_T_DAC_REGIONAL] = floor(dac[i * i_CYCLES_dim + j] * 10000. +0.5) / 10000.;
        }
      else{
        track[i].data[j].values[X_TRACK_SERIE_T_DAC_GLOBAL]   = track[i].mask;
        track[i].data[j].values[X_TRACK_SERIE_T_DAC_REGIONAL] = track[i].mask;
        }
      }
    }
  delete [] dac;
  }

/* *----------------------------------------------------------------------------
  get SLA */
  int sla_id;
  if( (status = nc_inq_varid(nc_id,"sla",&sla_id)) != NC_NOERR){
    check_error(status, "nc_inq_varid() failed", __LINE__,__FILE__, 1);
    }
  status = NC_NOERR + 13;
  if( (status = nc_get_att_float(nc_id, sla_id, "scale_factor", &scale)) != NC_NOERR){
    scale = 1;
    }
  status = NC_NOERR + 13;
  if( (status = nc_get_att_float(nc_id, sla_id, "offset", &offset)) != NC_NOERR){
    offset = 0;
    }
  status = NC_NOERR + 13;
  if( (status = nc_get_att_float(nc_id, sla_id, "_FillValue", &fillValue)) != NC_NOERR){
    check_error(status, "nc_get_att_float() failed", __LINE__,__FILE__, 1);
    }
  float *sla = new float[i_POINTS_dim * i_CYCLES_dim];
  if( (status = nc_get_var_float(nc_id,sla_id,sla)) != NC_NOERR){
    check_error(status,"nc_get_var_float() failed", __LINE__,__FILE__, 1);
    }
  for(size_t i = 0; i < track.size(); ++i){
    for(size_t j = 0; j < track[i].count; ++j){
      if(!is_equal(sla[i * i_CYCLES_dim + j], fillValue, epsilon) ){
        sla[i * i_CYCLES_dim + j] = sla[i * i_CYCLES_dim + j] * scale + offset;
        track[i].data[j].values[X_TRACK_SERIE_T_SSHA] = floor(sla[i * i_CYCLES_dim + j] * 10000. +0.5) / 10000.; // caution tides and DAC already applied
        }
      else{
        track[i].data[j].values[X_TRACK_SERIE_T_SSHA] = track[i].mask;
        }
      }
    }
  delete [] sla;

  if( (status = nc_close(nc_id)) != NC_NOERR){
    check_error(status, "nc_close() failed", __LINE__,__FILE__, 1);
    }

/* *----------------------------------------------------------------------------
  remove data with mask values */
  for(size_t r = 0; r < track.size(); ++r) {
    size_t n = track[r].count;
    data_t* tmp = new data_t[n]; //over alloc
    size_t count = 0;
    for(size_t i = 0; i < n; i++){
      float mask  = track[r].mask;
      double time = track[r].data[i].time;
      float ssha  = track[r].data[i].values[X_TRACK_SERIE_T_SSHA];
      float tide  = track[r].data[i].values[X_TRACK_SERIE_T_TIDE_GOT4_7];
      bool time_is_valid = (!is_equal(time, static_cast<double>(mask), double(1e-3)));
      bool ssha_is_valid = (!is_equal(ssha, mask, float(1e-3)));
      bool tide_is_valid = (!is_equal(tide, mask, float(1e-3)));
      if(ssha_is_valid != tide_is_valid) {
        printf("time serie anomaly record %, point %d: ssh flag=%d tide flag=%d\n", r, i, ssha_is_valid, tide_is_valid);
        }
      if(time_is_valid && ssha_is_valid && tide_is_valid) {
        tmp[count] = track[r].data[i];
        count++;
        }
      }
    track[r].count = count;
    if(count == 0){
      delete[] tmp;
      continue;
      }
    delete [] track[r].data;
    track[r].data = new data_t[count];
    for(size_t i = 0; i < count; ++i){
      track[r].data[i] = tmp[i];
      }
    delete[] tmp;
    }
 

/* *----------------------------------------------------------------------------
  set back tides and dac, already applied in SLA */
  for(size_t r = 0; r < track.size(); ++r) {
    if(track[r].count==0) continue;
//    printf("track %d: treating %d valid data\n",r,track[r].count);
    for(size_t i = 0; i < track[r].count; ++i){
      float ssha =track[r].data[i].values[X_TRACK_SERIE_T_SSHA];
      float dac  =track[r].data[i].values[X_TRACK_SERIE_T_DAC_GLOBAL];
      float tide =track[r].data[i].values[X_TRACK_SERIE_T_TIDE_GOT4_7];
      ssha+=dac+tide;
      track[r].data[i].values[X_TRACK_SERIE_T_SSHA] =ssha;
      }
   }
 
/* *----------------------------------------------------------------------------
  remove empty time series */
/**----------------------------------------------------------------------------
  erase totalement inefficace, il faut trouver une autre soluti !!! */
  for(size_t r = 0; r < track.size(); ++r) {
    if(track[r].count == 0){
//      track.erase(track.begin() + r);
//      delete[] track[r].data;
      }
    }
 
  return(track);


}

