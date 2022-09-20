

/*******************************************************************************

  T-UGO tools, 2006-2009

  Unstructured Ocean Grid initiative

Contributors:

  Florent Lyard      LEGOS/CNRS, Toulouse, France
  Yoann Le Bars      LEGOS, Toulouse, France (PhD)
  Laurent Roblou     LEGOS, Toulouse, France
  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
  Frédéric Dupont    Université de Laval à Québec, Canada

E-mail: florent.lyard@legos.obs-mip.fr

*******************************************************************************/

#define MAIN_SOURCE

#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "tools-define.h"
#include "tools-structures.h"

#include "functions.h"
#include "fe.h"
#include "bmg.h"
#include "map.h"
#include "ascii.h"
#include "archive.h"
#include "geo.h"
#include "sym-io.h"
#include "xtrack-io.h"


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int decode_trackfile_name(char *path, const char *name_template, int track, char **filename)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int  status=0;
  char dummy[256], *pointer,*prefix,*suffix;
  FILE *in;
  string s;
  size_t pos;

  /* *----------------------------------------------------------------------------
     Reconstruct the input file name*/
  (*filename) = new char[strlen(path) + strlen(name_template) + 256];

  /* *----------------------------------------------------------------------------
     former filename pattern*/
  sprintf(*filename,"%s/%s.%3.3d.ref.sla.nc",path,name_template,track);

  in=fopen(*filename,"r");

  if(in!=0) {
    fclose(in);
    status=0;
    return (status);
  }

  sprintf((*filename), "%s/%s", path, name_template);

  pointer = strstr((char *) name_template, "NUMTRACK");
  if(pointer == NULL) {
    status=-1;
    return (status);
  }

  //  s=string(filename);

  //  pos=strchr(name_template, "NUMTRACK");
  pos=(size_t) (pointer-name_template );

  prefix=new char[strlen(name_template)];
  suffix=new char[strlen(name_template)];

  strncpy(prefix,name_template,pos);
  prefix[pos]=0;
  strncpy(suffix,&(name_template[pos+8]),strlen(name_template)-8-pos);
  suffix[strlen(name_template)-8-pos]=0;

  sprintf(dummy, "%3.3d", track);
  sprintf((*filename), "%s/%s%s%s", path, prefix,dummy,suffix);

  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  track_data *read_data_CTOH(char *dir_input_alti, char *convention, int *nb_traces_ctoh, int choix_zone)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{

  int i,j,k,m,status;
  int nc_id,dim_cycle_id,dim_points_id;
  int lon_id,lat_id,mssh_id,cycle_id,point_id,time_id,sla_id,tide_id,mog2d_id;
  //  char filename[100];
  char *filename=0;
//   char filename2[100];
  int nb_data;
  int traces[80];

  int cpt=0;
  track_data *liste=NULL;

  double *temps=NULL;
  float *sla=NULL,*tide=NULL,*mog2d=NULL;
  FILE *f_data=NULL/*, *out=NULL*/;

  size_t i_CYCLES_dim,i_POINTS_dim;

  char mesg[70]= "Problème d'allocation mémoire - Programme READ TPJ";

  /*Ouverture du fichier contenant tous les numéros de trace de la zone*/

  switch(choix_zone) {
  case -1:
    f_data=fopen("traces.txt","r");
    if(f_data==0) {
      fprintf(stderr,"Probleme a l ouverture du fichier des traces\n");
      }
    while (feof(f_data) == 0) {
      cpt++;
      fscanf(f_data,"%d ",&traces[cpt-1]);
      }
    fclose(f_data);
    break;

  case 0:
    f_data=fopen("traces_ctoh_NEA.txt","r");
    if(f_data==0) {
      fprintf(stderr,"Probleme a l ouverture du fichier des traces NEA\n");
      }
    while (feof(f_data) == 0) {
      cpt++;
      fscanf(f_data,"%d ",&traces[cpt-1]);
      }
    fclose(f_data);
    break;

  case 1:
    f_data=fopen("traces_ctoh_MED.txt","r");
    if(f_data==0) {
      fprintf(stderr,"Probleme a l ouverture du fichier des traces MED\n");
      }
    while (feof(f_data) == 0) {
      cpt++;
      fscanf(f_data,"%d ",&traces[cpt-1]);
      }
    fclose(f_data);
    break;

  case 2:
    f_data=fopen("traces_ctoh_ATLSUD.txt","r");
    if(f_data==0) {
      fprintf(stderr,"Probleme a l ouverture du fichier des traces ATLSUD\n");
      }
    while (feof(f_data) == 0) {
      cpt++;
      fscanf(f_data,"%d ",&traces[cpt-1]);
      }
    fclose(f_data);
    break;
    }

  *nb_traces_ctoh = cpt;

  /*Allocation de mémoire pour le tableau de structure*/
  if ((liste=(track_data *) malloc(cpt*sizeof(track_data))) == NULL) perror(mesg);

  /*Boucle sur toutes les traces*/
  for (k=0;k<cpt;k++) {
    liste[k].num_trace=traces[k];
    printf("Numero de la trace : %d \n",liste[k].num_trace);
    /* Ouverture du fichier NetCDF de la trace concernée*/

    //     if (liste[k].num_trace<100 && liste[k].num_trace>=10) {
    //       sprintf(filename,"%s/%s.0%d.ref.sla.nc",dir_input_alti,convention,liste[k].num_trace);
    //       }
    //     else if (liste[k].num_trace<10) {
    //       sprintf(filename,"%s/%s.00%d.ref.sla.nc",dir_input_alti,convention,liste[k].num_trace);
    //       }
    //     else {
    //       sprintf(filename,"%s/%s.%d.ref.sla.nc",dir_input_alti,convention,liste[k].num_trace);
    //       }
    if(filename!=0) delete[] filename;
    status= decode_trackfile_name(dir_input_alti,convention, liste[k].num_trace, &filename);

    printf("Filename : %s\n", filename);
    status=nc_open(filename,NC_NOWRITE,&nc_id);
    NC_TRAP_ERROR(wexit,status,1,"error");

    /* Lecture de la dimension de chaque variable du fichier de données*/

    /*Nombre de cycles*/
    status=nc_inq_dimid(nc_id,"nbcycles",&dim_cycle_id);
    NC_TRAP_ERROR(wexit,status,1,"error");
    status=nc_inq_dimlen(nc_id,dim_cycle_id,&i_CYCLES_dim);
    NC_TRAP_ERROR(wexit,status,1,"error");

    /*Nombre de points*/
    status=nc_inq_dimid(nc_id,"nbpoints",&dim_points_id);
    NC_TRAP_ERROR(wexit,status,1,"error");
    status=nc_inq_dimlen(nc_id,dim_points_id,&i_POINTS_dim);
    NC_TRAP_ERROR(wexit,status,1,"error");

    /*Allocation de l'espace necessaire pour le stockage des données*/

    liste[k].nbpoints=i_POINTS_dim;
    printf("Nombre de points : %d \n",liste[k].nbpoints);
    liste[k].nbcycles=i_CYCLES_dim;
    printf("Nombre de cycles : %d\n",liste[k].nbcycles);
    printf("\n");

    size_t start[2]={0,0};
    size_t count[2]={liste[k].nbpoints, liste[k].nbcycles};
    nb_data=liste[k].nbcycles*liste[k].nbpoints;

    /*penser à mettre des warning au cas ou pbm de place*/
    if ((liste[k].lon=(float *) malloc(liste[k].nbpoints*sizeof(float))) == NULL) perror(mesg);
    if ((liste[k].lat=(float *) malloc(liste[k].nbpoints*sizeof(float))) == NULL) perror(mesg);
    if ((liste[k].mssh=(float *) malloc(liste[k].nbpoints*sizeof(float))) == NULL) perror(mesg);

    if ((liste[k].cycle=(int *) malloc(liste[k].nbcycles*sizeof(int))) == NULL) perror(mesg);
    if ((liste[k].point=(int *) malloc(liste[k].nbpoints*sizeof(int))) == NULL) perror(mesg);

    if ((liste[k].time=(double **) malloc(liste[k].nbpoints*sizeof(double))) == NULL) perror(mesg);
    if ((temps=(double*)malloc(liste[k].nbpoints*liste[k].nbcycles*sizeof(double))) == NULL) perror(mesg);

    for(i=0;i<liste[k].nbpoints;i++) {
      if ((liste[k].time[i]=(double*)malloc(liste[k].nbcycles*sizeof(double))) == NULL) perror(mesg);
      }

    if ((liste[k].sla=(float **) malloc(liste[k].nbpoints*sizeof(float *))) == NULL) perror(mesg);
    if ((sla=(float*) malloc(liste[k].nbpoints*liste[k].nbcycles*sizeof(float))) == NULL) perror(mesg);

    for(i=0;i<liste[k].nbpoints;i++) {
      if ((liste[k].sla[i]=(float*)malloc(liste[k].nbcycles*sizeof(float))) == NULL) perror(mesg);
      }

    if ((liste[k].tide=(float **) malloc(liste[k].nbpoints*sizeof(float *))) == NULL) perror(mesg);
    if ((tide=(float*)malloc(liste[k].nbpoints*liste[k].nbcycles*sizeof(float))) == NULL) perror(mesg);

    for(i=0;i<liste[k].nbpoints;i++) {
      if ((liste[k].tide[i]=(float*)malloc(liste[k].nbcycles*sizeof(float))) == NULL) perror(mesg);
      }


    if ((liste[k].mog2d=(float **) malloc(liste[k].nbpoints*sizeof(float *))) == NULL) perror(mesg);
    if ((mog2d=(float*)malloc(liste[k].nbpoints*liste[k].nbcycles*sizeof(float))) == NULL) perror(mesg);

    for(i=0;i<liste[k].nbpoints;i++) {
      if ((liste[k].mog2d[i]=(float*)malloc(liste[k].nbcycles*sizeof(float))) == NULL) perror(mesg);
      }

    /*Lecture du contenu des variables et envoi dans la structure de données*/

    /*  vecteur "lon" */
    status=nc_inq_varid(nc_id,"lon",&lon_id);
    NC_TRAP_ERROR(wexit,status,1,"error");

    status=nc_get_var_float(nc_id,lon_id,liste[k].lon);
    NC_TRAP_ERROR(wexit,status,1,"error");

    /*  vecteur "lat" */
    status=nc_inq_varid(nc_id,"lat",&lat_id);
    NC_TRAP_ERROR(wexit,status,1,"error");

    status=nc_get_var_float(nc_id,lat_id,liste[k].lat);
    NC_TRAP_ERROR(wexit,status,1,"error");

    /*  vecteur "mssh" */
    status=nc_inq_varid(nc_id,"mssh",&mssh_id);
    NC_TRAP_ERROR(wexit,status,1,"error");

    status=nc_get_var_float(nc_id,mssh_id,liste[k].mssh);
    NC_TRAP_ERROR(wexit,status,1,"error");

    /*  vecteur "cycles" */
    status=nc_inq_varid(nc_id,"cycle",&cycle_id);
    NC_TRAP_ERROR(wexit,status,1,"error");

    status=nc_get_var_int(nc_id,cycle_id,liste[k].cycle);
    NC_TRAP_ERROR(wexit,status,1,"error");

    /*  vecteur "points" */
    status=nc_inq_varid(nc_id,"point",&point_id);
    NC_TRAP_ERROR(wexit,status,1,"error");

    status=nc_get_var_int(nc_id,point_id,liste[k].point);
    NC_TRAP_ERROR(wexit,status,1,"error");

    /*  vecteur "time" */

    status=nc_inq_varid(nc_id,"time",&time_id);
    NC_TRAP_ERROR(wexit,status,1,"error");

    status=nc_get_vara_double(nc_id,time_id,start,count,temps);
    NC_TRAP_ERROR(wexit,status,1,"error");

    for (j=0;j<liste[k].nbcycles;j++) {
      for (i=0;i<liste[k].nbpoints;i++) {
//         m=nb_data-(i+1)*liste[k].nbcycles+j;
        m=i*liste[k].nbcycles+j;
        liste[k].time[i][j]=temps[m];
        }
      }

    /*  vecteur "SLA" */

    status=nc_inq_varid(nc_id,"slafilt",&sla_id);
    //NC_CHKERR_BASE_LINE(status,"error");
    if(status!= NC_NOERR) status=nc_inq_varid(nc_id,"sla",&sla_id);
    NC_TRAP_ERROR(wexit,status,1,"error");

    status=nc_get_vara_float(nc_id,sla_id,start,count,sla);
    NC_TRAP_ERROR(wexit,status,1,"error");

    for (j=0;j<liste[k].nbcycles;j++) {
      for (i=0;i<liste[k].nbpoints;i++) {
//         m=nb_data-(i+1)*liste[k].nbcycles+j;
        m=i*liste[k].nbcycles+j;
        liste[k].sla[i][j]=sla[m];
        }
      }

    /*  vecteur "Tide" */

    status=nc_inq_varid(nc_id,"tidefilt",&tide_id);
    //NC_CHKERR_BASE_LINE(status,"error");
    if(status!= NC_NOERR) status=nc_inq_varid(nc_id,"tide",&tide_id);
    NC_TRAP_ERROR(wexit,status,1,"error");

    status=nc_get_vara_float(nc_id,tide_id,start,count,tide);
    NC_TRAP_ERROR(wexit,status,1,"error");


    for (j=0;j<liste[k].nbcycles;j++) {
      for (i=0;i<liste[k].nbpoints;i++) {
//         m=nb_data-(i+1)*liste[k].nbcycles+j;
        m=i*liste[k].nbcycles+j;
        liste[k].tide[i][j]=tide[m];
        }
      }

    /*  vecteur "Mog2D" */

    status=nc_inq_varid(nc_id,"mog2dfilt",&mog2d_id);
    //NC_CHKERR_BASE_LINE(status,"error");
    if(status!= NC_NOERR) status=nc_inq_varid(nc_id,"mog2d",&mog2d_id);
    NC_TRAP_ERROR(wexit,status,1,"error");

    status=nc_get_vara_float(nc_id,mog2d_id,start,count,mog2d);
    NC_TRAP_ERROR(wexit,status,1,"error");


    for (j=0;j<liste[k].nbcycles;j++) {
      for (i=0;i<liste[k].nbpoints;i++) {
        //         m=nb_data-(i+1)*liste[k].nbcycles+j;
        m=i*liste[k].nbcycles+j;
        liste[k].mog2d[i][j]=mog2d[m];
        }
      }

    /*Formation de fichiers contenant les coordonn�es des points le long de la trace => sert ensuite pour calculer les coordonn�es des points de croisement (cf. Gwen)*/

    /*    if (liste[k].num_trace<100 && liste[k].num_trace>=10) */
    /*     { */
    /*       sprintf(filename2,"%s/tracemoy_Ps0%dtp.dat",dir_input_alti,liste[k].num_trace); */
    /*     } */
    /*   else if (liste[k].num_trace <10) */
    /*     { */
    /*      sprintf(filename2,"%s/tracemoy_Ps00%dtp.dat",dir_input_alti,liste[k].num_trace); */
    /*     } */
    /*   else  */
    /*     { */
    /*       sprintf(filename2,"%s/tracemoy_Ps%dtp.dat",dir_input_alti,liste[k].num_trace); */
    /*     } */


    /*      out=fopen(filename2,"w"); */


    /*        for (i=0;i<liste[k].nbpoints;i++) */
    /* 	 { */
    /* 	   fprintf(out,"%f %f %f %d \n",liste[k].lon[i],liste[k].lat[i],liste[k].mssh[i],0);   */
    /* 	 } */

    /*        fclose(out); */

    /* Fermeture du fichier NetCDF*/

    status=nc_close(nc_id);
    NC_TRAP_ERROR(wexit,status,1,"error");
    }

  printf("fin du programme de lecture CTOH \n\n");

  free(temps);
  free(sla);
  free(tide);
  free(mog2d);

  return liste;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int n,status;
  char *path=NULL,output[1024];
  char *keyword,*sfile=NULL,*convention=NULL;
  string cmd;
  track_data *data;
  int nb_traces_ctoh;
  int choix_zone=-1;

  fct_echo( argc, argv);

  n=1;
  while (n < argc) {
    keyword=strdup(argv[n]);
    switch (keyword[0]) {
      case '-':
        switch (keyword[1]) {

        case 'c' :
          convention= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'o' :
          strcpy(output,argv[n+1]);
          n++;
          n++;
          break;

        case 'p' :
          path= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'z' :
          sscanf(argv[n+1],"%d",&choix_zone);
          n++;
          n++;
          break;

        default:
          STDOUT_BASE_LINE("unknown option %s\n",keyword);
          exit(-1);
        }
        break;

      default:
        if(sfile==NULL) {
          sfile= strdup(argv[n]);
          printf("input file=%s\n",sfile);
          n++;
          }
        else {
          STDOUT_BASE_LINE("unknown option %s\n",keyword);
          exit(-1);
          }
        break;
      }
      free(keyword);
    }

  if(convention==0) {
    convention=strdup("track-raw.TPJ1");
    printf("use %s as default rootname\n",keyword);
    }

//   int track=62;
//   char *filename;
//   status=decode_trackfile_name(path, convention, track, &filename);

  data=read_data_CTOH(path, convention, &nb_traces_ctoh, choix_zone);

  convention=strdup("track-ref.TP+J1");
  status= formatage_CTOH (nb_traces_ctoh, data, "./", convention);

  printf("end of xtrack-format ... \n");
}
