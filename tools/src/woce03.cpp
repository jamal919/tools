/* ----------------------------------------------------------------------
 * Copyright (c) 2001-2009 LEGOS
 * All rights reserved.
 *
 * Redistribution and use  in source  and binary  forms,  with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 *   1. Redistributions of  source  code must retain the  above copyright
 *      notice, this list of conditions and the following disclaimer.
 *   2. Redistributions in binary form must reproduce the above copyright
 *      notice,  this list of  conditions and the following disclaimer in
 *      the  documentation  and/or  other   materials provided  with  the
 *      distribution.
 *
 * THIS  SOFTWARE IS PROVIDED BY  THE  COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND  ANY  EXPRESS OR IMPLIED  WARRANTIES,  INCLUDING,  BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES  OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR  PURPOSE ARE DISCLAIMED. IN  NO EVENT SHALL THE COPYRIGHT
 * HOLDERS OR      CONTRIBUTORS  BE LIABLE FOR   ANY    DIRECT, INDIRECT,
 * INCIDENTAL,  SPECIAL,  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 * BUT NOT LIMITED TO, PROCUREMENT OF  SUBSTITUTE GOODS OR SERVICES; LOSS
 * OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN  CONTRACT, STRICT LIABILITY, OR
 * TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
 * USE   OF THIS SOFTWARE, EVEN   IF ADVISED OF   THE POSSIBILITY OF SUCH
 * DAMAGE.
 *
 * DATE OF CREATION: 2001
 * AUTHOR: Thierry LETELLIER, Laurent ROBLOU
 * CONTACT: thierry.letellier@free.fr,laurent.roblou@legos.obs-mip.fr
 *
 * MESSAGE FROM THE AUTHORS:
 * In consideration of this free of charge access, users shall communicate
 * the results obtained directly  from the software. They shall inform the
 * authors of such publications.
 * Thanks!
 *
 ---------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <alloca.h>

#include "analyse_alti.h"
#include "aktarus.h"


/* Functions list : */
/*  question                                                                  */
/*  map_interpolation: interpolates data from a buffer                        */
/*  foud_z: reads complex buffer                                              */
/*  bmg_readheader: reads bimg file header                                    */
/*  bmg_loadc1_analyse_alti: load complex buffer                              */
/*  loadingEC: loads either loading either tides from models                  */
/*  Wave_select: selects a wave                                               */
/*  load_serie: loads informations on sealevel time series from the data file */
/*  header: read the data file header and returns the number of Xover points  */
/*  system                                                                    */


/*----------------------------------------------------------------------------*/
void question(char inifilename[60],conf_struct *conf)
/*----------------------------------------------------------------------------*/
{
  FILE *ini;
  char a,ofilename[256],line[256],*search;
  int nsearch,i,strlong;

  ini=fopen(inifilename,"r");

  
  fgets(line,sizeof(line),ini);
  sscanf(line,"%d",&conf->nbond);
  fgets(line,sizeof(line),ini);
  sscanf(line,"%d",&conf->nbdata);
  fgets(line,sizeof(line),ini);
  sscanf(line,"%lf",&conf->Trepet);
  fgets(line,sizeof(line),ini);
  sscanf(line,"%d",&conf->filter_type);
   if (conf->filter_type==1) {
      conf->geo=calloc(4,sizeof(double));
      sscanf(line,"%d  %lf  %lf %lf %lf",&conf->filter_type,&conf->geo[0],&conf->geo[1],&conf->geo[2],&conf->geo[3]);
    }
  if (conf->filter_type==2) {
      sscanf(line,"%d  %s ",&conf->filter_type,&conf->poly_filename);
      
    }
 fgets(line,sizeof(line),ini);
  sscanf(line,"%d",&conf->inv_bar);
  fgets(line,sizeof(line),ini);
  sscanf(line,"%d",&conf->load_effect);
  fgets(line,sizeof(line),ini);
  sscanf(line,"%d",&conf->solid_effect);
  fgets(line,sizeof(line),ini);
  sscanf(line,"%d",&conf->analys_type);
  fgets(line,sizeof(line),ini);
  sscanf(line,"%d",&conf->ios_output);
  fgets(line,sizeof(line),ini);
  sscanf(line,"%d",&conf->mgr_output);
  fgets(line,sizeof(line),ini);
  sscanf(line,"%c",&conf->format);
  fgets(line,sizeof(line),ini);
  sscanf(line,"%d",&conf->time_unit);
  fgets(line,sizeof(line),ini);
  sscanf(line,"%d",&conf->time_ref);
  switch (conf->time_ref)
    {
    case(1950): {conf->CTO=0;break;}
    case(1958): {conf->CTO=1;break;}
    case(1985): {conf->CTO=2;break;}
    case(2000): {conf->CTO=3;break;}
    default: quit("error in time reference");
    }
  if(conf->format=='h')conf->CTO=0;
  fgets(line,sizeof(line),ini);
  sscanf(line,"%d",&conf->fft);
  fgets(line,sizeof(line),ini);
  sscanf(line,"%d",&conf->drift);
  fgets(line,sizeof(line),ini);
  sscanf(line,"%d",&conf->verbose);
  fgets(line,sizeof(line),ini);
  sscanf(line,"%d",&conf->CPU);
  if(conf->CPU<1)quit("error NB CPU init file");
  

 printf("\n");

 printf(">>>>>>>>>>  starting analysis >>>>>>>>>>>>\n\n");
 printf("data file: %s\n\n",conf->data_filename);
 printf("number of waves: %d,\nfiltering option: %d, \nIB correction option: %d,\n",conf->nbond,conf->filter_type,conf->inv_bar);
 printf("loading correction option: %d, \nanalyse type: %d,\n\n",conf->load_effect,conf->analys_type);
 if (conf->fft==1) printf("FFT computation,\n\n");
 if (conf->drift==1) printf("sea level trend computation,\n\n");
 fclose(ini);

 search=strrchr(conf->program_name,'/');
 if(search!=NULL) sprintf(conf->program_name,"%s",search+sizeof(char));

 conf->justfilename=strrchr(conf->data_filename,'/');
  if(conf->justfilename==NULL) conf->justfilename=conf->data_filename;
  else conf->justfilename+=sizeof(char);

  search=strrchr(conf->justfilename,'.');
  nsearch=strlen(conf->justfilename)-strlen(search);
  strncpy(conf->rootname,conf->justfilename,nsearch);
  conf->rootname[nsearch]='\0';

  
}



/*----------------------------------------------------------------------------*/

int WaveSetHBuffer(char *name, char *h_models, char *ext, int longperiod,grid_t *grid_OT, fcomplex **buffer_OT, fcomplex *mask_OT)

/*----------------------------------------------------------------------------*/
{
  int na_XX,status_OT;
  char input[260];
  fcomplex *buffer,mask;

  status_OT=-1;
  if (!(longperiod))  sprintf(input,"%s/%s.%s",h_models,name,ext);
  else sprintf(input,"%s/%s.ext",h_models,name,ext);

  bmg_loadc1_analyse_alti(input,1,1,1,grid_OT,&na_XX,buffer_OT,&mask,&status_OT);
/*   found_z(grid_OT[i].nx,grid_OT[i],lat,lon,buffer_OT[i],mask_XX,&z_XX); */


  *mask_OT=mask;

  return(status_OT);
  
}



/*----------------------------------------------------------------------------*/

void loadingEC(char *input,double lat, double lon, fcomplex *Z,fcomplex *mask)

/*----------------------------------------------------------------------------*/
{

grid_t grid_EC;
fcomplex *buffer_EC=NULL,z_EC,mask_EC;
int na_EC,status;



bmg_loadc1_analyse_alti(input,1,1,1,&grid_EC,&na_EC,&buffer_EC,mask,&status);

 if(status !=0) {
     *Z=*mask;
     return;
   }

mask_EC=*mask;

 if(lon>180) lon-=360;
/*faire un test sur lon entre -180 et +180*/

found_z(na_EC,grid_EC,lat,lon,buffer_EC,mask_EC,&z_EC);

/*  *mask=mask_EC; */

 free(buffer_EC);

*Z=z_EC;


}/*end*/



/*----------------------------------------------------------------------------*/

int Wave_select(int i, tidal_wave *carac,int detided)

/*----------------------------------------------------------------------------*/
{

  int verbose;


verbose=1;

switch (i)
{

 case 1 : {*carac=wM2; break;}
 case 2 : {*carac=wS2; break;}
 case 3 : {*carac=wN2; break;}
 case 4 : {*carac=wK2; break;}
 case 5 : {*carac=wK1; break;}
 case 6 : {*carac=wO1; break;}
 case 7 : {*carac=wQ1; break;}
 case 8 : {*carac=w2N2; break;}
/*----------------------------------------------------------------------------*/
 case 9 : {*carac=wMu2; break;}
 case 10 : { *carac=wNu2;break;}
 case 11 : {*carac=wL2; break;}
 case 12 : {*carac=wT2; break;}
 case 13 : {*carac=wLa2; break;}
 case 14 :
   {
     if (detided==2) {*carac=wKJ2; break;}
     else {verbose=0;break;}
   }
 case 15 :
   {
     if (detided==2) {*carac=wR2; break;}
     else {verbose=0;break;}
   }
/*----------------------------------------------------------------------------*/
 case 16 : {*carac=wP1; break;}
 case 17 : {*carac=wOO1; break;}
 case 18 : {*carac=wJ1; break;}
 case 19 : {*carac=wPhi1; break;}
 case 20 : {*carac=wPi1; break;}
 case 21 :
  {
    if (detided==0) {verbose=0;break; }
    else {*carac=wPsi1; break;}
  }
 case 22 : {*carac=wRo1; break;}
 case 23 : {*carac=wSig1; break;}
 case 24 : {*carac=wTta1; break;}
 case 25 : {*carac=w2Q1; break;}
 case 26 : {*carac=wKi1; break;}
/*----------------------------------------------------------------------------*/
 case 27 : {*carac=wMf; break;}
 case 28 : {*carac=wMm; break;}
 case 29 : {*carac=wMtm; break;}
 case 30 :
   {
     if (detided==2) {*carac=wMSqm; break;}
     else {verbose=0;break;}
   }
 case 31 : {*carac=wSsa; break;}
 case 32 : {*carac=wSa; break;}
 case 33 :
   {
     if (detided==2) {*carac=wMSm; break;}
     else {verbose=0;break;}
   }
 case 34 : {*carac=wMSf; break;}
 case 35 : {*carac=wMqm; break;}
 case 36 : {*carac=wMStm; break;}

/*----------------------------------------------------------------------------*/
 case 37 : {*carac=wM4; break;}
 case 38 : {*carac=wMS4; break;}
 case 39 : {*carac=wMN4; break;}
 case 40 : {*carac=wS4; break;}
 case 41 : {*carac=wN4; break;}

/*----------------------------------------------------------------------------*/
 case 42 : {*carac=wS1; break;}
 case 50 : {*carac=wM1; break;}

 default : verbose=0;}/*switch*/


return(verbose);
}


/*-------------------------------------------------------------------------------*/

int questionmixte(int *position , wavelist S  , typeA *t , wavelist *n)

/*-------------------------------------------------------------------------------*/
{

  int i,type,cnt;
  int *keep;
  int finA=0,finb=0;
  int ctA=0,ctB=0,tmp=0;

  keep=(int *)malloc(S.n*sizeof(int));
  cnt=0;

      for (i=0;i<S.n;i++) {
          if (S.wave[i].Ap!=0) {
              printf("choix du type d'analyse de l'onde %s (1=ortho   0=harmo)\n",S.wave[i].name);
              scanf("%d",&type);
              if (type==0) keep[i]=0;
              else keep[i]=1;
            }
          else keep[i]=0;
          
        }/*for i*/

      for (i=0;i<S.n;i++) if (keep[i]==1){if (S.wave[i].nT==1)ctA++;if (S.wave[i].nT==2)ctB++;}
      if (ctA>0) {tmp++;}
      if (ctB>0) {tmp++;}
          


  for (i=0;i<S.n;i++) {
      if (keep[i]==0) {
          position[cnt]=i;
          n->wave[cnt]=S.wave[i];
          cnt++;
        }
   }

  t->h=cnt;

  for (i=0;i<S.n;i++) {
      if (keep[i]!=0) {
          position[cnt]=i;
          n->wave[cnt]=S.wave[i];
          cnt++;
        }
   }
  t->o=cnt-t->h;
  free(keep);

return(tmp);

}


/*----------------------------------------------------------------------------*/

void Path4Model(int i,char *path,char *ext)

/*----------------------------------------------------------------------------*/
{
  char FES99[80], FES02[80],GOT[80],MOG2D_archive[80],MOG2D[80];
  char ext_FES99[80],ext_FES02[80],ext_lp_FES02[80],ext_GOT[80],ext_MOG2D_archive[80],ext_MOG2D[80];
 /*
 sprintf(MOG2D_new,"/data/nivmer/travail_en_cours/simulations/medsea/bench-tide/december-2002");
 sprintf(ext_MOG2D_new,"s2c.den");*/
 sprintf(MOG2D,"/calcul/maewo/roblou/tides-mog2d");
 sprintf(ext_MOG2D,"bmg");

 sprintf(MOG2D_archive,"/data/nivmer/produits/archives-medsea/test42");
 sprintf(ext_MOG2D_archive,"den.new");

 sprintf(FES02,"/data/nivmer/produits/tides-reorganized/legos/FES2002/elevation");
 sprintf(ext_FES02,"den.a");
 sprintf(ext_lp_FES02,"den.h");

 sprintf(FES99,"/data/nivmer/produits/tides-reorganized/legos/FES98-99/elevation");
 sprintf(ext_FES99,"den.3");

 sprintf(GOT,"/data/nivmer/produits/tides-reorganized/others/GOT00/elevation");
 sprintf(ext_GOT,"den");

 switch(i)
   {
   case 1:
     { sprintf(path,MOG2D);
     sprintf(ext,ext_MOG2D);
     break; };/* M2 */
   case 2:
     { sprintf(path,MOG2D);
     sprintf(ext,ext_MOG2D);
     break;};/* S2 */
   case 3:
     { sprintf(path,MOG2D);
     sprintf(ext,ext_MOG2D);
     break;}; /* N2 */
   case 4:
     { sprintf(path,MOG2D);
     sprintf(ext,ext_MOG2D);
     break;}; /* K2 */
   case 5:
     { sprintf(path,FES02);
     sprintf(ext,ext_FES02);
     break;}; /* K1 */
   case 6:
     { sprintf(path,GOT);
     sprintf(ext,ext_GOT);
     break;}; /* 01 */
   case 7:
     { sprintf(path,MOG2D);
     sprintf(ext,ext_MOG2D);
     break;}; /* Q1 */
   case 8:
     { sprintf(path,MOG2D);
     sprintf(ext,ext_MOG2D);
     break;}; /* 2N2 */
   case 16:
     { sprintf(path,FES02);
     sprintf(ext,ext_FES02);
     break;}; /* P1 */
   case 27:
     { sprintf(path,MOG2D);
     sprintf(ext,ext_MOG2D);
     break;}; /* Mf */
   case 28:
     { sprintf(path,MOG2D);
     sprintf(ext,ext_MOG2D);
     break;}; /* Mm */
   case 29:
     { sprintf(path,MOG2D);
     sprintf(ext,ext_MOG2D);
     break;}; /* Mtm */
   case 30:
     { sprintf(path,MOG2D);
     sprintf(ext,ext_MOG2D);
     break;}; /* MSqm */
   case 31:
     {sprintf(path,MOG2D);
     sprintf(ext,ext_MOG2D);
     break;}; /* Ssa */
   case 37:
     { sprintf(path,MOG2D);
     sprintf(ext,ext_MOG2D);
     break;}; /* M4 */
   case 42:
     { sprintf(path,MOG2D);
     sprintf(ext,ext_MOG2D);
     break;}; /* S1 */
   }/* end switch*/

}/*end void*/


