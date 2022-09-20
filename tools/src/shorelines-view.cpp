
/**************************************************************************

  T-UGO tools, 2006-2009

  Unstructured Ocean Grid initiative

Contributors:

  Florent Lyard      LEGOS/CNRS, Toulouse, France
  Yoann Le Bars      LEGOS, Toulouse, France (PhD)
  Laurent Roblou     LEGOS, Toulouse, France
  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
  Frédéric Dupont    Université de Laval à Québec, Canada

E-mail: florent.lyard@legos.obs-mip.fr

***************************************************************************/

#include "pocview_prototype.hpp"

/*les fichiers de cotes comme les polygones sont stockes sous la forme de chaine
  cette methode tres simple permet de passer outre les allocations */
/*cette methode permet de plus, de facilite les editions (add remove ...) de donnees*/
/*la gestion de chaine est toutefois tres lente : il est donc plus rapide d'ajouter un structure qui regroupe la chaine
  et une fonction qui regenere la structure a volonte en fonction de l'etat de la chaine*/
/*cette gestion separee pourra etre utilise pour l'ajout de la commande UNDO dans l'edition des polygones*/

/*j'ai trouve le moyen de gerer les listes plus vite
  je repasse donc tout en liste*/


#define swabi2(i2) (((i2) >> 8) + (((i2) & 255) << 8))
#define swabi4(i4) (((i4) >> 24) + (((i4) >> 8) & 65280) + (((i4) & 65280) << 8) + (((i4) & 255) << 24))

/* ############################################### */
/* #                                              #*/
/* #                                              #*/
/* # identifie le format des fichiers de cotes    #*/
/* #        et envoie le chargement               #*/
/* #                                              #*/
/* ############################################## */

int loading_coast_file(coast_struct *coast, char *filename)
{

  char *extention;
  int rien;

  //check the filename
  extention=(char *)calloc(52,sizeof(char));
  extention=strrchr(filename,'.');

  /*if(strcmp(extention,".cst")==0 )  loading_cst_format_coast_file(filename) ;*/
  if(strcmp(extention,".cst")==0 )  rien=load_coast_file_cst_format(coast,filename) ;
  /*if(strcmp(extention,".ngdc")==0 )  loading_ngdc_format_coast_file(filename) ;  */
  if(strcmp(extention,".ngdc")==0 ) rien=load_coast_file_ngdc_format(coast,filename);
  /*l'appel de cette function permet de mettre les donnes du fichier de cotes dans une liste et non pas dans une structure*/
  /*l'ordre des segments est inverse en raison de l'appel a la commande prepend */
  //free(extention);

  //  from_list_to_struct(coast->plg_list, &coast->plg, &rien);


  return(coast_min_max_detection(coast) );

}

/* ##################################### */
/* #                                   # */
/* #   This routine save coast line    #*/
/* #       in cst fortran format       #*/
/* # Thierry LETELLIER le 02May2005    # */
/* #                                   # */
/* ##################################### */

int save_coast_in_fortran_cst_format(GtkWidget *appel)
{

  int size_to_fortran;
  FILE *coast_file;
  int seg_ID;
  int pt_seg;
  int i,j;
  Fposition_struct point;
  char *filename;
  int rstatus;
  int value_3;
  plg_t *plg;
  int nb_plg;
  int record;
  int first,last;
  GList *plg_list;
  map_C *map;
  int doc_ID;
  double lon_mem;

  filename=(char *)calloc(256,sizeof(char));

  rstatus=file_save_selection(NULL,filename);
  coast_file=fopen(filename,"wb");

  doc_ID=map_identification(&map);

  if (coast_file==NULL) {printf("#Coastline File not found : %s",filename);return(0);}

  plg_list=g_list_nth(map->coast->plg_list,0);
  record=1;
  while(plg_list!=NULL) {
      plg=(plg_t *)g_list_nth_data(plg_list,0);
      plg_list=g_list_nth(plg_list,1);
          size_to_fortran=3*sizeof(int);
          seg_ID= i+1;
          value_3=record;
          pt_seg=plg->npt+record-1;

 #if NEED_SWAP == 1
         swap_bytes(&size_to_fortran,SIZEOF_INT) ;
         swap_bytes(&seg_ID,SIZEOF_INT) ;
         swap_bytes(&value_3,SIZEOF_INT) ;
         swap_bytes(&pt_seg,SIZEOF_INT) ;
 #endif
 
          fwrite(&size_to_fortran , sizeof(int) , 1 , coast_file);
          fwrite(&seg_ID, sizeof(int) , 1 , coast_file);
          fwrite(&value_3 , sizeof(int) , 1 , coast_file);
          fwrite(&pt_seg, sizeof(int) , 1 , coast_file);
          fwrite(&size_to_fortran , sizeof(int) , 1 , coast_file);

          size_to_fortran=plg->npt*sizeof(Fposition_struct);
          
 #ifdef NEED_SWAP	
         swap_bytes(&size_to_fortran,SIZEOF_INT) ;
 #endif
 
          fwrite(&size_to_fortran , sizeof(int) , 1 , coast_file);
    
          for(j=0;j<plg->npt;j++) {
              point.lon=(float)plg->x[j];
              lon_mem=  plg->x[j];
              point.lat=(float)plg->y[j];
              
 #ifdef NEED_SWAP	
         swap_bytes(&point.lon,SIZEOF_FLOAT) ;
         swap_bytes(&point.lat,SIZEOF_FLOAT) ;
 #endif
              fwrite(&point , sizeof(Fposition_struct) , 1 , coast_file);
            }

          fwrite(&size_to_fortran , sizeof(int) , 1 , coast_file);
          record+=plg->npt;
    }
  fclose(coast_file);

}







/* ############################################### */
/* #                                              #*/
/* #                                              #*/
/* #        detection du cadre des cotes         #*/
/* #       Thierry LETELLIER 26Jan2005           #*/
/* #                                              #*/
/* #                                              #*/
/* ############################################## */
int coast_min_max_detection(coast_struct *coast)
{
  int i,j,n;
  int nb_plg;
  GList *accel;
  plg_t *plg;
  double x_max,x_min,y_max,y_min;

  if( coast->GeoCadre==NULL) coast->GeoCadre=new cadre_C;

  y_max=-100;
  y_min=+100;
  x_max=-200;
  x_min=+380;

  nb_plg=g_list_length(coast->plg_list);
  accel=g_list_nth(coast->plg_list,0);

  while(accel!=NULL)
    {
      plg=(plg_t *)g_list_nth_data(accel,0);
      accel=g_list_nth(accel,1);
          for(n=0;n<plg->npt;n++) {
              y_max=MAX(y_max,plg->y[n]);
              y_min=MIN(y_min,plg->y[n]);
              x_max=MAX(x_max,plg->x[n]);
              x_min=MIN(x_min,plg->x[n]);
              coast->max_nb_pt++;
            }
    }
  coast->GeoCadre->implement(x_min,x_max,y_min,y_max);
  return(1);
}

/* ############################################### */
/* #                                              #*/
/* #                                              #*/
/* # lecture du fichier cotes NGDC dans une liste #*/
/* #       Thierry LETELLIER 09May2005            #*/
/* #                                              #*/
/* #                                              #*/
/* ############################################## */
int load_coast_file_ngdc_format(coast_struct *coast,char *filename)
{

/*cette function devait permettre de lire le fichier de cote directement dans une liste
et ainsi permettre une recherche, un ajout, retrait de donnes plus rapide
l'appel a la commande prepend inverse l'ordre des segments dans la liste*/


struct GSHHS {	/* Global Self-consistant Hierarchical High-resolution Shorelines */
        int id;				/* Unique polygon id number, starting at 0 */
        int n;				/* Number of points in this polygon */
        int level;			/* 1 land, 2 lake, 3 island_in_lake, 4 pond_in_island_in_lake */
        int west, east, south, north;	/* min/max extent in micro-degrees */
        int area;			/* Area of polygon in 1/10 km^2 */
        short int greenwich;		/* Greenwich is 1 if Greenwich is crossed */
        short int source;		/* 0 = CIA WDBII, 1 = WVS */
};

 double w, e, s, n, area;
 char source;
 FILE	*fp;
 int	k, max_east = 270000000;
 GdkPoint p;
 struct GSHHS h;
 plg_t *plg;
 int ok;

 coast->plg_list==NULL;

 if ((fp = fopen (filename, "rb")) == NULL ) {
   __ERR_BASE_LINE__( "gshhs:  Could not find file %s.\n", filename);
   exit (EXIT_FAILURE);
 }

 while (fread(&h, sizeof (struct GSHHS), 1, fp) == 1)
   {
#ifdef NEED_SWAP
   h.id = swabi4 ((unsigned int)h.id);
   h.n = swabi4 ((unsigned int)h.n);
   h.level = swabi4 ((unsigned int)h.level);
   h.west = swabi4 ((unsigned int)h.west);
   h.east = swabi4 ((unsigned int)h.east);
   h.south = swabi4 ((unsigned int)h.south);
   h.north = swabi4 ((unsigned int)h.north);
   h.area = swabi4 ((unsigned int)h.area);
   h.greenwich = swabi2 ((unsigned int)h.greenwich);
   h.source = swabi2 ((unsigned int)h.source);
#endif
   w = h.west  * 1.0e-6;
   e = h.east  * 1.0e-6;
   s = h.south * 1.0e-6;
   n = h.north * 1.0e-6;
   source = (h.source == 1) ? 'W' : 'C';
   area = 0.1 * h.area;

   plg=(plg_t *)calloc(1,sizeof(plg_t));
   plg->npt=h.n;
   plg->x=(double *)calloc(h.n,sizeof(double));
   plg->y=(double *)calloc(h.n,sizeof(double));

   ok=0;

   for (k = 0; k < h.n; k++) {
     if (fread (&p, sizeof(GdkPoint), 1, fp) != 1) {
       __ERR_BASE_LINE__( "gshhs:  Error reading file %s for polygon %d, point %d.\n", filename, h.id, k);
       exit (EXIT_FAILURE);
       }
#ifdef NEED_SWAP
     p.x = swabi4 ((unsigned int)p.x);
     p.y = swabi4 ((unsigned int)p.y);
#endif
     plg->id=k;
     plg->x[k] = (h.greenwich && p.x > max_east) ? p.x * 1.0e-6 - 360.0 : p.x * 1.0e-6;
     plg->y[k] = p.y * 1.0e-6;
     if(plg->x[k]>FRAME_xmax)plg->x[k]-=360;
     else if(plg->x[k]<FRAME_xmin)plg->x[k]+=360;
     }

   /*Ajout temporaire !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
   /*    if( max_east != 180000000)coast->plg_list=g_list_prepend(coast->plg_list,plg); */
   //if( plg->npt>120000)coast->plg_list=g_list_prepend(coast->plg_list,plg);
   //if((w<180)&&(n<60) &&(s>-65) &&(e>70)&&(plg->npt>10)&&(h.level==1) )coast->plg_list=g_list_prepend(coast->plg_list,plg);
   /*   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

    coast->plg_list=g_list_prepend(coast->plg_list,plg);
    plg=NULL;
    max_east = 180000000;	/* Only Eurasiafrica needs 270 */
    }

  fclose (fp);
}


/* ############################################### */
/* #                                              #*/
/* #                                              #*/
/* #        lecture du ficher coast              #*/
/* #          format SUN Fortran                  #*/
/* #       Thierry LETELLIER 26Jan2005           #*/
/* #                                              #*/
/* #                                              #*/
/* ############################################## */
int load_coast_file_cst_format(coast_struct *coast,char *filename)
{
  int n,rstatus;
  FILE *coast_file;
  float lat,lon;
  int i,j;
  int rien;
  float point[2];
  int sizetoread,seg_id,value_3,value_4,readsize,sizetoseek,seekedsize;

  plg_t *plg;

  coast_file=fopen(filename,"rb");
  if (coast_file==NULL) {
    printf("#Coastline File not found : %s",filename);
    return(0);
    }
  i=0;

  while(!feof(coast_file)) {
    plg=(plg_t *)calloc(1,sizeof(plg_t));
    n=fread(&sizetoread,sizeof(int),1,coast_file);/* not necesarry to swap*/
    if(n==1) {
      fread(&seg_id,sizeof(int),1,coast_file);
#ifdef NEED_SWAP	
      swap_bytes(&seg_id,SIZEOF_INT) ;
#endif
      plg->id=i;
      fread(&value_3,sizeof(int),1,coast_file);
#ifdef NEED_SWAP	
      swap_bytes(&value_3,SIZEOF_INT) ;
#endif
      n=fread(&value_4,sizeof(int),1,coast_file);/* not necesarry to swap*/
      fread(&readsize,sizeof(int),1,coast_file); /* not necesarry to swap*/
      fread(&sizetoseek,sizeof(int),1,coast_file);
#ifdef NEED_SWAP	
      swap_bytes(&sizetoseek,SIZEOF_INT) ;
#endif
      plg->npt=sizetoseek/sizeof(float)/2;
      plg->x=(double *)calloc(plg->npt,sizeof(double));
      plg->y=(double *)calloc(plg->npt,sizeof(double));
      for(j=0;j<plg->npt;j++) {
        fread(&point,2*sizeof(float),1,coast_file);
#ifdef NEED_SWAP	
        swap_bytes(&(point[0]),SIZEOF_FLOAT) ;
        swap_bytes(&(point[1]),SIZEOF_FLOAT) ;
// 	  #if SIZEOF_FLOAT == 8
// 	  point[0]=fswap_64(point[0]);
// 	  #endif
// 	  #if SIZEOF_FLOAT == 4
// 	  point[0]=fswap_32(point[0]);
// 	  #endif
// 	  #if SIZEOF_FLOAT == 2
// 	  point[0]=fswap_16(point[0]);
// 	  #endif

// 	  #if SIZEOF_FLOAT == 8
// 	  point[1]=fswap_64(point[1]);
// 	  #endif
// 	  #if SIZEOF_FLOAT == 4
// 	  point[1]=fswap_32(point[1]);
// 	  #endif
// 	  #if SIZEOF_FLOAT == 2
// 	  point[1]=fswap_16(point[1]);
// 	  #endif
#endif
        plg->x[j]=point[0];
        plg->y[j]=point[1];
        if(plg->x[j]>FRAME_xmax)plg->x[j]-=360;
        else if(plg->x[j]<FRAME_xmin)plg->x[j]+=360;
        }
      fread(&seekedsize,sizeof(int),1,coast_file);
      coast->plg_list=g_list_prepend(coast->plg_list,plg);
      i++;
      }
   }
   fclose(coast_file);
}
