
#define  LEGEND_MAIN
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
#include "legend.def"

#include "tools-structures.h"

#include "legend.h"

#include "geo.h"
#include "map.h"

/*----------------------------------------------------------------------------*/
/* Fonctions appelees depuis Fortran */
/*----------------------------------------------------------------------------*/
#ifdef _add_

extern void lgd_free_ ();
#pragma weak lgd_free_ = lgd_free

#endif
/*----------------------------------------------------------------------------*/
int DefaultColumn=0;

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void lgd_free(int *rstatus)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/*----------------------------------------------------------------------------*/
  int	error;
/*----------------------------------------------------------------------------*/

if(legendstab != NULL)
  free(legendstab);

lgd_nmax= 0;

return;

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void lgd_init(int *rstatus)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{

  if(legendstab != NULL) free(legendstab);
  lgd_n=0;
  lgd_nmax=0;
  lgd_lastedited=0;
  lgd_defaulttype=LGD_TEXT_SYMBOL;
}


/*-----------------------------------------------------------------------
 Allocate memory for legend
 -----------------------------------------------------------------------*/
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void lgd_alloc(int n,int *rstatus)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int error;

  lgd_free(rstatus);

  if(n != 0) {
    lgd_nmax=n;
    exitIfNull(
      legendstab=(legend_t *) malloc(n*sizeof(legend_t))
      );
  }

  rstatus=0;//LGD_STATUS_OK;
}


/*-----------------------------------------------------------------------
 ReAllocate memory for legend
 -----------------------------------------------------------------------*/
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void lgd_realloc(int n,int *rstatus)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
int	m,error,k;
legend_t *newtab=NULL;
bool ok;

if(n < lgd_nmax) return;
if(n == 0) return;

m=MAX(lgd_nmax+10,n);

if(legendstab == NULL)
{
  exitIfNull(
    legendstab=(legend_t *) malloc(m*sizeof(legend_t))
    );
  lgd_nmax=m;
  return;
}

exitIfNull(
  newtab=(legend_t *)  malloc(lgd_nmax*sizeof(legend_t))
  );
for(k=0;k<lgd_nmax;k++)
  newtab[k]=legendstab[k];

free(legendstab);

exitIfNull(
  legendstab=(legend_t *) malloc(m*sizeof(legend_t))
  );
for(k=0;k<lgd_nmax;k++)
  newtab[k]=legendstab[k];
lgd_nmax = m;

*rstatus=LGD_STATUS_OK;

return;
}

/*-----------------------------------------------------------------------
 Clear all
 -----------------------------------------------------------------------*/
void lgd_clearall(int *rstatus)
  {
  return;
  }

/*-----------------------------------------------------------------------
 Transform to Cartesian
-----------------------------------------------------------------------*/
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void lgd_tocartesian01(legend01_t *legend,int *rstatus)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
int m,dum,n;
int ok;
float t,p,x,y;

  t=legend->T;
  p=legend->P;
//  geo_recale(&t,&xenv_min,&xenv_max);
  legend->X=x;
  legend->Y=y;

  return;
}

/*-----------------------------------------------------------------------
 Transform to Cartesian
-----------------------------------------------------------------------*/
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void lgd_tocartesian(legend_t *tab,int *ntab, int *rstatus)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
int m,n;
int ok;
float t,p,x,y;
legend01_t *dum1=NULL;
legend02_t *dum2=NULL;

for(m=0;m<*ntab;m++) {
  switch (tab[m].Type) {
    case LGD_TEXT_SYMBOL:
    case LGD_SYMBOL:
    dum1=(legend01_t *) tab[m].ptr;
    t=dum1->T;
    p=dum1->P;
    tab[m].T=t;
    tab[m].P=p;
//    ok=geo_recale(&t,&xenv_min,&xenv_max);
    tab[m].X=x;
    tab[m].Y=y;
    dum1->X=x;
    dum1->Y=y;
    break;

    case LGD_GRAPH:
    break;
    }
  }


return;
}

/*-----------------------------------------------------------------------
 Save legend
 -----------------------------------------------------------------------*/
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void lgd_save01(FILE *fp, legend01_t *ptr)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int	 n, l;
  legend01_t legend;

  legend=*ptr;

  if(legend.Type == LGD_TEXT_SYMBOL) {
    fprintf(fp,"#TEXT&SYMBOL\n");
    fprintf(fp,"%f %f \n",legend.T, legend.P);
    fprintf(fp,"%f %f \n",legend.X, legend.Y);
    l=strlen(legend.Text);
    fprintf(fp,"%s\n",legend.Text);
    fprintf(fp,"%d %d %d \n",legend.Font,
                             legend.FontType,
                             legend.FontSize);
    fprintf(fp,"%d %d %d \n",legend.Orientation,
                             legend.Hjustify,
                             legend.Vjustify);
    fprintf(fp,"%d %d %d \n",legend.TextFG,
                             legend.TextBG,
                             legend.TextFramed);
    fprintf(fp,"%d \n",legend.Symbol);
    fprintf(fp,"%f \n",legend.SymbolSize);
    fprintf(fp,"%d %d \n",legend.PenColour, legend.FillColour);
    }
  else if(legend.Type == LGD_SYMBOL  ) {
    fprintf(fp,"#SYMBOL\n");
    fprintf(fp,"%f %f %d %f %d %d\n",legend.T,legend.P,
                                     legend.Symbol,
                                     legend.SymbolSize,
                                     legend.PenColour,
                                     legend.FillColour);
    }

return;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int lgd_save02(FILE *fp, legend02_t *legend,char *fmt, char *comments)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,l,n,status;

  fprintf(fp,"#GRAPH-XYZ\n");
  fprintf(fp,"%d %d %d",legend->ID,legend->np,legend->nz);
  fprintf(fp," %d\n",legend->pen1);
  for (n=0; n < legend->np;n++) {
    fprintf(fp,"%6d %12.3f %12.3f",legend->points[n].N,legend->points[n].t,legend->points[n].p);
    for(i=0;i<legend->nz;i++) fprintf(fp," %12f",legend->points[n].z[i]);
    fprintf(fp,"\n");
    }

  status=LGD_STATUS_OK;
  return(status);

  end2 :
    fclose(fp);
    status=LGD_STATUS_READ_ERROR;
    return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int lgd_save03(FILE *fp, legend03_t *legend,char *fmt, char *comments)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int	i, n, l, status;

  fprintf(fp,"#PROFILE\n");
  if(fprintf(fp,"%d %d %d\n",legend->ID,legend->ndepths,legend->ncolumns) ==0) goto end2;
//  if(fprintf(fp,"%d\n",legend->pen1) ==0) goto end2;
  if(fprintf(fp,"%f %f\n",legend->points[0].t,legend->points[0].p) ==0) goto end2;
  for (n=0; n < legend->ndepths;n++) {
    if(fprintf(fp,"%f",legend->points[0].depths[n]) ==0) goto end2;
    for(i=0;i<legend->ncolumns;i++) if(fprintf(fp," %f",legend->points[0].z[n][i]) ==0) goto end2;
    if(fprintf(fp,"\n") ==0) goto end2;
    }

  status=LGD_STATUS_OK;
  return(status);

  end2 :
    fclose(fp);
    status=LGD_STATUS_READ_ERROR;
    return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int lgd_save(const char *out, legend_t *legend, int nlegend, char *fmt, char *comments)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
int   n, l;
FILE *fp=NULL;
int   status;

 if ((fp = fopen(out,"w")) == NULL) {
   __ERR_BASE_LINE__("Probleme a l ouverture du fichier : %s \n", out);
   exit(-1);
   }

fprintf(fp,"%d\n",nlegend);

for(n=0;n<nlegend;n++) {
  if(legend[n].Type == LGD_TEXT_SYMBOL) {
    lgd_save01(fp,(legend01_t *) legend[n].ptr);
    }
  else if(legend[n].Type == LGD_SYMBOL) {
    lgd_save01(fp,(legend01_t *) legend[n].ptr);
    }
  else if(legend[n].Type == LGD_GRAPH) {
    status=lgd_save02(fp,(legend02_t *) legend[n].ptr,fmt,comments);
    }
  else if(legend[n].Type == LGD_PROFILE) {
    status=lgd_save03(fp,(legend03_t *) legend[n].ptr,fmt,comments);
    }
  }

if(comments!=0) {
  fprintf(fp,"%s",comments);
  }

fclose(fp);

return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void lgd_savedefault(char *fich, int *rstatus)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{

*rstatus=lgd_save(fich, legendstab, lgd_n, NULL,NULL);
return;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void lgd_copystyle01(legend01_t *legend)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
int n;
legend01_t *dum=NULL;

for(n=0;n<lgd_n;n++) {
  dum=(legend01_t *) legendstab[n].ptr;
  dum->Font=legend->Font;
  dum->FontType=legend->FontType;
  dum->FontSize=legend->FontSize;
  dum->TextFG=legend->TextFG;
  dum->TextBG=legend->TextBG;
  dum->TextFramed=legend->TextFramed;
  dum->Symbol=legend->Symbol;
  dum->SymbolSize=legend->SymbolSize;
  dum->PenColour=legend->PenColour;
  dum->FillColour=legend->FillColour;
  }
return;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void lgd_copystyle02(legend02_t *legend)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
int n;
legend02_t *dum=NULL;

for(n=0;n<lgd_n;n++) {
  dum=(legend02_t *) legendstab[n].ptr;
  dum->pen0=legend->pen0;
  dum->pen1=legend->pen1;
  dum->currentz=legend->currentz;
  }
return;

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void lgd_defaultlegend01(legend01_t *legend,int *rstatus)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/*----------------------------------------------------------------------------*/
int    length,i,j;
int    ok;
float  t,p,x,y,z;
/*----------------------------------------------------------------------------*/

  legend->Type=0;
  t=legend->T;
  p=legend->P;
  legend->X=x;
  legend->Y=y;
  sprintf(legend->Text,"NEW LEGEND");
  legend->Font=0;
  legend->FontType=0;
  legend->FontSize=12;
  legend->Orientation=0;
  legend->Hjustify=0;
  legend->Vjustify=0;
  legend->Symbol=1;
  legend->SymbolSize=1.0;
  legend->TextFG=1;
  legend->TextBG=0;
  legend->TextFramed=0;
  legend->PenColour=0;
  legend->FillColour=1;

  return;

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void lgd_defaultlegend(legend_t *legend,int *rstatus)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
int    length,i,j;
int    ok;
float  t,p,x,y,z;
legend01_t *dum1=NULL;
legend02_t *dum2=NULL;


switch (legend->Type)
  {
  case LGD_TEXT_SYMBOL:
  legend->ptr=(char *) new legend_t;
  dum1=(legend01_t *) legend->ptr;
  dum1->T=legend->T;
  dum1->P=legend->P;
  lgd_defaultlegend01((legend01_t *) legend->ptr,rstatus);
  break;

  case LGD_SYMBOL:
  legend->ptr=(char *) new legend_t;
/*  legend->ptr->T=legend->T;
  legend->ptr->P=legend->P;*/
  lgd_defaultlegend01((legend01_t *) legend->ptr,rstatus);
  break;

  case LGD_GRAPH:
  legend->ptr=(char *) new legend_t;
  dum2=(legend02_t *) legend->ptr;
  dum2->currentz=0;
  dum2->pen0=0;
  dum2->pen1=0;
  break;
  }
return;

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int lgd_read01(FILE *fp, legend01_t *legend)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int	n,status;
  int 	ok;
  float x,y;
  char *a=NULL,c;

  if (legend->Type==LGD_TEXT_SYMBOL) {
    if(fscanf(fp,"%f %f",&legend->T,&legend->P) != 2)
      goto end2;
    if(fscanf(fp,"%f %f",&legend->X,&legend->Y) != 2)
      goto end2;
    lgd_tocartesian01(legend, &status);
    do {c=fgetc(fp);}  while ((c != '\n') && !feof(fp));
    a=legend->Text;
    do {*a=fgetc(fp); a++;}  while ((*(a-1) != '\n') && !feof(fp));
    *(a-1)=0;
/*
    if(fscanf(fp,"%s",   legend->Text) != 1)
      goto end2;
*/
    legend->Font=0;
    legend->FontType=0;
    if(fscanf(fp,"%d %d %d",&legend->Font,
                            &legend->FontType,
                            &legend->FontSize) != 3)
      goto end2;
    if(fscanf(fp,"%d %d %d",&legend->Orientation,
                            &legend->Hjustify,
                            &legend->Vjustify) != 3)
      goto end2;
    if(fscanf(fp,"%d %d %d",&legend->TextFG,
                            &legend->TextBG,
                            &legend->TextFramed) != 3)
      goto end2;
    if(fscanf(fp,"%d ",  &legend->Symbol) != 1)
      goto end2;
    if(fscanf(fp,"%f ",  &legend->SymbolSize) != 1)
      goto end2;
    if(fscanf(fp,"%d %d",&legend->PenColour,&legend->FillColour) != 2)
      goto end2;
    }
  else if (legend->Type==LGD_SYMBOL) {
    if(fscanf(fp,"%f %f %d %f %d %d", &legend->T,&legend->P,
                                      &legend->Symbol,&legend->SymbolSize,
                                      &legend->PenColour,&legend->FillColour) != 6)
      goto end2;
                                
    legend->X=x;
    legend->Y=y;
    strcpy(legend->Text,"\0");
    }

  return(0);

  end2 :
    fclose(fp);
    return(LGD_STATUS_READ_ERROR);

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void lgd_read02(FILE *fp, legend02_t *legend, int *rstatus)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
int	i,k,n;
int 	ok;
float 	x,y;
int renumber=0;

  if(fscanf(fp,"%d %d %d",&legend->ID,&legend->np,&legend->nz) != 3)
    goto end2;
  if(legend->np < 0) goto end2;
  if(legend->nz < 0) goto end2;
  if(fscanf(fp,"%d",&legend->pen1) != 1)
    goto end2;
  legend->points=new point_t[legend->np*sizeof(point_t)];
  for (n=0; n < legend->np;n++) {
    if(fscanf(fp,"%d %f %f",&legend->points[n].N,&legend->points[n].t,
                                                 &legend->points[n].p) != 3)
      goto end2;
    legend->points[n].z=new float[legend->nz];
    for(i=0;i<legend->nz;i++) if(fscanf(fp,"%f",&legend->points[n].z[i]) != 1)
      goto end2;
    }

  legend->pen0=0;
//  legend->framed=0;
  legend->currentz=DefaultColumn;
  legend->text=0;

/* *----------------------------------------------------------------------------
  check for weird-formed node number*/
  for (n=1; n < legend->np;n++) {
    if(legend->points[n].N==legend->points[n-1].N) {
      renumber=1;
      break;
      }
    }
/* *----------------------------------------------------------------------------
  if necessary, apply re-numbering*/
  if(renumber==1) {
    printf("warning: node numbering is deficient, xscan will apply patch procedure\n");
    for (n=0; n < legend->np;n++) {
      legend->points[n].N+=n;
      }
    }

  *rstatus=LGD_STATUS_OK;
  return;

  end2 :
    fclose(fp);
    *rstatus=LGD_STATUS_READ_ERROR;
    return;

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void lgd_read03(FILE *fp, legend02_t *legend, int *rstatus)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
int	i,k,n;
int 	ok;
float 	x,y;

  if(fscanf(fp,"%d %d %d",&legend->ID,&legend->np,&legend->nz) != 3)
    goto end2;
  if(legend->np < 0) goto end2;
  if(legend->nz < 0) goto end2;
  if(fscanf(fp,"%d",&legend->pen1) != 1)
    goto end2;
  legend->points=new point_t[legend->np];
  for (n=0; n < legend->np;n++) {
    if(fscanf(fp,"%d %f %f",&legend->points[n].N,&legend->points[n].p,
                                                 &legend->points[n].t) != 3)
      goto end2;
    legend->points[n].z=new float[legend->nz];
    for(i=0;i<legend->nz;i++) if(fscanf(fp,"%f",&legend->points[n].z[i]) != 1)
      goto end2;
    }

  legend->pen0=0;
  legend->framed=0;
  legend->currentz=DefaultColumn;
  legend->text=0;

  *rstatus=LGD_STATUS_OK;
  return;

end2 :
  fclose(fp);
  *rstatus=LGD_STATUS_READ_ERROR;
  return;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void lgd_read04(FILE *fp, legend02_t *legend, int *rstatus)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,k,n,no_z=0;
  int ok;
  float x,y;

  if(fscanf(fp,"%d %d %d",&legend->ID,&legend->np,&legend->nz) != 3)
    goto end2;
  if(legend->np < 0) goto end2;
  if(legend->nz == 0) no_z=1;
//  if(legend->nz < 0) goto end2;
  if(fscanf(fp,"%d",&legend->pen1) != 1) goto end2;
  legend->points=new point_t[legend->np];
  for (n=0; n < legend->np;n++) {
    if(fscanf(fp,"%f %f",&legend->points[n].t, &legend->points[n].p) != 2)
      goto end2;
    legend->points[n].N=n+1;
    if(no_z==1) {
      legend->nz=1;
      legend->points[n].z=new float[legend->nz];
      legend->points[n].z[0]=n;
      }
    else {
      legend->points[n].z=new float[legend->nz];
      for(i=0;i<legend->nz;i++) if(fscanf(fp,"%f",&legend->points[n].z[i]) != 1)
        goto end2;
      }
    }

  legend->pen0=0;
  legend->framed=0;
  legend->currentz=DefaultColumn;
  legend->text=0;

  *rstatus=LGD_STATUS_OK;
  return;

end2 :
  fclose(fp);
  *rstatus=LGD_STATUS_READ_ERROR;
  return;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void lgd_read05(FILE *fp, legend03_t *legend, int *rstatus)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int   i,k,n;
  int   ok;
  float x,y;
  int renumber=0;
  legend01_t symbol; //=legend01_default;

  if(fscanf(fp,"%d %d %d",&legend->ID,&legend->ndepths,&legend->ncolumns) != 3)
    goto end2;
  if(legend->ndepths  <= 0) goto end2;
  if(legend->ncolumns <= 0) goto end2;

  legend->np=1;
  legend->points=new point3D_t[legend->np];
  if(fscanf(fp,"%f %f",&legend->points[0].t,&legend->points[0].p) != 2) goto end2;
  
  legend->points[0].depths=new float[legend->ndepths];
  legend->points[0].z=new float *[legend->ndepths];
  
  for (n=0; n < legend->ndepths;n++) {
    if(fscanf(fp,"%f",&legend->points[0].depths[n]) != 1) goto end2;
    legend->points[0].z[n]=new float[legend->ncolumns];
    for(i=0;i<legend->ncolumns;i++) if(fscanf(fp,"%f",&legend->points[0].z[n][i]) != 1)
      goto end2;
    }

  legend->pen0=0;
  legend->framed=0;
  legend->currentz=DefaultColumn;
  legend->text=new legend01_t;
  
  symbol.T=legend->points[0].t;
  symbol.P=legend->points[0].p;
  
  *(legend->text)=symbol;
  
/**----------------------------------------------------------------------------
  check for weird-formed node number*/
  for (n=1; n < legend->np;n++) {
    if(legend->points[n].N==legend->points[n-1].N) {
      renumber=1;
      break;
      }
    }
/**----------------------------------------------------------------------------
  if necessary, apply re-numbering*/
  if(renumber==1) {
    printf("warning: node numbering is deficient, xscan will apply patch procedure\n");
    for (n=0; n < legend->np;n++) {
      legend->points[n].N+=n;
      }
    }

  *rstatus=LGD_STATUS_OK;
  return;

  end2 :
    fclose(fp);
    *rstatus=LGD_STATUS_READ_ERROR;
    return;

}


/*-----------------------------------------------------------------------
 Load legend
-------------------------------------------------------------------------
     The symbol definition is in kernel.def file.
     3 types of legend are already defined:
     * Text & Symbold
     * Symbol only
     * Graphs
 -----------------------------------------------------------------------*/
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void lgd_load(const char *in,  int *rstatus)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/*----------------------------------------------------------------------------*/
int	n,status;
char	NAME[255],keyword;
int 	ok;
float 	x,y;
FILE   *fp=NULL;
legend01_t *ptr;
/*----------------------------------------------------------------------------*/

if((fp = fopen(in,"r")) == NULL)
  goto end3;

if(fscanf(fp,"%d",&lgd_n) != 1)
  goto end2;

lgd_alloc(MAX(1000,lgd_n),rstatus);

lgd_n=0;
for(n=0;n<lgd_nmax;n++) {
  if(fscanf(fp,"%s",NAME) != 1)
    goto end1;
  if(strcmp(NAME,"#TEXT&SYMBOL") == 0) {
    legendstab[n].Type=LGD_TEXT_SYMBOL;
    exitIfNull(
      legendstab[n].ptr=(char *) malloc(sizeof(legend01_t))
      );
    ptr=(legend01_t *) legendstab[n].ptr;
    ptr->Type=LGD_TEXT_SYMBOL;
    status=lgd_read01(fp,(legend01_t*) legendstab[n].ptr);
    if(status ==0) lgd_n++;
    }
  else if(strcmp(NAME,"#SYMBOL") == 0) {
    legendstab[n].Type=LGD_SYMBOL;
    exitIfNull(
      legendstab[n].ptr=(char *) malloc(sizeof(legend01_t))
      );
    ptr=(legend01_t *) legendstab[n].ptr;
    ptr->Type=LGD_SYMBOL;
    status=lgd_read01(fp,(legend01_t*) legendstab[n].ptr);
    if(status ==0) lgd_n++;
    }
  else if(strcmp(NAME,"#GRAPH-XYZ") == 0) {
    legendstab[n].Type=LGD_GRAPH;
    exitIfNull(
      legendstab[n].ptr=(char *) malloc(sizeof(legend02_t))
      );
    lgd_read02(fp,(legend02_t*) legendstab[n].ptr,&status);
    if(status ==0) lgd_n++;
    }
  else if(strcmp(NAME,"#GRAPH-YXZ") == 0) {
    legendstab[n].Type=LGD_GRAPH;
    exitIfNull(
      legendstab[n].ptr=(char *) malloc(sizeof(legend02_t))
      );
    lgd_read03(fp,(legend02_t*) legendstab[n].ptr,&status);
    if(status ==0) lgd_n++;
    }
  else if(strcmp(NAME,"#TRACK-XY") == 0) {
    legendstab[n].Type=LGD_GRAPH;
    legendstab[n].ptr=new char[sizeof(legend02_t)];
    lgd_read04(fp,(legend02_t*) legendstab[n].ptr,&status);
    if(status ==0) lgd_n++;
    }
  else if(strcmp(NAME,"#PROFILE") == 0) {
        legendstab[n].Type=LGD_PROFILE;
        legendstab[n].ptr=new char[sizeof(legend03_t)];
        lgd_read05(fp,(legend03_t*) legendstab[n].ptr,&status);
        if(status ==0) lgd_n++;
        }
  else {
    printf("#Should be a LGD type, is not... %s\n",NAME);
    goto end2;
    }
  if(status !=0) goto end3;

  end05:  continue;
}
end1 :   fclose(fp);
         rstatus=LGD_STATUS_OK;
         lgd_tocartesian(legendstab,&lgd_n,rstatus);
         return;

end2 :   fclose(fp);
         *rstatus=LGD_STATUS_READ_ERROR;
         lgd_tocartesian(legendstab,&lgd_n,rstatus);
         return;

end3 :   *rstatus=LGD_STATUS_UNABLE_TO_OPEN;
         return;

}


#undef LEGEND_MAIN
