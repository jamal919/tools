
/*******************************************************************************
T-UGO tools, 2006-2009

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

VERSION :

*******************************************************************************/

#include "../config.h"

#include "swap.h"

#include "gshhs.h"
#include "polygones.h"

#if HAVE_LIBSHP
#if HAVE_SHAPEFIL_H
#include <shapefil.h>
#elif HAVE_LIBSHP_SHAPEFIL_H
#include <libshp/shapefil.h>
#endif
#endif

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int plg_inquire_shp(const char *filename, int *np, int *nblocks, int *shapetype)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// Get number of shapes and shape types of an SHP file
/**
\param *filename
\param *np pointer to the total number of points.
\param *nblocks pointer to the number of polygons.
\param *shapetype pointer to the type of shapes.
\return \c errno if compiled, \c ENOEXEC otherwise
See http://shapelib.maptools.org/shp_api.html for more info
*/
/*----------------------------------------------------------------------------*/
{
#if HAVE_SHAPEFIL_H || HAVE_LIBSHP
  SHPHandle shpfile;
  int i;
  SHPObject *shp;
  
  shpfile=SHPOpen(filename,"rb");
  if(!shpfile)TRAP_ERR_RETURN(errno,1,"SHPOpen(\"%s\",,) error (%d %s)\n",filename,errno,strerror(errno));
  
  SHPGetInfo(shpfile,nblocks,shapetype,NULL,NULL);
  for(i=0;i<*nblocks;i++){
    shp=SHPReadObject(shpfile,i);
    switch(shp->nSHPType){
    case SHPT_POINT:
      *np+=1;
      break;
    case SHPT_POLYGON:
    case SHPT_POLYGONZ:
    case SHPT_ARC:
    case SHPT_ARCZ:
      *np+=shp->nParts;
      break;
    default:
      (*nblocks)--;
      }
    SHPDestroyObject(shp);
    }
  
  SHPClose(shpfile);
  
  return 0;
#else
  return ENOEXEC;
#endif
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_read_shp (const char *filename, plg_t *polygones, int np, int mode)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// Read shapes from an SHP file
/**
\param *filename
\param *polygones pointer to the polygones
\param np number of shapes
\returns STUFF
*/
// http://shapelib.maptools.org/shp_api.html
/*----------------------------------------------------------------------------*/
{
#if HAVE_SHAPEFIL_H || HAVE_LIBSHP
  SHPHandle shpfile;
  int k,l,oi,p;//object and polygone indexes
  SHPObject *shp;
  plg_t *pp;
  
  shpfile=SHPOpen(filename,"rb");
  
  p=0;
  for(oi=0;oi<np;oi++){
    shp=SHPReadObject(shpfile,oi);
    
    switch(shp->nSHPType){
    
    case SHPT_POINT:
      pp=&polygones[p];
      
      pp->init(1,(plg_init_mode) mode);
      
//       STDERR_BASE_LINE_FUNC("[%d]:([%g;%g],[%g;%g])\n",oi,shp->dfXMin,shp->dfXMax,shp->dfYMin,shp->dfYMax);
      pp->t[0]=shp->dfXMin;
      pp->p[0]=shp->dfYMin;
      
      p++;
      break;
    
    case SHPT_POLYGON:
    case SHPT_POLYGONZ:
    case SHPT_ARC:
    case SHPT_ARCZ:
      for(l=0;l<shp->nParts;l++) {
        pp=&polygones[p];
        if(l==shp->nParts-1)
          pp->npt=shp->nVertices-shp->panPartStart[l];
        else
          pp->npt=shp->panPartStart[l+1]-shp->panPartStart[l];
        pp->init((plg_init_mode) mode);
        for (k = 0; k < pp->npt; k++) {
          pp->t[k]=shp->padfX[shp->panPartStart[l]+k];
          pp->p[k]=shp->padfY[shp->panPartStart[l]+k];
          pp->x[k]=pp->t[k];
          pp->y[k]=pp->p[k];
          }
        p++;
        }
      break;
      }
    
    SHPDestroyObject(shp);
    }
  
  SHPClose(shpfile);
  
  return 0;
#else
  return -1;
#endif
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 int plg_write_shp (const char *filename, const plg_t *polygones, int np)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// Read shapes from an SHP file
/**
\param *filename
\param *polygones pointer to the polygones
\param np number of shapes
\returns STUFF
*/
// http://shapelib.maptools.org/shp_api.html
/*----------------------------------------------------------------------------*/
{
#if HAVE_SHAPEFIL_H || HAVE_LIBSHP
  SHPHandle shpfile;
  int status, p;//polygone index
  SHPObject *shp;
  
  
  shpfile=SHPCreate(filename,SHPT_POLYGON);
  for(p=0;p<np;p++){
    int nSHPType=SHPT_POLYGON;
    int nVertices=polygones[p].npt;
    double *padfX=polygones[p].t;
    double *padfY=polygones[p].p;
    double *padfZ=0;
    shp=SHPCreateSimpleObject( nSHPType,  nVertices, padfX, padfY, padfZ);
    status=SHPWriteObject( shpfile, -1, shp );
    SHPDestroyObject(shp);
    }

  SHPClose(shpfile);
  return(0);

#else
  return(-1);
#endif
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int plg_read_gshhs_1_6 (const char *filename, plg_t *polygones)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double w, e, s, n, area, lon, lat;
  int np,ns;
  char source;
  FILE  *fp=NULL;
  int  k, max_east = 270000000, n_read, flip;
  struct  POINT p;
  struct GSHHS_1_6 h;

  fp = fopen (filename, "rb");
  if(fp==0) TRAP_ERR_RETURN(errno,1,"Could not open %s (%s)\n",filename,strerror(errno));

  n_read = fread (&h, sizeof (struct GSHHS_1_6), 1, fp);
  flip = (! (h.level > 0 && h.level < 5));  /* Take as sign that byte-swabbing is needed */

  ns=0;
  np=0;
  while (n_read == 1) {
    if (flip) {
      h.id = swabi4 ((unsigned int)h.id);
      h.n  = swabi4 ((unsigned int)h.n);
      h.level = swabi4 ((unsigned int)h.level);
      h.west  = swabi4 ((unsigned int)h.west);
      h.east  = swabi4 ((unsigned int)h.east);
      h.south = swabi4 ((unsigned int)h.south);
      h.north = swabi4 ((unsigned int)h.north);
      h.area  = swabi4 ((unsigned int)h.area);
      h.version  = swabi4 ((unsigned int)h.version);
      h.greenwich = swabi2 ((unsigned int)h.greenwich);
      h.source = swabi2 ((unsigned int)h.source);
      }
    w = h.west  * 1.0e-6;  /* Convert from microdegrees to degrees */
    e = h.east  * 1.0e-6;
    s = h.south * 1.0e-6;
    n = h.north * 1.0e-6;
    source = (h.source == 1) ? 'W' : 'C';  /* Either WVS or CIA (WDBII) pedigree */
    area = 0.1 * h.area;      /* Now im km^2 */

//    printf ("P %6d%8d%2d%2c%13.3f%10.5f%10.5f%10.5f%10.5f\n", h.id, h.n, h.level, source, area, w, e, s, n);
    ns++;
    np=h.n;
    polygones[ns].npt=np;
    polygones[ns].t= new double[polygones[ns].npt];
    polygones[ns].p= new double[polygones[ns].npt];
    polygones[ns].x= new double[polygones[ns].npt];
    polygones[ns].y= new double[polygones[ns].npt];

    for (k = 0; k < h.n; k++) {
      if (fread (&p, sizeof(struct POINT), 1, fp) != 1) {
        fprintf (stderr, "gshhs:  Error reading file %s for polygon %d, point %d.\n", filename, h.id, k);
        return(-1);
        }
      if (flip) {
        p.x = swabi4 ((unsigned int)p.x);
        p.y = swabi4 ((unsigned int)p.y);
        }
      lon = (h.greenwich && p.x > max_east) ? p.x * 1.0e-6 - 360.0 : p.x * 1.0e-6;
      lat = p.y * 1.0e-6;
      polygones[ns].t[k]=lon;
      polygones[ns].p[k]=lat;
      polygones[ns].x[k]=lon;
      polygones[ns].y[k]=lat;
      }
    max_east = 180000000;  /* Only Eurasiafrica needs 270 */
    n_read = fread(&h, sizeof (struct GSHHS_1_6), 1, fp);
    }

  fclose (fp);

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int plg_inquire_gshhs_1_6 (const char *filename, int *ns, int *np)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double w, e, s, n, area;
  char source;
  FILE  *fp=NULL;
  int  max_east = 270000000, n_read, flip;
  struct GSHHS_1_6 h;
  int status=0;

  fp = fopen (filename, "rb");
  if(fp==0) TRAP_ERR_RETURN(errno,1,"Could not open %s (%s)\n",filename,strerror(errno));

  n_read = fread (&h, sizeof (struct GSHHS_1_6), 1, fp);
  flip = (! (h.level > 0 && h.level < 5));  /* Take as sign that byte-swabbing is needed */

  *ns=0;
  *np=0;
  while (n_read == 1) {
    if (flip) {
      h.id = swabi4 ((unsigned int)h.id);
      h.n  = swabi4 ((unsigned int)h.n);
      h.level = swabi4 ((unsigned int)h.level);
      h.west  = swabi4 ((unsigned int)h.west);
      h.east  = swabi4 ((unsigned int)h.east);
      h.south = swabi4 ((unsigned int)h.south);
      h.north = swabi4 ((unsigned int)h.north);
      h.area  = swabi4 ((unsigned int)h.area);
      h.version  = swabi4 ((unsigned int)h.version);
      h.greenwich = swabi2 ((unsigned int)h.greenwich);
      h.source = swabi2 ((unsigned int)h.source);
      }
    if((h.source!=0) && (h.source!=1)) {
      status=-1;
      goto finished;
      }
    (*ns)++;
    *np+=h.n;
    w = h.west  * 1.0e-6;  /* Convert from microdegrees to degrees */
    e = h.east  * 1.0e-6;
    s = h.south * 1.0e-6;
    n = h.north * 1.0e-6;
    source = (h.source == 1) ? 'W' : 'C';  /* Either WVS or CIA (WDBII) pedigree */
    area = 0.1 * h.area;      /* Now im km^2 */

//    printf ("P %6d%8d%2d%2c%13.3f%10.5f%10.5f%10.5f%10.5f\n", h.id, h.n, h.level, source, area, w, e, s, n);
    fseek (fp, (long)(h.n * sizeof(struct POINT)), SEEK_CUR);
    max_east = 180000000;  /* Only Eurasiafrica needs 270 */
    n_read = fread(&h, sizeof (struct GSHHS_1_6), 1, fp);
    }

finished:
  fclose (fp);

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int plg_read_gshhs_2_0 (const char *filename, plg_t *polygones)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double w, e, s, n, area, lon, lat;
  int np,ns;
  char /*source,*/greenwich;
  FILE  *fp;
  int  k, max_east = 270000000, n_read, flip;
  struct  POINT p;
  struct GSHHS_2_0 h;

  fp = fopen (filename, "rb");
  if(fp==0) TRAP_ERR_RETURN(errno,1,"Could not open %s (%s)\n",filename,strerror(errno));

  n_read = fread ((void *)&h, (size_t)sizeof (struct GSHHS_2_0), (size_t)1, fp);
  flip = (! (h.flag > 0 && h.flag < 5));  /* Take as sign that byte-swabbing is needed */

  ns=0;
  np=0;
  while (n_read == 1) {
    if (flip) {
      h.id = swabi4 ((unsigned int)h.id);
      h.n  = swabi4 ((unsigned int)h.n);
      h.flag = swabi4 ((unsigned int)h.flag);
      h.west  = swabi4 ((unsigned int)h.west);
      h.east  = swabi4 ((unsigned int)h.east);
      h.south = swabi4 ((unsigned int)h.south);
      h.north = swabi4 ((unsigned int)h.north);
      h.area  = swabi4 ((unsigned int)h.area);
//      h.version  = swabi4 ((unsigned int)h.version);
//      h.greenwich = swabi2 ((unsigned int)h.greenwich);
//      h.source = swabi2 ((unsigned int)h.source);
      }
    w = h.west  * 1.0e-6;  /* Convert from microdegrees to degrees */
    e = h.east  * 1.0e-6;
    s = h.south * 1.0e-6;
    n = h.north * 1.0e-6;
//    source = (h.source == 1) ? 'W' : 'C';  /* Either WVS or CIA (WDBII) pedigree */
    area = 0.1 * h.area;      /* Now im km^2 */

//    printf ("P %6d%8d%2d%2c%13.3f%10.5f%10.5f%10.5f%10.5f\n", h.id, h.n, h.level, source, area, w, e, s, n);
    np+=h.n;
    polygones[ns].npt=h.n;
    polygones[ns].t= new double[polygones[ns].npt];
    polygones[ns].p= new double[polygones[ns].npt];
    polygones[ns].x= new double[polygones[ns].npt];
    polygones[ns].y= new double[polygones[ns].npt];
    greenwich = (h.flag >> 16) & 1;
    for (k = 0; k < h.n; k++) {
      if (fread ((void *)&p, (size_t)sizeof(struct POINT), (size_t)1, fp) != 1) {
        fprintf (stderr, "gshhs:  Error reading file %s for polygon %d, point %d.\n", filename, h.id, k);
        return(-1);
        }
      if (flip) {
        p.x = swabi4 ((unsigned int)p.x);
        p.y = swabi4 ((unsigned int)p.y);
        }
      lon = (greenwich && p.x > max_east) ? p.x * 1.0e-6 - 360.0 : p.x * 1.0e-6;
      lat = p.y * 1.0e-6;
      polygones[ns].t[k]=lon;
      polygones[ns].p[k]=lat;
      polygones[ns].x[k]=lon;
      polygones[ns].y[k]=lat;
//      printf ("%10.5f%9.5f\n", lon, lat);
      }
    ns++;
    max_east = 180000000;  /* Only Eurasiafrica needs 270 */
    n_read = fread((void *)&h, (size_t)sizeof (struct GSHHS_2_0), (size_t)1, fp);
    }

  fclose (fp);

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int plg_inquire_gshhs_2_0 (const char *filename, int *ns, int* np)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double w, e, s, n, area;
  FILE  *fp;
  int  max_east = 270000000, n_read, flip;
  struct GSHHS_2_0 h;

  fp = fopen (filename, "rb");
  if(fp==0) TRAP_ERR_RETURN(errno,1,"Could not open %s (%s)\n",filename,strerror(errno));

  n_read = fread ((void *)&h, (size_t)sizeof (struct GSHHS_2_0), (size_t)1, fp);
  flip = (! (h.flag > 0 && h.flag < 5));  /* Take as sign that byte-swabbing is needed */

  *ns=0;
  *np=0;
  while (n_read == 1) {
    if (flip) {
      h.id = swabi4 ((unsigned int)h.id);
      h.n  = swabi4 ((unsigned int)h.n);
      h.flag = swabi4 ((unsigned int)h.flag);
      h.west  = swabi4 ((unsigned int)h.west);
      h.east  = swabi4 ((unsigned int)h.east);
      h.south = swabi4 ((unsigned int)h.south);
      h.north = swabi4 ((unsigned int)h.north);
      h.area  = swabi4 ((unsigned int)h.area);
//      h.version  = swabi4 ((unsigned int)h.version);
//      h.greenwich = swabi2 ((unsigned int)h.greenwich);
//      h.source = swabi2 ((unsigned int)h.source);
      }
    *np+=h.n;
    w = h.west  * 1.0e-6;  /* Convert from microdegrees to degrees */
    e = h.east  * 1.0e-6;
    s = h.south * 1.0e-6;
    n = h.north * 1.0e-6;
//    source = (h.source == 1) ? 'W' : 'C';  /* Either WVS or CIA (WDBII) pedigree */
    area = 0.1 * h.area;      /* Now im km^2 */

//    printf ("P %6d%8d%2d%2c%13.3f%10.5f%10.5f%10.5f%10.5f\n", h.id, h.n, h.level, source, area, w, e, s, n);
    fseek (fp, (long)(h.n * sizeof(struct POINT)), SEEK_CUR);
    (*ns)++;
    max_east = 180000000;  /* Only Eurasiafrica needs 270 */
    n_read = fread((void *)&h, (size_t)sizeof (struct GSHHS_2_0), (size_t)1, fp);
    }

  fclose (fp);

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int plg_inquire_gshhs (const char *filename, int *ns, int* np, int *fmt)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;

  status=plg_inquire_gshhs_1_6 (filename, ns, np);
  
  if(status!=0) {
    status=plg_inquire_gshhs_2_0 (filename, ns, np);
    if(status==0) *fmt=PLG_FORMAT_GSHHS_2_0;
    }
  else {
    *fmt=PLG_FORMAT_GSHHS_1_6;
    }
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int plg_read_gshhs (const char *filename,  plg_t *polygones, int fmt)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=-1;

  switch(fmt) {
    case PLG_FORMAT_GSHHS_1_6:
      status=plg_read_gshhs_1_6 (filename, polygones);
      break;
    case PLG_FORMAT_GSHHS_2_0:
      status=plg_read_gshhs_2_0 (filename, polygones);
      break;
    }
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int plg_read_mif (const char *filename,float *rlon, float *rlat, int *first, int *last, int *ns)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  FILE  *fp=NULL;
  int  k, n, nitems, np;
  char keyword[256];
  double lon,lat,lon0,lat0;

  *ns=0;
  np=0;
  fp = fopen (filename, "r");
  if(fp==0) TRAP_ERR_RETURN(errno,1,"Could not open %s (%s)\n",filename,strerror(errno));

  lon0=-1.e+10;
  lat0=-1.e+10;


  while (strcmp(keyword,"Data")!=0) {
    nitems=fscanf(fp,"%s",keyword);
    }
  while(!feof(fp)) {
    while (strcmp(keyword,"Pline")!=0) {
      nitems=fscanf(fp,"%s",keyword);
      if(nitems<=0) break;
      }
    if(nitems<=0) break;
    nitems=fscanf(fp,"%d",&n);
//    printf("segment %d, npoints %d\n",*ns,n);
    nitems=fscanf(fp,"%lf %lf",&lon,&lat);
    if((lon==lon0) && (lat==lat0)) {
      np+=n;
      last[*ns]=np;
      }
    else {
      (*ns)++;
      np+=n;
      first[*ns]=np-n+1;
      last[*ns]=np;
      }
    rlon[last[*ns]-n]=lon;
    rlat[last[*ns]-n]=lat;
    for(k=1;k<n;k++) {
      nitems=fscanf(fp,"%lf %lf",&lon,&lat);
      rlon[last[*ns]-n+k]=lon;
      rlat[last[*ns]-n+k]=lat;
      }
    lon0=lon;
    lat0=lat;
    strcpy(keyword,"");
    }

  fclose (fp);
  printf("#segments %d, #points %d\n",*ns,np);

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int plg_inquire_mif (const char *filename, int *ns, int* np)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  FILE  *fp=NULL;
  int  n, nitems;
  char keyword[256];

  *ns=0;
  *np=0;
  fp = fopen (filename, "r");
  if(fp==0) TRAP_ERR_RETURN(errno,1,"Could not open %s (%s)\n",filename,strerror(errno));

  while (strcmp(keyword,"Data")!=0) {
    nitems=fscanf(fp,"%s",keyword);
    }
  while(!feof(fp)) {
    while (strcmp(keyword,"Pline")!=0) {
      nitems=fscanf(fp,"%s",keyword);
      if(nitems<=0) break;
      }
    if(nitems<=0) break;
    nitems=fscanf(fp,"%d",&n);
//    printf("segment %d, npoints %d\n",*ns,n);
    (*ns)++;
    *np+=n;
    strcpy(keyword,"");
    }

  fclose (fp);
  printf("#segments %d, #points %d\n",*ns,*np);

  return(0);
}

typedef struct {
  double  xmin,xmax,ymin,ymax;
  }vframe_t;

typedef struct {
  vframe_t  frame;
  short    type;
  size_t   pos;
  size_t   ndata;
} vct00_header_t;

typedef struct {
  short   pen;
  double  lon,lat;
} vct00_data_t;

typedef struct {
  vct00_header_t *headers;
  vct00_data_t   **data;
  size_t   nblocks;
} vct00_t;

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int plg_read_vct00 (const char *filename,float *rlon, float *rlat, int *first, int *last, int *ns)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int  np;
  FILE *in=NULL;
  int k,m,n;
  int  type,val[32];
  short sval[10];
  vct00_t buffer;

  *ns=0;
  np=0;

  in=fopen(filename,"rb");
  if (in == NULL) {
    return(-1);
    }

/*------------------------------------------------------------------------------
  scan numbre of blocks */
  buffer.nblocks=0;
  type=0;
  while (type != -1) {
    for (k=0;k<2;k++) {
      fread(&(val[k]),sizeof(int),1,in);
      }
    fread(&(sval[0]),sizeof(short),1,in);

    for (k=2;k<4;k++) {
      fread(&(val[k]),sizeof(int),1,in);
      }
    fread(&(sval[1]),sizeof(short),1,in);

    for (k=4;k<6;k++) {
      fread(&(val[k]),sizeof(int),1,in);
      }
    fread(&(sval[2]),sizeof(short),1,in);

    for (k=6;k<8;k++) {
      fread(&(val[k]),sizeof(int),1,in);
      }
    fread(&(sval[3]),sizeof(short),1,in);

    type=sval[0];
    buffer.nblocks++;
    }

  *ns=buffer.nblocks-1;

  if ((buffer.headers=(vct00_header_t *) malloc(buffer.nblocks*sizeof(vct00_header_t)))==NULL) {
    TRAP_ERR_EXIT(-1,"");perror("buffer.headers");
  }

  rewind(in);

  for (n=0;n<buffer.nblocks;n++) {
    for (k=0;k<2;k++) {
      fread(&(val[k]),sizeof(int),1,in);
      }
    fread(&(sval[0]),sizeof(short),1,in);

    buffer.headers[n].pos  =val[0];
    buffer.headers[n].ndata=val[1];

    for (k=2;k<4;k++) {
      fread(&(val[k]),sizeof(int),1,in);
      }
    fread(&(sval[1]),sizeof(short),1,in);

    for (k=4;k<6;k++) {
      fread(&(val[k]),sizeof(int),1,in);
      }
    fread(&(sval[2]),sizeof(short),1,in);

    for (k=6;k<8;k++) {
      fread(&(val[k]),sizeof(int),1,in);
      }
    fread(&(sval[3]),sizeof(short),1,in);

    buffer.headers[n].frame.ymax=(double) val[4]/1000000.0;
    buffer.headers[n].frame.xmin=(double) val[5]/1000000.0;

    buffer.headers[n].frame.ymin=(double) val[6]/1000000.0;
    buffer.headers[n].frame.xmax=(double) val[7]/1000000.0;

    buffer.headers[n].type=sval[0];
    }

  *ns=-1;
  np=0;
  for (n=0;n<buffer.nblocks-1;n++) {
    for (m=0;m<buffer.headers[n].ndata;m++) {
      for (k=8;k<10;k++) {
        fread(&(val[k]),sizeof(int),1,in);
        }
      fread(&(sval[3]),sizeof(short),1,in);
      if(sval[3]==0) {
        *ns=*ns+1;
        first[*ns]=np+1;
        if(*ns>1) last[*ns-1]=np;
        }
      rlon[np]=(float) val[9]/1000000.0;
      rlat[np]=(float) val[8]/1000000.0;
      np++;
      }
    }

  fclose(in);
  printf("#segments %d, #points %d\n",*ns,np);

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int plg_inquire_vct00 (const char *filename, int *ns, int* np)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  FILE *in=NULL;
  int k,m,n;
  int  type,val[32];
  short sval[10];
  vct00_t buffer;

  *ns=0;
  *np=0;

  in=fopen(filename,"rb");
  if (in == NULL) {
    return(-1);
    }

/*------------------------------------------------------------------------------
  scan number of blocks */
  buffer.nblocks=0;
  type=0;
  while (type != -1) {
    for (k=0;k<2;k++) {
      fread(&(val[k]),sizeof(int),1,in);
      }
    fread(&(sval[0]),sizeof(short),1,in);

    for (k=2;k<4;k++) {
      fread(&(val[k]),sizeof(int),1,in);
      }
    fread(&(sval[1]),sizeof(short),1,in);

    for (k=4;k<6;k++) {
      fread(&(val[k]),sizeof(int),1,in);
      }
    fread(&(sval[2]),sizeof(short),1,in);

    for (k=6;k<8;k++) {
      fread(&(val[k]),sizeof(int),1,in);
      }
    fread(&(sval[3]),sizeof(short),1,in);

    type=sval[0];
    buffer.nblocks++;
    }

  if ((buffer.headers=(vct00_header_t *) malloc(buffer.nblocks*sizeof(vct00_header_t)))==NULL) {
    TRAP_ERR_EXIT(-1,"");perror("buffer.headers");
  }

  rewind(in);

  for (n=0;n<buffer.nblocks;n++) {
    for (k=0;k<2;k++) {
      fread(&(val[k]),sizeof(int),1,in);
      }
    fread(&(sval[0]),sizeof(short),1,in);

    buffer.headers[n].pos  =val[0];
    buffer.headers[n].ndata=val[1];
    *np+=buffer.headers[n].ndata;

    for (k=2;k<4;k++) {
      fread(&(val[k]),sizeof(int),1,in);
      }
    fread(&(sval[1]),sizeof(short),1,in);

    for (k=4;k<6;k++) {
      fread(&(val[k]),sizeof(int),1,in);
      }
    fread(&(sval[2]),sizeof(short),1,in);

    for (k=6;k<8;k++) {
      fread(&(val[k]),sizeof(int),1,in);
      }
    fread(&(sval[3]),sizeof(short),1,in);

    buffer.headers[n].frame.ymax=(double) val[4]/1000000.0;
    buffer.headers[n].frame.xmin=(double) val[5]/1000000.0;

    buffer.headers[n].frame.ymin=(double) val[6]/1000000.0;
    buffer.headers[n].frame.xmax=(double) val[7]/1000000.0;

    buffer.headers[n].type=sval[0];
    }

  *ns=-1;
  for (n=0;n<buffer.nblocks-1;n++) {
    for (m=0;m<buffer.headers[n].ndata;m++) {
      for (k=8;k<10;k++) {
        fread(&(val[k]),sizeof(int),1,in);
        }
      fread(&(sval[3]),sizeof(short),1,in);
      if(sval[3]==0) *ns=*ns+1;
      }
    }

  fclose(in);
  printf("#segments %d, #points %d\n",*ns,*np);

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int plg_inquire_vct00_fortran (const char *filename, int *ns, int* np)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
int status;

status=plg_inquire_vct00_fortran (filename, ns, np);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int plg_read_cst (const char *filename, plg_t *polygones)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  FILE  *fp=NULL;
  int  k, np, ns, n_read, flip;

  unsigned int values[8],blocksize;
  float *buffer;

  fp = fopen (filename, "rb");
  if(fp==0) TRAP_ERR_RETURN(errno,1,"Could not open %s (%s)\n",filename,strerror(errno));

  n_read = fread (values, (size_t) 5*sizeof (int), 1, fp);

  if(values[0]!=values[4]) goto error;

  flip = (values[2]!=1);  /* Take as sign that byte-swabbing is needed */
  ns=0;
  np=0;
  while (n_read == 1) {
    if (flip) {
      for(k=0;k<5;k++) {
        values[k] = swabi4 ((unsigned int) values[k]);
        }
      }
    np=values[3]-values[2]+1;
    polygones[ns].init(np,PLG_INIT_SEPARATE);
    
    buffer=new float[2*np];
    n_read = fread (&blocksize, sizeof (int), 1, fp);
    if (flip)  blocksize= swabi4 ((unsigned int) blocksize);
    n_read = fread (buffer, (size_t) 2*np*sizeof (float), 1, fp);
    n_read = fread (&blocksize, sizeof (int), 1, fp);
    if (flip)  blocksize= swabi4 ((unsigned int) blocksize);
    if (flip) {
      for(k=0;k<2*np;k++) {
        buffer[k] = swap(buffer[k]);
        }
      }
    for(k=0;k<np;k++) {
      polygones[ns].t[k]=buffer[2*k];
      polygones[ns].p[k]=buffer[2*k+1];
      polygones[ns].x[k]=buffer[2*k];
      polygones[ns].y[k]=buffer[2*k+1];
      }
    n_read = fread (values, (size_t) 5*sizeof (int), 1, fp);
    ns++;
    delete[] buffer;
    }

  fclose (fp);
  return(0);

error:
  fclose (fp);
  return(-1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int plg_write_cst (const char *filename, const plg_t *polygones, int npolygones)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  FILE  *fp=NULL;
  int  k, np, ns, n_items, flip;

  unsigned int values[8],blocksize;
  float *buffer=NULL;

  fp = fopen (filename, "wb");
  if(fp==0) TRAP_ERR_RETURN(errno,1,"Could not open %s (%s)\n",filename,strerror(errno));
#define SWAPIO
#ifdef SWAPIO
  flip=1;
#else
  flip=0;
#endif
  np=1;
  for(ns=0;ns<npolygones;ns++) {
    values[0]=12;
    values[1]=ns;
    values[2]=np;
    values[3]=np+polygones[ns].npt-1;
    values[4]=12;
    if (flip) {
      for(k=0;k<5;k++) {
        values[k] = swabi4 ((unsigned int) values[k]);
        }
      }
    n_items = fwrite (values, (size_t) 5*sizeof (int), 1, fp);
    buffer=new float[2*polygones[ns].npt];
    for(k=0;k<polygones[ns].npt;k++) {
      buffer[2*k]  =polygones[ns].t[k];
      buffer[2*k+1]=polygones[ns].p[k];
      }
    if (flip) {
      for(k=0;k<2*polygones[ns].npt;k++) {
        buffer[k] = swap(buffer[k]);
        }
      }
    blocksize=4*2*polygones[ns].npt;
    if (flip)  blocksize= swabi4 ((unsigned int) blocksize);
    n_items = fwrite (&blocksize, sizeof (int), 1, fp);
    n_items = fwrite (buffer, (size_t) 2*polygones[ns].npt*sizeof (float), 1, fp);
    n_items = fwrite (&blocksize, sizeof (int), 1, fp);
    delete[] buffer;
    }

  fclose (fp);

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int plg_inquire_cst (const char *filename, int *ns, int *np)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double n;
  FILE  *fp=NULL;
  int  k, n_read, flip;

  unsigned int values[8];

  fp = fopen (filename, "r");
  if(fp==0) TRAP_ERR_RETURN(errno,1,"Could not open %s (%s)\n",filename,strerror(errno));

  n_read = fread (values, (size_t) 5*sizeof (int), 1, fp);
  flip = (values[2]!=1);  /* Take as sign that byte-swabbing is needed */
  *ns=0;
  *np=0;
  while (n_read == 1) {
    if (flip) {
      for(k=0;k<5;k++) {
        values[k] = swabi4 ((unsigned int) values[k]);
        }
      }
    (*ns)++;
    *np+=values[3]-values[2]+1;
    n=values[3]-values[2]+1;
    fseek (fp, (long)(2 *n * sizeof(float)+2*sizeof (int)), SEEK_CUR);
    n_read = fread (values, (size_t) 5*sizeof (int), 1, fp);
    }

  fclose (fp);

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int plg_read_xiso (const char *filename, plg_t *polygones)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  FILE  *fp=NULL;
  int  np, ns;
  int  k,n, nitems, recordsize, flip;
  float z;

  unsigned int values[8],blocksize;
  double *buffer=NULL;

  fp = fopen (filename, "rb");
  if(fp==0) TRAP_ERR_RETURN(errno,1,"Could not open %s (%s)\n",filename,strerror(errno));

  nitems = fread (&recordsize, sizeof (int), 1, fp);
  flip = (recordsize!=2);  /* Take as sign that byte-swabbing is needed */
  nitems = fread (&n, sizeof (int), 1, fp);
  nitems = fread (&z, sizeof (float), 1, fp);
  nitems = fread (&recordsize, sizeof (int), 1, fp);

  nitems = fread (&recordsize, sizeof (int), 1, fp);
  nitems = fread (values, (size_t) 3*sizeof (int), 1, fp);
  nitems = fread (&recordsize, sizeof (int), 1, fp);

  ns=0;
  np=0;
  while (nitems == 1) {
    if (flip) {
      n = swabi4 ((unsigned int) n);
      }
    np=n;
    polygones[ns].npt=np;
    polygones[ns].t= new double[polygones[ns].npt];
    polygones[ns].p= new double[polygones[ns].npt];
    polygones[ns].x= new double[polygones[ns].npt];
    polygones[ns].y= new double[polygones[ns].npt];
    buffer=new double[2*np];
    nitems = fread (&blocksize, sizeof (int), 1, fp);
    if (flip)  blocksize= swabi4 ((unsigned int) blocksize);
    nitems = fread (buffer, (size_t) 2*np*sizeof (double), 1, fp);
    nitems = fread (&blocksize, sizeof (int), 1, fp);
    if (flip)  blocksize= swabi4 ((unsigned int) blocksize);
    if (flip) {
      for(k=0;k<2*np;k++) {
        buffer[k] = swap(buffer[k]);
        }
      }
    for(k=0;k<np;k++) {
      polygones[ns].t[k]=buffer[2*k];
      polygones[ns].p[k]=buffer[2*k+1];
      polygones[ns].x[k]=buffer[2*k];
      polygones[ns].y[k]=buffer[2*k+1];
      }
    nitems = fread (&recordsize, sizeof (int), 1, fp);
    nitems = fread (&n, sizeof (int), 1, fp);
    nitems = fread (&z, sizeof (float), 1, fp);
    nitems = fread (&recordsize, sizeof (int), 1, fp);

    nitems = fread (&recordsize, sizeof (int), 1, fp);
    nitems = fread (values, (size_t) 3*sizeof (int), 1, fp);
    nitems = fread (&recordsize, sizeof (int), 1, fp);
    ns++;
    delete[] buffer;
    }

  fclose (fp);

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int plg_inquire_xiso (const char *filename, int *ns, int *np)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  FILE  *fp=NULL;
  int  n, nitems, recordsize, flip;
  float z;

  unsigned int values[8];

  fp = fopen (filename, "rb");
  if(fp==0) TRAP_ERR_RETURN(errno,1,"Could not open %s (%s)\n",filename,strerror(errno));

  nitems = fread (&recordsize, sizeof (int), 1, fp);
  flip = (recordsize!=2);  /* Take as sign that byte-swabbing is needed */
  nitems = fread (&n, sizeof (int), 1, fp);
  nitems = fread (&z, sizeof (float), 1, fp);
  nitems = fread (&recordsize, sizeof (int), 1, fp);

  nitems = fread (&recordsize, sizeof (int), 1, fp);
  nitems = fread (values, (size_t) 3*sizeof (int), 1, fp);
  nitems = fread (&recordsize, sizeof (int), 1, fp);

  *ns=0;
  *np=0;
  while (nitems == 1) {
    if (flip) {
        n = swabi4 ((unsigned int) n);
       }
    (*ns)++;
    *np+=n;

    fseek (fp, (long)(2 *n * sizeof(double)+2*sizeof (int)), SEEK_CUR);

    nitems = fread (&recordsize, sizeof (int), 1, fp);
    nitems = fread (&n, sizeof (int), 1, fp);
    nitems = fread (&z, sizeof (float), 1, fp);
    nitems = fread (&recordsize, sizeof (int), 1, fp);

    nitems = fread (&recordsize, sizeof (int), 1, fp);
    nitems = fread (values, (size_t) 3*sizeof (int), 1, fp);
    nitems = fread (&recordsize, sizeof (int), 1, fp);
    }

  fclose (fp);

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int plg_write_boundaries (const char *filename, const plg_t *polygones, int npolygones)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  FILE  *fp=NULL;
  int  k,n;
  float depth=0.0;

  fp = fopen (filename, "w");
  if(fp==0) TRAP_ERR_RETURN(errno,1,"Could not open %s (%s)\n",filename,strerror(errno));

  for(n=0;n<npolygones;n++) {
    fprintf(fp,"%d %f\n",n+1,depth);
    for(k=0;k<polygones[n].npt-1;k++) {
      fprintf(fp,"%lf %lf\n",polygones[n].t[k],polygones[n].p[k]);
      }
    fprintf(fp,"%d %d\n",polygones[n].npt,-9999);
    }
  fclose (fp);

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int plg_read_boundaries (const char *filename, plg_t *polygones)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  FILE  *fp=NULL;
  int  k,m,n, dum, np, ns;
  char line[1024],*read;

  fp = fopen (filename, "r");
  if(fp==0) TRAP_ERR_RETURN(errno,1,"Could not open %s (%s)\n",filename,strerror(errno));

  ns=0;
  np=0;
  while(!feof(fp)) {
    fgets(line,1024,fp);
    if(feof(fp)) goto finished;
    sscanf(line,"%d %d",&n,&dum);
    m=0;
    while(dum!=-9999) {
      read=fgets(line,1024,fp);
      if(read==0) {m++;break;}
      sscanf(line,"%d %d",&n,&dum);
      m++;
      }
    polygones[ns].npt=m-1;
    polygones[ns].t= new double[polygones[ns].npt];
    polygones[ns].p= new double[polygones[ns].npt];
    polygones[ns].x= new double[polygones[ns].npt];
    polygones[ns].y= new double[polygones[ns].npt];
    ns++;
    }

finished:
  rewind(fp);

  for(n=0;n<ns;n++) {
    fgets(line,1024,fp);
    for(k=0;k<polygones[n].npt;k++) {
      fgets(line,1024,fp);
      line[strlen(line)-1]='\0';
      sscanf(line,"%lf %lf",&(polygones[n].t[k]),&(polygones[n].p[k]));
      }
    fgets(line,1024,fp);
    }

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int plg_inquire_boundaries (const char *filename, int *ns, int *np)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  FILE  *fp=NULL;
  int  m,n, dum;
  char line[1024],*read;

  fp = fopen (filename, "r");
  if(fp==0) TRAP_ERR_RETURN(errno,1,"Could not open %s (%s)\n",filename,strerror(errno));
  
  *ns=0;
  *np=0;
  while(!feof(fp)) {
    fgets(line,1024,fp);
    if(feof(fp)) goto finished;
    sscanf(line,"%d %d",&n,&dum);
    m=0;
    while(dum!=-9999) {
      read=fgets(line,1024,fp);
      if(read==0) break;
      sscanf(line,"%d %d",&n,&dum);
      m++;
      }
    (*np)+=m;
    (*ns)++;
    }

finished:
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_read_descriptor (const char *filename, vector<plg_t> & boundaries, bool check, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  FILE  *fp=NULL;
  int  m,n,p;
  int  id,numblock,nitems,status;
  char line[1024],key[1024],scanfile[1024];
  plg_t *polygones=NULL;
  int npolygones;
  plg_desc_t descriptor;

  fp = fopen (filename, "r");
  if(fp==0) TRAP_ERR_RETURN(errno,1,"Could not open %s (%s)\n",filename,strerror(errno));

/*------------------------------------------------------------------------------
  count the number of blocks*/
  descriptor.nblock=0;
  fgets(line,1024,fp);
  while(!feof(fp)) {
    if(strcmp(line,"End_of_Block")==0) descriptor.nblock++;
    nitems=fscanf(fp,"%s",&line);
//    if(strcmp(line,"End_of_Description")==0) goto finished;
    }
  rewind(fp);

/*------------------------------------------------------------------------------
  count the number of segments in blocks*/
  fgets(line,1024,fp);
  sscanf(line,"%s",scanfile);
  status=plg_load_scan(scanfile, &polygones, &npolygones);
  if(status!=0) TRAP_ERR_RETURN(status,1,"plg_load_scan(\"%s\",,) error (%d: %s)\n",scanfile,status,strerror(status));

  descriptor.b=new plg_block_t[descriptor.nblock];

  for(n=0;n<descriptor.nblock;n++) {
    fgets(line,1024,fp);
    m=0;
    nitems=fscanf(fp,"%s",&key);
    fgets(line,1024,fp);
    while(strcmp(key,"End_of_Block")!=0) {
      fgets(line,1024,fp);
      if(m==0) {
        fgets(line,1024,fp);
        fgets(line,1024,fp);
        }
      m++;
      nitems=fscanf(fp,"%s",&key);
      fgets(line,1024,fp);
      }
    descriptor.b[n].segments.n=m;
    }

  rewind(fp);

/*------------------------------------------------------------------------------
  get descriptor intormations*/
  descriptor.scanfile=new char[strlen(scanfile)+1];
  fgets(line,1024,fp);
  sscanf(line,"%s",descriptor.scanfile);

  for(n=0;n<descriptor.nblock;n++) {
    fgets(line,1024,fp);
//    if(feof(fp)) goto finished;
    nitems=sscanf(line,"%d",&numblock);
    descriptor.b[n].segments.p=new plg_t[descriptor.b[n].segments.n];
    descriptor.b[n].code=new char[descriptor.b[n].segments.n];
    for(m=0;m<descriptor.b[n].segments.n;m++) {
      fgets(line,1024,fp);
      sscanf(line,"%d",&(descriptor.b[n].segments.p[m].id));
      fgets(line,1024,fp);
      sscanf(line,"%c",&(descriptor.b[n].code[m]));
      if( (descriptor.b[n].code[m]!='M') and (descriptor.b[n].code[m]!='T') and (descriptor.b[n].code[m]!='X') ) {
        printf("boundary block %d, segement %d: code %c is not valid (should be M/X for maritime boundary or T for terrestrial boundary), abandon\n",n,m);
        return(-1);
        }
      if(m==0) {
        fgets(line,1024,fp);
        sscanf(line,"%d",&(descriptor.b[n].orientation));
        fgets(line,1024,fp);
        sscanf(line,"%c",&(descriptor.b[n].type));
        if( (descriptor.b[n].type!='E') && (descriptor.b[n].type!='I') ) {
          printf("boundary block %d: type %c is not valid (should be E for external boundary or I for internal boundary), abandon\n",n);
          return(-1);
          }
        }
      }
    fgets(line,1024,fp);
    }

/*------------------------------------------------------------------------------
  process descriptor intormations*/
  
  status=plg_load_scan(scanfile, &polygones, &npolygones);
  
  if(check) {
    status=plg_checkAutoSecant(polygones, npolygones, PLG_SPHERICAL);
    if(status!=0) {
      printf("polygon set (%s) is not properly built, abandon\n",scanfile);
      return(-1);
      }
    }

  for(n=0;n<descriptor.nblock;n++) {
    int target[2];
    plg_t *tmp=new plg_t;
/*------------------------------------------------------------------------------
    agregate block segments*/
    int first=plg_identify(descriptor.b[n].segments.p[0].id, polygones, npolygones);
    if(first==-1) {
      printf("boundary block %d is not consistent (segment %d), abandon\n",n,0);
      return(-1);
      }
    target[0]=first;
    
    polygones[first].SetFlag(descriptor.b[n].code[0]);
      
    for(m=1;m<descriptor.b[n].segments.n;m++) {
      id=descriptor.b[n].segments.p[m].id;
      p=plg_identify( id, polygones, npolygones);
      if(p==-1) {
        printf("boundary block %d is not consistent (segment %d does not exist), abandon\n",n,p);
        return(-1);
        }
      target[1]=p;
      polygones[p].SetFlag(descriptor.b[n].code[m]);
      status=plg_concat (target, polygones, npolygones);
      if(status==-1) {
        printf("boundary block %d is not consistent (segment %d not connected), abandon\n",n,p);
        return(-1);
        }
      }
      
    tmp->duplicate(polygones[first]);
    boundaries.push_back(*tmp);
    }

  if(check) {
    status=plg_checkAutoSecant(boundaries, PLG_SPHERICAL);
    if(status!=0) {
      printf("boundary set is not safe, abandon\n");
      return(-1);
      }
    }

  fclose(fp);
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_read_descriptor(const char *filename, plg_t* & boundaries, int & nboundaries, bool check, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  vector<plg_t> tmp;
  
  status=plg_read_descriptor(filename, tmp, check, debug);
  
  if(status!=0){
    plg_destroy_entries(tmp);
    TRAP_ERR_RETURN(status,1,"plg_read_descriptor(\"%s\",,,) error : see above for messages\n",filename);
    }
  
  plg_vector2array(tmp, boundaries, nboundaries);
  
  status=plg_destroy_entries(tmp);
  
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_split(plg_t target, vector<plg_t> & out, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  reorganize polygon with in multiple polygons of similar code
  
  codes:
  
  - 'T' : standard coastal (rigid) limit
  - 'M' : standard ocean   (open)  limit
  - 'X' : locked   ocean   (open)  limit

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  size_t npassed=0, newlines=0, nlines=target.npt-1;
  
  size_t k=0;
  while (npassed!=nlines) {
    if(npassed>nlines) TRAP_ERR_EXIT(-1,"bug \n",__func__);
    if(target.flag[k]=='M') {
      npassed++;
      size_t l=(k+1) % (nlines);
      if (npassed<nlines) {
        while (target.flag[l]=='M') {
          l=(l+1) % (nlines);
          npassed++;
          if (npassed==nlines) break;
          }
        }
/*------------------------------------------------------------------------------
      duplicate open sections */
      plg_t p;
      p.duplicate(target,k,l);
      out.push_back(p);
      k=l;
      newlines++;
      }
    else if(target.flag[k]=='X') {
      npassed++;
      size_t l=(k+1) % (nlines);
      if (npassed<nlines) {
        while (target.flag[l]=='X') {
          l=(l+1) % (nlines);
          npassed++;
          if (npassed==nlines) break;
          }
        }
/*------------------------------------------------------------------------------
      duplicate open sections */
      plg_t p;
      p.duplicate(target,k,l);
      out.push_back(p);
      k=l;
      newlines++;
      }
    else {
      npassed++;
      size_t l=(k+1) % (nlines);
      if (npassed<nlines) {
        while (target.flag[l]=='T') {
          l=(l+1) % (nlines);
          npassed++;
          if (npassed==nlines) break;
          }
        }

/*------------------------------------------------------------------------------
      extract rigid section*/
      plg_t p;
      p.duplicate(target,k,l);
      out.push_back(p);
      newlines++;
      if(npassed>nlines) TRAP_ERR_EXIT(-1,"bug \n",__func__);
      k=l;
      }
    }
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  vector<plg_t> plg_split(const vector<plg_t> & p, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------
  reorganize polygon with in multiple polygons of similar code */
{
  int status;

  vector<plg_t> out;
  
  for(size_t s=0;s<p.size();s++) status=plg_split(p[s], out, debug);
  
  return(out);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_write_descriptor (string rootname, vector<plg_t> & boundaries, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  FILE  *fp=NULL;
  int  m,n;
  int  status;
  vector<plg_t> polygons, rigids;
  int first;
  plg_desc_t descriptor;
  string DescriptorName=rootname+".desc";
  string PolygonesName =rootname+"-splitted.scan";
  string RigidsName    =rootname+"-rigid.scan";

  fp = fopen (DescriptorName.c_str(), "w");
  if(fp==0) TRAP_ERR_RETURN(errno,1,"Could not open "+DescriptorName+" (%s)\n",strerror(errno));

/*------------------------------------------------------------------------------
  number of blocks*/
  descriptor.nblock=boundaries.size();

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  split rigid and open segments in separate polygons
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  polygons=plg_split(boundaries, debug);
  
/*------------------------------------------------------------------------------
  identify rigid pieces and store for further assembly if required */
  for(int s=0;s<polygons.size();s++) {
    if(polygons[s].flag[0]=='T') rigids.push_back(polygons[s]);
    }
  status=plg_save(RigidsName.c_str(), PLG_FORMAT_SCAN, rigids);

  rigids.clear();
  
/*------------------------------------------------------------------------------
  count the number of segments in blocks*/
  fprintf(fp,"%s\n",PolygonesName.c_str());
  status=plg_save(PolygonesName.c_str(), PLG_FORMAT_SCAN, polygons);

  descriptor.b=new plg_block_t[descriptor.nblock];

  first=0;
  for(n=0;n<descriptor.nblock;n++) {
    vector<plg_t> p;
    status=plg_split(boundaries[n],p,false);
    descriptor.b[n].segments.n=p.size();
    fprintf(fp,"%6d\n",n+1);
    for(m=0;m<descriptor.b[n].segments.n;m++) {
      p[m].id=first+m+1;
      fprintf(fp,"%6d\n",p[m].id);
      fprintf(fp,"%c\n",p[m].flag[0]);
      if(m==0) {
        int orientation=sign(boundaries[n].area());
        fprintf(fp,"%6d\n",orientation);
        if(n==0) fprintf(fp,"%c\n",'E');
        else fprintf(fp,"%c\n",'I');
        }
      }
    fprintf(fp,"%s\n","End_of_Block");
    first+=p.size();
    }

  for(int s=0;s<polygons.size();s++) polygons[s].destroy();
  polygons.clear();
  
  
  fclose(fp);
  return(0);
}

