#include <config.h>

#include <stdio.h>
#include <fstream>

#include "tools-structures.h"

#include "polygones.h"
#include "statistic.h"
#include "geo.h"
#include "xyz.h"

// #if HAVE_SHAPEFIL_H || HAVE_LIBSHP
// #include <shapefil.h>
// #endif

#if HAVE_LIBSHP
#if HAVE_SHAPEFIL_H
#include <shapefil.h>
#elif HAVE_LIBSHP_SHAPEFIL_H
#include <libshp/shapefil.h>
#endif
#endif

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int plg_load_histolitt(const char *filename, plg_t **polygones, int *npolygones)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  plg_t *plg_polygones;
  int k,i,npoly=0,nitems;
  FILE *polyfile;
  char line[1024],*s;
  vector<int> npt;
  int count;
  
  *polygones=0;
  *npolygones=0;

  polyfile=fopen(filename,"r");
  if(polyfile==NULL) {
    TRAP_ERR_RETURN(errno,1,"Can not open %s (%d: %s)\n",filename,errno,strerror(errno));
    }

  count=0;
  npoly=0;
  do {
    s=fgets(line, 1024, polyfile);
    if (s == 0) break;
    if(s[0]=='>' or s[0]=='\r') {
      if(npoly!=0) npt.push_back(count);
      npoly++;
      count=0;
      }
    else {
      count++;
      }
    } while(!feof(polyfile));

  if(count!=0) npt.push_back(count);
  rewind(polyfile);

  plg_polygones=new plg_t[npoly];

  count=0;
  size_t nlines=0;
  for(k=0;k<npoly;k++) {
    s=fgets(line, 1024, polyfile);
    if (s == 0) break;
    nlines++;
    if(s[0]=='>' or s[0]=='\n') count++;
    plg_t *p=&plg_polygones[k];
    p->npt=npt[k];
    p->init(PLG_INIT_SEPARATE);
    p->id=k;
    for(i=0;i<p->npt;i++) {
      fgets(line, 1024, polyfile);
      nitems=sscanf(line,"%lf %lf",&(p->x[i]),&(p->y[i]));
      if(nitems!=2) {
        p->npt--;
        continue;
        }
      nlines++;
      if(p->x[i]==0. and p->y[i]==0.0) {
        printf("read error at line=%d\n",nlines);
        }
      p->t[i]=p->x[i];
      p->p[i]=p->y[i];
      }
    }

  fclose(polyfile);

  *polygones=plg_polygones;
  *npolygones=npoly;

  return(-0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int plg_load_KML(const char *filename, plg_t **polygones, int *npolygones)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  plg_t *plg_polygones;
  int run,npoly,npt,nitems;
  ifstream file;
  string word;
  vector<int> npts;
  plg_t *p=0;
  bool inPlg;

  file.open(filename);
  if(file.fail()) {
    TRAP_ERR_RETURN(errno,1,"Can not open %s (%d: %s)\n",filename,errno,strerror(errno));
    }
  
/*------------------------------------------------------------------------------
  scan twice */
  for(run=0;true;run++){
    
    inPlg=false;
    npoly=0;
    
    while(true) {
      file>>word;
      
      if(not file.good())
        break;
      
      if(not inPlg){
        
/*------------------------------------------------------------------------------
        entering a polygone zone */
        if(word.find("<coordinates>")!=string::npos){
          inPlg=true;
          npt=0;
          
          if(run>0){
            p=&plg_polygones[npoly];
            p->init(npts[npoly],PLG_INIT_SHARED);
            }
          
          }
        
        continue;
        }
      
/*------------------------------------------------------------------------------
      leaving a polygone zone */
      if(word.find("</coordinates>")!=string::npos){
        inPlg=false;
        
        if(run<1){
          npts.push_back(npt);
          }
        
        npoly++;
        continue;
        }
      
/*------------------------------------------------------------------------------
      inside a polygone zone */
      if(run>0){
        nitems=sscanf(word.c_str(),"%lf%*[ ,]%lf",&(p->x[npt]),&(p->y[npt]));
        }
      
      npt++;
      
      }
    
    if(run<0 and inPlg){
      STDERR_BASE_LINE_FUNC("EOF reached while within a <coordinates> XML element in %s\n",filename);
      npts.push_back(npt);
      npoly++;
      }
    
    if(run>0)
      break;
    
/*------------------------------------------------------------------------------
    scanned once */
    file.clear();
    file.seekg(0,file.beg);
    
    plg_polygones=new plg_t[npoly];
    }
  
  *polygones=plg_polygones;
  *npolygones=npoly;
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int plg_load_scan(const char *filename, plg_t **polygones, int *npolygones)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  plg_t *plg_polygones;
  int k,tmp,i,npt,npoly=0,nitems;
  FILE *polyfile;
  char line[1024],a;

  polyfile=fopen(filename,"r");
  if(polyfile==NULL) {
    TRAP_ERR_RETURN(errno,1,"Can not open %s (%d: %s)\n",filename,errno,strerror(errno));
    }

  do {
    nitems=fscanf(polyfile,"%d",&tmp);
    if (nitems !=1 ) break;
    fscanf(polyfile,"%d",&npt);
    do { a=fgetc(polyfile);   }  while ((a != '\n') && !feof(polyfile));
    for(i=0;i<npt;i++) {
//      fscanf(polyfile,"%d %lf %lf",&tmp,&x,&y);
      fgets(line, 1024, polyfile);
      }
    npoly++;
    } while(!feof(polyfile));

  rewind(polyfile);

  plg_polygones=new plg_t[npoly];

  for(k=0;k<npoly;k++) {
    if ( fscanf(polyfile,"%d",&tmp) !=1 ) break;
    plg_t *p=&plg_polygones[k];
    fscanf(polyfile,"%d",&p->npt);
    do { a=fgetc(polyfile);   }  while ((a != '\n') && !feof(polyfile));
    p->init(PLG_INIT_SEPARATE);
    p->id=tmp;
    for(i=0;i<p->npt;i++) {
      fgets(line, 1024, polyfile);
      nitems=sscanf(line,"%d %lf %lf",&tmp,&(p->x[i]),&(p->y[i]));
      p->t[i]=p->x[i];
      p->p[i]=p->y[i];
      }
    }

  fclose(polyfile);

  *polygones=plg_polygones;
  *npolygones=npoly;

  return(-0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_write_scan (const char *filename, const plg_t *polygones, int npolygones)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  FILE  *fp;
  int  k, ns;

  fp = fopen (filename, "w");

  if (fp == 0)
    TRAP_ERR_RETURN(errno,"fopen(\"%s\",\"w\") error (%d %s)\n",filename,errno,strerror(errno));

  for(ns=0;ns<npolygones;ns++) {
    if(polygones[ns].npt==0) continue;
    fprintf(fp,"%6d %d\n",ns+1,polygones[ns].npt);
    if(polygones[ns].flag==0) {
      for(k=0;k<polygones[ns].npt;k++) {
        fprintf(fp,"%6d %13.7lf %13.7lf\n",k+1,polygones[ns].t[k],polygones[ns].p[k]);
        }
      }
    else {
      for(k=0;k<polygones[ns].npt;k++) {
 //       fprintf(fp,"%6d %lf %lf %d\n",k+1,polygones[ns].t[k],polygones[ns].p[k],polygones[ns].flag[k]);
        fprintf(fp,"%6d %13.7lf %13.7lf\n",k+1,polygones[ns].t[k],polygones[ns].p[k]);
        }
      }
    }

  fclose (fp);

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_write_pocformula (const char *filename, const plg_t *polygones, int npolygones)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  FILE  *fp;
  int  k, ns;

  fp = fopen (filename, "w");

  if (fp == 0)
    TRAP_ERR_RETURN(errno,"fopen(\"%s\",\"w\") error (%d %s)\n",filename,errno,strerror(errno));

  for(ns=0;ns<npolygones;ns++){
    if(polygones[ns].npt==0) continue;
    
    fprintf(fp,"plg%d =",ns+1);
    
    for(k=0;k<polygones[ns].npt;k++) {
      if(k>0)
        fprintf(fp," &");
      fprintf(fp," %.7f+%.7fj",k+1,polygones[ns].t[k],polygones[ns].p[k]);
      }
    
    fprintf(fp,"\n");
    }

  fclose (fp);

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_load_gshhs(const char *filename, plg_t **polygones, int *npolygones)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  plg_t *plg_polygones;
  int status,npt,npoly=0, fmt;

  status=plg_inquire_gshhs (filename, &npoly,&npt,&fmt);

  plg_polygones=new plg_t[npoly];

  status=plg_read_gshhs (filename, plg_polygones, fmt);

  *polygones=plg_polygones;
  *npolygones=npoly;

  return(-0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int plg_load_cst(const char *filename, plg_t **polygones, int *npolygones)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  plg_t *plg_polygones;
  int status,npt,npoly=0;

  status=plg_inquire_cst (filename, &npoly,&npt);

  plg_polygones=new plg_t[npoly];

  status=plg_read_cst (filename, plg_polygones);

  *polygones=plg_polygones;
  *npolygones=npoly;

  return(-0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int plg_load_xiso(const char *filename, plg_t **polygones, int *npolygones)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  plg_t *plg_polygones;
  int status,npt,npoly=0;

  status=plg_inquire_xiso (filename, &npoly,&npt);

  plg_polygones=new plg_t[npoly];

  status=plg_read_xiso (filename, plg_polygones);

  *polygones=plg_polygones;
  *npolygones=npoly;

  return(-0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int plg_load_boundaries(const char *filename, plg_t **polygones, int *npolygones)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  plg_t *plg_polygones;
  int status,npt,npoly=0;

  status=plg_inquire_boundaries (filename, &npoly,&npt);

  plg_polygones=new plg_t[npoly];

  status=plg_read_boundaries (filename, plg_polygones);

  *polygones=plg_polygones;
  *npolygones=npoly;

  return(-0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_load_shp(const char *filename, plg_t **polygones, int *npolygones, int mode)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  plg_t *plg_polygones;
  int status,npoly=0,nblocks=0,shapetype;
  
#if HAVE_SHAPEFIL_H || HAVE_LIBSHP
/**----------------------------------------------------------------------------
  one block may contain more than one segement */
  status=plg_inquire_shp (filename, &npoly, &nblocks, &shapetype);
  if(status!=0){
    TRAP_ERR_RETURN(status,1,"plg_inquire_shp(\"%s\",,,) error %d: %s\n",filename,status,strerror(status));
    }
  
  if(shapetype!=SHPT_POLYGON && shapetype!=SHPT_POLYGONZ && shapetype!=SHPT_POINT && shapetype!=SHPT_ARC && shapetype!=SHPT_ARCZ && shapetype!=SHPT_NULL)
    TRAP_ERR_RETURN(-1,1,"Shape %d not supported\n",shapetype);

  if(npoly==0) TRAP_ERR_RETURN(-1,1,"no polygons found\n");
  
  plg_polygones=new plg_t[npoly];

  status=plg_read_shp(filename, plg_polygones, nblocks, mode);

  *polygones=plg_polygones;
  *npolygones=npoly;
  return status;
#else
  status=ENOEXEC;

  STDERR_BASE_LINE("please install developpment files for libshp\n");
  
  return(status);
#endif
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int plg_load_netcdf(const char *filename, plg_t **polygones, int *npolygones)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  plg_t *plg_polygones;
  int status,npoly=0;

  *polygones=0;
  *npolygones=0;
  
  status=plg_inquire_netcdf (filename, &npoly);
  if(status!=0) return(-1);

  plg_polygones=new plg_t[npoly];

  plg_read_netcdf (filename, plg_polygones,npoly);

  *polygones=plg_polygones;
  *npolygones=npoly;

  return(-0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_load_XY(const char *filename, plg_t **polygons, int *npolygons)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  plg_t *p;
  int status;
  double *x, *y, *z, mask;
  int k,ndata;
  char *proj4_options=0;
  string header="";
  
  *polygons=0;
  *npolygons=0;
  
  status=xyz_loadraw (filename, header, proj4_options, x, y, z, &mask, ndata);
  if(status!=0) return(-1);

  p=new plg_t[1];
  p[0].npt=ndata;

  for(k=0;k<p[0].npt;k++) {
    p[0].t=x;
    p[0].p=y;
    }
  
  delete[] z;
  
  *polygons=p;
  *npolygons=1;

  return(-0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int plg_load_MAPINFO(const char *filename, plg_t **polygons, int *npolygons)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  vector<plg_t> p;
  int nitems;
  vector<double> x, y;
  double xx,yy;
  int k;
  FILE *in;
  char *line;
  
  in=fopen(filename,"r");
  
  *polygons=0;
  *npolygons=0;
  
  line=new char[1024];
  
  while(!feof(in)) {
    plg_t q;
    nitems=0;
    fgets(line,1024,in);
    if(feof(in)) break;
    if(strstr(line,"C element")==0) continue;
    fgets(line,1024,in);
    while(strstr(line,"C")==0) {
      sscanf(line,"%lf %lf",&xx,&yy);
      x.push_back(xx);
      y.push_back(yy);
      fgets(line,1024,in);
      if(feof(in)) break;
      }
    q.init(x.size(),PLG_INIT_SEPARATE);
    for(k=0;k<q.npt;k++) {
      q.t[k]=x[k];
      q.p[k]=y[k];
      q.x[k]=x[k];
      q.y[k]=y[k];
      }
    p.push_back(q);
    x.clear();
    y.clear();
    }
    
  fclose(in);
  
  delete[] line;
  
  plg_array_t a=plg_vector2array(p);
  
  *polygons=a.p;
  *npolygons=a.n;

  return(-0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_print_formats()

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  const char *c;

  #define PLG_FORMAT(x) c=#x;printf("  %d %s\n",x,&c[11])
  PLG_FORMAT(PLG_FORMAT_BIN);
  PLG_FORMAT(PLG_FORMAT_SCAN);
  PLG_FORMAT(PLG_FORMAT_GSHHS);
  PLG_FORMAT(PLG_FORMAT_XISO);
  PLG_FORMAT(PLG_FORMAT_BOUNDARIES);
  PLG_FORMAT(PLG_FORMAT_SHP);
  PLG_FORMAT(PLG_FORMAT_NETCDF);
  PLG_FORMAT(PLG_FORMAT_XY);
  PLG_FORMAT(PLG_FORMAT_MAPINFO);
  PLG_FORMAT(PLG_FORMAT_HISTOLITT);
  PLG_FORMAT(PLG_FORMAT_POCFORMULA);
  #undef PLG_FORMAT
  
  return 0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_print(const plg_t & p)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  for(int k=0;k<p.npt;k++) {
    printf("%d %lfE %lfN %lfm %lfm\n",p.t[k],p.p[k],p.x[k],p.y[k]);
    }
  
  return 0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_format_from_name(const string &format_name)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  const char *c;

  #define PLG_FORMAT(x) c=#x;if(format_name==&c[11])return x
  PLG_FORMAT(PLG_FORMAT_BIN);
  PLG_FORMAT(PLG_FORMAT_ASCII);
  PLG_FORMAT(PLG_FORMAT_SCAN);
  PLG_FORMAT(PLG_FORMAT_GSHHS);
  PLG_FORMAT(PLG_FORMAT_MIF);
  PLG_FORMAT(PLG_FORMAT_VCT00);
  PLG_FORMAT(PLG_FORMAT_XISO);
  PLG_FORMAT(PLG_FORMAT_BOUNDARIES);
  PLG_FORMAT(PLG_FORMAT_SHP);
  PLG_FORMAT(PLG_FORMAT_NETCDF);
  PLG_FORMAT(PLG_FORMAT_XY);
  PLG_FORMAT(PLG_FORMAT_MAPINFO);
  PLG_FORMAT(PLG_FORMAT_HISTOLITT);
  #undef PLG_FORMAT
  
  return PLG_FORMAT_UNKNOWN;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

string plg_format_extension(int format)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  string s;
  
  switch(format){
    case PLG_FORMAT_BIN:
      s=".cst";
      break;
    case PLG_FORMAT_SHP:
      s=".shp";
      break;
    case PLG_FORMAT_NETCDF:
      s=".nc";
      break;
    }
  
  return s;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int plg_find_format(const string &filename)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  if(strrncmp(filename,".cst")==0)
    return PLG_FORMAT_BIN;
  
  if(strrncmp(filename,".shp")==0)
    return PLG_FORMAT_SHP;
  
  if(strrncmp(filename,".nc")==0)
    return PLG_FORMAT_NETCDF;
  
  if(strrncmp(filename,".fr")==0)
    return PLG_FORMAT_BOUNDARIES;
  
  if(strrncmp(filename,".plg")==0 ||
     strrncmp(filename,".scan")==0)
    return PLG_FORMAT_SCAN;
  
  if(strrncasecmp(filename,".kml")==0)
    return PLG_FORMAT_KML;
  
  if(strrncmp(filename,".pocformula")==0)
    return PLG_FORMAT_POCFORMULA;
  
  return PLG_FORMAT_UNKNOWN;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_load(const string &filename, int format,plg_t **polygones,int *npolygones)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// loads a polygone
/**
\param filename
\param format if #PLG_FORMAT_UNKNOWN, call plg_find_format()
\param **polygones
\param *npolygones number of polygones loaded
\returns 0 or a NetCDF error code if error
*/
/*----------------------------------------------------------------------------*/
{
  const char *f=filename.c_str();
  int status;
  int mode;
  
  if(format==PLG_FORMAT_UNKNOWN){
    format=plg_find_format(filename);
    }

  switch (format) {
    case PLG_FORMAT_BIN:
      status=plg_load_cst(f, polygones, npolygones);
      break;

    case PLG_FORMAT_SCAN:
      status=plg_load_scan(f, polygones, npolygones);
      break;

    case PLG_FORMAT_GSHHS:
      status=plg_load_gshhs(f, polygones, npolygones);
      break;

    case PLG_FORMAT_MIF:
//      status=plg_load_mif(f, polygones, npolygones);
      break;

    case PLG_FORMAT_XISO:
      status=plg_load_xiso(f, polygones, npolygones);
      break;

    case PLG_FORMAT_BOUNDARIES:
      status=plg_load_boundaries(f, polygones, npolygones);
      break;

    case PLG_FORMAT_SHP:
      mode=PLG_INIT_SEPARATE;
      status=plg_load_shp(f, polygones, npolygones, mode);
      break;

    case PLG_FORMAT_NETCDF:
      status=plg_load_netcdf(f, polygones, npolygones);
      break;

    case PLG_FORMAT_XY:
      status=plg_load_XY(f, polygones, npolygones);
      break;

    case PLG_FORMAT_MAPINFO:
      status=plg_load_MAPINFO(f, polygones, npolygones);
      break;

    case PLG_FORMAT_HISTOLITT:
      status=plg_load_histolitt(f, polygones, npolygones);
      break;

    case PLG_FORMAT_KML:
      status=plg_load_KML(f, polygones, npolygones);
      break;

   default:
      return NC_EBADTYPE;/* data format error */

    }

  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_load(const string &filename, int format, vector<plg_t> & polygons)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
{
  int status;
  plg_t *tmp=0;
  int npolygones=0;
  
  status=plg_load(filename, format, &tmp, &npolygones);

  if(status==0) copy(&polygons,tmp, npolygones);
  
  deletep(&tmp);
  
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_load(const string &filename, vector<plg_t> & polygons)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// load a polygone, detecting the format
/**
\param filename
\param **polygones
\param *npolygones number of polygones loaded
\returns 0 or a NetCDF error code if error
*/
/*----------------------------------------------------------------------------*/
{
  int status;
  
  status=plg_load(filename,PLG_FORMAT_UNKNOWN,polygons);

  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_load(const string &filename,plg_t **polygones,int *npolygones)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// load a polygone, detecting the format
/**
\param filename
\param **polygones
\param *npolygones number of polygones loaded
\returns 0 or a NetCDF error code if error
*/
/*----------------------------------------------------------------------------*/
{
  int status;
  
  status=plg_load(filename,PLG_FORMAT_UNKNOWN,polygones,npolygones);

  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int plgF_load(const char *filename, int & format, double *rlon, double *rlat, int *first, int *last, int & ns)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  plg_t *polygones;
  int npolygones;
  
  status=plg_load(filename, format, &polygones, &npolygones);
  
  if(status!=0) return(status);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int plg_save(const char *filename, int format, const plg_t *polygones, int npolygones)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;

  if(format==PLG_FORMAT_UNKNOWN){
    format=plg_find_format(filename);
    }

  switch (format) {
    case PLG_FORMAT_BIN:
      status=plg_write_cst(filename, polygones, npolygones);
      break;

    case PLG_FORMAT_SCAN:
      status=plg_write_scan(filename, polygones, npolygones);
      break;

    case PLG_FORMAT_POCFORMULA:
      status=plg_write_pocformula(filename, polygones, npolygones);
      break;

    case PLG_FORMAT_GSHHS:
      status=-1;
      break;
      
    case PLG_FORMAT_SHP:
      status=plg_write_shp(filename, polygones, npolygones);
      break;

    case PLG_FORMAT_NETCDF:
      status=plg_write_netcdf(filename, polygones, npolygones);
      break;

    default:
      return NC_EBADTYPE;

    }

  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int plg_save(const char *filename, const plg_t & polygone)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status,format;
  
  format=plg_find_format(filename);

  status=plg_save(filename,format,&polygone,1);

  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int plg_save(const char *filename, int format, const vector<plg_t> & polygones)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;

  plg_t *tmp=copy(polygones);
  
  status=plg_save(filename, format, tmp, polygones.size());
  
  delete[] tmp;

  return status;
}
