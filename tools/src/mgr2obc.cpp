


/***************************************************************************
T-UGO tools, 2006-2011

  Unstructured Ocean Grid initiative
  
E-mail: florent.lyard@legos.obs-mip.fr
***************************************************************************/
/**
\file

\author  Florent Lyard      LEGOS/CNRS, Toulouse, France
\author  Yoann Le Bars      LEGOS, Toulouse, France (PhD)
\author  Laurent Roblou     LEGOS/CNRS, Toulouse, France
\author  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
\author  Frédéric Dupont    Université de Laval à Québec, Canada
\author  Clement MAYET      LEGOS, Toulouse, France (PhD)

VERSION :

\date  27/04/2011 : Clement MAYET: change header (doxygen like)

\brief


***************************************************************************/

#define MAIN_SOURCE

#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cmath>

#include "tools-structures.h"

#include "functions.h"
#include "map.h"
#include "fe.h"
#include "geo.h"
#include "mgr.h"

extern int get_openlimits(const char *,const char *, double **, double **, int *);
extern int load_HamonicOBCs(const char *filename, spectrum_t & spectrum, vector<tidaldata_t> & TidalData);

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int get_positions(const char *positions, double **x, double **y, int *count)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int nitems, status=0;
  FILE *in;
  
  vector<double> xtmp,ytmp;
  
  double xx,yy;
  
  in=fopen(positions,"r");
  
  if(in==0) return(-1);
  
  while (nitems!=0) {
    nitems=fscanf(in,"%lf %lf",&xx,&yy);
    if(nitems==0) break;
    xtmp.push_back(xx);
    ytmp.push_back(yy);
    }
  
  *count=xtmp.size();
  
  exitIfNull(
    (*x)=new double[*count]
    );
  exitIfNull(
    (*y)=new double[*count]
    );

  for(size_t n=0;n<xtmp.size();n++) {
    (*x)[n]=xtmp[n];
    (*y)[n]=ytmp[n];
    }
  
  fclose(in);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int get_positions(const char *belfile,const char *meshfile, double **x, double **y, int *count)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  if(meshfile!=0) {
    status=get_openlimits(belfile, meshfile, x, y, count);
    }
  else {
    status=get_positions(belfile, x, y, count);
    }
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int parse_GridOptions(const  string options, char* & filename, metagrid_t & meta)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  vector<string> tokens,keys,values;
  int k;
  string delimiter=" ";
  
  filename=0;
  
  delimiter=" ";
  tokens=string_split(options, delimiter);
  
  delimiter="=";
  for(k=0;k<tokens.size();k++) {
    vector<string> tmp=string_split(tokens[k], delimiter);
    keys.push_back(tmp[0]);
    values.push_back(tmp[1]);
    }
    
  for(k=0;k<keys.size();k++) {
    if(keys[k]=="file") {
      filename=strdup(values[k].c_str());
      continue;
      }
    if(keys[k]=="vlon_t") {
      meta.z_gnames.vlon=strdup(values[k].c_str());
      continue;
      }
    if(keys[k]=="vlat_t") {
      meta.z_gnames.vlat=strdup(values[k].c_str());
      continue;
      }
    if(keys[k]=="vmask_t") {
      meta.z_gnames.vmask=strdup(values[k].c_str());
      continue;
      }
    if(keys[k]=="vtopo_t") {
      meta.z_gnames.vtopo=strdup(values[k].c_str());
      continue;
      }
    if(keys[k]=="vlon_f") {
      meta.f_gnames.vlon=strdup(values[k].c_str());
      continue;
      }
    if(keys[k]=="vlat_f") {
      meta.f_gnames.vlat=strdup(values[k].c_str());
      continue;
      }
    if(keys[k]=="vmask_f") {
      meta.f_gnames.vmask=strdup(values[k].c_str());
      continue;
      }
    if(keys[k]=="vtopo_f") {
      meta.f_gnames.vtopo=strdup(values[k].c_str());
      continue;
      }
    }
    
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double *lon,*lat;
  double d, dmin, dmax=4.0;
  int count,nndes,status,*elts,found;
  int i,j,k,l,m,n;
  FILE *file;
  FILE *out;
  char *keyword;
  char *meshfile=NULL,*belfile=NULL,*output=NULL,*path=NULL;
  vector<string> mgrfiles;
  char *wave[1024],filename[1024];
  fcomplex *z,*u,*v;
  float a[3000],G[3000];  /// HERE !!!
  float **Ha=NULL,**Hg=NULL,**Ua=NULL,**Ug=NULL,**Va=NULL,**Vg=NULL,mask;
  mesh_t mesh;
  int nwave=0,analysis;
  spectrum_t spectrum;
  string s,echofile;
  vector<mgr_t> mgr;
  int nmgr;
  int nkeep,*keep,w;
  FILE *echo;
  int element=FE_TRIANGLE;
  string grid_options="";
  bool debug=false;
  double tag=NAN;
  bool append=false;
  spectrum_t          imported_spectrum;
  vector<tidaldata_t> imported_TidalData;
  

  fct_echo( argc, argv);

  n=1;
  while (n < argc) {
    keyword=strdup(argv[n]);
    if(strstr(argv[n],"--quadrangle")!=0){
      element=FE_QUADRANGLE;
      n++;
      continue;
      }
    if(strstr(argv[n],"--grid")!=0){
      grid_options=strstr(argv[n],"--grid")+7;
      n++;
      continue;
      }
    if(strstr(argv[n],"-dmax")!=0){
      sscanf(argv[n+1],"%lf",&dmax);
      n++;
      n++;
      continue;
      }
    if(strstr(argv[n],"--append")!=0){
      append=true;
      n++;
      continue;
      }
    if(strstr(argv[n],"--debug")!=0){
      debug=true;
      n++;
      continue;
      }
    switch (keyword[0]) {
      case '-':
        switch (keyword[1]) {
        case 'm' :
          meshfile= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'b' :
          belfile= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'd' :
        case 'g' :
          mgrfiles.push_back(argv[n+1]);
          n++;
          n++;
          break;

        case 'o' :
          output= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'p' :
          path= strdup(argv[n+1]);
          n++;
          n++;
          break;

        default:
          __OUT_BASE_LINE__("unknown option %s\n",keyword);
          exit(-1);
        }
        break;

      default:
        wave[nwave]= strdup(argv[n]);
        printf("input wave=%s\n",wave[nwave]);
        spectrum.n=nwave+1;
        nwave++;
        n++;
        break;
      }
      free(keyword);
    }
    
  if (output==0) {
    output=strdup("tides.obc");
    }
  if(append) status=load_HamonicOBCs(output, imported_spectrum, imported_TidalData);

/* *----------------------------------------------------------------------------
  Load tide gauge database */
  if(mgrfiles.size()==NULL) goto error;
  
  for(int k=0; k<mgrfiles.size(); k++) {
    nmgr=mgr_load(mgrfiles[k].c_str(), mgr);
    }
    
  spectrum.waves=new tidal_wave[spectrum.n];
  for (k=0;k<nwave;k++) {
    strcpy(spectrum.waves[k].name,wave[k]);
    }

  if(path ==NULL) path=strdup(".");

  if(grid_options!="") {
    char *filename;
    metagrid_t meta;
    grid_t grid;
    mesh_t mesh;
    string rootname="mgr2obc";
    bool debug=true;
    
    status= parse_GridOptions(grid_options, filename, meta);
       
    status=quadrangle_ImportStructured(filename, 0, mesh, meta, tag, rootname, "", debug);
    if(status !=0) goto error;
   
    meshfile=strdup("atlas2obc-quadrangle.nei");
    status=fe_savemesh(meshfile,MESH_FILE_FORMAT_TRIGRID, mesh);
    if(status !=0) goto error;
    }
  
//  status=get_openlimits(belfile,meshfile,&lon,&lat,&count);
  
  if(imported_TidalData.size()==0) {
    status=get_positions(belfile,meshfile,&lon,&lat,&count);
    if(status!=0) goto error;
    }
  else {
    count=imported_TidalData.size();
    lon=new double[count];
    lat=new double[count];
    for(n=0;n<count;n++) {
      lon[n]=imported_TidalData[n].lon;
      lat[n]=imported_TidalData[n].lat;
      }
    }

  elts=(int *) malloc(count*sizeof(int));
  keep= new int[nmgr];

  nkeep=0;
  for (m=0;m<nmgr;m++) {
    dmin=1.e+10;
    for(k=0; k<count; k++) {
      d=geo_km(lon[k],lat[k],mgr[m].loc.lon,mgr[m].loc.lat);
      if(d<dmin) {
        elts[k]=m;
        dmin=d;
        }
      }
    if(dmin<dmax) {
      keep[nkeep]=m;
      nkeep++; 
      printf("TG %s kept (d=%lf, dmax=%lf)\n",mgr[m].name ,dmin,dmax);
      }
    else {
      if(debug) printf("TG %d too far (d=%lf, dmax=%lf)\n",m,dmin,dmax);
      }
    }
    
  count=nkeep;
    
  Ha=new float*[nwave];
  Hg=new float*[nwave];
  Ua=new float*[nwave];
  Ug=new float*[nwave];
  Va=new float*[nwave];
  Vg=new float*[nwave];

  for (k=0;k<nwave;k++) {
    Ha[k]=new float[count];
    Hg[k]=new float[count];
    Ua[k]=new float[count];
    Ug[k]=new float[count];
    Va[k]=new float[count];
    Vg[k]=new float[count];
    }

  z= new fcomplex[nkeep];
  u= new fcomplex[nkeep];
  v= new fcomplex[nkeep];

  for (k=0;k<nwave;k++) {
/*------------------------------------------------------------------------------
    Tidal heights and currents at open limits */
    sprintf(filename,"%s.obc",wave[k]);
    out=fopen(filename,"w");
    fprintf(out,"%s\n",wave[k]);
    for(i=0; i<nkeep; i++) {
      m=keep[i];
      w=mgr[m].wave_index(wave[k]);
      if(w==-1) {
        __TRAP_ERR_EXIT__(-1,"cannot find %s wave in %d\n",wave[k],m);
        }
      z[i]=fct_polar2complex(mgr[m].data[w].amp,360.0-mgr[m].data[w].phi);
/*------------------------------------------------------------------------------
      convert elevation amplitude in m*/
      a[k]=abs(z[i])*1.0;
      G[k]=360.0-arg(z[i])*r2d;
      if(G[k]>360.0) G[k]-=360.0;
      if(G[k]<0.0)   G[k]+=360.0;
      Ha[k][i]=a[k];
      Hg[k][i]=G[k];
      fprintf(out,"%8f %8f %f %f\n",mgr[m].loc.lon,mgr[m].loc.lat,a[k],G[k]);
      }
    fclose(out);
    }

  if ((out=fopen(output,"w")) ==NULL) {
    __ERR_BASE_LINE__("cannot open obc file : %s \n",output);
    exit(-1);
    }
    
//   fprintf(out,"%d %d %s\n",count,nwave,"M");
  fprintf(out,"%d %d %s\n",count,nwave+imported_spectrum.n,"M");
 
  for (k=0;k<imported_spectrum.n;k++) {
/* *-----------------------------------------------------------------------------
    tidal heights and currents at open limits*/
    fprintf(out,"%s\n",imported_spectrum.waves[k].name);
    for(i=0; i<count; i++) {
/*---------------------------------------------------------------------
      elevation*/
      a[0]=imported_TidalData[i].Ha[k];
      G[0]=imported_TidalData[i].Hg[k];
      fprintf(out,"%8f %8f %f %f\n",lon[i],lat[i],a[0],G[0]);
      }
    }
    
  for (k=0;k<nwave;k++) {
/* *-----------------------------------------------------------------------------
    tidal heights and currents at open limits*/
    fprintf(out,"%s\n",wave[k]);
    for(i=0; i<nkeep; i++) {
      m=keep[i];
      w=mgr[m].wave_index(wave[k]);
      if(w==-1) {
        __TRAP_ERR_EXIT__(-1,"cannot find %s wave in %d\n",wave[k],m);
        }
      a[k]=Ha[k][i];
      G[k]=Hg[k][i];
      fprintf(out,"%8f %8f %f %f\n",mgr[m].loc.lon,mgr[m].loc.lat,a[k],G[k]);
      }
    }
  fclose(out);

  delete[] z;
  delete[] u;
  delete[] v;

end: 
  printf("end of %s ... \n",argv[0]);
  free(elts);
  __ERR_BASE_LINE__("exiting\n");
  exit(0);
  
error:
 __OUT_BASE_LINE__("error detected, quit ... \n");
  exit(-1);
}
