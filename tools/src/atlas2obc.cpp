
/*******************************************************************************

  T-UGO tools, 2006-2013

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

<!-- USE FIND AND REPLACE TO EDIT THIS LINE, SO THAT print_help IS ALSO UPDATED -->
\brief Creates an open boundary file for TUGOm.

<!-- A LINK TO main() or print_help() WILL NOT LINK TO THE RIGHT SOURCE ! -->
See the main function for how this works
and the print_help function for how to use this.
*/
/*----------------------------------------------------------------------------*/

#define MAIN_SOURCE

#include "version-macros.def" //for VERSION and REVISION

#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cmath>

#include "tools-structures.h"

#include "functions.h"
#include "map.h"
#include "fe.h"
#include "tides.h"
#include "archive.h"
#include "poc-netcdf.hpp"

extern int get_openlimits(const char *,const char *, double **, double **, int *);
extern int load_HamonicOBCs(const char *filename, spectrum_t & spectrum, vector<tidaldata_t> & TidalData);

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int get_positions(const char *positions, double **x, double **y, int *count)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int nitems, status=0;
  FILE *in;
  char line[1024];
  
  vector<double> xtmp,ytmp;
  
  double xx,yy;
  
  in=fopen(positions,"r");
  
  if(in==0) return(-1);
  
  fgets(line, 1024, in);

  bool latlon=(strstr(line, "YX")!=0);
  bool lonlat=(strstr(line, "YX")!=0);
  
  if(not latlon and not lonlat) {
    rewind(in);
    printf("no header XY or YX found, assume lat, lon file\n");
    latlon=true;
    rewind(in);
    }

  while (nitems!=0) {
    nitems=fscanf(in,"%lf %lf",&xx,&yy);
    if(nitems!=2) break;
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
  
  if(latlon) {
    for(size_t n=0;n<xtmp.size();n++) {
      (*y)[n]=xtmp[n];
      (*x)[n]=ytmp[n];
      }
    }
  else {
    for(size_t n=0;n<xtmp.size();n++) {
      (*x)[n]=xtmp[n];
      (*y)[n]=ytmp[n];
      }
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

  int parse_GridOptions(const string options, char* & gridfile, char* & maskfile, metagrid_t & meta)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  vector<string> tokens,keys,values;
  int k;
  string delimiter=" ";
  
  gridfile=0;
  
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
      meta.gridfile=strdup(values[k].c_str());
      continue;
      }
    if(keys[k]=="maskfile") {
      meta.maskfile=strdup(values[k].c_str());
      continue;
      }
    if(keys[k]=="topofile") {
      meta.topofile=values[k];
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
    printf("\n%s: keyword <%s> not reckognized \n\n", __func__, keys[k].c_str());
    TRAP_ERR_EXIT(-1,"exiting\n");
    }
  
  if(meta.maskfile=="") meta.maskfile=meta.gridfile.c_str();
  
  maskfile=strdup(meta.maskfile.c_str());
  gridfile=strdup(meta.gridfile.c_str());

//   if(meta.topofile=="") meta.topofile=gridfile;
    
  return(0);
  
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void print_help(const char *prog_name)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/** \brief prints help of this programme.

\param prog_name name of the program to be printed by this help function : argv[0]
*/
/*----------------------------------------------------------------------------*/
{
/** The code of the body of this function is :
    \code /**/ // COMPILED CODE BELOW !!!
  printf("\n"
    "NAME AND VERSION\n"
    "  %s version " VERSION " " REVISION "\n", prog_name);
  printf("\n"
    "USE\n"
    "  %s OPTIONS wave1 [ wave2 ... ]\n",prog_name);
  printf("\n"
    "DESCRIPTION\n"
/* USE FIND AND REPLACE TO EDIT THE LINE BELOW, SO THAT THE FILE HEAD IS ALSO UPDATED */
    "  Creates an open boundary file for TUGOm.\n"
    "\n"
    "OPTIONS :\n"
    "  --help,-h : Show this help and exit.\n"
    "  --quadrangle : Use quadrangles.\n"
    "  --zero : Set 0 depth as masked.\n"
    "  --append : Append points to existing .obc\n"
    "  -m : followed by path to mesh file in .nei format\n"
    "  -b : followed by path to boundary element file in .bel format\n"
    "  -p : followed by path to atlas files directory. Default: ./\n"
    "  -c : followed by atlas file name convention.\n"
    "      There can be up to 3 uses of this option. All extra uses are ignored.\n"
    "      See below for convention.\n"
    "  -o : followed by path to created open boundary file. Default: tides.obc\n"
    "  -v : followed by an empty-terminated list (empty character \"\") of up to 6 variable names. Default: Ha Hg Ua Ug Va Vg\n"
    "  -g : followed by a string such as `file=<value of rootname> <input keys of mesh_meta>' from the .intg for TUGOm\n"
    "\n"
    "TIPS\n"
    "  To get all available atlases :\n"
    "dir=...;ext=.FES2014.nc;f=($(cd $dir;echo *$ext));%s ... -p $dir -c WAVE$ext ${f[@]/$ext}\n",prog_name);
  
  print_tide_decode_atlasname_help();
  /** \endcode */
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double *lon=NULL,*lat=NULL;
  int count,status;
  int status_u,status_v;
  int i,k,n;
  FILE *out=NULL;
  char *keyword=NULL;
  char *meshfile=NULL,*belfile=NULL,*output=NULL,*path=NULL,*atlas_convention[3]={0,0,0};
  char *atlas_directory=NULL,*atlas_file=NULL;
  char filename[1024];
  char *wave[1024];
  float **Ha=NULL,**Hg=NULL,**Ua=NULL,**Ug=NULL,**Va=NULL,**Vg=NULL,mask=NC_FILL_DOUBLE;
  float a[3],G[3];
  int nwave=0;
  spectrum_t spectrum;
  string grid_options="";
  int element=FE_TRIANGLE;
  double tag=NAN;
  
  const int nvars=6;
  char *varnames[nvars]={0,0,0,0,0,0};
  string ampName,phaName;
  bool append=false;

  spectrum_t          imported_spectrum;
  vector<tidaldata_t> imported_TidalData;
  
  spectrum_t reference=initialize_tide();

  if(argc<=1){
    print_help(argv[0]);
    exit(0);
    }
  
  fct_echo( argc, argv);

  n=1;
  while (n < argc) {
    keyword=strdup(argv[n]);

    if(strcmp(argv[n],"--quadrangle")==0){
      element=FE_QUADRANGLE;
      n++;
      continue;
      }
    if(strcmp(argv[n],"--zero")==0){
      tag=0;
      n++;
      continue;
      }
    if(strcmp(argv[n],"-g")==0){
      grid_options=strdup(argv[n+1]);
      n++;
      n++;
      continue;
      }
    if(strcmp(argv[n],"--append")==0){
      append=true;
      n++;
      continue;
      }
    if(strcmp(argv[n],"--help")==0){
      print_help(argv[0]);
      exit(0);
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

/* *----------------------------------------------------------------------
        naming convention for tidal atlases*/
        case 'c' :
          if(atlas_convention[0]==0) {
            atlas_convention[0]= strdup(argv[n+1]);
            }
          else if(atlas_convention[1]==0) {
            atlas_convention[1]= strdup(argv[n+1]);
            }
          else if(atlas_convention[2]==0) {
            atlas_convention[2]= strdup(argv[n+1]);
            }
          n++;
          n++;
          break;

        case 'o' :
          output= strdup(argv[n+1]);
          n++;
          n++;
          break;

/* *----------------------------------------------------------------------
        path for tidal atlases*/
        case 'p' :
          path= strdup(argv[n+1]);
          n++;
          n++;
          break;

/* *----------------------------------------------------------------------
        variable names */
        case 'v' :
          if(varnames[0]){
            __FILE_LINE__(stdout,"multiple use of -v : the last one overrides the previous.");
            for(k=0;k<nvars;k++){
              deletep((void**)(&varnames[k]),free);
              }
            }
/* ----------------------------------------------------------------------------
          empty-terminated list of up to 6 variable names */
          for(k=0;k<nvars;k++){
            if(argv[n+1+k][0]==0){
              k++;
              break;
              }
            varnames[k]=strdup(argv[n+1+k]);
            }
          n++;
          n+=k;
          break;

        case 'h' :
          print_help(argv[0]);
          exit(0);

        default:
          printf("unknown option %s\n",keyword);
          print_help(argv[0]);
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

  if(path!=NULL) atlas_directory=strdup(path);
  
  if(atlas_convention[1]==0) {
    atlas_convention[1]= strdup(atlas_convention[0]);
    }
  if(atlas_convention[2]==0) {
    atlas_convention[2]= strdup(atlas_convention[0]);
    }
  
  printf("number of waves: %d\n",nwave);
  if(nwave==0){
    printf("*** No wave given ***\n");
    print_help(argv[0]);
    wexit(-1);
    }

  exitIfNull(spectrum.waves=new tidal_wave[spectrum.n]);
  
  if (output==0) {
    output=strdup("tides.obc");
    }
  
  if(append) status=load_HamonicOBCs(output, imported_spectrum, imported_TidalData);
  else if(belfile==0){
    printf("*** You must either append points to existing .obc or give boundary element file ***\n");
    print_help(argv[0]);
    wexit(-1);
    }
  
  for (k=0;k<nwave;k++) {
    spectrum.add(reference, wave[k],0);
    }

  if(path ==NULL) path=strdup(".");
  
  if(grid_options!="") {
    char *gridfile=0, *maskfile=0;
    metagrid_t meta;
    grid_t grid;
    mesh_t mesh;
    string rootname="atlas2obc";
    bool debug=true;
    
    status= parse_GridOptions(grid_options, gridfile, maskfile, meta);
    if(status!=0) TRAP_ERR_EXIT(status,"parse_GridOptions(\""+grid_options+"\",\"%s\",\"%s\",) error\n", gridfile, maskfile);
    status=quadrangle_ImportStructured(gridfile, maskfile, mesh, meta, tag, rootname, "", debug);
    if(status!=0) TRAP_ERR_EXIT(status,"quadrangle_ImportStructured(\"%s\",\"%s\",,,%g,\""+rootname+"\",%d) error\n", gridfile, maskfile, tag, debug);
    meshfile=strdup("atlas2obc-quadrangle.nei");
    status=fe_savemesh(meshfile,MESH_FILE_FORMAT_TRIGRID, mesh);
    }
  
  if(imported_TidalData.size()==0) {
    status=get_positions(belfile,meshfile,&lon,&lat,&count);
    if(status!=0) TRAP_ERR_EXIT(status,"get_positions(\"%s\",\"%s\",,,) error\n", belfile,meshfile);
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
 
  for(k=0;k<nvars;k++) {
    if(varnames[k]!=0) continue;
    const char *defVarName;
    switch(k){
      case 0:defVarName="Ha";break;
      case 1:defVarName="Hg";break;
      case 2:defVarName="Ua";break;
      case 3:defVarName="Ug";break;
      case 4:defVarName="Va";break;
      case 5:defVarName="Vg";break;
      default: TRAP_ERR_EXIT(ENOEXEC,"coding error\n");
      }
    varnames[k]=strdup(defVarName);
    }
 
  for (k=0;k<nwave;k++) {
/* *-----------------------------------------------------------------------------
    read tidal atlas database */
    status=tide_decode_atlasname(atlas_directory, atlas_convention[0], wave[k], 0, &atlas_file);
    if(status!=0) TRAP_ERR_EXIT(status,"tide_decode_atlasname(\"%s\",\"%s\",\"%s\",,) error\n", atlas_directory, atlas_convention[0], wave[k]);
    printf("treating %s wave, from %s\n",wave[k],atlas_file);
    
    atlas_grid_or_mesh_t gm;
    int verbose=0;
    const char *wave_null=NULL;
    int strict=0;
    
    ampName=replace(varnames[0],"WAVE",wave[k]);
    phaName=replace(varnames[1],"WAVE",wave[k]);
    status=tide_atlas2positions(atlas_file, ampName.c_str(), phaName.c_str(), lon, lat, count, Ha[k], Hg[k], mask, &gm, verbose, wave_null, strict);
    if(status!=0) TRAP_ERR_EXIT(status,"tide_atlas2positions(\"%s\",\"%s\",\"%s\",...) error\n", atlas_file,varnames[0],varnames[1]);
    
    status=tide_decode_atlasname(atlas_directory, atlas_convention[1], wave[k], 0, &atlas_file);
    if(status!=0) TRAP_ERR_EXIT(status,"tide_decode_atlasname(\"%s\",\"%s\",\"%s\",,) error\n", atlas_directory, atlas_convention[1], wave[k]);
    ampName=replace(varnames[2],"WAVE",wave[k]);
    phaName=replace(varnames[3],"WAVE",wave[k]);
    status_u=tide_atlas2positions(atlas_file,ampName.c_str(),phaName.c_str(),lon,lat,count,Ua[k],Ug[k],mask, &gm, verbose, wave_null, strict);
    
    status=tide_decode_atlasname(atlas_directory, atlas_convention[2], wave[k], 0, &atlas_file);
    if(status!=0) TRAP_ERR_EXIT(status,"tide_decode_atlasname(\"%s\",\"%s\",\"%s\",,) error\n", atlas_directory, atlas_convention[2], wave[k]);
    ampName=replace(varnames[4],"WAVE",wave[k]);
    phaName=replace(varnames[5],"WAVE",wave[k]);
    status_v=tide_atlas2positions(atlas_file,ampName.c_str(),phaName.c_str(),lon,lat,count,Va[k],Vg[k],mask, &gm, verbose, wave_null, strict);
/* *-----------------------------------------------------------------------------
    tidal heights and currents at open limits*/
    sprintf(filename,"%s.obc",wave[k]);
    out=fopen(filename,"w");
    if(out == NULL) TRAP_ERR_EXIT(errno,"Could not open %s (%s)\n",filename,strerror(errno));
//     fprintf(out,"%s\n",wave[k]);
    fprintf(out,"%s %d M\n",wave[k], count);
    for(i=0; i<count; i++) {
      fprintf(out,"%12.6f %12.6f",lon[i],lat[i]);
      if(Hg[k][i]<0.0) Hg[k][i]+=360.0;
      a[0]=Ha[k][i];
      G[0]=Hg[k][i];
      
/*------------------------------------------------------------------------------
      elevation*/
      if(a[0]!=mask) {
        fprintf(out," %9.4f %9.2f",a[0],G[0]);
        }
      if((status_u==0) and (status_v==0)) {
/*------------------------------------------------------------------------------
        barotropic currents*/
        if(Ug[k][i]<0.0) Ug[k][i]+=360.0;
        if(Vg[k][i]<0.0) Vg[k][i]+=360.0;
        if(a[0]!=mask) fprintf(out," %9.4f %9.2f %9.4f %9.2f",Ua[k][i],Ug[k][i],Va[k][i],Vg[k][i]);
        }
      fprintf(out,"\n");
      }
    fclose(out);
    }

/* *-----------------------------------------------------------------------------
  archive T-UGOm expected file*/
  out=fopen(output,"w");
  if (out ==NULL) {
    TRAP_ERR_EXIT(-1,"cannot open obc file : %s \n",output);
    }
  
  fprintf(out,"%d %d %s\n",count,nwave+imported_spectrum.n,"M");
  for (k=0;k<imported_spectrum.n;k++) {
/* *-----------------------------------------------------------------------------
    tidal heights and currents at open limits*/
    fprintf(out,"%s\n",imported_spectrum.waves[k].name);
    for(i=0; i<count; i++) {
/*------------------------------------------------------------------------------
      elevation*/
      fprintf(out,"%12.6f %12.6f",lon[i],lat[i]);
      a[0]=imported_TidalData[i].Ha[k];
      G[0]=imported_TidalData[i].Hg[k];
      if(a[0]!=mask) {
        fprintf(out," %9.4f %9.2f",a[0],G[0]);
        }
      if((status_u==0) and (status_v==0)) {
        if(a[0]!=mask) fprintf(out," %9.4f %9.2f %9.4f %9.2f",imported_TidalData[i].Ua[k],imported_TidalData[i].Ug[k],
                                   imported_TidalData[i].Va[k],imported_TidalData[i].Vg[k]);
        }
      fprintf(out,"\n");
      }
    }
  
  for (k=0;k<nwave;k++) {
/* *-----------------------------------------------------------------------------
    tidal heights and currents at open limits*/
    fprintf(out,"%s\n",wave[k]);
    for(i=0; i<count; i++) {
/*------------------------------------------------------------------------------
      elevation*/
      fprintf(out,"%12.6f %12.6f",lon[i],lat[i]);
      a[0]=Ha[k][i];
      G[0]=Hg[k][i];
      if(a[0]!=mask) {
        fprintf(out," %9.4f %9.2f",a[0],G[0]);
        }
      if((status_u==0) and (status_v==0)) {
        if(a[0]!=mask) fprintf(out," %9.4f %9.2f %9.4f %9.2f",Ua[k][i],Ug[k][i],Va[k][i],Vg[k][i]);
        }
      fprintf(out,"\n");
      }
    }
  fclose(out);

  STDOUT_BASE_LINE("end of %s ... \n",argv[0]);
  exit(0);
}
