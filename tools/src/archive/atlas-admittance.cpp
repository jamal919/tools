 #define MAIN_SOURCE
/**************************************************************************

  T-UGO tools, 2006-2010

  Unstructured Ocean Grid initiative

Contributors:

  Florent Lyard      LEGOS/CNRS, Toulouse, France
  Yoann Le Bars      LEGOS, Toulouse, France (PhD)
  Laurent Roblou     LEGOS/CNRS, Toulouse, France
  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
  Frédéric Dupont    Université de Laval à Québec, Canada

E-mail: florent.lyard@legos.obs-mip.fr

***************************************************************************/

#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

#include "tools-structures.h"

#include "tides.h"
#include "admittance.h"
#include "tides.def"
#include "functions.h"
#include "bmg.h"
#include "map.h"
#include "fe.h"
#include "archive.h"

#include "netcdf-proto.h"
#include "poc-netcdf-namespace.h"

#ifndef VERBOSE
#define VERBOSE 1
#endif

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int nc_loadc1(std::string fileName, grid_t *grid, fcomplex **cbuffer, fcomplex *cmask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  using namespace pocnetcdf;
/*-----------------------------------------------------------------------
  load netcdf grid */

  int status = -1;
  pocnetcdf::cdfgbl_t global;cdf_globalinfo(fileName.c_str(),&global,0);//pocnetcdf::cdfgbl_t global(fileName);

  int variable = cdf_identify(global, "Ha");
  if( (status = cdf_loadvargrid_2d((char *) fileName.c_str(), variable, grid)) !=0) {
    return(1);
  }

/*-----------------------------------------------------------------------
  load netcdf variable */

  variable = cdf_identify(global, "Ha");

  variable_t varinfo;
  float *abuf = new float[grid->nx * grid->ny];
  float rmask[2] = {static_cast<float>(1.e+35), static_cast<float>(1.e+35)};
  if( (status = cdf_loadvar_r1_2d ((char *) fileName.c_str(), variable, 0, 0, *grid, grid->nx, abuf, &(rmask[0]), &varinfo)) != NC_NOERR) {
    varinfo.reset();
    global.reset();
    return(1);
    }

  varinfo.reset();
  variable = cdf_identify(global, "Hg");

  float *gbuf = new float[grid->nx * grid->ny];
  if( (status = cdf_loadvar_r1_2d ((char *) fileName.c_str(), variable, 0, 0, *grid, grid->nx, gbuf, &(rmask[1]), &varinfo)) != NC_NOERR) {
    varinfo.reset();
    global.reset();
    return(1);
    }

  /*-----------------------------------------------------------------------
  build c1 buffer */

  *cmask = fcomplex(rmask[0], rmask[0]);
  *cbuffer = new fcomplex[grid->nx * grid->ny];
  for(size_t j = 0; j < grid->ny; j++){
    for(size_t i = 0; i < grid->nx; i++){
      size_t n = i + grid->nx * j;
      if( (abuf[n] != rmask[0]) && (gbuf[n] != rmask[1]) ) {
        (*cbuffer)[n] = fcomplex(abuf[n] * cos(gbuf[n] *d2r),
                             -abuf[n] * sin(gbuf[n] *d2r));
        }
      else {
        (*cbuffer)[n] = *cmask;
        }
      }
    }

  delete [] abuf;
  delete [] gbuf;

  varinfo.reset();
  global.reset();

  //clean grid3d TO DO

  return(0);

}// end nc_loadc1()


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int nc_savec1(string fileName, const grid_t grid, const fcomplex *cbuffer, const fcomplex cmask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  using namespace pocnetcdf;

  // convert complex field in amplitude and phase lag fields

  float *abuf = new float[grid.nx * grid.ny];
  float *gbuf = new float[grid.nx * grid.ny];

  for(size_t i = 0; i < grid.nx * grid.ny; i++){
    if( cbuffer[i] != cmask ){
      float x = real(cbuffer[i]);
      float y = imag(cbuffer[i]);
      abuf[i] = sqrt(x * x + y * y);
      gbuf[i] = atan2(-y, x) *r2d;
      }
    else {
      abuf[i] = real(cmask);
      gbuf[i] = real(cmask);
      }
    }

  // write separate fields
  (void) poc_createfile((char *) fileName.c_str());

  pocgrd_t z_ncgrid;
  (void) poc_sphericalgrid_xy((char *) fileName.c_str(), (char *) "", grid, &z_ncgrid);

  pocnetcdf::cdfvar_t amplitude((char *) "Ha",
                                real(cmask),
                                (char *) "m",
                                1.,
                                0.,
                                (char *) "tidal_elevation_amplitude",
                                (char *) "tidal elevation amplitude",
                                (char *) "Ha",
                                z_ncgrid);
  (void) create_ncvariable((char *) fileName.c_str(), &amplitude);

  int status = NC_NOERR + 1;
  if( (status = poc_write_xy((char *) fileName.c_str(), grid, amplitude.id, abuf)) != NC_NOERR){
#if VERBOSE > 0
    cerr << "ERROR poc_write_xy(): unknown format "<<endl;
#endif
    return(1);
    }
  amplitude.reset();

  pocnetcdf::cdfvar_t phase_lag((char *) "Hg",
                                real(cmask),
                                (char *) "degree",
                                1.,
                                0.,
                                (char *) "tidal_elevation_phase_lag",
                                (char *) "tidal elevation phase lag",
                                (char *) "Hg",
                               z_ncgrid);
  (void) create_ncvariable((char *) fileName.c_str(), &phase_lag);

  if( (status = poc_write_xy((char *) fileName.c_str(), grid, phase_lag.id, gbuf)) != NC_NOERR){
#if VERBOSE > 0
    cerr << "ERROR poc_write_xy(): unknown format "<<endl;
#endif
    return(1);
    }
  phase_lag.reset();

// clean z_ncgrid (to do)
  //z_ncgrid.reset();

  delete [] abuf;
  delete [] gbuf;

  return(0);

}// end nc_savec1()


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int s2c_savec1(string fileName, const mesh_t mesh, string msg, const fcomplex *cbuffer, const fcomplex cmask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  ofstream out(fileName.c_str(), ios::out | ios::trunc);

  if(out){
    // write header
    out << msg;
    // convert complex field in amplitude and phase lag fields
    float a = real(cmask);
    float g = real(cmask);

    for(size_t k = 0; k < mesh.nvtxs; k++){
      if( cbuffer[k] != cmask ){
        float x = real(cbuffer[k]);
        float y = imag(cbuffer[k]);
        a = sqrt(x * x + y * y);
        g = atan2(-y, x) *r2d;
        }
      else {
        a = real(cmask);
        g = real(cmask);
        }
      out /*<< setw(7) */<< k + 1 <<
          " " << setprecision (6) << a <<
          " " << setprecision (9) << g << endl;
      }
    out.close();
    return(0);
  }

  else{
    return(1);
    }

}// end s2c_savec1()

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int s2c_savec2(string fileName, const mesh_t mesh, string msg, const fcomplex *u, const fcomplex *v, const fcomplex cmask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  float x,y,a[2],g[2];
  ofstream out(fileName.c_str(), ios::out | ios::trunc);

  if(out){
    // write header
    out << msg;
    // convert complex field in amplitude and phase lag fields

    for(size_t k = 0; k < mesh.nvtxs; k++){
      if(( u[k] != cmask )&&( v[k] != cmask )){
        x = real(u[k]);
        y = imag(u[k]);
        a[0] = sqrt(x * x + y * y);
        g[0] = atan2(-y, x) *r2d;
        x = real(v[k]);
        y = imag(v[k]);
        a[1] = sqrt(x * x + y * y);
        g[1] = atan2(-y, x) *r2d;
        }
      else {
        a[0] = real(cmask);
        g[0] = real(cmask);
        a[1] = real(cmask);
        g[1] = real(cmask);
        }
      out /*<< setw(7) */<< k + 1 <<
          " " << setprecision (6) << a[0] <<
          " " << setprecision (9) << g[0] <<
          " " << setprecision (6) << a[1] <<
          " " << setprecision (9) << g[1] << endl;
      }
    out.close();
    return(0);
  }

  else{
    return(1);
    }

}


// static bool expired=expire(20130408,20130508);


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  string atlasDirectory = ".";
  string atlasConvention = "";
  string meshFile = "";
  string outpath = "";
  string input = "";

  bool __known_format = false;
  
//#if VERBOSE > 1
  fct_echo(argc,argv);
//#endif
  printf("\n ^^^^^^^^^^^^^ start computation ^^^^^^^^^^^^^\n");

  int n = 1;
  while (n < argc) {
    char *keyword = strdup(argv[n]);
    switch (keyword[0]) {
      case '-':
        switch (keyword[1]) {

            /*----------------------------------------------------------------------
            naming convention for tidal atlases*/
          case 'c' :
            atlasConvention = argv[n+1];
            n++;
            n++;
            break;

          /*----------------------------------------------------------------------
            path for tidal atlases*/
          case 'p' :
            atlasDirectory = argv[n+1];
            n++;
            n++;
            break;

            /*----------------------------------------------------------------------
            FE mesh (full) name (if needed)*/
          case 'm' :
            meshFile = argv[n+1];
            n++;
            n++;
            break;

          /*----------------------------------------------------------------------
            path for output tidal map*/
          case 'o' :
            outpath = argv[n+1];
            n++;
            n++;
            break;

        default:
#if VERBOSE > 0
          cerr << "ERROR : unknown option " << keyword << endl;
#endif
          __OUT_BASE_LINE__("\n ^^^^^^^^^^^^^ computation aborted ^^^^^^^^^^^^^\n");
          exit(-1);
        }
        break;

        /*----------------------------------------------------------------------
          input file (NETCDF, BMG or S2C) */
      default:
          input = argv[n];
          n++;
        break;
      }
      free(keyword);
      }

    tidal_wave wave = wave_from_name((char *) input.c_str());
    wave.init();

    // define basis constituents for admittance functions
    spectrum_t basisList;
    basisList.n = basisList.nmax = 3;
    basisList.waves = new tidal_wave[basisList.n];
    switch (wave.nT) {
      case 0:{
        basisList.waves[0] = wMf;
        basisList.waves[1] = wMm;
        basisList.waves[2] = wMtm;
        (void) initialize_omega(&basisList);
        break;
      }
      case 1:{
        basisList.waves[0] = wK1;
        basisList.waves[1] = wO1;
        basisList.waves[2] = wQ1;
        (void) initialize_omega(&basisList);
        break;
        }
      case 2:{
        basisList.waves[0] = wM2;
        basisList.waves[1] = wN2;
        basisList.waves[2] = wK2;        // pb grid K2 (<= FES2002!!)
        (void) initialize_omega(&basisList);
        break;
        }
      default:{
#if VERBOSE > 0
        cerr << "ERROR : admittance is not possible for "<< input << endl;
#endif
        __OUT_BASE_LINE__("\n ^^^^^^^^^^^^^ computation aborted ^^^^^^^^^^^^^\n");
        exit(1);
        }
      }// end switch()

    // default configuration initialization
    size_t pos = string::npos;
    string rootname(atlasConvention);
    string ext = "";
    if( (pos = rootname.find_last_of(".")) != string::npos){
      ext = rootname.substr(pos + 1);
      rootname = rootname.substr(0, pos);
      if( (pos = rootname.find_last_of("/")) != string::npos){
        rootname = rootname.substr(pos);
        }
      if( (pos = rootname.find_first_of("WAVE")) != string::npos){
        rootname = rootname.replace(pos, strlen("WAVE"), input);
        }
      }

    string outFormat = "";
    string extension = "";
    if( ext.compare("nc") == 0 ){
      outFormat = "NETCDF";
      extension = ".nc";
      __known_format = true;
      }

    if( (ext.compare("bmg") == 0) || (ext.compare(0, 3, "den") == 0) ){ // for den.a, den.h, den.something_else
      outFormat = "BMG";
      extension = ".bmg";
      __known_format = true;
      }

    if( ext.compare("s2c") == 0 ){
      outFormat = "S2C";
      extension = ".s2c";
      __known_format = true;
      if(meshFile.empty()){
#if VERBOSE > 0
        cerr << "ERROR : missing FE mesh" << endl;
#endif
        __OUT_BASE_LINE__("\n ^^^^^^^^^^^^^ computation aborted ^^^^^^^^^^^^^\n");
        exit(1);
        }
      }

    if( ext.compare("v2c") == 0 ){
      outFormat = "V2C";
      extension = ".v2c";
      __known_format = true;
      if(meshFile.empty()){
#if VERBOSE > 0
        cerr << "ERROR : missing FE mesh" << endl;
#endif
      __OUT_BASE_LINE__("\n ^^^^^^^^^^^^^ computation aborted ^^^^^^^^^^^^^\n");
      exit(1);
      }
    }

    if(outpath.empty()){
      outpath = strdup("./");
      }
    outpath += '/';

    string outfile = "";
    if( __known_format){
      outfile = outpath + rootname + extension;
      }
    else{
#if VERBOSE > 0
      cerr << "ERROR : unknown atlas format" << endl;
#endif
      __OUT_BASE_LINE__("\n ^^^^^^^^^^^^^ computation aborted ^^^^^^^^^^^^^\n");
      exit(1);
      }

  // read basis constituents

  mesh_t mesh;
  bool mesh_read = false;
  grid_t grid[3];
  fcomplex *cbuffer[3] = {0, 0, 0}, cmask[3] = {0, 0, 0};
  int status = -1;
  char *atlas_file = 0;
  int nvalues=-1;

  for(size_t w = 0; w < basisList.n; w++){
    // build the basis file name
    (void) tide_decode_atlasname((char *) atlasDirectory.c_str(),
            (char *) atlasConvention.c_str(), basisList.waves[w].name, &atlas_file);
#if VERBOSE > 1
    printf("loading %s wave from %s \n",
             basisList.waves[w].name, atlas_file);
#endif

    // read file wrt format
    // z = ( A cos -G (in rad) ; A sin -G (in rad) )
    if(outFormat.compare("BMG") == 0){
      float time = 0;
      (void) bmg_loadc1(atlas_file, 1, 1, 1, &(grid[w]), &(cbuffer[w]), &(cmask[w]), &time);
      (void) map_completegridaxis(&(grid[w]), 2);
      } // end load bmg

    if(outFormat.compare("NETCDF") == 0){
      if( (status = nc_loadc1(atlas_file, &(grid[w]), &(cbuffer[w]), &(cmask[w]))) != 0){
#if VERBOSE > 0
        cerr << "ERROR : NETCDF format not readable" << endl;
#endif
        __OUT_BASE_LINE__("\n ^^^^^^^^^^^^^ computation aborted ^^^^^^^^^^^^^\n");
        exit(1);
        }
      }// end load netcdf

    if(outFormat.compare("S2C") == 0){
      if( ( ! meshFile.empty()) && (! mesh_read) ){
        (void) fe_init(&mesh);
        if( (status = fe_readmesh((char *) meshFile.c_str(), MESH_FILE_FORMAT_TRIGRID, &mesh)) != 0){
#if VERBOSE > 0
          cerr << "ERROR : mesh format is not TRIGRID" << endl;
#endif
          __OUT_BASE_LINE__("\n ^^^^^^^^^^^^^ computation aborted ^^^^^^^^^^^^^\n");
          exit(1);
          }
        if( (status = fe_list(&mesh)) != 0){
#if VERBOSE > 0
          cerr << "ERROR : fe_list() failed"<<endl;
#endif
          __OUT_BASE_LINE__("\n ^^^^^^^^^^^^^ computation aborted ^^^^^^^^^^^^^\n");
          exit(1);
          }
        mesh_read = true;
        }// end load mesh()

      if( (cbuffer[w] = new fcomplex[mesh.nvtxs]) == 0){
#if VERBOSE > 0
        cerr << "ERROR : new[] allocation failed" << endl;
#endif
        __OUT_BASE_LINE__("\n ^^^^^^^^^^^^^ computation aborted ^^^^^^^^^^^^^\n");
        exit(1);
        }
      status=quoddy_loadc1(atlas_file, mesh.nvtxs, cbuffer[w]);
      if( status!= 0){
        __OUT_BASE_LINE__("\n ^^^^^^^^^^^^^ computation aborted ^^^^^^^^^^^^^\n");
        exit(1);
        }
      nvalues=mesh.nvtxs;
      } //end load s2c

    if(outFormat.compare("V2C") == 0){
      if( ( ! meshFile.empty()) && (! mesh_read) ){
        (void) fe_init(&mesh);
        if( (status = fe_readmesh((char *) meshFile.c_str(), MESH_FILE_FORMAT_TRIGRID, &mesh)) != 0){
#if VERBOSE > 0
          cerr << "ERROR : mesh format is not TRIGRID" << endl;
#endif
          __OUT_BASE_LINE__("\n ^^^^^^^^^^^^^ computation aborted ^^^^^^^^^^^^^\n");
          exit(1);
          }
        if( (status = fe_list(&mesh)) != 0){
#if VERBOSE > 0
          cerr << "ERROR : fe_list() failed"<<endl;
#endif
          __OUT_BASE_LINE__("\n ^^^^^^^^^^^^^ computation aborted ^^^^^^^^^^^^^\n");
          exit(1);
          }
        mesh_read = true;
        }// end load mesh()

      if( (cbuffer[w] = new fcomplex[2*mesh.nvtxs]) == 0){
#if VERBOSE > 0
        cerr << "ERROR : new[] allocation failed" << endl;
#endif
        __OUT_BASE_LINE__("\n ^^^^^^^^^^^^^ computation aborted ^^^^^^^^^^^^^\n");
        exit(1);
        }
      status=quoddy_loadc2(atlas_file, mesh.nvtxs, cbuffer[w], &(cbuffer[w][mesh.nvtxs]));
      if( status!= 0){
        __OUT_BASE_LINE__("\n ^^^^^^^^^^^^^ computation aborted ^^^^^^^^^^^^^\n");
        exit(1);
        }
      nvalues=2*mesh.nvtxs;
      } //end load s2c

    // clean
    delete [] atlas_file;
    } // end load basis constituents

  // compute atlas expansion by admittance
  int *method = admittance_check(basisList, wave);
  double coef[3] = {static_cast<double>(0), static_cast<double>(0), static_cast<double>(0)};
  switch(method[wave.nT]){
    case 2:{
      (void) admittance_lweightP(wave, coef);
      for(size_t k = 0; k < 2; k++){
/* *----------------------------------------------------------------------------
        Bug !! : 28/05/2010 */
//        coef[k] *= basisList.waves[k].Ap;
        coef[k] *= wave.Ap;
        }
      break;
      }
    case 1:{
      (void) admittance_sweightP(wave, coef);
      for(size_t k = 0; k < 3; k++){
/* *----------------------------------------------------------------------------
        Bug !! : 28/05/2010*/
//        coef[k] *= basisList.waves[k].Ap;
        coef[k] *= wave.Ap;
        }
      break;
      }
    default:{
#if VERBOSE > 0
      cerr << "ERROR : can not compute admittance coefficients" << endl;
#endif
      __OUT_BASE_LINE__("\n ^^^^^^^^^^^^^ computation aborted ^^^^^^^^^^^^^\n");
      exit(1);
      }
    }
  delete [] basisList.waves;

  fcomplex *outBuffer = 0;
  if((outFormat.compare("S2C") != 0) && (outFormat.compare("V2C") != 0)){ // NETCDF and BMG cases
    if( (outBuffer = new fcomplex[grid[0].nx * grid[0].ny]) == 0){
#if VERBOSE > 0
      cerr << "ERROR : new[] allocation failed" << endl;
#endif
      __OUT_BASE_LINE__("\n ^^^^^^^^^^^^^ computation aborted ^^^^^^^^^^^^^\n");
      exit(1);
      }

    for(size_t j = 0; j < grid[0].ny; j++){
      for(size_t i = 0; i < grid[0].nx; i++){
        size_t n = i + grid[0].nx * j;
        if( (cbuffer[0][n] != cmask[0])
         && (cbuffer[1][n] != cmask[1])
         && (cbuffer[2][n] != cmask[2]) ){
          outBuffer[n].real() = coef[0] * cbuffer[0][n].real()
            + coef[1] * cbuffer[1][n].real()
            + coef[2] * cbuffer[2][n].real();
          outBuffer[n].imag() = coef[0] * cbuffer[0][n].imag()
            + coef[1] * cbuffer[1][n].imag()
            + coef[2] * cbuffer[2][n].imag();
              }
        else{
          outBuffer[n] = cmask[0];
          }
        }
      }// end loop
    }// end NETCDF and BMG cases

  else{ // S2C case...
    if( (outBuffer = new fcomplex[nvalues]) == NULL){
#if VERBOSE > 0
      cerr << "ERROR : new[] allocation failed" << endl;
#endif
      __OUT_BASE_LINE__("\n ^^^^^^^^^^^^^ computation aborted ^^^^^^^^^^^^^\n");
      exit(1);
      }
    else{
      for(size_t n = 0; n < nvalues; n++){
        outBuffer[n].real() = coef[0] * cbuffer[0][n].real()
            + coef[1] * cbuffer[1][n].real()
            + coef[2] * cbuffer[2][n].real();
        outBuffer[n].imag() = coef[0] * cbuffer[0][n].imag()
            + coef[1] * cbuffer[1][n].imag()
            + coef[2] * cbuffer[2][n].imag();
        }
      }
    }// end S2C case

  delete [] cbuffer[0];
  delete [] cbuffer[1];
  delete [] cbuffer[2];

  admittance_terminate();

  // write file
  if(outFormat.compare("BMG") == 0){
    float time = 0;
    (void) bmg_savec1((char*) outfile.c_str(), 1, 1, 1, grid[0], outBuffer, time, cmask[0]);
    }

  if(outFormat.compare("NETCDF") == 0){
    if( (status = nc_savec1(outfile, grid[0], outBuffer, cmask[0])) != 0){
#if VERBOSE > 0
      cerr << "ERROR : nc_savec1() failed" << endl;
#endif
      __OUT_BASE_LINE__("\n ^^^^^^^^^^^^^ computation aborted ^^^^^^^^^^^^^\n");
      exit(1);
      }
    }

  if(outFormat.compare("S2C") == 0){
    ostringstream oss;
    oss << "This file is generated with " << argv[0] <<
        "\nAssociated mesh file is " << meshFile <<
        "\nBasis constituents are take from files in " << atlasDirectory << endl;
    if( (status = s2c_savec1(outfile, mesh, oss.str(), outBuffer, cmask[0])) != 0){
#if VERBOSE > 0
      cerr << "ERROR : s2c_savec1() failed" << endl;
#endif
      __OUT_BASE_LINE__("\n ^^^^^^^^^^^^^ computation aborted ^^^^^^^^^^^^^\n");
      exit(1);
      }
    for (size_t i = 0; i < mesh.nvtxs; i++) {
      delete[] mesh.vertices[i].ngh;
      }
    delete [] mesh.vertices;
    delete [] mesh.triangles;
    }

  if(outFormat.compare("V2C") == 0){
    ostringstream oss;
    oss << "This file is generated with " << argv[0] <<
        "\nAssociated mesh file is " << meshFile <<
        "\nBasis constituents are take from files in " << atlasDirectory << endl;
    if( (status = s2c_savec2(outfile, mesh, oss.str(), outBuffer, &(outBuffer[mesh.nvtxs]), cmask[0])) != 0){
#if VERBOSE > 0
      cerr << "ERROR : s2c_savec1() failed" << endl;
#endif
      __OUT_BASE_LINE__("\n ^^^^^^^^^^^^^ computation aborted ^^^^^^^^^^^^^\n");
      exit(1);
      }
    for (size_t i = 0; i < mesh.nvtxs; i++) {
      delete[] mesh.vertices[i].ngh;
      }
    delete [] mesh.vertices;
    delete [] mesh.triangles;
    }

  delete [] outBuffer;
  atlasConvention.clear();
  atlasDirectory.clear();

  printf("\n ^^^^^^^^^^^^^ end of computation ^^^^^^^^^^^^^\n");
  return(0);

abort:
  __OUT_BASE_LINE__("\n ^^^^^^^^^^^^^ computation aborted ^^^^^^^^^^^^^\n");
  exit(1);

}
