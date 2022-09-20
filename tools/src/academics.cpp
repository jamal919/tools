
/*******************************************************************************
  
  T-UGO tools, 2006-2009

  Unstructured Ocean Grid initiative
  
*******************************************************************************/
/**
\file

<!-- USE FIND AND REPLACE TO EDIT THIS LINE, SO THAT print_help IS ALSO UPDATED -->
\brief 
**/
/*----------------------------------------------------------------------------*/

#define MAIN_SOURCE

#include <unistd.h>

#include <iostream>
#include <vector>

#include "config.h"

#include "tools-define.h"
#include "tools-structures.h"
#include "tides.h"
#include "xtrack-io.h"

#ifndef VERBOSE
#define VERBOSE 1
#endif

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  char option = 'a';
  double sampling = 0, duration = 0, lag = 0;
  double lon = 0, lat = 0;
  vector<double> amplitude(0), phase(0);
  double noise_level = 0;
  bool gaussian_noise = false;
  vector<string> wave(0);
  string str = "";
  string output = "prediction.dat";
  string help =  "\nSYNTAX academics -s <sampling> (days) -d <duration> (days)\n\
     -l <lag> (1st sample in CNES days) -a <[amplitude1,amplitude2,...] | amplitude> (m)\n\
      -g <[phase1, phase2,...] | phase> (degrees)\n\
      -t <lon> (degrees) -p <lat> (degrees) -n <noise> (%)\n\
      -w <[wave1,wave2,...] | wave> [-o <output>]\n\
     [-h]\n\n";
  
  fct_echo(argc,argv);
 
  while ((option = getopt (argc, argv, "hd:s:f:l:a:g:w:o:t:p:n:")) != -1) {
    switch (option) {
    case 'h':
      cerr << help << endl;
      TRAP_ERR_EXIT(1,"exiting\n");
      break;
      
    case 's':
      if (optarg) str = optarg;
      if(!str.empty()) {
        sampling = atof(str.c_str());
        str.clear();
        }
      break;

    case 'd':
      if (optarg) str = optarg;
      if(!str.empty()) {
        duration = atof(str.c_str());
        str.clear();
        }
      break;

    case 'l':{
      if (optarg) str = optarg;
      if(!str.empty()) {
        lag = atof(str.c_str());
        str.clear();
        }
      break;

    case 'a':
      if (optarg) {
        str = optarg;
        size_t pos = string::npos;
        while( (pos = str.find_first_of("[,]")) != string::npos){
          if( (pos == 0) ||  (pos == str.size()) ){
            str = str.substr(pos + 1, string::npos - pos);
            continue;
            }
          else {
            amplitude.push_back( atof(str.substr(0, pos).c_str()) );
            str = str.substr(pos + 1, string::npos - pos);
            pos = string::npos;
            }
          }
        if(amplitude.empty()) {
          amplitude.push_back(atof(str.c_str()));
          }
        str.clear();
      }
      break;

    case 'g':
      if (optarg) {
        str = optarg;
        size_t pos = string::npos;
        while( (pos = str.find_first_of("[,]")) != string::npos){
          if( (pos == 0) ||  (pos == str.size()) ){
            str = str.substr(pos + 1, string::npos - pos);
            continue;
            }
          else {
            phase.push_back( atof(str.substr(0, pos).c_str()) );
            str = str.substr(pos + 1, string::npos - pos);
            pos = string::npos;
            }
          }
        if(phase.empty()) {
          phase.push_back(atof(str.c_str()));
          }
        str.clear();
        }
      break;

    case 't':
      if (optarg) str = optarg;
      if(!str.empty()) {
        lon = atof(str.c_str());
        str.clear();
        }
      break;

    case 'p':
      if (optarg) str = optarg;
      if(!str.empty()) {
        lat = atof(str.c_str());
        str.clear();
        }
      break;
      
    case 'n':
      if (optarg){
        str = optarg;
        noise_level = atof(str.c_str());
        gaussian_noise = true;
        str.clear();
        }
        break;

    case 'w':
    if (optarg){
      str = optarg;
      size_t pos = string::npos;
      while( (pos = str.find_first_of("[,]")) != string::npos){
        if( (pos == 0) ||  (pos == str.size()) ){
          str = str.substr(pos + 1, string::npos - pos);
          continue;
          }
        else{
          wave.push_back(str.substr(0, pos));
          str = str.substr(pos + 1, string::npos - pos);
          pos = string::npos;
          }
        }
      if(wave.empty()) {
        wave.push_back(str.c_str());
        }
      str.clear();
    }
    break;

    case 'o':
      if (optarg) output = optarg;
      break;
    }
  
    case '?':
      cerr << "unknown option" << endl ;
      cerr << help << endl;
      return(1);
      
    }
  }
  
  
#if VERBOSE > 1
  printf("\n\n-------------- beginning computation --------------\n\n\n");
  fflush(stdin);
#endif
  
  int n = (int)(duration / sampling) + 1;
  double *time = new double[n];
  for(size_t t = 0; t < n; ++t){
    time[t] = sampling * t ;
    }
  double *tides = new double[n];

  spectrum_t PredictionList;
  if(wave.empty()) {
    cerr << "syntax error: wave list missing" << endl << help << endl;
    }
  PredictionList.n = PredictionList.nmax = wave.size();
  PredictionList.waves = new tidal_wave[PredictionList.n];
  for(size_t w = 0; w != PredictionList.n; w++){
    PredictionList.waves[w] = wave_from_name((char*) wave[w].c_str());
    }
  
  if(amplitude.empty()){
    cerr << "syntax error: amplitude list missing" << endl << help << endl;
    }
  else {
    if(amplitude.size() < wave.size()){
      amplitude.resize(wave.size(), amplitude[0]);
#if VERBOSE > 0
    cout << "WARNING: missing amplitude values set to " << amplitude[0] << endl;
#endif
      }
    }
  
  if(phase.empty()){
    cerr << "ERROR: amplitude list missing" << endl << help << endl;
    }
  else{
    if(phase.size() < wave.size()){
      phase.resize(wave.size(), phase[0]);
#if VERBOSE > 0
      cout << "WARNING: missing phase lag values set to " << phase[0] << endl;
#endif
      }
    }
 
  hconstant_t seed;
  seed.a = new float[PredictionList.n];
  for(size_t w = 0; w != PredictionList.n; w++){
    seed.a[w] = amplitude[w];
    }
  seed.G = new float[PredictionList.n];
  for(size_t w = 0; w != PredictionList.n; w++){
    seed.G[w] = phase[w];
    }

  FILE *outFile = 0;
  if( (outFile = fopen((char *)output.c_str(),"w")) == 0 ){
#if VERBOSE > 1
    cout << "ERROR : can not create " << output << <<endl;
#endif
    return(1);
    }
  else {
#if VERBOSE > 1
    cout << "results are written in "<< output << endl << endl;
#endif
    }
  int ns = ( (lag > 0 ) ? ((int)(sampling / lag) + 1) : (1));
  XTRACK_ref_write_header(outFile, ns);
 
  std::vector<int> indexRecord(0);
  for(size_t i = 0; i < ns; ++i){
    indexRecord.push_back(i);
    }

  serie_t currentRecord;
  currentRecord.mask = 99.9999;
  currentRecord.count = n;
  currentRecord.data = new data_t[n];
  for(size_t t = 0; t < n; ++t){
    for(size_t c = 1; c < 30; ++c)
      currentRecord.data[t].values[c] = currentRecord.mask;
    }
  for(size_t t = 0; t < n; ++t){
    currentRecord.data[t].values[1]  = 0;
    currentRecord.data[t].values[10] = 0;
    currentRecord.data[t].values[11] = 0;
    currentRecord.data[t].values[15] = 0;
    }
  
  const gsl_rng_type *T = gsl_rng_default;
  gsl_rng *r = gsl_rng_alloc(T);
  gsl_rng_env_setup();
  int nodal = 1;

  for(size_t i = 0; i < ns; ++i){
    for(size_t t = 0; t < n; ++t){
      time[t] = sampling * t + i * lag;
      currentRecord.data[t].time = time[t];
      }
    harmonic_predictionTS(tides, time, n, PredictionList, seed.a, seed.G, nodal);
    currentRecord.data[i].lon = lon;
    currentRecord.data[i].lat = lat;
    if(gaussian_noise == true){
      for(size_t t = 0; t < n; ++t){
        #warning should be from the IFT of a random-phase white spectrum
        currentRecord.data[t].values[0] = tides[t] + gsl_ran_gaussian(r, noise_level);
        }
      }
    else{
      for(size_t t = 0; t < n; ++t){
        currentRecord.data[t].values[0] = tides[t];
        }
      }
    XTRACK_WriteAsciiRecords(outFile, &currentRecord, &(indexRecord[i]), 1, 0);
    }

  delete [] seed.a;
  delete [] seed.G;
  indexRecord.clear();
  delete [] time;
  delete [] tides;
  delete [] currentRecord.data;
  delete [] PredictionList.waves;
  
  fclose(outFile);

  printf(" ^^^^^^^^^^^^^ end of computation ^^^^^^^^^^^^^\n");
  return(0);

}
