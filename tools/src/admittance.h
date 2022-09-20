#ifndef ADMITTANCE_H
#define ADMITTANCE_H 1

#include "admittance-mts.h"

extern  void admittance_test(spectrum_t WaveList,spectrum_t AdmittanceList, int nby);
extern  int admittance_init(const spectrum_t &basisList, const tidal_wave &w);
extern  void admittance_terminate();
extern  int *admittance_check(const spectrum_t &WaveList, int option);
extern  int *admittance_check(const spectrum_t &WaveList, const tidal_wave &w);
extern  void admittance_sweightP(tidal_wave w, double coef[3]);
extern  int admittance_verify(tidal_wave w);
extern  void admittance_lweightP(tidal_wave w, double coef[2]);
extern  void admittance_compute(spectrum_t WaveList, hconstant_t ztide, int *available);
extern  int admittance_error_unsafe(const spectrum_t& s, const int *deduce, vector<double>& error);

#endif