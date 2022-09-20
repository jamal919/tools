#ifndef ADMITTANCE_MTS_H
#define ADMITTANCE_MTS_H 1

#define nspecies 3
#define ADMITTANCE_UNAVAILABLE 0
#define ADMITTANCE_SPLINE 1
#define ADMITTANCE_LINEAR 2

extern void admittance_scoefficient(tidal_wave wave[3],double ri[3][3]);
extern void admittance_lcoefficient(tidal_wave wave[2],double ri[2][2]);
extern void admittance_sweight(const tidal_wave &w,const tidal_wave wave[3],const double ri[3][3],double coef[3]);
extern void admittance_lweight(const tidal_wave &w,const tidal_wave wave[2],const double ri[2][2],double coef[2]);
extern void admittance_scompute(const tidal_wave &w,const tidal_wave wave[3],const complex<float> elevation[3],const double ri[3][3],float *amp,float *pha);
extern void admittance_lcompute(const tidal_wave &w,const tidal_wave wave[2],const complex<float> elevation[2],const double ri[2][2],float *amp,float *pha);

#include "tides.h"

class admittance_t {
private :
public :

  double response_s[nspecies][3][3];
  tidal_wave basewave_s[nspecies][3];
  int windex_s[nspecies][3];
  int count_s[nspecies];

  double response_l[nspecies][2][2];
  tidal_wave basewave_l[nspecies][2];
  int windex_l[nspecies][2];
  int count_l[nspecies];

  int method[nspecies];
  
  admittance_t() {
    int i,j;
    for(j=0;j<nspecies;j++){
      for(i=0;i<3;i++){
        this->windex_s[j][i]=-1;
        this->basewave_s[j][i].reset();
        }
      for(i=0;i<2;i++){
        this->windex_l[j][i]=-1;
        this->basewave_l[j][i].reset();
        }
      }
    }
  
};


extern int  admittance_mts_init (admittance_t *admittance, int option=0);
extern int *admittance_mts_check(admittance_t *admittance,const spectrum_t &WaveList, int *available=NULL);

extern void admittance_mts_compute(const admittance_t &admittance,spectrum_t WaveList, hconstant_t ztide,int *available, float mask, bool force, int verbose=1);
extern void admittance_mts_terminate(admittance_t);

extern void admittance_mts_sweightP(const admittance_t &admittance,const tidal_wave &w,double coef[3]);
extern void admittance_mts_lweightP(const admittance_t &admittance,const tidal_wave &w,double coef[2]);

extern int  admittance_mts_verify(const tidal_wave &w, const int *method);
extern int  admittance_mts_verify(const admittance_t &admittance,const spectrum_t & WaveList,int *available,int verbose);

extern int  admittance_mts_error(const admittance_t &admittance,const spectrum_t&,const int *,vector<double>&);

#endif