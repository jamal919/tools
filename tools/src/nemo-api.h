
#if NEMO_API_H == 0
#define NEMO_API_H 1

#define TEM_ID    0
#define SAL_ID    1
#define RHO_ID    2
#define U_ID      3
#define V_ID      4
#define W_ID      5
#define SSH_ID    6
#define G_ID      7
#define M_ID      8

#include "constants.h"

#ifdef MAIN_SOURCE
#define GLOBAL
#define INITIALVALUE(x) =x
#else
#define GLOBAL extern
#define INITIALVALUE(x)
#endif

GLOBAL double ctrlLon INITIALVALUE(NAN);
GLOBAL double ctrlLat INITIALVALUE(NAN);
GLOBAL int ctrlI INITIALVALUE(-9);
GLOBAL int ctrlJ INITIALVALUE(-9);
GLOBAL int64_t ctrlM INITIALVALUE(-1);
GLOBAL FILE *ctrlF INITIALVALUE(0);

#undef GLOBAL
#undef INITIALVALUE

extern int SG_Wmodes_v1(const grid_t & grid, float *R, float mask, grid_t & topogrid, float *topo, float topomask, int i, int j, float *Nbar, float **celerity, float **modes, int nmodes, double dz);
extern int SG_Wmodes_v2(const grid_t & grid, bool ascending, double factor, float *T, float *S, float *R, float mask, grid_t & topogrid, float *topo, float topomask, int i, int j, float *Nbar, float **celerity, float **modes, int *nmodes, int nkeep,bool rigid_lid);

extern  int BruntVassala(grid_t & grid, float *R, float mask, grid_t & topogrid, float *topo, float topomask, int i, int j, float *N, double *z);
extern  int BruntVassala(grid_t & grid,grid_t & w_grid, float *T, float *S, float *R, float fmask, grid_t & topogrid, float *topo, float topomask, int i, int j, float *N, double *z);

extern int compute_Wmode_v2(int nlevels,double *z, double *rho, double g, double *N, double *c, double **modes);

extern double water_density(double t, double s, double z, int option=2);
extern double water_density_check(int option);

extern int map_loadfield3D(const char *datafile, const char *gridfile, const char *var, grid_t & grid, float * &buffer, float & mask, int verbose);
extern int map_loadfield3D(const char *datafile, const char *gridfile, const char *var, int frame, grid_t & grid, float * &buffer, float & mask, int verbose);


extern int save_SGXYZT_C(const char *output, grid_t grid, complex<double> *buffer, complex<double> mask,
                         const char *varnameA, const char *varnameG, const char *unit, const char *standard, int createfile, int creategrid, const char * vlabel, int frame);

extern int save_SGXYT_C(const char *output, grid_t grid, complex<double> **buffer, complex<double> mask,
                        const char *varnameA, const char *varnameG, const char *unit, const char *standard, int createfile, int creategrid, int range[2]);

extern int save_SGXYZ_C(const char *output, grid_t & grid, complex<double> *buffer, complex<double> mask,
                         const char *varnameA, const char *varnameG, const char *unit, const char *standard, int createfile, int creategrid, const char *hlabel, const char *vlabel);

extern int save_SGXY_C(const char *output, grid_t & grid, complex<double> *buffer, complex<double> mask,
                         const char *varnameA, const char *varnameG, const char *unit, const char *standard, int createfile, int creategrid);

extern int save_SGXYZ(const char *output, grid_t grid, float *buffer, float mask, const char *varname,const char *unit, const char *standard,
                      int createfile, int creategrid, const char *hlabel, const char *vlabel);
extern int save_SGXYZ(const char *output, grid_t grid, double *buffer, double mask, const char *varname,const char *unit, const char *standard,
                      int createfile, int creategrid, const char *hlabel, const char *vlabel);

extern int save_SGXYT(const char *output, grid_t & grid, float *buffer,  float mask,  const char *varname,const char *unit, const char *standard,  int createfile, int creategrid, int frame);
extern int save_SGXYT(const char *output, grid_t & grid, double *buffer, double mask, const char *varname,const char *unit, const char *standard,  int createfile, int creategrid, int frame);

extern int save_SGXYT(const char *output, grid_t & grid, float **buffer,  float mask,  const char *varname,const char *unit, const char *standard,  int createfile, int creategrid, int range[2]);
extern int save_SGXYT(const char *output, grid_t & grid, double **buffer, double mask, const char *varname,const char *unit, const char *standard,  int createfile, int creategrid, int range[2]);

extern  int save_SGXYZT(const char *output, grid_t grid, float **buffer, float mask, const char *varname,const char *unit, const char *standard,
                        int createfile, int creategrid, const char *hlabel, const char *vlabel, int range[2]);
extern  int save_SGXYZT(const char *output, grid_t grid, float *buffer, float mask, const char *varname, const char *unit, const char *standard,
                        int createfile, int creategrid, const char *hlabel, const char *vlabel, int frame);

extern  int VerticalModes_SpectralDecomposition(const string & rootname, int maxmodes, string & wave,const string varnames[5], int nrecomposed, bool debug);
extern  int VerticalModes_SpectralDecomposition(const string & filename,const string & varname, int maxmodes,
                                          double ** &modes, double & mask, complex<double> **p, complex<double> umask, complex<double> ** &p_decomposition, int & MinUmodes, bool debug);

extern  int compute_vertical_mode(const char *bathymetry, int nRequestedProcs, grid_t & t_grid, grid_t & w_grid, float *T, float *S, float *R, float mask,const char *rootname, int maxmodes);


extern int VerticalModes_SpectralDecomposition_NEMO(const string & rootname, const string maskFile, int maxmodes, string & wave, string filenames[4], string varnames[4], int nrecomposed, bool debug, bool recomposition);
extern int VerticalModes_SpectralEnergy(const string maskFile, string ModesFile, int maxmodes, string & wave, int nrecomposed, int verbose, bool debug);
extern int VerticalModes_SpectralDiagnostics(const string maskFile, string ModesFile, int maxmodes, string & wave, int nrecomposed, int verbose, bool debug);

extern int NEMO(const string *input, const string maskFile, const string *grid_variables, const string *variables,const char *rootname,const char *bathymetry, int nRequestedProcs, const char *units, int maxmodes, int verbose);
extern int NEMO_FV(const string *input,const string maskFile,const string *variables,const char *rootname,const char *bathymetry, int nRequestedProcs, const char *units, int maxmodes, int verbose);
extern int ORCA_4(const string *input,const char *rootname,const char *bathymetry, int nRequestedProcs, const char *units, int maxmodes, int verbose);
extern  int SpectralEnergy_RAW(const string maskFile, grid_t & grid, complex<double> **u, complex<double> **v, complex<double> **w, complex<double> **p, string & wave, int verbose, bool debug);

#endif
