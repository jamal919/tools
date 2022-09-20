/** \file
\brief declaration of topo functions
\date reviewed 1 Aug 2011
\author  Damien Allain

\todo library : compile with topo01.cpp
*/

#ifndef __TOPO_PROTOTYPES
#define __TOPO_PROTOTYPES

#include "tools-structures.h"
#include "polygones.h"

extern size_t topo_PolygonsSelection(grid_t grid, plg_t *polygones, int npolygones, signed char *selected, bool spherical=true, int initialise=1);
extern size_t topo_PolygonsSelection(grid_t grid, const char *poly, signed char *selected, bool spherical=true, int initialise=1);
extern size_t topo_PolygonsSelection(grid_t grid, vector<plg_t> & polygons, signed char *selected, bool spherical, int initialise);

extern int topo_RangeSelection(const grid_t & grid, float *topo, float mask, range_t<float> range, signed char *selected, int initialise);

extern int topo_import(char *bathymetry, const char *varname, const grid_t & topogrid, float *topo, float topomask, float zmin, float zmax, char *poly, int masked_only, bool debug);
extern int topo_import(const char *bathymetry, const char *varname, const grid_t & topogrid, float *topo, float topomask, float zmin, float zmax, char *poly, int masked_only, short *tag, short value, const char *input_format, const char *proj4, bool debug);

extern int topo_import(char *bathymetry, const char *varname, const grid_t & topogrid, float *topo, float topomask);

extern int topo_operation(const char *bathymetry, int signus, grid_t topogrid, float *topo, float topomask, float zmin, float zmax, char *poly, int masked_only, int positive_only, int persistence);
extern int topo_set(const grid_t & topogrid, float *topo, float topomask, float zmin, float zmax, const char *poly, int masked_only, int valid_only, int exclusion, float value, short *tag, short tagvalue, bool debug);
extern int topo_shift(const grid_t & topogrid, float *topo, float topomask, float zmin, float zmax, const char *poly, int masked_only, float value, short *tag, short tagvalue, bool debug);

extern int topo_checks(const char *rootname, grid_t topogrid, float * topo, float topomask);

extern int topo_loadfield_cdf(const char *input,const char *varname, grid_t *grid, signed char **buffer, signed char *mask, bool debug);
extern int topo_loadfield_cdf(const char *input,const char *varname, grid_t *grid, short **buffer, short *mask, bool debug);
extern int topo_loadfield_cdf(const char *input,const char *varname, grid_t *grid, float **buffer, float *mask, bool debug);

extern int topo_loadfield(const char *input, const char *varname, grid_t *grid, float **buffer, float *mask, bool debug);
extern int topo_loadfield(const char *input, grid_t *grid, float **buffer, float *mask, bool debug);
extern int topo_loadfield(const char *filename, const char *varname, grid_t *grid, short **buffer, short *mask, bool debug);
extern int topo_loadfield(const char *input, grid_t *grid, short **buffer, short *mask, bool debug);

extern int topo_savefield(const char *filename,grid_t grid, float *buffer, float mask);
extern int topo_save(const char *rootname, char **filename, const char *format, grid_t topogrid, float *ftopo, float ftopomask, bool debug);

extern int topo_save(const string & output, const char *format, grid_t topogrid, float *ftopo, float ftopomask, bool debug);
extern int topo_save(string rootname, string & filename, const char *format, grid_t topogrid, float *ftopo, float ftopomask, bool debug);


extern int topo_hydro(char *hydro, int signus, grid_t topogrid, float *topo, float topomask, int persistence);

extern int loess_filter_init(const grid_t & grid, float scale, float* & weight, int & nx, int & ny);
extern int loess_filter(const grid_t & grid, float *buf, float mask, float *weight, int nx, int ny, int k, int l, float & out);
extern int loess_filter(const grid_t & grid, float *buf, float mask, float scale, int k, int l, float & out);

extern int topo_smooth_spherical(grid_t grid, float *topo, float topomask, float maxscale, float zmin, float zmax, const char *poly, const char *trusted);

extern int topo_smooth_cartesian(const grid_t & grid, float *topo, float topomask, float scale, range_t<float> range, bool keep_masked, const char *poly, const char *trusted, int filter);
extern int topo_smooth_cartesian(const grid_t & grid, float *topo, float topomask, float scale, range_t<float> range, bool keep_masked, vector<plg_t> & polygons, const char *trusted, int filter);


#endif
