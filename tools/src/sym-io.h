#if SYM_IO == 0
#define SYM_IO 1

#include "geo.h"

class ntbk_grid_t
{
  /** see http://sirocco.omp.obs-mip.fr/outils/Symphonie/Documentation/SymphonieDocNotebook.htm#grid */
 private:

 public:
  double dXb,dYb;                 /* cartesian steps */
  double dLon,dLat;               /* spherical steps: both unused here */
  double Phi_0;                   /* mercator central latitude */
  double Longi_0,Latit_0;         /* point at i=I0,j=J0 */
  double Angle_0;                 /* rotation (deg) */
  double I0,J0;                   /* see above */
  double RayonTerre;              /* Earth Radius (m) */
  int    TypeGrid,I1D;            /* see load_notebook(), unused here */
  int    MECO, NECO, NR;          /* ni,nj,nk */
  double Pole_LON,Pole_LAT;       /* polar grid pole projection rotation */

  ntbk_grid_t(){
    dXb=dYb=0.0;
    dLon=dLat=0.0;
    Phi_0=0.0;
    Longi_0=Latit_0=0.0;
    Angle_0=0.0;
    I0=J0=-1;
    RayonTerre=6370949. ;
    TypeGrid=I1D=1;
    MECO=NECO=NR=0;
    Pole_LON=0.;
    Pole_LAT=90.;
    }
  ntbk_grid_t &operator = (const ntbk_grid_t &source) {
    return(*this);
    }

};

extern int symio_createfile(char *filename,size_t x_len, size_t y_len, const grid_t & grid);
extern int symio_writefile(char *filename, float *, float *, float *, float *);

extern int savenodes(char *filename, serie_t metadata);
extern grid_t map_get_spherical(geo_t projection,grid_t cgrid);

extern int read_notebook(const char *input, ntbk_grid_t *notebook);
//extern int load_notebook_obsolete(char *input, grid_t *cgrid, geo_t *projection);
extern int load_notebook(const char *input, grid_t *cgrid, grid_t *sgrid, geo_t *projection=NULL);
extern int notebook2grids(const ntbk_grid_t & notebook, grid_t *cgrid, grid_t *sgrid, bool extend=true, geo_t *projection=NULL);

extern int save_notebook(const char *input, const ntbk_grid_t & notebook);

extern int symphonie_loadmask(char *,grid_t,float *,float *);
extern int symphonie_savemask(char *,grid_t,float *,float *);

#endif
