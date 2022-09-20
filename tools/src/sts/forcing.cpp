
/**************************************************************************

  T-UGOm hydrodynamic ocean model, 2006-2012

  Unstructured Ocean Grid initiative

Contributors:

  Florent Lyard      LEGOS/CNRS, Toulouse, France
  Cyril Nguyen       LA/CNRS,    Toulouse, France
  Laurent Roblou     LEGOS/CNRS, Toulouse, France
  Damien Allain      LEGOS/CNRS, Toulouse, France
  
  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
  Frédéric Dupont    Université de Laval à Québec, Canada
  
  Yoann Le Bars      PhD, LEGOS, Toulouse, France
  Yves Soufflet      Post-doctorant, LEGOS, Toulouse, France
  Clement Mayet      PhD, LEGOS, Toulouse, France

***************************************************************************/

#ifdef HAVE_MPI
#include "mpi.h"
#endif

#include "tugo-prototypes.h"
#include "tides.h"
#include "poc-netcdf-data.hpp"
#include "map.h"

#ifdef PERTURB
#   include "perturb.h"
#endif

extern int tide_complex(double, spectrum_t, float *, float *, int);

extern void   meteo_forcing2D(double, mesh_t &, state2d_t *,parameter_t , action_t *);


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int read_loading01(int nndes, parameter_t data, int *found)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------
 Loading and self-attraction from diagnostic fields
 Read the amplitude (meter) and phase lag, convert into
 real and imaginary part 
----------------------------------------------------------------------*/
{
  int i, k, m, n, items, status;
  float zr, zi, a, G;
  float *amp, *pha, mask;
  FILE *in;
  float dummy;
  double x, y, dum;
  char char_dum[1024];
  char *filename;
  size_t count[3], start[3], var;

  grid_t LSAgrid;

  filename = new char[1024];

  for(k = 0; k < WaveList.n; k++) {
    strcpy(filename, LSA_directory);
    filename = strcat(filename, "/");
    filename = strcat(filename, WaveList.waves[k].name);
    filename = strcat(filename, ".load");
    if((in = fopen(filename, "r")) == NULL) {
      continue;
    } else
      printf("opening loading file: %s \n", filename);

    m = 999;
    for(n = 1; n <= 4; n++) {
      fgets(char_dum, m, in);
      }

    for(n = 0; n < nndes; n++) {
      items = fscanf(in, "%f %f %f %f", &zr, &zi, &a, &G);
      if(items != 4)
        check_error(-1,"read error", __LINE__, __FILE__, 1);
      G *= d2r;
      data.LSA[n].z[k]= a*fcomplex(cos(-G),sin(-G));
      }
    found[k] = 1;
    fclose(in);
    }
  zaparr(filename);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int decode_loadingname(char *wave, char **filename)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------
  Atmospheric pressure and wind from AM fields; load the field at or
  just before t
 ----------------------------------------------------------------------*/
{
  int k, i, nt, nframe, status, hour;
  char dummy[256], *pointer;
  char *tmp,*p2;
  FILE *out;

  size_t *start, *count;
  cdfvar_t info;

/*----------------------------------------------------------------------
  build the atlas file name*/
  (*filename) = new char[strlen(LSA_directory) + 256];
  sprintf((*filename), "%s/%s", LSA_directory, LSA_convention);

/*----------------------------------------------------------------------
  check existence name*/
  out = fopen(*filename, "r");
  status = (out == NULL);

  switch (status) {
    case 0:
/*----------------------------------------------------------------------
      file exists, do nothing more*/
      fclose(out);
      break;

    default:
/*----------------------------------------------------------------------
      use format information*/
      pointer = strstr((*filename), "WAVE");
      if(pointer != NULL) {
        tmp=strdup(*filename);
        p2 = strstr(tmp, "WAVE");
        sprintf(dummy, "%s", wave);
        strcpy(pointer, dummy);
        strcpy(pointer+strlen(dummy), p2+4);
        status=0;
        }
      pointer = strstr((*filename), "wave");
      if(pointer != NULL) {
        tmp=strdup(*filename);
        p2 = strstr(tmp, "wave");
        sprintf(dummy, "%s", wave);
        int k;
        for(k=0;k<strlen(dummy);k++) {
          dummy[k]=tolower(dummy[k]);
          }
        strcpy(pointer, dummy);
        strcpy(pointer+strlen(dummy), p2+4);
        status=0;
        }
    }

  return (status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int read_loading02(int nndes, parameter_t data, int *found)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------
 Loading and self-attraction from diagnostic fields
 Read the amplitude (meter) and phase lag, convert into
 real and imaginary part 
----------------------------------------------------------------------*/
{
  int k, m, n, status, verbose=0;
  float factor;
  fcomplex *z,cmask;
  double x, y;
  char *filename=0;
  size_t length;
  grid_t LSAgrid;
  decoded_t decoded;

  for(k = 0; k < WaveList.n; k++) {
    poc_global_t global;
    const poc_var_t *ampVar,*phaVar;
    const poc_att_t *ampUnitAtt;
    string ampUnit;
    
    status=decode_loadingname(WaveList.waves[k].name,&filename);
    if(gCPU_ID==gCPU_MASTER) {
      printf("seeking tidal LSA for %s in %s\n",WaveList.waves[k].name,filename);
      }
    
    status=poc_inq(filename, &global, verbose);
    if(status!=0) {
      if(gCPU_ID==gCPU_MASTER) printf("File not found, assumes null LSA for %s\n",WaveList.waves[k].name);
      continue;
      }
    
/**----------------------------------------------------------------------------
    identigy LSA variables (amplitude and phase lag) */
    ampVar=global.variables.findP("LSAa");
    phaVar=global.variables.findP("LSAg");

/**----------------------------------------------------------------------------
    patches for old formats */
    if(!ampVar) {
      ampVar=global.variables.findP("Ha");
      }
    if(!phaVar) {
      phaVar=global.variables.findP("Hg");
      }

    if(!ampVar) {
      printf("could not identify LSA amplitude variable in file\n");
      break;
      }
    if(!phaVar) {
      printf("could not identify LSA phase lag variable in file\n");
      break;
      }
    
    status = poc_get_grid(filename, *ampVar, &LSAgrid);
    if(status != 0) {
      delete[] filename;
      continue;
      }
    length=LSAgrid.nx * LSAgrid.ny;
    
    if(gCPU_ID==gCPU_MASTER) printf("file %s, amplitude id:%d phase id:%d\n",filename, ampVar->id, phaVar->id);
    
    ampUnitAtt=ampVar->attributes.findP("units");
    if(ampUnitAtt){
      ampUnit=ampUnitAtt->as_string();
      }
    
    if(ampUnit=="cm")
      factor = 1.0e-02;
    else if(ampUnit=="m")
      factor = 1.0;
    else{
      if(gCPU_ID==gCPU_MASTER) printf("no known units given in %s for "+ampVar->name+", assume loading in meters\n", filename);
      factor = 1.0;
      }
    
    z = new fcomplex[length];
    
    status=poc_get_cvara(filename, ampVar->name, phaVar->name, -1, z, verbose);
    cmask=NC_FILL_DOUBLE;
    if(status != 0)  goto error;

    if(factor!=1.)
    for(m = 0; m < length; m++) {
      if(z[m] != cmask) {
       z[m] *= factor;
       }
      }
    
    for(n = 0; n < nndes; n++) {
      x = data.lon[n] * 180. / M_PI;
      y = data.lat[n] * 180. / M_PI;
      if(x < LSAgrid.xmin)
        x = x + 360.0;
      if(x > LSAgrid.xmax)
        x = x - 360.0;
      if(x < LSAgrid.xmin - LSAgrid.dx / 2.)
        x = x + 360.0;
      if(x > LSAgrid.xmax + LSAgrid.dx / 2.)
        x = x - 360.0;
      status = map_interpolation(LSAgrid, z, cmask, x, y, &(data.LSA[n].z[k]));
      if(status != 0) {
        map_printgrid(LSAgrid);
        printf("interpolation error: node %d lon=%lf lat=%lf \n", n, x, y);
        goto error;
        }
      }
    found[k] = 1;
    zaparr(z);

    LSAgrid.free();
    
    delete[] filename;
    }

  status=0;
  return (status);

error:
  return (status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int initialize_LSA(mesh_t & mesh, int paire, parameter_t *data)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**----------------------------------------------------------------------
  Loading and self-attraction from pre-computed fields*/
{
  int i, k, m, n, items, status;
  int zdim,udim;
  int z_discretisation, u_discretisation;
  int *found;
  complex<float> *buffer, mask;
  char *comment[3];
  char filename[1024];
//   int fmt=MESH_FAMILY_QUODDY; /// HERE !!!
  
  if(gCPU_ID==gCPU_MASTER) {
    printf(SEPARATOR_1);
    printf("initialise loading and self-attraction terms\n");
    }

  status=paire_dimension(mesh,paire,&zdim,&udim);
  status=paire_discretisation_id(paire, &z_discretisation, &u_discretisation);

  data->nLSA=zdim;
  for(n = 0; n < zdim; n++) {
    data->LSA[n].z = new fcomplex[WaveList.n];
    }

  for(k = 0; k < WaveList.n; k++) {
    for(n = 0; n < zdim; n++) {
      data->LSA[n].z[k] = fcomplex(0.,0.);
      }
    }

  found = new int[WaveList.n];
  for(k = 0; k < WaveList.n; k++)
    found[k] = 0;

/*----------------------------------------------------------------------
  comment:

  there is a significant danger if loading files and wave list are
  not consistent, also admittance added waves are poorly treated

  something need to be done to avoid future problems

----------------------------------------------------------------------*/

/**----------------------------------------------------------------------------
  first try the old ascii format (given at LGP1 nodes) */
  status = read_loading01(zdim, *data, found);

/**----------------------------------------------------------------------------
  then try the netcf format (structured grid database) */
  status = read_loading02(zdim, *data, found);
  check_error(status, "loading failed", __LINE__, __FILE__, 1);

/**----------------------------------------------------------------------------
  if allowed, use admittance to extend LSA database to minor astronomical tides */
//   if(admittance)
//     status = admittance_loading(zdim, *data, found); /// HERE !!!

  for(k = 0; k < WaveList.n; k++) {
    if(found[k] == 0) {
      if(gCPU_ID==gCPU_MASTER) __OUT_BASE_LINE__("Warning : wave %10s has no valid loading : %gcm\n", WaveList.waves[k].name, WaveList.waves[k].Ap);
      }
    }

  data->have_LSA=found;
  
/**----------------------------------------------------------------------------
  some ouputs */
  comment[0]=new char[256];
  comment[1]=new char[256];

  buffer=new complex<float> [zdim];
  for(k = 0; k < WaveList.n; k++) {
    if(found[k] == 0) continue;
    for(n = 0; n < zdim; n++) {
      buffer[n]=data->LSA[n].z[k];
      }
    sprintf(filename, "%s/%s.load.s2c", gOutputPath,WaveList.waves[k].name);
    sprintf(comment[0], "snapshot, %s", tugo_version);
    sprintf(comment[1], "Loading and self-attraction (m)\n");

//     status=fe_ascii_savec1((const char*) filename, mesh, buffer, mask, z_discretisation, fmt, comment); /// HERE !!!

    
    }
  delete[] buffer;

  delete[] comment[0];
  delete[] comment[1];

  if(gCPU_ID==gCPU_MASTER) {
    printf(SEPARATOR_2);
    }
    
  return (status);
}
