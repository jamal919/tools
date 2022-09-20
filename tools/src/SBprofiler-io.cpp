
/**************************************************************************

  T-UGO tools, 2006-2009

  Unstructured Ocean Grid initiative

Contributors:

  Florent Lyard      LEGOS/CNRS, Toulouse, France
  Yoann Le Bars      LEGOS, Toulouse, France (PhD)
  Laurent Roblou     LEGOS/CNRS, Toulouse, France
  Cyril Nguen        LA, Toulouse, France
  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
  Frédéric Dupont    Université de Laval à Québec, Canada

E-mail: florent.lyard@legos.obs-mip.fr

***************************************************************************/

#include "config.h"

#include <stdio.h>
#include <string.h>
#include <fstream>
#include <string>

#include "tools-structures.h"

#include "tides.h"
#include "fe.h"
#include "archive.h"
#include "netcdf-proto.h"
#include "poc-time.h"
#include "mgr.h"
#include "functions.h"

using namespace std;

/* *----------------------------------------------------------------------------
Seabird profiler ascii format
...
** Station: 10A
** Date(UTC): 12-Apr-2011
** Time(UTC): 22:33
** Lat: 37-44.99N
** Long: 141-05.03E
...
# nquan = 14
# nvalues = 20
# units = specified
# name 0 = scan: Scan Count
# name 1 = prDM: Pressure, Digiquartz [db]
# name 2 = depSM: Depth [salt water, m]
# name 3 = t090C: Temperature [ITS-90, deg C]
# name 4 = c0S/m: Conductivity [S/m]
# name 5 = altM: Altimeter [m]
# name 6 = dz/dtM: Descent Rate [m/s]
# name 7 = modError: Modulo Error Count
# name 8 = pumps: Pump Status
# name 9 = nbin: number of scans per bin
# name 10 = sal00: Salinity [PSU]
# name 11 = sigma-é00: Density [sigma-theta, Kg/m^3]
# name 12 = potemp090C: Potential Temperature [ITS-90, deg C]
# name 13 = flag: flag
# span 0 =      16258,      20359
# span 1 =      4.000,     23.000
# span 2 =      3.972,     22.831
# span 3 =     8.1214,     8.2987
# span 4 =   3.474641,   3.485433
# span 5 =       9.72,      28.94
# span 6 =      0.085,      0.385
# span 7 =          0,          0
# span 8 =          1,          1
# span 9 =         51,        191
# span 10 =    33.2431,    33.3024
# span 11 =    25.8509,    25.9232
# span 12 =     8.1191,     8.2969
# span 13 = 0.0000e+00, 0.0000e+00
# interval = decibars: 1
# start_time = Apr 12 2011 22:33:19
# bad_flag = -9.990e-29
...
*END*
*/

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int SBprofiler_info(char *filename, double **time)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int ncid,status;
  int time_id,reference_id;
  int vdim,tdim;
  int u_id,v_id,w_id;
  cdfvar_t u_info,v_info,w_info,reference_info;
  cdfgbl_t global;
  char *s;
  date_t *origin;
  int hour,minute,seconds;
  FILE *in;

  origin=new date_t;

  int nlevels, nsamples;

  std::ifstream input(filename);

  int cnt = 0;
  std::string line;
  do{
    std::getline(input, line);
    cnt++;
    }
  while( line.compare(0, 5, "*End*") != 0);

  status = nc_inq_varid(ncid, "JULD", &time_id);
  status = nc_inq_varid(ncid, "REFERENCE_DATE_TIME", &reference_id);

  status = nc_inq_varid(ncid, "UVEL_ADCP", &u_id);
  status = nc_inq_varid(ncid, "VVEL_ADCP", &v_id);
  status = nc_inq_varid(ncid, "WVEL_ADCP", &w_id);

  vdim=cdf_identify_dimension(global,"N_LEVEL");
  tdim=cdf_identify_dimension(global,"N_DATE_TIME");

  status= cdf_varinfo(ncid, u_id, &u_info);

  s=new char[32];
  status=nc_get_var_text(ncid,reference_id,s);
  nc_check_error(status,__LINE__,__FILE__);

  sscanf(s,"%4d%2d%2d%2d%2d%2d",&(origin->year),&(origin->month),&(origin->day),&hour,&minute,&seconds);
  origin->second=hour*3600.+minute*60.+seconds;

  nlevels=global.dimension[vdim].length;
  nsamples=global.dimension[tdim].length;

  *time=new double[nsamples];
  status=nc_get_var_double(ncid,time_id,*time);
  nc_check_error(status,__LINE__,__FILE__);

//  initial = time_getcnesdate(first, 's');
  status = nc_close(ncid);

}

