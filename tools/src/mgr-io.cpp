

/*******************************************************************************

  T-UGO tools, 2006-2009

  Unstructured Ocean Grid initiative

Contributors:

  Florent Lyard      LEGOS/CNRS, Toulouse, France
  Yoann Le Bars      LEGOS, Toulouse, France (PhD)
  Laurent Roblou     LEGOS/CNRS, Toulouse, France
  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
  Frédéric Dupont    Université de Laval à Québec, Canada

E-mail: florent.lyard@legos.obs-mip.fr

*******************************************************************************/

#include <config.h>

#include <vector>
#include <iomanip>
#include <fstream>

#include "tools-structures.h"
#include "netcdf-proto.h"
#include "mgr.h"
#include "tides.h"
#include "tides.def"
#include "list.h"
#include "zapper.h"


inline bool wcompare(tidal_wave first, tidal_wave second){
  return(first.omega < second.omega);
}

inline bool compare_nocase (tidal_wave first, tidal_wave second)
{
  return -strcasecmp(first.name,second.name);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int mgr_save(const char *filename, vector<mgr_t> mgr, const string & format, const char *wave)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int nmgr=mgr.size();
  
  if( mgr_compact(mgr, nmgr) != 0){
    return(1);
    }
  
  if(format=="ASCII"){
    return(mgr_save_ascii(filename , mgr));
    }
  else if(format=="LEGOS-ASCII"){
    return(mgr_save_ascii(filename , mgr));
    }
  else if(format=="NETCDF"){
    return(mgr_save_netcdf(filename, mgr));
    }
  else if(format=="LEGOS-NETCDF"){
    return(mgr_save_netcdf(filename, mgr));
    }
  else if(format=="LEGOS-OBS"){
    return(mgr_save_obs(filename, mgr,wave));
    }
  else if(format=="BIO"){
    return(mgr_save_BIO(filename, mgr));
    }
  else {
    printf("invalid format for mgr file: %s\n",format.c_str());
    }

  return(1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

string shellGetDate()

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  time_t current;

  time(&current);
  return(ctime(&current));
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int get_nWave(vector<mgr_t>mgr)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int nWave = 76; // max
  const int nmgr=mgr.size();
  
  if(nmgr<1)
    return 0;
  
  size_t idx = 0;
  for(size_t m = 0; m < nmgr; m++){
    if(mgr[m].nwave > mgr[idx].nwave){
      idx = m;
      }
    }
  nWave = mgr[idx].nwave;
  
  return(nWave);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void createFile(const char *filename, vector<mgr_t> mgr_serie)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
    int status = NC_NOERR + 13;
    
    int nRecords=mgr_serie.size();
    
    // create file
    int  ncid = NC_EBADID;
    status = nc_create(filename, NC_CLOBBER, &ncid);
    NC_TRAP_ERROR(wexit,status,1,"");
   
    // define dimensions
    int spectrum_dim = -1;
    status = NC_NOERR + 13;
    int nWaves = get_nWave(mgr_serie);
    status = nc_def_dim(ncid, "constituent", nWaves, &spectrum_dim);
    NC_TRAP_ERROR(wexit,status,1,"");
    int strlen_dim = -1;
    status = NC_NOERR + 13;
    status = nc_def_dim(ncid, "namelength", 10, &strlen_dim);
    NC_TRAP_ERROR(wexit,status,1,"");
    int records_dim = -1;
    status = NC_NOERR + 13;
    status = nc_def_dim(ncid, "records", nRecords, &records_dim);
    NC_TRAP_ERROR(wexit,status,1,"");
    
    // set global attributes
    string field = "CF 1.4";
    status = NC_NOERR + 13;
    status = nc_put_att_text(ncid, NC_GLOBAL, "Conventions", field.size()+1,field.c_str());
    NC_TRAP_ERROR(wexit,status,1,"");
    status = NC_NOERR + 13;
    if( (status = nc_put_att_text(ncid, NC_GLOBAL, "file_name", strlen(filename), filename)) != NC_NOERR)
      NC_TRAP_ERROR(wexit,status,1,"");
    string keyword = shellGetDate();
    status = NC_NOERR + 13;
    status = nc_put_att_text(ncid, NC_GLOBAL, "date", keyword.size()-1,keyword.c_str());
    NC_TRAP_ERROR(wexit,status,1,"");
    field = "altimetry-detidor";
    status = NC_NOERR + 13;
    status = nc_put_att_text(ncid, NC_GLOBAL, "software", field.size()+1,field.c_str());
    NC_TRAP_ERROR(wexit,status,1,"");
    field = "2.0";
    status = NC_NOERR + 13;
    status = nc_put_att_text(ncid, NC_GLOBAL, "version", field.size()+1,field.c_str());
    NC_TRAP_ERROR(wexit,status,1,"");
    field = "CTOH/LEGOS, http://ctoh.legos.obs-mip.fr";
    status = NC_NOERR + 13;
    status = nc_put_att_text(ncid, NC_GLOBAL, "production",field.size()+1,field.c_str());
    NC_TRAP_ERROR(wexit,status,1,"");
    field = "CTOH/LEGOS, http://ctoh.legos.obs-mip.fr";
    status = NC_NOERR + 13;
    status = nc_put_att_text(ncid, NC_GLOBAL, "institution",field.size()+1,field.c_str());
    NC_TRAP_ERROR(wexit,status,1,"");
    field = "product creation " + keyword;
    status = NC_NOERR + 13;
    status = nc_put_att_text(ncid, NC_GLOBAL, "history", field.size()-1,field.c_str());
    NC_TRAP_ERROR(wexit,status,1,"");
    field = "Tidal harmonic constants along the long-term TOPEX/Poseidon, Jason-1/2 ground track";
    status = NC_NOERR + 13;
    status = nc_put_att_text(ncid, NC_GLOBAL, "title", field.size()+1,field.c_str());
    NC_TRAP_ERROR(wexit,status,1,"");
    field = "ctoh_products@legos.obs-mip.fr";
    status = NC_NOERR + 13;
    status = nc_put_att_text(ncid, NC_GLOBAL, "contact", field.size()+1,field.c_str());
    NC_TRAP_ERROR(wexit,status,1,"");
    
    field = mgr_serie[0].debut;
    size_t pos = 0;
    while ( (pos = field.find('/')) != string::npos) {
      field.replace(pos, 1, "-");
    }
    status = NC_NOERR + 13;
    status = nc_put_att_text(ncid, NC_GLOBAL, "time_series_first_date_used", field.size()+1,field.c_str());
    NC_TRAP_ERROR(wexit,status,1,"");
      
    field = mgr_serie[0].fin;
    while ( (pos = field.find('/')) != string::npos) {
      field.replace(pos, 1, "-");
    }
    status = NC_NOERR + 13;
    status = nc_put_att_text(ncid, NC_GLOBAL, "time_series_last_date_used", field.size()+1,field.c_str());
    NC_TRAP_ERROR(wexit,status,1,"");
    field = mgr_serie[0].name;
    status = NC_NOERR + 13;
    status = nc_put_att_text(ncid, NC_GLOBAL, "time_series_file_name", field.size()+1,field.c_str());
    NC_TRAP_ERROR(wexit,status,1,"");
    field = mgr_serie[0].origine;
    status = NC_NOERR + 13;
    status = nc_put_att_text(ncid, NC_GLOBAL, "time_series_origin", field.size()+1,field.c_str());
    NC_TRAP_ERROR(wexit,status,1,"");
    field = mgr_serie[0].validation;
    status = NC_NOERR + 13;
    status = nc_put_att_text(ncid, NC_GLOBAL, "time_series_validation", field.size()+1,field.c_str());
    NC_TRAP_ERROR(wexit,status,1,"");
    
    
    // complete file
    status = NC_NOERR + 13;
    status = nc_enddef (ncid);
    NC_TRAP_ERROR(wexit,status,1,"");
    status = NC_NOERR + 13;
    status = nc_close(ncid);
    NC_TRAP_ERROR(wexit,status,1,"");

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void putVariables(const char *filename, vector<mgr_t> mgr_serie)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
    string attribute="";
    double mask = static_cast<double>(1.0e35);
    double value = mask;
    short ivalue = static_cast<int>(32767);

    int nmgr=mgr_serie.size();
    int status = NC_NOERR + 13;

    int ncid = NC_EBADID;
    status = nc_open(filename, NC_WRITE, &ncid);
    NC_TRAP_ERROR(wexit,status,1,"");

    // enter define mode
    status = nc_redef(ncid);
    NC_TRAP_ERROR(wexit,status,1,"");


    // get dimension ids
    int spectrum_dimid = -1;
    status = NC_NOERR + 13;
    status = nc_inq_dimid(ncid,"constituent",&spectrum_dimid);
    NC_TRAP_ERROR(wexit,status,1,"");

    int strlen_dimid = -1;
    status = NC_NOERR + 13;
    status = nc_inq_dimid(ncid,"namelength",&strlen_dimid);
    NC_TRAP_ERROR(wexit,status,1,"");

    int records_dimid = -1;
    status = NC_NOERR + 13;
    status = nc_inq_dimid(ncid,"records",&records_dimid);
    NC_TRAP_ERROR(wexit,status,1,"");


    // get dimension length
    size_t spectrum_dimlen;
    status = NC_NOERR + 13;
    status = nc_inq_dimlen(ncid,spectrum_dimid,&spectrum_dimlen);
    NC_TRAP_ERROR(wexit,status,1,"");

    size_t strlen_dimlen;
    status = NC_NOERR + 13;
    status = nc_inq_dimlen(ncid,strlen_dimid,&strlen_dimlen);
    NC_TRAP_ERROR(wexit,status,1,"");

    size_t records_dimlen;
    status = NC_NOERR + 13;
    status = nc_inq_dimlen(ncid,records_dimid,&records_dimlen);
    NC_TRAP_ERROR(wexit,status,1,"");


    // define variables
    int constituent_id = NC_EBADID;
    int constituent_dims[2] = {spectrum_dimid,strlen_dimid};
    status = NC_NOERR + 13;
    status = nc_def_var(ncid, "constituentname", NC_CHAR, 2, constituent_dims, &constituent_id);
    NC_TRAP_ERROR(wexit,status,1,"");
    attribute = "tidal_spectrum";
    status = nc_put_att_text(ncid,constituent_id,"long_name",attribute.size(),attribute.c_str());
    NC_TRAP_ERROR(wexit,status,1,"");
    attribute = "tidal_spectrum";
    status = nc_put_att_text(ncid,constituent_id,"standard_name",attribute.size(),attribute.c_str());
    NC_TRAP_ERROR(wexit,status,1,"");
    attribute = "spectrum";
    status = nc_put_att_text(ncid,constituent_id,"short_name",attribute.size(),attribute.c_str());
    NC_TRAP_ERROR(wexit,status,1,"");

    int lat_id = NC_EBADID;
    int lat_dims[1] = {records_dimid};
    status = NC_NOERR + 13;
    status = nc_def_var(ncid, "lat", NC_DOUBLE, 1, lat_dims, &lat_id);
    NC_TRAP_ERROR(wexit,status,1,"");
    attribute = "latitude";
    status = nc_put_att_text(ncid,lat_id,"long_name",attribute.size(),attribute.c_str());
    NC_TRAP_ERROR(wexit,status,1,"");
    attribute = "latitude";
    status = nc_put_att_text(ncid,lat_id,"standard_name",attribute.size(),attribute.c_str());
    NC_TRAP_ERROR(wexit,status,1,"");
    attribute = "lat";
    status = nc_put_att_text(ncid,lat_id,"short_name",attribute.size(),attribute.c_str());
    NC_TRAP_ERROR(wexit,status,1,"");
    attribute = "degrees_north";
    status = nc_put_att_text(ncid,lat_id,"units",attribute.size(),attribute.c_str());
    NC_TRAP_ERROR(wexit,status,1,"");
    value = static_cast<double>(-90.0);
    status = nc_put_att_double(ncid,lat_id,"valid_min", NC_DOUBLE, 1, &value);
    NC_TRAP_ERROR(wexit,status,1,"");
    value = 90.0;
    status = nc_put_att_double(ncid,lat_id,"valid_max", NC_DOUBLE, 1, &value);
    NC_TRAP_ERROR(wexit,status,1,"");
    value = 1.0;
    status = nc_put_att_double(ncid,lat_id,"scale_factor", NC_DOUBLE, 1, &value);
    NC_TRAP_ERROR(wexit,status,1,"");
    value = 0.0;
    status = nc_put_att_double(ncid,lat_id,"add_offset", NC_DOUBLE, 1, &value);
    NC_TRAP_ERROR(wexit,status,1,"");
    value = mask;
    status = nc_put_att_double(ncid,lat_id,"_FillValue", NC_DOUBLE, 1, &value);
    NC_TRAP_ERROR(wexit,status,1,"");

    int lon_id = NC_EBADID;
    int lon_dims[1] = {records_dimid};
    status = NC_NOERR + 13;
    status = nc_def_var(ncid, "lon", NC_DOUBLE, 1, lon_dims, &lon_id);
    NC_TRAP_ERROR(wexit,status,1,"");
    attribute = "longitude";
    status = nc_put_att_text(ncid,lon_id,"long_name",attribute.size(),attribute.c_str());
    NC_TRAP_ERROR(wexit,status,1,"");
    attribute = "longitude";
    status = nc_put_att_text(ncid,lon_id,"standard_name",attribute.size(),attribute.c_str());
    NC_TRAP_ERROR(wexit,status,1,"");
    attribute = "lon";
    status = nc_put_att_text(ncid,lon_id,"short_name",attribute.size(),attribute.c_str());
    NC_TRAP_ERROR(wexit,status,1,"");
    attribute = "degrees_east";
    status = nc_put_att_text(ncid,lon_id,"units",attribute.size(),attribute.c_str());
    NC_TRAP_ERROR(wexit,status,1,"");
    value = static_cast<double>(-180.0);
    status = nc_put_att_double(ncid,lon_id,"valid_min", NC_DOUBLE, 1, &value);
    NC_TRAP_ERROR(wexit,status,1,"");
    value = 180.0;
    status = nc_put_att_double(ncid,lon_id,"valid_max", NC_DOUBLE, 1, &value);
    NC_TRAP_ERROR(wexit,status,1,"");
    value = 1.0;
    status = nc_put_att_double(ncid,lon_id,"scale_factor", NC_DOUBLE, 1, &value);
    NC_TRAP_ERROR(wexit,status,1,"");
    value = 0.0;
    status = nc_put_att_double(ncid,lon_id,"add_offset", NC_DOUBLE, 1, &value);
    NC_TRAP_ERROR(wexit,status,1,"");
    value = mask;
    status = nc_put_att_double(ncid,lon_id,"_FillValue", NC_DOUBLE, 1, &value);
    NC_TRAP_ERROR(wexit,status,1,"");


    int amplitude_id = NC_EBADID;
    int amplitude_dims[2] = {records_dimid, spectrum_dimid};
    status = NC_NOERR + 13;
    status = nc_def_var(ncid, "amplitude", NC_DOUBLE, 2, amplitude_dims, &amplitude_id);
    NC_TRAP_ERROR(wexit,status,1,"");
    attribute = "tidal_elevation_amplitude";
    status = nc_put_att_text(ncid,amplitude_id,"long_name",attribute.size(),attribute.c_str());
    NC_TRAP_ERROR(wexit,status,1,"");
    attribute = "amplitude";
    status = nc_put_att_text(ncid,amplitude_id,"standard_name",attribute.size(),attribute.c_str());
    NC_TRAP_ERROR(wexit,status,1,"");
    attribute = "Ha";
    status = nc_put_att_text(ncid,amplitude_id,"short_name",attribute.size(),attribute.c_str());
    NC_TRAP_ERROR(wexit,status,1,"");
    attribute = "m";
    status = nc_put_att_text(ncid,amplitude_id,"units",attribute.size(),attribute.c_str());
    NC_TRAP_ERROR(wexit,status,1,"");
    value = 0.0;
    status = nc_put_att_double(ncid,amplitude_id,"valid_min", NC_DOUBLE, 1, &value);
    NC_TRAP_ERROR(wexit,status,1,"");
/*    value = 300.0;//!!
    status = nc_put_att_double(ncid,amplitude_id,"valid_max", NC_DOUBLE, 1, &value);
    NC_TRAP_ERROR(wexit,status,1,"");*/
    value = 1.0;
    status = nc_put_att_double(ncid,amplitude_id,"scale_factor", NC_DOUBLE, 1, &value);
    NC_TRAP_ERROR(wexit,status,1,"");
    value = 0.;
    status = nc_put_att_double(ncid,amplitude_id,"add_offset", NC_DOUBLE, 1, &value);
    NC_TRAP_ERROR(wexit,status,1,"");
    value = mask;
    status = nc_put_att_double(ncid,amplitude_id,"_FillValue", NC_DOUBLE, 1, &value);
    NC_TRAP_ERROR(wexit,status,1,"");


    int phase_id = NC_EBADID;
    int phase_dims[2] = {records_dimid, spectrum_dimid};
    status = NC_NOERR + 13;
    status = nc_def_var(ncid, "phase_lag", NC_DOUBLE, 2, phase_dims, &phase_id);
    NC_TRAP_ERROR(wexit,status,1,"");
    attribute = "tidal_elevation_phase_lag";
    status = nc_put_att_text(ncid,phase_id,"long_name",attribute.size(),attribute.c_str());
    NC_TRAP_ERROR(wexit,status,1,"");
    attribute = "phase_lag";
    status = nc_put_att_text(ncid,phase_id,"standard_name",attribute.size(),attribute.c_str());
    NC_TRAP_ERROR(wexit,status,1,"");
    attribute = "Hg";
    status = nc_put_att_text(ncid,phase_id,"short_name",attribute.size(),attribute.c_str());
    NC_TRAP_ERROR(wexit,status,1,"");
    attribute = "degrees";
    status = nc_put_att_text(ncid,phase_id,"units",attribute.size(),attribute.c_str());
    NC_TRAP_ERROR(wexit,status,1,"");
    value = 0.0;
    status = nc_put_att_double(ncid,phase_id,"valid_min", NC_DOUBLE, 1, &value);
    NC_TRAP_ERROR(wexit,status,1,"");
    value = 360.0;
    status = nc_put_att_double(ncid,phase_id,"valid_max", NC_DOUBLE, 1, &value);
    NC_TRAP_ERROR(wexit,status,1,"");
    value = 1.0;
    status = nc_put_att_double(ncid,phase_id,"scale_factor", NC_DOUBLE, 1, &value);
    NC_TRAP_ERROR(wexit,status,1,"");
    value = 0.0;
    status = nc_put_att_double(ncid,phase_id,"add_offset", NC_DOUBLE, 1, &value);
    NC_TRAP_ERROR(wexit,status,1,"");
    value = mask;
    status = nc_put_att_double(ncid,phase_id,"_FillValue", NC_DOUBLE, 1, &value);
    NC_TRAP_ERROR(wexit,status,1,"");

    int error_id = NC_EBADID;
    int error_dims[2] = {records_dimid, spectrum_dimid};
    status = NC_NOERR + 13;
    status = nc_def_var(ncid, "mean_square_error", NC_DOUBLE, 2, error_dims, &error_id);
    NC_TRAP_ERROR(wexit,status,1,"");
    attribute = "tidal_analysis_mean_square_error";
    status = nc_put_att_text(ncid,error_id,"long_name",attribute.size(),attribute.c_str());
    NC_TRAP_ERROR(wexit,status,1,"");
    attribute = "error";
    status = nc_put_att_text(ncid,error_id,"standard_name",attribute.size(),attribute.c_str());
    NC_TRAP_ERROR(wexit,status,1,"");
    attribute = "mse";
    status = nc_put_att_text(ncid,error_id,"short_name",attribute.size(),attribute.c_str());
    NC_TRAP_ERROR(wexit,status,1,"");
    attribute = "m";
    status = nc_put_att_text(ncid,error_id,"units",attribute.size(),attribute.c_str());
    NC_TRAP_ERROR(wexit,status,1,"");
    value = 0.0;
    status = nc_put_att_double(ncid,error_id,"valid_min", NC_DOUBLE, 1, &value);
    NC_TRAP_ERROR(wexit,status,1,"");
/*    value = 300.0;//!!
    status = nc_put_att_double(ncid,error_id,"valid_max", NC_DOUBLE, 1, &value);
    NC_TRAP_ERROR(wexit,status,1,"");*/
    value = 1.0;
    status = nc_put_att_double(ncid,error_id,"scale_factor", NC_DOUBLE, 1, &value);
    NC_TRAP_ERROR(wexit,status,1,"");
    value = 0.;
    status = nc_put_att_double(ncid,error_id,"add_offset", NC_DOUBLE, 1, &value);
    NC_TRAP_ERROR(wexit,status,1,"");
    value = mask;
    status = nc_put_att_double(ncid,error_id,"_FillValue", NC_DOUBLE, 1, &value);
    NC_TRAP_ERROR(wexit,status,1,"");

    int bg_contamination_error_id = NC_EBADID;
    int bg_contamination_error_dims[2] = {records_dimid, spectrum_dimid};
    status = NC_NOERR + 13;
    status = nc_def_var(ncid, "bg_contamination_error", NC_DOUBLE, 2, bg_contamination_error_dims, &bg_contamination_error_id);
    NC_TRAP_ERROR(wexit,status,1,"");
    attribute = "tidal analysis error from background contamination";
    status = nc_put_att_text(ncid, bg_contamination_error_id,"long_name",attribute.size(),attribute.c_str());
    NC_TRAP_ERROR(wexit,status,1,"");
    attribute = "background noise error";
    status = nc_put_att_text(ncid, bg_contamination_error_id,"standard_name",attribute.size(),attribute.c_str());
    NC_TRAP_ERROR(wexit,status,1,"");
    attribute = "background error";
    status = nc_put_att_text(ncid, bg_contamination_error_id,"short_name",attribute.size(),attribute.c_str());
    NC_TRAP_ERROR(wexit,status,1,"");
    attribute = "m";
    status = nc_put_att_text(ncid, bg_contamination_error_id,"units",attribute.size(),attribute.c_str());
    NC_TRAP_ERROR(wexit,status,1,"");
    value = 0.0;
    status = nc_put_att_double(ncid, bg_contamination_error_id,"valid_min", NC_DOUBLE, 1, &value);
    NC_TRAP_ERROR(wexit,status,1,"");
/*    value = 300.0;//!!
    status = nc_put_att_double(ncid, bg_contamination_error_id,"valid_max", NC_DOUBLE, 1, &value);
    NC_TRAP_ERROR(wexit,status,1,"");*/
    value = 1.0;
    status = nc_put_att_double(ncid, bg_contamination_error_id,"scale_factor", NC_DOUBLE, 1, &value);
    NC_TRAP_ERROR(wexit,status,1,"");
    value = 0.;
    status = nc_put_att_double(ncid, bg_contamination_error_id,"add_offset", NC_DOUBLE, 1, &value);
    NC_TRAP_ERROR(wexit,status,1,"");
    value = mask;
    status = nc_put_att_double(ncid, bg_contamination_error_id,"_FillValue", NC_DOUBLE, 1, &value);
    NC_TRAP_ERROR(wexit,status,1,"");

    int nmes_id = NC_EBADID;
    int nmes_dims[1] = {records_dimid};
    status = NC_NOERR + 13;
    status = nc_def_var(ncid, "sample", NC_SHORT, 1, nmes_dims, &nmes_id);
    NC_TRAP_ERROR(wexit,status,1,"");
    attribute = "number_of_samples";
    status = nc_put_att_text(ncid,nmes_id,"long_name",attribute.size(),attribute.c_str());
    NC_TRAP_ERROR(wexit,status,1,"");
    attribute = "samples";
    status = nc_put_att_text(ncid,nmes_id,"standard_name",attribute.size(),attribute.c_str());
    NC_TRAP_ERROR(wexit,status,1,"");
    attribute = "samples";
    status = nc_put_att_text(ncid,nmes_id,"short_name",attribute.size(),attribute.c_str());
    NC_TRAP_ERROR(wexit,status,1,"");
    ivalue = static_cast<int>(0);
    status = nc_put_att_short(ncid,nmes_id,"valid_min", NC_SHORT, 1, &ivalue);
    NC_TRAP_ERROR(wexit,status,1,"");
    ivalue = static_cast<int>(1000); //!!!
    status = nc_put_att_short(ncid,nmes_id,"valid_max", NC_SHORT, 1, &ivalue);
    NC_TRAP_ERROR(wexit,status,1,"");
    ivalue = static_cast<int>(1);
    status = nc_put_att_short(ncid,nmes_id,"scale_factor", NC_SHORT, 1, &ivalue);
    NC_TRAP_ERROR(wexit,status,1,"");
    ivalue = static_cast<int>(0);
    status = nc_put_att_short(ncid,nmes_id,"add_offset", NC_SHORT, 1, &ivalue);
    NC_TRAP_ERROR(wexit,status,1,"");
    ivalue = static_cast<int>(32767);
    status = nc_put_att_short(ncid,nmes_id,"_FillValue", NC_SHORT, 1, &ivalue);
    NC_TRAP_ERROR(wexit,status,1,"");

    // end definitions
    status = NC_NOERR + 13;
    status = nc_enddef(ncid);
    NC_TRAP_ERROR(wexit,status,1,"");

    // provide values
    size_t idx = 0;
    for(size_t mgr = 0; mgr < nmgr; mgr++){
      if(mgr_serie[mgr].nwave > mgr_serie[idx].nwave){
         idx = mgr;
        }
      }

    // constituent name
    spectrum_t fullSpectrum;
    fullSpectrum.init(mgr_serie[idx].nwave);

    size_t start[2] = {0,0};
    size_t count[2] = {1,strlen_dimlen};
    char wName[strlen_dimlen];
    for(size_t w = 0; w < mgr_serie[idx].nwave; w++){
      start[0] = w;
      strcpy(wName, mgr_serie[idx].data[w].constituent.name);
      count[1] = strlen(wName);
      status = nc_put_vara_text(ncid, constituent_id, start, count, wName);
      NC_TRAP_ERROR(wexit,status,1,"");
      tidal_wave wave = wave_from_name(wName);
      fullSpectrum.add(wave);
      wave.reset();
    }

    // latitude
    double *lat = new double[records_dimlen];
    for(size_t mgr = 0; mgr < nmgr; mgr++){
      lat[mgr] = mgr_serie[mgr].loc.lat;
    }
    status = nc_put_var_double(ncid,lat_id,lat);
    NC_TRAP_ERROR(wexit,status,1,"");
    delete [] lat;

    // longitude
    double *lon = new double[records_dimlen];
    for(size_t mgr = 0; mgr < nmgr; mgr++){
      lon[mgr] = mgr_serie[mgr].loc.lon;
    }
    status = nc_put_var_double(ncid,lon_id,lon);
    NC_TRAP_ERROR(wexit,status,1,"");
    delete [] lon;


    // amplitude
      double *a = new double[records_dimlen * spectrum_dimlen];
      for(size_t i = 0; i < records_dimlen * spectrum_dimlen; ++i){
        a[i] = mask;
      }
      
      for(size_t mgr = 0; mgr < nmgr; mgr++){
        for(size_t w = 0; w < mgr_serie[mgr].nwave; w++){
          int k = fullSpectrum.wave_index(mgr_serie[mgr].data[w].constituent.name);
          if(k != -1){
            a[mgr * spectrum_dimlen + k] = mgr_serie[mgr].data[w].amp;
          }
        }
      }
    status = nc_put_var_double(ncid,amplitude_id,a);
    NC_TRAP_ERROR(wexit,status,1,"");
    delete [] a;


    // phase lag
    double *g = new double[records_dimlen * spectrum_dimlen];
    for(size_t i = 0; i < records_dimlen * spectrum_dimlen; ++i){
      g[i] = mask;
    }
    for(size_t mgr = 0; mgr < nmgr; mgr++){
      for(size_t w = 0; w < mgr_serie[mgr].nwave; w++){
        int k = fullSpectrum.wave_index(mgr_serie[mgr].data[w].constituent.name);
        if(k != -1){
          g[mgr * spectrum_dimlen + k] = mgr_serie[mgr].data[w].phi;
        }
      }
    }
    status = nc_put_var_double(ncid,phase_id,g);
    NC_TRAP_ERROR(wexit,status,1,"");
    delete [] g;


    // LSE formal error
    double *e = new double[records_dimlen * spectrum_dimlen];
    for(size_t i = 0; i < records_dimlen * spectrum_dimlen; ++i){
      e[i] = mask;
    }
    for(size_t mgr = 0; mgr < nmgr; mgr++){
      for(size_t w = 0; w < mgr_serie[mgr].nwave; w++){
        int k = fullSpectrum.wave_index(mgr_serie[mgr].data[w].constituent.name);
        if(k != -1){
          e[mgr * spectrum_dimlen + k] = mgr_serie[mgr].data[w].error;
        }
      }
    }
    status = nc_put_var_double(ncid, error_id, e);
    NC_TRAP_ERROR(wexit,status,1,"");


    // ocean back-ground contamination error
    for(size_t i = 0; i < records_dimlen * spectrum_dimlen; ++i){
      e[i] = mask;
      }
    for(size_t mgr = 0; mgr < nmgr; mgr++){
      for(size_t w = 0; w < mgr_serie[mgr].nwave; w++){
        int k = fullSpectrum.wave_index(mgr_serie[mgr].data[w].constituent.name);
        if(k != -1){
          e[mgr * spectrum_dimlen + k] = mgr_serie[mgr].data[w].bg_contamination_error;
          }
        }
      }
    status = nc_put_var_double(ncid, bg_contamination_error_id, e);
    NC_TRAP_ERROR(wexit,status,1,"");
    delete [] e;

    // samples used
    double *nmes = new double[records_dimlen];
    for(size_t mgr = 0; mgr < nmgr; mgr++){
      nmes[mgr] = mgr_serie[mgr].duree;
      }
    status = nc_put_var_double(ncid,nmes_id,nmes);
    NC_TRAP_ERROR(wexit,status,1,"");
    delete [] nmes;

    // complete file
    status = NC_NOERR + 13;
    status = nc_close(ncid);
    NC_TRAP_ERROR(wexit,status,1,"");


}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int mgr_save_netcdf(const char *filename, vector<mgr_t> mgr_serie)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  if(mgr_serie.size()<1)
    return -1;
  
  createFile(filename, mgr_serie);
  putVariables(filename, mgr_serie);

  return(0);

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int mgr_load_netcdf(const char *filename, vector<mgr_t> & mgr)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status = NC_NOERR + 13;
  
  int ncid = NC_EBADID;
  status = nc_open(filename, NC_NOWRITE, &ncid);
  NC_TRAP_ERROR(wexit,status,1,"");
  
  // get number of stations
  int dimid = NC_NOERR + 25;
  status = NC_NOERR + 13;
  status = nc_inq_dimid(ncid, "records", &dimid);
  NC_TRAP_ERROR(wexit,status,1,"nc_inq_dimid((\"%s\"),\"records\",) error",filename);
  size_t nmgr = NC_NOERR + 25;
  status = NC_NOERR + 13;
  status = nc_inq_dimlen(ncid, dimid, &nmgr);
  NC_TRAP_ERROR(wexit,status,1,"");


  // get number of constituents
  dimid = NC_NOERR + 25;
  status = NC_NOERR + 13;
  status = nc_inq_dimid(ncid, "constituent", &dimid);
  NC_TRAP_ERROR(wexit,status,1,"");
  size_t nwave =  NC_NOERR + 25;
  status = NC_NOERR + 13;
  status = nc_inq_dimlen(ncid, dimid, &nwave);
  NC_TRAP_ERROR(wexit,status,1,"");
  

  // get global attributes for common metadata
  char *sptr = new char[200];
  status = nc_get_att_text(ncid, NC_GLOBAL, "time_series_file_name", sptr);
  NC_TRAP_ERROR(wexit,status,1,"");
  string time_series_file_name = sptr;
  delete [] sptr;

  sptr = new char[200];
  status = nc_get_att_text(ncid, NC_GLOBAL, "time_series_origin", sptr);
  NC_TRAP_ERROR(wexit,status,1,"");
  string time_series_origin = sptr;
  delete [] sptr;

  sptr = new char[200];
  status = nc_get_att_text(ncid, NC_GLOBAL, "time_series_validation", sptr);
  NC_TRAP_ERROR(wexit,status,1,"");
  string time_series_validation = sptr;
  delete [] sptr;

  sptr = new char[200];
  status = nc_get_att_text(ncid, NC_GLOBAL, "time_series_first_date_used", sptr);
  NC_TRAP_ERROR(wexit,status,1,"");
  string time_series_first_date_used = sptr;
  delete [] sptr;


  sptr = new char[200];
  status = nc_get_att_text(ncid, NC_GLOBAL, "time_series_last_date_used", sptr);
  NC_TRAP_ERROR(wexit,status,1,"");
  string time_series_last_date_used = sptr;
  delete [] sptr;



  // get longitudes of stations
  string loc_units = "degrees";

  int varid = NC_NOERR + 25;
  status = NC_NOERR + 13;
  status = nc_inq_varid(ncid, "lon", &varid);
  NC_TRAP_ERROR(wexit,status,1,"");

  float *lon = new float[nmgr];
  status = NC_NOERR + 13;
  status = nc_get_var_float(ncid, varid, lon);
  NC_TRAP_ERROR(wexit,status,1,"");


  // get latitudes of stations
  varid = NC_NOERR + 25;
  status = NC_NOERR + 13;
  status = nc_inq_varid(ncid, "lat", &varid);
  NC_TRAP_ERROR(wexit,status,1,"");

  float *lat = new float[nmgr];
  status = NC_NOERR + 13;
  status = nc_get_var_float(ncid, varid, lat);
  NC_TRAP_ERROR(wexit,status,1,"");
    

  // get number of samples used for each stations
  varid = NC_NOERR + 25;
  status = NC_NOERR + 13;
  status = nc_inq_varid(ncid, "sample", &varid);
  NC_TRAP_ERROR(wexit,status,1,"");

  short *sample = new short[nmgr];
  status = NC_NOERR + 13;
  status = nc_get_var_short(ncid, varid, sample);
  NC_TRAP_ERROR(wexit,status,1,"");
    

  // get constituents list
  varid = NC_NOERR + 25;
  status = NC_NOERR + 13;
  status = nc_inq_varid(ncid, "constituentname", &varid);
  NC_TRAP_ERROR(wexit,status,1,"");

  dimid = NC_NOERR + 25;
  status = NC_NOERR + 13;
  status = nc_inq_dimid(ncid, "namelength", &dimid);
  NC_TRAP_ERROR(wexit,status,1,"");
  size_t namelength =  NC_NOERR + 25;
  status = NC_NOERR + 13;
  status = nc_inq_dimlen(ncid, dimid, &namelength);
  NC_TRAP_ERROR(wexit,status,1,"");

  status = NC_NOERR + 13;
  sptr = new char[nwave * namelength];
  status = nc_get_var_text(ncid, varid, sptr);
  NC_TRAP_ERROR(wexit,status,1,"");

  char **constituentname = new char*[nwave];
  for(size_t w = 0; w < nwave; ++w){
    constituentname[w] = new char[namelength];
    for(size_t c = 0; c < namelength; c++){
      constituentname[w][c] = sptr[w * namelength + c];
      }
    }
      

  // get values for amplitude
  int a_varid = NC_NOERR + 25;
  status = NC_NOERR + 13;
  status = nc_inq_varid(ncid, "amplitude", &a_varid);
  NC_TRAP_ERROR(wexit,status,1,"");

  double *a = new double[nmgr * nwave];
  status = NC_NOERR + 13;
  status = nc_get_var_double(ncid, a_varid, a);
  NC_TRAP_ERROR(wexit,status,1,"");

  float fillValue;
  status = NC_NOERR + 13;
  status = nc_get_att_float(ncid, a_varid, "_FillValue", &fillValue);
  NC_TRAP_ERROR(wexit,status,1,"");


  // get values for phase lag
  int g_varid = NC_NOERR + 25;
  status = NC_NOERR + 13;
  status = nc_inq_varid(ncid, "phase_lag", &g_varid);
  NC_TRAP_ERROR(wexit,status,1,"");

  double *g = new double[nmgr * nwave];
  status = NC_NOERR + 13;
  status = nc_get_var_double(ncid, g_varid, g);
  NC_TRAP_ERROR(wexit,status,1,"");


  // get values for mean square error
  int e_varid = NC_NOERR + 25;
  status = NC_NOERR + 13;
  status = nc_inq_varid(ncid, "mean_square_error", &e_varid);
  NC_TRAP_ERROR(wexit,status,1,"");

  double *e = new double[nmgr * nwave];
  status = NC_NOERR + 13;
  status = nc_get_var_double(ncid, e_varid, e);
  NC_TRAP_ERROR(wexit,status,1,"");


  // get variable ID for background contamination error (might not exist so request may fail)
  double *bg = 0;
  bool got_bg = false;

  int bg_varid = NC_NOERR + 25;
  status = NC_NOERR + 13;
  status = nc_inq_varid(ncid, "bg_contamination_error", &bg_varid);
  if(status!=0){
    nc_check_error(status, __LINE__, __FILE__, 0);
    }
  else{
    got_bg = true;
    bg = new double[nmgr * nwave];
    status = NC_NOERR + 13;
    status = nc_get_var_double(ncid, bg_varid, bg);
    NC_TRAP_ERROR(wexit,status,1,"");
    }

  // fill array of tidal constants
//   *mgr_serie = new mgr_t*[nmgr];
  spectrum_t reference_spectrum = initialize_tide();
  
  for(size_t m = 0; m < nmgr; m++){
    mgr_t gauge;
//    (*mgr_serie)[m] = new mgr_t();
    // index is assigned to loop index
    gauge.number = m;
    // name is assigned to filename_index
    std::stringstream tmp;
    int characteristic = floor(log10(static_cast<float>(nmgr)));
    tmp.width (characteristic + 1);
    tmp << setfill('0') << m;
    string name = time_series_file_name + tmp.str();
    tmp.str().clear();
    strcpy(gauge.name, name.c_str());
    name.clear();
    
    // location as read from lon/lat variables
    gauge.loc.lon = lon[m];
    gauge.loc.lat = lat[m];
    gauge.loc.depth = 0;
    
    // everything else is read from file global attributes
    strcpy(gauge.origine, time_series_origin.c_str());
    strcpy(gauge.validation, time_series_validation.c_str());
    gauge.loc.units = strdup(loc_units.c_str());
    gauge.mindex = 1;
    gauge.duree = sample[m];
    strcpy(gauge.debut, time_series_first_date_used.c_str());
    strcpy(gauge.fin, time_series_last_date_used.c_str());

    gauge.data = new mgr_data_t[nwave];

    // tidal constants et al.
    size_t  count = 0;
    for(size_t w = 0; w < nwave; w++){
      tidal_wave wave=reference_spectrum.wave(constituentname[w]);
//       gauge.data[w].code=wave.code;
      gauge.data[w].amp   = a[m * nwave + w];
      gauge.data[w].phi   = g[m * nwave + w];
      gauge.data[w].error = e[m * nwave + w];
      if(got_bg == true){
        gauge.data[w].bg_contamination_error = bg[m * nwave + w];
        gauge.data[w].error = bg[m * nwave + w];
        }
      strcpy(gauge.data[w].constituent.name, constituentname[w]);
      //gauge.data[w]->constituent =;
      //gauge.data[w]->code =;
      //gauge.data[w].index = 1;
      if( gauge.data[w].amp != fillValue){
        count++;
        }
      }
    gauge.nwave = count;
    mgr.push_back(gauge);
    }
  
  loc_units.clear();

  status = NC_NOERR + 13;
  status = nc_close(ncid);
  NC_TRAP_ERROR(wexit,status,1,"");
  
  delete [] lon;
  delete [] lat;
  delete [] sample;

  delete [] a;
  delete [] g;
  delete [] e;
  delete [] bg;


  time_series_file_name.clear();
  time_series_origin.clear();
  time_series_validation.clear();
  time_series_first_date_used.clear();
  time_series_last_date_used.clear();

  return(nmgr);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int mgr_loadRAY(const char *filename,vector<mgr_t> & mgr_serie)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* Few comments (LR)
     This program converts Richard's Ray TG data set into mgr_t**.
     Amplitude units are assumed to be cm (so scaled to be set in meters)
*/
{
  int j,k;
  int nwaves;
  int nitems;
  int nmgr;
  double lon, lat, amp, phi;
  double scale_factor = 0.01;
  char sdum[50];
  tidal_wave wave;

  FILE *infile=NULL;
  assert( (infile=fopen(filename,"r")) != NULL );

  string origin = string("R. Ray - ") + string(filename);	

  char *line = NULL;
  assert( (line = (char *) malloc(200 * sizeof(char))) != NULL );
 
  int nst=0;
  while( fgetc(infile) != EOF) {
    mgr_t mgr_station;
    // allocate memory and initialisation
    // first line: station informations... later (too ugly!)
    if(fgets(line, 200, infile)==0) break;
    strncpy(mgr_station.name, (char *) line, 20);
    mgr_station.name[20]=0;
    mgr_station.mindex = 1;
    mgr_station.number = nst++;
    strcpy(mgr_station.origine, (char *) origin.c_str());
    // 2nd line: coordinates and nwaves
    if(fgets(line,200,infile)==0) break;
    nitems=sscanf(line,"%lf %lf %d",&lat,&lon,&nwaves);
    if(nitems!=3) {
      printf("mgr_loadRAY: troubles at line %s \n",line);
      }
//    assert( (nitems=sscanf(line,"%lf %lf %d",&lat,&lon,&nwaves)) ==3);

    mgr_station.loc.lat = lat;
    mgr_station.loc.lon = lon;
    mgr_station.loc.depth=0.;
    mgr_station.loc.units = strdup("degrees");

    // next lines: harmonic constants
    mgr_station.data = new mgr_data_t[nwaves];
    for(k = 0, j = 0; k<nwaves;k++) {
      assert( (fgets(line,200,infile)) != NULL );
      assert( (nitems=sscanf(line,"%s %lf %lf",&sdum,&amp,&phi)) ==3);
      wave = wave_from_name(sdum);
      if(strcmp(wave.name, "wNUL") != 0) {
        sprintf(mgr_station.data[j].constituent.name,"%s",wave.name);
        mgr_station.data[j].amp = amp * scale_factor;
        mgr_station.data[j].phi = phi;
//         mgr_station.data[j].code = wave.code;
        //mgr_station.data[j].index = 1;// ???
        j++;
        }
      }
    mgr_station.nwave = j;
    nwaves = mgr_station.nwave;
    mgr_serie.push_back(mgr_station);
    }

  nmgr = mgr_serie.size();

  // free, close, ...
  assert( fclose(infile)==0);
  free(line);

  return(nmgr);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int mgr_loadATG_DATABASE(const char *filename,vector<mgr_t> & mgr_serie)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* Few comments (LR)
     This program converts Richard's Ray TG data set into mgr_t**.
     Amplitude units are assumed to be cm (so scaled to be set in meters)
*/
{
  int j,k;
  int num, nwaves;
  int nitems;
  int nmgr;
  double lon, lat, amp, phi;
  double scale_factor = 0.01;
  char sdum[50],instrument[128];
  tidal_wave wave;

  const char *waves[12]={"Q1","O1","P1","K1","N2","M2","S2","K2"};

  FILE *infile=NULL;
  assert( (infile=fopen(filename,"r")) != NULL );

  string origin = string("ATG_DATABASE - ") + string(filename);	

  char *line = NULL;
  assert( (line = (char *) malloc(200 * sizeof(char))) != NULL );
 
  
  for(k=0;k<12;k++) fgets(line,200,infile);
  
  int nst=0;
  while( fgetc(infile) != EOF) {
    mgr_t mgr_station;
    if(fgets(line,200,infile)==0) break;
    nitems=sscanf(line,"%d %42c %lf %lf %d",&num,mgr_station.name,&lat,&lon,sdum);
    if(nitems!=5) {
      printf("mgr_loadATG_DATABASE: troubles at line %s \n",line);
      }
    mgr_station.number = nst++;
//    strcpy(mgr_station.origine, (char *) origin.c_str());
    mgr_station.loc.lat = lat;
    mgr_station.loc.lon = lon;
    mgr_station.loc.depth=0.;
    mgr_station.loc.units = strdup("degrees");
/**----------------------------------------------------------------------------
    skeep observation type line */
    fgets(line,200,infile);
    nitems=sscanf(line,"%s",instrument);
/**----------------------------------------------------------------------------
    skeep observation type line */
    fgets(mgr_station.origine,200,infile);
    mgr_station.origine[strlen(mgr_station.origine)-2]=0;
/**----------------------------------------------------------------------------
    Q1,O1,P1,K1,N2,M2,S2,K2 */
    nwaves=8;
    mgr_station.data = new mgr_data_t[nwaves];
    for(k = 0, j = 0; k<nwaves;k++) {
      wave = wave_from_name(waves[k]);
      if(strcmp(wave.name, "wNUL") != 0) {
        sprintf(mgr_station.data[j].constituent.name,"%s",wave.name);
//         mgr_station.data[j].code = wave.code;
        //mgr_station.data[j].index = 1;// ???
        j++;
        }
      }
    for(k = 0, j = 0; k<nwaves;k++) {
      nitems=fscanf(infile,"%lf",&amp);
      if(nitems==0) amp=0;
      wave = wave_from_name(waves[k]);
      if(strcmp(wave.name, "wNUL") != 0) {
        mgr_station.data[j].amp = amp * scale_factor;
        j++;
        }
      }
    for(k = 0, j = 0; k<nwaves;k++) {
      nitems=fscanf(infile,"%lf",&phi);
      if(nitems==0) phi=0;
      wave = wave_from_name(waves[k]);
      if(strcmp(wave.name, "wNUL") != 0) {
        mgr_station.data[j].phi = phi;
        j++;
        }
      }
    mgr_station.nwave = j;
    nwaves = mgr_station.nwave;
    if(strstr(instrument,"BPR")!=0)
      mgr_serie.push_back(mgr_station);
/**----------------------------------------------------------------------------
    skeep empty line */
    fgets(line,200,infile);
    fgets(line,200,infile);
    }

  nmgr = mgr_serie.size();

  assert( fclose(infile)==0);
  free(line);

  return(nmgr);
}
