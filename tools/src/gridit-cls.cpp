
/**************************************************************************

  T-UGO tools, 2006-2009

  Unstructured Ocean Grid initiative

Contributors:

  Florent Lyard      LEGOS/CNRS, Toulouse, France
  Yoann Le Bars      LEGOS, Toulouse, France (PhD)
  Laurent Roblou     LEGOS, Toulouse, France
  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
  Frédéric Dupont    Université de Laval à Québec, Canada

E-mail: florent.lyard@legos.obs-mip.fr

***************************************************************************/

#define MAIN_SOURCE

#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

using namespace std;


 
#include "tools-structures.h"
#include "constants.h"

#include "functions.h"
#include "map.h"
#include "fe.h"
#include "netcdf-proto.h"
#include "archive.h"
#include "geo.h"
#include "sym-io.h"

#include "rutin.h"     /*  rutin.h contains common utility routines  */


/*----------------------------------------------------------------------------*/

int cls_write_xy(char *filename, grid_t grid, int id, short *buffer)

/*----------------------------------------------------------------------------*/
{
  int  i,j,m,n,ncid;
  short *tmp;

  size_t start[2];
  size_t count[2];
  int status=0;

  start[0]=0;
  start[1]=0;
  count[0]=grid.nx;
  count[1]=grid.ny;

  tmp=new short[grid.nx*grid.ny];
  for (j=0;j<grid.ny;j++) {
    for (i=0;i<grid.nx;i++) {
      m=grid.nx*j+i;
      n=grid.ny*i+j;
      tmp[n]=buffer[m];
      }
    }

  status=nc_open(filename,NC_WRITE,&ncid);
  nc_check_error(status,__LINE__,__FILE__);
  if(status !=0) goto error;
  status=nc_put_vara_short(ncid,id,start,count,tmp);
  nc_check_error(status,__LINE__,__FILE__);
  if(status !=0) goto error;
  status=nc_close(ncid);
  nc_check_error(status,__LINE__,__FILE__);
  if(status !=0) goto error;

error:
  delete[] tmp;
  return(status);
}

/*----------------------------------------------------------------------------*/

int cls_write_xy(char *filename, grid_t grid, int id, float *buffer)

/*----------------------------------------------------------------------------*/
{
  int  i,j,m,n,ncid;
  float *tmp;

  size_t start[2];
  size_t count[2];
  int status=0;

  start[0]=0;
  start[1]=0;
  count[0]=grid.nx;
  count[1]=grid.ny;

  tmp=new float[grid.nx*grid.ny];
  for (j=0;j<grid.ny;j++) {
    for (i=0;i<grid.nx;i++) {
      m=grid.nx*j+i;
      n=grid.ny*i+j;
      tmp[n]=buffer[m];
    }
  }
  printf("VERIF : %lf %lf \n",tmp[100],tmp[1000]);
  status=nc_open(filename,NC_WRITE,&ncid);
  nc_check_error(status,__LINE__,__FILE__);
  if(status !=0) goto error;
  status=nc_put_vara_float(ncid,id,start,count,tmp);
  nc_check_error(status,__LINE__,__FILE__);
  if(status !=0) goto error;
  status=nc_close(ncid);
  nc_check_error(status,__LINE__,__FILE__);
  if(status !=0) goto error;

error:
  delete[] tmp;
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int cls_standardaxis_xy(const char *vname,const char *units,const char *lname, double mask, char **dname, cdfvar_t *variable)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,l,status;
  char *tmpname;

  variable->id=-1;
  variable->type=NC_DOUBLE;

  variable->ndim=1;

  variable->initatt(3);

  variable->name=strdup(vname);

  k=0;
  variable->dim[k].name=strdup(dname[k]);

  k=0;
  variable->att[k].type=NC_CHAR;
  variable->att[k].name=strdup("units");
  variable->att[k].data=strdup(units);
  variable->att[k].length=strlen( variable->att[k].data);

  k++;
  variable->att[k].type=NC_CHAR;
  variable->att[k].name=strdup("long_name");
  variable->att[k].data=strdup(lname);
  variable->att[k].length=strlen( variable->att[k].data);


  k++;
  variable->att[k].type=NC_DOUBLE;
  variable->att[k].name=strdup("_FillValue");
  variable->att[k].data=(char *) malloc(sizeof(double));
  for(l=0;l<sizeof(double);l++) variable->att[k].data[l]=(char) *((char *) &mask +l);
  variable->att[k].length=1;

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int cls_standardaxis_xy(const char *vname,const char *units,const char *lname, int mask,const char **dname, cdfvar_t *variable)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,l,status;
  char *tmpname;

  variable->id=-1;
  variable->type=NC_INT;

  variable->ndim=1;

  variable->initatt(3);

  variable->name=strdup(vname);

  k=0;
  variable->dim[k].name=strdup(dname[k]);

  k=0;
  variable->att[k].type=NC_CHAR;
  variable->att[k].name=strdup("units");
  variable->att[k].data=strdup(units);
  variable->att[k].length=strlen( variable->att[k].data);

  k++;
  variable->att[k].type=NC_CHAR;
  variable->att[k].name=strdup("long_name");
  variable->att[k].data=strdup(lname);
  variable->att[k].length=strlen( variable->att[k].data);


  k++;
  variable->att[k].type=NC_INT;
  variable->att[k].name=strdup("_FillValue");
  variable->att[k].data=(char *) malloc(sizeof(int));
  for(l=0;l<sizeof(int);l++) variable->att[k].data[l]=(char) *((char *) &mask +l);
  variable->att[k].length=1;

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

cdfvar_t cls_floatvariable_xy(const char *name,float mask,const char *units,float scale,const char *longname,pocgrd_t grid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  cdfvar_t variable;
  int k,l,status;
  char *tmpname;
  float jd;

  variable.id=-1;
  variable.type=NC_FLOAT;

  variable.initdim(2);

  variable.initatt(5);

  variable.name=strdup(name);

  k=0;
  variable.dim[k].name=strdup("NbLongitudes");

  k++;
  variable.dim[k].name=strdup("NbLatitudes");

  k=0;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("units");
  variable.att[k].data=strdup(units);
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("long_name");
  variable.att[k].data=strdup(longname);
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_FLOAT;
  variable.att[k].name=strdup("_FillValue");
  variable.att[k].data=(char *) malloc(sizeof(float));
  for(l=0;l<sizeof(float);l++) variable.att[k].data[l]=(char) *((char *) &mask +l);
  variable.att[k].length=1;
/*
  k++;
  variable.att[k].type=NC_FLOAT;
  variable.att[k].name=strdup("scale_factor");
  variable.att[k].data=(char *) malloc(sizeof(float));
  for(l=0;l<sizeof(float);l++) variable.att[k].data[l]=(char) *((char *) &scale +l);
  memcpy(variable.att[k].data, &scale,sizeof(float));
  variable.att[k].length=1;

  k++;
  variable.att[k].type=NC_FLOAT;
  variable.att[k].name=strdup("Date_CNES_JD");
  jd=-1;
  variable.att[k].data=(char *) malloc(sizeof(float));
  for(l=0;l<sizeof(float);l++) variable.att[k].data[l]=(char) *((char *) &jd +l);
  memcpy(variable.att[k].data, &scale,sizeof(float));
  variable.att[k].length=1;
*/
  return(variable);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

cdfvar_t cls_shortvariable_xy(const char *name,short mask,char *units,float scale,char *longname,pocgrd_t grid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  cdfvar_t variable;
  int k,l,status;
  char *tmpname;
  float jd;

  variable.id=-1;
  variable.type=NC_SHORT;

  variable.initdim(2);

  variable.initatt(6);

  variable.name=strdup(name);

  k=0;
  variable.dim[k].name=strdup("NbLongitudes");

  k++;
  variable.dim[k].name=strdup("NbLatitudes");

  k=0;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("units");
  variable.att[k].data=strdup(units);
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("long_name");
  variable.att[k].data=strdup(longname);
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_SHORT;
  variable.att[k].name=strdup("_FillValue");
  variable.att[k].data=(char *) malloc(sizeof(short));
  for(l=0;l<sizeof(short);l++) variable.att[k].data[l]=(char) *((char *) &mask +l);
  variable.att[k].length=1;

  k++;
  variable.att[k].type=NC_FLOAT;
  variable.att[k].name=strdup("scale_factor");
  variable.att[k].data=(char *) malloc(sizeof(float));
  for(l=0;l<sizeof(float);l++) variable.att[k].data[l]=(char) *((char *) &scale +l);
  memcpy(variable.att[k].data, &scale,sizeof(float));
  variable.att[k].length=1;

  k++;
  variable.att[k].type=NC_FLOAT;
  variable.att[k].name=strdup("Date_CNES_JD");
  jd=0;
  variable.att[k].data=(char *) malloc(sizeof(float));
  for(l=0;l<sizeof(float);l++) variable.att[k].data[l]=(char) *((char *) &jd +l);
  memcpy(variable.att[k].data, &jd,sizeof(float));
  variable.att[k].length=1;

  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("date");
  variable.att[k].data=strdup("1950-01-01 00:00:00.000000 UTC");
  variable.att[k].length=strlen( variable.att[k].data);

   return(variable);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int cls_sphericalgrid_xy(const char *input,const char *name, grid_t grid, pocgrd_t *cdfgrid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{

  int k,l,m,status;
  cdfvar_t variable[10];
  char     *dname[10];
  char     *tmp[10];
  double mask,buffer[2];
  int dimlgth[10];
  int ncid,stat,dim_id,imask;

  cdfgbl_t global;
  cdfgrid->id=-1;
  char *lonname, *latname, *levelsname;
  char *xname, *yname;

/*------------------------------------------------------------------------

  Defines a grid to which variable will be fully or partially connected.

  1) defines dimensions
     beware the unlimited dimension

  2) creates standard horizontal and vertical axis


  3) returns informations to be passed to variable definition and creation routines

*------------------------------------------------------------------------*/

  dname[0]=strdup("NbLongitudes");
  dname[1]=strdup("NbLatitudes");
  dname[2]=strdup("LatLon");

/*------------------------------------------------------------------------
  enter define mode */
  stat=nc_open(input, NC_WRITE, &ncid);
  nc_check_error(stat,__LINE__,__FILE__);
  stat=nc_redef(ncid);
  nc_check_error(stat,__LINE__,__FILE__);

/*------------------------------------------------------------------------
  create dimensions */
  dimlgth[0]=grid.nx;
  dimlgth[1]=grid.ny;
  for(k=0;k<2;k++) {
    stat=nc_def_dim(ncid,dname[k], dimlgth[k], &dim_id);
    nc_check_error(stat,__LINE__,__FILE__);
    }

  dimlgth[2]=2;
  stat=nc_def_dim(ncid,dname[k], dimlgth[k], &dim_id);
  nc_check_error(stat,__LINE__,__FILE__);

  stat=nc_enddef (ncid);
  nc_check_error(stat,__LINE__,__FILE__);
  stat=nc_close(ncid);
  nc_check_error(stat,__LINE__,__FILE__);

/*------------------------------------------------------------------------
  create axis variable */
  tmp[0]=dname[2];

  cdfgrid->x=new cdfvar_t;
  cdfgrid->lon=new cdfvar_t;
  cdfgrid->lat=new cdfvar_t;

  imask=2147483647;

  status=cls_standardaxis_xy("LatLon", (char *)("count"), (char *)("No sense but necessary for some automatic tools"),imask,tmp,cdfgrid->x);
  status=create_ncvariable(input, cdfgrid->x);

  mask=1.84467440737096e+19;

  status=cls_standardaxis_xy("LatLonMin", (char *)("degree"), (char *)("Latitude/Longitude of south/west corner"),mask,tmp,cdfgrid->lon);
  status=create_ncvariable(input, cdfgrid->lon);
  status=cls_standardaxis_xy("LatLonStep", (char *)("degree"), (char *)("latitude/longitude steps"),mask,tmp,cdfgrid->lat);
  status=create_ncvariable(input, cdfgrid->lat);

  stat=nc_open(input,NC_WRITE,&ncid);
  nc_check_error(stat,__LINE__,__FILE__);

  buffer[0]=grid.ymin;
  buffer[1]=grid.xmin;

  stat=nc_put_var_double(ncid,cdfgrid->lon->id,buffer);
  nc_check_error(stat,__LINE__,__FILE__);

  buffer[0]=grid.dy;
  buffer[1]=grid.dy;

  stat=nc_put_var_double(ncid,cdfgrid->lat->id,buffer);
  nc_check_error(stat,__LINE__,__FILE__);

  stat=nc_close(ncid);
  nc_check_error(stat,__LINE__,__FILE__);

/*------------------------------------------------------------------------
  get file information */
  status=cdf_globalinfo(input,&global,0);

  cdfgrid->name=strdup(name);
  cdfgrid->filename=strdup(input);

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int cls_createfile(char *filename,char *title)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int  ncid; /* netCDF id */
  char *name,*value=NULL;
  int stat,id, lgth;

  stat=nc_create(filename, NC_CLOBBER, &ncid);
  nc_check_error(stat,__LINE__,__FILE__);

  value=strdup("CF1.X, POC revision 11/2009");
  lgth=strlen(value);
  stat=nc_put_att_text(ncid, NC_GLOBAL, "Conventions", lgth,value);
  nc_check_error(stat,__LINE__,__FILE__);
  if(value!=NULL) free(value);
  value=NULL;

  value=strdup("GRID_DOTS");
  lgth=strlen(value);
  stat=nc_put_att_text(ncid, NC_GLOBAL, "FileType", lgth, value);
  nc_check_error(stat,__LINE__,__FILE__);
  if(value!=NULL) free(value);
  value=NULL;

  lgth=strlen(title);
  stat=nc_put_att_text(ncid, NC_GLOBAL,"title", lgth, title);
  nc_check_error(stat,__LINE__,__FILE__);
  if(value!=NULL) free(value);
  value=NULL;

  lgth=strlen(filename);
  stat=nc_put_att_text(ncid, NC_GLOBAL, "file_name", lgth, filename);
  nc_check_error(stat,__LINE__,__FILE__);

  value=strdup("POC/Noveltis NetCDF interface (CLS conventions) produced this file");
  lgth=strlen(value);
  stat=nc_put_att_text(ncid, NC_GLOBAL, "production", lgth, value);
  nc_check_error(stat,__LINE__,__FILE__);
  if(value!=NULL) free(value);
  value=NULL;

  time_t creation_time;
  stat=time(&creation_time);

  value=strdup(ctime(&creation_time));
  lgth=strlen(value);
  value[lgth-1]=0;
  lgth--;
  stat=nc_put_att_text(ncid, NC_GLOBAL, "creation_date", lgth, value);
  nc_check_error(stat,__LINE__,__FILE__);
  if(value!=NULL) free(value);
  value=NULL;

   /* leave define mode */
  stat=nc_enddef(ncid);
  nc_check_error(stat,__LINE__,__FILE__);
  stat=nc_close(ncid);
  nc_check_error(stat,__LINE__,__FILE__);

   return 0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double t0,dT,pulsation;
  double *serie[500],a1,p1,a2,p2,zr,zi,d;

  float x,y;
  float dummy,*rbuffer;
  float *rbufferx,*rbuffery;
  float spec=-9999,scale;
  short smask,*sbuf;
  float rmask,*rbuf;
  float *rbufx,*rbufy;

  int nndes,status,*elts,fmt,fefmt;
  int i,j,k,l,L,m,n;
  int nitems;
  FILE *file;
  FILE *out;
  char *keyword,*zone=NULL,*gridfile=NULL,*notebook=NULL,*title=NULL;
  char *meshfile=NULL, *depthfile=NULL,*output=NULL,*path=NULL;
  char *zgridfile=NULL,*ugridfile=NULL,*vgridfile=NULL;
  char hfile[1024],ufile[1024];
  char *wave[1024];
  grid_t grid,cgrid;
  fcomplex *cbuffer,*cbuf,cmask;
  fcomplex *cbufferx,*cbuffery,*cbufx,*cbufy;
  mesh_t mesh;
  int nwave=0,analysis,ncid;
  spectrum_t spectrum;
  geo_t projection;
  string cmd;
  cdfvar_t variable;
  pocgrd_t ncgrid;
  string ErrorMessage;

  fct_echo( argc, argv);

  n=1;
  while (n < argc) {
    keyword=strdup(argv[n]);
    switch (keyword[0]) {
      case '-':
        switch (keyword[1]) {
        case 'z' :
          zone= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'g' :
          gridfile= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'n' :
          notebook= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'h' :
          zgridfile= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'u' :
          ugridfile= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'v' :
          vgridfile= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'm' :
          meshfile= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'd' :
          depthfile= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'a' :
          sscanf(argv[n+1],"%d",&analysis);
          n++;
          n++;
          break;

        case 's' :
          sscanf(argv[n+1],"%f",&spec);
          n++;
          n++;
          break;

        case 'o' :
          output= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'p' :
          path= strdup(argv[n+1]);
          n++;
          n++;
          break;

        default:
          __OUT_BASE_LINE__("unknown option %s\n",keyword);
          exit(-1);
        }
        break;

      default:
          wave[nwave]= strdup(argv[n]);
//        printf("input wave=%s\n",wave[nwave]);
          spectrum.n=nwave+1;
          nwave++;
          n++;
        break;
      }
      free(keyword);
    }

  spectrum.waves=new tidal_wave[spectrum.n];
  for (k=0;k<nwave;k++) {
    strcpy (spectrum.waves[k].name,wave[k]);
    }

  rmask=spec;
  cmask=fcomplex(spec,spec);

  if((zone == NULL) && (gridfile == NULL) && (notebook == NULL) && (zgridfile == NULL)) goto error;

  if(output == NULL) output=new char[1024];


  if(output ==NULL){
    ErrorMessage=(string) "missing output filename";
    check_error(-1, ErrorMessage, __LINE__, __FILE__, 1);
    }
  if(path == NULL) path=strdup(".");

  if(zone != NULL) {
    grid=get_zonegrid(zone);
    status=map_completegridaxis_2(&grid);
    }

  if(gridfile != NULL) {
    status=cdf_loadvargrid (gridfile,0,&grid);
//    cmd=string("ncdum -h ")+string(gridfile)+string(" > gridit.cdl");
//    system(cmd.c_str());
    }

 if(notebook != NULL) {
    status=load_notebook(notebook, &cgrid, &grid, &projection);
    printf("%s (notebook file) processed\n",notebook);
    }

  if(meshfile != NULL) {
    status=fe_readmesh(meshfile,MESH_FILE_FORMAT_TRIGRID,&mesh);
    if(status !=0) goto error;
    status=fe_list(&mesh);
    if(status !=0) goto error;
    }
  else {
    ErrorMessage="no mesh file specified; abort...";
    check_error(-1, ErrorMessage, __LINE__, __FILE__, 1);
    goto error;
    }

  nndes=mesh.nvtxs;
  //fe_initaffine(&mesh);

  elts=fe_scan_elements(mesh,grid,0);
  if(elts==NULL) goto error;

  grid.nz=1;
  grid.z=NULL;

  cbuffer=new complex<float>[nndes];
  rbuffer=new float[nndes];

  cbuf  =new complex<float>[grid.nx*grid.ny];
  sbuf  =new short[grid.nx*grid.ny];
  rbuf  =new float[grid.nx*grid.ny];
  rbufx =new float[grid.nx*grid.ny];
  rbufy =new float[grid.nx*grid.ny];

//   if(depthfile!=NULL) {
//     status=quoddy_loadr1(depthfile, nndes, rbuffer);
//     if(status!=0) {
//       ErrorMessage="failure in reading depth file " + (string) depthfile;
//       check_error((const int) -1, ErrorMessage, __LINE__, __FILE__, 1);
//       }
//     status=fe_map(mesh,rbuffer,grid,elts,rbuf,rmask);
//     if(status!=0) goto error;
//     poc_standardvariable_xy(&variable,"bathymetry",rmask,"m",1., 0.,"bathymetry","bathymetry","bathymetry",ncgrid);
//     status=create_ncvariable(output, &variable);
//     status=poc_write_xy(output, grid, variable.id,rbuf);
//     status=free_ncvariable(&variable);
//     }

  title=strdup("T-UGO/SpEnOI tidal atlas");
  for (k=0;k<nwave;k++) {
    sprintf(output,"%s_ocean.nc",wave[k]);
    status=cls_createfile(output,title);
    status=cls_sphericalgrid_xy(output,"",grid,&ncgrid);
    sprintf(hfile,"%s/%s.ele.%2.2d.s2c",path,wave[k],analysis);
    if(status!=0) goto error;
    printf("gridit: treating %s wave from harmonic analysis #%d \n",wave[k],analysis);
    cbuf=(fcomplex *) malloc(grid.nx*grid.ny*sizeof(fcomplex));
    status=quoddy_loadc1(hfile, nndes, cbuffer);
    if(status!=0) {
      printf("cannot read %s\n",hfile);
      goto error;
      }
    status=fe_map(mesh,cbuffer,grid,elts,cbuf,cmask);
    for (i=0;i<grid.nx*grid.ny;i++) {
      x=real(cbuf[i]);
      y=imag(cbuf[i]);
      if(x!=rmask) {
/*------------------------------------------------------------------------
        convert fcomplex h in amplitude (meters) and phase lag (degrees)*/
        rbufx[i]=sqrt(x*x+y*y);
        rbufy[i]=atan2(-y,x)*r2d;
        }
      else {
        rbufx[i]=rmask;
        rbufy[i]=rmask;
        }
      }

    smask=32767;
    rmask= 1.844674e19;

//    scale=0.001;
    scale=0.01; // Passage des metres aux centimètres !
    for (i=0;i<grid.nx*grid.ny;i++) {
      if(rbufx[i]!=rmask) {
        rbuf[i]=rbufx[i]/scale;
        }
      else {
        rbuf[i]=rmask;
        }
      }

    scale=1.;
    variable=cls_floatvariable_xy("Grid_0001",rmask,"cm",scale,"ocean_tide_amplitude",ncgrid);
    status=create_ncvariable(output, &variable);
    status=cls_write_xy(output, grid, variable.id,rbuf);
    variable.destroy();

    scale=1.;
    for (i=0;i<grid.nx*grid.ny;i++) {
      if(rbufy[i]!=rmask) {
        rbuf[i]=rbufy[i]/scale;
        }
      else {
        rbuf[i]=smask;
        }
      }
    variable=cls_floatvariable_xy("Grid_0002",rmask,"degree",scale,"ocean_tide_phase_lag",ncgrid);
    status=create_ncvariable(output, &variable);
    status=cls_write_xy(output, grid, variable.id, rbuf);
    variable.destroy();
    }

  delete[] rbuffer;
  delete[] cbuffer;
  delete[] rbuf;
  delete[] rbufx;
  delete[] rbufy;
  delete[] sbuf;

end:
  free(elts);
  __ERR_BASE_LINE__("exiting\n");exit(0);
error:
  __OUT_BASE_LINE__("gridit: error detected, quit ... \n");
  exit(-1);
}
