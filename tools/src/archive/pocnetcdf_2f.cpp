
/***************************************************************************
 *   Copyright (C) 2006 by Cyril NGUYEN                                    *
 *   nguc@aeropc30                                                         *
  ***************************************************************************/


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <cstdlib>
#include "tools-structures.h"     
#include "poc-netcdf.h"     
#include "cfortran.h"

using namespace std;

   pocgrd_t t_ncgrid[20];
   grid_t t_grid[20];   
   cdfvar_t variable[1000];
 
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int cdf_put_global_att_float(char *filename,char *name,nc_type xtype, size_t len, float *fp)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int option,status, ncid;
  
  printf("cdf_put_global_att_float %s %d %s\n",name,len,filename);
  
  status = nc_open(filename, NC_WRITE, &ncid);
  ncredef(ncid);
  nc_check_error(status,__LINE__,__FILE__);
  if(status != NC_NOERR) goto error;
  
  status = nc_put_att_float(ncid,NC_GLOBAL,name,xtype,len,fp);
  nc_check_error(status,__LINE__,__FILE__);
  if(status != NC_NOERR) goto error;
  
  status = nc_close(ncid);
  nc_check_error(status,__LINE__,__FILE__);
  
  return(0);
error:
  printf("ERROR cdf_put_global_att_float\n");
  nc_check_error(status,__LINE__,__FILE__);
  return(-1);
  }


/*------------------------------------------------------------------------------------------------------------------*/
void poccreatefile(char *filename,int status) 
  {
    printf("Creation de %s\n",filename);
    status= poc_createfile(filename); 
  };

FCALLSCSUB2(poccreatefile,POCCREATEFILE,poccreatefile,STRING,INT)
/*-------------------------------------------------------------------------------------------------------------------*/


/*-------------------------------------------------------------------------------------------------------------------*/
  void pocsphericalgridxyzt(char *filename, int grid_id,char *grid_ref,int nx, int ny,int nz, 
                            signed char *mask,double *lon,double *lat, double *levels,int status)
  {
    float maskv=1.e+35;
    pocgrd_t ncgrid;
    grid_t grid;   

    int tt,i,j,k;

    /* define grid */	
    grid.nx=nx;
    grid.ny=ny;
    grid.nz=nz; 
    /* alloc mask */
    /*grid.mask=(signed char *) malloc(grid.nx*grid.ny*sizeof(char));*/
    /*grid.x=(double *) malloc(grid.nx*grid.ny*sizeof(double));*/
     /*grid.y=(double *) malloc(grid.nx*grid.ny*sizeof(double));  */
    /*grid.z=(double *) malloc(grid.nx*grid.ny*grid.nz*sizeof(double));*/

    grid.mask=mask;
    grid.x=lon;
    grid.y=lat;
    grid.z=levels;
    cout << grid_ref << endl;
    status=poc_sphericalgrid_xyzt(filename,grid_ref,"",grid,&ncgrid);
    t_grid[grid_id]=grid;
    t_ncgrid[grid_id]=ncgrid;
  }
FCALLSCSUB11(pocsphericalgridxyzt,PCSHERICALGRIDXYZT,pocsphericalgridxyzt,STRING,INT,STRING,INT,INT,INT,BYTEV,
             DOUBLEV,DOUBLEV,DOUBLEV,INT)

/*-------------------------------------------------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------------------------------------------------*/
  void pocsphericalgridxietast(char *filename, int grid_id,char *grid_ref,int nx, int ny,int nz, 
                            signed char *mask,double *lon,double *lat, double *levels,int status)
  {
    float maskv=1.e+35;
    pocgrd_t ncgrid;
    grid_t grid;   

    int tt,i,j,k;

    /* define grid */	
    grid.nx=nx;
    grid.ny=ny;
    grid.nz=nz; 
    /* alloc mask */
    /*grid.mask=(signed char *) malloc(grid.nx*grid.ny*sizeof(char));*/
    /*grid.x=(double *) malloc(grid.nx*grid.ny*sizeof(double));*/
     /*grid.y=(double *) malloc(grid.nx*grid.ny*sizeof(double));  */
    /*grid.z=(double *) malloc(grid.nx*grid.ny*grid.nz*sizeof(double));*/

    grid.mask=mask;
    grid.x=lon;
    grid.y=lat;
    grid.z=levels;
    cout << grid_ref << endl;
    status=poc_sphericalgrid_xietast(filename,grid_ref,grid,&ncgrid);
    t_grid[grid_id]=grid;
    t_ncgrid[grid_id]=ncgrid;
  }
FCALLSCSUB11(pocsphericalgridxietast,PCSHERICALGRIDXIETAST,pocsphericalgridxietast,STRING,INT,STRING,INT,INT,INT,BYTEV,
             DOUBLEV,DOUBLEV,DOUBLEV,INT)

/*-------------------------------------------------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------------------------------------------------*/
  void pocsphericalgridxyzwt(char *filename, int grid_id,char *grid_ref,int nx, int ny,int nz, 
                             signed char *mask,double *lon,double *lat, double *levels,int nw, int status)
  {
    float maskv=1.e+35;
    pocgrd_t ncgrid;
    grid_t grid;   

    int tt,i,j,k;

    /* define grid */	
    grid.nx=nx;
    grid.ny=ny;
    grid.nz=nz; 
    /* alloc mask */
    /*grid.mask=(signed char *) malloc(grid.nx*grid.ny*sizeof(char));*/
    /*grid.x=(double *) malloc(grid.nx*grid.ny*sizeof(double));*/
     /*grid.y=(double *) malloc(grid.nx*grid.ny*sizeof(double));  */
    /*grid.z=(double *) malloc(grid.nx*grid.ny*grid.nz*sizeof(double));*/

    grid.mask=mask;
    grid.x=lon;
    grid.y=lat;
    grid.z=levels;
    cout << grid_ref << endl;
    status=poc_sphericalgrid_xyzwt(filename,grid_ref,nw,grid,&ncgrid);
    t_grid[grid_id]=grid;
    t_ncgrid[grid_id]=ncgrid;
  }
FCALLSCSUB12(pocsphericalgridxyzwt,PCSHERICALGRIDXYZWT,pocsphericalgridxyzwt,STRING,INT,STRING,INT,INT,INT,BYTEV,
             DOUBLEV,DOUBLEV,DOUBLEV,INT,INT)

/*-------------------------------------------------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------------------------------------------------*/
  void pocsphericalgridxietaswt(char *filename, int grid_id,char *grid_ref,int nx, int ny,int nz, 
                             signed char *mask,double *lon,double *lat, double *levels,int nw, int status)
  {
    float maskv=1.e+35;
    pocgrd_t ncgrid;
    grid_t grid;   

    int tt,i,j,k;

    /* define grid */	
    grid.nx=nx;
    grid.ny=ny;
    grid.nz=nz; 
    /* alloc mask */
    /*grid.mask=(signed char *) malloc(grid.nx*grid.ny*sizeof(char));*/
    /*grid.x=(double *) malloc(grid.nx*grid.ny*sizeof(double));*/
     /*grid.y=(double *) malloc(grid.nx*grid.ny*sizeof(double));  */
    /*grid.z=(double *) malloc(grid.nx*grid.ny*grid.nz*sizeof(double));*/

    grid.mask=mask;
    grid.x=lon;
    grid.y=lat;
    grid.z=levels;
    cout << grid_ref << endl;
    status=poc_sphericalgrid_xietaswt(filename,grid_ref,nw,grid,&ncgrid);
    t_grid[grid_id]=grid;
    t_ncgrid[grid_id]=ncgrid;
  }
FCALLSCSUB12(pocsphericalgridxietaswt,PCSHERICALGRIDXIETASWT,pocsphericalgridxietaswt,STRING,INT,STRING,INT,INT,INT,BYTEV,
             DOUBLEV,DOUBLEV,DOUBLEV,INT,INT)

/*-------------------------------------------------------------------------------------------------------------------*/


/*-------------------------------------------------------------------------------------------------------------------*/
  void pocstandardvariablexy(char *file,float filv,int grid_id,int var_id,char *name,char *units,float scale, float offset,
                      char *standardname, char *longname, char *shortname, int status)
  {    
    char *gridfilename;
    /* test d'existence d'une grille correspondant a l'identifiant*/
    poc_standardvariable_xy(&variable[var_id],name,filv, units,scale, offset, standardname, longname,
                                               shortname,t_ncgrid[grid_id]);
    
    gridfilename = strdup(t_ncgrid[grid_id].filename);
    printf("gridfilename= %s \n",gridfilename);
    if (!strcmp(gridfilename,file)) {
      status=create_ncvariable(file, &(variable[var_id])); 
    } else {
      status=create_ncvariable_indfile(file, &(variable[var_id]),gridfilename); 
    }

  }

FCALLSCSUB12(pocstandardvariablexy,POCSTANDARDVARIABLEXY,pocstandardvariablexy,STRING,FLOAT,INT,INT,STRING,STRING,FLOAT,
             FLOAT,STRING,STRING,STRING,INT)
/*-------------------------------------------------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------------------------------------------------*/
  void pocstandardvariablexieta(char *file,float filv,int grid_id,int var_id,char *name,char *units,float scale, float offset,
                      char *standardname, char *longname, char *shortname, int status)
  {    
    char *gridfilename;
    /* test d'existence d'une grille correspondant a l'identifiant*/
    variable[var_id]=poc_standardvariable_xieta(name,filv, units,scale, offset, standardname, longname,
                                               shortname,t_ncgrid[grid_id]);
    
    gridfilename = strdup(t_ncgrid[grid_id].filename);
    printf("gridfilename= %s \n",gridfilename);
    if (!strcmp(gridfilename,file)) {
      status=create_ncvariable(file, &(variable[var_id])); 
    } else {
      status=create_ncvariable_indfile(file, &(variable[var_id]),gridfilename); 
    }

  }

FCALLSCSUB12(pocstandardvariablexieta,POCSTANDARDVARIABLEXIETA,pocstandardvariablexieta,STRING,FLOAT,INT,INT,STRING,STRING,FLOAT,
             FLOAT,STRING,STRING,STRING,INT)
/*-------------------------------------------------------------------------------------------------------------------*/


/*-------------------------------------------------------------------------------------------------------------------*/
  void pocstandardvariablexietaw(char *file,float filv,int grid_id,int var_id,char *name,char *units,float scale, float offset,
                    char *standardname, char *longname, char *shortname, int status)
  {    

    char *gridfilename;
    /* test d'existence d'une grille correspondant a l'identifiant*/
    variable[var_id]=poc_standardvariable_xietaw(name,filv, units,scale, offset, standardname, longname,
      shortname,t_ncgrid[grid_id]); 
    gridfilename = strdup(t_ncgrid[grid_id].filename);
    if (!strcmp(gridfilename,file)) {
      status=create_ncvariable(file, &(variable[var_id])); 
    } else {
      status=create_ncvariable_indfile(file, &(variable[var_id]),gridfilename); 
    }

  }

FCALLSCSUB12(pocstandardvariablexietaw,POCSTANDARDVARIABLEXIETAW,pocstandardvariablexietaw,STRING,FLOAT,INT,INT,STRING,STRING,FLOAT,
             FLOAT,STRING,STRING,STRING,INT)
/*-------------------------------------------------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------------------------------------------------*/
  void pocstandardvariablexyz(char *file,float filv,int grid_id,int var_id,char *name,char *units,float scale, float offset,
                    char *standardname, char *longname, char *shortname, int status)
  {
    char *gridfilename;
    

    /* test d'existence d'une grille correspondant a l'identifiant*/
    status=poc_standardvariable_xyz(name,filv, units,scale, offset, standardname, longname,
      shortname,t_ncgrid[grid_id],variable[var_id]); 
    gridfilename = strdup(t_ncgrid[grid_id].filename);
    if (!strcmp(gridfilename,file)) {
      status=create_ncvariable(file, &(variable[var_id])); 
    } else {
      status=create_ncvariable_indfile(file, &(variable[var_id]),gridfilename); 
    }

  }

FCALLSCSUB12(pocstandardvariablexyz,POCSTANDARDVARIABLEXYZ,pocstandardvariablexyz,STRING,FLOAT,INT,INT,STRING,STRING,FLOAT,
             FLOAT,STRING,STRING,STRING,INT)
/*-------------------------------------------------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------------------------------------------------*/
  void pocstandardvariablexietas(char *file,float filv,int grid_id,int var_id,char *name,char *units,float scale, float offset,
                    char *standardname, char *longname, char *shortname, int status)
  {
    char *gridfilename;
    

    /* test d'existence d'une grille correspondant a l'identifiant*/
    variable[var_id]=poc_standardvariable_xietas(name,filv, units,scale, offset, standardname, longname,
      shortname,t_ncgrid[grid_id]); 
    gridfilename = strdup(t_ncgrid[grid_id].filename);
    if (!strcmp(gridfilename,file)) {
      status=create_ncvariable(file, &(variable[var_id])); 
    } else {
      status=create_ncvariable_indfile(file, &(variable[var_id]),gridfilename); 
    }

  }

FCALLSCSUB12(pocstandardvariablexietas,POCSTANDARDVARIABLEXIETAS,pocstandardvariablexietas,STRING,FLOAT,INT,INT,STRING,STRING,FLOAT,
             FLOAT,STRING,STRING,STRING,INT)
/*-------------------------------------------------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------------------------------------------------*/
void pocstandardvariablexietasw(char *file,float filv,int grid_id,int var_id,char *name,char *units,float scale, float offset,
                    char *standardname, char *longname, char *shortname, int status)
  {
    char *gridfilename;
    

    /* test d'existence d'une grille correspondant a l'identifiant*/
    variable[var_id]=poc_standardvariable_xietasw(name,filv, units,scale, offset, standardname, longname,
      shortname,t_ncgrid[grid_id]); 
    gridfilename = strdup(t_ncgrid[grid_id].filename);
    if (!strcmp(gridfilename,file)) {
      status=create_ncvariable(file, &(variable[var_id])); 
    } else {
      status=create_ncvariable_indfile(file, &(variable[var_id]),gridfilename); 
    }

  }

FCALLSCSUB12(pocstandardvariablexietasw,POCSTANDARDVARIABLEXIETASW,pocstandardvariablexietasw,STRING,FLOAT,INT,INT,STRING,STRING,FLOAT,
             FLOAT,STRING,STRING,STRING,INT)
/*-------------------------------------------------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------------------------------------------------*/
void pocstandardvariablexywt(char *file,float filv,int grid_id,int var_id,char *name,char *units,float scale, float offset,
                    char *standardname, char *longname, char *shortname, int status)
  {
    char *gridfilename;
    

    /* test d'existence d'une grille correspondant a l'identifiant*/
    variable[var_id]=poc_standardvariable_xywt(name,filv, units,scale, offset, standardname, longname,
      shortname,t_ncgrid[grid_id]); 
    gridfilename = strdup(t_ncgrid[grid_id].filename);
    if (!strcmp(gridfilename,file)) {
      status=create_ncvariable(file, &(variable[var_id])); 
    } else {
      status=create_ncvariable_indfile(file, &(variable[var_id]),gridfilename); 
    }
  }

FCALLSCSUB12(pocstandardvariablexywt,POCSTANDARDVARIABLEXYWT,pocstandardvariablexywt,STRING,FLOAT,INT,INT,STRING,STRING,FLOAT,
             FLOAT,STRING,STRING,STRING,INT)
/*-------------------------------------------------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------------------------------------------------*/
void pocstandardvariablexietawt(char *file,float filv,int grid_id,int var_id,char *name,char *units,float scale, float offset,
                    char *standardname, char *longname, char *shortname, int status)
  {
    char *gridfilename;
    

    /* test d'existence d'une grille correspondant a l'identifiant*/
    variable[var_id]=poc_standardvariable_xietawt(name,filv, units,scale, offset, standardname, longname,
      shortname,t_ncgrid[grid_id]); 
    gridfilename = strdup(t_ncgrid[grid_id].filename);
    if (!strcmp(gridfilename,file)) {
      status=create_ncvariable(file, &(variable[var_id])); 
    } else {
      status=create_ncvariable_indfile(file, &(variable[var_id]),gridfilename); 
    }
  }

FCALLSCSUB12(pocstandardvariablexietawt,POCSTANDARDVARIABLEXIETAWT,pocstandardvariablexietawt,STRING,FLOAT,INT,INT,STRING,STRING,FLOAT,
             FLOAT,STRING,STRING,STRING,INT)
/*-------------------------------------------------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------------------------------------------------*/
void pocstandardvariablexyzt(char *file,float filv,int grid_id,int var_id,char *name,char *units,float scale, float offset,
                      char *standardname, char *longname, char *shortname, int status)
  {
    char *gridfilename;
    

    /* test d'existence d'une grille correspondant a l'identifiant*/
    status=poc_standardvariable_xyzt(name,filv, units,scale, offset, standardname, longname,
                                               shortname,t_ncgrid[grid_id], variable[var_id]);
    gridfilename = strdup(t_ncgrid[grid_id].filename);
    if (!strcmp(gridfilename,file)) {
      status=create_ncvariable(file, &(variable[var_id])); 
    } else {
      status=create_ncvariable_indfile(file, &(variable[var_id]),gridfilename); 
    }
  } 

FCALLSCSUB12(pocstandardvariablexyzt,POCSTANDARDVARIABLEXYZT,pocstandardvariablexyzt,STRING,FLOAT,INT,INT,STRING,STRING,FLOAT,
             FLOAT,STRING,STRING,STRING,INT)
/*-------------------------------------------------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------------------------------------------------*/
void pocstandardvariablexietast(char *file,float filv,int grid_id,int var_id,char *name,char *units,float scale, float offset,
                      char *standardname, char *longname, char *shortname, int status)
  {
    char *gridfilename;
    

    /* test d'existence d'une grille correspondant a l'identifiant*/
    variable[var_id]=poc_standardvariable_xietast(name,filv, units,scale, offset, standardname, longname,
                                               shortname,t_ncgrid[grid_id]);
    gridfilename = strdup(t_ncgrid[grid_id].filename);
    if (!strcmp(gridfilename,file)) {
      status=create_ncvariable(file, &(variable[var_id])); 
    } else {
      status=create_ncvariable_indfile(file, &(variable[var_id]),gridfilename); 
    }
  } 

FCALLSCSUB12(pocstandardvariablexietast,POCSTANDARDVARIABLEXIETAST,pocstandardvariablexietast,STRING,FLOAT,INT,INT,STRING,STRING,FLOAT,
             FLOAT,STRING,STRING,STRING,INT)
/*-------------------------------------------------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------------------------------------------------*/
  void pocstandardtime(char *file,int var_id,char *vname,char *dname,char *units,int year,int month,int day,float second, int status)
    {
    date_t origin;
    
    origin.year=year;
    origin.month=month;
    origin.day=day;
    origin.second=second;
    
    variable[var_id]=poc_standardtime(vname,dname,units,origin);
    status=create_ncvariable(file, &(variable[var_id]));      
    }
FCALLSCSUB10(pocstandardtime,POCSTANDARDTIME,pocstandardtime,STRING,INT,STRING,STRING,STRING,INT,INT,INT,FLOAT, INT)
/*-------------------------------------------------------------------------------------------------------------------*/


/*-------------------------------------------------------------------------------------------------------------------*/
  void pocstandardvariablexyzwt(char *file,float filv,int grid_id,int var_id,char *name,char *units,float scale, float offset,
                      char *standardname, char *longname, char *shortname, int status)
  {
    char *gridfilename;
    

    /* test d'existence d'une grille correspondant a l'identifiant*/
    variable[var_id]=poc_standardvariable_xyzwt(name,filv, units,scale, offset, standardname, longname,
                                               shortname,t_ncgrid[grid_id]);
    gridfilename = strdup(t_ncgrid[grid_id].filename);
    if (!strcmp(gridfilename,file)) {
      status=create_ncvariable(file, &(variable[var_id])); 
    } else {
      status=create_ncvariable_indfile(file, &(variable[var_id]),gridfilename); 
    }

  } 

FCALLSCSUB12(pocstandardvariablexyzwt,POCSTANDARDVARIABLEXYZWT,pocstandardvariablexyzwt,STRING,FLOAT,INT,INT,STRING,STRING,FLOAT,
             FLOAT,STRING,STRING,STRING,INT)
/*-------------------------------------------------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------------------------------------------------*/
  void pocstandardvariablexietaswt(char *file,float filv,int grid_id,int var_id,char *name,char *units,float scale, float offset,
                      char *standardname, char *longname, char *shortname, int status)
  {
    char *gridfilename;
    

    /* test d'existence d'une grille correspondant a l'identifiant*/
    variable[var_id]=poc_standardvariable_xietaswt(name,filv, units,scale, offset, standardname, longname,
                                               shortname,t_ncgrid[grid_id]);
    gridfilename = strdup(t_ncgrid[grid_id].filename);
    if (!strcmp(gridfilename,file)) {
      status=create_ncvariable(file, &(variable[var_id])); 
    } else {
      status=create_ncvariable_indfile(file, &(variable[var_id]),gridfilename); 
    }

  } 

FCALLSCSUB12(pocstandardvariablexietaswt,POCSTANDARDVARIABLEXIETASWT,pocstandardvariablexietaswt,STRING,FLOAT,INT,INT,STRING,STRING,FLOAT,
             FLOAT,STRING,STRING,STRING,INT)
/*-------------------------------------------------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------------------------------------------------*/
  void pocwritevar3dframe(char *file,int grid_id,int var_id,int frame,float *buf3d, int status)
{
  status=poc_write_xyzt(file, t_grid[grid_id], frame, variable[var_id].id, buf3d);
}
FCALLSCSUB6(pocwritevar3dframe,POCWRITEVAR3DFRAME,pocwritevar3dframe,STRING,INT,INT,INT,FLOATV, INT)
/*-------------------------------------------------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------------------------------------------------*/
  void pocwritevar3d(char *file,int grid_id,int var_id,int ndim, int *start,float *buf3d, int status)
{
  status=write_ncfile3d(file, t_grid[grid_id], ndim, start, variable[var_id].id, buf3d);
}
FCALLSCSUB7(pocwritevar3d,POCWRITEVAR3D,pocwritevar3d,STRING,INT,INT,INT,INTV,FLOATV, INT)
/*-------------------------------------------------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------------------------------------------------*/
  void pocstandardvariablexyt(char *file,float filv,int grid_id,int var_id,char *name,char *units,float scale, float offset,
                    char *standardname, char *longname, char *shortname, int status)
{
  char *gridfilename;

  /* test d'existence d'une grille correspondant a l'identifiant*/
  variable[var_id]=poc_standardvariable_xyt(name,filv, units,scale, offset, standardname, longname,
                                               shortname,t_ncgrid[grid_id]);
  gridfilename = strdup(t_ncgrid[grid_id].filename);
  if (!strcmp(gridfilename,file)) {
    status=create_ncvariable(file, &(variable[var_id])); 
    }
  else {
    status=create_ncvariable_indfile(file, &(variable[var_id]),gridfilename); 
    }
}

FCALLSCSUB12(pocstandardvariablexyt,POCSTANDARDVARIABLEXYT,pocstandardvariablexyt,STRING,FLOAT,INT,INT,STRING,STRING,FLOAT,
             FLOAT,STRING,STRING,STRING,INT)
/*-------------------------------------------------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------------------------------------------------*/
  void pocstandardvariablexietat(char *file,float filv,int grid_id,int var_id,char *name,char *units,float scale, float offset,
                    char *standardname, char *longname, char *shortname, int status)
  {
    char *gridfilename;
    

    /* test d'existence d'une grille correspondant a l'identifiant*/
    variable[var_id]=poc_standardvariable_xietat(name,filv, units,scale, offset, standardname, longname,
                                               shortname,t_ncgrid[grid_id]);
    gridfilename = strdup(t_ncgrid[grid_id].filename);
    if (!strcmp(gridfilename,file)) {
      status=create_ncvariable(file, &(variable[var_id])); 
    } else {
      status=create_ncvariable_indfile(file, &(variable[var_id]),gridfilename); 
    }

  }

FCALLSCSUB12(pocstandardvariablexietat,POCSTANDARDVARIABLEXIETAT,pocstandardvariablexietat,STRING,FLOAT,INT,INT,STRING,STRING,FLOAT,
             FLOAT,STRING,STRING,STRING,INT)
/*-------------------------------------------------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------------------------------------------------*/
  void pocwritevar2dframe(char *file,int grid_id,int var_id,int frame,float *buf2d, int status)
    {
      
          status=write_ncfile2dframe(file, t_grid[grid_id], frame, variable[var_id].id, buf2d);
    }
FCALLSCSUB6(pocwritevar2dframe,POCWRITEVAR2DFRAME,pocwritevar2dframe,STRING,INT,INT,INT,FLOATV, INT)
/*-------------------------------------------------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------------------------------------------------*/
  void pocwritevar2d(char *file,int grid_id,int var_id,int ndim, int *start,float *buf3d, int status)
    {
      status=write_ncfile2d(file, t_grid[grid_id], ndim, start, variable[var_id].id, buf3d);
    }
FCALLSCSUB7(pocwritevar2d,POCWRITEVAR2D,pocwritevar2d,STRING,INT,INT,INT,INTV,FLOATV, INT)
/*-------------------------------------------------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------------------------------------------------*/
  void pocwritexy(char *file,int grid_id,int var_id,float *buf2d, int status)
    {
      status=poc_write_xy(file,t_grid[grid_id],variable[var_id].id, buf2d);
    }
FCALLSCSUB5(pocwritexy,POCWRITEXY,pocwritexy,STRING,INT,INT,FLOATV, INT)
/*-------------------------------------------------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------------------------------------------------*/
  void pocwritetime(char *filename, int frame, int id, double time, int status)
    {
      status=poc_writetime(filename,frame,variable[id].id,time); 
    }
FCALLSCSUB5(pocwritetime,POCWRITETIME,pocwritetime,STRING,INT,INT,DOUBLE, INT)
/*-------------------------------------------------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------------------------------------------------*/
/*   GLOBAL ATTRIBUTES */
/*-------------------------------------------------------------------------------------------------------------------*/
  void cdfputglobalattfloat(char *filename,char *name,
        int nbelt, float *fp,int status)
    {
    size_t len;
    len = nbelt;    

    /*printf("appel a cdf_put_global_att_float \n"); */
    status = cdf_put_global_att_float(filename, name, NC_FLOAT, len, fp);
    } 
FCALLSCSUB5(cdfputglobalattfloat,CDFPUTGLOBALATTFLOAT,cdfputglobalattfloat,STRING,STRING,INT,FLOATV,INT)

/*-------------------------------------------------------------------------------------------------------------------*/
  void pocstandardvariablexyw(char *file,float filv,int grid_id,int var_id,char *name,char *units,float scale, float offset,
                    char *standardname, char *longname, char *shortname, int status)
  {    

    char *gridfilename;
    /* test d'existence d'une grille correspondant a l'identifiant*/
    variable[var_id]=poc_standardvariable_xyw(name,filv, units,scale, offset, standardname, longname,
      shortname,t_ncgrid[grid_id]); 
    gridfilename = strdup(t_ncgrid[grid_id].filename);
    if (!strcmp(gridfilename,file)) {
      status=create_ncvariable(file, &(variable[var_id])); 
    } else {
      status=create_ncvariable_indfile(file, &(variable[var_id]),gridfilename); 
    }

  }

FCALLSCSUB12(pocstandardvariablexyw,POCSTANDARDVARIABLEXYW,pocstandardvariablexyw,STRING,FLOAT,INT,INT,STRING,STRING,FLOAT,
             FLOAT,STRING,STRING,STRING,INT)
/*-------------------------------------------------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------------------------------------------------*/
void pocstandardvariablexyzw(char *file,float filv,int grid_id,int var_id,char *name,char *units,float scale, float offset,
                    char *standardname, char *longname, char *shortname, int status)
  {
    char *gridfilename;
    

    /* test d'existence d'une grille correspondant a l'identifiant*/
    variable[var_id]=poc_standardvariable_xyzw(name,filv, units,scale, offset, standardname, longname,
      shortname,t_ncgrid[grid_id]); 
    gridfilename = strdup(t_ncgrid[grid_id].filename);
    if (!strcmp(gridfilename,file)) {
      status=create_ncvariable(file, &(variable[var_id])); 
    } else {
      status=create_ncvariable_indfile(file, &(variable[var_id]),gridfilename); 
    }

  }

FCALLSCSUB12(pocstandardvariablexyzw,POCSTANDARDVARIABLEXYZW,pocstandardvariablexyzw,STRING,FLOAT,INT,INT,STRING,STRING,FLOAT,
             FLOAT,STRING,STRING,STRING,INT)
/*-------------------------------------------------------------------------------------------------------------------*/
