
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "tools-define.h"
#include "tools-structures.h"

#include "fe.h"
#include "map.h"
#include "archive.h"
#include "netcdf-proto.h"


   /* variable ids */
extern   int lon_id;
extern   int lat_id;
extern   int time_id;
extern   int h_id;
extern   int u_id;
extern   int v_id;
extern   int p_id;
extern   int ibd_id;
extern tidal_wave wS2,wS1,wSa;

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  grid_t get_zonegrid_fatal(char *zone)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{

  grid_t grid;

  if (strcmp(zone,"indien")==0) {
    /* Frame indien */
    map_set2Dgrid(&grid,+20,-60,+150,+20,1.,1.);
    }
  else if(strcmp(zone,"medsea")==0) {
    /* Frame medsea */
    map_set2Dgrid(&grid,-10,+30,+37.5,+47.5,.125,.125);
    }
  else if(strcmp(zone,"gibraltar")==0) {
    /* Zoom Gibraltar */
    map_set2Dgrid(&grid,-6.5,+35.5,-4.5,+36.5,.01,.01);
    }
  else if(strcmp(zone,"kerguelen")==0) {
    /* Frame kerguelen */
    map_set2Dgrid(&grid,45,-55,80,-35,.25,.25);
    }
  else if(strcmp(zone,"global")==0) {
    /* Frame global */
    map_set2Dgrid(&grid,0.,-90.,360.,90.,1.,1.);
    }
  else if(strcmp(zone,"global-0.5")==0) {
    /* Frame global */
    map_set2Dgrid(&grid,0.,-90.,360.,90.,0.5,0.5);
    }
  else if(strcmp(zone,"anglet")==0) {
    /* Frame global */
    map_set2Dgrid(&grid,-2.25,+43.25,-1.25,+44.25,0.025,0.025);
    }
  else if(strcmp(zone,"SWAtl1")==0) {
    /* Frame global */
    map_set2Dgrid(&grid,300.,-50.,335.,-30.,0.5,0.5);
    }
  else if(strcmp(zone,"caspian")==0) {
    /* Frame global */
    map_set2Dgrid(&grid,46.0,36.0,55.0,47.5,0.05,0.05);
    }
  else if(strcmp(zone,"sicily")==0) {
    /* Frame strait of Sicily */
    map_set2Dgrid(&grid,9.0,30.0,20.0,40.0,0.05,0.05);
    }
  else if(strcmp(zone,"albicocca")==0) {
    /* Frame global */
    map_set2Dgrid(&grid,5.0,40.0,12.0,45.0,0.02,0.02);
    }
  else if(strcmp(zone,"cataluna")==0) {
    /* Frame global */
    map_set2Dgrid(&grid,-0.5,38.5,6.5,44.0,0.02,0.02);
    }
  else if(strcmp(zone,"caledonie")==0) {
    /* Frame  */
    map_set2Dgrid(&grid,140,-35,185.,-5,.1,.1);
    }
  else if(strcmp(zone,"test")==0) {
    map_set2Dgrid(&grid,-67.0,+45.0,-66.25,+45.5,0.005,0.005);
    }
  else
    {
    __ERR_BASE_LINE__("exiting\n");exit(-1);
    }

  grid.modeH=0;
  
  grid.x=(double *) malloc(grid.nx*sizeof(double));
  grid.y=(double *) malloc(grid.ny*sizeof(double));

  return(grid);
}

  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  float *buf2d,mask=1.e+35;
  float *buf3d;

  double time,start_time,end_time;
  int day,month,year,first_day;
  int frame,ncid,id,input_id;
  int *element,m;
  size_t *start,*count;

  int i,k,minstep,maxstep,nndes,nlayers,stride=1,n,status;
  size_t size;
  char *h_root=NULL,*keyword,*zone=NULL,*fmt=NULL;
  char hfile[1024],*pathname=NULL,*outpath=NULL;
  char *gridfile=NULL;
  char rootname[1024]="\0",filename[1024]="\0",tidefile[1024]="\0",*s;
  grid_t grid;
  grid_t grid2d;
  date_t date,model_date,start_date,end_date,actual;
  meta_archive_t hinfo;
  int options[10]={0,0,0,1,0,0,0,0,0,0};
  int detide_s2=0,detide_s1=0,detide_sa=0;
  fcomplex **constants=NULL;
  double  *tide=NULL,cnestime;
  float *buffer[2][3]={{NULL,NULL,NULL},{NULL,NULL,NULL}};
  float *fv2d=NULL;
  float *fv3d=NULL;
  cdfvar_t variable[10];
  pocgrd_t ncgrid;

  fct_echo(argc,argv);
 
  n=1;
  while (n < argc)
    {
    keyword=strdup(argv[n]);
    switch (keyword[0])
      {
      case '-':
        switch (keyword[1])
        {
/*------------------------------------------------------------------------
        rootname for  model simulations sea state files */
        case 'r' :
          h_root= strdup(argv[n+1]);
          n++;
          n++;
          break;

/*------------------------------------------------------------------------
        pathname for model simulations files */
        case 'd' :
          pathname= strdup(argv[n+1]);
          n++;
          n++;
          break;

/*------------------------------------------------------------------------
        region of extraction [1] */
        case 'z' :
          zone= strdup(argv[n+1]);
          n++;
          n++;
          break;

/*------------------------------------------------------------------------
        grid of extraction   [2]*/
        case 'g' :
          gridfile= strdup(argv[n+1]);
          n++;
          n++;
          break;

/*------------------------------------------------------------------------
        starting date */
        case 's' :
          s= strdup(argv[n+1]);
          n++;
          n++;
          sscanf(s,"%d/%d",&start_date.month,&start_date.year);
          free(s);
          break;

/*------------------------------------------------------------------------
        ending date */
        case 'e' :
          s= strdup(argv[n+1]);
          n++;
          n++;
          sscanf(s,"%d/%d",&end_date.month,&end_date.year);
          free(s);
          break;
 
/*------------------------------------------------------------------------
        pathname for extraction files */
        case 'o' :
          outpath= strdup(argv[n+1]);
          n++;
          n++;
          break;

/*------------------------------------------------------------------------
        extraction file format */
        case 'f' :
          fmt= strdup(argv[n+1]);
          n++;
          n++;
          break;

/*------------------------------------------------------------------------
        extraction interval (2 = 1frame over 2)*/
        case 'i' :
          s= strdup(argv[n+1]);
          n++;
          n++;
          sscanf(s,"%d",&stride);
         free(s);
          break;

        default:
          __OUT_BASE_LINE__("unknown option %s\n",keyword);
          exit(-1);
        }
        break;

      default:
        __OUT_BASE_LINE__("unknown option %s\n",keyword);
        exit(-1);
        break;
      }
      free(keyword);
    }

  if(fmt==NULL) fmt= strdup("cdf");
  if((zone == NULL) && (gridfile == NULL)) goto error;

  if(zone != NULL) {
    grid=get_zonegrid_fatal(zone);
    status=map3d_completegridaxis(&grid);
    }
  else
    zone=strdup("grid");

  if(gridfile != NULL) {
    status=cdf_loadvargrid(gridfile,0,&grid2d);
    }
  else
    {
    grid2d= map_getgrid2d(grid);
    }

  if(h_root==NULL) {
    h_root= strdup("archive");
    printf("use <archive> as root name for sea state archive files\n");
    }

  if(pathname==NULL) {
    pathname= strdup(".");
    printf("use <.> as default path for archive files\n");
    }

 if(outpath==NULL) {
    outpath = strdup(".");
    printf("use <.> as default path for output files\n");
    }

  for (year=start_date.year;year<=end_date.year;year++) {
    for (month=1;month<=12;month++) {
      if((year==start_date.year) && (month < start_date.month)) continue;
      if((year==end_date.year)  &&  (month > end_date.month))  break;
      sprintf(hfile,"%s/%s-%4.4d.%2.2d.nc",pathname,h_root,year,month);

      status=cdf_archiveinfo(hfile,&hinfo);

      status=fe_list(&hinfo.mesh);
      element=(int *) fe_scan_elements(hinfo.mesh,grid2d,0);
      if (element == NULL) goto error;

      nndes=hinfo.mesh.nvtxs;
      nlayers=hinfo.mesh.nlayers;

      minstep=0;
      maxstep=hinfo.nframe;

      fv2d=(float *) malloc(nndes*hinfo.nframe*sizeof(float));
      if(fv2d == NULL) {
        printf("#memory allocation error for buffer N= %d \n",nndes);
        goto error;
        }

      fv3d=(float *) malloc(nndes*hinfo.mesh.nlayers*sizeof(float));
      if(fv3d == NULL) {
        printf("#memory allocation error for buffer N= %d \n",nndes);
        goto error;
        }

      grid.nz=hinfo.mesh.nlayers;

      buf2d=(float *) malloc(grid.nx*grid.ny*sizeof(float));
      if (buf2d == NULL) goto error;
      buf3d=(float *) malloc(grid.nx*grid.ny*grid.nz*sizeof(float));
      if (buf3d == NULL) goto error;

      status=nc_open(hfile,NC_NOWRITE,&input_id);

/*------------------------------------------------------------------------
      compute bathymetry*/
      grid.bottom=(double *) malloc(grid.nx*grid.ny*sizeof(double));
      status=nc_inq_varid(input_id,"depth",&id);
      if(status!=0) goto error;
      status=nc_get_var_float(input_id,id,fv2d);
      if(status!=0) goto error;
      status=fe_map(hinfo.mesh,fv2d,grid2d,element,buf2d,mask);
      for(m=0;m<grid.nx*grid.ny;m++) {
        grid.bottom[m]=(double) buf2d[m];
       }

/*------------------------------------------------------------------------
      compute sigma-levels*/
      grid.sigma=(double *) malloc(grid.nx*grid.ny*grid.nz*sizeof(double));
      for(k=0;k<nlayers;k++) {
        for(n=0;n<hinfo.mesh.nvtxs;n++)  fv2d[n]=hinfo.mesh.vertices[n].sigma[k];
        status=fe_map(hinfo.mesh,fv2d,grid2d,element,buf2d,mask);
        for(m=0;m<grid.nx*grid.ny;m++) grid.sigma[k*grid.nx*grid.ny+m]=buf2d[m];
       }

/*------------------------------------------------------------------------
      compute z-levels*/
      grid.z=(double *) malloc(grid.nx*grid.ny*grid.nz*sizeof(double));
      for(k=0;k<nlayers;k++) {
        for(m=0;m<grid.nx*grid.ny;m++)
          if(buf2d[m]!=mask) {
            grid.z[k*grid.nx*grid.ny+m]=grid.sigma[k*grid.nx*grid.ny+m]*grid.bottom[m];
            }
          else
            grid.z[k*grid.nx*grid.ny+m]=mask;
       }

/*------------------------------------------------------------------------
      compute mask*/
      grid.mask=(signed char *) malloc(grid.nx*grid.ny*sizeof(char));
/*
      for(n=0;n<hinfo.mesh.nvtxs;n++)  fv2d[m]=1.;
      status=fe_map(hinfo.mesh,fv2d,grid2d,element,buf2d,mask);
*/
      for(m=0;m<grid.nx*grid.ny;m++)
        if(buf2d[m]!=mask) grid.mask[m]=1;
        else  grid.mask[m]=0;

      printf("create the gridded file: %2.2d/%4d\n",month,year);
      sprintf(filename,"%s/%s-%4.4d.%2.2d.huv.nc",outpath,zone,year,month);

      status= poc_createfile(filename);
      status=poc_sphericalgrid_xyzt(filename,"","",grid,&ncgrid);

      printf("geometry done\n");
      mask=1.e+35;

/* *------------------------------------------------------------------------
      treat elevations */
      if(options[0]==1) {
        status=poc_standardvariable_xyt("surfelev",mask,"m",1.f, 0.f,"sea_surface_elevation","sea surface elevation","h",ncgrid,variable[0]);
        status=create_ncvariable(filename, &(variable[0]));
        status=nc_inq_varid(input_id,"surfelev",&id);
        if(status!=0) goto error;
        status=poc_getbuffer(input_id,id,fv2d);
        if(status!=0) goto error;
        for(frame=minstep;frame<maxstep;frame+=stride) {
          status=fe_map(hinfo.mesh,&(fv2d[nndes*frame]),grid2d,element,buf2d,mask);
          if(status!=0) goto error;
/*------------------------------------------------------------------------
          convert elevation in meters
          for(m=0;m<grid.nx*grid.ny;m++) if(buf2d[m]!=mask) buf2d[m]/=100.;*/
          status=write_ncfile2d_obsolete(filename,  grid2d, frame, variable[0].id, buf2d);
          if(status!=0) goto error;
          }
        variable[0].destroy();
        }

/*------------------------------------------------------------------------
      treat barotropic currents */
      if(options[1]==1) {
        status=poc_standardvariable_xyt("ubar",mask,"m/s",1.f, 0.f,"mean_eastward_velocity","mean_eastward_velocity","ubar",ncgrid,variable[0]);
        status=poc_standardvariable_xyt("vbar",mask,"m/s",1.f, 0.f,"mean_northward_velocity","mean_northward_velocity","vbar",ncgrid,variable[1]);
        status=create_ncvariable(filename, &(variable[0]));
        status=create_ncvariable(filename, &(variable[1]));
        status=nc_inq_varid(input_id,"ubar",&id);
        if(status!=0) goto error;
        status=nc_get_var_float(input_id,id,fv2d);
        if(status!=0) goto error;
        for(frame=minstep;frame<maxstep;frame+=stride) {
          status=fe_map(hinfo.mesh,&(fv2d[nndes*frame]),grid2d,element,buf2d,mask);
          if(status!=0) goto error;
/*------------------------------------------------------------------------
          convert currents in meters/s */
          for(m=0;m<grid.nx*grid.ny;m++) if(buf2d[m]!=mask) buf2d[m]/=100.;
          status=write_ncfile2d_obsolete(filename,  grid2d, frame, variable[0].id, buf2d);
          if(status!=0) goto error;
          }
        status=nc_inq_varid(input_id,"vbar",&id);
        if(status!=0) goto error;
        status=nc_get_var_float(input_id,id,fv2d);
        if(status!=0) goto error;
        for(frame=minstep;frame<maxstep;frame+=stride) {
          status=fe_map(hinfo.mesh,&(fv2d[nndes*frame]),grid2d,element,buf2d,mask);
          if(status!=0) goto error;
/*------------------------------------------------------------------------
          convert currents in meters/s */
          for(m=0;m<grid.nx*grid.ny;m++) if(buf2d[m]!=mask) buf2d[m]/=100.;
          status=write_ncfile2d_obsolete(filename,  grid2d, frame, variable[1].id, buf2d);
          if(status!=0) goto error;
          }
        variable[0].destroy();
        variable[1].destroy();
        }

/* *------------------------------------------------------------------------
      treat temperature */
      if(options[2]==1) {
        status=poc_standardvariable_xyzt("T",mask,"C",1., 0.,"sea_water_temperature","sea water temperature","T",ncgrid,variable[0]);
        status=create_ncvariable(filename, &(variable[0]));
        status=nc_inq_varid(input_id,"T",&id);
        if(status!=0) goto error;
        start=(size_t *) malloc(3*sizeof(size_t));
        count=(size_t *) malloc(3*sizeof(size_t));
        for(frame=minstep;frame<maxstep;frame+=stride) {
          printf("temperature %d\n",frame);
          start[0]=frame;
          start[1]=0;
          start[2]=0;
          count[0]=1;
          count[1]=nndes;
          count[2]=nlayers;
          status=fe_extract(input_id,id,start,count,fv3d);
          for(k=0;k<nlayers;k++) {
            for(m=0;m<nndes;m++) fv2d[m]=fv3d[nlayers*m+k];
            status=fe_map(hinfo.mesh,fv2d,grid2d,element,&(buf3d[grid.nx*grid.ny*k]),mask);
            if(status!=0) goto error;
            }
          status=write_ncfile3d_obsolete(filename, grid, frame, variable[0].id, buf3d);
          if(status!=0) goto error;
          }
        variable[0].destroy();
        }


/*------------------------------------------------------------------------
      treat salinity */
      if(options[3]==1) {
        status=poc_standardvariable_xyzt("S",mask,"PSU",1., 0.,"sea_water_salinity","sea water salinity","S",ncgrid,variable[0]);
        status=create_ncvariable(filename, &(variable[0]));
        status=nc_inq_varid(input_id,"S",&id);
        if(status!=0) goto error;
        start=(size_t *) malloc(3*sizeof(size_t));
        count=(size_t *) malloc(3*sizeof(size_t));
        for(frame=minstep;frame<maxstep;frame+=stride) {
          printf("salinity %d\n",frame);
          start[0]=frame;
          start[1]=0;
          start[2]=0;
          count[0]=1;
          count[1]=nndes;
          count[2]=nlayers;
          status=fe_extract(input_id,id,start,count,fv3d);
          for(k=0;k<nlayers;k++) {
            for(m=0;m<nndes;m++) fv2d[m]=fv3d[nlayers*m+k];
            status=fe_map(hinfo.mesh,fv2d,grid2d,element,&(buf3d[grid.nx*grid.ny*k]),mask);
            if(status!=0) goto error;
           }
          status=write_ncfile3d_obsolete(filename, grid, frame, variable[0].id, buf3d);
          if(status!=0) goto error;
          }
        variable[0].destroy();
      }


/*------------------------------------------------------------------------
      treat density */
      if(options[4]==1) {
        status=poc_standardvariable_xyzt("RHO",mask,"NONE",1., 0.,"sea_water_density","sea water density","RHO",ncgrid,variable[0]);
        status=create_ncvariable(filename, &(variable[0]));
        status=nc_inq_varid(input_id,"RHO",&id);
        if(status!=0) goto error;
        start=(size_t *) malloc(3*sizeof(size_t));
        count=(size_t *) malloc(3*sizeof(size_t));
        for(frame=minstep;frame<maxstep;frame+=stride) {
          printf("density %d\n",frame);
          start[0]=frame;
          start[1]=0;
          start[2]=0;
          count[0]=1;
          count[1]=nndes;
          count[2]=nlayers;
          status=fe_extract(input_id,id,start,count,fv3d);
          for(k=0;k<nlayers;k++) {
            for(m=0;m<nndes;m++) fv2d[m]=fv3d[nlayers*m+k];
            status=fe_map(hinfo.mesh,fv2d,grid2d,element,&(buf3d[grid.nx*grid.ny*k]),mask);
            if(status!=0) goto error;
            }
          status=write_ncfile3d_obsolete(filename, grid, frame, variable[0].id, buf3d);
          if(status!=0) goto error;
          }
        variable[0].destroy();
        }

/*------------------------------------------------------------------------
      treat U */
      if(options[5]==1) {
        status=poc_standardvariable_xyzt("U",mask,"m/s",1., 0.,"eastward_sea_water_velocity","eastward sea water velocity","U",ncgrid,variable[0]);
        status=create_ncvariable(filename, &(variable[0]));
        status=nc_inq_varid(input_id,"U",&id);
        if(status!=0) goto error;
        start=(size_t *) malloc(3*sizeof(size_t));
        count=(size_t *) malloc(3*sizeof(size_t));
        for(frame=minstep;frame<maxstep;frame+=stride) {
          printf("U %d\n",frame);
          start[0]=frame;
          start[1]=0;
          start[2]=0;
          count[0]=1;
          count[1]=nndes;
          count[2]=nlayers;
          status=fe_extract(input_id,id,start,count,fv3d);
          for(k=0;k<nlayers;k++) {
            for(m=0;m<nndes;m++) fv2d[m]=fv3d[nlayers*m+k];
            status=fe_map(hinfo.mesh,fv2d,grid2d,element,&(buf3d[grid.nx*grid.ny*k]),mask);
            if(status!=0) goto error;
            }
          status=write_ncfile3d_obsolete(filename, grid, frame, variable[0].id, buf3d);
          if(status!=0) goto error;
          }
        variable[0].destroy();
        }

      status=nc_close(input_id);

      archive_freeinfo(&hinfo);
      for(k=0;k<3;k++) {
        for(i=0;i<2;i++) if(buffer[i][k]!=NULL) free(buffer[i][k]);
        }
      }
    }

end: __OUT_BASE_LINE__("end of convert ... \n");

  exit(0);

error:

  __ERR_BASE_LINE__("exiting\n");exit(-1);
 
 
}
