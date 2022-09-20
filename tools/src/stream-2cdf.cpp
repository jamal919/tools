
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
  else if(strcmp(zone,"global-0.125")==0) {
    /* Frame global */
    map_set2Dgrid(&grid,0.,-90.,360.,90.,0.125,0.125);
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

  grid.nz=1;
  grid.z=NULL;

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
  char *meshfile=NULL,*keyword,*zone=NULL,*fmt=NULL;
  char hfile[1024],*pathname=NULL,*outpath=NULL;
  char *gridfile=NULL;
  char rootname[1024]="\0",filename[1024]="\0",tidefile[1024]="\0",*s;
  grid_t grid;

  fcomplex **constants=NULL;
  float *buffer[2][3]={{NULL,NULL,NULL},{NULL,NULL,NULL}};
  float *fv2d=NULL,*fv2dx=NULL,*fv2dy=NULL;
  float *fv3d=NULL;
  cdfvar_t variable[10];
  pocgrd_t ncgrid;

  mesh_t mesh;

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
        extraction file format */
        case 'm' :
          meshfile= strdup(argv[n+1]);
          n++;
          n++;
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
    status=cdf_loadvargrid(gridfile,0,&grid);
    }

  if(pathname==NULL) {
    pathname= strdup(".");
    printf("use <.> as default path for archive files\n");
    }
 
 if(outpath==NULL) {
    outpath = strdup(".");
    printf("use <.> as default path for output files\n");
    }

  if(meshfile != NULL) {
    status=fe_readmesh(meshfile,MESH_FILE_FORMAT_TRIGRID,&mesh);
    if(status !=0) goto error;
    status=fe_list(&mesh);
    if(status !=0) goto error;
    }
  else  {
    printf("no mesh file specified; abort...\n");
    goto error;
    }
 
  element=(int *) fe_scan_elements(mesh,grid,0);
  if (element == NULL) goto error;

  nndes=mesh.nvtxs;
 
  fv2d=(float *) malloc(nndes*sizeof(float));
  if(fv2d == NULL) goto error;
  fv2dx=(float *) malloc(nndes*sizeof(float));
  if(fv2dx == NULL) goto error;
  fv2dy=(float *) malloc(nndes*sizeof(float));
  if(fv2dy == NULL) goto error;

  buf2d=(float *) malloc(grid.nx*grid.ny*sizeof(float));
  if (buf2d == NULL) goto error;

/*------------------------------------------------------------------------
  compute mask*/
  grid.mask=(signed char *) malloc(grid.nx*grid.ny*sizeof(char));

  for(m=0;m<grid.nx*grid.ny;m++)
    if(buf2d[m]!=mask) grid.mask[m]=1;
    else  grid.mask[m]=0;

  sprintf(filename,"%s/residual-%s.nc",outpath,zone);

  mask=1.e+35;

  status= poc_createfile(filename);

  status=poc_sphericalgrid_xyzt(filename,"","",grid,&ncgrid);

/*------------------------------------------------------------------------
  treat mean elevation */
  poc_standardvariable_xy(& variable[0],"z0",mask,"m",1., 0.,"mean_elevation","mean elevation","z0",ncgrid);

  status=create_ncvariable(filename, &(variable[0]));
  if(status!=0) goto error;

  sprintf(hfile,"%s/%s",pathname,"ele.mean.s2r");

  status=quoddy_loadr1(hfile, nndes, fv2d);
  if(status!=0) {
    printf("cannot read %s\n",hfile);
    goto error;
    }

  status=fe_map(mesh,fv2d,grid,element,buf2d,mask);
  if(status!=0) goto error;

/*------------------------------------------------------------------------
  convert to meters */
  for(m=0;m<grid.nx*grid.ny;m++)  if(buf2d[m]!=mask) buf2d[m]/=100.;

  status= poc_write_xy(filename, grid,variable[0].id, buf2d);

    variable[0].destroy();

/*------------------------------------------------------------------------
  treat mean current u */
  poc_standardvariable_xy(& variable[0],"u",mask,"m/s",1., 0.,"eastward_mean_current","eastward mean current","u",ncgrid);

  status=create_ncvariable(filename, &(variable[0]));
  if(status!=0) goto error;

  poc_standardvariable_xy(& variable[1],"v",mask,"m/s",1., 0.,"northward_mean_current","northward mean current","v",ncgrid);

  status=create_ncvariable(filename, &(variable[1]));
  if(status!=0) goto error;

  sprintf(hfile,"%s/%s",pathname,"uv.mean.v2r");

  status=quoddy_loadr2(hfile, nndes, fv2dx, fv2dy);
  if(status!=0) {
    printf("cannot read %s\n",hfile);
    goto error;
    }

  status=fe_map(mesh,fv2dx,grid,element,buf2d,mask);
  if(status!=0) goto error;

/*------------------------------------------------------------------------
  convert to m/s */
  for(m=0;m<grid.nx*grid.ny;m++)  if(buf2d[m]!=mask) buf2d[m]/=100.;

  status= poc_write_xy(filename, grid,variable[0].id, buf2d);

  status=fe_map(mesh,fv2dy,grid,element,buf2d,mask);
  if(status!=0) goto error;

/*------------------------------------------------------------------------
  convert to m/s */
  for(m=0;m<grid.nx*grid.ny;m++)  if(buf2d[m]!=mask) buf2d[m]/=100.;

  status= poc_write_xy(filename, grid,variable[1].id, buf2d);

  variable[0].destroy();
  variable[1].destroy();

/*------------------------------------------------------------------------
  treat mean transport U */
  poc_standardvariable_xy(& variable[0],"U",mask,"m^2/s",1., 0.,"eastward_mean_transport","eastward mean transport","U",ncgrid);

  status=create_ncvariable(filename, &(variable[0]));
  if(status!=0) goto error;

  poc_standardvariable_xy(& variable[1],"V",mask,"m^2/s",1., 0.,"northward_mean_transport","northward mean transport","V",ncgrid);

  status=create_ncvariable(filename, &(variable[1]));
  if(status!=0) goto error;

  sprintf(hfile,"%s/%s",pathname,"T.mean.v2r");

  status=quoddy_loadr2(hfile, nndes, fv2dx, fv2dy);
  if(status!=0) {
    printf("cannot read %s\n",hfile);
    goto error;
    }

  status=fe_map(mesh,fv2dx,grid,element,buf2d,mask);
  if(status!=0) goto error;

  status= poc_write_xy(filename, grid,variable[0].id, buf2d);

  status=fe_map(mesh,fv2dy,grid,element,buf2d,mask);
  if(status!=0) goto error;

  status= poc_write_xy(filename, grid,variable[1].id, buf2d);

  variable[0].destroy();
  variable[1].destroy();

  sprintf(filename,"%s/streamf-%s.nc",outpath,zone);

  mask=1.e+35;

  status= poc_createfile(filename);

  status=poc_sphericalgrid_xyzt(filename,"","",grid,&ncgrid);

/*------------------------------------------------------------------------
  treat stream function */
  poc_standardvariable_xy(& variable[0],"sf",mask,"Sv",1., 0.,"stream_function","stream_function","sf",ncgrid);

  status=create_ncvariable(filename, &(variable[0]));
  if(status!=0) goto error;

  sprintf(hfile,"%s/%s",pathname,"stream.s2r");

  status=quoddy_loadr1(hfile, nndes, fv2d);
  if(status!=0) {
    printf("cannot read %s\n",hfile);
    goto error;
    }

  status=fe_map(mesh,fv2d,grid,element,buf2d,mask);
  if(status!=0) goto error;

  status= poc_write_xy(filename, grid,variable[0].id, buf2d);

  variable[0].destroy();


/*------------------------------------------------------------------------
  treat potential function */
  poc_standardvariable_xy(& variable[0],"pf",mask,"Sv",1., 0.,"potential_function","potential_function","sf",ncgrid);

  status=create_ncvariable(filename, &(variable[0]));
  if(status!=0) goto error;

  sprintf(hfile,"%s/%s",pathname,"potential.s2r");

  status=quoddy_loadr1(hfile, nndes, fv2d);
  if(status!=0) {
    printf("cannot read %s\n",hfile);
    goto error;
    }

  status=fe_map(mesh,fv2d,grid,element,buf2d,mask);
  if(status!=0) goto error;

  status= poc_write_xy(filename, grid,variable[0].id, buf2d);

  variable[0].destroy();


/*------------------------------------------------------------------------
  treat mean non-divergent transport U */
  poc_standardvariable_xy(& variable[0],"Us",mask,"m^2/s",1., 0.,"eastward_mean_non_divergent_transport","eastward mean non divergent transport","Us",ncgrid);

  status=create_ncvariable(filename, &(variable[0]));
  if(status!=0) goto error;

  poc_standardvariable_xy(& variable[1],"Vs",mask,"m^2/s",1., 0.,"northward_mean_non_divergent_transport","northward mean non divergent transport","Vs",ncgrid);

  status=create_ncvariable(filename, &(variable[1]));
  if(status!=0) goto error;

  sprintf(hfile,"%s/%s",pathname,"non-divergent.v2r");

  status=quoddy_loadr2(hfile, nndes, fv2dx, fv2dy);
  if(status!=0) {
    printf("cannot read %s\n",hfile);
    goto error;
    }

  status=fe_map(mesh,fv2dx,grid,element,buf2d,mask);
  if(status!=0) goto error;

  status= poc_write_xy(filename, grid,variable[0].id, buf2d);

  status=fe_map(mesh,fv2dy,grid,element,buf2d,mask);
  if(status!=0) goto error;

  status= poc_write_xy(filename, grid,variable[1].id, buf2d);

  variable[0].destroy();
  variable[1].destroy();

/*------------------------------------------------------------------------
  treat mean divergent transport U */
  poc_standardvariable_xy(& variable[0],"Up",mask,"m^2/s",1., 0.,"eastward_mean_divergent_transport","eastward mean divergent transport","Up",ncgrid);

  status=create_ncvariable(filename, &(variable[0]));
  if(status!=0) goto error;

  poc_standardvariable_xy(& variable[1],"Vp",mask,"m^2/s",1., 0.,"northward_mean_divergent_transport","northward mean divergent transport","Vp",ncgrid);

  status=create_ncvariable(filename, &(variable[1]));
  if(status!=0) goto error;

  sprintf(hfile,"%s/%s",pathname,"divergent.v2r");

  status=quoddy_loadr2(hfile, nndes, fv2dx, fv2dy);
  if(status!=0) {
    printf("cannot read %s\n",hfile);
    goto error;
    }

  status=fe_map(mesh,fv2dx,grid,element,buf2d,mask);
  if(status!=0) goto error;

  status= poc_write_xy(filename, grid,variable[0].id, buf2d);

  status=fe_map(mesh,fv2dy,grid,element,buf2d,mask);
  if(status!=0) goto error;

  status= poc_write_xy(filename, grid,variable[1].id, buf2d);

  variable[0].destroy();
  variable[1].destroy();

end: __OUT_BASE_LINE__("end of convert ... \n");

  exit(0);

error:

  __ERR_BASE_LINE__("exiting\n");exit(-1);
 
 
}
