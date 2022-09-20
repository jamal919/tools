
/*******************************************************************************

  T-UGO tools, 2006-2009

  Unstructured Ocean Grid initiative

Contributors:

  Florent Lyard      LEGOS/CNRS, Toulouse, France
  Yoann Le Bars      LEGOS, Toulouse, France (PhD)
  Laurent Roblou     LEGOS, Toulouse, France
  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
  Frédéric Dupont    Université de Laval à Québec, Canada

E-mail: florent.lyard@legos.obs-mip.fr

*******************************************************************************/

/**-************************************************************************

  2D interpolation from structured to unstructured grid

*******************************************************************************/

#define MAIN_SOURCE

#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "tools-structures.h"

#include "fe.h"
#include "bmg.h"
#include "netcdf-proto.h"
#include "xyz.h"
#include "grd.h"
#include "map.h"
#include "geo.h"
#include "archive.h"
#include "ascii.h"
#include "functions.h"

#define XYZ 0
#define YXZ 1

#define FMT_CST 0
#define FMT_XYZ 1
#define FMT_BMG 2
#define FMT_CDF 3
#define FMT_ASC 4
#define FMT_GRD 5


#define ROW    0
#define COLUMN 1
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int cross_reference(mesh_t mesh, mesh_t node, int *translation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int   m,n;
  double x1,y1,x2,y2,d,*r;

  if(mesh.nvtxs!=node.nvtxs) goto error;

  r=new double[mesh.nvtxs];

  for(n=0;n<mesh.nvtxs;n++) {
    translation[n]=-1;
    }

  for(n=0;n<mesh.nvtxs;n++) {
    x1=mesh.vertices[n].lon;
    y1=mesh.vertices[n].lat;
    r[n]=1.e+10;
    for(m=0;m<node.nvtxs;m++) {
      x2=node.vertices[m].lon;
      y2=node.vertices[m].lat;
      d=(x2-x1)*(x2-x1)+(y2-y1)*(y2-y1);
      if(d<r[n]) {
        r[n]=d;
        translation[n]=m;
        }
      }
    }

  for(n=0;n<mesh.nvtxs;n++) {
    if(r[n]>1.e-03) goto error;
    }

  for(n=0;n<mesh.nvtxs;n++) {
    if(translation[n]==-1) goto error;
    }

  delete[] r;
  return (0);

 error:
  return (-1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int surfref_saver1(char *name, mesh_t mesh, float *buffer, char **comment)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int   n,nitems;
  FILE *in;

  in=fopen(name,"w");
  if(in == NULL) {
    printf("unable to open %s\n",name);
    return (-1);
    }

  nitems=fprintf(in,"%s\n",comment[0]);
  if(nitems ==0) goto error;
  nitems=fprintf(in,"%s\n",comment[1]);
  if(nitems ==0) goto error;

  for(n=0;n<mesh.nvtxs;n++) {
    nitems=fprintf(in,"%6d %12.6lf %12.6lf %12.6f\n",n+1,mesh.vertices[n].lon, mesh.vertices[n].lat, buffer[n]);
    if(nitems ==0) goto error;
    }

  fclose(in);
  return (0);

 error:
  fclose(in);
  return (-1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int surfref_loadr1(char *name, mesh_t mesh, float *buffer, char **comment)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int   n,k,nitems;
//   char c;
  FILE *in;
  double x,y;

  in=fopen(name,"r");
  if(in == NULL) {
    printf("unable to open %s\n",name);
    return (-1);
    }

  fgets(comment[0],1024,in);
  comment[0][strlen(comment[0])-1]=0;
  fgets(comment[1],1024,in);
  comment[1][strlen(comment[1])-1]=0;

  for(n=0;n<mesh.nvtxs;n++) {
    nitems=fscanf(in,"%d %lf %lf %f\n",&k,&x,&y, &(buffer[n]));
//    do { c=fgetc(in);   }  while ((c != '\n') && !feof(in));
    if(fabs(x-mesh.vertices[n].lon)>1.e-04) {
      printf("position mismatch %lf %lf\n",x,mesh.vertices[n].lon);
      goto error;
      }
    if(fabs(y-mesh.vertices[n].lat)>1.e-04) {
      printf("position mismatch %lf %lf\n",y,mesh.vertices[n].lat);
      goto error;
      }
    if(nitems ==0) goto error;
    }

  fclose(in);
  return (0);

 error:
  fclose(in);
  return (-1);
}



/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int loadgrid (const char* filename, char *proj4,grid_t *grid, int *ncolumn, int fmt)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;

  switch(fmt) {
    case FMT_CST:
      status= shomcst_loadgrid (filename, grid, 0);
      *ncolumn=1;
      break;

    case FMT_XYZ:
      status=xyz_loadgrid(filename, proj4, grid, ncolumn);
      break;

    case FMT_GRD:
      status=grd_loadgrid(filename, grid);
      *ncolumn=1;
      break;

    case FMT_BMG:
      status=bmg_loadgrid(filename, grid);
      break;

    }

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int loadr1 (char* filename,grid_t grid, int pos, float *buf, float *mask, int fmt )

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{

  int status;
  char output[1024],vname[1024];
  pocgrd_t ncgrid;
  cdfvar_t variable;
  float time;

  switch(fmt) {
    case FMT_CST:
      status= shomcst_read(filename, grid, pos,buf,mask );
      break;

    case FMT_XYZ:
      status=xyz_loadr1(filename,grid,pos,buf,mask);
      break;

    case FMT_GRD:
      status=grd_loadr1(filename,grid,buf,mask);
      break;

    case FMT_BMG:
      status= bmg_loadr1(filename,1,1,1,grid,buf,mask,&time);
      break;

    }

//  if((grid.modeH==2) && (fmt==FMT_CST)) {
  if(grid.modeH==2) {
    if(strrchr(filename,'/')==0)
      sprintf(output,"%s.nc",filename);
     else
      sprintf(output,"%s.nc",(char *) (strrchr(filename,'/')+1));
    status= poc_createfile(output);
    status=poc_sphericalgrid_xyzt(output,"","",grid,&ncgrid);
    sprintf(vname,"varaiable_%3.3d",pos);
    poc_standardvariable_xy(&variable,vname,*mask,"no_units",1., 0.,"unspecified","unspecified","unspecified",ncgrid);
    status=create_ncvariable(output, &variable);
    status=poc_write_xy(output,  grid, variable.id,buf);
    variable.destroy();
    }

  return(status);

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  float *rbuf[2],*tmp[4];
  float *rbufx[2],*rbufy[2],time;

  double  t,p;
  float  rmask;
  float  spec,value;
  int status,fmt=FMT_XYZ,fefmt,ncolumn;
  int k,n;
  int append=0;
  char *keyword,*zone,*from=NULL;
  char *meshfile=NULL,*input=NULL,*fmt_in=NULL,*fmt_out=NULL,*nodefile=NULL;
  char *proj4=0;
  char output[1024]="",filename[1024]="",*comment[2];
  grid_t grid;
  fcomplex *cbuf[2],cmask;
  fcomplex *cbufx[2],*cbufy[2];
  mesh_t mesh,nodes;
  bmgheader_t header;

  fct_echo( argc, argv);

  spec=999.9;

  comment[0]=(char *)malloc(1024);
  comment[1]=(char *)malloc(1024);
  sprintf(comment[0],"The C->CPP transformation");
  sprintf(comment[1],"The C->CPP transformation");

  n=1;
  while (n < argc) {
    keyword=strdup(argv[n]);
    if(strcmp(keyword,"-append")==0) {
      append=1;
      n++;
      continue;
      }
    switch (keyword[0]) {
      case '-':
        switch (keyword[1]) {
        case 'n' :
          nodefile= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'm' :
          meshfile= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'z' :
          zone= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'f' :
          sscanf(argv[n+1],"%s",fmt_out);
          n++;
          n++;
          break;

        case 's' :
          sscanf(argv[n+1],"%f",&spec);
          n++;
          n++;
          break;

        case 't' :
          fmt_in=strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'o' :
          sscanf(argv[n+1],"%s",output);
          n++;
          n++;
          break;

        default:
          STDOUT_BASE_LINE("unknown option %s\n",keyword);
          exit(-1);
        }
        break;

      default:
        if(input==NULL) {
          input= strdup(argv[n]);
          printf("input file=%s\n",input);
          n++;
          }
        else {
          STDOUT_BASE_LINE("unknown option %s\n",keyword);
          exit(-1);
          }
        break;
      }
      free(keyword);
    }

//   status= shomcst_loadgrid ("/home/models/SurfRef/dev-3/data/ModelesMaree/cstMORBIHAN_12042005", &grid, 0);
//   rbuf[0]=new float[grid.nx*grid.ny];
//   status= shomcst_read   ("/home/models/SurfRef/dev-3/data/ModelesMaree/cstMORBIHAN_12042005", grid, 0, rbuf[0], &rmask);

  rmask=spec;
  cmask=fcomplex(spec,spec);

  if(fmt_in==0) fmt_in=strdup("xyz");

  if(strcmp(fmt_in,"xyz")==0) {
    fmt=FMT_XYZ;
    fefmt=S2R;
    }
  else if(strcmp(fmt_in,"ascii")==0) {
    fmt=FMT_ASC;
    fefmt=S2R;
    }
  else if(strcmp(fmt_in,"grd")==0) {
    fmt=FMT_GRD;
    fefmt=S2R;
    }
  else if(strcmp(fmt_in,"bmg")==0){
    fmt=FMT_BMG;
    status=bmg_getinfo (input,&header);
    switch (header.nv) {
      case 1:
        if(header.nd==1) fefmt=S2R;
        if(header.nd==2) fefmt=V2R;
        break;

      case 2:
        if(header.nd==1) fefmt=S2C;
        if(header.nd==2) fefmt=V2C;
        break;

      default:
        goto error;
      }
    }
  else if(strcmp(fmt_in,"cst")==0) {
    fmt=FMT_CST;
    fefmt=S2R;
    }
  else if(strcmp(fmt_in,"netcdf")==0) {
    fmt=FMT_CDF;
    fefmt=S2R;
    }

  if(nodefile!=0) {
    status=fe_readnodes(nodefile,NODE_FILE_FORMAT_TRIGRID,&mesh);
    from=strdup("node-array");
    if(status !=0) goto error;
    }

  if(meshfile!=0) {
    status=fe_readmesh(meshfile,MESH_FILE_FORMAT_TRIGRID,&mesh);
    from=strdup("mesh-array");
    if(status !=0) goto error;
    }

//   if(strrchr(input,'/')==0)
//     sprintf(output,"%s.extraction",input);
//   else sprintf(output,"%s.extraction",(char *) (strrchr(input,'/')+1));

  if(strlen(output)==0) sprintf(output,"%s",input);

  switch (fefmt) {
    case S2R:
/* *-----------------------------------------------------------------------------
      read xyz file and interpret grid */
      status=loadgrid(input,proj4,&grid,&ncolumn,fmt);
      if(status !=0) {
        STDOUT_BASE_LINE("cannot load grid in structured file=%s\n",input);
        exit(-1);
        }
      rbuf[0]=new float[grid.nx*grid.ny];
      rbuf[1]=new float[mesh.nvtxs];
/* *-----------------------------------------------------------------------------
      read xyz database */
      for(k=0;k<ncolumn;k++) {
        status=loadr1(input,grid,k,rbuf[0],&rmask,fmt);
        if(status !=0) {
          STDOUT_BASE_LINE("cannot load data in structured file=%s\n",input);
          exit(-1);
          }
        if(strrchr(output,'/')==0)
          sprintf(filename,"%s.%s.%2.2d.surfref",output,from,k);
        else
          sprintf(filename,"%s.%s.%2.2d.surfref",(char *) (strrchr(output,'/')+1),from,k);
        for (n=0;n<mesh.nvtxs;n++) rbuf[1][n]=rmask;
        if(append==1) {
/* *-----------------------------------------------------------------------------
          load output file if exists (update mode) */
          status=surfref_loadr1(filename, mesh, rbuf[1],comment);
          }
        for (n=0;n<mesh.nvtxs;n++) {
          t=mesh.vertices[n].lon;
          p=mesh.vertices[n].lat;
          status=map_interpolation(grid, rbuf[0], rmask, t, p, &value);
          if(status==0){
            if(value!=rmask) {
              rbuf[1][n]=value;
              }
            }
          }
        status=surfref_saver1(filename, mesh, rbuf[1],comment);
        if(strrchr(output,'/')==0)
          sprintf(filename,"%s.%s.%2.2d.s2r",output,from,k);
        else
          sprintf(filename,"%s.%s.%2.2d.s2r",(char *) (strrchr(output,'/')+1),from,k);
        status=quoddy_saver1(filename, mesh.nvtxs, rbuf[1],comment);
        if(status!=0) goto error;
        }
      delete[] rbuf[0];
      delete[] rbuf[1];
      break;
    case V2R:
      rbufx[0]=(float *)malloc(mesh.nvtxs*sizeof(float));
      rbufx[1]=(float *)malloc(nodes.nvtxs*sizeof(float));
      rbufy[0]=(float *)malloc(mesh.nvtxs*sizeof(float));
      rbufy[1]=(float *)malloc(nodes.nvtxs*sizeof(float));
      if(status!=0) goto error;
      for (n=0;n<nodes.nvtxs;n++) {
        t=nodes.vertices[n].lon;
        p=nodes.vertices[n].lat;
        if(status!=0) goto error;
        else {
          rbufx[1][n]=rmask;
          rbufy[1][n]=rmask;
          }
        }
      status=quoddy_saver2(output, nodes.nvtxs, rbufx[1], rbufy[1],comment);
      free(rbufx[0]);
      free(rbufx[1]);
      free(rbufy[0]);
      free(rbufy[1]);
      break;
    case S2C:
      cbuf[0]=NULL;
      status= bmg_loadc1(input,1,1,1,&grid,&cbuf[0],&cmask,&time);
      cbuf[1]=(fcomplex *)malloc(nodes.nvtxs*sizeof(fcomplex));
      for (n=0;n<nodes.nvtxs;n++) {
        t=nodes.vertices[n].lon;
        p=nodes.vertices[n].lat;
        status=map_interpolation(grid, cbuf[0], cmask, t, p, &cbuf[1][n]);
        if(status!=0) goto error;
        }
      status=quoddy_savec1(output, nodes.nvtxs, cbuf[1],comment);
      free(cbuf[0]);
      free(cbuf[1]);
      break;
    case V2C:
      cmask=fcomplex(999.9,999.9);
      cbufx[0]=(fcomplex *)malloc(mesh.nvtxs*sizeof(fcomplex));
      cbufx[1]=(fcomplex *)malloc(nodes.nvtxs*sizeof(fcomplex));
      cbufy[0]=(fcomplex *)malloc(mesh.nvtxs*sizeof(fcomplex));
      cbufy[1]=(fcomplex *)malloc(nodes.nvtxs*sizeof(fcomplex));
      tmp[0]=(float *)malloc(mesh.nvtxs*sizeof(float));
      tmp[1]=(float *)malloc(mesh.nvtxs*sizeof(float));
      tmp[2]=(float *)malloc(mesh.nvtxs*sizeof(float));
      tmp[3]=(float *)malloc(mesh.nvtxs*sizeof(float));
      for (n=0;n<mesh.nvtxs;n++) tmp[0][n]=real(cbufx[0][n]);
      for (n=0;n<mesh.nvtxs;n++) tmp[1][n]=imag(cbufx[0][n]);
      for (n=0;n<mesh.nvtxs;n++) tmp[2][n]=real(cbufy[0][n]);
      for (n=0;n<mesh.nvtxs;n++) tmp[3][n]=imag(cbufy[0][n]);
      for (n=0;n<nodes.nvtxs;n++) {
        t=nodes.vertices[n].lon;
        p=nodes.vertices[n].lat;
        if(status!=0) goto error;
        else {
          cbufx[1][n]=cmask;
          cbufy[1][n]=cmask;
          }
        if(status!=0) goto error;
        }
      if(status!=0) goto error;
      status=quoddy_savec2(output, nodes.nvtxs, cbufx[1], cbufy[1],comment);
      free(cbufx[0]);
      free(cbufx[1]);
      free(cbufy[0]);
      free(cbufy[1]);
      free(tmp[0]);
      free(tmp[1]);
      free(tmp[2]);
      free(tmp[3]);
      break;

    default:
      break;
    }


  printf("end of regular2node ... \n");

error:
  TRAP_ERR_EXIT(-1,"exiting\n");
}
