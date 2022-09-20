#define MAIN_SOURCE

#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
 

#include "tools-structures.h"

#include "fe.h"
#include "map.h"
#include "archive.h"
#include "functions.h"


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  float *rbuf[2],*tmp[4];
  float *rbufx[2],*rbufy[2];

  float  t,p;
  float  dummy,rmask;
  float  spec,waiting_real_float,waiting_imag_float;
  int nndes,status,*elts,fmt,fefmt;
  int i,j,k,l,L,m,n;
  int nitems,element;
  FILE *file;
  FILE *out;
  char *keyword,*zone;
  char *meshfile[2]={NULL,NULL},*sfile=NULL,*cfile=NULL,*format=NULL,*nodefile=NULL;
  char *comment[2];
  char *root,output[1024],rootname[1024];
  grid_t grid;
  fcomplex *cbuf[2],cmask;
  fcomplex *cbufx[2],*cbufy[2];
  mesh_t mesh[2];
  list_t list;
  spec=999.9;

  fct_echo( argc, argv);

  n=1;
  while (n < argc) {
    keyword=strdup(argv[n]);
    switch (keyword[0]) {
      case '-':
        switch (keyword[1]) {
        case 'i' :
          meshfile[0]= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'o' :
          meshfile[1]= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 's' :
          sscanf(argv[n+1],"%f",&spec);
          n++;
          n++;
          break;

        case 'f' :
          format=strdup(argv[n]);
          n++;
          n++;
          break;

        default:
          __OUT_BASE_LINE__("unknown option %s\n",keyword);
          exit(-1);
        }
        break;

      default:
        if(sfile==NULL) {
          sfile= strdup(argv[n]);
          printf("input file=%s\n",sfile);
          n++;
          }
        else {
          __OUT_BASE_LINE__("unknown option %s\n",keyword);
          exit(-1);
          }
        break;
      }
      free(keyword);
    }

  rmask=spec;
  cmask=fcomplex(spec,spec);

  if(strstr(sfile,".s2r") != NULL)    fefmt=S2R;
  if(strstr(sfile,".s2c") != NULL)    fefmt=S2C;
  if(strstr(sfile,".v2r") != NULL)    fefmt=V2R;
  if(strstr(sfile,".v2c") != NULL)    fefmt=V2C;

  for(k=0;k<2;k++) {
    if(meshfile[k] != NULL) {
      status=fe_readmesh(meshfile[k],MESH_FILE_FORMAT_TRIGRID,&mesh[k]);
      if(status !=0) goto error;
      status=fe_list(&(mesh[k]));
      if(status !=0) goto error;
      }
    else {
      printf("no mesh file specified; abort...\n");
      goto error;
      }
    }

  fe_AllocateList(mesh[1], &grid,&list);
  fe_CreateTriangleList(grid, mesh[0], &list);

  //fe_initaffine(&mesh);

  if(strrchr(sfile,(int) '/') != NULL) {
    l=strlen(strrchr(sfile,(int) '/'));
    strncpy(rootname,strrchr(sfile,(int) '/')+1,l);
    }
  else  strcpy(rootname,sfile);

  comment[0]=(char *)malloc(256);
  comment[1]=(char *)malloc(256);

  switch (fefmt){
    case S2R:
      rbuf[0]=(float *)malloc(mesh[0].nvtxs*sizeof(float));
      rbuf[1]=(float *)malloc(mesh[1].nvtxs*sizeof(float));
      status=quoddy_loadr1(sfile, mesh[0].nvtxs, rbuf[0]);
      if(status!=0) goto error;
      for (n=0;n<mesh[1].nvtxs;n++) {
        t=mesh[1].vertices[n].lon;
        p=mesh[1].vertices[n].lat;
/*         element=fe_whichelement(mesh,t,p); */
        t=map_recale(grid,t);
        k=(int)( floor((t-grid.xmin)/grid.dx) );
        l=(int)( floor((p-grid.ymin)/grid.dy) );
        element=fe_whichelement_inlist( mesh[0],list.elements[k][l], (double) t,(double) p);
        if(element>0) {
          status= fe_intpl_LGP1(mesh[0], rbuf[0], t,p,element,&rbuf[1][n]);
          if(status!=0) goto error;
          }
//        else rbuf[1][n]=rmask;
        else {
          m=fe_nearest_vertex (mesh[0],t,p,0,0);
          rbuf[1][n]=rbuf[0][m];
          }
        }
      sprintf(output,"%s.new.s2r",rootname);
      sprintf(comment[0],"%c POC/Noveltis, extracted from %s",169,rootname);
      sprintf(comment[1],"what?");
      status=quoddy_saver1(output, mesh[1].nvtxs, rbuf[1],comment);
      if(status!=0) goto error;
      free(rbuf[0]);
      free(rbuf[1]);
      break;
    case V2R:
      rbufx[0]=(float *)malloc(mesh[0].nvtxs*sizeof(float));
      rbufx[1]=(float *)malloc(mesh[1].nvtxs*sizeof(float));
      rbufy[0]=(float *)malloc(mesh[0].nvtxs*sizeof(float));
      rbufy[1]=(float *)malloc(mesh[1].nvtxs*sizeof(float));
      sprintf(comment[0],"Mog2D simulation, %c POC/Noveltis, extracted from %s",169,rootname);
      sprintf(comment[1],"model bathymetry (m)");
      status=quoddy_loadr2(sfile, mesh[0].nvtxs, rbufx[0], rbufy[0]);
      if(status!=0) goto error;
      for (n=0;n<mesh[1].nvtxs;n++) {
        t=mesh[1].vertices[n].lon;
        p=mesh[1].vertices[n].lat;
/*         element=fe_whichelement(mesh,t,p); */
        t=map_recale(grid,t);
        k=(int)( floor((t-grid.xmin)/grid.dx) );
        l=(int)( floor((p-grid.ymin)/grid.dy) );
        element=fe_whichelement_inlist( mesh[0],list.elements[k][l], (double) t,(double) p);
        if(element>0) {
          status= fe_intpl_LGP1(mesh[0], rbufx[0], t,p,element,&rbufx[1][n]);
          status= fe_intpl_LGP1(mesh[0], rbufy[0], t,p,element,&rbufy[1][n]);
          if(status!=0) goto error;
          }
        else {
          rbufx[1][n]=rmask;
          rbufy[1][n]=rmask;
          }
        }
      status=quoddy_saver2(output, mesh[1].nvtxs, rbufx[1], rbufy[1],comment);
      free(rbufx[0]);
      free(rbufx[1]);
      free(rbufy[0]);
      free(rbufy[1]);
      break;
    case S2C:
      cbuf[0]=(fcomplex *)malloc(mesh[0].nvtxs*sizeof(fcomplex));
      cbuf[1]=(fcomplex *)malloc(mesh[1].nvtxs*sizeof(fcomplex));
      tmp[0]=(float *)malloc(mesh[0].nvtxs*sizeof(float));
      tmp[1]=(float *)malloc(mesh[0].nvtxs*sizeof(float));
      status=quoddy_loadc1(sfile,mesh[0].nvtxs , cbuf[0]);
      for (n=0;n<mesh[0].nvtxs;n++) tmp[0][n]=real(cbuf[0][n]);
      for (n=0;n<mesh[0].nvtxs;n++) tmp[1][n]=imag(cbuf[0][n]);
      for (n=0;n<mesh[1].nvtxs;n++) {
        t=mesh[1].vertices[n].lon;
        p=mesh[1].vertices[n].lat;
/*         element=fe_whichelement(mesh,t,p); */
        t=map_recale(grid,t);
        k=(int)( floor((t-grid.xmin)/grid.dx) );
        l=(int)( floor((p-grid.ymin)/grid.dy) );
        element=fe_whichelement_inlist( mesh[0],list.elements[k][l], (double) t,(double) p);
        if(element>0) {
          status= fe_intpl_LGP1(mesh[0], tmp[0], t,p,element,&waiting_real_float);
          status= fe_intpl_LGP1(mesh[0], tmp[1], t,p,element,&waiting_imag_float);
          cbuf[1][n]=fcomplex(waiting_real_float,waiting_imag_float);
          if(status!=0) goto error;
          }
        else cbuf[1][n]=cmask;
        if(status!=0) goto error;
        }
      if(status!=0) goto error;
      sprintf(comment[0],"Mog2D simulation, %c POC/Noveltis, extracted from %s",169,rootname);
      sprintf(comment[1],"tidal elevation, harmonic constants: amplitude and phase lag ref. Greenwich (m,)");
      status=quoddy_savec1(output, mesh[1].nvtxs, cbuf[1],comment);
      free(cbuf[0]);
      free(cbuf[1]);
      free(tmp[0]);
      free(tmp[1]);
      break;
    case V2C:
      //cmask.r=999.9;
      //cmask.i=999.9; CPP MODIF Thierry
      cmask=fcomplex(999.9,999.9);
      cbufx[0]=(fcomplex *)malloc(mesh[0].nvtxs*sizeof(fcomplex));
      cbufx[1]=(fcomplex *)malloc(mesh[1].nvtxs*sizeof(fcomplex));
      cbufy[0]=(fcomplex *)malloc(mesh[0].nvtxs*sizeof(fcomplex));
      cbufy[1]=(fcomplex *)malloc(mesh[1].nvtxs*sizeof(fcomplex));
      tmp[0]=(float *)malloc(mesh[0].nvtxs*sizeof(float));
      tmp[1]=(float *)malloc(mesh[0].nvtxs*sizeof(float));
      tmp[2]=(float *)malloc(mesh[0].nvtxs*sizeof(float));
      tmp[3]=(float *)malloc(mesh[0].nvtxs*sizeof(float));
      status=quoddy_loadc2(sfile,mesh[0].nvtxs , cbufx[0], cbufy[0]);
      for (n=0;n<mesh[0].nvtxs;n++) tmp[0][n]=real(cbufx[0][n]);
      for (n=0;n<mesh[0].nvtxs;n++) tmp[1][n]=imag(cbufx[0][n]);
      for (n=0;n<mesh[0].nvtxs;n++) tmp[2][n]=real(cbufy[0][n]);
      for (n=0;n<mesh[0].nvtxs;n++) tmp[3][n]=imag(cbufy[0][n]);
      for (n=0;n<mesh[1].nvtxs;n++) {
        t=mesh[1].vertices[n].lon;
        p=mesh[1].vertices[n].lat;
/*         element=fe_whichelement(mesh,t,p); */
        t=map_recale(grid,t);
        k=(int)( floor((t-grid.xmin)/grid.dx) );
        l=(int)( floor((p-grid.ymin)/grid.dy) );
        element=fe_whichelement_inlist( mesh[0],list.elements[k][l], (double) t,(double) p);
        if(element>0) {
          status= fe_intpl_LGP1(mesh[0], tmp[0], t,p,element,&waiting_real_float);
          status= fe_intpl_LGP1(mesh[0], tmp[1], t,p,element,&waiting_imag_float);
          cbufx[1][n]=fcomplex(waiting_real_float,waiting_imag_float);
          status= fe_intpl_LGP1(mesh[0], tmp[2], t,p,element,&waiting_real_float);
          status= fe_intpl_LGP1(mesh[0], tmp[3], t,p,element,&waiting_imag_float);
          cbufy[1][n]=fcomplex(waiting_real_float,waiting_imag_float);

          if(status!=0) goto error;
          }
        else {
          cbufx[1][n]=cmask;
          cbufy[1][n]=cmask;
          }
        if(status!=0) goto error;
        }
      if(status!=0) goto error;
      sprintf(comment[0],"Mog2D simulation, %c POC/Noveltis, extracted from %s",169,rootname);
      sprintf(comment[1],"tidal currents (east,north), harmonic constants: amplitude and phase lag ref. Greenwich (m/s,)");
      status=quoddy_savec2((const char *) output, mesh[1].nvtxs, cbufx[1], cbufy[1],comment);
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


end: printf("end of mesh2mesh ... \n");
error:
  __ERR_BASE_LINE__("exiting\n");exit(-1);
}
