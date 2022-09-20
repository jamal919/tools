
/**************************************************************************

  Double Kelvin Wave for canal test

  
Contributors:

  Florent Lyard      LEGOS/CNRS, Toulouse, France
  Damien Allain      LEGOS/CNRS, Toulouse, France
  Yves Soufflet      LEGOS/CNRS, Toulouse, France

Date: 26/01/2011

***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <cmath>
//#include "constants.h"
#include <complex>

#include "config.h"


#include "fe.h"

#include "geo.h"
#include "functions.h"
#include "tools-structures.h"
#include "tides.h"
#include "archive.h"

#include "rutin.h"     /*  rutin.h contains common utility routines  */

using namespace std;

 
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int get_openlimits_obsolete(char *belfile,char *meshfile, double **x, double **y, int *count)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------

----------------------------------------------------------------------*/
{
  mesh_t mesh;
  int *code;
  int k,l,m,n,n1,n2,status;
  int nopen=0;

  if(meshfile != NULL) {
    status=fe_readmesh(meshfile,MESH_FILE_FORMAT_TRIGRID,&mesh);
    if(status !=0) goto error;
    status=fe_list(&mesh);
    if(status !=0) goto error;
    }
  else {
    printf("no mesh file specified; abort...\n");
    goto error;
    }

  status= fe_edgetable(&mesh,0,0);

  status= fe_vertex_element_tables(&mesh);

  status= fe_codetable2(&mesh,0,1,0);

  status=fe_read_boundarycode(belfile,mesh,2);

  code=(int *) malloc(mesh.nvtxs*sizeof(int));

  for (n=0;n<mesh.nvtxs;n++) code[n]=0;

  for (n=0;n<mesh.nedges;n++) {
    if(mesh.edges[n].code==5) {
      n1=mesh.edges[n].extremity[0];
      n2=mesh.edges[n].extremity[1];
      code[n1]=1;
      code[n2]=1;
      nopen++;
      }
    }

  *count=0;
  for (n=0;n<mesh.nvtxs;n++)
    if(code[n]!=0) (*count)++;

  (*x)=(double *) malloc(*count*sizeof(double));
  (*y)=(double *) malloc(*count*sizeof(double));

/* *----------------------------------------------------------------------------
  check all limits (in most cases not necessary, but who knows users can do?)*/
  m=0;
  for(l=0;l<mesh.nlimits;l++) {
//    printf("treat limit %d...\n",mesh.limits[l].code);
    for(k=0;k<mesh.limits[l].nvertex;k++) {
      n=mesh.limits[l].vertex[k];
      if(code[n]!=0) {
        (*x)[m]=mesh.vertices[n].lon;
        (*y)[m]=mesh.vertices[n].lat;
        m++;
        }
      }
    }

//   m=0;
//   for (n=0;n<mesh.nvtxs;n++)
//     if(code[n]!=0) {
//       (*x)[m]=mesh.vertices[n].lon;
//       (*y)[m]=mesh.vertices[n].lat;
//       m++;
//       }

  free(code);

  return(0);

  error:
  return(-1);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int get_openlimits(const char *belfile, const char *meshfile, double **x, double **y, int *count)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  mesh_t mesh;
  int *code=NULL;
  int k,l,m,n,n1,n2,status;
  int nopen=0;
  int stopon_EdgeError=1, stopon_PinchError=0;

  if(meshfile != NULL) {
    status=fe_readmesh(meshfile,MESH_FILE_FORMAT_TRIGRID,&mesh);
    if(status !=0) goto error;
/**----------------------------------------------------------------------------
  test - test - test - test - test - test - test - test - test - test - test */
//    status=fe_list(&mesh);
    status=quadrangle_list (mesh,3,4);
    if(status !=0) goto error;
    }
  else {
    printf("no mesh file specified; abort...\n");
    goto error;
    }

  if(mesh.ntriangles !=0) {
    status=fe_edgetable(&mesh,0,0);
    }
  else {
    status=fe_edgetable_Q(&mesh,0);
    for(m=0;m<mesh.nquadrangles;m++) {
      status=fe_initaffine_spherical(mesh, &(mesh.quadrangles[m]), m);
      }
    }

  status=fe_vertex_crosstables02(&mesh);

  status=fe_codetable1(&mesh, 0, stopon_EdgeError, stopon_PinchError);
  if(status!=0) goto error;

  status=fe_read_boundarycode(belfile,mesh,2);
  if(status!=0) goto error;

  exitIfNull(
    code=new int[mesh.nvtxs]
    );

  for (n=0;n<mesh.nvtxs;n++) code[n]=0;

  for (n=0;n<mesh.nedges;n++) {
    if(mesh.edges[n].code==5) {
      n1=mesh.edges[n].extremity[0];
      n2=mesh.edges[n].extremity[1];
      code[n1]=1;
      code[n2]=1;
      nopen++;
      }
    }

  *count=0;
  for (n=0;n<mesh.nvtxs;n++)
    if(code[n]!=0) (*count)++;

  exitIfNull(
    (*x)=new double[*count]
    );
  exitIfNull(
    (*y)=new double[*count]
    );


/*------------------------------------------------------------------------------
  check all limits (in most cases not necessary, but who knows users can do?)*/
  m=0;
  for(l=0;l<mesh.nlimits;l++) {
//    printf("treat limit %d...\n",mesh.limits[l].code);
    for(k=0;k<mesh.limits[l].nvertex;k++) {
      n=mesh.limits[l].vertex[k];
      if(code[n]!=0) {
        (*x)[m]=mesh.vertices[n].lon;
        (*y)[m]=mesh.vertices[n].lat;
        m++;
        }
      }
    }

  delete[] code;

  return(0);

  error:
  return(-1);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int load_HamonicOBCs(const char *filename, spectrum_t & spectrum, vector<tidaldata_t> & TidalData)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, idum, j, k, nitems,status=0;
  float fdum;
  char name[1024],units[256],msg[1024],char_dum[1024];
  char *s;
  double factor;
  FILE *in1;
  tidal_wave *dum;
  int NTidalData, NTidalWave;
//   tidaldata_t *TidalData;
  
  spectrum_t reference=initialize_tide();
  
  in1 = fopen(filename, "r");
  
  if(in1 == NULL) {
    status = -1;
    return(status);
    }

  s=fgets(char_dum, 1024, in1);
  nitems=sscanf(char_dum,"%d %d %s", &NTidalData, &NTidalWave, units);

  if(nitems==2) {
    printf("\nno units given in %s, assumes cm\n",filename);
    factor=0.01;
    }
  else {
    if((strcmp(units,"cm")==0) || (strcmp(units,"CM")==0)) {
      factor=0.01;
      }
    if((strcmp(units,"m")==0)||(strcmp(units,"M")==0)) {
      factor=1.0;
      }
    printf("\nunits given in %s: %s \n",filename,units);
    }

  fprintf(stdout, "%d tidal data points. %d tidal input waves\n", NTidalData, NTidalWave);

  if(NTidalData != 0) {
    for(i = 0; i < NTidalData; i++) {
      tidaldata_t tmp;
      tmp.Ha = new float[NTidalWave] ();
      tmp.Hg = new float[NTidalWave] ();
      tmp.Ua = new float[NTidalWave] ();
      tmp.Ug = new float[NTidalWave] ();
      tmp.Va = new float[NTidalWave] ();
      tmp.Vg = new float[NTidalWave] ();
      TidalData.push_back(tmp);
      }
    }

//   WaveList.n     = NTidalWave;
//   WaveList.nmax  = NTidalWave;
//   WaveList.waves = new tidal_wave[WaveList.nmax];
//   WaveList.prescribed = new int[WaveList.nmax];
//   for(i = 0; i < WaveList.nmax; i++) {
//     WaveList.prescribed[i]=-1;
//     }
  spectrum.init(NTidalWave);
  
  j = -1;
  printf("Tidal boundary conditions found for waves : ");
  while(!feof(in1) && j < NTidalWave - 1) {
    s=fgets(char_dum, 1024, in1);
    if(s==0) {
      check_error(-1, "missing wave in OBC tidal file", __LINE__, __FILE__, 1);
      }
    sscanf(char_dum,"%s", name);
/**----------------------------------------------------------------------------
    Patch for name mismatch in FES atlas */
    if(strcmp("Msqm", name)==0) strcpy(name,"MSqm"); 
/**----------------------------------------------------------------------------
    indentify tidal wave from spectrum list */
    k=spectrum.add(reference, name,0);
    j++;
    if(k<0){
      nitems=sscanf(name,"%lg",&spectrum.waves[j].omega);
      if(nitems==1){
        //this is safer than a strcpy
        snprintf(spectrum.waves[j].name,TIDAL_WAVE_NAME_LEN,"%gdph",spectrum.waves[j].omega);
        k=spectrum.n;
        }
      }
    if(k<0) {
      fprintf(stdout, "->>>%s<<<<- not identified %d %d/%d\n", name,nitems,k,spectrum.n);
      fflush(stdout);
      status=-1;
      check_error(status, "unknown tidal wave in OBC's file", __LINE__, __FILE__, 1);
      }

    printf("%s ", spectrum.waves[j].name);

    for(i = 0; i < NTidalData; i++) {
      s=fgets(char_dum, 1024, in1);
      nitems=sscanf(char_dum,"%lf %lf %f %f %f %f %f %f", &TidalData[i].lon, &TidalData[i].lat,
                                                     &(TidalData[i].Ha[j]), &(TidalData[i].Hg[j]),
                                                     &(TidalData[i].Ua[j]), &(TidalData[i].Ug[j]),
                                                     &(TidalData[i].Va[j]), &(TidalData[i].Vg[j]));
//       TidalData[i].lon   *= d2r;
//       TidalData[i].lat   *= d2r;
      TidalData[i].Ha[j] *= factor;
//       TidalData[i].Hg[j] *= d2r;
      switch (nitems) {
        case 4:
          TidalData[i].Ua[j] = 999.9;
          TidalData[i].Ug[j] = 999.9;
          TidalData[i].Va[j] = 999.9;
          TidalData[i].Vg[j] = 999.9;
          break;
        case 8:
          TidalData[i].Ua[j] *= factor;
//           TidalData[i].Ug[j] *= d2r;
          TidalData[i].Va[j] *= factor;
//           TidalData[i].Vg[j] *= d2r;
          break;
        default:
          check_error(-1, "illegal format", __LINE__, __FILE__, 1);
          break;
        }
      }
    }

  printf("\n\n");
  
  fclose(in1);
    
  return(status);

}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int save_TidalOBCs(spectrum_t s, tidalOBC_t *data, int ndata)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int i,j,k,l,m,n;
  int nitems;
  FILE *out=NULL;
  char filename[1024];
  float mask;
  float a[3],G[3];
  mesh_t mesh;
  
  for(k=0;k<s.n;k++) {
    sprintf(filename,"%s.obc",s.waves[k].name);
    out=fopen(filename,"w");
    if ((out=fopen(filename,"w")) ==NULL) {
      fprintf(stderr,"cannot open obc file : %s \n",filename);
      return(-1);
      }
    fprintf(out,"%s\n",s.waves[k].name);
    for(i=0; i<ndata; i++) {
/*---------------------------------------------------------------------
      elevation*/
      if(data[i].Htide.G[k]<0.0) data[i].Htide.G[k]+=360.0;
      a[0]=data[i].Htide.a[k];
      G[0]=data[i].Htide.G[k];
      if((data[i].Utide.a==0) || (data[i].Vtide.a==0)) {
        nitems=fprintf(out,"%12.4f %12.4f %9.4f %9.2f\n",data[i].lon,data[i].lat,a[0],G[0]);
        }
      else {
        a[1]=data[i].Utide.a[k];
        G[1]=data[i].Utide.G[k];
        a[2]=data[i].Vtide.a[k];
        G[2]=data[i].Vtide.G[k];
        if(G[1]<0.0) G[1]+=360.0;
        if(G[2]<0.0) G[2]+=360.0;
        nitems=fprintf(out,"%8f %8f %f %f %f %f %f %f\n",data[i].lon,data[i].lat,a[0],G[0],a[1],G[1],a[2],G[2]);
        }
      }
    }
  fclose(out);

  sprintf(filename,"tides.obc");
  out=fopen(filename,"w");
  if (out ==NULL) {
    fprintf(stderr,"cannot open obc file : %s \n",filename);
    return(-1);
    }

  fprintf(out,"%d %d %s\n",ndata,s.n,"M");
  for(k=0;k<s.n;k++) {
/* *-----------------------------------------------------------------------------
    tidal heights and currents at open limits*/
    fprintf(out,"%s\n",s.waves[k].name);
    for(i=0; i<ndata; i++) {
/*---------------------------------------------------------------------
      elevation*/
      if(data[i].Htide.G[k]<0.0) data[i].Htide.G[k]+=360.0;
      a[0]=data[i].Htide.a[k];
      G[0]=data[i].Htide.G[k];
      if((data[i].Utide.a==0) || (data[i].Vtide.a==0)) {
        nitems=fprintf(out,"%12.4f %12.4f %9.4f %9.2f\n",data[i].lon,data[i].lat,a[0],G[0]);
        }
      else {
        a[1]=data[i].Utide.a[k];
        G[1]=data[i].Utide.G[k];
        a[2]=data[i].Vtide.a[k];
        G[2]=data[i].Vtide.G[k];
        if(G[1]<0.0) G[1]+=360.0;
        if(G[2]<0.0) G[2]+=360.0;
        nitems=fprintf(out,"%8f %8f %f %f %f %f %f %f\n",data[i].lon,data[i].lat,a[0],G[0],a[1],G[1],a[2],G[2]);
        }
      }
    }
  fclose(out);
  
  return(0);
}
