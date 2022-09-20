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

#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "tools-structures.h"

#include "rutin.h"
#include "netcdf.h"
#include "geo.h"
#include "swap.h"
#include "sym-io.h"
#include "netcdf-proto.h"
#include "polygones.h"

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int read_graph2d(char *input, grid_t grid, float *buf, char *mask, float spec,int nrec)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

{
  int k,nitems,status=0;
  long offset;
  FILE *in;

  in=fopen(input,"r");

  offset=grid.nx*grid.ny*sizeof(float)*(nrec-1);
  fseek(in,offset,SEEK_SET);
  nitems=fread(buf,4,grid.nx*grid.ny,in);

#if NEED_SWAP == 1
  for(k=0;k<nitems;k++) f3swap(&buf[k]);
#endif

  for (k=0;k<nitems;k++) {
    if(mask[k]==0) buf[k]=spec;
    }

  fclose(in);

  if(nitems==grid.nx*grid.ny) status=0;
  else status=-1;

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int read_graphinfo(char *input, char *wave, int **records,int *n)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,nitems,status=0;
  long offset;
  FILE *in;
  char line[1024],*pointer;
  char dummy[1024];
  const char *identifier[]={"amp maree:","phase maree:","u_A:","u_G:","v_A:","v_G:"};

  int nid=6;

  in=fopen(input,"r");

  *records=(int *)malloc(nid*sizeof(int));
  for (k=0;k<nid;k++) {
    (*records)[k]=-1;
    }

  while(!feof(in)){
    fgets(line,sizeof(line),in);
    pointer=strstr(line,wave);
    if(pointer != NULL) {
      for (k=0;k<nid;k++){
        pointer=strstr(line,identifier[k]);
        if(pointer != NULL) {
          sscanf(pointer,"%s %s %d",dummy,dummy,&(*records)[k]);
          printf("%s %s: %d \n",identifier[k],wave,(*records)[k]);
          }
        }
      }
    }

  fclose(in);

  for (k=0;k<nid;k++) {
    if((*records)[k]==-1) status=-1;
    }
  *n=nid;
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
#define rsize 1000

  float  z;
  geo_t projection;
  float spec=1.e+10;
  double  x,y;
  int i,j,k,l,m,n,ncid,status;
  int nid,*records;
  size_t size;
  FILE *in,*out;
  char *wave,*mask;
  char *keyword,*zone,*s;
  char *rootname=NULL,*output=NULL,*input=NULL,*pathname=NULL,*info=NULL;
  char *bathy=NULL,*notebook=NULL,*poly=NULL,*maskfile=NULL;
  grid_t cgrid,grid;
  float  *buf,*tmp;
  spectrum_t spectrum;
  size_t nitems,offset;


  plg_t *polygones=NULL;
  int npolygones=0;

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
        case 'r' :
          rootname= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'd' :
          pathname= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'i' :
          info= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'm' :
          maskfile= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'n' :
          notebook= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'w' :
          wave= strdup(argv[n+1]);
          n++;
          n++;
          break;

        default:
          __OUT_BASE_LINE__("unknown option %s\n",keyword);
          exit(-1);
        }
        break;

      default:
        if(input==NULL) {
          input= strdup(argv[n]);
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

  output=(char *)malloc(1024);

  if(rootname==NULL) rootname= strdup("symphonie");
  if(wave==NULL) {
    printf("wave missing ( -w [wavename] )...\n");
    goto error;
    }

  if(info==NULL) {
   printf("information file missing ( -i [filename] )...\n");
    goto error;
    }

  if(mask==NULL) {
    printf("mask file missing ( -m [filename] )...\n");
    goto error;
    }

  if(input==NULL) {
    printf("input file missing ( [filename] )...\n");
    goto error;
    }

/*-----------------------------------------------------------------------------
  Read notebook data */
  if(notebook!=0) {
    status=load_notebook(notebook, &cgrid, &grid, &projection);
    printf("%s (notebook file) processed\n",notebook);
    }
  else {
    printf("notebook file missing (-n [filename] )...\n");
    goto error;
    }

  sprintf(output,"%s.barotropic-tide.nc",wave);
  
/**----------------------------------------------------------------------------
  obsolete call, will be suppressed */
  status=create_headershort(output, grid, spectrum);
  if(status !=0) goto error;

  size=grid.nx*grid.ny;
  buf=(float *)malloc(size*sizeof(float));
  mask=(char *)malloc(size*sizeof(char));
  for (k=0;k<size;k++) {
    mask[k]=1;
    }

  spectrum.n=1;
  spectrum.waves=new tidal_wave[spectrum.n];

  for (k=0;k<spectrum.n;k++) {
    strcpy(spectrum.waves[k].name,wave);
    }

  status=read_graphinfo(info, wave, &records,&nid);
  if(status!=0) goto error;

  status= read_graph2d(maskfile, grid,buf, mask, spec, 2);
  if(status!=0) goto error;
  for (k=0;k<size;k++) {
    if(buf[k]==0.) mask[k]=0;
    else             mask[k]=1;
    }

  k=0;
  status=nc_open(output,NC_WRITE,&ncid);
  if(status !=0) goto error;

/**----------------------------------------------------------------------------
  obsolete call, will be suppressed */
  status=nc_writeaxis(ncid,grid);
  if(status !=0) goto error;

  status=nc_writespectrum(ncid, spectrum);
  if(status !=0) goto error;

  for (k=0;k<nid;k++) {
    status= read_graph2d(input, grid,buf, mask, spec, records[k]);
    if(status!=0) goto error;
/**----------------------------------------------------------------------------
    obsolete call, will be suppressed */
    status=nc_write_r1(ncid, grid, 4+k, 0, buf, spec);
    if(status!=0) goto error;
    }

  status = nc_close(ncid);
  __ERR_BASE_LINE__("exiting\n");exit(0);

 error:
  __ERR_BASE_LINE__("exiting\n");exit(-1);

}
