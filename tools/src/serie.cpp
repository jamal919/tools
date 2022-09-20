#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "tools-structures.h"
     
#include "map.h"
#include "bmg.h"
#include "poc-time.h"

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 void serie_r1(FILE *in, FILE *out, bmgheader_t header, double *t, double *p)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  grid_t grid;
  float *buf;
  float  h[10],time,mask;
  int m,status;

  bmg_h2g (header, &grid);
  mask=header.spec;
  buf=(float *)malloc(header.ni*header.nj*sizeof(float));

  for (m=0;m<header.nt+1;m++) {
      /* FHL, 18/11/2001: time in NASA days */
    status=bmg_read_r(in, header, 1,1,m+1, buf, &time);
    if(status != 0) break;
    status=map_interpolation(grid,buf,mask,t[0],p[0],&h[0]);
    if(status != 0) break;
    fprintf(out, "%f %f %f %f \n",time,h[0]);
    }

  free(buf);
  return;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 void serie_r2(FILE *in, FILE *out, bmgheader_t header, double *t, double *p)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  grid_t grid;
  float *buf;
  float  h[10],time,mask;
  int m,status;

  bmg_h2g (header, &grid);
  mask=header.spec;
  buf=(float *)malloc(header.ni*header.nj*sizeof(float));

  for (m=0;m<header.nt+1;m++) {
      /* FHL, 18/11/2001: time in NASA days */
    status=bmg_read_r(in, header, 1,1,m+1, buf, &time);
    printf("%s\n",sgetnasadate(time*24.));

    if(status != 0) break;
    status=map_interpolation(grid,buf,mask,t[0],p[0],&h[0]);
    if(status != 0) break;
    status=bmg_read_r(in, header, 1,2,m+1, buf, &time);
    if(status != 0) break;
    status=map_interpolation(grid,buf,mask,t[0],p[0],&h[1]);
    if(status != 0) break;
    fprintf(out, "%f %f %f %f \n",time,h[0],h[1]);
    fflush(out);
    }

  free(buf);
  return;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
 
  FILE *in,*out;
  int i,l,L,m,n,status;
  char *bmgfile,*outfile=NULL;
  char *keyword;
  char rootname[256]="\0",tmp[256];
  double t[10],p[10];
  float h[10],time,mask;
  float *buf,*tmp1,*tmp2,tmp3;
  grid_t grid;
  bmgheader_t header;

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
        case 'b' :
          bmgfile= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'o' :
          outfile= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'p' :
          n++;
          sscanf(argv[n],"%f",&t[0]);
          n++;
          sscanf(argv[n],"%f",&p[0]);
          n++;
          break;

        default:
          printf("unknown option %s\n",keyword);
          __OUT_BASE_LINE__("usage: serie.exe -b {input file} -p {lon lat} -o {output file} \n");
          exit(-1);
        }
        break;

      case '?':
        __OUT_BASE_LINE__("usage: serie.exe -b {input file} -p {lon lat} -o {output file} \n");
        exit(-1);
        break;

      default:
        printf("unknown option %s\n",keyword);
        __OUT_BASE_LINE__("usage: serie.exe -b {input file} -p {lon lat} -o {output file} \n");
        exit(-1);
        break;
      }
      free(keyword);
    }

  l=strlen(strrchr(bmgfile,(int) '/'));
  L=strlen(bmgfile);
  strncpy(rootname,bmgfile,L-l+1);
  printf("rootname: %s \n",rootname);

  if(!outfile) {
    outfile=(char *)malloc(256);
    strncpy(outfile,bmgfile,L-l);
    strcat (outfile,".out");
    }


/*
  p[0]= 43.400;
  t[0]=  3.683;
 
  p[1]= 42.5000;
  t[1]=  4.6435;
*/
/*
  p= 42.875;
  t=  5.875;
1-    3.142    42.532
2-    3.060    42.920
3-    3.693    43.398
4-    5.233    43.450
5-    6.228    43.000*/

  if ((in = fopen(bmgfile, "rb")) == NULL) {
    __ERR_BASE_LINE__("exiting\n");exit(-1);
    }
  if ((out = fopen(outfile, "w")) == NULL) {
    __ERR_BASE_LINE__("exiting\n");exit(-1);
    }

  status=bmg_readheader (in, &header);
  if(header.nd==1) {
    serie_r1(in,out,header,t,p);
    }
  else
    {
    serie_r2(in,out,header,t,p);
    }

  fclose(in);
  fclose(out);
  
}
