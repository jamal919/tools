
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

#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "tools-structures.h"
#include "constants.h"

#include "fe.h"
#include "map.h"
#include "geo.h"
#include "rutin.h"     /*  rutin.h contains common utility routines  */

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

float SquareDst(float lon1, float lat1, float lon2, float lat2)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  float SquareDstq;
  float dX,dY;
  double  R = 6.3675e8;

  lon2=geo_recale(lon2,lon1);

  dX = R*cos(0.5*(lat1+lat2))*(lon1-lon2);
  dY = R*(lat1-lat2);
  SquareDstq = dX*dX + dY*dY;

  return(SquareDstq);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  float zero=0,dummy,min_ds,ds_tst,plon,plat;

  int i,j,k,l,L,n,ndum,ntrack;
  int status,nst,*list,first;
  mesh_t mesh;

  brozouf *track;

  size_t size,offset;
  FILE *in,*out;

  char fout[100][256], outname[256]="\0",rootname[256]="\0",sfile[256];
  char *option,*keyword,*s,*meshfile=NULL,*input=NULL,output[256];
  
  fct_echo(argc,argv);
 
  if(argc < 2) goto abort;
  n=1;
  while (n < argc)
    {
    keyword=strdup(argv[n]);
    switch (keyword[0])
      {
      case '-':
        switch (keyword[1])
        {
        case 'i' :
          input= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'h' :
          goto abort;
          break;

        case 'm' :
          meshfile= strdup(argv[n+1]);
          n++;
          n++;
          break;

        default:
          printf("unknown option %s\n",keyword);
          goto abort;
          break;
        }
        break;

      default:
        n++;
        break;
      }
      free(keyword);
    }

  if(meshfile==NULL) goto abort;
  if(input==NULL) goto abort;
  
  status=fe_readmesh(meshfile,MESH_FILE_FORMAT_TRIGRID,&mesh);

  status=fe_list(&mesh);


  if ((in=fopen(input,"r"))== NULL)
    gmerror("Cannot open input T/P file \n");
    
  sprintf(sfile,"%s",input);
  printf("#sfile='%s'\n",sfile);
  strncpy(output,sfile,strlen(sfile)-3);
  printf("#output='%s'\n",output);
  strcat(output,"beta");
  printf("#output='%s'\n",output);
  
  out=fopen(output,"w");
  
  while(!feof(in))
    {
    fscanf(in, "%d %d ",&ntrack, &nst);
    fprintf(out, "%6d %6d \n",ntrack, nst);
    printf("#ntrak=%d, nst=%d \n",ntrack, nst);
    track=(brozouf *)malloc(nst*sizeof(brozouf));
    for(k=0; k<nst;k++) {
      fscanf(in,"%f %f %d",&track[k].t,&track[k].p,&n);
      printf("lon=%f lat=%f \n",track[k].t, track[k].p);
    /*  track[k].t/=1.e+06;
      track[k].p/=1.e+06;*/
      for(i=0;i<6;i++) track[k].beta[i]=0.0;
      for(i=0;i<6;i++) track[k].node[i]=0;
      status=fe_beta(mesh,track[k].t,track[k].p,&track[k].elt,&(track[k].node[0]),&(track[k].beta[0]));
/*       printf("# ef_beta status=%d \n",status); */
      if(status == 0) {
        track[k].h=0;
        for(i=0;i<3;i++) {
          track[k].h+=track[k].beta[i]*mesh.vertices[track[k].node[i]-1].h;
          }
        fprintf(out,"%12.6f %12.6f %8.1f %6d %6d %6d %6d %8.6f %8.6f %8.6f \n",
                track[k].t, track[k].p,track[k].h/100.,track[k].elt+1,
                track[k].node[0],track[k].node[1],track[k].node[2],
                track[k].beta[0],track[k].beta[1],track[k].beta[2]);
  
/*        printf("%12.6f %12.6f %8.1f %6d %6d %6d %6d %8.6f %8.6f %8.6f \n",  */
/*               track[k].t, track[k].p,track[k].h,track[k].elt, */
/*               track[k].node[0],track[k].node[1],track[k].node[2], */
/*               track[k].beta[0],track[k].beta[1],track[k].beta[2]); */
       }
      else
      {
       min_ds = 1.e30;
       plon = track[k].t * d2r;
       plat = track[k].p * d2r;
       for(j=0; j<mesh.nvtxs; j++) {
         ds_tst = SquareDst(plon,plat,mesh.vertices[j].lon*d2r,mesh.vertices[j].lat*d2r);
         if(ds_tst < min_ds) {
           min_ds = ds_tst;
           track[k].node[0] = j+1;
           }
         }
       /*si distance inf�ieure �10 km */
       if(min_ds/1.e+10 < 100.) {
       track[k].node[1]=1;
       track[k].node[2]=1;
    
       track[k].beta[0]=1.0;
       track[k].beta[1]=0.0;
       track[k].beta[2]=0.0;
    
       track[k].elt=0;
       track[k].h=0;
       for(i=0;i<3;i++) {
         track[k].h+=track[k].beta[i]*mesh.vertices[track[k].node[i]-1].h;
         }
       }
       else
       {
       track[k].node[0]=0;
       track[k].node[1]=0;
       track[k].node[2]=0;
    
       track[k].beta[0]=0.0;
       track[k].beta[1]=0.0;
       track[k].beta[2]=0.0;
    
       track[k].elt=0;
       track[k].h=0;
       }

       fprintf(out,"%12.6f %12.6f %8.1f %6d %6d %6d %6d %8.6f %8.6f %8.6f \n",
               track[k].t, track[k].p,track[k].h/100.,track[k].elt,
               track[k].node[0],track[k].node[1],track[k].node[2],
               track[k].beta[0],track[k].beta[1],track[k].beta[2]);

        /*printf("%12.6f %12.6f %8.1f %6d %6d %6d %6d %8.6f %8.6f %8.6f \n",
               track[k].t, track[k].p,track[k].h,track[k].elt,
               track[k].node[0],track[k].node[1],track[k].node[2],
               track[k].beta[0],track[k].beta[1],track[k].beta[2]); */
    
      }
     
     }
     
    free(track);
    }

  endup: printf("# end of file ...\n");

  fclose(out);
  fclose(in);
  free (list);
  __ERR_BASE_LINE__("exiting\n");exit(0);
  abort:
     printf("Usage: interpolator.exe -i [satellite location file] -m [finite element mesh file] \n");
     __OUT_BASE_LINE__("Output will be saved into [satellite location file].beta\n");
  exit(-1);
  }
