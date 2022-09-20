#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
 
#include "tools-structures.h"

#include "map.h"
#include "bmg.h"
#include "fe.h"

#include "rutin.h"     /*  rutin.h contains common utility routines  */


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int count=0,z;

  double  t,p;
  float  dummy,rmask;
  float  spec;
  int nndes,status,*elts,fmt,fefmt;
  int i,j,k,l,L,m,n;
  int nitems,element;
  FILE *file;
  FILE *out;
  char *keyword,*zone;
  char *meshfile=NULL,*sfile=NULL,*cfile=NULL,*format=NULL,*nodefile=NULL;
  char file1[256],file2[256];
  char *root,output[1024];
  mesh_t mesh,nodes;

  fct_echo(argc,argv);
 
  spec=999.9;

  n=1;
  while (n < argc)
    {
    keyword=strdup(argv[n]);
    switch (keyword[0])
      {
      case '-':
        switch (keyword[1])
        {
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

        case 's' :
          sscanf(argv[n+1],"%f",&spec);
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
        else
          {
          __OUT_BASE_LINE__("unknown option %s\n",keyword);
          exit(-1);
          }
        break;
      }
      free(keyword);
    }

  out=fopen("/home/models/issyk-kul/topo/sounding.xyz","w");
  
  nodes.vertices=new vertex_t[100000];//1E5 elements
  nodes.nvtxs=0;
  
  file=fopen("/home/models/issyk-kul/topo/iki.vec","r");

  n=-1;
  while(n!=0)
    {
    fscanf(file,"%d %d",&z,&n);
    if(z!=0) {
      for(i=0;i<n;i++) {
        fscanf(file,"%lf %lf",&(nodes.vertices[i].lon),&(nodes.vertices[i].lat));
        }
      continue;
      }
    printf("%d %d %d\n",z,n,nodes.nvtxs);
    count ++;
    for(i=0;i<n;i++) {
      fscanf(file,"%lf %lf",&(nodes.vertices[i].lon),&(nodes.vertices[i].lat));
      nodes.vertices[i].h=z;
      fprintf(out,"%lf %lf %d\n",nodes.vertices[i].lon,nodes.vertices[i].lat,z);
      }
    nodes.nvtxs+=n;
    }
    
  fclose(file);

  file=fopen("/home/models/issyk-kul/topo/iki.vec","r");

  n=-1;
  while(n!=0)
    {
    fscanf(file,"%d %d",&z,&n);
    if(z==0) {
      for(i=0;i<n;i++) {
        fscanf(file,"%lf %lf",&(nodes.vertices[i].lon),&(nodes.vertices[i].lat));
        }
      continue;
      }
    printf("%d %d %d\n",z,n,nodes.nvtxs);
    count ++;
    for(i=0;i<n;i++) {
      fscanf(file,"%lf %lf",&(nodes.vertices[i].lon),&(nodes.vertices[i].lat));
      nodes.vertices[i].h=z;
      fprintf(out,"%lf %lf %d\n",nodes.vertices[i].lon,nodes.vertices[i].lat,z);
      }
    nodes.nvtxs+=n;
    }
    
  fclose(file);


  file=fopen("/home/models/issyk-kul/topo/ikb.pnv","r");

  n=-1;
  while(n!=0)
    {
    fscanf(file,"%d %d",&z,&n);
    printf("%d %d %d\n",z,n,nodes.nvtxs);
    count ++;
    for(i=0;i<n;i++) {
      fscanf(file,"%lf %lf",&(nodes.vertices[i].lon),&(nodes.vertices[i].lat));
      nodes.vertices[i].h=z;
      if(z>100)
        fprintf(out,"%lf %lf %d\n",nodes.vertices[i].lon,nodes.vertices[i].lat,z);
      }
    nodes.nvtxs+=n;
    }
    


  fclose(out);
  
  __ERR_BASE_LINE__("exiting\n");exit(0);
 
  file=fopen("/home/models/issyk-kul/topo/iki.vec","r");

  out=fopen("/home/models/issyk-kul/topo/shoreline.scan","w");
  
  
  nodes.vertices=new vertex_t[100000];//1E5 elements
  nodes.nvtxs=0;
  
  while(n!=0)
    {
    fscanf(file,"%d %d",&z,&n);
    if(z!=0) {
      for(i=0;i<n;i++) {
        fscanf(file,"%lf %lf",&(nodes.vertices[i].lon),&(nodes.vertices[i].lat));
        }
      continue;
      }
    printf("%d %d %d\n",z,n,nodes.nvtxs);
    count ++;
    fprintf(out,"%d %d\n",count,n);
    for(i=0;i<n;i++) {
      fscanf(file,"%lf %lf",&(nodes.vertices[i].lon),&(nodes.vertices[i].lat));
/*      printf("%d %lf %lf\n",i+1,nodes.vertices[i].lon,nodes.vertices[i].lat);*/
      nodes.vertices[i].h=z;
      fprintf(out,"%d %lf %lf\n",i+1,nodes.vertices[i].lon,nodes.vertices[i].lat);
      }
    nodes.nvtxs+=n;
    }
    
  fclose(file);
  fclose(out);
 


end: printf("end of delaunay ... \n");
  free(elts);
error:
  __ERR_BASE_LINE__("exiting\n");exit(-1);
}
