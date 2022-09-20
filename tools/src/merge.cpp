
/**************************************************************************

  T-UGO tools, 2006-2011

  Unstructured Ocean Grid initiative

***************************************************************************/
/** \file

\author  Florent Lyard      LEGOS/CNRS, Toulouse, France. florent.lyard@legos.obs-mip.fr
\author  Yoann Le Bars      LEGOS, Toulouse, France (PhD)
\author  Laurent Roblou     LEGOS/CNRS, Toulouse, France
\author  Damien Allain      LEGOS/CNRS, Toulouse, France
\author  Sara Fleury        LEGOS/CNRS, Toulouse, France
\author  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
\author  Frédéric Dupont    Université de Laval à Québec, Canada

<!-- USE FIND AND REPLACE TO EDIT THIS LINE, SO THAT print_help IS ALSO UPDATED -->
\brief USED BY CTOH ?

<!-- A LINK TO main() or print_help() WILL NOT LINK TO THE RIGHT SOURCE ! -->
See the main function for how this works
and the print_help function for how to use this.
*/
/*----------------------------------------------------------------------------*/

#define MAIN_SOURCE

#include "version-macros.def" //for VERSION and REVISION

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "rutin.h"
#include "functions.h"
#include "poc-assertions.h"


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void print_help(const char *prog_name)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/** \brief prints help of this programme.

\param prog_name name of the program to be printed by this help function : argv[0]
*/
/*----------------------------------------------------------------------------*/
{
/** The code of the body of this function is :
    \code /**/ // COMPILED CODE BELOW !!!
  printf("\n"
    "NAME AND VERSION\n"
    "  %s version " VERSION " " REVISION "\n", prog_name);
  printf("\n"
    "USE\n"
    "  %s ...\n",prog_name);
  printf("\n"
    "DESCRIPTION\n"
/* USE FIND AND REPLACE TO EDIT THE LINE BELOW, SO THAT THE FILE HEAD IS ALSO UPDATED */
    "  ... USED BY CTOH ?\n"
    "  ...\n"
    "\n"
    "OPTIONS :\n"
    "  -h  Show this help and exit.\n"
    "  -m  \n"
    "  -o  \n"); /** \endcode */
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  float z,dum[3],velocity[2],elevation,dummy,*buffer;
  double time,time2,tlim,dT;
  int fmt,i,j,k,l,ndum,rstatus,shift;
  int n,status,nst;
  int *sample,hr,mn,sc;
  size_t size,offset;
  FILE *out,*in1,*in2, *in;
  char *file1=NULL,*file2=NULL,*file3=NULL, *mesh, *keyword;
 
 //correction d'erruers le 15Nov2005 Thierry
    FILE *in5;
    int Nnodes,max_nb;
  fprintf(stderr,"%s  -starting computation *********\n",argv[0]);

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

        case 'm' :
          mesh= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'o' :
          file3= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'h':
          print_help(argv[0]);
          exit(0);

        default:
          __OUT_BASE_LINE__("unknown option %s\n",keyword);
          print_help(argv[0]);
          exit(-1);
        }
        break;

      default:
        if(file1 == NULL) file1= strdup(argv[n]);
        else file2= strdup(argv[n]);
        n++;
        break;
      }
      free(keyword);
    }


  if(mesh == NULL) {__ERR_BASE_LINE__("exiting\n");exit(-1);}
  if(file1 == NULL) {__ERR_BASE_LINE__("exiting\n");exit(-1);}
  if(file2 == NULL) {__ERR_BASE_LINE__("exiting\n");exit(-1);}
  if(file3 == NULL) {__ERR_BASE_LINE__("exiting\n");exit(-1);}

  if ((in5 = fopen(mesh, "r")) == NULL) {
    gmerror("Cannot open NeighbourFile in input_and_Mallocate()");
    }

  fscanf(in5, "%d ", &Nnodes);
//  printf("#Number of nodes      : %d \n",Nnodes);
  fscanf(in5, "%d ", &max_nb);
  fscanf(in5, "%f %f %f %f", &dummy,&dummy,&dummy,&dummy);

  fclose(in5);

  printf("merge %s with %s in %s\n",file1,file2,file3);

  in1=fopen(file1,"rb");
  in2=fopen(file2,"rb");
  out=fopen(file3,"wb");

  buffer=(float *)malloc(Nnodes*3*sizeof(float));

  size=fread(&time,sizeof(double),1,in1);
  offset=3*Nnodes*sizeof(float);
  fseek(in1,offset,SEEK_CUR);
  size=fread(&time2,sizeof(double),1,in1);

  dT=time2-time;

  size=fread(&tlim,sizeof(double),1,in2);
  rewind(in1);
  rewind(in2);
  
  i=0;
  while (size != 0)
    {
    size=0;
    i++;
    if(time < tlim-dT) in=in1;
    else in=in2;
    if(!(size=fread(&time,sizeof(double),1,in))) break;
    fwrite(&time,sizeof(double),1,out);
    if(fmod(i,50.0) ==0)
      printf("t= %6.1f hours, t= %6.1f hours\n",time/3600.0,tlim/3600.0);
    size=fread (buffer, sizeof(float),3*Nnodes,in);
    size=fwrite(buffer, sizeof(float),3*Nnodes,out);
    }
  
  fclose(in1);
  fclose(in2);

  fclose(out);


  __ERR_BASE_LINE__("%s -computation complete ^^^^^^^^^\n",argv[0]);
  exit(0);
}
