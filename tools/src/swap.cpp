#include <config.h>

#include <stdio.h>
#include <string.h>

#include "tools-structures.h"

#include "swap.h"


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T> T swap_template(T v)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  const size_t n=sizeof(T);
  union{
    T v;
    char b[n];//bytes
  } in, out;
  size_t i;
  
  in.v=v;
  
  for(i=0;i<n;i++)
    out.b[i]=in.b[n-1-i];
  
  return out.v;
}

/*----------------------------------------------------------------------------*/
/// swaps endian
/**
\date created 9 Aug 2011
\author Damien Allain

\param v value to swap
\returns swapped value
*/
/*----------------------------------------------------------------------------*/
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

short swap(short v)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  short result;
  
  result=swap_template(v);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int swap(int v)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=swap_template(v);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

unsigned swap(unsigned v)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  unsigned result;
  
  result=swap_template(v);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

float swap(float v)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  float result;
  
  result=swap_template(v);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

double swap(double v)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double result;
  
  result=swap_template(v);
  
  return result;
}


/*------------------------------------------------------------------------
  
  Name        :  shswap

--------------------------------------------------------------------------*/
short shswap( short shrtint )
{
  union{
    short l;
    unsigned char byte[2];
  } into, outof;
  
  into.l = shrtint;
  outof.byte[0]=into.byte[1];
  outof.byte[1]=into.byte[0];
  return( outof.l );
}

/*------------------------------------------------------------------------
  
  Name        :  lswap

--------------------------------------------------------------------------*/
int lswap( int long_word )
{
  union{
    int l;
    char byte[4];
  } inword, outword;
  
  inword.l=long_word;
  outword.byte[0] = inword.byte[3];
  outword.byte[1] = inword.byte[2];
  outword.byte[2] = inword.byte[1];
  outword.byte[3] = inword.byte[0];
  return( outword.l );
}

/*------------------------------------------------------------------------
  
  Name        :  fswap

--------------------------------------------------------------------------*/
float fswap( float long_word )
{
  union{
    float l;
    char byte[4];
  } inword, outword;
  
  inword.l=long_word;
  outword.byte[0] = inword.byte[3];
  outword.byte[1] = inword.byte[2];
  outword.byte[2] = inword.byte[1];
  outword.byte[3] = inword.byte[0];
  return( outword.l );
}

/*------------------------------------------------------------------------
  
  Name        :  f3swap

--------------------------------------------------------------------------*/
void f3swap(float *long_word)
{

 
  union {
    float l;
    unsigned char byte[4];
    } inword, outword;

  union {
    float *l;
    unsigned char *byte;
    } buffer;
  float dum;
  int i;
  
/*
  buffer= (union{float l; unsigned char byte[4];} *) long_word;
*/
  buffer.l= long_word;
  for (i=0;i<4;i++) inword.byte[i]=buffer.byte[i];
  outword.byte[0] = inword.byte[3];
  outword.byte[1] = inword.byte[2];
  outword.byte[2] = inword.byte[1];
  outword.byte[3] = inword.byte[0];
/*
  printf("%d %d %d %d\n",inword.byte[0],inword.byte[1],
                         inword.byte[2],inword.byte[3]);
  printf("%d %d %d %d\n",outword.byte[0],outword.byte[1],
                         outword.byte[2],outword.byte[3]);
*/
  for (i=0;i<4;i++) buffer.byte[i]=outword.byte[i];
                         
}


/*------------------------------------------------------------------------
  
  Name        :  dswap
ATTENTION: dswap a ete renomm� dswap2, sinon il y avait un pb avec dswap
definie dans kernel/
Cette fonction n'intervient que dans UserInOut(), ligne 184
--------------------------------------------------------------------------*/
double dswap2( double dword )
{
  union {
    double d;
    unsigned char bytes[8];
    } inval, outval,test;

  inval.d = dword;

/* Earlier version, may be working with DEC and pgcc, but NOT with cc*/

  outval.bytes[0] = inval.bytes[7];
  outval.bytes[1] = inval.bytes[6];
  outval.bytes[2] = inval.bytes[5];
  outval.bytes[3] = inval.bytes[4];
  outval.bytes[4] = inval.bytes[3];
  outval.bytes[5] = inval.bytes[2];
  outval.bytes[6] = inval.bytes[1];
  outval.bytes[7] = inval.bytes[0];

  /* Version, may be working with cc, but not pgcc  */
/*
  outval.bytes[0] = inval.bytes[3];
  outval.bytes[1] = inval.bytes[2];
  outval.bytes[2] = inval.bytes[1];
  outval.bytes[3] = inval.bytes[0];
  outval.bytes[4] = inval.bytes[7];
  outval.bytes[5] = inval.bytes[6];
  outval.bytes[6] = inval.bytes[5];
  outval.bytes[7] = inval.bytes[4];
*/


/*
  printf("%d,%d,%d,%d,%d,%d,%d,%d \n",outval.bytes[0],outval.bytes[1],
                                      outval.bytes[2],outval.bytes[3],
                                      outval.bytes[4],outval.bytes[5],
                                      outval.bytes[6],outval.bytes[7]);

  test.d=1.0;
  printf("%d,%d,%d,%d,%d,%d,%d,%d \n",test.bytes[0],test.bytes[1],
                                      test.bytes[2],test.bytes[3],
                                      test.bytes[4],test.bytes[5],
                                      test.bytes[6],test.bytes[7]);

*/
  return( outval.d );
}


/*------------------------------------------------------------------------
  
  Name        :  cswap

--------------------------------------------------------------------------*/
fcomplex cswap( fcomplex dword )
{
  union
  {
    float f;
    unsigned char bytes[4];
  } inval,outval[2];
    
  inval.f = real(dword);
  outval[0].bytes[0] = inval.bytes[3];
  outval[0].bytes[1] = inval.bytes[2];
  outval[0].bytes[2] = inval.bytes[1];
  outval[0].bytes[3] = inval.bytes[0];

  inval.f = imag(dword);
  outval[1].bytes[0] = inval.bytes[3];
  outval[1].bytes[1] = inval.bytes[2];
  outval[1].bytes[2] = inval.bytes[1];
  outval[1].bytes[3] = inval.bytes[0];

  return(fcomplex(outval[0].f,outval[1].f));
}

/*------------------------------------------------------------------------
  
  Name        :  gswap

--------------------------------------------------------------------------*/
grid_t gswap( grid_t grid)
{
  grid_t tmp;

  tmp.nx=lswap(grid.nx);
  tmp.ny=lswap(grid.ny);
  tmp.modeH=lswap(grid.modeH);
  tmp.xmin=dswap2(grid.xmin);
  tmp.ymin=dswap2(grid.ymin);
  tmp.xmax=dswap2(grid.xmax);
  tmp.ymax=dswap2(grid.ymax);
  tmp.dx=dswap2(grid.dx);
  tmp.dy=dswap2(grid.dy);

  return(tmp);
}

/*------------------------------------------------------------------------
  
  Name        :  date_swap

--------------------------------------------------------------------------*/
date_t date_swap( date_t date)
{
  date_t tmp;
  float dum;

  tmp.day=lswap(date.day);
  tmp.month=lswap(date.month);
  tmp.year=lswap(date.year);
    
/*  tmp.second=fswap(date.second); */
  
  dum=date.second;
  f3swap(&dum);
  
  tmp.second=dum;

/*   printf("date_swap: %f %f \n",tmp.second,date.second); */

  return( tmp );
}

/*------------------------------------------------------------------------
  
  Name        :  swap_header

--------------------------------------------------------------------------*/
void swap_header( bmgheader_t *header )
{
 int i,nlevel=1;
  header->ni	= lswap(header->ni);
  header->nj	= lswap(header->nj);
  header->nk	= lswap(header->nk);
  header->nt	= lswap(header->nt);
  header->nd	= lswap(header->nd);
  header->code	= lswap(header->code);

  f3swap(&header->xmin);
  f3swap(&header->ymin);
  f3swap(&header->dx);
  f3swap(&header->dy);

/*
  header->xmin=dswap2(header->xmin);
  header->ymin=dswap2(header->ymin);
  header->dx=dswap2(header->dx);
  header->dy=dswap2(header->dy);
*/

  f3swap(&header->spec);
  for(i=0;i<nlevel;i++)
    f3swap(&header->levels[i]);
}

