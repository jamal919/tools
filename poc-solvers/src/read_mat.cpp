
#include <stdio.h>
#include <string.h> 
#include <malloc.h>
#include "matrix_io.h"


/* ========================================================================== */
/* === is_blank_line ======================================================== */
/* ========================================================================== */

static int is_blank_line	/* TRUE if s is a blank line, FALSE otherwise */
(
    char *s
)
{
    int c, k ;
    for (k = 0 ; k <= MAXLINE ; k++)
    {
	c = s [k] ;
	if (c == '\0') 
	{
	    /* end of line */
	    break ;
	}
	if (!isspace (c))
	{
	    /* non-space character */
	    return (FALSE) ;
	}
    }
    return (TRUE) ;
}

/* ========================================================================== */
/* === get_line ============================================================= */
/* ========================================================================== */

/* Get the next input line, discarding comments. */

static int get_line	/* returns # of items read, or -1 if error */
(
    FILE *f,		/* file to read from */
    long *i,		/* row index */
    long *j,		/* column index */
    double *x,		/* real part */
    double *z,		/* imaginary part */
    int *stype)		/* stype, as determined from Matrix Market header,
			 * but with additional options for the skew-symmetric
	* and complex symmetric cases:
	* 1: symmetric, with upper part stored (not in Matrix Market format
	* 0: unsymmetric (Matrix Market "general")
	* -1: real symmetric or complex Hermitian, with lower part stored
	*	(Matrix Market "real symmetric" or "complex hermitian")
	* -2: real or complex skew symmetric
	* -3: complex symmetric
	*/
{
    char *p, s [MAXLINE+1] ;
    int c, c2, is_complex ;
    *i = 0 ;
    *j = 0 ;
    *x = 0 ;
    *z = 0 ;
    for ( ; ; )
    {
	s [0] = '\0' ;
	s [1] = '\0' ;
	s [MAXLINE] = '\0' ;
	if (fgets (s, MAXLINE, f) == NULL)
	{
	    /* end of file */
	    return (EMPTY) ;
	}
	if (s [0] == '%')
	{
	    /* a comment line */
	    if (strncmp (s, "%%MatrixMarket", 14) == 0)
	    {
		/* this is a Matrix Market header, with the format:
		 * %%MatrixMarket matrix coord type storage */
		p = s ;

		/* get "matrix" token */
		while (*p && !isspace (*p)) p++ ;
		while (*p &&  isspace (*p)) p++ ;
		c = tolower (*p) ;
		if (c != 'm')
		{
		    /* bad format */
		    return (EMPTY) ;
		}

		/* get "coord" token */
		while (*p && !isspace (*p)) p++ ;
		while (*p &&  isspace (*p)) p++ ;
		c = tolower (*p) ;
		if (c != 'c')
		{
		    /* bad format, only "coordinate" is supported;
		     * "array" is not supported */
		    return (EMPTY) ;
		}

		/* get type token (real, pattern, complex, integer) */
		while (*p && !isspace (*p)) p++ ;
		while (*p &&  isspace (*p)) p++ ;
		c = tolower (*p) ;
		if (!(c == 'r' || c == 'p' || c == 'c' || c == 'i'))
		{
		    /* bad format */
		    return (EMPTY) ;
		}

		is_complex = (c == 'c') ;

		/* get storage token (general, hermitian, symmetric, or
		 * skew-symmetric) */
		while (*p && !isspace (*p)) p++ ;
		while (*p &&  isspace (*p)) p++ ;
		c = tolower (*p) ;
		c2 = tolower (*(p+1)) ;
		if (c == 'g')
		{
		    /* "general" storage (unsymmetric matrix) */
		    *stype = 0 ;
		}
		else if (c == 's' && c2 == 'y')
		{
		    /* "symmetric" */
		    if (is_complex)
		    {
			/* complex symmetric, lower triangular part present */
			*stype = -3 ;
		    }
		    else
		    {
			/* real symmetric, lower triangular part present */
			*stype = -1 ;
		    }
		}
		else if (c == 'h')
		{
		    /* "hermitian" matrix, lower triangular part present */
		    *stype = -1 ;
		}
		else if (c == 's' && c2 == 'k')
		{
		    /* "skew-symmetric" (real or complex) */
		    *stype = -2 ;
		}
		else
		{
		    /* bad format */
		    return (EMPTY) ;
		}
	    }
	}
	else
	{
	    /* an entry, or a blank line */
	    if (is_blank_line (s))
	    {
		/* the line is blank, continue and get the next line */
		continue ;
	    }
	    /* this line contains an entry */
	    return (sscanf (s, "%ld %ld %lg %lg\n", i, j, x, z)) ;
	}
    }
}

triplet *read_triplet(char *file,int wanaOne_based)

{
  FILE *f;
  triplet *T ;
  double *Tx ;
  int *Ti, *Tj;
  long l1, l2 ;
  int nitems, nrow, ncol, nnz, stype,ignore,one_based,k;
  int i,j,imax,jmax;
  double x, z ;

  /* Open file*/
  f = fopen (file, "r") ;
  stype = 999 ;
  nitems = get_line (f, &l1, &l2, &x, &z, &stype) ;
  printf("nitems=%d, l12,l2,x,z,stype=%d;%d,%e,%e,%d \n",nitems,l1,l2,x,z,stype);
  ncol = l1;
  nrow=l2;
  nnz = x ;
  if (nrow != ncol)
    {
      stype = 0 ;
    }
  else if (nitems == 4)
    {
      /* first line contains: m n nnz stype */
      if (z < 0)
	{
	  stype = -1 ;
	}
      else if (z > 0)
	{
	  stype = 1 ;
	}
      else
	{
	  stype = 0 ;
	}
    }
    one_based = TRUE ;
    Tx = NULL ;
    Ti = NULL ;
    Tj = NULL ;
    Ti = (int *) malloc(sizeof(int) * nnz); 
    Tj = (int *) malloc(sizeof(int) * nnz); 
    Tx = (double *) malloc(sizeof(double) * nnz); 
    printf("Apres allocate \n");
    for (k = 0 ; k < nnz ; k++)
    {
	nitems = get_line (f, &l1, &l2, &x, &z, &ignore) ;

	i = l1 ;
	j = l2 ;
	Ti [k] = i ;
	Tj [k] = j ;
	Tx [k] = x ;
	if (i == 0 || j == 0)
	{
	    one_based = FALSE ;
	}

	imax = MAX (i, imax) ;
	jmax = MAX (j, jmax) ;
    }
    
    fclose(f);
    printf("lecture ok \n");
    if ( (one_based) && (wanaOne_based==0))
      {
	/* input matrix is one-based; convert matrix to zero-based */
	for (k = 0 ; k < nnz ; k++)
	  {
	    Ti [k]-- ;
	    Tj [k]-- ;
	  }
      }
    if ( (one_based) && (wanaOne_based==0))
      {
	/* input matrix is one-based; convert matrix to zero-based */
	for (k = 0 ; k < nnz ; k++)
	  {
	    Ti [k]-- ;
	    Tj [k]-- ;
	  }
      }
    if ( (!one_based) && (wanaOne_based==1))
      {
	/* input matrix is one-based; convert matrix to zero-based */
	for (k = 0 ; k < nnz ; k++)
	  {
	    Ti [k]++ ;
	    Tj [k]++ ;
	  }
      }
    T = NULL;
    T = (triplet *) malloc(sizeof(triplet));
    printf("lecture ok \n");
    printf("Avant affectation \n");
    printf("T->nnz=%d nnz=,%d\n",T->nnz,nnz);
    T->nnz = nnz;
    printf("T->nnz=%d \n",T->nnz);
    T->nrow = nrow;
    T->ncol = ncol;
    printf("Apres 3 affectation \n");
    T->i = Ti;
    T->j = Tj;
    T->x = Tx;
    T->stype=stype;

    printf("sortie read_mat \n");
    return(T);
}

/*-------------------------------------------------------*/
/*-------------------------------------------------------*/
int read_vect(char *file,long sz,double *Tx)
{
   FILE *f; 
   float tmp;

   int k;
   f = fopen (file, "r") ;
   for (k = 0 ; k < sz ; k++)
     {
       fscanf(f,"%f\n",&tmp);
       Tx[k]=tmp;
    }
   fclose(f);   
 }


/*******************************************************/
/********* Complex routines ****************************/
/*******************************************************/

tripletz *read_tripletz(char *file,int wanaOne_based)

{
  FILE *f;
  tripletz *T ;
  complex<double> *Tx ;
  int *Ti, *Tj;
  long l1, l2 ;
  int nitems, nrow, ncol, nnz, stype,ignore,one_based,k;
  int i,j,imax,jmax;
  double  x, z ;

  /* Open file*/
  f = fopen (file, "r") ;
  stype = 999 ;
  nitems = get_line (f, &l1, &l2, &x, &z, &stype) ;
  printf("nitems=%d, l12,l2,x,z,stype=%d;%d,%e,%e,%d \n",nitems,l1,l2,x,z,stype);
  ncol = l1;
  nrow=l2;
  nnz = x ;
  if (nrow != ncol)
    {
      stype = 0 ;
    }
  else if (nitems == 4)
    {
      /* first line contains: m n nnz stype */
      if (z < 0)
	{
	  stype = -1 ;
	}
      else if (z > 0)
	{
	  stype = 1 ;
	}
      else
	{
	  stype = 0 ;
	}
    }
    one_based = TRUE ;
    Tx = NULL ;
    Ti = NULL ;
    Tj = NULL ;
    Ti = (int *) malloc(sizeof(int) * nnz); 
    Tj = (int *) malloc(sizeof(int) * nnz); 
    Tx = (complex<double> *) malloc(sizeof(complex<double>) * nnz); 
    printf("Apres allocate \n");
    for (k = 0 ; k < nnz ; k++)
    {
	nitems = get_line (f, &l1, &l2, &x, &z, &ignore) ;

	i = l1 ;
	j = l2 ;
	Ti [k] = i ;
	Tj [k] = j ;
	Tx [k] = x ;
	if (i == 0 || j == 0)
	{
	    one_based = FALSE ;
	}

	imax = MAX (i, imax) ;
	jmax = MAX (j, jmax) ;
    }
    
    fclose(f);
    printf("lecture ok \n");
    if ( (one_based) && (wanaOne_based==0))
      {
	/* input matrix is one-based; convert matrix to zero-based */
	for (k = 0 ; k < nnz ; k++)
	  {
	    Ti [k]-- ;
	    Tj [k]-- ;
	  }
      }
    if ( (one_based) && (wanaOne_based==0))
      {
	/* input matrix is one-based; convert matrix to zero-based */
	for (k = 0 ; k < nnz ; k++)
	  {
	    Ti [k]-- ;
	    Tj [k]-- ;
	  }
      }
    if ( (!one_based) && (wanaOne_based==1))
      {
	/* input matrix is one-based; convert matrix to zero-based */
	for (k = 0 ; k < nnz ; k++)
	  {
	    Ti [k]++ ;
	    Tj [k]++ ;
	  }
      }
    T = NULL;
    T = (tripletz *) malloc(sizeof(tripletz));
    printf("lecture ok \n");
    printf("Avant affectation \n");
    printf("T->nnz=%d nnz=,%d\n",T->nnz,nnz);
    T->nnz = nnz;
    printf("T->nnz=%d \n",T->nnz);
    T->nrow = nrow;
    T->ncol = ncol;
    printf("Apres 3 affectation \n");
    T->i = Ti;
    T->j = Tj;
    T->x = Tx;
    T->stype=stype;

    printf("sortie read_mat \n");
    return(T);
}

/*-------------------------------------------------------*/
/*-------------------------------------------------------*/
int read_vectz(char *file,long sz,complex<double> *Tx)
{
   FILE *f; 
   float tmp;
   
   int k;
   f = fopen (file, "r") ;
   for (k = 0 ; k < sz ; k++)
     {
       fscanf(f,"%f\n",&tmp);
       Tx[k]= tmp;
    }
   fclose(f);   
 }

