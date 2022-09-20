
/**************************************************************************

  POC-SOLVERS interface, 2006-2012

  Part of the Unstructured Ocean Grid initiative and SYMPHONIE suite

Contributors:

  Florent Lyard      LEGOS/CNRS, Toulouse, France
  Cyril Nguyen       LA/CNRS,    Toulouse, France
  Damien Allain      LEGOS/CNRS, Toulouse, France
  Yoann Le Bars      PhD, LEGOS, Toulouse, France

***************************************************************************/

#include <malloc.h>

#include "poc-solvers.h"
#include "matrix_io.h"

#include "string.h"

// #ifdef UMFPACK
// #include "umfpack.h"
// #endif

/* convert triplet between each form coo, csr, csc */
/* coo = 0, csc = 1, csr =2*/

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int convert(triplet *Mat, int destform, int base, int verbose) {

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  int offset;
  int status=-1;
  
  /*printf("base=%d\n",base);*/
  offset = base - Mat->base;
  
  if (Mat->comptype == destform && offset==0) {
/**----------------------------------------------------------------------------
    nothing to do, assumes it was properly done by the user */
    if(verbose==1) printf("matrix format conversion from %d toward %d with offset=%d, nothing to be done\n",Mat->comptype,destform,offset);
    return(0);
    }

  if(verbose==1)
    printf("matrix format conversion from %d toward %d with offset=%d (source=%d target=%d)\n",Mat->comptype,destform,offset,Mat->base, base);
  
  switch(Mat->comptype) {
    
    case COO:
      switch(destform) {
        case CSC: 
          coo2csc(Mat,offset);
          Mat->comptype = destform;
          status=0;
          break;
        case CSR: /* To csr */
          printf("(double matrix) conversion from coo to csr not yet implemented\n");
          break;
        case BND: /* To bnd */
          printf("(double matrix) conversion from coo to bnd not yet implemented\n");
          break;
        case COO: /* To csr */
          printf("(double matrix) conversion from coo to coo not yet implemented\n");
          break;
        case COO_COL_ARRANGED:
          printf("(double matrix) conversion from coo to coo (COL arranged) not yet implemented, offset=%d\n",offset);
          break;
        default:
          printf("(double matrix) conversion from coo to %d not yet implemented\n",destform);
          break;
        }
      break;
      
    case CSC:
      switch(destform) {
        case COO:
          if (Mat->stype >  0) csc2coosym(Mat, offset);
          if (Mat->stype == 0) status=csc2coo(Mat, offset);
          if(status!=0) return(status);
          Mat->comptype = destform;
          status=0;
          break;
        case COO_COL_ARRANGED: /* To coo for pastix */
          if(verbose) printf("complex matrix conversion from csc to coo\n");
/*------------------------------------------------------------------------------
          poc-solvers refurbishing : bug found */
//           csr2cooz(Mat, offset);
          csr2coo(Mat, offset);
          Mat->comptype = destform;
          status=0;
          break;
        case BND:
          status = csc2bnd(Mat, offset);
          Mat->comptype = destform;
          status=0;
          break;
        default:
          printf("(double matrix) conversion from csc to %d not yet implemented\n",destform);
          break;
        }
      break;
      
    case CSR: /* *** From csr *** */
        switch(destform) {   
          case COO: /* To coo */
            if(verbose) printf("complex matrix conversion from csr to coo\n");
            if (Mat->stype >  0) csc2coosym(Mat, offset);
            else if (Mat->stype == 0) csr2coo(Mat, offset);
            Mat->comptype = destform;
            status=0;
            break;
          default:
            printf("(complex) matrix conversion from csc to %d not yet implemented\n",destform);
            break;
          }
      break;
      
    case BND: /* *** From bnd *** */
      printf("(double matrix) conversion from BND not yet implemented\n",destform);
      break;
      
    default : /*Not yet*/
      printf("(double matrix) matrix format not yet supported \n");
      break;
    }
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int csc2coosym(triplet *Mat, int offset) {
 
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int status=0;
  int k,m,count,count2;
  int *lig,*col;
  double *val;

  count = 0;
  count2 =0;
  /* find nnz sym values*/
  for (m=0;m<Mat->nrow;m++) {
    for (k=Mat->i[m];k<Mat->i[m+1];k++) {
      if (Mat->j[count] <= m)  {
        count2++;
        }
      count ++;
      }
    }
  
  lig = (int *) malloc (count2 * sizeof (int)) ;
  col = (int *) malloc (count2 * sizeof (int)) ;
  val = (double *) malloc (count2 * sizeof (double)) ;
  count = 0;
  count2 =0;
  for (m=0;m<Mat->nrow;m++) {
    for (k=Mat->i[m];k<Mat->i[m+1];k++) {
      if (Mat->j[count] <= m)  {
	lig[count2] = m+offset;
	col[count2]=Mat->j[count]+offset;
	val[count2] = Mat->x[count];
	count2++;
        }
      count ++;
      }
    }
  Mat->nnz = count2;

  free(Mat->j);
  free(Mat->i);
  free(Mat->x);

  Mat->i = lig;
  Mat->j = col;
  Mat->x = val;

  return(status);
}

// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   int csc2coo(triplet *Mat, int offset)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
// /**@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//     
//   ordering array conversion, matrix untouched
//   
//   connexion with initial matrix structure (as in T-UGOm and tools) needs
//   to be carefuly handled
//   
//   Documentation:
//     http://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.coo_matrix.html
//     
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
// {
//   int status=0;
//   int i,j,k,l,m,n;
//   size_t count,count2;
//   int *lig,*col;
//   char *c;
// 
//   count  = 0;
//   count2 = Mat->nnz;
//   
//   lig = (int *)    malloc (count2 * sizeof (int));
//   col = (int *)    malloc (count2 * sizeof (int));
//   
//   count  = 0;
//   count2 = 0;
//   
//   for (m=0;m<Mat->nrow;m++) {
//     for (k=Mat->i[m];k<Mat->i[m+1];k++) {
//       lig[count2] = m+offset;          /// row incidence array
//       col[count2] = Mat->j[k]+offset;  /// column incidence array
//       count2++;
//       count ++;
//       }
//     }
// 
//   free(Mat->j);
// /**----------------------------------------------------------------------------
//   now well handled, re-activated */
//   free(Mat->i);
//     
//   Mat->i = lig;
//   Mat->j = col;
//   
//   Mat->base+=offset;
// 
//   return(status);
// }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int csr2coo(triplet *Mat, int offset)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

/**@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    
  ordering array conversion, matrix untouched

  csr is close to coo, matrix arrangement are the same (thus untouched)

  connexion with initial matrix structure (as in T-UGOm and tools) needs
  to be carefuly handled
  
  Documentation:
    http://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.coo_matrix.html
    
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
{
  int status=0;
  int i,j,k,l,m,n,count,count2;
  int *lig,*col;
  complex<double> *val;
  char *c;

  count = 0;
  count2 = Mat->nnz;
  
  int base=Mat->base;
  
  lig = (int *) malloc (count2 * sizeof (int)) ;
  col = (int *) malloc (count2 * sizeof (int)) ;
  
  count = 0;
  count2 =0;
  for (m=0;m<Mat->nrow;m++) {
    for (k=Mat->i[m]-base;k<Mat->i[m+1]-base;k++) {
      lig[count2] = m+base+offset;
      col[count2] = Mat->j[k]+offset;
      count2++;
      count ++;
      }
    }
  free(Mat->j);
  free(Mat->i);

  Mat->i = lig;
  Mat->j = col;

  Mat->base+=offset;

  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int csc2coo(triplet *Mat, int offset)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

/**@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    
  ordering array conversion:
  
  csc is pretty different from coo, matrix re-arrangement necessary

  connexion with initial matrix structure (as in T-UGOm and tools) needs
  to be carefuly handled
  
  Documentation:
    http://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.coo_matrix.html
    
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
{
  int status=0;
  int i,j,l,n,row,col,count;
  int *line,*column,*cardinal;
  size_t *pointer, pos, k, nnz;
  double *values;
  char *c;

  count = 0;
  nnz = Mat->nnz;
  
  int base=Mat->base;
  
  line   = (int *) malloc (nnz * sizeof (int)) ;
  column = (int *) malloc (nnz * sizeof (int)) ;
  
  for(k=0;k<nnz;k++) {
    line[k]=column[k]=-1;
    }
  
  values = (double *) malloc (nnz * sizeof (double)) ;
  
  cardinal= (int *) malloc (Mat->nrow * sizeof (int)) ;
  
  for (row=0;row<Mat->nrow;row++) cardinal[row]=0;
    
  count = 0;
  for (col=0;col<Mat->ncol;col++) {
    for (k=Mat->i[col]-base;k<Mat->i[col+1]-base;k++) {
      row=Mat->j[k]-base;
      cardinal[row]++;
      }
    }
    
  pointer= (size_t *) malloc ((Mat->nrow+1) * sizeof (size_t));
  
  pointer[0]=0;
  for (row=0;row<Mat->nrow;row++) pointer[row+1]=pointer[row]+cardinal[row];
  
  for (col=0;col<Mat->ncol;col++) {
    for (k=Mat->i[col]-base;k<Mat->i[col+1]-base;k++) {
      row=Mat->j[k]-base;
      if(row>Mat->nrow-1) {
        printf("#%s: bad!\n",__FUNCTION__);
        return(-1);
        }
      pos=pointer[row];
      if(pos>Mat->nnz-1) {
        printf("#%s: bad!\n",__FUNCTION__);
        return(-1);
        }
      line[pos]   = row+base+offset;
      column[pos] = col+base+offset;
      values[pos] = Mat->x[k];
      pointer[row]++;
      count ++;
      }
    }
    
  for(k=0;k<nnz;k++) {
    if(line[k]==-1) {
      printf("#%s: bad!\n",__FUNCTION__);
      return(-1);
      }
    if(column[k]==-1) {
      printf("#%s: bad!\n",__FUNCTION__);
      return(-1);
      }
    }
    
  free(Mat->j);
  free(Mat->i);
  free(Mat->x);
  
  Mat->i = line;
  Mat->j = column;
  Mat->x = values;

  free(cardinal);
  free(pointer);
  
  Mat->base+=offset;

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int coo2csc(triplet *Mat,int offset) {

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  int status=0,bcl,k,p;
  int *Ap, *Ai, *w;
  double *Ax;
  char *c;
  /* to check diff with umfpack triplet2col */
  /*int *Ap2, *Ai2;*/
  /*double *Ax2;*/

  Ap = (int *) malloc ((Mat->nrow+1) * sizeof (int)) ;
  w  = (int *) malloc ((Mat->nrow+1) * sizeof (int)) ;
  Ai = (int *) malloc (Mat->nnz * sizeof (int)) ;
  Ax = (double *) malloc (Mat->nnz * sizeof (double)) ; 
  
  for (bcl=0;bcl <  Mat->nnz;bcl++) {
    /*printf("Mat->i[%d]= %d,Mat->j[%d]= %d \n",bcl,Mat->i[bcl],bcl,Mat->j[bcl]);*/
    Mat->i[bcl]= Mat->i[bcl]+offset;
    Mat->j[bcl]= Mat->j[bcl]+offset;
    /*printf("Mat->i[%d]= %d,Mat->j[%d]= %d \n",bcl,Mat->i[bcl],bcl,Mat->j[bcl]);*/
    /*c = getchar();*/
  }
  for (k = 0 ; k < Mat->nrow+1 ; k++) w[k]=0.0;
  /* (k = 0 ; k < Mat->nnz ; k++) printf("i[%d]=%d\n",k,Mat->i[k]) ;  */
  for (k = 0 ; k < Mat->nnz ; k++) w[Mat->i[k]]++ ;  
  /* (k = 0 ; k < Mat->nrow ; k++) printf("w[%d]=%d\n",k,w[k]) ;  */
  cs_cumsum (Ap, w, Mat->nrow) ;
  p=0;
  for (k = 0 ; k < Mat->nnz ; k++)
    if(Mat->x[k] == 0) p++;
  for (k = 0 ; k < Mat->nnz ; k++)
    {
      Ai [p = w [Mat->i[k]]++] = Mat->j[k] ;    /* A(i,j) is the pth entry in C */
      Ax [p] = Mat->x[k] ;
    }
  /* Few difference between umf compress and cs */
/*   Ap2= (int *) malloc ((Mat->nrow+1) * sizeof (int)) ; */
/*   Ai2 = (int *) malloc (Mat->nnz * sizeof (int)) ; */
/*   Ax2 = (double *) malloc (Mat->nnz * sizeof (double)) ;  */
/*   status = umfpack_di_triplet_to_col(Mat->nrow,Mat->nrow , Mat->nnz,   */
/*  				      Mat->i,  Mat->j, Mat->x,  */
/*  				     Ap2, Ai2, Ax2, (int *) NULL) ;  */
/*   printf("umfpack_di_triplet_to_col status=%d\n",status);  */
/*   printf("bcl, Ap , Ap2, diff\n"); */
/*   for (k=0;k<Mat->nrow+1 ; k++) { */
/*     if ((Ap[k]-Ap2[k]) != 0) { */
/*       printf("%d : %d  %d %d\n",k,Ap[k],Ap2[k],Ap[k]-Ap2[k]); */
/*       printf("   CS:%d, %e  ||  UM:%d, %e \n",Ai[Ap[k]],Ax[Ap[k]],Ai2[Ap2[k]],Ax2[Ap2[k]]); */
/*       status++; */
/*       c= getchar(); */
/*     } */
/*   } */
/*   printf("nbdiff =%d\n",status); */
  status=0;
  free(Mat->i);
  free(Mat->j);
  free(Mat->x);
  Mat->i=Ap;
  Mat->j=Ai;
  Mat->x=Ax; 
  Mat->comptype=1;
  Mat->base = Mat->base+offset;

  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int getbwd(triplet *Mat,int *ml, int *mu) {

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  int i,k,ldist,status=0;
  *ml = - Mat->nrow;
  *mu = - Mat->nrow;
  /*for(i=0;i < Mat->nrow; i++) printf("Ap==i = %d \n",Mat->i[i]);*/
  /*for(i=0;i < Mat->nnz; i++) printf("Ai==j = %d \n",Mat->j[i]);*/
  for(i = 0;i < Mat->nrow;i++) { /* i == colonne*/
    for (k = Mat->i[i];k< Mat->i[i+1];k++) { 
      /*Mat->j[k] == Ligne */
      ldist=i-Mat->j[k]; 
      *mu = MAX(*mu,ldist);
      *ml = MAX(*ml,-ldist);    
    }
  }
  /*printf("ml=%d,mu=%d \n",*ml,*mu);*/
  return(status);
}

/*-------------------------------------------------------------*/
/*------------convert csc matric to lapack banded format ------*/
/*-------------------------------------------------------------*/

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int csc2bnd(triplet *Mat,int offset) {

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  int m,ml,mu,mdiag,lowd,n,dg;
  int k,lda,count=0,status=0;
  status=getbwd(Mat,&ml,&mu);
  m = 2*ml+mu+1;
  dg =  ml+mu; /*pb if ml != mu*/
  lda= 2*ml+mu+1;
  Mat->a = (double *) malloc ((Mat->nrow * lda) * sizeof (double)) ;
  
  for(n=0; n < (Mat->nrow * lda) ;n++) {
    Mat->a[n]=0;
  }
  for(m = 0; m < Mat->nrow; m++) {  /*row index */
    for(k = Mat->i[m]; k < Mat->i[m+1]; k++) {
      n = Mat->j[k];     /*column index */
      /*Mat->a[lda * n + m - n + dg] = Mat->x[k];*/
      /* Mat->a[lda * m + n - m + dg] = Mat->x[k];*/
      Mat->a[lda * n + m - n + dg] = Mat->x[k];
    }
  }
  free(Mat->x);
  free(Mat->j);
  free(Mat->i);
  count=0;
  Mat->comptype=3;
  Mat->lda=lda;
  Mat->ml=ml;
  Mat->mu=mu;

  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int coo_issym(triplet *Mat) {

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  /* methode brutale */
  /* return 1 si Mat est symetrique 0 sinon*/
  int bcl1,bcl2,li,co,count,ok=1,issym=1;
  double val;
  
  for(bcl1=0;bcl1 < Mat->nnz;bcl1++) {
    li  = Mat->i[bcl1];
    co  = Mat->j[bcl1];
    val = Mat->x[bcl1];
    if (li != co) {
      count = bcl1;
      ok=0;
      for(bcl2=0;bcl2 <Mat->nnz;bcl2++) {
	/*printf("      bcl2=%d    %d,%d,%e %e\n",
	  bcl2, Mat->i[bcl2], Mat->j[bcl2], Mat->x[bcl2],Mat->x[bcl2]-val);*/
	if (li == Mat->j[bcl2]) {
	  if (co == Mat->i[bcl2]) {
	    if(val == Mat->x[bcl2]) {
	      ok++;
	      bcl2 = Mat->nnz;
	      /*printf("      elt sym ok\n");*/
	    }
	  }
	}
      }
    }
    if (ok == 0) bcl1= Mat->nnz;
  }
  return(ok);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int csc_issym(triplet *Mat) {

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  /* methode brutale */
  /* return 1 si Mat est symetrique 0 sinon*/
  int bcl1,bcl2,li,co,count,ok=1,issym=1;
  double val;
  
}
