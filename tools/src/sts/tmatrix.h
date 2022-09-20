//-----------------------tmatrix prototype and inline function------------

#ifndef TMATRIX_H
#define TMATRIX_H

// Note that the  argument "Ttype" is not used in the subroutine
// It is necessary because C++ only figures out the type from arguments.

/* Matrix arithmetic is very difficult to implement. We want to use an
   efficient and portable library, for instance Boost. */

template < typename T > inline T ** tmatrix(T ** Ttype, const int nrows, const int ncols);

template < typename T > inline T ** t1Darray(T ** matrix, const int ndim);

//--------------------------------------------------------------------
//--------------------------------------------------------------------
template < typename T > inline T ** tmatrix(T **Ttype, const T value, const int nrows,  const int ncols)
{
  T **matrix;

  matrix = new T *[nrows];
  if (!matrix)
    check_error(-1, "memory allocation error", __LINE__, __FILE__, 1);
  for(int i = 0; i < nrows; i++) {
    matrix[i] = new T[ncols];
    if (!matrix[i])
      check_error(-1, "memory allocation error", __LINE__, __FILE__, 1);
    for(int n=0;n<ncols;n++) matrix[i][n]=value;
    }
  return matrix;
}

//--------------------------------------------------------------------   
template < typename T > inline T * t1Darray(T * Ttype, const T value, const int ndim)
{
  T *array;
  int n;
  array = new T[ndim];
  for(n=0;n<ndim;n++) array[n]=value;
  return array;
}

#endif /* TMATRIX_H */
