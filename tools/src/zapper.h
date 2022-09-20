
/**************************************************************************

  T-UGOm hydrodynamic ocean model, 2006

  Unstructured Ocean Grid initiative

Contributors:

  Florent Lyard      LEGOS/CNRS, Toulouse, France
  Laurent Roblou     LEGOS/CNRS, Toulouse, France
  Yoann Le Bars      LEGOS, Toulouse, France
  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
  Frédéric Dupont    Université de Laval à Québec, Canada

***************************************************************************/

#ifndef ZAPPER_H
#   define ZAPPER_H
#   include <assert.h>
template < typename T > inline void zapsc(T & x);
template < typename T > inline void zaparr(T & x);
// Put an assert to check if x is NULL, this is to catch
// program "logic" errors early. Even though delete works
// fine with NULL by using assert you are actually catching
// "bad code" very early

// Defining Zap using templates
// Use zap instead of delete as this will be very clean

template < class T > inline void zapsc(T & x)
{
  assert(x != NULL);
  delete x;
  x = NULL;
}

// In C++ the reason there are 2 forms of the delete operator is - because
// there is no way for C++ to tell the difference between a pointer to
// an object and a pointer to an array of objects. The delete operator
// relies on the programmer using "[]" to tell the two apart.
// Hence, we need to define zaparr function below.
// To delete array of pointers

template < class T > inline void zaparr(T & x)
{
  if(x==NULL) {
/**----------------------------------------------------------------------------
    error*/
    printf("zappar fatal error\n");
    exit(0);
    }
  assert(x != NULL);
  delete[]x;
  x = NULL;
}

//The zap() function will delete the pointer and set it NULL. This will
// ensure that even if multiple zap()'s are called
// on the same deleted pointer then the program will not crash.
//        zap(pFirstname);
//        //zap(pFirstname); // no core dumps.  Because pFirstname is NULL now
//        //zap(pFirstname); // no core dumps.  Because pFirstname is NULL now
//
//        zap(pLastname);
//        zap(pJobDescription);
//
//        int *iiarray = new int[10];
//        zaparr(iiarray);

#endif /* ZAPPER_H */
