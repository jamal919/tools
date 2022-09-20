
/**************************************************************************

  T-UGOm hydrodynamic ocean model & T-UGO tools, 2006-2012

  Unstructured Ocean Grid initiative

Contributors:

  Florent Lyard      LEGOS/CNRS, Toulouse, France
  Cyril Nguyen       LA/CNRS,    Toulouse, France
  Laurent Roblou     LEGOS/CNRS, Toulouse, France
  Damien Allain      LEGOS/CNRS, Toulouse, France
  
  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
  Frédéric Dupont    Université de Laval à Québec, Canada
  
  Yoann Le Bars      PhD, LEGOS, Toulouse, France
  Yves Soufflet      Post-doctorant, LEGOS, Toulouse, France
  Clement Mayet      PhD, LEGOS, Toulouse, France

***************************************************************************/

#include "tools-structures.h"
#include "fe.h"
#include "matrix.h"
#include "zapper.h"

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int init_packed_CQP0xCQP0(mesh_t & mesh, discretisation_t & descriptor, ordering_t *ordering)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**----------------------------------------------------------------------
  Compressed Row storage/Compressed Sparse Row
  Packed matrix with (reordered) neighbours list as row index

  !!! in our case (geometryic matrix symetry), equivalent to

  Compressed Column storage/Compressed Sparse Column/Harwell-Boing Sparse Matrix
  Packed matrix with (reordered) neighbours list as col index

  QP0xQP0, diagonal

  ----------------------------------------------------------------------*/
{
  int m;
  int nelts;

  nelts = mesh.nquadrangles;

/*----------------------------------------------------------------------
  allocate arrays*/
  ordering->cardinal  = new int[nelts];
  ordering->pointer   = new int[nelts + 1];
  ordering->incidence = new int[nelts];

/*----------------------------------------------------------------------
  */
  for(m = 0; m < nelts; m++) {
    (ordering->cardinal)[m]  = 1;
    (ordering->pointer)[m]   = m;
    (ordering->incidence)[m] = m;
    }
  (ordering->pointer)[m] = nelts;

  ordering->nrows = nelts;
  return (0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int init_packed_CQP0xCQP0_01(mesh_t & mesh, discretisation_t & descriptor, ordering_t *ordering)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**----------------------------------------------------------------------
  Compressed Row storage/Compressed Sparse Row
  Packed matrix with (reordered) neighbours list as row index

  !!! in our case (geometryic matrix symetry), equivalent to

  Compressed Column storage/Compressed Sparse Column/Harwell-Boing Sparse Matrix
  Packed matrix with (reordered) neighbours list as col index

  QP0xQP0, with 1 

  ----------------------------------------------------------------------*/
{
  int k,m,status;
  int ncols, ndof;
  int *tmp;
  size_t n, count;

  ndof = mesh.nquadrangles;
  ncols = 5;

/*----------------------------------------------------------------------
  allocate arrays*/
  ordering->cardinal  = new int[ndof];
  ordering->pointer   = new int[ndof + 1];

  tmp = new int[ndof * ncols];
  
  for(n = 0; n < ndof * ncols; n++)
    tmp[n] = -1;

/*----------------------------------------------------------------------
  first level of neighbours*/
  count = 0;
  for(n = 0; n < ndof; n++) {
    ordering->pointer[n]  = count;
    tmp[n*ncols]=n;
    ordering->cardinal[n]=1;
    count+= 1;
    for(k = 0; k < 4; k++) {
      size_t nn=mesh.quadrangles[n].edges[k];
      if (mesh.edges[nn].nshared==1) continue;
      size_t m1=mesh.edges[nn].shared[0];
      size_t m2=mesh.edges[nn].shared[1];
      if (m1==n) {
        tmp[n*ncols + ordering->cardinal[n]] = m2;
        }
      else {
        tmp[n*ncols + ordering->cardinal[n]] = m1;
        }
      ordering->cardinal[n]++;
      count++;
      }
    }
    
  ordering->pointer[n] = count;
  ordering->incidence = new int[count];

  for(n = 0; n < ndof; n++) {
    for(k = 0; k < ordering->cardinal[n]; k++) {
      size_t pos=ordering->pointer[n] + k;
      ordering->incidence[pos] = tmp[n * ncols + k];
      }
    }
  delete[] tmp;
  
  ordering->nrows = ndof;
  status= matrix_reorder(ordering);

  status=matrix_check_consistency(ordering->incidence,ordering->pointer,ordering->cardinal,ndof);
  return (0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int init_packed_CQP0xCQP0_02(mesh_t & mesh, discretisation_t & descriptor, ordering_t *ordering)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**----------------------------------------------------------------------
  Compressed Row storage/Compressed Sparse Row
  Packed matrix with (reordered) neighbours list as row index

  !!! in our case (geometryic matrix symetry), equivalent to

  Compressed Column storage/Compressed Sparse Column/Harwell-Boing Sparse Matrix
  Packed matrix with (reordered) neighbours list as col index

  QP0xQP0, with 1 

  ----------------------------------------------------------------------*/
{
  int k,m,status;
  int ncols, ndof;
  int *tmp;
  size_t n, count;

  ndof = mesh.nquadrangles;
  ncols = 9;

/*----------------------------------------------------------------------
  allocate arrays*/
  ordering->cardinal  = new int[ndof];
  ordering->pointer   = new int[ndof + 1];

  tmp = new int[ndof * ncols];
  
  for(n = 0; n < ndof * ncols; n++)
    tmp[n] = -1;

/*----------------------------------------------------------------------
  first level of neighbours*/
  count = 0;
  for(n = 0; n < ndof; n++) {
    ordering->pointer[n]  = count;
    tmp[n*ncols]=n;
    ordering->cardinal[n]=1;
    count+= 1;
    for(k = 0; k < 4; k++) {
      size_t nn=mesh.quadrangles[n].edges[k];
      if (mesh.edges[nn].nshared==1) continue;
      size_t m1=mesh.edges[nn].shared[0];
      size_t m2=mesh.edges[nn].shared[1];
      if (m1==n) {
        tmp[n*ncols + ordering->cardinal[n]] = m2;
        }
      else {
        tmp[n*ncols + ordering->cardinal[n]] = m1;
        }
      ordering->cardinal[n]++;
      count++;
      }
    }
    
  ordering->pointer[n] = count;
  ordering->incidence = new int[count];

  for(n = 0; n < ndof; n++) {
    for(k = 0; k < ordering->cardinal[n]; k++) {
      size_t pos=ordering->pointer[n] + k;
      ordering->incidence[pos] = tmp[n * ncols + k];
      }
    }
  delete[] tmp;
  
  ordering->nrows = ndof;
  status= matrix_reorder(ordering);

  status=matrix_check_consistency(ordering->incidence,ordering->pointer,ordering->cardinal,ndof);
  return (0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int init_packed_CQN1xCQP0_01(mesh_t & mesh, discretisation_t & descriptor, ordering_t *ordering)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**----------------------------------------------------------------------------
  Compressed Row storage/Compressed Sparse Row
  Packed matrix with (reordered) neighbours list as row index

  !!! in our case (geometryic matrix symetry), equivalent to

  Compressed Column storage/Compressed Sparse Column/Harwell-Boing Sparse Matrix
  Packed matrix with (reordered) neighbours list as col index

  CQN1xCQP0

-----------------------------------------------------------------------------*/
{
  int k,m,status;
  int ncols, ndof;
  int *tmp;
  size_t n, count;

  ndof = descriptor.nnodes;
  ncols = 2;

/*----------------------------------------------------------------------
  allocate arrays*/
  ordering->cardinal  = new int[ndof];
  ordering->pointer   = new int[ndof + 1];

  tmp = new int[ndof * ncols];
  
  for(n = 0; n < ndof * ncols; n++)
    tmp[n] = -1;

/*----------------------------------------------------------------------
  first level of neighbours*/
  count = 0;
  for(n = 0; n < ndof; n++) {
    ordering->pointer[n]  = count;
    ordering->cardinal[n] = mesh.edges[n].nshared;
    count += ordering->cardinal[n];
    for(k = 0; k < mesh.edges[n].nshared; k++) {
      tmp[n*ncols + k] = mesh.edges[n].shared[k];
      }
    }
    
  ordering->pointer[n] = count;
  ordering->incidence = new int[count];

  for(n = 0; n < ndof; n++) {
    for(k = 0; k < ordering->cardinal[n]; k++) {
      size_t pos=ordering->pointer[n] + k;
      ordering->incidence[pos] = tmp[n * ncols + k];
//       if(ordering->incidence[pos] == -1) {
//         check_error(-1, "matrix index anomaly", __LINE__, __FILE__, 1);
//         }
      }
    }
  delete[] tmp;
  
  ordering->nrows = ndof;
  status= matrix_reorder(ordering);

  status=matrix_check_consistency(ordering->incidence,ordering->pointer,ordering->cardinal,ndof);
  return (0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int init_packed_CQN1xCQP0_02(mesh_t & mesh, discretisation_t & descriptor, ordering_t *ordering)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**----------------------------------------------------------------------------
  Compressed Row storage/Compressed Sparse Row
  Packed matrix with (reordered) neighbours list as row index

  !!! in our case (geometryic matrix symetry), equivalent to

  Compressed Column storage/Compressed Sparse Column/Harwell-Boing Sparse Matrix
  Packed matrix with (reordered) neighbours list as col index

  CQN1xCQP0

-----------------------------------------------------------------------------*/
{
  int k,l,m,status;
  int n1,n2;
  int ncols, ndof;
  int *tmp;
  size_t n, count;

  ndof  = descriptor.nnodes;
  ncols = 6;

/*----------------------------------------------------------------------
  allocate arrays*/
  ordering->cardinal  = new int[ndof];
  ordering->pointer   = new int[ndof + 1];

  tmp = new int[ndof * ncols];
  
  for(n = 0; n < ndof * ncols; n++)
    tmp[n] = -1;

/*----------------------------------------------------------------------
  first level of neighbours*/
  count = 0;
  for(n = 0; n < ndof; n++) {
    ordering->pointer[n]  = count;
    n1= mesh.edges[n].extremity[0];
    n2= mesh.edges[n].extremity[1];
    ordering->cardinal[n] = mesh.vertices[n1].nelmts+mesh.vertices[n2].nelmts-mesh.edges[n].nshared;
    l=0;
    for(k = 0; k < mesh.vertices[n1].nelmts; k++) {
      tmp[n*ncols + l] = mesh.vertices[n1].elmts[k];
      l++;
      }
    for(k = 0; k < mesh.vertices[n2].nelmts; k++) {
      m=mesh.vertices[n2].elmts[k];
      if( vpos(m,&tmp[n*ncols],l)==-1) {
        tmp[n*ncols + l] = m;
        l++;
        }
      }
    if(l!=ordering->cardinal[n]) {
      printf("init_packed_CQN1xCQP0, warning\n");
      ordering->cardinal[n]=l;
      }
    count += ordering->cardinal[n];
    }
    
  ordering->pointer[n] = count;
  ordering->incidence = new int[count];

  for(n = 0; n < ndof; n++) {
    for(k = 0; k < ordering->cardinal[n]; k++) {
      size_t pos=ordering->pointer[n] + k;
      ordering->incidence[pos] = tmp[n * ncols + k];
//       if(ordering->incidence[pos] == -1) {
//         check_error(-1, "matrix index anomaly", __LINE__, __FILE__, 1);
//         }
      }
    }
  delete[] tmp;
  
  ordering->nrows = ndof;
  status= matrix_reorder(ordering);

  status=matrix_check_consistency(ordering->incidence,ordering->pointer,ordering->cardinal,ndof);
  return (0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int init_packed_CQN1xCQP0_02(mesh_t & mesh, discretisation_t & descriptor, ordering_t *ordering, int depth)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**----------------------------------------------------------------------
  Compressed Row storage/Compressed Sparse Row
  Packed matrix with (reordered) neighbours list as row index

  !!! in our case (geometryic matrix symetry), equivalent to

  Compressed Column storage/Compressed Sparse Column/Harwell-Boing Sparse Matrix
  Packed matrix with (reordered) neighbours list as col index

  ncP1 (rows) x LGP0 (columns) based on 2 level neighbouring relationships

  ----------------------------------------------------------------------*/
{
  int status;
  int i, m, mm, n, j, k, l;
  int row, nrows, col, ncols, loop, first, last, additional;
  int *tmp;
  int *incidence, *cardinal, *pointer;

  size_t count;

  nrows = 2;
  ncols = 6;
  additional=4;
  for(loop=0;loop<depth;loop++) {
    ncols+=additional;
    additional=additional*2;
    }

/*----------------------------------------------------------------------
  allocate arrays*/
  cardinal = new int[mesh.nedges];
  pointer  = new int[mesh.nedges + 1];
  tmp      = new int[mesh.nedges * ncols];

/*----------------------------------------------------------------------
  node index in packed matrix*/
  for(n = 0; n < mesh.nedges * ncols; n++)
    tmp[n] = -1;

/*----------------------------------------------------------------------
  first level of neighbours*/
  count = 0;
  for(n = 0; n < mesh.nedges; n++) {
    pointer[n]  = count;      
    cardinal[n] = mesh.edges[n].nshared;
    for(k = 0; k < mesh.edges[n].nshared; k++) {
      tmp[n*ncols + k] = mesh.edges[n].shared[k];
      }
    first=0;
    for(loop=0;loop<depth;loop++) {
      last=cardinal[n];
      for(k = first; k < last; k++) {
        m=tmp[n*ncols + k];
        for(l = 0; l < 4; l++) {
          int nn=mesh.quadrangles[m].edges[l];
	  for(i=0;i<mesh.edges[nn].nshared;i++) {
            mm=mesh.edges[nn].shared[i];
            if(vpos(mm,&(tmp[n*ncols]),cardinal[n])==-1) {
              tmp[n*ncols + cardinal[n]] = mm;
              if(cardinal[n]==ncols) {
                check_error(-1, "matrix index overflow", __LINE__, __FILE__, 1);
               }
              cardinal[n]++;
              }
            }
          }
        }
      first=last;
      }
    count += cardinal[n];
    }
  pointer[n] = count;

  incidence = new int[count];

  for(n = 0; n < mesh.nedges; n++) {
    for(k = 0; k < cardinal[n]; k++) {
      incidence[pointer[n] + k] = tmp[n * ncols + k];
      if(incidence[pointer[n] + k] == -1) {
        check_error(-1, "matrix index anomaly", __LINE__, __FILE__, 1);
        }
      }
    }

  zaparr(tmp);

  ordering->pointer   = pointer;
  ordering->incidence = incidence;
  ordering->cardinal  = cardinal;
  ordering->nrows     = mesh.nedges;
  
  status=matrix_reorder(ordering);

  status=matrix_check_consistency(ordering->incidence,ordering->pointer,ordering->cardinal,mesh.nedges);


  return (0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int init_packed_CQN1xCQP0(mesh_t & mesh, discretisation_t & descriptor, ordering_t *ordering)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
//  status=init_packed_CQN1xCQP0_02(mesh, descriptor, ordering);
  status=init_packed_CQN1xCQP0_02(mesh, descriptor, ordering,1);
  
  return (status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int init_packed_CQP1xCQP1(mesh_t & mesh, discretisation_t & descriptor, ordering_t *ordering)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**----------------------------------------------------------------------------
  Compressed Row storage/Compressed Sparse Row
  Packed matrix with (reordered) neighbours list as row index

  !!! in our case (geometryic matrix symetry), equivalent to

  Compressed Column storage/Compressed Sparse Column/Harwell-Boing Sparse Matrix
  Packed matrix with (reordered) neighbours list as col index

  QP1xQP0

-----------------------------------------------------------------------------*/
{
  int k,m,status;
  int ncols, ndof;
  int *tmp;
  size_t n, count;

  ndof  = descriptor.nnodes;
  ncols = 9;

/*----------------------------------------------------------------------
  allocate arrays*/
  ordering->cardinal  = new int[ndof];
  ordering->pointer   = new int[ndof + 1];

  tmp = new int[ndof * ncols];
  
  for(n = 0; n < ndof * ncols; n++)
    tmp[n] = -1;

/*----------------------------------------------------------------------
  first level of neighbours*/
  count = 0;
  for(n = 0; n < ndof; n++) {
    ordering->pointer[n]  = count;
    ordering->cardinal[n] = descriptor.nodes[n].nnghbs+1;
    count += ordering->cardinal[n];
    tmp[n*ncols] = n;
    for(k = 0; k < descriptor.nodes[n].nnghbs; k++) {
      tmp[n*ncols + k + 1] = descriptor.nodes[n].nghbs[k];
      }
    }
    
  ordering->pointer[n] = count;
  ordering->incidence = new int[count];

  for(n = 0; n < ndof; n++) {
    for(k = 0; k < ordering->cardinal[n]; k++) {
      size_t pos=ordering->pointer[n] + k;
      ordering->incidence[pos] = tmp[n * ncols + k];
      if(ordering->incidence[pos] == -1) {
        check_error(-1, "matrix index anomaly", __LINE__, __FILE__, 1);
        }
      }
    }
  delete[] tmp;

  ordering->nrows = ndof;
  
  status= matrix_reorder(ordering);

  status=matrix_check_consistency(ordering->incidence,ordering->pointer,ordering->cardinal,ndof);
  status=matrix_check_symmetry(ordering->incidence,ordering->pointer,ordering->cardinal,ndof);

  return (0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int init_packed_CQN1xCQN1(mesh_t & mesh, discretisation_t & descriptor, ordering_t *ordering)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**----------------------------------------------------------------------------
  Compressed Row storage/Compressed Sparse Row
  Packed matrix with (reordered) neighbours list as row index

  !!! in our case (geometryic matrix symetry), equivalent to

  Compressed Column storage/Compressed Sparse Column/Harwell-Boing Sparse Matrix
  Packed matrix with (reordered) neighbours list as col index

  QP1xQP0

-----------------------------------------------------------------------------*/
{
  int k,m,status;
  int ncols, ndof;
  int *tmp;
  size_t n, count;

  ndof  = descriptor.nnodes;
  ncols = 1;

/*----------------------------------------------------------------------
  allocate arrays*/
  ordering->cardinal  = new int[ndof];
  ordering->pointer   = new int[ndof + 1];

  tmp = new int[ndof * ncols];
  
  for(n = 0; n < ndof * ncols; n++)
    tmp[n] = -1;

/*----------------------------------------------------------------------
  first level of neighbours*/
  count = 0;
  for(n = 0; n < ndof; n++) {
    ordering->pointer[n]  = count;
    ordering->cardinal[n] = 1;
    count += ordering->cardinal[n];
    tmp[n*ncols] = n;
    }
    
  ordering->pointer[n] = count;
  ordering->incidence = new int[count];

  for(n = 0; n < ndof; n++) {
    for(k = 0; k < ordering->cardinal[n]; k++) {
      size_t pos=ordering->pointer[n] + k;
      ordering->incidence[pos] = tmp[n * ncols + k];
      if(ordering->incidence[pos] == -1) {
        check_error(-1, "matrix index anomaly", __LINE__, __FILE__, 1);
        }
      }
    }
  delete[] tmp;

  ordering->nrows = ndof;
  
  status= matrix_reorder(ordering);

  status=matrix_check_consistency(ordering->incidence,ordering->pointer,ordering->cardinal,ndof);
  status=matrix_check_symmetry(ordering->incidence,ordering->pointer,ordering->cardinal,ndof);

  return (0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int init_packed_CQP0xCQP1(mesh_t & mesh, discretisation_t & descriptor, ordering_t *ordering)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**----------------------------------------------------------------------
  Compressed Row storage/Compressed Sparse Row
  Packed matrix with (reordered) neighbours list as row index

  !!! in our case (geometryic matrix symetry), equivalent to

  Compressed Column storage/Compressed Sparse Column/Harwell-Boing Sparse Matrix
  Packed matrix with (reordered) neighbours list as col index

  LGP0 (rows) x LGP1 (column) based on first level neighbouring relationships

  ----------------------------------------------------------------------*/
{
  int i, m, n, n1, n2, n3, j, k, l;
  int status;
  int opposite[2], row, col, ncols;
  double h;
  quadrangle_t *elt;
  int *tmp;
  int *incidence, *cardinal, *pointer;

  size_t count;

  ncols = 4;

  elt = mesh.quadrangles;

/*----------------------------------------------------------------------
  allocate arrays*/
  cardinal = new int[mesh.nquadrangles];
  pointer  = new int[mesh.nquadrangles + 1];
  tmp      = new int[mesh.nquadrangles * ncols];

/*----------------------------------------------------------------------
  node index in packed matrix*/
  for(n = 0; n < mesh.nquadrangles * ncols; n++)
    tmp[n] = -1;

/*----------------------------------------------------------------------
  first level of neighbours*/
  count = 0;
  for(n = 0; n < mesh.nquadrangles; n++) {
    pointer[n] = count;      /*pointer= column start index */
    cardinal[n] = 4;
    count += cardinal[n];
    for(k = 0; k <  4; k++) {
      tmp[n*ncols + k] = mesh.quadrangles[n].vertex[k];
      }
    }
  pointer[n] = count;

  incidence = new int[count];

  for(n = 0; n < mesh.nquadrangles; n++) {
    for(k = 0; k < cardinal[n]; k++) {
       /*incidence= row index */
      incidence[pointer[n] + k] = tmp[n * ncols + k];
      if(incidence[pointer[n] + k] == -1) {
        check_error(-1, "matrix index anomaly", __LINE__, __FILE__, 1);
        }
      }
    }

  zaparr(tmp);

  ordering->pointer   = pointer;
  ordering->incidence = incidence;
  ordering->cardinal  = cardinal;
  ordering->nrows     = mesh.nquadrangles;

  status= matrix_reorder(ordering);

  status=matrix_check_consistency(ordering->incidence,ordering->pointer,ordering->cardinal,mesh.nquadrangles);
  
  return (status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int init_packed_CQP0xCQN1(mesh_t & mesh, discretisation_t & descriptor, ordering_t *ordering)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**----------------------------------------------------------------------
  Compressed Row storage/Compressed Sparse Row
  Packed matrix with (reordered) neighbours list as row index

  !!! in our case (geometryic matrix symetry), equivalent to

  Compressed Column storage/Compressed Sparse Column/Harwell-Boing Sparse Matrix
  Packed matrix with (reordered) neighbours list as col index

  LGP0 (rows) x LGP1 (column) based on first level neighbouring relationships

  ----------------------------------------------------------------------*/
{
  int i, m, n, n1, n2, n3, j, k, l;
  int status;
  int opposite[2], row, col, ncols;
  double h;
  quadrangle_t *elt;
  int *tmp;
  int *incidence, *cardinal, *pointer;

  size_t count;

  ncols = 4;

  elt = mesh.quadrangles;

/*----------------------------------------------------------------------
  allocate arrays*/
  cardinal = new int[mesh.nquadrangles];
  pointer  = new int[mesh.nquadrangles + 1];
  tmp      = new int[mesh.nquadrangles * ncols];

/*----------------------------------------------------------------------
  node index in packed matrix*/
  for(n = 0; n < mesh.nquadrangles * ncols; n++)
    tmp[n] = -1;

/*----------------------------------------------------------------------
  first level of neighbours*/
  count = 0;
  for(n = 0; n < mesh.nquadrangles; n++) {
    pointer[n] = count;      /*pointer= column start index */
    cardinal[n] = 4;
    count += cardinal[n];
    for(k = 0; k <  4; k++) {
      tmp[n*ncols + k] = mesh.quadrangles[n].edges[k];
      }
    }
  pointer[n] = count;

  incidence = new int[count];

  for(n = 0; n < mesh.nquadrangles; n++) {
    for(k = 0; k < cardinal[n]; k++) {
       /*incidence= row index */
      incidence[pointer[n] + k] = tmp[n * ncols + k];
      if(incidence[pointer[n] + k] == -1) {
        check_error(-1, "matrix index anomaly", __LINE__, __FILE__, 1);
        }
      }
    }

  zaparr(tmp);

  ordering->pointer   = pointer;
  ordering->incidence = incidence;
  ordering->cardinal  = cardinal;
  ordering->nrows     = mesh.nquadrangles;

  status= matrix_reorder(ordering);

  status=matrix_check_consistency(ordering->incidence,ordering->pointer,ordering->cardinal,mesh.nquadrangles);
  
  return (status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int init_packed_CQP1xCQP0(mesh_t & mesh, discretisation_t & descriptor, ordering_t *ordering)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**----------------------------------------------------------------------
  Compressed Row storage/Compressed Sparse Row
  Packed matrix with (reordered) neighbours list as row index

  !!! in our case (geometryic matrix symetry), equivalent to

  Compressed Column storage/Compressed Sparse Column/Harwell-Boing Sparse Matrix
  Packed matrix with (reordered) neighbours list as col index

  LGP0 (rows) x LGP1 (column) based on first level neighbouring relationships

  ----------------------------------------------------------------------*/
{
  int i, m, n, n1, n2, n3, j, k, l;
  int status;
  int opposite[2], row, col, ncols;
  double h;
  quadrangle_t *elt;
  int *tmp;
  int *incidence, *cardinal, *pointer;
  size_t ndof;

  size_t count;

  ncols = 4;

  elt = mesh.quadrangles;
  
  ndof=descriptor.nnodes;

/*----------------------------------------------------------------------
  allocate arrays*/
  cardinal = new int[ndof];
  pointer  = new int[ndof + 1];
  tmp      = new int[ndof * ncols];

/*----------------------------------------------------------------------
  node index in packed matrix*/
  for(n = 0; n < ndof * ncols; n++)
    tmp[n] = -1;

/*----------------------------------------------------------------------
  first level of neighbours*/
  count = 0;
  for(n = 0; n < ndof; n++) {
    pointer[n] = count;
    cardinal[n] = 4;
    count += cardinal[n];
    for(k = 0; k <  4; k++) {
      tmp[n*ncols + k] = mesh.edges[n].shared[k];
      }
    }
  pointer[n] = count;

  incidence = new int[count];

  for(n = 0; n < mesh.nquadrangles; n++) {
    for(k = 0; k < cardinal[n]; k++) {
      incidence[pointer[n] + k] = tmp[n * ncols + k];
      if(incidence[pointer[n] + k] == -1) {
        check_error(-1, "matrix index anomaly", __LINE__, __FILE__, 1);
        }
      }
    }

  zaparr(tmp);

  ordering->pointer   = pointer;
  ordering->incidence = incidence;
  ordering->cardinal  = cardinal;
  ordering->nrows     = ndof;

  status= matrix_reorder(ordering);

  status=matrix_check_consistency(ordering->incidence,ordering->pointer,ordering->cardinal,mesh.nquadrangles);
  
  return (status);
}
