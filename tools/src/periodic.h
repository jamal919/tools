/***************************************************************************
T-UGO hydrodynamic ocean model, 2006-2012

  Unstructured Ocean Grid initiative

\file

Contributors:

\author  Florent Lyard      LEGOS/CNRS, Toulouse, France
\author  Cyril Nguyen       LA/CNRS,    Toulouse, France
\author  Laurent Roblou     LEGOS, Toulouse, France
\author  Damien Allain      LEGOS/CNRS, Toulouse, France

\author  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
\author  Frédéric Dupont    Université de Laval à Québec, Canada

\author  Yoann Le Bars      PhD, LEGOS, Toulouse, France
\author  Yves Soufflet      Post-doctorant, LEGOS, Toulouse, France
\author  Clement MAYET      PhD, LEGOS, Toulouse, France

VERSION :

\date Jan 24, 2012 : Clement MAYET : first documentation

\brief defines useful tools for periodic boundary conditions implementation

***************************************************************************/

#ifndef __PERIODIC_H
#define __PERIODIC_H

#include <list>

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  Associates two sections of open boundary to specify periodic boundary conditions
 
  Periodic boundary conditions are defined between two sections of open boundary. 
  
  This class permits to handle and store informations about the associated edges, 
  vertices and elements.
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

class periodic_t {
private :
public :
  
/*------------------------------------------------------------------------------
  geometric periodicity data */

  list<int> side_1;         /* list of edge numbers of the first  part of open boundary */
  list<int> side_2;         /* list of edge numbers of the second part of open boundary */

  paire_t *vertex_paire;    /* table of paire_t defining the corresponding vertices of side_1 and side_2 */
  paire_t *edge_paire;      /* table of paire_t defining the corresponding edges    of side_1 and side_2 */
  paire_t *element_paire;   /* table of paire_t defining the corresponding elements of side_1 and side_2 */
  
  int  nedge_paire, nvertex_paire, nelement_paire;
  
/*------------------------------------------------------------------------------
  numerical periodicity data*/
  vector<paire_t> LGP0_paire, LGP1_paire, NCP1_paire, LGP2_paire, CQN1_paire;

  bool activated;


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  periodic_t()

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  {
    edge_paire=NULL;
    vertex_paire=NULL;
    element_paire=NULL;
    nedge_paire=nvertex_paire=nelement_paire=0;
    activated=false;
  }
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void destroy()
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  {
    if(edge_paire!=0) {
      delete[] edge_paire;
      edge_paire=0;
      };
    if(vertex_paire!=0) {
      delete[] vertex_paire;
      vertex_paire=0;
      };
    if(element_paire!=0) {
      delete[] element_paire;
      element_paire=0;
      };
    nedge_paire=nvertex_paire=nelement_paire=0;
    LGP0_paire.clear();
    LGP1_paire.clear();
    LGP2_paire.clear();
    NCP1_paire.clear();
    CQN1_paire.clear();
    side_1.clear();
    side_2.clear();
    }

/** Constructor : get the two sides to associate from  mesh_t and code and store into side_1 and side_2
*
* @param mesh
* @param code the code corresponding to periodic boundary conditions (code in bel file )
*/

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  periodic_t(mesh_t & mesh, int code) {
    
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
    edge_paire=NULL;
    vertex_paire=NULL;
    element_paire=NULL;
    nedge_paire=nvertex_paire=nelement_paire=0;
    activated=false;

    int l, n;
    int start;
    int found=0;

/*----------------------------------------------------------------------------
    seek first edge with non-periodic code */
    for(l=0; l<mesh.limits[0].nedges;l++) {
      n=mesh.limits[0].edges[l];
      if(mesh.edges[n].code != code) {
        start=l;
        break;
        }
      }
/*----------------------------------------------------------------------------
    seek first edge with periodic code */
    for(l=start; l<mesh.limits[0].nedges;l++) {
      n=mesh.limits[0].edges[l];
      if(mesh.edges[n].code == code) {
        start=l;
        found=1;
        break;
        }
      }

    if (found == 1) {                       //There is at least one edge with periodic code
/*----------------------------------------------------------------------------
      collect edge with periodic code */
      while(mesh.edges[n].code == code) {
        this->side_1.push_back(n);
        l++;
        l=l%mesh.limits[0].nedges;
        n=mesh.limits[0].edges[l];
        }
/*----------------------------------------------------------------------------
      skip edge with non-periodic code */
      while(mesh.edges[n].code != code) {
        l++;
        l=l%mesh.limits[0].nedges;
        n=mesh.limits[0].edges[l];
        }
/*----------------------------------------------------------------------------
      collect edge with periodic code */
      while(mesh.edges[n].code == code) {
        this->side_2.push_back(n);
        l++;
        l=l%mesh.limits[0].nedges;
        n=mesh.limits[0].edges[l];
        }
      activated=true;
      side_2.reverse();
      }
    }

  /** Associate paires of edges and vertex
   *
   * @param mesh
   * @return int status (0 if ok, -1 if not)
   */
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int create_paires(mesh_t & mesh, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
   int nedges=0, i=0;
/*-----------------------------------------------------------------------------
    test list length : */
    if (side_1.size() == side_2.size()) {
      nedges=side_1.size();
      }
    else{
      printf("Error building periodic boundary conditions, the two sides don't have the same size. Check the .bel file\n");
      exit(-1);
      }

    nedge_paire = nedges;
    if(nedges==0) {
      nvertex_paire=0;
      return(0);
      }
    nvertex_paire = nedges+1;
    edge_paire   = new paire_t[nedge_paire];
    vertex_paire = new paire_t[nvertex_paire];

    std::list<int>::iterator it1(side_1.begin());
    std::list<int>::iterator it1_end(side_1.end());
    std::list<int>::iterator it2(side_2.begin());

    for (;it1 != it1_end; it1++){
      if(*it1 < 0 or *it1 > mesh.nedges) return(-1);
      if(*it2 < 0 or *it2 > mesh.nedges) return(-1);
/*----------------------------------------------------------------------------
      associate edge paires*/
      edge_paire[i].value[0]=*it1;
      edge_paire[i].value[1]=*it2;
/*----------------------------------------------------------------------------
      associate vertex paires (first extremity of side_1 with last of side_2)*/
      vertex_paire[i].value[0]=mesh.edges[*it1].extremity[0];
      vertex_paire[i].value[1]=mesh.edges[*it2].extremity[1];
      if(debug) printf("vertex paire num %d : %d   %d \n", i,vertex_paire[i].value[0],vertex_paire[i].value[1]);

/*----------------------------------------------------------------------------
      this is the last edge, we have to associate the last point of edges (nvertex = nedges+1) */
      if (i == nedges-1) { // if this is the last edge
        vertex_paire[i+1].value[0]=mesh.edges[*it1].extremity[1];
        vertex_paire[i+1].value[1]=mesh.edges[*it2].extremity[0];
        if(debug) printf("vertex paire num %d : %d   %d \n", i+1,vertex_paire[i+1].value[0],vertex_paire[i+1].value[1]);
        }
      it2++;
      i++;
      }
    return(0);
  }
  
  /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  /** This is the general alias function
   *
   * looks in the paire_t table of type "type" (vertex, edge, element) for the num of the element associated
   * to the one having number "target"
   * @param target
   * @param type
   * @return associated element's number
   */
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int alias(int target, int type, int mode)  const {

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
    paire_t *paire;
    int paire_size=0;
/*------------------------------------------------------------------------------------
    mode bi-direcionnal (0), right alias (1), left alias (2) */

    switch (type) {
      case 0:
        paire=vertex_paire;
        paire_size=nvertex_paire;
        break;
      case 1:
        paire=edge_paire;
        paire_size=nedge_paire;
        break;
      case 2:
        paire=element_paire;
        paire_size=nelement_paire;
      }
      
    int n=0;
    switch (mode) {
      case 0:
        for (n=0;n<paire_size;n++) {
          if(paire[n].value[0]==target) {
            return(paire[n].value[1]);
            }
          if(paire[n].value[1]==target) {
            return(paire[n].value[0]);
            }
          }
        break;
      case 1:
        for (n=0;n<paire_size;n++) {
          if(paire[n].value[0]==target) {
            return(paire[n].value[1]);
            }
          }
      case 2:
        for (n=0;n<paire_size;n++) {
          if(paire[n].value[1]==target) {
            return(paire[n].value[0]);
            }
          }
      }
    return(-1);
  }

  /** Get the vertex associated to the target vertex (vertex having number target)
   *
   * @param target
   * @return num of associated vertex
   */
  int vertex_alias(int target, int mode) const {
    return alias(target,0,mode);
  }

  /** Get the edge associated to the target edge (edge having number target)
   *
   * @param target
   * @return num of associated edge
   */
  int edge_alias(int target, int mode) const {
    return alias(target,1,mode);
  }

  /** Get the element associated to the target element (element having number target)
   *
   * @param target
   * @return num of associated element
   */
  int element_alias(int target, int mode) const {
    return alias(target,2,mode);
  }

};


extern int rhs_periodic2D(periodic_t & periodic, double *rhs, int discretisation);
extern int rhs_periodic2D(periodic_t & periodic, complex<double> *rhs, int discretisation);

// extern int rhs_periodic3D(mesh_t & mesh, periodic_t & periodic, complex<double> *rhs, ordering_t *);

// extern int rhs_periodic01(periodic_t *periodic, complex<double> *rhs, ordering_t *, int discretisation);
// extern int rhs_periodic01(periodic_t *periodic, double          *rhs, ordering_t *, int discretisation);

extern int rhs_periodic01(const vector<paire_t> *paires, complex<double> *rhs, ordering_t *ordering, int discretisation);
extern int rhs_periodic01(const vector<paire_t> *paires, double          *rhs, ordering_t *ordering, int discretisation);

extern int rhs_periodic3D(mesh_t & mesh, periodic_t & periodic, complex<double> **rhs, int vdim, int discretisation);

extern int PackedMatrix_periodic_generic(const mesh_t &  mesh, const vector<paire_t> & paires, hypermatrix_t  & M);
extern int PackedMatrix_periodic_generic(const mesh_t &  mesh, const vector<paire_t> & paires, hyperzmatrix_t & M);

extern int PackedMatrix_periodic(const periodic_t & periodic, hypermatrix_t & M,int);
extern int PackedMatrix_periodic(const periodic_t & periodic, hyperzmatrix_t & M,int);

extern int DiagonalMatrix_periodic_LGP1(mesh_t & mesh, periodic_t & periodic, double *M);
extern int force_periodic(const mesh_t & mesh, periodic_t & periodic, double *rhs);

// extern int sum_periodic(mesh_t mesh, periodic_t periodic, double *rhs);

extern int sum_periodic2D(mesh_t & mesh, periodic_t & periodic, double *rhs, int discretisation);
extern int sum_periodic2D(mesh_t & mesh, periodic_t & periodic, complex<double> *rhs, int discretisation);

extern int periodic_ordering(ordering_t *ordering, vector<paire_t> paires, bool debug);

extern int PackedMatrix_periodic02(const mesh_t &  mesh, const vector<paire_t> & paires, hypermatrix_t  & M);
extern int PackedMatrix_periodic02(const mesh_t &  mesh, const vector<paire_t> & paires, hyperzmatrix_t & M);

extern int periodic_neighbours(discretisation_t & descriptor, vector<paire_t> paires, bool debug);
 
extern int mesh_periodic_LGP2(mesh_t  & mesh, periodic_t *periodic, bool debug);
extern int mesh_periodic_LGP1(mesh_t  & mesh, periodic_t *periodic, bool debug);
extern int mesh_periodic_LGP0(mesh_t  & mesh, periodic_t *periodic, bool debug);
extern int mesh_periodic_NCP1(mesh_t  & mesh, periodic_t *periodic, bool debug);
extern int mesh_periodic_CQN1(mesh_t  & mesh, periodic_t *periodic, bool debug);

extern int check_periodic2D(mesh_t & mesh, vector<paire_t> & paires, complex<double> *rhs);
extern int check_periodic2D(mesh_t & mesh, vector<paire_t> & paires, double *rhs);

extern int MomentumSystem_periodic(mesh_t mesh, periodic_t periodic, double **Au, double **Av, double *Bu, double *Bv);

#endif
