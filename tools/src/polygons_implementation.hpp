
/**************************************************************************

  T-UGOm hydrodynamic ocean model, 2006-2017

  Unstructured Ocean Grid initiative

Developers:

  Florent Lyard      LEGOS/CNRS, Toulouse, France
  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
  Cyril Nguyen       LA/CNRS,    Toulouse, France
  Damien Allain      LEGOS/CNRS, Toulouse, France
  
Contributors:
  
  Yoann Le Bars      PhD, LEGOS, Toulouse, France
  Yves Soufflet      Post-doctorant, LEGOS, Toulouse, France
  Clement Mayet      PhD, LEGOS, Toulouse, France
  Frédéric Dupont    Canada environement, Canada

***************************************************************************/

#ifndef POLYGONS_IMPLEMENTATION_HPP
#define POLYGONS_IMPLEMENTATION_HPP

/**
 * \file polygons_implementation.hpp
 * \brief Elementary cycles seeking.
 * \author Yoann LE BARS
 * \version 2.1
 * \date February 7th, 2007
 * \date February 18th, 2007
 * \date June 6th, 2007
 * \date April 15th, 2008
 * \date April 17th, 2008
 * \date May 26th, 2008
 * \date June 17th, 2008
 * \date August 11th, 2008
 *
 * Elementary cycles seeking, which is used to find polygons into drag map.
 * Cf. internal reports for more details.
 *
 * Copyright: see COPYING file.
 */

#include <stack>
#include <limits>
#include <algorithm>
#include <cmath>
#include <functional>
#include <fstream>

#include "polygons.hpp"

namespace Polygons {
// --- Classes and typedefs for cycles gestion. -----------------------------

//    typedef std::deque<size_t>::iterator StackIterator;
  /// Maximum number of vertices a graph can contain.
  const size_t maxVertices = numeric_limits<size_t>::max() - 1;

  /// Type that represents an elementary cycle.
  typedef std::list<Adjacency> ElementaryCycles;
  /// Type that represents the adjacencies of a vertex of a graph.
  typedef std::vector<Adjacency> Adjacencies;

  /**
   * \class Stockage
   * \brief Class for easily remove a point of an adjacency list.
   *
   * Constains the label of an element and a flag that indicates wether the
   * element is present or not.
   */
  class Stockage {
  private:
    /// Label of the element.
    size_t elt;
    /// The presence flag.
    bool pres;

  public:
    // Constructors
    /**
     * \brief Constructs a Stockage with a label and a flag.
     * \param _elt The label.
     * \param _pres The flag.
     */
    Stockage (size_t _elt, bool _pres): elt (_elt), pres (_pres) {}
    /**
     * \brief Copy constructor.
     * \param s The Stockage to be copied.
     */
    Stockage (const Stockage &s): elt (s.elt), pres (s.pres) {}

    /**
     * \brief Controlled access to the label of the element.
     * \return The value.
     */
    size_t value (void) const {return elt;}
    /**
     * \brief Indicates wether the element is present or not.
     * \return True if the element is present, else false.
     */
    bool isPresent (void) const {return pres;}
     /// Set the element as not present.
    void notPresent (void) {pres = false;}
    /// Set the element as present.
    void present (void) {pres = true;}
  };  // class Stockage






/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  /**
   * \brief Equality operator overloading.
   * \param s1 First Stockage to be compared.
   * \param s2 Second Stockage to be compared.
   * \return True if s1 and s2 are equal, else false.
   */
  inline bool operator == (const Stockage &s1, const Stockage &s2) {
    return (s1.value() == s2.value()) && (s1.isPresent() == s2.isPresent());
  }  // bool operator == (const Stockage &, const Stockage &)

  /// List of adjacencies in an easy to remove manner.
  typedef list<Stockage> PresentAdjacency;
  /// Representation of all the adjacencies, in an easy to remove way.
  typedef vector<PresentAdjacency> RemovableAdjacencies;

  // --- Procedures and functions for elementary cycles seeking. --------------

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  /**
   * \fn void nocycle (size_t x, size_t y, RemovableAdjacencies &A, Adjacencies
   * &B)
   * \brief Algorithm 1.
   * \param x Label of the first point.
   * \param y Label of the second point.
   * \param A Adjacencies list of each vertex. Will be modified by the
   * procedure.
   * \param B A set of list, one for each vertex. Will be modified by the
   * procedure
   *
   * See internal report.
   */
  inline void nocycle (size_t x, size_t y, RemovableAdjacencies &A,
                       Adjacencies &B) {
    B.at(y).push_front(x);
    if (x < A.size()) {
      // The vertex to be treated.
      PresentAdjacency::iterator p = find(A[x].begin(), A[x].end(),Stockage (y, true));
      if (p != A[x].end()) p->notPresent();
    }
  }

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  /**
   * \brief Algorithm 2.
   * \param x Label of a vertex in the graph.
   * \param mark For each vertex, indicates whether it is mark or not. Will be
   * modified by the procedure.
   * \param A Adjacencies list of each vertices. Will be modified by the
   * procedure.
   * \param B A set of list, one for each vertex. Will be modified by the
   * procedure.
   *
   * See internal report.
   */
  static void unmark (size_t x, std::vector<bool> &mark,
                      RemovableAdjacencies &A, Adjacencies &B) {
    mark.at(x) = false;
    if (x < B.size()) {
      // p is the vertex being treated.
      for (Adjacency::iterator p = B[x].begin(); p != B[x].end(); p++) {
        // Label of the current vertex.
        const size_t y = *p;
        if (y < A.size()) {
          // The vertex labelled x if it is adjacent to p.
          PresentAdjacency::iterator q = find(A[y].begin(), A[y].end(),
                                              Stockage (x, false));
          if (q != A[y].end()) q->present();
          if (mark.at(y)) unmark(y, mark, A, B);
        }
      }
      B[x].clear();
    }
  }

  /**
   * \brief Algorithm 3.
   * \param v Label of the node.
   * \param q1 Position of v in the stack internal to this function.
   * \param mark Indicates if a node is marked or not. Will be modified by the
   * function.
   * \param N Number of nodes in the graph.
   * \param reach Indicates if a node had been treated or not. Will be modified
   * by the function.
   * \param A Adjacencies list of each vertices. Will be modified by the
   * function.
   * \param B A set of list, one for each vertices. Will be modified by the
   * function.
   * \param elementaryCycles Contains the elementary cycles that has been
   * found. Will be modified by the function.
   * \return True if an elementary cycle has been found, else false.
   *
   * See internal report.
   */
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  static bool seek_cycle (size_t v, size_t q1, std::vector<bool> &mark, size_t N,
		     std::vector<bool> &reach, RemovableAdjacencies &A,
		     Adjacencies &B, ElementaryCycles &elementaryCycles, std::vector<size_t> & position, std::deque<size_t> & s) {
    
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
    // Simplifies the gestion of the stack.
//    typedef std::deque<size_t>::iterator StackIterator;
    
/* Indicates the position of every vertex in the stack. Conserve its value
   from one call of the function to another. */
//    static std::vector<size_t> position (N);
    
/* The internal stack. Stay the same from one call of the function to
    another. */
//    static std::deque<size_t> s;
    s.push_front(v);
    // Number of nodes in the stack.
    const size_t t = s.size();
    
    // First position of v in the stack.
    const size_t q = reach.at(v)? q1: t;
    // Return value of the function.
    bool f = false;
    mark.at(v) = true;
    position.at(v) = t;

    // p is a reference on the node being treated.
    for (PresentAdjacency::const_iterator p = A.at(v).begin(); p != A.at(v).end(); p++) {
      if (p->isPresent()) {
      // Current node.
        const size_t w = p->value();
        if (!mark.at(w)) {
//          printf("try cycle from %d with highest degree node %d\n",v+1,w+1);
	  if (seek_cycle(w, q, mark, N, reach, A, B, elementaryCycles,position,s)) {
           f = true;
            }
	  else {
            nocycle(v, w, A, B);
            }
          }
        else if (position.at(w) <= q) {
	// Stock the elementary cycle that has been found.
	// Reference to v.
	  const deque<size_t>::const_reverse_iterator i = find(s.rbegin(), s.rend(), v);
	// Reference to w.
	  const deque<size_t>::const_reverse_iterator j = find(s.rbegin(), s.rend(), w);
	// The first one in the stack between v and w.
	  const deque<size_t>::const_reverse_iterator mini = min(i, j);
	// The last one.
	  const deque<size_t>::const_reverse_iterator maxi = max(i, j);
	// The cycle itself.
	  Adjacency theCycle;
	  for (deque<size_t>::const_reverse_iterator it = mini; it != maxi; it++)
	    theCycle.push_back(*it);
	  theCycle.push_back(*maxi);

/*------------------------------------------------------------------------------
         Determine if the elementary cycle has to be retain.*/
	// Size of the cycle.
	  const size_t cycleSize = theCycle.size();
//          printf("found cycle from %d with %d nodes\n",v+1,cycleSize);
	  if (cycleSize > 2) {
	  // Indicates whether the cycle is acceptable or not.
	    bool acceptable = true;
//            printf("found cycle %d from %d with %d nodes\n",elementaryCycles.size(),v+1,cycleSize);
	  // i is a reference to a cycle.
	    for (ElementaryCycles::const_iterator i = elementaryCycles.begin();
	       (i != elementaryCycles.end()) && acceptable; i++) {
	    // The tested cycle.
	      const Adjacency testedCycle = *i;
	    // Indicates whether the cycle has been found or not.
	      bool temporary = true;
	      if (testedCycle.size() == cycleSize) {
	        temporary = false;
	      // j is a reference to a node of the tested cycle.
	        for (Adjacency::const_iterator j = testedCycle.begin(); (!temporary) && (j != testedCycle.end()); j++) {
		// The node.
		  const size_t current = *j;
		// Reference to the node if it had been found.
		  const Adjacency::const_iterator k = find(theCycle.begin(), theCycle.end(), current);
		  if (k == theCycle.end()) temporary = true;
	          }
	        }
	      if (!temporary) acceptable = false;
	      }
	    if (acceptable) {
              elementaryCycles.push_back(theCycle);
//              printf("found cycle %d to be eligible\n",elementaryCycles.size());
              }
	    }

	// Back to the third algorithm.
	  f = true;
          } 
        else nocycle(v, w, A, B);
        }
      }
    s.pop_front();
    if (f) unmark(v, mark, A, B);
    reach.at(v) = true;
    position.at(v) = N + 1;

    return f;
  }
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  static bool seek_cycle (size_t v, size_t q1, std::vector<bool> &mark, size_t N,
                     std::vector<bool> &reach, RemovableAdjacencies &A,
                     Adjacencies &B, ElementaryCycles &elementaryCycles, int maxsize)
    
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
 {
   // Simplifies the gestion of the stack.
//     typedef std::deque<size_t>::iterator StackIterator;
    
/* *----------------------------------------------------------------------------
    Indicates the position of every vertex in the stack. Conserve its value
    from one call of the function to another. */
    static std::vector<size_t> position (N);
    
/* *----------------------------------------------------------------------------
    The internal stack. Stay the same from one call of the function to another. */
    static std::deque<size_t> stack;
    size_t test;
    bool debug=false;

/* *----------------------------------------------------------------------------
    add v to stack*/
    stack.push_front(v);

    // Number of nodes in the stack.
    const size_t t = stack.size();
    // First position of v in the stack.
    const size_t q = reach.at(v)? q1: t;
    // Return value of the function.
    bool f = false;

    mark.at(v) = true;
    position.at(v) = t;
//   if (t>maxsize) return(f);

    // p is a reference on the node being treated.
    for (PresentAdjacency::const_iterator p = A.at(v).begin(); p != A.at(v).end(); p++) {
      if (p->isPresent()) {
      // Current node.
        const size_t w = p->value();
         if (!mark.at(w) ) {
/* *----------------------------------------------------------------------------
           w not yet passed, continue the process*/
//         printf("try cycle from %d with highest degree node %d\n",v+1,w+1);
//         t=position.at(w) - position.at(v);
           if (seek_cycle(w, q, mark, N, reach, A, B, elementaryCycles,maxsize)) {
           f = true;
//            printf("successful\n",v+1,w+1);
            }
          else {
            nocycle(v, w, A, B);
            }
         test=position.at(w) - position.at(v);
          }
       else if (position.at(w) <= q) {
/* *----------------------------------------------------------------------------
           w already passed*/
        // Stock the elementary cycle that has been found.
        // Reference to v.
          const deque<size_t>::const_reverse_iterator i = find(stack.rbegin(), stack.rend(), v);
        // Reference to w.
          const deque<size_t>::const_reverse_iterator j = find(stack.rbegin(), stack.rend(), w);
        // The first one in the stack between v and w.
          const deque<size_t>::const_reverse_iterator mini = min(i, j);
        // The last one.
          const deque<size_t>::const_reverse_iterator maxi = max(i, j);
        // The cycle itself.
          Adjacency theCycle;
          for (deque<size_t>::const_reverse_iterator it = mini; it != maxi; it++)
            theCycle.push_back(*it);
          theCycle.push_back(*maxi);

/* *----------------------------------------------------------------------------
         Determine if the elementary cycle has to be retain.*/
        // Size of the cycle.
          size_t cycleSize = theCycle.size();
//          printf("found cycle from %d with %d nodes\n",v+1,cycleSize);
/* *----------------------------------------------------------------------------
          Eliminate over-sized cycles if prescribed (-1 means no limit)*/
//          printf("found cycle %d to be eligible, size=%d\n",elementaryCycles.size(),cycleSize);
          if(maxsize!=-1) {
             if (cycleSize > maxsize) cycleSize=0;
             }
          if (cycleSize > 2) {
          // Indicates whether the cycle is acceptable or not.
            bool acceptable = true;
//            printf("found cycle %d from %d with %d nodes\n",elementaryCycles.size(),v+1,cycleSize);
          // i is a reference to a cycle.
            for (ElementaryCycles::const_iterator i = elementaryCycles.begin();
               (i != elementaryCycles.end()) && acceptable; i++) {
            // The tested cycle.
              const Adjacency testedCycle = *i;
            // Indicates whether the cycle has been found or not.
              bool temporary = true;
              if (testedCycle.size() == cycleSize) {
                temporary = false;
              // j is a reference to a node of the tested cycle.
                for (Adjacency::const_iterator j = testedCycle.begin(); (!temporary) && (j != testedCycle.end()); j++) {
                // The node.
                  const size_t current = *j;
                // Reference to the node if it had been found.
                  const Adjacency::const_iterator k = find(theCycle.begin(), theCycle.end(), current);
                  if (k == theCycle.end()) temporary = true;
                  }
                }
              if (!temporary) acceptable = false;
              }
            if (acceptable) {
              elementaryCycles.push_back(theCycle);
              if(debug) printf("found cycle %d to be acceptable\n",elementaryCycles.size());
              }
            }

        // Back to the third algorithm.
          f = true;
          }
        else nocycle(v, w, A, B);
        }
      }
    stack.pop_front();
    if (f) unmark(v, mark, A, B);
    reach.at(v) = true;
    position.at(v) = N + 1;

    return f;
  }

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  /**
   * \brief Algorithm 5.
   * \param v Node that indicates the convex component.
   * \param treated Indicates whether the node has been treated or not. Will
   * be modified by the function.
   * \param D The graph.
   * \return The node with the highest degree in the convex component.
   *
   * See internal report.
   */
  template <class T> static size_t seek (size_t v, std::vector<bool> &treated, const GraphBase<T> &D) {
    assert(v < D.size());

    // List of adjacencies of the node v.
    const Adjacency adjacency = D[v].adjacency();
    // Inspected vertice with the highest degree.
    size_t vertice = v;
    // Max degree for a node in the convex component.
    size_t maxDegree = adjacency.size();
    treated.at(v) = true;

    // p is a reference to a node.
    for (Adjacency::const_iterator p = adjacency.begin(); p != adjacency.end(); p++) {
      // The current node.
      const size_t w = *p;
      if (!treated.at(w)) {
        // Adjacent node to w with the highest degree.
        const size_t n = seek(w, treated, D);
        // Local maximal degree.
        const size_t localMax = D.at(n).adjacency().size();
        if (localMax > maxDegree) {
          maxDegree = localMax;
          vertice = n;
        }
      }
    }

    return vertice;
  }

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  /**
   * \fn size_t successor (size_t v, const Adjacency &A)
   * \brief Function which returns the successor of a vertice of a cycle.
   * \param v The vertex.
   * \param A Adjacencies list of each vertices.
   * \return The successor of v in the cycle.
   */
  inline size_t successor (size_t v, const Adjacency &A) {
    // Reference to a vertex.
    Adjacency::const_iterator p = find(A.begin(), A.end(), v);
    if (++p == A.end()) p = A.begin();
    return *p;
  }

  /**
   * \fn PolygonSet<Coordinate, Z0> base (const ElementaryCycles &set, const
   * Graph<Coordinate> &G)
   * \brief Algorithm 9.
   * \param set Set of all polygons.
   * \param G The graph
   * \return Minimal fundamental base of the set of polygons.
   *
   * See internal report.
   */
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <class Coordinate, class Z0>
  static PolygonSetBase<Coordinate, Z0> base (const ElementaryCycles &set,
					      const GraphBase<Coordinate> &G) {
    
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
    // Make the set of all find polygons
    // Contains the base of the polygons set.
    PolygonSetBase<Coordinate, Z0> base_;
    bool debug=false;
    
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
    store all cycles as elligible polygons
 
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
    
    // i is a reference to a cycle.
    for (ElementaryCycles::const_iterator i = set.begin(); i != set.end(); i++) {
      // The current cycle.
      const Adjacency cycle = *i;
      // The polygon associate to the cycle.
      PolygonBase<Coordinate, Z0> polygon;
      // place is a reference to a node of the cycle.
      size_t count=0;
      for (Adjacency::const_iterator place = cycle.begin(); place != cycle.end(); place++) {
        polygon.push_back(G[*place].coordinates());
        count++;
        }
//      printf("cycle has %d nodes\n",count);
      base_.push_back(polygon);
      }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
    eliminate polygons made of 2 smaller ones concatenation
 
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
    
redo:
/*------------------------------------------------------------------------------
    Sort the base in decreasing area order*/
    base_.sort(greater< PolygonBase<Coordinate, Z0> >());
/*------------------------------------------------------------------------------
    Dump full set of cycles*/
//    base_.dump();

    for (typename PolygonSetBase<Coordinate, Z0>::iterator i = base_.begin(); i != base_.end();i++) {
      const PolygonBase<Coordinate, Z0> p = *i;
      double global=p.area2();
      if(debug) printf("testing area %lf\n",global);
      for (typename PolygonSetBase<Coordinate, Z0>::const_iterator j = i; j != base_.end(); j++) {
        const PolygonBase<Coordinate, Z0> q = *j;
        for (typename PolygonSetBase<Coordinate, Z0>::const_iterator k = j; k != base_.end(); k++) {
          const PolygonBase<Coordinate, Z0> r = *k;
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
          accuracy issue: fixed limit (1.e-6) will fail in tiny polygon cases
 
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
//           if(fabs(q.area2()+r.area2()-global) < 1.e-6) {
          if(fabs(q.area2()+r.area2()-global) < 1.e-6*global) {
/*------------------------------------------------------------------------------
            polygon is the concatenation of j and k, thus not elementary */
             if(debug) printf("we have a paire in %d cycles; areas %lf %lf %lf  %lf\n",base_.size(),q.area2(),r.area2(),q.area2()+r.area2(),global);
             base_.erase(i);
             goto redo;
             }
           }
        }
     }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
    Construct the minimum basis
 
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
    // i is a reference to a polygon.
    for (typename PolygonSetBase<Coordinate, Z0>::iterator i = base_.begin(); i != base_.end();) {
     /** i is not incremented here because I need to save it before */
     /* Indicates whether or not the current polygon must be remove from the base. */
      bool remove = false;
      // j is a reference to a vertex of the current polygon.
      for (typename PolygonSetBase<Coordinate, Z0>::const_iterator j = base_.begin(); (j != base_.end()) && !remove; j++) {
        if (j != i) {
        // The polygon to be tested.
          const PolygonBase<Coordinate, Z0> testPolygon = *i;
        // Current polygon.
          const PolygonBase<Coordinate, Z0> polygon = *j;
        // Indicates wether the polygon to be tested is suspect or not.
          bool suspect = true;
        // Number of vertex in the polygon.
          const size_t n = polygon.size();
        // k is the label of a vertex of the polygon.
/* *----------------------------------------------------------------------------
          look if one point of j at least is out of i*/
          for (size_t k = 0; suspect && (k < n); k++)
            if (testPolygon.where(polygon[k]) == out) suspect = false;
/* *----------------------------------------------------------------------------
          if not so, j fully included in i...*/
          if (suspect) {
            suspect = false;
/* *@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Development notes:

  Check : MANDATORY !!!

  Note:

  25/10/2008:

    the following test is not correct !!!


    for (size_t k = 0; !suspect && (k < n) ; k++) {
      if (testPolygon.where(polygon[k]) == vertex) suspect = true;
      }

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */

//             size_t count=0;
//
// 	    for (size_t k = 0; k < n; k++) {
// /*            k is the label of a vertex of the polygon */
// 	      if (testPolygon.where(polygon[k]) == vertex) count++;
//               }
//             if(count!=n) {
//               suspect = true;
//               }
// 	    if (suspect) {
//               remove = true;
//               printf("remove cycle %d\n");
//               }
            }
          }
        }
      // saving i into p
      typename PolygonSetBase<Coordinate, Z0>::iterator p = i;
      /* iterating i before removing, otherwise no one can knows what will
         happend */
      i++;
      // remove p if it has to be remove
      if (remove) base_.erase(p);
    }

/*------------------------------------------------------------------------------
    Sort the base in decreasing area order*/
    base_.sort(greater< PolygonBase<Coordinate, Z0> >());

    return base_;
  }

  /**
   * \brief Algorithm 7.
   * \param D The graph.
   * \return The minimum fundamental basis of the graph.
   *
   * See internal report.
   */
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <class Coordinate, class Z0>
   PolygonSetBase<Coordinate, Z0> searchPolygons (const GraphBase<Coordinate> &D) {
    
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
    // Number of nodes in the graph.
    bool debug=false;
    const size_t N = D.size();
    if (N >= maxVertices) throw CapacityOverflow ();

    // Adjacencies list of each vertex.
    RemovableAdjacencies A (N);
    // Set of list, one for each vertex.
    Adjacencies B (N);
    // Indicates if a node is marked or not.
    vector<bool> mark (N);
    // Indicates if a node had been reached or not.
    vector<bool> reach (N);
    // Indicates if a node had been treated or not.
    vector<bool> treated (N);
    // List of elementary cycles.
    ElementaryCycles elementaryCycles;
    std::vector<size_t> position (N);
    std::deque<size_t> s;

    // Initialization
    // j is the label of a node of the graph.
    for (size_t j = 0; j < N; j++) {
      // List of adjacency of the current node.
      const Adjacency adjacency = D[j].adjacency();
      // p is the reference to a node.
      for (Adjacency::const_iterator p = adjacency.begin(); p!= adjacency.end(); p++) {
        A[j].push_back(Stockage (*p, true));
        }
      mark[j] = false;
      reach[j] = false;
      treated[j] = false;
      }

    // seeking
    // v is the label of a vertex of the graph.
    for (size_t v = 0; v < N; v++) { 
      if (!treated[v]) {
        size_t m;
        m=seek(v, treated, D);
        if(debug) printf("try cycle from %d with highest degree node %d\n",v+1,m+1);
        seek_cycle(m, 0, mark, N, reach, A, B, elementaryCycles,position,s);
        }
      }
    if(debug) printf("found %d cycles\n", elementaryCycles.size());
    return base<Coordinate, Z0>(elementaryCycles, D);
  } 

  /**
   * \brief Algorithm 7.
   * \param D The graph.
   * \return The minimum fundamental basis of the graph.
   *
   * See internal report.
   */
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <class Coordinate, class Z0>
    static PolygonSetBase<Coordinate, Z0> searchPolygons (const GraphBase<Coordinate> &D,int maxsize) {
    
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
    // Number of nodes in the graph.
    const size_t N = D.size();
    bool debug=false;

    if (N >= maxVertices) throw CapacityOverflow ();

    // Adjacencies list of each vertex.
    RemovableAdjacencies A (N);
    // Set of list, one for each vertex.
    Adjacencies B (N);
    // Indicates if a node is marked or not.
    vector<bool> mark (N);
    // Indicates if a node had been reached or not.
    vector<bool> reach (N);
    // Indicates if a node had been treated or not.
    vector<bool> treated (N);
    // List of elementary cycles.
    ElementaryCycles elementaryCycles;

    // Initialization
    // j is the label of a node of the graph.
    for (size_t j = 0; j < N; j++) {
      // List of adjacency of the current node.
      const Adjacency adjacency = D[j].adjacency();
      // p is the reference to a node.
      for (Adjacency::const_iterator p = adjacency.begin(); p!= adjacency.end(); p++) {
        A[j].push_back(Stockage (*p, true));
        }
      mark[j] = false;
      reach[j] = false;
      treated[j] = false;
      }


    // seeking
    // v is the label of a vertex of the graph.
    for (size_t v = 0; v < N; v++) {
      if (!treated[v]) {
        size_t m;
        m=seek(v, treated, D);
        if(debug) printf("try cycle from %d with highest degree node %d\n",v+1,m+1);
        seek_cycle(m, 0, mark, N, reach, A, B, elementaryCycles,maxsize);
//        cycle(seek(v, treated, D), 0, mark, N, reach, A, B, elementaryCycles);
        }
      }
    if(debug) printf("found %d cycles\n", elementaryCycles.size());
    return base<Coordinate, Z0>(elementaryCycles, D);
  }

  /**
   * \brief Algorithm 7-bis.
   * \param D The graph.
   * \return The minimum fundamental basis of the graph.
   *
   * See internal report.
   */
  template <class Coordinate>
  static ElementaryCycles searchCycles (const GraphBase<Coordinate> &D,int maxsize) {
    // Number of nodes in the graph.
    const size_t N = D.size();
    bool debug=false;
    
    if (N >= maxVertices) throw CapacityOverflow ();

    // Adjacencies list of each vertex.
    RemovableAdjacencies A (N);
    // Set of list, one for each vertex.
    Adjacencies B (N);
    // Indicates if a node is marked or not.
    vector<bool> mark (N);
    // Indicates if a node had been reached or not.
    vector<bool> reach (N);
    // Indicates if a node had been treated or not.
    vector<bool> treated (N);
    // List of elementary cycles.
    ElementaryCycles elementaryCycles;

    // Initialization
    // j is the label of a node of the graph.
    for (size_t j = 0; j < N; j++) {
      // List of adjacency of the current node.
      const Adjacency adjacency = D[j].adjacency();
      // p is the reference to a node.
      for (Adjacency::const_iterator p = adjacency.begin(); p!= adjacency.end(); p++) {
        A[j].push_back(Stockage (*p, true));
        }
      mark[j] = false;
      reach[j] = false;
      treated[j] = false;
      }

    // seeking
    // v is the label of a vertex of the graph.
    for (size_t v = 0; v < N; v++) {
      if (!treated[v]) {
        size_t m;
        m=seek(v, treated, D);
        if(debug) printf("try cycle from %d with highest degree node %d\n",v+1,m+1);
        seek_cycle(m, 0, mark, N, reach, A, B, elementaryCycles,maxsize);
//        cycle(seek(v, treated, D), 0, mark, N, reach, A, B, elementaryCycles);
        }
      }
    if(debug) printf("found %d cycles\n", elementaryCycles.size());
    return (elementaryCycles);
  }

  /**
   * \brief Algorithm 10: computes the area of a triangle times 2.
   * \param a First vertex of the triangle.
   * \param b Second vertex of the triangle.
   * \param c Third vertex of the triangle.
   * \return The area.
   *
   * See [O'Rourke, 1998] and internal report.
   */
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <class T> inline T triangleArea2 (const PointBase<T> &a,
                                           const PointBase<T> &b,
                                           const PointBase<T> &c) {
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  return (b.x() - a.x()) * (c.y() - a.y())  - (c.x() - a.x()) * (b.y() - a.y());
  }

  /**
   * \brief Read the reference points in a file.
   * \param fileName Name of the file containing the reference points.
   * \return A vector containing all the reference points.
   */

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <class Coordinate, class Z0>
  static vector< ReferencePoint<Coordinate, Z0> > readReferencePoints (const std::string &fileName)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/* *----------------------------------------------------------------------------
  Reads data associated with polygones */
  // Input flow on the file.
  ifstream file (fileName.c_str());
  if (!file) throw ReadError (fileName);

  // Vector containing the reference points.
  vector< ReferencePoint<Coordinate, Z0> > v;
  
  // Abscissa, ordinate of the point.
  Coordinate x,y;
  // Value associated to the point.
  Z0 z0;

  file >> x >> y >> z0; /// HERE !!!
  while (!file.eof()) {
    if (file.bad()) throw ReadError (fileName);
    v.push_back(ReferencePoint<Coordinate, Z0> (x, y, z0));
    file >> x >> y >> z0;
    }

  return v;
  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <class T> GraphBase<T>::GraphBase (const std::string &fileName)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
 {
  mesh_t mesh;
  int status;
  
/*------------------------------------------------------------------------------
  Reads a graph in a trigrid format file */
  status=fe_readmesh(fileName.c_str(), MESH_FILE_FORMAT_TRIGRID, &mesh);
  
  if (status != 0)
    throw ReadError (fileName);

/* *----------------------------------------------------------------------------
    i is the label of a node of the graph. */
  for (size_t i = 0; i < static_cast<size_t>(mesh.nvtxs); i++) {
/* *----------------------------------------------------------------------------
      Current point. */
      const PointBase<T> p (static_cast<T>(mesh.vertices[i].lon),
                            static_cast<T>(mesh.vertices[i].lat));
/* *----------------------------------------------------------------------------
      List of adjacency of the point. */
    Adjacency adj;
/* *----------------------------------------------------------------------------
      j is the label of an adjacency point. */
      for (size_t j = 0; j < static_cast<size_t>(mesh.vertices[i].nngh); j++)
        adj.push_back(static_cast<size_t>(mesh.vertices[i].ngh[j]));
      this->push_back(ElementBase<T> (p, adj));
      }
    }

template <class T> GraphBase<T>::GraphBase (const mesh_t & mesh) {
/* *----------------------------------------------------------------------------
    Reads a graph in a file. The graph stocked in trigrid format. */
    int status;
    
/* *----------------------------------------------------------------------------
    i is the label of a node of the graph. */
  for (size_t i = 0; i < static_cast<size_t>(mesh.nvtxs); i++) {
/* *----------------------------------------------------------------------------
      Current point. */
      const PointBase<T> p (static_cast<T>(mesh.vertices[i].lon),
                            static_cast<T>(mesh.vertices[i].lat));
/* *----------------------------------------------------------------------------
      List of adjacency of the point. */
      Adjacency adj;
/* *----------------------------------------------------------------------------
      j is the label of an adjacency point. */
      for (size_t j = 0; j < static_cast<size_t>(mesh.vertices[i].nngh); j++)
        adj.push_back(static_cast<size_t>(mesh.vertices[i].ngh[j]));
      this->push_back(ElementBase<T> (p, adj));
      }
    }

  /**
   * \fn void PolygonBase<Coordinate, Z0>::computeArea2 (void) const
   * \brief Algorithm 11: computes the area of the polygon times 2 and stock it
   * in a2.
   *
   * See [O'Rourke, 1998].
   *
   * In [O'Rourke, 1998], the polygon has to be described anticlockwise, but if
   * the absolute value of the sum is transmitted, this is useless.
   */
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <class Coordinate, class Z0>
  void PolygonBase<Coordinate, Z0>::computeArea2 (void) const {

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  areaComputed = true;

  // Starting point of the polygon.
  const PointBase<Coordinate> p = at(0);
  // Sum of the area of all the triangles.
  Coordinate sum = 0;
  // Number of vertices in the polygon.
  const size_t n = this->size();
  // i is the label of the current vertex.
  for (size_t i = 1; i < n; i++) {
    // Current vertex.
    const PointBase<Coordinate> a = operator[](i);
    // Next vertex.
    const PointBase<Coordinate> next = operator[]((i + 1) % n);
    sum += triangleArea2(p, a, next);
    }

  areax2 = abs(sum);
  }

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  /**
   * \fn Position PolygonBase<Coordinate, Z0>::where
   * (const PointBase<Coordinate> &q) const
   * \brief Algorithm 8: gives the position of a point relatively to the
   * polygon.
   * \param q The point to be tested.
   * \return The position as a Position (type).
   *
   * See [O'Rourke, 1998].
   */
  template <class Coordinate, class Z0>
  Position PolygonBase<Coordinate, Z0>::where (const PointBase<Coordinate> &q)
    const {
    // Number of vertices in the polygon.
    const size_t n = this->size();
    assert (n != 0);

/*------------------------------------------------------------------------------
    How many times the edges of the polygon are crossed on the right.*/
    size_t rCross = 0;
/*------------------------------------------------------------------------------
    How many times the edges of the polygon are crossed on the left. */
    size_t lCross = 0;
    // i is the current vertex.
    for (size_t i = 0; i < n; i++) {
      // Current polygon vertex.
      const PointBase<Coordinate> P1 = vector< PointBase<Coordinate> >::operator[](i) - q;
      if (P1 == PointBase<Coordinate> (0, 0)) return vertex;
      // Label of the next vertex.
      const size_t i1 = (i + n - 1) % n;
      // Next vertex.
      const PointBase<Coordinate> P2 = vector< PointBase<Coordinate> >::operator[](i1) - q;
      // Are the edge of the polygon crossed on the right ?
      const bool rStrad = (P1.y() > 0) != (P2.y() > 0);
      // Are the edge of the polygon crossed on the left ?
      const bool lStrad = (P1.y() < 0) != (P2.y() < 0);
      if (rStrad || lStrad) {
        // Is the crossing really done?
        const Coordinate x = (P1.x() * P2.y() - P2.x() * P1.y()) / (P2.y() - P1.y());
        if (rStrad && (x > 0)) rCross++;
        if (lStrad && (x < 0)) lCross++;
        }
      }
    // Precise where is the cross done.
    const size_t dummy = rCross % 2;
    if (dummy != (lCross % 2)) return edge;
    if (dummy == 1) return in;
    else return out;
  }

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <class Coordinate, class Z0> int PolygonSetBase<Coordinate, Z0>::setz0
    (const vector< ReferencePoint<Coordinate, Z0> > &v)
    
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  /**
   * \brief Algorithm 12.
   * \param v Set of ReferencePoint.
   *
   * See internal report.
   */
{     
    int status=0;
/*------------------------------------------------------------------------------
    initialisation; q is a reference to a polygon */
    for (iterator q = this->begin(); q != this->end(); q++)
      q->setz0(defaultz0_);
    
/*------------------------------------------------------------------------------
    i is the label of the current ReferencePoint */
    for (size_t i = 0; i < v.size(); i++) {
      // Reference to a polygon.
      iterator p;
      int count=0;
      int plg;
      ReferencePoint<Coordinate, Z0> P=v[i];
      const Coordinate x0 = P.x();
      this->recale(x0);
/*------------------------------------------------------------------------------
      q is a reference to the polygon being tested */
      int index=0;
      for (iterator q = this->begin(); q != this->end(); q++) {
/*------------------------------------------------------------------------------
        Position of the current ReferencePoint relatively to the polygon being tested */
        const Position pos = q->where(v[i]);
        if ((pos == vertex) || (pos == edge))
          throw OnBoundaryPoint (static_cast<double>(v[i].x()), static_cast<double>(v[i].y()), i);
        if (pos == in) {
          p = q;
          plg=index;
          count++;
          }
        index++;
        }
      switch (count) {
        case 0:
          printf("no polygons found for %lf %lf\n",v[i].x(),v[i].y());
          status=-1;
          break;
        case 1:
          printf("polygon prescription %2d : lon=%lf lat=%lf value=%lf plg=%d\n",i+1,v[i].x(),v[i].y(),v[i].z0(),plg+1);
          p->setz0(v[i].z0());
          break;
        default:
          printf("Warning: too many polygons (%d) found for %lf %lf\n",count,v[i].x(),v[i].y());
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    assume last found polygon is the right one (as polygons are given in
    descending area order)

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
          printf("polygon prescription %2d : lon=%lf lat=%lf value=%lf plg=%d\n",i+1,v[i].x(),v[i].y(),v[i].z0(),plg+1);
          p->setz0(v[i].z0());
          break;
      }
    }
  return(status);
  }

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <class Coordinate, class Z0> void PolygonBase<Coordinate, Z0>::dumpCoord(FILE *out, size_t s) 

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
    size_t p;

    const size_t n = this->size();
    p=0;
    for (size_t i = 0; i < n; i++) {
      p++;
      const PointBase<Coordinate> P1 = vector< PointBase <Coordinate> >::operator[](i);
//      const PointBase<Coordinate> P1 = PointBase<Coordinate> (0, 0);
      double x=P1.x();
      double y=P1.y();
      fprintf(out,"%d %lf %lf\n",p,x,y);
      }
    for (size_t i = 0; i < 1; i++) {
      p++;
      const PointBase<Coordinate> P1 = vector< PointBase <Coordinate> >::operator[](i);
//      const PointBase<Coordinate> P1 = PointBase<Coordinate> (0, 0);
      double x=P1.x();
      double y=P1.y();
      fprintf(out,"%d %lf %lf\n",p,x,y);
      }
  }

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <class Coordinate, class Z0> void PolygonSetBase<Coordinate, Z0>::dump(const char *filename)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
    // q is a reference to a polygon.
  FILE *out;
  
  if(filename!="") out=fopen(filename,"w");
  else out=stdout;
                                          
    size_t s=0;
    for (iterator q = this->begin(); q != this->end(); q++) {
      s++;
      fprintf(out,"%d %d %lf\n",s,q->size()+1,q->z0());
      q->dumpCoord(out, s);
      };
  
  if(out!=stdout) fclose(out);

  }

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <class Coordinate, class Z0> void PolygonBase<Coordinate, Z0>::recale(Coordinate x0, size_t s) 

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
    size_t p;

    const size_t n = this->size();
    if(isnan(x0)) {
      p=0;
      PointBase<Coordinate> P0 = vector< PointBase <Coordinate> >::operator[](0);
      x0 = P0.x();
      }
    for (size_t i = 0; i < n; i++) {
      p++;
      PointBase<Coordinate> & P1 = vector< PointBase <Coordinate> >::operator[](i);
//       double x=P1.x();
//       double y=P1.y();
//       printf("%d %lf %lf\n",p,x,y);
//       double test=this->x(i);
      Coordinate & xx = P1.x();
      if(xx > x0 +180.0) xx-=360.0;
      if(xx < x0 -180.0) xx+=360.0;
      x0=xx;
      }
  }

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <class Coordinate, class Z0> void PolygonSetBase<Coordinate, Z0>::recale(Coordinate x0)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
    // q is a reference to a polygon.
    size_t s=0;
    for (iterator q = this->begin(); q != this->end(); q++) {
      s++;
//       printf("%d %d\n",s,q->size()+1);
      q->recale(x0,s);
      };

  }

  /**
   * \brief Determine the value of z0 for every node of the mesh.
   * \param mesh The mesh.
   * \return A vector containing the value of z0 associated to each point of
   * the mesh.
   */

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <class Coordinate, class Z0>
  std::vector<Z0> PolygonSetBase<Coordinate, Z0>::determineValue (const mesh_t &mesh) const

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/*------------------------------------------------------------------------------
    Number of vertices in the mesh. */
    const size_t N = mesh.nvtxs;
    
/*------------------------------------------------------------------------------
    Contains the value of z0 associated to each point of the mesh. */
    std::vector<Z0> vz0 (N);

    for (size_t i = 0; i < N; i++) {
      vz0[i] = defaultz0_;
      double t=mesh.vertices[i].lon;
      double p=mesh.vertices[i].lat;
      for (const_iterator polygon = this->begin(); polygon != this->end(); polygon++) {
        const PointBase<Coordinate> p0 = polygon->at(0);
        if(t<p0.x()-180.) t+=360.;
        if(t>p0.x()+180.) t-=360.;

        if (polygon->where(PointBase<Coordinate> (static_cast<Coordinate>(t), 
                                                  static_cast<Coordinate>(p))) != out) {
          vz0[i] = polygon->z0();
          }
        }
      }

    return vz0;
  }

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <class Coordinate, class Z0>
  std::vector<Z0> PolygonSetBase<Coordinate, Z0>::determineValue (const grid_t & grid) const

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/*------------------------------------------------------------------------------
    Number of vertices in the mesh. */
    const size_t N = grid.Hsize();
/*------------------------------------------------------------------------------
    Contains the value of z0 associated to each point of the mesh. */
    std::vector<Z0> vz0 (N);

    for (size_t j = 0; j < grid.ny; j++) {
      for (size_t i = 0; i < grid.nx; i++) {
	int m=grid.nx*j+i;
        vz0[m] = defaultz0_;
        double t,p;
        grid.xy(i,j,t,p);
        for (const_iterator polygon = this->begin(); polygon != this->end(); polygon++) {
          const PointBase<Coordinate> p0 = polygon->at(0);
          if(t<p0.x()-180.) t+=360.;
          if(t>p0.x()+180.) t-=360.;
          if (polygon->where(PointBase<Coordinate> (static_cast<Coordinate>(t), 
                                                    static_cast<Coordinate>(p))) != out) {
            vz0[m] = polygon->z0();
            }
          }
        }
      }

    return vz0;
  }

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <class Coordinate, class Z0>
  std::vector<Z0> PolygonSetBase<Coordinate, Z0>::determineFlag (const grid_t &grid) const

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/*------------------------------------------------------------------------------
    Number of vertices in the mesh. */
    const size_t N = grid.Hsize();
/*------------------------------------------------------------------------------
    Contains the value of z0 associated to each point of the mesh. */
    std::vector<Z0> vz0 (N);

    for (size_t j = 0; j < grid.ny; j++) {
      for (size_t i = 0; i < grid.nx; i++) {
	int m=grid.nx*j+i;
	int flag=0;
        vz0[m] = defaultz0_;
        double t,p;
        grid.xy(i,j,t,p);
        for (const_iterator polygon = this->begin(); polygon != this->end(); polygon++) {
          const PointBase<Coordinate> p0 = polygon->at(0);
	  flag++;
          if(t<p0.x()-180.) t+=360.;
          if(t>p0.x()+180.) t-=360.;
          if (polygon->where(PointBase<Coordinate> (static_cast<Coordinate>(t), 
                                                    static_cast<Coordinate>(p))) != out) {
            vz0[m] = flag;
            }
          }
        }
      }

    return vz0;
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <class Coordinate, class Z0>
  std::vector<Z0> PolygonSetBase<Coordinate, Z0>::determineFlag (const mesh_t &mesh) const

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/*------------------------------------------------------------------------------
    Number of vertices in the mesh. */
    const size_t N = mesh.nvtxs;
/*------------------------------------------------------------------------------
    Contains the value of z0 associated to each point of the mesh. */
    std::vector<Z0> vz0 (N);

    for (size_t i = 0; i < N; i++) {
//      vz0[i] = defaultz0_;
      int flag=0;
      vz0[i] = 0;
      double t=mesh.vertices[i].lon;
      double p=mesh.vertices[i].lat;
      for (const_iterator polygon = this->begin(); polygon != this->end(); polygon++) {
        const PointBase<Coordinate> p0 = polygon->at(0);
        flag++;
        if(t<p0.x()-180.) t+=360.;
        if(t>p0.x()+180.) t-=360.;

	if (polygon->where(PointBase<Coordinate> (static_cast<Coordinate>(t), static_cast<Coordinate>(p))) != out)
	  vz0[i] = 1;
// 	  vz0[i] = flag;
          }
      }

    return vz0;
  }

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <class Coordinate, class Z0>
  std::vector<Z0> PolygonSetBase<Coordinate, Z0>::determineFlagCartesian (const mesh_t &mesh, projPJ PJ) const

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/*------------------------------------------------------------------------------
    Number of vertices in the mesh. */
    const size_t N = mesh.nvtxs;
/*------------------------------------------------------------------------------
    Contains the value of z0 associated to each point of the mesh. */
    std::vector<Z0> vz0 (N);
    
/*------------------------------------------------------------------------------
    temporary polygon set in cartesian coordinates */
    PolygonSetBase<Coordinate, Z0> tmp;
    
    for (const_iterator polygon = this->begin(); polygon != this->end(); polygon++) {
      PolygonBase<Coordinate, Z0> p=*polygon;
      for (int k=0; k<p.size(); k++) {
        PointBase<Coordinate> *point = &(p[k]);
        double tt=point->x();
        double pp=point->y();
        double x,y;
        geo_to_projection(PJ, pp, tt, &x, &y);
        const PointBase<Coordinate> cartesian(x,y);
        *point = cartesian;
        }
      tmp.push_back(p);
      }

    for (size_t i = 0; i < N; i++) {
      vz0[i] = 0;
      double t=mesh.vertices[i].lon;
      double p=mesh.vertices[i].lat;

      for (const_iterator polygon = tmp.begin(); polygon != tmp.end(); polygon++) {
        const PointBase<Coordinate> p0 = polygon->at(0);
//         if(t<p0.x()-180.) t+=360.;
//         if(t>p0.x()+180.) t-=360.;
        double x,y;
        geo_to_projection(PJ, p, t, &x, &y);
        if (polygon->where(PointBase<Coordinate> (static_cast<Coordinate>(x), static_cast<Coordinate>(y))) != out)
          vz0[i] = 1;
          }
      }

    tmp.clear();
    
    return vz0;
  }

  /**
   * \brief Set the value of z0 in the mesh.
   * \param polygonName Name of the file containing the descriptions of the
   * polygons.
   * \param pointName Name of the file containing the reference points.
   * \param mesh The mesh itself.
   * \param defaultZ0 The default value of z0 in the mesh.
   * \return The value of z0 for each node of the mesh.
   *
   * Take a mesh file name, a reference point file name and a mesh and set the
   * values of z0 in the mesh.
   */

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <class Coordinate, class Z0>
  static std::vector<Z0> initialiseValue (const std::string &polygonName,
                                          const std::string &pointName,
                                          const mesh_t &mesh, const Z0 &defaultZ0, 
                                          const char *proj4_options, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{ 
  int status;
  bool debug=false;
  int maxsize=-1;
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    
    Read the file named polygonName, which contains a graph, and search the
    polygons in
    
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

    GraphBase<Coordinate> graph=GraphBase<Coordinate>(polygonName);
    if(debug) printf("found %d nodes in %s\n",graph.size(),polygonName.c_str());
    
/*------------------------------------------------------------------------------
    next call is problematic when re-entering the polygons seek procedures */
//    PolygonSetBase<Coordinate, Z0> polygons =searchPolygons<Coordinate, Z0>(graph,maxsize);

/*------------------------------------------------------------------------------
    this one is safer, to be fixed */
    PolygonSetBase<Coordinate, Z0> polygons=searchPolygons<Coordinate, Z0>(graph);
    if(verbose==1) printf("found %d polygons in %s\n",polygons.size(),polygonName.c_str());
    polygons.recale(nan(""));

/*------------------------------------------------------------------------------
    set the default z0 in the polygon set */
    polygons.setDefaultValue(defaultZ0);

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    
    assign value/flag following polygons partition
    
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

    if(pointName!="NONE") {
/*------------------------------------------------------------------------------
      read the file named pointName, which contains reference points */
      static vector< ReferencePoint<Coordinate, Z0> > prescribed=Polygons::readReferencePoints<Coordinate, Z0> (pointName);
/*------------------------------------------------------------------------------
      set the value of z0 in the associated polygons */
      status=polygons.setz0(prescribed);
      string tmp="prescribed-values.plg";
      polygons.dump(tmp.c_str());
      if(status!=0) {
/*------------------------------------------------------------------------------
        dump final set of cycles*/
        polygons.dump("failed.plg");
        check_error(-1, "polygon method, value setting failed", __LINE__, __FILE__, 1) ;
        }
/*------------------------------------------------------------------------------
      compute the z0 for each vertex in the mesh */
      std::vector<Coordinate> vz0 = polygons.determineValue(mesh);
      return vz0;
      }
    else {
      std::vector<Coordinate> vz0;
      if(proj4_options==0) {
        vz0 = polygons.determineFlag(mesh);
        }
      else {
        projPJ PJ=init_projection(proj4_options, false, 0);
        vz0 = polygons.determineFlagCartesian(mesh, PJ);
        }
      return vz0;
      }
    graph.clear();
  }

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <class Coordinate, class Z0>
  static std::vector<Z0> initialiseValue (const std::string &polygonName,
                                          const std::string &pointName,
                                          const mesh_t &mesh, const Z0 &defaultZ0)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{ 
  int status;
    int maxsize=-1;
/**----------------------------------------------------------------------------
    Read the file named polygonName, which contains a graph, and search the
    polygons in */
    GraphBase<Coordinate> graph=GraphBase<Coordinate>(polygonName);
    printf("found %d nodes in %s\n",graph.size(),polygonName.c_str());
    
/**----------------------------------------------------------------------------
    next call is problematic when re-entering the plogyns seek procedures */
//    PolygonSetBase<Coordinate, Z0> polygons =searchPolygons<Coordinate, Z0>(graph,maxsize);
/**----------------------------------------------------------------------------
    this one is safer, to be fixed */
    PolygonSetBase<Coordinate, Z0> polygons =searchPolygons<Coordinate, Z0>(graph);
    printf("found %d polygons in %s\n",polygons.size(),polygonName.c_str());

/**----------------------------------------------------------------------------
    Set the default z0 in the polygonSet */
    polygons.setDefaultValue(defaultZ0);

/**----------------------------------------------------------------------------
    Dump final set of cycles*/
//    polygons.dump();

/**----------------------------------------------------------------------------
    Read the file named pointName, which contains reference points, and set
    the value of z0 in the associated polygons. */
    if(pointName!="NONE") {
      status=polygons.setz0(Polygons::readReferencePoints<Coordinate, Z0> (pointName));
      if(status!=0) check_error(-1, "polygon method, value setting failed", __LINE__, __FILE__, 1) ;
/**----------------------------------------------------------------------------
      compute the z0 for each vertex in the mesh */
      std::vector<Coordinate> vz0 = polygons.determineValue(mesh);
      return vz0;
      }
    else {
      std::vector<Coordinate> vz0 = polygons.determineFlag(mesh);
      return vz0;
      }
    graph.clear();
  }
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <class Coordinate, class Z0> 
  static std::vector<Z0> initialiseValue (const std::string &polygonFilename,
                                          const std::string &pointFilename,
                                          const grid_t &grid, const Z0 &defaultZ0)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{ 
  int status;
//     int maxsize=-1;
    bool debug=false;
    
/*------------------------------------------------------------------------------
    Read the file named polygonName, which contains a graph, and search the
    polygons in */
    GraphBase<Coordinate> graph=GraphBase<Coordinate>(polygonFilename);
    if(debug) printf("found %d nodes in %s\n",graph.size(),polygonFilename.c_str());
    
/**----------------------------------------------------------------------------
    next call is problematic when re-entering the polygons seek procedures */
//    PolygonSetBase<Coordinate, Z0> polygons =searchPolygons<Coordinate, Z0>(graph,maxsize);
/**----------------------------------------------------------------------------
    this one is safer, to be fixed */
    PolygonSetBase<Coordinate, Z0> polygons =searchPolygons<Coordinate, Z0>(graph);
    if(debug) printf("found %d polygons in %s\n",polygons.size(),polygonFilename.c_str());

/*------------------------------------------------------------------------------
    Set the default z0 in the polygonSet */
    polygons.setDefaultValue(defaultZ0);

/*------------------------------------------------------------------------------
    Dump final set of cycles*/
//    polygons.dump();

/*------------------------------------------------------------------------------
    Read the file named pointName, which contains reference points, and set
    the value of z0 in the associated polygons. */
    if(pointFilename!="NONE") {
//       status=polygons.setz0(Polygons::readReferencePoints<Coordinate, Z0> (pointFilename));
      static vector< ReferencePoint<Coordinate, Z0> > prescribed=Polygons::readReferencePoints<Coordinate, Z0> (pointFilename);
/*------------------------------------------------------------------------------
      set the value of z0 in the associated polygons */
      status=polygons.setz0(prescribed);
/**----------------------------------------------------------------------------
      compute the z0 for each vertex in the mesh */
      std::vector<Z0> vz0 = polygons.determineValue(grid);
      string tmp="prescribed-values.plg";
      polygons.dump(tmp.c_str());
      if(status!=0) check_error(-1, "polygon method, value setting failed", __LINE__, __FILE__, 1) ;
      return vz0;
      }
    else {
      std::vector<Z0> vz0 = polygons.determineFlag(grid);
      return vz0;
      }
    graph.clear();
  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <class Coordinate, class Z0>
  PolygonSetBase<Coordinate, Z0> load (const std::string &polygonFilename, const std::string &pointFilename)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{ 
  int status;
  int maxsize=-1;
/*------------------------------------------------------------------------------
  Read the file named polygonName, which contains a graph, and search the
  polygons in */
  GraphBase<Coordinate> graph=GraphBase<Coordinate>(polygonFilename);
  printf("found %d nodes in %s\n",graph.size(),polygonFilename.c_str());
    
  PolygonSetBase<Coordinate, Z0> polygons =searchPolygons<Coordinate, Z0>(graph,maxsize);
  printf("found %d polygons in %s\n",polygons.size(),polygonFilename.c_str());

/*------------------------------------------------------------------------------
    Dump final set of cycles*/
//    polygons.dump();

/*------------------------------------------------------------------------------
  Read the file named pointName, which contains reference points, and set
  the value of z0 in the associated polygons. */
  if(pointFilename!="NONE") {
    status=polygons.setz0(Polygons::readReferencePoints<Coordinate, Z0> (pointFilename));
    if(status!=0) check_error(-1, "polygon method, value setting failed", __LINE__, __FILE__, 1) ;
    }
  graph.clear();
  return(polygons);
}

}

#endif
