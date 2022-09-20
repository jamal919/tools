
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

/******************************************************************************
* Header for polygons seeking.                                                *
* --------------------------------------------------------------------------- *
* Created: march 19th, 2007                                                   *
* Last modification: august 9th, 2007                                         *
* Contribuers: Yoann LE BARS                                                  *
******************************************************************************/

#ifndef POLYGONS_HPP
#define POLYGONS_HPP

/**
 * \file
 * \brief Header for polygons seeking.
 * \author Yoann LE BARS
 * \version 2.2
 * \date February 18th, 2007
 * \date March 19th, 2007
 * \date June 6th, 2007
 * \date April 9th, 2008
 * \date April 10th, 2008
 * \date April 14th, 2008
 * \date May 15th, 2008
 * \date May 26th, 2008
 * \date June 6th, 2008
 * \date August 11th, 2008
 * \date August 28th, 2008
 *
 * Copyright: see COPYING file.
 */

#include <list>
#include <vector>
#include <string>
#include <assert.h>
#include <iostream>

#include "polygons_exceptions.hpp"
#include "fe.h"
#include "constants.h"

#include <proj_api.h>
/**
 * \namespace Polygons
 * \brief Namespace for polygons gestion.
 */
namespace Polygons {
  using PolygonsExceptions::GraphError;
  using PolygonsExceptions::CapacityOverflow;
  using PolygonsExceptions::ReadError;
  using PolygonsExceptions::OnBoundaryPoint;
  using PolygonsExceptions::Negative;

  /**
   * \typedef std::list<size_t> Adjacency
   * \brief Describes the points a given point of a graph is adjacent with.
   */
  typedef std::list<size_t> Adjacency;

  /**
   * \enum Position
   * \brief Indicates the position of a point relatively to a polygon.
   */
  enum Position {
    /// The point is strictly out of the polygon.
    out,
    /// The point is on an edge of the polygon, but not on a vertex.
    edge,
    /// The point is a vertex of the polygon.
    vertex,
    /// The point is strictly in the polygon.
    in
  };

  /**
   * \class PointBase
   * \brief Representation of a two dimensional point.
   *
   * The point coordinates can be of any type.
   */
template <class T> class PointBase {
  private:
    /// Abscissa of the point.
    T abs;
    /// Ordinate of the point.
    T ord;

  public:
    /**
     * \brief Default constructor.
     * \param x_ Abscissa of the point.
     * \param y_ Ordinate of the point.
     *
     * If abscissa and ordinate are note given, then their values are the one
     * of the default constructor of the type T.
     */
    explicit PointBase (const T &x_ = T (), const T &y_ = T ()):
      abs (x_), ord (y_) {}  // PointBase (const T & = T (), const T & = T ())
    /**
     * \brief Copy constructor.
     * \param p The point to be copied.
     */
    PointBase (const PointBase<T> &p): abs (p.abs), ord (p.ord) {}

    /**
     * \brief Controlled access to the abscissa.
     * \return The value of the abscissa.
     *
     * To modify the abscissa is possible.
     */
    T &x (void) {return abs;}
    /**
     * \brief Controlled access to the ordinate.
     * \return The value of the ordinate.
     *
     * To Modify the ordinate is possible.
     */
    T &y (void) {return ord;}
    /**
     * \brief Constant controlled access to the abscissa.
     * \return The value of the abscissa.
     *
     * To modify the abscissa is not possible.
     */
    T x (void) const {return abs;}
    /**
     * \brief Constant controlled access to the ordinate.
     * \return The value of the ordinate.
     *
     * To Modify the ordinate is not possible.
     */
    T y (void) const {return ord;}

    /**
     * \brief Copy operator.
     * \param p The point to be copied.
     * \return The current object.
     *
     * The current object receive the value to be copied.
     */
    PointBase<T> &operator = (const PointBase<T> &p) {
      abs = p.abs;
      ord = p.ord;
      return *this;
      } 
      
    /**
     * \brief Tabular access operator.
     * \param i Label of the asked element.
     * \return The value of the asked element.
     *
     * 0 is the abscissa, 1 the ordinate.
     * To modify the value is possible.
     */
    T &operator [] (size_t i) {
      assert (i < 2);
      switch (i) {
        case 0: return abs;
        case 1: return ord;
      }
      } 
    
    /**
     * \brief Constant tabular access operator.
     * \param i Label of the asked element.
     * \return The value of the asked element.
     *
     * 0 is the abscissa, 1 the ordinate.
     * To modify the value is not possible.
     */
    T operator [] (size_t i) const {
      assert (i < 2);
      switch (i) {
        case 0: return abs;
        case 1: return ord;
      }
      } 
  };

  /**
   * \fn bool operator == (const PointBase<T> &p1, const PointBase<T> &p2)
   * \brief Equality operator overloading.
   * \param p1 First PointBase to be compared.
   * \param p2 Second PointBase to be compared.
   * \return True if the two point are identical, else false.
   */
template <class T>
  inline bool operator == (const PointBase<T> &p1, const PointBase<T> &p2) {
    return (p1.x() == p2.x()) && (p1.y() == p2.y());
    }
    
  /**
   * \fn bool operator != (const PointBase<T> &p1, const PointBase<T> &p2)
   * \brief Non-equality operator overloading.
   * \param p1 First PointBase to be compared.
   * \param p2 Second PointBase to be compared.
   * \return True if the two point base are different, else false.
   */
  template <class T>
  inline bool operator != (const PointBase<T> &p1, const PointBase<T> &p2) {
    return (p1.x() != p2.x()) || (p1.y() != p2.y());
    }  
    
  /**
   * \fn PointBase<T> operator + (const PointBase<T> &p1, const PointBase<T>
   * &p2)
   * \brief Addition operator overloading.
   * \param p1 First PointBase.
   * \param p2 PointBase to be added to the first one.
   * \return A new PointBase, which contains the results of the addition.
   */
  template <class T>
  inline PointBase<T> operator + (const PointBase<T> &p1,
                                  const PointBase<T> &p2) {
    return PointBase<T> (p1.x() + p2.x(), p1.y() + p2.y());
    }
    
  /**
   * \fn PointBase<T> &operator += (PointBase<T> &p1, const PointBase<T>
   * &p2)
   * \brief Inplace addition operator overloading.
   * \param p1 First PointBase.
   * \param p2 PointBase to be added to the first one.
   * \return A new PointBase, which contains the results of the addition.
   */
  template <class T>
  inline PointBase<T> &operator += (PointBase<T> &p1, const PointBase<T> &p2) {
    p1.x() += p2.x();
    p1.y() += p2.y();
    return p1;
    }  
  
  /**
   * \fn PointBase<T> operator - (const PointBase<T> &p1, const PointBase<T> &p2)
   * \brief Substraction operator overloading.
   * \param p1 First PointBase.
   * \param p2 PointBase to be substract to the first one.
   * \return A new PointBase, which contains the results of the substraction.
   */
  template <class T>
  inline PointBase<T> operator - (const PointBase<T> &p1,
                                  const PointBase<T> &p2) {
    return PointBase<T> (p1.x() - p2.x(), p1.y() - p2.y());
    }  
    
  /**
   * \fn PointBase<T> &operator -= (PointBase<T> &p1, const PointBase<T>
   * &p2)
   * \brief Inplace substraction operator overloading.
   * \param p1 First PointBase.
   * \param p2 PointBase to be added to the first one.
   * \return A new PointBase, which contains the results of the addition.
   */
  template <class T>
  inline PointBase<T> &operator -= (PointBase<T> &p1, const PointBase<T> &p2) {
    p1.x() -= p2.x();
    p1.y() -= p2.y();
    return p1;
    }
  
  /**
   * \brief Multiplication operator overloading.
   * \param a The scalar to multiply with.
   * \param p The point to be multiply.
   * \return p times a.
   */
  template <class T>
  inline PointBase<T> operator * (const T &a, const PointBase<T> &p) {
    return PointBase<T> (a * p.x(), a * p.y());
    }
  
  /**
   * \brief Inplace multiplication operator overloading.
   * \param a The scalar to multiply with.
   * \param p The point to be multiply.
   * \return p after multiplication.
   */
  template <class T>
  inline PointBase<T> &operator *= (PointBase<T> &p, const T &a) {
    p.x() *= a;
    p.y() *= a;
    return p;
    }
    
  /**
   * \brief Division operator overloading.
   * \param p The point to be divided.
   * \param a The scalar to divide with.
   * \return p divided by a.
   */
  template <class T>
  inline PointBase<T> operator / (const PointBase<T> &p, const T &a) {
    return PointBase<T> (p.x() / a, p.y() / a);
    }
    
  /**
   * \brief Inplace division operator overloading.
   * \param a The scalar to divide with.
   * \param p The point to be divided.
   * \return The point after the division.
   */
  template <class T>
  inline PointBase<T> &operator /= (PointBase<T> &p, const T &a) {
    p.x() /= a;
    p.y() /= a;
    return p;
    } 
    
  /**
   * \brief Operator << overloading.
   * \param os The output stream.
   * \param point The Point for the output.
   * \return The modified stream.
   */
  template <class T>
  std::ostream &operator << (std::ostream &os, const PointBase<T> &point) {
    os << '(' << point.x() << ", " << point.y() << ')';
    return os;
    }
    
  /**
   * \brief Operator >> overloading.
   * \param is The input stream.
   * \param point The Point for the input.
   * \return The modified stream.
   */
  template <class T>
  std::istream &operator >> (std::istream &is, PointBase<T> &point) {
    is >> point.x() >> point.y();
    return is;
    }

  /**
   * \class ElementBase
   * \brief Class to represent an element of a graph.
   *
   * An element of a graph is a point which has some adjacent points.
   */
  template <class T> class ElementBase {
  private:
    /// The point itself.
    PointBase<T> p;
    /// List of its adjacent points.
    Adjacency adj;

  public:
    /**
     * \brief Default constructor.
     * \param _p Coordinate of the element.
     * \param _adj List of the points adjacent to the element.
     *
     * If no point is given, then an empty element is created.
     */
    explicit ElementBase (const PointBase<T> &_p = PointBase<T> (),
                          const Adjacency &_adj = Adjacency ()):
      p (_p), adj (_adj) {
      }
      
    /**
     * \brief Copy constructor.
     * \param el The element to be copied.
     */
    ElementBase (const ElementBase<T> &el): p (el.p), adj (el.adj) {}

    /**
     * \brief Controlled access to the coordinates of the element.
     * \return The coordinates.
     *
     * To modify the coordinates is possible.
     */
    PointBase<T> &coordinates (void) {return p;}
    /**
     * \brief Controlled access to the adjacency list of the element.
     * \return The List.
     *
     * To modify the adjacency list is possible.
     */
    Adjacency &adjacency (void) {return adj;}
    /**
     * \brief Controlled access to the coordinates of the element.
     * \return The coordinates.
     *
     * To modify the coordinates is not possible.
     */
    PointBase<T> coordinates (void) const {return p;}
    /**
     * \brief Controlled access to the adjacency list of the element.
     * \return The List.
     *
     * To modify the adjacency list is not possible.
     */
    Adjacency adjacency (void) const {return adj;}

    /**
     * \brief Copy operator overloading.
     * \param el The element to be copied.
     * \return A reference to the Element.
     */
    ElementBase<T> &operator = (const ElementBase<T> &el) {
      p = el.p;
      adj = el.adj;
      return *this;
      } 
    }; 

  /**
   * \class GraphBase
   * \brief Class for graph representation.
   */
  template <class T> class GraphBase: public std::vector< ElementBase<T> > {
  public:
    /**
     * \brief Constructor that reads a graph in a file containing the graph.
     * \param fileName The name of the file.
     */
    explicit GraphBase (const std::string &);
    explicit GraphBase (const mesh_t &);
    };

/* *----------------------------------------------------------------------------
   Definition of PolygonBase */

  /**
   * \class PolygonBase
   * \brief Class for polygons gestion.
   *
   * A polygon is composed of points and a value for z0.
   */
  template <class Coordinate, class Z0> class PolygonBase: public std::vector< PointBase<Coordinate> > 
{
  private:
    /// Value of z0 in the polygon.
    Z0 _z0;
    /// The area of the polygon times 2.
    mutable Coordinate areax2;
    /// True if the area of the polygon had been computed, else false.
    mutable bool areaComputed;

    // Compute the area of the polygon times 2 and stock it in a2.
    void computeArea2 (void) const;

    /**
     * \typedef std::vector< PointBase<Coordinate> > ElementType
     * \brief A vector containing points.
     *
     * Simplify the code: it is more simple to read and understand and the
     * template gestion is facilitated.
     */
    typedef std::vector< PointBase<Coordinate> > ElementType;

  public:
    /**
     * \brief Default constructor.
     * \param _vert Vector containing the vetices of the polygon.
     * \param z0_ Value of z0 in the polygon.
     *
     * If no vector is given, then an empty polygon is created.
     */
    explicit PolygonBase (const ElementType &_vert = ElementType (),
                          const Z0 &z0_ = Z0 ()):
      ElementType (_vert), _z0 (z0_), areaComputed (false) {
      }
      
    /**
     * \brief Copy constructor.
     * \param p The polygon to be copied.
     */
    PolygonBase (const PolygonBase<Coordinate, Z0> &p):
      ElementType (p), _z0 (p._z0), areax2 (p.areax2), areaComputed (p.areaComputed) {
      }

    /**
     * \brief Controlled access to z0.
     * \return The value of z0.
     */
    Z0 z0 (void) const {return _z0;}
    
    /**
     * \brief Changes the value of z0 in the polygon.
     * \param value The new value of z0.
     */
    void setz0 (const Z0 &value) {
//      if (value < 0) throw Negative ();
      _z0 = value;
      }

    // Gives the position of a point relatively to the polygon.
    Position where (const PointBase<Coordinate> &) const;

    /**
     * \brief Gives the area of the polygon times 2.
     * \return The area times 2.
     */
    Coordinate area2 (void) const {
      if (!areaComputed) computeArea2();
      return areax2;
      }

    /**
     * \brief Controlled access to the vertices of the polygon.
     * \param i The label of the vertex.
     * \return The vertex as a Point.
     *
     * Modification of the vertex is possible.
     */
    PointBase<Coordinate> &at (size_t i) {
      areaComputed = false;
      return ElementType::at(i);
      }
    
    /**
     * \brief Controlled access to the vertices of the polygon.
     * \param i The label of the vertex.
     * \return The vertex as a Point.
     *
     * Modification of the vertex is not possible.
     */
    PointBase<Coordinate> at (size_t i) const {return ElementType::at(i);}
    
    /**
     * \brief Add a vertex to the polygon, at the tail.
     * \param p The vertex to add.
     */
    void push_back (const PointBase<Coordinate> &p) throw () {
      ElementType::push_back(p);
      areaComputed = false;
      } 

    /**
     * \brief Uncontrolled access operator overloading.
     * \param i Label of the vertex.
     * \return The vertex as a Point.
     *
     * Modification of the vertex is possible.
     */
    PointBase<Coordinate> &operator [] (size_t i) {
      areaComputed = false;
      return ElementType::operator[](i);
      }
    
    /**
     * \brief Uncontrolled access operator overloading.
     * \param i Label of the vertex.
     * \return The vertex as a Point.
     *
     * Modification of the vertex is not possible.
     */
    PointBase<Coordinate> operator [] (size_t i) const {
      return ElementType::operator[](i);
      }
      
    /**
     * \brief Copy operator overloading.
     * \param p The polygon to copy.
     * \return The current polygon, after modification.
     */
    PolygonBase &operator = (const PolygonBase<Coordinate, Z0> &p) {
      areax2 = p.areax2;
      _z0 = p._z0;
      areaComputed = p.areaComputed;
      ElementType::operator=(p);
      return *this;
      }
      
    void dumpCoord (FILE *, size_t);
    void recale(Coordinate, size_t);
  }; 

  /**
   * \brief Equality operator overloading.
   * \param p1 First polygon to be compared.
   * \param p2 Second polygon to be compared.
   * \return True if the two polygons are identical, else false.
   */
  template <class Coordinate, class Z0>
  inline bool operator == (const PolygonBase<Coordinate, Z0> &p1,
                           const PolygonBase<Coordinate, Z0> &p2) {
    return (p1.area2() == p2.area2())
           && (static_cast< std::vector< PointBase<Coordinate> > >(p1)
               == static_cast< std::vector< PointBase<Coordinate> > >(p1));
  }
  
  /**
   * \brief Unequality operator overloading.
   * \param p1 First polygon to be compared.
   * \param p2 Second polygon to be compared.
   * \return True if the two polygons are different, else false.
   */
  template <class Coordinate, class Z0>
  inline bool operator != (const PolygonBase<Coordinate, Z0> &p1,
                           const PolygonBase<Coordinate, Z0> &p2) {
    return (p1.area2() != p2.area2())
           || (static_cast< std::vector< PointBase<Coordinate> > >(p1)
               == static_cast< std::vector< PointBase<Coordinate> > >(p1));
  }
  
  /**
   * \brief Inferior operator overloading.
   * \param p1 First polygon to be compared.
   * \param p2 Second polygon to be compared.
   * \return True if area(p1) < area(p2), else false.
   */
  template <class Coordinate, class Z0>
  inline bool operator < (const PolygonBase<Coordinate, Z0> &p1,
                          const PolygonBase<Coordinate, Z0> &p2) {
    return p1.area2() < p2.area2();
  }
  
  /**
   * \brief Inferior or equal operator overloading.
   * \param p1 First polygon to be compared.
   * \param p2 Second polygon to be compared.
   * \return True if area(p1) <= area(p2), else false.
   */
  template <class Coordinate, class Z0>
  inline bool operator <= (const PolygonBase<Coordinate, Z0> &p1,
                           const PolygonBase<Coordinate, Z0> &p2) {
    return p1.area2() <= p2.area2();
  }
  
  /**
   * \brief Superior operator overloading.
   * \param p1 First polygon to be compared.
   * \param p2 Second polygon to be compared.
   * \return True if area(p1) > area(p2), else false.
   */
  template <class Coordinate, class Z0>
  inline bool operator > (const PolygonBase<Coordinate, Z0> &p1,
                          const PolygonBase<Coordinate, Z0> &p2) {
    return p1.area2() > p2.area2();
  }
  
  /**
   * \brief Superior or equal operator overloading.
   * \param p1 First polygon to be compared.
   * \param p2 Second polygon to be compared.
   * \return True if area(p1) >= area(p2), else false.
   */
  template <class Coordinate, class Z0>
  inline bool operator >= (const PolygonBase<Coordinate, Z0> &p1,
                           const PolygonBase<Coordinate, Z0> &p2) {
    return p1.area2() >= p2.area2();
  }

  /**
   * \class ReferencePoint
   * \brief Class for reference point gestion.
   *
   * A ReferencePoint is a point with an associated z0, which indicates the
   * value of z0 in the polygon in which the ReferencePoint lay.
   */
  template <class Coordinate, class Z0>
  class ReferencePoint: public PointBase<Coordinate> {
  private:
    /// Value of z0 in the localization given by the coordinates of the point.
    Z0 z0_;

  public:
    /// Default constructor.
    ReferencePoint (void) {z0_ = static_cast<Z0>(z0_def);}
    /**
     * \brief Constructs a ReferencePoint with a Point and a value of z0.
     * \param p The point.
     * \param _z0 The value for z0.
     */
    ReferencePoint (const PointBase<Coordinate> &p, const Z0 &_z0):
      PointBase<Coordinate> (p), z0_ (_z0) {
//      if (z0_ < 0) throw Negative ();
    }
    
    /**
     * \brief Constructs a ReferencePoint with an abscissa, an ordinate and a
     * value of z0.
     * \param _x abscissa
     * \param _y ordinate
     * \param _z0 Value of z0.
     */
    ReferencePoint (const Coordinate &_x, const Coordinate &_y,
                    const Z0 &_z0): PointBase<Coordinate> (_x, _y), z0_ (_z0) {
//      if (z0_ < 0) throw Negative ();
    }
    
    /**
     * \brief Copy constructor.
     * \param p The ReferencePoint to be copied.
     */
    ReferencePoint (const ReferencePoint<Coordinate, Z0> &p):
      PointBase<Coordinate> (p), z0_ (p.z0_) {
//      if (z0_ < 0) throw Negative ();
      }

    /**
     * \brief Changes the value of z0.
     * \param _z0 The new value for z0.
     */
    void setz0 (const Z0 &_z0) {
//      if (_z0 < 0) throw Negative ();
      z0_ = _z0;
      }
      
    /**
     * \brief Gives the value of z0 associated to the ReferencePoint.
     * \return The value.
     */
    Z0 z0 (void) const {return z0_;}

    /**
     * \brief Copy operator overloading.
     * \param p The ReferencePoint to be copied.
     * \return The current ReferencePoint, after being modified.
     */
    ReferencePoint<Coordinate, Z0> operator =
    (const ReferencePoint<Coordinate, Z0> &p) {
      z0_ = p.z0_;
      PointBase<Coordinate>::operator=(p);
      return *this;
      }
  };

/* *----------------------------------------------------------------------------
   Definition of PolygonSetBase */

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  /**
   * \class PolygonSetBase
   * \brief Data structure to stock all the polygons.
   */
  template <class Coordinate, class Z0>
  class PolygonSetBase: public std::list< PolygonBase<Coordinate, Z0> > {
  private:
    /// Default value of z0 in the PolygonSet.
    Z0 defaultz0_;
    /**
     * \typedef std::list< Polygon<Coordinate, Z0> > ElementType
     * \brief List containing the polygons of the PolygonSet.
     *
     * Simplify the code: it is more simple to read and understand and the
     * template gestion is facilitated.
     */
    typedef std::list< PolygonBase<Coordinate, Z0> > ElementType;

  public:
    /// Define the iterator type of PolygonSet.
    typedef typename ElementType::iterator iterator;
    /// Define the constant iterator type of PolygonSet.
    typedef typename ElementType::const_iterator const_iterator;

    /// Default constructor.
    PolygonSetBase (void) throw (): ElementType () {  // Default constructor
      defaultz0_ = static_cast<Z0>(z0_def);
      }
    /**
     * \brief Copy constructor.
     * \param p The PolygonSet to be copied.
     */
    PolygonSetBase (const PolygonSetBase<Coordinate, Z0> &p):
      ElementType (p), defaultz0_ (p.defaultz0_) {
      }

    // Determine the value of z0 in every polygons
    int setz0 (const std::vector< ReferencePoint<Coordinate, Z0> > &);

    // Determine the value of z0 for every points of a mesh
    std::vector<Z0> determineValue (const mesh_t &) const;

    // Determine the interior flag for every points of a mesh
    std::vector<Z0> determineFlag (const mesh_t &) const;

    // Determine the interior flag for every points of a mesh
    std::vector<Z0> determineFlagCartesian (const mesh_t &, projPJ PJ) const;

    // Determine the value of z0 for every points of a grid
    std::vector<Z0> determineValue (const grid_t &) const;

    // Determine the interior flag for every points of a mesh
    std::vector<Z0> determineFlag (const grid_t &) const;

    /**
     * \brief Set the default value of z0 in the PolygonSet.
     * \param _defaultZ0 The value to became the default z0.
     */
    void setDefaultValue (const Z0 &_defaultZ0) {
/**----------------------------------------------------------------------------
      commented, as used for various purposes */
//      if (_defaultZ0 < 0) throw Negative ();
      defaultz0_ = _defaultZ0;
      }
    /**
     * \brief Gives the default value of z0 in the PolygonSet.
     * \return The value.
     */
    Z0 defaultValue (void) const {return defaultz0_;}

    /**
     * \brief Gives the default value of z0 in the PolygonSet.
     * \return The value.
     */
    void dump(const char *);
    void recale(Coordinate);

    /**
     * \brief Equality operator overloading.
     * \param p The PolygonSet to be compared to the current one.
     * \return True if the two PolygonSet are equal, else false.
     */
    bool operator == (const PolygonSetBase<Coordinate, Z0> &p) const {
      return (defaultz0_ == p.defaultz0_) && ElementType::operator==(p);
      }
    /**
     * \brief Copy operator overloading
     * \param p The PolygonSet to be copied.
     * \return The current PolygonSet, after being modified.
     */
    PolygonSetBase &operator = (const PolygonSetBase<Coordinate, Z0> &p) {
      defaultz0_ = p.defaultz0_;
      ElementType::operator=(p);
      return *this;
      }
  } ;  /* template <class Coordinate = double, class Z0 = double>
  class PolygonSet */

  // Read the reference points in a file.
  template <class Coordinate, class Z0>
  std::vector< ReferencePoint<Coordinate, Z0> > readReferencePoints
  (const string &);

  // Search polygons in a graph.
  template <class Coordinate, class Z0>
  PolygonSetBase<Coordinate, Z0> searchPolygons (const GraphBase<Coordinate> &);

  template <class Coordinate, class Z0>
  PolygonSetBase<Coordinate, Z0> searchPolygons (const GraphBase<Coordinate> &, int maxsize);

  /* Take a mesh file name, a reference point file name and a mesh and set the
    values of z0 in the mesh. */
  template <class Coordinate, class Z0>
  std::vector<Z0> initialiseValue (const string &, const string &, const mesh_t &, const Z0 &, const char *proj4_options, int verbose=0);
  
  template <class Coordinate, class Z0>
  std::vector<Z0> initialiseValue (const string &, const string &, const mesh_t &, const Z0 &);

  template <class Coordinate, class Z0>
  PolygonSetBase<Coordinate, Z0> load (const std::string &polygonFilename, const std::string &pointFilename);

  /**
   * \typedef PointBase<double> Point
   * \brief Specialization of PointBase with type double for coordinates
   * representation.
   */
  typedef PointBase<double> Point;
  /**
   * \typedef ElementBase<double> Element
   * \brief Specialization of ElementBase with type double for coordinates
   * representation.
   */
  typedef ElementBase<double> Element;
  /**
   * \typedef GraphBase<double> Graph
   * \brief Specialization of PolygonBase with type double for coordinates
   * representation.
   */
  typedef GraphBase<double> Graph;
  /**
   * \typedef PolygonBase<double, double> Polygon
   * \brief Specialization of PolygonBase with type double for coordinates
   * and z0 representation.
   */
  typedef PolygonBase<double, double> Polygon;
  /**
   * \typedef PolygonSetBase<double, double> PolygonSet
   * \brief Specialization of PolygonSetBase with type double for coordinates
   * and z0 representation.
   */
  typedef PolygonSetBase<double, double> PolygonSet;
}

#include "polygons_implementation.hpp"
  
  extern void initializeLocalValue(const string &polygonsFileName, const string &pointsFileName,
                    const mesh_t &mesh, float *z0, double defaultZ0, const char*, int);
  

#endif  /* POLYGONS_HPP */
