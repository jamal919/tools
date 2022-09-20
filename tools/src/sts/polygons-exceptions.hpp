/******************************************************************************
*               P O L Y G O N S -  E X C E P T I O N S . H P P                *
* --------------------------------------------------------------------------- *
* Header which contains classes for exceptions gestion for polygons.          *
* --------------------------------------------------------------------------- *
* Created: august 9th, 2007                                                   *
* Last modification: august 9th, 2007                                         *
* Contribuers: Yoann LE BARS                                                  *
******************************************************************************/

#ifndef POLYGONS_EXCEPTIONS_HPP
#define POLYGONS_EXCEPTIONS_HPP

#include <string>

namespace PolygonsExceptions {
  class GraphError {};
  class CapacityOverflow: public GraphError {};
  class Negative:         public GraphError {};
  class ReadError:        public GraphError {
  private:
    std::string _fileName;

  public:
    ReadError (const std::string &fileName_) throw (): _fileName (fileName_) {}

    std::string fileName (void) const throw () {return _fileName;}
  };  // class ReadError
  class OnBoundaryPoint:  public GraphError {
  private:
    double x_, y_;
    size_t ind_;

  public:
    OnBoundaryPoint (double _x, double _y, size_t _ind) throw ():
      x_ (_x), y_ (_y), ind_ (_ind) {
    }  // OnBoundaryEdge (double, double, size_t)
    double x (void) const throw () {return x_;}
    double y (void) const throw () {return y_;}
    size_t ind (void) const throw () {return ind_;}
  };  //class OnBoundaryPoint
}  // namespace PolygonsExceptions

#endif  // #ifndef POLYGONS_EXCEPTIONS_HPP
