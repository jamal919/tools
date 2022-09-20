#ifndef POLYGONS_EXCEPTIONS_HPP
#define POLYGONS_EXCEPTIONS_HPP

/**
 * \file
 * \brief Header which contains all classes for exceptions gestion for
 * polygons.
 * \author Yoann LE BARS
 * \version 2.0
 * \date August 9th, 2007
 * \date April 10th, 2008
 * \date May 26th, 2008
 * \date June 17th, 2008
 *
 * Copyright: see COPYING file.
 */

#include <exception>
#include <string>
#include <sstream>

/**
 * \namespace PolygonsExceptions
 *
 * \brief Namespace for all the exceptions that can happend with polygons.
 */
namespace PolygonsExceptions {
  /**
   * \class GraphError
   * \brief Base class for all exceptions with graph.
   */
  class GraphError: public std::exception {
  public:
    const char *what (void) const throw () = 0;
  };  // class GraphError
  /**
   * \class CapacityOverflow
   * \brief Generated when trying to acces out of range data.
   */
  class CapacityOverflow: public GraphError {
  public:
    const char *what (void) const throw () {
      return "capacity overflow";
    }  // const char *what (void) const
  };  // class CapacityOverflow
  /**
   * \class Negative
   * \brief Generated when a value that should be positive is negative.
   */
  class Negative: public GraphError {
  public:
    const char *what (void) const throw () {
      return "a value that should be positive is negative";
    }  // const char *what (void) const
  };  // class Negative
  /**
   * \class ReadError
   * \brief Generated when an error occurs during the reading of a file.
   */
  class ReadError: public GraphError {
  private:
    /// Name of the file that generates the error.
    std::string _fileName;

  public:
    /**
     * \brief Constructor
     * \param fileName_ Name of the invalid file.
     */
    explicit ReadError (const std::string &fileName_) throw ():
      _fileName (fileName_) {}  // ReadError (const std::string &)
    /**
     * \brief Destructor
     */
    ~ReadError (void) throw () {}

    /**
     * \brief Indicates which file is wrong.
     * \return The name of the invalid file.
     */
    std::string fileName (void) const throw () {return _fileName;}

    const char *what (void) const throw () {
      std::ostringstream oss;
      oss << "impossible to read in file: \"" << _fileName << '\"';
      return oss.str().c_str();
    }  // const char *what (void) const
  };  // class ReadError
  /**
   * \class OnBoundaryPoint
   * \brief Generated when a point in on the boundary of a graph.
   */
  class OnBoundaryPoint: public GraphError {
  private:
    /// Abscissa of the point.
    double x_;
    /// Ordinate of the point.
    double y_;
    /// Label of the point.
    size_t ind_;

  public:
    /**
     * \brief Constructor
     * \param _x Abscissa of the point.
     * \param _y Ordinate of the point.
     * \param _ind Label of the point.
     */
    OnBoundaryPoint (double _x, double _y, size_t _ind) throw ():
      x_ (_x), y_ (_y), ind_ (_ind) {
    }  // OnBoundaryEdge (double, double, size_t)
    /**
     * \brief Destructor
     */
    ~OnBoundaryPoint (void) throw () {}

    /**
     * \brief Element access.
     * \return Abscissa of the point.
     */
    double x (void) const throw () {return x_;}
    /**
     * \brief Element access.
     * \return Ordinate of the point.
     */
    double y (void) const throw () {return y_;}
    /**
     * \brief Element access.
     * \return Label of the point.
     */
    size_t ind (void) const throw () {return ind_;}

    const char *what (void) const throw () {
      std::ostringstream oss;
      oss << "the point " << ind_ << " (" << x_ << ", " << y_
          << ") is on the boundary of a polygon";
      return oss.str().c_str();
    }  // const char *what (void) const
  };  //class OnBoundaryPoint
}  // namespace PolygonsExceptions

#endif  // #ifndef POLYGONS_EXCEPTIONS_HPP
