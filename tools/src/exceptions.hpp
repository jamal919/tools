/*****************************************************************************
*                      E X C E P T I O N S . H P P                           *
* -------------------------------------------------------------------------- *
* Header which refer to all classes for exceptions gestion.                  *
* -------------------------------------------------------------------------- *
* Created: december 7th, 2006                                                *
* Last modification: june 11th, 2007                                         *
* Contribuers: Yoann LE BARS                                                 *
*****************************************************************************/

#ifndef EXCEPTIONS_HPP
#define EXCEPTIONS_HPP

/**
 * \file
 * \brief Header which refer to all classes for exceptions gestion.
 * \author Yoann LE BARS
 * \version 1.0
 * \date December 7th, 2006
 * \date June 11th, 2007
 * \date August 13th, 2008
 */

#include <exception>

#include "polygons_exceptions.hpp"

/**
 * \namespace TugoExceptions
 * \brief Namespace containing classes for general exceptions.
 */
namespace TugoExceptions {
  // Namespace composition
  using namespace PolygonsExceptions;

  // --- Exceptions that can occur during the initialization process. --------

  /**
   * \brief Base class for exceptions that can occur during the
   * initialization process.
   *
   * Purely virtual: no instance of this class.
   */
  class InitExceptions: public std::exception {
  public:
    /**
     * \brief Describe the exception.
     * \return A C style char string.
     *
     * Purely virtual
     */
    virtual const char *what (void) throw () = 0;
  };  // class InitExceptions

  /// \brief Exception that occurs when two contents are incompatible.
  class MismatchingContent: public InitExceptions {
  public:
    /**
     * \brief Describe the exception.
     * \return A C style char string.
     */
    virtual const char *what (void) throw () {
      return "mismatching contents";
    }  // const char *what (void)
  };  // class MismatchingContent
}  // namespace TugoExceptions

#endif  //#ifndef EXCEPTIONS_HPP
