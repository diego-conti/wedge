/***************************************************************************
 *   Copyright (C) 2007, 2008, 2009 by Diego Conti, diego.conti@unimib.it  *
 *                                                                         *
 *   This file is part of Wedge.                                           *
 *                                                                         *
 *   Wedge is free software; you can redistribute it and/or modify         *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 3 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   Wedge is distributed in the hope that it will be useful,              *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with Wedge; if not, write to the                                *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#ifndef WEXCEPTION_H
#define WEXCEPTION_H
/** @ingroup Base */
  
/** @{ 
 * @file wexception.h
 * @brief Exception template classes
 */

#include <stdexcept>
#include <string>
#include <ginac/ex.h>

namespace Wedge {

namespace internal {
std::string ToString(int n);
}
	
/** @brief Base template class for exceptions in Wedge
 * @param T An exception class, derived from std::exception
 */
template<typename T> class WedgeException : public T {
public:
/** @brief Create an exception object
 * @param what_arg A string describing the error
 * @param in_file The name of the file, as returned by __FILE__ 
 * @param at_line The line number, as returned by __LINE__  
 */
	WedgeException(const std::string& what_arg, const char* in_file, int at_line) throw() :
		T (std::string("Error in wedge0.0/")+in_file+", line "+internal::ToString(at_line)+" : "+ what_arg) {}
	virtual ~WedgeException() throw() {}
};


/** @brief Exception for out of range errors 
 */
class OutOfRange : public WedgeException<std::out_of_range> {
public:
/** @brief Create an exception object
 * @param number The out of range value
 * @param in_file The name of the file, as returned by __FILE__ 
 * @param at_line The line number, as returned by __LINE__  
 */
	OutOfRange(const char* in_file, int at_line, int number) throw() :
		WedgeException<std::out_of_range>("Value " + internal::ToString(number) + " is out of range",in_file,at_line) {}
/** @overload
 */ 
	OutOfRange(const char* in_file, int at_line) throw() :
		WedgeException<std::out_of_range>("Out of range",in_file,at_line) {} 
};


/** @brief Exception thrown by functions presently lacking an implementation 
 */
class NotImplemented : public WedgeException<std::logic_error> {
public:
/** @brief Create an exception object
 * @param function_name The name of the function not implemented
 * @param in_file The name of the file, as returned by __FILE__ 
 * @param at_line The line number, as returned by __LINE__  
 */
	NotImplemented(const char* in_file, int at_line,const char* function_name) : 
		WedgeException<std::logic_error>(std::string("Function ")+function_name+" not implemented",in_file, at_line) {} 
/** @overload
 */ 
	NotImplemented(const char* in_file, int at_line) : 
		WedgeException<std::logic_error>("Function not implemented",in_file, at_line) {}
};


/** @brief Exception for invalid argument errors 
 */

class InvalidArgument : public WedgeException<std::invalid_argument> {
public:
	GiNaC::ex arg;	///< The invalid argument
/** @brief Create an exception object
 * @param arg The invalid argument
 * @param in_file The name of the file, as returned by __FILE__ 
 * @param at_line The line number, as returned by __LINE__  
 */
	InvalidArgument(const char* in_file, int at_line,GiNaC::ex arg) throw();		
/** @overload
 */
	InvalidArgument(const char* in_file, int at_line) throw() :
		WedgeException<std::invalid_argument>("Invalid argument passed to function ",in_file,at_line) {}
	virtual ~InvalidArgument() throw() {}
};


/** @brief Exception thrown when a declaration is found inconsistent 
 */
class InconsistentDeclaration : public WedgeException<std::runtime_error> {
public:
/** @brief Create an exception object
 * @param declaring_what A description of the inconsistent declaration
 * @param in_file The name of the file, as returned by __FILE__ 
 * @param at_line The line number, as returned by __LINE__  
 */
	InconsistentDeclaration(const char* in_file, int at_line,const std::string& declaring_what) : 
		WedgeException<runtime_error>("Declaration of "+declaring_what+" is inconsistent with preexisting conditions", in_file, at_line) {}
};

/** @brief Exception thrown by ConcreteManifold when the dimension is less than one*/
struct DiscreteManifold : public WedgeException<std::runtime_error> {
	DiscreteManifold(const char* in_file, int at_line, int dimension) :   
		WedgeException<std::runtime_error>("Cannot create a manifold of dimension"+ internal::ToString(dimension), in_file, at_line) {}
};


} /** @} */

#endif
