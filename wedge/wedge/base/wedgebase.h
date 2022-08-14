/*******************************************************************************
 *  Copyright (C) 2007-2022 by Diego Conti, diego.conti@unimib.it 
 *  This file is part of Wedge.                                           
 *  Wedge is free software; you can redistribute it and/or modify         
 *  it under the terms of the GNU General Public License as published by  
 *  the Free Software Foundation; either version 3 of the License, or     
 *  (at your option) any later version.                                   
 *                                                                          
 *  Wedge is distributed in the hope that it will be useful,              
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of        
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         
 *  GNU General Public License for more details.                          
 *                                                                           
 *  You should have received a copy of the GNU General Public License     
 *  along with Wedge; if not, write to the                                
 *   Free Software Foundation, Inc.,                                       
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             
 *  
 *******************************************************************************/
 #ifndef WEDGEBASE_H
#define WEDGEBASE_H

/** @defgroup Base Base of the %Wedge framework
 * 
 * This module implements functionality essential to %Wedge, though mathematically uninteresting. 
 * */

 /** @{ 
 * @file wedgebase.h
 * @brief Basic objects of %Wedge
 */
 
//standard library includes
#include <string>
#include <set>
#include <vector>
#include <list>
#include <algorithm>

#include <cassert>
#include <limits.h>
#include <stdexcept>

#include <sstream>
#include <fstream> 
#include <iostream>

#include  <memory>


//ginac
#include <ginac/ginac.h>
//other includes in wedge
#include "wexception.h"
#include "wedge/wedgeconfig.h"
#include "wedge/convenience/printcontext.h"

/** @brief Main namespace */
namespace Wedge {
using namespace GiNaC;
using namespace std;
using GiNaC::function;

typedef int ZeroBased; ///< Type used for zero-based indices (e.g. C-style arrays)
typedef int OneBased; ///< Type used for one-based indices (e.g. maple-style arrays)

/** @brief This class implements one-based index semantics for vectors. 
 * 
 * As a rule of thumb, in %Wedge this class is  used for those vectors in which the index has some mathematical meaning (for instance,
 * one usually wants a basis of a vector space to have a one-based index). An extensive use might (sligthly) impact performance.  
 *
 * @note  ExVector inherits the standard zero-based operator[] from std::vector. So if v is an ExVector, v[0] and v(1) are equivalent expressions.
 * This convention is shared by other unrelated classes, e.g. LinearCombinations.
 */
struct ExVector : public exvector {
	ex& operator() (OneBased i) {return operator[](i-1);} ///<Return the i-th element of this vector as an lvalue
	ex operator() (OneBased i) const {return operator[](i-1);}///< Return the i-th element of this vector
	using exvector::exvector;
	using exvector::operator=;
	ExVector& operator=(ExVector&&)=default;
	ExVector& operator=(const ExVector&)=default;
	ExVector(const ExVector&)=default;
	ExVector(ExVector&&)=default;
};



/** @brief %Function object implementing a total ordering (needed by the STL)
 * 
 * @warning The ordinary < operator is not a total ordering for classes derived from basis. When using the STL,
 * the compare class basic_is_less should be passed as a template argument when a total ordering is needed (e.g. for sort, set...) 
 */
struct basic_is_less: public std::binary_function<basic, basic, bool>
{
	bool operator() (const basic& lh, const basic& rh) const { return lh.compare(rh) < 0; }		
};


////////////////////////////////////////////////////////////////////////////////
// 								Implementation						  	
////////////////////////////////////////////////////////////////////////////////


/** @brief Helper template functions to add a range of elements to a container
 * 
 * This function is a workaround for the fact that different container classes have different interfaces, consider e.g. std::set and GiNaC::lst in contrast to  std::vector and std::list
 */
template<typename Container, typename Iterator> Container& Insert(Container& container, Iterator begin, Iterator end)
{
	container.insert(container.end(),begin,end);
	return container;
}

/** @overload
 */
template<class Key, class Compare, class Alloc, typename Iterator> set<Key,Compare,Alloc>& Insert(set<Key,Compare,Alloc>& container, Iterator begin, Iterator end)
{
	container.insert(begin,end);
	return container;
}

/** @overload
 */
template<typename Iterator> lst& Insert(lst& container, ex x, Iterator begin,Iterator end)
{
	while (begin!=end) container.append(*begin++);
	return container;
}

/** @brief Helper template functions to add an element to a container
 * 
 * This function is a workaround for the fact that different container classes have different interfaces, consider e.g. std::set and GiNaC::lst in contrast to  std::vector and std::list
 */
template<typename Container> Container& Insert(Container& container, typename Container::value_type obj)
{
	container.push_back(obj);
	return container;
}

/** @overload
 */
template<class Key, class Compare, class Alloc> set<Key,Compare,Alloc>& Insert(set<Key,Compare,Alloc>& container, Key object )
{
	container.insert(object);
	return container;
}

/** @overload
 */
inline lst& Insert(lst& container, ex x)
{
	container.append(x);
	return container;
}



/** @deprecated This is made obsolete by Ginac 1.4.0
 */
void SplitIntoRealAndImaginary(ex,ex& realpart,ex& imaginarypart);

/** @deprecated This is made obsolete by Ginac 1.4.0
 */
inline ex RealPart(ex e)
{
	ex re,im;
	SplitIntoRealAndImaginary(e,re,im);
	return re;
}

/** @deprecated This is made obsolete by Ginac 1.4.0
 */
inline ex ImaginaryPart(ex e)
{
	ex re,im;
	SplitIntoRealAndImaginary(e,re,im);
	return im;
}






} /** @} */

#endif /*WEDGEBASE_H*/
