/*******************************************************************************
 *  Copyright (C) 2007-2023 by Diego Conti, diego.conti@unipi.it 
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
 #ifndef LINEARCOMBINATIONS_H_
#define LINEARCOMBINATIONS_H_

/** @ingroup LinearAlgebra */ 

/** @{ 
 * @file linearcombinations.h
 * @brief Containers of linear combinations with automatic reduction
 */
 
#include "wedge/base/wedgealgebraic.h"
#include "wedge/base/expressions.h"
#include "wedge/linearalgebra/bilinear.h"
#include "wedge/linearalgebra/ginaclinalg.h"

namespace Wedge {
using namespace GiNaC;

/** @brief Abstract base class for containers of linear combinations that perform automatic reduction on element insertion
 *
 * - A LinearCombinations object is constructed from a list of degree one polynomials, say \f$P_i=\sum_j a_{ij}v_j+c_i\f$, where the \f$v_j\f$ are simple elements of some type V and \f$a_{ij},c_i\f$ are coefficients. 
* - The polynomials are stored in expanded form, meaning that ex::expand()) is invoked on each element.
* - Elements can be subsequently added, STL-style, but not modified or removed
* - LinearCombinations provides the framework for automatic reduction, meaning that when elements are inserted, a minimal set of linearly independent elements in \f$\mathbb{R}[V]\f$ is chosen so as to preserve the filtration \f$\{ Span\{P_1,\dots,P_n\}, 1\leq n\leq size()\}\f$.
* - The actual algorithms that perform the reduction are implemented in the derived class templates Basis and AffineBasis
*/
class LinearCombinations {
public:
/** @brief Construct an empty basis
 */
	LinearCombinations() {size_=0;}
	
	virtual ~LinearCombinations() {}
	
/** @brief Construct a basis from an exvector of linearly independent elements
 *  @param basis a vector of linearly independent elements
 * 
 *  For efficiency, no check is made that the elements are linearly independent
 */
	LinearCombinations(const exvector& basis) 
	{
		e.reserve(basis.size());
		for (exvector::const_iterator i=basis.begin();i!=basis.end();++i)
			e.push_back(i->expand());
		size_=e.size();
	}

/** @brief Construct a basis from a range of linearly independent elements
 *  @param [from,to) A range of elements to add   
 * 
 *   For efficiency, no check is made that the elements are linearly independent
 */
	template<typename InputIterator>  LinearCombinations(InputIterator from, InputIterator to)
	{
		while (from!=to)
			e.push_back(from++->expand());
		size_=e.size();
	}

/** @brief Construct basis from an exvector of linearly independent elements
 *  @param basis A vector of linearly independent elements
 * 
 *  For efficiency, no check is made that the elements are linearly independent
 */
	void operator=(const exvector& basis) {
		e.clear();
		e.reserve(basis.size());
		for (exvector::const_iterator i=basis.begin();i!=basis.end();++i)
			e.push_back(i->expand());
		size_=e.size();
		MarkNotInitialized();		
	}

//ExVector semantics
	typedef exvector::value_type value_type; 			///<Implements container semantics, in analogy with std::vector
	typedef exvector::iterator iterator;				///<Implements container semantics, in analogy with std::vector
	typedef exvector::const_iterator const_iterator; 		///<Implements container semantics, in analogy with std::vector
	operator ExVector() const {return ExVector(begin(),end());}	///<Convert to type ExVector FIXME should be explicit but this breaks a number of things

	inline const_iterator begin() const {return e.begin();} 	///<Implements container semantics, in analogy with std::vector 
	inline const_iterator end() const {return e.begin()+size();} ///<Implements container semantics, in analogy with std::vector

/** @brief Implements container semantics, in analogy with std::vector

@warning This member is meant to be used in conjunction with LinearCombinations::insert.  Behaviour is undefined if non-const iterators
are used to modify the elements of LinearCombinations directly.
*/
	iterator begin()  {return e.begin();}
/** @brief Implements container semantics, in analogy with std::vector

@warning This member is meant to be used in conjunction with LinearCombinations::insert.  Behaviour is undefined if non-const iterators
are used to modify the elements of LinearCombinations directly.
*/
	iterator end()  {return e.begin()+size();}
	int size() const {assert(size_>=0 && size_<=e.size()); return size_;} ///<Return the number of elements of this basis
	
/** @brief Return the n-th element of this basis
 * @param n A zero-based index
 * @return The n-th element of this basis
 */
	ex operator[] (ZeroBased n)const {return e[n];}	

/** @brief Return the n-th element of this basis
 * @param n A one-based index
 * @return The n-th element of this basis
 */
	ex operator() (OneBased n) const {return e(n);}


/** @brief Add a range of elements to this basis
 *  @param at The position at which the range is to be inserted (e.g. this->end() to append to the existing generators)
 *  @param [from,to) A range of elements to add
 * 
 *  The elements in the range are not required to be linearly independent. After adding the elements,
 *  this function clears all redundant elements, in a flag-preserving way. 
 */
 
	template<typename InputIterator> void insert(iterator at, InputIterator from, InputIterator to)
	{
		assert(begin()<=at && at<=end());
		if (end()!=e.end())
			e.erase(end(),e.end()); //remove all the ''extra'' elements
		list<ex> expanded;
		while (from!=to)
			expanded.push_back(from++->expand());
		e.insert(at,expanded.begin(),expanded.end());
		size_=e.size();		
		Prune();	//remove redundant generators from basis
		MarkNotInitialized();
	}

/** @brief Remove all elements from this basis
 */
	void clear() {
		e.clear();
		size_=0;
		MarkNotInitialized();
	}

/** @brief Append a vector at the end of the basis 
 *  @param v A vector
 * 
 * If v lies in LinearCombinations, nothing happens.
 */
	virtual void push_back(ex v)
	{
		v=v.expand();
	    	insert(end(),&v,(&v)+1);
	}
 /** @brief Test whether the basis is empty
  *  @return true if the basis contains no generators 
 */	
	bool empty() const {
		return size_==0;	
	}

	

/** @brief Comparison operator, implemented to make mathematical sense
 * 
 * Two bases are equal if they consist of the same elements, arranged in the same order. The enlarged bases
 * are not compared. 
*/
	bool operator==(const LinearCombinations& o) const
	{
		return size()==o.size() && equal(begin(),end(),e.begin());
	}
protected:
	mutable ExVector e;		///<Enlarged basis. The first _size elements are the generators, and they are not modified by const members.
		
	int size_;	///<Size of this basis. 

/** @brief Mark this basis as not initialized
 * 
 *  Subclasses should call this function after modifying the elements of this->e directly. 
 */
	virtual void MarkNotInitialized() {}

/** @brief eliminates redundant generators, preserving order (\sa Basis).
 * 
 * Like do_const_prune(), but also checks that the first _size elements are independent.
*/	
	void CheckIndependenceAndPrune() const
	{
		if (ConstPrune()!=size_) throw WedgeException<std::logic_error>("A const element altered the basis, or Basis erroneously marked as initialized",__FILE__,__LINE__);	
	}
	
private:	
/** @brief eliminates redundant generators, preserving order (\sa Basis).
 * 
 * Like ConstPrune(), except that prune() updates the _size member in such a way that the span of the first _size
 * elements stays the same
*/
	void Prune() {
		size_=ConstPrune(); 
	}
/** @brief eliminates redundant generators, preserving order (\sa Basis).
 * 
 * The k-th element of the resulting basis will be the first generator \f$e_i\f$ such that \f$\langle e_1,\dots,e_i\rangle\f$ has dimension k.
 * Returns the dimension of the vector space spanned by the first _size elements
*/
	virtual int ConstPrune() const =0;

};

} /** @} */
#endif


