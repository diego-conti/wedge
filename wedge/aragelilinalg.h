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
#ifndef ARAGELILINEARTRAITS_H_
#define ARAGELILINEARTRAITS_H_

/** @ingroup ExternalAlgorithms */ 

/** @{ 
 * @file aragelilinalg.h
 * @brief Alternative implementation of linear algebra algorithms, using %Arageli
 * 
 * Tested with Arageli 2.2.1.191
 *
 * @note The default, %GiNaC-based algorithms appear to be faster
 */
#include "anylinalg.h"
#include "ginaclinalg.h"
#include "wexception.h"
#include "wedgealgebraic.h"
#include "linear.h"


namespace Arageli {
	GiNaC::ex abs(GiNaC::ex a)
	{
		return GiNaC::abs(a);	
	}
	bool is_null(GiNaC::ex a)
	{
		return a.expand().is_zero();
	}
}

#include <arageli/arageli.hpp> 

namespace Arageli {	
	template<> struct type_traits< GiNaC::ex > : public type_traits_default<GiNaC::ex>
	{
		static const bool is_specialized = true;

		/// True iff T is fraction (but not necessary a rational number, for last see below).
		static const bool is_rational = true;
		
		static const bool is_number = false;
		static const bool is_integer_number = false;
		static const bool is_polynom = false;
		static const bool is_real_number = false;
		static const bool is_rational_number = false;
		static const bool is_complex_number = false;
		static const bool is_ring = true;	// +, -, *, null, unit
		static const bool is_field = true;	// +, -, *, /, null, unit
		static const bool is_finite = false;
		static const bool is_additive_group = true;	// +, -, null
		static const bool is_multiplicative_group =	true;	// *, /, unit
		static const bool is_integer_modulo_ring = false;
		static const bool is_matrix = false;
		static const bool is_vector = false;
		
		static const bool has_zero_divisor = false;
		static const bool has_commutative_multiplication = true;
		static const bool has_commutative_addition = true;;
		static const bool has_associative_multiplication = true;
		static const bool has_associative_addition = true;
		static const bool has_distributive_multiplication = true;
		static const bool has_distributive_addition = true;
		static const bool has_null = true;
		static const bool has_unit = true;
		static const bool has_opposite_unit = true;
		static const bool is_entire_ring = true;
	
		/// True iff type is composite type consists another elements.
		static const bool is_aggregate = false;
	
		/// Type of each element if T is composite type.
		typedef GiNaC::ex element_type;
	
		template <typename T1, bool REFCNT2>
		struct other_element_type_refcnt;
	
		typedef type_category::type category_type;
		static const category_type category_value;
			
	};	
}

namespace Wedge {
using namespace GiNaC;


/** @brief Alternative implementation of linear algebra algorithms (uses %Arageli)
 */
struct ArageliLinAlgAlgorithms {
	/** @brief A class to compute a minimal set of independent rows in a matrix
	 */	
	class IndependenceMatrix {
		Arageli::matrix<ex> m;
		vector<int> independentRows;		
	public:	
		typedef vector<int>::const_iterator const_iterator; ///< Iterator type for the container of the indipendent rows
	/** @brief Construct a zero matrix
	 *  @param r Number of rows
	 *  @param c Number of columns
	 */
		IndependenceMatrix(int r, int c) :m (c,r,Arageli::fromsize_t()) {}

	/** @brief Compute a minimal set of linearly independent rows
	 * 
	 * The set of independent rows can subsequently be retrieved as the range
	 * [IndependentRowsBegin(),IndependentRowsEnd())
	 */
		void ChooseLinearlyIndependentRows () {
			Arageli::matrix<ex> b,q;
			//reduce to echelon form, obtaining the set of independent columns						
			Arageli::rref(m,b,q,independentRows); 
		}
		
	/** @brief Return the entry (r,c) as an lvalue
	 *  @param r A zero-based index
	 *  @param c A zero-based index
	 */		
		ex& M(ZeroBased r, ZeroBased c)  {
			//since we want the set of independent rows, and Arageli::rref gives the columns,
			//we effectively transpose the matrix by swapping r and c
			return m.at(c,r);
		}
		
	/** @brief Return the entry (r,c)
	 *  @param r A zero-based index
	 *  @param c A zero-based index
	 */				
		ex M(ZeroBased r, ZeroBased c) const {return m.at(c,r);}	
		
	/** @brief The range  [IndependentRowsBegin(),IndependentRowsEnd()) is a list of integers, each 
	 * representing a row in the matrix.
	 *  @return The beginning of the range
	 */
	 	const_iterator IndependentRowsBegin() const {return independentRows.begin();}
	/** @brief The range  [IndependentRowsBegin(),IndependentRowsEnd()) is a list of integers, each 
	 * representing a row in the matrix.
	 * @return The end of the range
	 */
		const_iterator IndependentRowsEnd() const {return independentRows.end();}
	
	/** @brief The number of rows
	 */
		inline int rows() const {return m.nrows();}
	/** @brief The number of columns
	 */
		inline int cols() const {return m.ncols();}				
	};
	
	/** @brief A class to compute the inverse of a matrix
	 * 
	 *  The interface is the same as that of GiNaC::matrix \sa GiNaC::matrix 
	 */
	class InverseMatrix {
		Arageli::matrix<ex> m;		
	public:
		InverseMatrix(int r, int c) : m(r,c,Arageli::fromsize_t()) 
		{
			assert(r==c);
		}
	
	/** @brief Return the entry (r,c) as an lvalue
	 *  @param r A zero-based index
	 *  @param c A zero-based index
	 */		
		ex& operator()(ZeroBased r, ZeroBased c)  {return m.at(r,c);}
	
	/** @brief Return the entry (r,c)
	 *  @param r A zero-based index
	 *  @param c A zero-based index
	 */		
		ex operator()(ZeroBased r, ZeroBased c) const {return m.at(r,c);}
		
		/** @brief Returns the inverse matrix
		 * 
		 * @warning This implementation alters *this
		 * */
		matrix inverse()
		{
			m.inverse();
			exvector v(m.begin(),m.end());
			assert(m.nrows()==m.ncols());
			assert(v.size()==m.nrows()*m.ncols());
			return matrix(m.nrows(),m.ncols(),v);
		}		
		
		int cols() const {return m.ncols();} ///< The number of columns
		int rows() const {return m.nrows();} ///< The number of rows
	}; 
	
/** @brief Solve a linear system of equations
* 
* @sa GiNaC::lsolve()
*/ 
	inline static lst lsolve(lst eqns,lst unknowns)
	{ 
		return GinacLinAlgAlgorithms::lsolve(eqns,unknowns);
	}
};


/** @brief Overloaded output operator for matrices
 *  @param os An output stream
 *  @param m A matrix
 *  @return A reference to os
 */

template<typename ostream>
ostream& operator<<(ostream& os,  const ArageliLinAlgAlgorithms::IndependenceMatrix& m)
{
	return internal::Output<ArageliLinAlgAlgorithms>(os,m);
}

/** @brief Overloaded output operator for matrices
 *  @param os An output stream
 *  @param m A matrix
 *  @return A reference to os
 */

template<typename ostream>
ostream& operator<<(ostream& os,  const ArageliLinAlgAlgorithms::InverseMatrix& m)
{
	return internal::Output<ArageliLinAlgAlgorithms>(os,m);
}

} /** @} */
#endif /*ARAGELILINEARTRAITS_H_*/
