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
 #ifndef GINACLINEARTRAITS_H_
#define GINACLINEARTRAITS_H_

/** @ingroup ExternalAlgorithms  */

/** @{ 
 * @file ginaclinalg.h
 * @brief Standard implementation of linear algebra algorithms, using %GiNaC
 */

#include "wedge/base/wedgealgebraic.h"
#include "wedge/base/expressions.h"
#include "wedge/linearalgebra/anylinalg.h"
#include "wedge/linearalgebra/linear.h"

namespace Wedge {
using namespace GiNaC;

/** @brief Default implementation of linear algebra algorithms (uses %GiNaC)
 */
struct GinacLinAlgAlgorithms {
	/** @brief A class to compute a minimal set of independent rows in a matrix
	 */	
	class IndependenceMatrix  {
	public:
		typedef list<int>::const_iterator const_iterator; ///< Iterator type for the container of the indipendent rows
	
	/** @brief Construct a zero matrix
	 *  @param r Number of rows
	 *  @param c Number of columns
	 */	
		IndependenceMatrix(int r, int c) : m(r*c,ex(0)),col(c) {
			n_rows=r; n_cols=c;
			for (int i=0;i<c;i++) col[i]=i;	
		}
		
	/** @brief Return the entry (r,c) as an lvalue
	 *  @param r A zero-based index
	 *  @param c A zero-based index
	 */		
		inline ex& M(ZeroBased r, ZeroBased c)
		{
			return m[r*n_cols+col[c]]; 
		}
		
	/** @brief Return the entry (r,c)
	 *  @param r A zero-based index
	 *  @param c A zero-based index
	 */				
		inline ex M(ZeroBased r, ZeroBased c) const
		{
			return m[r*n_cols+col[c]]; 
		}
	
	/** @brief Compute a minimal set of linearly independent rows
	 * 
	 * The set of independent rows can subsequently be retrieved as the range
	 * [IndependentRowsBegin(),IndependentRowsEnd())
	 */
		void ChooseLinearlyIndependentRows();
	/** @brief The range  [IndependentRowsBegin(),IndependentRowsEnd()) is a list of integers, each 
	 * representing a row in the matrix.
	 *  @return The beginning of the range
	 */	
		const_iterator IndependentRowsBegin() const {return independentRows.begin();}
	/** @brief The range  [IndependentRowsBegin(),IndependentRowsEnd()) is a list of integers, each 
	 * representing a row in the matrix.
	 *  @return The end of the range
	 */	
		const_iterator IndependentRowsEnd() const {return independentRows.end();}
		
		inline int rows() const {return n_rows;} ///< The number of rows
		inline int cols() const {return n_cols;} ///< The number of columns
	protected:
		bool pivot(unsigned row, unsigned col);
		inline void SwapColumns(int c1,int c2)
		{
			int swap=col[c1];
			col[c1]=col[c2];
			col[c2]=swap;
		}
		
		list<int> independentRows;
	private:
		exvector m;
		int n_rows,n_cols;
		vector<int> col;				
	};
		
	//typedef matrix InverseMatrix;	///< Use the class GiNaC::matrix to compute the inverse of a matrix  
	

/** @brief A class to compute the inverse of a matrix
 *
 * The code of this class is copied from GiNaC::matrix, with the difference that if the entries of the matrix are not rational functions, Gauss elimination
 * is used as opposed to Bareiss.
 */	
	class InverseMatrix : public matrix {
	public:
		InverseMatrix(const matrix& m) : matrix(m) {}
		using matrix::matrix;
		using matrix::operator=;
//		matrix_init<ex, exvector::iterator> operator=(const ex & x) {
//			return matrix::operator=(x);
//		}

		matrix inverse() {
			if (row != col)
			throw (std::logic_error("matrix::inverse(): matrix not square"));
	
			matrix identity(row,col);
			for (unsigned i=0; i<row; ++i)
				identity(i,i) = 1;
			matrix vars(row,col);
			for (unsigned r=0; r<row; ++r)
				for (unsigned c=0; c<col; ++c)
					vars(r,c) = symbol();
	
			matrix sol(row,col);
			try {
				exmap repl;
				int i=0;
				while (i<nops() && repl.empty()) 
					op(i++).to_rational(repl);
				LOG_DEBUG(op(i-1));
				LOG_DEBUG(repl);
				sol = this->solve(vars,identity,
					repl.empty() ?  solve_algo::automatic : solve_algo::gauss
				);
			} 
			catch (const std::runtime_error & e) {
				if (e.what()==std::string("matrix::solve(): inconsistent linear system"))
					throw (std::runtime_error("matrix::inverse(): singular matrix"));
				else
					throw;
			}
			return sol;		
		}	
	};
		
	inline static matrix MatrixInverse(const matrix& m) {
		return InverseMatrix{m}.inverse();
	}
	
/** @brief Solve a linear system of equations
 * 
 * A wrapper about GiNaC::lsolve(eqns,unknowns) (throwing more informative exceptions)
*/ 
	inline static lst lsolve(lst eqns,lst unknowns)
	{
		try {
			ex result=GiNaC::lsolve(eqns,unknowns);
			assert(is_a<lst>(result));
			return ex_to<lst>(result);
		}
// ginac's lsolve may throw either logic_error or invalid_argument which also derives from logic_error.
		catch (const std::logic_error& c)
		{
			LOG_ERROR(eqns);
			LOG_ERROR(unknowns);
			throw;
		}
		/*{ //paranoid check commented out
			ex result2=GiNaC::lsolve(eqns,unknowns,solve_algo::gauss);
			int free=0, free2=0;
			for (int i=0;i<result2.nops();++i)
				if (result2.op(i).lhs()==result2.op(i).rhs())
					++free2;
			for (int i=0;i<result.nops();++i)
				if (result.op(i).lhs()==result.op(i).rhs())
					++free;
			if(free!=free2)
			{
				LOG_ERROR(eqns);
				LOG_ERROR(result);
				LOG_ERROR(result2);
			}
			for (int i=0;i<eqns.nops();++i)
			{
				ex h=eqns.op(i);
				if (!(h.lhs()-h.rhs()).subs(result).normal().is_zero())
				{
					LOG_ERROR(eqns);
					LOG_ERROR(eqns.op(i));
					LOG_ERROR(result);
					LOG_ERROR(result2);
					break;
				}
			}
		}*/
	}	
};

typedef GinacLinAlgAlgorithms DefaultLinAlgAlgorithms;

/** @brief Overloaded output operator for matrices
 *  @param os An output stream
 *  @param m A matrix
 *  @return A reference to os
 */

inline ostream& operator<<(ostream& os,  const GinacLinAlgAlgorithms::IndependenceMatrix& m)
{
	return internal::Output<GinacLinAlgAlgorithms>(os,m);
}

/** @brief Overloaded output operator for matrices
 *  @param os An output stream
 *  @param m A matrix
 *  @return A reference to os
 */

//inline ostream& operator<<(ostream& os,  const GinacLinAlgAlgorithms::InverseMatrix& m)
//{
//	return internal::Output<GinacLinAlgAlgorithms>(os,m);
//}

} /** @} */
#endif /*GINACLINEARTRAITS_H_*/
