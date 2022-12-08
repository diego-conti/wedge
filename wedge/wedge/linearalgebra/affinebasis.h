/*******************************************************************************
 *  Copyright (C) 2007-2022 by Diego Conti, diego.conti@unipi.it 
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
 #ifndef AFFINEBASIS_H_
#define AFFINEBASIS_H_

/** @ingroup LinearAlgebra */ 

/** @{ 
 * @file affinebasis.h
 * @brief Affine spaces as bases of independent degree one polynomials
 */
 
#include "linearcombinations.h"

namespace Wedge {
using namespace GiNaC;

/** @brief A class to represent a basis of generators for some subspace of the space of degree one polynomials in T
 *
 * This class is analogous to Basis<T>, with the following differences:
 * - Its elements are allowed to have a scalar term, e.g. \f$2v+1\f$ where \f$v\f$ has type T.
 * - There are no functions Components() and dual()
 *
 * @note The type T need not derive from Vector. For instance, one can describe an affine subspace of a VSpace<V> using an AffineBasis<VSpace<V>::Coordinate> object.
 */

template<typename T, typename LinAlgAlgorithms=DefaultLinAlgAlgorithms> class AffineBasis : public LinearCombinations, public LinAlgAlgorithms {
public:
/** @brief Construct a basis from a range of linearly independent elements
 *  @param [from,to) A range of elements to add   
 * 
 *   For efficiency, no check is made that the elements are linearly independent
 */
	template<typename InputIterator>  AffineBasis(InputIterator from, InputIterator to) : LinearCombinations(from,to)  {}
/** @overload
*/
	AffineBasis(const exvector& basis) : LinearCombinations(basis) {}
/** @brief Construct an empty basis
*/
	AffineBasis() {}
private:
	int ConstPrune() const {
		exvector symbols;
		GetSimple<T>(symbols,this->e.begin(),this->e.end());
		typename LinAlgAlgorithms::IndependenceMatrix m(this->e.size(),symbols.size()+1);
		for (int i=0;i<this->e.size();++i)
		{
			internal::NormalFormHelper<Poly1<T> > v;
			this->e[i].accept(v);
			int j;
			for (j=0;j<symbols.size();++j)
			{
				exmap::const_iterator k=v.coeffs.find(symbols[j]);
				if (k!=v.coeffs.end()) 
					m.M(i,j)=k->second;
			}
			exmap::const_iterator k=v.coeffs.find(1);
			if (k!=v.coeffs.end()) 
				m.M(i,j)=k->second;
		}
				
		m.ChooseLinearlyIndependentRows();
		int new_size=0;
		int j=0;
		for (typename LinAlgAlgorithms::IndependenceMatrix::const_iterator i=m.IndependentRowsBegin();i!=m.IndependentRowsEnd();i++,j++)
		{
			if (*i<this->size_) ++new_size;
			if (*i!=j) this->e[j]=this->e[*i];			
		}
		this->e.resize(j);
		return new_size;		
	}

};

} /** @} */
#endif

