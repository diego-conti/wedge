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
#ifndef NORMALFORM_H
#define NORMALFORM_H

#include "wedgebase.h"

namespace Wedge {

/** @brief The normal form of a linear combination of vectors (either of type Vector or LambdaVector)
*
* This class allows iterating through nonzero coefficients. The order is determined by ex_is_less, evaluated on vectors.
* Dereferencing the iterator gives a pair<ex,ex> where the first element is a VectorBase wrapped in an ex, and the second element a scalar expression
*
* test:is_a<VectorBase>(*i);
*/

class VectorNormalForm {
	friend class LambdaVectorNormalForm;
	exmap coefficients;
public:
	VectorNormalForm() {};
	VectorNormalForm(ex linear_combination);

	using const_iterator = exmap::const_iterator;
	const_iterator begin() const {return coefficients.begin();}
	const_iterator end() const {return coefficients.end();}
/** @brief Return a normal form in which coefficients are collected
*
* e.g. a*e(1)-a*e(2)+a*e(3) -> a*(e(1)-e(2)+e(3))
*/
	ex CollectCoefficients() const;
};


struct LambdaVectorNormalForm {
	ex scalar_part;
	VectorNormalForm vector_part, lambda_vector_part;
	LambdaVectorNormalForm(ex linear_combination);
	using const_iterator = exmap::const_iterator;
};

}
#endif
