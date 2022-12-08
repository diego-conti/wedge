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
#ifndef SIMPLIFIER_H
#define SIMPLIFIER_H

#include "wedge/base/expressions.h"
#include "ginac/matrix.h"
namespace Wedge {

struct Simplifier {
/** @brief Simplify a scalar expression
*
* Override in a derived class to define a custom simplifier*/
	virtual ex SimplifyScalar(ex e) const {return e.expand().normal();}
/** @brief Simplify a vector of expressions, either scalar or not
*
* Can be overridden in derived class to improve performance*/
	virtual ExVector Simplify(const exvector& e) const {
		ExVector result; result.reserve(e.size());
		for (exvector::const_iterator i=e.begin();i!=e.end();++i)
			result.push_back(Simplify(*i));
		return result;
	}

/** @brief Simplify a scalar equation
*
* TODO what about vector equations?
*/
	virtual ex ScalarEquationSimplify(ex e) const {return SimplifyScalar(e).numer();}

/** @brief simplify a generic expression 
*
* TODO is lst allowed?
*/
	ex Simplify(ex v) const; 

	set<ex,ex_is_less> Simplify(const set<ex,ex_is_less>& v) const {
		exvector x{v.begin(),v.end()};
		ExVector y=Simplify(exvector(x.begin(),x.end()));
		return set<ex,ex_is_less>{y.begin(),y.end()};
	}
	matrix Simplify(const matrix& m) const {
		exvector as_vector; as_vector.reserve(m.nops());
		for (int i=0;i<m.nops();++i) as_vector.push_back(m.op(i));
		auto simplified_vector=Simplify(as_vector);
		return matrix{m.rows(),m.cols(),lst{simplified_vector.begin(),simplified_vector.end()}};
	}

};


extern Simplifier default_simplifier;

//TODO test and document
struct TrigonometricSimplifier : public Simplifier {
	ex SimplifyScalar(ex e) const override;
};

}
#endif

