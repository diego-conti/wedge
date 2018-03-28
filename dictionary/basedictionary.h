/***************************************************************************
 *   Copyright (C) 2008 by Diego Conti, diego.conti@unimib.it              *
 *                                                                         *
 *   This file is part of Dictionary.                                      *
 *                                                                         *
 *   Dictionary is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 3 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   Dictionary is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with Dictionary; if not, write to the                           *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#ifndef BASEDICTIONARY_H_
#define BASEDICTIONARY_H_
/** @file basedictionary.h
 * @brief Abstract base class for class Dictionary and d operator semantics
 */

#include <wedge/manifold.h>
#include <wedge/liegroup.h>
#include <wedge/function.h>
#include <wedge/logging.h>
#include <wedge/leibniz.h>
#include <wedge/lambda.h>
#include "invariantforms.h"
#include "vvaluedforms.h"

using namespace Wedge;

class BaseDictionary;

class dOperator : public Leibniz<CompositeElement, Function>
{
	const BaseDictionary& d;
public:
	dOperator(const BaseDictionary& m) : d(m) {}
	void visit(const ncmul& alpha) {	
		Result()=0;
		bool sign=true; 	//true stands for +
		exvector v(alpha.begin(),alpha.end());
		
		for (unsigned i=0;i<alpha.nops();i++) {
			v[i]=this->RecursiveVisit(v[i]);
			if (sign) 
				this->Result()+=ncmul(v);
			else 
				this->Result()-=ncmul(v);
			v[i]=alpha.op(i);
			sign=IsOdd<CompositeElement>(v[i]) xor sign;
		}
	}
	void visit(const Wedge::Function& f); 
	void visit(const CompositeElement& alpha);
};

class BaseDictionary {
	exmap table;
	mutable dOperator thedOperator;
public:
	BaseDictionary() : thedOperator(*this) {}
	virtual ex d(const CompositeElement& alpha) const=0;
	virtual ex d(const Function& f) const
	{
		exmap::const_iterator it=table.find(f);
		if (it!=table.end()) return it->second;
		else {
			throw Manifold::dException(__FILE__,__LINE__,f);
			//ex dr=1/(4*RadiusSquared())*d(RadiusSquared()); 
			//return LieDerivative(dr,f)*d(RadiusSquared());
		}
	}
/** @brief Take d of a form expressed as a CompositeElement
 * @param alpha A CompositeElement with functions of r as coefficients
 */
	ex d(ex alpha) const {
		alpha.expand().accept(thedOperator);
		return thedOperator.GetResult().expand();
	}

	void Declare_d(ex f, ex df)
	{
		assert(is_a<Function>(f));
		table[f]=df;
	}
/** @brief Returns a symbolic function that represents the radial coordinate, may equal   \f$r=\sqrt{aa}\f$
 */
	virtual ex RadialCoordinate()  const=0;


	/** @brief Evaluate a form or function at the principal point, meaning a fixed point in the sphere bundle
	 * @param form A DifferentialForm with functions of the a's as coefficients
	 * @return form A DifferentialForm whose coefficients do not depend on either the a's or r
	 */
	virtual ex AtPrincipalPoint(ex form) const=0;

	/** @brief Evaluate a form or function at a principal point, meaning a fixed point in the sphere bundle \f$\{aa=r^2\}\f$
	 * @param form A DifferentialForm with functions of the a's as coefficients
	 * @param r The radius of the sphere bundle
	 * @return form A DifferentialForm whose coefficients do not depend on either the a's or r
	 */
	virtual ex AtPrincipalPoint(ex form, ex r) const=0;

	/** @brief Evaluate a form or function at the generic principal point, meaning a point on the line (slice).
	 * @param form A DifferentialForm with functions of r and the a's as coefficients
	 * @return form A DifferentialForm with functions of r as coefficients
	 */
	virtual ex AtGenericPoint(ex form) const=0;

	/** @brief Returns a pointer to a manifold representing the total space of the principal bundle
	 */
	virtual const Manifold* TotalSpace() const=0;
};

namespace Wedge {
namespace internal {
template<typename DegreeType> class ComputeDegree<CompositeElement,DegreeType > : public RecursiveVisitor<DegreeType>,public CompositeElement::visitor,public ncmul::visitor,public mul::visitor,public add::visitor,public basic::visitor
{
	void visit(const mul & x) {
		this->Result()=0;
		for (unsigned i=0;i<x.nops();i++)
			this->Result()+=this->RecursiveVisit(x.op(i));
	}
	void visit(const ncmul & x) {
		this->Result()=0;
		for (unsigned i=0;i<x.nops();i++)
			this->Result()+=this->RecursiveVisit(x.op(i));
	}
	void visit(const CompositeElement& x) {
		this->Result()=x.degree();
	}

	void visit(const add & x) {
		unsigned i=0;
		this->Result()=this->RecursiveVisit(x.op(i++));
		while (i<x.nops())
			if (this->RecursiveVisit(x.op(i++))!=this->Result()) throw InhomogeneousExpression(__FILE__,__LINE__);
	}
	void visit(const basic&) {
		this->Result()=0;
	}
};
}}


#endif
