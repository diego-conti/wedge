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

#include "differentialform.h"
#include "ginac/ginac.h"
#include "ginac/container.h"
#include "ginac/archive.h"
#include <iostream>
#include "utilities.h"
#include <sstream>
#include "wedgealgebraic.h"
#include "normalform.h"

namespace Wedge {
 using namespace  GiNaC;

VectorField::VectorField() :  RegClass(Name("Form","F")(ID)) {}

VectorField::VectorField(const Name& name) : RegClass(name) {}

namespace internal {

class HookOperator : public IBilinearOperator<AssociativeOperator<DifferentialForm>,Derivation<DifferentialForm>  > {
public:
	ex Apply(const VectorField& left,const VectorField& right) const;
};

ex HookOperator::Apply(const VectorField& left,const VectorField& right) const
{
	//Notice that I cannot use basic::compare here, because left and right might have different tinfo's
	//e.g. if left is a VectorField and right is a DifferentialOneForm.
	return left==right? 1 : 0;
}

}


ex Hook1(ex left, ex right)
{
	static internal::HookOperator hook;
	ex result=internal::HookOperator::BilinearOperator(left.expand(),right.expand(),&hook).expand();
		//adjust sign so that e^{12..n}\hook e^{12..n}=1
	return (Degree<DifferentialForm>(left) % 4<2)? result : -result;
}	


ex TrivialPairing(const VectorNormalForm& x, const VectorNormalForm& y)
{
	ex result;
	auto i=x.begin(), j=y.begin();
	while (i!=x.end() && j!=y.end()) {
		int cmp=i->first.compare(j->first);
		if (cmp>0) ++j;
		else if (cmp<0) ++i;
		else result+=(i++)->second * (j++)->second;
	}
	return result;
}

ex TrivialPairingImpl(ex left, ex right)
{
	LambdaVectorNormalForm x{left}, y{right};
	ex result=x.scalar_part*y.scalar_part;
	result+=TrivialPairing(x.vector_part,y.vector_part);
	result+=TrivialPairing(x.lambda_vector_part,y.lambda_vector_part);
	return result;
}


ex Hook(const LambdaVector& v, const Vector& w) {return 0;}

ex Hook(const Vector& v, const Vector& w) {
	ex result;
	if (v==w) result=1;
	return result;
}

ex Hook(const Vector& v, const LambdaVector& w)
{
	exvector result;
	int sign=1;
	auto j=w.begin();
	while(!(v==ex_to<Vector>(*j))) {	//why not !=?
		result.push_back(*j);
		sign=-sign;
		if (++j==w.end()) return 0;
	}
	while(++j!=w.end()) result.push_back(*j);
	return sign*LambdaVector{move(result)};
}


ex Hook(const LambdaVector& v, const LambdaVector& w)
{
	auto i=v.begin(), j=w.begin();
	exvector result;
	int sign=1;
	int v_parity=(v.nops()%2==0)? 1: -1;
	while(i!=v.end() && j!=w.end()) {
		if (*i==*j) {
			++i; 
			++j;
			v_parity=-v_parity;
		}
		else {
			sign*=v_parity;
			result.push_back(*j++);
		}
	}
	if (i==v.end()) {
		result.insert(result.end(),j,w.end());
		return sign*LambdaVector{move(result)};
	}
	else return 0;
}


ex Hook(ex v, ex w)
{
	LambdaVectorNormalForm left{v}, right{w};

	ex result=left.scalar_part*w;
	for (const auto& x : left.vector_part)
	for (const auto& y : right.vector_part) {
		ex hook=Hook(ex_to<Vector>(x.first),ex_to<Vector>(y.first));
		LOG_DEBUG(x.first<<"\\hook "<<y.first<<"="<< hook);
		if (!hook.is_zero()) result+=(x.second*y.second).expand()*hook;
	}
	for (const auto& x : left.vector_part)
	for (const auto& y : right.lambda_vector_part) {
		ex hook=Hook(ex_to<Vector>(x.first),ex_to<LambdaVector>(y.first));
		LOG_DEBUG(x.first<<"\\hook "<<y.first<<"="<< hook);
		if (!hook.is_zero()) result+=(x.second*y.second).expand()*hook;
	}
	for (const auto& x : left.lambda_vector_part)
	for (const auto& y : right.lambda_vector_part) {
		ex hook=Hook(ex_to<LambdaVector>(x.first),ex_to<LambdaVector>(y.first));
		LOG_DEBUG(x.first<<"\\hook "<<y.first<<"="<< hook);
		if (!hook.is_zero()) result+=(x.second*y.second).expand()*hook;
	}
	return result;
}
}
