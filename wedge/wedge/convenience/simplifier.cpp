/***************************************************************************
 *   Copyright (C) 2015 by Diego Conti, diego.conti@unimib.it              *
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

#include "wedge/convenience/simplifier.h"
#include <ginac/ginac.h>
#include "wedge/base/normalform.h"
namespace Wedge {
using namespace GiNaC;

Simplifier default_simplifier;

//TODO use Simplify (exvector) to increase performance
ex Simplifier::Simplify(ex v) const {
	ex result;
	LambdaVectorNormalForm normal_form{v};
	for (auto coeff : normal_form.vector_part)
		result+=coeff.first * SimplifyScalar(coeff.second);
	for (auto coeff : normal_form.lambda_vector_part)
		result+=coeff.first * SimplifyScalar(coeff.second);
	result+=SimplifyScalar(normal_form.scalar_part);
	return result;
}

ex TrigonometricSimplifier::SimplifyScalar(ex e) const {
		e=e.expand(expand_options::expand_function_args).subs(cos(Pi/2-wild())==sin(wild())).subs(sin(Pi/2-wild())==cos(wild()));
		e=e.subs(cos(2*wild())==2*cos(wild())*cos(wild())-1).subs(sin(2*wild())==2*sin(wild())*cos(wild())).expand();

		exset matches;
		e.find(pow(sin(wild(1)),wild(2)),matches);
		lst subs;
		for (exset::const_iterator i=matches.begin();i!=matches.end();++i) {
			ex exponent=i->op(1);
			if (is_a<numeric>(exponent)) {
				const numeric& as_numeric = ex_to<numeric>(exponent);
				if (as_numeric.is_even())
					subs.append(pow(sin(wild()),as_numeric)==pow(1-pow(cos(wild()),2),as_numeric/2));
				else if (as_numeric.is_odd())
					subs.append(pow(sin(wild()),as_numeric)==sin(wild())*pow(1-pow(cos(wild()),2),(as_numeric-1)/2));
			}
		}
		return e.subs(subs).expand().normal();
	}

}
