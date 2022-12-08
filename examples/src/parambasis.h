/***************************************************************************
 *   Copyright (C) 2009 by Diego Conti, diego.conti@unipi.it              *
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

#ifndef PARAMBASIS_H
#define PARAMBASIS_H
 
#include <wedge/wedge.h>

namespace Wedge {
struct ContainsIf {
	enum Flag {Generically, NotGenerically, Always, Never} flag;
	ex moredata;
	ContainsIf(Flag f, ex data) : flag(f), moredata(data) {}
	ContainsIf(Flag f) : flag(f) {assert(f==Always || f==Never);}
};

ostream& operator<<(ostream& os, ContainsIf c)
{
	switch (c.flag) {
		case ContainsIf::Generically :
			os<<"The element belongs generically to the space, provided the following elements are well-defined"<<endl<<c.moredata<<endl; break;
		case ContainsIf::NotGenerically :
			os<<"The element does not generically belong to the space, because the following equations have generically no solution"<<endl<<c.moredata; break;
		case ContainsIf::Always :
			os<<"The element belongs to the space"<<endl; break;
		case ContainsIf::Never :
			os<<"The element does not belong to the space"<<endl; break;
		default:
			assert(false);
	}
	return os;
}



/** @brief Container representing a flag of vector spaces depending on parameters
 *
 * Analogous to Basis<T>, except that:
 * - basis elements depend on parameters.
 * - no automatic reduction is performed
 * Should the parameters have a type?
 */
template<typename T> class ParamBasis : public LinearCombinations {
public:
/////////////////////////////////////////////////////////////////////
// 		Constructors
/////////////////////////////////////////////////////////////////////
	ParamBasis() {}
	ParamBasis(const exvector& generators) : LinearCombinations(generators) {}
 	template<typename InputIterator> ParamBasis(InputIterator from,InputIterator to) : LinearCombinations(from,to) { }
	ContainsIf Contains(ex v) const {
		exvector symbols;
		GetSimple<T>(symbols,begin(),end());
		GetSimple<T>(symbols,size_);

		lst unknowns;
		for (int i=1;i<=size();++i)	
			unknowns.append(symbol("x"+ToString(i)));

		//let e_i = a_{ij} v_j, v = b_j v_j
		//we must set up the linear system a_{ij} x_i = b_j
		lst eqns;
		for (int j=0;j<symbols.size();++j)
		{
			ex lhs, rhs;
			for (int i=0;i<size();++i)
			{
				internal::NormalFormHelper<T> v;
				e[i].accept(v);	
				exmap::const_iterator k=v.coeffs.find(symbols[j]);
				if (k!=v.coeffs.end()) //then the component a_{ij} of e[i] along symbols[j] is k->second
					lhs+=k->second * unknowns[i];
			}
			internal::NormalFormHelper<T> v;
			exmap::const_iterator k=v.coeffs.find(symbols[j]);
			if (k!=v.coeffs.end()) //then the component b_j of v along symbols[j] is k->second
				rhs=k->second;
			eqns.append(lhs==rhs);
		}
	
		LOG_INFO(eqns);
		ex sol=	lsolve(eqns,unknowns);
		LOG_INFO(sol);
		if (sol==lst()) return ContainsIf(ContainsIf::NotGenerically, eqns);
		else return ContainsIf(ContainsIf::Generically, sol);
	}
private:
	int ConstPrune() const {return size_;}

};
}
#endif

