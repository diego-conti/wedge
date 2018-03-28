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

 #ifndef LEIBNIZ_H
#define LEIBNIZ_H
/** @ingroup LinearAlgebra */
/** @{ 
 * @file leibniz.h
 * @brief Operators satisfying the Leibniz rule
 */
 
#include <ginac/ginac.h>
#include "utilities.h"
#include "lambda.h"
#include "linear.h"

namespace Wedge {
using namespace GiNaC;


/** @brief Linear operator over an \f$R\f$-module \f$R\otimes_{\mathbb{R}}V\f$ satisfying the %Leibniz rule
 * @param V The type of a vector space
 * @param R A type derived from symbol, representing the ring of coefficients
 */  

template<typename V, typename R> 
class Leibniz : 
	public function::visitor, public R::visitor,
	public AdditiveOperator<V>, public mul::visitor, public power::visitor, public ncmul::visitor
{
public:
	enum {ExpectsExpandedInput=1};
protected:
	void visit(const GiNaC::function& alpha) {		
		vector<R> variables;
		GetSymbols<R>(variables,alpha);
		this->Result()=0;
		for (typename vector<R>::const_iterator i=variables.begin();i!=variables.end();i++)
		{
			this->Result()+=this->RecursiveVisit(*i)*alpha.diff(*i);
		}
	}	
	void visit(const mul& alpha)
	{
		//a mul represents a commutative product x_1 ... x_n, so we implement the Leibniz rule by
		//(Dx_1)x_2...x_n + ... + (Dx_n)x_1...x_{n-1}
		this->Result()=0;
		for (unsigned i=0;i<alpha.nops();++i) {
			exvector v; v.reserve(alpha.nops());
			int j=0;
			v.push_back(this->RecursiveVisit(alpha.op(i)));
			while (j<i)
				v.push_back(alpha.op(j++));
			++j;
			while (j<alpha.nops())
				v.push_back(alpha.op(j++));
			this->Result()+=ncmul(v);
		}
	}
	void visit(const power& alpha)
	{
		ex basis=alpha.op(0);
		ex exponent=alpha.op(1);		
		ex dexponent=this->RecursiveVisit(exponent).expand();
		this->Result()=0;
		if (!dexponent.is_zero())
		{
			assert(is_a<numeric>(basis));
			ex logbasis=log(ex_to<numeric>(basis));
			this->Result()=alpha* dexponent*logbasis;
		}
		ex dbasis=this->RecursiveVisit(basis);			
		this->Result()+=power(basis,exponent-1) * exponent * dbasis;
		
	}	
	void visit(const ncmul& alpha)
	{
		LOG_WARN(alpha);
		throw WedgeException<std::invalid_argument>("Cannot apply Leibniz rule to unexpanded noncommutative product",__FILE__,__LINE__);
	}
};

namespace internal {

/** @internal
@brief Template specialization \see DerivationOver
*/

template<typename LeftType,typename W, typename R, typename Bilinear> 
class MyRightOperator<LeftType,Leibniz<W,R>,Bilinear> :
 		public Leibniz<W,R>
{
	const Bilinear* bil;
	const LeftType& pv;
	typedef Leibniz<W,R> RightOperator;
public:
	MyRightOperator(const LeftType& v,const Bilinear* bil) 	: pv(v)
	{
		this->bil=bil; 
	} 
protected:
	void visit(const W& w) {
		this->Result()=bil->Apply(pv,w);
	}
	void visit(const R& w) {
		this->Result()=bil->Apply(pv,w);
	}
};

}




} /** @} */
#endif
