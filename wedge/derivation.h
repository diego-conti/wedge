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

 #ifndef WEDGEDERIVATION_H
#define WEDGEDERIVATION_H
/** @ingroup LinearAlgebra */
/** @{ 
 * @file derivation.h
 * @brief Operators on algebras; derivations
 */
 
#include <ginac/ginac.h>
#include "utilities.h"
#include "lambda.h"
#include "linear.h"
#include "leibniz.h"

namespace Wedge {
using namespace GiNaC;

template<typename T> class AssociativeOperator;


/** @brief Abstract base class for a linear associative operator on the exterior algebra over V
 * @param T The type Lambda<V> of decomposable elements in the exterior algebra
 *
 * This class can be used to construct a left-associative bilinear operator, meaning that
 *  \f$ B(v_1\wedge\dots\wedge v_k, w)=B(v_1,B(\dots,B(v_k,b)))\f$
 *
 * Notice that in particular this implies that \f$B(v,B(v,w))=0\f$ for all \f$v,w\f$.
 */  
template<typename V> class AssociativeOperator<Lambda<V> > : public LinearOperator<Lambda<V> > {
};

namespace internal {

/** @internal
@brief Template specialization \see AssociativeOperator
*/
  template<typename V,typename RightOperator, typename Bilinear> 
class MyLeftOperator <AssociativeOperator<Lambda<V> >,RightOperator, Bilinear > :
	public LinearOperator<Lambda<V> > 
{
	typedef LinearOperator<Lambda<V> > LeftOperator;	
	ex w;
	const Bilinear* bil;

public:
	MyLeftOperator(ex w,const Bilinear* bil) {this->w=w; this->bil=bil;}

protected:
	void visit(const V& v) {
		MyRightOperator<V,RightOperator,Bilinear> oper(v,bil);			
		w.accept(oper);
		this->Result()=oper.GetResult();
	}
	void visit(const Lambda<V>& v) {
		this->Result()=w;
		for (ncmul::const_reverse_iterator i=v.rbegin();i!=v.rend();i++)
		{
			MyRightOperator<V,RightOperator,Bilinear > oper(ex_to<V>(*i),bil);			
			this->Result().accept(oper);
			this->Result()=oper.GetResult();
		}		
	}	
};

}


template<typename T,  bool SKEW=true> class Derivation;

/** @brief Abstract base class for a linear derivation on the exterior algebra Lambda<V> over the field of real numbers
 * @param T The type Lambda<V> of decomposable elements in the exterior algebra
 * @param SKEW Determines whether the derivation is skew-symmetric or symmetric  
 */  

  template<typename T, bool SKEW> 
class Derivation <Lambda<T>, SKEW > : 
	public LinearOperator<Lambda<T> >, public ncmul::visitor
{
public:
	enum {ExpectsExpandedInput=1};
	typedef T OperatesOn;
protected:
	void visit(const Lambda<T>& alpha)
	{
		this->Result()=0;
		bool sign=true; 	//true stands for +
		exvector v(alpha.begin(),alpha.end());
		
		for (unsigned i=0;i<alpha.nops();i++) {			
			v[i]=this->RecursiveVisit(v[i]);
			if (sign) 
				this->Result()+=ncmul(v);
			else 
				this->Result()-=ncmul(v);
			v[i]=alpha.op(i);
			sign=(not SKEW) or (IsOdd<Lambda<T> >(v[i]) xor sign);
		}
	}
	void visit(const ncmul& alpha)
	{
		LOG_WARN(alpha);
		throw WedgeException<std::invalid_argument>("Cannot take derivation of unexpanded noncommutative product",__FILE__,__LINE__);
	}

};



/** @brief Abstract base class for a linear derivation on the exterior algebra Lambda<V> over the ring R
 * @param T The type Lambda<V> of decomposable elements in the exterior algebra
 * @param R The type of simple elements in the ring R
 * @param SKEW Determines whether the derivation is skew-symmetric or symmetric
 *
 * Expects an expanded input
 *
 * @remark The derivation must be linear over \f$\mathbb{R}\f$, but not necessarily over \f$R\f$.
 * @warning R must derive from symbol, because GiNaC::ex::diff requires so
 */  

  template<typename LambdaT, typename R,bool SKEW=true> 
class DerivationOver : 
	public Leibniz<LambdaT,R>
{
protected:
	void visit(const LambdaT& alpha)
	{
		this->Result()=0;
		bool sign=true; 	//true stands for +
		exvector v(alpha.begin(),alpha.end());
		
		for (unsigned i=0;i<alpha.nops();i++) {			
			v[i]=this->RecursiveVisit(v[i]);
			if (sign) 
				this->Result()+=ncmul(v);
			else 
				this->Result()-=ncmul(v);
			v[i]=alpha.op(i);
			sign=(not SKEW) or (IsOdd<LambdaT >(v[i]) xor sign);
		}
	}
}
;

namespace internal {

/** @internal
@brief Template specialization \see DerivationOver
*/
  template<typename LeftType,typename W,typename R,typename Bilinear, bool SKEW> 
class MyRightOperator<LeftType,DerivationOver<Lambda<W>,R,SKEW>,Bilinear> :
 		public DerivationOver<Lambda<W>,R,SKEW >	
{
	const Bilinear* bil;
	const LeftType& pv;
	typedef DerivationOver<Lambda<W>,R,SKEW > RightOperator;
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
