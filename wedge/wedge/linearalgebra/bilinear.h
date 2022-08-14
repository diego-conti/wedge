/*******************************************************************************
 *  Copyright (C) 2007-2022 by Diego Conti, diego.conti@unimib.it 
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
 #ifndef BILINEAR_H_
#define BILINEAR_H_

/** @ingroup LinearAlgebra */
/** @{ 
 * @file bilinear.h
 * @brief Bilinear operators on vector spaces
 */
#include "wedge/linearalgebra/linear.h"
#include "wedge/base/expressions.h"

namespace Wedge {
 using namespace  GiNaC;
 
// ****************************************************************************
// *							BilinearOperator						  *
// ****************************************************************************
namespace internal {
/** @brief Used internally by IBilinearOperator
 *
 * @warning Since version 0.2.3, a pointer to the constructor parameter v is stored internally. This means that the MyRightOperator object should not outlive the argument v.
 */
template<typename LeftType,typename RightOperator, typename Bilinear> class MyRightOperator : public RightOperator {
public:
	MyRightOperator(const LeftType& v,const Bilinear* bil) : pv(v)
	{
		this->bil=bil;
	} 
protected:
	void visit(const typename RightOperator::OperatesOn& w) {
		this->Result()=bil->Apply(pv,w);
	}
private:
	const Bilinear* bil;
	const LeftType& pv;
};

/** @brief Used internally by IBilinearOperator
 */
template<typename LeftOperator, typename RightOperator, typename Bilinear> class MyLeftOperator :
	public LeftOperator 
{
	typedef typename LeftOperator::OperatesOn LeftType;
public:
	MyLeftOperator(ex w,const Bilinear* bil) {this->w=w; this->bil=bil;}
	void visit(const LeftType& v) {
		MyRightOperator<LeftType,RightOperator,Bilinear> oper(v,bil);			
		w.accept(oper);
		this->Result()=oper.GetResult();
	}
private:
	ex w;
	const Bilinear* bil;	
};

}

/** @internal 
 * @brief Abstract base class for bilinear operators
 * @param LeftOperator A type describing the operator \f$ v\to B(v,w)\f$
 * @param RightOperator A type describing the operator \f$ w\to B(v,w)\f$
 * 
 * To define a bilinear operator, subclass from IBilinearOperator and define
 * 
 * 		ex Apply(Right v, Left w) const
 * 
 * in your subclass. This function may be virtual; it is not declared here because
 * the types of its arguments cannot be determined uniquely (consider for instance linear operators on Lambda<V>)
 *
 * @todo Consider changing the implementation. For instance, one could write first \f$v,w\f$ in terms of fixed bases of \f$V,W\f$, and then
 * invoke Apply. This would be faster in the case that the basis consists of simple elements, and could be generalized to use orthonormal bases, potentially fixing the CliffordDot bug and avoiding the non-optimal call to Components() in Connection::Nabla. The problem is fitting in the Leibniz, Derivation and Associative variants.
 */
template<typename LeftOperator, typename RightOperator> struct IBilinearOperator {	
	/** @brief Evaluate the bilinear operator
	 * @param v,w The operands
	 * @param _this The object that implements the bilinear operator 
	 * @param This A class derived from IBilinearOperator that defines ex Apply(Right v, Left w) const
	 */ 
	template<typename This> static inline ex BilinearOperator(ex v, ex w, const This* _this) {
		if (LeftOperator::ExpectsExpandedInput!=0) v=v.expand();
		if (RightOperator::ExpectsExpandedInput!=0) w=w.expand();
		internal::MyLeftOperator<LeftOperator,RightOperator,This> oper(w,_this);
		v.accept(oper);
		return oper.GetResult();
	}
};

// ****************************************************************************
// *							TrivialPairingOperator						  *
// ****************************************************************************
/** @brief Class implementing the trivial pairing, namely the scalar product for which simple (or decomposable) elements are orthogonal
 *  @param V The type of simple elements in the relevant vector space
 *
 *  @remark Some functions, like VSpace::Subspace(), accept a bilinear operator as an argument; this class can be used to that effect.
 * 
 *  @sa TrivialPairing
 *  @note To work with the exterior algebra over V, use TrivialPairingOperator<Lambda<V> >
 */
template<typename V> class TrivialPairingOperator : public IBilinearOperator<LinearOperator<V>,LinearOperator<V> > {
public:	
	ex Apply(const V& v, const V& w) const
	{
		//return v.compare_same_type(w)==0? 1 : 0;
		return v==w? 1:0;
	}
};

ex TrivialPairingImpl(ex v,ex w);


/** @brief The trivial pairing, namely the scalar product for which simple (or decomposable) elements are orthogonal
 *  @param v,w Elements in the vector space whose simple elements have type V 
 *  @return The expression \f$\sum a_ib_i\f$, where \f$ v=\sum a_iv_i, w=\sum b_iv_i\f$, and the \f$v_i\f$ are simple. 
 
 *  @note To work with the exterior algebra over V, use TrivialPairingOperator<Lambda<V> >
 */
template<typename V> ex TrivialPairing(ex v,ex w)
{
//	return TrivialPairingImpl(v,w);

	static TrivialPairingOperator<V> oper;
	return TrivialPairingOperator<V>::BilinearOperator(v.expand(),w.expand(),&oper).expand();

			
}

}

/** @} */

#endif /*BILINEAR_H_*/
