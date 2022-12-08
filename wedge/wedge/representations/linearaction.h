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
 #ifndef LINEARACTION_H_
#define LINEARACTION_H_
/** @defgroup Representations Representations of Lie groups and Lie algebras  */
/** @{ 
  * @file linearaction.h 
  * @brief Base class for linear actions (representations)
*/

#include "wedge/manifolds/concretemanifold.h"
#include "wedge/linearalgebra/affinebasis.h"
#include "wedge/liealgebras/liesubgroup.h"
#include "wedge/representations/gl.h"
#include "wedge/linearalgebra/lambda.h"
#include "wedge/linearalgebra/linear.h"

using namespace  GiNaC;
namespace Wedge {

namespace internal { class LinearActionVisitor; }

/** @brief Non-template base class of template class LinearAction<T>
 *
 * An object of type LinearActionBase represents a linear endomorphism. This is an algebraic class: one can add and multiply LinearActionBase's,
 * obtaining the natural sum and composition of endomorphisms.
 * @sa LinearAction
*/
class LinearActionBase : public Register<LinearActionBase, Vector>::Algebraic {
	friend class internal::LinearActionVisitor;
public:
	LinearActionBase() {}	
	LinearActionBase(const lst& subs) : linear_subs(subs) {}
	/* Functions used internally by the GiNaC framework: */
	static const char* static_class_name() {return "LinearAction";}
	unsigned precedence() const {return 50;}
	int compare_same_type(const GiNaC::basic & other) const {
		const LinearActionBase& o =static_cast<const LinearActionBase&>(other);
		return linear_subs.compare(o.linear_subs);
	}
protected:
	lst linear_subs;
	ex eval_ncmul(const exvector & v) const; ///< Implements the composition of endomorphisms
	unsigned return_type() const {return return_types::noncommutative;}
};

/** @brief An endomorphism of a vector space T, representing the linear action of an element of a Lie algebra or group
 *
 * An object of type LinearAction<T> represents a linear endomorphism. This is an algebraic class: one can add and multiply LinearAction<T>'s,
 * obtaining the natural sum and composition of endomorphisms.
 *
 * To be used in conjunction with AlgebraAction.
*/
template<typename T> class LinearAction : public LinearActionBase {
public:
	LinearAction() {}
/** @brief Construct a LinearAction element from a list of substitutions
 * @param subs A list of substitutions of the form v==w, where v and w are vectors of type T
*/
	LinearAction(const lst& subs) : LinearActionBase(subs) {}
/** @brief Construct a LinearAction element from its action in terms of a basis
 * @param from A basis of independent elements spanned by T
 * @param to The images of the elements of the specified basis
*/
	LinearAction(const exvector& from, const exvector& to) : LinearActionBase(LinearMapToSubstitutions<T>(from,to)) {}
/** @brief Initialize the LinearAction element from a list of substitutions
 * @param subs A list of substitutions of the form v==w, where v and w are vectors of type T
*/
	void InitializeFromSubs(const lst& subs);
/** @brief Initialize the LinearAction element from its action in terms of a basis
 * @param from A basis of independent elements spanned by T
 * @param to The images of the elements of the specified basis
*/
	void InitializeFromBasis(const exvector& from, const exvector& to) {linear_subs=LinearMapToSubstitutions<T>(from,to);}
};

template<typename T, typename R> ex AlgebraAction(ex g, ex w);

/** @brief Compute the space of invariants in the exterior algebra over a representation
 * @param representation A class defining a representation, e.g. SO3Representation<T>
 * @param frame A basis of elements of type Lambda1<T>
 * @param p A positive integer
 * @return The space of invariants in \f$Lambda^p T\f$.
*/
template<typename Representation> Subspace<Lambda<typename Representation::ActsOnType> > Invariant_pForms(const Representation& representation, const IBasis<Lambda1<typename Representation::ActsOnType> >& frame, int p) 
{
	VectorSpace<Lambda<typename Representation::ActsOnType> > pforms=pForms(frame,p);
	list<ex> eqns;
	representation.template GetEquationsTrivialAction<Lambda<typename Representation::ActsOnType> >(eqns,pforms.GenericElement());
	return pforms.SubspaceFromEquations(eqns.begin(),eqns.end());
}


/** @brief Compute the space of invariants in the exterior algebra over a representation
 * @param representation A class defining a representation, e.g. SO3Representation<T>
 * @param e A basis of elements of type T, spanning a representation \f$V\f$
 * @return The space of invariant symmetric tensors in \f$V\otimes V\f$.
*/
template<typename Representation> Subspace<Tensor<typename Representation::ActsOnType,typename Representation::ActsOnType> > InvariantMetrics(const Representation& representation, const ExVector& e) 
{
	VectorSpace<Tensor<typename Representation::ActsOnType,typename Representation::ActsOnType> > metrics;
	for (int i=1;i<=e.size();++i)
	{
		metrics.AddGenerator(TensorProduct<typename Representation::ActsOnType, typename Representation::ActsOnType> (e(i),e(i)));
		for (int j=i+1;j<=e.size();++j)
			metrics.AddGenerator(TensorProduct<typename Representation::ActsOnType, typename Representation::ActsOnType> (e(i),e(j))+TensorProduct<typename Representation::ActsOnType, typename Representation::ActsOnType> (e(j),e(i)));
	}
	list<ex> eqns;
	representation.template GetEquationsTrivialAction<Tensor<typename Representation::ActsOnType,typename Representation::ActsOnType> >(eqns,metrics.GenericElement());
	return metrics.SubspaceFromEquations(eqns.begin(),eqns.end());
};

namespace internal {

class LinearActionVisitor : public LinearOperator<LinearActionBase> {
	ex w;
public:
	LinearActionVisitor(ex w) {this->w=w;}
	void visit(const LinearActionBase& g)
	{
		Result()=w.subs(g.linear_subs);
	}
};

template<typename T,typename W> class AlgebraActionVisitor;

template<typename T> class AlgebraActionVisitor<T,T> {
	ex g;
public:
	AlgebraActionVisitor(ex g){this->g=g;}
	ex RecursiveVisit(ex w) {
		LinearActionVisitor v(w);
		return v.RecursiveVisit(g);
	}
};

template<typename T,typename W> class AlgebraActionVisitor<T,Lambda<W> > : public Derivation<Lambda<W>,false>
{
	ex g;
public:
	AlgebraActionVisitor(ex g) {this->g=g;}
	void visit(const W& w)
	{
		this->Result()=AlgebraAction<T,W>(g,w);
	};
};

template<typename T,typename V,typename W> class AlgebraActionVisitor<T,Tensor<V,W> >: public LinearOperator<Tensor<V,W> >
{	
	ex g;
public:
	AlgebraActionVisitor(ex g) {this->g=g;}
	void visit(const Tensor<V,W>& p)
	{
		ex left=TensorProduct<V,W>(AlgebraAction<T,V>(g,p.v),p.w);
		ex right=TensorProduct<V,W>(p.v,AlgebraAction<T,W>(g,p.w));
		this->Result()=TensorProduct<V,W>(AlgebraAction<T,V>(g,p.v),p.w)+
			TensorProduct<V,W>(p.v,AlgebraAction<T,W>(g,p.w));
	}
};


}


/** @brief Compute the action of an element of a Lie group on a vector space
 * @param g An element of a Lie group together with a representation; combination of type LinearAction<T> for some T
 * @param w A vector in the representation T, or any induced representation
 * @result The action of g on w.
 *
 * Implemented by means of a substitution, so that the vector w can be in any induced representation, e.g. Lambda<T> or Tensor<T,T>. Thus the action of g on w may be polynomial in g. 
 * @sa AlgebraAction
*/
ex GroupAction(ex g, ex w);

/** @brief Compute the action of an element of a Lie Algebra on a vector space
 * @param g A vector of type LinearAction<T>, representing an element of a Lie algebra together with a representation
 * @param w A vector of the representation R
 * @result The (infinitesimal) action of g on w.
 *
 * The type R can be anything like Lambda<T>, Tensor<T,Lambda<T> >  and so on. Notice that this is the  induced action that is linear in g.
 * @sa GroupAction
*/
template<typename T, typename R> ex AlgebraAction(ex g, ex w)
{
	internal::AlgebraActionVisitor<T,R> v(g);
	return v.RecursiveVisit(w.expand()).expand();
}

/** @brief An endomorphism of a vector space T, representing the linear action of an element of a Lie algebra or group
 *
 * An object of type LinearAction<T> represents a linear endomorphism. This is an algebraic class: one can add and multiply LinearAction<T>'s,
 * obtaining the natural sum and composition of endomorphisms.
 *
 * To be used in conjunction with AlgebraAction.
 */
template<typename T>
inline void LinearAction<T>::InitializeFromSubs(const lst& subs) {
	linear_subs = subs;
}

}



/** @} */

#endif /*LINEARACTION_H_*/
