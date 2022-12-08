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
#ifndef MANIFOLD_H
#define MANIFOLD_H

/** @defgroup Manifolds Manifolds, differential forms, Lie groups */
/**  @{ */
/**  @file manifold.h 
 * @brief Generic manifolds 
 */

#include "wedge/linearalgebra/derivation.h"
#include "wedge/linearalgebra/vectorspace.h"
#include "wedge/linearalgebra/linear.h"
#include "wedge/linearalgebra/pforms.h"
#include "wedge/manifolds/differentialform.h"
#include "wedge/manifolds/function.h"

namespace Wedge {
using namespace  GiNaC;

/** @brief Type representing a possibly non-standard basis of differential one-forms
 * 
 * Non-standard means that the elements may not be simple.
 */ 
 
typedef Basis<DifferentialOneForm,DefaultLinAlgAlgorithms> Frame;

/**
 * @brief Represents a parallelizable manifold
 *
 *  A manifold is described by a basis of 1-forms, which one might think of as a local frame. This is an abstract class; a concrete manifold is defined
 * by subclassing Manifold, and setting up a basis of 1-forms in the constructor. The exterior derivative of a differential form can be computed
 * with Manifold::d
 * 
 * @warning The return value of most members of Manifold (and its derived classes) are only valid in the scope where the Manifold object 
 * itself is valid.
 * 
 * \todo Allow for non-parallelizable manifolds arising as quotients, e.g. homogeneous spaces
*/

class Manifold : public IBilinearOperator<LinearOperator<VectorField>,DerivationOver<DifferentialForm,Function,false>  > {
public:
	/** @brief Exception thrown by d operator when the Manifold object does not know how to take d of some differential form
	 *  
	 */
	struct dException : public WedgeException<std::runtime_error> {
		ex alpha;
		dException(const char* in_file,int at_line, ex _alpha) throw() : 
			WedgeException<std::runtime_error> (std::string("Don't know how to take d of ")+ToString(_alpha),in_file,at_line), alpha(_alpha) {}
		~dException() throw() {}
	}; 	

	struct ddNotZeroException  : public WedgeException<std::runtime_error> {
		ex alpha;
		string structure_constants;
		ddNotZeroException(const char* in_file,int at_line, ex _alpha, const Manifold& M) ;
		const char* what() const throw();
		~ddNotZeroException() throw() {}
	};
	Manifold();
	virtual ~Manifold() {}

  /**
   * @brief Returns the dimension of the manifold
   * @returns The dimension of the manifold
  */
	virtual int Dimension() const {return e().size();}

  /**
   * @brief Returns the k-th element of a global basis of 1-forms
   * @param k An integer in the range [1,dimension]
   * @return An ex object representing a differential 1-form
  */
	ex e(OneBased k) const {return e()(k);}

  /**
   * @brief Returns a frame for this manifold
   * @return A Frame object, consisting of a basis of one-forms
   * 
   * The elements of the frame need not be simple.
   * @note Since the frame is returned by reference, the caller is allowed to iterate through the frame STL-style
  */
	virtual	const Frame& e() const=0;

  /**
   * @brief The Hodge star operator induced by the identifications of forms with vector fields
   * @param alpha A differential form on this manifold
   * @returns The interior product \f$\alpha\lrcorner e^1\wedge\dotsb \wedge e^n\f$
   * 
   * FIXME alpha is not allowed to have degree zero at the moment
   * This operator does not depend on the choice of a metric or frame
   * @sa Hook
  */
	ex HodgeStar(ex alpha) const;

  /**
   * @brief The exterior derivative operator \f$d\f$
   * @param alpha A differential form on this manifold
   * @returns The exterior derivative \f$d\alpha\f$ of \f$\alpha\f$ 
   * @exception Manifold::dException This Manifold object does not know how to take d of alpha
   * 
   * The default implementation gives zero on constants and throws an exception otherwise.
   * 
  */
	virtual ex d(ex alpha) const;
  /**
   * @brief Return the vector space of p-forms
   * @param p An integer in the range [1,dimension]
   * @returns The vector space consisting of forms of degree p
   * 
   * @remark p has to be greater than zero because functions are not DifferentialForm's.
  */
	VectorSpace<DifferentialForm> pForms(int p) const;

  /**@brief Verify that \f$d^2\f$ is zero
   * 
   * @exception WedgeException<ddNotZeroException> \f$d^2\f$ is not zero
  */
	virtual void Check_ddZero();

  /**@brief Return true if member function d() can be used to compute the exterior derivative
   * @return true if this manifold object knows how to compute d of its forms   
   *
   * This function returns false for those subclasses that do not reimplement d()
  */
	virtual bool KnowsHowToCompute_d() const {return false;} 

/** @brief Compute the Lie Derivative 
  * @param X An ex representing a vector field
  * @param f An ex representing a function or a differential form
  * @returns The Lie derivative \f$\mathcal{L}_Xf\f$
  *
  * @note The Lie derivative does not respect our identification of DifferentialOneForm's with VectorField's. Thus, this function only works for differential forms; for vector fields, use LieBracket instead.
 */
	virtual ex LieDerivative(ex X, ex f) const;

/** @brief Compute the Lie bracket of two vector fields
 * @param X,Y Vector fields on the manifold
 * @return The vector field \f$[X,Y]\f$
 */
	virtual ex LieBracket(ex X, ex Y) const;

/** @brief For %internal use
*/
	virtual ex Apply(const VectorField& X,const Function& f) const
	{
		return f.Derive(X,*this);
	}
/** @brief For %internal use
*/
	ex Apply(const VectorField& X,const DifferentialForm& alpha) const
	{
		return Hook(X,d(alpha))+d(Hook(X,alpha));
	}
/** @brief For %internal use
*/
	ex Apply(const VectorField& X,const VectorField& alpha) const
	{
		//due to the way DerivationOver works, alpha is a DifferentialOneForm in disguise
		assert(is_a<DifferentialOneForm>(alpha));
		const DifferentialOneForm& beta=static_cast<const DifferentialOneForm&>(alpha);
		return Hook(X,d(beta))+d(Hook(X,beta));
	}
	
protected:
	ex constant_function;	///< The constant function \f$f\equiv 1\f$ on this manifold
private:
	bool operator==(const Manifold&) const;	///< Not defined
};


} /** @} */
#endif
