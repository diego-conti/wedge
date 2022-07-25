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

#ifndef RIEMANNIANSTRUCTURE_H
#define RIEMANNIANSTRUCTURE_H
/** @defgroup RiemannianGeometry Riemannian geometry, special geometries, connections, curvature*/ 

/** @{ 
 * @file riemannianstructure.h
 * @brief Riemannian metrics as G-structures
 */
 
#include "../structures/gstructure.h"
#include "../structures/spinor.h"
#include "wedge/manifolds/manifold.h"
#include "wedge/linearalgebra/vectorspace.h"

namespace Wedge {

class RiemannianStructure;

namespace internal { 
class RiemannianHookOperator : public IBilinearOperator<AssociativeOperator<DifferentialForm>,Derivation<DifferentialForm>  > 
{
	const RiemannianStructure* structure;
public:
	RiemannianHookOperator(const RiemannianStructure* s) {structure=s;} 
	ex Apply(const VectorField& left,const VectorField& right) const;
};

}

/** @brief Riemannian metric on a manifold, viewed as an O(n)-structure
 */
 
class RiemannianStructure : public GStructure,
	public IBilinearOperator<LinearOperator<Spinor>,LinearOperator<Spinor> >, 	//scalar product 
	public IBilinearOperator<LinearOperator<DifferentialForm>,LinearOperator<DifferentialForm> >, 	//scalar product
	public IBilinearOperator<AssociativeOperator<DifferentialForm>,LinearOperator<Spinor> > 	//Clifford multiplication
{
	friend class internal::RiemannianHookOperator;
public:
/** @brief Define a Riemannian structure on a manifold, represented by the choice of a frame
 *  @param manifold The manifold on which the structure is defined.
 *  @param orthonormalframe The orthonormal frame defining the metric
 * 
 *  @warning Caller must ensure that the pointer remains valid.
 */ 
	RiemannianStructure(const Manifold* manifold, const Frame& orthonormalframe) :
		GStructure(manifold,orthonormalframe), hookOperator(this) {}

/** @brief Overloaded copy constructor
 *  
 * @warning Every RiemannianStructure object contains a pointer to a Manifold object, representing the underlying manifold. The pointer is copied by
 * this constructor, so the caller must ensure that it remains valid.
 *
 * @todo Replace all pointers with boost:shared_ptr
 */ 
	RiemannianStructure(const RiemannianStructure& o) : GStructure(o), hookOperator(this) {
		/* write a message to the log, since use of copy constructor may lead to subtle bugs, i.e.
		 * RiemannianStructure f() {
		 * 	ConcreteManifold M(5);
		 *  return RiemannianStructure(&M,M.e()); //crash
		 * }
		 */
		LOG_WARN("RiemannianStructure copy constructor invoked; beware the pointer.");
	}

/** @brief Overloaded assignment operator
 *  
 * @warning Every RiemannianStructure object contains a pointer to a Manifold object, representing the underlying manifold. The pointer is copied by
 * this operator, so the caller must ensure that it remains valid.
 */ 
	RiemannianStructure& operator=(const RiemannianStructure& o) {
		//notice that this assignment operator does not alter hookOperator.
		GStructure::operator=(o);		
		return *this;
	}

/**
   @brief Returns the k-th element of a global basis of complex spinors
   @param k An index in the range [0,DimensionOfSpinorRepresentation())
   @return A section of the complex spinor bundle
*/
	ex u(ZeroBased k) const;	

/** @brief Compute the Clifford action of a form on a spinor
 * @param alpha A differential form
 * @param psi A spinor
 * @return The spinor \f$\alpha\cdot\psi\f$.
	 
 * @remark This is the Clifford multiplication for differential forms. If the adapted frame
 * does not consist of simple elements, Clifford multiplication for vector fields can be obtained as
 * \f[X\cdot\psi =\sum_i e^i(X)\cdot\psi\f]
 * where \f$e^i(X)\f$ is computed using Wedge::TrivialPairing
	
 * @bug If the adapted frame does not consist of simple elements and \f$\alpha\f$ has degree greater than one, the result is incorrect.
 * This is because the Clifford multiplication is only "associative" when one decomposes in orthonormal factors.	   

 * @remark Clifford multiplication is implemented using the formulae of 
 * [Baum; Friedrich; Grunewald; Kath: Twistor and Killing spinors on Riemannian manifolds. Seminarberichte [Seminar Reports], 108. Humboldt Universität, Sektion Mathematik, Berlin, 1990. 179 pp.] 
 * with the basis \f$ u(-\epsilon_1,\dots,-\epsilon_m)=u_{\epsilon_m+2\epsilon_{m-1}+\dots+2^{m-1}\epsilon_1}\f$, 
 *	
 * @sa [Conti-Fino: Calabi-Yau cones from contact reduction, arXiv:0710.4441]
*/
	ex CliffordDot(ex alpha, ex psi) const {
			return 
			IBilinearOperator<AssociativeOperator<DifferentialForm>,LinearOperator<Spinor> >::BilinearOperator(alpha,psi,this);
	}
/**
   @brief Returns the rank of the complex spinor bundle
   @return The number \f$2^{[n/2]}\f$, where \f$n\f$ is the manifold's dimension
 */
	int DimensionOfSpinorRepresentation() const;

/**
   @brief Compute the scalar product of two vector fields, differential forms or spinors
   @param op1,op2 Differential forms, vector fields or spinors
   @return Their scalar product as a real number (function)
 */
	template<typename T> ex ScalarProduct (ex op1, ex op2) const {
		return RealPart(
			IBilinearOperator<LinearOperator<T>,LinearOperator<T> >::BilinearOperator
				(op1.expand(),op2.expand().conjugate(),this)
		);
	}

/**
   @brief Compute the norm of a differential form or spinor
   @param op A form or spinor
   @return Its norm as a real number (or function)
*/
	template<typename T> ex SquareNorm (ex op) const {
			return ScalarProduct<T>(op,op);
	}
 /** @brief The interior product, computed using the metric
  * @param alpha A differential form
  * @param beta A differential form
  * @returns The interior product \f$\alpha\lrcorner\beta\f$
  * 
  * @sa Wedge::Hook
 */
 	ex Hook(ex alpha, ex beta) const;
 	
 	/** @brief The Hodge star, computed using the metric
 	 *  @param alpha A differential form on the manifold to which this RiemannianStructure refers
 	 *  @return The Hodge ''dual'' \f$ *\alpha\f$ of \f$\alpha\f$
 	 * 
 	 * @sa Manifold::HodgeStar
 	 */
 	ex HodgeStar(ex alpha) const;
 	
protected:
 	RiemannianStructure(const Manifold* manifold) : GStructure(manifold), hookOperator(this) {}
private:
 	internal::RiemannianHookOperator hookOperator;
public:
 	ex Apply (const DifferentialForm& op1, const DifferentialForm& op2) const;		///< For %internal use
 	ex Apply (const DifferentialForm&, const VectorField&) const {return 0;}		///< For %internal use
 	ex Apply (const VectorField&, const DifferentialForm&) const {return 0;} 		///< For %internal use
 	ex Apply (const VectorField& op1, const VectorField& op2) const;				///< For %internal use
 	ex Apply (const Spinor& op1, const Spinor& op2) const;							///< For %internal use
 	ex Apply (const DifferentialForm& alpha, const Spinor& spinor) const;			///< For %internal use
 	ex Apply (const VectorField& alpha, const Spinor& spinor) const;				///< For %internal use
};

template<> ex RiemannianStructure::ScalarProduct<VectorField> (ex op1, ex op2) const;


} /** @} */

#endif
