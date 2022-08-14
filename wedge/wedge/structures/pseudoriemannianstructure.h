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
#ifndef PSEUDORIEMANNIANSTRUCTURE_H
#define PSEUDORIEMANNIANSTRUCTURE_H
/** @ingroup RiemannianGeometry

/** @{ 
 * @file pseudoriemannianstructure.h
 * @brief Pseudoriemannian metrics as G-structures
 */
 
#include "../structures/gstructure.h"
#include "wedge/manifolds/manifold.h"
#include "wedge/linearalgebra/vectorspace.h"
#include "wedge/linearalgebra/bilinearform.h"

namespace Wedge {


/** @brief Pseudoriemannian metric on a manifold, 
 * 
 * This is an abstract base class
 */
 
class PseudoRiemannianStructure : public virtual GStructure
{
	ex Hook(ex) const = delete;	//Not implemented. Use ScalarProduct().Interior() instead.
public:
	virtual const BilinearForm& ScalarProduct() const=0;

	ex HodgeStar(ex form) const {return ScalarProduct().Interior(form,ncmul(e()));}

/** @brief Define a Pseudoriemannian structure on a manifold, represented by the choice of a frame
 *  @param manifold The manifold on which the structure is defined.
 *  @param orthonormalframe A coframe with respect to which the metric is defined; could be orthonormal or not.
 * 
 *  @warning Caller must ensure that the pointer remains valid.
 */ 
	PseudoRiemannianStructure(const Manifold* manifold, const Frame& frame) :
		GStructure(manifold,frame) {}

/** @brief Overloaded copy constructor
 *  
 * @warning Every PseudoriemannianStructure object contains a pointer to a Manifold object, representing the underlying manifold. The pointer is copied by
 * this constructor, so the caller must ensure that it remains valid.
 *
 * @todo Replace all pointers with shared_ptr
 */ 
	PseudoRiemannianStructure(const PseudoRiemannianStructure& o) : GStructure(o) {
		/* write a message to the log, since use of copy constructor may lead to subtle bugs, i.e.
		 * PseudoriemannianStructure f() {
		 * 	ConcreteManifold M(5);
		 *  return PseudoriemannianStructure(&M,M.e()); //crash
		 * }
		 */
		LOG_WARN("PseudoRiemannianStructure copy constructor invoked; beware the pointer.");
	}

	virtual pair<ex,matrix> DecomposeRicci(matrix ricci) const=0;
};

/** @brief Pseudoriemannian metric on a manifold, represented by an orthonormal coframe
 */
 
class StandardPseudoRiemannianStructure : public PseudoRiemannianStructure {
	const int p_;
	const StandardScalarProduct scalar_product;
public:
/** @brief 
 *  @param manifold The manifold on which the structure is defined.
 *  @param orthonormal_frame A coframe with respect to which the metric is defined; could be orthonormal or not.
 *  @param p An integer that determines the signature
 *
 * The orthonormal frame e_1,..., e_n is assumed to satisfy <e_i,e_j>=\delta_{ij} for i\leq p, -\delta_{ij} for i>p
*/
	StandardPseudoRiemannianStructure(const Manifold* manifold, const Frame& orthonormal_frame, int p) : 
		GStructure{manifold, orthonormal_frame}, PseudoRiemannianStructure(manifold,orthonormal_frame), p_{p}, scalar_product{orthonormal_frame,p} {}
	const BilinearForm& ScalarProduct() const {return scalar_product;}

	pair<ex,matrix> DecomposeRicci(matrix ricci) const override;
};


/** @brief Pseudoriemannian metric on a manifold, represented by a matrix relative to some frame
 */

class PseudoRiemannianStructureByMatrix : public PseudoRiemannianStructure {
	const ScalarProductDefinedByMatrix scalar_product;
public:
/** @brief
 *  @param manifold The manifold on which the structure is defined.
 *  @param frame A coframe with respect to which the metric is defined; could be orthonormal or not.
 *  @param m The metric relative to the coframe
*/
	PseudoRiemannianStructureByMatrix(const Manifold* manifold, const Frame& frame, const matrix& m) :
		GStructure{manifold, frame}, PseudoRiemannianStructure(manifold,frame), scalar_product{frame, m} {}
	const BilinearForm& ScalarProduct() const {return scalar_product;}

	pair<ex,matrix> DecomposeRicci(matrix ricci) const override {throw NotImplemented(__FILE__,__LINE__,"PseudoRiemannianStructureByMatrix::DecomposeRicci");}
};

} /** @} */

#endif
