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
#include "pseudoriemannianstructure.h"

namespace Wedge {

class RiemannianStructure;


/** @brief Riemannian metric on a manifold, viewed as an O(n)-structure
 */
 
class RiemannianStructure : public PseudoRiemannianStructureByOrthonormalFrame {
public:
/** @brief Define a Riemannian structure on a manifold, represented by the choice of a frame
 *  @param manifold The manifold on which the structure is defined.
 *  @param orthonormalframe The orthonormal frame defining the metric
 * 
 *  @warning Caller must ensure that the pointer remains valid.
 */ 
	RiemannianStructure(const Manifold* manifold, const Frame& orthonormal_frame);

/**
   @brief Compute the scalar product of two vector fields, differential forms or spinors
   @param op1,op2 Differential forms, vector fields or spinors
   @return Their scalar product as a real number (function)
 */
	template<typename T> ex ScalarProduct (ex op1, ex op2) const;
	using PseudoRiemannianStructureByOrthonormalFrame::ScalarProduct;

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
 	
private:
	class RiemannianHookOperator;
	struct Deleter {
  		void operator()(RiemannianHookOperator* r);    
	};
 	unique_ptr<RiemannianHookOperator,Deleter> hookOperator;
};

template<> ex RiemannianStructure::ScalarProduct<VectorField> (ex op1, ex op2) const;
template<> ex RiemannianStructure::ScalarProduct<DifferentialForm> (ex op1, ex op2) const;
template<> ex RiemannianStructure::ScalarProduct<Spinor> (ex op1, ex op2) const;



} /** @} */

#endif
