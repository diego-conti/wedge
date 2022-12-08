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
#ifndef MANIFOLDWITH_H_
#define MANIFOLDWITH_H_


/** @ingroup RiemannianGeometry */ 

#include "../structures/structures.h"

/** @{ 
 * @file manifoldwith.h 
 * @brief Generic manifolds with a given G-structure
 */
 
namespace Wedge {
	
/** @brief A generic Riemannian manifold with a G-structure
 * 
 * A ManifoldWith object represents the generic manifold with a G-structure satisfying certain
 * intrinsic torsion conditions, which can be imposed by using Declare_d and DeclareNabla.
 * 
 * In contrast to Has_dTable, the action of d is not determined from a dTable, but from 
 * the Levi Civita connection.
 * 
 * @remark Instances of this class are not required to satisfy \f$d^2=0\f$ for all values in the parameters. The conditions 
 * in the parameters corresponding to \f$d^2=0\f$ can be obtained with a call to TorsionFreeConnection<false>::GetEquations_ddZero()
*/

template<class Structure> class ManifoldWith : public ConcreteManifold, public Structure {
	LeviCivitaConnection<false> connection;	
public:
/** @brief Constructor for GStructures with a fixed dimension
 */  
	ManifoldWith() : ConcreteManifold(Structure::dimension), Structure(this,ConcreteManifold::e()),		
		connection(this,*this)
		{}
		
/** @brief Constructor for GStructures with variable dimension
 *  @param dimension The dimension of the manifold
 */ 
	ManifoldWith(int dimension) : ConcreteManifold(dimension), Structure(this,ConcreteManifold::e()),
		connection(this,*this)
		{}

/** @brief Return the Levi Civita connection
 *
 * To declare conditions on the Levi Civita connection, use ManifoldWith::DeclareZero
 */
	const LeviCivitaConnection<false>& LeviCivita() const {
		return connection;
	}
	
	const Frame& e() const {return ConcreteManifold::e();}	
	ex e(OneBased k) const {return ConcreteManifold::e(k);}

/** @brief Impose conditions on the Christoffel symbols
 *  @param alpha A differential form
 *  @param dalpha d of alpha
 */
	void Declare_d(ex alpha, ex dalpha)
	{
		connection.Declare_d(alpha,dalpha);
	}

/** @brief Impose conditions on the Christoffel symbols
 *  @param X A vector field 
 *  @param alpha A section of the bundle associated to the type T
 *  @param nabla_Xalpha The covariant derivative \f$\nabla_X\alpha\f$
 */	
	template<typename T> void DeclareNabla(ex X, ex alpha, ex nabla_Xalpha)
	{
		connection.DeclareNabla<T>(X, alpha,nabla_Xalpha);	
	}
/** @brief Eliminate some connection parameters by imposing linear conditions on an expression
 * @param alpha An expression depending linearly on the parameters
 * @exception InconsistentDeclaration Thrown if alpha cannot be zero for any choice of the parameters
 * @exception std::invalid_argument Thrown by %GiNaC framework if the dependence on the parameters is not linear
 * 
 * This function solves \f$\alpha=0\f$ and updates the parameters accordingly.
 * 
 * @sa HasParameters::DeclareZero
 */
	void DeclareZero(ex alpha)
	{
		connection.DeclareZero(alpha);
	}
	
/** @overload
 */
	template<typename Iterator> void DeclareZero(Iterator begin, Iterator end)
	{
		connection.DeclareZero(begin,end);
	}		
/** @brief Compute covariant derivative
 *  @param X A vector field 
 *  @param alpha A section of the bundle associated to the type T
 * 	@return The covariant derivative \f$\nabla_X\alpha\f$
**/
	template<typename T> ex Nabla(ex X, ex alpha) const
	{
		return connection.Nabla<T>(X, alpha);
	}
	
	ex d(ex alpha) const 
	{
		return connection.d(alpha);
	}

	bool KnowsHowToCompute_d() const {return true;} 

	ex LieBracket(ex X, ex Y) const {
		return Nabla<VectorField>(X,Y)-Nabla<VectorField>(Y,X);
	}

//Overloading is needed since HodgeStar is defined both in RiemannianStructure and in Manifold. Since the frame is the standard frame, we can safely use the one in Manifold
	ex HodgeStar(ex alpha) const {return Manifold::HodgeStar(alpha);}
};
} /** @} */
#endif /*MANIFOLDWITH_H_*/
