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
#ifndef RIEMANNIANCONNECTION_H_
#define RIEMANNIANCONNECTION_H_
#include "wedge/structures/riemannianstructure.h"
#include "wedge/connections/connection.h"
/** @ingroup RiemannianGeometry */ 

/** @{ 
 * @file riemannianconnection.h
 * @brief Riemannian connections
 */


namespace Wedge {

class RiemannianConnection;

namespace internal {
template<> class CovariantDerivative<Spinor> :
	public IBilinearOperator<LinearOperator<VectorField>, Leibniz<Spinor,Function> >
{
	const RiemannianConnection& connection;
public:
	CovariantDerivative(const RiemannianConnection& c) : connection(c) {}
	ex Apply(const VectorField& vfield, const Spinor& Alpha) const; 
	ex Apply(const VectorField& vfield, const Function& Alpha) const;	
};
}

/** @brief Specialization of Connection for Riemannian connections 
*/  
class RiemannianConnection: 
	public virtual Connection
{
	friend class internal::CovariantDerivative<Spinor>;
protected:
	RiemannianStructure g;
public:
/** @brief Construct a Riemannian connection on a Riemannian manifold
 *  @param manifold The manifold on which the connection is defined
 *  @param g A riemannian structure on manifold
 *  @param christoffel The symbol to use for the connection parameter
 * 
 * In Wedge, Riemannian structures are represented as orthonormal frames. The frame underlying
 * g is used to represent this connection as a matrix.
 * 
 * @warning The pointer to manifold is stored internally. Caller is responsible for making sure
 * that it remains valid until the Connection objects is destroyed.
*/
	RiemannianConnection(const Manifold* manifold, const RiemannianStructure& g, const Name& christoffel=N.Gamma);

	RiemannianConnection(const Manifold* manifold, const RiemannianStructure& g, bool DoNotInitialize); ///< Overloaded constructor for subclasses that initialize the connection matrix themselves

/** @copydoc Connection::Nabla
 * 
 * This function overloads Connection::Nabla, which only works for differential forms and vector fields. For RiemannianConnection's, Section may also be Spinor, in which case the formula is \f[\nabla \psi = -\frac12\sum_{i<j} \omega_{ij} \otimes e_i\cdot e_j\cdot\psi\f].
**/
	template<typename Section> ex Nabla(ex X, ex alpha) const
	{
		internal::CovariantDerivative<Section> d(*this);
		return internal::CovariantDerivative<Section>::BilinearOperator(X,alpha,&d).expand();
	}

/** @copydoc Connection::DeclareNabla
 *
 * This function overloads Connection::DeclareNabla, which only works for differential forms and vector fields. For RiemannianConnection's, Section may also be Spinor
**/
	template<typename Section> void DeclareNabla(ex X, ex alpha, ex beta)
	{	
		alpha=alpha.expand().normal();
		beta=beta.expand().normal();
		if (alpha.is_zero()) {
			if (beta.is_zero()) return;
			else throw InconsistentDeclaration(__FILE__,__LINE__,"Christoffel symbols");
		}
		ex GenericNablaAlpha=Nabla<Section>(X,alpha);
		DeclareZero(GenericNablaAlpha-beta);
	}
};



} /** @} */
#endif /*RIEMANNIANCONNECTION_H_*/
