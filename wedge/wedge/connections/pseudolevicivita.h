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
#ifndef PSEUDOLEVICIVITA_H_
#define PSEUDOLEVICIVITA_H_

/** @ingroup RiemannianGeometry */ 

/** @{ 
 * @file pseudolevicivita.h
 * @brief Levi-Civita connection on a pseudoriemannian manifold
 */
 
#include "wedge/structures/pseudoriemannianstructure.h"
#include "wedge/connections/torsionfreeconnection.h"

namespace Wedge {

class CovariantDerivativeSpinor;

/** @brief The Levi-Civita connection  on a pseudo riemannian manifold that (already) knows how to 
 * take d of its forms.*/  
class PseudoLeviCivitaConnection : public virtual Connection, public TorsionFreeConnection<true> {
	struct Deleter {	//necessary to use unique_ptr with forward reference
		void operator() (CovariantDerivativeSpinor* p);
	};
	unique_ptr<CovariantDerivativeSpinor,Deleter> covariant_derivative_spinor;
public:
/** @brief Construct the Levi-Civita connection
 *  @param manifold The manifold on which the connection is defined
 *  @param structure A pseudoriemannian structure on manifold
 *  @param christoffel The symbol to use for the connection parameters
 * 
 * @warning The pointer argument is stored internally. Caller is responsible for making sure 
 * that it remains valid.
 * @exception WedgeException<std::logic_error> Thrown if manifold does not know how to take d of its forms
*/
	PseudoLeviCivitaConnection(const Manifold* manifold, const PseudoRiemannianStructure& structure, const Name& christoffel=N.Gamma);
	
	template<typename Section> ex Nabla(ex X, ex alpha) const {
		return Connection::Nabla<Section>(X,alpha);
	}
};

template<> 
ex PseudoLeviCivitaConnection::Nabla<Spinor>(ex X, ex psi) const;

}

#endif /*PSEUDOLEVICIVITA_H_*/

