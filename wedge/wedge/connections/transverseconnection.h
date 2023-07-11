/*******************************************************************************
 *  Copyright (C) 2007-2023 by Diego Conti, diego.conti@unipi.it 
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
#ifndef TRANSVERSECONNECTION_H
#define TRANSVERSECONNECTION_H
/** @ingroup RiemannianGeometry 
 *  @{ */

/** @file oneill.h 
 *  @brief Levi Civita connection on tranverse Riemannian structures
 */

#include "wedge/structures/transversestructure.h"
#include "wedge/connections/pseudolevicivita.h"
#include "wedge/convenience/simplifier.h"

namespace Wedge {

/** @brief The Levi Civita connection on a Riemannian foliation/submersion.

 Implements O'Neill formulas.
*/
class TransverseLeviCivitaConnection : public PseudoLeviCivitaConnection {
	TransverseRiemannianStructure g;	//horrible. we are keeping a double copy!
protected:
	matrix RicciAsMatrix() const {return PseudoLeviCivitaConnection::RicciAsMatrix();}
	matrix CurvatureForm() const {return PseudoLeviCivitaConnection::CurvatureForm();}
public:
	TransverseLeviCivitaConnection(const Manifold* manifold, const TransverseRiemannianStructure& _g,const Name& christoffel=N.Gamma) : 
		Connection(manifold,_g.e(),true), PseudoLeviCivitaConnection(manifold, _g, christoffel) , g(_g) {}
	
/** @brief Project a form or vector field on the vertical part
*/
	template<typename Section> ex V(ex X) const {return g.V<Section>(X);}

/** @brief Project a form or vector field on the horizontal part
*/
	template<typename Section> ex H(ex X) const {return g.H<Section>(X);}

/** @brief Compute the tensor T, \f$T(X,Y)=(\nabla_{X_V} Y_H)_V+(\nabla_{X_V} Y_V)_H\f$
 * @param X,Y %Vector fields on the total space 
*/
	ex T(ex X, ex Y) const;
/** @brief Compute the tensor A, \f$A(X,Y)=(\nabla_{X_H} Y_H)_V+(\nabla_{X_H} Y_V)_H\f$
 * @param X,Y %Vector fields on the total space 
*/
	ex A(ex X, ex Y) const;

/** @brief Compute the Ricci tensor of the Levi Civita connection on \f$M\f$
 * @param useRicciFormula If true, compute the ricci on the total space and apply the formula, otherwise compute the curvature first and apply the definition
 * @return The Ricci tensor as a matrix
 * 
 * The result represents the Ricci tensor as a matrix.
 * @sa [Besse, Einstein manifolds]
**/
	virtual matrix BaseRicciAsMatrix(bool useRicciFormula=false,const Simplifier& simplifier=default_simplifier) const;

/** @brief Compute the curvature 2-form of the connection on the base \f$M\f$
 *  @return The curvature 2-form as a GiNaC::matrix (whose indices are zero-based)
 * @exception WedgeException<std::logic_error> Thrown if the manifold does not know how to take d of its forms
 * 
 * Computes the curvature form  \f$\Omega^M\f$ in terms of the curvature of the total space \f$\Omega^M\f$ by the O'Neill formula
 * \f$2\Omega^X_{ij}(X,Y)=2\Omega^M_{ij}(X,Y)-2\langle A(X,Y),A_{ij}\rangle +\langle A(Y,e_i),A(X,e_j)\rangle -\langle A(X,e_i),A(Y,e_j)\rangle\f$. The
 * curvature \b tensor \f$R\f$ can be obtained by \f$\langle R(X,Y)e_i,e_j\rangle =2\Omega_{ji}(X,Y)\f$ 
 * 
 * @sa [O'Neill, The fundamental equation of a submersion], but notice O'Neill defines the curvature with the opposite sign
 *
 * @note We cannot just overload CurvatureForm because it's virtual
 */
	virtual matrix BaseCurvatureForm (const Simplifier& simplifier=default_simplifier) const;
};

}

#endif

