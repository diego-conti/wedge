/***************************************************************************
 *   Copyright (C) 2008, 2009 by Diego Conti				   *
 *   diego.conti@unimib.it                                                 *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 3 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#ifndef STABILIZER_H_
#define STABILIZER_H_

#include "linearaction.h"
namespace Wedge {

/** @ingroup Representations
/** @{ 
  * @file stabilizer.h 
  * @brief Stabilizers of Lie algebra actions
*/

/** @brief Return the dimension of the \f$G\f$-orbit of an element in a representation \f$R\f$
 * @param G The object representing the Lie group
 * @param representation The representation inducing the representation \f$R\f$
 * @param v An element in the representation \f$R\f$
 * @return The dimension of the \f$ G\f$-orbit of \f$v\f$ in \f$R\f$
 */

   template<typename R, typename Representation> 
int OrbitDimension(const LieGroupWithoutParameters& G,const Representation& representation, ex v)
{
	exvector image;
	for (exvector::const_iterator i=G.e().begin();i!=G.e().end();++i)
		image.push_back(representation.template Action<R>(*i,v));
	LOG_DEBUG(image);
	return VectorSpace<R>(image).Dimension();
}

/** @brief Return the Lie algebra of the stabilizer of an element in a representation \f$R\f$ with respect to the action of a group \f$G\f$
 * @param G The object representing the Lie group
 * @param representation The representation inducing the representation \f$R\f$
 * @param v An element in the representation \f$R\f$
 * @return The subspace of \f$\mathfrak{g}\f$ that fixes \f$v\f$
 */

   template<typename R, typename Representation> 
Subspace<DifferentialForm > StabilizerAlgebra(const LieGroupWithoutParameters& G,const Representation& representation, ex v)
{
	VectorSpace<DifferentialForm> V=G.pForms(1);
	ex image=representation.template Action<R>(V.GenericElement(),v);
	list<ex> eqns;
	GetCoefficients<DifferentialForm>(eqns,image);
	return V.SubspaceFromEquations(eqns.begin(),eqns.end());
}

/** @brief Return the stabilizer of an element in a representation \f$R\f$ with respect to the action of a group \f$G\f$
 * @param G The object representing the Lie group
 * @param representation The representation inducing the representation \f$R\f$
 * @param v An element in the representation \f$R\f$
 * @return The stabilizer of \f$\alpha\f$ in \f$G\f$ as an abstract Lie group (as opposed to a subgroup of G)
 * @throw DiscreteManifold if the stabilizer algebra is zero, since AbstractLieSubgroup can only represent Lie groups of positive dimension
 */

   template<typename R, typename Representation> 
AbstractLieSubgroup<false> Stabilizer(const LieGroupWithoutParameters& G,const Representation& representation, ex v)		
{
	return  AbstractLieSubgroup<false>(G,StabilizerAlgebra<R>(G,representation,v).e());
}


/** @} */
}
#endif /*STABILIZER_H_*/

