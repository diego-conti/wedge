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
#ifndef ADJOINT_H
#define ADJOINT_H_
/** @ingroup Representations */
/** @{
  * @file adjoin.h
  * @brief Explicit representations of the Lie algebra \f$\mathfrak{gl}(n,\mathbb{R})\f$
*/

#include "linearaction.h"
#include "liegroup.h"
namespace Wedge {

/** @brief The standard representation of \f$\mathfrak{gl}(n)\f$ on \f$\mathbb{R}^n\f$
*/
class AdjointRepresentation {
	const LieGroup* G;
	lst subs;	//the substitutions mapping an element of so(n) to the corresponding LinearAction object
public:
	typedef VectorField ActsOnType;

/** @brief Construct the adjoint representation of \f$G\f$
 * @param G A Lie group
 *
 * @warning Caller must ensure that the pointer G remains valid.
*/
	AdjointRepresentation(const LieGroup* G) {
		this->G=G;
	}
/** @brief Action of an element of \f$\mathfrak{g}\f$ on a vector
 * @param A An element of \f$\mathfrak{g}\f$
 * @param v An element of a representation induced by the adjoint representation
 * @return The element \f$ A\cdot v\f$
 */
	template<typename R> ex Action(ex A, ex v) const {
		ExVector ei_goes_to(G->Dimension());
		for (int i=1;i<=ei_goes_to.size();++i)
		{
			ei_goes_to(i)=G->LieBracket(A,G->e(i));
		}
		LinearAction<VectorField> a(G->e(),ei_goes_to);
		return AlgebraAction<VectorField,R>(a,v);
	}
/** @brief Computes the conditions on a generic element in the representation, for the algebra to act trivially on it.
 * @param container A container where the equation are to be stored
 * @param v A (generic) element of a representation of \f$\mathfrak{g}\f$
 * @return A reference to container
 */
	template<typename R, typename Container>
		Container& GetEquationsTrivialAction(Container& container, ex v) const
	{
		for (int i=1;i<=G->Dimension();++i)
			GetCoefficients<R>(container,Action<R>(G->e(i),v));
		return container;
	}
};




}
/** @} */

#endif /*ADJOINT_H_*/
