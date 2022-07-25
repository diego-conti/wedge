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
#ifndef REPGL_H_
#define REPGL_H_
/** @ingroup Representations */
/** @{ 
  * @file repgl.h 
  * @brief Explicit representations of the Lie algebra \f$\mathfrak{gl}(n,\mathbb{R})\f$
*/

#include "wedge/representations/linearaction.h"
#include "wedge/representations/gl.h"

namespace Wedge {

/** @brief The standard representation of \f$\mathfrak{gl}(n)\f$ on \f$\mathbb{R}^n\f$
*/
template<typename T>  class GLRepresentation {
	const GL* G;
	ExVector e;
	lst subs;	//the substitutions mapping an element of so(n) to the corresponding LinearAction object
public:	
	typedef T ActsOnType;

/** @brief Construct the standard real representation of \f$\mathfrak{gl}(n)\f$ on \f$\mathbb{R}^n\f$ corresponding to the given choice of basis
 * @param G A point to the object representing the Lie group GL(n)
 * @param frame A basis of n vectors of type T
 *
 * @warning Caller must ensure that the pointer G remains valid.
*/
	GLRepresentation(const GL* G, const ExVector& frame) : e(frame) {
		this->G=G;
		assert(e.size()==G->n());

	}
/** @brief Action of an element of \f$\mathfrak{gl}(n)\f$ on a vector
 * @param A An element of \f$\mathfrak{gl}(n)\f$
 * @param v An element of a representation of \f$\mathfrak{gl}(n)\f$
 * @return The element \f$ A\cdot v\f$
 */
	template<typename R> ex Action(ex A, ex v) const {
		ExVector comps=G->e().Components(A);
		ExVector eigoesto(G->n());
		for (int i=1;i<=G->n();++i)
		for (int j=1;j<=G->n();++j)
		{
			ex comp=G->GetDoubleIndexed(comps, i,j);
			eigoesto(j)+=comp*e(i);
		}
		LinearAction<T> a(e,eigoesto);
		return AlgebraAction<T,R>(a,v);
	}
/** @brief Computes the conditions on a generic element in the representation, for the algebra to act trivially on it.
 * @param container A container where the equation are to be stored
 * @param v A (generic) element of a representation of \f$\mathfrak{gl}(n)\f$
 * @return A reference to container
 */
	template<typename R, typename Container>
		Container& GetEquationsTrivialAction(Container& container, ex v) const
	{
		for (int i=1;i<=G->n();++i)
		for (int j=1;j<=G->n();++j)
			GetCoefficients<R>(container,Action<R>(G->A(i,j),v));
		return container;
	}
};




}
/** @} */

#endif /*REPSO_H_*/
