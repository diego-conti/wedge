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
#ifndef FUNCTION_H_
#define FUNCTION_H_
/** @ingroup Manifolds  
*  @{ */

/**  @file wedge/function.h
 * @brief Functions on manifolds
 */
#include <ginac/basic.h>
#include "wedgealgebraic.h"
#include "differentialform.h"

namespace Wedge 
{
using namespace GiNaC;
class Manifold;

/** @brief Class representing a function on a manifold
 * 
 * A function is treated essentially a zero-degree form, but the types Function and DifferentialForm
 * are not related. Polymorphic mechanisms must be employed by functions that take a differential form
 * as a parameter to ensure that objects of type Function are handled correctly.
 *  
 * @note Function is derived from symbol because coordinates on a manifold are functions. This enables one
 * to take the derivative with respect to a coordinate (since GiNaC can only take derivatives with respect
 * to objects whose type is either GiNaC::symbol or a derived class).
 */

class Function : public Register<Function,realsymbol>::Algebraic
{
public:
	static const char* static_class_name() {return "Function";}
 /** @brief Create a function object
  * @param name The name of this function (only used in output).
  * 
  * Use this constructor to define a new function
  * @note Distinct instances of Function are regarded as independent.  
 */	
	Function(const Name& name) : RegClass(name) {;}	
 /** @overload  
 */
	Function() : RegClass(N.f(serial)) {}
 /** @brief Compute a Lie derivative of this function \f$f\f$
 * @param Y A vector field on a manifold \f$M\f$
 * @param M A Manifold object
 * @return The Lie derivative \f$\mathcal{L}_Y f\f$
 */
	virtual ex Derive(const VectorField& Y, const Manifold& M) const;

/* @brief Compute the derivative of this function \f$f\f$ with respect to a (coordinate) function
 * @param t A (coordinate) function on a manifold \f$M\f$
 * @param M A Manifold object
 * @return The derivative \f$\frac{\partial f}{\partial t}\f$
 *
	virtual ex Derive(const Function& t, const Manifold& M) const;
	*/
};

} /** @} */

#endif /*FUNCTION_H_*/
