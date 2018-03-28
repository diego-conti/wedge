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

#ifndef DIFFERENTIALFORM_H
#define DIFFERENTIALFORM_H
/** @ingroup Manifolds 
 *  @{ */
/**  @file differentialform.h 
 * @brief Vector fields, differential forms
 */
#include <ginac/basic.h>
#include <ginac/ex.h>
#include <string>
#include <cassert>
#include "lambda.h"
#include "linear.h"
#include "derivation.h"
#include "wedgealgebraic.h"

namespace Wedge {
 using namespace  GiNaC;
 using namespace std;

/** @brief The class representing a vector field on a manifold
 * 
 * In %Wedge, objects are ''identified with their duals''. This means that if X,Y,Z are instances of VectorField, then the same objects
 * X,Y,Z can also represent the dual basis of 1-forms.
 *
 * @note The class VectorField is also used to represent differential one-forms. However, the product of vectorField's does not represent
 * a differential form of higher degree (in fact, this product is commutative). \sa DifferentialOneForm 
 */

class VectorField : public Register<VectorField,Vector>::NamedNumbered
{
public:
 /** @brief Create a @em simple vector field
  * @param name The name of this vector field (only used in output).
  * 
  * Use this constructor to define a new vector field.
  * @note Distinct instances of VectorField are regarded as linearly independent.  
 */	
	VectorField(const Name& name);

 /** @overload  
 */
	VectorField();
	
	static const char* static_class_name() {return "VectorField";}
};

/** @brief Type of a differential one-form on a manifold; derived from VectorField.
 * 
 * Differs from VectorField in that the product of differentialOneForm's is skew-symmetric.
 * Differential forms of higher degree are represented as (wedge) products of differentialOneForm's. \sa DifferentialForm
 * 
 * @note This type implements the wedge product. The exterior derivative operator is implemented in Wedge::Manifold and its derived classes
 * 
 * @note  Wedge identifies differential one-forms with vector fields; hence, functions 
 * that take a differential one-form as a parameter should declare that parameter as a VectorField.
*/
typedef Lambda1<VectorField> DifferentialOneForm;

/** @brief Type representing a differential form on a manifold
 * 
 * The wedge product of differentialOneForm's is a differentialForm.  
*/
typedef Lambda<VectorField> DifferentialForm;


/** @brief The interior product, or "hook" operator
  * @param left A differential form (viewed as a multi-vector field)
  * @param right A differential form
  * @returns The interior product of left into right
  * 
  * @note Because of the identification of vector fields with differential one-forms, the result does not
  * depend on the choice of a metric. Equivalently, the interior product is computed with respect
  * to the metric for which the simple elements are an orthonormal basis.  @sa RiemannianStructure::Hook
 */

ex Hook(ex left, ex right);

} /** @} */

#endif
