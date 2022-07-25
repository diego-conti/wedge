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

/** @ingroup RiemannianGeometry */ 

/** @{ 
 * @file spinor.h
 * @brief Spinor algebra
 */
 
#ifndef SPINOR_H
#define SPINOR_H

#include "wedge/base/wedgealgebraic.h"
#include "wedge/base/expressions.h"

namespace Wedge {
using namespace GiNaC;
using namespace std;

/** @brief The class representing a spinor field on a Riemannian spin manifold
 * 
 * Since a Riemannian structure is represented by a choice of an orthonormal frame, spinors can be represented
 * accordingly in terms of a compatible basis of the complex spin representation. Objects of type Spinor represent
 * elements of this basis.
 * 
 * Spinors are implemented as strings of bits, but their %internal structure is not directly accessible. Objects of type spinor 
 * are constructed by calling RiemannianStrucure::u; the class RiemannianStructure also implements Clifford multiplication. 
 */
class Spinor : public Register<Spinor,Vector>::Algebraic
{
	friend class RiemannianStructure;
  	Spinor(unsigned index, int dimension);
	Spinor(vector<bool> index);				
public:
	int compare_same_type(const basic &other) const;	
	Spinor() {};
	unsigned return_type() const {return return_types::commutative;}
	static const char* static_class_name() {return "Spinor";}
private:
	void print(const print_context &c, unsigned level) const;	///< Overloaded from basic::print
	vector<bool> a;
};


} /** @} */

#endif
