/*******************************************************************************
 *  Copyright (C) 2007-2022 by Diego Conti, diego.conti@unimib.it 
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
 *******************************************************************************//** @ingroup RiemannianGeometry */ 

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

/** @brief The class representing a spinor field on a pseudo-Riemannian spin manifold
 * 
 * This class works with pseudo-Riemanannian structures defined in terms of an orthonormal frame, @sa PseudoRiemannianStructureByOrthonormalFrame
 * Spinors can be represented accordingly in terms of a compatible basis of the complex spin representation. Objects of type Spinor represent
 * elements of this basis.
 * 
 * Spinor objects should not be manipulated directly. Objects of type spinor should be constructed by calling RiemannianStrucure::u or PseudoRiemannianStructureByOrthonormalFrame::u.
 * The classes RiemannianStructure and PseudoRiemannianStructureByOrthonormalFrame also implement Clifford multiplication. 
 */
class Spinor : public Register<Spinor,Vector>::Algebraic
{
	friend class RiemannianStructure;
  	Spinor(unsigned index, int dimension);
	Spinor(vector<bool> index);				
public:
	int compare_same_type(const basic &other) const;	
	Spinor()=default;

/** @brief Construct a spinor
 * @param signs the sequence (\epsilon_1,\dotsc, \epsilon_m)
 * @returns the spinor u(\epsilon_m,\dotsc, \epsilon_1) in the notation of 
 * [Baum, H. and Kath, I. Parallel Spinors and Holonomy Groups onPseudo-Riemannian Spin Manifolds. Annals of Global Analysis and Geometry 17: 1–17, 1999.]
 * 
 * Note the inversion in the order.
 */
	static Spinor from_epsilons(const vector<int>& signs) {
		vector<bool> as_bools;
		transform(signs.begin(),signs.end(), back_inserter(as_bools), [] (int x) {return x<0;});
		return Spinor{as_bools};
	}

	unsigned return_type() const {return return_types::commutative;}
	static const char* static_class_name() {return "Spinor";}
/** @brief Invert an index
 * @param j an index in the interval [1,m]
 *  @return the spinor u(\epsilon_m,... , -\epsilon_j, ..., \epsilon_1)
 */
	Spinor reflect(OneBased j) const {
		vector<bool> reversed=a;
		reversed[j-1]=!reversed[j-1];
		return Spinor{reversed};
	}
/** @brief Return the product of the first j signs
 * @param j an index in the interval [1,m]
 *  @return the product \epsilon_1,..., \epsilon_j
 */
	int product_up_to(OneBased j) const {
		 return std::count(a.begin(),a.begin()+j,true)%2==0? 1: -1;
	}
private:
	void print(const print_context &c, unsigned level) const;	///< Overloaded from basic::print
	vector<bool> a;	//true represents -1, false 1
};


} /** @} */

#endif
