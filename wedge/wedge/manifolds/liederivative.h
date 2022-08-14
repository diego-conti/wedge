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
 *******************************************************************************/
#ifndef LIEDERIVATIVE_H
#define LIEDERIVATIVE_H
/** @ingroup Manifolds 
 *  @{ */
/**  @file liederivative.h 
 * @brief Lie derivatives of a differential form
 */
#include <list>
#include "wedge/base/wedgealgebraic.h"
#include "wedge/manifolds/manifold.h"

namespace Wedge {
 using namespace  GiNaC;
 using namespace std;

/** @brief The class representing the Lie derivative of a function w.r.t. a vector field
 * 
 */
class LieDerivative : public Register<LieDerivative,Function>::Named
{
	list<VectorField> X;
	Function f;
	static NameWithNoIndices CreateName(const list<VectorField>& v, const Function& u)
	{
		string name("("), texname;
		assert(!v.empty());
		for (list<VectorField>::const_iterator i=v.begin();i!=v.end();++i)
		{
			name+=i->get_name();
			texname+="\\mathcal{L}_{"+i->get_tex_name()+"}";
		}
		name+=u.get_name()+")";
		texname+=u.get_tex_name();
		return NameWithNoIndices(name,texname);
	}
	static NameWithNoIndices CreateName(const VectorField& v, const Function& u)
	{
		return NameWithNoIndices("("+v.get_name()+u.get_name()+")","\\mathcal{L}_{"+v.get_tex_name()+"}"+u.get_tex_name());
	}
public:
 /** @brief Create a Lie Derivative object, representing the function \f$\mathcal{L}_Xu\f$
  * @param X A vector field
  * @param u A function
  * 
  * Use this constructor to define a new Lie derivative object
 */	
	LieDerivative(const VectorField& X, const Function& u) : RegClass(CreateName(X,u)), f(u) {
		this->X.push_back(X);
	}

 /** @overload  
 */
	LieDerivative(const list<VectorField>& v, const Function& u) : RegClass(CreateName(v,u)), X(v),f(u)
	{
	}

 /** @overload  
 */
	LieDerivative() {}
	
	static const char* static_class_name() {return "LieDerivative";}
/** @brief Compute the Lie derivative of this LieDerivative object with respect to a vector field
 * @param Y A vector field on M
 * @param M A manifold
 * @return The covariant derivative \f$\mathcal{L}_Y f\f$, where \f$f\f$ stands for this LieDerivative object
 */
	ex Derive(const VectorField& Y, const Manifold& M) const
	{
		list<VectorField> v(X);
		v.push_front(Y);
		return Canonicalize(v,f,v.begin(),M);
	}

	bool is_equal_same_type(const GiNaC::basic& o) const
	{		
		return compare_same_type(o)==0; 
	}
	
	int compare_same_type(const GiNaC::basic& o) const
	{
		assert(is_a<LieDerivative>(o));
		const LieDerivative& other=static_cast<const LieDerivative&>(o);
		assert(!X.empty() && !other.X.empty());
		int cmp=f.compare(other.f);
		if (cmp!=0) return cmp;
		list<VectorField>::const_iterator i=X.begin(),j=other.X.begin();
		while (true) {
			cmp=i->compare_same_type(*j);
			if (cmp!=0) return cmp;
			++i; ++j;
			if (i==X.end())
			{
				if (j==other.X.end())
					return 0;
				else
					return -1;
			}
			else if (j==other.X.end())
				return 1;
		}
	}
protected:
	//override symbol::calchash so that the serial number is not used
	unsigned calchash() const { 
 		hashvalue = f.gethash();
		for (list<VectorField>::const_iterator i=X.begin();i!=X.end();++i)
		{
			hashvalue =  internal::rotate_left(hashvalue);
			hashvalue ^= i->gethash();
		}
		setflag(status_flags::hash_calculated);
		return hashvalue;
	}
	virtual void print(const print_context &c, unsigned level= 0) const {
		RegClass::print(c,level);
		if (dynamic_cast<const print_tree*>(&c)!=NULL)
		{
			for (list<VectorField>::const_iterator i=X.begin();i!=X.end();++i)
				i->print(c,level+1);
			f.print(c,level+1);
		}
	}
private:
/** @brief Construct a LieDerivative object from a sequence of vector fields
 * @param v A sequence of vector fields \f$X_1,\dotsc,X_n\f$ on a manifold \f$M\f$
 * @param f A function on \f$M\f$
 * @param at An iterator pointing to the unique element of v that may not be lesser than or equal to the following element
 * @param M The Manifold object
 * @return A sum of LieDerivative objects equivalent to \f$X_1\dotsb X_nf\f$, each of which is canonical.
 *
 * @internal There is no eval() memeber for LieDerivative's, as they are canonicalized on construction.
 */
	static ex Canonicalize(list<VectorField>& v, const Function& f, list<VectorField>::iterator at, const Manifold& M);
};


} /** @} */

#endif
