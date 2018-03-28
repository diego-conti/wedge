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
#ifndef TENSOR_H_
#define TENSOR_H_
/** @ingroup LinearAlgebra */ 

/** @{ 
 * @file wedge/tensor.h
 * @brief Tensor products
 */
 
#include <ginac/ex.h>
#include "bilinear.h"
#include "wedgealgebraic.h"

namespace Wedge
{
using namespace GiNaC;

/** @brief Simple elements of a tensor product \f$ V\otimes W\f$
 * @param V Type of simple elements in V
 * @param W Type of simple elements in W 
 * 
 * Instances of Tensor<V,W> represent objects \f$ v\otimes w\f$ where v is a simple element of V and w is a simple element of W.
 */
template<typename V,typename W> class Tensor : public Register<Tensor<V,W>,Vector>::Algebraic
{
public:
	ex v;
	ex w;
	static const char* static_class_name() {return "Tensor";}

	unsigned return_type() const {return return_types::noncommutative;}
#if RTTI_METHOD <= 1
 	tinfo_t return_type_tinfo() const {
 		return Registered<Tensor<V,W>,basic>::get_class_info_static().options.get_id();
	}
#endif 
	/** @brief Construct the tensor product of two simple elements */
	Tensor(ex _v, ex _w) : v(_v), w(_w) {}
	Tensor() {}

	int compare_same_type(const basic& other) const	
	{
		const Tensor<V,W>& o=static_cast<const Tensor<V,W>&>(other);
		int cmp=v.compare(o.v);
		if (cmp==0) return w.compare(o.w);
		else return cmp;
	}
/** @brief Overloaded substitution function 

@note There are limitations to using subs to change the type of an object. In fact if v has type V and w has type W, the substitution v==w will not map \f$v\otimes v\f$ to \f$w\otimes w\f$. However, the substitution TensorProduct<V,V>(v,v)==w will map \f$v\otimes v\f$ to \f$w\f$. A similar limitation holds for objects of type Lambda<V>.
*/
	virtual ex subs (const exmap &m, unsigned options=0) const;

protected:
	/** @brief Overloaded print function for nice \f$\text{\TeX}\f$ output
	 */	
	void print(const print_context &c, unsigned level = 0) const
	{
		if (dynamic_cast<const print_tree*>(&c)!=NULL) 
		{
			c.s << std::string(level, ' ') <<  this->class_name() << " @" << this
				<< std::hex << ", hash=0x" << this->hashvalue<< std::dec<<std::endl;
			v.print(c,level+1);
			w.print(c,level+1);
		}
		else 
			c.s<<v<<((dynamic_cast<const print_latex*>(&c)==NULL)? " X " : "\\otimes ")<<w;
	}
};

namespace internal {
	
template<typename V,typename W> class TensorProductOperator : public IBilinearOperator<LinearOperator<V>,LinearOperator<W> > {
public:
	ex Apply(const V& v, const W& w) const {
		return Tensor<V,W>(v,w);
	}
};

//helper class to expand selectively differential forms and the like
template<typename V> struct SelectiveExpand {
	static ex expand (ex v) {return v;}
};

template<typename V> struct SelectiveExpand<Lambda<V> > {
	static ex expand (ex v) {return v.expand();}
};

}

/** @brief Take the tensor product of a vector in V and a vector in W
 * @param v An element of a vector space whose simple elements have type V
 * @param w An element of a vector space whose simple elements have type W
 * @return The tensor product of v and w, as a linear combination of objects of type Tensor<V,W>
 */
template<typename V,typename W> ex TensorProduct(ex v, ex w)
{
	static internal::TensorProductOperator<V,W> t;	
	return internal::TensorProductOperator<V,W>::BilinearOperator(internal::SelectiveExpand<V>::expand(v),internal::SelectiveExpand<W>::expand(w),&t);
}


template<typename V,typename W> ex Tensor<V,W>::subs (const exmap &m, unsigned options) const 
	{
		exmap::const_iterator it;
		ex result=TensorProduct<V,W>(v.subs(m,options),w.subs(m,options));
		if (result!=*this) return result;
		if (options & subs_options::no_pattern) {
			it = m.find(result);
			if (it != m.end()) return it->second;
			else return result;
		}
		else 
			for (it = m.begin(); it != m.end(); ++it) {
#if (GINACLIB_MAJOR_VERSION<=1) && (GINACLIB_MINOR_VERSION<=4)
				lst repl_lst;
#else
				exmap repl_lst;
#endif
				if (this->match(ex_to<basic>(it->first), repl_lst))
					return it->second.subs(repl_lst, options | subs_options::no_pattern); 

			}
		return result;
	}

} /** @} */

#include "tensorlambda.h"

#endif /*TENSOR_H_*/
