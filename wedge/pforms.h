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
#ifndef PFORMS_H_
#define PFORMS_H_

#include "lambda.h"
#include "utilities.h"
/** @ingroup LinearAlgebra */
/** @{ 
 * @file pforms.h
 * @brief Compute exterior algebra over a given vector space
 */
namespace Wedge { 

namespace internal {
	template<typename V> class PFormsHelper : public IterateOverSubsets {
		const IBasis<Lambda1<V> >& space;
		list<ex> result;
		int k;
	public:
		PFormsHelper(const IBasis<Lambda1<V> >& _space) :
			space(_space) {}
		list<ex> GetGenerators(unsigned degree) {
			k=0;
			Iterate(degree,space.size());
			return result;
		}
	protected:
		virtual bool Apply(const vector<int>& form) {
			vector<int>::const_iterator i=form.begin();
			ex alpha=space[*i++];
			while (i!=form.end())
				alpha*=space[*i++];
			result.push_back(alpha.expand());
			return true;
		}
	};
}

/** @brief Compute the exterior algebra over a given vector space 
 * @param space A basis of a vector space \f$W\subset\Lambda^1(V)\f$
 * @param p The degree
 * @return The vector space \f$\Lambda^p(W)\f$
 */
template<typename V> VectorSpace<Lambda<V> > pForms(const IBasis<Lambda1<V> >& space, int p)
{
	if (p<0) throw OutOfRange(__FILE__,__LINE__,p);
	else if (p==0) 
	{
		lst l(1);
		return VectorSpace<Lambda<V> > (l.begin(),l.end());		
	}	
	else if (p>space.size())
	{
		return VectorSpace<Lambda<V> >();
	}
	
	internal::PFormsHelper<V> h(space);
	list<ex> gens=h.GetGenerators(p);
	return VectorSpace<Lambda<V> >(gens.begin(),gens.end());
}

/** @brief Compute the degree 2 exterior algebra over a subspace of some vector space
 * @param subspace A subspace  \f$W\subset W'\subset\Lambda^1(V)\f$ 
 * @return The subspace \f$\Lambda^2(W)\subset\Lambda^2(W')\f$  
 */
template<typename V> Subspace<Lambda<V> > TwoForms(const SubBasis<Lambda1<V> >& subspace)
{
	exvector Lambda2_h;
	for (typename SubBasis<Lambda1<V> >::const_iterator i=subspace.begin();i!=subspace.end();i++)
		for (typename SubBasis<Lambda1<V> >::const_iterator j=i+1;j!=subspace.end();j++)
		{
			Lambda2_h.push_back(*i * *j);			
		}
		
	exvector complement ;
	for (typename SubBasis<Lambda1<V> >::const_iterator i=subspace.full_begin();i!=subspace.full_end();i++)
		for (typename SubBasis<Lambda1<V> >::const_iterator j=subspace.complement_begin();j!=subspace.complement_end();j++)
		{
			if (i<j) complement.push_back(*i * *j);			
		}
	SubBasis<Lambda<V> > basis (Lambda2_h.begin(),Lambda2_h.end(),complement.begin(),complement.end());
	return 	basis;
}

}  /** @} */
#endif /*PFORMS_H_*/
