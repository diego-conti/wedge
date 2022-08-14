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
#ifndef SO_H_
#define SO_H_

#include "wedge/liealgebras/liegroup.h"
#include "wedge/liealgebras/liesubgroup.h"
namespace Wedge {

/** @ingroup Manifolds
 *  @{
 *  @file so.h
 *  @brief The special orthogonal group \f$ SO(n,\mathbb{R})\f$
*/

/** @brief The special orthogonal group \f$ SO(n,\mathbb{R})\f$
 */
class SO : public LieGroupWithoutParameters, public ConcreteManifold, public virtual Has_dTable
{
public:
/** @brief Define the linear group \f$ SO(n,\mathbb{R})\f$
 * 
 * Defines the structure constants according to \f$ dA^{ij}=-\frac12\sum_k A^{ik}\wedge A^{kj} \f$
 */  
	SO(int n) : ConcreteManifold(CreateFrame(n))
	{
		n_=n;
		for (int i=1;i<=n;i++)
			for (int j=i+1;j<=n;j++)
			{
				ex r=0;
				for (int k=1;k<=n;k++)
					r-=A(i,k)*A(k,j);
				Declare_d(A(i,j),r);
			}
	}
			
/** @brief Return the element of the canonical basis with 1 at (i,j) and -1 at (j,i)
 */
	ex A(OneBased i, OneBased j) const
	{
		return GetDoubleIndexed(e(),i,j);
	}

/** @brief Return the element of the canonical basis with 1 at (i,j) and -1 at (j,i)
 *  @param i,j One-based indices
 * @param v The vector to be used as a basis
 *
 * The point of this function is that v can be any vector of the correct length
*/
	ex GetDoubleIndexed(const ExVector& v, OneBased i, OneBased j) const
	{
		assert(i>0 && i<=n_);
		assert(j>0 && j<=n_);
		assert(v.size()==Dimension());
		if (i<j)
			return v(j+(i-1)*n_-(1+i)*i/2);
		else if (j==i)
			return 0;
		else 
			return -v(i+(j-1)*n_-(1+j)*j/2);
	}

	int n() const {return n_;}
private:
	int n_;
	exvector CreateFrame(int n)
	{
		exvector v;	v.reserve((n*(n-1))/2);
		for (int i=1;i<=n;i++)
			for (int j=i+1;j<=n;j++)
				v.push_back(DifferentialOneForm(N.e(i,j)));
		return v;
	}
};


/** @} */
}
#endif /*SO_H_*/
