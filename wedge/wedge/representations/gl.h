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
 #ifndef GL_H_
#define GL_H_

#include "wedge/liealgebras/liegroup.h"
#include "wedge/liealgebras/liesubgroup.h"
namespace Wedge {

/** @ingroup Manifolds
 *  @{
 *  @file gl.h
 *  @brief The general linear group \f$ GL(n,\mathbb{R})\f$
*/

/** @brief The linear group \f$ GL(n,\mathbb{R})\f$
 */
class GL : public LieGroupWithoutParameters, public ConcreteManifold, public virtual Has_dTable
{
public:
/** @brief Define the linear group \f$ GL(n,\mathbb{R})\f$
 * 
 * Defines the structure constants according to \f$ de^{ij}=-\sum_k e^{ik}\wedge e^{kj} \f$
 */  
	GL(int n) : ConcreteManifold(CreateFrame(n))
	{
		_n=n;
		for (int i=1;i<=n;i++)
			for (int j=1;j<=n;j++)
			{
				ex r=0;
				for (int k=1;k<=n;k++)
					r-=A(i,k)*A(k,j);
				Declare_d(A(i,j),r);
			}
	}
			
/** @brief Return the element of the canonical basis with 1 at (i,j)
 */
	ex A(OneBased i, OneBased j) const
	{		
		return GetDoubleIndexed(e(),i,j);
	}	

/** @brief Return the element of the canonical basis with 1 at (i,j)
 *  @param i,j One-based indices
 * @param v The vector to be used as a basis
 *
 * The point of this function is that v can be any vector of the correct length
*/	
	ex GetDoubleIndexed(const ExVector& v, OneBased i, OneBased j) const
	{
		assert(i>0 && i<=_n);
		assert(j>0 && j<=_n);
		assert(v.size()==Dimension());
		return v((i-1)*_n+j);		
	}

	int n() const {return _n;}

	ex MatrixTo_gl(const matrix& m) const
	{
		if (n()!=m.cols() || n()!=m.rows()) throw InvalidArgument(__FILE__,__LINE__,m);
		ex result;
		for (int i=0;i<n();++i)
			for (int j=0;j<n();++j)
				result+=m(i,j)*A(i+1,j+1);
		return result;
	}

	matrix glToMatrix(ex x) const
	{
		matrix result(n(),n());
		for (int i=0;i<n();++i)
			for (int j=0;j<n();++j)
				result(i,j)=Hook(A(i+1,j+1),x);
		return result;
	}
private:
	int _n;
	exvector CreateFrame(int n)
	{
		exvector v;	v.reserve(n*n);
		for (int i=1;i<=n;i++)
			for (int j=1;j<=n;j++)
				v.push_back(DifferentialOneForm(N.e(i,j)));
		return v;
	}
};


/** @} */
}
#endif /*GL_H_*/
