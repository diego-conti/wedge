/*******************************************************************************
 *  Copyright (C) 2007-2022 by Diego Conti, diego.conti@unipi.it 
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
#ifndef SU_H_
#define SU_H_

#include "wedge/liealgebras/liegroup.h"
#include "wedge/liealgebras/liesubgroup.h"
namespace Wedge {

/** @ingroup Manifolds
 *  @{
 *  @file su.h
 *  @brief The special unitary group \f$ SU(n)\f$
*/

/**  @brief The special unitary group \f$ SU(n)\f$
 */
class SU : public LieGroupWithoutParameters, public ConcreteManifold, public virtual Has_dTable
{
public:
/** @brief Define the special unitary group \f$ SU(n)\f$
 *
 * Uses the basis \f$\{E_{ij}, F_{ij}, G_i\}\f$ related to the standard basis \f$e_{ij}\f$ of \f$\mathfrak{gl}(n,\mathbb{C})\f$ by
 * \f[E_{ij}=e_{ij}-e_{ji},\quad F_{ij}=i(e_{ij}+e_{ji}), \quad G_{i}=i(e_{ii}-e_{i+1,i+1}).\f]
 * The structure constants are given by 
\f[
dE^{pq}=-E^{pi}\wedge E^{iq} + \sum_{k\neq p,q} F^{pk}\wedge F^{kq} - F^{pq}\wedge (G^p-G^{p+1}-G^{q}+G^{q-1}),\quad
dF^{pq}=-\sum_{r\neq p,q} (E^{pr}\wedge F^{rq} +E^{qr}\wedge F^{rp}) + E^{pq}\wedge (G^p-G^{p+1}-G^{q}+G^{q-1}),\quad
dG^p=\sum_{1\leq m\leq n,1\leq r\leq p} E^{mr}\wedge F^{rm}-\sum_{1\leq m\leq n,p<r\leq n} E^{mr}\wedge F^{rm}.
\f]
 */  
	SU(int n) : ConcreteManifold(n*n-1)
	{
		_n=n;
		for (int p=1;p<=n;++p) 			
			for (int q=p+1;q<=n;++q)
			{
				ex r;
				for (int k=1;k<=n;k++)
				{
					r-=E(p,k)*E(k,q);
					r+=F(p,k)*F(k,q);	//notice that F(p,p)=F(q,q)=0
				}
				r-=F(p,q)*(G(p)-G(p-1)-G(q)+G(q-1));
				Declare_d(E(p,q),r);

				r=0;
				for (int k=1;k<=n;k++)
				{
					r-=E(p,k)*F(k,q);
					r-=E(q,k)*F(k,p);	//notice that F(p,p)=F(q,q)=0
				}
				r+=E(p,q)*(G(p)-G(p-1)-G(q)+G(q-1));
				Declare_d(F(p,q),r);
			}
		for (int i=1;i<n;++i)
		{
			ex r;
			for (int m=1;m<=n;++m)  
			{
				for (int k=1;k<=i;++k)
					r+=E(m,k)*F(k,m);
				for (int k=i+1;k<=n;++k)
					r-=E(m,k)*F(k,m);
			}
			Declare_d(G(i), r);
		}	
		Check_ddZero();
	}
			
/** @brief Return the element E(i,j)
 */
	ex E(OneBased i, OneBased j) const
	{
		if (i<j) return e(DoubleIndexToOneBased(i,j));
		else if (j<i) return -e(DoubleIndexToOneBased(j,i)); 
		else return 0;
	}

/** @brief Return the element F(i,j)
 */
	ex F(OneBased i, OneBased j) const
	{		
		if (i<j) return e(DoubleIndexToOneBased(i,j)+(_n*(_n-1))/2);
		else if (j<i) return e(DoubleIndexToOneBased(j,i)+(_n*(_n-1))/2); 
		else return 0;
	}
/** @brief Return the element G(i)
 */
	ex G(OneBased i) const
	{
		assert(i>=0 && i<=_n);
		if (i==0 || i==_n) return 0;
		else return e(_n*(_n-1)+i);
	}	

/** @brief Convert a double skew symmetric index to a single index
 *  @param i,j One-based indices with i<j.
*/	
	int DoubleIndexToOneBased(OneBased i, OneBased j) const
	{
		assert(i>0 && i<=_n);
		assert(i<j && j<=_n);
		return (j+(i-1)*_n-(1+i)*i/2);
	}

	int n() const {return _n;} ///< Returns the value of n, for an object representing the group \f$SU(n)\f$ 
private:
	int _n;
	exvector CreateFrame(int n)
	{
		exvector v;
		v.reserve(n*n-1);
		for (int i=1;i<=n;i++)
			for (int j=i+1;j<=n;j++)
				v.push_back(DifferentialOneForm(N.E(i,j)));
		for (int i=1;i<=n;i++)
			for (int j=i+1;j<=n;j++)
				v.push_back(DifferentialOneForm(N.F(i,j)));
		for (int i=1;i<n;i++)
			v.push_back(DifferentialOneForm(N.G(i)));
		return v;
	}
};


/** @} */
}
#endif /*SU_H_*/
