/***************************************************************************
 *   Copyright (C) 2008, 2009 by Diego Conti				   *
 *   diego.conti@unimib.it                                                 *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 3 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#ifndef REPSL2_H_
#define REPSL2_H_
/** @ingroup Representations */
/** @{ 
  * @file repsl2.h 
  * @brief Explicit representations of the Lie algebras \f$\mathfrak{sl}(2,\mathbb{C})\f$ and \f$\mathfrak{su}(2)\f$
*/


#include "linearaction.h"

namespace Wedge {

/** @brief A real irreducible representation of \f$SO(3)\f$

We view \f$so(3)\f$ in terms of the generators satisfying the relations
\f[[H,X]=2Y, \quad [H,Y]=-2X, \quad [X,Y]=2H.\f]
*/
template<typename T> class SO3Representation {
	LinearAction<T> aH, aX,aY;
	struct Init {
		exvector e;
		exvector from, x_to, y_to,h_to;
		int n;
		ex w(ZeroBased k) const
		{
			assert(k<=n/2+1);
			if (k<0) return 0;
			if (k==n/2) {
				if (n%4==0) 	// if n/2 is even..
					return e[n]*sqrt(ex(2));	// ... w_k is non-zero and the last element
					else
					return 0;
			}
				else if (k==n/2+1)
				return (n%4==0) ? -w(n/2-1)*(n/2+1)*n/2 : w(n/2-1)*(n/2+1)*n/2;
			else
				return e[k]*NormalizeRoots(sqrt(factorial(k)/factorial(n-k)));
		}
		ex u(ZeroBased k) const
		{		
			assert(k<=n/2+1);
			if (k<0) return 0;
			if (k==n/2) {
				if (n%4!=0) 	// if n/2 is odd..
					return e[n]*sqrt(ex(2));	// ... u_k is non-zero and the last element
				else
					return 0;
			}
			else if (k==n/2+1)
				return (n%4==0) ? u(n/2-1)*(n/2+1)*n/2 : -u(n/2-1)*(n/2+1)*n/2;
			else
				return e[k+n/2]*NormalizeRoots(sqrt(factorial(k)/factorial(n-k)));
		}
		Init(const exvector& frame) : e(frame) {
			n=frame.size()-1;
			if (n%2!=0) throw InvalidArgument(__FILE__,__LINE__,n);
			for (int k=0;k<n/2;++k)
			{	
				from.push_back(w(k));
				from.push_back(u(k));
				x_to.push_back(k*(n-k+1)*w(k-1)-w(k+1));
				x_to.push_back(k*(n-k+1)*u(k-1)-u(k+1));
				y_to.push_back(k*(n-k+1)*u(k-1)+u(k+1));
				y_to.push_back(-k*(n-k+1)*w(k-1)-w(k+1));
				h_to.push_back((n-2*k)*u(k));
					h_to.push_back((2*k-n)*w(k));
			}
			int k=n/2;
			if (u(k).is_zero())
			{
				from.push_back(w(k));
				x_to.push_back(k*(n-k+1)*w(k-1)-w(k+1));
				y_to.push_back(k*(n-k+1)*u(k-1)+u(k+1));
				h_to.push_back((n-2*k)*u(k));
				}
			else
			{
				from.push_back(u(k));
				x_to.push_back(k*(n-k+1)*u(k-1)-u(k+1));
				y_to.push_back(-k*(n-k+1)*w(k-1)-w(k+1));
				h_to.push_back((2*k-n)*w(k));
			}
		}
	
	};
public:
	typedef T ActsOnType;	
/** @brief Construct the irreducible real representation of \f$SO(3)\f$ corresponding to the given choice of basis
 *
 * @param frame A vector of odd length of vectors of type T.
*/
	SO3Representation(const exvector& frame) {
		Init init(frame);		
		aH.InitializeFromBasis (init.from,init.h_to);
		aX.InitializeFromBasis(init.from,init.x_to);
		aY.InitializeFromBasis(init.from,init.y_to);
	}
	template<typename R> ex H(ex v) const {
		return AlgebraAction<T,R>(aH,v);
	}
	template<typename R> ex X(ex v) const {
		return AlgebraAction<T,R>(aX,v);
	}
	template<typename R> ex Y(ex v) const {
		return AlgebraAction<T,R>(aY,v);
	}

/** @brief Computes the conditions on a generic element in the representation, for the algebra to act trivially on it.
 * @param container A container where the equation are to be stored
 * @param v A (generic) element of a representation of \f$\mathfrak{so}(3)\f$
 * @return A reference to container
 */
	template<typename R, typename Container>
		Container& GetEquationsTrivialAction(Container& container, ex element) const
	{
		GetCoefficients<R>(container,H<R>(element));
		GetCoefficients<R>(container,X<R>(element));
		GetCoefficients<R>(container,Y<R>(element));
		return container;
	}
};


/** @brief An irreducible representation of \f$SL(2,\mathbb{C})\f$

We view \f$sl(2,\mathbb{C})\f$ in terms of the generators satisfying the relations
\f[[H,X]=2X, [H,Y]=-2Y, [X,Y]=H.\f]

@sa [Fulton-Harris, Representation Theory]
*/
template<typename T> class SL2Representation  {
	LinearAction<T> aH, aX,aY;
public:	
	typedef T ActsOnType;
/** @brief The representation \f$Sym(n,\mathbb{C}^2)\f$
*/
	SL2Representation(const ExVector& e) {
		int n=e.size();
		lst subs_H, subs_X, subs_Y;
		for (int i=1;i<=n;++i)
			subs_Y.append(e(i)==e(i+1));
		subs_Y.append(e(n+1)==0);
		for (int i=2;i<=n+1;++i)
			subs_X.append(e(i)==(i-1)*(n-i+2)*e(i-1));
		subs_X.append(e(1)==0);
		for (int i=1;i<=n+1;++i)
			subs_H.append(e(i)==(n-2*(i-1))*e(i));
		aH.set(subs_H);
		aX.set(subs_X);
		aY.set(subs_Y);
	}
	template<typename R> ex H(ex form) const {
		return AlgebraAction<T,R>(aH,form);
	}
	template<typename R> ex X(ex form) const {
		return AlgebraAction<T,R>(aX,form);
	}
	template<typename R> ex Y(ex form) const {
		return AlgebraAction<T,R>(aY,form);
	}
/** @brief Computes the conditions on a generic element in the representation, for the algebra to act trivially on it.
 * @param container A container where the equation are to be stored
 * @param v A (generic) element of a representation of \f$\mathfrak{sl}(2,\mathbb{C})\f$
 * @return A reference to container
 */
	template<typename R, typename Container>
		Container& GetEquationsTrivialAction(Container& container, ex element) const
	{
		GetCoefficients<R>(container,H<R>(element));
		GetCoefficients<R>(container,X<R>(element));
		GetCoefficients<R>(container,Y<R>(element));
	}
};


}
/** @} */

#endif /*REPSL2_H_*/
