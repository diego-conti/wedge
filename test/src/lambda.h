/*******************************************************************************
 *  Copyright (C) 2007-2023 by Diego Conti, diego.conti@unipi.it 
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
#ifndef TEST_LAMBDA_H_
#define TEST_LAMBDA_H_

#include "wedge/base/wedgebase.h"
#include "wedge/base/wedgealgebraic.h"
#include "wedge/base/expressions.h"
#include "wedge/base/normalform.h"
#include <cxxtest/TestSuite.h>
#include "test.h"
#include "wedge/linearalgebra/lambda.h"
#include "wedge/convenience/canonicalprint.h"
WEDGE_DECLARE_NAMED_ALGEBRAIC(V,Vector)
WEDGE_DECLARE_NAMED_ALGEBRAIC(W,Vector)

using namespace Wedge;
using namespace GiNaC;
/**
	@author Diego Conti <diego.conti@unipi.it>
	@brief Test suite for the Lambda class
		
*/

class LambdaTestSuite : public CxxTest::TestSuite  {
public:
	void testNormalForm()
	{
		V x(N.x),y(N.y),z(N.z);
		W a(N.a),b(N.b),c(N.c);
		exvector p;
		p.push_back((x+y)*(1+b)*(a+c));
		p.push_back((x+y)*(a+b)*(1+c)+z*c);
		p.push_back(((x-y)*(a+b)*(1+c)+z*c).expand());
		p.push_back(((x+y)*(a+b)*(1+c)+z*c).expand());
		exvector q;
		q.push_back((x+y)*(a+c+a*b+c*b));
		q.push_back((x+y)*(a+b+a*c+b*c)+z*c);
		q.push_back((x-y)*(a+b+a*c+b*c)+z*c);
		q.push_back((x+y)*(a+b+a*c+b*c)+z*c);
		TS_ASSERT_EQUALS(NormalForm<V>(p[0]),q[0]);	
		TS_ASSERT_EQUALS(NormalForm<V>(p[1]),q[1]);
		TS_ASSERT_EQUALS(NormalForm<V>(p[2]),q[2]);
		TS_ASSERT_EQUALS(NormalForm<V>(p[3]),q[3]);
		exvector r=NormalForm<V>(p);
		TS_ASSERT_EQUALS(q,r);		

		symbol t;				
		TS_ASSERT_EQUALS(NormalForm<Lambda<V> >(x),x);
		TS_ASSERT_EQUALS(NormalForm<Lambda<V> >(x*y+x*z),x*y+x*z);
		TS_ASSERT_EQUALS(NormalForm<Lambda<V> >(x*y+z+1),x*y+z+1);
		Lambda1<V> x1(N.x), y1(N.y); 
		TS_ASSERT_EQUALS(NormalForm<Lambda<V> >(x1*y1),x1*y1);		
		TS_ASSERT_EQUALS(NormalForm<Lambda<V> >(x1*y1+sin(t)*x1*y1),(1+sin(t))*x1*y1);

	}

	void testVectorNormalForm()
	{
		V x(N.x),y(N.y),z(N.z);
		W a(N.a),b(N.b),c(N.c);

		LambdaVectorNormalForm v{-x+2*y}, w{-z+2*b+sin(ex(3))};
		TS_ASSERT_EQUALS(v.scalar_part,0);
		TS_ASSERT_EQUALS(w.scalar_part,sin(ex(3)));

		TS_ASSERT_THROWS(LambdaVectorNormalForm{a*b},WedgeException<std::runtime_error>);
		TS_ASSERT_THROWS(LambdaVectorNormalForm{x+a*b},WedgeException<std::runtime_error>);
		TS_ASSERT_THROWS(LambdaVectorNormalForm{x+x*b},WedgeException<std::runtime_error>);
		TS_ASSERT_THROWS(LambdaVectorNormalForm{2*x*b},WedgeException<std::runtime_error>);
	}
	void testCanonicalPrint() {
		V x(N.x),y(N.y),z(N.z(1));
		symbol a("a1","a_1"),b("b");
		ex linear(a*x+a*b*b*y+(a+b)*z);
		TS_ASSERT_EQUALS(to_canonical_string(linear),"(a1+b)*z1+a1*b^2*y+a1*x");
		TS_ASSERT_EQUALS(to_latex_canonical_string(linear),"(a_1+b) z_1+a_1 b^2 y+a_1 x");
	}
	void testCanonicalPrintMinus() {
		V x(N.x),y(N.y),z(N.z(1));
		symbol a("a1","a_1"),b("b");
		ex linear(-a*x-y-(a+b)*z);
		TS_ASSERT_EQUALS(to_canonical_string(linear),"-y-(a1+b)*z1-a1*x");
		TS_ASSERT_EQUALS(to_latex_canonical_string(linear),"-y-(a_1+b) z_1-a_1 x");
	}
	void testLambdaVectorNormalForm()
	{
		ex sqrt3=sqrt(ex(3));
		symbol t;
		Lambda1<V> x(N.x),y(N.y),z(N.z);
		Lambda1<W> a(N.a),b(N.b),c(N.c);
		exvector p;
		p.push_back(x*(z-y));
		p.push_back(x*(y+z)-a);
		p.push_back(x*(sqrt3*y+z)-a);
		p.push_back(-2*exp(t)*a*b+3*exp(-t)*x+z*(x-y));
		exvector q;
		q.push_back(x*z-x*y);
		q.push_back(x*y+x*z-a);
		q.push_back(sqrt3*x*y+x*z-a);
		q.push_back(-2*exp(t)*a*b+3*exp(-t)*x+z*x-z*y);
		for (auto j: p)
		for (auto i: VectorNormalForm {j})
			TS_ASSERT(is_a<VectorBase>(i.first));

		for (int i=0;i<p.size();++i)
			TS_ASSERT_EQUALS(VectorNormalForm(p[i]).CollectCoefficients(),q[i]);
		for (auto x: p)
			TS_ASSERT_EQUALS(VectorNormalForm(x).CollectCoefficients(),x);

	}

	
	void testLambda() {
		Lambda1<V> x,y,z;
		set<ex,ex_is_less> simple;
		GetSimple<Lambda<V> > (simple,3*x*y+2*z);
		TS_ASSERT(simple.find(x*y)!=simple.end());
		TS_ASSERT(simple.find(z)!=simple.end());
		TS_ASSERT_EQUALS(simple.size(),2);
		
		set<ex,ex_is_less> coeffs;		
		GetCoefficients<Lambda<V> >(coeffs, 3*x*y+2*z);
		TS_ASSERT(coeffs.find(3)!=coeffs.end());
		TS_ASSERT(coeffs.find(2)!=coeffs.end());
		TS_ASSERT_EQUALS(coeffs.size(),2);

		V t;
		coeffs.clear();
		GetSymbols<Lambda<V> > (coeffs,3*x*y+2*z+t);
		TS_ASSERT(coeffs.find(x*y)!=coeffs.end());
		TS_ASSERT(coeffs.find(z)!=coeffs.end());
		TS_ASSERT(coeffs.find(t)!=coeffs.end());
	}

};

#endif
