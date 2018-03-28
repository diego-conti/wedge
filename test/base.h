/***************************************************************************
 *   Copyright (C) 2007, 2008, 2009 by Diego Conti			   *
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
#ifndef BASE_H_
#define BASE_H_

#include <list>
#include <vector>
#include <cxxtest/TestSuite.h>
#include <ginac/ginac.h>
#include "test.h"
#include "wedge/vectorspace.h"
#include "wedge/utilities.h"
#include "wedge/differentialform.h"
#include "wedge/tensor.h"
#include "wedge/lambda.h"
#include "wedge/wedgealgebraic.h"
#include "wedge/parameters.h"
#include "wedge/concretemanifold.h"
#include "wedge/normalform.h"


#ifdef HAVE_ARAGELI
#include "wedge/aragelilinalg.h"
#endif

using namespace Wedge;
using namespace GiNaC;
/**
	@author Diego Conti <diego.conti@unimib.it>
	@brief Test suite for the Base module in %Wedge
		
*/


WEDGE_DECLARE_NAMED_ALGEBRAIC(A,symbol);
WEDGE_DECLARE_ALGEBRAIC(B,A)
WEDGE_DECLARE_NAMED_ALGEBRAIC(C,realsymbol);


//test expressions.h
class BaseTestSuite : public CxxTest::TestSuite 
{
public:
	void testGetSymbols()
	{
		A a1(N.a(1)),a2(N.a(2));
		B b1,b2;
		list<ex> symbols;
		GetSymbols<A>(symbols,a1*a1*b1+log(a2)*pow(b2,3));
		list<ex>::const_iterator i=find(symbols.begin(),symbols.end(),b1);
		TS_ASSERT_DIFFERS(i,symbols.end());
		TS_ASSERT(is_a<B>(*i));
		i=find(symbols.begin(),symbols.end(),b2);
		TS_ASSERT_DIFFERS(i,symbols.end());
		TS_ASSERT(is_a<B>(*i));
		i=find(symbols.begin(),symbols.end(),a1);
		TS_ASSERT_DIFFERS(i,symbols.end());
		TS_ASSERT(is_a<A>(*i) && !is_a<B>(*i));
		TS_ASSERT(is_exactly_a<A>(*i));
		i=find(symbols.begin(),symbols.end(),a2);
		TS_ASSERT_DIFFERS(i,symbols.end());
		TS_ASSERT(is_a<A>(*i) && !is_a<B>(*i));
		TS_ASSERT(is_exactly_a<A>(*i));

		C c1(N.c(1));
		GetSymbols<C>(symbols, c1*a1);
		i=find(symbols.begin(),symbols.end(),c1);
		TS_ASSERT_DIFFERS(i,symbols.end());
		symbols.clear();
		GetSymbols<realsymbol>(symbols, c1*a1);
		i=find(symbols.begin(),symbols.end(),c1);
		TS_ASSERT_DIFFERS(i,symbols.end());
	}
	void testGetSymbols_old()	//this was written with the old version...
	{
		vector<V> v;
		v.push_back(V(N.v(0)));v.push_back(V(N.v(1)));v.push_back(V(N.v(2)));
		vector<W> w;
		w.push_back(W(N.w(0)));w.push_back(W(N.w(1)));w.push_back(W(N.w(2)));
		symbol s1,s2,s3;
		exvector a;
		a.push_back(v[0]*sin(w[0])*s1+v[1]*pow(2,w[1])*s2+v[2]*w[2]*s3);
		a.push_back(v[0]*v[1]*s1+w[0]*w[1]*s2);
			
		vector<V> sV;		
		GetSymbols<V>(sV,a.begin(),a.end());
		LOG_INFO(a);LOG_INFO(sV);
		TS_ASSERT(includes(sV.begin(),sV.end(),v.begin(),v.end()));
		GetSymbols<V>(sV,a[0]);
		LOG_INFO(a);LOG_INFO(sV);
		TS_ASSERT(includes(sV.begin(),sV.end(),v.begin(),v.end()));
		GetSymbols<V>(sV,a[1]);
		LOG_INFO(a);LOG_INFO(sV);
		TS_ASSERT(includes(sV.begin(),sV.end(),v.begin(),--v.end()));
		
		list<W> sW;
		GetSymbols<W>(sW,a.begin(),a.end());
		TS_ASSERT(includes(sW.begin(),sW.end(),w.begin(),w.end()));
		GetSymbols<W>(sW,a[0]);
		TS_ASSERT(includes(sW.begin(),sW.end(),w.begin(),w.end()));
		GetSymbols<W>(sW,a[1]);
		TS_ASSERT(includes(sW.begin(),sW.end(),w.begin(),--w.end()));
					
		list<ex> ssw;
		GetSymbols<Lambda<W> >(ssw,a.begin(),a.end());
		TS_ASSERT(includes(ssw.begin(),ssw.end(),w.begin(),--w.end()));
	}
	
	void testGetCoefficients()
	{
		vector<V> v;
		v.push_back(V(N.v(0)));v.push_back(V(N.v(1)));v.push_back(V(N.v(2)));
		vector<W> w;
		w.push_back(W(N.w(0)));w.push_back(W(N.w(1)));w.push_back(W(N.w(2)));
		realsymbol s1,s2,s3;
		exvector a;
		a.push_back(v[0]*sin(w[0])*s1+v[1]*pow(2,w[1])*s2+v[2]*w[2]*s3);
		a.push_back(v[0]*w[1]*s1+w[0]*w[1]*s2);
		
		exvector sV;
		TS_ASSERT_THROWS(GetSimple<V>(sV,v[0]*v[1]*s1+w[0]*w[1]*s2),WedgeException<std::runtime_error>);
		TS_ASSERT_THROWS(GetSimple<W>(sV,v[0]*v[1]*s1+w[0]*w[1]*s2),WedgeException<std::runtime_error>);
		
		GetSimple<V>(sV,a.begin(),a.end());
		TS_ASSERT(includes(sV.begin(),sV.end(),v.begin(),v.end()));
		GetSimple<V>(sV,a[0]);
		TS_ASSERT(includes(sV.begin(),sV.end(),v.begin(),v.end()));
		TS_ASSERT_THROWS(GetCoefficients<V> (sV,a),WedgeException<std::runtime_error>);
		
		list<ex> l;
		GetCoefficients<V> (l,a[0],withoutRHS);
		exvector eqns;
		eqns.push_back(sin(w[0])*s1);
		eqns.push_back(pow(2,w[1])*s2);
		eqns.push_back(s3);
		TS_ASSERT(includes(l.begin(),l.end(),eqns.begin(),eqns.end()));				
		
		l.clear();
		GetCoefficientsComplex<V> (l,a[0]+I*w[1]*s2*v[1],withoutRHS);
		eqns.push_back(w[1]*s2);
		TS_ASSERT(includes(l.begin(),l.end(),eqns.begin(),eqns.end()));
		
		l.clear();
		GetCoefficients<V> (l,a[0]+I*w[1]*s2*v[1],withoutRHS);
		eqns.clear();
		eqns.push_back(sin(w[0])*s1);
		eqns.push_back(pow(2,w[1])*s2+I*w[1]*s2);
		eqns.push_back(w[2]*s3);
		TS_ASSERT(includes(l.begin(),l.end(),eqns.begin(),eqns.end()));	
		
		l.clear();
		eqns.clear();
		TS_ASSERT_THROWS(GetCoefficients<V> (l,v[0]*v[1]+w[0]),WedgeException<std::runtime_error>);
		TS_ASSERT_THROWS(GetCoefficients<Poly1<V> > (l,v[0]*v[1]+w[0]),WedgeException<std::runtime_error>);
		GetCoefficients<Poly<V> > (l,v[0]*v[1]+w[0]);
		eqns.push_back(1);
		eqns.push_back(w[0]);
		TS_ASSERT(includes(l.begin(),l.end(),eqns.begin(),eqns.end()));		

		l.clear();
		TS_ASSERT_THROWS(GetCoefficients<V> (l,v[0]*w[0]+1),WedgeException<std::runtime_error>);
		GetCoefficients<Poly1<V> > (l,v[0]*w[0]+1);
		TS_ASSERT(includes(l.begin(),l.end(),eqns.begin(),eqns.end()));		
		l.clear();
		GetCoefficients<Poly<V> > (l,v[0]*w[0]+1);
		TS_ASSERT(includes(l.begin(),l.end(),eqns.begin(),eqns.end()));
	}
	
	void testStreams()
	{
		cout<<latex;
		TS_ASSERT(IsStreamTeX(cout));
		cout<<dflt;
		TS_ASSERT(!IsStreamTeX(cout));
		
		//only checks that no exception is thrown
		stringstream f; 
		const int max_size=5;
		exvector elements(max_size);
		for (int n=0;n<max_size;++n) {
			f<<dflt; f<<vector<V>(n);
			f<<latex; f<<vector<V>(n);
			f<<dflt; f<<exvector(n);
			f<<latex; f<<exvector(n);
			set<ex,ex_is_less> b(elements.begin(),elements.begin()+n);
			f<<dflt<<b;
			f<<latex<<b;
			f<<dflt<<exmap();
			f<<latex<<exmap();
		}

		list<V> a; 
		a.push_back(V(N.v));a.push_back(V(N.v));a.push_back(V(N.v));

		f<<dflt<<a;
		f<<latex<<a;
		list<ex> b;
		b.push_back(V(N.v));b.push_back(V(N.v));b.push_back(V(N.v));
		f<<dflt<<b;
		f<<latex<<b;
	}
	
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

	void testVariableNames() 
	{
		TS_ASSERT_EQUALS(N.x.plain(),"x");
		TS_ASSERT_EQUALS(N.T.plain(),"T");
		TS_ASSERT_EQUALS(N.psi.plain(),"psi");
		TS_ASSERT_EQUALS(N.Gamma.plain(),"Gamma");
		TS_ASSERT_EQUALS(N.x(1).plain(),"x1");
		TS_ASSERT_EQUALS(N.x(10).plain(),"x10");
		TS_ASSERT_EQUALS(N.x(1,2).plain(),"x1$2");
		TS_ASSERT_EQUALS(N.x(1).Tilde().plain(),"x~1");
		TS_ASSERT_EQUALS(N.x(1,2).Plus().plain(),"x+1$2");
		TS_ASSERT_EQUALS(N.x(1,2,3).plain(),"x1$2$3");
		TS_ASSERT_EQUALS(N.x(1,2,3)(4).plain(),"x1$2$3$4");
		TS_ASSERT_EQUALS(N.x(1,2)(3)(4).plain(),"x1$2$3$4");
	
	}

};

//test the parsing functions
class ParseTestSuite : public CxxTest::TestSuite 
{
public:
	void testParseMaple()
	{
		symbol a("a"),b("b");
		TS_ASSERT_THROWS(ParseMapleExpression(lst(),"a*2"), std::invalid_argument);
		TS_ASSERT_EQUALS(ParseMapleExpression(a,"a*2"), a*2);
		TS_ASSERT_EQUALS(ParseMapleExpression(a,"a*2=a"), a*2==a);
		ExVector expressions=ParseMapleExpressions(lst(a,b),"a^2,b-sqrt(a)");
		TS_ASSERT_EQUALS(expressions[0], a*a);
		TS_ASSERT_EQUALS(expressions[1], b-sqrt(a));

	}


	void testParseCocoa()
	{
#ifdef HAVE_COCOA
		logging_level=LOGGING_LEVEL_DEBUG;
		char s[]="x[1],x[1]+x[2], -x[1]+x[2],-2*x[1]-x[2],-2*x[1]^2+x[2]^2*x[0],3*x[1]*(x[1]-x[2]),1/2*x[1],x[0]^2 +1/2*x[0]";
		exvector x(10);
		for (int i=0;i<x.size();i++) x[i]=symbol("x["+ToString(i)+"]");
		exvector v=ParseCocoaExpressions(x,s);
		TS_ASSERT_EQUALS(v[0],x[1]);
		TS_ASSERT_EQUALS(v[1],x[1]+x[2]);
		TS_ASSERT_EQUALS(v[2],-x[1]+x[2]);
		TS_ASSERT_EQUALS(v[3],-2*x[1]-x[2]);
		TS_ASSERT_EQUALS(v[4],-2*pow(x[1],2)+pow(x[2],2)*x[0]);
		TS_ASSERT_EQUALS(v[5],3*x[1]*(x[1]-x[2]));
		TS_ASSERT_EQUALS(v[6],x[1]/2);
		TS_ASSERT_EQUALS(v[7],pow(x[0],2)+x[0]/2);
#endif
	}
	void testParseForms()
	{
		ConcreteManifold M(13);		
		TS_ASSERT_EQUALS(ParseDifferentialForm(M.e(),"[sqrt(2)]*1"),sqrt(ex(2))*M.e(1));		
		TS_ASSERT_EQUALS(ParseDifferentialForm(M.e(),"0"),0);
		
		exvector forms=ParseDifferentialForms(M.e(),"-1,2+[sqrt(3)]*34+0,123+7*4,-3-(1-2),1-2/3*(-3+45)");
		TS_ASSERT_EQUALS(-M.e(1),forms[0]);
		TS_ASSERT_EQUALS(M.e(2)+sqrt(ex(3))*M.e(3)*M.e(4),forms[1]);
		TS_ASSERT_EQUALS(M.e(1)*M.e(2)*M.e(3)+7*M.e(4),forms[2]);
		TS_ASSERT_EQUALS(-M.e(3)-M.e(1)+M.e(2),forms[3]);
		TS_ASSERT_EQUALS(M.e(1)-numeric(2,3)*(-M.e(3)+M.e(4)*M.e(5)),forms[4]);

		symbol a("a");
		forms=ParseDifferentialForms(M.e(),"1+[a]*a4,1/3*52^[1+sqrt(2)]*3,[pow(2,-1/3)]*(1+34),[a^(-2)]*(1+2)3,2*3+4*[sqrt(3)]*1",lst(a));
		TS_ASSERT_EQUALS(M.e(1)+a*M.e(10)*M.e(4),forms[0]);
		TS_ASSERT_EQUALS(1/ex(3)*M.e(5)*M.e(2)*(1+sqrt(ex(2)))*M.e(3),forms[1]);
		TS_ASSERT_EQUALS(pow(2,-1/ex(3))*(M.e(1)+M.e(3)*M.e(4)),forms[2]);
		TS_ASSERT_EQUALS(1/(a*a)*(M.e(1)+M.e(2))*M.e(3),forms[3]);
		TS_ASSERT_EQUALS(2*M.e(3)+4*sqrt(ex(3))*M.e(1),forms[4]);
		TS_ASSERT_THROWS(ParseDifferentialForm(M.e(),"[a]*2"), std::invalid_argument);

	}
};

#endif /*BASE_H_*/
