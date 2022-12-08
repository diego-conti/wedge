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
#ifndef BASE_H_
#define BASE_H_

#include "wedge/base/wedgebase.h"
#include "wedge/base/wedgealgebraic.h"
#include "wedge/base/expressions.h"
#include "wedge/convenience/latex.h"
#include "wedge/base/utilities.h"
#include <cxxtest/TestSuite.h>
#include "test.h"

using namespace Wedge;
using namespace GiNaC;
/**
	@author Diego Conti <diego.conti@unipi.it>
	@brief Test suite for the Base module in %Wedge
		
*/


WEDGE_DECLARE_NAMED_ALGEBRAIC(A,symbol);
WEDGE_DECLARE_ALGEBRAIC(B,A)
WEDGE_DECLARE_NAMED_ALGEBRAIC(C,realsymbol);

WEDGE_DECLARE_NAMED_ALGEBRAIC(V,Vector)

WEDGE_DECLARE_NAMED_ALGEBRAIC(W,Vector)


//test expressions.h
class ExpressionsTestSuite : public CxxTest::TestSuite 
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
	

	void testStreamsOutput() {
		cout<<latex;
		TS_ASSERT(IsStreamTeX(cout));
		cout<<dflt;
		TS_ASSERT(!IsStreamTeX(cout));

		ex a=A{N.a}, phi=C{N.phi};
		stringstream f;
		f<<a<<phi;
		TS_ASSERT_EQUALS(f.str(),"aphi"); f.str("");
		f<<latex<<a<<phi;
		TS_ASSERT_EQUALS(f.str(),"a\\phi"); f.str("");
		ex a1=A{N.a(1)}, phi1=C{N.phi(1)};
		f<<latex<<a1<<phi1;
		TS_ASSERT_EQUALS(f.str(),"a_1\\phi_1"); f.str("");
		f<<dflt<<a1<<phi1;
		TS_ASSERT_EQUALS(f.str(),"a1phi1"); f.str("");
	}

	void testStreamsDontThrow()
	{
		
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


class UtilitiesTestSuite : public CxxTest::TestSuite
{
public:
	  class IterateOverSubsetsTest: public IterateOverSubsets
	  {
		bool shouldStopEarly;
	  public:
		list<vector<int> > sequence;
		IterateOverSubsetsTest(bool shouldStopEarly) {
			this->shouldStopEarly=shouldStopEarly;
			vector<int> v(2);
			v[0]=0;v[1]=1; sequence.push_front(v);
			v[0]=0;v[1]=2; sequence.push_front(v);
			v[0]=1;v[1]=2; sequence.push_front(v);
			Iterate(2,3);
		}
		bool Apply(const vector<int>& v) {
			TS_ASSERT_EQUALS(v,sequence.back());
			sequence.pop_back();
			return (!shouldStopEarly || sequence.size()>1);
		}
	  };
	void testIterate() {
		IterateOverSubsetsTest t(false);
		TS_ASSERT(t.sequence.empty());
		IterateOverSubsetsTest w(true);
		TS_ASSERT(w.sequence.size()==1);
	}
};

#endif /*BASE_H_*/
