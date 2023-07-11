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
#ifndef POLY_H_
#define POLY_H_
#include <cxxtest/TestSuite.h>
#include "test.h"
#include "wedge/polynomialalgebra/cocoapolyalg.h"
#include "wedge/polynomialalgebra/polybasis.h"
#include "wedge/manifolds/coordinates.h"


WEDGE_DECLARE_NAMED_ALGEBRAIC(V,realsymbol)
WEDGE_DECLARE_NAMED_ALGEBRAIC(U,realsymbol)
class PolyTestSuite : public  CxxTest::TestSuite
{
	typedef PolyBasisImplementation<V,DefaultPolyAlgorithms> PolyBasisImpl;
public:
	void testPolyBasis()
	{
		V x(N.x),y(N.y),z(N.z),t(N.t),s(N.s);
		PolyBasis<V> I;
		TS_ASSERT(I.IdealIsZero());
		TS_ASSERT(I.empty());
		I.push_back(0);
		TS_ASSERT(I.IdealIsZero());
		TS_ASSERT(!I.empty());
		TS_ASSERT_EQUALS(I.Variables().size(),0);
		
		I.push_back(y+pow(x,3));
		TS_ASSERT(!I.IdealIsZero());
		TS_ASSERT(!I.empty());
		
		TS_ASSERT_EQUALS(count(I.Variables().begin(),I.Variables().end(),x),1);
		TS_ASSERT_EQUALS(count(I.Variables().begin(),I.Variables().end(),y),1);
		TS_ASSERT_EQUALS(count(I.Variables().begin(),I.Variables().end(),z),0);
		TS_ASSERT_EQUALS(count(I.Variables().begin(),I.Variables().end(),t),0);
		TS_ASSERT_EQUALS(I.VariableByName("x"),x);
		TS_ASSERT_EQUALS(I.VariableByName("y"),y);
		TS_ASSERT_THROWS(I.VariableByName("z"),InvalidArgument);
		TS_ASSERT_THROWS(I.VariableByName("t"),InvalidArgument);
		
		I.push_back(x*y+x*t);
		TS_ASSERT_EQUALS(count(I.Variables().begin(),I.Variables().end(),x),1);
		TS_ASSERT_EQUALS(count(I.Variables().begin(),I.Variables().end(),y),1);
		TS_ASSERT_EQUALS(count(I.Variables().begin(),I.Variables().end(),z),0);
		TS_ASSERT_EQUALS(count(I.Variables().begin(),I.Variables().end(),t),1);		
		TS_ASSERT_EQUALS(I.VariableByName("x"),x);
		TS_ASSERT_EQUALS(I.VariableByName("y"),y);
		TS_ASSERT_THROWS(I.VariableByName("z"),InvalidArgument);
		TS_ASSERT_EQUALS(I.VariableByName("t"),t);
		
		I.push_back(x*z+y*pow(t,4));
		TS_ASSERT_EQUALS(*(I.begin()),0); 		
		TS_ASSERT_EQUALS(*(I.begin()+1),y+pow(x,3));
		TS_ASSERT_EQUALS(*(I.begin()+2),x*y+x*t);
		TS_ASSERT_EQUALS(*(I.begin()+3),x*z+y*pow(t,4));
		
		ex pol=x*s;
		I.insert(I.begin(),&pol,&pol+1);
		TS_ASSERT_EQUALS(count(I.Variables().begin(),I.Variables().end(),x),1);
		TS_ASSERT_EQUALS(count(I.Variables().begin(),I.Variables().end(),y),1);
		TS_ASSERT_EQUALS(count(I.Variables().begin(),I.Variables().end(),z),1);
		TS_ASSERT_EQUALS(count(I.Variables().begin(),I.Variables().end(),t),1);		
		TS_ASSERT_EQUALS(count(I.Variables().begin(),I.Variables().end(),s),1);
		TS_ASSERT_EQUALS(I.VariableByName("x"),x);
		TS_ASSERT_EQUALS(I.VariableByName("y"),y);
		TS_ASSERT_EQUALS(I.VariableByName("z"),z);
		TS_ASSERT_EQUALS(I.VariableByName("t"),t);
		TS_ASSERT_EQUALS(I.VariableByName("s"),s);
		TS_ASSERT_EQUALS(*(I.begin()),x*s);		
		TS_ASSERT_EQUALS(*(I.begin()+1),0); 		
		TS_ASSERT_EQUALS(*(I.begin()+2),y+pow(x,3));
		TS_ASSERT_EQUALS(*(I.begin()+3),x*y+x*t);
		TS_ASSERT_EQUALS(*(I.begin()+4),x*z+y*pow(t,4));
		I.TakeVariableToFront(s);
		I.TakeVariableToFront(z);
		I.TakeVariableToBack(x);
		I.TakeVariableToBack(x);
		I.TakeVariableToFront(x);				
		TS_ASSERT_EQUALS(count(I.Variables().begin(),I.Variables().end(),x),1);
		TS_ASSERT_EQUALS(count(I.Variables().begin(),I.Variables().end(),y),1);
		TS_ASSERT_EQUALS(count(I.Variables().begin(),I.Variables().end(),z),1);
		TS_ASSERT_EQUALS(count(I.Variables().begin(),I.Variables().end(),t),1);		
		TS_ASSERT_EQUALS(count(I.Variables().begin(),I.Variables().end(),s),1);
		TS_ASSERT_EQUALS(I.VariableByName("x"),x);
		TS_ASSERT_EQUALS(I.VariableByName("y"),y);
		TS_ASSERT_EQUALS(I.VariableByName("z"),z);
		TS_ASSERT_EQUALS(I.VariableByName("t"),t);
		TS_ASSERT_EQUALS(I.VariableByName("s"),s);
		TS_ASSERT_EQUALS(I.Variables()[0],x);
		
		I.clear();
		TS_ASSERT(I.IdealIsZero());
		TS_ASSERT(I.empty());
		TS_ASSERT_EQUALS(I.Variables().size(),0);
		TS_ASSERT_THROWS(I.VariableByName("x"),InvalidArgument);
		TS_ASSERT_THROWS(I.VariableByName("y"),InvalidArgument);
		TS_ASSERT_THROWS(I.VariableByName("z"),InvalidArgument);
		TS_ASSERT_THROWS(I.VariableByName("t"),InvalidArgument);
		TS_ASSERT_THROWS(I.VariableByName("s"),InvalidArgument);
		TS_ASSERT_THROWS(I.TakeVariableToFront(s),InvalidArgument);
		TS_ASSERT_THROWS(I.TakeVariableToBack(s),InvalidArgument);
	}
	
	void testGroebner()
	{
		V x(N.x),y(N.y),z(N.z),t(N.t),s(N.s);
		PolyBasisImpl I;
		TS_ASSERT(I.IdealIsZero());
		TS_ASSERT(!I.IdealIsOne());
		TS_ASSERT(I.empty());
		I.push_back(0);
		TS_ASSERT(I.IdealIsZero());
		TS_ASSERT(!I.IdealIsOne());
		TS_ASSERT(!I.empty());
		TS_ASSERT_EQUALS(I.Variables().size(),0);
		
		I.push_back(y+pow(x,3));
		TS_ASSERT(!I.IdealIsZero());
		TS_ASSERT(!I.empty());
		TS_ASSERT(!I.IdealIsOne());

		I.push_back(x*y+pow(x,4)+1);
		TS_ASSERT(!I.IdealIsZero());
		TS_ASSERT(!I.empty());
		TS_ASSERT(I.IdealIsOne());
		TS_ASSERT(I.IdealContains(x));
		TS_ASSERT(I.IdealContains(z));
		
		I.Reduce();
		TS_ASSERT(!I.IdealIsZero());
		TS_ASSERT(!I.empty());
		TS_ASSERT(I.IdealIsOne());
		TS_ASSERT(I.IdealContains(x));
				
		I.clear();
		TS_ASSERT(I.IdealIsZero());
		TS_ASSERT(I.empty());
		TS_ASSERT(!I.IdealIsOne());
		TS_ASSERT(!I.IdealContains(x));
		I.Reduce();
		TS_ASSERT(I.IdealIsZero());
		TS_ASSERT(I.empty());
		TS_ASSERT(!I.IdealIsOne());
		TS_ASSERT(!I.IdealContains(x));
		I.push_back(0);
		I.Reduce();
		TS_ASSERT(I.IdealIsZero());
		TS_ASSERT(I.empty());
		TS_ASSERT(!I.IdealIsOne());
		TS_ASSERT(!I.IdealContains(x));
		
		I.push_back(x*y+z);
		I.push_back(t*s);
		TS_ASSERT(I.IdealContains(y*(x*y+z)+z*t*s));
		TS_ASSERT_EQUALS(I.ReduceModuloIdeal(y*(x*y+z)+z*t*s),0);		
		TS_ASSERT_EQUALS(I.ReduceModuloIdeal(s*s*y*x),I.ReduceModuloIdeal(-s*s*z));
		I.clear();
		TS_ASSERT_EQUALS(I.ReduceModuloIdeal(s*s*y*x),s*s*y*x);
		
		I.push_back(-3+z*z);
		I.push_back(x*x*x+z*y);
		I.push_back(x*x*x*x*z+3*x*y+z);
		TS_ASSERT(I.IdealIsOne());
		I.Reduce();
		TS_ASSERT(I.IdealIsOne());

		I.clear();
		I.push_back(-x+z*z);
		TS_ASSERT(!I.RadicalContains(z));
		TS_ASSERT(!I.RadicalContains(x));
		I.push_back(x*x);
		TS_ASSERT(I.RadicalContains(z));
		TS_ASSERT(I.RadicalContains(x));

		PolyBasisImpl J;
		TS_ASSERT(I.Intersected(J).IdealIsZero());
		TS_ASSERT(J.Intersected(I).IdealIsZero());
		J.push_back(0);
		TS_ASSERT(I.Intersected(J).IdealIsZero());
		TS_ASSERT(J.Intersected(I).IdealIsZero());
		J.push_back(x);
		TS_ASSERT(J.Intersected(I).IdealContains(x*x));
		TS_ASSERT(I.Intersected(J).IdealContains(x*x));
		TS_ASSERT(!J.Intersected(I).IdealContains(x));
		TS_ASSERT(!I.Intersected(J).IdealContains(x));
		J.clear();
		J.push_back(0);
		TS_ASSERT(!I.Intersected(J).IdealContains(x*x));
		TS_ASSERT(!J.Intersected(I).IdealContains(x*x));
		J.push_back(-2);
		TS_ASSERT(I.Intersected(J).IdealContains(x*x));
		TS_ASSERT(J.Intersected(I).IdealContains(x*x));
		
		I.clear();
		J.clear();
		I.push_back(x);
		J.push_back(x+1);
		TS_ASSERT(I.Intersected(J).IdealContains(x*(x+1)));
		TS_ASSERT(!I.Intersected(J).IdealContains(x*x));
	}

	void testGroebner_R_multipleroots()
	{
		V x(N.x),y(N.y),z(N.z),t(N.t),s(N.s);
		PolyBasisImplementation<V,CocoaPolyAlgorithms_R> I;
		exvector roots;
		roots.push_back(sqrt(ex(2)));
		roots.push_back(sqrt(ex(3)));
		roots.push_back(sqrt(ex(6)));
		TS_ASSERT(I.IdealIsZero());
		TS_ASSERT(!I.IdealIsOne());
		I.push_back(roots[0]*x+roots[1]*y+roots[2]);
		TS_ASSERT(!I.IdealIsZero());
		TS_ASSERT(!I.IdealIsOne());
		TS_ASSERT(!I.IdealContains(roots[1]));
		I.push_back(2*x+roots[2]*y);
		TS_ASSERT(!I.IdealIsZero());
		TS_ASSERT(I.IdealIsOne());
		TS_ASSERT(I.IdealContains(roots[1]));
		I.clear();
		I.push_back(roots[0]*x+roots[1]*y+roots[2]);
		I.push_back(roots[1]*x+roots[2]*y);
		I.push_back(roots[2]*x+roots[0]*y);
		TS_ASSERT(!I.IdealIsZero());
		TS_ASSERT(I.IdealIsOne());

		roots.clear();
		I.clear();
		roots.push_back(pow(10,ex(1)/2));
		roots.push_back(pow(10,ex(1)/3));
		roots.push_back(pow(10,ex(1)/6));
		TS_ASSERT(I.IdealIsZero());
		TS_ASSERT(!I.IdealIsOne());
		I.push_back(roots[0]*x+roots[1]*y+roots[2]);
		TS_ASSERT(!I.IdealIsZero());
		TS_ASSERT(!I.IdealIsOne());
		TS_ASSERT(!I.IdealContains(roots[1]));
		I.push_back(roots[1]*x+roots[2]*y);
		TS_ASSERT(!I.IdealIsZero());
		TS_ASSERT(I.IdealIsOne());
		TS_ASSERT(I.IdealContains(roots[1]));
		I.clear();
		I.push_back(roots[0]*x+roots[1]*y+roots[2]);
		I.push_back(roots[1]*x+roots[2]*y);
		I.push_back(roots[2]*x+roots[0]*y);
		TS_ASSERT(!I.IdealIsZero());
		TS_ASSERT(I.IdealIsOne());

		roots.clear();
		I.clear();
		roots.push_back(pow(ex(2)/3,ex(1)/2));
		roots.push_back(pow(ex(4)/5,ex(1)/3));
		roots.push_back(pow(ex(15)/2,ex(1)/6));
		I.push_back(roots[0]*x+roots[1]*y+roots[2]*z);
		I.push_back(roots[1]*x+roots[2]*y+roots[0]*z);
		TS_ASSERT(!I.IdealContains(x));
		TS_ASSERT(!I.IdealContains(y));
		TS_ASSERT(!I.IdealContains(z));
		TS_ASSERT(!I.IdealContains(roots[0]*z));
		I.push_back(roots[2]*x+roots[0]*y+roots[1]*z);
		TS_ASSERT(I.IdealContains(x));
		TS_ASSERT(I.IdealContains(y));
		TS_ASSERT(I.IdealContains(z));
		TS_ASSERT(I.IdealContains(roots[0]*z));		

		I.clear();
		I.push_back(-roots[0]*x+roots[1]*z*z);
		TS_ASSERT(!I.RadicalContains(z));
		TS_ASSERT(!I.RadicalContains(x));
		I.push_back(roots[2]*x*x);
		TS_ASSERT(I.RadicalContains(z));
		TS_ASSERT(I.RadicalContains(x));
	}

	void testGroebner_R()
	{
		V x(N.x),y(N.y),z(N.z),t(N.t),s(N.s);
		PolyBasisImplementation<V,CocoaPolyAlgorithms_R> I;
		TS_ASSERT(I.IdealIsZero());
		TS_ASSERT(!I.IdealIsOne());
		TS_ASSERT(I.empty());
		I.push_back(0);
		TS_ASSERT(I.IdealIsZero());
		TS_ASSERT(!I.IdealIsOne());
		TS_ASSERT(!I.empty());
		TS_ASSERT_EQUALS(I.Variables().size(),0);
		
		I.push_back(sqrt(ex(3))*y+pow(x,3));
		TS_ASSERT(!I.IdealIsZero());
		TS_ASSERT(!I.empty());
		TS_ASSERT(!I.IdealIsOne());

		I.push_back(3*x*y+sqrt(ex(3))*pow(x,4)+sqrt(ex(3)));
		TS_ASSERT(!I.IdealIsZero());
		TS_ASSERT(!I.empty());
		TS_ASSERT(I.IdealIsOne());
		TS_ASSERT(I.IdealContains(1));
		TS_ASSERT(I.IdealContains(x));
		TS_ASSERT(I.IdealContains(z));
		
		I.Reduce();
		TS_ASSERT(!I.IdealIsZero());
		TS_ASSERT(!I.empty());
		TS_ASSERT(I.IdealIsOne());
		TS_ASSERT(I.IdealContains(x));
				
		I.clear();
		TS_ASSERT(I.IdealIsZero());
		TS_ASSERT(I.empty());
		TS_ASSERT(!I.IdealIsOne());
		TS_ASSERT(!I.IdealContains(x));
		I.Reduce();
		TS_ASSERT(I.IdealIsZero());
		TS_ASSERT(I.empty());
		TS_ASSERT(!I.IdealIsOne());
		TS_ASSERT(!I.IdealContains(x));
		I.push_back(0);
		I.Reduce();
		TS_ASSERT(I.IdealIsZero());
		TS_ASSERT(I.empty());
		TS_ASSERT(!I.IdealIsOne());
		TS_ASSERT(!I.IdealContains(x));
		
		I.push_back(t*s);
		I.push_back(x*y+z*pow(2,ex(1)/3));		
		TS_ASSERT_EQUALS(I.ReduceModuloIdeal(s*s*y*x),I.ReduceModuloIdeal(-s*s*z*pow(2,ex(1)/3)));
		I.clear();
		TS_ASSERT_EQUALS(I.ReduceModuloIdeal(s*s*y*x),s*s*y*x);
	}

	//verify that everything works as expected with symbols of varying type
	void testVariantSymbols()
	{
		ManifoldWithCoordinates M(3);
		Function f(N.f);
		PolyBasisImplementation<Function,DefaultPolyAlgorithms> I;
		I.push_back(M.LieDerivative(M.e(1),M.LieDerivative(M.e(1),f))+f);
		I.push_back(M.LieDerivative(M.e(1),M.LieDerivative(M.e(1),f))+f+M.LieDerivative(M.e(2),f));
		I.push_back(M.LieDerivative(M.e(2),f)*f);
		TS_ASSERT_EQUALS(I.size(),3);
		TS_ASSERT_EQUALS(count(I.Variables().begin(),I.Variables().end(),f),1);
		TS_ASSERT_EQUALS(count(I.Variables().begin(),I.Variables().end(),M.LieDerivative(M.e(2),f)),1);
		TS_ASSERT_EQUALS(count(I.Variables().begin(),I.Variables().end(),M.LieDerivative(M.e(1),M.LieDerivative(M.e(1),f))),1);
		TS_ASSERT_EQUALS(count(I.Variables().begin(),I.Variables().end(),M.LieDerivative(M.e(1),f)),0);
		TS_ASSERT(I.IdealContains(M.LieDerivative(M.e(2),f)));
		TS_ASSERT(!I.IdealContains(f));
		I.Reduce();
		TS_ASSERT(I.IdealContains(M.LieDerivative(M.e(2),f)));
		TS_ASSERT(!I.IdealContains(f));
	}

	void testEliminate() 
	{
		U x(N.x),y(N.y),z(N.z),t(N.t),s(N.s);
		PolyBasisImplementation<U,DefaultPolyAlgorithms> basis;
		basis.push_back(x-y*y);
		basis.push_back(t-y*y);
		basis.push_back(x-x*z);
		lst subs=basis.Eliminate();
		cout<<subs<<endl;
		cout<<basis<<endl;
		cout<<(x-z*y*y).subs(subs)<<endl;
		TS_ASSERT((x-t).subs(subs).is_zero());
		TS_ASSERT(basis.IdealContains((x-z*y*y).subs(subs)));
	}
};

#endif
