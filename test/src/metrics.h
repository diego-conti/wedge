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


#ifndef SRC_METRICS_H_
#define SRC_METRICS_H_
#include <cxxtest/TestSuite.h>
#include "test.h"
#include "wedge/linearalgebra/bilinearform.h"
#include "wedge/structures/pseudoriemannianstructure.h"

using namespace Wedge;

class MetricsTestSuite : public CxxTest::TestSuite  {
public:
	void testStandardScalarProduct() {
		ConcreteManifold M(3);
		Frame e=ParseDifferentialForms(M.e(),"1+2*2, 2*2, 3");
		auto dual=e.dual();
		auto g =StandardScalarProduct(e,1);
		TS_ASSERT_EQUALS(g.OnForms(e(1),e(1)),1);
		TS_ASSERT_EQUALS(g.OnForms(e(1),e(2)),0);
		TS_ASSERT_EQUALS(g.OnForms(e(1),e(3)),0);
		TS_ASSERT_EQUALS(g.OnForms(e(2),e(2)),-1);
		TS_ASSERT_EQUALS(g.OnForms(e(3),e(3)),-1);
		TS_ASSERT_EQUALS(g.OnForms(e(1)*e(2),e(1)*e(2)),-1);
		TS_ASSERT_EQUALS(g.OnForms(e(1)*e(2),e(1)*e(3)),0);

		TS_ASSERT_EQUALS(g.OnVectors(dual(1),dual(1)),1);
		TS_ASSERT_EQUALS(g.OnVectors(dual(1),dual(2)),0);
		TS_ASSERT_EQUALS(g.OnVectors(dual(1),dual(3)),0);
		TS_ASSERT_EQUALS(g.OnVectors(dual(2),dual(2)),-1);
		TS_ASSERT_EQUALS(g.OnVectors(dual(3),dual(3)),-1);
	}
	void testScalarProductByOrthonormalFrame() {
		ConcreteManifold M(3);
		Frame e=ParseDifferentialForms(M.e(),"1+2*2, 2*2, 3");
		auto dual=e.dual();
		auto g =ScalarProductByOrthonormalFrame::FromSequenceOfSigns(e,{1,-1,-1});
		TS_ASSERT_EQUALS(g.OnForms(e(1),e(1)),1);
		TS_ASSERT_EQUALS(g.OnForms(e(1),e(2)),0);
		TS_ASSERT_EQUALS(g.OnForms(e(1),e(3)),0);
		TS_ASSERT_EQUALS(g.OnForms(e(2),e(2)),-1);
		TS_ASSERT_EQUALS(g.OnForms(e(3),e(3)),-1);
		TS_ASSERT_EQUALS(g.OnForms(e(1)*e(2),e(1)*e(2)),-1);
		TS_ASSERT_EQUALS(g.OnForms(e(1)*e(2),e(1)*e(3)),0);

		TS_ASSERT_EQUALS(g.OnVectors(dual(1),dual(1)),1);
		TS_ASSERT_EQUALS(g.OnVectors(dual(1),dual(2)),0);
		TS_ASSERT_EQUALS(g.OnVectors(dual(1),dual(3)),0);
		TS_ASSERT_EQUALS(g.OnVectors(dual(2),dual(2)),-1);
		TS_ASSERT_EQUALS(g.OnVectors(dual(3),dual(3)),-1);
	}

	void testScalarProductByOrthonormalFrame2() {
		ConcreteManifold M(3);
		Frame e=ParseDifferentialForms(M.e(),"1+2*2, 2*2, 3");
		auto dual=e.dual();
		auto g =ScalarProductByOrthonormalFrame::FromTimelikeIndices(e,{1,3});
		TS_ASSERT_EQUALS(g.OnForms(e(1),e(1)),-1);
		TS_ASSERT_EQUALS(g.OnForms(e(1),e(2)),0);
		TS_ASSERT_EQUALS(g.OnForms(e(1),e(3)),0);
		TS_ASSERT_EQUALS(g.OnForms(e(2),e(2)),1);
		TS_ASSERT_EQUALS(g.OnForms(e(3),e(3)),-1);
		TS_ASSERT_EQUALS(g.OnForms(e(1)*e(2),e(1)*e(2)),-1);
		TS_ASSERT_EQUALS(g.OnForms(e(1)*e(2),e(1)*e(3)),0);

		TS_ASSERT_EQUALS(g.OnVectors(dual(1),dual(1)),-1);
		TS_ASSERT_EQUALS(g.OnVectors(dual(1),dual(2)),0);
		TS_ASSERT_EQUALS(g.OnVectors(dual(1),dual(3)),0);
		TS_ASSERT_EQUALS(g.OnVectors(dual(2),dual(2)),1);
		TS_ASSERT_EQUALS(g.OnVectors(dual(3),dual(3)),-1);
	}

	void testScalarProductByOrtohogonalFrame() {
		ConcreteManifold M(3);
		Frame e=ParseDifferentialForms(M.e(),"1+2*2, 2*2, 3");
		auto dual=e.dual();
		auto g =ScalarProductByOrthogonalFrame::FromVectorSquareNorms(e, {1,2,3});
		TS_ASSERT_EQUALS(g.OnVectors(dual(1),dual(2)),0);
		TS_ASSERT_EQUALS(g.OnVectors(dual(1),dual(3)),0);
		TS_ASSERT_EQUALS(g.OnVectors(dual(2),dual(3)),0);
		TS_ASSERT_EQUALS(g.OnVectors(dual(1),dual(1)),1);
		TS_ASSERT_EQUALS(g.OnVectors(dual(2),dual(2)),2);
		TS_ASSERT_EQUALS(g.OnVectors(dual(3),dual(3)),3);

		TS_ASSERT_EQUALS(g.OnForms(e(1),e(2)),0);
		TS_ASSERT_EQUALS(g.OnForms(e(1),e(3)),0);
		TS_ASSERT_EQUALS(g.OnForms(e(2),e(3)),0);
		TS_ASSERT_EQUALS(g.OnForms(e(1),e(1)),1);
		TS_ASSERT_EQUALS(g.OnForms(e(2),e(2)),1/ex(2));
		TS_ASSERT_EQUALS(g.OnForms(e(3),e(3)),1/ex(3));		
	}

	void testScalarProductByOrtohogonalFrame2() {
		ConcreteManifold M(3);
		Frame e=ParseDifferentialForms(M.e(),"1+2*2, 2*2, 3");
		auto dual=e.dual();
		auto g =ScalarProductByOrthogonalFrame::FromFormSquareNorms(e, {1,1/ex(2),1/ex(3)});
		TS_ASSERT_EQUALS(g.OnVectors(dual(1),dual(2)),0);
		TS_ASSERT_EQUALS(g.OnVectors(dual(1),dual(3)),0);
		TS_ASSERT_EQUALS(g.OnVectors(dual(2),dual(3)),0);
		TS_ASSERT_EQUALS(g.OnVectors(dual(1),dual(1)),1);
		TS_ASSERT_EQUALS(g.OnVectors(dual(2),dual(2)),2);
		TS_ASSERT_EQUALS(g.OnVectors(dual(3),dual(3)),3);

		TS_ASSERT_EQUALS(g.OnForms(e(1),e(2)),0);
		TS_ASSERT_EQUALS(g.OnForms(e(1),e(3)),0);
		TS_ASSERT_EQUALS(g.OnForms(e(2),e(3)),0);
		TS_ASSERT_EQUALS(g.OnForms(e(1),e(1)),1);
		TS_ASSERT_EQUALS(g.OnForms(e(2),e(2)),1/ex(2));
		TS_ASSERT_EQUALS(g.OnForms(e(3),e(3)),1/ex(3));		
	}


	void testCliffordRiemannianOdd() {
		ConcreteManifold M(3);
		Frame e=ParseDifferentialForms(M.e(),"1+2*2, 2*2, 3");
		auto dual=e.dual();
		ex d1=dual(1), d2=dual(2), d3=dual(3);
		auto g =PseudoRiemannianStructureByOrthonormalFrame::FromTimelikeIndices(&M,e,{},CliffordConvention::BAUM_KATH);
		for (int i=0;i<g.DimensionOfSpinorRepresentation();++i) {
			ex u=g.u(i);
			TS_ASSERT_EQUALS(g.CliffordDot(d1,g.CliffordDot(d1,u)),-u);
			TS_ASSERT_EQUALS(g.CliffordDot(d2,g.CliffordDot(d2,u)),-u);
			TS_ASSERT_EQUALS(g.CliffordDot(d3,g.CliffordDot(d3,u)),-u);
			for (auto X : dual)
			for (auto Y : dual) {
				if (X!=Y) TS_ASSERT_EQUALS(g.CliffordDot(X,g.CliffordDot(Y,u))+g.CliffordDot(Y,g.CliffordDot(X,u)),0);
			}
		}
		auto timelikeindices=vector<int>{};
		TS_ASSERT_EQUALS(g.ScalarProduct().TimelikeIndices(),timelikeindices);

	}
	void testCliffordRiemannianEven() {
		ConcreteManifold M(4);		
		ex d1=M.e(1), d2=M.e(2),d3=M.e(3),d4=M.e(4);
		auto g =PseudoRiemannianStructureByOrthonormalFrame::FromTimelikeIndices(&M,M.e(),{},CliffordConvention::BAUM_KATH);
		for (int i=0;i<g.DimensionOfSpinorRepresentation();++i) {
			ex u=g.u(i);
			TS_ASSERT_EQUALS(g.CliffordDot(d1,g.CliffordDot(d1,u)),-u);
			TS_ASSERT_EQUALS(g.CliffordDot(d2,g.CliffordDot(d2,u)),-u);
			TS_ASSERT_EQUALS(g.CliffordDot(d3,g.CliffordDot(d3,u)),-u);
			TS_ASSERT_EQUALS(g.CliffordDot(d4,g.CliffordDot(d4,u)),-u);
			for (auto X : M.e())
			for (auto Y : M.e())
				if (X!=Y) TS_ASSERT_EQUALS(g.CliffordDot(X,g.CliffordDot(Y,u))+g.CliffordDot(Y,g.CliffordDot(X,u)),0);
		}
		auto timelikeindices=vector<int>{};
		TS_ASSERT_EQUALS(g.ScalarProduct().TimelikeIndices(),timelikeindices);

	}
	void testCliffordIndefiniteEven() {
		ConcreteManifold M(4);		
		ex d1=M.e(1), d2=M.e(2),d3=M.e(3),d4=M.e(4);
		auto g =PseudoRiemannianStructureByOrthonormalFrame::FromTimelikeIndices(&M,M.e(),{1,3},CliffordConvention::BAUM_KATH);
		for (int i=0;i<g.DimensionOfSpinorRepresentation();++i) {
			ex u=g.u(i);
			TS_ASSERT_EQUALS(g.CliffordDot(d1,g.CliffordDot(d1,u)),u);
			TS_ASSERT_EQUALS(g.CliffordDot(d2,g.CliffordDot(d2,u)),-u);
			TS_ASSERT_EQUALS(g.CliffordDot(d3,g.CliffordDot(d3,u)),u);
			TS_ASSERT_EQUALS(g.CliffordDot(d4,g.CliffordDot(d4,u)),-u);
			for (auto X : M.e())
			for (auto Y : M.e())
				if (X!=Y) TS_ASSERT_EQUALS(g.CliffordDot(X,g.CliffordDot(Y,u))+g.CliffordDot(Y,g.CliffordDot(X,u)),0);
		}
		auto timelikeindices=vector<int>{1,3};
		TS_ASSERT_EQUALS(g.ScalarProduct().TimelikeIndices(),timelikeindices);
		}


	void testCliffordIndefiniteOdd() {
		ConcreteManifold M(3);
		Frame e=ParseDifferentialForms(M.e(),"1+2*2, 2*2, 3");
		auto dual=e.dual();
		ex d1=dual(1), d2=dual(2), d3=dual(3);
		auto g =PseudoRiemannianStructureByOrthonormalFrame::FromTimelikeIndices(&M,e,{1,3},CliffordConvention::BAUM_KATH);
		for (int i=0;i<g.DimensionOfSpinorRepresentation();++i) {
			ex u=g.u(i);
			TS_ASSERT_EQUALS(g.CliffordDot(d1,g.CliffordDot(d1,u)),u);
			TS_ASSERT_EQUALS(g.CliffordDot(d2,g.CliffordDot(d2,u)),-u);
			TS_ASSERT_EQUALS(g.CliffordDot(d3,g.CliffordDot(d3,u)),u);
			for (auto X : dual)
			for (auto Y : dual)
				if (X!=Y) TS_ASSERT_EQUALS(g.CliffordDot(X,g.CliffordDot(Y,u))+g.CliffordDot(Y,g.CliffordDot(X,u)),0);
		}
		auto timelikeindices=vector<int>{1,3};
		TS_ASSERT_EQUALS(g.ScalarProduct().TimelikeIndices(),timelikeindices);

	}

	void testDiagonalCliffordIndefiniteEven() {
		ConcreteManifold M(4);		
		ex d1=M.e(1), d2=M.e(2),d3=M.e(3),d4=M.e(4);
		auto g =PseudoRiemannianStructureByOrthogonalFrame(&M,M.e(),{-1,2,-3,4},CliffordConvention::BAUM_KATH);
		for (int i=0;i<g.DimensionOfSpinorRepresentation();++i) {
			ex u=g.u(i);
			TS_ASSERT_EQUALS(g.CliffordDot(d1,g.CliffordDot(d1,u)).expand(),u);
			TS_ASSERT_EQUALS(g.CliffordDot(d2,g.CliffordDot(d2,u)).expand(),-2*u);
			TS_ASSERT_EQUALS(g.CliffordDot(d3,g.CliffordDot(d3,u)).expand(),3*u);
			TS_ASSERT_EQUALS(g.CliffordDot(d4,g.CliffordDot(d4,u)).expand(),-4*u);
			for (auto X : M.e())
			for (auto Y : M.e())
				if (X!=Y) TS_ASSERT_EQUALS((g.CliffordDot(X,g.CliffordDot(Y,u))+g.CliffordDot(Y,g.CliffordDot(X,u))).expand(),0);
		}
		auto timelikeindices=vector<int>{1,3};
		TS_ASSERT_EQUALS(g.ScalarProduct().TimelikeIndices(),timelikeindices);
	}

	void testDiagonalCliffordIndefiniteOdd() {
		ConcreteManifold M(3);
		Frame e=ParseDifferentialForms(M.e(),"1+2*2, 2*2, 3");
		auto dual=e.dual();
		ex d1=dual(1), d2=dual(2), d3=dual(3);
		possymbol a("a");
		auto g =PseudoRiemannianStructureByOrthogonalFrame(&M,e,{-a,a,-4},CliffordConvention::BAUM_KATH);
		for (int i=0;i<g.DimensionOfSpinorRepresentation();++i) {
			ex u=g.u(i);
			TS_ASSERT_EQUALS(g.CliffordDot(d1,g.CliffordDot(d1,u)).expand(),a*u);
			TS_ASSERT_EQUALS(g.CliffordDot(d2,g.CliffordDot(d2,u)).expand(),-a*u);
			TS_ASSERT_EQUALS(g.CliffordDot(d3,g.CliffordDot(d3,u)).expand(),4*u);
			for (auto X : dual)
			for (auto Y : dual)
				if (X!=Y) TS_ASSERT_EQUALS((g.CliffordDot(X,g.CliffordDot(Y,u))+g.CliffordDot(Y,g.CliffordDot(X,u))).expand(),0);
		}
		auto timelikeindices=vector<int>{1,3};
		TS_ASSERT_EQUALS(g.ScalarProduct().TimelikeIndices(),timelikeindices);

	}


};


class SpinorTestSuite : public CxxTest::TestSuite {
public:
	void testSpinor4() {
		ConcreteManifold M(4);
		auto g =PseudoRiemannianStructureByOrthonormalFrame::FromTimelikeIndices(&M,M.e(),{},CliffordConvention::BAUM_KATH);
		TS_ASSERT_EQUALS(g.u(0),g.u({1,1}));
		TS_ASSERT_EQUALS(g.u(1),g.u({-1,1}));
		TS_ASSERT_EQUALS(g.u(2),g.u({1,-1}));
		TS_ASSERT_EQUALS(g.u(3),g.u({-1,-1}));
		TS_ASSERT_EQUALS(g.DimensionOfSpinorRepresentation(),4);
		TS_ASSERT_THROWS(g.u(-1),OutOfRange);
		TS_ASSERT_THROWS(g.u(4),OutOfRange);
		TS_ASSERT_THROWS(g.u(vector<int>{1}),InvalidArgument);
		TS_ASSERT_THROWS(g.u({1,1,1}),InvalidArgument);
	}
	void testSpinor5() {
		ConcreteManifold M(5);
		auto g =PseudoRiemannianStructureByOrthonormalFrame::FromTimelikeIndices(&M,M.e(),{},CliffordConvention::BAUM_KATH);
		TS_ASSERT_EQUALS(g.u(0),g.u({1,1}));
		TS_ASSERT_EQUALS(g.u(1),g.u({-1,1}));
		TS_ASSERT_EQUALS(g.u(2),g.u({1,-1}));
		TS_ASSERT_EQUALS(g.u(3),g.u({-1,-1}));
		TS_ASSERT_EQUALS(g.DimensionOfSpinorRepresentation(),4);
		TS_ASSERT_THROWS(g.u(-1),OutOfRange);
		TS_ASSERT_THROWS(g.u(4),OutOfRange);
		TS_ASSERT_THROWS(g.u(vector<int>{1}),InvalidArgument);
		TS_ASSERT_THROWS(g.u({1,1,1}),InvalidArgument);
	}
	void testSpinor6() {
		ConcreteManifold M(6);
		auto g =PseudoRiemannianStructureByOrthonormalFrame::FromTimelikeIndices(&M,M.e(),{},CliffordConvention::BAUM_KATH);
		TS_ASSERT_EQUALS(g.u(0),g.u({1,1,1}));
		TS_ASSERT_EQUALS(g.u(1),g.u({-1,1,1}));
		TS_ASSERT_EQUALS(g.u(2),g.u({1,-1,1}));
		TS_ASSERT_EQUALS(g.u(3),g.u({-1,-1,1}));
		TS_ASSERT_EQUALS(g.DimensionOfSpinorRepresentation(),8);
		TS_ASSERT_THROWS(g.u(-1),OutOfRange);
		TS_ASSERT_THROWS(g.u(8),OutOfRange);
		TS_ASSERT_THROWS(g.u(vector<int>{1}),InvalidArgument);
		TS_ASSERT_THROWS(g.u({1,1}),InvalidArgument);
	}

	void testVolume3() {
		ConcreteManifold M3(3), M5(5),M7(7);
		
		do_test_volume(M3);
		do_test_volume(M5);
		do_test_volume(M7);
	}
private:
	//product by a sequence of vectors
	ex clifford(exvector e, ex u, const PseudoRiemannianStructureByOrthonormalFrame& g) {
		for (auto i=e.rbegin();i!=e.rend();++i)
			u=g.CliffordDot(*i,u);
		return u;
	}
	
	void do_test_volume(const Manifold& M) {		
		for (int r=0;r<=M.Dimension();++r) {
			int n=M.Dimension();
			int s=n-r;
			vector<int> timelike;
			for (int k=1;k<=s;++k) timelike.push_back(k);
			auto g_bk =PseudoRiemannianStructureByOrthonormalFrame::FromTimelikeIndices(&M,M.e(),timelike,CliffordConvention::BAUM_KATH);
			for (auto u : {g_bk.u(0),g_bk.u(1)})
				TS_ASSERT_EQUALS(clifford(M.e(),u,g_bk),pow(I,(r-s+1)/2)*u);			
			
			auto g_std =PseudoRiemannianStructureByOrthonormalFrame::FromTimelikeIndices(&M,M.e(),timelike,CliffordConvention::STANDARD);
			for (auto u : {g_std.u(0),g_std.u(1)}) 
				TS_ASSERT_EQUALS(clifford(M.e(),u,g_std),pow(I,(s+(n*(n+1))/2))*u);			
		}
	}
};


#endif /* SRC_STRUCTURES_H_ */
