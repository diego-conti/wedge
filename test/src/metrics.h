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



	void testCliffordRiemannian() {
		ConcreteManifold M(3);
		Frame e=ParseDifferentialForms(M.e(),"1+2*2, 2*2, 3");
		auto dual=e.dual();
		ex d1=dual(1), d2=dual(2), d3=dual(3);
		auto g =PseudoRiemannianStructureByOrthonormalFrame::FromTimelikeIndices(&M,e,{});
		for (int i=0;i<g.DimensionOfSpinorRepresentation();++i) {
			ex u=g.u(i);
			TS_ASSERT_EQUALS(g.CliffordDot(d1,g.CliffordDot(d1,u)),-u);
			TS_ASSERT_EQUALS(g.CliffordDot(d2,g.CliffordDot(d2,u)),-u);
			TS_ASSERT_EQUALS(g.CliffordDot(d3,g.CliffordDot(d3,u)),-u);
			TS_ASSERT_EQUALS(g.CliffordDot(d1,g.CliffordDot(d2,u))+g.CliffordDot(d2,g.CliffordDot(d1,u)),0);
			TS_ASSERT_EQUALS(g.CliffordDot(d1,g.CliffordDot(d3,u))+g.CliffordDot(d3,g.CliffordDot(d1,u)),0);
			TS_ASSERT_EQUALS(g.CliffordDot(d2,g.CliffordDot(d3,u))+g.CliffordDot(d3,g.CliffordDot(d2,u)),0);
		}
		auto timelikeindices=vector<int>{};
		TS_ASSERT_EQUALS(g.ScalarProduct().TimelikeIndices(),timelikeindices);

	}
	void testCliffordRiemannianEven() {
		ConcreteManifold M(4);		
		ex d1=M.e(1), d2=M.e(2),d3=M.e(3),d4=M.e(4);
		auto g =PseudoRiemannianStructureByOrthonormalFrame::FromTimelikeIndices(&M,M.e(),{});
		for (int i=0;i<g.DimensionOfSpinorRepresentation();++i) {
			ex u=g.u(i);
			TS_ASSERT_EQUALS(g.CliffordDot(d1,g.CliffordDot(d1,u)),-u);
			TS_ASSERT_EQUALS(g.CliffordDot(d2,g.CliffordDot(d2,u)),-u);
			TS_ASSERT_EQUALS(g.CliffordDot(d3,g.CliffordDot(d3,u)),-u);
			TS_ASSERT_EQUALS(g.CliffordDot(d4,g.CliffordDot(d4,u)),-u);
			TS_ASSERT_EQUALS(g.CliffordDot(d1,g.CliffordDot(d2,u))+g.CliffordDot(d2,g.CliffordDot(d1,u)),0);
			TS_ASSERT_EQUALS(g.CliffordDot(d1,g.CliffordDot(d3,u))+g.CliffordDot(d3,g.CliffordDot(d1,u)),0);
			TS_ASSERT_EQUALS(g.CliffordDot(d2,g.CliffordDot(d3,u))+g.CliffordDot(d3,g.CliffordDot(d2,u)),0);
		}
		auto timelikeindices=vector<int>{};
		TS_ASSERT_EQUALS(g.ScalarProduct().TimelikeIndices(),timelikeindices);

	}
	void testCliffordIndefiniteEven() {
		ConcreteManifold M(4);		
		ex d1=M.e(1), d2=M.e(2),d3=M.e(3),d4=M.e(4);
		auto g =PseudoRiemannianStructureByOrthonormalFrame::FromTimelikeIndices(&M,M.e(),{1,3});
		for (int i=0;i<g.DimensionOfSpinorRepresentation();++i) {
			ex u=g.u(i);
			TS_ASSERT_EQUALS(g.CliffordDot(d1,g.CliffordDot(d1,u)),u);
			TS_ASSERT_EQUALS(g.CliffordDot(d2,g.CliffordDot(d2,u)),-u);
			TS_ASSERT_EQUALS(g.CliffordDot(d3,g.CliffordDot(d3,u)),u);
			TS_ASSERT_EQUALS(g.CliffordDot(d4,g.CliffordDot(d4,u)),-u);
			TS_ASSERT_EQUALS(g.CliffordDot(d1,g.CliffordDot(d2,u))+g.CliffordDot(d2,g.CliffordDot(d1,u)),0);
			TS_ASSERT_EQUALS(g.CliffordDot(d1,g.CliffordDot(d3,u))+g.CliffordDot(d3,g.CliffordDot(d1,u)),0);
			TS_ASSERT_EQUALS(g.CliffordDot(d2,g.CliffordDot(d3,u))+g.CliffordDot(d3,g.CliffordDot(d2,u)),0);
		}
		auto timelikeindices=vector<int>{1,3};
		TS_ASSERT_EQUALS(g.ScalarProduct().TimelikeIndices(),timelikeindices);

	}


	void testCliffordIndefiniteOdd() {
		ConcreteManifold M(3);
		Frame e=ParseDifferentialForms(M.e(),"1+2*2, 2*2, 3");
		auto dual=e.dual();
		ex d1=dual(1), d2=dual(2), d3=dual(3);
		auto g =PseudoRiemannianStructureByOrthonormalFrame::FromTimelikeIndices(&M,e,{1,3});
		for (int i=0;i<g.DimensionOfSpinorRepresentation();++i) {
			ex u=g.u(i);
			TS_ASSERT_EQUALS(g.CliffordDot(d1,g.CliffordDot(d1,u)),u);
			TS_ASSERT_EQUALS(g.CliffordDot(d2,g.CliffordDot(d2,u)),-u);
			TS_ASSERT_EQUALS(g.CliffordDot(d3,g.CliffordDot(d3,u)),u);
			TS_ASSERT_EQUALS(g.CliffordDot(d1,g.CliffordDot(d2,u))+g.CliffordDot(d2,g.CliffordDot(d1,u)),0);
			TS_ASSERT_EQUALS(g.CliffordDot(d1,g.CliffordDot(d3,u))+g.CliffordDot(d3,g.CliffordDot(d1,u)),0);
			TS_ASSERT_EQUALS(g.CliffordDot(d2,g.CliffordDot(d3,u))+g.CliffordDot(d3,g.CliffordDot(d2,u)),0);
		}
		auto timelikeindices=vector<int>{1,3};
		TS_ASSERT_EQUALS(g.ScalarProduct().TimelikeIndices(),timelikeindices);

	}

};


class SpinorTestSuite : public CxxTest::TestSuite {
public:
	void testSpinor4() {
		ConcreteManifold M(4);
		auto g =PseudoRiemannianStructureByOrthonormalFrame::FromTimelikeIndices(&M,M.e(),{});
		TS_ASSERT_EQUALS(g.u(0),g.u({1,1}));
		TS_ASSERT_EQUALS(g.u(1),g.u({1,-1}));
		TS_ASSERT_EQUALS(g.u(2),g.u({-1,1}));
		TS_ASSERT_EQUALS(g.u(3),g.u({-1,-1}));
		TS_ASSERT_EQUALS(g.DimensionOfSpinorRepresentation(),4);
		TS_ASSERT_THROWS(g.u(-1),OutOfRange);
		TS_ASSERT_THROWS(g.u(4),OutOfRange);
		TS_ASSERT_THROWS(g.u(vector<int>{1}),InvalidArgument);
		TS_ASSERT_THROWS(g.u({1,1,1}),InvalidArgument);
	}
	void testSpinor5() {
		ConcreteManifold M(5);
		auto g =PseudoRiemannianStructureByOrthonormalFrame::FromTimelikeIndices(&M,M.e(),{});
		TS_ASSERT_EQUALS(g.u(0),g.u({1,1}));
		TS_ASSERT_EQUALS(g.u(1),g.u({1,-1}));
		TS_ASSERT_EQUALS(g.u(2),g.u({-1,1}));
		TS_ASSERT_EQUALS(g.u(3),g.u({-1,-1}));
		TS_ASSERT_EQUALS(g.DimensionOfSpinorRepresentation(),4);
		TS_ASSERT_THROWS(g.u(-1),OutOfRange);
		TS_ASSERT_THROWS(g.u(4),OutOfRange);
		TS_ASSERT_THROWS(g.u(vector<int>{1}),InvalidArgument);
		TS_ASSERT_THROWS(g.u({1,1,1}),InvalidArgument);
	}
	void testSpinor6() {
		ConcreteManifold M(6);
		auto g =PseudoRiemannianStructureByOrthonormalFrame::FromTimelikeIndices(&M,M.e(),{});
		TS_ASSERT_EQUALS(g.u(0),g.u({1,1,1}));
		TS_ASSERT_EQUALS(g.u(1),g.u({1,1,-1}));
		TS_ASSERT_EQUALS(g.u(2),g.u({1,-1,1}));
		TS_ASSERT_EQUALS(g.u(3),g.u({1,-1,-1}));
		TS_ASSERT_EQUALS(g.DimensionOfSpinorRepresentation(),8);
		TS_ASSERT_THROWS(g.u(-1),OutOfRange);
		TS_ASSERT_THROWS(g.u(8),OutOfRange);
		TS_ASSERT_THROWS(g.u(vector<int>{1}),InvalidArgument);
		TS_ASSERT_THROWS(g.u({1,1}),InvalidArgument);
	}
};


#endif /* SRC_STRUCTURES_H_ */