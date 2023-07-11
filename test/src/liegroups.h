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

#ifndef SRC_LIEGROUPS_H_
#define SRC_LIEGROUPS_H_


#include <cxxtest/TestSuite.h>
#include "test.h"
#include "wedge/liealgebras/liegroup.h"
#include "wedge/liealgebras/liegroupextension.h"
#include "wedge/liealgebras/liesubgroup.h"
#include "wedge/liealgebras/so.h"
#include "wedge/liealgebras/su.h"
#include "wedge/liealgebras/liegrouptostring.h"
#include "wedge/liealgebras/derivations.h"
#include "wedge/polynomialalgebra/cocoapolyalg.h"
#include "wedge/polynomialalgebra/polybasis.h"
#include "wedge/representations/gl.h"

//test liegroup.h, subgroup.h, gl.h, so.h, su.h
class LieGroupTestSuite : public CxxTest::TestSuite  {
public:
	void testStructureConstant() {
		ex x=StructureConstant{N.x};
		TS_ASSERT(is_a<symbol>(x));
		TS_ASSERT(x.info(info_flags::symbol));
	}

	void testAbstractGroup2() {
		AbstractLieGroup<> G("0,0,[sqrt(3)]*12,13");
		VectorSpace<DifferentialForm> twoforms=G.pForms(2);
		exvector eqns;
		GetCoefficients<DifferentialForm>(eqns,G.d(twoforms.GenericElement()));
		exvector sol;
		twoforms.GetSolutions(sol,eqns.begin(),eqns.end());
		TS_ASSERT_EQUALS(sol.size(),4);
		Subspace<DifferentialForm> V=twoforms.Subspace(sol.begin(),sol.end());
		TS_ASSERT_EQUALS(V.Dimension(),4);
		for (int i=1;i<=V.Dimension();++i)
		{
			TS_ASSERT_EQUALS(G.d(V.e(i)),0);
		}
	}
	void testAbstractGroupParameters() {
		ex a=symbol("a");
		auto G=AbstractLieGroup<true>{"0,0,[a]*12,[b1+c]*12+[b2]*13",a,NameRange(N.b,1,3),N.c};
		stringstream s;
		G.canonical_print(s);
		cout<<s.str();
		TS_ASSERT_EQUALS(s.str(),"(0,0,a*(e1*e2),(b1+c)*(e1*e2)+b2*(e1*e3))")
	}
	void testAbstractGroupParametersLst() {
		symbol a("a");
		auto G=AbstractLieGroup<true>{"0,0,[a]*12,[b1+c]*12+[b2]*13",N.c,lst{a},NameRange(N.b,1,3)};
		stringstream s;
		G.canonical_print(s);
		cout<<s.str();
		TS_ASSERT_EQUALS(s.str(),"(0,0,a*(e1*e2),(b1+c)*(e1*e2)+b2*(e1*e3))")
	}
	void testAbstractGroup() {
		AbstractLieGroup<> SU3("-23-45+2*67,13+46-57-[sqrt(3)]*58,-12-47+[sqrt(3)]*48-56,"
					"15-26+37-[sqrt(3)]*38,"
					"-14+27+36+[sqrt(3)]*28,-2*17+24-35,2*16-25-34,-[sqrt(3)]*25+[sqrt(3)]*34");
		TS_ASSERT_EQUALS(Hook(SU3.e(2),SU3.d(SU3.e(1))),-SU3.e(3));
		TS_ASSERT_EQUALS(Hook(SU3.e(2),SU3.d(SU3.e(2))),0);
		TS_ASSERT_EQUALS(Hook(SU3.e(2),SU3.d(SU3.e(3))),SU3.e(1));
		TS_ASSERT_EQUALS(Hook(SU3.e(2),SU3.d(SU3.e(4))),-SU3.e(6));
		TS_ASSERT_EQUALS(Hook(SU3.e(2),SU3.d(SU3.e(5))),SU3.e(7)+sqrt(ex(3))*SU3.e(8));
		TS_ASSERT_EQUALS(Hook(SU3.e(2),SU3.d(SU3.e(6))),SU3.e(4));
		TS_ASSERT_EQUALS(Hook(SU3.e(2),SU3.d(SU3.e(7))),-SU3.e(5));
		TS_ASSERT_EQUALS(Hook(SU3.e(2),SU3.d(SU3.e(8))),-sqrt(ex(3))*SU3.e(5));

		TS_ASSERT_EQUALS(SU3.LieBracket(SU3.e(1),SU3.e(2)),SU3.e(3));

		Subspace<DifferentialForm> Closed2forms=SU3.ClosedForms(2);
		for (int i=1;i<=Closed2forms.Dimension();++i)
			TS_ASSERT_EQUALS(SU3.d(Closed2forms.e(i)),0);
		TS_ASSERT_EQUALS(SU3.d(Closed2forms.GenericElement()),0);
		TS_ASSERT_EQUALS(Closed2forms.Dimension(),8);
		LOG_INFO(Closed2forms);
		TS_ASSERT_EQUALS(SU3.ExactForms(2),SU3.ClosedForms(2));
		TS_ASSERT(SU3.IsUnimodular());
	}

	void testCanonicalPrint() {
		AbstractLieGroup<> G("0,0,[sqrt(3)]*12,13+23");
		stringstream s;
		G.canonical_print(s);
		TS_ASSERT_EQUALS(s.str(),"(0,0,sqrt(3)*(e1*e2),e1*e3+e2*e3)")
		s.str("");
		s<<latex;
		G.canonical_print(s);
		TS_ASSERT_EQUALS(s.str(),"(0,0,\\sqrt{3} e^{12},e^{13}+e^{23})")

		AbstractLieGroup<> H("0,0,-12+[sqrt(3)]*12");
		s.str("");
		s<<dflt;
		H.canonical_print(s);
		TS_ASSERT_EQUALS(s.str(),"(0,0,(-1+sqrt(3))*(e1*e2))")
		s.str("");
		s<<latex;
		H.canonical_print(s);
		TS_ASSERT_EQUALS(s.str(),"(0,0,(-1+\\sqrt{3}) e^{12})")
	}


	void testGeneric() {

		GenericLieGroup G(3);
		PolyBasisImplementation<StructureConstant,DefaultPolyAlgorithms> pols;
		G.GetEquations_ddZero(pols);
		pols.Reduce();

		exvector h;
		h.push_back(G.e(1)+G.e(2));
		{
			LieSubgroup<true> H(G,h);
			H.Check_ddZero();
		}
		{
			AbstractLieSubgroup<true> H(G,h);
			H.Check_ddZero();
		}

		h.push_back(G.e(2)+G.e(3));
		{
			LieSubgroup<true> H2(G,h);
			TS_ASSERT_EQUALS(H2.LieBracket(H2.e(1),H2.e(1)),0);
			TS_ASSERT_EQUALS(H2.LieBracket(H2.e(2),H2.e(2)),0);
			TS_ASSERT(H2.pForms(1).Contains(H2.LieBracket(H2.e(1),H2.e(2))));
			TS_ASSERT(H2.pForms(2).Contains(H2.d(H2.e(1))));
			TS_ASSERT(H2.pForms(2).Contains(H2.d(H2.e(2))));

			PolyBasis<StructureConstant> pols2;
			H2.GetEquations_ddZero(pols2);
			for (PolyBasis<StructureConstant>::const_iterator i=pols2.begin();i!=pols2.end();i++)
				TS_ASSERT(pols.IdealContains(*i));
		}
		TS_ASSERT_THROWS(AbstractLieSubgroup<true> H2(G,h),WedgeException<std::runtime_error>);

	}

	void testSubgroup() {
		AbstractLieGroup<> G("23,31,12,0,0");
		TS_ASSERT_EQUALS(G.d(G.e(1)),G.e(2)*G.e(3));
		TS_ASSERT_EQUALS(G.d(G.e(2)),G.e(3)*G.e(1));
		TS_ASSERT_EQUALS(G.d(G.e(3)),G.e(1)*G.e(2));
		TS_ASSERT_EQUALS(G.d(G.e(4)),0);
		TS_ASSERT_EQUALS(G.d(G.e(5)),0);
		TS_ASSERT_EQUALS(G.LieBracket(G.e(1),G.e(2)),-G.e(3));
		TS_ASSERT_EQUALS(G.LieBracket(G.e(1)+2*G.e(3),G.e(2)),2*G.e(1)-G.e(3));

		exvector h; h.push_back(G.e(1)+G.e(2)); h.push_back(G.e(4));
		LieSubgroup<false> H(G, h);
		TS_ASSERT_EQUALS(H.LieBracket(H.e(1),H.e(2)),0);
		TS_ASSERT_EQUALS(H.LieBracket(H.e(1),H.e(1)),0);
		TS_ASSERT_EQUALS(H.LieBracket(H.e(2),H.e(2)),0);
		TS_ASSERT_EQUALS(H.Dimension(),2);
		TS_ASSERT_EQUALS(H.d(H.e(1)),0);
		TS_ASSERT_EQUALS(H.d(H.e(2)),0);
		TS_ASSERT(H.IsUnimodular());

		h.push_back(G.e(2)); h.push_back(G.e(3));
		LieSubgroup<false> H2(G, h);
		TS_ASSERT_EQUALS(H2.Dimension(),4);
		TS_ASSERT(H2.IsUnimodular());
		TS_ASSERT_EQUALS(H2.LieBracket(H2.e(1),H2.e(3)),-H2.e(4));
		TS_ASSERT_EQUALS(H2.LieBracket(H2.e(1),H2.e(4)),2*H2.e(3)-H2.e(1));

		vector<int> betti=H2.BettiNumbers();
		TS_ASSERT_EQUALS(betti.size(),5);
		TS_ASSERT_EQUALS(betti[0],1);
		TS_ASSERT_EQUALS(betti[1],1);
		TS_ASSERT_EQUALS(betti[2],0);
		TS_ASSERT_EQUALS(betti[3],1);
		TS_ASSERT_EQUALS(betti[4],1);

		betti=G.BettiNumbers();
		TS_ASSERT_EQUALS(betti[2],1);
		TS_ASSERT_EQUALS(betti[3],1);
		TS_ASSERT_EQUALS(betti[4],2);
		TS_ASSERT_EQUALS(betti[5],1);
	}

	void testAbstractSubgroup() {
		AbstractLieGroup<> G("23,31,12,0,0");
		TS_ASSERT_EQUALS(G.d(G.e(1)),G.e(2)*G.e(3));
		TS_ASSERT_EQUALS(G.d(G.e(2)),G.e(3)*G.e(1));
		TS_ASSERT_EQUALS(G.d(G.e(3)),G.e(1)*G.e(2));
		TS_ASSERT_EQUALS(G.d(G.e(4)),0);
		TS_ASSERT_EQUALS(G.d(G.e(5)),0);
		TS_ASSERT_EQUALS(G.LieBracket(G.e(1),G.e(2)),-G.e(3));
		TS_ASSERT_EQUALS(G.LieBracket(G.e(1)+2*G.e(3),G.e(2)),2*G.e(1)-G.e(3));

		exvector h; h.push_back(G.e(1)+G.e(2)); h.push_back(G.e(4));
		AbstractLieSubgroup<false> H(G, h);
		TS_ASSERT_EQUALS(H.LieBracket(H.e(1),H.e(2)),0);
		TS_ASSERT_EQUALS(H.LieBracket(H.e(1),H.e(1)),0);
		TS_ASSERT_EQUALS(H.LieBracket(H.e(2),H.e(2)),0);
		TS_ASSERT_EQUALS(H.Dimension(),2);
		TS_ASSERT_EQUALS(H.d(H.e(1)),0);
		TS_ASSERT_EQUALS(H.d(H.e(2)),0);
		TS_ASSERT(H.IsUnimodular());

		h.push_back(G.e(2)); h.push_back(G.e(3));
		AbstractLieSubgroup<false> H2(G, h);
		TS_ASSERT_EQUALS(H2.Dimension(),4);
		TS_ASSERT(H2.IsUnimodular());
		TS_ASSERT_EQUALS(H2.LieBracket(H2.e(1),H2.e(3)),-H2.e(4));
		TS_ASSERT_EQUALS(H2.LieBracket(H2.e(1),H2.e(4)),2*H2.e(3)-H2.e(1));

		vector<int> betti=H2.BettiNumbers();
		TS_ASSERT_EQUALS(betti.size(),5);
		TS_ASSERT_EQUALS(betti[0],1);
		TS_ASSERT_EQUALS(betti[1],1);
		TS_ASSERT_EQUALS(betti[2],0);
		TS_ASSERT_EQUALS(betti[3],1);
		TS_ASSERT_EQUALS(betti[4],1);

		betti=G.BettiNumbers();
		TS_ASSERT_EQUALS(betti[2],1);
		TS_ASSERT_EQUALS(betti[3],1);
		TS_ASSERT_EQUALS(betti[4],2);
		TS_ASSERT_EQUALS(betti[5],1);
	}


	//test so.h
	void testSO()
	{
		SO G(3);
		TS_ASSERT_EQUALS(G.Dimension(),3);
		TS_ASSERT_EQUALS(G.A(1,1),0);
		TS_ASSERT_EQUALS(G.A(2,2),0);
		TS_ASSERT_EQUALS(G.A(3,3),0);
		TS_ASSERT_EQUALS(G.A(1,2)+G.A(2,1),0);
		TS_ASSERT_EQUALS(G.LieBracket(G.A(1,2),G.A(1,3)),-G.A(2,3));
		TS_ASSERT_EQUALS(G.LieBracket(G.e(1),G.e(2)),-G.e(3));
		TS_ASSERT_EQUALS(G.LieBracket(G.e(3),G.e(1)),-G.e(2));
		TS_ASSERT_EQUALS(G.LieBracket(G.e(2),G.e(3)),-G.e(1));
	}

	static int delta(int i, int j) {return i==j? 1 : 0;}
	//test su.h
	void testSU()
	{
		{
			SU G(2);
			TS_ASSERT_EQUALS(G.Dimension(),3);
			TS_ASSERT_EQUALS(G.e(1),G.E(1,2));
			TS_ASSERT_EQUALS(G.e(2),G.F(1,2));
			TS_ASSERT_EQUALS(G.e(3),G.G(1));
			TS_ASSERT_EQUALS(G.LieBracket(G.e(1),G.e(2)),2*G.e(3));
			TS_ASSERT_EQUALS(G.LieBracket(G.e(1),G.e(3)),-2*G.e(2));
			TS_ASSERT_EQUALS(G.LieBracket(G.e(2),G.e(3)),2*G.e(1));

			TS_ASSERT_EQUALS(G.d(G.e(1)),-2*G.e(2)*G.e(3));
			TS_ASSERT_EQUALS(G.d(G.e(2)),2*G.e(1)*G.e(3));
			TS_ASSERT_EQUALS(G.d(G.e(3)),-2*G.e(1)*G.e(2));
		}

		{
			SU G(3);
			TS_ASSERT_EQUALS(G.Dimension(),8);
			TS_ASSERT_EQUALS(G.e(1),G.E(1,2));
			TS_ASSERT_EQUALS(G.e(2),G.E(1,3));
			TS_ASSERT_EQUALS(G.e(3),G.E(2,3));
			TS_ASSERT_EQUALS(G.e(4),G.F(1,2));
			TS_ASSERT_EQUALS(G.e(5),G.F(1,3));
			TS_ASSERT_EQUALS(G.e(6),G.F(2,3));
			TS_ASSERT_EQUALS(G.e(7),G.G(1));
			TS_ASSERT_EQUALS(G.e(8),G.G(2));
			for (int i=1;i<=3;++i)
			for (int j=i+1;j<=3;++j)
			for (int k=1;k<=3;++k)
			for (int l=k+1;l<=3;++l)
				TS_ASSERT_EQUALS(G.LieBracket(G.E(i,j),G.E(k,l)),-delta(i,l)*G.E(k,j)+delta(j,l)*G.E(k,i)+delta(i,k)*G.E(l,j)-delta(j,k)*G.E(l,i));

			TS_ASSERT_EQUALS(G.d(G.E(1,2)), -G.E(1,3)*G.E(3,2)+G.F(1,3)*G.F(3,2)-(G.F(1,2)*(2*G.G(1)-G.G(2))).expand());
			TS_ASSERT_EQUALS(G.d(G.E(1,3)), -G.E(1,2)*G.E(2,3)+G.F(1,2)*G.F(2,3)-(G.F(1,3)*(G.G(1)+G.G(2))).expand());
			TS_ASSERT_EQUALS(G.d(G.E(2,3)), -G.E(2,1)*G.E(1,3)+G.F(2,1)*G.F(1,3)-(G.F(2,3)*(2*G.G(2)-G.G(1))).expand());

			//verify invariance of ThreeForm under change of frame
			Frame l=ParseDifferentialForms(G.e(),"7,1,4,2,5,3,6,[sqrt(1/3)]*(7+2*8)");	//the Gambioli frame
			AbstractLieSubgroup<false> H(G,l);
			lst subs=LinearMapToSubstitutions<DifferentialForm>(l.dual(),H.e());
			TS_ASSERT_EQUALS(H.ThreeForm(),G.ThreeForm().subs(subs));
		}
	}



	void testExtension()
	{
		AbstractLieGroup<> G("23,31,12");
		TS_ASSERT_THROWS(LieGroupExtension H(G,2),InvalidArgument);
		LieGroupExtension H(G,3);
		H.Check_ddZero();

		LieGroupExtension H2(G,4);
		H2.Declare_d(H2.e(4),0);
		PolyBasisImplementation<StructureConstant, DefaultPolyAlgorithms> pol;
		H2.GetEquations_ddZero(pol);
		LieGroupExtension H3(H2,5);
		PolyBasisImplementation<StructureConstant, DefaultPolyAlgorithms> pol2;
		H3.GetEquations_ddZero(pol2);
		//for (PolyBasis_impl<StructureConstant, DefaultPolyAlgorithms>::const_iterator i=pol.begin();i!=pol.end();i++)
//			TS_ASSERT(pol2.IdealContains(*i));


	}

	void testDerivations() {
		AbstractLieGroup<> G("23,31,12");
		GL gl(3);
		auto der= derivations(G,gl);
		TS_ASSERT_EQUALS(der.Dimension(),3);
		TS_ASSERT(der.Contains(gl.A(1,2)-gl.A(2,1)));
		TS_ASSERT(der.Contains(gl.A(3,2)-gl.A(2,3)));
		TS_ASSERT(der.Contains(gl.A(1,3)-gl.A(3,1)));
	}

	void test_lie_group_to_string() {
		AbstractLieGroup<> G("23,31,12");
		TS_ASSERT_EQUALS(lie_group_to_string(G),"23,-13,12");
		AbstractLieGroup<> G2("0,12,13,14,15,16,17,18,19,1a");
		TS_ASSERT_EQUALS(lie_group_to_string(G2),"0,12,13,14,15,16,17,18,19,1a");
		AbstractLieGroup<> G3("0,-12,-2*13+12");
		TS_ASSERT_EQUALS(lie_group_to_string(G3),"0,-12,12-2*13");
		AbstractLieGroup<true> G4("0,[a]*12,[sqrt(2)]*13",N.a);
		TS_ASSERT_EQUALS(lie_group_to_string(G4),"0,[a]*12,[sqrt(2)]*13");
	}
};





#endif /* SRC_LIEGROUPS_H_ */
