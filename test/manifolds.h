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
#ifndef MANIFOLDS_
#define MANIFOLDS_

#include <list>
#include <vector>
#include <cxxtest/TestSuite.h>
#include <ginac/ginac.h>
#include "test.h"
#include "wedge/riemannianstructure.h"
#include "wedge/connection.h"
#include "wedge/liegroup.h"
#include "wedge/liegroupextension.h"
#include "wedge/gstructure.h"
#include "wedge/parse.h"
#include "wedge/structures.h"
#include "wedge/function.h"
#include "wedge/coordinates.h"
#include "wedge/polybasis.h"
#include "wedge/liesubgroup.h"
#include "wedge/gl.h"
#include "wedge/su.h"
#include "wedge/so.h"
#include "wedge/manifoldwith.h"
#include "wedge/liederivative.h"
#include "wedge/normalform.h"
#ifdef HAVE_COCOA
#include "wedge/cocoapolyalg.h"
#endif

using namespace GiNaC;
using namespace Wedge;
/**
	@author Diego Conti <diego.conti@unimib.it>
*/


//test manifold.h, differentialform.h, concretemanifold.h
class DifferentialFormsTestSuite : public CxxTest::TestSuite  {
public:
	class TestManifold : public Manifold 
	{
		Frame frame;
		exvector CreateFrame() 
		{
			exvector result;
			result.push_back(DifferentialOneForm(N.x));			
			result.push_back(DifferentialOneForm(N.y)+result[0]);			
			result.push_back(DifferentialOneForm(N.z)-result[1]);			
			result.push_back(DifferentialOneForm(N.t)+result[1]);
			return result;
		}
	public:
		TestManifold() : frame(CreateFrame()) {
			Check_ddZero();
		}
		ex d(ex alpha) const {return 0;}
		const Frame& e() const {return frame;}
	};

	class TestConcreteManifold  : public ConcreteManifold, public virtual Has_dTable 
	{
	public:
		TestConcreteManifold() : ConcreteManifold(3){
			Declare_d(e(1),e(2)*e(3));				  
			Declare_d(e(2),e(1)*e(2));
			Declare_d(e(3),0);
			TS_ASSERT_THROWS(Check_ddZero(),WedgeException<runtime_error>);
			Declare_d(e(2),e(1)*e(3));
			Check_ddZero();
			  
			TS_ASSERT_EQUALS(d(e(1)+e(2)),e(2)*e(3)+e(1)*e(3));
			TS_ASSERT_EQUALS(d(e(1)-e(2)),e(2)*e(3)-e(1)*e(3));
			TS_ASSERT_EQUALS(d(e(1)*e(2)),0);
			  
			VectorSpace<DifferentialForm> V=pForms(1);
			TS_ASSERT_EQUALS(V.Dimension(),3);
			TS_ASSERT(V.Contains(e(1)));
			TS_ASSERT(V.Contains(e(2)));
			TS_ASSERT(V.Contains(e(3)));
			V=pForms(2);
			TS_ASSERT_EQUALS(V.Dimension(),3);
			TS_ASSERT(V.Contains(e(1)*e(2)));
			TS_ASSERT(V.Contains(e(2)*e(3)));
			TS_ASSERT(V.Contains(e(3)*e(1)));
			V=pForms(3);
			TS_ASSERT_EQUALS(V.Dimension(),1);
			TS_ASSERT(V.Contains(e(1)*e(2)*e(3)));
			
			TS_ASSERT_EQUALS(HodgeStar(e(1)),e(2)*e(3));
			TS_ASSERT_EQUALS(d(e(1)-e(2)),e(2)*e(3)-e(1)*e(3));
			TS_ASSERT_EQUALS(d(e(1)*e(2)),0);

// FIXME this is a known bug: HodgeStar only works for positive degree forms at the moment
			TS_ASSERT_EQUALS(HodgeStar(1),e(1)*e(2)*e(3));


	  }
  };




	void testManifold()
	{
		TestManifold M;
		ex volume=1;
		for (Frame::const_iterator i=M.e().begin();i!=M.e().end();i++)
			volume*=*i;
		TS_ASSERT_EQUALS(M.HodgeStar(volume),1);

		VectorSpace<DifferentialForm> oneforms=M.pForms(1);		
		VectorSpace<DifferentialForm> twoforms=M.pForms(2);
		VectorSpace<DifferentialForm> threeforms=M.pForms(3);
		TS_ASSERT_EQUALS(M.Dimension(),4);
		TS_ASSERT_EQUALS(oneforms.Dimension(),4);
		TS_ASSERT_EQUALS(twoforms.Dimension(),6);
		TS_ASSERT_EQUALS(threeforms.Dimension(),4);
		
		for (IBasis<DifferentialForm>::const_iterator i=twoforms.e().begin();i!=twoforms.e().end();i++)
		for (IBasis<DifferentialForm>::const_iterator j=twoforms.e().begin();j!=twoforms.e().end();j++)
		{
			ex x=*i*M.HodgeStar(*j);
			TS_ASSERT_EQUALS(x.expand(),TrivialPairing<DifferentialForm>(*i,*j)*volume);
		}

		for (IBasis<DifferentialForm>::const_iterator i=oneforms.e().begin();i!=oneforms.e().end();i++)
		for (IBasis<DifferentialForm>::const_iterator j=oneforms.e().begin();j!=oneforms.e().end();j++)		
			TS_ASSERT_EQUALS((*i*M.HodgeStar(*j)).expand(),TrivialPairing<DifferentialForm>(*i,*j)*volume);
		
		//verify that a Manifold can correctly determine that d of a constant is zero		
		TS_ASSERT_EQUALS(M.d(0),0);
		TS_ASSERT_EQUALS(M.d(-2),0);
		TS_ASSERT_EQUALS(M.d(-sqrt(ex(3))),0);
		TS_ASSERT_EQUALS(M.d(-Pi+sin(ex(1))),0);
		TestConcreteManifold N;
	}
	

	void testForms()
	{
		DifferentialOneForm a,b,c;
		TS_ASSERT(a.is_equal(a));
		TS_ASSERT_EQUALS(a,a);
		ex d=a*(b+c)-a*b;
		TS_ASSERT(is_a<DifferentialForm>(Hook(a,a*b*c)));
		TS_ASSERT(is_a<DifferentialForm>(Hook(c,a*b*c)));
		TS_ASSERT_EQUALS(d.expand(),-c*a);		
		TS_ASSERT_EQUALS(Hook(a,d),c);
		TS_ASSERT_EQUALS(Hook(c,d),-a);
		TS_ASSERT_EQUALS(Hook(c*b,c*b),1);
		TS_ASSERT_EQUALS(Hook(c*b,a*c*b),a);
	}
	
	void testFrame() {
		ConcreteManifold M(5);
		ExVector e=ExVector(M.e());
		for (int i=0;i<e.size();i++)
				for (int j=i+1;i<e.size();i++)
					e[i]+=j*e[j];
		Frame frame(e);
		
		for (int i=0;i<frame.size();i++)
			for (int j=0;j<frame.size();j++) {
				TS_ASSERT_EQUALS(Hook(frame.dual()[i],frame[j]),i==j ? 1 : 0);
			}
// FIXME this is a known bug: Hook only works for positive degree forms at the moment
		for (int i=0;i<frame.size();i++)		
			TS_ASSERT_EQUALS(Hook(1,frame[i]),frame[i]);
		
		for (int i=0;i<frame.size();i++)
		{
			TS_ASSERT_EQUALS(frame[i],e[i]);
			exvector comps=frame.Components(frame[i]);
			for (int j=0;j<frame.size();j++)
				TS_ASSERT_EQUALS(comps[j],i==j? 1:0);

			comps=frame.Components(M.e()[i]);
			ex sum=0;
			for (int j=0;j<comps.size();j++)
				sum+=comps[j]*frame[j];
			TS_ASSERT_EQUALS(M.e()[i],sum);

		}
	}
	
	void testExpand() {
	  DifferentialOneForm a,b,c,d;
    ex p=2*a*b, q=2*c*d;
    TS_ASSERT_EQUALS((p*q).expand(), 4*a*b*c*d);
    symbol x;
    p=2*a*b*x;
    TS_ASSERT_EQUALS((p*q).expand(), 4*x*a*b*c*d);
    
	}
};

//test function.h, coordinates.h, liederivative.h
class FunctionTestSuite : public CxxTest::TestSuite  {
	class TwoSphere : public ManifoldWithCoordinates {		
	public:	
		ExVector OrthFrame()
		{			
			ExVector e(2);
			e(1)=2/(pow(x(1),2)+pow(x(2),2)+1)*d(x(1));
			e(2)=2/(pow(x(1),2)+pow(x(2),2)+1)*d(x(2));
			LOG_INFO(e);
			Frame x(e.begin(),e.end());
			LOG_INFO(x.dual());
			return e;
		}

		RiemannianStructure g;	
		LeviCivitaConnection<true> omega;
		
		TwoSphere() : ManifoldWithCoordinates(2), g(this,OrthFrame()), omega(this,g)
		{			
		}
	};
	
	class S3ByS2 : public ManifoldWithCoordinates {
		exvector CreateFrame()
		{			
			exvector e(3);
			for (int i=0;i<3;++i) e[i]=DifferentialOneForm(N.e(i+1));
			return e;
		}
		exvector OrthFrame()	//riemannian product of S^2 by S^3
		{			
			ExVector g=e();
			g(4)=2/(pow(x(1),2)+pow(x(2),2)+1)*d(x(1));
			g(5)=2/(pow(x(1),2)+pow(x(2),2)+1)*d(x(2));
			return g;
		}		
		RiemannianStructure g;	
	public:
		S3ByS2() : ManifoldWithCoordinates(CreateFrame(),2), g(this,OrthFrame())
		{
			Declare_d(e(1),e(2)*e(3));
			Declare_d(e(2),e(3)*e(1));
			Declare_d(e(3),e(1)*e(2));

			test();
		}
		void test()
		{
			Check_ddZero();
			TS_ASSERT_EQUALS(Dimension(),5);
			TS_ASSERT_EQUALS(d(x(1)),e(4));
			TS_ASSERT_EQUALS(d(x(2)),e(5));
			TS_ASSERT_EQUALS(d(sin(x(1))),cos(x(1))*e(4));
			TS_ASSERT_EQUALS(d(1/(x(1))),-1/pow(x(1),2)*e(4));						
		}
	};
public:
//test liederivative.h
	void testLieBracket() {
		LieDerivative a;
		ManifoldWithCoordinates M(3);
		ex dx1=M.d(M.x(1)), dx2=M.d(M.x(2)), dx3=M.d(M.x(3));
		ex f=Function(N.f);
		ex Xf=M.LieDerivative(dx1,f);
		ex Yf=M.LieDerivative(dx2,f);
		TS_ASSERT_DIFFERS(Xf,Yf);
		TS_ASSERT_DIFFERS(Xf,f);
		TS_ASSERT_EQUALS(M.LieDerivative(dx2,Xf),M.LieDerivative(dx1,Yf));
		TS_ASSERT_EQUALS(M.LieDerivative(dx1,M.x(1)),1);
		TS_ASSERT_EQUALS(M.LieDerivative(dx1,M.x(1)*M.x(2)),M.x(2));
		TS_ASSERT_EQUALS(M.LieDerivative(dx1,log(M.x(1))*M.x(2)),M.x(2)/M.x(1));

		S3ByS2 X;
		TS_ASSERT_EQUALS(X.LieBracket(X.e(1),X.e(2)),-X.e(3));
		TS_ASSERT_EQUALS(X.LieBracket(ex_to<VectorField>(X.e(1)),ex_to<VectorField>(X.e(2))),-X.e(3));
		TS_ASSERT_EQUALS(X.LieBracket(X.e(1),X.e(2)),-X.e(3));
		VectorField e1 {ex_to<VectorField>(X.e(1))}, e2 {ex_to<VectorField>(X.e(2))};
		TS_ASSERT_EQUALS(X.LieBracket(e1,e2),-X.e(3));
		TS_ASSERT_EQUALS(X.LieBracket(X.x(1)*X.e(1),X.e(2)),-X.x(1)*X.e(3));

		Xf=X.LieDerivative(X.e(1),f);
		Yf=X.LieDerivative(X.e(4),f);
		TS_ASSERT_EQUALS(Xf+Yf,X.LieDerivative(X.e(1)+X.e(4),f));
		TS_ASSERT_EQUALS(sqrt(3)*Xf,X.LieDerivative(sqrt(3)*X.e(1),f));
		TS_ASSERT_EQUALS(X.x(1)*Xf,X.LieDerivative(X.x(1)*X.e(1),f));
		TS_ASSERT_DIFFERS(Xf,Yf);
		TS_ASSERT_DIFFERS(Xf,f);
		TS_ASSERT_EQUALS(X.LieDerivative(X.e(4),Xf),X.LieDerivative(X.e(1),Yf));
		Xf=X.LieDerivative(X.e(1),f);
		Yf=X.LieDerivative(X.e(2),f);
		TS_ASSERT_DIFFERS(Xf,Yf);
		TS_ASSERT_DIFFERS(Xf,f);
		TS_ASSERT_EQUALS(X.LieDerivative(X.e(2),Xf)-X.LieDerivative(X.e(1),Yf),X.LieDerivative(X.e(3),f));


		f=pow(X.x(1),3);
		Xf=X.LieDerivative(X.e(1),f);
		Yf=X.LieDerivative(X.e(4),f);
		TS_ASSERT_DIFFERS(Xf,Yf);
		TS_ASSERT_DIFFERS(Xf,f);
		TS_ASSERT_EQUALS(Xf,0);
		TS_ASSERT_EQUALS(Yf,diff(f,ex_to<symbol>(X.x(1))));
		TS_ASSERT_EQUALS(X.LieDerivative(X.e(4),Xf),X.LieDerivative(X.e(1),Yf));
		TS_ASSERT_EQUALS(X.LieDerivative(X.e(1),X.e(2)*X.e(3)),X.LieDerivative(X.e(1),X.e(2))*X.e(3)+X.e(2)*X.LieDerivative(X.e(1),X.e(3)));		
		TS_ASSERT_EQUALS(X.LieDerivative(X.e(1),X.e(2)*X.e(4)),Hook(X.e(1),X.d(X.e(2)*X.e(4)))+X.d(Hook(X.e(1),X.e(2)*X.e(4))));
	}

	void testCoordinates() {
		ManifoldWithCoordinates M(3);
		ex dx1=M.d(M.x(1)), dx2=M.d(M.x(2)), dx3=M.d(M.x(3));
		TS_ASSERT_EQUALS(dx1,M.e(1));
		TS_ASSERT_EQUALS(dx2,M.e(2));
		TS_ASSERT_EQUALS(dx3,M.e(3));
		TS_ASSERT_EQUALS(M.Dimension(),3);
				
		TS_ASSERT_EQUALS(M.d(dx1),0);
		TS_ASSERT_EQUALS(M.d(M.x(1)*dx2),dx1*dx2);
		TS_ASSERT_EQUALS(M.d(M.x(1)*dx1),0);
		TS_ASSERT_EQUALS(M.d(M.x(1)*M.x(2)),M.x(1)*dx2+M.x(2)*dx1);
		TS_ASSERT_EQUALS(M.d(M.x(1)*M.x(1)),2*M.x(1)*dx1);
		TS_ASSERT_EQUALS(M.d(sin(M.x(1))),cos(M.x(1))*dx1);
		TS_ASSERT_EQUALS(M.d(1/(M.x(1))),-1/pow(M.x(1),2)*dx1);
		
		S3ByS2 N;
	}
	//(this also tests connection.h)
	void testS2() {
		TwoSphere S2;
		TS_ASSERT_EQUALS(S2.d(S2.x(1)*S2.e(2)),S2.e(1)*S2.e(2));		

		TS_ASSERT_EQUALS(S2.LieBracket(S2.g.e().dual()(1),S2.g.e().dual()(2)).expand(),(S2.x(1)*S2.g.e().dual()(2)-S2.x(2)*S2.g.e().dual()(1)).expand());
		TS_ASSERT_EQUALS(S2.g.ScalarProduct<VectorField>(S2.x(1)*S2.g.e(2)-S2.x(2)*S2.g.e(1),S2.g.e(1)),-16*S2.x(2)*power(1+S2.x(1)*S2.x(1)+S2.x(2)*S2.x(2),-4));

		LOG_INFO(S2.omega);
		matrix R=S2.omega.CurvatureForm();		
		LOG_INFO(R);
		TS_ASSERT_EQUALS(R(0,0).normal(),0);
		TS_ASSERT_EQUALS(R(1,1).normal(),0);
		TS_ASSERT_EQUALS(R(0,1).normal(),S2.g.e(1)*S2.g.e(2));		
		TS_ASSERT_EQUALS(R(1,0).normal(),-S2.g.e(1)*S2.g.e(2));
		
		ExVector e=S2.OrthFrame();
		Connection k(&S2,e);
		for (int i=1;i<=2;++i)
		for (int j=1;j<=2;++j)
			k.DeclareNabla<DifferentialForm>(S2.e(i),S2.e(j),S2.omega.Nabla<DifferentialForm>(S2.e(i),S2.e(j)));
		TS_ASSERT_EQUALS(ex(R).normal(),ex(k.CurvatureForm()).normal());
	}
};

//test liegroup.h, subgroup.h, gl.h, so.h, su.h
class LieGroupTestSuite : public CxxTest::TestSuite  {
public:	 
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

	void testGeneric() {
#ifdef HAVE_COCOA
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
#endif
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
#ifdef HAVE_COCOA		
		LieGroupExtension H2(G,4);
		H2.Declare_d(H2.e(4),0);		
		PolyBasisImplementation<StructureConstant, DefaultPolyAlgorithms> pol;
		H2.GetEquations_ddZero(pol);
		LieGroupExtension H3(H2,5);
		PolyBasisImplementation<StructureConstant, DefaultPolyAlgorithms> pol2;
		H3.GetEquations_ddZero(pol2);
		//for (PolyBasis_impl<StructureConstant, DefaultPolyAlgorithms>::const_iterator i=pol.begin();i!=pol.end();i++)
//			TS_ASSERT(pol2.IdealContains(*i));
#endif

	}

};


#endif /*MANIFOLDS_*/

