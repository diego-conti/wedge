/*
 * structures.h
 *
 *  Created on: Jul 21, 2022
 *      Author: diego
 */

#ifndef SRC_STRUCTURES_H_
#define SRC_STRUCTURES_H_

#include "wedge/manifolds/manifoldwith.h"
#include "wedge/liealgebras/liegroup.h"
#include "wedge/structures/structures.h"

using namespace Wedge;

class StructuresTestSuite : public CxxTest::TestSuite  {
public:
	class SU3 : public AbstractLieGroup<>  {
	public:
		SU3() : AbstractLieGroup<>("-23-45+2*67,13+46-57-[sqrt(3)]*58,-12-47+[sqrt(3)]*48-56,"
						"15-26+37-[sqrt(3)]*38,"
						"-14+27+36+[sqrt(3)]*28,-2*17+24-35,2*16-25-34,-[sqrt(3)]*25+[sqrt(3)]*34")
		{}
	};
//structures.h: SU3Structure, PSU3Structure
	void testStructures()
	{
#ifndef LONG_TEST
		return;
#endif
		{
			ConcreteManifold M(6);
			Frame e=M.e();
			SU3Structure P(&M,e);
			TS_ASSERT_EQUALS((P.psiplus()*P.omega()).expand(),0);
			TS_ASSERT_EQUALS((P.psiminus()*P.omega()).expand(),0);
			TS_ASSERT_EQUALS((P.psiplus()*P.psiminus()/4).expand(),(P.omega()*P.omega()*P.omega()/6).expand());
			TS_ASSERT_EQUALS((P.psiplus()*P.psiminus()/4).expand(),ParseDifferentialForm(e,"123456"));
			TS_ASSERT_EQUALS(P.SquareNorm<DifferentialForm>(P.psiplus()),4);
			TS_ASSERT_EQUALS(P.SquareNorm<DifferentialForm>(P.psiminus()),4);
			TS_ASSERT_EQUALS(P.SquareNorm<DifferentialForm>(P.omega()),3);

			e=ParseDifferentialForms(e,"1+2,1+3,1+4,1+5,1+6,6");
			P=SU3Structure (&M,e);
			TS_ASSERT_EQUALS((P.psiplus()*P.omega()).expand(),0);
			TS_ASSERT_EQUALS((P.psiminus()*P.omega()).expand(),0);
			TS_ASSERT_EQUALS((P.psiplus()*P.psiminus()/4).expand(),(P.omega()*P.omega()*P.omega()/6).expand());
			TS_ASSERT_EQUALS((P.psiplus()*P.psiminus()/4).expand(),ParseDifferentialForm(e,"123456").expand());
			TS_ASSERT_EQUALS(P.SquareNorm<DifferentialForm>(P.psiplus()),4);
			TS_ASSERT_EQUALS(P.SquareNorm<DifferentialForm>(P.psiminus()),4);
			TS_ASSERT_EQUALS(P.SquareNorm<DifferentialForm>(P.omega()),3);
		}
		{
			ConcreteManifold M(8);
			Frame e=M.e();
			PSU3Structure P(&M,e);
			VectorSpace<DifferentialForm> V=P.LambdaComponent(PSU3Structure::Lambda_2_8);
			TS_ASSERT_EQUALS(VectorSpace<DifferentialForm>(V.e_begin(),V.e_end()).Dimension(),8);
			V=P.LambdaComponent(PSU3Structure::Lambda_3_8);
			TS_ASSERT_EQUALS(VectorSpace<DifferentialForm>(V.e_begin(),V.e_end()).Dimension(),8);
			V=P.LambdaComponent(PSU3Structure::Lambda_4_8Plus);
			TS_ASSERT_EQUALS(VectorSpace<DifferentialForm>(V.e_begin(),V.e_end()).Dimension(),8);
			V=P.LambdaComponent(PSU3Structure::Lambda_4_8Minus);
			TS_ASSERT_EQUALS(VectorSpace<DifferentialForm>(V.e_begin(),V.e_end()).Dimension(),8);
			V=P.LambdaComponent(PSU3Structure::Lambda_6_8);
			TS_ASSERT_EQUALS(VectorSpace<DifferentialForm>(V.e_begin(),V.e_end()).Dimension(),8);

			V=P.LambdaComponent(PSU3Structure::Lambda_2_20);
			TS_ASSERT_EQUALS(VectorSpace<DifferentialForm>(V.e_begin(),V.e_end()).Dimension(),20);
			V=P.LambdaComponent(PSU3Structure::Lambda_3_20);
			TS_ASSERT_EQUALS(VectorSpace<DifferentialForm>(V.e_begin(),V.e_end()).Dimension(),20);
	//verify that the map \[[T\otimes\Lambda^3_{20}\to\Lambda^2\] induced by interior product is surjective.
			exvector Two;
			for (int j=1;j<=V.Dimension();j++)
				for (int i=1;i<=8;i++)
				{
					Two.push_back(Hook(M.e(i),V.e(j)));
				}
			TS_ASSERT_EQUALS(VectorSpace<DifferentialForm>(Two).Dimension(),28);

			V=P.LambdaComponent(PSU3Structure::Lambda_6_20);
			TS_ASSERT_EQUALS(VectorSpace<DifferentialForm>(V.e_begin(),V.e_end()).Dimension(),20);

			V=P.LambdaComponent(PSU3Structure::Lambda_4_27Plus);
			TS_ASSERT_EQUALS(VectorSpace<DifferentialForm>(V.e_begin(),V.e_end()).Dimension(),27);
			V=P.LambdaComponent(PSU3Structure::Lambda_4_27Minus);
			TS_ASSERT_EQUALS(VectorSpace<DifferentialForm>(V.e_begin(),V.e_end()).Dimension(),27);


			e=ParseDifferentialForms(e,"1+2,1+3,1+4,1+5,1+6,6,7,8");
			P=PSU3Structure (&M,e);

			V=P.LambdaComponent(PSU3Structure::Lambda_2_8);
			TS_ASSERT_EQUALS(VectorSpace<DifferentialForm>(V.e_begin(),V.e_end()).Dimension(),8);
			V=P.LambdaComponent(PSU3Structure::Lambda_3_8);
			TS_ASSERT_EQUALS(VectorSpace<DifferentialForm>(V.e_begin(),V.e_end()).Dimension(),8);
			V=P.LambdaComponent(PSU3Structure::Lambda_4_8Plus);
			TS_ASSERT_EQUALS(VectorSpace<DifferentialForm>(V.e_begin(),V.e_end()).Dimension(),8);
			V=P.LambdaComponent(PSU3Structure::Lambda_4_8Minus);
			TS_ASSERT_EQUALS(VectorSpace<DifferentialForm>(V.e_begin(),V.e_end()).Dimension(),8);
			V=P.LambdaComponent(PSU3Structure::Lambda_6_8);
			TS_ASSERT_EQUALS(VectorSpace<DifferentialForm>(V.e_begin(),V.e_end()).Dimension(),8);

			V=P.LambdaComponent(PSU3Structure::Lambda_2_20);
			TS_ASSERT_EQUALS(VectorSpace<DifferentialForm>(V.e_begin(),V.e_end()).Dimension(),20);
			V=P.LambdaComponent(PSU3Structure::Lambda_3_20);
			TS_ASSERT_EQUALS(VectorSpace<DifferentialForm>(V.e_begin(),V.e_end()).Dimension(),20);
				//verify that the map \[[T\otimes\Lambda^3_{20}\to\Lambda^2\] induced by interior product is surjective.
			Two.clear();
			for (int j=1;j<=V.Dimension();j++)
				for (int i=1;i<=8;i++)
				{
					Two.push_back(Hook(M.e(i),V.e(j)));
				}
			TS_ASSERT_EQUALS(VectorSpace<DifferentialForm>(Two).Dimension(),28);

			V=P.LambdaComponent(PSU3Structure::Lambda_6_20);
			TS_ASSERT_EQUALS(VectorSpace<DifferentialForm>(V.e_begin(),V.e_end()).Dimension(),20);

			V=P.LambdaComponent(PSU3Structure::Lambda_4_27Plus);
			TS_ASSERT_EQUALS(VectorSpace<DifferentialForm>(V.e_begin(),V.e_end()).Dimension(),27);
			V=P.LambdaComponent(PSU3Structure::Lambda_4_27Minus);
			TS_ASSERT_EQUALS(VectorSpace<DifferentialForm>(V.e_begin(),V.e_end()).Dimension(),27);

			matrix g=PSU3Structure::Metric(M,P.phi());
			Frame x=P.e().dual();
			TS_ASSERT_EQUALS(g.rows(),8);
			TS_ASSERT_EQUALS(g.cols(),8);
			ex conformalfactor;
			for (int i=0;i<8;i++)
				for (int j=0;j<8;j++)
				{
					exvector v=x.Components(M.e()[i]),w=x.Components(M.e()[j]);
					ex gij;
					for (int k=0;k<8;k++) gij+=v[k]*w[k];
					if (g(i,j).is_zero()) {
						TS_ASSERT(gij.is_zero());
					}
					else {
						if (conformalfactor.is_zero()) conformalfactor=gij/g(i,j);
						else TS_ASSERT_EQUALS(conformalfactor,gij/g(i,j));
					}
				}

			for (int i=0;i<8;i++)
			{
				ex ei=M.e()[i];
			//verify that the sigma's are right inverses of the pi's
				TS_ASSERT_EQUALS(P.PiZero(P.SigmaZero(ei)),ei);
				TS_ASSERT_EQUALS(P.PiPlus(P.SigmaPlus(ei)),ei);
				TS_ASSERT_EQUALS(P.PiMinus(P.SigmaMinus(ei)),ei);
				//verify that the sigma's take values in the spaces they are supposed to
				LOG_INFO(P.LambdaComponent(PSU3Structure::Lambda_4_8Minus));
				LOG_INFO(P.SigmaMinus(ei));
				TS_ASSERT(P.LambdaComponent(PSU3Structure::Lambda_4_8Minus).Contains(P.SigmaMinus(ei)));
				TS_ASSERT(P.LambdaComponent(PSU3Structure::Lambda_4_8Plus).Contains(P.SigmaPlus(ei)));
				TS_ASSERT(P.LambdaComponent(PSU3Structure::Lambda_6_8).Contains(P.SigmaZero(ei)));
			}
		}
	}

	//helper function for testliegroupstructures
	void helperTestPSU3(const Manifold& G, const GStructureHasParameters<PSU3Structure,false>& P)
	{
		IntrinsicTorsion it=P.GetIntrinsicTorsion();
		TS_ASSERT(P.LambdaComponent(PSU3Structure::Lambda_4_27Plus).Contains(it[PSU3Structure::Xi27Plus]));
		TS_ASSERT(P.LambdaComponent(PSU3Structure::Lambda_4_27Minus).Contains(it[PSU3Structure::Xi27Minus]));
		TS_ASSERT(P.LambdaComponent(PSU3Structure::Lambda_6_20).Contains(it[PSU3Structure::Xi20]));
		TS_ASSERT_EQUALS(G.d(P.phi()),it[PSU3Structure::Xi27Plus]+it[PSU3Structure::Xi27Minus]+P.SigmaPlus(it[PSU3Structure::Xi8Plus])+P.SigmaMinus(it[PSU3Structure::Xi8Minus]));
		TS_ASSERT(P.LambdaComponent(PSU3Structure::Lambda_6_8).Contains(G.d(P.starphi())-it[PSU3Structure::Xi20]));
	}
//liegroupstructures.h; intrinsic torsion of structures with no parameters
	void testliegroupstructures()
	{
		{
			LieGroupWith<SU3Structure> G("0,0,0,0,12+34,13+42","1,2,3,4,5,6");
			IntrinsicTorsion it=G.P.GetIntrinsicTorsion();
			TS_ASSERT_EQUALS(it[SU3Structure::W1Plus],0);
			TS_ASSERT_EQUALS(it[SU3Structure::W1Minus],ex(1)/3);
			TS_ASSERT_EQUALS(it[SU3Structure::W2Plus],0);
			TS_ASSERT_EQUALS(it[SU3Structure::W2Minus],ParseDifferentialForm(G.e(),"2/3*12+2/3*34-4/3*56"));
			TS_ASSERT_EQUALS(it[SU3Structure::W3],ParseDifferentialForm(G.e(),"-1/2*135-1/2*146-1/2*236+1/2*245"));
			TS_ASSERT_EQUALS(it[SU3Structure::W4],G.e(6));
			TS_ASSERT_EQUALS(it[SU3Structure::W5],0);
			TS_ASSERT_EQUALS(it.Type(),SU3Structure::W1Minus+SU3Structure::W2Minus + SU3Structure::W3+SU3Structure::W4);
			TS_ASSERT_EQUALS(it.GenericType()-it.Type(),SU3Structure::W1Plus+SU3Structure::W2Plus + SU3Structure::W5);
		}
		{
			LieGroupWith<PSU3Structure> G("-23-45+2*67,13+46-57-[sqrt(3)]*58,-12-47+[sqrt(3)]*48-56,"
						"15-26+37-[sqrt(3)]*38,"
						"-14+27+36+[sqrt(3)]*28,-2*17+24-35,2*16-25-34,-[sqrt(3)]*25+[sqrt(3)]*34",
						"1,6,7,3,2,4,5,8");
#ifdef LONG_TEST
			helperTestPSU3(G,G.P);
			helperTestPSU3(G,GStructureHasParameters<PSU3Structure,false>(&G,ParseDifferentialForms(G.e(),"1,3,7,5,4,2,6,8")));
#endif
		}

	}

	class StructureWithParameters : public SU3, public GStructureHasParameters<PSU3Structure,true> {
	public:
		StructureWithParameters() :  GStructureHasParameters<PSU3Structure,true>(this)
		{
			ExVector e=ParseDifferentialForms(SU3::e(),"1,6,7,5,4,2,3,8");
			for (int i=2;i<=8;i++)
				e(1)+=GStructureParameter("x"+ToString(i))*e(i);
			SetFrame(e);
		}
	};
//structures.h; test PSU3-structures with parameters
	void testStructureWParams()
	{
		StructureWithParameters P;
 		P.DeclareZero(P.d(P.starphi()).expand());
 		TS_ASSERT_EQUALS(P.d(P.starphi()),0);
 		TS_ASSERT_EQUALS(P.d(P.phi()),0);

		StructureWithParameters P2;
 		P2.DeclareZero(P2.d(P2.phi()).expand());
 		TS_ASSERT_EQUALS(P2.d(P2.starphi()),0);
 		TS_ASSERT_EQUALS(P2.d(P2.phi()),0);
 #ifdef HAVE_COCOA
 		{
			StructureWithParameters P;
 			PolyBasisImplementation<GStructureParameter, CocoaPolyAlgorithms_R> I;
 			P.GetEquationsForTorsionIn(I,PSU3Structure::Xi14+ PSU3Structure::Xi8Plus+PSU3Structure::Xi8Minus);
			I.Reduce();
			LOG_INFO(I);
 		}
#endif

	}


	void testcheckorthogonalityrelations()
	{
		ManifoldWith<PSU3Structure> M;
		VectorSpace<DifferentialForm> V1=M.LambdaComponent(PSU3Structure::Lambda_4_8Plus);
		VectorSpace<DifferentialForm> V2=M.LambdaComponent(PSU3Structure::Lambda_4_8Minus);
		VectorSpace<DifferentialForm> V3=M.LambdaComponent(PSU3Structure::Lambda_4_27Plus);
		VectorSpace<DifferentialForm> V4=M.LambdaComponent(PSU3Structure::Lambda_4_27Minus);

		for (IBasis<DifferentialForm>::const_iterator i=V1.e_begin();i!=V1.e_end();++i)
		{
			for (IBasis<DifferentialForm>::const_iterator j=V2.e_begin();j!=V2.e_end();++j)
				TS_ASSERT((*i * *j).expand().is_zero());
			for (IBasis<DifferentialForm>::const_iterator j=V3.e_begin();j!=V3.e_end();++j)
				TS_ASSERT((*i * *j).expand().is_zero());
			for (IBasis<DifferentialForm>::const_iterator j=V4.e_begin();j!=V4.e_end();++j)
				TS_ASSERT((*i * *j).expand().is_zero());
		}
		for (IBasis<DifferentialForm>::const_iterator i=V2.e_begin();i!=V2.e_end();++i)
		{
			for (IBasis<DifferentialForm>::const_iterator j=V3.e_begin();j!=V3.e_end();++j)
				TS_ASSERT((*i * *j).expand().is_zero());
			for (IBasis<DifferentialForm>::const_iterator j=V4.e_begin();j!=V4.e_end();++j)
				TS_ASSERT((*i * *j).expand().is_zero());
		}
		for (IBasis<DifferentialForm>::const_iterator i=V3.e_begin();i!=V3.e_end();++i)
			for (IBasis<DifferentialForm>::const_iterator j=V4.e_begin();j!=V4.e_end();++j)
				TS_ASSERT((*i * *j).expand().is_zero());
	}


};






#endif /* SRC_STRUCTURES_H_ */
