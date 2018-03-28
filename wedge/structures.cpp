/***************************************************************************
 *   Copyright (C) 2007, 2008, 2009 by Diego Conti, diego.conti@unimib.it  *
 *                                                                         *
 *   This file is part of Wedge.                                           *
 *                                                                         *
 *   Wedge is free software; you can redistribute it and/or modify         *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 3 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   Wedge is distributed in the hope that it will be useful,              *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with Wedge; if not, write to the                                *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#include "structures.h"
#include "parse.h"
#include "tensor.h" 
namespace Wedge {

//Defining forms
ex SU2Structure::psi() const
{
	return u(0);
}

ex SU2Structure::alpha() const
{
	return e(5);
}

ex SU2Structure::omega1() const
{
	return e(1)*e(2)+e(3)*e(4);
}

ex SU2Structure::omega2() const
{
	return e(1)*e(3)+e(4)*e(2);
}

ex SU2Structure::omega3() const
{
	return e(1)*e(4)+e(2)*e(3);
}


ex SU3Structure::omega() const
{
	return e(1)*e(2)+e(3)*e(4)+e(5)*e(6);
}

ex SU3Structure::psiplus() const
{
	return e(1)*e(3)*e(5)-e(1)*e(4)*e(6)-e(2)*e(3)*e(6)-e(2)*e(4)*e(5);		
}

ex SU3Structure::psiminus() const
{
	return e(1)*e(3)*e(6)+e(1)*e(4)*e(5)+e(2)*e(3)*e(5)-e(2)*e(4)*e(6);
}

ex SU3StructureDim7::psi() const
{
	return u(0);
}

ex SU3StructureDim7::alpha() const
{
	return e(7);
}

ex SU3StructureDim7::F() const
{
	return e(1)*e(2)+e(3)*e(4)+e(5)*e(6);
}

ex SU3StructureDim7::OmegaPlus() const
{
	return e(1)*e(3)*e(5)-e(1)*e(4)*e(6)-e(2)*e(3)*e(6)-e(2)*e(4)*e(5);
}

ex SU3StructureDim7::OmegaMinus() const
{
	return e(1)*e(3)*e(6)+e(1)*e(4)*e(5)+e(2)*e(3)*e(5)-e(2)*e(4)*e(6);
}


ex G2Structure::phi() const
{
	return e(1)*e(3)*e(5)-e(1)*e(4)*e(6)-e(2)*e(3)*e(6)-e(2)*e(4)*e(5)+
		 e(1)*e(2)*e(7)+e(3)*e(4)*e(7)+e(5)*e(6)*e(7);
	;
}

ex G2Structure::starphi() const
{
	ex omega=e(1)*e(2)+e(3)*e(4)+e(5)*e(6);
	ex psiminus=e(1)*e(3)*e(6)+e(1)*e(4)*e(5)+e(2)*e(3)*e(5)-e(2)*e(4)*e(6);
	return omega*omega/2+psiminus*e(7);
}


ex PSU3Structure::phi() const
{
	ExVector Omega=ParseDifferentialForms(e(),"45-67,47-56,46-75,46+75");
	return e(1)*e(2)*e(3)+e(1)*Omega(1)/2+e(2)*Omega(2)/2+e(3)*Omega(3)/2+sqrt(ex(3))/2*e(8)*Omega(4);
}

ex PSU3Structure::starphi() const
{
	ExVector Omega=ParseDifferentialForms(e(),"45-67,47-56,46-75,46+75");
	return (
		Omega(4)*Omega(4)*e(8)-e(2)*e(3)*e(8)*Omega(1)+e(1)*e(3)*e(8)*Omega(2)-e(1)*e(2)*e(8)*Omega(3)-
		sqrt(ex(3))*e(1)*e(2)*e(3)*Omega(4)
	).expand()/2;
}


//Intrinsic torsion
template<> IntrinsicTorsion  GStructureHasParameters<SU3Structure,false>::GetIntrinsicTorsion() const
{
	if (it.empty())
	{
		ex dpsiplus=M()->d(psiplus());
		ex dpsiminus=M()->d(psiminus());
		ex domega=M()->d(omega());	
		it[W5]=Hook(psiplus(),dpsiplus);
		it[W1Plus]=HodgeStar(omega()*dpsiplus)/6;
		it[W2Plus]=Hook(omega(),dpsiplus-it[W1Plus]*omega()*omega());
		it[W1Minus]=HodgeStar(omega()*dpsiminus)/6;
		it[W2Minus]=Hook(omega(),dpsiminus-it[W1Minus]*omega()*omega());	
		it[W4]=Hook(omega(),domega)/2;
		it[W3]=domega-it[W4]*omega()+ex(3)/2*it[W1Minus]*psiplus()-ex(3)/2*it[W1Plus]*psiminus();
		assert(dpsiplus.is_equal((psiplus()*it[W5]+it[W2Plus]*omega()+it[W1Plus]*omega()*omega()).expand()));
		assert(dpsiminus.is_equal((psiminus()*it[W5]+it[W2Minus]*omega()+it[W1Minus]*omega()*omega()).expand()));	
	}
	return it;
}

template<> IntrinsicTorsion GStructureHasParameters<PSU3Structure,false>::GetIntrinsicTorsion() const
{
	if (it.empty())
	{
		ex dphi=M()->d(phi());
		ex stardphi=HodgeStar(dphi);
		ex dstarphi=M()->d(starphi());
		it[Xi8Plus]=PiPlus(dphi);
		it[Xi8Minus]=PiMinus(dphi);
		it[Xi20]=dstarphi-SigmaZero(PiZero(dstarphi));
		it[Xi27Plus]=(dphi+stardphi)/2-SigmaPlus(it[Xi8Plus]);
		it[Xi27Minus]=(dphi-stardphi)/2-SigmaMinus(it[Xi8Minus]);
		
//FIXME: What happens if M is a submersion?
		LeviCivitaConnection<true> leviCivita(M(),*this);
		it[Xi14]=0;		
		for (int i=1;i<=8;i++)
		{
			ex diff=leviCivita.Nabla<DifferentialForm>(e().dual()(i),phi());
			diff-=e(i)*HodgeStar(dstarphi)/6+Hook(e(i),dphi)/4;
			it[Xi14]+=TensorProduct<DifferentialForm,DifferentialForm>(e(i),diff.expand());
		}
	}
	return it;
}


namespace internal {
/** @internal @brief Helper class for PSU3Structure::LambdaComponent 
 */
struct ReferencePSU3Structure {
	ConcreteManifold M;
	PSU3Structure P;
	vector<VectorSpace<DifferentialForm> > spaces;
	ReferencePSU3Structure() : M(8),P(&M,M.e()){}

	void init() {
		if (!spaces.empty()) return;
		spaces.resize(PSU3Structure::nLambdaComponents);
		
		ex phi=P.phi(),starphi=P.starphi();
		//2-forms	
		exvector L_2_8(8);
		for (int i=0;i<8;i++)
			L_2_8[i]=Hook(M.e()[i],phi);
		spaces[PSU3Structure::Lambda_2_8].SetBasis(L_2_8);
		
		VectorSpace<DifferentialForm> TwoForms=M.pForms(2);		
		exvector eqns(8);
		for (int i=0;i<8;i++)
			eqns[i]=TrivialPairing<DifferentialForm>(M.e()[i]*TwoForms.GenericElement(),phi);
		spaces[PSU3Structure::Lambda_2_20]=TwoForms.SubspaceFromEquations(eqns.begin(),eqns.end());
		assert(spaces[PSU3Structure::Lambda_2_20].Dimension()==20);
		
		//3-forms
		exvector L_3_20(20);		
		for (int i=0;i<20;i++)		 
			L_3_20[i]=M.HodgeStar(spaces[PSU3Structure::Lambda_2_20].e()[i]*phi);
		spaces[PSU3Structure::Lambda_3_20].SetBasis(L_3_20);
		
		exvector L_3_8(8);		
		for (int i=0;i<8;i++)
			L_3_8[i]=M.HodgeStar(L_2_8[i]*phi);
		spaces[PSU3Structure::Lambda_3_8].SetBasis(L_3_8);
		//4-forms
		exvector L_4_8Plus(8),L_4_8Minus(8);	
		for (int i=0;i<8;i++)
		{
			ex a=phi*M.e()[i].expand();
			ex b=Hook(M.e()[i],starphi).expand();
			L_4_8Plus[i]=a+b;
			L_4_8Minus[i]=a-b;
		}
		spaces[PSU3Structure::Lambda_4_8Plus].SetBasis(L_4_8Plus);
		spaces[PSU3Structure::Lambda_4_8Minus].SetBasis(L_4_8Minus);

		exvector vPlus,vMinus;
		for (int i=0;i<8;i++)
			for (int j=0;j<8;j++)
			{
				ex x=(M.e()[i]*L_3_8[j]).expand();
				ex starx=M.HodgeStar(x);
				vPlus.push_back(x+starx);
				vMinus.push_back(x-starx);
			}
		VectorSpace<DifferentialForm> VPlus(vPlus),VMinus(vMinus);
		Subspace<DifferentialForm> ZPlus=VPlus.Subspace(L_4_8Plus.begin(),L_4_8Plus.end(),TrivialPairingOperator<DifferentialForm>());
		Subspace<DifferentialForm> ZMinus=VMinus.Subspace(L_4_8Minus.begin(),L_4_8Minus.end(),TrivialPairingOperator<DifferentialForm>());
		spaces[PSU3Structure::Lambda_4_27Plus].SetBasis(Basis<DifferentialForm>(ZPlus.complement_begin(),ZPlus.complement_end()));
		spaces[PSU3Structure::Lambda_4_27Minus].SetBasis(Basis<DifferentialForm>(ZMinus.complement_begin(),ZMinus.complement_end()));

		//6-forms
		exvector L_6_8(8);		
		for (int i=0;i<8;i++)
			L_6_8[i]=M.HodgeStar(L_2_8[i]);
		spaces[PSU3Structure::Lambda_6_8].SetBasis(L_6_8);

		exvector L_6_20(20);		
		for (int i=0;i<20;i++)
			L_6_20[i]=M.HodgeStar(spaces[PSU3Structure::Lambda_2_20].e()[i]);
		spaces[PSU3Structure::Lambda_6_20].SetBasis(L_6_20);
	}
} referencePSU3Structure;
}

VectorSpace<DifferentialForm> PSU3Structure::LambdaComponent(LambdaComponentType type) const
{
	internal::referencePSU3Structure.init();	 	
	ex subs=LinearMapToSubstitutions<DifferentialForm>(internal::referencePSU3Structure.P.e(),e());
	exvector v; 
	v.reserve(internal::referencePSU3Structure.spaces[type].Dimension());
	for (Basis<DifferentialForm>::const_iterator i=internal::referencePSU3Structure.spaces[type].e_begin();
		i!=internal::referencePSU3Structure.spaces[type].e_end();
			++i)
		v.push_back(i->subs(subs).expand());
	return v;	 
}

bool PSU3Structure::IsGamma30Zero() const {
	ex dstarphi=M()->d(starphi());
	VectorSpace<DifferentialForm> L_2_20=LambdaComponent(Lambda_2_20);
	return (dstarphi*L_2_20.GenericElement()).expand().is_zero();
} 

ex PSU3Structure::PiPlus(ex fourform) const
{
	ex proj=(fourform+HodgeStar(fourform))/2;
	return HodgeStar(proj*phi());	
}

ex PSU3Structure::PiMinus(ex fourform) const
{
	ex proj=(fourform-HodgeStar(fourform))/2;
	return HodgeStar(proj*phi());	
}

ex PSU3Structure::PiZero(ex sixform) const
{	
	return Hook(starphi(),sixform);	
}

ex PSU3Structure::SigmaPlus(ex oneform) const
{
	ex fourform=oneform*phi();
	return -2*(fourform+HodgeStar(fourform))/5;
}

ex PSU3Structure::SigmaMinus(ex oneform) const
{
	ex fourform=oneform*phi();
	return 2*(fourform-HodgeStar(fourform))/5;
}

ex PSU3Structure::SigmaZero(ex oneform) const
{
	return 2*(starphi()*oneform).expand()/3;	
}		

matrix PSU3Structure::Metric(const Manifold& M,ex phi) 
{	
	matrix d(9,9);
	matrix g(8,8);
	for (int i=1;i<=8;i++)
	for (int j=1;j<=8;j++)
		d(i,j)=(Wedge::Hook(M.e(i),phi)*Wedge::Hook(M.e(j),phi)*phi).expand();
	for (int i=1;i<=8;i++)
	for (int j=1;j<=8;j++)
	{
		g(i-1,j-1)=0;
		for (int h=1;h<=8;h++)
		for (int k=1;k<=8;k++)
		{
			g(i-1,j-1)+=M.HodgeStar(d(h,j)*M.e(k))*M.HodgeStar(d(i,k)*M.e(h));
		}
	}
//	ex det=pow(g.determinant(),-ex(1)/9).expand();
//	for (int i=0;i<8;i++)
//	for (int j=0;j<8;j++)
// 		g(i,j)*=det;
	return g;	
}

//torsion classes
IntrinsicTorsionClass PSU3Structure::Xi8Plus(N.Xi(8).Plus());
IntrinsicTorsionClass PSU3Structure::Xi8Minus(N.Xi(8).Minus());
IntrinsicTorsionClass PSU3Structure::Xi27Plus(N.Xi(27).Plus());
IntrinsicTorsionClass PSU3Structure::Xi27Minus(N.Xi(27).Minus());
IntrinsicTorsionClass PSU3Structure::Xi20(N.Xi(20));
IntrinsicTorsionClass PSU3Structure::Xi14(N.Xi(14));

IntrinsicTorsionClass SU3Structure::W1Plus(N.W(1).Plus());
IntrinsicTorsionClass SU3Structure::W1Minus(N.W(1).Minus());
IntrinsicTorsionClass SU3Structure::W2Plus(N.W(2).Plus());
IntrinsicTorsionClass SU3Structure::W2Minus(N.W(2).Minus());
IntrinsicTorsionClass SU3Structure::W3(N.W(3));
IntrinsicTorsionClass SU3Structure::W4(N.W(4));
IntrinsicTorsionClass SU3Structure::W5(N.W(5));

}
