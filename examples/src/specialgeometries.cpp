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
#include <wedge/wedge.h>
using namespace Wedge;
using namespace std;


class MyGStructure : public GStructureHasParameters<SU3Structure, true> {
	Frame CreateFrame(const Frame& e)
	{
		Frame _e=e;
		GStructureParameter a(N.a),b(N.b),c(N.c);
		_e(1)=e(1)+a*e(2)+b*e(3)+c*e(4);
		return _e;
	}
public:
	MyGStructure(const Manifold* manifold) : GStructureHasParameters<SU3Structure, true>(manifold, CreateFrame(manifold->e())) {}
};

//test whether a family of SU(3)-structures on a fixed nilpotent Lie group contains a half-flat structure
void HalfFlatLieGroup()
{
	AbstractLieGroup<> G("0,0,12,13,14+23,34+52");
	cout<<"Lie group with structure constants"<<endl<<G<<endl;
	MyGStructure P(&G);
	cout<<"The SU(3)-structure: "<<P<<endl;
	PolyBasisImplementation<GStructureParameter,DefaultPolyAlgorithms> V;
	P.GetEquationsForTorsionIn(V,SU3Structure::W3+SU3Structure::W2Minus+SU3Structure::W1Minus);
	cout<<"is half-flat for the parameters in the ideal "<<endl<<V<<endl;
	V.Reduce();
	cout<<"or equivalently"<<endl<<V<<endl;	
}


//SU(3)<GL(6,R) as the group that fixes a symplectic form and a complex volume
void SU3AsAStructureGroup()
{
	ManifoldWith<SU3Structure> M;
	GL gl(6);
	GLRepresentation<VectorField> V(&gl,M.e());
	if (OrbitDimension<DifferentialForm>(gl,V,M.omega()),M.pForms(2).Dimension())
		cout<<"omega is stable"<<endl;	
	if (OrbitDimension<DifferentialForm>(gl,V,M.psiplus()),M.pForms(3).Dimension())
		cout<<"psi+ is stable"<<endl;	 
	AbstractLieSubgroup<false> G=Stabilizer<DifferentialForm>(gl,V,M.omega()+M.psiplus());
	cout<<"Joint stabilizer has dimension "<<G.Dimension();
	cout<<" and Betti numbers "<<G.BettiNumbers()<<endl;
}	

class SymplecticSixManifold : public ConcreteManifold {
public:
	SymplecticSixManifold() : ConcreteManifold(6) {
		cout<<"Symplectic 6-manifold"<<endl;
		ex omega=ParseDifferentialForm(e(),"12+34+56");
		Basis<DifferentialForm> b1,b2;
		for (int i=1;i<=6;++i)
		{
			b1.push_back(omega*e(i));
			b2.push_back(Hook(e(i),omega*omega));
		}
		cout<<"\\{\\alpha\\wedge\\omega,\\alpha\\in T^*X\\}="<<b1;
		cout<<"\\{v \\hook(\\omega\\wedge\\omega),v\\in TX\\}="<<b2;
		VectorSpace<DifferentialForm> ThreeForms=pForms(3);
		exvector eqns;
		GetCoefficients<DifferentialForm>(eqns,ThreeForms.GenericElement()*omega);
		cout<<"\\{\\phi\\in\\Lambda^3(X)\\mid\\omega\\wedge\\phi=0\\}="<<
			ThreeForms.SubspaceFromEquations(eqns.begin(),eqns.end()).e();
	}
};

int main()
{
	SymplecticSixManifold Y;
	HalfFlatLieGroup();
	SU3AsAStructureGroup();
	return 0;
}
