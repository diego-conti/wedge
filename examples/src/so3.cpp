/***************************************************************************
 *   Copyright (C) 2008, 2009 by Diego Conti				   *
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
#include <stdlib.h>

using namespace Wedge;
using namespace GiNaC;

//Irreducible representation of SO(3) on R^9
class SO3Structure : public RiemannianStructure {
public:
	enum {dimension=9};
	SO3Structure(const Manifold* manifold, const Frame& frame) : RiemannianStructure(manifold,frame) {}
	ex Gamma() const {
		return ParseDifferentialForm(e(), "7/3*3478+3456+5678-2468-2457+1/3 *2367-2358-1467-4/3 *1458-1368-8/3 *1357+1278+1234-8/3 *1256+2/3*sqrt(5)*4569-2/3*sqrt(5)*1689 -2/3*sqrt(5)*2589+2/3*sqrt(5)*1249+1/3*sqrt(35)*2349-1/3*sqrt(35)*3689-1/3*sqrt(35)*2789+1/3*sqrt(35)*4679");	
	}
	ex StarGamma() const {
		return Hook(Gamma(),ParseDifferentialForm(e(),"123456789"));
	}
};

//Irreducible representation of SO(3) on R^9 w.r.t a non-orthonormal-frame, so that no roots appear in the expression of the invariant forms.
class SO3StructureNoRoots : public RiemannianStructure {
public:
	enum {dimension=9};
	SO3StructureNoRoots(const Manifold* manifold, const Frame& frame) : RiemannianStructure(manifold,frame) {}
	ex Gamma() const {
		return Hook(StarGamma(),ParseDifferentialForm(e(),"123456789"));
	}
	ex StarGamma() const {
		return ParseDifferentialForm(e(),"2520 *35678+63 *34569+63 *56789-34789+2520*12378+24696*12569+35280*12358+63*12349-35280*12457+35280*15678-126*23679-63*23589-2520*13467-14*24689-35280*13456-63*24579-63*13689+63*12789-2520*23457-756*13579-63*14679+14*14589");
	}
};


void DisplaySO3Representation(int n) 
{
	cout<<latex<<"Irreducible action of SO(3) on R^"<<n+1<<" (with [H,X]=\\sqrt2 Y)"<<endl;
	ConcreteManifold M(n+1);
	SO3Representation<VectorField> V(M.e());
	matrix mH(n+1,n+1), mY(n+1,n+1), mX(n+1,n+1);

	symbol X("b"), Y("c"), H("a");
	for (int i=1;i<=n+1;++i)
	for (int j=1;j<=n+1;++j)
	{
		mH(i-1,j-1)=NormalizeRoots(TrivialPairing<DifferentialOneForm>(V.H<VectorField>(M.e(j)),M.e(i))/sqrt(ex(2)));
		mX(i-1,j-1)=NormalizeRoots(TrivialPairing<DifferentialOneForm>(V.X<VectorField>(M.e(j)),M.e(i))/sqrt(ex(2)));
		mY(i-1,j-1)=NormalizeRoots(TrivialPairing<DifferentialOneForm>(V.Y<VectorField>(M.e(j)),M.e(i))/sqrt(ex(2)));
	}
	cout<<"H acts by "<<mH<<endl;
	cout<<"X acts by "<<mX<<endl;
	cout<<"Y acts by "<<mY<<endl;
	cout<<"Generic element aH + bX + cY acts as: "<<endl<<(X*mX+H*mH+Y*mY).expand().evalm()<<endl;

	cout<<"The SO(3) action preserves the metric "<<InvariantMetrics(V,M.e()).e()<<endl;
	cout<<"and the forms "<<endl;
	for (int k=1;k<=n+1;++k)
	{	
		VectorSpace<DifferentialForm> inv=Invariant_pForms(V,M.e(),k);
		if (inv.Dimension()!=0) cout<<inv.e()<<endl;
	}
}

class Grassmannian : public ConcreteManifold {
	SO3Structure g;
public:
	Grassmannian () : ConcreteManifold(9) , g(this,e()) {				
		cout<<latex<<"*\\gamma="<<g.StarGamma()<<endl;
		Hessian(e(1)*e(3)*e(5)*e(7));
		Hessian(e(1)*e(2)*e(5)*e(6));
	}

	void Hessian(ex fourform)
	{
		cout<<"Computing Hessian of the map f(\\alpha)= *(\\alpha\\wedge *\\gamma)"<<endl;
		cout<<"at \\alpha="<<fourform<<" \\in Gr(4,9)"<<endl;
		SO so(9);
		SORepresentation<VectorField> V(&so,e());
		VectorSpace<DifferentialForm> oneforms=so.pForms(1);
		ex act=V.Action<DifferentialForm>(oneforms.GenericElement(),fourform);
		list<ex> eqns;
		GetCoefficients<DifferentialForm>(eqns,act);

		Subspace<DifferentialForm> m=oneforms.SubspaceFromEquations(eqns.begin(),eqns.end()).Complement();
		cout<<"If h\\subset so(n) is the stabilizer of \\alpha, its complement is spanned by = "<<m.e();
		matrix A(m.Dimension(),m.Dimension());
		for (int i=1;i<=m.Dimension();++i)
		for (int j=1;j<=m.Dimension();++j)
		{
			ex act=V.Action<DifferentialForm>(m.e(i),fourform);
			act=V.Action<DifferentialForm>(m.e(j),act);
			A(i-1,j-1)=NormalizeRoots(HodgeStar(act*g.StarGamma()));
		}
		cout<<"In this basis, the Hessian is given by H="<<A<<endl;

		ex v=-m.e(4)+m.e(8)*(-1/ex(6)+1/ex(12)*(17-3*sqrt(ex(41))))+m.e(10)*(1/ex(6)+1/ex(12)*(-17+3*sqrt(ex(41))))+m.e(13);
		cout<<"If v="<<v<<endl;
		ex fv;
		for (int j=1;j<=m.Dimension();++j)
		{
			ex act=V.Action<DifferentialForm>(v,fourform);
			act=V.Action<DifferentialForm>(m.e(j),act);
			fv+=(m.e(j)*NormalizeRoots(HodgeStar(act*g.StarGamma()))).expand();
		}
		cout<<"then Hv = "<<fv<<endl;

		VectorSpace<DifferentialForm> Spanv(&v,&v+1);
		if (Spanv.Contains(fv))
			cout<<"so v is a eigenvector with eigenvalue"<<Spanv.Components(fv)[0]<<endl;
		else
			cout<<"not an eigenvector!"<<endl;
		cout<<endl;
	}
};

void PrintStabilizer() {
	ManifoldWith<SO3StructureNoRoots> M;
	GL G(9);
	GLRepresentation<VectorField> V(&G,M.e());
	cout<<"Stabilizer of "<<M.StarGamma()<<":"<<endl;
	AbstractLieSubgroup<false> H=Stabilizer<DifferentialForm>(G,V,M.StarGamma());
	cout<<H<<endl;
	cout<<"Killing form:"<<endl;
	matrix kill=H.KillingForm();
	cout<<kill<<endl;
}


int main()
{
	cout<<latex;
	//Grassmannian grass;
	
	DisplaySO3Representation(8);
	PrintStabilizer();

	cout<<"The invariant 4-form under the irreducible action of SO(3) on R^9 is"<<endl;
	ConcreteManifold M(9);
	SO3Representation<VectorField> V(M.e());
	ex gamma=NormalizeRoots(3*Invariant_pForms(V,M.e(),4).e(1)/8);
	cout<<gamma<<endl;
	return 0;
}
