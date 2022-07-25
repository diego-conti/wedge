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

//Check whether gamma is E-stable wrt a "standard" codimension one subspace
//See [D. Conti, Embedding into manifolds with torsion, http://arxiv.org/abs/0812.4186v1 ]
void EStable(const Manifold& M, ex gamma) {
	cout<<latex<<"Verifying whether the form "<<gamma<<" on M^"<<M.Dimension()<<" is E-stable"<<endl;
	gamma=gamma.expand();
	for (int k=1;k<=M.Dimension();++k)
		{
			VectorSpace<DifferentialForm> b;
			for (int i=1;i<=M.Dimension();++i)
			for (int j=1;j<=M.Dimension();++j)
			{
				LinearAction<VectorField> a(lst{M.e(i)==M.e(j)});
				b.AddGenerator(AlgebraAction<VectorField,DifferentialForm>(a,gamma).subs(M.e(k)==0));
			}
			//dimension of Lambda^pE
			ex Lambda_pE=(binomial(M.Dimension()-1,Degree<DifferentialForm>(gamma)));
			if (b.Dimension()==Lambda_pE)
				cout<<"is E-stable wrt E="<<M.e(k)<<"^\\perp : \\dim\\Lambda^pE="<<b.Dimension()<<endl;
			else
				cout<<"not E-stable wrt E="<<M.e(k)<<"^\\perp : \\dim\\Lambda^pE="<<Lambda_pE<<">"<<b.Dimension()<<endl;
		}
}

//A structure in 10 dimensions defined by a "democratic" 6-form
//See [Chandrashekar Devchand, Jean Nuyts, Gregor Weingart, Matryoshka of Special Democratic Forms, http://arxiv.org/abs/0812.3012v1 ]
class SU4U1Structure : public RiemannianStructure {
public:
	enum {dimension=10};	///< The dimension of the manifold
/** @brief Define an \f$\rm SU(4)U(1)\f$-structure on a 10-manifold, represented by the choice of a frame
 *  @param manifold The manifold on which the \f$\rm SU(4)U(1)\f$-structure is defined.
 *  @param frame The adapted frame defining the structure.
 * 
 *  @warning Caller must ensure that the pointer remains valid. 
 */
	SU4U1Structure(const Manifold* manifold, const Frame& frame): RiemannianStructure(manifold, frame) {}
	ex Omega() const;					///< The defining 6-form
	ex StarOmega() const {return HodgeStar(Omega());} 	///< The defining 4-form
protected:
	SU4U1Structure(const Manifold* manifold) : RiemannianStructure(manifold) {}
};

ex SU4U1Structure::Omega() const {
	ExVector sigma(dimension);
	ex result;
	ex spin7=ParseDifferentialForm(e(),"1234+1256+1278+1357+1386+1485+1476+2385+2376+2475+2468+3456+3478+5678")*e(9)*e(10);
	LOG_INFO(spin7);
	for (int k=0;k<=8;k+=2)
	{
		for (int i=0;i<dimension;++i)
			sigma[i]=e()[(i+k)%10];
		LOG_INFO(sigma);
		LinearAction<VectorField> g(e(),sigma);
		result+=GroupAction(g,spin7);
		LOG_INFO(result);
	}
	return result;
}

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


int main()
{
	{
	SU G(4);
	cout<<latex<<G.KillingForm()<<endl;
	EStable(G,G.HodgeStar(G.ThreeForm()));
	}



	{
	SU G(3);
	cout<<latex<<G.KillingForm()<<endl;
	cout<<G<<endl;
	cout<<G.ThreeForm()<<endl;
	EStable(G,G.ThreeForm());

	GL gl(8);
	GLRepresentation<VectorField> V(&gl,G.e());
	AbstractLieSubgroup<false> H=Stabilizer<DifferentialForm>(gl,V,G.ThreeForm());
	cout<<H<<endl;
	cout<<"Killing form:"<<endl;
	cout<<H.KillingForm()<<endl;
	}

	ConcreteManifold X(9);
	EStable(X,SO3StructureNoRoots(&X,X.e()).StarGamma());

	ManifoldWith<SU4U1Structure> M;
	EStable(M,M.Omega());

	ConcreteManifold Sp2Sp1(8);
	ExVector omega=ParseDifferentialForms(Sp2Sp1.e(),"12+34+56+78,13+42+57+86,14+23+58+67");
	EStable(Sp2Sp1,(omega(1)*omega(1)+omega(2)*omega(2)+omega(3)*omega(3)));
}
