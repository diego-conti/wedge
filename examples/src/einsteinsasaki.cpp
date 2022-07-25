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

//Einstein-Sasaki manifolds in dimension 5
class ES5Manifold : public ManifoldWith<SU2Structure> {
public:
	ES5Manifold()
	{
		for (int i=1;i<=5;i++)			
			DeclareNabla<DifferentialForm>(e(i),alpha(),-Hook(e(i),omega1()));
		Declare_d(alpha(),-2*omega1());
		Declare_d(omega2(),3*alpha()*omega3());
		Declare_d(omega3(),-3*alpha()*omega2());

		DoCalculations();
	}

	void DoCalculations() {
		cout<<"Einstein-Sasaki 5-manifolds have a Killing spinor:"<<endl;
		for (int i=1;i<=5;++i)
		{
			cout<<"\\nabla_{"<<e(i)<<"}"<<psi()<<"=";
			cout<<Nabla<Spinor>(e(i),psi()).expand()<<endl;
			cout<<"= -1/2"<<e(i)<<"\\cdot "<<psi()<<"=";
			cout<<-CliffordDot(e(i),psi())/2<<endl;
		}
	}
//In theory the following code should print the Ricci tensor of a ES 5-manifold, which is 4 times the identity. In practice it will probably eat up all memory and crash your system...
	void printRicci() {
		cout<<"Ricci tensor of generic Einstein-Sasaki 5-manifold:"<<endl;
		ex Ric=LeviCivita().Ricci();
		PolyBasisImplementation<Function,DefaultPolyAlgorithms> I;
		LeviCivita().GetEquations_ddZero(I);
		I.Reduce();
		cout<< I.ReduceModuloIdeal<Tensor<DifferentialOneForm,DifferentialOneForm> >(Ric)<<endl;
	}

};


//Hyperkaehler 4-manifolds
struct HyperKaehlerFourManifold : public ConcreteManifold {
	HyperKaehlerFourManifold() : ConcreteManifold(4) {
		DoCalculations();
	}
	
	//same as ES5Manifold::printRicci applies here
	void printRicci() {
		cout<<"Generic hyperkaehler 4-manifold:"<<endl;
		RiemannianStructure g(this,e());
		ExVector omega=ParseDifferentialForms(e(),"12+34,13+42,14+23");
		LeviCivitaConnection<false> leviCivita(this,g);
		leviCivita.Declare_d(omega(1),0);
		leviCivita.Declare_d(omega(2),0);
		leviCivita.Declare_d(omega(3),0);
		cout<<"Connection form: "<<endl;
		cout<<leviCivita.AsMatrix()<<endl;
		for (int j=0;j<3;++j)
		for (int i=1;i<=4;++i)
			cout<<"\\nabla_{"<<e(i)<<"}"<<g.u(j)<<"="<<leviCivita.Nabla<Spinor>(e(i),g.u(j))<<endl;
		cout<<endl;
		PolyBasisImplementation<Function,DefaultPolyAlgorithms> I;
		leviCivita.GetEquations_ddZero(I);
		I.Reduce();
		cout<<"Ricci tensor: "<<endl;
		cout<<I.ReduceModuloIdeal<Tensor<DifferentialOneForm,DifferentialOneForm> >(leviCivita.Ricci())<<endl<<endl;
	}

	void DoCalculations()
	{
		cout<<"Generic 4-manifold with a parallel spinor (hyperkahler):"<<endl;
		RiemannianStructure g(this,e());
		LeviCivitaConnection<false> leviCivita(this,g);
		for (int i=1;i<=4;++i)
			leviCivita.DeclareNabla<Spinor>(e(i),g.u(0),0);
		cout<<":"<<leviCivita.d(leviCivita.d(e(1)))<<endl;
		cout<<"Connection form: "<<endl;
		cout<<leviCivita.AsMatrix()<<endl;
		ExVector omega=ParseDifferentialForms(e(),"12+34,13+42,14+23");
		for (int i=1;i<=3;++i)
			cout<<"d"<<omega(i)<<"="<<leviCivita.d(omega(1))<<endl;


	}
};


int main()
{
	cout<<latex;	//switch to latex output
	HyperKaehlerFourManifold X;
	ES5Manifold M;
	return 0;
}
