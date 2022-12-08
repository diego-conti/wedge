/***************************************************************************
 *   Copyright (C) 2008, 2009 by Diego Conti				   *
 *   diego.conti@unipi.it                                                 *
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

#include "wedge/wedge.h"

using namespace Wedge;
using namespace std;

struct X : public ManifoldWith<RiemannianStructure> {
	X() : ManifoldWith<RiemannianStructure>(4) {
		for (int i=1;i<=4;++i)
			DeclareNabla<Spinor>(e(i),u(0),0);
		cout<<latex<<"Generic 4-manifold with a parallel spinor (hyperkahler)"<<endl;
		cout<<"de^{12}+de^{34}="<<d(e(1)*e(2)+e(3)*e(4))<<endl;
		cout<<"de^{13}+de^{42}="<<d(e(1)*e(3)+e(4)*e(2))<<endl;
		cout<<"de^{14}+de^{23}="<<d(e(1)*e(4)+e(2)*e(3))<<endl;
		cout<<LeviCivita();

		TorsionZeroImpliesInvolutive();
		InvolutiveImpliesTorsionZero();
	}
	X(bool symplectic) : ManifoldWith<RiemannianStructure>(4) {
		Declare_d(e(1)*e(2)+e(3)*e(4),0);
		Declare_d(e(1)*e(3)+e(4)*e(2),0);
		Declare_d(e(1)*e(4)+e(2)*e(3),0);
		cout<<"Generic symplectic 4-manifold"<<endl;
		cout<<"de^{12}+de^{34}="<<d(e(1)*e(2)+e(3)*e(4))<<endl;
		cout<<"de^{13}+de^{42}="<<d(e(1)*e(3)+e(4)*e(2))<<endl;
		cout<<"de^{14}+de^{23}="<<d(e(1)*e(4)+e(2)*e(3))<<endl;
		TorsionZeroImpliesInvolutive();
		InvolutiveImpliesTorsionZero();
	}

	void InvolutiveImpliesTorsionZero() {
		cout<<"Assuming the lagrangian distributions <e1,e3>, <e2,e4> are involutive:"<<endl;
		DeclareZero(Hook(e(2),LieBracket(e(1),e(3))));
		DeclareZero(Hook(e(4),LieBracket(e(1),e(3))));
		DeclareZero(Hook(e(1),LieBracket(e(2),e(4))));
		DeclareZero(Hook(e(3),LieBracket(e(2),e(4))));
		for (int i=1;i<=4;++i)
			for (int j=i+1;j<=4;++j)
				cout<<"["<<e(i)<<","<<e(j)<<"]="<<LieBracket(e(i),e(j))<<endl;
		
		Connection omega(this,e(),N.Gamma.Prime());	
		for (int k=1;k<=4;++k) {
			omega.DeclareNabla<DifferentialForm>(e(k),e(1)*e(2)+e(3)*e(4),0);
			for (int i=1;i<=4;++i)
			for (int j=i%2+1;j<=4;j+=2)
				omega.DeclareZero(Hook(e(j),omega.Nabla<DifferentialForm>(e(k),e(i))));
		}
		for (int k=1;k<=4;++k)
		for (int i=k%2+1;i<=4;i+=2)
		for (int j=k%2+1;j<=4;j+=2)
			omega.DeclareZero(Hook(e(j),omega.Nabla<VectorField>(e(k),e(i))-LieBracket(e(k),e(i))));
		cout<<"Then the canonical symplectic connection"<<endl;	
		cout<<omega<<endl;
		exvector T=omega.Torsion();
		cout<<"has torsion"<<endl<<T<<endl<<endl;		
		exvector zerovector(4);
		assert(T==zerovector);
	}
	void TorsionZeroImpliesInvolutive() {
		Connection omega(this,e(),N.Gamma.Prime());	
		for (int k=1;k<=4;++k) {
			omega.DeclareNabla<DifferentialForm>(e(k),e(1)*e(2)+e(3)*e(4),0);
			for (int i=1;i<=4;++i)
			for (int j=i%2+1;j<=4;j+=2)
				omega.DeclareZero(Hook(e(j),omega.Nabla<DifferentialForm>(e(k),e(i))));
		}
		for (int k=1;k<=4;++k)
		for (int i=k%2+1;i<=4;i+=2)
		for (int j=k%2+1;j<=4;j+=2)
			omega.DeclareZero(Hook(e(j),omega.Nabla<VectorField>(e(k),e(i))-LieBracket(e(k),e(i))));


		cout<<"The bilagrangian splitting <e1,e3> + <e2,e4> makes it into a bilagrangian manifold iff canonical symplectic connection"<<endl;	
		cout<<omega<<endl;
		cout<<"has zero torsion, i.e."<<endl;
		exvector T=omega.Torsion();
		cout<<"0="<<T<<endl;
		cout<<"Imposing that the torsion is indeed zero, we obtain"<<endl;
		DeclareZero(T.begin(),T.end());
		for (int i=1;i<=4;++i)
		for (int j=i+1;j<=4;++j)
			cout<<"["<<e(i)<<","<<e(j)<<"]="<<LieBracket(e(i),e(j))<<endl;
		assert(Hook(e(2),LieBracket(e(1),e(3)))==0);
		assert(Hook(e(4),LieBracket(e(1),e(3)))==0);
		assert(Hook(e(1),LieBracket(e(2),e(4)))==0);
		assert(Hook(e(3),LieBracket(e(2),e(4)))==0);
		cout<<"implying that the distributions <e1,e3> and <e2,e4> are involutive"<<endl<<endl;
	}
};

int main()
{
	X M;
	X N(true);
	return 0;
}
