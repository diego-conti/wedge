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
#include <list>
#include <vector>
#include <ginac/ginac.h>
#include "wedge/torsionfreeconnection.h"
#include "wedge/coordinates.h"
#include "wedge/riemannianstructure.h"
#include "wedge/liegroup.h"
#include "wedge/tensor.h"
#include <iostream>

using namespace Wedge;
using namespace std;

struct Kodaira : public ConcreteManifold , public Has_dTable {
	Kodaira() : ConcreteManifold(4) {
		Declare_d(e(1),-e(4)*e(2));
		Declare_d(e(2),e(4)*e(1));
		Declare_d(e(3),e(1)*e(2));
		Declare_d(e(4),0);
	}
};

struct X : public ConcreteManifold , public Has_dTable {
	X() : ConcreteManifold(4) {
		Declare_d(e(1),0);
		Declare_d(e(2),0);
		Declare_d(e(3),e(1)*e(2));
		Declare_d(e(4),e(1)*e(3));
	}
};


int main()
{
	struct AlmostComplex : X {
		TorsionFreeConnection<true> omega;
		ex J(ex Y) {return Hook(Y,e(1)*e(2)+e(3)*e(4));}
		ex A(ex X,ex Y)	{
			return omega.Nabla<VectorField>(X,J(Y))-J(omega.Nabla<VectorField>(X,Y));
		}
		ex Q(ex X, ex Y) {
			return (A(J(Y),X)+J(A(Y,X))+2*J(A(X,Y)))/4;
		}
		AlmostComplex() : omega(this,e()) {//RiemannianStructure(this,e())) {
			Connection k(this,e());
			for (int i=1;i<=4;++i)
			for (int j=1;j<=4;++j)
				k.DeclareNabla<VectorField>(e(i),e(j),
					omega.Nabla<VectorField>(e(i),e(j))-Q(e(i),e(j)));
				cout<<k<<endl<<k.Torsion()<<endl;
	
			matrix nablaj(4,4);
			for (int i=1;i<=4;++i)
			for (int j=1;j<=4;++j)
				nablaj(i-1,j-1)=k.Nabla<VectorField>(e(i),J(e(j)))-J(k.Nabla<VectorField>(e(i),e(j)));
			cout<<nablaj<<endl;
		}
	};

	cout<<latex;
	AlmostComplex M;
	cout<<endl;
	return 0;
}
