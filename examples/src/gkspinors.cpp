/***************************************************************************
 *   Copyright (C) 2007, 2008, 2009 by Diego Conti			   *
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
#include <wedge/wedge.h>

using namespace Wedge;
using namespace std;

void CheckHypo()
{
	const int dimension=7;
	ConcreteManifold M(dimension);
	RiemannianStructure g(&M, M.e());
	RiemannianConnection omega(&M,g);
	cout<<"Almost-contact SU(3)-structure defined by a generalized Killing spinor"<<endl;
	matrix A(dimension,dimension);
	for (int i=0;i<dimension;i++)
		for (int j=i;j<dimension;j++)
			A(j,i)=A(i,j)=realsymbol("A"+ToString(i+1,j+1));
	for (int i=1;i<=dimension;i++)
	{
		ex nablapsi;
		for (int j=1;j<=dimension;j++)
		{
			nablapsi+=g.CliffordDot(M.e(j)*A(i-1,j-1),g.u(0));
		}
		omega.DeclareNabla<Spinor>(M.e(i),g.u(0),nablapsi/2);	//declare that psi is a generalized Killing spinor
	}
	ex alpha=M.e(dimension);
	ex F=ParseDifferentialForm(M.e(),"12+34+56");
	ex OmegaPlus=ParseDifferentialForm(M.e(),"135-146-236-245");
	ex OmegaMinus=ParseDifferentialForm(M.e(),"136+235+145-246");
	cout<<"Defining forms:"<<endl;
	cout<<"\\alpha="<<alpha<<endl;
	cout<<"F="<<F<<endl;
	cout<<"Re \\Omega="<<OmegaPlus<<endl;
	cout<<"Im \\Omega="<<OmegaMinus<<endl;

	cout<<"Covariant derivatives in terms of intrinsic torsion tensor A:"<<endl;
	for (int i=1;i<=dimension;i++)
	{
		cout<<"\\nabla_"<<i<<"\\alpha="<<omega.Nabla<DifferentialForm>(M.e(i),alpha)<<endl;
		cout<<"Re \\nabla_"<<i<<"\\Omega="<< omega.Nabla<DifferentialForm>(M.e(i),OmegaPlus)<<endl;
		cout<<"Im \\nabla_"<<i<<"\\Omega="<< omega.Nabla<DifferentialForm>(M.e(i),OmegaMinus)<<endl;
	}
}

int main()
{
	CheckHypo();
}
