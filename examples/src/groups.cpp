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

template<bool HasParams> void PrintStructureConstantsWrtFrame(const LieGroupHasParameters<HasParams>& G, const exvector& newframe)
{
	AbstractLieSubgroup<HasParams> H(G,newframe);
	cout<<"Structure constants with respect to new frame"<<endl;
	cout<<H<<endl;
	cout<<"Structure constants with respect to new frame (alternative format)"<<endl;
	cout<<H.StructureConstants()<<endl;
}

class ALieGroupWithParameters : public ConcreteManifold, public LieGroupHasParameters<true> {
	public:
	ALieGroupWithParameters() : ConcreteManifold(7)
	{
		ex a=StructureConstant(N.a);								
		Has_dTable::Declare_d(e(1),a*e(1)*e(4));		 
		Has_dTable::Declare_d(e(2),a*e(2)*e(4));
		Has_dTable::Declare_d(e(3),2*a*e(1)*e(2)+2*a*a*e(4)*e(7)-2*a*e(5)*e(6));
		Has_dTable::Declare_d(e(4),0);
		Has_dTable::Declare_d(e(5),a*e(4)*e(5));
		Has_dTable::Declare_d(e(6),a*e(4)*e(6));
		Has_dTable::Declare_d(e(7),-2*(e(1)*e(2)+e(3)*e(4)+e(5)*e(6)));

		exvector eqns;
		GetEquations_ddZero(eqns);
		ex zero[1];	//a vector with one element, which is equal to zero
		if (!(includes(zero,zero+1,eqns.begin(),eqns.end())))
			throw std::logic_error("d^2!=0");
	}
};


class Iwasawa : public AbstractLieGroup<> {
public:
	Iwasawa() : AbstractLieGroup<>("0,0,0,0,13+42,14+23")
	{
		Basis<DifferentialForm> b;
		for (int i=1;i<=6;++i)
		for (int j=i+1;j<=6;++j)
			b.push_back(d(e(i)*e(j)));
		cout<<latex<<b;
		cout<<b.Components(d(e(4)*e(5)))<<endl;
	}
};

int main()
{

	Iwasawa M;
	ALieGroupWithParameters G;	//define a certain family of Lie groups
	ExVector f=ParseDifferentialForms(G.e(),"-4,3,-2,1,-5,-6,7");	//we don't like the frame, so we choose another
	PrintStructureConstantsWrtFrame(G,f);		
	GStructureHasParameters<SU3StructureDim7,true> P(&G,f);	//the SU(3)-structure determined by the frame

	cout<<"dF = "<<G.d(P.F())<<endl;
	cout<<"d(alpha\\wedge\\Omega) = "<<G.d(P.alpha()*P.Omega())<<endl;
}
