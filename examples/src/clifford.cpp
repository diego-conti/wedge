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
#include "wedge/wedge.h"

using namespace Wedge;
using namespace std;

void printCliffordMultiplicationTable(int dim)
{
	cout<<dflt;	//normal output mode
	ConcreteManifold M(dim);
	RiemannianStructure g(&M,M.e());
	cout<<"Clifford multiplication table for dimension "<<dim<<endl;
	exvector v(g.DimensionOfSpinorRepresentation());
	for (int i=1;i<=dim;i++)
	{
		for (int k=0;k<g.DimensionOfSpinorRepresentation();k++)
		{
			v[k]+=M.e(i)*g.CliffordDot(M.e(i),g.u(k));
			cout<<M.e(i)<<"."<<g.u(k)<<"="<<g.CliffordDot(M.e(i),g.u(k))<<endl;
		}
	}
	cout<<latex<<v<<endl;	//latex output mode
}

int main()
{
	printCliffordMultiplicationTable(2);
	printCliffordMultiplicationTable(3);
	return 0;
}
