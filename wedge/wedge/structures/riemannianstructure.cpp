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

#include "../structures/riemannianstructure.h"

#include "../structures/spinor.h"

namespace Wedge {

ex RiemannianStructure::Apply(const DifferentialForm& v, const DifferentialForm& w) const
{
	if (v.nops()!=w.nops()) return 0;
	else {
		matrix m(w.nops(),w.nops());
		for (int i=0;i<w.nops();i++)
			for (int j=0;j<w.nops();j++)
				m(i,j)=	Apply(ex_to<DifferentialOneForm>(v.op(i)),ex_to<DifferentialOneForm>(w.op(j)));
		return m.determinant();
	}
}


ex RiemannianStructure::Apply (const VectorField& X, const Spinor& spinor) const
{
		const DifferentialOneForm& alpha=static_cast<const DifferentialOneForm&>(X);
		int dimension=e().size();
		ex result=0;
		ExVector components=e().Components(alpha);
		for (int i=1;i<=dimension;i++)
		{			
			if (components(i)!=0)
			{
				vector<bool> a=spinor.a;
				ex coeff;
				if (dimension%2==0 || i!=dimension) {
					int j=i/2;
					if (i%2==0) {
						bool sign=false;
						for (int h=0;h<j;h++) sign=sign xor a[h];
						a[j-1]=!a[j-1];
						coeff=sign ? 1 : -1;
					}
					else {
						bool sign=false;
						for (int h=0;h<j;h++) sign=sign xor a[h];						
						a[j]=!a[j];
						coeff=sign? -I : I;
					}
					if (dimension%2!=0) coeff=-coeff;
				}
				else coeff=1;

				if (dimension%2!=0) {
					//Clifford-multiply by e(dimension)
					bool sign=false;
					for (int h=0;h<a.size();h++) sign=sign xor a[h];
//At this point, we can take as coefficient \pm I sign.
//My thesis had the minus sign, in order that the complex volume as defined in Michelson-Lawson act identically.
//The choice to follow is consistent with
//	 Friedrich - Kim : The Einstein-Dirac equation on Riemannian spin manifolds					
					coeff*=(sign xor (dimension%4==3))? -I : I;					
				}
				result+=components(i)*coeff*ex(Spinor(a));
			}
		}
		return result;					
}



ex RiemannianStructure::u(ZeroBased n) const
{
	if (n<0) throw OutOfRange(__FILE__,__LINE__,n);
	return Spinor(n,M()->Dimension());
}

int RiemannianStructure::DimensionOfSpinorRepresentation() const
{
	return 1<<(M()->Dimension()/2);
}

ex RiemannianStructure::Apply (const VectorField& left,const VectorField& right) const
{
	ex scalarproduct;
	exvector comps1=e().Components(static_cast<const DifferentialOneForm&>(left));
	exvector comps2=e().Components(static_cast<const DifferentialOneForm&>(right));
	assert(comps1.size()==comps2.size());
	for (exvector::const_iterator i=comps1.begin(),j=comps2.begin();i!=comps1.end();i++,j++)
	{
		scalarproduct+=*i * *j;
	}
	return scalarproduct;	
}

ex RiemannianStructure::Apply (const Spinor& left,const Spinor& right) const
{
	return left==right? 1 : 0;
}

ex RiemannianStructure::HodgeStar(ex form1) const
{
	ex form=e()[0];
	for (int i=1;i<e().size();i++)
		form*=e()[i];
	return Hook(form1,form);
}

ex RiemannianStructure::Hook(ex left, ex right) const
{
	ex result=internal::RiemannianHookOperator::BilinearOperator(left.expand(),right.expand(),&hookOperator).expand();
		//adjust sign so that e^{12..n}\hook e^{12..n}=1
	return (Degree<DifferentialForm>(left) % 4<2)? result : -result;
}	

template<> ex RiemannianStructure::ScalarProduct<VectorField> (ex op1, ex op2) const
{
	ex result;
	for (exvector::const_iterator i=e().begin();i!=e().end();++i)
		result+=Wedge::Hook(op1,*i)*Wedge::Hook(op2,*i);
	return result;
}

namespace internal {
ex RiemannianHookOperator::Apply(const VectorField& left,const VectorField& right) const
{
	return structure->Apply(left,right);
}
}

}
