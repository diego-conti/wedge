/* (C) 2008 by Diego Conti, diego.conti@unimib.it              *
 *                                                                         *
 *   This file is part of Dictionary.                                      *
 *                                                                         *
 *   Dictionary is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 3 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   Dictionary is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with Dictionary; if not, write to the                           *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include "dictionary/dictionaries/SU3SO3Dictionary.h"
#include "dictionary/parsedictionary.h"
#include "parambasis.h"
#include <wedge/wedge.h>

TrigonometricSimplifier trig;

class SO3R2 : public ManifoldWithCoordinates {
	exvector CreateFrame() {
		exvector l;
		for (int i=1;i<=3;++i) l.push_back(DifferentialOneForm(N.e(i)));
		return l;		
	}
	ex r;
public:
	SO3R2 () : ManifoldWithCoordinates(CreateFrame(),2) {
		Declare_d(e(1),e(2)*e(3));
		Declare_d(e(2),e(3)*e(1));
		Declare_d(e(3),e(1)*e(2));
		r=Function(N.r);
		Declare_d(r,(x(1)*d(x(1))+x(2)*d(x(2)))/r);

		Frame f;
		f.push_back(r*(d(x(1))+x(2)*e(1)));
		f.push_back(r*(d(x(2))-x(1)*e(1)));
		f.push_back(cos(r)*e(2));
		f.push_back(cos(r)*e(3));
		f.push_back(e(1));
		TransverseRiemannianStructure P(this, f,4);
		

		
		
		
		TransverseLeviCivitaConnection omega(this,P);
		cout<<omega.AsMatrix()<<endl;
		cout<<omega.Nabla<DifferentialForm>(e(1),e(2))<<endl;
		cout<<static_cast<PseudoLeviCivitaConnection>(omega).CurvatureForm()<<endl;
		matrix r=omega.BaseRicciAsMatrix(true);
		matrix r2=omega.BaseRicciAsMatrix(false);
		cout<<ex(r).expand().normal()<<endl;
		cout<<ex(r2).expand().normal()<<endl;
		for (int i=0;i<3;++i)
		for (int j=0;j<3;++j)
			if (!(r(i,j)-r2(i,j)).expand().numer().is_zero())
				cout<<i<<j<<endl;

	}
	
};

void DictionaryAsLocal()
{
	SU3SO3Dictionary Y(false,1);	
	Y.CreateAlgebra();
	ex r=Y.RadialCoordinate();
	
	lst subs;

	ExVector pullbackframe=Y.PullBackFrame();
	exvector f;
	f.push_back(pullbackframe(6));
	f.push_back(sqrt(2-2*cos(r))/r*pullbackframe(7));
	f.push_back(sqrt(2-2*cos(r))/r*pullbackframe(8));
	f.push_back(sqrt(2+2*cos(r))*pullbackframe(1));
	f.push_back(sqrt(2+2*cos(r))*pullbackframe(2));
	f.push_back(sqrt(ex(2))*cos(r)*pullbackframe(3));
	f.push_back(sqrt(ex(2))*cos(r)*pullbackframe(4));
	f.push_back(2*pullbackframe(5));
	for (exvector::iterator i=f.begin();i!=f.end();++i) *i=i->subs(subs);

	Frame frame(f.begin(),f.end());	
	LOG_INFO(MaxLength(10000)<<frame);

/*	coefficients.push_back(sqrt(2*F(r))/r);
	coefficients.push_back(sqrt(2*J(r)));
	coefficients.push_back(2*M(r));
	coefficients.push_back(2*(2/(2*J(r)*F(r)+M(r)*M(r)*(F(r)+J(r)))));
*/
	TransverseRiemannianStructure g=Y.InvariantStructure(frame);

	TransverseLeviCivitaConnection P(&Y.GtimesR,g);
	cout<<latex<<OmitArgument(r);
	matrix r2=P.BaseRicciAsMatrix(true);
	matrix r3=P.BaseRicciAsMatrix(false);
	assert(ex(r2-r3).expand().is_zero());
	for (int i=0;i<=7;++i)
	for (int j=0;j<=7;++j)
	{
		ex Ric_ij=trig.Simplify((r2(i,j)));
		if (!Ric_ij.is_zero())
			cout<<"\\ric_{"<<i+1<<j+1<<"}="<<(Ric_ij).normal().expand().normal()<<"\\"<<endl;
		ex Ric_ij2=trig.Simplify((r3(i,j)));
		if (!(Ric_ij2-Ric_ij).expand().is_zero())
			cout<<"BUT \\ric_{"<<i+1<<j+1<<"}="<<(Ric_ij2).normal().expand().normal()<<"\\"<<endl;
		

	}

	for (int i=0;i<=7;++i)
	for (int j=i+1;j<=7;++j)
		if ((r2(i,i))==(r2(j,j))) cout<<"\\ric_{"<<i+1<<i+1<<"}="<<"\\ric_{"<<j+1<<j+1<<"}"<<endl;

	ExVector diag;
	for (int i=0;i<=7;++i)
		diag.push_back(r2(i,i).normal().expand().normal());
	cout<<"\\pi(Ric)_1="<<diag(1)-diag(8)+diag(2)+diag(3)-diag(7)-diag(5)<<endl;
	cout<<"\\pi(Ric)_2="<<diag(6)+diag(4)-diag(7)-diag(5)<<endl;

}

int main()
{
	cout<<latex;
	DictionaryAsLocal();
	return 0;
}
