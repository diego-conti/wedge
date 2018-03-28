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
#include <list>
#include <vector>
#include <ginac/ginac.h>
#include "wedge/torsionfreeconnection.h"
#include "wedge/coordinates.h"
#include "wedge/riemannianstructure.h"
#include "wedge/liegroup.h"
#include "wedge/tensor.h"
#include "wedge/transversestructure.h"
#include "wedge/transverseconnection.h"
#include <iostream>

using namespace Wedge;
using namespace std;

struct S3TimesS1 : public ConcreteManifold , public Has_dTable {
	S3TimesS1() : ConcreteManifold(4) {
		Declare_d(e(1),e(2)*e(3));
		Declare_d(e(2),e(3)*e(1));
		Declare_d(e(3),e(1)*e(2));
		Declare_d(e(4),0);
	}
};

struct SomeCalculations : public S3TimesS1 {
 	SomeCalculations() {	
		cout<<"Some calculations on S^3\\times S^1"<<endl;
 		cout<<e(1)*(e(1)+e(2))<<endl;
 		cout<<d(e(1)+e(4))<<endl;
 		cout<<d(e(1)*e(2))<<endl;
		cout<<LieDerivative(e(1),e(2)*e(4))<<endl;	
 	}
};

struct HermitianConnection : public S3TimesS1 {

	HermitianConnection() {
		cout<<"A flat connection on S^3\\times S^1"<<endl;
		{
			Connection omega(this,e());
			for (int i=1;i<=4;++i)
				for (int j=1;j<=4;++j)
					omega.DeclareNabla<DifferentialForm>(e(i),e(j),0);	
			cout<<omega.Torsion();
			cout<<omega.CurvatureForm();
		}
		cout<<"A hermitian connection on S^3\\times S^1"<<endl;
		RiemannianStructure g(this,e());
		{
			RiemannianConnection omega(this,g);
			for (int i=1;i<=4;++i)
			{
				omega.DeclareNabla<DifferentialForm>(e(i),e(1)*e(2)+e(3)*e(4),0);
				omega.DeclareNabla<DifferentialForm>(e(i),e(4),0);
			}
			cout<<omega.Torsion();
			cout<<omega.CurvatureForm();
		}
		cout<<"The Levi-Civita connection on S^3\\times S^1"<<endl;
		{
			LeviCivitaConnection<true> omega(this,g);
			cout<<omega.Torsion();
			cout<<omega.CurvatureForm();

	
		}
	}
};

//the cotangent bundle over the two-sphere (insanely slow)
struct TS2 : public ManifoldWithCoordinates {
	TS2() : ManifoldWithCoordinates(4) {
		cout<<"Cotangent bundle over S^2"<<endl;
		ExVector e(4);
		e(1)=2/(pow(x(1),2)+pow(x(2),2)+1)*d(x(1));
		e(2)=2/(pow(x(1),2)+pow(x(2),2)+1)*d(x(2));

		ex omega=x(2)*e(1)-x(1)*e(2);
		e(3)=d(x(3))-x(4)*omega;
		e(4)=d(x(4))+x(3)*omega;

		RiemannianStructure flatMetric(this,e);	
		cout<<"Curvature of flat metric"<<endl;	
		cout<<LeviCivitaConnection<true>(this,flatMetric).CurvatureForm()<<endl;
		
		ex rsquared=pow(x(3),2)+pow(x(4),2);
		e(1)*=pow((1+rsquared),ex(1)/4);
		e(2)*=pow((1+rsquared),ex(1)/4);
		e(3)*=pow((1+rsquared),-ex(1)/4);
		e(4)*=pow((1+rsquared),-ex(1)/4);

		RiemannianStructure EguchiHanson(this,e);		
		cout<<"Curvature of Eguchi-Hanson metric"<<endl;	
		cout<<LeviCivitaConnection<true>(this,EguchiHanson).CurvatureForm()<<endl<<endl;		
	}
}; 

struct S3 : public ConcreteManifold , public Has_dTable {
	S3() : ConcreteManifold(3) {				
		Declare_d(e(1),e(2)*e(3));
		Declare_d(e(2),e(3)*e(1));
		Declare_d(e(3),e(1)*e(2));
	
		printCurvature();
	}
		
	void printCurvature() {
		cout<<"S^3 as a Lie group, standard basis, using dtable:"<<endl;
		LeviCivitaConnection<true> leviCivita(this,RiemannianStructure(this,e()));
		cout<<"Curvature tensor of Levi-Civita connection: "<<endl;
		cout<<leviCivita.CurvatureForm()<<endl;
		cout<<"Ricci tensor of Levi-Civita connection: "<<endl;
		cout<<leviCivita.Ricci()<<endl<<endl;
	}	
};

struct ThreeManifold : public ConcreteManifold {
	ThreeManifold() : ConcreteManifold(3) {
		printCurvature();
	}
	
	void printCurvature() {	
		cout<<"S^3 as a Lie group, non-standard basis, no dtable"<<endl;
		exvector v;
		v.push_back(e(1)+e(2));
		v.push_back(e(2));
		v.push_back(e(3)+2*e(1));
		Frame e(v);
		cout<<"In the basis"<<e<<endl;

//		TorsionFreeConnection<false> connection(this,e);
//		connection.Declare_d(e(1),e(2)*e(3));
//		connection.Declare_d(e(2),e(3)*e(1));
//		connection.Declare_d(e(3),e(1)*e(2));
//		cout<<"Ricci tensor of a generic torsion-free connection: "<<endl;
//		cout<<connection.Ricci()<<endl;	
	
		LeviCivitaConnection<false> leviCivita(this,RiemannianStructure(this,e));
		leviCivita.Declare_d(e(1),e(2)*e(3));
		leviCivita.Declare_d(e(2),e(3)*e(1));
		leviCivita.Declare_d(e(3),e(1)*e(2));
		cout<<"Curvature tensor of Levi-Civita connection: "<<endl;
		cout<<leviCivita.CurvatureForm()<<endl;
		cout<<"Ricci tensor of Levi-Civita connection: "<<endl;
		cout<<leviCivita.Ricci()<<endl<<endl;
	}	
};



//compute the Levi-Civita connection for a left-invariant metric on a Lie group
matrix LeviCivitaOfLieGroup(const char* struct_constants)
{
	AbstractLieGroup<> G(struct_constants);
	RiemannianStructure g(&G,G.e());
	LeviCivitaConnection<true> omega(&G,g);
	return omega.AsMatrix();
}


class CP2 : public AbstractLieGroup<> {
public:
	CP2 () : AbstractLieGroup<>("-23-45+2*67,13+46-57-[sqrt(3)]*58,-12-47+[sqrt(3)]*48-56,"
					"15-26+37-[sqrt(3)]*38,"
					"-14+27+36+[sqrt(3)]*28,-2*17+24-35,2*16-25-34,-[sqrt(3)]*25+[sqrt(3)]*34") {
		cout<<"CP^2 as the symmetric space SU(3)/U(2)"<<endl;
		lst hforms,vforms;
		hforms=e(2),e(3),e(4),e(5);
		vforms=e(1),e(6),e(7),e(8);
		TransverseRiemannianStructure g(this,hforms.begin(),hforms.end(),vforms.begin(),vforms.end());
		TransverseLeviCivitaConnection p(this,g);
		cout<<"Ricci tensor"<<endl;
		cout<<p.BaseRicciAsMatrix(true);
		cout<<p.BaseRicciAsMatrix(false);
		cout<<"Curvature form"<<endl;
		cout<<p.BaseCurvatureForm();		
	}
};


int main()
{
	cout<<latex;	//latex output
	CP2 cp2;
	S3TimesS1 M;
	cout<<"Some calculations on S^3\\times S^1"<<endl;
	cout<<M.e(1)*(M.e(1)+M.e(2))<<endl;
	cout<<M.d(M.e(1)+M.e(4))<<endl;
	cout<<M.d(M.e(1)*M.e(2))<<endl;
	cout<<M.LieDerivative(M.e(1),M.e(2)*M.e(4))<<endl;

	SomeCalculations();
	HermitianConnection();
	S3 M1;
	ThreeManifold M2;
	cout<<"Levi-Civita connection form of the Heisenberg group:"<<endl;
	cout<<LeviCivitaOfLieGroup("0,12,13")<<endl;
	//TS2 X;
	return 0;
}

