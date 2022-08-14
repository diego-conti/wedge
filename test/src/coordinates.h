/*******************************************************************************
 *  Copyright (C) 2007-2022 by Diego Conti, diego.conti@unimib.it 
 *  This file is part of Wedge.                                           
 *  Wedge is free software; you can redistribute it and/or modify         
 *  it under the terms of the GNU General Public License as published by  
 *  the Free Software Foundation; either version 3 of the License, or     
 *  (at your option) any later version.                                   
 *                                                                          
 *  Wedge is distributed in the hope that it will be useful,              
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of        
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         
 *  GNU General Public License for more details.                          
 *                                                                           
 *  You should have received a copy of the GNU General Public License     
 *  along with Wedge; if not, write to the                                
 *   Free Software Foundation, Inc.,                                       
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             
 *  
 *******************************************************************************/


#ifndef SRC_COORDINATES_H_
#define SRC_COORDINATES_H_

#include <cxxtest/TestSuite.h>
#include "test.h"

#include "wedge/base/normalform.h"

#include "wedge/connections/connection.h"
#include "wedge/liealgebras/liegroup.h"
#include "wedge/liealgebras/liegroupextension.h"
#include "wedge/liealgebras/liesubgroup.h"

#include "wedge/structures/riemannianstructure.h"
#include "wedge/structures/gstructure.h"
#include "wedge/structures/structures.h"
#include "wedge/manifolds/manifoldwith.h"

#include "wedge/manifolds/function.h"
#include "wedge/manifolds/coordinates.h"
#include "wedge/manifolds/liederivative.h"


using namespace GiNaC;
using namespace Wedge;
/**
	@author Diego Conti <diego.conti@unimib.it>
*/




//test function.h, coordinates.h, liederivative.h connection.h
class FunctionTestSuite : public CxxTest::TestSuite  {
	class TwoSphere : public ManifoldWithCoordinates {
	public:
		ExVector OrthFrame()
		{
			ExVector e(2);
			e(1)=2/(pow(x(1),2)+pow(x(2),2)+1)*d(x(1));
			e(2)=2/(pow(x(1),2)+pow(x(2),2)+1)*d(x(2));
			LOG_INFO(e);
			Frame x(e.begin(),e.end());
			LOG_INFO(x.dual());
			return e;
		}

		RiemannianStructure g;
		LeviCivitaConnection<true> omega;

		TwoSphere() : ManifoldWithCoordinates(2), g(this,OrthFrame()), omega(this,g)
		{
		}
	};

	class S3ByS2 : public ManifoldWithCoordinates {
		exvector CreateFrame()
		{
			exvector e(3);
			for (int i=0;i<3;++i) e[i]=DifferentialOneForm(N.e(i+1));
			return e;
		}
		exvector OrthFrame()	//riemannian product of S^2 by S^3
		{
			ExVector g=e();
			g(4)=2/(pow(x(1),2)+pow(x(2),2)+1)*d(x(1));
			g(5)=2/(pow(x(1),2)+pow(x(2),2)+1)*d(x(2));
			return g;
		}
		RiemannianStructure g;
	public:
		S3ByS2() : ManifoldWithCoordinates(CreateFrame(),2), g(this,OrthFrame())
		{
			Declare_d(e(1),e(2)*e(3));
			Declare_d(e(2),e(3)*e(1));
			Declare_d(e(3),e(1)*e(2));

			test();
		}
		void test()
		{
			Check_ddZero();
			TS_ASSERT_EQUALS(Dimension(),5);
			TS_ASSERT_EQUALS(d(x(1)),e(4));
			TS_ASSERT_EQUALS(d(x(2)),e(5));
			TS_ASSERT_EQUALS(d(sin(x(1))),cos(x(1))*e(4));
			TS_ASSERT_EQUALS(d(1/(x(1))),-1/pow(x(1),2)*e(4));
		}
	};
public:
//test liederivative.h
	void testLieBracket() {
		LieDerivative a;
		ManifoldWithCoordinates M(3);
		ex dx1=M.d(M.x(1)), dx2=M.d(M.x(2)), dx3=M.d(M.x(3));
		ex f=Function(N.f);
		ex Xf=M.LieDerivative(dx1,f);
		ex Yf=M.LieDerivative(dx2,f);
		TS_ASSERT_DIFFERS(Xf,Yf);
		TS_ASSERT_DIFFERS(Xf,f);
		TS_ASSERT_EQUALS(M.LieDerivative(dx2,Xf),M.LieDerivative(dx1,Yf));
		TS_ASSERT_EQUALS(M.LieDerivative(dx1,M.x(1)),1);
		TS_ASSERT_EQUALS(M.LieDerivative(dx1,M.x(1)*M.x(2)),M.x(2));
		TS_ASSERT_EQUALS(M.LieDerivative(dx1,log(M.x(1))*M.x(2)),M.x(2)/M.x(1));

		S3ByS2 X;
		TS_ASSERT_EQUALS(X.LieBracket(X.e(1),X.e(2)),-X.e(3));
		TS_ASSERT_EQUALS(X.LieBracket(ex_to<VectorField>(X.e(1)),ex_to<VectorField>(X.e(2))),-X.e(3));
		TS_ASSERT_EQUALS(X.LieBracket(X.e(1),X.e(2)),-X.e(3));
		VectorField e1 {ex_to<VectorField>(X.e(1))}, e2 {ex_to<VectorField>(X.e(2))};
		TS_ASSERT_EQUALS(X.LieBracket(e1,e2),-X.e(3));
		TS_ASSERT_EQUALS(X.LieBracket(X.x(1)*X.e(1),X.e(2)),-X.x(1)*X.e(3));

		Xf=X.LieDerivative(X.e(1),f);
		Yf=X.LieDerivative(X.e(4),f);
		TS_ASSERT_EQUALS(Xf+Yf,X.LieDerivative(X.e(1)+X.e(4),f));
		TS_ASSERT_EQUALS(sqrt(3)*Xf,X.LieDerivative(sqrt(3)*X.e(1),f));
		TS_ASSERT_EQUALS(X.x(1)*Xf,X.LieDerivative(X.x(1)*X.e(1),f));
		TS_ASSERT_DIFFERS(Xf,Yf);
		TS_ASSERT_DIFFERS(Xf,f);
		TS_ASSERT_EQUALS(X.LieDerivative(X.e(4),Xf),X.LieDerivative(X.e(1),Yf));
		Xf=X.LieDerivative(X.e(1),f);
		Yf=X.LieDerivative(X.e(2),f);
		TS_ASSERT_DIFFERS(Xf,Yf);
		TS_ASSERT_DIFFERS(Xf,f);
		TS_ASSERT_EQUALS(X.LieDerivative(X.e(2),Xf)-X.LieDerivative(X.e(1),Yf),X.LieDerivative(X.e(3),f));


		f=pow(X.x(1),3);
		Xf=X.LieDerivative(X.e(1),f);
		Yf=X.LieDerivative(X.e(4),f);
		TS_ASSERT_DIFFERS(Xf,Yf);
		TS_ASSERT_DIFFERS(Xf,f);
		TS_ASSERT_EQUALS(Xf,0);
		TS_ASSERT_EQUALS(Yf,diff(f,ex_to<symbol>(X.x(1))));
		TS_ASSERT_EQUALS(X.LieDerivative(X.e(4),Xf),X.LieDerivative(X.e(1),Yf));
		TS_ASSERT_EQUALS(X.LieDerivative(X.e(1),X.e(2)*X.e(3)),X.LieDerivative(X.e(1),X.e(2))*X.e(3)+X.e(2)*X.LieDerivative(X.e(1),X.e(3)));
		TS_ASSERT_EQUALS(X.LieDerivative(X.e(1),X.e(2)*X.e(4)),Hook(X.e(1),X.d(X.e(2)*X.e(4)))+X.d(Hook(X.e(1),X.e(2)*X.e(4))));
	}

	void testCoordinates() {
		ManifoldWithCoordinates M(3);
		ex dx1=M.d(M.x(1)), dx2=M.d(M.x(2)), dx3=M.d(M.x(3));
		TS_ASSERT_EQUALS(dx1,M.e(1));
		TS_ASSERT_EQUALS(dx2,M.e(2));
		TS_ASSERT_EQUALS(dx3,M.e(3));
		TS_ASSERT_EQUALS(M.Dimension(),3);

		TS_ASSERT_EQUALS(M.d(dx1),0);
		TS_ASSERT_EQUALS(M.d(M.x(1)*dx2),dx1*dx2);
		TS_ASSERT_EQUALS(M.d(M.x(1)*dx1),0);
		TS_ASSERT_EQUALS(M.d(M.x(1)*M.x(2)),M.x(1)*dx2+M.x(2)*dx1);
		TS_ASSERT_EQUALS(M.d(M.x(1)*M.x(1)),2*M.x(1)*dx1);
		TS_ASSERT_EQUALS(M.d(sin(M.x(1))),cos(M.x(1))*dx1);
		TS_ASSERT_EQUALS(M.d(1/(M.x(1))),-1/pow(M.x(1),2)*dx1);

		S3ByS2 N;
	}
	//(this also tests connection.h)
	void testS2() {
		TwoSphere S2;
		TS_ASSERT_EQUALS(S2.d(S2.x(1)*S2.e(2)),S2.e(1)*S2.e(2));

		TS_ASSERT_EQUALS(S2.LieBracket(S2.g.e().dual()(1),S2.g.e().dual()(2)).expand(),(S2.x(1)*S2.g.e().dual()(2)-S2.x(2)*S2.g.e().dual()(1)).expand());
		TS_ASSERT_EQUALS(S2.g.ScalarProduct<VectorField>(S2.x(1)*S2.g.e(2)-S2.x(2)*S2.g.e(1),S2.g.e(1)),-16*S2.x(2)*power(1+S2.x(1)*S2.x(1)+S2.x(2)*S2.x(2),-4));

		LOG_INFO(S2.omega);
		matrix R=S2.omega.CurvatureForm();
		LOG_INFO(R);
		TS_ASSERT_EQUALS(R(0,0).normal(),0);
		TS_ASSERT_EQUALS(R(1,1).normal(),0);
		TS_ASSERT_EQUALS(R(0,1).normal(),S2.g.e(1)*S2.g.e(2));
		TS_ASSERT_EQUALS(R(1,0).normal(),-S2.g.e(1)*S2.g.e(2));

		ExVector e=S2.OrthFrame();
		Connection k(&S2,e);
		for (int i=1;i<=2;++i)
		for (int j=1;j<=2;++j)
			k.DeclareNabla<DifferentialForm>(S2.e(i),S2.e(j),S2.omega.Nabla<DifferentialForm>(S2.e(i),S2.e(j)));
		TS_ASSERT_EQUALS(ex(R).normal(),ex(k.CurvatureForm()).normal());
	}
};



#endif /* SRC_COORDINATES_H_ */
