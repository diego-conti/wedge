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
#ifndef CONNECTIONS_H_
#define CONNECTIONS_H_

#include <cxxtest/TestSuite.h>
#include "test.h"
#include "wedge/structures/riemannianstructure.h"
#include "wedge/structures/transversestructure.h"
#include "wedge/structures/structures.h"
#include "wedge/structures/pseudoriemannianstructure.h"
#include "wedge/connections/pseudolevicivita.h"
#include "wedge/connections/transverseconnection.h"
#include "wedge/liealgebras/liegroup.h"
#include "wedge/liealgebras/liesubgroup.h"
#include "wedge/structures/gstructure.h"
#include "wedge/manifolds/manifoldwith.h"
#include "wedge/manifolds/coordinates.h"
#include "wedge/structures/submersion.h"

using namespace GiNaC;
using namespace Wedge;

//test riemannianstructure.h
class RiemannianTestSuite : public CxxTest::TestSuite 
{
public:
	class RiemannianManifold: public ConcreteManifold
	{
	public:
		RiemannianStructure structure;
		RiemannianManifold() : ConcreteManifold(5), structure(this,ExVector(e())) {}
		void MetricTest() {			
			TS_ASSERT_EQUALS(structure.DimensionOfSpinorRepresentation(),4);
			TS_ASSERT_EQUALS(structure.ScalarProduct<DifferentialForm>(structure.e(1),structure.e(1)),1);
						
			TS_ASSERT_EQUALS(structure.ScalarProduct<DifferentialForm>(
					structure.e(1)+3*structure.e(2)+5*structure.e(3),
					structure.e(2)-2*structure.e(3)-4*structure.e(4)
				), 3-10);
			TS_ASSERT_EQUALS(structure.ScalarProduct<DifferentialForm>(
					structure.e(1)*structure.e(2),
					structure.e(2)*structure.e(1)
				), -1);
				
			TS_ASSERT_EQUALS(structure.SquareNorm<DifferentialForm>(structure.e(1)*structure.e(2)),1);
			TS_ASSERT_EQUALS(structure.SquareNorm<DifferentialForm>(ParseDifferentialForm(structure.e(),"123+145")),2);		
				
			TS_ASSERT_EQUALS(structure.ScalarProduct<Spinor>(
					structure.u(1)+3*structure.u(2)+I*5*structure.u(3),
					I*structure.u(2)+2*I*structure.u(3)-4*structure.u(4)
				), 10);			

			TS_ASSERT_EQUALS(structure.HodgeStar(structure.e(1)*structure.e(4)*structure.e(5)),structure.e(2)*structure.e(3));
			TS_ASSERT_EQUALS(structure.HodgeStar(structure.e(3)*structure.e(5)),-structure.e(1)*structure.e(2)*structure.e(4));
			TS_ASSERT_EQUALS(structure.ScalarProduct<Spinor>(structure.u(0),structure.u(1)),0);
			TS_ASSERT_EQUALS(structure.ScalarProduct<Spinor>(structure.u(0)-structure.u(1),structure.u(1)),-1);
			TS_ASSERT_EQUALS(structure.SquareNorm<Spinor>(structure.u(0)+I*structure.u(1)),2);
			TS_ASSERT_EQUALS(structure.SquareNorm<DifferentialForm>(structure.e(1)-2*structure.e(2)),5);
			TS_ASSERT_EQUALS(structure.Hook(structure.e(1)+structure.e(2),structure.e(1)*structure.e(2)),structure.e(2)-structure.e(1));

			exvector dual=structure.e().dual();	//vector fields
			for (exvector::const_iterator i=dual.begin();i!=dual.end();++i)
			for (exvector::const_iterator j=dual.begin();j!=dual.end();++j)
				if (i==j) {
					TS_ASSERT_EQUALS(structure.ScalarProduct<VectorField>(*i,*j),1);
				}
				else {
					TS_ASSERT_EQUALS(structure.ScalarProduct<VectorField>(*i,*j),0);
				}
		}

		void CliffordTest() {
			TS_ASSERT_EQUALS(structure.CliffordDot(structure.e(1),
								structure.CliffordDot(structure.e(2),structure.u(0))),
								structure.CliffordDot(structure.e(1)*structure.e(2),structure.u(0)));
			
			for (int i=1;i<=Dimension();i++)
				for (int j=1;j<=Dimension();j++)
					for (int k=0;k<structure.DimensionOfSpinorRepresentation();k++) {
				LOG_DEBUG(e(i)<<"*"<<structure.u(k)<<"="<<structure.CliffordDot(e(i),structure.u(k)));
				ex ijk=structure.CliffordDot(structure.e(i),structure.CliffordDot(e(j),I*structure.u(k)));
				ex ijk2=structure.CliffordDot(structure.e(i)*structure.e(j),I*structure.u(k));
				if (i!=j) TS_ASSERT_EQUALS(ijk,ijk2);
				ex jik=structure.CliffordDot(structure.e(j),structure.CliffordDot(structure.e(i),I*structure.u(k)));
				TS_ASSERT_EQUALS(ijk+jik,i==j?-2*I*structure.u(k):0);
					}
		}
	};		
	 
	void testRiemannian( void )
	{
		RiemannianManifold M;
		M.MetricTest();
		M.CliffordTest();
		ExVector e(5);
		e(1)=M.e(1);
		e(2)=M.e(2)+e(1);
		e(3)=M.e(3)+e(2);
		e(4)=M.e(4)+e(3);
		e(5)=M.e(5)+e(4);
		RiemannianStructure g(&M,e);
		M.structure=g;
		TS_ASSERT_EQUALS(M.structure.e(1),e(1));
		TS_ASSERT_EQUALS(M.structure.e(5),e(5));
		M.MetricTest();
	      //M.CliffordTest();	//fails due to a bug in CliffordDot (see documentation of CliffordDot)
	}
};


class ConnectionTestSuite : public CxxTest::TestSuite  {
public:
	class S3 : public ConcreteManifold, public virtual Has_dTable {
	public:
		S3() : ConcreteManifold(3) {				
			Declare_d(e(1),e(2)*e(3));
			Declare_d(e(2),e(3)*e(1));
			Declare_d(e(3),e(1)*e(2));
		}
	};
	
	class ThreeManifold : public virtual Has_dTable {
		Frame frame; 
		Function x,y,z;
	public:
		ThreeManifold() {
			for (int i=1;i<=3;++i)
				frame.push_back(DifferentialOneForm(N.e(1)));
		}
		const Frame& e() const {return frame;}
		ex e(OneBased i) const {return frame(i);}
	};
	
	//test connection.h
	void testConnection()
	{
		{
			ThreeManifold M;
			exvector v;
			v.push_back(M.e(1)+M.e(2));
			v.push_back(M.e(2));
			v.push_back(M.e(3)+2*M.e(1));
			Frame e(v);
	
			Connection connection(&M,e);
			connection.DeclareNabla<DifferentialForm>(e.dual()(1),e(1),e(2));
			connection.DeclareNabla<DifferentialForm>(e.dual()(1),e(2),e(2));
			connection.DeclareNabla<DifferentialForm>(e.dual()(1),e(3),e(3));
			TS_ASSERT_EQUALS(connection.Nabla<DifferentialForm>(e.dual()(1),e(1)*e(2)),e(1)*e(2));
			TS_ASSERT_EQUALS(connection.Nabla<DifferentialForm>(e.dual()(1),e(1)*e(3)),e(2)*e(3)+e(1)*e(3));
			connection.DeclareNabla<DifferentialForm>(e.dual()(2),e(1),0);
			connection.DeclareNabla<DifferentialForm>(e.dual()(2),e(2),0);
			connection.DeclareNabla<DifferentialForm>(e.dual()(2),e(3),0);
			TS_ASSERT_EQUALS(connection.Nabla<DifferentialForm>(e.dual()(1)+e.dual()(2),e(1)*e(2)),e(1)*e(2));
			TS_ASSERT_EQUALS(connection.Nabla<DifferentialForm>(e.dual()(1)+e.dual()(2),e(1)*e(3)),e(2)*e(3)+e(1)*e(3));
			TS_ASSERT_THROWS(connection.CurvatureForm(),WedgeException<std::logic_error>);
			TS_ASSERT_THROWS(connection.Torsion(),WedgeException<std::logic_error>);

			//check that the contraction commutes with the covariant derivative	
			ex X=symbol("x")*e.dual()(1)+symbol("y")*e.dual()(2)+symbol("z")*e.dual()(3);
			for (int i=1;i<=3;++i)
			for (int j=1;j<=3;++j)
				 TS_ASSERT_EQUALS(Hook(e.dual()(j),connection.Nabla<DifferentialForm>(X,e(i))),-Hook(connection.Nabla<VectorField>(X,e.dual()(j)),e(i)));				
		//verify the Leibnitz rule
			Function f(N.f),g(N.g);
			ex F=f;	
			TS_ASSERT_EQUALS(connection.Nabla<DifferentialForm>(M.e(1),F),M.LieDerivative(M.e(1),F));
			TS_ASSERT_EQUALS(connection.Nabla<DifferentialForm>(M.e(1),F),connection.Nabla<VectorField>(M.e(1),F));
			TS_ASSERT_EQUALS(connection.Nabla<VectorField>(e(1),F),M.LieDerivative(e(1),F));
			F=pow(f,2)+sin(g);
			TS_ASSERT_EQUALS(connection.Nabla<DifferentialForm>(e(1),F),M.LieDerivative(e(1),F));
			TS_ASSERT_EQUALS(connection.Nabla<VectorField>(e(1),F),M.LieDerivative(e(1),F));
			F=pow(f,2)*log(g);
			TS_ASSERT_EQUALS(connection.Nabla<DifferentialForm>(e(1),F),M.LieDerivative(e(1),F));
			TS_ASSERT_EQUALS(connection.Nabla<VectorField>(e(1),F),M.LieDerivative(e(1),F));

			TS_ASSERT_EQUALS(connection.Nabla<DifferentialForm>(X,F),M.LieDerivative(X,F));
			TS_ASSERT_EQUALS(connection.Nabla<VectorField>(X,F),M.LieDerivative(X,F));
			F=pow(f,2)+sin(g);
			TS_ASSERT_EQUALS(connection.Nabla<DifferentialForm>(X,F),M.LieDerivative(X,F));
			TS_ASSERT_EQUALS(connection.Nabla<VectorField>(X,F),M.LieDerivative(X,F));
			F=pow(f,2)*log(g);
			TS_ASSERT_EQUALS(connection.Nabla<DifferentialForm>(X,F),M.LieDerivative(X,F));
			TS_ASSERT_EQUALS(connection.Nabla<VectorField>(X,F),M.LieDerivative(X,F));

			for (int i=1;i<=3;++i)
			{
				ex F=f;		
				ex nonlinearterm=(connection.Nabla<DifferentialForm>(X,F)*e(i)).expand();
				TS_ASSERT_DIFFERS(connection.Nabla<DifferentialForm>(X,F*e(i))-(F*connection.Nabla<DifferentialForm>(X,e(i))).expand(),0);
				TS_ASSERT_EQUALS(connection.Nabla<DifferentialForm>(X,F*e(i))-(F*connection.Nabla<DifferentialForm>(X,e(i))).expand(),nonlinearterm);
				F=f*g;
				nonlinearterm=(connection.Nabla<DifferentialForm>(X,F)*e(i)).expand();
				TS_ASSERT_DIFFERS(connection.Nabla<DifferentialForm>(X,F*e(i))-(F*connection.Nabla<DifferentialForm>(X,e(i))).expand(),0);
				TS_ASSERT_EQUALS(connection.Nabla<DifferentialForm>(X,F*e(i))-(F*connection.Nabla<DifferentialForm>(X,e(i))).expand(),nonlinearterm);
				F=pow(f,2)+sin(g);
				nonlinearterm=(connection.Nabla<DifferentialForm>(X,F)*e(i)).expand();
				TS_ASSERT_DIFFERS(connection.Nabla<DifferentialForm>(X,F*e(i))-(F*connection.Nabla<DifferentialForm>(X,e(i))).expand(),0);
				TS_ASSERT_EQUALS(connection.Nabla<DifferentialForm>(X,F*e(i))-(F*connection.Nabla<DifferentialForm>(X,e(i))).expand(),nonlinearterm);
				F=pow(f,2)*log(g);
				nonlinearterm=(connection.Nabla<DifferentialForm>(X,F)*e(i)).expand();
				TS_ASSERT_DIFFERS(connection.Nabla<DifferentialForm>(X,F*e(i))-(F*connection.Nabla<DifferentialForm>(X,e(i))).expand(),0);
				TS_ASSERT_EQUALS(connection.Nabla<DifferentialForm>(X,F*e(i))-(F*connection.Nabla<DifferentialForm>(X,e(i))).expand(),nonlinearterm);
			}
		}
	
		{
			S3 M;
			exvector v;
			v.push_back(M.e(1)+M.e(2));
			v.push_back(M.e(2));
			v.push_back(M.e(3)+2*M.e(1));
			Frame e(v);
	
			Connection connection(&M,e);
			connection.DeclareNabla<DifferentialForm>(e.dual()(1),e(1),e(2));
			connection.DeclareNabla<DifferentialForm>(e.dual()(1),e(2),e(2));
			connection.DeclareNabla<DifferentialForm>(e.dual()(1),e(3),e(3));
			TS_ASSERT_EQUALS(connection.Nabla<DifferentialForm>(e.dual()(1),e(1)*e(2)),e(1)*e(2));
			TS_ASSERT_EQUALS(connection.Nabla<DifferentialForm>(e.dual()(1),e(1)*e(3)),e(2)*e(3)+e(1)*e(3));
			connection.DeclareNabla<DifferentialForm>(e.dual()(2),e(1),0);
			connection.DeclareNabla<DifferentialForm>(e.dual()(2),e(2),0);
			connection.DeclareNabla<DifferentialForm>(e.dual()(2),e(3),0);
			TS_ASSERT_EQUALS(connection.Nabla<DifferentialForm>(e.dual()(1)+e.dual()(2),e(1)*e(2)),e(1)*e(2));
			TS_ASSERT_EQUALS(connection.Nabla<DifferentialForm>(e.dual()(1)+e.dual()(2),e(1)*e(3)),e(2)*e(3)+e(1)*e(3));
			connection.DeclareNabla<DifferentialForm>(e.dual()(3),e(1),e(1));
			connection.DeclareNabla<DifferentialForm>(e.dual()(3),e(2),0);
			connection.DeclareNabla<DifferentialForm>(e.dual()(3),e(3),0);
		//check that the covariant derivative works as expected on vector fields
			TS_ASSERT_EQUALS(connection.Nabla<VectorField>(e.dual()(1),e.dual()(1)),0);
			TS_ASSERT_EQUALS(connection.Nabla<VectorField>(e.dual()(2),e.dual()(1)),0);
			TS_ASSERT_EQUALS(connection.Nabla<VectorField>(e.dual()(3),e.dual()(1)),-e.dual()(1));
			TS_ASSERT_EQUALS(connection.Nabla<VectorField>(e.dual()(1),e.dual()(2)),-e.dual()(1)-e.dual()(2));
			TS_ASSERT_EQUALS(connection.Nabla<VectorField>(e.dual()(2),e.dual()(2)),0);
			TS_ASSERT_EQUALS(connection.Nabla<VectorField>(e.dual()(3),e.dual()(2)),0);
			TS_ASSERT_EQUALS(connection.Nabla<VectorField>(e.dual()(1),e.dual()(3)),-e.dual()(3));
			TS_ASSERT_EQUALS(connection.Nabla<VectorField>(e.dual()(2),e.dual()(3)),0);
			TS_ASSERT_EQUALS(connection.Nabla<VectorField>(e.dual()(3),e.dual()(3)),0);
		//check that the connection form components match with the above
			matrix connectionform(3,3);
			connectionform={{-e(3),-e(1),0},
				{0,-e(1),0},
				{0,0,-e(1)}};
			TS_ASSERT_EQUALS(connection.AsMatrix(),connectionform);

			matrix curvature(3,3);
			curvature={{-M.d(e(3)),-M.d(e(1))+e(3)*e(1),0},
				{0,-M.d(e(1)),0},
				{0,0,-M.d(e(1))}};
			TS_ASSERT_EQUALS(connection.CurvatureForm(),curvature);
			exvector torsion;
			torsion.push_back(M.d(e(1))-e(3)*e(1)-e(1)*e(2));
			torsion.push_back(M.d(e(2))-e(1)*e(2));
			torsion.push_back(M.d(e(3))-e(1)*e(3));
			TS_ASSERT_EQUALS(connection.Torsion(),torsion);	

		//verify the Leibnitz rule
			ex X=symbol("x")*e.dual()(1)+symbol("y")*e.dual()(2)+symbol("z")*e.dual()(3);
			Function f(N.f),g(N.g);
			ex F=f;	
			TS_ASSERT_EQUALS(connection.Nabla<DifferentialForm>(X,F),M.LieDerivative(X,F));
			TS_ASSERT_EQUALS(connection.Nabla<VectorField>(X,F),M.LieDerivative(X,F));
			F=pow(f,2)+sin(g);
			TS_ASSERT_EQUALS(connection.Nabla<DifferentialForm>(X,F),M.LieDerivative(X,F));
			TS_ASSERT_EQUALS(connection.Nabla<VectorField>(X,F),M.LieDerivative(X,F));
			F=pow(f,2)*log(g);
			TS_ASSERT_EQUALS(connection.Nabla<DifferentialForm>(X,F),M.LieDerivative(X,F));
			TS_ASSERT_EQUALS(connection.Nabla<VectorField>(X,F),M.LieDerivative(X,F));
			for (int i=1;i<=3;++i)
			{
				ex F=f;		
				ex nonlinearterm=(connection.Nabla<DifferentialForm>(X,F)*e(i)).expand();
				TS_ASSERT_DIFFERS(connection.Nabla<DifferentialForm>(X,F*e(i))-(F*connection.Nabla<DifferentialForm>(X,e(i))).expand(),0);
				TS_ASSERT_EQUALS(connection.Nabla<DifferentialForm>(X,F*e(i))-(F*connection.Nabla<DifferentialForm>(X,e(i))).expand(),nonlinearterm);
				F=f*g;
				nonlinearterm=(connection.Nabla<DifferentialForm>(X,F)*e(i)).expand();
				TS_ASSERT_DIFFERS(connection.Nabla<DifferentialForm>(X,F*e(i))-(F*connection.Nabla<DifferentialForm>(X,e(i))).expand(),0);
				TS_ASSERT_EQUALS(connection.Nabla<DifferentialForm>(X,F*e(i))-(F*connection.Nabla<DifferentialForm>(X,e(i))).expand(),nonlinearterm);
				F=pow(f,2)+sin(g);
				nonlinearterm=(connection.Nabla<DifferentialForm>(X,F)*e(i)).expand();
				TS_ASSERT_DIFFERS(connection.Nabla<DifferentialForm>(X,F*e(i))-(F*connection.Nabla<DifferentialForm>(X,e(i))).expand(),0);
				TS_ASSERT_EQUALS(connection.Nabla<DifferentialForm>(X,F*e(i))-(F*connection.Nabla<DifferentialForm>(X,e(i))).expand(),nonlinearterm);
				F=pow(f,2)*log(g);
				nonlinearterm=(connection.Nabla<DifferentialForm>(X,F)*e(i)).expand();
				TS_ASSERT_DIFFERS(connection.Nabla<DifferentialForm>(X,F*e(i))-(F*connection.Nabla<DifferentialForm>(X,e(i))).expand(),0);
				TS_ASSERT_EQUALS(connection.Nabla<DifferentialForm>(X,F*e(i))-(F*connection.Nabla<DifferentialForm>(X,e(i))).expand(),nonlinearterm);
			}
		}	
	}

	//test riemannianconnection.h
	void testRiemannianConnection()
	{
		{
			ThreeManifold M;
	
			exvector v;
			v.push_back(M.e(1)+M.e(2));
			v.push_back(M.e(2));
			v.push_back(M.e(3)+2*M.e(1));
			Frame e(v);
	
			RiemannianConnection connection(&M,RiemannianStructure(&M,e));
			connection.DeclareNabla<DifferentialForm>(e.dual()(1),e(1),e(2));
			TS_ASSERT_THROWS(connection.DeclareNabla<DifferentialForm>(e.dual()(1),e(2),e(2)),InconsistentDeclaration);
			connection.DeclareNabla<DifferentialForm>(e.dual()(1),e(2),-e(1));
			TS_ASSERT_EQUALS(connection.Nabla<DifferentialForm>(e.dual()(1),e(3)),0);	//because it's riemannian		
			TS_ASSERT_EQUALS(connection.Nabla<DifferentialForm>(e.dual()(1),e(1)*e(2)),0);
			TS_ASSERT_EQUALS(connection.Nabla<DifferentialForm>(e.dual()(1),e(1)*e(3)),e(2)*e(3));		
			connection.DeclareNabla<DifferentialForm>(e.dual()(2),e(1),0);
			connection.DeclareNabla<DifferentialForm>(e.dual()(2),e(2),0);
			connection.DeclareNabla<DifferentialForm>(e.dual()(3),e(1),0);
			connection.DeclareNabla<DifferentialForm>(e.dual()(3),e(2),0);
			
			TS_ASSERT_EQUALS(connection.Nabla<DifferentialForm>(e.dual()(1)+e.dual()(2),e(3)),0);	//because it's riemannian		
			TS_ASSERT_EQUALS(connection.Nabla<DifferentialForm>(e.dual()(1)+e.dual()(2),e(1)*e(2)),0);
			TS_ASSERT_EQUALS(connection.Nabla<DifferentialForm>(e.dual()(1)+e.dual()(2),e(1)*e(3)),e(2)*e(3));
			TS_ASSERT_THROWS(connection.CurvatureForm(),WedgeException<std::logic_error>);
			TS_ASSERT_THROWS(connection.Torsion(),WedgeException<std::logic_error>);

			//check that the contraction commutes with the covariant derivative
			ex X=symbol("x")*e.dual()(1)+symbol("y")*e.dual()(2)+symbol("z")*e.dual()(3);
			for (int i=1;i<=3;++i)
			for (int j=1;j<=3;++j)
				TS_ASSERT_EQUALS(Hook(e.dual()(j),connection.Nabla<DifferentialForm>(X,e(i))),-Hook(connection.Nabla<VectorField>(X,e.dual()(j)),e(i)));				

		//verify the Leibnitz rule
			Function f(N.f),g(N.g);
			ex F=f;	
			TS_ASSERT_EQUALS(connection.Nabla<DifferentialForm>(X,F),M.LieDerivative(X,F));
			TS_ASSERT_EQUALS(connection.Nabla<VectorField>(X,F),M.LieDerivative(X,F));
			F=pow(f,2)+sin(g);
			TS_ASSERT_EQUALS(connection.Nabla<DifferentialForm>(X,F),M.LieDerivative(X,F));
			TS_ASSERT_EQUALS(connection.Nabla<VectorField>(X,F),M.LieDerivative(X,F));
			F=pow(f,2)*log(g);
			TS_ASSERT_EQUALS(connection.Nabla<DifferentialForm>(X,F),M.LieDerivative(X,F));
			TS_ASSERT_EQUALS(connection.Nabla<VectorField>(X,F),M.LieDerivative(X,F));
			for (int i=1;i<=3;++i)
			{
				ex F=f;		
				ex nonlinearterm=(connection.Nabla<DifferentialForm>(X,F)*e(i)).expand();
				TS_ASSERT_DIFFERS(connection.Nabla<DifferentialForm>(X,F*e(i))-(F*connection.Nabla<DifferentialForm>(X,e(i))).expand(),0);
				TS_ASSERT_EQUALS(connection.Nabla<DifferentialForm>(X,F*e(i))-(F*connection.Nabla<DifferentialForm>(X,e(i))).expand(),nonlinearterm);
				F=f*g;
				nonlinearterm=(connection.Nabla<DifferentialForm>(X,F)*e(i)).expand();
				TS_ASSERT_DIFFERS(connection.Nabla<DifferentialForm>(X,F*e(i))-(F*connection.Nabla<DifferentialForm>(X,e(i))).expand(),0);
				TS_ASSERT_EQUALS(connection.Nabla<DifferentialForm>(X,F*e(i))-(F*connection.Nabla<DifferentialForm>(X,e(i))).expand(),nonlinearterm);
				F=pow(f,2)+sin(g);
				nonlinearterm=(connection.Nabla<DifferentialForm>(X,F)*e(i)).expand();
				TS_ASSERT_DIFFERS(connection.Nabla<DifferentialForm>(X,F*e(i))-(F*connection.Nabla<DifferentialForm>(X,e(i))).expand(),0);
				TS_ASSERT_EQUALS(connection.Nabla<DifferentialForm>(X,F*e(i))-(F*connection.Nabla<DifferentialForm>(X,e(i))).expand(),nonlinearterm);
				F=pow(f,2)*log(g);
				nonlinearterm=(connection.Nabla<DifferentialForm>(X,F)*e(i)).expand();
				TS_ASSERT_DIFFERS(connection.Nabla<DifferentialForm>(X,F*e(i))-(F*connection.Nabla<DifferentialForm>(X,e(i))).expand(),0);
				TS_ASSERT_EQUALS(connection.Nabla<DifferentialForm>(X,F*e(i))-(F*connection.Nabla<DifferentialForm>(X,e(i))).expand(),nonlinearterm);
			}
		}
		{
			S3 M;
	
			exvector v;
			v.push_back(M.e(1)+M.e(2));
			v.push_back(M.e(2));
			v.push_back(M.e(3)+2*M.e(1));
			Frame e(v);
	
			RiemannianConnection connection(&M,RiemannianStructure(&M,e));
			connection.DeclareNabla<DifferentialForm>(e.dual()(1),e(1),e(2));
			TS_ASSERT_THROWS(connection.DeclareNabla<DifferentialForm>(e.dual()(1),e(2),e(2)),InconsistentDeclaration);
			connection.DeclareNabla<DifferentialForm>(e.dual()(1),e(2),-e(1));
			TS_ASSERT_EQUALS(connection.Nabla<DifferentialForm>(e.dual()(1),e(3)),0);	//because it's riemannian		
			TS_ASSERT_EQUALS(connection.Nabla<DifferentialForm>(e.dual()(1),e(1)*e(2)),0);
			TS_ASSERT_EQUALS(connection.Nabla<DifferentialForm>(e.dual()(1),e(1)*e(3)),e(2)*e(3));		
			connection.DeclareNabla<DifferentialForm>(e.dual()(2),e(1),0);
			connection.DeclareNabla<DifferentialForm>(e.dual()(2),e(2),0);
			connection.DeclareNabla<DifferentialForm>(e.dual()(3),e(1),0);
			connection.DeclareNabla<DifferentialForm>(e.dual()(3),e(2),0);
			
			TS_ASSERT_EQUALS(connection.Nabla<DifferentialForm>(e.dual()(1)+e.dual()(2),e(3)),0);		
			TS_ASSERT_EQUALS(connection.Nabla<DifferentialForm>(e.dual()(1)+e.dual()(2),e(1)*e(2)),0);
			TS_ASSERT_EQUALS(connection.Nabla<DifferentialForm>(e.dual()(1)+e.dual()(2),e(1)*e(3)),e(2)*e(3));
			
			matrix curvature(3,3);
			curvature={{0,-M.d(e(1)),0},
				{M.d(e(1)),0,0},
				{0,0,0}};
			TS_ASSERT_EQUALS(connection.CurvatureForm(),curvature);
			exvector torsion;
			torsion.push_back(M.d(e(1))-e(1)*e(2));
			torsion.push_back(M.d(e(2)));
			torsion.push_back(M.d(e(3)));
			TS_ASSERT_EQUALS(connection.Torsion(),torsion);
			
			//check that the contraction commutes with the covariant derivative
			ex X=symbol("x")*e.dual()(1)+symbol("y")*e.dual()(2)+symbol("z")*e.dual()(3);
			for (int i=1;i<=3;++i)
			for (int j=1;j<=3;++j)
				TS_ASSERT_EQUALS(Hook(e.dual()(j),connection.Nabla<DifferentialForm>(X,e(i))),-Hook(connection.Nabla<VectorField>(X,e.dual()(j)),e(i)));				

		//verify the Leibnitz rule
			Function f(N.f),g(N.g);
			ex F=f;	
			TS_ASSERT_EQUALS(connection.Nabla<DifferentialForm>(X,F),M.LieDerivative(X,F));
			TS_ASSERT_EQUALS(connection.Nabla<VectorField>(X,F),M.LieDerivative(X,F));
			F=pow(f,2)+sin(g);
			TS_ASSERT_EQUALS(connection.Nabla<DifferentialForm>(X,F),M.LieDerivative(X,F));
			TS_ASSERT_EQUALS(connection.Nabla<VectorField>(X,F),M.LieDerivative(X,F));
			F=pow(f,2)*log(g);
			TS_ASSERT_EQUALS(connection.Nabla<DifferentialForm>(X,F),M.LieDerivative(X,F));
			TS_ASSERT_EQUALS(connection.Nabla<VectorField>(X,F),M.LieDerivative(X,F));
			for (int i=1;i<=3;++i)
			{
				ex F=f;		
				ex nonlinearterm=(connection.Nabla<DifferentialForm>(X,F)*e(i)).expand();
				TS_ASSERT_DIFFERS(connection.Nabla<DifferentialForm>(X,F*e(i))-(F*connection.Nabla<DifferentialForm>(X,e(i))).expand(),0);
				TS_ASSERT_EQUALS(connection.Nabla<DifferentialForm>(X,F*e(i))-(F*connection.Nabla<DifferentialForm>(X,e(i))).expand(),nonlinearterm);
				F=f*g;
				nonlinearterm=(connection.Nabla<DifferentialForm>(X,F)*e(i)).expand();
				TS_ASSERT_DIFFERS(connection.Nabla<DifferentialForm>(X,F*e(i))-(F*connection.Nabla<DifferentialForm>(X,e(i))).expand(),0);
				TS_ASSERT_EQUALS(connection.Nabla<DifferentialForm>(X,F*e(i))-(F*connection.Nabla<DifferentialForm>(X,e(i))).expand(),nonlinearterm);
				F=pow(f,2)+sin(g);
				nonlinearterm=(connection.Nabla<DifferentialForm>(X,F)*e(i)).expand();
				TS_ASSERT_DIFFERS(connection.Nabla<DifferentialForm>(X,F*e(i))-(F*connection.Nabla<DifferentialForm>(X,e(i))).expand(),0);
				TS_ASSERT_EQUALS(connection.Nabla<DifferentialForm>(X,F*e(i))-(F*connection.Nabla<DifferentialForm>(X,e(i))).expand(),nonlinearterm);
				F=pow(f,2)*log(g);
				nonlinearterm=(connection.Nabla<DifferentialForm>(X,F)*e(i)).expand();
				TS_ASSERT_DIFFERS(connection.Nabla<DifferentialForm>(X,F*e(i))-(F*connection.Nabla<DifferentialForm>(X,e(i))).expand(),0);
				TS_ASSERT_EQUALS(connection.Nabla<DifferentialForm>(X,F*e(i))-(F*connection.Nabla<DifferentialForm>(X,e(i))).expand(),nonlinearterm);
			}

			ex Ric; matrix RicM(3,3);
			for (int i=1;i<=3;++i)
			for (int j=1;j<=3;++j)
			{
				ex X=e.dual()(i);
				ex Y=e.dual()(j);
				ex Ric_XY;
				for (int k=1;k<=3;++k)
				{
					Ric_XY+=Hook(connection.Nabla<VectorField>(e.dual()(k),connection.Nabla<VectorField>(X,Y)),e(k));
					Ric_XY-=Hook(connection.Nabla<VectorField>(X,connection.Nabla<VectorField>(e.dual()(k),Y)),e(k));
					Ric_XY-=Hook(connection.Nabla<VectorField>(M.LieBracket(e.dual()(k),X),Y),e(k));
				}
				Ric+=Ric_XY*TensorProduct<DifferentialOneForm,DifferentialOneForm>(e(i),e(j));
				RicM(i-1,j-1)=Ric_XY;
			}
			TS_ASSERT_EQUALS(Ric,connection.Ricci());
			TS_ASSERT_EQUALS(RicM,connection.RicciAsMatrix());
		}	
	}		

	//test torsionfreeconnection.h (no dtable)
	void testTorsionFreeConnection()
	{
		ConcreteManifold M(3);

		exvector v;
		v.push_back(M.e(1)+M.e(2));
		v.push_back(M.e(2));
		v.push_back(M.e(3)+2*M.e(1));
		Frame e(v);

		TorsionFreeConnection<false> connection(&M,e);
		connection.Declare_d(e(1),e(2)*e(3));
		LOG_INFO(connection.AsMatrix());
		TS_ASSERT_EQUALS(connection.d(e(1)).expand(),(e(2)*e(3)).expand());
		TS_ASSERT_EQUALS(connection.d(e(2)*e(3)).expand(),0);
		connection.Declare_d(e(2),e(3)*e(1));
		connection.Declare_d(e(3),e(1)*e(2));
		TS_ASSERT(connection.Nabla<DifferentialForm>(e.dual()(1),e(1))!=0);
		TS_ASSERT(connection.Nabla<DifferentialForm>(e.dual()(1),e(2))!=-e(3)/2);
		TS_ASSERT(connection.Nabla<DifferentialForm>(e.dual()(1),e(3))!=e(2)/2);

		//check that the contraction commutes with the covariant derivative
		ex X=symbol("x")*e.dual()(1)+symbol("y")*e.dual()(2)+symbol("z")*e.dual()(3);
		for (int i=1;i<=3;++i)
		for (int j=1;j<=3;++j)
			TS_ASSERT_EQUALS(Hook(e.dual()(j),connection.Nabla<DifferentialForm>(X,e(i))),-Hook(connection.Nabla<VectorField>(X,e.dual()(j)),e(i)));				
		
		LeviCivitaConnection<false> leviCivita(&M,RiemannianStructure(&M,e));
		leviCivita.Declare_d(e(1),e(2)*e(3));
		TS_ASSERT_EQUALS(leviCivita.d(e(1).expand()),(e(2)*e(3)).expand());
		TS_ASSERT_EQUALS(leviCivita.d(e(2)*e(3)),0);
		leviCivita.Declare_d(e(2),e(3)*e(1));
		leviCivita.Declare_d(e(3),e(1)*e(2));
		TS_ASSERT_EQUALS(leviCivita.Nabla<DifferentialForm>(e.dual()(1),e(1)),0);
		TS_ASSERT_EQUALS(leviCivita.Nabla<DifferentialForm>(e.dual()(1),e(2)),-e(3)/2);
		TS_ASSERT_EQUALS(leviCivita.Nabla<DifferentialForm>(e.dual()(1),e(3)),e(2)/2);
		
		matrix R=leviCivita.CurvatureForm();		
		TS_ASSERT_EQUALS(R(0,0),0);
		TS_ASSERT_EQUALS(R(1,1),0);
		TS_ASSERT_EQUALS(R(2,2),0);
		TS_ASSERT_EQUALS(R(0,1),-R(1,0));
		TS_ASSERT_EQUALS(R(0,2),-R(2,0));
		TS_ASSERT_EQUALS(R(1,2),-R(2,1));
		TS_ASSERT_EQUALS(R(0,1),(e(1)*e(2)/4).expand());
		TS_ASSERT_EQUALS(R(0,2),(e(1)*e(3)/4).expand());
		TS_ASSERT_EQUALS(R(1,2),(e(2)*e(3)/4).expand());
		
		ex ric=TensorProduct<DifferentialOneForm,DifferentialOneForm>(e[0],e[0])/2+
			TensorProduct<DifferentialOneForm,DifferentialOneForm>(e[1],e[1])/2+
			TensorProduct<DifferentialOneForm,DifferentialOneForm>(e[2],e[2])/2;
		TS_ASSERT_EQUALS(leviCivita.Ricci(),ric);

		//verify the Leibnitz rule
			Function f(N.f),g(N.g);
			ex F=f;	
			TS_ASSERT_EQUALS(connection.Nabla<DifferentialForm>(X,F),M.LieDerivative(X,F));
			TS_ASSERT_EQUALS(connection.Nabla<VectorField>(X,F),M.LieDerivative(X,F));
			F=pow(f,2)+sin(g);
			TS_ASSERT_EQUALS(connection.Nabla<DifferentialForm>(X,F),M.LieDerivative(X,F));
			TS_ASSERT_EQUALS(connection.Nabla<VectorField>(X,F),M.LieDerivative(X,F));
			F=pow(f,2)*log(g);
			TS_ASSERT_EQUALS(connection.Nabla<DifferentialForm>(X,F),M.LieDerivative(X,F));
			TS_ASSERT_EQUALS(connection.Nabla<VectorField>(X,F),M.LieDerivative(X,F));
			for (int i=1;i<=3;++i)
			{
				ex F=f;		
				ex nonlinearterm=(connection.Nabla<DifferentialForm>(X,F)*e(i)).expand();
				TS_ASSERT_DIFFERS(connection.Nabla<DifferentialForm>(X,F*e(i))-(F*connection.Nabla<DifferentialForm>(X,e(i))).expand(),0);
				TS_ASSERT_EQUALS(connection.Nabla<DifferentialForm>(X,F*e(i))-(F*connection.Nabla<DifferentialForm>(X,e(i))).expand(),nonlinearterm);
				F=f*g;
				nonlinearterm=(connection.Nabla<DifferentialForm>(X,F)*e(i)).expand();
				TS_ASSERT_DIFFERS(connection.Nabla<DifferentialForm>(X,F*e(i))-(F*connection.Nabla<DifferentialForm>(X,e(i))).expand(),0);
				TS_ASSERT_EQUALS(connection.Nabla<DifferentialForm>(X,F*e(i))-(F*connection.Nabla<DifferentialForm>(X,e(i))).expand(),nonlinearterm);
				F=pow(f,2)+sin(g);
				nonlinearterm=(connection.Nabla<DifferentialForm>(X,F)*e(i)).expand();
				TS_ASSERT_DIFFERS(connection.Nabla<DifferentialForm>(X,F*e(i))-(F*connection.Nabla<DifferentialForm>(X,e(i))).expand(),0);
				TS_ASSERT_EQUALS(connection.Nabla<DifferentialForm>(X,F*e(i))-(F*connection.Nabla<DifferentialForm>(X,e(i))).expand(),nonlinearterm);
				F=pow(f,2)*log(g);
				nonlinearterm=(connection.Nabla<DifferentialForm>(X,F)*e(i)).expand();
				TS_ASSERT_DIFFERS(connection.Nabla<DifferentialForm>(X,F*e(i))-(F*connection.Nabla<DifferentialForm>(X,e(i))).expand(),0);
				TS_ASSERT_EQUALS(connection.Nabla<DifferentialForm>(X,F*e(i))-(F*connection.Nabla<DifferentialForm>(X,e(i))).expand(),nonlinearterm);
			}
	}		

	//test torsionfreeconnection.h (with dtable)
	void testConnection2() {
			S3 M;
		{		//standard frame
			ExVector e=M.e();
			LeviCivitaConnection<true> leviCivita(&M,RiemannianStructure(&M,e));
			TS_ASSERT_EQUALS(leviCivita.Torsion(),ExVector(M.Dimension()));
			
			TS_ASSERT_EQUALS(leviCivita.Nabla<DifferentialForm>(e(1),e(1)),0);
			TS_ASSERT_EQUALS(leviCivita.Nabla<DifferentialForm>(e(1),e(2)),-e(3)/2);
			TS_ASSERT_EQUALS(leviCivita.Nabla<DifferentialForm>(e(1),e(3)),e(2)/2);
						
			matrix R=leviCivita.CurvatureForm();		
			TS_ASSERT_EQUALS(R(0,0),0);
			TS_ASSERT_EQUALS(R(1,1),0);
			TS_ASSERT_EQUALS(R(2,2),0);
			TS_ASSERT_EQUALS(R(0,1),-R(1,0));
			TS_ASSERT_EQUALS(R(0,2),-R(2,0));
			TS_ASSERT_EQUALS(R(1,2),-R(2,1));
			TS_ASSERT_EQUALS(R(0,1),(e(1)*e(2)/4).expand());
			TS_ASSERT_EQUALS(R(0,2),(e(1)*e(3)/4).expand());
			TS_ASSERT_EQUALS(R(1,2),(e(2)*e(3)/4).expand());
			
			ex ric=TensorProduct<DifferentialOneForm,DifferentialOneForm>(e[0],e[0])/2+
			TensorProduct<DifferentialOneForm,DifferentialOneForm>(e[1],e[1])/2+
			TensorProduct<DifferentialOneForm,DifferentialOneForm>(e[2],e[2])/2;
			TS_ASSERT_EQUALS(leviCivita.Ricci(),ric);
		}
		{	//non-standard frame
			ExVector v;
			ex sqrt2=sqrt(ex(2));
			v.push_back(1/sqrt2*M.e(1)+M.e(2)/2+M.e(3)/2);
			v.push_back(-1/sqrt2*M.e(1)+M.e(2)/2+M.e(3)/2);
			v.push_back(-M.e(2)/sqrt2+M.e(3)/sqrt2);		
			LeviCivitaConnection<true> leviCivita(&M,RiemannianStructure(&M,v));
		
			TS_ASSERT_EQUALS(leviCivita.Nabla<DifferentialForm>(M.e().dual()(1),M.e(1)),0);
			TS_ASSERT_EQUALS(leviCivita.Nabla<DifferentialForm>(M.e().dual()(1),M.e(2)),-M.e(3)/2);
			TS_ASSERT_EQUALS(leviCivita.Nabla<DifferentialForm>(M.e().dual()(1),M.e(3)),M.e(2)/2);
		
			matrix R=leviCivita.CurvatureForm();		
			TS_ASSERT_EQUALS(R(0,0),0);
			TS_ASSERT_EQUALS(R(1,1),0);
			TS_ASSERT_EQUALS(R(2,2),0);
			TS_ASSERT_EQUALS(R(0,1),-R(1,0));
			TS_ASSERT_EQUALS(R(0,2),-R(2,0));
			TS_ASSERT_EQUALS(R(1,2),-R(2,1));
			TS_ASSERT_EQUALS(R(0,1),(v(1)*v(2)/4).expand());
			TS_ASSERT_EQUALS(R(0,2),(v(1)*v(3)/4).expand());
			TS_ASSERT_EQUALS(R(1,2),(v(2)*v(3)/4).expand());
			
			ex ric=TensorProduct<DifferentialOneForm,DifferentialOneForm>(v[0],v[0])/2+
			TensorProduct<DifferentialOneForm,DifferentialOneForm>(v[1],v[1])/2+
			TensorProduct<DifferentialOneForm,DifferentialOneForm>(v[2],v[2])/2;
			TS_ASSERT_EQUALS(leviCivita.Ricci(),ric);			
		}
	}


	void testExpandBug() 
	{
		AbstractLieGroup<> G("0,0,12,13");
		StructureConstant a(N.a), c(N.c);
		RiemannianStructure g(&G,ParseDifferentialForms(G.e(),"-3+4,-4,2,1"));
		LeviCivitaConnection<true> omega(&G,g);
		ex x=(G.e(2)*G.e(3)-G.e(2)*G.e(4))*G.e(1);
		for (int i=1;i<=G.Dimension();++i) {
			ex l=omega.Nabla<DifferentialForm>(G.e(i),x);
			ex m=omega.Nabla<DifferentialForm>(G.e(i),x.expand());
			TS_ASSERT_EQUALS(l,m);
		}
	}

	void testPseudoLeviCivita()
	{
		AbstractLieGroup<> M("0,12,13");
		exvector v;
		v.push_back(M.e(1)+M.e(2));
		v.push_back(M.e(2));
		v.push_back(M.e(3)+2*M.e(1));
		Frame e(v);

		StandardPseudoRiemannianStructure g(&M,e,2);
		PseudoLeviCivitaConnection leviCivita(&M,g);
		for (int i=1;i<=3;++i)
		for (int j=1;j<=3;++j) {
			ex tij=	leviCivita.Nabla<VectorField>(M.e(i),M.e(j))-
				leviCivita.Nabla<VectorField>(M.e(j),M.e(i))-
				M.LieBracket(M.e(i),M.e(j));
			TS_ASSERT(tij.expand().is_zero());
		}
		for (int i=1;i<=3;++i)
		for (int j=i+1;j<=3;++j) {
			ex tij=	leviCivita.Nabla<VectorField>(e(i),e(j))-
				leviCivita.Nabla<VectorField>(e(j),e(i))-
				M.LieBracket(e(i),e(j));
			TS_ASSERT(tij.expand().is_zero());
		}
		for (int i=1;i<=3;++i)
		for (int j=1;j<=3;++j)
		for (int k=j;k<=3;++k)
	 	{
			ex nabla_ej=leviCivita.Nabla<DifferentialForm>(e(i),e(j));
			ex nabla_ek=leviCivita.Nabla<DifferentialForm>(e(i),e(k));
			ex tijk=g.ScalarProduct().OnForms(nabla_ej,e(k)) +g.ScalarProduct().OnForms(nabla_ek,e(j));
			TS_ASSERT(tijk.expand().is_zero());
		}
		for (int i=1;i<=3;++i)
		for (int j=1;j<=3;++j)
		for (int k=j;k<=3;++k)
	 	{
			ex nabla_ej=leviCivita.Nabla<DifferentialForm>(M.e(i),M.e(j));
			ex nabla_ek=leviCivita.Nabla<DifferentialForm>(M.e(i),M.e(k));
			ex tijk=g.ScalarProduct().OnForms(nabla_ej,M.e(k)) +g.ScalarProduct().OnForms(nabla_ek,M.e(j));
			TS_ASSERT(tijk.expand().is_zero());
		}

		matrix ricci=leviCivita.RicciAsMatrix();
		
		AbstractLieSubgroup<false> M_other_frame(M,e.dual());
		
		StandardPseudoRiemannianStructure g_other_frame(&M_other_frame,M_other_frame.e(),2);
		PseudoLeviCivitaConnection leviCivita_other_frame(&M_other_frame,g_other_frame);

		matrix ricci_other_frame=leviCivita_other_frame.RicciAsMatrix();
		TS_ASSERT(!ricci.is_zero_matrix());
		TS_ASSERT_EQUALS(ricci, ricci_other_frame);

	}
};



//manifoldwith.h
class GStructureTestSuite : public CxxTest::TestSuite  {
public:	
	class ES5Manifold : public ManifoldWith<SU2Structure> {
	public:
		ES5Manifold()
		{
			for (int i=1;i<=5;i++)			
				DeclareNabla<DifferentialForm>(e(i),alpha(),-Hook(e(i),omega1()));
			//Declare_d(alpha(),-2*omega1());c
			Declare_d(omega1(),0);
			Declare_d(omega2(),3*alpha()*omega3());
			Declare_d(omega3(),-3*alpha()*omega2());
		}
	
		void DoTest() {
			TS_ASSERT_EQUALS(d(alpha()),-2*omega1());			
			TS_ASSERT_EQUALS(d(omega2()),(3*alpha()*omega3()).expand());
			TS_ASSERT_EQUALS(d(omega3()),-(3*alpha()*omega2()).expand());
			TS_ASSERT_EQUALS(d(omega1()),0);
			TS_ASSERT_EQUALS(d(3*alpha()*omega3()),0);
			TS_ASSERT_EQUALS(d(-3*alpha()*omega2()),0);
#ifdef HAVE_COCOA
			ex Ric=LeviCivita().Ricci();
			LOG_INFO((NormalForm<Tensor<DifferentialOneForm,DifferentialOneForm> >(Ric)));
			PolyBasisImplementation<ConnectionParameter,DefaultPolyAlgorithms> I;
			LeviCivita().GetEquations_ddZero(I);
#ifdef LONG_TEST
			I.Reduce();	//speeds up by ~25% the call to ReduceModuloIdeal			
			Ric=I.ReduceModuloIdeal<Tensor<DifferentialOneForm,DifferentialOneForm> >(Ric);
			ex g;
			for (int i=1;i<=5;i++)
				g+=TensorProduct<DifferentialOneForm,DifferentialOneForm>(e(i),e(i));
			TS_ASSERT_EQUALS(Ric,4*g);
#endif
#endif
		}
	};

	void testIT() {
		ES5Manifold M;
		M.DoTest();	
	}

};


class U2R3 : public ManifoldWithCoordinates {
	exvector CreateFrame() {
		exvector l;
		for (int i=1;i<=4;++i) l.push_back(DifferentialOneForm(N.e(i)));
		return l;		
	}
	ex r;
	ex AtPrincipalPoint(ex alpha, ex r) const {
		return alpha.subs(lst{x(1)==r,x(2)==0,x(3)==0});
	}
	ex approximate(ex alpha, ex r) const {
		ex res=AtPrincipalPoint(alpha,r);
		ExVector delta(3);
		for (int i=1;i<=3;++i)
		{
			delta(i)=x(i)-AtPrincipalPoint(x(i),r);
			res+=AtPrincipalPoint(alpha.diff(ex_to<symbol>(x(i))),r)*delta(i);
		}
		for (int i=1;i<=3;++i)
		for (int j=1;j<=3;++j)
			res+=AtPrincipalPoint(alpha.diff(ex_to<symbol>(x(i))).diff(ex_to<symbol>(x(j))),r)*delta(i)*delta(j)/2;
		res=res.expand();
		return res;
	}

public:
	U2R3 () : ManifoldWithCoordinates(CreateFrame(),3) {
		const int dimension=4;
		ExVector de=ParseDifferentialForms(e(),"23, 31, 12, 0");
		for (int i=1;i<=dimension;++i)
			Declare_d(e(i),de(i));
		r=Function(N.r);
		Declare_d(r,(x(1)*d(x(1))+x(2)*d(x(2))+x(3)*d(x(3)))/r);

		ex b1=d(x(1))+x(2)*e(3)-x(3)*e(2);
		ex b2=d(x(2))-x(1)*e(3)+x(3)*e(1);
		ex b3=d(x(3))+x(1)*e(2)-x(2)*e(1);
		
		Frame f;
		f.push_back(r*b1);
		f.push_back(b2);
		f.push_back(b3);
		f.push_back(e(4));
		f.push_back(e(1));
		f.push_back(e(2));
		f.push_back(e(3));
		{
			TransverseRiemannianStructure P(this, f,dimension);
			TS_ASSERT_THROWS(P.AssertIntegrable(),WedgeException<std::invalid_argument>);
		}
	
		ExVector y;
		y.push_back((x(1)*b1+x(2)*b2+x(3)*b3)/r);
		ex s=pow(x(1),2)+pow(x(2),2);
		y.push_back((-x(2)*b1+x(1)*b2)/sqrt(s));
		y.push_back((x(1)*x(3)*b1+x(2)*x(3)*b2-s*b3)/(r*sqrt(s)));
		y.push_back(e(4));
		y.push_back(e(1));
		y.push_back(e(2));
		y.push_back(e(3));		
		f.clear();
		for (exvector::const_iterator i=y.begin();i!=y.end();++i)
			f.push_back(approximate(i->subs(r==sqrt(pow(x(1),2)+pow(x(2),2)+pow(x(3),2))),1));
		TransverseRiemannianStructure P(this,f,dimension);
		lst eqns;
		P.GetEquationsIntegrable(eqns);
//		lst subs(x(3)==sqrt(pow(r,2)-pow(x(1),2)-pow(x(2),2)));
		lst subs{x(1)==1,x(2)==0,x(3)==0,r==1};//there are too many square roots for simplications to work generically, so we work at a point.
		for (lst::const_iterator i=eqns.begin();i!=eqns.end();++i)
			TS_ASSERT_EQUALS(i->subs(subs).expand().normal(),0);
		TransverseLeviCivitaConnection omega(this,P);
		for (int i=1;i<=dimension;++i) 
		{
			LOG_INFO(f(i));
			LOG_INFO(omega.H<DifferentialForm>(f(i)));
			TS_ASSERT_EQUALS((omega.H<VectorField>(f.dual()(i))-f.dual()(i)).expand().subs(subs).expand().normal(),0);
			TS_ASSERT_EQUALS((omega.H<DifferentialForm>(f(i))-f(i)).expand().subs(subs).expand().normal(),0);
			TS_ASSERT_EQUALS(omega.V<VectorField>(f.dual()(i)).expand().subs(subs).expand().normal(),0);	
			TS_ASSERT_EQUALS(omega.V<DifferentialForm>(f(i)).expand().subs(subs).expand().normal(),0);
		}	
		LOG_MSG("I got here");
		for (int i=dimension+1;i<=f.size();++i) 
		{
			TS_ASSERT_EQUALS((omega.V<VectorField>(f.dual()(i))-f.dual()(i)).subs(subs).expand().normal(),0);	
			TS_ASSERT_EQUALS((omega.V<DifferentialForm>(f(i))-f(i)).subs(subs).expand().normal(),0);
			TS_ASSERT_EQUALS(omega.H<VectorField>(f.dual()(i)).subs(subs).expand().normal(),0);
			TS_ASSERT_EQUALS(omega.H<DifferentialForm>(f(i)).subs(subs).expand().normal(),0);
		}	
		LOG_MSG("I got here");
		for (int i=1;i<=f.size();++i) 
		for (int j=i+1;j<=f.size();++j) {
			ex X=f.dual()(i), Y=f.dual()(j);
			TS_ASSERT_EQUALS((omega.Nabla<VectorField>(X,Y)-omega.Nabla<VectorField>(Y,X)-LieBracket(X,Y)).subs(subs).expand().normal(),0);
		}
		LOG_MSG("I got here");
		for (int i=1;i<=dimension;++i) {
			for (int j=1;j<=dimension;++j) {
				ex X=f.dual()(i), Y=f.dual()(j);
				ex AXY=omega.A(X,Y);
				TS_ASSERT_EQUALS((AXY+omega.A(Y,X)).subs(subs).expand().normal(),0);
				ex XY=LieBracket(X,Y)/2;
				TS_ASSERT_EQUALS((omega.V<VectorField>(XY)-AXY).subs(subs).expand().normal(),0);
			}
		}
		LOG_MSG("I got here");
		matrix r=omega.BaseRicciAsMatrix(true);
		LOG_MSG("I got here");
		matrix r2=omega.BaseRicciAsMatrix(false);
		for (int i=0;i<dimension;++i)
		for (int j=0;j<dimension;++j)
			TS_ASSERT_EQUALS((r(i,j)-r2(i,j)).subs(subs).expand().normal(),0);
	}
	
};

class SubmersionTestSuite  : public CxxTest::TestSuite  {
public:	
	//test submersion.h
	void testSubmersion() {
		AbstractLieGroup<> S3("23,31,12");
		exvector hforms,vforms;
		hforms.push_back(S3.e(1));
		hforms.push_back(S3.e(2));
		vforms.push_back(S3.e(3));
		Submersion p(S3,hforms,vforms);
		TS_ASSERT(p.IsBasic<DifferentialForm>(S3.e(1)*S3.e(2)));
		TS_ASSERT(p.IsBasic<DifferentialForm>(S3.e(1)*S3.e(2)));
		TS_ASSERT(!p.IsBasic<DifferentialForm>(S3.e(1)));
		TS_ASSERT(!p.IsBasic<DifferentialForm>(S3.e(3)));

		TS_ASSERT_EQUALS(p.H<DifferentialForm>(S3.e(1)+S3.e(3)),S3.e(1));
		TS_ASSERT_EQUALS(p.V<DifferentialForm>(S3.e(1)+S3.e(3)),S3.e(3));
		TS_ASSERT_EQUALS(p.H<VectorField>(S3.e(1)+S3.e(3)),S3.e(1));
		TS_ASSERT_EQUALS(p.V<VectorField>(S3.e(1)+S3.e(3)),S3.e(3));

		TransverseRiemannianStructure P(&p, p.e());
		P.AssertIntegrable();
		TS_ASSERT_EQUALS(P.ScalarProduct().Interior(S3.e(1),S3.e(1)*S3.e(2)),S3.e(2));
		TS_ASSERT_EQUALS(P.HodgeStar(S3.e(1)),S3.e(2));
		TS_ASSERT_EQUALS(P.ScalarProduct().OnForms(S3.e(1),S3.e(1)+S3.e(2)),1);
		TS_ASSERT_EQUALS(P.ScalarProduct().OnVectors(S3.e(1),S3.e(1)+S3.e(2)),1);
		
		logging_level=LOGGING_LEVEL_DEBUG;
		TransverseLeviCivitaConnection E(&S3,P);
		matrix T(3,3), A(3,3),T2(3,3),A2(3,3);
		A=	{{0,-S3.e(3)/2,S3.e(2)/2},
			{S3.e(3)/2,0,-S3.e(1)/2},
			{0,0,0}};
		for (int i=1;i<=3;++i)
		for (int j=1;j<=3;++j)
		{
			T2(i-1,j-1)=E.T(S3.e(i),S3.e(j));
			A2(i-1,j-1)=E.A(S3.e(i),S3.e(j));
		}
		TS_ASSERT_EQUALS(A,A2);
		TS_ASSERT_EQUALS(T,T2);

		TS_ASSERT_EQUALS(E.V<VectorField>(p.LieBracket(S3.e(1),S3.e(2)))/2,E.A(S3.e(1),S3.e(2)));
		matrix Omega=E.BaseCurvatureForm();
		TS_ASSERT_EQUALS(Omega(0,1),S3.e(1)*S3.e(2));
		TS_ASSERT_EQUALS(Omega(1,0),-S3.e(1)*S3.e(2));
		TS_ASSERT_EQUALS(Omega(0,0),0);
		TS_ASSERT_EQUALS(Omega(1,1),0);
		
		matrix R=E.BaseRicciAsMatrix(true);
		TS_ASSERT_EQUALS(R(0,0),1);
		TS_ASSERT_EQUALS(R(1,1),1);
		TS_ASSERT_EQUALS(R(0,1),0);
		TS_ASSERT_EQUALS(R(1,0),0);
		TS_ASSERT_EQUALS(R,E.BaseRicciAsMatrix(false));
	}	
	void testSubmersion2() {
		AbstractLieGroup<> G("0,12,13,14,15,16");
		ExVector e=G.e();
		e(1)=e(1)+sqrt(symbol("a"))*e(2);
		//e(4)=e(4)-log(symbol("b"))*e(2);
		TransverseRiemannianStructure P(&G, e,3);
		P.AssertIntegrable();

		for (int i=1;i<=3;++i)
			TS_ASSERT_EQUALS(P.H<DifferentialForm>(e(i)),e(i));
		for (int i=4;i<=6;++i)
			TS_ASSERT_EQUALS(P.H<DifferentialForm>(e(i)),0);
		for (int i=1;i<=3;++i)
			TS_ASSERT_EQUALS(P.H<VectorField>(P.e().dual()(i)),P.e().dual()(i));
		for (int i=4;i<=6;++i)
			TS_ASSERT_EQUALS(P.H<VectorField>(P.e().dual()(i)),0);
		for (int i=1;i<=3;++i)
			TS_ASSERT_EQUALS(P.V<DifferentialForm>(e(i)),0);
		for (int i=4;i<=6;++i)
			TS_ASSERT_EQUALS(P.V<DifferentialForm>(e(i)),e(i));
		for (int i=1;i<=3;++i)
			TS_ASSERT_EQUALS(P.V<VectorField>(P.e().dual()(i)),0);
		for (int i=4;i<=6;++i)
			TS_ASSERT_EQUALS(P.V<VectorField>(P.e().dual()(i)),P.e().dual()(i));
			
		TransverseLeviCivitaConnection omega(&G, P);
		matrix r=omega.BaseRicciAsMatrix(true);
		matrix r2=omega.BaseRicciAsMatrix(false);
		TS_ASSERT_EQUALS(r,r2);
	}
	void testSubmersion3() {
	//this fails horribly but I suspect there is something wrong conceptually...
		//cout<<"skipping testSubmersion3"<<endl;
//		U2R3 g;		
	}
	
};
#endif /*LIEGROUPS_H_*/

