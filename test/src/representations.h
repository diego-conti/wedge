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
#ifndef ALGEBRA_H_
#define ALGEBRA_H_
#include <cxxtest/TestSuite.h>
#include "test.h"
#include "wedge/representations/linearaction.h"
#include "wedge/representations/repsl2.h"
#include "wedge/liealgebras/so.h"
#include "wedge/representations/repso.h"
#include "wedge/representations/repgl.h"
#include "wedge/representations/stabilizer.h"
#include "wedge/manifolds/manifoldwith.h"

using namespace GiNaC;
using namespace Wedge;

class RepresentationsTestSuite : public CxxTest::TestSuite 
{		
public:
	void testSO3() {
		for (int k=3;k<=9;k+=2)
		{
			ConcreteManifold M(k);
			SO3Representation<VectorField> V(M.e());

			for (int i=1;i<=M.Dimension();++i)
			{
				ex test=M.e(i);
				TS_ASSERT_EQUALS(V.X<DifferentialForm>(V.Y<DifferentialForm>(test))-V.Y<DifferentialForm>(V.X<DifferentialForm>(test)),2*V.H<DifferentialForm>(test));
				TS_ASSERT_EQUALS(V.H<DifferentialForm>(V.X<DifferentialForm>(test))-V.X<DifferentialForm>(V.H<DifferentialForm>(test)),2*V.Y<DifferentialForm>(test));
				TS_ASSERT_EQUALS(V.H<DifferentialForm>(V.Y<DifferentialForm>(test))-V.Y<DifferentialForm>(V.H<DifferentialForm>(test)),-2*V.X<DifferentialForm>(test));
			}
			if (M.Dimension()>1) {
				ex test=M.e(1)*M.e(2);
				TS_ASSERT_EQUALS(V.X<DifferentialForm>(V.Y<DifferentialForm>(test))-V.Y<DifferentialForm>(V.X<DifferentialForm>(test)),2*V.H<DifferentialForm>(test));
				TS_ASSERT_EQUALS(V.H<DifferentialForm>(V.X<DifferentialForm>(test))-V.X<DifferentialForm>(V.H<DifferentialForm>(test)),2*V.Y<DifferentialForm>(test));
				TS_ASSERT_EQUALS(V.H<DifferentialForm>(V.Y<DifferentialForm>(test))-V.Y<DifferentialForm>(V.H<DifferentialForm>(test)),-2*V.X<DifferentialForm>(test));
			}
			if (M.Dimension()>2) {
				ex test=M.e(1)*M.e(2)*M.e(3);
				TS_ASSERT_EQUALS(V.X<DifferentialForm>(V.Y<DifferentialForm>(test))-V.Y<DifferentialForm>(V.X<DifferentialForm>(test)),2*V.H<DifferentialForm>(test));
				TS_ASSERT_EQUALS(V.H<DifferentialForm>(V.X<DifferentialForm>(test))-V.X<DifferentialForm>(V.H<DifferentialForm>(test)),2*V.Y<DifferentialForm>(test));
				TS_ASSERT_EQUALS(V.H<DifferentialForm>(V.Y<DifferentialForm>(test))-V.Y<DifferentialForm>(V.H<DifferentialForm>(test)),-2*V.X<DifferentialForm>(test));
			}
		}
	}

	void testAction() {
		ConcreteManifold M(5);
		SO3Representation<VectorField> V(M.e());
		ex x,left,right;
		x=TensorProduct<VectorField,VectorField>(M.e(1),M.e(2));
		left=V.H<Tensor<VectorField,VectorField> > (x);
		right=TensorProduct<VectorField,VectorField> (V.H<VectorField>(M.e(1)),M.e(2))+
			TensorProduct<VectorField,VectorField> (M.e(1),V.H<VectorField>(M.e(2)));
		TS_ASSERT_EQUALS(left,right);
		left=V.X<Tensor<VectorField,VectorField> > (x);
		right=TensorProduct<VectorField,VectorField> (V.X<VectorField>(M.e(1)),M.e(2))+
			TensorProduct<VectorField,VectorField> (M.e(1),V.X<VectorField>(M.e(2)));
		TS_ASSERT_EQUALS(left,right);
		left=V.Y<Tensor<VectorField,VectorField> > (x);
		right=TensorProduct<VectorField,VectorField> (V.Y<VectorField>(M.e(1)),M.e(2))+
			TensorProduct<VectorField,VectorField> (M.e(1),V.Y<VectorField>(M.e(2)));
		TS_ASSERT_EQUALS(left,right);

		x=TensorProduct<VectorField,DifferentialForm>(M.e(1),M.e(2)*M.e(3));
		left=V.H<Tensor<VectorField,DifferentialForm> > (x);
		right=TensorProduct<VectorField,DifferentialForm> (V.H<VectorField>(M.e(1)),M.e(2)*M.e(3))+
			TensorProduct<VectorField,DifferentialForm> (M.e(1),V.H<DifferentialForm>(M.e(2)*M.e(3)));
		TS_ASSERT_EQUALS(left,right);
		left=V.X<Tensor<VectorField,DifferentialForm> > (x);
		right=TensorProduct<VectorField,DifferentialForm> (V.X<VectorField>(M.e(1)),M.e(2)*M.e(3))+
			TensorProduct<VectorField,DifferentialForm> (M.e(1),V.X<DifferentialForm>(M.e(2)*M.e(3)));
		TS_ASSERT_EQUALS(left,right);
		left=V.Y<Tensor<VectorField,DifferentialForm> > (x);
		right=TensorProduct<VectorField,DifferentialForm> (V.Y<VectorField>(M.e(1)),M.e(2)*M.e(3))+
			TensorProduct<VectorField,DifferentialForm> (M.e(1),V.Y<DifferentialForm>(M.e(2)*M.e(3)));
		TS_ASSERT_EQUALS(left,right);

		x=M.e(2)*M.e(3);
		left=V.H<DifferentialForm> (x);
		right=V.H<DifferentialForm>(M.e(2))*M.e(3)+M.e(2)*V.H<DifferentialForm>(M.e(3));
		TS_ASSERT_EQUALS(left,right);
		left=V.X<DifferentialForm> (x);
		right=V.X<DifferentialForm>(M.e(2))*M.e(3)+M.e(2)*V.X<DifferentialForm>(M.e(3));
		TS_ASSERT_EQUALS(left,right);
		left=V.Y<DifferentialForm> (x);
		right=V.Y<DifferentialForm>(M.e(2))*M.e(3)+M.e(2)*V.Y<DifferentialForm>(M.e(3));
		TS_ASSERT_EQUALS(left,right);
	}

	void testInvariant()
	{
		{
		ConcreteManifold M(5);
		SO3Representation<VectorField> V(M.e());
		TS_ASSERT_EQUALS(Invariant_pForms(V,M.e(),1).Dimension(),0);
		LOG_INFO(1);

		TS_ASSERT_EQUALS(Invariant_pForms(V,M.e(),2).Dimension(),0);
		LOG_INFO(1);

		TS_ASSERT_EQUALS(Invariant_pForms(V,M.e(),3).Dimension(),0);
		LOG_INFO(1);

		TS_ASSERT_EQUALS(Invariant_pForms(V,M.e(),4).Dimension(),0);
		LOG_INFO(1);

		TS_ASSERT_EQUALS(Invariant_pForms(V,M.e(),5).Dimension(),1);
		LOG_INFO(1);

		TS_ASSERT(Invariant_pForms(V,M.e(),5).Contains(M.pForms(5).e(1)));
		VectorSpace<Tensor<VectorField,VectorField> > metrics=InvariantMetrics(V,M.e());
		TS_ASSERT_EQUALS(metrics.Dimension(),1);
		ex metric;
		for (int i=1;i<=M.Dimension();++i)
			metric+=TensorProduct<VectorField,VectorField>(M.e(i),M.e(i));
		TS_ASSERT(metrics.Contains(metric));
		}
		{
		ConcreteManifold M(7);
		SO3Representation<VectorField> V(M.e());
		TS_ASSERT_EQUALS(Invariant_pForms(V,M.e(),1).Dimension(),0);
		LOG_INFO(1);
		TS_ASSERT_EQUALS(Invariant_pForms(V,M.e(),7).Dimension(),1);
		LOG_INFO(1);
		/*VectorSpace<Tensor<VectorField,VectorField> > metrics=InvariantMetrics(V,M.e());
		TS_ASSERT_EQUALS(metrics.Dimension(),1);
		ex metric;
		for (int i=1;i<=M.Dimension();++i)
			metric+=TensorProduct<VectorField,VectorField>(M.e(i),M.e(i));
		TS_ASSERT(metrics.Contains(metric));*/
		}
	}

	void testSO() {
		ConcreteManifold M(5);
		SO G(5);
		SORepresentation<VectorField> V(&G, M.e());

		TS_ASSERT_EQUALS(V.Action<DifferentialForm>(G.A(1,2), M.e(1)), -M.e(2));
		TS_ASSERT_EQUALS(V.Action<DifferentialForm>(G.A(1,2), M.e(2)), M.e(1));
		TS_ASSERT_EQUALS(V.Action<DifferentialForm>(G.A(1,2), M.e(3)), 0);
		TS_ASSERT_EQUALS(V.Action<DifferentialForm>(G.A(1,2), M.e(4)), 0);
		TS_ASSERT_EQUALS(V.Action<DifferentialForm>(G.A(1,2), M.e(5)), 0);

		TS_ASSERT_EQUALS(V.Action<DifferentialForm>(G.A(1,2), M.e(1)*M.e(2)), 0);
		TS_ASSERT_EQUALS(V.Action<DifferentialForm>(G.A(1,2), M.e(2)*M.e(3)), M.e(1)*M.e(3));
		TS_ASSERT_EQUALS(V.Action<DifferentialForm>(G.A(1,2), M.e(3)*M.e(4)), 0);

		TS_ASSERT_EQUALS(V.Action<VectorField>(G.A(1,2), M.e(1)), -M.e(2));
		TS_ASSERT_EQUALS(V.Action<VectorField>(G.A(1,2), M.e(2)), M.e(1));
		TS_ASSERT_EQUALS(V.Action<VectorField>(G.A(1,2), M.e(3)), 0);
		TS_ASSERT_EQUALS(V.Action<VectorField>(G.A(1,2), M.e(4)), 0);
		TS_ASSERT_EQUALS(V.Action<VectorField>(G.A(1,2), M.e(5)), 0);

		for (int i=1;i<=5;++i)
		for (int j=i+1;j<=5;++j)
		{
			TS_ASSERT_EQUALS(V.Action<VectorField>(G.A(i,i),M.e(j)),0);
			TS_ASSERT_EQUALS(V.Action<VectorField>(G.A(i,j),M.e(i)),-M.e(j));
			TS_ASSERT_EQUALS(V.Action<VectorField>(G.A(j,i),M.e(i)),M.e(j));
		}

		VectorSpace<Tensor<VectorField,VectorField> > metrics=InvariantMetrics(V,M.e());
		TS_ASSERT_EQUALS(metrics.Dimension(),1);
		ex metric;
		for (int i=1;i<=M.Dimension();++i)
			metric+=TensorProduct<VectorField,VectorField>(M.e(i),M.e(i));
		TS_ASSERT(metrics.Contains(metric));		
	}	


	//test stabilizer.h
	void testStabilizer()
	{
		{
			ConcreteManifold N(3);
			GL gl(3);
			GLRepresentation<VectorField> V(&gl, N.e());
			TS_ASSERT_EQUALS(OrbitDimension<DifferentialForm>(gl,V,N.e(1)),3);
			TS_ASSERT_EQUALS(Stabilizer<DifferentialForm>(gl,V,N.e(1)).Dimension(),6);
			VectorSpace<DifferentialForm> W=StabilizerAlgebra<DifferentialForm>(gl,V,N.e(1));
			for (int i=1;i<=3;++i)
			{
				TS_ASSERT(W.Contains(gl.A(i,2)));
				TS_ASSERT(W.Contains(gl.A(i,3)));
			}

			TS_ASSERT_EQUALS(Stabilizer<DifferentialForm>(gl,V,N.e(1)+N.e(2)).Dimension(),6);
			TS_ASSERT_EQUALS(OrbitDimension<DifferentialForm>(gl,V,N.e(1)+N.e(2)),3);
			W=StabilizerAlgebra<DifferentialForm>(gl,V,N.e(1)+N.e(2));
			for (int i=1;i<=3;++i)
			{
				TS_ASSERT(W.Contains(gl.A(i,1)-gl.A(i,2)));
				TS_ASSERT(W.Contains(gl.A(i,3)));
			}

			TS_ASSERT_EQUALS(Stabilizer<DifferentialForm>(gl,V,N.e(1)*N.e(2)).Dimension(),6);
			TS_ASSERT_EQUALS(OrbitDimension<DifferentialForm>(gl,V,N.e(1)*N.e(2)),3);
		}		
		{
			ManifoldWith<SU3Structure> M;
			GL gl(M.Dimension());
			GLRepresentation<VectorField> V(&gl, M.e());
			TS_ASSERT_EQUALS(OrbitDimension<DifferentialForm>(gl,V,M.omega()),M.pForms(2).Dimension()); 
			TS_ASSERT_EQUALS(OrbitDimension<DifferentialForm>(gl,V,M.psiplus()),M.pForms(3).Dimension()); 
			TS_ASSERT_EQUALS(OrbitDimension<DifferentialForm>(gl,V,M.omega()+M.psiplus()),36-8); 

			AbstractLieSubgroup<false> SU3=Stabilizer<DifferentialForm>(gl,V,M.omega()+M.psiplus());
			TS_ASSERT_EQUALS(SU3.Dimension(),8);
			TS_ASSERT_EQUALS(SU3.ClosedForms(2).Dimension(),8);
			TS_ASSERT_EQUALS(SU3.ExactForms(2),SU3.ClosedForms(2));
			TS_ASSERT(SU3.IsUnimodular());		
			
			AbstractLieSubgroup<false> SU2=Stabilizer<DifferentialForm>(gl,V,M.e(1)+M.omega()+M.psiplus());
			TS_ASSERT_EQUALS(SU2.Dimension(),3);		
			TS_ASSERT_EQUALS(SU2.BettiNumbers(),AbstractLieGroup<>("23,31,12").BettiNumbers());
			
			Subspace<DifferentialForm> su3=StabilizerAlgebra<DifferentialForm>(gl,V,M.omega()+M.psiplus());
			TS_ASSERT(su3.Contains(gl.A(1,2)-gl.A(2,1)+gl.A(3,4)-gl.A(4,3)-2*gl.A(5,6)+2*gl.A(6,5)));
			TS_ASSERT(!su3.Contains(gl.A(1,2)-gl.A(2,1)));
			TS_ASSERT(!su3.Contains(gl.A(3,4)-gl.A(4,3)));
			TS_ASSERT(!su3.Contains(gl.A(5,6)-gl.A(6,5)));
		}
		{
			ManifoldWith<SU3Structure> M;
			SO so(M.Dimension());
			SORepresentation<VectorField> V(&so, M.e());
			TS_ASSERT_EQUALS(OrbitDimension<DifferentialForm>(so,V,M.omega()),15-9); 
			TS_ASSERT_EQUALS(OrbitDimension<DifferentialForm>(so,V,M.psiplus()),15-8); 
			TS_ASSERT_EQUALS(OrbitDimension<DifferentialForm>(so,V,M.omega()+M.psiplus()),15-8); 

			AbstractLieSubgroup<false> SU3=Stabilizer<DifferentialForm>(so,V,M.psiplus());
			TS_ASSERT_EQUALS(SU3.Dimension(),8);
			TS_ASSERT_EQUALS(SU3.ClosedForms(2).Dimension(),8);
			TS_ASSERT_EQUALS(SU3.ExactForms(2),SU3.ClosedForms(2));
			TS_ASSERT(SU3.IsUnimodular());		
			
			AbstractLieSubgroup<false> SU2=Stabilizer<DifferentialForm>(so,V,M.e(1)+M.psiplus());
			TS_ASSERT_EQUALS(SU2.Dimension(),3);		
			TS_ASSERT_EQUALS(SU2.BettiNumbers(),AbstractLieGroup<>("23,31,12").BettiNumbers());
			
			Subspace<DifferentialForm> su3=StabilizerAlgebra<DifferentialForm>(so,V,M.psiplus());
			TS_ASSERT(su3.Contains(so.A(1,2)+so.A(3,4)-2*so.A(5,6)));
			TS_ASSERT(!su3.Contains(so.A(1,2)));
			TS_ASSERT(!su3.Contains(so.A(3,4)));
			TS_ASSERT(!su3.Contains(so.A(5,6)));
		}
	}
};
#endif /*ALGEBRA_H_*/

