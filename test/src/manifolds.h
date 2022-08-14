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
#ifndef MANIFOLDS_
#define MANIFOLDS_

#include <cxxtest/TestSuite.h>
#include "test.h"

#include "wedge/base/normalform.h"

#include "wedge/manifolds/function.h"
#include "wedge/manifolds/coordinates.h"
#include "wedge/manifolds/liederivative.h"



using namespace GiNaC;
using namespace Wedge;
/**
	@author Diego Conti <diego.conti@unimib.it>
*/


//test manifold.h, differentialform.h, concretemanifold.h
class DifferentialFormsTestSuite : public CxxTest::TestSuite  {
public:
	class TestManifold : public Manifold 
	{
		Frame frame;
		exvector CreateFrame() 
		{
			exvector result;
			result.push_back(DifferentialOneForm(N.x));			
			result.push_back(DifferentialOneForm(N.y)+result[0]);			
			result.push_back(DifferentialOneForm(N.z)-result[1]);			
			result.push_back(DifferentialOneForm(N.t)+result[1]);
			return result;
		}
	public:
		TestManifold() : frame(CreateFrame()) {
			Check_ddZero();
		}
		ex d(ex alpha) const {return 0;}
		const Frame& e() const {return frame;}
	};

	class TestConcreteManifold  : public ConcreteManifold, public virtual Has_dTable 
	{
	public:
		TestConcreteManifold() : ConcreteManifold(3){
			Declare_d(e(1),e(2)*e(3));				  
			Declare_d(e(2),e(1)*e(2));
			Declare_d(e(3),0);
			TS_ASSERT_THROWS(Check_ddZero(),WedgeException<runtime_error>);
			Declare_d(e(2),e(1)*e(3));
			Check_ddZero();
			  
			TS_ASSERT_EQUALS(d(e(1)+e(2)),e(2)*e(3)+e(1)*e(3));
			TS_ASSERT_EQUALS(d(e(1)-e(2)),e(2)*e(3)-e(1)*e(3));
			TS_ASSERT_EQUALS(d(e(1)*e(2)),0);
			  
			VectorSpace<DifferentialForm> V=pForms(1);
			TS_ASSERT_EQUALS(V.Dimension(),3);
			TS_ASSERT(V.Contains(e(1)));
			TS_ASSERT(V.Contains(e(2)));
			TS_ASSERT(V.Contains(e(3)));
			V=pForms(2);
			TS_ASSERT_EQUALS(V.Dimension(),3);
			TS_ASSERT(V.Contains(e(1)*e(2)));
			TS_ASSERT(V.Contains(e(2)*e(3)));
			TS_ASSERT(V.Contains(e(3)*e(1)));
			V=pForms(3);
			TS_ASSERT_EQUALS(V.Dimension(),1);
			TS_ASSERT(V.Contains(e(1)*e(2)*e(3)));
			
			TS_ASSERT_EQUALS(HodgeStar(e(1)),e(2)*e(3));
			TS_ASSERT_EQUALS(d(e(1)-e(2)),e(2)*e(3)-e(1)*e(3));
			TS_ASSERT_EQUALS(d(e(1)*e(2)),0);

// FIXME this is a known bug: HodgeStar only works for positive degree forms at the moment
			TS_ASSERT_EQUALS(HodgeStar(1),e(1)*e(2)*e(3));


	  }
  };




	void testManifold()
	{
		TestManifold M;
		ex volume=1;
		for (Frame::const_iterator i=M.e().begin();i!=M.e().end();i++)
			volume*=*i;
		TS_ASSERT_EQUALS(M.HodgeStar(volume),1);

		VectorSpace<DifferentialForm> oneforms=M.pForms(1);		
		VectorSpace<DifferentialForm> twoforms=M.pForms(2);
		VectorSpace<DifferentialForm> threeforms=M.pForms(3);
		TS_ASSERT_EQUALS(M.Dimension(),4);
		TS_ASSERT_EQUALS(oneforms.Dimension(),4);
		TS_ASSERT_EQUALS(twoforms.Dimension(),6);
		TS_ASSERT_EQUALS(threeforms.Dimension(),4);
		
		for (IBasis<DifferentialForm>::const_iterator i=twoforms.e().begin();i!=twoforms.e().end();i++)
		for (IBasis<DifferentialForm>::const_iterator j=twoforms.e().begin();j!=twoforms.e().end();j++)
		{
			ex x=*i*M.HodgeStar(*j);
			TS_ASSERT_EQUALS(x.expand(),TrivialPairing<DifferentialForm>(*i,*j)*volume);
		}

		for (IBasis<DifferentialForm>::const_iterator i=oneforms.e().begin();i!=oneforms.e().end();i++)
		for (IBasis<DifferentialForm>::const_iterator j=oneforms.e().begin();j!=oneforms.e().end();j++)		
			TS_ASSERT_EQUALS((*i*M.HodgeStar(*j)).expand(),TrivialPairing<DifferentialForm>(*i,*j)*volume);
		
		//verify that a Manifold can correctly determine that d of a constant is zero		
		TS_ASSERT_EQUALS(M.d(0),0);
		TS_ASSERT_EQUALS(M.d(-2),0);
		TS_ASSERT_EQUALS(M.d(-sqrt(ex(3))),0);
		TS_ASSERT_EQUALS(M.d(-Pi+sin(ex(1))),0);
		TestConcreteManifold N;
	}
	

	void testForms()
	{
		DifferentialOneForm a,b,c;
		TS_ASSERT(a.is_equal(a));
		TS_ASSERT_EQUALS(a,a);
		ex d=a*(b+c)-a*b;
		TS_ASSERT(is_a<DifferentialForm>(Hook(a,a*b*c)));
		TS_ASSERT(is_a<DifferentialForm>(Hook(c,a*b*c)));
		TS_ASSERT_EQUALS(d.expand(),-c*a);		
		TS_ASSERT_EQUALS(Hook(a,d),c);
		TS_ASSERT_EQUALS(Hook(c,d),-a);
		TS_ASSERT_EQUALS(Hook(c*b,c*b),1);
		TS_ASSERT_EQUALS(Hook(c*b,a*c*b),a);
	}
	
	void testFrame() {
		ConcreteManifold M(5);
		ExVector e=ExVector(M.e());
		for (int i=0;i<e.size();i++)
				for (int j=i+1;i<e.size();i++)
					e[i]+=j*e[j];
		Frame frame(e);
		
		for (int i=0;i<frame.size();i++)
			for (int j=0;j<frame.size();j++) {
				TS_ASSERT_EQUALS(Hook(frame.dual()[i],frame[j]),i==j ? 1 : 0);
			}
// FIXME this is a known bug: Hook only works for positive degree forms at the moment
		for (int i=0;i<frame.size();i++)		
			TS_ASSERT_EQUALS(Hook(1,frame[i]),frame[i]);
		
		for (int i=0;i<frame.size();i++)
		{
			TS_ASSERT_EQUALS(frame[i],e[i]);
			exvector comps=frame.Components(frame[i]);
			for (int j=0;j<frame.size();j++)
				TS_ASSERT_EQUALS(comps[j],i==j? 1:0);

			comps=frame.Components(M.e()[i]);
			ex sum=0;
			for (int j=0;j<comps.size();j++)
				sum+=comps[j]*frame[j];
			TS_ASSERT_EQUALS(M.e()[i],sum);

		}
	}
	
};

#endif /*MANIFOLDS_*/

