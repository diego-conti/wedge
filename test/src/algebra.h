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
#ifndef ALGEBRA_H_
#define ALGEBRA_H_

#include "wedge/base/wedgebase.h"
#include "wedge/base/wedgealgebraic.h"
#include "wedge/base/expressions.h"

#include <cxxtest/TestSuite.h>
#include "test.h"

#include "wedge/linearalgebra/vectorspace.h"
#include "wedge/linearalgebra/lambda.h"
#include "wedge/linearalgebra/derivation.h"
#include "wedge/linearalgebra/linear.h"
#include "wedge/linearalgebra/bilinearform.h"
#include "wedge/linearalgebra/affinebasis.h"
#include "wedge/manifolds/concretemanifold.h"

using namespace GiNaC;
using namespace Wedge;


WEDGE_DECLARE_NAMED_ALGEBRAIC(V,Vector)
WEDGE_DECLARE_NAMED_ALGEBRAIC(W,Vector)

class BilinearOperatorsTestSuite : public CxxTest::TestSuite 
{		
	class ABilinearOperator : public IBilinearOperator<LinearOperator<V>,LinearOperator<V> >
	{
	public:
		ex Apply(const V& v, const V& w) const
		{
			return Tensor<V,V>(v,w);		
		}
	};
	
	class ALeftAssociativeBilinearOperator : public IBilinearOperator<AssociativeOperator<Lambda<V> >,LinearOperator<Lambda<V> > >
	{
	public:
		ex Apply(const V& v, const V& w) const
		{
			return static_cast<const Lambda1<V>&>(v)*static_cast<const Lambda1<V>&>(w);		
		}
		ex Apply(const V& v, const Lambda<V>& w) const
		{
			return  static_cast<const Lambda1<V>&>(v)*w;		
		}
	};
public:
	void testBilinearOperators()
	{
		Lambda1<V> v1,v2,v3,v4,v5;
		ABilinearOperator op1;		
		ALeftAssociativeBilinearOperator op2;
		TS_ASSERT_EQUALS(ABilinearOperator::BilinearOperator(v1,v2,&op1),(TensorProduct<V,V>(v1,v2)));
		TS_ASSERT_EQUALS(ABilinearOperator::BilinearOperator(2*v1+3*v2,4*v2,&op1),(TensorProduct<V,V>(8*v1+12*v2,v2)));
		
		TS_ASSERT_EQUALS(ALeftAssociativeBilinearOperator::BilinearOperator(v1,v2,&op2),v1*v2);		
		TS_ASSERT_EQUALS(ALeftAssociativeBilinearOperator::BilinearOperator(v1*v3,v2,&op2),v1*v3*v2);		
		TS_ASSERT_EQUALS(ALeftAssociativeBilinearOperator::BilinearOperator(v1+v4*v3,v2,&op2),v1*v2+v4*v3*v2);		
		TS_ASSERT_EQUALS(ALeftAssociativeBilinearOperator::BilinearOperator(v1*v2*v3,v4*v5,&op2),v1*v2*v3*v4*v5);		
		TS_ASSERT_EQUALS(ALeftAssociativeBilinearOperator::BilinearOperator(v1,v2*v3*v4*v5,&op2),v1*v2*v3*v4*v5);	
	}
	
	//this really tests wedgealgebraic.h
	void testTrivialPairing()
	{
		V v1(N.v(1)),v2(N.v(2));
		Lambda1<V> v3 (N.v(3));
		TS_ASSERT(!is_exactly_a<basic>(ex(v1)));
		TS_ASSERT(is_exactly_a<V>(ex(v1)));
		TS_ASSERT(is_a<basic>(ex(v1)));
		TS_ASSERT_EQUALS(V::static_class_name(),v1.class_name());
//		const GiNaC::tinfo_static_t* q1=&V::tinfo_static;
//		const GiNaC::tinfo_static_t* q2=&basic::tinfo_static;
//		TS_ASSERT_EQUALS(q1,v1.tinfo());
//		TS_ASSERT_EQUALS(q1,v2.tinfo());
//		TS_ASSERT_DIFFERS(q1,q2);		
		
		TS_ASSERT_EQUALS(TrivialPairing<V>(v1,v2),0);		
		TS_ASSERT_EQUALS(TrivialPairing<V>(v1,v3),0);
		TS_ASSERT_EQUALS(TrivialPairing<V>(v2,v2),1);
		TS_ASSERT_EQUALS(TrivialPairing<V>(v1+v2+v3,v1),1);		
		TS_ASSERT_EQUALS(TrivialPairing<V>(v1+v2+v3,v2),1);		
		TS_ASSERT_EQUALS(TrivialPairing<V>(v1+v2+v3,v3),1);		
		
		DifferentialOneForm x;
		TS_ASSERT_EQUALS(TrivialPairing<VectorField>(x,x),1);
		const VectorField& y=static_cast<const VectorField&>(x);
		TS_ASSERT_EQUALS(TrivialPairing<VectorField>(y,x),1);
		TS_ASSERT_EQUALS(TrivialPairing<VectorField>(x,y),1);
		TS_ASSERT_EQUALS(TrivialPairing<VectorField>(y,y),1);

	}
};

class BilinearFormTestSuite : public CxxTest::TestSuite  {
public:

	void testNeutral() {
		ConcreteManifold M(8);
		NeutralProduct scalar_product(M.e());
		for (int i=1;i<=4;++i) {
			TS_ASSERT_EQUALS(scalar_product.OnOneForms(M.e(i),M.e(1)),0);
			TS_ASSERT_EQUALS(scalar_product.OnVectors(M.e(i),M.e(1)),0);
			TS_ASSERT_EQUALS(scalar_product.OnForms(M.e(i),M.e(1)),0);
			TS_ASSERT_EQUALS(scalar_product.OnOneForms(M.e(i),M.e(i+4)),1);
			TS_ASSERT_EQUALS(scalar_product.OnVectors(M.e(i),M.e(i+4)),1);
			TS_ASSERT_EQUALS(scalar_product.OnForms(M.e(i),M.e(i+4)),1);
			TS_ASSERT_EQUALS(scalar_product.Sharp(M.e(i)),M.e(i+4));
			TS_ASSERT_EQUALS(scalar_product.Flat(M.e(i)),M.e(i+4));
			TS_ASSERT_EQUALS(scalar_product.Sharp(M.e(i+4)),M.e(i));
			TS_ASSERT_EQUALS(scalar_product.Flat(M.e(i+4)),M.e(i));
		}
		TS_ASSERT_EQUALS(scalar_product.Interior(M.e(5),M.e(1)*M.e(2)),M.e(2));
	}

	void testNeutral2()  {
		ConcreteManifold M(8); 
		NeutralProduct g(M.e());
		ExVector e=M.e();
		TS_ASSERT_EQUALS(g.OnVectors(e(1),e(5)),1);
		TS_ASSERT_EQUALS(g.OnVectors(e(5),e(1)+e(5)),1);
		TS_ASSERT_EQUALS(g.OnForms(e(5),e(1)+e(5)),1);
		TS_ASSERT_EQUALS(g.OnOneForms(e(5),e(1)+e(5)),1);
		TS_ASSERT_EQUALS(g.OnForms(e(1)*e(5),e(1)*e(5)),-1);
		TS_ASSERT_EQUALS(g.OnForms(e(1)*e(2),e(1)*e(5)+e(2)*e(6)),0);
		TS_ASSERT_EQUALS(g.OnForms(e(1)*e(2),e(6)*e(5)+e(2)*e(6)),-1);
		TS_ASSERT_EQUALS(g.OnForms(e(1)*e(2)*e(5)*e(8),-e(5)*e(6)*e(1)*e(4)),-1);

		ConcreteManifold N(6);
	}

	void testPositiveDefinite()  {
		ConcreteManifold M(8); 
		PositiveDefiniteScalarProduct g(M.e());
		ExVector e=M.e();
		assert(g.OnVectors(e(1),e(5))==0);
		assert(g.OnVectors(e(5),e(1)+e(5))==1);
		assert(g.OnForms(e(5),e(1)+e(5))==1);
		assert(g.OnOneForms(e(5),e(1)+e(5))==1);
		assert(g.OnForms(e(1)*e(5),e(1)*e(5))==1);
		assert(g.OnForms(e(1)*e(2),e(1)*e(5)+e(2)*e(6))==0);
		assert(g.OnForms(e(1)*e(2),e(6)*e(5)+e(2)*e(6))==0);
		assert(g.OnForms(e(1)*e(2)*e(5)*e(8),-e(1)*e(2)*e(5)*e(8))==-1);
	}
	void testStandardScalarProduct()  {
		ConcreteManifold M(8); 
		StandardScalarProduct g(M.e(),3);
		ExVector e=M.e();
		assert(g.OnVectors(e(1),e(5))==0);
		assert(g.OnVectors(e(3),e(3))==1);
		assert(g.OnVectors(e(4),e(4))==-1);
		assert(g.OnVectors(e(8),e(8))==-1);
		assert(g.OnVectors(e(5),e(1)+e(5))==-1);
		assert(g.OnForms(e(5),e(1)+e(5))==-1);
		assert(g.OnOneForms(e(5),e(1)+e(5))==-1);
		assert(g.OnForms(e(1)*e(5),e(1)*e(5))==-1);
		assert(g.OnForms(e(1)*e(2),e(1)*e(5)+e(2)*e(6))==0);
		assert(g.OnForms(e(1)*e(2),e(6)*e(5)+e(2)*e(6))==0);
		assert(g.OnForms(e(1)*e(2)*e(5)*e(8),-e(1)*e(2)*e(5)*e(8))==-1);
	}

	void testNonStandardFrame() {
		ConcreteManifold M(6); 
		Frame e=ParseDifferentialForms(M.e(),"1+2,2-2*1,3+2*2, 4+2*1,5-6,6");
		StandardScalarProduct g(e,3);
		assert(g.OnVectors(e.dual()(1),e.dual()(5))==0);
		assert(g.OnVectors(e.dual()(3),e.dual()(3))==1);
		assert(g.OnVectors(e.dual()(4),e.dual()(4))==-1);
		assert(g.OnVectors(e.dual()(5),e.dual()(1)+e.dual()(5))==-1);
		assert(g.OnForms(e(5),e(1)+e(5))==-1);
		assert(g.OnOneForms(e(5),e(1)+e(5))==-1);
		assert(g.OnForms(e(1)*e(5),e(1)*e(5))==-1);
		assert(g.OnForms(e(1)*e(2),e(1)*e(5)+e(2)*e(6))==0);
		assert(g.OnForms(e(1)*e(2),e(6)*e(5)+e(2)*e(6))==0);
		assert(g.OnForms(e(1)*e(2),e(1)*e(2)+e(2)*e(6))==1);
	}
};


class BasisTestSuite : public CxxTest::TestSuite 
{
public:
	//test continer semantics in Basis
	void testBasisAsContainer() 
	{
		V v1(N.v(1)),v2(N.v(2)),v3(N.v(3)),v4(N.v(4)),v5(N.v(5));
		exvector l;
		l.push_back(v1);l.push_back(v1+v2);l.push_back(v2+v3+v4);
		Basis<V> basis(l.begin(),l.end());
		LOG_INFO(basis);
		TS_ASSERT(!basis.empty());
		TS_ASSERT_EQUALS(basis.end()-basis.begin(),basis.size());		
		TS_ASSERT_EQUALS(l,ExVector(basis));
		IBasis<V>::const_iterator i;
		int i0,i1;
		for (i0=0,i1=1,i=basis.begin();i0<basis.size() && i1<=basis.size() && i!=basis.end();i0++, i1++, i++)
		{
			TS_ASSERT_EQUALS(basis[i0],*i);
			TS_ASSERT_EQUALS(basis(i1),*i);
		}
		TS_ASSERT(i0==basis.size() && i1>basis.size() && i==basis.end());
		
		list<ex> m;
		m.push_back(v3+v5);
		m.push_back(v4+v2);
		basis.insert(basis.begin()++,m.begin(),m.end());
		l.insert(l.begin()++,m.begin(),m.end());
		TS_ASSERT_EQUALS(basis.end()-basis.begin(),basis.size());		
		TS_ASSERT_EQUALS(l,ExVector(basis));
		basis.push_back(m.front()+m.back());
		TS_ASSERT_EQUALS(l,ExVector(basis));
		basis.push_back(V(N.v(6)));
		TS_ASSERT_EQUALS(l.size(),basis.size()-1);
		
		basis.clear();
		TS_ASSERT_EQUALS(basis.size(),0);
		TS_ASSERT_EQUALS(basis.end()-basis.begin(),basis.size());
		TS_ASSERT(basis.empty());
		basis=exvector(m.begin(),m.end());
		TS_ASSERT_EQUALS(basis.end()-basis.begin(),basis.size());		
		TS_ASSERT_EQUALS(exvector(m.begin(),m.end()),ExVector(basis));		
	}

	//test the behaviour of Basis::Components() using one fixed example 
	void testBasisComponents() {
		V v1(N.v(1)),v2(N.v(2)),v3(N.v(3)),v4(N.v(4)),v5(N.v(5));
		exvector l;
		l.push_back(v1);l.push_back(v1+v2);l.push_back(v1+v3+v4);
		Basis<V> basis(l.begin(),l.end());
		TS_ASSERT_EQUALS(basis.size(),3);
		TS_ASSERT_EQUALS(basis[0],l[0]);
		TS_ASSERT_EQUALS(basis[1],l[1]);
		TS_ASSERT_EQUALS(basis[2],l[2]);

		TS_ASSERT(basis.NotInitialized());
		ExVector components=basis.Components(v1-v2);
		TS_ASSERT(!basis.NotInitialized());
		TS_ASSERT_EQUALS(components.size(),3);
		TS_ASSERT_EQUALS(components[0],2);
		TS_ASSERT_EQUALS(components[1],-1);
		TS_ASSERT_EQUALS(components[2],0);		
		TS_ASSERT_THROWS(basis.Components(v4),NotInSpan);
		components = basis.AllComponents(v4);
		TS_ASSERT_EQUALS(components.size(),4);		
		TS_ASSERT_EQUALS(v4,components[0]*basis[0]+components[1]*basis[1]+components[2]*basis[2]+components[3]*basis[3]);
		
		
		components=basis.Components(v3+v4);
		TS_ASSERT_EQUALS(components.size(),3);
		TS_ASSERT_EQUALS(components[0],-1);
		TS_ASSERT_EQUALS(components[1],0);
		TS_ASSERT_EQUALS(components[2],1);		
	}

	Basis<V> ConstructBasis()
	{
		V v1(N.v(1)),v2(N.v(2)),v3(N.v(3)),v4(N.v(4)),v5(N.v(5));
		exvector l;
		l.push_back(v1);
		l.push_back(v2+rand()*v1);
		l.push_back(v4+rand()*v1+rand()*v2+rand()*v3);
		LOG_INFO(l);
		return Basis<V>(l.begin(),l.end());
	}
	
	void testBasisDual() 
	{
		Basis<V> basis=ConstructBasis();		
		for (int i=0;i<basis.size();i++)
			for (int j=0;j<basis.size();j++)
				TS_ASSERT_EQUALS(TrivialPairing<V>(basis.dual()[i],basis[j]),i==j ? 1 : 0)
		TS_ASSERT_EQUALS(basis.size(),basis.dual().size());
	}
	
	void testSubBasis()
	{
		exvector e=ConstructBasis();
		exvector::iterator i=e.begin();
		exvector f; f.push_back(*i);
		SubBasis<V> b(f.begin(),f.end(),++i,e.end());
		TS_ASSERT_EQUALS(exvector(f.begin(),f.end()),exvector(b.begin(),b.end()));
		TS_ASSERT_EQUALS(exvector(i,e.end()),exvector(b.complement_begin(),b.complement_end()));
		
		TS_ASSERT_EQUALS(b.dual().size(),e.size());
		for (int i=0;i<e.size();i++)
			for (int j=0;j<e.size();j++)
				TS_ASSERT_EQUALS(TrivialPairing<V>(b.dual()[i],e[j]),i==j ? 1 : 0)
		
		TS_ASSERT_EQUALS(exvector(e.begin(),e.end()),exvector(b.full_begin(),b.full_end()));
		TS_ASSERT_EQUALS(b.DimensionOfContainingSpace(),e.size());
		V v(N.v);
		e.push_back(v);
		b.insert(b.end(),e.begin(),e.end());
		TS_ASSERT_EQUALS(exvector(e.begin(),e.end()),exvector(b.begin(),b.end()));
		TS_ASSERT_EQUALS(b.complement_begin(),b.complement_end());
		TS_ASSERT_EQUALS(exvector(e.begin(),e.end()),exvector(b.full_begin(),b.full_end()));
		TS_ASSERT_EQUALS(b.DimensionOfContainingSpace(),e.size());

		V w(N.w);
		exvector g; g.push_back(w);				
		b.AddToComplement(g.begin(),g.end());

		TS_ASSERT_EQUALS(exvector(e.begin(),e.end()),exvector(b.begin(),b.end()));
		TS_ASSERT_EQUALS(g,exvector(b.complement_begin(),b.complement_end()));
		TS_ASSERT_EQUALS(b.DimensionOfContainingSpace(),e.size()+1);
	}
//test LinearMapAsSubstitutions
	void testLinearMapAsSubstitutions()
	{
		V x,y,z;
		exvector e,f=ConstructBasis(); lst subs;
		e.push_back(x+y);
		e.push_back(2*x+y+z);
		e.push_back(z-x);
		subs=LinearMapToSubstitutions<V>(e,f,solve_algo::gauss);
		for (int i=0;i<e.size();++i)
			TS_ASSERT_EQUALS(e[i].subs(subs),f[i]);
		subs=LinearMapToSubstitutions<V>(e,f,solve_algo::bareiss);
		for (int i=0;i<e.size();++i)
			TS_ASSERT_EQUALS(e[i].subs(subs),f[i]);
		subs=LinearMapToSubstitutions<V>(e,f,solve_algo::automatic);
		for (int i=0;i<e.size();++i)
			TS_ASSERT_EQUALS(e[i].subs(subs),f[i]);

		subs=LinearMapToSubstitutions<V>(e.rbegin(),e.rend(),f.begin(),f.end(),solve_algo::gauss);
		for (int i=0;i<e.size();++i)
			TS_ASSERT_EQUALS(e[e.size()-i-1].subs(subs),f[i]);
		subs=LinearMapToSubstitutions<V>(e.rbegin(),e.rend(),f.begin(),f.end(),solve_algo::bareiss);
		for (int i=0;i<e.size();++i)
			TS_ASSERT_EQUALS(e[e.size()-i-1].subs(subs),f[i]);
		subs=LinearMapToSubstitutions<V>(e.rbegin(),e.rend(),f.begin(),f.end(),solve_algo::automatic);
		for (int i=0;i<e.size();++i)
			TS_ASSERT_EQUALS(e[e.size()-i-1].subs(subs),f[i]);

		e=ConstructBasis();
		subs=LinearMapToSubstitutions<V>(e,f,solve_algo::gauss);
		for (int i=0;i<e.size();++i)
			TS_ASSERT_EQUALS(e[i].subs(subs),f[i]);
		subs=LinearMapToSubstitutions<V>(e,f,solve_algo::bareiss);
		for (int i=0;i<e.size();++i)
			TS_ASSERT_EQUALS(e[i].subs(subs),f[i]);
		subs=LinearMapToSubstitutions<V>(e,f,solve_algo::automatic);
		for (int i=0;i<e.size();++i)
			TS_ASSERT_EQUALS(e[i].subs(subs),f[i]);

		subs=LinearMapToSubstitutions<V>(e.rbegin(),e.rend(),f.begin(),f.end(),solve_algo::gauss);
		for (int i=0;i<e.size();++i)
			TS_ASSERT_EQUALS(e[e.size()-i-1].subs(subs),f[i]);
		subs=LinearMapToSubstitutions<V>(e.rbegin(),e.rend(),f.begin(),f.end(),solve_algo::bareiss);
		for (int i=0;i<e.size();++i)
			TS_ASSERT_EQUALS(e[e.size()-i-1].subs(subs),f[i]);
		subs=LinearMapToSubstitutions<V>(e.rbegin(),e.rend(),f.begin(),f.end(),solve_algo::automatic);
		for (int i=0;i<e.size();++i)
			TS_ASSERT_EQUALS(e[e.size()-i-1].subs(subs),f[i]);
	}
};

class VectorSpaceTestSuite : public CxxTest::TestSuite
{
	ExVector ConstructBasis()
	{
		V v1(N.v(1)),v2(N.v(2)),v3(N.v(3)),v4(N.v(4)),v5(N.v(5));
		ExVector l;
		l.push_back(v1);
		l.push_back(v2+rand()*v1);
		l.push_back(v4+rand()*v1+rand()*v2+rand()*v3);
		LOG_INFO(l);
		return l;
	}
public:		
	void testVectorSpace() {
		ExVector r=ConstructBasis();
		VectorSpace<V> S1(r);
		TS_ASSERT_EQUALS(ExVector(S1.e()),r);
		TS_ASSERT_EQUALS(S1.Dimension(),r.size());
		for (int i=1;i<=r.size();i++)
			TS_ASSERT_EQUALS(S1.e(i),r(i));
		
		list<ex> basis(r.begin(),r.end());			
		VectorSpace<V> S2(basis.begin(),basis.end(),"c");
		TS_ASSERT_EQUALS(ExVector(S2.e()),r);
		TS_ASSERT_EQUALS(S2.Dimension(),r.size());
		TS_ASSERT_EQUALS(S1.e(),S2.e());
		
		V w(N.w);				
		S2.AddGenerator(w);
		TS_ASSERT_DIFFERS(S1.e(),S2.e());
		TS_ASSERT_EQUALS(*--(S2.e_end()),w);
		TS_ASSERT_EQUALS(S2.Dimension(),S1.Dimension()+1);
		TS_ASSERT(S2.Contains(w));
		TS_ASSERT(S2.Contains(rand()*w+rand()*S2.e(1)));
		TS_ASSERT(!S1.Contains(w));
		TS_ASSERT(!S1.Contains(w+rand()*S1.e(1)));
		TS_ASSERT(!S1.Contains(sqrt(ex(2345))*w+rand()*S1.e(1)));
		
		basis.push_back(w);			
		S1.SetBasis(Basis<V>(basis.begin(),basis.end()));
		TS_ASSERT_EQUALS(S1.e(),S2.e());

		S1.SetBasis(Basis<V>(r.begin(),r.end()));
		TS_ASSERT_DIFFERS(S1.e(),S2.e());
		TS_ASSERT_EQUALS(*--(S2.e_end()),w);
		TS_ASSERT_EQUALS(S2.Dimension(),S1.Dimension()+1);

		for (exvector::iterator i=r.begin();i!=r.end();i++)
			*i+=w;
		S1.AddGenerators(r.begin(),r.end());
		TS_ASSERT_EQUALS(S1.Dimension(),S2.Dimension());
		TS_ASSERT(S1.Contains(w));
		TS_ASSERT(S1.Contains(rand()*w+rand()*S1.e(1)));		
	 }
	
	void testVectorSpaceEquality() {
		ExVector r=ConstructBasis();
		VectorSpace<V> S1(r);
		r.front()+=r.back();
		swap(r[0],r[1]);
		VectorSpace<V> S2(r);
		LOG_INFO(S1);
		LOG_INFO(S2);
		TS_ASSERT_EQUALS(S1,S2);
		TS_ASSERT_EQUALS(S2,S1);

		ExVector r2=ConstructBasis();
		S2.AddGenerators(r2.begin(),r2.end());
		TS_ASSERT_DIFFERS(S1,S2);
		TS_ASSERT_DIFFERS(S2,S1);
		r2[0]=r2[1]-2*r2.back();
		S1.AddGenerators(r2.begin(),r2.end());
		TS_ASSERT_DIFFERS(S1,S2);
		TS_ASSERT_DIFFERS(S2,S1);
	}
	 
	void testSubspace()
	{
		ExVector r=ConstructBasis();
		exvector v; v.push_back(r(1));
		
		Subspace<V> S1(r.begin(),++r.begin(),r.begin(),r.end());
		Subspace<V> S2(r.begin(),++r.begin(),++r.begin(),r.end());
		Subspace<V> S3(v,r);
		Subspace<V> S4(v,r,TrivialPairingOperator<V>());
		
		TS_ASSERT_EQUALS(S1.Dimension(),v.size());
		TS_ASSERT_EQUALS(S2.Dimension(),v.size());
		TS_ASSERT_EQUALS(S3.Dimension(),v.size());
		TS_ASSERT_EQUALS(S4.Dimension(),v.size());

		TS_ASSERT(equal(S1.e_begin(),S1.e_end(),v.begin()));
		TS_ASSERT(equal(S2.e_begin(),S2.e_end(),v.begin()));
		TS_ASSERT(equal(S3.e_begin(),S3.e_end(),v.begin()));
		TS_ASSERT(equal(S4.e_begin(),S4.e_end(),v.begin()));
		TS_ASSERT(equal(S1.complement_begin(),S1.complement_end(),++r.begin()));
		TS_ASSERT(equal(S2.complement_begin(),S2.complement_end(),++r.begin()));
		TS_ASSERT(equal(S3.complement_begin(),S3.complement_end(),++r.begin()));		

		TS_ASSERT_EQUALS(S1.Complement().Dimension(),r.size()-v.size());
		TS_ASSERT_EQUALS(S2.Complement().Dimension(),r.size()-v.size());
		TS_ASSERT_EQUALS(S3.Complement().Dimension(),r.size()-v.size());
		TS_ASSERT_EQUALS(S4.Complement().Dimension(),r.size()-v.size());
		
		TS_ASSERT(equal(S1.full_begin(),S1.full_end(),S2.full_begin()));		
		TS_ASSERT(equal(S1.full_begin(),S1.full_end(),S3.full_begin()));		
		TS_ASSERT(equal(S1.e_begin(),S1.e_end(),S2.e_begin()));		
		TS_ASSERT(equal(S1.e_begin(),S1.e_end(),S3.e_begin()));
		TS_ASSERT(equal(S1.e_begin(),S1.e_end(),S4.e_begin()));				
		
		TS_ASSERT(S1.Contains(r(1)));
		TS_ASSERT(S4.Contains(r(1)));
		TS_ASSERT(!S1.Contains(r(2)));
		TS_ASSERT(!S4.Contains(r(2)));
		
		S1.SetBasis(S1.Complement().e());
		S4.SetBasis(S4.Complement().e());

		TS_ASSERT(!S1.Contains(r(1)));
		TS_ASSERT(!S4.Contains(r(1)));
	
		TS_ASSERT(S1.Contains(S1.Project(r(1))));
		TS_ASSERT(S4.Contains(S4.Project(r(1))));
		TS_ASSERT(S2.Contains(S1.ProjectOnComplement(r(1))));
		TS_ASSERT(S2.Contains(S4.ProjectOnComplement(r(1))));
		TS_ASSERT(S1.Contains(S1.Project(S1.GenericElement())));
		TS_ASSERT(S4.Contains(S4.Project(S4.GenericElement())));
		
		TS_ASSERT_EQUALS(r(1),S1.Project(r(1))+S1.ProjectOnComplement(r(1)));
		TS_ASSERT_EQUALS(r(1),S4.Project(r(1))+S4.ProjectOnComplement(r(1)));

		TS_ASSERT_EQUALS(S1.GenericElement(),S1.Project(S1.GenericElement())+S1.ProjectOnComplement(S1.GenericElement()));
		TS_ASSERT_EQUALS(S4.GenericElement(),S4.Project(S4.GenericElement())+S4.ProjectOnComplement(S4.GenericElement()));

		SubBasis<V> b;
		b.insert(b.end(),r.begin(),++r.begin());
		S1.SetBasis(b);
		TS_ASSERT_EQUALS(S1.Dimension(),1);
		TS_ASSERT_EQUALS(S1.Complement().Dimension(),0);
		S1.SetComplement(++r.begin(),r.end());
		TS_ASSERT(S1.Contains(r(1)));
		TS_ASSERT(!S1.Contains(r(2)));
		TS_ASSERT_EQUALS(S1.Project(r(2)),0);
		S1.SetComplement(r.end(),r.end());
		TS_ASSERT(S1.Contains(r(1)));
		//TS_ASSERT_THROWS(S1.Project(r(2)),InvalidArgument);
	}
	 
	void testVSpace()
	{
		ExVector r=ConstructBasis();
		VectorSpace<V> S(r);		
		TS_ASSERT_EQUALS(S.SubspaceFromEquations(static_cast<ex*>(NULL),static_cast<ex*>(NULL)).Dimension(),S.Dimension());

		LOG_INFO(r);		
		ex eqn=TrivialPairing<V>(S.GenericElement(),r(1));		
		Subspace<V> S1=S.SubspaceFromEquations(&eqn,(&eqn)+1);
		Subspace<V> S2=S.SubspaceFromEquations(&eqn,(&eqn)+1,TrivialPairingOperator<V>());
		Subspace<V> S3=S.Subspace(r.begin(),++r.begin());
		Subspace<V> S4=S.Subspace(r.begin(),++r.begin(),TrivialPairingOperator<V>());
		
		TS_ASSERT_EQUALS(S1.Dimension(),S.Dimension()-1);
		TS_ASSERT_EQUALS(S2.Dimension(),S.Dimension()-1);
		TS_ASSERT_EQUALS(S3.Dimension(),1);
		TS_ASSERT_EQUALS(S4.Dimension(),1);
		
		TS_ASSERT(equal(S1.e_begin(),S1.e_end(),S2.e_begin()));
		TS_ASSERT(equal(S3.e_begin(),S3.e_end(),S4.e_begin()));
		TS_ASSERT(equal(S2.complement_begin(),S2.complement_end(),S4.e_begin()));
		TS_ASSERT(equal(S2.e_begin(),S2.e_end(),S4.complement_begin()));
		
		TS_ASSERT_EQUALS(TrivialPairing<V>(S1.GenericElement(),r(1)),0);
		TS_ASSERT_EQUALS(TrivialPairing<V>(S2.GenericElement(),r(1)),0);
		
	}

	void testAffineSpace() {
		ExVector r=ConstructBasis();
		VectorSpace<V> V1(r.begin(),r.end()),V2(r.begin(),--r.end());
		AffineBasis<VectorSpace<V>::Coordinate> x;
		V1.GetOrthogonalEquations(x,V2,TrivialPairingOperator<V>());
		TS_ASSERT_EQUALS(x.size(),V2.Dimension());
		Basis<V> orth;
		ex v;
		V1.GetSolutions(orth,x.begin(),x.end(),&v);
		TS_ASSERT_EQUALS(v,0);
		TS_ASSERT_EQUALS(orth.size(),V1.Dimension()-V2.Dimension());
		TS_ASSERT_EQUALS(TrivialPairing<V>(V2.GenericElement(),orth(1)),0);

		exvector y;
		V1.GetOrthogonalEquations(y,V2,TrivialPairingOperator<V>());
		y.back()+=1;
		orth.clear();
		TS_ASSERT_THROWS(V1.GetSolutions(orth,y.begin(),y.end()),WedgeException<std::runtime_error>);
		V1.GetSolutions(orth,y.begin(),y.end(),&v);
		TS_ASSERT_DIFFERS(v,0);
		TS_ASSERT_EQUALS(orth.size(),V1.Dimension()-x.size());
		TS_ASSERT_EQUALS(TrivialPairing<V>(V2.GenericElement(),orth(1)),0);

		ex eqn=V1.coordinate(1)+V1.coordinate(2)+1;
		exvector sol;
		V1.GetSolutions(sol,&eqn,&eqn+1,&v);
		TS_ASSERT_EQUALS(sol.size(),V1.Dimension()-1);
		TS_ASSERT_EQUALS(V1.e().Components(v)(1)+V1.e().Components(v)(2)+1,0);
		TS_ASSERT_THROWS(V1.SubspaceFromEquations(&eqn,&eqn+1),WedgeException<std::runtime_error>);
	}

//affinebasis.h
	void testAffineBasis()
	{
		AffineBasis<V> A;
		V v1,v2,v3;
		A.push_back(v1+1);
		A.push_back(v1);
		TS_ASSERT_EQUALS(A.size(),2);
		A.push_back(1);
		TS_ASSERT_EQUALS(A.size(),2);
		A.push_back(v2+v3);
		TS_ASSERT_EQUALS(A.size(),3);

		A.clear();
		ex x=symbol()*(1+v1)+symbol()*(v2+v3);
		GetCoefficients<symbol>(A,x);
		TS_ASSERT_EQUALS(count(A.begin(),A.end(),1+v1),1);		
		TS_ASSERT_EQUALS(count(A.begin(),A.end(),v2+v3),1);
		TS_ASSERT_EQUALS(A.size(),2);

		ExVector r=ConstructBasis();
		VectorSpace<V> V1(r.begin(),r.end());

//		AffineSpace<V> B;
//		B.AddEquation(V1.coordinate(1)+1);
//		B.AddEquation(V1.coordinate(1));
//		TS_ASSERT_EQUALS(B.codimension(),-1);
		
		AffineBasis<VectorSpace<V>::Coordinate> b;
		b.push_back(V1.coordinate(1)+1);
		b.push_back(V1.coordinate(2)+V1.coordinate(1)+2);
		TS_ASSERT_EQUALS(b.size(),2);
		b.push_back(V1.coordinate(2));
		TS_ASSERT_EQUALS(b.size(),3);
		b.push_back(1);
		TS_ASSERT_EQUALS(b.size(),3);

	}

	void testSubspaceFrom() {
		ExVector r=ConstructBasis();
		VectorSpace<V> V1(r.begin(),r.end());
		ex eqn=V1.coordinate(1);
		Subspace<V> S=V1.SubspaceFromEquations(&eqn,&eqn+1);
		TS_ASSERT_EQUALS(S.Dimension(),V1.Dimension()-1);
		TS_ASSERT(S.Contains(V1.e(2)));

		list<ex> eqns;
		eqns.push_back(V1.coordinate(1)+sqrt(ex(3))*V1.coordinate(2));
		eqns.push_back(V1.coordinate(1)+V1.coordinate(3));
		S=V1.SubspaceFromEquations(eqns.begin(),eqns.end());
		TS_ASSERT_EQUALS(S.Dimension(),V1.Dimension()-2);
		TS_ASSERT(S.Contains(V1.e(1)-1/sqrt(ex(3))*V1.e(2)-V1.e(3)));
		exvector a;
		V1.GetSolutions(a,eqns.begin(),eqns.end());
		TS_ASSERT_EQUALS(a.size(),V1.Dimension()-2);
		for (exvector::const_iterator i=a.begin();i!=a.end();++i)
		{
		 	ExVector comps=V1.e().Components(*i);
			TS_ASSERT_EQUALS(comps(1)+sqrt(ex(3))*comps(2),0);
			TS_ASSERT_EQUALS(comps(1)+comps(3),0);
		}
	}

	void testOperators()
	{
		ex r=symbol("r");
		V v,w;
		TS_ASSERT_EQUALS(TrivialPairing<V>(v,w),0);
		TS_ASSERT_EQUALS(TrivialPairing<V>(r*v,v),r);
		TS_ASSERT_EQUALS(TrivialPairing<V>(r*v,v+v*r),r+r*r);
		TS_ASSERT_EQUALS(TrivialPairing<V>(r*v,r*r*v),r*r*r);
	}
};

class MultilinearTestSuite : public CxxTest::TestSuite
{	
	public:
	
	void testpForms()
	{		
		ex vv[5];
		for (int i=0;i<5;++i)
			vv[i]=Lambda1<V>(N.v(i));
		Basis<Lambda1<V> > v(vv,vv+5);
		
		VectorSpace<Lambda<V> > S=pForms(v,1);
		TS_ASSERT_EQUALS(S.Dimension(),5);
		TS_ASSERT(S.Contains(v[0]));
		TS_ASSERT(S.Contains(v[1]));
		TS_ASSERT(S.Contains(v[2]));
		TS_ASSERT(S.Contains(v[3]));
		TS_ASSERT(S.Contains(v[4]));

		S=pForms(v,2);
		TS_ASSERT_EQUALS(S.Dimension(),10);
		for (Basis<Lambda1<V> >::const_iterator i=v.begin();i!=v.end();i++)
			for (Basis<Lambda1<V> >::const_iterator j=i+1;j!=v.end();j++)
				TS_ASSERT(S.Contains(*i * *j));

		S=pForms(v,3);
		TS_ASSERT_EQUALS(S.Dimension(),10);

		S=pForms(v,4);
		TS_ASSERT_EQUALS(S.Dimension(),5);
		TS_ASSERT(S.Contains(v[1]*v[2]*v[3]*v[4]));
		TS_ASSERT(S.Contains(v[0]*v[2]*v[3]*v[4]));
		TS_ASSERT(S.Contains(v[0]*v[1]*v[3]*v[4]));
		TS_ASSERT(S.Contains(v[0]*v[1]*v[2]*v[4]));
		TS_ASSERT(S.Contains(v[0]*v[1]*v[2]*v[3]));

		S=pForms(v,5);
		TS_ASSERT_EQUALS(S.Dimension(),1);
		TS_ASSERT(S.Contains(v[0]*v[1]*v[2]*v[3]*v[4]));
		
		SubBasis<Lambda1<V> > a(vv,vv+3,vv+3,vv+5);
		Subspace<Lambda<V> > W=TwoForms(a);
		TS_ASSERT_EQUALS(W.Dimension(),3);
		TS_ASSERT(W.Contains(v[0]*v[1]));
		TS_ASSERT(W.Contains(v[0]*v[2]));
		TS_ASSERT(W.Contains(v[1]*v[2]));
		TS_ASSERT_EQUALS(W.Complement().Dimension(),7);
	}
	
	void testLambda() 
	{
		V a(N.a);
//		TS_ASSERT_EQUALS(a,Lambda1<V>(a));
//		TS_ASSERT_EQUALS(Lambda1<V>(a),a);
		TS_ASSERT_EQUALS(a,V(a));
		TS_ASSERT_EQUALS(ex(a),ex(V(a)));
//		TS_ASSERT_DIFFERS(ex(a),ex(Lambda1<V>(a)));
		
		Lambda1<V> v,w;			//elements of V\subset\Lambda(V)
		TS_ASSERT_EQUALS(v*w+w*v,0);
		TS_ASSERT_EQUALS((v*w*v).expand(),0);
		TS_ASSERT_EQUALS((v*w*v).expand(),0);

		TS_ASSERT(v.is_equal(v));		
		TS_ASSERT((v*w).is_equal(v*w));
		TS_ASSERT_EQUALS(TrivialPairing<V>(v,v),1);
		TS_ASSERT_EQUALS(TrivialPairing<Lambda<V> >(v*w,v*w),1);
		TS_ASSERT_EQUALS(TrivialPairing<Lambda<V> >(v,v),1);
		
		exvector vector;
		vector.push_back(v);vector.push_back(w);
		ex alpha=Lambda<V>(vector);		
		if (!is_a<Lambda<V> >(alpha)) alpha=-alpha;
		TS_ASSERT(is_a<Lambda<V> >(alpha));
		TS_ASSERT(is_a<ncmul>(alpha));		
	}
	
	void testTensor()
	{
		V v(N.v);
		W w(N.w);
	
		list<ex> vw; vw.push_back(TensorProduct<V,W>(v,w));
		VectorSpace<Tensor<V,W> > VW(vw.begin(),vw.end());
		TS_ASSERT_EQUALS(VW.Dimension(),1);
		
//		TS_ASSERT_THROWS_ANYTHING((VectorSpace<Tensor<V,V> >(vw)));	//bisogna riscrivere VectorSpace
//		TS_ASSERT_THROWS_ANYTHING((TensorProduct<V,W>(v,v)));
//		TS_ASSERT_THROWS_ANYTHING((TensorProduct<V,W>(w,v)));
//		TS_ASSERT_THROWS_ANYTHING((TensorProduct<V,V>(v,w)));
//		TS_ASSERT_THROWS_ANYTHING((TensorProduct<V,V>(w,v)));
	
		V v2(N.v(2));
		list<ex> vv;
		vv.push_back(TensorProduct<V,V>((v+v2),(v+v2)));
		vv.push_back(TensorProduct<V,V>(v,v));
		vv.push_back(TensorProduct<V,V>(v,v2));
		vv.push_back(TensorProduct<V,V>(v2,v));
		vv.push_back(TensorProduct<V,V>(v2,v2));	
		VectorSpace<Tensor<V,V> > VV(vv.begin(),vv.end());
		TS_ASSERT_EQUALS(VV.Dimension(),4);

		ex a=TensorProduct<V,V>((v+v2),(v+v2));
		ex b=TensorProduct<V,V>(-v,v);
		ex c=TensorProduct<V,V>(v,-v2);
		LOG_INFO(a);
		LOG_INFO(b);
		LOG_INFO(c);
		TS_ASSERT_EQUALS((TrivialPairing<Tensor<V,V> >)(a,a),4);
		TS_ASSERT_EQUALS((TrivialPairing<Tensor<V,V> >)(a,b),-1);
		TS_ASSERT_EQUALS((TrivialPairing<Tensor<V,V> >)(a,c),-1);
		TS_ASSERT_EQUALS((TrivialPairing<Tensor<V,V> >)(b,b),1);
		TS_ASSERT_EQUALS((TrivialPairing<Tensor<V,V> >)(b,c),0);
		TS_ASSERT_EQUALS((TrivialPairing<Tensor<V,V> >)(c,c),1);
		Lambda1<V> e(N.e),f(N.f),g(N.g);
		a=TensorProduct<Lambda<V>,V>(e*f,(v+v2));
		b=TensorProduct<Lambda<V>,V>((e+g)*f,v2);
		c=TensorProduct<Lambda<V>,V>(e+f*g,v);
		LOG_INFO(a);
		LOG_INFO(b);
		LOG_INFO(c);
		TS_ASSERT_EQUALS((TrivialPairing<Tensor<Lambda<V>,V> >)(a,a),2);
		TS_ASSERT_EQUALS((TrivialPairing<Tensor<Lambda<V>,V> >)(a,b),1);
		TS_ASSERT_EQUALS((TrivialPairing<Tensor<Lambda<V>,V> >)(a,c),0);
		TS_ASSERT_EQUALS((TrivialPairing<Tensor<Lambda<V>,V> >)(b,b),2);
		TS_ASSERT_EQUALS((TrivialPairing<Tensor<Lambda<V>,V> >)(b,c),0);
		TS_ASSERT_EQUALS((TrivialPairing<Tensor<Lambda<V>,V> >)(c,c),2);

		c=TensorProduct<V,Lambda<V> >(v,e);
		b=TensorProduct<V,Lambda<V> >(v,e);
		TS_ASSERT_EQUALS(c,b);
		TS_ASSERT_DIFFERS(c,0);
		TS_ASSERT_EQUALS((TrivialPairing<Tensor<V,Lambda<V> > >)(c,c),1);

		a=TensorProduct<V,Lambda<V> >((v+v2),e*f);
		b=TensorProduct<V,Lambda<V> >(v2,(e+g)*f);
		c=TensorProduct<V,Lambda<V> >(v,e+f*g);
		LOG_INFO(a);
		LOG_INFO(b);
		LOG_INFO(c);
		TS_ASSERT_EQUALS((TrivialPairing<Tensor<V,Lambda<V> > >)(a,a),2);
		TS_ASSERT_EQUALS((TrivialPairing<Tensor<V,Lambda<V> > >)(a,b),1);
		TS_ASSERT_EQUALS((TrivialPairing<Tensor<V,Lambda<V> > >)(a,c),0);
		TS_ASSERT_EQUALS((TrivialPairing<Tensor<V,Lambda<V> > >)(b,b),2);
		TS_ASSERT_EQUALS((TrivialPairing<Tensor<V,Lambda<V> > >)(b,c),0);
		TS_ASSERT_EQUALS((TrivialPairing<Tensor<V,Lambda<V> > >)(c,c),2);
		
		a=TensorProduct<Lambda<V>,Lambda<V> >(e*f,f*g);
		b=TensorProduct<Lambda<V>,Lambda<V> >(e*f*g,f*g);
		c=TensorProduct<Lambda<V>,Lambda<V> >(f*g,e*f);
		TS_ASSERT_EQUALS((TrivialPairing<Tensor<Lambda<V>,Lambda<V> > >)(a,a),1);
		TS_ASSERT_EQUALS((TrivialPairing<Tensor<Lambda<V>,Lambda<V> > >)(a,b),0);
		TS_ASSERT_EQUALS((TrivialPairing<Tensor<Lambda<V>,Lambda<V> > >)(a,c),0);
		TS_ASSERT_EQUALS((TrivialPairing<Tensor<Lambda<V>,Lambda<V> > >)(b,b),1);
		TS_ASSERT_EQUALS((TrivialPairing<Tensor<Lambda<V>,Lambda<V> > >)(b,c),0);
		TS_ASSERT_EQUALS((TrivialPairing<Tensor<Lambda<V>,Lambda<V> > >)(c,c),1);
		
		a=TensorProduct<V,Lambda<V> >(e,f);
		const Tensor<V,Lambda<V> >& s=ex_to<Tensor<V,Lambda<V> > >(a);
		TS_ASSERT(is_a<Lambda1<V> >(s.v));
		TS_ASSERT(is_a<Lambda1<V> >(s.w));
		TS_ASSERT_EQUALS(a.subs(ex(e)==0),0);

		a=TensorProduct<V,V>(e,f);
		b=TensorProduct<Lambda1<V>,V>(e,f);
		c=TensorProduct<V,Lambda1<V> >(e,f);
		const Tensor<V,V>& p=ex_to<Tensor<V,V> >(a);
		TS_ASSERT(is_a<Lambda1<V> >(p.v));
		TS_ASSERT(is_a<Lambda1<V> >(p.w));
		const Tensor<V,Lambda1<V> >& q=ex_to<Tensor<V,Lambda1<V> > >(c);
		TS_ASSERT(is_a<Lambda1<V> >(q.v));
		TS_ASSERT(is_a<Lambda1<V> >(q.w));
		ex left,right;
		left=ex_to<Tensor<V,V> >(a).v; right=ex_to<Tensor<V,V> >(a).w;
		TS_ASSERT_EQUALS(left,e);
		TS_ASSERT_EQUALS(right,f);
		left=ex_to<Tensor<Lambda1<V>,V> >(b).v; right=ex_to<Tensor<Lambda1<V>,V> >(b).w;
		TS_ASSERT_EQUALS(left,e);
		TS_ASSERT_EQUALS(right,f);
		left=ex_to<Tensor<V,Lambda1<V> > >(c).v; right=ex_to<Tensor<V,Lambda1<V> > >(c).w;
		TS_ASSERT_EQUALS(left,e);
		TS_ASSERT_EQUALS(right,f);

		a=TensorProduct<V,V>(e,f);
		b=TensorProduct<Lambda<V>,Lambda<V> >(e*g,f*g);
		c=TensorProduct<Tensor<V,V>,V>(a,e);
		right=TensorProduct<V,V>(f,f);
		TS_ASSERT_EQUALS(a.subs(ex(e)==f),right);
		TS_ASSERT_EQUALS(b.subs(ex(e)==g),0);
		right=TensorProduct<Lambda<V>,Lambda<V> >(f*g,f*g);
		TS_ASSERT_EQUALS(b.subs(ex(e)==f),right);
		TS_ASSERT_EQUALS(c.subs(ex(f)==0),0);
		TS_ASSERT_EQUALS(c.subs(ex(e)==0),0);
		right=TensorProduct<Tensor<V,V>,V>(a.subs(ex(e)==g),g);
		TS_ASSERT_EQUALS(c.subs(ex(e)==g),right);
		TS_ASSERT_EQUALS((f+c).subs(c==g),f+g);

		a=DifferentialOneForm(N.a), b=DifferentialOneForm(N.b), c=DifferentialOneForm(N.c);
		ex h=DifferentialOneForm(N.h), k=DifferentialOneForm(N.k);
		TS_ASSERT_EQUALS((TensorProduct<DifferentialForm,DifferentialForm>)(a*(a+b), a),(TensorProduct<DifferentialForm,DifferentialForm>)(a*b,a));
		TS_ASSERT_EQUALS((TensorProduct<DifferentialForm,DifferentialForm>)(a,a*(a+b)),(TensorProduct<DifferentialForm,DifferentialForm>)(a,a*b));
		TS_ASSERT_EQUALS((TensorProduct<DifferentialForm,DifferentialForm>)(a*(a+b), a*a),0);
		TS_ASSERT_EQUALS((TensorProduct<DifferentialForm,DifferentialForm>)((c+h)*(a+b)*k,a), (TensorProduct<DifferentialForm,DifferentialForm>)(c*a*k+h*a*k+c*b*k+h*b*k,a));
		TS_ASSERT_EQUALS((TensorProduct<DifferentialForm,DifferentialForm>)((c*a+h*a)*k,a),(TensorProduct<DifferentialForm,DifferentialForm>)(c*a*k+h*a*k,a));
		TS_ASSERT_EQUALS((TensorProduct<DifferentialForm,DifferentialOneForm>)((c*a+h*a)*k,a),(TensorProduct<DifferentialForm,DifferentialOneForm>)(c*a*k+h*a*k,a));
		TS_ASSERT_EQUALS((TensorProduct<DifferentialForm,DifferentialForm>)(a,(c*a+h*a)*k),(TensorProduct<DifferentialForm,DifferentialForm>)(a,c*a*k+h*a*k));
		TS_ASSERT_EQUALS((TensorProduct<DifferentialOneForm,DifferentialForm>)(a,(c*a+h*a)*k),(TensorProduct<DifferentialOneForm,DifferentialForm>)(a,c*a*k+h*a*k));
	}
};

class LinearOperatorsTestSuite : public CxxTest::TestSuite 
{
	template<typename V> class AnAdditiveOperator : public AdditiveOperator<V>, public mul::visitor
	{
	protected:
		void visit(const basic& v)
		{
			this->Result()=0;
		}
		void visit(const V& v)
		{
			this->Result()= 2*v;
		}
		void visit(const mul& v)
		{
			this->Result()=  2*v;
		}
	};
	
	class AnAdditiveOperator_Lambda : public AdditiveOperator<Lambda<V> >, public mul::visitor
	{
	protected:
		void visit(const basic& v)
		{
			this->Result()=0;
		}
		void visit(const V& v)
		{
			this->Result()= 2*v;
		}
		void visit(const Lambda<V>& v)
		{
			this->Result()= 2*v;
		}
		void visit(const mul& v)
		{
			this->Result()= 2*v;
		}
	};
	
	class ANotSoAdditiveOperator : public AdditiveOperator<V>, public mul::visitor
	{
	protected:
		void visit(const basic& v)
		{
			this->Result()=0;
		}
		void visit(const V& v)
		{
			this->Result()= v*v;
		}
		void visit(const mul& v)
		{
			this->Result()=0;
		}
	};

	class ANotSoLinearOperator_Lambda : public LinearOperator<Lambda<V> >
	{
	protected:
		void visit(const basic& v)
		{
			this->Result()=0;
		}
		void visit(const V& v)
		{
			this->Result()= v*v;
		}
		void visit(const Lambda<V>& v)
		{
			this->Result()= v*v;
		}
	};
	
	class ANotSoLinearOperator : public LinearOperator<V>
	{
	protected:
		void visit(const V& v)
		{
			this->Result()= v*v;
		}
	};
	
	template<bool SKEW> class ADerivation : public Derivation<Lambda<V>,SKEW >
	{
	public:
		ex alpha;
		ADerivation() {
			alpha=Lambda1<V>(N.alpha);
			if (SKEW==false) alpha*=Lambda1<V>(N.beta);
		}
	protected:
		void visit(const V& v)
		{			
			this->Result()= alpha*static_cast<const Lambda1<V>&>(v);
		}
	};

	template<bool SKEW> class ADerivationOverSymbol : public DerivationOver<Lambda<V>,symbol,SKEW >
	{
	public:
		ex alpha;
		ADerivationOverSymbol() {
			alpha=Lambda1<V>(N.alpha);
			if (SKEW==false) alpha*=Lambda1<V>(N.beta);
		}	
	protected:
		void visit(const symbol&)
		{
			this->Result()=alpha;
		}
		void visit(const V& v)
		{
			this->Result()= alpha*static_cast<const Lambda1<V>&> (v);
		}
	};
		
public:
	void testLinearOperators()
	{
		V v(N.v),w(N.w);
		AnAdditiveOperator<V> op0;
		TS_ASSERT_EQUALS(op0.RecursiveVisit(v),2*v);
		TS_ASSERT_EQUALS(op0.RecursiveVisit(v+w),2*v+2*w);
		
		ANotSoAdditiveOperator op1;
		TS_ASSERT_EQUALS(op1.RecursiveVisit(v),v*v);
		TS_ASSERT_EQUALS(op1.RecursiveVisit(v+w),v*v+w*w);
		
		ANotSoLinearOperator op2;
		symbol x;
		TS_ASSERT_EQUALS(op2.RecursiveVisit(v),v*v);
		TS_ASSERT_EQUALS(op2.RecursiveVisit(v+w),v*v+w*w);
		TS_ASSERT_EQUALS(op2.RecursiveVisit(v-w),v*v-w*w);		
		TS_ASSERT_EQUALS(op2.RecursiveVisit(x*v+2*w),x*v*v+2*w*w);
		TS_ASSERT_EQUALS(op2.RecursiveVisit(v+sqrt(5)*w),v*v+sqrt(5)*w*w);
		
		Lambda1<V> e(N.e),f(N.f),g(N.g),h(N.h);
		ANotSoLinearOperator_Lambda op3;
		AnAdditiveOperator_Lambda op4;
		TS_ASSERT_EQUALS(op3.RecursiveVisit(e),0);
		TS_ASSERT_EQUALS(op3.RecursiveVisit(e*f),0);
		TS_ASSERT_EQUALS(op3.RecursiveVisit(e+f),0);
		TS_ASSERT_EQUALS(op3.RecursiveVisit(e*f+g*h),0);		
		TS_ASSERT_EQUALS(op4.RecursiveVisit(e),2*e);
		TS_ASSERT_EQUALS(op4.RecursiveVisit(e*f),2*e*f);
		TS_ASSERT_EQUALS(op4.RecursiveVisit(e+f),2*e+2*f);
		TS_ASSERT_EQUALS(op4.RecursiveVisit(e*f+g*h),2*e*f+2*g*h);

		AnAdditiveOperator<Tensor<Lambda<V>,V> > op5;
		ex a=TensorProduct<Lambda<V>,V>(e,f);
		TS_ASSERT_EQUALS(op5.RecursiveVisit(a),2*a);
		a=TensorProduct<Lambda<V>,V>(e*g*h,f);
		TS_ASSERT_EQUALS(op5.RecursiveVisit(a),2*a);		
	}
	
	void testDerivations()
	{
		Lambda1<V> e(N.e),f(N.f),g(N.g),h(N.h);
		ADerivation<false> op1;
		ADerivation<true> op2;
		ADerivationOverSymbol<false> op3;
		ADerivationOverSymbol<true> op4;
		TS_ASSERT_EQUALS(op1.RecursiveVisit(e*f*g).expand(),3*op1.alpha*e*f*g);
		TS_ASSERT_EQUALS(op2.RecursiveVisit(e*f*g).expand(),3*op2.alpha*e*f*g);
		TS_ASSERT_EQUALS(op3.RecursiveVisit(e*f*g).expand(),3*op3.alpha*e*f*g);
		TS_ASSERT_EQUALS(op4.RecursiveVisit(e*f*g).expand(),3*op4.alpha*e*f*g);
		symbol x("x");
		TS_ASSERT_EQUALS(op1.RecursiveVisit(x*e*f*g).expand(),3*x*op1.alpha*e*f*g);
		TS_ASSERT_EQUALS(op2.RecursiveVisit(x*e*f*g).expand(),3*x*op2.alpha*e*f*g);
		TS_ASSERT_EQUALS(op3.RecursiveVisit(x*e*f*g).expand(),op3.alpha*(3*x*e*f*g+e*f*g));
		TS_ASSERT_EQUALS(op4.RecursiveVisit(x*e*f*g).expand(),op4.alpha*(3*x*e*f*g+e*f*g));
		TS_ASSERT_EQUALS(op1.RecursiveVisit(pow(x,3)*e*f*g).expand(),pow(x,3)*3*op1.alpha*e*f*g);
		TS_ASSERT_EQUALS(op2.RecursiveVisit(pow(x,3)*e*f*g).expand(),pow(x,3)*3*op2.alpha*e*f*g);
		TS_ASSERT_EQUALS(op3.RecursiveVisit(pow(x,3)*e*f*g).expand(),op3.alpha*(pow(x,3)*3*e*f*g+3*x*x*e*f*g));
		TS_ASSERT_EQUALS(op4.RecursiveVisit(pow(x,3)*e*f*g).expand(),op4.alpha*(pow(x,3)*3*e*f*g+3*x*x*e*f*g));
		TS_ASSERT_EQUALS(op1.RecursiveVisit(sin(x)*e*f*g).expand(),sin(x)*3*op1.alpha*e*f*g);
		TS_ASSERT_EQUALS(op2.RecursiveVisit(sin(x)*e*f*g).expand(),sin(x)*3*op2.alpha*e*f*g);
		TS_ASSERT_EQUALS(op3.RecursiveVisit(sin(x)*e*f*g).expand(),op3.alpha*(sin(x)*3*e*f*g+cos(x)*e*f*g));
		TS_ASSERT_EQUALS(op4.RecursiveVisit(sin(x)*e*f*g).expand(),op4.alpha*(sin(x)*3*e*f*g+cos(x)*e*f*g));		
		TS_ASSERT_EQUALS(op4.RecursiveVisit(sin(x)),cos(x)*op4.alpha);		
		TS_ASSERT_EQUALS(op4.RecursiveVisit(sin(x)*e).expand(),op4.alpha*(sin(x)*e+cos(x)*e));
		TS_ASSERT_EQUALS(op4.RecursiveVisit(sin(x)*e*f).expand(),op4.alpha*(sin(x)*2*e*f+cos(x)*e*f));
	}
	
};

//expressions.h
class ExpressionsTestStuite : public CxxTest::TestSuite 
{
public:
	void testGetCoefficients() 
	{
		Lambda1<V> x,y,z;
		symbol a,b,c;
		set<ex,ex_is_less> coeffs;
		GetCoefficients<Lambda<V> >(coeffs,a*x*y+b*z+c);
		TS_ASSERT_DIFFERS(coeffs.find(a),coeffs.end());
		TS_ASSERT_DIFFERS(coeffs.find(b),coeffs.end());
		TS_ASSERT_DIFFERS(coeffs.find(c),coeffs.end());

		TS_ASSERT_THROWS(GetCoefficients<V>(coeffs,a*x*y+b*z+c),WedgeException<std::runtime_error>);
		TS_ASSERT_THROWS(GetCoefficients<V>(coeffs,a*x*y),WedgeException<std::runtime_error>);
		TS_ASSERT_THROWS(GetCoefficients<V>(coeffs,b*z+c),WedgeException<std::runtime_error>);
		coeffs.clear();
		GetCoefficients<V>(coeffs,a*(x+z)+b*z);
		TS_ASSERT_DIFFERS(coeffs.find(a+b),coeffs.end());
		TS_ASSERT_DIFFERS(coeffs.find(a),coeffs.end());		
	}
};
//lambda.h
class DegreeTestStuite : public CxxTest::TestSuite 
{
public:
	void testIsOdd()
	{
		Lambda1<V> a,y,z;
		ex x(a);
		TS_ASSERT_EQUALS(Degree<Lambda<V> >(x*y*z),3);
		TS_ASSERT_EQUALS(Degree<Lambda<V> >(x*y*x),0);
		TS_ASSERT_EQUALS(Degree<Lambda<V> >(x*y+symbol("r")*x*z),2);
		TS_ASSERT_THROWS(Degree<Lambda<V> >(x*y+x),InhomogeneousExpression);
		TS_ASSERT_THROWS(Degree<Lambda<V> >(x*y*z+x),InhomogeneousExpression);

		x=sqrt(127)*3*a;
		TS_ASSERT_EQUALS(Degree<Lambda<V> >(x*y*z),3);
		TS_ASSERT_EQUALS(Degree<Lambda<V> >(x*y*x),0);
		TS_ASSERT_EQUALS(Degree<Lambda<V> >(x*y+symbol("r")*x*z),2);
		TS_ASSERT_THROWS(Degree<Lambda<V> >(x*y+x),InhomogeneousExpression);
		TS_ASSERT_THROWS(Degree<Lambda<V> >(x*y*z+x),InhomogeneousExpression);

		x=sqrt(127)*3*a+y;
		TS_ASSERT_EQUALS(Degree<Lambda<V> >(x*y*z),3);
		TS_ASSERT_EQUALS(Degree<Lambda<V> >(x*y*x),0);
		TS_ASSERT_EQUALS(Degree<Lambda<V> >(x*y+symbol("r")*x*z),2);
		TS_ASSERT_THROWS(Degree<Lambda<V> >(x*y+x),InhomogeneousExpression);
		TS_ASSERT_THROWS(Degree<Lambda<V> >(x*y*z+x),InhomogeneousExpression);
	}
	void testDegree()
	{
		Lambda1<V> a,y,z;
		ex x(a);
		TS_ASSERT_EQUALS(IsOdd<Lambda<V> >(x*y*z),true);
		TS_ASSERT_EQUALS(IsOdd<Lambda<V> >(x*y*x),false);
		TS_ASSERT_EQUALS(IsOdd<Lambda<V> >(x*y+symbol("r")*x*z),false);
		TS_ASSERT_THROWS(IsOdd<Lambda<V> >(x*y+x),InhomogeneousExpression);
		TS_ASSERT_EQUALS(IsOdd<Lambda<V> >(x*y*z+x),true);

		x=sqrt(127)*3*a;
		TS_ASSERT_EQUALS(IsOdd<Lambda<V> >(x*y*z),true);
		TS_ASSERT_EQUALS(IsOdd<Lambda<V> >(x*y*x),false);
		TS_ASSERT_EQUALS(IsOdd<Lambda<V> >(x*y+symbol("r")*x*z),false);
		TS_ASSERT_THROWS(IsOdd<Lambda<V> >(x*y+x),InhomogeneousExpression);
		TS_ASSERT_EQUALS(IsOdd<Lambda<V> >(x*y*z+x),true);

		x=sqrt(127)*3*a+y;
		TS_ASSERT_EQUALS(IsOdd<Lambda<V> >(x*y*z),true);
		TS_ASSERT_EQUALS(IsOdd<Lambda<V> >(x*y*x),false);
		TS_ASSERT_EQUALS(IsOdd<Lambda<V> >(x*y+symbol("r")*x*z),false);
		TS_ASSERT_THROWS(IsOdd<Lambda<V> >(x*y+x),InhomogeneousExpression);
		TS_ASSERT_EQUALS(IsOdd<Lambda<V> >(x*y*z+x),true);
	}
};



#endif /*ALGEBRA_H_*/
