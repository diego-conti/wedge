/*******************************************************************************
 *  Copyright (C) 2007-2023 by Diego Conti, diego.conti@unipi.it 
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
#include "wedge/manifolds/concretemanifold.h"
#include "wedge/linearalgebra/vectorspace.h"


namespace Wedge {
 using namespace GiNaC;
	


Manifold::ddNotZeroException::ddNotZeroException(const char* in_file,int at_line, ex _alpha, const Manifold& M) 
	: WedgeException<std::runtime_error> (std::string("d^2 not zero ")+ToString(_alpha),in_file,at_line), alpha(_alpha) 
{
	assert(M.Dimension()>0);	//in principle we cannot reach this code if the dimension is zero...
	structure_constants="("+ToString(M.d(M.e(1)));
	for (int i=2;i<=M.Dimension();++i)
		structure_constants+=","+ToString(M.d(M.e(i)));
	structure_constants+=")";
}	

const char* Manifold::ddNotZeroException::what() const throw()
{
	string result=WedgeException<std::runtime_error>::what();
	result+=", Manifold: \n";
	result+= structure_constants;
	return result.c_str();
}
	

////////////////////////////////////////////////////////////////
//					implement manifold.h
////////////////////////////////////////////////////////////////

ex Manifold::HodgeStar(ex form1) const
{
	ex form=e()[0];
	for (int i=1;i<Dimension();i++)
		form*=e()[i];
	return Hook(form1,form);
}

Manifold::Manifold() : constant_function(Function(Name("One","\\mathbf1")))
{
}


VectorSpace<DifferentialForm> Manifold::pForms(int p) const {
	if (p<=0 || p>Dimension()) throw OutOfRange(__FILE__,__LINE__,p);
	return Wedge::pForms(e(),p);
}

void Manifold::Check_ddZero()
{
	for (int i=0;i<Dimension();i++)
	{
		ex ddei=d(d(e()[i])).expand();
		list<ex> eqns;
		GetCoefficients<DifferentialForm>(eqns,ddei);
		if (!eqns.empty()) {
			LOG_ERROR(d(e()[i]));
			LOG_ERROR(ddei);			
			throw ddNotZeroException(__FILE__,__LINE__,e()[i],*this);
		}
	}
}

ex Manifold::LieDerivative(ex X, ex f) const
{
	return BilinearOperator(X,f,this).expand();	
}

ex Manifold::LieBracket(ex X, ex Y) const {
	ex XY;
	for (int k=1;k<=Dimension();k++)
	{
		ex XYk=LieDerivative(X,TrivialPairing<VectorField>(e(k),Y))-LieDerivative(Y,TrivialPairing<VectorField>(e(k),X));
		XYk-=Hook(Y,Hook(X,d(e(k))));
		XY+=XYk*e().dual()(k);
	}
	return XY;
}

////////////////////////////////////////////////////////////////
//					implement concretemanifold.h
////////////////////////////////////////////////////////////////
ex Has_dTable::Apply(const VectorField& X,const Function& f) const
{
	exmap::const_iterator it=table.find(f);
	if (it!=table.end()) return Hook(X,it->second);
	else {
//		for (it=table.begin();it!=table.end();++it)
//			if (it->second==X)
//				return f.Derive(it->first,*this);
		return f.Derive(X,*this);
	}
}

class Trivial_dOperator : public DerivationOver<DifferentialForm,Function>
{
	ex constant_function_;
public:
	Trivial_dOperator(ex constant_function ) : constant_function_(constant_function) {}
	void visit(const Wedge::Function& f) {
		if (f==constant_function_) Result()=0;
		else throw Manifold::dException(__FILE__,__LINE__,f);
	}
	void visit(const VectorField& alpha) {
		throw Manifold::dException(__FILE__,__LINE__,alpha);
	}
};

ex Manifold::d(ex alpha) const {
	Trivial_dOperator the_d_operator(constant_function);
	alpha.expand().accept(the_d_operator);
	return the_d_operator.GetResult().expand();
}

class Has_dTable::dOperator : 
	public DerivationOver<DifferentialForm,Function>
{
	const Has_dTable& manifold;
public:
	dOperator(const Has_dTable& m) : manifold(m) {}
	void visit(const Wedge::Function& f) {
		exmap::const_iterator it=manifold.table.find(f);
		if (it!=manifold.table.end()) Result()=it->second;
		else {
			Result()=0;
			for (int i=0;i<manifold.Dimension();++i)
				Result()+=manifold.e()[i]*manifold.LieDerivative(manifold.e().dual()[i],f);
		}
	}
	void visit(const VectorField& alpha) {
		exmap::const_iterator it=manifold.table.find(alpha);
		if (it!=manifold.table.end()) Result()=it->second;
		else {
			LOG_ERROR(alpha);
			LOG_ERROR(manifold.table);
			throw dException(__FILE__,__LINE__,alpha);
		}
	}
};

Has_dTable::Has_dTable(const Has_dTable& o) : Manifold(o), the_d_operator(new dOperator(*this))
{
	table=o.table;	
}

Has_dTable::Has_dTable(): the_d_operator(new dOperator(*this))
{	
	Declare_d(constant_function,0);	
}

ex Has_dTable::d(ex alpha) const {
	alpha.expand().accept(*the_d_operator);
	return the_d_operator->GetResult().expand();
}


Frame ConcreteManifold::CreateFrame(int dimension)
{
	if (dimension<1) throw DiscreteManifold(__FILE__,__LINE__,dimension);
	exvector basis;
	for (int i=1; i<=dimension;i++)	
		basis.push_back(DifferentialOneForm(N.e(i)));
	return basis;
}

ConcreteManifold::ConcreteManifold(exvector frame) : frame(frame)
{
}

ConcreteManifold::ConcreteManifold(int dimension) : frame(CreateFrame(dimension))
{
}


}
