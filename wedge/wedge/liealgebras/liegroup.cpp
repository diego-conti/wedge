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
#include "wedge/liealgebras/liegroup.h"
#include "wedge/convenience/parse.h"
#include "wedge/convenience/canonicalprint.h"

namespace Wedge {

void LieGroup::canonical_print(ostream& os) const {
	os<<"(";
	auto i=e().begin();
	Wedge::canonical_print(os,d(*i));
	while (++i!=e().end()) {
		os<<",";
		Wedge::canonical_print(os,NormalForm<DifferentialForm>(d(*i)));
	}
	os<<")";
}


namespace internal {

int GetFrameLength(string structureConstants) {
	return count(structureConstants.begin(),structureConstants.end(),',')+1;
}

}

AbstractLieGroup<false>::AbstractLieGroup(const char* structure_constants)  : ConcreteManifold(internal::GetFrameLength(structure_constants))
{
	exvector forms=ParseDifferentialForms(e(),structure_constants);
	assert(forms.size()==Dimension());
	exvector::const_iterator i=forms.begin(),j=e().begin();
	while (i!=forms.end())
	{		
		LOG_DEBUG(*j);
		LOG_DEBUG(*i);
		Declare_d(*j++,*i++);				
	}
	Check_ddZero();
}

list<ex> collate() {
	return {};
}

AbstractLieGroup<true>::AbstractLieGroup(DelegatedConstructor,const char* structure_constants,const lst& parameters)
	: ConcreteManifold(internal::GetFrameLength(structure_constants)) {
	auto de=ParseDifferentialForms(e(),structure_constants,parameters);
	assert(de.size()==ConcreteManifold::Dimension());
	exvector::const_iterator i=de.begin(),j=e().begin();
	while (i!=de.end())
	{		
		LOG_DEBUG(*j);
		LOG_DEBUG(*i);
		Has_dTable::Declare_d(*j++,*i++);				
	}
	Check_ddZero();
}



Subspace<DifferentialForm> LieGroupHasParameters<false>::ClosedForms(int degree) const
{
	if (degree<=0 || degree>Dimension()) throw OutOfRange(__FILE__,__LINE__,degree);
	VectorSpace<DifferentialForm> forms=pForms(degree);
	list<ex> equations;
	GetCoefficients<DifferentialForm> (equations,d(forms.GenericElement()));
	exvector sol;
	forms.GetSolutions(sol,equations.begin(),equations.end());
	return forms.SubspaceFromEquations(equations.begin(),equations.end());
}
VectorSpace<DifferentialForm> LieGroupHasParameters<false>::ExactForms(int degree) const
{	 
	if (degree<0 || degree>Dimension()) throw OutOfRange(__FILE__,__LINE__,degree);
	if (degree<=1) return VectorSpace<DifferentialForm>(); //trivial vector space
	VectorSpace<DifferentialForm> forms=pForms(degree-1);
	exvector basis;
	basis.reserve(forms.Dimension());
	for (int i=1;i<=forms.Dimension();i++)
		basis.push_back(d(forms.e(i)));
	return VectorSpace<DifferentialForm>(basis.begin(),basis.end());
}

bool LieGroupHasParameters<false>::IsUnimodular() const
{
	ex trace;
	for (int i=1;i<=Dimension();i++)	
		trace+=Hook(e().dual()(i),d(e(i)));
	trace=NormalForm<DifferentialForm>(trace);
	return trace.is_zero();
}
vector<int> LieGroupHasParameters<false>::BettiNumbers() const
{
	vector<int> v(Dimension()+1);
	v[0]=1;
	for (int i=1;i<=Dimension();i++)
		v[i]=ClosedForms(i).Dimension()-ExactForms(i).Dimension();
	return v;
}

matrix LieGroup::KillingForm() const
{
	matrix R(Dimension(),Dimension());
	for (int i=0;i<Dimension();++i)
	for (int j=i;j<Dimension();++j)
	{
		for (int k=0;k<Dimension();++k)
			R(i,j)+=TrivialPairing<VectorField>(e()[k],LieBracket(e()[i],LieBracket(e()[j],e()[k])));
		R(j,i)=R(i,j);
	}
	return R;
}

ex LieGroup::ThreeForm() const
{
	matrix Killing=KillingForm();
	ex res;
	for (int i=1;i<=Dimension();++i)
	for (int j=i+1;j<=Dimension();++j)
	for (int k=j+1;k<=Dimension();++k)
	for (int l=1;l<=Dimension();++l)
		res+=-Killing(l-1,k-1)*Hook(e(i)*e(j),d(e(l)))*e(i)*e(j)*e(k);
	return res;
}

}

