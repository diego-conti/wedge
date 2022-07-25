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

#include "wedge/liealgebras/liegroup.h"
#include "wedge/convenience/parse.h"

namespace Wedge {

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

void AbstractLieGroup<true>::Initialize(const exvector& forms) 
{
	assert(forms.size()==Dimension());
	exvector::const_iterator i=forms.begin(),j=e().begin();
	while (i!=forms.end())
	{		
		LOG_DEBUG(*j);
		LOG_DEBUG(*i);
		Has_dTable::Declare_d(*j++,*i++);				
	}
	Check_ddZero();
}



AbstractLieGroup<true>::AbstractLieGroup(const char* structureConstants, const Name& n1) : ConcreteManifold(internal::GetFrameLength(structureConstants))
{
	ex symbols=StructureConstant(n1);
	Initialize(ParseDifferentialForms(e(),structureConstants,symbols));
}


AbstractLieGroup<true>::AbstractLieGroup(const char* structureConstants, const Name& n1, const Name& n2) : ConcreteManifold(internal::GetFrameLength(structureConstants))
{
	Initialize(ParseDifferentialForms(e(),structureConstants, lst{StructureConstant(n1),StructureConstant(n2)}));
}


AbstractLieGroup<true>::AbstractLieGroup(const char* structureConstants, const Name& n1, const Name& n2, const Name& n3) : ConcreteManifold(internal::GetFrameLength(structureConstants))
{
	;
	Initialize(ParseDifferentialForms(e(),structureConstants,lst {StructureConstant(n1),StructureConstant(n2),StructureConstant(n3)}));
}


AbstractLieGroup<true>::AbstractLieGroup(const char* structureConstants, const Name& n1, const Name& n2, const Name& n3, const Name& n4) : ConcreteManifold(internal::GetFrameLength(structureConstants))
{
	Initialize(ParseDifferentialForms(e(),structureConstants,lst{StructureConstant(n1),StructureConstant(n2),StructureConstant(n3),StructureConstant(n4)}));
}


lst& GetNewStructureConstants(lst& variables, const NameRange& n)
{
	for (NameRange::const_iterator i=n.begin();i!=n.end();++i)
		variables.append(StructureConstant(*i));
	return variables;
}

AbstractLieGroup<true>::AbstractLieGroup(const char* structureConstants, const NameRange& n1) : ConcreteManifold
(internal::GetFrameLength(structureConstants))
{
	lst variables;
	Initialize(ParseDifferentialForms(e(),structureConstants,GetNewStructureConstants(variables,n1)));
}

AbstractLieGroup<true>::AbstractLieGroup(const char* structureConstants, const NameRange& n1, const NameRange& n2) : ConcreteManifold
(internal::GetFrameLength(structureConstants))
{
	lst variables;
	GetNewStructureConstants(variables,n1); GetNewStructureConstants(variables,n2);
	Initialize(ParseDifferentialForms(e(),structureConstants,variables));
}

AbstractLieGroup<true>::AbstractLieGroup(const char* structureConstants, const NameRange& n1, const NameRange& n2, const NameRange& n3) : ConcreteManifold
(internal::GetFrameLength(structureConstants))
{
	lst variables;
	GetNewStructureConstants(variables,n1); GetNewStructureConstants(variables,n2); GetNewStructureConstants(variables,n3);
	Initialize(ParseDifferentialForms(e(),structureConstants,variables));
}


AbstractLieGroup<true>::AbstractLieGroup(const char* structureConstants, const Name& n1, const Name& n2, const Name& n3, const Name& n4, const Name& n5) : ConcreteManifold(internal::GetFrameLength(structureConstants))
{
	Initialize(ParseDifferentialForms(e(),structureConstants,lst{StructureConstant(n1),StructureConstant(n2),StructureConstant(n3),StructureConstant(n4),StructureConstant(n5)}));
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

