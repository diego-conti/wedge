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
#include "../structures/transversestructure.h"

namespace Wedge 
{

void ReduceSubstitutions(lst& substitutions, const Simplifier& simplifier)
{
	for (int i=0;i<substitutions.nops();++i) {
		ex& op = substitutions.let_op(i);
		if (!is_a<relational>(op)) throw InvalidArgument(__FILE__,__LINE__,op);
		op = (op.lhs() == simplifier.Simplify(op.rhs()));
	}
}


void TransverseRiemannianStructure::Reduce(const Simplifier& simplifier) {
	ReduceSubstitutions(horizontal,simplifier);
	ReduceSubstitutions(horizontalf,simplifier);
	ReduceSubstitutions(vertical,simplifier);
	ReduceSubstitutions(verticalf,simplifier);
}

TransverseRiemannianStructure::TransverseRiemannianStructure(const Manifold* M, const Frame& coframe, int _k) : StandardPseudoRiemannianStructure{M,coframe,M->Dimension()}, k{_k} {
	LOG_INFO(k);
	LOG_INFO(e());
	Frame frame=e().dual();
	LOG_INFO(frame);
	exvector symbols;
	GetSimple<VectorField>(symbols, e());

	for (exvector::const_iterator i=symbols.begin();i!=symbols.end();++i)
	{
		ex fi;
		for (int j=1;j<=k;++j)
			fi+=Wedge::Hook(*i,e(j))*frame(j);
		horizontal.append(*i==fi.expand());
	}
	LOG_INFO(horizontal);
	for (exvector::const_iterator i=symbols.begin();i!=symbols.end();++i)
	{
		ex fi;
		for (int j=k+1;j<=e().size();++j)
			fi+=Wedge::Hook(*i,e(j))*frame(j);
		vertical.append(*i==fi.expand());
	}
	LOG_INFO(vertical);
	{
		exvector images(e().begin(),e().begin()+k);
		images.resize(e().size());
		horizontalf=LinearMapToSubstitutions<DifferentialOneForm>(
			e(),images);
	}
	LOG_INFO(horizontalf);
	{
		exvector images(k);
		images.insert(images.end(),e().begin()+k,e().end());
		verticalf=LinearMapToSubstitutions<DifferentialOneForm>(e(),images);
	}
	LOG_INFO(verticalf);
}

template<> ex TransverseRiemannianStructure::V<VectorField>(ex X) const
{
	return X.subs(vertical).expand();
}

template<> ex TransverseRiemannianStructure::H<VectorField>(ex X) const
{
	return X.subs(horizontal).expand();
}

template<> ex TransverseRiemannianStructure::V<DifferentialForm>(ex X) const
{
	return X.subs(verticalf).expand();
}

template<> ex TransverseRiemannianStructure::H<DifferentialForm>(ex X) const
{
	return X.subs(horizontalf).expand();
}
}

