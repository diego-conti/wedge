/***************************************************************************
 *   Copyright (C) 2008 by Diego Conti, diego.conti@unimib.it              *
 *                                                                         *
 *   This file is part of Dictionary.                                      *
 *                                                                         *
 *   Dictionary is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 3 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   Dictionary is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with Dictionary; if not, write to the                           *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#include "invariantforms.h"
#include <sstream>
#include "latex.h"
	
void CohomogeneityOneInvariantForms::setupAlgebra(ex orbitType)
{
	vector<list<ex> > found_so_far;			
	int l=0;
	while (l<compositeGenerators.size() &&  !compositeGenerators[l].empty())
	{
		for (list<CompositeElement>::const_iterator j=compositeGenerators[l].begin();j!= compositeGenerators[l].end();j++)
		{
			CompositeElement elem=*j;				
			ex e=CompositeElementToEx(elem).subs(orbitType).expand();
			int degree=Wedge::Degree<DifferentialForm>(e);
			if (found_so_far.size()<degree+1) found_so_far.resize(degree+1);
			found_so_far[degree].push_back(e);
		}
		l++;
	}
	
	algebra.clear();
	for (int degree=0;degree<found_so_far.size();++degree)
	{
		VectorSpace<DifferentialForm> forms_of_this_degree(found_so_far[degree].begin(),found_so_far[degree].end());
		algebra.push_back(forms_of_this_degree);
	}	

}

bool CohomogeneityOneInvariantForms::Eligible(CompositeElement elem) const
{
	//check whether already a composite generator.
	if (find(compositeGenerators[elem.size()].begin(),compositeGenerators[elem.size()].end(),elem)!=
		compositeGenerators[elem.size()].end())
		return false;
	
	CompositeElement::iterator i=elem.begin();
	int l=elem.size()-1;
	if (l==0) return true;
	const list<CompositeElement>& generatorsDegreeL=compositeGenerators[l];
	while (true)
	{
		SimpleElement toRemove=*i;
		CompositeElement::iterator j=elem.erase(i);
		if (j==elem.end()) break;		
		if (find(generatorsDegreeL.begin(),generatorsDegreeL.end(),elem)==generatorsDegreeL.end()) 
		{
//			for (list<CompositeElement>::const_iterator i=generatorsDegreeL.begin();i!=generatorsDegreeL.end();++i)
//			{
//				LOG_DEBUG(ex(elem));
//				LOG_DEBUG(ex(*i));
//				LOG_DEBUG(i->size());
//				LOG_DEBUG(elem.size());
//				LOG_DEBUG(elem.compare_same_type(*i));
//			}
//			LOG_DEBUG(ex(elem));
			return false;			
		}		
		elem.insert(i=j,toRemove);
	}
	return true;
}

/** @brief Used internally by CohomogeneityOneInvariantForms
 */
struct EligibleElements {
	vector<CompositeElement> elems;
	exvector forms;
	void push_back(const CompositeElement& elem, ex form) {elems.push_back(elem);forms.push_back(form);}
	const CompositeElement& FormToElem(ex form)
	{
		exvector::const_iterator it=find(forms.begin(),forms.end(),form);		
		assert (it!=forms.end());
		int n=it-forms.begin();
		assert(n>=0 && n<elems.size());
		return elems[n];		
	}				
};

/** @brief Overloaded output operator, used for logging
 */
std::ostream& operator<<(std::ostream& os, const EligibleElements& x)
{	
	if (!x.elems.empty()) {
		vector<CompositeElement>::const_iterator i=x.elems.begin();
		os<<ex(*i++);
		while (i!=x.elems.end())
			os<<", " <<ex(*i++);
	}
	return os<<endl;
}

void CohomogeneityOneInvariantForms::AddGeneratorsRelativeToOrbitType(ex orbitType)
{
	setupAlgebra(orbitType);
	compositeGenerators.resize(10);
	//setup compositeGenerators[1]. This is because simpleGenerators need not be independent.
	for (SimpleElement i=simpleGenerators.begin();i!=simpleGenerators.end();i++)
	{
		CompositeElement elem;
		elem.push_back(i);		
		ex new_elem=CompositeElementToEx(elem).subs(orbitType).expand();
		int degree=Wedge::Degree<DifferentialForm>(new_elem);
		if(algebra.size()<degree+1) algebra.resize(degree+1);
		assert(degree>=0);
		if (degree==0 && !new_elem.is_zero())	//add the radial function
			compositeGeneratorsSpace[0].AddGenerator(elem);
		if (algebra[degree].AddGenerator(new_elem))
		{			
			LOG_DEBUG(ex(elem));
			assert (compositeGenerators.size()>=2);
			compositeGenerators[1].push_back(elem);
			if (compositeGeneratorsSpace.size()<degree+1) compositeGeneratorsSpace.resize(degree+1);
			compositeGeneratorsSpace[degree].AddGenerator(elem);											
		}
	}
	

	int l=1;	
	//setup compositeGenerators[l], l>1
	while (l<compositeGenerators.size() && !compositeGenerators[l].empty())				 
	{
		vector<EligibleElements> eligibleElementsByDegree(algebra.size());
		for (list<CompositeElement>::const_iterator j=compositeGenerators[l].begin();j!= compositeGenerators[l].end();++j)
		{				
			CompositeElement elem=*j;
			LOG_DEBUG(ex(elem));
			assert(elem.size()==l);
			ex e=CompositeElementToEx(elem);
								
			for (list<CompositeElement>::const_iterator k=compositeGenerators[1].begin();k!=compositeGenerators[1].end();++k)
			{				
				SimpleElement i=k->front();
				if (i<elem.back()) continue;				
				elem.push_back(i);
				if (Eligible(elem)) {										
					LOG_DEBUG("eligible element"<<ex(elem));
					ex new_elem=(e * *i).subs(orbitType).expand();
					LOG_DEBUG(new_elem);
					int degree=Wedge::Degree<DifferentialForm>(new_elem);
					assert(degree>=0);																				
					if(algebra.size()<degree+1) algebra.resize(degree+1);
					//check if already in algebra; otherwise, mark as eligible 
					if (degree<=max_degree && !algebra[degree].Contains(new_elem))
					{
						if (eligibleElementsByDegree.size()<degree+1) eligibleElementsByDegree.resize(degree+1);
						eligibleElementsByDegree[degree].push_back(elem,new_elem);
					}
				}
				elem.pop_back();
			}
		}

		for (int degree=1;degree<eligibleElementsByDegree.size();degree++)
		{
			LOG_INFO(eligibleElementsByDegree[degree]);
			int i=algebra[degree].Dimension();
			algebra[degree].AddGenerators(eligibleElementsByDegree[degree].forms.begin(),eligibleElementsByDegree[degree].forms.end());
			while (i<algebra[degree].Dimension())
			{				
				CompositeElement elem=eligibleElementsByDegree[degree].FormToElem(algebra[degree].e()[i]);
				LOG_DEBUG(algebra[degree].e()[i]);
				LOG_DEBUG(ex(elem));
				if (compositeGenerators.size()<=l+1)
					compositeGenerators.resize(2*(l+2));
				compositeGenerators[l+1].push_back(elem);	//notice that this procedure produces a sorted list
					
				if (compositeGeneratorsSpace.size()<degree+1) compositeGeneratorsSpace.resize(degree+1);
				compositeGeneratorsSpace[degree].AddGenerator(elem);	
				
				++i;
			}			
		}		
		
		++l;
	}	
	
}

ex CohomogeneityOneInvariantForms::ExToCompositeElement(ex form) const
{
	form=form.expand();
	if (form.is_zero()) return 0;
	ex result;
	int degree;
	try {
		degree=Wedge::Degree<DifferentialForm>(form);
	}
	catch (const InhomogeneousExpression&)
	{
		Wedge::internal::NormalFormHelper<DifferentialForm> v;
		form.accept(v);
		exvector bydegree(algebra.size());
		for (exmap::const_iterator i=v.coeffs.begin();i!=v.coeffs.end();++i)
		{
			int degree=Wedge::Degree<DifferentialForm>(i->first);
			assert(degree>0 || i->second==0);
			assert(degree<bydegree.size());
			bydegree[degree]+=i->second*i->first;
		}		
		for (int i=1;i<bydegree.size();++i)
			result+=ExToCompositeElement(bydegree[i]);
		return result;
	}
	if (degree==0) {
		return form.normal();
	}
	if (degree<0 || degree>=algebra.size())
	{
		LOG_ERROR(degree);
		LOG_ERROR(form);
		throw InvalidArgument(__FILE__,__LINE__,form);
	}
	try {
		ExVector lambda_i=algebra[degree].Components(form,solve_algo::gauss);
		assert(lambda_i.size()==p_forms(degree).Dimension());
		assert(algebra[degree].Dimension()==p_forms(degree).Dimension());
		int j=1;
		LOG_DEBUG(form);	
		LOG_DEBUG(algebra[degree].e());	
		LOG_DEBUG(lambda_i);
		for (int i=0;i<lambda_i.size();i++,j++)
			result+=lambda_i[i].normal()* p_forms(degree).e(j);
		LOG_DEBUG(result);
		assert(j==p_forms(degree).Dimension()+1);
	}
	catch (const NotInSpan&)
	{
		LOG_ERROR(Equation("form",NormalForm<DifferentialForm>(form)));
		ExVector algebradegree;
		for (int i=1;i<=algebra[degree].Dimension();++i)
			algebradegree.push_back(NormalForm<DifferentialForm>(algebra[degree].e(i)));
		LOG_ERROR(algebradegree);
		LOG_ERROR(algebra[degree].e().AllComponents(form));
		throw;
	}
	return result;	
} 

namespace Dictionary_internal {
class CompositeToExOperator: public AdditiveOperator<CompositeElement>, public ncmul::visitor, public power::visitor, public mul::visitor
{
	const CohomogeneityOneInvariantForms* forms;
public:
	CompositeToExOperator(const CohomogeneityOneInvariantForms* forms) {this->forms=forms;}
	void visit(const CompositeElement& c) 	
	{
		Result()=forms->CompositeElementToEx(c);
	}
	void visit(const ncmul& c)
	{
		exvector v;
		v.reserve(c.nops());		
		for (unsigned i=0;i<c.nops();i++)
		{			
			ex image=this->RecursiveVisit(c.op(i));
			v.push_back(image);
		}
		this->Result()=ncmul(v);
	}
	void visit(const power& p)
	{
		this->Result()=power(this->RecursiveVisit(p.op(0)),p.op(1));
	}
	void visit(const basic& m)
	{
		this->Result()=m;
	}
	void visit(const mul& m)
	{
		exvector v; v.reserve(m.nops());
		for (int i=0;i<m.nops();++i)
			v.push_back(this->RecursiveVisit(m.op(i)));
		this->Result()=mul(v);
	}
};

}

ex CohomogeneityOneInvariantForms::CompositeElementToEx(ex form) const
{
	Dictionary_internal::CompositeToExOperator v(this);
	form.accept(v);
	return v.GetResult();
}


