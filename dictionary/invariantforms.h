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
#ifndef INVARIANTFORMS_H_
#define INVARIANTFORMS_H_
/** @file invariantforms.h
 *  @brief The main algorithm, independent of Lie group and representation data
 */

#include <wedge/liegroup.h>
#include <wedge/function.h>
#include <wedge/logging.h>
using namespace GiNaC;
using namespace Wedge;

/** @brief Base class for objects that have a "name", to be used in output.
*/
class Named {
	string _name,_texname;
public:
	Named(string name, string texname) : _name(name),_texname(texname) {}
	Named() {}	
	Named(const Name& name)  {_name=name.plain(); _texname=name.tex();}

	string name() const {return _name;}
	string texname() const {return _texname;}
	void set_name(string name, string texname) {_name=name; _texname=texname;}
	void set_name(const Name& name) {_name=name.plain(); _texname=name.tex();}
};

/** @brief Type of a formal contraction \f$b(\alpha_1,\dotsc,\alpha_k)\f$
*/
class NamedForm : public ex, public Named
{
public:
	NamedForm() {}
	NamedForm(const Name& name, ex form) : ex(form), Named(name) {}
	NamedForm(string name, string texname, ex form) : ex(form), Named(name,texname) {}
};

/** @brief The type of a "simple" generator, i.e. a formal contraction 
  @internal Represented by an iterator pointing to an element of CohomogeneityOneInvariantForms::simpleGenerators
*/
typedef vector<NamedForm>::const_iterator SimpleElement;

/** @brief Type of a product of formal contractions, i.e. an element of \f$\mathcal{C}\f$

Composite elements of the algebra are by definitions (wedge products of) lists of simple generators. So a composite element is a list \f$(a_1,...a_l)\f$ of (pointers to) simple generators. We say that \f$l\f$ is the length of the element.

Length and lexicographic order among elements of same length define a "canonical" choice for a basis, consisting of the
''smallest'' non-zero element \f$\alpha\f$, the smallest element independent of \f$\alpha\f$, and so on. 

Given a composite generator \f$(a_1,...,a_l)\f$, we consider elements \f$(a_1,...,a_l,b)\f$ where 
- \f$b\geq a_l\f$
- \f$(a_1, ...,a_{k-1},a_{k+1},...,a_l,b)\f$ is in the basis of generators of length \f$l\f$

Applying this criteria one obtains from the set of ''smallest'' generators of length \f$l\f$
a set that is guaranteed to contain the set of smallest generators of length \f$l+1\f$. 
*/
class CompositeElement : public list<SimpleElement>, public Register<CompositeElement,Vector>::Algebraic {
	friend class CohomogeneityOneInvariantForms;
public:
	static const char* static_class_name() {return "CompositeElement";}
	unsigned return_type() const {
		return return_types::noncommutative;
		for (const_iterator i=begin();i!=end();++i)
			if (Wedge::Degree<DifferentialForm>(**i)>0)
				return return_types::noncommutative;
		return return_types::commutative;
	}
#if RTTI_METHOD <= 1
	tinfo_t return_type_tinfo() const {return this->get_class_info_static().options.get_id();}
#endif
	int degree() const {
		int degree=0;
		for (const_iterator i=begin();i!=end();++i)
			degree+=Wedge::Degree<DifferentialForm>(**i);
	return degree;
	}		

	void print(const print_context &c, unsigned = 0) const {
		if(empty()) return;
		bool LaTeX=dynamic_cast<const print_latex*>(&c)!=NULL;
		const_iterator i=begin();
		c.s<<	(LaTeX? (*i)->texname(): (*i)->name());
		while (++i!=end())
			c.s<<( LaTeX? "\\," + (*i)->texname(): " " + (*i)->name());
	}		

	bool operator==(const CompositeElement& other) const
	{
		return compare_same_type(other)==0;
	}
	bool operator!=(const CompositeElement& other) const
	{
		return compare_same_type(other)!=0;
	}
	int compare_same_type(const GiNaC::basic & o) const
	{
		assert(is_a<CompositeElement>(o));
		const CompositeElement& other=static_cast<const CompositeElement&>(o);			
		const_iterator i=begin(),j=other.begin();
		while (true)			
		{
			if (i==end())
				return (j==other.end())? 0 : -1;
			else if (j==other.end())
				return 1;
			else if (*i<*j)
				return -1;
			else if (*i>*j)
				return 1;				
			++i,++j;				
		}
		assert(false);			
	}
};


/** @brief A class to compute and represent the space of invariant forms on a cohomogeneity one
 * associated bundle \f$G\times_H V\f$, where H acts sphere-transitively on V.
 * 
 * @remark The algorithm is taken from  [D. Conti, Invariant forms, associated bundles and Calabi-Yau metrics, J. Geom. Phys. 57 (2007) 2483-2508]. The documentation uses notation from the paper.
*/
class CohomogeneityOneInvariantForms : public virtual Logging {
public:	
/** @brief Construct the CohomogeneityOneInvariantForms object, without actually starting the calculations
 @param max_degree If specified, forms of degree higher than max_degree are ignored (which may speed things up)
 @internal Debug output can be enabled by passing LOGGING_LEVEL_DEBUG to the constructor of base class Logging
*/	
	CohomogeneityOneInvariantForms(int max_degree=1000) : Logging (LOGGING_LEVEL_INFO) {
		this->max_degree=max_degree;
	}	
//	typedef list<CompositeElement> CompositeGeneratorList;	///<A list of composite generators

private:
	int max_degree;	//only forms up to degree max_degree are computed
	vector<NamedForm> simpleGenerators;	
	vector< list<CompositeElement> >  compositeGenerators; 		//compositeGenerators[l] represents a sorted list of generators of length l	
	vector<VectorSpace<DifferentialForm> > algebra;	
	vector<VectorSpace<CompositeElement> > compositeGeneratorsSpace;
	
	void setupAlgebra(ex orbitType); ///< Setup the vector algebra using compositeGenerators, and evaluating at the point orbitType

	void AddGeneratorsRelativeToOrbitType(ex orbitType); ///< Add to compositeGenerators the generators that arise considering the given orbitType
	
// check whether \f$(a_1, ...,a_{k-1},a_{k+1},...,a_l,b)\f$ is in the basis of generators of length \f$l\f$
	bool Eligible (CompositeElement elem) const;
public:
/** @brief Conversion routine
 * @param alpha An invariant form in this algebra (as a linear combination of DifferentialForm's) or a function
 * @return The form \f$\alpha\f$ as a linear combination of CompositeElement's
 * 
 * @note The form can be of mixed degree, but in this case it cannot have a scalar part (because GiNaC does not allow a numeric scalar term)
 */
	ex ExToCompositeElement(ex alpha) const;
/** @brief Conversion routine
 * @param alpha An invariant form in this algebra (as a linear combination of CompositeElement's)
 * @return The form \f$\alpha\f$ as a linear combination of DifferentialForm's
 */
	ex CompositeElementToEx(ex alpha) const;

/** @overload
*/
	ex CompositeElementToEx(const CompositeElement& alpha) const
	{
		ex e=1;
		for (CompositeElement::const_iterator i=alpha.begin();i!=alpha.end();i++)
			e*=**i;
		return e;	
	}
	
/** @brief Return a generator of the space of invariant functions (zero-forms) 
*/
	ex InvariantFunction() const
	{
		LOG_DEBUG(compositeGeneratorsSpace[0]);
		assert(!compositeGeneratorsSpace.empty() && compositeGeneratorsSpace[0].Dimension()>0); 
		assert(!compositeGeneratorsSpace[0].e(1).is_zero());
		return compositeGeneratorsSpace[0].e(1);	
	}

/** @brief Information about a specific representation.
 * 
 * Each data member is a lst of equations in the coordinates that characterize the relevant point.
 */  
	struct RepresentationInfo {
		ex specialPoint;
		ex principalPoint;
		ex genericPoint;
		void AppendCondition(ex condition)
		{
			assert(is_a<lst>(specialPoint)); 
			lst l=ex_to<lst>(specialPoint);
			l.append(condition);
			specialPoint=l;
			assert(is_a<lst>(principalPoint)); 
			l=ex_to<lst>(principalPoint);
			l.append(condition);
			principalPoint=l;
			assert(is_a<lst>(genericPoint)); 
			l=ex_to<lst>(genericPoint);
			l.append(condition);
			genericPoint=l;
		}
	};

/** @brief Create the algebra of invariant forms generated by a given list of "simple" elements 
 * @param [from, to) A range of NamedForm representing the generators. Need not be independent. 
 * @param representationInfo Contains information about the fibre \f$V\f$ and the action of \f$H\f$ on it.
*/
	template<typename InputIterator> void CreateFromGenerators(InputIterator from, InputIterator to, RepresentationInfo representationInfo)
	{
		simpleGenerators=vector<NamedForm>(from,to);
		list<CompositeElement> ZeroLengthGenerators;
		ZeroLengthGenerators.push_back(CompositeElement());
		compositeGenerators.resize(1);
		compositeGenerators[0]=ZeroLengthGenerators;
		compositeGeneratorsSpace.resize(1);	//compositeGeneratorsSpace[0] always contains an invariant function
		AddGeneratorsRelativeToOrbitType(representationInfo.specialPoint);		
		LOG_INFO(compositeGeneratorsSpace);
		AddGeneratorsRelativeToOrbitType(representationInfo.principalPoint);
		LOG_INFO(compositeGeneratorsSpace);
		
//Having computed a basis, recompute the generators at generic point. This would
//not be necessary if one is only interested in working at principal point
		for (int i=1;i< compositeGeneratorsSpace.size();i++)
		{
				ExVector e=compositeGeneratorsSpace[i].e();
				ExVector f; f.reserve(e.size());
				for (ExVector::const_iterator j=e.begin();j!=e.end();j++)
					f.push_back(CompositeElementToEx(*j).subs(representationInfo.genericPoint).expand());
				assert(algebra.size()>i);
				assert(algebra[i].Dimension()==f.size());
				algebra[i].SetBasis(Basis<DifferentialForm>(f.begin(),f.end()));
		}
	}

/** @brief Return the space of invariant \f$p\f$-forms
*/
	VectorSpace<CompositeElement> p_forms(int degree) const
	{
		assert(degree>=0);
		if (degree>=compositeGeneratorsSpace.size()) {
			LOG_WARN(degree);
			return VectorSpace<CompositeElement>();
		}
		return compositeGeneratorsSpace[degree];		
	}
};


#endif /*INVARIANTFORMS_H_*/
