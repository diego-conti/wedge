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
#ifndef LIEGROUP_H
#define LIEGROUP_H
#include "wedge/manifolds/concretemanifold.h"
#include "wedge/convenience/parse.h"
#include "wedge/linearalgebra/vectorspace.h"
#include "wedge/linearalgebra/lambda.h"
#include "wedge/base/parameters.h"
#include "wedge/polynomialalgebra/polybasis.h"

/** @ingroup Manifolds
 *  @{
 * 
 *  @file liegroup.h
 *  @brief Lie groups 
*/


namespace Wedge {

/** @brief Class template with explicit specializations.
 */
template<bool WithParams> class LieGroupHasParameters {};

/////////////////////////////////////////////////////////////////////////////////
//				Abstract base class for Lie groups
/////////////////////////////////////////////////////////////////////////////////

/** @brief A (connected) Lie group, represented in %Wedge by its Lie algebra.
 * 
 * @note The frame associated to a LieGroup object is assumed to be left-invariant.
*/
class LieGroup : public virtual Manifold {
	//ensure noone derives from LieGroup except through LieGroupHasParameters
	template<bool WithParams> friend class LieGroupHasParameters;	
	LieGroup() {}	
public:
/** @brief Returns the structure constants (useful for output)
 * @return The vector \f$de^1,\dotsc,de^n\f$
 */ 
	ExVector StructureConstants() const
	{
		ExVector dei(Dimension());
		for (int i=1;i<=Dimension();i++)			
			dei(i)=NormalForm<DifferentialForm>(d(e(i)));
		return dei;	
	}		
/** @brief Returns the Killing form of the Lie group
 *  @return The Killing form as a symmetric matrix
 */ 
	matrix KillingForm() const;

/** @brief Returns the standard three form of the Lie group
 * @return The three-form \f$\phi(X,Y,Z)=\frac16\langle[X,Y],Z\rangle\f$
 * @todo Optimize
 */ 
	ex ThreeForm() const;

	/** @brief Print out the structure constants in the form (de^1,...,de^n).
	 * Wedge's canonical printing functions are used for consistent output through different runs.
	 */

	ostream& canonical_print(ostream& os) const;
};

/** @brief Overloaded output operator
 * 
 * Display the structure constants relative to a lieGroup in terms of Lie brackets
 * 
 * @todo define a manipulator that enables one to output the forms rather than the vector fields
 */
template<class charT, class traits> std::basic_ostream<charT,traits>& operator<<(std::basic_ostream<charT,traits>& os,const LieGroup& G)
{
	for (int i=1;i<=G.Dimension();i++)
		for (int j=i+1;j<=G.Dimension();j++)
		{
			ex XiXj;
			for (int k=1;k<=G.Dimension();k++)
				XiXj-=TrivialPairing<DifferentialForm>(G.d(G.e(k)),G.e(i)*G.e(j))*G.e().dual()(k);
			if (!XiXj.is_zero()) os<<"["<<G.e(i)<<","<<G.e(j)<<"]="<<XiXj<<endl;
		}
	return os;
}

/////////////////////////////////////////////////////////////////////////////////
//			Template specializations for groups with/without parameters
/////////////////////////////////////////////////////////////////////////////////

/** @brief Abstract base class for Lie groups whose structure constants do not depend on any parameter
 * 
 * @remark Instances of this class are required to satisfy \f$d^2=0\f$. 
*/
template<> class LieGroupHasParameters<false>  : public LieGroup {
public:
/** @brief Compute the space of closed left-invariant forms
 * @param degree An integer
 * @return The space of closed left-invariant forms of the indicated degree
 */
	Subspace<DifferentialForm> ClosedForms(int degree) const;
/** @brief Compute the space of exact left-invariant forms
 * @param degree An integer
 * @return The space of exact left-invariant forms of the indicated degree
 */
	VectorSpace<DifferentialForm> ExactForms(int degree) const;	
/** @brief Compute the Betti numbers of the Lie algebra 
 * @return The Betti numbers as a zero-based vector
 * 
 * @remark If the group is compact, the  Betti numbers of the Lie algebra coincide with the Betti numbers of the
 * group. Notice, however, that this does not hold in general.
 */
	vector<int> BettiNumbers() const;
	
 /** @brief Test whether the group is unimodular
  *  
 * @remark A Lie group is unimodular if \f$\operatorname{tr} \operatorname{ad}(X)=0 \f$ for all \f$X\f$ in the Lie algebra.
 * 
 * @todo Add functions IsNilpotent, IsSolvable, IsCompact, and analogous functions for Lie groups with parameters.
 * A possible implementation would be to: 
 * - define "attribute" classes Nilpotent, Unimodular, and so on
 * - declare template functions Is<class Attribute>, for Lie groups without parameters, and GetEquationFor<class Attribute>, 
 * for Lie groups with parameters
 * - define a specialization for each attribute, outside of the definition of the class (I mean, non-inline)
 * - define a template class that enables one to define, say, an attribute And<Nilpotent, Unimodular> 
 * - define a default partial specialization (or something legal to that effect!) of Is<And<X,Y> > that invokes Is<X> and Is<Y> and then takes the logical and.
 * - when it makes sense, define a full specialization of Is<And<X,Y> >. This could be advantageous in some cases, where 
 * computing two things at once may be better than computing them in sequence.  
 * - find some way of handling combinations of more than two attributes (possibly using boost::mpl, or macros).
 * - apply same technique to GStructure, for intrinsic torsion computations  
 */
	bool IsUnimodular() const;	
};


WEDGE_DECLARE_NAMED_ALGEBRAIC(StructureConstant,realsymbol)

/** @brief Abstract base class for Lie groups with structure constants depending on parameters
 * 
 * The parameters must have type StructureConstant
 * 
 * @remark Instances of this class are not required to satisfy \f$d^2=0\f$ for all values in the parameters. 
 */
template<> class LieGroupHasParameters<true> : 
	public LieGroup, public HasParameters<StructureConstant>, public virtual Has_dTable {
public:
/** @brief Compute the equations in the parameters corresponding to \f$d^2=0\f$
 * @param I A container where the equations are to be stored
 * @return A reference to I 
 * 
 * @sa PolyBasis, PolyBasis_impl 
 */
	template<typename Container> Container& GetEquations_ddZero(Container& I) const
	{			
		for (int i=1;i<=Dimension();i++)
		{
			GetCoefficients<DifferentialForm>(I,d(d(e(i))));			
		}		
		return I;
	}
	
/** @brief Impose conditions on the parameters so that \f$d\alpha=\beta\f$
 * @param alpha,beta Differential forms
 * 
 * The assumption is that the equations \f$d\alpha=\beta\f$, \f$d\beta=0\f$ depend linearly on the parameters.
 */
	void Declare_d(ex alpha, ex beta)
	{
		DeclareZero(d(beta));
		DeclareZero(d(alpha)-beta);
	}	

/** @brief Update dtable by imposing \f$d\alpha=\beta\f$
 * @param alpha,beta Differential forms as in Has_dTable::Declare_d
 */

	void ReplaceIn_dTable(ex alpha, ex beta) {
		Has_dTable::Declare_d(alpha,beta);
	}
private:
	void DeclareConditions(const lst& list_of_equations) {
		for (exmap::const_iterator i=dTable().begin();i!=dTable().end();i++)					
			Has_dTable::Declare_d(i->first,i->second.subs(list_of_equations));
	}
};


typedef LieGroupHasParameters<true> LieGroupWithParameters;
typedef LieGroupHasParameters<false> LieGroupWithoutParameters;


/////////////////////////////////////////////////////////////////////////////////+
//					Ready-for-use Lie group classes
/////////////////////////////////////////////////////////////////////////////////

/** @brief Instances of AbstractLieGroup represent groups defined by their structure constants,
 * instead of being defined as subgroups of \f$GL(n,R)\f$ or another group.
 * 
 * @note This class is not abstract in the C++ sense. Abstract only means ``not (presented as) a subgroup of \f$GL(n,R)\f$''. 
*/
template<bool WithParams=false> class AbstractLieGroup;

template<>
class AbstractLieGroup<false> : public LieGroupHasParameters<false>, public ConcreteManifold, public virtual Has_dTable {
public:
  /** @brief Define a Lie group by the structure constants given in Salamon's notation
   *  @param structure_constants A string of the form, e.g., "0,0,12,-13+2*42"
   * @sa ParseDifferentialForms(const exvector& frame, const char* to_parse)
  */
	AbstractLieGroup(const char* structure_constants);
	AbstractLieGroup(const string& structure_constants): AbstractLieGroup{structure_constants.c_str()} {}
};


template<> 
class AbstractLieGroup<true> : public LieGroupHasParameters<true>, public ConcreteManifold, public virtual Has_dTable {
	struct DelegatedConstructor {} delegate_constructor_tag;
	AbstractLieGroup(DelegatedConstructor,const char* structure_constants,const lst& parameters);
public:
  /** @brief Define a Lie group by the structure constants given in Salamon's notation
   *  @param structure_constants A string of the form, e.g., "0,0,12,-13+2*42"; bracket notation such as [sqrt(3)] is allowed
   *  @params parameters a sequence of symbols, lists of symbols, Name's, or NameIndex's corresponding to parameters used in the definition of the Lie algebra
   * In order to use a parameter, say a, use the bracket notation as in [a]*12. The parameter needs to also be passed in the constructor, either as an ex, or an element of a lst of ex,
   * or by name, in which case a StructureConstant object with that name is created.
   *
   * @sa ParseDifferentialForms(const exvector& frame, const char* to_parse)
  */
	template<typename... Parameters>
	AbstractLieGroup(const char* structure_constants, Parameters&&... parameters);
	template<typename... Parameters>
	AbstractLieGroup(const string& structure_constants, Parameters&&... parameters);
};

list<ex> collate();	//trivial specialization returning empty list
template<typename... Parameters>
list<ex> collate(list<ex> parameter_lst, Parameters&&... parameters) {
	parameter_lst.splice(parameter_lst.end(),collate(std::forward<Parameters>(parameters)...));
	return parameter_lst;
}
template<typename... Parameters>
list<ex> collate(ex parameter, Parameters&&... parameters) {
	list<ex> as_list=is_a<lst>(parameter)? list<ex>{parameter.begin(),parameter.end()}: list<ex>{parameter};
	return collate(move(as_list),std::forward<Parameters>(parameters)...);
}
template<typename... Parameters>
list<ex> collate(const exvector& parameter_lst, Parameters&&... parameters) {
	return collate(list<ex>{parameter_lst.begin(),parameter_lst.end()},std::forward<Parameters>(parameters)...);
}
template<typename... Parameters>
list<ex> collate(const Name& parameter_name, Parameters&&... parameters) {
	ex parameter=StructureConstant{parameter_name};
	return collate(parameter,std::forward<Parameters>(parameters)...);
}
template<typename... Parameters>
list<ex> collate(const NameRange& range, Parameters&&... parameters) {
	list<ex> structure_constants;
	std::transform(range.begin(),range.end(),back_inserter(structure_constants),
			[] (NameAndIndex&& name) {return StructureConstant{move(name)};});
	return collate(move(structure_constants),std::forward<Parameters>(parameters)...);
}

template<typename... Parameters>
AbstractLieGroup<true>::AbstractLieGroup(const char* structure_constants, Parameters&&... parameters) :
	AbstractLieGroup{delegate_constructor_tag,structure_constants,collate(std::forward<Parameters>(parameters)...)} {}
template<typename... Parameters>
AbstractLieGroup<true>::AbstractLieGroup(const string& structure_constants, Parameters&&... parameters) :
	AbstractLieGroup{delegate_constructor_tag,structure_constants.c_str(),std::forward<Parameters>(parameters)...} {}


/** @brief A generic Lie group, whose structure constants depend on parameters
 */ 

class GenericLieGroup : public ConcreteManifold, public LieGroupWithParameters {
public:	
/** @brief Construct a generic Lie group, all of whose structure constants are parameters
 *  @param dimension The dimension of the group
 */
	GenericLieGroup(int dimension) : ConcreteManifold(dimension)
	{
		VectorSpace<DifferentialForm> twoforms=pForms(2);
		for (int i=1;i<=dimension;i++)
		{
			ex de;
			for (int j=1;j<=twoforms.Dimension();j++)
				de+=StructureConstant(N.a(i,j))*twoforms.e(j);
			Has_dTable::Declare_d(e(i),de);			
		}
	}
};

} /** @} */
#include "liegroupstructures.h"
#endif
