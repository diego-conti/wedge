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

#ifndef POLYBASIS_H
#define POLYBASIS_H

/** @ingroup LinearAlgebra */ 

/** @{ 
 * @file polybasis.h
 * @brief Bases of polynomials
 */


#include "wedge/base/wedgealgebraic.h"
#include "wedge/base/expressions.h"
#include "wedge/linearalgebra/lambda.h"
#include <type_traits>

namespace Wedge {
 using namespace  GiNaC;
 using namespace std;

/** @brief A basis of polynomials generating an ideal
 * 
 * This class is essentially a container with interface similar to std::vector. 
 * To apply Groebner bases methods, use the derived class PolyBasis_impl.
 * 
 * @remark This class represents a list of polynomials, but also the ideal they generate. Members which refer to the ideal start
 * with the word Ideal.
 */
template <typename Variable> class PolyBasis {
public:
	typedef exvector::value_type value_type; 			///<Implements container semantics, in analogy with std::vector
	typedef exvector::iterator iterator;				///<Implements container semantics, in analogy with std::vector
	typedef exvector::const_iterator const_iterator; 	///<Implements container semantics, in analogy with std::vector

	inline const_iterator begin() const {return polynomials.begin();} 	///<Implements container semantics, in analogy with std::vector 
	inline const_iterator end() const {return polynomials.end();} 		///<Implements container semantics, in analogy with std::vector

/** @brief Implements container semantics, in analogy with std::vector

@warning This member is meant to be used in conjunction with PolyBasis::insert.  Behaviour is undefined if non-const iterators
are used to modify the elements of PolyBasis directly.
*/
	inline iterator begin() {return polynomials.begin();} 		
/** @brief Implements container semantics, in analogy with std::vector

@warning This member is meant to be used in conjunction with PolyBasis::insert.  Behaviour is undefined if non-const iterators
are used to modify the elements of PolyBasis directly.
*/ 
	inline iterator end() {return polynomials.end();} 			
 
   /** @brief Create an ideal of polynomials from a range of generators
    * @param [gens_begin,gens_end) Polynomials with indeterminate type Variable
  */
	template<typename Iterator> PolyBasis(Iterator gens_begin, Iterator gens_end) : polynomials(gens_begin,gens_end) {
		SanityCheck(gens_begin,gens_end); 
		Update();
	} 
   /** @overload
  */
	PolyBasis() {}

	virtual ~PolyBasis() {}
   /** @brief Add a generator to the ideal
    * @param polynomial A generator to be appended to the list of generators
  */	
	void push_back(ex polynomial)
	{
		SanityCheck(&polynomial, &polynomial+1);
		polynomials.push_back(polynomial);		
		Update();
	}
   /** @brief Remove a generator from the ideal
    * @param iterator An iterator pointing to a generator
  */	
	
	iterator erase(iterator pos) {
		iterator result=polynomials.erase(pos);	
		Update();
		return result;
	}

  /** @brief Add some generators
   *  @param at An iterator at which the range is to be added
   *  @param [begin,end) A range of polynomials to be added
   */
	template<typename Iterator> void insert(iterator at, Iterator begin, Iterator end)
	{
		SanityCheck(begin,end);
		polynomials.insert(at,begin,end);		
		Update();
	}

  /** @brief Remove all generators
  */
	void clear() {
		polynomials.clear();
		Update();
	}

  /** @brief Test whether the basis is empty
   *  @return true if the basis contains no generators (either zero or non-zero)
   *  @sa IdealIsZero 
  */	
	bool empty() const
	{
		return polynomials.empty();
	}
	
	int size() const {return polynomials.size();}	///< Return the number of polynomials in the basis

  /** @brief Test whether this is the zero ideal
   *  @return true if all the elements of the basis are zero 
  */	
	bool IdealIsZero() const
	{
		const_iterator i=begin();
		while (i!=end())
			if (!i++->is_zero()) return false;
		return true;		
	}

	const exvector& Variables() const {return variables;}	///< Return the variables  appearing in the polynomials of this basis
//	const exvector& e() const {return polynomials;}
	
  /** @brief Return a variable by name    
  */	
	ex VariableByName(const char* name) const
	{
		for (typename exvector::const_iterator i=variables.begin();i!=variables.end();++i)
		{
			assert(is_a<Variable>(*i));
			if (ex_to<Variable>(*i).get_name()==name) return *i;
		}
		throw InvalidArgument(__FILE__,__LINE__);	
	}
	
	
  /** @brief Change the ordering of the variables appearing in the ideal
   *  @param x A variable to be taken to the front position
   *  
   * The ordering is relevant for Groebner basis reduction. 
  */	
	void TakeVariableToFront(ex x)
	{
		exvector::iterator i=find(variables.begin(),variables.end(),x);
		if (i==variables.end()) throw InvalidArgument(__FILE__,__LINE__,x);
		variables.erase(i);		
		variables.insert(variables.begin(),x);
	}

  /** @brief Change the ordering of the variables appearing in the ideal
   *  @param x A variable to be taken to the back position
   *  
   * The ordering is relevant for Groebner basis reduction. 
  */
	void TakeVariableToBack(ex x)
	{
		exvector::iterator i=find(variables.begin(),variables.end(),x);
		if (i==variables.end()) throw InvalidArgument(__FILE__,__LINE__,x);
		variables.erase(i);
		variables.push_back(x);		
	}
	
   /** @brief Eliminate variables by taking square-free factorization and solving linear equations
   *  
   * Returns a list of substitutions that have been applied
  */
	lst Eliminate();
protected:
	exvector variables;
	exvector polynomials;
	
/** @brief Invoked after an element is inserted or modified
 */
	virtual void Update() {
		variables.clear();
		GetSymbols<Variable>(variables,polynomials.begin(),polynomials.end());
	}
private:
	template<typename Iterator>
	static void SanityCheck(Iterator begin, Iterator end);	///< Test whether a range of polynomials is valid
};

/** @brief A basis of polynomials generating an ideal, enabling reduction by Groebner basis methods
 * 
 * The template parameter PolyAlgorithms determines the implementation to use; only CoCoA is supported at the moment.
 * 
 * @warning Use PolyBasis for polynomials for which the Groebner basis implementation is not applicable. Recall that CoCoA only works with rational numbers. 
 */
template<typename Variable, typename PolyAlgorithms> class PolyBasisImplementation : public PolyBasis<Variable>, public PolyAlgorithms::Initializer {
public:	
	/** @brief Construct a basis from a range of generators
	 * @param [gens_begin,gens_end) A range of polynomials
	 */
	template<typename Iterator> PolyBasisImplementation(Iterator gens_begin, Iterator gens_end) : PolyBasis<Variable>(gens_begin,gens_end) 
	{
		basis_is_reduced=false;
	}
	
	/** @overload
	 */
	PolyBasisImplementation() {basis_is_reduced=false;} 
	
	/** @brief Reduce the number of generators by Groebner basis methods
	 * @return A reference to this
	 */
	PolyBasisImplementation& Reduce()
	{
		if (!basis_is_reduced) {
			exvector reduced= PolyAlgorithms::template IdealReduce<Variable>(this->variables, this->polynomials.begin(),this->polynomials.end());
			this->polynomials=reduced;
			PolyBasis<Variable>::Update();
			basis_is_reduced=true;
		}
		return *this;
	}

/** @brief Test whether the ideal contains a polynomial
 * @param p A polynomial
 * @return true if p is in this ideal
 * 
 * @todo Consider defining a class Ideal, in analogy with VectorSpace.
 */
	bool IdealContains(ex p) const {		
		return PolyAlgorithms::template IdealContains<Variable>(this->variables, this->polynomials.begin(),this->polynomials.end(),p);
	}

/** @brief Test whether the radical of this ideal contains a polynomial
 * @param p A polynomial
 * @return true if some power of p is in this ideal
 */
	bool RadicalContains(ex p) const {		
		return PolyAlgorithms::template RadicalContains<Variable>(this->variables, this->polynomials.begin(),this->polynomials.end(),p);
	}

/** @brief Compute the intersection with another ideal */
	PolyBasisImplementation Intersected(const PolyBasisImplementation& I) const {
		ExVector intersection =  PolyAlgorithms::template IdealIntersection<Variable>(this->polynomials.begin(),this->polynomials.end(),I.polynomials.begin(),I.polynomials.end());
		return PolyBasisImplementation(intersection.begin(),intersection.end());
	}

/** @brief Compute a reduction of a polynomial modulo this ideal
 * @param p A polynomial
 * @return A polynomial representing p modulo this ideal
 * 
 * @todo Consider defining a class Ideal, in analogy with VectorSpace.
 */
	ex ReduceModuloIdeal(ex p) const {		
		
		return PolyAlgorithms::template ElementModuloIdeal<Variable>(this->variables, this->polynomials.begin(),this->polynomials.end(),p);
	}
	
/** @brief Test whether the ideal is the full ring, i.e. it contains one
 * @return true if the ideal contains one
 *	 
 * @todo It would be nice to have a function that says whether V(I) is empty, over R or Q...  
 */
	bool IdealIsOne() const {
		
		return PolyAlgorithms::template IdealIsOne<Variable>(this->variables, this->polynomials.begin(),this->polynomials.end());
	}
	
	void Update()
	{
		PolyBasis<Variable>::Update();
		MarkAsNotInitialized();			
	}
	
/** @brief Reduce the coefficients of an expression modulo this ideal
 * @param e A linear combination of objects of type T, with polynomials in Variable as coefficients
 * @return An expression obtained from e by reducing the coefficients modulo this ideal 
 */
	template<typename T> ex ReduceModuloIdeal(ex e) const
	{			
		e=e.expand();	
		internal::NormalFormHelper<T> v;
		e.accept(v);
		
		ex result;
		for (exmap::const_iterator i=v.coeffs.begin();i!=v.coeffs.end();i++)
		{	
			 result+=ReduceModuloIdeal(i->second)*i->first; 
		}
		return result;
	}
	
private:	
	mutable bool basis_is_reduced;
	
	void MarkAsNotInitialized()
	{		
		basis_is_reduced=false;
	}
};

	
template<typename Variable> template<typename Iterator> void PolyBasis<Variable>::SanityCheck(Iterator begin, Iterator end)
{
	for (;begin!=end;++begin)
		if (is_a<relational>(*begin)) throw WedgeException<invalid_argument>("Polynomial expected, but an equation was passed instead",__FILE__,__LINE__);
}

/** @brief Overloaded output operator
 */
template<typename Variable>  std::ostream& operator<<(std::ostream& out, const PolyBasis<Variable>& I)
{
	out<<I.Variables().size()<<" coordinates:"<<endl<<I.Variables()<<endl;
	out<<I.size()<<" polynomials:"<<endl;
	return out<<Range(I.begin(),I.end());		
}

namespace internal {
//helper class to single out the equation the form kx + p=0, with p not depending on x and of the lowest degree
struct EqnToSolve {
	exvector::iterator lowest;			//points to equation to solve
	int lowestdegree=numeric_limits<int>::max();	//degree of lowest
	ex var;						//variable to solve in, or zero if no variable found
	//find equation to solve, starting from first. If nonzero constant is found, return the 1 ideal.
	template<typename Parameter> EqnToSolve(PolyBasis<Parameter>& basis) {
		const exvector& variables{basis.Variables()};
		//we single out the equation of lowest degree that is degree one in at least one variable
		for (exvector::iterator i=basis.begin();i!=basis.end();++i)
		{
			int deg=Degree<Poly<Parameter> >(*i);
			if (deg==0 && !i->expand().is_zero()) {		//zeroes will be removed later
				LOG_DEBUG(*i);
				basis.clear(); basis.push_back(1);	//make this the 1 ideal
				var=0;
				return;
			}
			else if (deg<lowestdegree) {
				exvector::const_iterator x=variables.begin();
				//look for a variable x such that the equation has the form kx + p=0, with p not depending on x.
				while (x!=variables.end() && Degree<Poly<Parameter> >(*i - i->subs(*x==0))!=1) ++x;
				if (x!=variables.end()) {
					var=*x;
					lowest=i;
					lowestdegree=deg;
				}
			}
		}
	}
};

//simplify a polynomial equation using square free factorization
ex SimplifyPolyEqn(ex x,lst variables=lst());

}

template<typename Parameter> lst PolyBasis<Parameter>::Eliminate() {
	static_assert(is_base_of<symbol,Parameter>::value,"PolyBasis::Eliminate requires that Parameter derive from symbol");
	exvector eliminated_vars;
	exvector replaced_with;
	//the first thing to do would be taking the radical of the zero ideal, but I can't see how to get CoCoALib to do that.

	for (auto& x: *this) x=internal::SimplifyPolyEqn(x.expand(),lst{variables.begin(),variables.end()});
	while (true) {
		internal::EqnToSolve eqn{*this}; 
		
		if (!eqn.var.is_zero()) {
			ex sol=lsolve(*eqn.lowest==0,eqn.var);
			assert(sol!=lst());
			*eqn.lowest=0;
			eliminated_vars.push_back(eqn.var);
			replaced_with.push_back(sol);

			ex sub=(eqn.var==sol);
			for (auto& eq: replaced_with) {
				eq= eq.subs(sub);	//eliminate var in subtitution
			}
			for (auto& x : polynomials)
				if (!x.is_zero())	//eliminate var in equation
				try {
					 x=internal::SimplifyPolyEqn(x,lst{variables.begin(),variables.end()});
				}
				catch (...)
				{
					LOG_ERROR(x);
					LOG_ERROR(sub);
					LOG_ERROR(x.subs(sub));
					LOG_ERROR(variables);
					throw;
				}
			Update();
		}
		else break;
	}
	//remove zeroes
	iterator i=polynomials.begin();
	while (i!=polynomials.end()) {
		if (i->is_zero()) i=polynomials.erase(i);
		else ++i;
	}
	lst result;
	assert(eliminated_vars.size()==replaced_with.size());
	for (int i=0;i<eliminated_vars.size();++i)
		result.append(eliminated_vars[i]==replaced_with[i]);
	return result;
}



}
#endif
