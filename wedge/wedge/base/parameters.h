/*******************************************************************************
 *  Copyright (C) 2007-2022 by Diego Conti, diego.conti@unipi.it 
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
#ifndef PARAMETERS_H_
#define PARAMETERS_H_
#include "wedgebase.h"
#include "../linearalgebra/ginaclinalg.h"

/** @ingroup Base
 *  @{
 * 
 *  @file parameters.h
 *  @brief Objects depending on parameters
*/
 


namespace Wedge {

namespace internal {
	/** @brief Helper class representing a subset of a set of ex
	*/
	class Subset : public lst
	{
		const set<ex,ex_is_less>& superset;
	public:
		Subset(const set<ex,ex_is_less>& contains_this) : superset(contains_this) {}

		typedef ex value_type;
		void push_back(ex x) {
			if (superset.find(x)!=superset.end()) lst::append(x);
		}
	};

	template<typename Iterator> Subset& Insert(Subset& container, Iterator begin, Iterator end)
	{
		while (begin!=end) container.push_back(*begin++);
		return container;
	}
}


/** @brief Abstract base class for mathematical objects that depend on a list of parameters
 * 
 * Deriving from this class permits to distinguish among different sets of parameters corresponding to the same type
 * - This class stores internally a list of the relevant parameters. 
 * - The parameters must either be created by the Parameter member function, or "registered" by invoking the member function StoreParameters
 * \sa HasParameters
 */
template <typename TypeOfParameter> class HasParameterList {
public:
/** @brief Eliminate some parameters by imposing linear conditions on an expression
 * @param alpha An expression depending linearly on the parameters
 * @exception InconsistentDeclaration Thrown if alpha cannot be zero for any choice of the parameters
 * @exception std::invalid_argument Thrown by %GiNaC framework if the dependence on the parameters is not linear
 * 
 * This function solves \f$\alpha=0\f$ and updates the parameters accordingly.
 * 
 * @sa PolyBasis, PolyBasis_impl 
 */
	void DeclareZero(ex alpha)
	{
		DeclareZero(&alpha, (&alpha)+1);
	}
	
/** @overload
 */
	template<typename Iterator> void DeclareZero(Iterator begin, Iterator end)
	{		
		list<ex> eqns;
		for (Iterator i=begin;i!=end;++i)		
			GetCoefficientsComplex(eqns,*i,withRHS);
		if (eqns.empty()) return;		
		internal::Subset unknowns(parameters);
		GetSymbols<TypeOfParameter>(unknowns,eqns.begin(),eqns.end());
		LOG_DEBUG(unknowns);
		LOG_DEBUG(eqns);
		ex sol=DefaultLinAlgAlgorithms::lsolve(lst(eqns),unknowns);
		LOG_DEBUG(sol);
		if (sol==lst()) {
			LOG_ERROR(Range(begin,end));
			LOG_ERROR(eqns);			
			throw InconsistentDeclaration(__FILE__,__LINE__,"parameters");
		}
		assert(is_a<lst>(sol));
		DeclareConditions(ex_to<lst> (sol));
	}
protected:
/** @brief Declare that all the symbols of type TypeOfParameter are to be viewed as parameters
 * @param [from,to] A range of ex
 */
	template<typename Iterator> void StoreParameters(Iterator from, Iterator to)
	{
		GetSymbols<TypeOfParameter>(parameters,from,to);			
	}	
/** @brief Construct a parameter
 *
 * Subclasses should either construct the parameter objects directly and add them with StoreParameters, or use those returned by this function.
*/
	ex Parameter() {
		//parameters.append(TypeOfParameter());
		//return parameters.back();
		return *parameters.insert(TypeOfParameter()).first;
	}
/** @overload
*/
	template<typename T> ex Parameter(T t) {
		//parameters.append(TypeOfParameter(t));
		//return parameters.back();
		return *parameters.insert(TypeOfParameter(t)).first;
	}
private:
/** @brief Declare conditions in the parameters
 *  @param list_of_equations Equations of the form symbol==value
 * 
 *  The list of equations has the form of the return value of a call to GiNaC::lsolve 
 */	
	virtual void DeclareConditions(const lst& list_of_equations)=0;

	set<ex,ex_is_less> parameters;
};

/** @brief Abstract base class for mathematical objects that depend on parameters
 * 
 * This class treats all objects of type TypeOfParameter as parameters. \sa HasParameterList
 */
template <typename TypeOfParameter> class HasParameters {	
public:
/** @brief Eliminate some parameters by imposing linear conditions on an expression
 * @param alpha An expression depending linearly on the parameters
 * @exception InconsistentDeclaration Thrown if alpha cannot be zero for any choice of the parameters
 * @exception std::invalid_argument Thrown by %GiNaC framework if the dependence on the parameters is not linear
 * 
 * This function solves \f$\alpha=0\f$ and updates the parameters accordingly.
 * 
 * @sa PolyBasis, PolyBasis_impl 
 */
	void DeclareZero(ex alpha)
	{
		DeclareZero(&alpha, (&alpha)+1);
	}
	
/** @overload
 */
	template<typename Iterator> void DeclareZero(Iterator begin, Iterator end)
	{		
		list<ex> eqns;
		for (Iterator i=begin;i!=end;++i)		
			GetCoefficientsComplex(eqns,*i,withRHS);
		if (eqns.empty()) return;		
		list<ex> unknowns;
		GetSymbols<TypeOfParameter>(unknowns,eqns.begin(),eqns.end());
		LOG_DEBUG(unknowns);
		LOG_DEBUG(eqns);
		ex sol=DefaultLinAlgAlgorithms::lsolve(lst(eqns),lst(unknowns));
		if (sol==lst()) {
			LOG_ERROR(Range(begin,end));
			LOG_ERROR(eqns);			
			throw InconsistentDeclaration(__FILE__,__LINE__,"parameters");
		}
		assert(is_a<lst>(sol));
		DeclareConditions(ex_to<lst> (sol));
	}

	template<typename PolyContainer> PolyContainer DeclareZero(PolyContainer&& equations)
	{
		DeclareConditions(equations.Eliminate());
		return equations;
	}
private:
/** @brief Declare conditions in the parameters
 *  @param list_of_equations Equations of the form symbol==value
 * 
 *  The list of equations has the form of the return value of a call to GiNaC::lsolve 
 */	
	virtual void DeclareConditions(const lst& list_of_equations)=0;
};

} /** @} */
#endif /*PARAMETERS_H_*/
