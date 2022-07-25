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
#ifndef BASIS_H_
#define BASIS_H_

/** @ingroup LinearAlgebra */ 

/** @{ 
 * @file basis.h
 * @brief Vector spaces as bases of independent vectors
 */
 
#include "linearcombinations.h"

namespace Wedge {
using namespace GiNaC;

/** @brief Exception thrown by IBasis::Components, VSpace::Components etc. when the element is not in the vector space
*/
class NotInSpan: public WedgeException<runtime_error> {
public:
	NotInSpan(const char* in_file, int at_line, ex v) : WedgeException<runtime_error>("Not in vector space: "+ToString(v),in_file,at_line) {}
};


/** @brief Abstract class for a basis of generators for some vector space.
 * 
 * In Wedge, a vector space is defined by a basis of generators. An object of type derived from Basis<T> is a containter
 * with semantics similar to std::vector, and it represents a list of independent vectors
 * in the vector space corresponding to the type T. 
 * 
 * Elements can be added to an existing IBasis with the insert(iterator, InputIterator, InputIterator) function; in this case, 
 * a minimal set of linearly independent elements is chosen so as to preserve the filtration 
 * \f$\{ Span\{e_1,\dots,e_n\}, 1\leq n\leq size()\}\f$.
 * 
 * IBasis also provides a function Components() to determine the components of a vector in terms of the generators. 
 * This function, as well as dual(), requires computing the inverse of a  matrix, and is therefore only suitable for low dimensions.
 * This computation is referred to as ``initialization'' of the basis, and is only performed when needed.
 * 
 * @remark An IBasis is essentially an improved ExVector. In fact, if no elements are added after construction, IBasis is equivalent to Wedge::ExVector,
 * with a minimal overhead.   
 */

template<typename T> class IBasis : public LinearCombinations  {	
public:	
/** @brief Construct an empty basis
*/
	IBasis() {}
/** @brief Construct a basis from a range of linearly independent elements
 *  @param [from,to) A range of elements to add   
 * 
 *   For efficiency, no check is made that the elements are linearly independent
 */
	template<typename InputIterator>  IBasis(InputIterator from, InputIterator to) : LinearCombinations(from,to)  {}
/** @overload
*/
	IBasis(const exvector& basis) : LinearCombinations(basis) {}

/** @brief Compute the components of a vector with respect to the basis
 *  @param v A vector (in a vector space of type T)
 *  @param lies_in_span (out) If specified, the referenced bool indicates whether v lies in the vector space spanned by the basis  
 *  @return The components of the vector in terms of the basis
 * 
 *  This function requires initializazion, and is therefore slow for high dimensions. 
 *
 *  In other words, when this function is called the first time, Basis computes the components of all 
 *  simple vectors of type T.  
 * 
 * @exception NotInSpan Thrown if v does not lie in the space spanned by this basis and lies_in_span is NULL
 */
	ExVector Components(ex v, bool* lies_in_span=NULL) const
	{
		v=v.expand();
		if (NotInitialized()) SetBasis();
		map<ex,ExVector,ex_is_less>::const_iterator it=inverse.find(v);		
		if (it!=inverse.end()) {
			ExVector result= it->second;
			for (int i=size();i<result.size();i++)
				if (result[i]!=0)	//v is not in the vector space spanned by this basis!
					throw NotInSpan(__FILE__,__LINE__,v);
			result.resize(size());
			return result;			
		}
		else {	
			ExVector result(size());
			internal::NormalFormHelper<T> helper;
			v.accept(helper);
			for (it=inverse.begin();it!=inverse.end();it++)
			{
				exmap::const_iterator k=helper.coeffs.find(it->first);
				if (k!=helper.coeffs.end()) 
					for (int i=0;i<size();++i)					
						result[i]+=k->second*it->second[i];
			}

			ex test;
			for (int i=0;i<size();i++)
			{
				result[i]=result[i].expand();
				test+=result[i]*e[i];
			}			
			//check whether v=test
			list<ex> eqns;
			GetCoefficients<T>(eqns,v-test);
			
			if (lies_in_span==NULL) {
				for (list<ex>::iterator i=eqns.begin();i!=eqns.end();i++)
					if (!i->numer().expand().is_zero()) {
						LOG_WARN(v);
						LOG_WARN(test);
						LOG_WARN((v-test).expand());
						LOG_WARN(ExVector(*this));
						LOG_WARN(size_);
						throw InvalidArgument(__FILE__,__LINE__,v);
					}
			}
			else {
				*lies_in_span=true;
				for (list<ex>::iterator i=eqns.begin();i!=eqns.end();i++)
					if (!i->normal().is_zero()) {					
						*lies_in_span=false;
					}
			}
			return result;
		}
	}
	
/** @brief Like Components(ex), but computes the components of a vector with respect to the enlarged basis
 *  @param v A vector in the span of the enlarged basis
 *  @return an ExVector containing the components 
 * 
 *  The enlarged basis is normally only used internally; it is obtained from the basis
 *  by appending simple elements until every simple element appearing in any generator
 *  lies in the span of this enlarged basis.
 * 
 * It is assumed that v is in the span of the enlarged basis. No exception is thrown if this condition does not hold. 
 */
	ExVector AllComponents(ex v) const
	{
		v=v.expand();
		if (NotInitialized()) SetBasis();
		map<ex,ExVector,ex_is_less>::const_iterator it=inverse.find(v);
		if (it!=inverse.end()) {
			return  it->second;
		}
		else {
			ExVector result(e.size());
			for (it=inverse.begin();it!=inverse.end();it++)
			{
				ex coeff=TrivialPairing<T>(it->first,v);	//this works because it->first is simple
				if (coeff!=0)  {
					for (int i=0;i<e.size();i++)
						result[i]+=coeff*it->second[i];						
				}
			}
			return result;
		}				
	}	
	ex e_dual(OneBased i) const {if (NotInitialized()) SetBasis(); return dual_(i);}
	exvector::const_iterator dual_begin() const {if (NotInitialized()) SetBasis(); return dual_.begin();}
	exvector::const_iterator dual_end() const {return dual_.end();}
	
/** @brief Compute (if needed) and return the dual basis  
 *  @return an ExVector containing the dual basis
 *
 * Recall that vector spaces are consistently identified with their duals by the existence of a standard
 * basis, consisting of all simple elements. So, the dual basis is still a vector of elements of type T 
 *
 *  This function requires initializazion, and is therefore slow for high dimensions.
 *
 * @note The dual vector may be longer than size() (if the object is of type SubBasis)
 */
	const ExVector& dual() const {if (NotInitialized()) SetBasis(); return dual_;}

/** @brief Check whether this basis has been initialized
 *  @return true if this basis requires initializing
 * 
 *  The basis initializes itself when needed; the caller should not worry about that, except
 *  that initialization might take a long time.
 */
	bool NotInitialized() const {return dual_.empty() && !e.empty();}

	void push_back(ex v)
	{
	 	// if basis is initialized, then it is optimal to call Components(); otherwise, better call insert()	 	
	    	if (!NotInitialized())
    		{
    			bool contained;
    			static_cast<void>(Components(v,&contained));
    			if (contained) return;
	    	}
	    	insert(end(),&v,(&v)+1);
	}

protected:
/** @brief Mark this basis as not initialized
 * 
 *  Call this function after modifying the elements of this->e directly. 
 */
	void MarkNotInitialized() {
		inverse.clear();
		dual_.clear();
	}

/** @brief Initialize this basis
 *  
 * It is assumed that the elements of this->e are linearly independent when this function is called
 */ 
	virtual void SetBasis()	const=0;

/** @brief Associates to a simple vector its components in this basis.
 * 
 *  Computed when needed.
 */
	mutable std::map<ex,ExVector,ex_is_less> inverse;

/** @brief The dual basis as an ExVector.
 * 
 *  Computed when needed.
 */	
	mutable ExVector dual_;									
};


/** @brief Implements IBasis<T> using the specified linear algebra algorithms
 * @param T An algebraic class, representing simple elements of a vector space
 * @param LinAlgAlgorithms A type implementing the linear algebra algorithms that Basis is to use
 */
template<typename T, typename LinAlgAlgorithms=DefaultLinAlgAlgorithms> class Basis  : public LinAlgAlgorithms, public IBasis<T> {
public:
/** @brief Construct an empty basis
 */
	Basis() {}

/** @brief Construct a basis from an exvector of linearly independent elements
 *  @param basis a vector of linearly independent elements
 * 
 *  For efficiency, no check is made that the elements are linearly independent
 */	
	Basis(const exvector& basis) : IBasis<T>(basis) {}

/** @brief Construct a basis from a range of linearly independent elements
 *  @param [from,to) a range of linearly independent elements
 * 
 *  For efficiency, no check is made that the elements are linearly independent
 */	
	template<typename InputIterator>  Basis(InputIterator from, InputIterator to) : IBasis<T>(from,to) 
	{}

protected:
	void SetBasis()	const //assume prune() has been called
	{		
		exvector symbols;
		GetSimple<T>(symbols,this->e.begin(),this->e.end());
		int oldsize=this->e.size();
		if (symbols.size()>this->e.size())	//not a basis of the full space,
		{										//so complete it
			this->e.insert(this->e.end(),symbols.begin(),symbols.end());
			IBasis<T>::CheckIndependenceAndPrune();
		}
		assert(symbols.size()==this->e.size());
		int dimension=symbols.size();
		
		exvector symbols_not_in_basis;
		symbols_not_in_basis.reserve(oldsize);
		exvector::const_iterator i=this->begin()+oldsize, j=symbols.begin();
		while (j!=symbols.end())
		{
			if (i!=this->e.end() && *i==*j)
				i++,j++;
			else 
			{	
				symbols_not_in_basis.push_back(*j);
				j++;
			}			
		}
		assert(symbols_not_in_basis.size()==oldsize);
		
	//setup the inverse vector, containing the components of ''simple'' elements w.r.t basis
			
//write e_i=a_{ij}symbols_not_in_basis_j + a'_{ir}e_r with i,j in [0,oldsize), r in [oldsize,dimension)
//We must solve symbols_not_in_basis_k=b_{kl}e_l + b'_{ks}e_s, with k,l in [0,oldsize), s in [oldsize,dimension) 
//this means that e_i =a_{ij}b_{jl}e_l + a'_{ir}e_r + a_{ij}b'_{js}e_s, and in particular b is the inverse of a

		typename LinAlgAlgorithms::InverseMatrix a(oldsize,oldsize);
		for (int i=0;i<oldsize;i++)
			for (int j=0;j<oldsize;j++)
				a(i,j)=TrivialPairing<T>(symbols_not_in_basis[j],(*this)[i]);
		
		LOG_DEBUG(*this);
		LOG_DEBUG(symbols);
		LOG_DEBUG(a);
		matrix b;
		try {
			b=a.inverse();	//destroys a
		}
		catch (std::runtime_error)
		{
			LOG_ERROR(symbols);
			LOG_ERROR(*this);
			ex e=symbols[0]-(*this)[0];
			e.print(print_tree(cout));
			LOG_ERROR(a);
			assert(false);	
		}		
		
//the dual basis must of course satisfy e^k(e_i)=\delta_{ik}, i.e.
//a_{ij}e^k(symbols_not_in_basis_j) + a'_{ir} \delta_{kr} = \delta_{ik}
//so e^k=b_{lk}symbols_not_in_basis_l when k in [0,oldsize)

		this->dual_.clear();
		this->dual_.reserve(oldsize);
		for (int k=0;k<oldsize;k++)
		{
			ex e=0;
			for (int l=0;l<oldsize;l++)
				e+=b(l,k)*symbols_not_in_basis[l];
			this->dual_.push_back(e.expand());
		}
		
		this->inverse.clear();
		for (int i=0;i<oldsize;i++)
		{	
			ex e=0;
			this->inverse[symbols_not_in_basis[i]]=ExVector(dimension);
			for (int j=0;j<oldsize;j++)
			{
				this->inverse[symbols_not_in_basis[i]][j]=b(i,j);
				e+=b(i,j)* (*this)[j];	
			}
			for (int j=oldsize;j<dimension;j++)
				this->inverse[symbols_not_in_basis[i]][j]=-TrivialPairing<T>(e,(*this)[j]);
		}
		for (int i=oldsize;i<dimension;i++)
		{
			ExVector v(dimension);
			v[i]=1;
			this->inverse[this->e[i]]=v;			
		}		
	}
private:
	int ConstPrune() const {
		exvector symbols;
		GetSimple<T>(symbols,this->e.begin(),this->e.end());
		typename LinAlgAlgorithms::IndependenceMatrix m(this->e.size(),symbols.size());
		for (int i=0;i<this->e.size();++i)
		{
			internal::NormalFormHelper<T> v;
			this->e[i].accept(v);
			for (int j=0;j<symbols.size();++j)
			{
				exmap::const_iterator k=v.coeffs.find(symbols[j]);
				if (k!=v.coeffs.end()) 
					m.M(i,j)=k->second;
			}
		}
				
		m.ChooseLinearlyIndependentRows();
		int new_size=0;
		int j=0;
		for (typename LinAlgAlgorithms::IndependenceMatrix::const_iterator i=m.IndependentRowsBegin();i!=m.IndependentRowsEnd();i++,j++)
		{
			LOG_DEBUG(*i);
			if (*i<this->size_) ++new_size;
			if (*i!=j) this->e[j]=this->e[*i];			
		}
		LOG_DEBUG(*this);		
		this->e.resize(j);
		return new_size;		
	}
	
};


/** @brief %Basis of a subspace of some bigger space
 * @param T An algebraic class, representing simple elements of a vector space
 * @param LinAlgAlgorithms A type implementing the linear algebra algorithms that SubBasis is to use \sa Basis
 * 
 * A SubBasis object is a basis of some vector space reflecting the direct sum \f$ Space = Subspace \oplus Complement\f$    
 */
template<typename T, typename LinAlgAlgorithms=DefaultLinAlgAlgorithms> class SubBasis : public Basis<T, LinAlgAlgorithms> {
public:
/** @brief Construct an empty basis
 */
	SubBasis()
	{
		comprehensive_dimension=0;
	}

/** @brief Construct a sub-basis object from a basis of generators and a basis of the complement
 *  @param [basis_from,basis_to) A basis of the subspace
 *  @param [not_in_basis_from,not_in_basis_to) A basis of the complement
 * 
 * Defines a sub-basis object, corresponding to a basis of \f$Span [from,to)\f$ as a subspace of 
 * \f[Span[basis\_from,basis\_to)\oplus Span[not\_in\_basis\_from,not\_in\_basis\_to).\f]
 * 
 * For efficiency, no check is made that the elements are linearly independent and the sum is direct. 
 */
	template<typename InputIterator1, typename InputIterator2>  
	SubBasis(InputIterator1 basis_from, InputIterator1 basis_to,
		InputIterator2 not_in_basis_from, InputIterator2 not_in_basis_to)	  
	{
		while (basis_from!=basis_to)
			this->e.push_back(basis_from++->expand());
		this->size_=this->e.size();
		while (not_in_basis_from!=not_in_basis_to)
			this->e.push_back(not_in_basis_from++->expand());
		comprehensive_dimension=this->e.size();
		this->MarkNotInitialized();
	}

/** @brief Add a range of elements to current basis (of the subspace)
 *  @param at The position at which the range is to be inserted (e.g. this->end() to append to the existing generators)
 *  @param [from,to) A range of elements to add
 *  
 *  The elements in the range are not required to be linearly independent. After adding the elements,
 *  this function clears all redundant elements, in a flag-preserving way. 
 *  
 * \sa Basis 
 */
	template<typename InputIterator> void insert(typename IBasis<T>::iterator at, InputIterator from, InputIterator to)
	{
		assert(at<=this->end());
		list<ex> expanded;
		while (from!=to)
		{
			expanded.push_back(from++->expand());
			this->size_++;
			comprehensive_dimension++;	
		}
		this->e.insert(at,expanded.begin(),expanded.end());
		Prune();	//remove redundant generators from basis
		this->MarkNotInitialized();
	}

/** @brief Return the dimension of the containing space 
 */
	int DimensionOfContainingSpace() const {return comprehensive_dimension;}
	 
/** @brief Add elements to the the complement until the space contains the given range of generators
 * @param [from,to) Elements to be added to the complement
 * 
 * Assumes the current basis is in fact a basis (i.e. prune() has been called); this condition is
 * always satisfied unless the elements of the subbasis have been modified directly.
 * 
 * The complement is enlarged so that
 * \f$ Subspace \oplus Complement \ni \{from,\dots, to)\}\f$
*/
	template<typename InputIterator> void AddToComplement(InputIterator from, InputIterator to)
	{
		//this->_size=this->e.size();
		//
		//prune();	//remove redundant generators from basis
		//if(this->_size!=this->e.size())
		//			throw WedgeException<std::runtime_error>("Not a basis", "Basis");
		exvector expanded;
		while (from!=to)
			expanded.push_back(from++->expand());
		this->e.insert(this->begin()+DimensionOfContainingSpace(),expanded.begin(),expanded.end());
		comprehensive_dimension +=expanded.size();
		Prune();
		this->MarkNotInitialized();
	}

	typename IBasis<T>::const_iterator inline complement_begin() const {return this->end();}///< Returns an iterator pointing to the first element of the basis of the complement
	typename IBasis<T>::const_iterator inline complement_end() const {return this->begin()+DimensionOfContainingSpace();}///< Returns an iterator pointing after the last element of the basis of the complement
	typename IBasis<T>::const_iterator inline full_begin() const {return this->begin();} ///< Returns an iterator pointing to the first element of the basis of the containing space
	typename IBasis<T>::const_iterator inline full_end() const {return this->begin()+DimensionOfContainingSpace();} ///< Returns an iterator pointing after the last element of the basis of the containing space

private:	
	void Prune() {
		exvector symbols;
		GetSimple<T>(symbols,this->e.begin(),this->e.end());
		typename LinAlgAlgorithms::IndependenceMatrix m(this->e.size(),symbols.size());
		for (int i=0;i<this->e.size();++i)
		{
			internal::NormalFormHelper<T> v;
			this->e[i].accept(v);
			for (int j=0;j<symbols.size();++j)
			{
				exmap::const_iterator k=v.coeffs.find(symbols[j]);
				if (k!=v.coeffs.end()) 
					m.M(i,j)=k->second;
			}
		}
				
		m.ChooseLinearlyIndependentRows();
		int new_size=0;
		int new_comprehensive_dimension=0;
		int j=0;
		for (typename LinAlgAlgorithms::IndependenceMatrix::const_iterator i=m.IndependentRowsBegin();i!=m.IndependentRowsEnd();i++,j++)
		{
			LOG_DEBUG(*i);
			if (*i<this->size_) ++new_size;
			if (*i<comprehensive_dimension) ++new_comprehensive_dimension;
			if (*i!=j) this->e[j]=this->e[*i];			
		}
		LOG_DEBUG(*this);
		this->size_=new_size;
		comprehensive_dimension=new_comprehensive_dimension;
		this->e.resize(j);
	}
	
	int comprehensive_dimension;	///< Dimension of the containing space	
};



/** @brief Overloaded output operator */
template<typename T, typename LinAlgAlgorithms> inline std::ostream& operator<<(std::ostream& os, const SubBasis<T,LinAlgAlgorithms>& basis)
{
	os<<"Subbasis:"<<endl;
	os<<Range(basis.begin(),basis.end());
	os<<"Complement:"<<endl;
	return os<<Range(basis.complement_begin(),basis.complement_end());
}


} /** @} */
#endif /*BASIS_H_*/
