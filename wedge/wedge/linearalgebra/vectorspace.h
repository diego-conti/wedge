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
#ifndef WEDGEVECTORSPACE_H
#define WEDGEVECTORSPACE_H
/** @ingroup LinearAlgebra */ 

/** @{ 
 * @file vectorspace.h
 * @brief Vector spaces, affine spaces, subspaces
 */
#include "wedge/base/wedgealgebraic.h"
#include "wedge/linearalgebra/linear.h"
#include "wedge/linearalgebra/basis.h"
#include <ostream>

namespace Wedge {

template<typename T> class Subspace;

/** @brief Exception thrown by members of VSpace when the set of affine equations is inconsistent
*/
class EmptyAffineSpace: public WedgeException<runtime_error> {
public:
	EmptyAffineSpace(const char* in_file, int at_line) : WedgeException<runtime_error>("Affine space is empty",in_file,at_line) {}
};
	

// ****************************************************************************
// *							VSpace										  *
// ****************************************************************************
/** @brief Abstract base class for vector spaces spanned by elements of type T
 *  @param T By definition, the type of "simple" elements
 *
 * A vector space is represented by a basis of generators, which in turn are linear combinations 
 * of elements of a certain type T, derived from GiNaC::basic. Distinct elements of type T are assumed to be independent.
 * 
 * Standard coordinates are given on each vector space, since the basis is part of the definition. These coordinates
 * can be used to identify linear or affine subspaces.
 * 
 * Computationally intensive methods depend on the template parameter LinAlgAlgorithms.
 */
template<typename T> class VSpace {
public:	
/** @brief A coordinate on this vector space
* 
*  Mathematically, coordinates on a vectors space are elements of its dual. The generic element of a vector space
*  spanned by \f$v_1, \dots, v_n\f$ is a linear combination \f$\sum_i a_i v_i\f$, where the \f$v_i\f$ 
*  are of type V and the \f$a_i\f$ are of type VSpace<V>::Coordinate.
*
* @note Coordinates are complex by default.
*/ 
	class Coordinate : public Register<Coordinate,symbol>::Algebraic
	{
	public:
		static const char* static_class_name() {return "Coordinate";}
	 /** @brief Constructor
	  *  @param name A name to be given to this Coordinate (for output only).
	  *  @param tex_name \f$\text{\TeX}\f$ version of the name of this coordinate
	  */
		Coordinate(const Name& name) : Coordinate::RegClass(name) 	{}
		Coordinate() {}
	};

/** @brief Create a trivial vector space (i.e., the set {0}). 
 * @param coordinate_name The letter or string associated to coordinates on this vector space (for output only)
 */	
	VSpace<T>(const Name& coordinate=N.lambda) : coordinate_name(coordinate.PlaceHolder()) {
	}

	virtual ~VSpace() {}

////////////////////////////////////////////////////////////
// 			Basic linear algebra
////////////////////////////////////////////////////////////
/**  
 * @return The dimension of this vector space
 */
	int Dimension() const {return e().size();}

/** @brief Returns the generic element of this vector space 
 *  @returns  The generic element of this vector space as a linear combination \f$\sum_i a_i e_i\f$, where the \f$e_i\f$
 *  are generators and the \f$a_i\f$ are coordinates
 */
	ex GenericElement() const;

/** @brief Check whether this space contains a given vector
 *  @param v A linear combination of elements of type T
 *  @return true if v lies in this space
 */
	bool Contains(ex v) const
	{
		if (e().NotInitialized())
		{
			Basis<T> b(e());
			b.insert(b.end(),&v, (&v)+1);
			int delta=b.size()-e().size();
			if (delta!=0 && delta!=1)
			{
				LOG_ERROR(delta);
				LOG_ERROR(e());
				LOG_ERROR(b);
				LOG_ERROR(v);
			}
			return delta==0;	
		}
		else  
		{
			bool lies_in_span;
			static_cast<void> (e().Components(v,&lies_in_span));
			return lies_in_span;
		}
	}


/** @brief Compute the components of a vector with respect to the basis
 *  @param v A vector in this vector space
 * @param options Specifies the options to be passed to GiNaC::lsolve
 *  @return The components of the vector in terms of the basis
 * @exception NotInSpan Thrown if v does not lie in this vector space
 *
 * @note One can alternatively call e().Components(), which behaves differently. Performance may vary.
 */
	ExVector Components(ex v,unsigned options = solve_algo::automatic) const
	{
		if (e().NotInitialized())
		{
			lst eqns;
			ex c=lst(coordinates.begin(),coordinates.end());
			
			GetCoefficients<T>(eqns,GenericElement()-v,withRHS);
			ex sol=lsolve(eqns,c,options);
			if (sol==lst()) throw NotInSpan(__FILE__,__LINE__,NormalForm<T>(v));

			c=c.subs(sol).expand();
			return ExVector(c.begin(),c.end());
		}
		else return e().Components(v);
	}
/** @brief Check whether two VSpace objects represent the same vector space
 *  @param V A VSpace<T> object to be compared with *this
 *  @return true if the vector spaces have the same dimension and one contains a basis of the other, false ohterwise
 */
	bool operator==(const VSpace& V)
	{
		if (Dimension()!=V.Dimension()) return false;
		if (e().NotInitialized()) {
			for (int i=1;i<=Dimension();i++)
				if (!V.Contains(e(i))) return false;
		}
		else { 
			for (int i=1;i<=Dimension();i++)
				if (!Contains(V.e(i))) return false;
		}
		return true;
	}

/////////////////////////////////////////////////////////////////////
// 		Functions that determine a subspace
/////////////////////////////////////////////////////////////////////

/** @brief Compute the equations corresponding to the orthogonal complement of a given subspace 
 * @param container A container where the equations are to be stored
 *  @param subspace A subspace of this vector space
 *  @param metric An object defining the metric to use    
 * @returns A reference to container
 */
	template<typename Container, typename Metric> Container& GetOrthogonalEquations(Container& container, const VSpace<T>& subspace, const Metric& metric) const
	{
		return GetOrthogonalEquations(container,subspace.e_begin(),subspace.e_end(),metric);
	}
/** @overload
*/
	template<typename Container, typename Iterator, typename Metric> Container& GetOrthogonalEquations(Container& container, Iterator subbasis_begin, Iterator subbasis_end, const Metric& metric) const;

/** @brief Compute a basis of solutions for the given linear or affine equations on this vector space
 * @param container The container where the basis is to be stored
 * @param [first,last) A range of affine equations in the coordinates of this vector space
 * @param vector (out) If specified, will contain an element of the affine space.
 * @returns A reference to container
 * @throw EmptyAffineSpace is thrown if there is no solution
 * @throw WedgeException<std::runtime_error> is thrown if vector is not specified and the equations are not linear
*/	
	template<typename LinAlgAlgorithms,typename Container,typename Iterator> Container& GetSolutions(Container& container, Iterator first, Iterator last, ex* vector=NULL) const;

/** @overload
*/
	template<typename Container,typename Iterator> Container& GetSolutions(Container& container, Iterator first, Iterator last, ex* vector=NULL) const
	{
		return GetSolutions<DefaultLinAlgAlgorithms>(container,first,last,vector);
	}


/** @brief Split this space into a given subspace and its ortogonal
  *  @param [subbasis_begin, subbasis_end) A range of generators whose span is the subspace to return
  *  @param metric The metric to be used to compute the (orthogonal) complement
  * @returns A Subspace object
 */
	template<typename LinAlgAlgorithms,typename Iterator, typename Metric>
		Wedge::Subspace<T> Subspace(Iterator subbasis_begin, Iterator subbasis_end, const Metric& metric) const;

/** @overload
 */
	template<typename Iterator, typename Metric>
		Wedge::Subspace<T> Subspace(Iterator subbasis_begin, Iterator subbasis_end, const Metric& metric) const
	{
		return Subspace<DefaultLinAlgAlgorithms>(subbasis_begin,subbasis_end,metric);
	}

/** @brief Construct a Subspace object
*  @param [subbasis_begin, subbasis_end) A range of generators whose span is the subspace to return
* @returns A Subspace object; the complement is chosen arbitrarily
* 
* Much faster than Subspace(Iterator, Iterator, const Metric&)
*/
   template<typename Iterator> 
	Wedge::Subspace<T> Subspace(Iterator subbasis_begin, Iterator subbasis_end) const;

/** @brief Construct a Subspace object from a range of linear equations
 *  @param [eqns_begin,eqns_end) A range of equations defining the subspace
 *  @param metric The metric to be used to compute the (orthogonal) complementt
 * @returns A Subspace object
 */
   template<typename LinAlgAlgorithms,typename InputIterator,typename Metric>
   		Wedge::Subspace<T> SubspaceFromEquations(InputIterator eqns_begin,InputIterator eqns_end, const Metric& metric) const;

/** @overload
 */
   template<typename InputIterator,typename Metric>
   		Wedge::Subspace<T> SubspaceFromEquations(InputIterator eqns_begin,InputIterator eqns_end, const Metric& metric) const
   	{
   		return 	SubspaceFromEquations<DefaultLinAlgAlgorithms,InputIterator,Metric>(eqns_begin,eqns_end,metric);
   	}

/** @brief Construct a Subspace object from a list of linear equations
 *  @param [eqns_begin,eqns_end) A range of equations defining the subspace
 * @returns A Subspace object; the complement is chosen arbitrarily
 * 
 * Much faster than SubspaceFromEquations(InputIterator,InputIterator, const Metric&)
 */
   template<typename LinAlgAlgorithms,typename InputIterator>
   		Wedge::Subspace<T> SubspaceFromEquations(InputIterator eqns_begin,InputIterator eqns_end) const;

/** @overload
 */
   template<typename InputIterator>
   		Wedge::Subspace<T> SubspaceFromEquations(InputIterator eqns_begin,InputIterator eqns_end) const
   	{
   		return SubspaceFromEquations<DefaultLinAlgAlgorithms,InputIterator>(eqns_begin,eqns_end);   		
   	}

/** @brief helper function to translate the result of a call to GiNaC::lsolve to a basis of solutions
*
* @sa GetSolutions
* @todo This should be protected. It is public now because the parambasis segment of the code is outside this class.
*/
	template<typename Container> Container& GetSolutionsFromGenericSolution(Container& container, lst gensol, ex* vector=NULL) const;

/////////////////////////////////////////////////////////////////////
// 		Access to elements of the basis
/////////////////////////////////////////////////////////////////////

/** @brief Return the k-th element of a basis
 *  @param k an integer in the range [1,Dimension()]
 *  @returns An ex representing the k-th element of a basis
*/
	ex e(OneBased k) const {return e()(k);}

/** @brief Returns the basis defining this vector space
 */
	virtual const IBasis<T>& e() const=0;
	
	typename IBasis<T>::const_iterator e_begin() const {return e().begin();} ///< Returns an iterator pointing to the first element of the basis of the subspace   
	typename IBasis<T>::const_iterator e_end() const {return e().end();}///< Returns an iterator pointing after the last element of the basis of the subspace
	
/** @brief Return the n-th coordinate
 *  @param n An integer in the range [1,Dimension()]
 *  @returns An ex representing the n-th coordinate corresponding to the basis of this space
 * 
 * Every VSpace has a fixed basis, which defines a standard set of coordinates.
*/
	ex coordinate(OneBased n) const
	{
		assert(n>0 && n<=coordinates.size());
		return coordinates[n-1];
	}
/** @brief Provide access to the vector of coordinates
*/
	exvector::const_iterator coordinate_begin() const {
		return coordinates.begin();
	}
/** @brief Provide access to the vector of coordinates
*/	
	exvector::const_iterator coordinate_end() const {
		return coordinates.begin()+e().size();
	}	
protected:	
/** @brief Setup a basis of standard coordinates on this space
 * 
 * If generators are added to the space, coordinates which are relative to the old elements of the basis remain valid.
*/
	void UpdateCoordinates()	
	{
		if (e().size()>coordinates.size())	
		coordinates.reserve(e().size());

		for (unsigned i=coordinates.size();i<e().size();i++)
			coordinates.push_back(Coordinate(coordinate_name(i+1)));
	}


private:
	NameAndIndex coordinate_name;
	exvector coordinates;
};

// ****************************************************************************
// *							VectorSpace										  *
// ****************************************************************************
/** @brief A vector space defined by a basis
 *  @param T Type of simple elements
 *  @param LinAlgAlgorithms Type determining the algorithms to use
 */
template<typename T, typename LinAlgAlgorithms=DefaultLinAlgAlgorithms> class VectorSpace : public VSpace<T> {	
public:
   /** @brief Create a trivial vector space (i.e., the set {0}).
     * @param coordinate_name The name associated to coordinates on this vector space (for output only)
    */	
	VectorSpace(string coordinate_name="lambda") : VSpace<T>(coordinate_name) {}

   /** @brief Create a vector space spanned by a family of elements     
     * @param generators A container of linear combinations of T, not necessarily linearly independent
     * @param coordinate_name The name associated to coordinates on this vector space (for output only)
     * 
     * This constructor selects a subset of linearly independent generators to use as a basis.
     * In the case that the generators are already known to be linearly independent, a faster alternative 
     * is to use the default constructor and then call setBasis(). 
    */
    VectorSpace(const exvector& generators,string coordinate_name="lambda"): VSpace<T>(coordinate_name) 
	{
 		basis.insert(basis.end(),generators.begin(),generators.end());
	    	this->UpdateCoordinates();
	}

   /** @brief Create a vector space spanned by a family of elements
     * @param [from,to) A range of linear combinations of T, not necessarily linearly independent
     * @param coordinate_name The name associated to coordinates on this vector space (for output only)
	 *
     * This constructor selects a subset of linearly independent generators to use as a basis.
     * In the case that the generators are already known to be linearly independent, a faster alternative 
     * is to use the default constructor and then call setBasis().  
    */
    template<typename InputIterator>
	VectorSpace(InputIterator from,InputIterator to,string coordinate_name="lambda") : VSpace<T>(coordinate_name) 
 	{
 		basis.insert(basis.end(),from,to);
	    	this->UpdateCoordinates();
	}

/** @brief Copy constructor
 * 
 * Overrides the default constructor, allowing one to initialize a VectorSpace object from a Subspace object.
 */ 
	VectorSpace(const VSpace<T>& V) : VSpace<T>(V) 
	{
 		basis.insert(basis.end(),V.e().begin(),V.e().end());
	    	this->UpdateCoordinates();
	}

	ex e(OneBased k) const {return e()(k);}
	typename IBasis<T>::const_iterator e_begin() const {return basis.begin();} ///< Returns an iterator pointing to the first element of the basis of the subspace   
	typename IBasis<T>::const_iterator e_end() const {return basis.end();}///< Returns an iterator pointing after the last element of the basis of the subspace


	const Basis<T,LinAlgAlgorithms>& e() const {return basis;}

	/** @brief Change the basis defining this vector space
	 * @param basis The new basis
	 */	
	void SetBasis(const Basis<T,LinAlgAlgorithms>& basis)
	{
		this->basis=basis;
	    	this->UpdateCoordinates();
	}
	
	/** @brief Add some generators to this vector space
	 * @param [from,to) A range of linear combinations of elements of type T
	 * @return The resulting variation in dimension (a non-negative integer)
	 */
  
    template<typename InputIterator> int AddGenerators(InputIterator from, InputIterator to)
    {
    	int olddim=this->Dimension();
    	basis.insert(basis.end(),from,to);
    	this->UpdateCoordinates();
    	return this->Dimension()-olddim;
    }

	/** @brief Add a generator to this vector space
	 * @param e A linear combination of elements of type T
	 * @return true If the dimension has increased (i.e. if v was not already in the space)
	 */ 
    bool AddGenerator(ex e)
    {
    	// if basis is initialized, then it is optimal to call Components(); otherwise, better call IBasis::insert()
    	if (basis.NotInitialized() || !this->Contains(e))
    	{
    		int olddim=this->Dimension();
//    		lst my_list(e);    		
//    		basis.insert(basis.end(),my_list.begin(),my_list.end());
    		basis.insert(basis.end(),&e,(&e)+1);    		
    		this->UpdateCoordinates();
    		return this->Dimension()>olddim;
    	}
    	return false;
    }
/** @brief Reset this object, so it represents the trivial vector space
 */

	void Reset() {basis.clear(); this->UpdateCoordinates();}
protected:
	Basis<T,LinAlgAlgorithms> basis; ///< The basis defining this vector space
};


// ****************************************************************************
// *							Subspace									  *
// ****************************************************************************
/** @brief A subspace of a vector space
 *  @param T The type of simple elements in the vector space
 * 
 *  A subspace represents a splitting
 *  \f$ containing space = subspace \oplus complement\f$.
 *  In other words, the complement is part of the definition of the subspace.
 */
template<typename T> class Subspace : public VSpace<T>
{	
public:
/** @brief Create an empty subspace object
 * 
 *  Use operator= or setBasis() to initialize the object
 */
	Subspace() {};

/** @brief Construct a Subspace object from a SubBasis object
 */	
	Subspace(const SubBasis<T>& basis) 
	{
		SetBasis(basis);	
	}

/** @brief Change the basis defining this vector space
 * @param basis The new basis
*/
	void SetBasis(const SubBasis<T>& basis)
	{
		this->basis=basis;
		this->UpdateCoordinates();
	}

/** @brief Reset the subspace to the trivial subspace of the trivial vector space
 */	
	void Reset() {basis.clear(); this->UpdateCoordinates();}
		
/** @brief Create a subspace of a given vector space
 * @param basisOfSubspace A basis of the desired subspace
 * @param complementaryGenerators Elements in the containing space
 * 
 * The containing space is computed as the sum \f$Span(basisOfSubspace)+Span(complementaryGenerators)\f$
 * Note that this need not be a direct sum
*/
	template<typename Container1, typename Container2>
	Subspace(const Container1& basisOfSubspace, const Container2& complementaryGenerators)
	{		
		basis.insert(basis.end(),basisOfSubspace.begin(),basisOfSubspace.end());
		basis.AddToComplement(complementaryGenerators.begin(),complementaryGenerators.end());
		this->UpdateCoordinates();
	}

/** @brief Create a subspace of a given vector space, using a given metric to determine the complement
 * @param basisOfSubspace A basis of the desired subspace
 * @param fullSpaceGenerators Generators of the containing space
 * @param metric The metric
 * 
 * The containing space is computed as the orthogonal sum 
 * \f$Span(basisOfSubspace)\oplus Span(basisOfSubspace)^\perp =  Span(fullSpaceGenerators)\f$
*/
	template<typename Container1, typename Container2, typename Metric>
	Subspace(const Container1& basisOfSubspace, const Container2& fullSpaceGenerators, const Metric& metric)
	{
		VectorSpace<T> V(fullSpaceGenerators.begin(),fullSpaceGenerators.end());
		*this=V.Subspace(basisOfSubspace.begin(),basisOfSubspace.end(),metric);
	}

/** @brief Create a subspace of a given vector space
 * @param [basisOfSubspace_begin,basisOfSubspace_end) a basis of the desired subspace
 * @param [complementaryGenerators_begin, complementaryGenerators_end) elements in the containing space
 * 
 * The containing space is computed as the sum  \f$Span(basisOfSubspace)+Span(complementaryGenerators)\f$
 * Note that this need not be a direct sum.
 * 
 * The complement is constructed as a subspace of \f$Span(complementaryGenerators)\f$. If the above sum
 * is direct, the complement is determined uniquely. 
*/
	template<typename Iterator1, typename Iterator2> Subspace(
		Iterator1 basisOfSubspace_begin,Iterator1 basisOfSubspace_end,
		Iterator2 complementaryGenerators_begin,Iterator2 complementaryGenerators_end)
	{		
		basis.insert(basis.end(),basisOfSubspace_begin,basisOfSubspace_end);
		basis.AddToComplement(complementaryGenerators_begin,complementaryGenerators_end);
		this->UpdateCoordinates();
	}

/** @brief Returns the complement of this subspace
 */
	Subspace Complement() const
	{			
		Subspace S;
		S.SetBasis(SubBasis<T>(complement_begin(),complement_end(),this->e_begin(),this->e_end()));
		return S;
	}

/** @brief Change the containing space, or the basis of the complement.
 * @param [complementaryGenerators_begin,complementaryGenerators_end) The new basis of the complement
 *
 * This function assumes that complementaryGenerators are independent vectors, and 
 * redefines the complement as \f$Span(complementaryGenerators)\f$. It is also assumed that the sum
 * \f$ subspace \oplus Span(complementaryGenerators\f$ is direct.
 */
	template<typename Iterator> void SetComplement(Iterator complementaryGenerators_begin, Iterator complementaryGenerators_end)
	{
		SubBasis<T> b(this->basis.begin(),this->basis.end(),complementaryGenerators_begin,complementaryGenerators_end);
		basis=b;
		this->UpdateCoordinates();
	}
	
/** @brief Project a vector on this subspace
 *  @param v A vector in the containing space
 *  @return The projection of v on the subspace
 */
	ex Project(ex v) const {
		ExVector components=this->e().AllComponents(v);
		ex result=0;
		for (int i=0;i<this->Dimension();i++)
			result+=this->e()[i]*components[i];
		return result;		
	}

/** @brief Project a vector on the complement
 *  @param v A vector in the containing space
 *  @return The projection of v on the complement  
 */
	ex ProjectOnComplement(ex v) const {
		ExVector components=this->e().AllComponents(v);
		ex result=0;
		for (int i=this->Dimension();i<basis.DimensionOfContainingSpace();i++)
			result+=this->e()[i]*components[i];
		return result;
	}
	
	ex e(OneBased k) const {return e()(k);}

	const SubBasis<T>& e() const {return basis;}
	typename IBasis<T>::const_iterator e_begin() const {return basis.begin();} ///< Returns an iterator pointing to the first element of the basis of the subspace   
	typename IBasis<T>::const_iterator e_end() const {return basis.end();}///< Returns an iterator pointing after the last element of the basis of the subspace	
	typename IBasis<T>::const_iterator complement_begin() const {return basis.complement_begin();}///< Returns an iterator pointing to the first element of the basis of the complement
	typename IBasis<T>::const_iterator complement_end() const {return basis.complement_end();}///< Returns an iterator pointing after the last element of the basis of the complement
	typename IBasis<T>::const_iterator full_begin() const {return basis.full_begin();} ///< Returns an iterator pointing to the first element of the basis of the containing space
	typename IBasis<T>::const_iterator full_end() const {return basis.full_end();} ///< Returns an iterator pointing after the last element of the basis of the containing space
private:
	SubBasis<T> basis;
};


/** @brief Implements a linear map as a list of substitutions
 *  @param [v_begin,v_end) A range of linearly independent elements in V
 *  @param [fv_begin,fv_end) The corresponding range of images under the linear map
 *  @param options Parameter to pass to GiNaC::lsolve (has no effect if the basis b is smaller than the set of simple elements it contains)
 *  @return A list of equations, to be used as an argument to ex::subs
 */
template<typename V, typename Iterator1, typename Iterator2> 
lst LinearMapToSubstitutions(Iterator1 v_begin,Iterator1 v_end, Iterator2 fv_begin, Iterator2 fv_end,unsigned options = solve_algo::automatic)
{
	Basis<V> source(v_begin,v_end);
	return LinearMapToSubstitutions(source, fv_begin, fv_end);
}

/** @brief Implements a linear map as a list of substitutions
 *  @param source A basis of linearly independent elements in V
 *  @param [fv_begin,fv_end) The corresponding range of images under the linear map
 *  @param options Parameter to pass to GiNaC::lsolve (has no effect if the basis b is smaller than the set of simple elements it contains)
 *  @return A list of equations, to be used as an argument to ex::subs
 */
template<typename V, typename LinAlgAlgorithms, typename Iterator> 
lst LinearMapToSubstitutions(const Basis<V,LinAlgAlgorithms>& source, Iterator fv_begin, Iterator fv_end,unsigned options = solve_algo::automatic)
{
	LOG_DEBUG(source);

	exvector symbols; GetSimple<V>(symbols,source.begin(),source.end());
	lst result;

	if (symbols.size()==source.size())	//we can use GiNaC::lsolve
	{
		VectorSpace<V> s; 
		s.SetBasis(source);
		LOG_DEBUG(s);
		for (exvector::const_iterator i=symbols.begin();i!=symbols.end();i++)
		{		
			exvector comps=s.Components(*i,options);
			LOG_DEBUG(comps);
			assert(comps.size()==symbols.size());		
			ex fi;
			Iterator j=fv_begin;
			exvector::const_iterator k=comps.begin();
			while (j!=fv_end)
			{
				assert(k!=comps.end());
				fi+=*k * *j;
				++j,++k;
			}
			assert(k==comps.end());
			result.append(*i == fi.expand());
		}
	}
	else 				//we have to use Basis
		for (exvector::const_iterator i=symbols.begin();i!=symbols.end();i++)
		{		
			exvector comps=source.AllComponents(*i);
			assert(comps.size()==symbols.size());
			LOG_DEBUG(comps);
			ex fi;
			Iterator j=fv_begin;
			exvector::const_iterator k=comps.begin();
			while (j!=fv_end)
			{
				assert(k!=comps.end());
				fi+=*k * *j;
				++j,++k;
			}			
			result.append(*i == fi.expand());
		}
	return result;
}



/** @brief Implements a linear map as a list of substitutions
 *  @param v Linearly independent elements in V
 *  @param fv The images of the elements in v
 *  @return A list of equations, to be used as an argument to ex::subs
 */
template<typename V, typename LinAlgAlgorithms> lst LinearMapToSubstitutions(const Basis<V,LinAlgAlgorithms>& v, const exvector& fv,unsigned options = solve_algo::automatic)
{
	return LinearMapToSubstitutions(v,fv.begin(),fv.end(),options);
}

/** @brief Implements a linear map as a list of substitutions
 *  @param v Linearly independent elements in V
 *  @param fv The images of the elements in v
 *  @return A list of equations, to be used as an argument to ex::subs
 */
template<typename V> lst LinearMapToSubstitutions(const exvector& v, const exvector& fv,unsigned options = solve_algo::automatic)
{
	return LinearMapToSubstitutions<V>(v.begin(),v.end(),fv.begin(),fv.end(),options);
/*	assert(v.size()==fv.size());
	Basis<V> source(v.begin(),v.end());
	VectorSpace<V> s(source);

	exvector symbols; GetSimple<V>(symbols,v.begin(),v.end());
	lst result;

	for (exvector::const_iterator i=symbols.begin();i!=symbols.end();i++)
	{		
		exvector comps=s.Components(*i,options);
		assert(comps.size()==symbols.size());		
		ex fi=0;
		for (int j=0;j<fv.size();j++)
		{
			fi+=comps[j] * fv[j];			
		}
		result.append(*i == fi);
	}		
	return result;
*/
}

///////////////////////////////////////////////////
//	 	Implementation
///////////////////////////////////////////////////

template<typename T> ex VSpace<T>::GenericElement() const
{
	ex gen_elem;
	exvector::const_iterator i=e().begin();
	exvector::const_iterator j=coordinates.begin();
	while (i!=e().end()) {
		assert(j!=coordinates.end());
		gen_elem+=*j++ * *i++;
	}
	return gen_elem;
}


template<typename T> 
template<typename Container> Container& VSpace<T>::GetSolutionsFromGenericSolution(Container& container, lst gensol, ex* vector) const
{
	if (gensol==lst() && Dimension()>0) 
		throw EmptyAffineSpace(__FILE__,__LINE__);
	
	list<ex> free; 
	lst substitutions;
	for (lst::const_iterator i=gensol.begin();i!=gensol.end();++i)
		if (i->lhs()==i->rhs()) {
			free.push_back(i->lhs());
			substitutions.append(i->lhs()==0);
		}

	ex gen_elem=GenericElement().subs(gensol);
//compute a solution, chosen by setting all free variables to zero
	ex onesolution=gen_elem.subs(substitutions);
	if (vector!=NULL) *vector=onesolution;
	if (!onesolution.is_zero())  	//if the space solution is an affine space...
	{
		if (vector==NULL) {	//...make sure the caller was expecting it...
			LOG_WARN(onesolution);
			LOG_WARN(*this);
			throw WedgeException<std::runtime_error>("The affine space of solutions is not a vector space",__FILE__,__LINE__);
		}
		gen_elem-=onesolution;	//...and translate to a vector space
	}
	for (list<ex>::const_iterator i=free.begin();i!=free.end();i++)
	{
		ex element=(gen_elem.subs(*i==1)).subs(substitutions);
		Insert(container,element);
	}
	return container;
}


template<typename T> 
template<typename LinAlgAlgorithms,typename Container,typename Iterator> Container& VSpace<T>::GetSolutions(Container& container, Iterator first, Iterator last, ex* vector) const
{	
	lst eqns;
	while (first!=last)
		eqns.append((*first++)==0);
	assert(coordinates.size()==Dimension());
	lst result=LinAlgAlgorithms::lsolve(eqns,lst(coordinates.begin(),coordinates.end()));
	try {
		GetSolutionsFromGenericSolution(container,result,vector);
	}
	catch (...) {
		LOG_WARN(eqns);
		LOG_WARN(coordinates);
		throw;
	}
	return container;
}


template<typename T>
template<typename LinAlgAlgorithms, typename Iterator, typename Metric>
Subspace<T> VSpace<T>::Subspace(Iterator subbasis_begin, Iterator subbasis_end, const Metric& metric) const
{
	static_cast<void>	//compile-time check
	(
		static_cast<const IBilinearOperator<LinearOperator<T>,LinearOperator<T> >& >(metric)
	);		
	
	exvector orthequations;
	GetOrthogonalEquations(orthequations,subbasis_begin,subbasis_end,metric);
	exvector orthogonal;
	GetSolutions<LinAlgAlgorithms>(orthogonal, orthequations.begin(),orthequations.end());
	return Wedge::Subspace<T>(subbasis_begin,subbasis_end,orthogonal.begin(),orthogonal.end()); 
}

template<typename T>
template<typename Iterator>
Subspace<T> VSpace<T>::Subspace(Iterator subbasis_begin, Iterator subbasis_end) const
{
	return Wedge::Subspace<T>(subbasis_begin,subbasis_end,e_begin(),e_end());
}

template<typename T>
template<typename Container, typename Iterator, typename Metric> Container& VSpace<T>::GetOrthogonalEquations(Container& container, Iterator subbasis_begin, Iterator subbasis_end, const Metric& metric) const
{
	static_cast<void>	//compile-time check
	(
		static_cast<const IBilinearOperator<LinearOperator<T>,LinearOperator<T> >& >(metric)
	);		
	
	ex gen_elem=GenericElement();
	while (subbasis_begin!=subbasis_end)
	{
		ex product=IBilinearOperator<LinearOperator<T>,LinearOperator<T> > ::BilinearOperator(gen_elem,*subbasis_begin++,&metric);
		Insert(container,product);
	}
	return container;
}

template<typename T>
template<typename LinAlgAlgorithms,typename InputIterator, typename Metric>
		Subspace<T> VSpace<T>::SubspaceFromEquations(InputIterator eqns_begin,InputIterator eqns_end, const Metric& metric) const
{
	exvector basis;
	GetSolutions<LinAlgAlgorithms>(basis, eqns_begin,eqns_end);
	return Subspace<LinAlgAlgorithms>(basis.begin(),basis.end(),metric);
}


template<typename T>
template<typename LinAlgAlgorithms,typename InputIterator>
		Subspace<T> VSpace<T>::SubspaceFromEquations(InputIterator eqns_begin,InputIterator eqns_end) const
{
	exvector basis;
	GetSolutions<LinAlgAlgorithms>(basis, eqns_begin,eqns_end);
	return Subspace(basis.begin(),basis.end());
}

/** @brief Overloaded output operator */
template<class charT, class traits, class T> std::basic_ostream<charT,traits>& operator<<(std::basic_ostream<charT,traits>& out, const VSpace<T>& V)
{
	if (V.Dimension()==0) return out<<"Trivial vector space"<<endl;
	else return out<<V.e();
}

/** @brief Overloaded output operator */
template<class charT, class traits, class T> std::basic_ostream<charT,traits>& operator<<(std::basic_ostream<charT,traits>& out, const Subspace<T>& V)
{
	out<<"Subspace:"<<Range(V.e_begin(),V.e_end());
	out<<"Complement:"<<Range(V.complement_begin(),V.complement_end());
	return out;
}


} /** @} */

#endif
