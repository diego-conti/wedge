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
#ifndef WEDGE_GSTRUCTURE_H
#define WEDGE_GSTRUCTURE_H

/** @ingroup RiemannianGeometry */ 

/** @{ 
 * @file gstructure.h
 * @brief Abstract G-structures
 */
 
#include "wedge/manifolds/concretemanifold.h"
#include "wedge/base/parameters.h"

namespace Wedge {

WEDGE_DECLARE_NAMED_ALGEBRAIC(IntrinsicTorsionClass,realsymbol)

/** @brief A class to represent the intrinsic torsion of a G-structure
 *
 * IntrinsicTorsion is an associative array, mapping each intrinsic torsion class to the relative component
 * 
 * @remark The intrinsic torsion of a  G-structure (with \f$G\subset O(n)\f$) is a section of a bundle with fibre \f$ T^*\otimes \mathfrak{g}^\perp\f$.
 * The intrinsic torsion classes are the irreducible components of this G-module. 
 */
struct IntrinsicTorsion : public map<IntrinsicTorsionClass,ex,basic_is_less>
{
	/** @brief Return the intrinsic torsion class
	 *  @return A formal sum of intrinsic torsion classes
	 * 
	 *  A null result means that the structure is integrable. A linear combination, say, \f$W_1+W_2\f$ means that
	 *  the non-zero components of the intrinsic torsion are exactly \f$W_1\f$ and \f$W_2\f$.
	 */
	ex Type() const {
		ex type;
		for (const_iterator i=begin();i!=end();i++)
			if (!i->second.expand().is_zero()) type+=i->first;
		return type;
	}

	/** @brief Return the generic intrinsic torsion class
	 *  @return A formal sum of intrinsic torsion classes
	 *  @exception WedgeException<std::logic_error> Thrown if the object is not the result of a call to GStructure::GetIntrinsicTorsion() 
	 */
	ex GenericType() const
	{
		if (empty()) throw WedgeException<std::logic_error>("Cannot obtain generic type from uninitialized IntrinsicTorsion object",__FILE__,__LINE__);
		ex type;
		for (const_iterator i=begin();i!=end();i++)
			type+=i->first;
		return type;		
	}

/** @brief Check whether all coefficients of a formal torsion class are non-negative
 * @param type A linear combination of objects of type IntrinsicTorsionClass
 *
 * @todo It would be nice to be able to  compare intrinsic torsion classes like, say, divisors.
 */	
	static bool IsNonNegative(ex type)
	{
		exvector v;
		GetSimple<IntrinsicTorsionClass>(v,type);
		for (int i=0;i<v.size();i++)
			if (type.coeff(v[i])<0) return false;
		return true;
	}
	
	
/** @brief Check whether two torsion classes have no component in common
 * @param type1,type2 Linear combinations of objects of type IntrinsicTorsionClass, with all coefficients equal to 0 or 1
 * @return true if all coefficients of type1+type2 are either 0 or 1  
 */	
	static bool AreTypesComplementary(ex type1, ex type2)
	{
		ex type=type1+type2;
		exvector v; GetSimple<IntrinsicTorsionClass>(v,type);
		for (int i=0;i<v.size();i++)
			if (type.coeff(v[i])>1) return false;
		return true;
	}		
};

WEDGE_DECLARE_NAMED_ALGEBRAIC(GStructureParameter,realsymbol)

class LieGroup;

/** @brief A class representing a G-structure, possibly with parameters
 * 
 * Mathematically, a G-structure is a reduction to \f$G\subset GL(n,\mathbb{R})\f$ of the bundle of frames. This class
 * represents a G-structure by the choice of a single frame, i.e. a fixed section of the G-reduction, which will be called 
 * \em the \em adapted \em frame.  
 * 
 * @note The parameters of a generic GStructure should be of type GStructureParameter; they are assumed to be constant
 * on the manifold.
 *   
 * @remark Specific G-structures, relative to a specific choice of G, are implemented as subclasses of GStructure 
 * respecting the natural hierarchy, in the sense that if H is a subgroup of K,  then HStructure derives from KStructure. In particular, Riemannian metrics are defined
 * as GStructure's, and so GStructure's corresponding to compact G should derive from RiemannianStructure.
 */

class GStructure {
public:
/** @brief An element of the adapted frame for this GStructure
 *  @param i An index
 *  @returns The i-th element of the adapted frame
 */ 
	ex e(OneBased i) const {
		return frame[i-1];
	};
/** @brief The adapted frame for this GStructure 
 */ 
	const Frame& e() const {
		return frame;	
	}

	virtual ~GStructure() {}
protected:
/** @brief Define a GStructure on a manifold, represented by the choice of a frame
 *  @param manifold The manifold on which the GStructure is defined.
 *  @param frame The adapted frame defining the structure. By definition, acting on the frame
 *  by an element of G would define the same structure, though with a different %internal representation
 * 
 *  @warning Caller must ensure that the pointer remains valid.
 * 
 *  @note This constructor is protected because it constructs an object with no information on the structure group G; hence,
 *  it is assumed  that this information is contained in a subclass.
 */
	GStructure(const Manifold* manifold, const Frame& frame) : manifold_{manifold} {
		SetFrame(frame);
	};


/** @brief Alter the GStructure by redefining the adapted frame
 *  @param frame The adapted frame. By definition, acting on the frame
 *  by an element of G would define the same structure, though with a different %internal representation.
 */
	void SetFrame(const Frame& frame) {
		assert(frame.size()==manifold_->Dimension());
		this->frame=frame;
	}

	const Manifold* M() const {return manifold_;} ///< The underlying manifold
	
/** @brief Construct a GStructure object without initializing the frame
 *  @param manifold The manifold on which the GStructure is defined.
 * 
 *  @warning Caller must ensure that the pointers remain valid.
 *  
 *  Used by those derived classes that initialize the frame in the constructor
 */	
	GStructure(const Manifold* manifold) : manifold_{manifold} {}

private:
	const Manifold* manifold_ {nullptr};
	Frame frame;
};

/** @brief Output operator
 * 
 * Should be reimplemented for subclasses
 */ 
template<class charT, class traits> std::basic_ostream<charT,traits>& operator<<(std::basic_ostream<charT,traits>& os, const GStructure& gStructure)
{
	return os<<"GStructure with adapted frame:"<<endl<<gStructure.e();
}


/** @brief Class template with explicit specializations.
 */
template<class Structure, bool WithParams> class GStructureHasParameters {};


/////////////////////////////////////////////////////////////////////////////////
//			Template specializations for GStructures with/without parameters
/////////////////////////////////////////////////////////////////////////////////

/** @brief Abstract base class for GStructure's whose intrinsic torsion does not depend on any parameter
 *  
 * @param Structure A class derived from GStructure, e.g. SU3Structure, G2Structure... 
*/
template<class Structure> class GStructureHasParameters<Structure, false> : public Structure {
public:
/** @brief Define a GStructure on a manifold, represented by the choice of a frame.
 *  @param manifold The manifold on which the GStructure is defined.
 *  @param frame The adapted frame defining the structure.
 *  
 *  @warning Caller must ensure that the pointer remains valid.
 * 
 * @remark Since the intrinsic torsion should not depend on parameters, the \f$d\f$ operator on the manifold argument should not depend on parameters.
 * @exception WedgeException<std::logic_error> Thrown if manifold does not know how to take d of its forms
 */
	GStructureHasParameters(const Manifold* manifold, const Frame& frame) : Structure(manifold,frame) {
		if (!manifold->KnowsHowToCompute_d())
		{
			LOG_MSG("Trying to construct a GStructureHasParameters<.,false> object from a Manifold object that does not know how to take d of its forms. Was TorsionFreeConnection<false> intended?");
			throw WedgeException<std::logic_error>("KnowsHowToCompute_d returned false in GStructureHasParameters<.,false>::GStructureHasParameters()",__FILE__,__LINE__);
		}
	}
	
/** @brief Compute the intrinsic torsion of this G-structure
 * 
 *  @note The result is cached in GStructureHasParameters::it
 *
 * @todo This function is declared but not defined, and a definition should be given for each choice of the template parameter, using template specialization. Some specializations are missing.
 */
	IntrinsicTorsion GetIntrinsicTorsion() const;
protected:
	void SetFrame(const Frame& frame) {Structure::SetFrame(frame); it.clear();}
	GStructureHasParameters(const Manifold* manifold) : Structure(manifold) {}
	mutable IntrinsicTorsion it;
}; 

/** @brief Abstract base class for GStructure's whose intrinsic torsion depends on parameters
 * @param Structure A class derived from GStructure, e.g. SU3Structure, G2Structure...  
 * 
 * There can be parameters in the adapted frame (which must have type GStructureParameter), or parameters that determine
 * the action of d (e.g. the manifold may be a LieGroupWithParameters, or even a ManifoldWith).
 */
template<class Structure> class GStructureHasParameters<Structure, true> : public Structure, public HasParameters<GStructureParameter>
{
public:
/** @brief Define a GStructure on a manifold, represented by the choice of a frame
 *  @param manifold The manifold on which the GStructure is defined.
 *  @param frame The adapted frame defining the structure, whose elements depend on parameters.
 * 
 * The parameters may have type GStructureParameter, in which case one can use the functions provided by HasParameters. 
 *  
 *  @warning Caller must ensure that the pointer remains valid.
 */
	GStructureHasParameters(const Manifold* manifold, const Frame& frame) : Structure(manifold,frame) {}
	
/** @brief Compute the equations in the parameters corresponding to a certain torsion class
 * @param container (out) A reference to a container object where the equations are to be stored.  
 * @param TorsionClass A formal sum of IntrinsicTorsionClass elements
 * @return A reference to the container
 * 
 * If, say, \f$ T^*\otimes \mathfrak{g}^\perp\f$ splits as \f$W1\oplus W2\oplus W_3\f$, we say a G-structure has torsion class W1+W2 
 * if the torsion is contained in \f$W1\oplus W2\subset T^*\otimes \mathfrak{g}^\perp\f$, i.e. the W3 component is zero. Invoking this function
 * with TorsionClass=W1+W2 returns the equations corresponding to this condition.
 * 
 * @note The type Container can be PolyBasis<GStructureParameter>, PolyBasis<LieGroupParameter> etc.
 *
 * @todo This function is declared but not defined, and a definition should be given for each choice of the template parameter, using template specialization. Some specializations are missing.   
 */	
	template<typename Container> Container& GetEquationsForTorsionIn(Container& container, ex TorsionClass) const;
protected:
/** @copydoc GStructure::GStructure(const Manifold*)
 */
	GStructureHasParameters(const Manifold* manifold) : Structure(manifold) {}

	void DeclareConditions(const lst& list_of_equations)
	{
		exvector frame; frame.reserve(this->e().size());
		for (Frame::const_iterator i=this->e().begin();i!=this->e().end();i++)
			frame.push_back(i->subs(list_of_equations));
		this->SetFrame(frame);		
	}
};

} /** @} */
#endif
