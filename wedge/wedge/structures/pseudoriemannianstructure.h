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
#ifndef PSEUDORIEMANNIANSTRUCTURE_H
#define PSEUDORIEMANNIANSTRUCTURE_H
/** @ingroup RiemannianGeometry

/** @{ 
 * @file pseudoriemannianstructure.h
 * @brief Pseudoriemannian metrics as G-structures
 */
 
#include "../structures/gstructure.h"
#include "wedge/manifolds/manifold.h"
#include "wedge/linearalgebra/vectorspace.h"
#include "wedge/linearalgebra/bilinearform.h"
#include "spinor.h"

namespace Wedge {


/** @brief Pseudoriemannian metric on a manifold, 
 * 
 * This is an abstract base class
 */
 
class PseudoRiemannianStructure : public GStructure
{
	ex Hook(ex) const = delete;	//Not implemented. Use ScalarProduct().Interior() instead.
public:
	virtual const BilinearForm& ScalarProduct() const=0;

	ex HodgeStar(ex form) const {return ScalarProduct().Interior(form,ncmul(e()));}

/** @brief Define a Pseudoriemannian structure on a manifold, represented by the choice of a frame
 *  @param manifold The manifold on which the structure is defined.
 *  @param orthonormalframe A coframe with respect to which the metric is defined; could be orthonormal or not.
 * 
 *  @warning Caller must ensure that the pointer remains valid.
 */ 
	PseudoRiemannianStructure(const Manifold* manifold, const Frame& frame) :
		GStructure(manifold,frame) {}

/** @brief Overloaded copy constructor
 *  
 * @warning Every PseudoriemannianStructure object contains a pointer to a Manifold object, representing the underlying manifold. The pointer is copied by
 * this constructor, so the caller must ensure that it remains valid.
 *
 * @todo Replace all pointers with shared_ptr
 */ 
	PseudoRiemannianStructure(const PseudoRiemannianStructure& o) : GStructure(o) {
		/* write a message to the log, since use of copy constructor may lead to subtle bugs, i.e.
		 * PseudoriemannianStructure f() {
		 * 	ConcreteManifold M(5);
		 *  return PseudoriemannianStructure(&M,M.e()); //crash
		 * }
		 */
		LOG_WARN("PseudoRiemannianStructure copy constructor invoked; beware the pointer.");
	}

	virtual pair<ex,matrix> DecomposeRicci(matrix ricci) const=0;
};

/** @brief Pseudoriemannian metric on a manifold, represented by an orthonormal coframe
 * 
 * @deprecated Use PseudoRiemannianStructureByFrame::FromSignature
 */
 
class StandardPseudoRiemannianStructure : public PseudoRiemannianStructure {
	const int p_;
	const StandardScalarProduct scalar_product;
public:
/** @brief 
 *  @param manifold The manifold on which the structure is defined.
 *  @param orthonormal_frame A coframe with respect to which the metric is defined; could be orthonormal or not.
 *  @param p An integer that determines the signature
 *
 * The orthonormal frame e_1,..., e_n is assumed to satisfy <e_i,e_j>=\delta_{ij} for i\leq p, -\delta_{ij} for i>p
*/
	StandardPseudoRiemannianStructure(const Manifold* manifold, const Frame& orthonormal_frame, int p) : 
		PseudoRiemannianStructure(manifold,orthonormal_frame), p_{p}, scalar_product{orthonormal_frame,p} {}
	const BilinearForm& ScalarProduct() const {return scalar_product;}

	pair<ex,matrix> DecomposeRicci(matrix ricci) const override;
};


/** @brief Pseudoriemannian metric on a manifold, represented by a matrix relative to some frame
 */

class PseudoRiemannianStructureByMatrix : public PseudoRiemannianStructure {
	const ScalarProductDefinedByMatrix scalar_product;
	PseudoRiemannianStructureByMatrix(const Manifold* manifold, const Frame& frame, ScalarProductDefinedByMatrix&& scalar_product) :
		PseudoRiemannianStructure(manifold,frame), scalar_product{std::move(scalar_product)} {}
public:
/** @brief
 *  @param manifold The manifold on which the structure is defined.
 *  @param frame A coframe with respect to which the metric is defined; could be orthonormal or not.
 *  @param m The metric relative to the coframe
*/
	static PseudoRiemannianStructureByMatrix FromMatrixOnFrame(const Manifold* manifold, const Frame& frame, const matrix& m) {
		return PseudoRiemannianStructureByMatrix{manifold,frame,ScalarProductDefinedByMatrix::OnFrame(frame,m)};
	}
	static PseudoRiemannianStructureByMatrix FromMatrixOnCoframe(const Manifold* manifold, const Frame& frame, const matrix& m) {
		return PseudoRiemannianStructureByMatrix{manifold,frame,ScalarProductDefinedByMatrix::OnCoframe(frame,m)};
	}

	const BilinearForm& ScalarProduct() const {return scalar_product;}

	pair<ex,matrix> DecomposeRicci(matrix ricci) const override {throw NotImplemented(__FILE__,__LINE__,"PseudoRiemannianStructureByMatrix::DecomposeRicci");}
};

/** @brief Enum class for choice of convention in the Clifford product.
 * 
 * The values are:
 * 
 * - BAUM_KATH for the convention of 
 * 		Helga Baum and Ines Kath. “Parallel Spinors and Holonomy Groups on Pseudo-Riemannian Spin Manifolds”. en. In: Annals of Global Analysis and
 * 		Geometry 17.1 (Feb. 1999), pp. 1–17. issn: 1572-9060. doi: 10.1023/A:1006556630988,	
 * 		for which to make the volume form acts as i^{(r-s+1)/2} in odd dimension 
 * - STANDARD to make the volume form act as i^{q+n(n+1)/2} on spinors of positive chirality (even dimension) or all spinors (odd dimension)
 */
enum class CliffordConvention {
	BAUM_KATH, 	/* convention of:
		Helga Baum and Ines Kath. “Parallel Spinors and Holonomy Groups on Pseudo-Riemannian Spin Manifolds”. en. In: Annals of Global Analysis and
		Geometry 17.1 (Feb. 1999), pp. 1–17. issn: 1572-9060. doi: 10.1023/A:1006556630988,	
	*/
	STANDARD	/* choose signs so that the volume form  acts as i^{q+n(n+1)/2} on spinors of positive chirality (even dimension) or all spinors (odd dimension)*/	
};

class CliffordProduct;
class CliffordProductForm;


/** @brief Base class for pseudoriemannian metrics on a manifold, represented by an orthogonal or orthonormal coframe
 */
 
 class PseudoRiemannianStructureByFrame : public PseudoRiemannianStructure {	
protected:
	PseudoRiemannianStructureByFrame(const Manifold* manifold, const Frame& frame, CliffordProduct* clifford_product);	
	pair<ex,matrix> DecomposeRicci(matrix ricci) const override;
public:
/**
   @brief Returns the k-th element of a global basis of complex spinors
   @param k An index in the range [0,DimensionOfSpinorRepresentation())
   @return The spinor \f$ u(\epsilon_m,\dots,\epsilon_1)\f$ where \epsilon_i=1 if the i-th least significant digit in base 2 of k is 0 and -1 otherwise
*/
	ex u(ZeroBased k) const;	

/**
   @brief Returns the complex spinor u(\epsilon_m,...,\epsilon_1)
   @param signs The sequence \epsilon_1,...,\epsilon_m
   @return A section of the complex spinor bundle
*/
	ex u(const vector<int>& signs) const;

/** @brief Compute the Clifford action of a vector field on a spinor
 * @param X A vector field
 * @param psi A spinor
 * @return The spinor \f$ X\cdot\psi\f$.
	 
  * @remark The implementation of Clifford multiplication is based on the formulae of 
 * [Baum, H. and Kath, I. Parallel Spinors and Holonomy Groups on Pseudo-Riemannian Spin Manifolds. Annals of Global Analysis and Geometry 17: 1–17, 1999.]
  *  
  * If the dimension n is odd, the Clifford multiplication by e_n is chosen with the sign that makes the volume form act as i^{(r-s+1)/2} or i^(s+n(n+1)/2) according to the value of the CliffordConvention parameter passed on initialization @sa CliffordConvention
*/
	ex CliffordDot(ex X, ex psi) const;

/** @brief Compute the Clifford action of a vector field on a spinor
 * @param alpha a differential form 
 * @param psi A spinor
 * @return The spinor \f$ alpha\cdot\psi\f$.
*/
	ex CliffordDotByForm(ex alpha, ex psi) const;
/**
   @brief Returns the rank of the complex spinor bundle
   @return The number \f$2^{[n/2]}\f$, where \f$n\f$ is the manifold's dimension
 */
	int DimensionOfSpinorRepresentation() const;
protected:
	struct Deleter {
  		void operator()(CliffordProduct* r);    
  		void operator()(CliffordProductForm* r);    
	};
 	unique_ptr<CliffordProduct,Deleter> clifford_product_operator;
	unique_ptr<CliffordProductForm,Deleter> clifford_product_form_operator;
};

/** @brief Pseudoriemannian metric on a manifold, represented by an orthonormal coframe
 */
 
class PseudoRiemannianStructureByOrthonormalFrame : public PseudoRiemannianStructureByFrame {	
	const ScalarProductByOrthonormalFrame scalar_product;
protected:
	PseudoRiemannianStructureByOrthonormalFrame(const Manifold* manifold, const Frame& frame, ScalarProductByOrthonormalFrame&& scalar_product, CliffordConvention clifford_convention);
public:
/** @brief 
 *  @param manifold The manifold on which the structure is defined.
 *  @param orthonormal_frame An orthonormal coframe with respect to which the metric is defined
 *  @param p A sequence of one-based indices corresponding to timelike elements in the orthonormal coframe; if empty, the metric is positive definite
 *
 * The orthonormal frame e_1,..., e_n is assumed to satisfy <e_i,e_j>=-\delta_ij or \delta_ij according to whether i is in timelike_indices
*/
	static PseudoRiemannianStructureByOrthonormalFrame FromTimelikeIndices(const Manifold* manifold, const Frame& orthonormal_frame, const list<int>& timelike_indices, CliffordConvention clifford_convention) {
		return PseudoRiemannianStructureByOrthonormalFrame{manifold,orthonormal_frame,ScalarProductByOrthonormalFrame::FromTimelikeIndices(orthonormal_frame,timelike_indices),clifford_convention};
	}
	static PseudoRiemannianStructureByOrthonormalFrame FromTimelikeIndices(const Manifold* manifold, const Frame& orthonormal_frame, const vector<int>& timelike_indices, CliffordConvention clifford_convention) {
		return FromTimelikeIndices(manifold,orthonormal_frame,list<int>(timelike_indices.begin(),timelike_indices.end()),clifford_convention);
	}
	static PseudoRiemannianStructureByOrthonormalFrame FromTimelikeIndices(const Manifold* manifold, const Frame& orthonormal_frame, std::initializer_list<int> timelike_indices, CliffordConvention clifford_convention) {
		return FromTimelikeIndices(manifold,orthonormal_frame,list<int>(timelike_indices.begin(),timelike_indices.end()),clifford_convention);
	}

/** @brief 
 *  @param manifold The manifold on which the structure is defined.
 *  @param orthonormal_frame An orthonormal coframe with respect to which the metric is defined
 *  @param signs A sequence \epsilon_1,...,\epsilon_n, where each element is either 1 or -1
 *
 * The metric is taken to be of the form \epsilon_1 e^1\otimes e^1+...+\epsilon_n e^n\otimes e^n
*/
	static PseudoRiemannianStructureByOrthonormalFrame FromSequenceOfSigns(const Manifold* manifold, const Frame& orthonormal_frame, const vector<int>& signs, CliffordConvention clifford_convention) {
		return PseudoRiemannianStructureByOrthonormalFrame{manifold,orthonormal_frame,ScalarProductByOrthonormalFrame::FromSequenceOfSigns(orthonormal_frame,signs),clifford_convention};
	}
/** @brief 
 *  @param manifold The manifold on which the structure is defined.
 *  @param orthonormal_frame An orthonormal coframe with respect to which the metric is defined
 *  @param signature The pair p,q defining the signature
 *
 * The metric is taken to be of the form e^1\otimes e^1+... e^p\otimes e^p - e^{p+1}\otimes e^{p+1}- ... - e^n\otimes e^n
*/
	static PseudoRiemannianStructureByOrthonormalFrame FromSignature(const Manifold* manifold, const Frame& orthonormal_frame, pair<int,int> signature, CliffordConvention clifford_convention) {
		int p=signature.first;
		vector<int> signs(manifold->Dimension());
		fill(signs.begin(),signs.begin()+p,1);
		fill(signs.begin()+p,signs.end(),-1);		
		return PseudoRiemannianStructureByOrthonormalFrame{manifold,orthonormal_frame,ScalarProductByOrthonormalFrame::FromSequenceOfSigns(orthonormal_frame,signs),clifford_convention};
	}
	const ScalarProductByOrthonormalFrame& ScalarProduct() const override {return scalar_product;}
};

/** @brief Pseudoriemannian metric on a manifold, represented by an orthogonal coframe and the diagonal entries of the metric tensor
 */
 
class PseudoRiemannianStructureByOrthogonalFrame : public PseudoRiemannianStructureByFrame {	
	const ScalarProductByOrthogonalFrame scalar_product;
public:
/** @brief 
 *  @param manifold The manifold on which the structure is defined.
 *  @param orthogonal_coframe An orthogonal coframe e^1,...,e^n with respect to which the metric is defined
 *  @param g A sequence g_1,...,g_n such that the metric takes the form g_1e^1\otimes e^1+... + g_ne^n\otimes e^n
*/
	PseudoRiemannianStructureByOrthogonalFrame(const Manifold* manifold, const Frame& orthogonal_coframe, const ExVector& g, CliffordConvention clifford_convention);

	const ScalarProductByOrthogonalFrame& ScalarProduct() const override {return scalar_product;}
};


} /** @} */

#endif
