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
#ifndef BILINEARFORM_H
#define BILINEARFORM_H
/** @ingroup LinearAlgebra*/ 

/** @{ 
 * @file bilinearform.h
 * @brief Scalar products on the tangent space
 */
 
#include "wedge/manifolds/manifold.h"
#include "wedge/linearalgebra/vectorspace.h"

namespace Wedge {
/** @brief A bilinear form on the tangent space
*
* Abstract base class
*/
class BilinearForm {
public:
	BilinearForm() {}
	virtual ex OnOneForms(ex v, ex w) const=0;
	virtual ex OnVectors(ex v, ex w) const =0;
/** @brief Computes the bilinear form on a pair of forms
 *
 * Default implementation puts v and w in normal form, i.e. expanded in linear combinations of decomposable forms; on each pair of decomposable forms, the scalar product is computed invoking OnOneForms
 */
	virtual ex OnForms(ex v, ex w) const;
/** @brief The musical isomorphism, \\alpha \\mapsto \\sum \langle \alpha, e^j\rangle e_j */
	virtual ex Sharp(ex form) const =0;
	virtual ex Flat(ex vector) const =0;
/** @brief Computes the interior product of a form into another
 */
	ex Interior(ex v, ex w) const;
	virtual ~BilinearForm() {}
private:
	ex OnSimpleForms(ex v, ex w) const;	//v and w are either DifferentialForm's or containers of the form [x,y,z] representing the form x\wedge y\wedge z
};

/** @brief A bilinear form on the tangent space, represented by a coframe and a matrix
*
* Abstract base class, which implements bilinearity. Derived class must define the matrix that determines the bilinear form
*/

class BilinearFormWithFrame : public BilinearForm {
	Frame frame_;
public:
/** @brief Initialize the coframe*/
	BilinearFormWithFrame(const Frame& frame) : frame_(frame) {};
/** @brief Computes the bilinear form on a pair of one-forms
 *
 * The one-forms are decomposed according to the reference frame; then the matrix provided by the subclass is used
 */
	ex OnOneForms(ex v, ex w) const {
		ex result;
		ExVector comps_v=frame_.Components(v);
		ExVector comps_w=frame_.Components(w);
		for (int i=1;i<=comps_v.size();++i)
		for (int j=1;j<=comps_w.size();++j)
			result+=MatrixEntry(i,j)*comps_v(i)*comps_w(j);
		return result.expand();
	}
/** @brief Computes the bilinear form on a pair of vectors
 *
 * The one-forms are decomposed according to the reference frame; then the matrix provided by the subclass is used
 */
	ex OnVectors(ex v, ex w) const {
		ex result;
		ExVector comps_v, comps_w; 
		for (exvector::const_iterator i=frame_.begin();i!=frame_.end();++i) {
			comps_v.push_back(Hook(v,*i));
			comps_w.push_back(Hook(w,*i));
		}
		for (int i=1;i<=comps_v.size();++i)
		for (int j=1;j<=comps_w.size();++j)
			result+=InverseMatrixEntry(i,j)*comps_v(i)*comps_w(j);
		return result.expand();
	}

/** @brief The musical isomorphism, \\alpha \\mapsto \\sum \langle \alpha, e^j\rangle e_j */
	virtual ex Sharp(ex form) const {
		if (sharp_subs==lst()) {
			exvector e_sharp;
			for (int i=1;i<=frame_.size();++i) {
				ex ei_sharp; 
				for (int j=1;j<=frame_.size();++j)
					ei_sharp+=MatrixEntry(i,j)*frame_.dual()(j);
				e_sharp.push_back(ei_sharp);
			}
			sharp_subs=LinearMapToSubstitutions<DifferentialOneForm>(frame_,e_sharp);
		}
		return form.subs(sharp_subs);
	}
/** @brief The musical isomorphism, v\\mapsto \\sum \langle v, e_j\rangle e^j */
	virtual ex Flat(ex vector) const {
		if (flat_subs==lst()) {
			exvector e_flat;
			for (int i=1;i<=frame_.size();++i) {
				ex ei_flat; 
				for (int j=1;j<=frame_.size();++j)
					ei_flat+=InverseMatrixEntry(i,j)*frame_(j);
				e_flat.push_back(ei_flat);
			}
			flat_subs=LinearMapToSubstitutions<DifferentialOneForm>(frame_.dual(),e_flat);
		}
		return vector.subs(flat_subs);
	}


protected:
/** @brief Returns the (i,j) entry of the matrix that represents the scalar product on one-forms wrt given frame*/
	virtual ex MatrixEntry(OneBased i, OneBased j) const=0;	///
/** @brief Returns the (i,j) entry of the matrix that represents the scalar product on vectors wrt given frame*/
	virtual ex InverseMatrixEntry(OneBased i, OneBased j) const=0;
	int DimensionOfSpace() const {return frame_.size();}
private:
	mutable lst sharp_subs, flat_subs;
};

/** @brief Riemannian scalar product represented by an orthonormal coframe
*/
class PositiveDefiniteScalarProduct : public BilinearFormWithFrame {
public:
/** @brief Define a Riemannian scalar product in terms of an orthonormal coframe
 @param frame An orthonormal coframe
*/
	PositiveDefiniteScalarProduct(const Frame& frame) : BilinearFormWithFrame(frame) {}
protected:
	ex MatrixEntry(OneBased i, OneBased j) const {
		return i==j ? 1 : 0;
	}
	ex InverseMatrixEntry(OneBased i, OneBased j) const {
		return i==j ? 1 : 0;
	}
};

/** @brief Pseudoriemannian scalar product represented by an orthonormal coframe
*/
class StandardScalarProduct : public BilinearFormWithFrame {
	int p_, q_;
public:
/** @brief Define a Pseudoriemannian scalar product in terms of an orthonormal coframe
* @param frame An orthonormal coframe
* @param p, q Integers defining the signature
*
* The orthonormal frame e_1,..., e_n is assumed to satisfy <e_i,e_j>=\delta_{ij} for i\leq p, -\delta_{ij} for p<i\leq p+q, and 0 otherwise.
*/
	StandardScalarProduct(const Frame& frame, int p, int q) : BilinearFormWithFrame(frame) {p_=p; q_=q;}
/** @brief Define a Pseudoriemannian scalar product in terms of an orthonormal coframe
* @param frame An orthonormal coframe
* @param p An integer defining the signature
*
* The orthonormal frame e_1,..., e_n is assumed to satisfy <e_i,e_j>=\delta_{ij} for i\leq p, -\delta_{ij} for i>p
*/
	StandardScalarProduct(const Frame& frame, int p) : BilinearFormWithFrame(frame) {p_=p; q_=frame.size()-p;}
protected:
	ex MatrixEntry(OneBased i, OneBased j) const {
		if (i<=p_) 
			return i==j? 1 : 0;
		else if (i<=p_+q_)
			return i==j? -1 : 0;
		else return 0;
	}
	ex InverseMatrixEntry(OneBased i, OneBased j) const {
		assert(p_+q_==DimensionOfSpace());	//only works if non degenerate
		return MatrixEntry(i,j);		//representative matrix coincides with its inverse
	}
};


/** @brief Neutral scalar product represented by a matrix of the form ((0 I) (I 0))
*/
class NeutralProduct : public BilinearFormWithFrame {
public:
	NeutralProduct(const Frame& frame) : BilinearFormWithFrame(frame) {assert(DimensionOfSpace()%2==0);}
protected:
	ex MatrixEntry(OneBased i, OneBased j) const {
		int n=DimensionOfSpace()/2;
		if (i-j==n || j-i==n) return 1; else return 0; 
	}
	ex InverseMatrixEntry(OneBased i, OneBased j) const {
		int n=DimensionOfSpace()/2;
		if (i-j==n || j-i==n) return 1; else return 0; 
	}

};


class ScalarProductDefinedByMatrix : public BilinearFormWithFrame {
	matrix m_, m_inverse_;
public:
	ScalarProductDefinedByMatrix(const Frame& frame, matrix m) : BilinearFormWithFrame(frame),m_{m}, m_inverse_{m.inverse()} {}
protected:
	ex MatrixEntry(OneBased i, OneBased j) const {
		return m_(i-1,j-1);
	}
	ex InverseMatrixEntry(OneBased i, OneBased j) const {
		return m_inverse_(i-1,j-1);
	}
};

}
#endif
