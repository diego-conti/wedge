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

/** @brief Pseudoriemannian scalar product represented by an orthonormal coframe, without assuming that the space-like elements come first in the orthonormal basis
*/
class ScalarProductByOrthonormalFrame : public BilinearFormWithFrame {
	ExVector signs; 
	ScalarProductByOrthonormalFrame(const Frame& frame, const ExVector& signs) : BilinearFormWithFrame(frame), signs{signs} {}
public:
/** @brief Define a Pseudoriemannian scalar product in terms of an orthonormal coframe
* @param frame An orthonormal coframe
* @param timelike_indices a sequence of one-based indices corresponding to timelike elements in the orthonormal coframe; if empty, the metric is positive definite
*
* The orthonormal frame e_1,..., e_n is assumed to satisfy <e_i,e_j>=-\delta_ij or \delta_ij according to whether i is in timelike_indices
*/
	static ScalarProductByOrthonormalFrame FromTimelikeIndices(const Frame& frame, const list<int>& timelike_indices) {
		ExVector signs(frame.size());
		std::fill(signs.begin(),signs.end(),1);
		for (int i: timelike_indices)
			signs(i)=-1;
		return ScalarProductByOrthonormalFrame{frame,signs};		
	}
/** @brief Define a Pseudoriemannian scalar product in terms of an orthonormal coframe
* @param frame An orthonormal coframe e_1,...,e_n
* @param signs A sequence \epsilon_1,...,\epsilon_n, where each element is either 1 or -1
*
* The metric is taken to be of the form \epsilon_1 e^1\otimes e^1+...+\epsilon_n e^n\otimes e^n
*/
	static ScalarProductByOrthonormalFrame FromSequenceOfSigns(const Frame& frame, const vector<int>& signs) {
		return ScalarProductByOrthonormalFrame{frame,ExVector{signs.begin(),signs.end()}};
	}
	/**
   @brief Returns the list of timelike elements in the orthonormal frame
   @return A list of one-based indices corresponding to timelike elements in the orthonormal frame
 */
	vector<int> TimelikeIndices() const {
		vector<int> result;
		for (int i=1;i<=signs.size();++i)
			if (signs(i)<0) result.push_back(i);
		return result;
	}

protected:
	ex MatrixEntry(OneBased i, OneBased j) const override {
		return i==j? signs(i) : 0;
	}
	ex InverseMatrixEntry(OneBased i, OneBased j) const override {
		return i==j? signs(i) : 0;
	}
};

/** @brief Pseudoriemannian scalar product represented by an orthogonal coframe
*/
class ScalarProductByOrthogonalFrame : public BilinearFormWithFrame {
	ExVector g; 	//the sequence <e_1,e_1>, ... , <e_n,e_n>
	ScalarProductByOrthogonalFrame(const Frame& frame, const ExVector& g) : BilinearFormWithFrame(frame), g{g} {}
public:
/** @brief Define a Pseudoriemannian scalar product in terms of an orthonormal coframe
*  @param orthogonal_coframe An orthogonal coframe e^1,...,e^n with respect to which the metric is defined
*  @param g A sequence g_1,...,g_n such that the metric takes the form g_1e^1\otimes e^1+... + g_ne^n\otimes e^n
*/
	static ScalarProductByOrthogonalFrame FromVectorSquareNorms(const Frame& frame, const ExVector& g) {
		return ScalarProductByOrthogonalFrame{frame,g};		
	}
/** @brief Define a Pseudoriemannian scalar product in terms of an orthonormal coframe
*  @param orthogonal_coframe An orthogonal coframe e^1,...,e^n with respect to which the metric is defined
*  @param g the sequence <e^1,e^1>, ... , <e^n,e^n>
*/
	static ScalarProductByOrthogonalFrame FromFormSquareNorms(const Frame& frame, const ExVector& g) {
		ExVector g_inverse;
		g_inverse.reserve(g.size());
		std::transform(g.begin(),g.end(),back_inserter(g_inverse),[] (ex x) {return 1/x;});
		return ScalarProductByOrthogonalFrame{frame,g_inverse};		
	}
	/**
   @brief Returns the list of timelike elements in the orthogonal frame
   @return A list of one-based indices corresponding to timelike elements in the orthogonal frame
 */
	vector<int> TimelikeIndices() const {
		vector<int> result;
		for (int i=1;i<=g.size();++i)
			if (g(i)<0) result.push_back(i);
		return result;
	}

protected:
	ex MatrixEntry(OneBased i, OneBased j) const override {
		return i==j? 1/g(i) : 0;
	}
	ex InverseMatrixEntry(OneBased i, OneBased j) const override {
		return i==j? g(i) : 0;
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
	matrix metric_on_coframe,metric_on_frame;
	ScalarProductDefinedByMatrix(const Frame& frame, const matrix& metric_on_coframe, const matrix& metric_on_frame)
		: BilinearFormWithFrame(frame),metric_on_coframe{metric_on_coframe}, metric_on_frame{metric_on_frame} {}
public:
	static ScalarProductDefinedByMatrix OnCoframe(const Frame& frame, const matrix& metric_on_coframe) {
		return ScalarProductDefinedByMatrix{frame, metric_on_coframe,metric_on_coframe.inverse()};
	}
	static ScalarProductDefinedByMatrix OnFrame(const Frame& frame, const matrix& metric_on_frame) {
		return ScalarProductDefinedByMatrix{frame, metric_on_frame.inverse(),metric_on_frame};
	}
protected:
	ex MatrixEntry(OneBased i, OneBased j) const {
		return metric_on_coframe(i-1,j-1);
	}
	ex InverseMatrixEntry(OneBased i, OneBased j) const {
		return metric_on_frame(i-1,j-1);
	}
};

}
#endif
