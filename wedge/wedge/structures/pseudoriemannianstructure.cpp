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
#include "../structures/pseudoriemannianstructure.h"
#include "pseudoriemannianstructure.h"

namespace Wedge {

pair<ex,matrix> StandardPseudoRiemannianStructure::DecomposeRicci(matrix ricci) const  {
	int dim=M()->Dimension();
	assert(ricci.rows()==dim);
	assert(ricci.cols()==dim);
	ex s;
	for (int i=0;i<p_;++i)
		s+=ricci(i,i);
	for (int i=p_;i<dim;++i)
		s-=ricci(i,i);
	ex normalized_s=s/dim;
	for (int i=0;i<p_;++i)
		ricci(i,i)-=normalized_s;
	for (int i=p_;i<dim;++i)
		ricci(i,i)+=normalized_s;
	return {normalized_s,ricci};
}


pair<ex,matrix> PseudoRiemannianStructureByFrame::DecomposeRicci(matrix ricci) const  {
	throw NotImplemented(__FILE__,__LINE__,"PseudoRiemannianStructureByOrthonormalFrame::DecomposeRicci");
}

ex PseudoRiemannianStructureByFrame::u(ZeroBased n) const {	
	return Spinor::from_index_and_dimension(n,M()->Dimension());
}

ex PseudoRiemannianStructureByFrame::u(const vector<int>& signs) const {
	int m=M()->Dimension()/2;
	if (signs.size()!=m) throw InvalidArgument(__FILE__,__LINE__,lst{list<ex>(signs.begin(),signs.end())});
	return Spinor::from_epsilons(signs);
}


int PseudoRiemannianStructureByFrame::DimensionOfSpinorRepresentation() const
{
	return 1<<(M()->Dimension()/2);
}

class PseudoRiemannianStructureByFrame::CliffordProduct : public IBilinearOperator<LinearOperator<VectorField>,LinearOperator<Spinor>> {
	Frame coframe;
	ExVector taus;	//tau_k=i if e_k is timelike and 1 if spacelike
	int s;	//the second element of the signature (r,s), i..e the number of taus that equal i

	ex alpha_j(OneBased j) const {
		int r=taus.size()-s;
		if (j%2==0) return 1;
		else if (j==coframe.size() && (s-r+1)%4) return -I;
		else return I;
	}

	ex dot(OneBased j, const Spinor& spinor) const {
		ex coeff=taus(j)*spinor.product_up_to(j/2)*alpha_j(j);
		if (((j-1)/2)%2) coeff=-coeff;
		if (j==taus.size() && j%2==1)  return coeff*spinor;
		else return coeff*spinor.reflect((j+1)/2);
	}
	CliffordProduct(const Frame& coframe, const ExVector& taus, int s) : coframe{coframe}, taus{taus},s{s} {}
public:
	static CliffordProduct* FromTimelikeIndices(const Frame& coframe, const vector<int>& timelike_indices) {
		ExVector taus(coframe.size(),1);
		for (int i: timelike_indices) 
			taus(i)=I;		
		return new CliffordProduct{coframe, taus,timelike_indices.size()};
	}
	static CliffordProduct* FromSquareNormsOfFrameVectors(const Frame& coframe, const ExVector& square_norms_of_frame_vectors) {
		ExVector taus;
		int s=0;
		for (int i=0;i<coframe.size();++i) {
			taus.push_back(sqrt(square_norms_of_frame_vectors[i]));
			if (square_norms_of_frame_vectors[i]<0) ++s;		}
		return new CliffordProduct(coframe,taus,s);
	}
	ex Apply (const VectorField& X, const Spinor& spinor) const
	{
		ex result;		
		for (int i=1;i<=coframe.size();++i) {
			ex component = TrivialPairing<VectorField>(X,coframe(i));			
			if (!component.is_zero())
				result+=component*dot(i,spinor);
		}
		return result;
	}	
};

class PseudoRiemannianStructureByFrame::CliffordProductForm : public IBilinearOperator<AssociativeOperator<DifferentialForm>,LinearOperator<Spinor>> {
	const PseudoRiemannianStructureByFrame::CliffordProduct& clifford;
	const PseudoRiemannianStructure& g;
public:
	CliffordProductForm(const PseudoRiemannianStructureByFrame::CliffordProduct& clifford, const PseudoRiemannianStructure& g) : clifford{clifford}, g{g} {}
	ex Apply (const VectorField& alpha, const Spinor& psi) const {
		auto X=g.ScalarProduct().Sharp(alpha);
		return PseudoRiemannianStructureByFrame::CliffordProduct::BilinearOperator(X,psi,&clifford);
	}
};


PseudoRiemannianStructureByFrame::PseudoRiemannianStructureByFrame(const Manifold *manifold, const Frame &frame, CliffordProduct *clifford_product)
	: PseudoRiemannianStructure(manifold,frame), clifford_product_operator{clifford_product}, clifford_product_form_operator{new CliffordProductForm(*clifford_product_operator,*this)} {}

void PseudoRiemannianStructureByFrame::Deleter::operator() (CliffordProduct* r) {
	delete r;
}
void PseudoRiemannianStructureByFrame::Deleter::operator() (CliffordProductForm* r) {
	delete r;
}

ex PseudoRiemannianStructureByFrame::CliffordDot(ex X, ex psi) const {
	return CliffordProduct::BilinearOperator(X,psi,clifford_product_operator.get());
}
ex PseudoRiemannianStructureByFrame::CliffordDotByForm(ex alpha, ex psi) const {
	return CliffordProductForm::BilinearOperator(alpha,psi,clifford_product_form_operator.get());
}

PseudoRiemannianStructureByOrthonormalFrame::PseudoRiemannianStructureByOrthonormalFrame(const Manifold* manifold, const Frame& frame, ScalarProductByOrthonormalFrame&& scalar_product) :
	PseudoRiemannianStructureByFrame(manifold,
		frame,
		CliffordProduct::FromTimelikeIndices(frame, scalar_product.TimelikeIndices())
	), 
	scalar_product{std::move(scalar_product)}
	{}

PseudoRiemannianStructureByOrthogonalFrame::PseudoRiemannianStructureByOrthogonalFrame(const Manifold *manifold, const Frame &orthogonal_coframe, const ExVector &g) :
	PseudoRiemannianStructureByFrame(
		manifold,
		orthogonal_coframe,
		CliffordProduct::FromSquareNormsOfFrameVectors(orthogonal_coframe, g)
	),
	scalar_product{ScalarProductByOrthogonalFrame::FromVectorSquareNorms(orthogonal_coframe,g)}
	 {}

}
