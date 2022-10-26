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
#include "../structures/pseudoriemannianstructure.h"

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

pair<ex,matrix> PseudoRiemannianStructureByOrthonormalFrame::DecomposeRicci(matrix ricci) const  {
	throw NotImplemented(__FILE__,__LINE__,"PseudoRiemannianStructureByOrthonormalFrame::DecomposeRicci");
}

ex PseudoRiemannianStructureByOrthonormalFrame::u(ZeroBased n) const {
	int m=M()->Dimension()/2;
	if (n<0) throw OutOfRange(__FILE__,__LINE__,n);
	vector<int> signs(m,1);
	signs.reserve(m);
	for (int i=0;i<m;++i, n/=2)
		if (n%2) signs[m-i-1]=-1;	//take into account reversal
	if (n!=0) throw OutOfRange(__FILE__,__LINE__,n);
	return Spinor::from_epsilons(signs);
}

ex PseudoRiemannianStructureByOrthonormalFrame::u(const vector<int>& signs) const {
	int m=M()->Dimension()/2;
	if (signs.size()!=m) throw InvalidArgument(__FILE__,__LINE__,lst{list<ex>(signs.begin(),signs.end())});
	return Spinor::from_epsilons(signs);
}


int PseudoRiemannianStructureByOrthonormalFrame::DimensionOfSpinorRepresentation() const
{
	return 1<<(M()->Dimension()/2);
}

class CliffordProduct : public IBilinearOperator<LinearOperator<VectorField>,LinearOperator<Spinor>> {
	Frame orthonormal_coframe;
	ExVector taus;	//tau_k=i if e_k is timelike and 1 if spacelike
	int s;	//the second element of the signature (r,s), i..e the number of taus that equal i

	ex alpha_j(OneBased j) const {
		if (j%2==0) return 1;
		else if (j==orthonormal_coframe.size() && s%2) return -I;
		else return I;
	}

	ex dot(OneBased j, const Spinor& spinor) const {
		ex coeff=taus(j)*spinor.product_up_to(j/2)*alpha_j(j);
		if (((j-1)/2)%2) coeff=-coeff;
		if (j==taus.size() && j%2==1)  return coeff*spinor;
		else return coeff*spinor.reflect((j+1)/2);
	}
public:
	CliffordProduct(const Frame& frame, const vector<int>& timelike_indices) : orthonormal_coframe{frame}, taus(frame.size(),1), s{timelike_indices.size()} {
		//std::fill(taus.begin(),taus.end(),ex{1});
		for (int i: timelike_indices) 
			taus(i)=I;		
	}
	ex Apply (const VectorField& X, const Spinor& spinor) const
	{
		ex result;		
		for (int i=1;i<=orthonormal_coframe.size();++i) {
			ex component = TrivialPairing<VectorField>(X,orthonormal_coframe(i));			
			if (!component.is_zero())
				result+=component*dot(i,spinor);
		}
		return result;
	}	
};


ex PseudoRiemannianStructureByOrthonormalFrame::CliffordDot(ex X, ex psi) const {
		CliffordProduct clifford{e(), scalar_product.TimelikeIndices()};
		return CliffordProduct::BilinearOperator(X,psi,&clifford);
}

}