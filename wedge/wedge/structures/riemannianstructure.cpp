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
#include "../structures/riemannianstructure.h"

#include "../structures/spinor.h"

namespace Wedge {

ex RiemannianStructure::HodgeStar(ex form1) const
{
	ex form=e()[0];
	for (int i=1;i<e().size();i++)
		form*=e()[i];
	return Hook(form1,form);
}

ex RiemannianStructure::Hook(ex left, ex right) const
{
	ex result=internal::RiemannianHookOperator::BilinearOperator(left.expand(),right.expand(),&hookOperator).expand();
		//adjust sign so that e^{12..n}\hook e^{12..n}=1
	return (Degree<DifferentialForm>(left) % 4<2)? result : -result;
}	

template<> ex RiemannianStructure::ScalarProduct<VectorField> (ex op1, ex op2) const {
	return ScalarProduct().OnVectors(op1,op2);
}
template<> ex RiemannianStructure::ScalarProduct<DifferentialForm> (ex op1, ex op2) const {
	return ScalarProduct().OnForms(op1,op2);
}

template<> ex RiemannianStructure::ScalarProduct<Spinor> (ex op1, ex op2) const {
	return RealPart(TrivialPairing<Spinor>(op1,op2.conjugate()));
}

namespace internal {
ex RiemannianHookOperator::Apply(const VectorField& left,const VectorField& right) const
{
	return structure->ScalarProduct().OnOneForms(left,right);
}
}

}
