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
#include "liederivative.h"


namespace Wedge {
 using namespace  GiNaC;
 using namespace std;

 ex LieDerivative::Canonicalize(list<VectorField>& v, const Function& f, list<VectorField>::iterator at, const Manifold& M)
 {
	 assert(!v.empty());
	 list<VectorField>::iterator next=at; ++next;
	 //if at is the last element, nothing to do
	 if (next==v.end()) return LieDerivative(v,f);

	 VectorField X=*at, Y=*next;
	 //if the sequence v is not sorted, use the Lie bracket to obtain a sum of sorted sequences
	 if (X.compare_same_type(Y)>0)
	 {
		 *at=Y; *next=X;
		 ex result=Canonicalize(v,f,next,M);
		 *at=X;	*next=Y;	//restore v to its original state

		 ex t;
		 if (++next==v.end())
			 t=f;
		 else
			 t=LieDerivative(list<VectorField>(next,v.end()), f);
		 LOG_INFO(X);
		 LOG_INFO(Y);
		 t=M.LieDerivative(M.LieBracket(X,Y),t);
		 LOG_INFO(t);
		 for (list<VectorField>::const_iterator i=at;i!=v.begin();)
			 t=M.LieDerivative(*--i,t);
		 return result+t;
	 }
	 //if the sequence is sorted, return it as it is
	 else
		 return LieDerivative(v,f);
 }

}
