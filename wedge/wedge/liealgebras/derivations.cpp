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

#include "derivations.h"
#include "wedge/representations/repgl.h"

namespace Wedge {
ex Xbracket(const LieGroup& G, const GLRepresentation<VectorField>& V, ex A, ex X, ex Y) {
	ex Ax=V.Action<VectorField>(A,X);
	ex Ay=V.Action<VectorField>(A,Y);
	ex Axy=V.Action<VectorField>(A,G.LieBracket(X,Y));
	return G.LieBracket(Ax,Y)+G.LieBracket(X,Ay)-Axy;
}

exvector Xbrackets(const LieGroup& G, const GLRepresentation<VectorField>& V, ex A) {
		exvector Xbrackets;
		for (int i=1;i<=G.Dimension();++i)
		for (int j=i+1;j<=G.Dimension();++j) 
			Xbrackets.push_back(Xbracket(G,V,A,G.e(i),G.e(j)));				
		return Xbrackets;
}


lst equations_such_that_linear_map_is_derivation(const LieGroup& G, const GL& gl, ex f) {
	lst eqns;
	auto X=Xbrackets(G,GLRepresentation<VectorField>(&gl,G.e()),f);		
	GetCoefficients<VectorField>(eqns,X);	
	return eqns;
}
	

VectorSpace<DifferentialForm> derivations(const LieGroup& G,const GL& Gl)  {
		auto gl=Gl.pForms(1);
		auto generic_matrix =gl.GenericElement();				
		auto eqns=equations_such_that_linear_map_is_derivation(G,Gl,generic_matrix);
		lst sol;		
		gl.GetSolutions(sol,eqns.begin(),eqns.end());		
		return {sol.begin(),sol.end()};
}

}