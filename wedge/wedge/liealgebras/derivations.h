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
#ifndef DERIVATIONS_H
#define DERIVATIONS_H
#include "liegroup.h"
#include "wedge/representations/gl.h"
#include "wedge/representations/linearaction.h"

namespace Wedge {

/** @brief Return a lst of equations that an element f of gl, acting on G with the action induced by G.e(), must satisfy in order to be a derivation of G */
lst equations_such_that_linear_map_is_derivation(const LieGroup& G, const GL& gl, ex f);

/** @brief Return the space of derivations of n n-dimensional Lie algebra without parameters as a subspace of gl(n,R) */
VectorSpace<DifferentialForm> derivations(const LieGroup& G,const GL& Gl);


}
#endif