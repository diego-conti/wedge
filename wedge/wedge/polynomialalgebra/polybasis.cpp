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
 #include "wedge/polynomialalgebra/polybasis.h"

namespace Wedge {

namespace internal {

ex SimplifyPowerEqn(ex x)
{
	return is_a<power>(x)? x.op(0) : x;
}

ex SimplifyPolyEqn(ex x,lst variables)
{
	x=sqrfree(x.expand(),variables);	//bug in ginac requires expand before sqrfree??
	if (is_a<mul>(x)) {
		exvector ops;
		for (ex y : x) ops.push_back(SimplifyPowerEqn(y));
		return mul(ops);
	}
	else return SimplifyPowerEqn(x);
}

}
}
