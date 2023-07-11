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
 /*
 * fderivative.h
 *
 *  Created on: Apr 27, 2015
 *      Author: diego
 */

#ifndef WEDGE_FDERIVATIVE_H_
#define WEDGE_FDERIVATIVE_H_

#include <ginac/fderivative.h>

namespace Wedge {
using namespace GiNaC;
using namespace std;

/** Return the multiindex associated to a fderivative object.
 *
 * @param f an fderivative object
 * @result a vector with a number of elements equal to the number of arguments of f whose values
 * are the order of derivation
 *
 * e.g. D[0,0]f(x,y,z) is mapped to the vector [2,0,0]
 */
vector<int> FunctionDerivativeToMultiIndex(const fderivative& f);

/** Create an fderivative object from a multiindex, a serial number, and a list of arguments */
ex FunctionDerivativeFromMultiIndex(unsigned serial, const vector<int>& orders, const exvector& args);


}
#endif /* WEDGE_FDERIVATIVE_H_ */
