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
#ifndef CANONICAL_PRINT_H
#define CANONICAL_PRINT_H
#include "wedge/base/wedgebase.h"

namespace Wedge {

/** @brief Canonical output of a Ginac ex, with terms sorted in a consistent order
 * Notice that no simplification is performed, since there is no "universal" simplification command in GiNaC.
 * Thus, (a+b)*x and a*x+b*x will give different outputs.
 * The caller should invoke the appropriate simplifiation command, e.g. canonical_print(x.normal());
 */
void canonical_print(ostream& os, ex x);
string to_canonical_string(ex x);
string to_latex_canonical_string(ex x);

string to_string_using(const print_context* pc, ex x, int level=0);
string to_string_using(ostream& os, ex x, int level=0);

}
#endif
