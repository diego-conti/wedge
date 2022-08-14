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
#ifndef PARSE_H
#define PARSE_H

/** @ingroup Base */ 
/** @file parse.h
 * @brief Parsing functions 
 */
namespace Wedge {
class ParseError : public WedgeException<runtime_error> {
public:
	ParseError(string parsing,const char* file, int line) : WedgeException<runtime_error> ("Error parsing: "+parsing,file,line) {}
	ParseError(int index,const char* file, int line) : WedgeException<runtime_error> ("Out of bounds in parsing: "+internal::ToString(index),file,line) {}
};

	/** @brief Convert a string into a list of forms  
	 * @param frame A reference frame, in the guise of a vector of one-forms 
	 * @param to_parse A string of comma-separated expressions, each representing a form in the given frame
	 * @return An ExVector of forms
	 * 
	 * @remark The notation generalizes that of 
	 * [S. Salamon :Complex structures on nilpotent Lie algebras. J. Pure Appl. Algebra 157 (2001), no. 2-3, 311--333]
	 * and is referred to as \em Salamon's \em notation in this documentation.
	 *  
	 *  A few examples (where \f$e^1\dotsc e^5\f$ denote the elements of the reference frame):
 	*  - 1   -> \f$e^1\f$
 	*  - 1+2  -> \f$e^1+e^2\f$
 	*  - 4*1+sqrt(2)*45  -> \f$4e^1+\sqrt2e^4\wedge e^5\f$
 	* 
 	*  @note This notation will not work for frames of length greater than 9.
	 */
	ExVector ParseDifferentialForms(const exvector& frame, const char* to_parse);
	/** @overload
	 */
	ExVector ParseDifferentialForms(const exvector& frame, const char* to_parse, ex symbols);
	
	/** @overload
	 */
	ex ParseDifferentialForm(const exvector& frame, const char* to_parse);
	
	/** @overload
	 */
	ex ParseDifferentialForm(const exvector& frame, const char* to_parse, ex symbols);
	
	/** @brief Convert a string into a list of expressions  
	 * @param symbols A lst of symbols that appear in the expression 
	 * @param to_parse A string of comma-separated expressions, each representing a maple expression
	 * @return An ExVector
	 */
	ExVector ParseMapleExpressions(ex symbols, const char* to_parse);
	
	/** @overload
	 */
	ex ParseMapleExpression(ex symbols, const char* to_parse);
	
	/** @brief Convert a string into a list of expressions  
	 * @param symbols A lst of symbols that appear in the expression 
	 * @param to_parse A string of comma-separated expressions, each representing a %CoCoA expression
	 * @return An ExVector
	 */
	ExVector ParseCocoaExpressions(const exvector& symbols, const char* to_parse);
	
	/** @overload
	 */	
	ex ParseCocoaExpression(const exvector& symbols, const char* to_parse);
}
#endif
