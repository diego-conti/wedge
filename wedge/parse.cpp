/***************************************************************************
 *   Copyright (C) 2014 by Diego Conti, diego.conti@unimib.it  *
 *                                                                         *
 *   This file is part of Wedge.                                           *
 *                                                                         *
 *   Wedge is free software; you can redistribute it and/or modify         *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 3 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   Wedge is distributed in the hope that it will be useful,              *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with Wedge; if not, write to the                                *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#include "parse.h"
#include <ginac/ginac.h>
#include <string>
#include "wexception.h"
#include "wedgealgebraic.h"
#include <iostream>
#include "spiritgrammars.hpp"

/**
	 @brief Convert a string into a list of forms  
	 * @param frame A reference frame, in the guise of a vector of one-forms 
	 * @param to_parse A string of comma-separated expressions, each representing a form in the given frame
	 * @return An ExVector of forms
	 * 
	 * @remark Suppose the frame is e^1,...,e^{15}. The following expressions are valid:
 	*  - 1   -> \f$e^1\f$
 	*  - 12   -> \f$e^1\wedge e^2\f$
	*  - 0	-> \f$0\f$
 	*  - -1+2 -> \f$-e^1+e^2\f$
	*  - 2*2 -> \f$2e^2\f$. The first two is interpreted as a constant because it is followed by a *.  
	*  - 2/3*2 -> \f$\frac23e^2\f$. The first coefficient is interpreted as a fraction.  
	*  - a + 1 -> \f$ e^{10}+e^1 \f$. a is hexadecimal 10 (the digits are 0, ..., 9, a, ..., f; 0 is not used). 
	*  - [sqrt(2)]*45  -> \f$\sqrt2e^4\wedge e^5\f$. The expression between square brackets is interpreted as a ginac expression
	*  - (1+2)^(3+2*4), or  (1+2)(3+2*4)  -> \f$(e^1+e^2)\wedge(e^3+2e^4)\f$
 	*  - 4*1+sqrt(2)*45  -> \f$4e^1+\sqrt2e^4\wedge e^5\f$
	*  - 1 2 -> \f$e^1 \wedge e^2\f$
	*  - [a^(1/3)]*2 -> \f$a^(1/3)e^2\f$. The symbol a must be passed as an argument for this to be handled correctly.
	*
	* More formally: a form is the sum of terms: if x, y are terms, then x, -x, x-y, x+y, are forms.
	* A term is either 0 or a product of factors. If x, y are factors, x y or x^y is a form representing \f$ x\wedge y \f$
	* If x is a form, (x) is a factor
	* If x is either a GiNaC expression in square brackets or a positive rational number followed by a *, then x is a constant
	* If x is a constant and y is a factor, x*y is a factor representing the product
	* If x is a sequence of digits between 1 and min(15,frame.size()) not followed by a *, x is a factor representing a form, namely the product of the corresponding one-forms in the frame
 	*  @note This notation will not work for frames of length greater than 15.

 * @remark The notation generalizes that of 
	 * [S. Salamon :Complex structures on nilpotent Lie algebras. J. Pure Appl. Algebra 157 (2001), no. 2-3, 311--333]
	 * and is referred to as \em Salamon's \em notation in this documentation.
*/



namespace Wedge {

//Namespace Parser

GiNaC::symtab SymtabFromSymbols(ex symbols) {
	symtab table;
	if (is_a<symbol>(symbols)) {
		string name=ex_to<symbol> (symbols).get_name(); 
		table[name]=symbols;
	}
	else {
		assert(is_a<lst>(symbols));
		for (int i=0;i<symbols.nops();++i) {
			assert(is_a<symbol>(symbols.op(i)));
			string name=ex_to<symbol> (symbols.op(i)).get_name(); 
			table[name]=symbols.op(i);
		}
	}
	return table;
}


const char* endOfString(const char* begin)
{
	const char* result=begin;
	while (*result) ++result; 
	return result;
}


ex ParseDifferentialForm(const exvector& reference_coframe, const string& to_parse, ex symbols) 
{
	return DifferentialFormParser{reference_coframe,SymtabFromSymbols(symbols)}.ParseForm(to_parse.begin(),to_parse.end());
}
ex ParseDifferentialForm(const exvector& reference_coframe, const char* to_parse, ex symbols) 
{
	return DifferentialFormParser{reference_coframe,SymtabFromSymbols(symbols)}.ParseForm(to_parse,endOfString(to_parse));
}
ex ParseDifferentialForm(const exvector& reference_coframe, const string& to_parse) 
{
	return DifferentialFormParser{reference_coframe}.ParseForm(to_parse.begin(),to_parse.end());
}
ex ParseDifferentialForm(const exvector& reference_coframe, const char* to_parse) 
{
	return DifferentialFormParser{reference_coframe}.ParseForm(to_parse,endOfString(to_parse));
}

ExVector ParseDifferentialForms(const exvector& reference_coframe, const string& to_parse, ex symbols) 
{
	return DifferentialFormParser{reference_coframe,SymtabFromSymbols(symbols)}.ParseForms(to_parse.begin(),to_parse.end());
}
ExVector ParseDifferentialForms(const exvector& reference_coframe, const char* to_parse, ex symbols) 
{
	return DifferentialFormParser{reference_coframe,SymtabFromSymbols(symbols)}.ParseForms(to_parse,endOfString(to_parse));
}
ExVector ParseDifferentialForms(const exvector& reference_coframe, const string& to_parse) 
{
	return DifferentialFormParser{reference_coframe}.ParseForms(to_parse.begin(),to_parse.end());
}
ExVector ParseDifferentialForms(const exvector& reference_coframe, const char* to_parse) 
{
	return DifferentialFormParser{reference_coframe}.ParseForms(to_parse,endOfString(to_parse));
}



ex ParseMapleExpression(ex symbols,const char* to_parse) {
	MapleParser p(SymtabFromSymbols(symbols));
	return p.ParseExpression(to_parse,endOfString(to_parse));
}

ExVector ParseMapleExpressions(ex symbols,const char* to_parse) {
	MapleParser p(SymtabFromSymbols(symbols));
	return p.ParseExpressions(to_parse,endOfString(to_parse));
}

ex ParseMapleExpression(ex symbols,string to_parse) {
	MapleParser p(SymtabFromSymbols(symbols));
	return p.ParseExpression(to_parse.begin(),to_parse.end());
}

ExVector ParseMapleExpressions(ex symbols, string to_parse) {
	MapleParser p(SymtabFromSymbols(symbols));
	return p.ParseExpressions(to_parse.begin(),to_parse.end());
}

ex ParseCocoaExpression(const exvector& symbols, string to_parse) {
	CocoaParser p(symbols);
	return p.ParsePolynomial(to_parse.begin(),to_parse.end());
}
ExVector ParseCocoaExpressions(const exvector& symbols, string to_parse) {
	CocoaParser p(symbols);
	return p.ParsePolynomials(to_parse.begin(),to_parse.end());
}
ex ParseCocoaExpression(const exvector& symbols,const char* to_parse) {
	CocoaParser p(symbols);
	return p.ParsePolynomial(to_parse,endOfString(to_parse));
}
ExVector ParseCocoaExpressions(const exvector& symbols,const char* to_parse) {
	CocoaParser p(symbols);
	return p.ParsePolynomials(to_parse,endOfString(to_parse));
}


}

