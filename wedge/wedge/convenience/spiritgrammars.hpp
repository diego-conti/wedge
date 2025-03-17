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

//the part of the parsing code that depends on boost::spirit.

#include <boost/spirit/home/x3.hpp>
#include <boost/config/warning_disable.hpp>
#include <boost/fusion/include/at.hpp>
#include "spiritsemanticactions.hpp"

namespace x3 = boost::spirit::x3;
namespace ascii = boost::spirit::x3::ascii;
namespace fusion = boost::fusion;
namespace Wedge {

GiNaC::symtab SymtabFromSymbols(ex symbols);

template<typename TypeToDisplay> class DisplayType;


/** Template class to parse a Maple expression or equation or a list of Maple expressions or equations, depending on whether template parameter Result
 is ex or ExVector

@throws std::invalid_argument if the string contains a symbol that does not appear in the symbols list
*/

class FormParser {
	static const exvector* reference_coframe_;
	ex OneFormFromIndex(int index) const {
		if (index<=0 || index>reference_coframe_->size()) throw ParseError(index,__FILE__,__LINE__);
		return reference_coframe_->at(index-1);
	}
	ex OneFormFromDigit(char digit) const {
		if (digit>'0' && digit<='9') return OneFormFromIndex(digit-'0');
		else if (digit>='a' && digit<='z') return OneFormFromIndex(digit-'a'+10);
		else if (digit>='A' && digit<='Z') return OneFormFromIndex(digit-'A'+10+26);
		else throw ParseError(digit,__FILE__,__LINE__);
	}
	template<typename Container>
	ex Parse(Container digit_sequence) const {
		auto i=digit_sequence.begin();
		ex result=OneFormFromDigit(*i);
		while (++i!=digit_sequence.end()) result*=OneFormFromDigit(*i);
		return result;
	}
public:
	static void SetReferenceCoframe(const exvector* reference_coframe) noexcept {
		 reference_coframe_=reference_coframe;
	}
	template<typename Context> void operator() (const Context& ctx) const {
		x3::_val(ctx)=Parse(x3::_attr(ctx));
	}
};

class ExpressionParser {
	static unique_ptr<GiNaC::parser> ginac_parser;
	static bool is_symtab_trivial;
	ex Parse(const string& char_sequence) const {return (*ginac_parser)(char_sequence);}
public:
	static void SetSymbolTable(const symtab& symbol_table) {
		is_symtab_trivial=false;
	 	ginac_parser = make_unique<GiNaC::parser>(symbol_table,true);
	}
	static void ResetSymbolTable() {
		if (!is_symtab_trivial) {
		 	ginac_parser = make_unique<GiNaC::parser>(symtab(),true);
			is_symtab_trivial=true;
		}
	}
	template<typename Context> void operator() (const Context& ctx) const {
		x3::_val(ctx)=Parse(x3::_attr(ctx));
	}
};

class CocoaVariableParser {
	static const exvector* variables_;
public:
	static void SetVariables(const exvector* variables) {
		variables_=variables;
	}
	template<typename Context> void operator() (const Context& ctx) const {
		auto index = x3::_attr(ctx);
		x3::_val(ctx)=variables_->at(index);
	}
};


const exvector* FormParser::reference_coframe_=nullptr;
unique_ptr<GiNaC::parser> ExpressionParser:: ginac_parser = make_unique<GiNaC::parser>(symtab(),true);
const exvector* CocoaVariableParser::variables_=nullptr;
bool ExpressionParser::is_symtab_trivial=true;


namespace DifferentialFormGrammar {
	using namespace SemanticActions;
	ExpressionParser expression_parser;
	FormParser form_parser;

	x3::rule<class integer_id, unsigned int>	//declare an incomplete type as template argument, used by X3 framework as an ID; notice that this will create problems if a class with that name is defined elsewhere
	        integer = "Integer";
	x3::rule<class ginac_expression_id, ex>
	        ginac_expression = "GinacExpression";
	x3::rule<class constant_id, ex>
	        constant = "Constant";
	x3::rule<class simple_form_id, ex>
	        simple_form = "SimpleForm";
	x3::rule<class factor_id, ex>
	        factor = "Factor";
	x3::rule<class product_id, ex>
	        product = "Product";
	x3::rule<class quotient_id, ex>
	        quotient = "Quotient";
	x3::rule<class wedge_product_id, ex>
	        wedge_product = "WedgeProduct";
	x3::rule<class zero, ex>
	        zero = "Zero";
	x3::rule<class term_id, ex>
	        term = "Term";
	x3::rule<class negative_term_id, ex>
	        negative_term = "NegativeTerm";
	x3::rule<class signed_term_id, ex>
	        signed_term = "SignedTerm";
	x3::rule<class form_id, ex>
	        form = "Form";
	x3::rule<class forms_id, ExVector>
	        forms = "Forms";

	auto integer_def = x3::uint_;
	auto ginac_expression_def='[' >> ((*(ascii::char_ - ']'))[expression_parser]) >> ']';
	auto quotient_def =  (integer >> '/'> integer)[read_quotient];
	auto constant_def = ginac_expression | quotient | integer;
	auto simple_form_def=(+x3::alnum)[form_parser];
	auto product_def = (constant >>  '*' > factor)[read_product];
	auto factor_def = product |  simple_form |  ('('> form > ')');
	auto wedge_product_def=	(factor % -x3::lit('^'))[read_wedge_product];	//the wedge is optional
	auto zero_def=x3::lit('0') [read_zero];
	auto term_def = zero_def | wedge_product_def;
	auto negative_term_def = ( '-' > term )[read_opposite];
	auto signed_term_def = negative_term | term;
	auto form_def =	(signed_term % -x3::lit('+'))[read_sum];
	auto forms_def = (form % ',') > x3::eoi;

	BOOST_SPIRIT_DEFINE(integer,ginac_expression,constant,simple_form,product, factor,term,signed_term,form,forms,zero,wedge_product,negative_term,quotient);	
}



class DifferentialFormParser {
public:
	DifferentialFormParser(const exvector& reference_coframe)
	{
		FormParser::SetReferenceCoframe(&reference_coframe);
	}
	DifferentialFormParser(const exvector& reference_coframe, const symtab& symbol_table)
	{
		ExpressionParser::SetSymbolTable(symbol_table);
		FormParser::SetReferenceCoframe(&reference_coframe);
	}
	~DifferentialFormParser() noexcept {
		ExpressionParser::ResetSymbolTable();
		FormParser::SetReferenceCoframe(nullptr);
	}

	template<typename Iterator> ex ParseForm(Iterator begin, Iterator end)
	{
		ex result;
		if (!x3::phrase_parse(begin,end,DifferentialFormGrammar::form,ascii::space,result))			
			throw ParseError("Error parsing differential form",__FILE__,__LINE__);
		return result;
	}
	template<typename Iterator> ExVector ParseForms(Iterator begin, Iterator end)
	{
		ExVector result;
		if (!x3::phrase_parse(begin,end,DifferentialFormGrammar::forms,ascii::space,result))			
			throw ParseError("Error parsing differential forms",__FILE__,__LINE__);
		return result;
	}
};


namespace CocoaParserGrammar {
	using namespace SemanticActions;
	CocoaVariableParser cocoa_variable_parser;

	x3::rule<class variable_id, ex>
		variable = "Variable";
	x3::rule<class power_id, ex>
		power = "Power";
	x3::rule<class fraction_id, ex>
		fraction = "Fraction";
	x3::rule<class monomial_id, ex>
		monomial = "Monomial";
	x3::rule<class number_id, ex>
		number = "Number";
	x3::rule<class signed_number_id, ex>
		signed_number = "SignedNumber";
	x3::rule<class elem_id, ex>
		elem = "Elem";
	x3::rule<class term_id, ex>
		term = "PositiveTerm";
	x3::rule<class signed_term_id, ex>
		signed_term = "SignedTerm";
	x3::rule<class polynomial_id, ex>
		polynomial = "Polynomial";
	x3::rule<class equation_id, ex>
		equation = "Equation";
	x3::rule<class equation_or_polynomial_id, ex>
		equation_or_polynomial = "EquationOrPolynomial";
	x3::rule<class list_id, ExVector>
		list = "List";
	x3::rule<class negative_term_id, ex>
	        negative_term = "NegativeTerm";

	auto variable_def = "x[" > x3::uint_[cocoa_variable_parser]   > ']';
	auto power_def = (variable >> ('^'> x3::uint_ ))[read_power];
	auto monomial_def = power | variable;
	auto fraction_def = (x3::uint_ >>'/' > x3::uint_)[read_quotient] ;
	auto number_def = fraction | x3::uint_;
	auto signed_number_def = number | (x3::lit('(') > '-' > number > ')')[read_opposite];
	auto elem_def = ('(' > polynomial > ')') | number | monomial;
	auto term_def  = (elem % '*')[read_mul];
	auto negative_term_def = ( '-' > term )[read_opposite];
	auto signed_term_def = negative_term | term;
	auto polynomial_def = -x3::lit('+') > (signed_term % -x3::lit('+'))[read_sum];
	auto equation_def = (polynomial >> '=' > polynomial)[read_equation];
	auto equation_or_polynomial_def = equation | polynomial;
	auto list_def = (equation_or_polynomial_def % ',') > x3::eoi; 

	BOOST_SPIRIT_DEFINE(variable, monomial,number,signed_number,elem,term,signed_term,polynomial,equation,equation_or_polynomial,list,power,fraction,negative_term);
}



class CocoaParser {
public:
	CocoaParser(const exvector& symbols){
		CocoaVariableParser::SetVariables(&symbols);
	}
	~CocoaParser() {
		CocoaVariableParser::SetVariables(nullptr);
	}

	template<typename Iterator> ex ParsePolynomial(Iterator begin, Iterator end)
	{
		ex result;
		try {
			if (!x3::phrase_parse(begin,end,CocoaParserGrammar::equation_or_polynomial,ascii::space,result))			
			throw ParseError("Error parsing Cocoa polynomial "+string{begin,end},__FILE__,__LINE__);
		}
		catch (...) {
			throw ParseError("Error parsing Cocoa polynomial "+string{begin,end},__FILE__,__LINE__);
		}
		return result;
	}

	template<typename Iterator> ExVector ParsePolynomials(Iterator begin, Iterator end) {
		ExVector result;
		try {
			if (!x3::phrase_parse(begin,end,CocoaParserGrammar::list,ascii::space,result))			
			throw ParseError("Error parsing Cocoa polynomial",__FILE__,__LINE__);
		}
		catch (...) {
			throw ParseError("Error parsing Cocoa polynomial "+string{begin,end},__FILE__,__LINE__);
		}
		return result;
	}
};



namespace MapleParserGrammar {
	using namespace SemanticActions;

	ExpressionParser expression_parser;

	x3::rule<class expression_id, ex>
	        expression = "Expression";
	x3::rule<class equation_id, ex>
		equation = "Equation";
	x3::rule<class expression_or_equation_id, ex>
	        expression_or_equation = "ExpressionOrEquation";
	x3::rule<class list_id, ExVector>
	        list = "List";

	auto expression_def = (+(ascii::char_ - ascii::char_(",="))) [expression_parser];
	auto equation_def = (expression >> '=' > expression)[read_equation];
	auto expression_or_equation_def = equation | expression;
	auto list_def = (expression_or_equation % ',') >x3::eoi;
	BOOST_SPIRIT_DEFINE(expression,equation,expression_or_equation,list);
}


class MapleParser {
public:
	MapleParser(const symtab& symbol_table)
	{
		ExpressionParser::SetSymbolTable(symbol_table);
	}

	~MapleParser() {
		ExpressionParser::ResetSymbolTable();
	}

	template<typename Iterator> ex ParseExpression(Iterator begin, Iterator end)
	{
		ex result;
		if (!x3::phrase_parse(begin,end,MapleParserGrammar::expression_or_equation,ascii::space,result))			
			throw ParseError("Error parsing Maple expression",__FILE__,__LINE__);
		return result;
	}

	template<typename Iterator> ExVector ParseExpressions(Iterator begin, Iterator end) {
		ExVector result;
		if (!x3::phrase_parse(begin,end,MapleParserGrammar::list,ascii::space,result))			
			throw ParseError("Error parsing Maple expressions",__FILE__,__LINE__);
		return result;
	}

};

}
