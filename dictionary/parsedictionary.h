/***************************************************************************
 *   Copyright (C) 2008 by Diego Conti, diego.conti@unimib.it              *
 *                                                                         *
 *   This file is part of Dictionary.                                      *
 *                                                                         *
 *   Dictionary is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 3 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   Dictionary is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with Dictionary; if not, write to the                           *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#ifndef PARSEDICTIONARY_H_
#define PARSEDICTIONARY_H_
/** @file parsedictionary.h
 * @brief Code to construct an element of a given dictionary directly from a string
 */

#include <boost/spirit/include/classic_core.hpp>
#include <boost/spirit/include/classic_dynamic.hpp>
#include <string>
#include <list>
#include <wedge/wexception.h>
#include <wedge/wedgealgebraic.h>
#include "dictionary.h"

using namespace std;
using namespace boost::spirit::classic;
using namespace Wedge;

template<typename Dictionary> ex 
	ParseDictionaryElement(const Dictionary& dictionary, const char* to_parse, lst symbols);
/** @overload
*/
template<typename Dictionary> ex 
	ParseDictionaryElement(const Dictionary& dictionary, const char* to_parse)
{
	return ParseDictionaryElement(dictionary,to_parse,lst());
}

namespace Dictionary_internal {
template<typename Dictionary>
	struct DictionaryParser : public grammar<DictionaryParser<Dictionary> >
{
	struct Data {
		const Dictionary& dictionary; 
		lst symbols;
		Data(const Dictionary& _dictionary,lst _symbols) : dictionary(_dictionary),symbols(_symbols) {}
		string contractionname;		
		list<typename Dictionary::NamedVValuedForm> forms;
		list<ex> parsing_stack;
	};

	Data* data;	//need a pointer here to work around the fact that the parser is passed as a const reference

	struct ContrAction  {
		Data& data;
		ContrAction(const DictionaryParser<Dictionary>& parser) : data(*parser.data) {}
		void operator() (char const& a) const {
			data.parsing_stack.push_back(data.dictionary.Contract(data.contractionname,data.forms));
 			data.forms.clear();			
 			data.contractionname.erase();
		}
	};

	struct ProductAction  {
		Data& data;
		ProductAction(const DictionaryParser& parser) : data(*parser.data) {}
		void operator() (char const* a, char const* b) const {
			ex last=data.parsing_stack.back();
			data.parsing_stack.pop_back();
			data.parsing_stack.back()*=last;
		}
	};

	struct ContractionName  {
		Data& data;
		ContractionName(const DictionaryParser& parser) : data(*parser.data) {}
		void operator() (char const* a, char const* b) const {
			data.contractionname=string(a,b-a);
		}
	};

	struct SumAction {
		Data& data;
		SumAction(const DictionaryParser& parser) : data(*parser.data) {}
		void operator() (char const* a, char const* b) const {
			ex op2=data.parsing_stack.back();
			data.parsing_stack.pop_back();
			data.parsing_stack.back()+=op2;
		}
	};

	struct DifferenceAction  {
		Data& data;
		DifferenceAction(const DictionaryParser& parser) : data(*parser.data) {}
		void operator() (char const* a, char const* b) const {
			ex op2=data.parsing_stack.back();
			data.parsing_stack.pop_back();
			data.parsing_stack.back()-=op2;
		}
	};

	struct NegateAction {
		Data& data;
		NegateAction(const DictionaryParser& parser) : data(*parser.data) {}
		void operator() (char const* a, char const* b) const {
			data.parsing_stack.back()=-data.parsing_stack.back();
		}
	};

	struct GinacAction  {
		Data& data;
		GinacAction(const DictionaryParser& parser) : data(*parser.data) {}
		void operator() (const char* a, const char* b) const {			
			data.parsing_stack.push_back(ex(string(a,b-a),data.symbols));
		}
	};


	struct PowerAction  {
		Data& data;
		PowerAction(const DictionaryParser& parser) : data(*parser.data) {}
		void operator() (const char* a, const char* b) const {
			ex exponent=data.parsing_stack.back();
			data.parsing_stack.pop_back();
			data.parsing_stack.back()=pow(data.parsing_stack.back(),exponent);
		}
	};

	struct LetterAction  {
		Data& data;
		typename Dictionary::NamedVValuedForm letter;
		LetterAction(const DictionaryParser& parser,const typename Dictionary::NamedVValuedForm& _letter) : data(*parser.data),letter(_letter) {}
		void operator() (const char* a, const char* b) const {			
			data.forms.push_back(letter);
		}
	};	

	struct FormAction  {
		Data& data;
		NamedForm form;
		FormAction(const DictionaryParser& parser,const NamedForm& _form) : data(*parser.data),form(_form) {}
		void operator() (const char* a, const char* b) const {			
			data.parsing_stack.push_back(form);
		}
	};	

	template <typename ScannerT>
			struct definition
	{
		
		definition(const DictionaryParser& self)
		{	
			expression =
					(('-' >> term)[NegateAction(self)] | term)
					>> *(   ('+' >> term)[SumAction(self)]
					|   ('-' >> term)[DifferenceAction(self)]
						)
					;

			term = ('['>>form>>']' | ginacexpression | power) >>  *('['>>form>>']' | ginacexpression | power)[ProductAction(self)];

			power = contraction>> !((ch_p('^')>>exponent)[PowerAction(self)]);

			exponent = (ginac_in_parentheses | (lexeme_d[+(alnum_p)]))[GinacAction(self)];
			
			ginacexpression=(
					+ (ginac_in_parentheses | (anychar_p-')'-'('-'{'-'['-'+'-'-'))
					)[GinacAction(self)];

			ginac_in_parentheses = (ch_p('(')>> + (ginac_in_parentheses | (anychar_p-')'-'('))>>ch_p(')')); 

			contraction=ch_p('{')>>(letter | (*~ch_p(','))[ContractionName(self)])>> ','>>(letter % ch_p(','))>>ch_p('}')[ContrAction(self)];                                                                                                
			
			list<typename Dictionary::NamedVValuedForm> letters=self.data->dictionary.Letters();
			for (typename list<typename Dictionary::NamedVValuedForm>::iterator i=letters.begin();i!=letters.end();i++)
			{				
				letter=letter.copy()|(str_p(i->name().c_str())[LetterAction(self,*i)]);			
			}
			
			list<NamedForm> f=self.data->dictionary.FormsOnBase();
			for (list<NamedForm>::iterator i=f.begin();i!=f.end();i++)
			{				
				form=form.copy()|(str_p(i->name().c_str())[FormAction(self,*i)]);			
			}
		}

		rule<ScannerT> ginacexpression,contraction,term,expression,exponent,power,ginac_in_parentheses;
		stored_rule<ScannerT> letter,form;

		rule<ScannerT> const& start() const { return expression; }
	};

	DictionaryParser(const Dictionary& dictionary,lst symbols) {
		this->data=new Data(dictionary,symbols);
	}
	
	~DictionaryParser() {delete data;}
};
}

/** @brief Convert a string representing an element of the dictionary to a differential form
 * @param dictionary The dictionary to use
 * @param to_parse The string to parse
 * @param symbols Symbols appearing in the expression; they must have different names, and none may be called r.
 * @return A CompositeElement corresponding to the specified dictionary element
 * 
 * The syntax is as following:
 * - If x is a form on the base, the string "[x]" represents that form
 * - If a and b are letters, the string {a,b} means the contraction of a and b
 * - If a and b are letters and sigma is a contraction taking two elements, {sigma,a,b} means the contraction of a and b through sigma
 * - The product of two forms is indicated by juxtaposition
 * - A form can be multiplied by a coefficient, which is parsed like the constructor ex::ex(const string&)
 * - Parentheses are allowed in the coefficients, but the forms cannot appear enclosed in parentheses. For instance, (1+sqrt(3))*2 {a,a} is a valid expression, but 2({a,a}+{a,b}) is not.
 * - r denotes the radial coordinate
 *
 * Symbols passed as arguments are converted to type symbol in the result.
*/
template<typename Dictionary> ex 
	ParseDictionaryElement(const Dictionary& dictionary, const char* to_parse, lst symbols)
{
	for (lst::const_iterator i=symbols.begin();i!=symbols.end();++i)
		if (!is_a<symbol>(*i)) throw InvalidArgument(__FILE__,__LINE__,*i);
		else {
			const symbol& s=ex_to<symbol>(*i);
			if (s.get_name()=="r") throw WedgeException<std::invalid_argument>("Symbol r is reserved and cannot be passed as an argument",__FILE__,__LINE__);
		}
	symbol r("r");
	symbols.append(r);
	Dictionary_internal::DictionaryParser<Dictionary> p(dictionary,symbols);
	bool result=parse(to_parse,p,space_p).full;
	if (!result) throw WedgeException<runtime_error>("Error parsing "+string(to_parse),__FILE__,__LINE__);
	return dictionary.MakeGlobal(p.data->parsing_stack.back().subs(r==dictionary.RadialCoordinate()));
}

#endif /*PARSEDICTIONARY_H_*/
