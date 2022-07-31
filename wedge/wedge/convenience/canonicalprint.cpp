/***************************************************************************
 *   Copyright (C) 2007-2022 by Diego Conti, diego.conti@unimib.it         *
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
#include "canonicalprint.h"
#include <regex>
namespace Wedge {
using namespace std;
using namespace GiNaC;
string to_string_using(const print_context* pc, ex x, int level) {
	stringstream s;
	unique_ptr<print_context> strpc{pc? pc->duplicate(s) : new print_dflt(s)};
	x.print(*strpc,level);
	return s.str();
}
string to_string_using(ostream& os, ex x, int level) {
	auto pc=get_print_context(os);
	return to_string_using(pc,x,level);
}
string to_canonical_string_using(const print_context& c, ex x, int level);

int precedence_for_printing(ex x) {
	return ex_to<basic>(x).precedence();
}
bool is_latex(const print_context& c) {
	return dynamic_cast<const print_latex*>(&c);
}

struct Term {
	string representation;
	string representation_for_comparison;	//a string extracted from the actual representation to be used for ordering
	string ncmul_factor;					//the factor of type ncmul inside the term, if it is a product, or ""
	int type_order;	//an integer depending on the expression's type, ensuring that muls and ncmuls come after other objects

	//if x is an ncmul or a product with an ncmul factor, return a representation of the ncmul. This ensures that
	//a term such as 2e^{12} comes before a term such as e^{34}.
	static string ncmul_within(ex x) {
		if (is_a<mul>(x)) {
			for (int i=0;i<x.nops();++i)
				if (is_a<ncmul>(x.op(i))) return to_canonical_string(x.op(i));
		}
		else if (is_a<ncmul>(x))
			return to_canonical_string(x);
		return{};
	}

	static string extract_essential(const string s) {
		int pos=0;
		if (s[pos]=='-') ++pos;
		if (s[pos]=='\\') ++pos;
		int power=s.find(pos,'^');
		int length=(power==string::npos)? power-pos : string::npos;
		return s.substr(pos,length);
	}

	static bool is_product(ex x) {
		return is_a<ncmul>(x) || is_a<mul>(x);
	}
	static int type_code(ex x) {
		return (is_product(x)&& is_product(-x))? 1: 0;
	}
	auto to_tuple_for_comparison() const {
		return std::tie(type_order,ncmul_factor,representation_for_comparison);
	}
public:
	Term(const print_context& c, ex x, int level) :
		representation {to_canonical_string_using(c,x,level)},
		representation_for_comparison{extract_essential(representation)},
		ncmul_factor{ncmul_within(x)},
		type_order{type_code(x)}
		{}
	bool operator<(const Term& other) const {
		return to_tuple_for_comparison()<other.to_tuple_for_comparison();
	}
	string rep() const {return representation;}
};

class Terms {
	set<Term> terms; //no two terms should have the same string representation
protected:
	virtual void print_first_term(const print_context& pc, const Term& term) const =0;
	virtual void print_term_after_first(const print_context& pc, const Term& term) const=0;
	int precedence;
public:
	Terms(const print_context& c, const expairseq& s) : precedence(s.precedence()) {
		for (int i=0;i<s.nops();++i)
			terms.emplace(c,s.op(i),s.precedence());
	}
	//TODO handle level
	void print(const print_context& pc,int level) const {
		if (precedence<level) pc.s<<"(";
		auto i=terms.begin();
		print_first_term(pc,*i);
		while (++i!=terms.end())
			print_term_after_first(pc,*i);
		if (precedence<level) pc.s<<")";
	}
};

class Product : public Terms {
protected:
	void print_first_term(const print_context& pc, const Term& term) const override {
		pc.s<<term.rep();
	}
	void print_term_after_first(const print_context& pc, const Term& term) const {
		char sep= is_latex(pc)? ' ' : '*';
		pc.s<<sep<<term.rep();
	}
public:
	using Terms::Terms;
};
class Sum : public Terms {
protected:
	void print_first_term(const print_context& pc, const Term& term) const override {
		pc.s<<term.rep();
	}
	void print_term_after_first(const print_context& pc, const Term& term) const {
		string s=term.rep();
		assert(s.length()>0);
		if (s[0]!='-') pc.s<<"+";
		pc.s<<term.representation;
	}
public:
	using Terms::Terms;
};

class CanonicalPrintVisitor : public visitor,public add::visitor, public mul::visitor, public power::visitor, public basic::visitor {
	const print_context& c;
	int level;

	static bool is_negative_numeric(ex x) {
		return  is_a<numeric>(x) && ex_to<numeric>(x).is_negative();
	}
	static bool is_positive_integer(ex x) {
		return  is_a<numeric>(x) && ex_to<numeric>(x).is_pos_integer();
	}
	static bool is_alpha_symbol(ex x) {
		static auto alpha=std::regex{"[a-zA-Z]*"};
		return is_a<symbol>(x) && regex_match(ex_to<symbol>(x).get_name(),alpha);
	}
	void print_latex_exponent(ex x) {
		if (is_positive_integer(x) || is_alpha_symbol(x)) c.s<<to_canonical_string_using(c,x,0);
		else c.s<<"{"<<to_canonical_string_using(c,x,0)<<"}";
	}
	void print_dflt_exponent(ex x) {
		static int power_precedence=power{2,3}.precedence();
		c.s<<to_canonical_string_using(c,x,power_precedence);
	}
	void print_square_root(ex argument) {
		if (is_latex(c)) c.s<<"\\sqrt{"<<to_canonical_string_using(c,argument,0)<<"}";
		else c.s<<"sqrt("<<to_canonical_string_using(c,argument,0)<<")";
	}
	void print_power(ex base, ex exponent) {
		static auto precedence=power(2,ex(1)/20).precedence();
		c.s<<to_canonical_string_using(c,base,precedence);
		c.s<<"^";
		if (is_latex(c)) print_latex_exponent(exponent);
		else print_dflt_exponent(exponent);
	}
public:
	CanonicalPrintVisitor(const print_context& c, int level): c{c},level{level} {}
	void visit(const mul& x) {
		if (is_negative_numeric(x.op(x.nops()-1))) {
			c.s<<"-";	//handle level
			(-x).normal().accept(*this);
		}
		else Product{c,x}.print(c,level);
	}

	void visit(const add& x) {
		Sum{c,x}.print(c,level);
	}
	void visit (const basic& x)	{
		int prev_level=level;
		level=x.precedence();
		x.print(c,prev_level);
	}
	void visit (const power& x) {
		if (x.op(1)==ex(1)/2) print_square_root(x.op(0));
		else print_power(x.op(0),x.op(1));
	}
};

const print_context& print_context_from_stream(ostream& os, unique_ptr<print_context>& store) {
	print_context* pc=get_print_context(os);
	if (pc) return *pc;
	store=make_unique<print_dflt>(os);
	return *store;
}

void canonical_print(ostream& os, ex x) {
	unique_ptr<print_context> newpc;
	auto& pc=print_context_from_stream(os,newpc);
	auto visitor=CanonicalPrintVisitor{pc,0};
	x.normal().accept(visitor);
}

string to_canonical_string_using(const print_context& c, ex x, int level) {
	stringstream s;
	unique_ptr<print_context> strpc{c.duplicate(s)};
	auto visitor=CanonicalPrintVisitor{*strpc,level};
	x.normal().accept(visitor);
	return s.str();
}

string to_canonical_string(ex x) {
	stringstream s;
	canonical_print(s,x);
	return s.str();
}
}
