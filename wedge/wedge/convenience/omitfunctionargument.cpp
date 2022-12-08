/*******************************************************************************
 *  Copyright (C) 2007-2022 by Diego Conti, diego.conti@unipi.it 
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
#include "wedge/base/wedgebase.h"
#include "wedge/manifolds/fderivative.h"
#include "wedge/convenience/omitfunctionargument.h"
#include "wedge/convenience/printcontext.h"

namespace Wedge {
using namespace std;

print_omit_function_argument::Init print_omit_function_argument::init;

GINAC_IMPLEMENT_PRINT_CONTEXT(print_omit_function_argument, print_latex)

print_omit_function_argument::print_omit_function_argument() : print_latex(std::cout) {}

//had to patch ginac/function.cpp to make this work. the problem is that ginac::function's behaviour cannot be overridden with set_print_func
void print_function_omit_argument(const GiNaC::function & p,
		const print_omit_function_argument & c,
		unsigned level)
{
	// get the precedence of the 'function' class
	unsigned function_prec = p.precedence();

	// if the parent operator has the same or a higher precedence
	// we need parentheses
	if (level >= function_prec)
		c.s << '(';

	c.s<<p.get_latex_name();

	if (!c.ShouldOmit(p)) {
		c.s<<'('<<p.op(0);
		for (int i=1;i<p.nops();++i)
			c.s<<","<<p.op(i);
		c.s<<')';
	}

	if (level >= function_prec)
		c.s << ')';

}


void print_overall_coeff(const mul& p, const print_omit_function_argument & c,unsigned level)
{
	ex overallcoeff=p.op(p.nops()-1);
	if (!is_a<numeric>(overallcoeff)) return;
	//overall coefficient is the last operand
	const numeric &coeff = ex_to<numeric>(overallcoeff);
	if (coeff.csgn() == -1)
		c.s << '-';
	if (!coeff.is_equal(1) &&
			!coeff.is_equal(-1)) {
		if (coeff.is_rational()) {
			if (coeff.is_negative())
				(-coeff).print(c);
			else
				coeff.print(c);
		} else {
			if (coeff.csgn() == -1)
				(-coeff).print(c, p.precedence());
			else
				coeff.print(c, p.precedence());
		}
		//overallcoefficient is a numeric so no omission here
	}
}




void print_mul_omit_argument(const mul & p,	const print_omit_function_argument & c,unsigned level)
{
	if (p.precedence() <= level)
		c.s << "{(";


	// Separate factors into those with negative numeric exponent
	// and all others
	exvector neg_powers, others;

	for (int i=0;i<p.nops()-1;++i) {
		ex op=p.op(i);
		if (is_a<power>(op) && op.op(1)<0) neg_powers.push_back(op);
		else others.push_back(op);
	}
	ex last_coefficient=p.op(p.nops()-1);
	//numeric factors and sign are absorbed in the last coefficient
	if (!is_a<numeric>(last_coefficient))  {
		if (is_a<power>(last_coefficient) && last_coefficient.op(1)<0) neg_powers.push_back(last_coefficient);
		else others.push_back(last_coefficient);
	}
	else {
		const numeric &coeff = ex_to<numeric>(last_coefficient);
		if (coeff.csgn() == -1) {c.s<<"-"; last_coefficient=-last_coefficient;}
		if (last_coefficient!=1) c.s<<last_coefficient;
	}

	if (!neg_powers.empty()) {
		// Factors with negative exponent are printed as a fraction
		c.s << "\\frac{";
		mul(others).eval().print(c);
		c.s << "}{";
		(1/mul(neg_powers)).eval().print(c);
		c.s << "}";
	} else {
		// All other factors are printed in the ordinary way
		c.s<<' ';
		others[0].print(c,p.precedence());
		for (int i=1;i<others.size();++i)
		{
			//if this term is a function with omitted argument and the next term is in parentheses, write \\cdot
			if (c.ShouldOmit(others[i-1])) {
				const basic& thisterm=ex_to<basic>(others[i]);
				if (thisterm.precedence()<=p.precedence()) c.s<<"\\cdot";
			}
			c.s<<' ';
			others[i].print(c, p.precedence());
		}
	}
	if (p.precedence() <= level)
		c.s << ")}";

}



void print_fderivative_omit_argument(const fderivative & p,
		const print_omit_function_argument & c,
		unsigned level)
{
	// get the precedence of the 'function' class
	unsigned function_prec = p.precedence();

	// if the parent operator has the same or a higher precedence
	// we need parentheses
	if (level >= function_prec)
		c.s << '(';

	if (level >= function_prec)
		c.s << '(';

	vector<int> multiindex=FunctionDerivativeToMultiIndex(p);

	if (c.ShouldOmit(p)) {

		c.s<<p.get_latex_name();

		assert(multiindex.size()==1);
		assert(multiindex[0]>0);
		switch (multiindex[0]) {
		case 0: assert(false);
		case 1: c.s<<"'"; break;
		case 2: c.s<<"''"; break;
		case 3: c.s<<"'''"; break;
		default:
			c.s<<"^{("<<multiindex[0]<<")}";
		}
	}
	else {
		int sum=0;
		for (vector<int>::const_iterator k=multiindex.begin();k!=multiindex.end();++k)
			sum+=*k;
		assert(sum>0);
		if (sum>1)
			c.s<<"\\frac{\\partial^{"<<sum<<"}"<<p.get_latex_name()<<"}";
		else
			c.s<<"\\frac{\\partial "<<p.get_latex_name()<<"}";
		assert(multiindex.size()==p.nops());

		c.s<<"{";
		for (int k=0;k<multiindex.size();++k)
			if (multiindex[k]==1) c.s<<"\\partial_"<<k+1;	//use one-based indices
			else if (multiindex[k]>1) c.s<<"\\partial_"<<k+1<<"^{"<<multiindex[k]<<"}";
		c.s<<"}";
		c.s<<'('<<p.op(0);
			for (int i=1;i<p.nops();++i)
				c.s<<","<<p.op(i);
			c.s<<')';
	}

	if (level >= function_prec)
		c.s << ')';
}

std::ostream& operator<<(std::ostream& os, OmitArgument x)
{
	set_print_context(os,print_omit_function_argument(os,x));
	return os;
}

print_omit_function_argument::Init::Init()
{
	// install a new print method for function objects
	set_print_func<GiNaC::function, print_omit_function_argument>(print_function_omit_argument);
	set_print_func<mul, print_omit_function_argument>(print_mul_omit_argument);
	set_print_func<fderivative, print_omit_function_argument>(print_fderivative_omit_argument);
}

}
