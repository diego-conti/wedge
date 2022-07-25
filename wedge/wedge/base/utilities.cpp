/***************************************************************************
 *   Copyright (C) 2007, 2008, 2009 by Diego Conti, diego.conti@unimib.it  *
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

#include "wedgebase.h"
#include "utilities.h"
//#include "vectorspace.h"

namespace Wedge {
 using namespace  GiNaC; 

bool IterateOverPermutations::Iterate(int n, int startFromPermutationNumber, int endAtPermutationNumber )
{
	assert(n>0);
	std::vector<int> Value(n);
	for (int i=0;i<n;i++) Value[i]=i;		

	int todo=endAtPermutationNumber-startFromPermutationNumber;
	while (startFromPermutationNumber-->0 && next_permutation(Value.begin(),Value.end())) ;

	if (endAtPermutationNumber<0)
	{
		bool result;
		do {
			result=Apply(Value);
		}
		while (result && next_permutation(Value.begin(),Value.end())) ;		
		return result;
	}
	else {
		while (Apply(Value) && todo-->0 && next_permutation(Value.begin(),Value.end())) ;
		return todo<0;
	}         
}
	
void IterateOverSubsets::Iterate(int m, int n)
{
	if (m<0 || m>n ) throw OutOfRange(__FILE__,__LINE__,m);
	else if (m==0) return;
	vector<int> choice(m);
	int l=0;
	choice[0]=0;
	while(choice[0]!=n) {
		if (l==m-1) {
			if (!Apply(choice)) return;
			++choice[l];
		}
		else {
			choice[l+1]=choice[l]+1;
			l++;
		}
		if (choice[l]==n) {
			while (l>0 && ++choice[--l]==n) ;
		}
	}
}

namespace internal {
class RealPartVisitor : public visitor, public  add::visitor, public mul::visitor, public ncmul::visitor, public basic::visitor, public symbol::visitor, public  numeric::visitor {
	ex real,imag;
	friend void Wedge::SplitIntoRealAndImaginary(ex e,ex& realpart,ex& imaginarypart) ;
public:	
	void visit (const add& e) {
		ex rsum=0,isum=0;
		for (int i=0;i<e.nops();i++)
		{
			e.op(i).accept(*this);
			rsum+=real; isum+=imag;
		}
		real=rsum; imag=isum;
	}
	void visit (const mul& e) {
		ex rprod=1,iprod=0;
		for (int i=0;i<e.nops();i++) {
			e.op(i).accept(*this);
			ex old_rprod=rprod;
			rprod=rprod*real-iprod*imag;
			iprod=old_rprod*imag+iprod*real;
		}
		real=rprod; imag=iprod;
	}
	void visit (const ncmul& e) {
		ex rprod=1,iprod=0;
		for (int i=0;i<e.nops();i++) {
			e.op(i).accept(*this);
			ex old_rprod=rprod;
			rprod=rprod*real-iprod*imag;
			iprod=old_rprod*imag+iprod*real;
		}
		real=rprod; imag=iprod;
	}
	void visit (const basic& e) {real=e; imag=0;}	///< @deprecated Unsafe and obsolete
	void visit (const symbol& e) {
		if (e.get_domain()!=domain::complex)  {
			real=e; 
			imag=0;
		}
		else throw InvalidArgument(__FILE__,__LINE__,e);
	}

	void visit (const numeric& e) {real=e.real(); imag=e.imag();}
};

}

void SplitIntoRealAndImaginary(ex e,ex& realpart,ex& imaginarypart) 
{
	internal::RealPartVisitor v;
	e.accept(v);
	realpart=v.real;
	imaginarypart=v.imag;
}

namespace internal {
ex DumbNumericRoot(numeric x, numeric n)
{
	numeric coeff=1;
	numeric i=2;
	numeric isq;
	while ((isq=pow(i,n))<x)
	{
		while (mod(x,isq)==0)
		{			
			x/=(isq);
			coeff*=i;
		}
		i++;
	}
	return coeff*pow(ex(x),ex(1)/n);
}
}

#if (GINACLIB_MAJOR_VERSION<=1) && (GINACLIB_MINOR_VERSION<=4)
typedef lst Matches;
#else
typedef exset Matches;
#endif

ex NormalizeRoots(ex expression)
{	
	//x^m * y^n -> x^
	expression=expression.expand();
	Matches matches;
	lst substitutions;
	if (expression.find(pow(wild(1),wild(2)),matches))
	{
		numeric exponent=0;		
		for (Matches::const_iterator i=matches.begin();i!=matches.end();i++)
		{
			assert(is_a<power>(*i));
			if (is_a<numeric>(i->op(1))) {				
				numeric n=ex_to<numeric>(i->op(1));				
				if (n.is_rational())
				{
					exponent =
						exponent.is_zero()? n :
							1/lcm(n.numer()*exponent.denom(),n.denom()*exponent.numer())*(n.numer()*exponent.numer());					
				}									
			}
		}
		for (Matches::const_iterator i=matches.begin();i!=matches.end();i++)
		{			
			if (is_a<numeric>(i->op(1))) {
				numeric thisexponent=ex_to<numeric>(i->op(1));
				assert((thisexponent/exponent).is_integer());	
				ex newbase=pow(i->op(0),thisexponent/exponent);
				substitutions.append(*i==pow(newbase,exponent));				
			}
		}		
	}
	expression=expression.subs(substitutions).expand();
	
	// x^z * y^z -> (xy)^z
	ex newexpr= expression.subs(pow(wild(1),wild(2))*pow(wild(3),wild(2))==pow(wild(1)*wild(3),wild(2)),subs_options::algebraic).expand();
	while (newexpr!=expression)
	{
		expression=newexpr;
	 	newexpr= expression.subs(pow(wild(1),wild(2))*pow(wild(3),wild(2))==pow(wild(1)*wild(3),wild(2)),subs_options::algebraic).expand();
	}
	
	//handle numeric roots
	expression=expression.expand();
	substitutions.remove_all();
	matches=Matches();
	if (expression.find(power(wild(1),wild(2)),matches))
	{
		for (Matches::const_iterator i=matches.begin();i!=matches.end();i++)
		{
			assert(is_a<power>(*i));
			if (is_a<numeric>(i->op(0)) && is_a<numeric>(i->op(1))) {
				numeric n=ex_to<numeric>(i->op(0));
				numeric m=ex_to<numeric>(i->op(1));				
				if (n.is_rational() && m.is_rational())
				{
					ex numer=power(internal::DumbNumericRoot(n.numer(),m.denom()),m.numer());
					ex denom=power(internal::DumbNumericRoot(n.denom(),m.denom()),m.numer());
					substitutions.append(*i==numer/denom);
				}									
			}
		}		
	}
	expression=expression.subs(substitutions).expand();

	// (x^y)^z -> x^{yz}
	newexpr= expression.subs(pow(pow(wild(1),wild(2)),wild(3))==pow(wild(1),wild(2)*wild(3)),subs_options::algebraic).expand();
	while (newexpr!=expression)
	{
		expression=newexpr;
		newexpr= expression.subs(pow(pow(wild(1),wild(2)),wild(3))==pow(wild(1),wild(2)*wild(3)),subs_options::algebraic).expand();
	}

	return expression;	
}



}
