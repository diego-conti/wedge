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
#include "wedge/base/wedgebase.h"
#include "wedge/base/wedgealgebraic.h"
#include "wedge/convenience/parse.h"
#include "wedge/convenience/simplifier.h"
#include "wedge/convenience/canonicalprint.h"
#include <cxxtest/TestSuite.h>
#include "wedge/linearalgebra/lambda.h"
#include "wedge/manifolds/concretemanifold.h"
#include "test.h"
#include <regex>

using namespace Wedge;
using namespace GiNaC;

WEDGE_DECLARE_NAMED_ALGEBRAIC(V,Vector)

WEDGE_DECLARE_NAMED_ALGEBRAIC(W,Vector)

WEDGE_DECLARE_NAMED_ALGEBRAIC(X,symbol)


//test the parsing functions
class ParseTestSuite : public CxxTest::TestSuite 
{
public:
	void testParseGinac() {
		//test if GiNaC's parser works correctly on square roots
		//(there is a subtle bug in parsers/default_reader.cpp which relies on functions being aligned, which is architecture-dependent)
		ex p("sqrt(3)"s,lst{});
		TS_ASSERT_EQUALS(p*p,3);
		symbol a("a");
		ex q("sqrt(a)",lst{a});
		TS_ASSERT_EQUALS(q*q,a);
	}
	void testParseMaple()
	{
		symbol a("a"),b("b");
		TS_ASSERT_THROWS(ParseMapleExpression(lst(),"a*2"), std::invalid_argument);
		TS_ASSERT_EQUALS(ParseMapleExpression(a,"a*2"), a*2);
		TS_ASSERT_EQUALS(ParseMapleExpression(a,"a*2=a"), a*2==a);
		ExVector expressions=ParseMapleExpressions(lst{a,b},"a^2,b-sqrt(a)");
		TS_ASSERT_EQUALS(expressions[0], a*a);
		TS_ASSERT_EQUALS(expressions[1], b-sqrt(a));

	}


	void testParseCocoa()
	{
		logging_level=LOGGING_LEVEL_DEBUG;
		char s[]="x[1],x[1]+x[2], -x[1]+x[2],-2*x[1]-x[2],-2*x[1]^2+x[2]^2*x[0],3*x[1]*(x[1]-x[2]),1/2*x[1],x[0]^2 +1/2*x[0]";
		exvector x(10);
		for (int i=0;i<x.size();i++) x[i]=symbol("x["+ToString(i)+"]");
		exvector v=ParseCocoaExpressions(x,s);
		TS_ASSERT_EQUALS(v[0],x[1]);
		TS_ASSERT_EQUALS(v[1],x[1]+x[2]);
		TS_ASSERT_EQUALS(v[2],-x[1]+x[2]);
		TS_ASSERT_EQUALS(v[3],-2*x[1]-x[2]);
		TS_ASSERT_EQUALS(v[4],-2*pow(x[1],2)+pow(x[2],2)*x[0]);
		TS_ASSERT_EQUALS(v[5],3*x[1]*(x[1]-x[2]));
		TS_ASSERT_EQUALS(v[6],x[1]/2);
		TS_ASSERT_EQUALS(v[7],pow(x[0],2)+x[0]/2);
	}
	void testParseForms()
	{
		ConcreteManifold M(13);		
		TS_ASSERT_EQUALS(ParseDifferentialForm(M.e(),"[sqrt(2)]*1"),sqrt(ex(2))*M.e(1));		
		TS_ASSERT_EQUALS(ParseDifferentialForm(M.e(),"0"),0);
		
		exvector forms=ParseDifferentialForms(M.e(),"-1,2+[sqrt(3)]*34+0,123+7*4,-3-(1-2),1-2/3*(-3+45)");
		TS_ASSERT_EQUALS(-M.e(1),forms[0]);
		TS_ASSERT_EQUALS(M.e(2)+sqrt(ex(3))*M.e(3)*M.e(4),forms[1]);
		TS_ASSERT_EQUALS(M.e(1)*M.e(2)*M.e(3)+7*M.e(4),forms[2]);
		TS_ASSERT_EQUALS(-M.e(3)-M.e(1)+M.e(2),forms[3]);
		TS_ASSERT_EQUALS(M.e(1)-numeric(2,3)*(-M.e(3)+M.e(4)*M.e(5)),forms[4]);

		symbol a("a");
		forms=ParseDifferentialForms(M.e(),"1+[a]*a4,1/3*52^[1+sqrt(2)]*3,[pow(2,-1/3)]*(1+34),[a^(-2)]*(1+2)3,2*3+4*[sqrt(3)]*1",lst{a});
		TS_ASSERT_EQUALS(M.e(1)+a*M.e(10)*M.e(4),forms[0]);
		TS_ASSERT_EQUALS(1/ex(3)*M.e(5)*M.e(2)*(1+sqrt(ex(2)))*M.e(3),forms[1]);
		TS_ASSERT_EQUALS(pow(2,-1/ex(3))*(M.e(1)+M.e(3)*M.e(4)),forms[2]);
		TS_ASSERT_EQUALS(1/(a*a)*(M.e(1)+M.e(2))*M.e(3),forms[3]);
		TS_ASSERT_EQUALS(2*M.e(3)+4*sqrt(ex(3))*M.e(1),forms[4]);
		TS_ASSERT_THROWS(ParseDifferentialForm(M.e(),"[a]*2"), std::invalid_argument);
		TS_ASSERT_EQUALS(ParseDifferentialForm(M.e(),"-12+[a]*34+5*(6+7)+8(9+a)",a), -M.e(1)*M.e(2)+a*M.e(3)*M.e(4)+5*(M.e(6)+M.e(7))+M.e(8)*(M.e(9)+M.e(10)));
	}
	void testParseFormsBounds() {
		ConcreteManifold N(3);
		TS_ASSERT_THROWS(ParseDifferentialForm(N.e(),"4"),ParseError);
		TS_ASSERT_THROWS(ParseDifferentialForm(N.e(),"1+4"),ParseError);
		TS_ASSERT_THROWS(ParseDifferentialForm(N.e(),"14"),ParseError);
		TS_ASSERT_THROWS(ParseDifferentialForm(N.e(),"2*4"),ParseError);
		ConcreteManifold M(13);
		TS_ASSERT_THROWS(ParseDifferentialForm(M.e(),"e"),ParseError);
		TS_ASSERT_THROWS(ParseDifferentialForm(M.e(),"e1"),ParseError);
		TS_ASSERT_THROWS(ParseDifferentialForm(M.e(),"1e"),ParseError);
		TS_ASSERT_THROWS(ParseDifferentialForm(M.e(),"cg"),ParseError);
		TS_ASSERT_THROWS(ParseDifferentialForm(M.e(),"gc"),ParseError);
		TS_ASSERT_THROWS(ParseDifferentialForm(M.e(),"1e"),ParseError);
		TS_ASSERT_THROWS(ParseDifferentialForm(M.e(),"1+e"),ParseError);
	}
	void testParseFormsLarge() {
		ConcreteManifold M(61);
		TS_ASSERT_EQUALS(ParseDifferentialForm(M.e(),"a"),M.e(10));
		TS_ASSERT_EQUALS(ParseDifferentialForm(M.e(),"z"),M.e(35));
		TS_ASSERT_EQUALS(ParseDifferentialForm(M.e(),"A"),M.e(36));
		TS_ASSERT_EQUALS(ParseDifferentialForm(M.e(),"Z"),M.e(61));
		TS_ASSERT_THROWS(ParseDifferentialForm(M.e(),"$"),ParseError);
		TS_ASSERT_THROWS(ParseDifferentialForm(M.e(),"."),ParseError);
		TS_ASSERT_THROWS(ParseDifferentialForm(M.e(),"\\"),ParseError);
	}
	void testToStringUsing() {
		symbol psi("psi");
		stringstream s;
		TS_ASSERT_EQUALS(get_print_context(s),nullptr);
		TS_ASSERT_EQUALS(to_string_using(s,psi),"psi");
		s<<latex;
		TS_ASSERT_EQUALS(to_string_using(s,psi),"\\psi");
		s<<dflt;
		TS_ASSERT_EQUALS(to_string_using(s,psi),"psi");
		power(psi,2).print(*make_unique<print_dflt>(s));
		TS_ASSERT_EQUALS(s.str(),"psi^2");
		TS_ASSERT_EQUALS(to_string_using(s,power(psi,2)),"psi^2");
	}
	void testCanonicalPrintDflt() {
		symbol a1("a1"), a2("a2"), psi("psi"),a("a");

		TS_ASSERT_EQUALS(to_canonical_string(-a),"-a");
		TS_ASSERT_EQUALS(to_canonical_string(-1/a),"-a^(-1)");
		TS_ASSERT_EQUALS(to_canonical_string(-a*psi),"-a*psi");
		TS_ASSERT_EQUALS(to_canonical_string(a2-a),"-a+a2");
		TS_ASSERT_EQUALS(to_canonical_string(-a2+a),"a-a2");
		TS_ASSERT_EQUALS(to_canonical_string(2*a2-2*a),"-2*a+2*a2");

		TS_ASSERT_EQUALS(to_canonical_string(-a+a1*a1-a2-psi),"-a+a1^2-a2-psi");
		TS_ASSERT_EQUALS(to_canonical_string(1/a1*a*psi*a/a1),"a1^(-2)*a^2*psi");
		TS_ASSERT_EQUALS(to_canonical_string(a1/a),"a1*a^(-1)");
		TS_ASSERT_EQUALS(to_canonical_string((a+a1)/a2),"(a+a1)*a2^(-1)");
		TS_ASSERT_EQUALS(to_canonical_string(a*a1),"a*a1");
		TS_ASSERT_EQUALS(to_canonical_string(a*a1+a),"a+a*a1");

		TS_ASSERT_EQUALS(to_canonical_string(sqrt(ex(2))*a1+a),"a+sqrt(2)*a1");

	}
	string to_latex_canonical_string(ex x) {
		stringstream s;
		s<<latex;
		canonical_print(s,x);
		return s.str();
	}

	void testCanonicalPrintLatex() {
		X a1(N.a(1)), a2(N.a(2)),psi(N.psi), a(N.a),z(N.z);

		TS_ASSERT_EQUALS(to_latex_canonical_string(-1/(a*a)),"-a^{-2}");
		TS_ASSERT_EQUALS(to_latex_canonical_string(-a),"-a");
		TS_ASSERT_EQUALS(to_latex_canonical_string(-a*psi),"-a \\psi");
		TS_ASSERT_EQUALS(to_latex_canonical_string(-z*psi),"-\\psi z");
		TS_ASSERT_EQUALS(to_latex_canonical_string(-a+a1*a1-a2-psi),"-a+a_1^2-a_2-\\psi");
		TS_ASSERT_EQUALS(to_latex_canonical_string(1/a1*a*psi*a/a1),"a^2 a_1^{-2} \\psi");
		TS_ASSERT_EQUALS(to_latex_canonical_string(a1/a),"a^{-1} a_1");
		TS_ASSERT_EQUALS(to_latex_canonical_string((a+a1)/a2),"(a+a_1) a_2^{-1}");
		TS_ASSERT_EQUALS(to_latex_canonical_string(a*a1),"a a_1");
		TS_ASSERT_EQUALS(to_latex_canonical_string(a*a1+a),"a+a a_1");
		TS_ASSERT_EQUALS(to_latex_canonical_string(pow(a,pow(a,a))),"a^{a^a}");
		TS_ASSERT_EQUALS(to_latex_canonical_string(pow(psi,pow(psi,psi))),"\\psi^{\\psi^\\psi}");
		TS_ASSERT_EQUALS(to_latex_canonical_string(1/ex(2)*(1+a1)*z),"\\frac{1}{2} (1+a_1) z");

	}

	void testCanonicalPrintForms() {
		ConcreteManifold M(6);
		symbol a("a");
		TS_ASSERT_EQUALS(to_latex_canonical_string(ParseDifferentialForm(M.e(),"2*12")),"2 e^{12}");
		TS_ASSERT_EQUALS(to_latex_canonical_string(ParseDifferentialForm(M.e(),"2*12+34")),"2 e^{12}+e^{34}");
		TS_ASSERT_EQUALS(to_latex_canonical_string(ParseDifferentialForm(M.e(),"-2*12-34")),"-2 e^{12}-e^{34}");
		TS_ASSERT_EQUALS(to_latex_canonical_string(ParseDifferentialForm(M.e(),"-[a]*12-34",a)),"-a e^{12}-e^{34}");

	}
};

class SimplifierTestSuite: public CxxTest::TestSuite  {
public:
	void test_default_simplifier() {
		realsymbol t("t");
		ex x=1+t, y=V{N.eta};
		TS_ASSERT_EQUALS(default_simplifier.Simplify(x),x);
		//TS_ASSERT_EQUALS(default_simplifier.Simplify(sin(t)*sin(t)+cos(t)*cos(t)),1);
		TS_ASSERT_EQUALS(default_simplifier.Simplify((t*t-1)/(t-1)),t+1);
		TS_ASSERT_EQUALS(default_simplifier.Simplify(x*y),x*y);
		TS_ASSERT_EQUALS(default_simplifier.Simplify(x+y),x+y);

		Lambda1<V> w,u,z;
		symbol a,b,c;
		ex p=w*u*sin(a)+b+z;
		TS_ASSERT_EQUALS(default_simplifier.Simplify(p),p);
	}
};

