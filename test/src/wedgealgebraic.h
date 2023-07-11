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
#include <cxxtest/TestSuite.h>
#include "test.h"

using namespace Wedge;
using namespace GiNaC;

WEDGE_DECLARE_NAMED_ALGEBRAIC(A,symbol);
WEDGE_DECLARE_ALGEBRAIC(B,A)
WEDGE_DECLARE_NAMED_ALGEBRAIC(C,realsymbol);
WEDGE_DECLARE_ALGEBRAIC(D,basic)

class WedgeAlgebraicInheritanceTestSuite : public CxxTest::TestSuite 
{
public:
	void testInheritanceFromBasic() {
		ex d=D{};
		TS_ASSERT(is_a<D>(d));
		TS_ASSERT(is_a<basic>(d));
		TS_ASSERT(typeid(D::inherited)==typeid(basic));		
		TS_ASSERT_EQUALS(D{}.get_class_info().options.get_name(),"D");
		TS_ASSERT_EQUALS(D{}.get_class_info().options.get_parent_name(),"basic");
	}

	void testRegisterRealsymbol() {
		C ac(N.c);
		auto& classinfo=ac.get_class_info();
		TS_ASSERT_EQUALS(classinfo.options.get_name(),"C");
		TS_ASSERT_EQUALS(classinfo.options.get_parent_name(),"symbol");
		auto parent=classinfo.get_parent();
		TS_ASSERT_EQUALS(parent->options.get_name(),"symbol");
	}

	void testInheritanceFromRealSymbol() {
		ex c=C{};
		TS_ASSERT_EQUALS(c.real_part(),c);
		TS_ASSERT_EQUALS(c.imag_part(),0);
		TS_ASSERT_DIFFERS(c.real_part(),0);
		TS_ASSERT(is_a<C>(c));
		TS_ASSERT(is_a<basic>(c));
		TS_ASSERT(is_a<realsymbol>(c));
		TS_ASSERT(typeid(C::inherited)==typeid(realsymbol));		
	}

		void testNamedAlgebraic() {
			A a1(N.a(1)),a2(N.a(2));			
			TS_ASSERT_DIFFERS(a1,a2);
			TS_ASSERT_DIFFERS(a1,2*a1);
			ex a1ex=a1, a2ex=a2;
			TS_ASSERT_EQUALS(a1ex,a1);
			TS_ASSERT_EQUALS(a1,a1ex);
			TS_ASSERT((a1-a1ex).is_zero());
			TS_ASSERT_DIFFERS(a1ex,a2ex);
			ex a1ex_symbol=static_cast<symbol>(a1),a2ex_symbol=static_cast<symbol>(a2); 
			TS_ASSERT(BASIC_DERIVED_DIFFERS(a1ex_symbol,a2ex_symbol));
			
			//the ones below fail. This is because GiNaC considers objects with different types different.
			//TS_ASSERT(BASIC_DERIVED_EQUALS(a1ex_symbol,a1));
			//TS_ASSERT(BASIC_DERIVED_EQUALS(a1,a1ex_symbol));
			//TS_ASSERT(BASIC_DERIVED_EQUALS(a1ex,a1ex_symbol));
		}
		
		void testInheritance() {
			ex a=A{};
			ex b=B{};
			TS_ASSERT(is_a<B>(b));
			TS_ASSERT(is_a<A>(b));
			TS_ASSERT(is_a<A>(a));
			TS_ASSERT(!is_a<B>(a));		
			ex basa=ex_to<A>(b);
			TS_ASSERT(is_a<A>(basa));
			TS_ASSERT_EQUALS(basa,b);			
		}
};

class MyAlgebraicClass : public Register<MyAlgebraicClass,basic>::Algebraic
{
	public:
	static const char* static_class_name() {return "MyAlgebraicClass";}
	virtual int compare_same_type(const GiNaC::basic & other) const {return basic::compare_same_type(other);}
};

class MyNumberedClass : public Register<MyNumberedClass,basic>::Numbered
{
public:
	static const char* static_class_name() {return "MyNumberedClass";}
};

class MyInheritedClass : public Register<MyInheritedClass,MyNumberedClass>::Algebraic
{
public:
	static const char* static_class_name() {return "MyInheritedClass";}
	int compare_same_type(const GiNaC::basic & other) const {
		return MyNumberedClass::compare_same_type(other);
	}
};

class MyNamedClass : public Register<MyNamedClass,basic>::Named
{
	public:
	MyNamedClass(const Name& name) : Register<MyNamedClass,basic>::Named(name) {}
	virtual int compare_same_type(const GiNaC::basic & other) const {return basic::compare_same_type(other);}
	static const char* static_class_name() {return "MyNamedClass";}
};

class MyNamedNumberedClass :
	public Register<MyNamedNumberedClass,basic>::NamedNumbered
{
public:
	MyNamedNumberedClass(const Name& name)  : Register<MyNamedNumberedClass,basic>::NamedNumbered (name) {}
	static const char* static_class_name() {return "MyNamedNumberedClass";}
};

WEDGE_DECLARE_NAMED_ALGEBRAIC(X,symbol)

WEDGE_DECLARE_NAMED_ALGEBRAIC(Y,realsymbol)

class AlgebraicTestSuite : public CxxTest::TestSuite  {
public:
	void testLst() {	//check for a problem that could arise when linking with ginac and wedge in the wrong order
		lst p;
		TS_ASSERT(p.info(info_flags::list));
		TS_ASSERT(ex(p).info(info_flags::list));
		p.append(1);
		TS_ASSERT(p.info(info_flags::list));
		TS_ASSERT(ex(p).info(info_flags::list));
		p.append(2);
		TS_ASSERT(p.info(info_flags::list));
		TS_ASSERT(ex(p).info(info_flags::list));
	}

	void testSymbol() {
		X x;
		TS_ASSERT(x.info(info_flags::symbol));
		TS_ASSERT(is_a<symbol>(x));
		Y y;
		TS_ASSERT(y.info(info_flags::symbol));
		TS_ASSERT(is_a<symbol>(y));
	}	

	void testAlgebraic() {
		ex myAlgebraicObject=MyAlgebraicClass();
		
		ex myNumberedObject=MyNumberedClass();
		ex myNumberedObject2=MyNumberedClass();
		ex e=myNumberedObject+myNumberedObject2;
		TS_ASSERT_EQUALS(myNumberedObject+myNumberedObject2,e);
		TS_ASSERT(myNumberedObject!=myNumberedObject2);
		{
			string name("A name");
			MyNamedClass myNamedObject(name);
			stringstream ss;
			ss<<myNamedObject;		
			TS_ASSERT_EQUALS(ss.str(),name);
		}
		{
			string name("Another name");
			MyNamedNumberedClass myNamedNumberedObject(name);
			stringstream ss;
			ss<<myNamedNumberedObject;		
			TS_ASSERT_EQUALS(ss.str(),name);
		}
	}
	
	void testDiff() {
		Y y;
		X x;
		TS_ASSERT_EQUALS(diff(y,x),0);	//this gives a wrong assertion if Y defines it's own compare_same_type rather than symbol::compare_same_type
	}
	
	
	class MyVisitor : public MyNumberedClass::visitor, public visitor
	{
	public:
		MyVisitor() {flag=false;}
		bool flag;
		void visit(const MyNumberedClass& e)
		{
			flag=true;	
		} 
	};
	
	void testDoubleInheritance()
	{
		ex e=MyInheritedClass();		
		TS_ASSERT(is_a<MyInheritedClass>(e))
		const MyInheritedClass& inherited=ex_to<MyInheritedClass>(e);		
		e=inherited;
		MyVisitor v;
		e.accept(v);
		TS_ASSERT(v.flag);
		
		TS_ASSERT(is_a<MyNumberedClass>(e));		
		const MyNumberedClass& numbered=ex_to<MyNumberedClass>(e);		
		e=numbered;
		e.accept(v);
		TS_ASSERT(v.flag);
		
		TS_ASSERT(is_a<MyNumberedClass>(e));
		const MyInheritedClass& in2=dynamic_cast<const MyInheritedClass&>(numbered);		
	}
	
};


