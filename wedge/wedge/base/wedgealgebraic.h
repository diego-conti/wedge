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
#ifndef WEDGEALGEBRAIC_H
#define WEDGEALGEBRAIC_H
/** @ingroup Base */
 
/** @{ 
 * @file wedgealgebraic.h 
 * @brief Template classes used in %Wedge to write classes supporting %GiNaC's RTTI system.

Many classes in %Wedge are derived from GiNaC::basic. In this documentation, they are referred to as ''algebraic'' classes, since one is allowed to perform algebraic operations on them. Early versions of %GiNaC implement a Runtime Type Information (RTTI) system, which is essential in the handling of algebraic expressions and requires some extra work in the definition  of each class. The template class Registered contains most of this extra code. GiNaC provides two macros serving the same function (viz. GINAC_DECLARE_REGISTERED_CLASS and GINAC_IMPLEMENT_REGISTERED_CLASS), but Registered has the advantage
of supporting template classes as well. More precisely, one can derive a template class from basic (or symbol, or matrix...) and use the Registered template to make %GiNaC recognize each instance of that template.

The template can also be instructed to implement additional "attributes". There are currently two attributes; the Named attribute implements a name for each object of that class, and the Numbered attribute implements a unique ID for each object, which is used for object comparison.

In order to handle attributes (considering the cases where one derives from a class which is either a Registered itself, or a GiNaC class which implements either attribute), the template instance to derive is given as a member of another template class call Register.

@note Registered does not provide a general method to define classes with archiving support (though one can always add archiving support ''by hand'' if needed).
 
To define an Algebraic class follow these steps:
 - derive your class from Registered, as in class YourClass : public Register<YourClass, basic>::Algebraic
 - use Named, Numbered or NameNumbered instead of Algebraic to apply specific attributes.
 - if the Numbered attribute is not specified, consider declaring and implementing compare_same_type(const GiNaC::basic & other), depending on what the superclass does.
 - define a public static function returning the classname, i.e. static const char* static_class_name() {return "YourClass";}

 @remark The above applies equally well to template classes.

As a shortcut to define a class which contains no data except the ID, %Wedge provides the macro WEDGE_DECLARE_ALGEBRAIC(classname, superclass),
which defines a Numbered algebraic class named classname and derived from superclass. The macro WEDGE_DECLARE_NAMED_ALGEBRAIC is similar, but also implements the Named feature.
  
@warning No two objects with different types can be considered equal by operator==(ex,ex). This means that if A derives from B and B is numbered,
 then an object of type A is equal to an object of type B if they have the same ID, but if the two objects are converted  to ex's, they will be no longer equal.  @see TrivialPairingOperator<Lambda<V> >    
 */

#include "wedge/base/wedgebase.h"
#include "wedge/base/classname.h"
#include "wedge/base/logging.h"
#include "wedge/convenience/parse.h"
#include "wedge/convenience/named.h"

/* @brief Defines a Numbered algebraic class named classname and derived from superclass.
 * 
 * This macro is a shortcut for the procedure described in the documentation of this file.
 */
#define WEDGE_DECLARE_ALGEBRAIC(classname,superclass)									\
class classname : public Register<classname,superclass>::Numbered							\
{															\
public:															\
	static const char* static_class_name() {return #classname;}							\
};

/* @brief Defines a Named, Numbered algebraic class named classname and derived from superclass.
 * 
 * This macro is a shortcut for the procedure described in the documentation of this file.
 */
#define WEDGE_DECLARE_NAMED_ALGEBRAIC(classname,superclass)								\
class classname : public Register<classname,superclass>::NamedNumbered							\
{															\
public:															\
	classname() {}													\
	classname(const Name& name) : Register<classname,superclass>::NamedNumbered(name) {}								\
	static const char* static_class_name() {return #classname;}							\
};
 

namespace Wedge {
using namespace std;
using namespace GiNaC;


//Wedge's attribute systems. Used internally by the registration system.
struct True {	static const bool value=true;	};
struct False {	static const bool value=false;	};
class NamedAttribute {};
class NumberedAttribute {};

//by default, a class's attributes are inherited from the immediate ancestor
template<typename T,typename Attribute> struct Has : public Has<typename T::inherited,Attribute> {};

//for basic, all attributes are false
template<typename Attribute> struct Has<basic,Attribute> : public False {};
//symbol is both Named and Numbered
template<> struct Has<symbol,NamedAttribute> : public True {};
template<> struct Has<symbol,NumberedAttribute> : public True {};

//a class has an attribute if it typedefs the attribute class (see template class Register)
template<typename T> struct Has<T,typename T::Named> : public True {};
template<typename T> struct Has<T,typename T::Numbered> : public True {};

//a class containing the typedefs of the base class to use for registering
template<typename subclass, typename superclass> struct Register;

//TODO what about ginac::function ?


/** @brief Template base class of all Wedge classes, replaces Ginac's RTTI macros

 A RegisteredNamed class, like GiNaC::symbol, contains both a name for default output and a name for \f$\text{\TeX}\f$ output.
 The \f$\text{\TeX}\f$ name, if not provided, is generated by a call to ParseVariableName.
*/
template<typename subclass,typename supername,bool makeNumbered,bool makeNamed> class Registered;

template<bool HasName> struct ForwardName;
template<> struct ForwardName<true> {
	static const Name& cast(const Name& name) {return name;} 
};
template<> struct ForwardName<false> {
	static void cast(const Name& name) {
		throw WedgeException<std::invalid_argument>("Trying to create a nameless object with name " + name.plain(),__FILE__,__LINE__); 
	}
};

template<typename subclass,typename supername> class Registered<subclass,supername,false,false> : public supername {
//constructors
public:
	Registered() {};
	template<typename T> Registered(T x) : supername(x) {}	//required by Lambda
	Registered(const Name& name) : supername(ForwardName<Has<supername, NamedAttribute>::value>::cast(name)) {}
//implements Registered
protected:
	typedef Registered<subclass,supername,false,false> RegClass;
private:
	static GiNaC::registered_class_info reg_info;		
public:
	typedef supername inherited; 
	static GiNaC::registered_class_info &get_class_info_static() { return reg_info; } 		
	virtual const GiNaC::registered_class_info &get_class_info() const { return subclass::get_class_info_static(); }
	virtual GiNaC::registered_class_info &get_class_info() { return subclass::get_class_info_static(); }
	virtual const char *class_name() const { return subclass::get_class_info_static().options.get_name(); }		
	class visitor { 																		
	public: 																				
		virtual void visit(const subclass &) = 0; 											
		virtual ~visitor() {}																
	};																						
	virtual  supername* duplicate() const { return new subclass(static_cast<const subclass&>(*this)); } 					
	virtual void accept(GiNaC::visitor & v) const {											
		if (visitor *p = dynamic_cast<visitor *>(&v)) p->visit(static_cast<const subclass&>(*this));
		else inherited::accept(v); 															
	}
	/** @brief Overloaded comparison operator
	 * 
	 * Unlike GiNaC::compare, this operator may return equality for objects with different types
	 * 
	 */
	bool operator==(const subclass& other) const
	{
		return this->compare_same_type(other)==0;
	}
	static const char* static_class_name() {throw WedgeException<std::logic_error>("Member static_class_name not defined",__FILE__,__LINE__);}
};

template<typename subclass,typename supername> class Registered<subclass,supername,true,false> : public supername {
//constructors
public:
	Registered(const Name& name) : supername(ForwardName<Has<supername, NamedAttribute>::value>::cast(name)) {}
	Registered() {ID=++lastID;}
//attributes
	typedef NumberedAttribute Numbered;
//implements Registered
protected:
	typedef Registered<subclass,supername,true,false> RegClass;
private:
	static GiNaC::registered_class_info reg_info;		
public:
	typedef supername inherited; 
	static GiNaC::registered_class_info &get_class_info_static() { return reg_info; } 		
	virtual const GiNaC::registered_class_info &get_class_info() const { return subclass::get_class_info_static(); }
	virtual GiNaC::registered_class_info &get_class_info() { return subclass::get_class_info_static(); }
	virtual const char *class_name() const { return subclass::get_class_info_static().options.get_name(); }		
	class visitor { 																		
	public: 																				
		virtual void visit(const subclass &) = 0; 											
		virtual ~visitor() {}																
	};																						
	virtual  supername* duplicate() const { return new subclass(static_cast<const subclass&>(*this)); } 					
	virtual void accept(GiNaC::visitor & v) const {											
		if (visitor *p = dynamic_cast<visitor *>(&v)) p->visit(static_cast<const subclass&>(*this));
		else inherited::accept(v); 															
	}
	/** @brief Overloaded comparison operator
	 * 
	 * Unlike GiNaC::compare, this operator may return equality for objects with different types
	 * 
	 */
	bool operator==(const subclass& other) const
	{
		return this->compare_same_type(other)==0;
	}
	static const char* static_class_name() {throw WedgeException<std::logic_error>("Member static_class_name not defined",__FILE__,__LINE__);}
//implements Numbered
protected:
	int ID;				//object ID
	static int lastID;		//last assigned object ID
	unsigned calchash() const {
		this->hashvalue = basic::calchash()+ ID;
		this->setflag(status_flags::hash_calculated);
		return this->hashvalue;
	}
public:
	int compare_same_type(const GiNaC::basic& o) const
	{
		assert(dynamic_cast<const subclass*>(&o)!=NULL);
		const subclass& other=static_cast<const subclass&>(o);
		if (ID<other.ID) return -1;
		else if (ID>other.ID) return 1;
		else return 0;
	}
};

template<typename subclass,typename supername> class Registered <subclass,supername,true,true>: public supername {
//constructors
public:
	Registered() {ID=++lastID;}
	Registered(const Name& n) : name(n.plain()), TeX_name(n.tex()) {ID=++lastID;}
//attributes
	typedef NumberedAttribute Numbered;
	typedef NamedAttribute Named;
//implements Registered
protected:
	typedef Registered<subclass,supername,true,true> RegClass;
private:
	static GiNaC::registered_class_info reg_info;		
public:
	typedef supername inherited; 
	static GiNaC::registered_class_info &get_class_info_static() { return reg_info; } 		
	virtual const GiNaC::registered_class_info &get_class_info() const { return subclass::get_class_info_static(); }
	virtual GiNaC::registered_class_info &get_class_info() { return subclass::get_class_info_static(); }
	virtual const char *class_name() const { return subclass::get_class_info_static().options.get_name(); }		
	class visitor { 																		
	public: 																				
		virtual void visit(const subclass &) = 0; 											
		virtual ~visitor() {}																
	};																						
	virtual  supername* duplicate() const { return new subclass(static_cast<const subclass&>(*this)); } 					
	virtual void accept(GiNaC::visitor & v) const {											
		if (visitor *p = dynamic_cast<visitor *>(&v)) p->visit(static_cast<const subclass&>(*this));
		else inherited::accept(v); 															
	}
	/** @brief Overloaded comparison operator
	 * 
	 * Unlike GiNaC::compare, this operator may return equality for objects with different types
	 * 
	 */
	bool operator==(const subclass& other) const
	{
		return compare_same_type(other)==0;
	}
	static const char* static_class_name() {throw WedgeException<std::logic_error>("Member static_class_name not defined",__FILE__,__LINE__);}
//implements Numbered
protected:
	int ID;				//object ID
	static int lastID;		//last assigned object ID
	unsigned calchash() const {
		this->hashvalue = basic::calchash()+ ID;
		this->setflag(status_flags::hash_calculated);
		return this->hashvalue;
	}
public:
	int compare_same_type(const GiNaC::basic& o) const
	{
		assert(dynamic_cast<const subclass*>(&o)!=NULL);
		const subclass& other=static_cast<const subclass&>(o);
		if (ID<other.ID) return -1;
		else if (ID>other.ID) return 1;
		else return 0;
	}

//implements Named	
	virtual void print(const print_context &c, unsigned level= 0) const {
		if (dynamic_cast<const print_latex*>(&c)!=NULL)
			c.s<<this->get_tex_name();
		else if (dynamic_cast<const print_tree*>(&c)!=NULL) 
			c.s << std::string(level, ' ') << this->get_name() << " (" << this->class_name() << ")" << " @" << this
				<<", ID="<<this->ID<< std::hex << ", hash=0x" << this->hashvalue<< std::dec<<std::endl;
		else c.s<<this->get_name();
	}
protected:
	string name;
	string TeX_name;
public:
	string get_name() const {return name;}
	string get_tex_name() const {return TeX_name;}	
};


template<typename subclass,typename supername> class Registered <subclass,supername,false,true>: public supername {
//constructors
public:
	Registered() {}
	Registered(const Name& n) : name(n.plain()), TeX_name(n.tex()) {}
//attributes
	typedef NamedAttribute Named;
//implements Registered
protected:
	typedef Registered<subclass,supername,false,true> RegClass;
private:
	static GiNaC::registered_class_info reg_info;		
public:
	typedef supername inherited; 
	static GiNaC::registered_class_info &get_class_info_static() { return reg_info; } 		
	virtual const GiNaC::registered_class_info &get_class_info() const { return subclass::get_class_info_static(); }
	virtual GiNaC::registered_class_info &get_class_info() { return subclass::get_class_info_static(); }
	virtual const char *class_name() const { return subclass::get_class_info_static().options.get_name(); }		
	class visitor { 																		
	public: 																				
		virtual void visit(const subclass &) = 0; 											
		virtual ~visitor() {}																
	};																						
	virtual  supername* duplicate() const { return new subclass(static_cast<const subclass&>(*this)); } 					
	virtual void accept(GiNaC::visitor & v) const {											
		if (visitor *p = dynamic_cast<visitor *>(&v)) p->visit(static_cast<const subclass&>(*this));
		else inherited::accept(v); 															
	}
	/** @brief Overloaded comparison operator
	 * 
	 * Unlike GiNaC::compare, this operator may return equality for objects with different types
	 * 
	 */
	bool operator==(const subclass& other) const
	{
		return compare_same_type(other)==0;
	}
	static const char* static_class_name() {throw WedgeException<std::logic_error>("Member static_class_name not defined",__FILE__,__LINE__);}
//implements Named	
	virtual void print(const print_context &c, unsigned level= 0) const {
		if (dynamic_cast<const print_latex*>(&c)!=NULL)
			c.s<<this->get_tex_name();
		else if (dynamic_cast<const print_tree*>(&c)!=NULL) 
			c.s << std::string(level, ' ') << this->get_name() << " (" << this->class_name() << ")" << " @" << this
				<< std::hex << ", hash=0x" << this->hashvalue<< std::dec<<std::endl;
		else c.s<<this->get_name();
	}
protected:
	string name;
	string TeX_name;
public:
	string get_name() const {return name;}
	string get_tex_name() const {return TeX_name;}	
};


//classes whose immediate superclass is symbol, possymbol or realsymbol must use a different template class to harmonize symbol's interface with
//the interface of Registered<.,.,true,true>
template<typename subclass,typename supername> class RegisteredSymbol : public supername {
//constructors
public:
	RegisteredSymbol() {}
	RegisteredSymbol(const Name& n) : supername(n.plain(), n.tex()) {}
//attributes
	typedef NumberedAttribute Numbered;
	typedef NamedAttribute Named;
//implements Registered
protected:
	typedef RegisteredSymbol<subclass,supername> RegClass;
private:
	static GiNaC::registered_class_info reg_info;		
public:
	typedef supername inherited; 
	static GiNaC::registered_class_info &get_class_info_static() { return reg_info; } 		
	virtual const GiNaC::registered_class_info &get_class_info() const { return subclass::get_class_info_static(); }
	virtual GiNaC::registered_class_info &get_class_info() { return subclass::get_class_info_static(); }
	virtual const char *class_name() const { return subclass::get_class_info_static().options.get_name(); }		
	class visitor { 																		
	public: 																				
		virtual void visit(const subclass &) = 0; 											
		virtual ~visitor() {}																
	};																						
	virtual  supername* duplicate() const { return new subclass(static_cast<const subclass&>(*this)); } 					
	virtual void accept(GiNaC::visitor & v) const {											
		if (visitor *p = dynamic_cast<visitor *>(&v)) p->visit(static_cast<const subclass&>(*this));
		else inherited::accept(v); 															
	}
	/** @brief Overloaded comparison operator
	 * 
	 * Unlike GiNaC::compare, this operator may return equality for objects with different types
	 * 
	 */
	bool operator==(const subclass& other) const
	{
		return this->compare_same_type(other)==0;
	}
	static const char* static_class_name() {throw WedgeException<std::logic_error>("Member static_class_name not defined",__FILE__,__LINE__);}
//implements Named	
	string get_tex_name() const {return this->TeX_name;}
};




template<typename subclass, typename superclass> struct Register {
	typedef Registered<subclass, superclass,false,false> Algebraic;
	typedef Registered<subclass, superclass,false,!Has<superclass,NamedAttribute>::value> Named;
	typedef Registered<subclass, superclass,!Has<superclass,NumberedAttribute>::value,false> Numbered;
	typedef Registered<subclass, superclass,!Has<superclass,NumberedAttribute>::value,!Has<superclass,NamedAttribute>::value> NamedNumbered;
};

template<typename subclass> struct Register<subclass,symbol> {
	typedef RegisteredSymbol<subclass, symbol> Algebraic;
	typedef RegisteredSymbol<subclass, symbol> Named;
	typedef RegisteredSymbol<subclass, symbol> Numbered;
	typedef RegisteredSymbol<subclass, symbol> NamedNumbered;
};
template<typename subclass> struct Register<subclass,realsymbol> {
	typedef RegisteredSymbol<subclass, realsymbol> Algebraic;
	typedef RegisteredSymbol<subclass, realsymbol> Named;
	typedef RegisteredSymbol<subclass, realsymbol> Numbered;
	typedef RegisteredSymbol<subclass, realsymbol> NamedNumbered;
};
template<typename subclass> struct Register<subclass,possymbol> {
	typedef RegisteredSymbol<subclass, possymbol> Algebraic;
	typedef RegisteredSymbol<subclass, possymbol> Named;
	typedef RegisteredSymbol<subclass, possymbol> Numbered;
	typedef RegisteredSymbol<subclass, possymbol> NamedNumbered;
};




namespace internal {

//used internally to allow for the fact that realsymbol::visitor accepts a symbol rather than a realsymbol;
template<typename T> struct VisitorClass
{
	typedef T type;
};

template<> struct VisitorClass<realsymbol> : VisitorClass<symbol> {};
template<> struct VisitorClass<possymbol> : VisitorClass<symbol> {};

}

template<typename subclass,typename supername> int Registered<subclass,supername,true,false>::lastID;
template<typename subclass,typename supername> int Registered<subclass,supername,true,true>::lastID;

// The following serves the same purpose as GINAC_IMPLEMENT_REGISTERED_CLASS(subclass,superclass);
template<typename subclass,typename superclass> 
  GiNaC::registered_class_info Registered<subclass,superclass,true,true>::reg_info =
	GiNaC::registered_class_info(GiNaC::registered_class_options(
		subclass::static_class_name(),
		className<superclass>,
		typeid(subclass)
	));
template<typename subclass,typename superclass> 
  GiNaC::registered_class_info Registered<subclass,superclass,true,false>::reg_info =
	GiNaC::registered_class_info(GiNaC::registered_class_options(
		subclass::static_class_name(),
		className<superclass>,
		typeid(subclass)
	));
template<typename subclass,typename superclass> 
  GiNaC::registered_class_info Registered<subclass,superclass,false,true>::reg_info =
	GiNaC::registered_class_info(GiNaC::registered_class_options(
		subclass::static_class_name(),
		className<superclass>,
		typeid(subclass)
	));
template<typename subclass,typename superclass> 
  GiNaC::registered_class_info Registered<subclass,superclass,false,false>::reg_info =
	GiNaC::registered_class_info(GiNaC::registered_class_options(
		subclass::static_class_name(),
		className<superclass>,
		typeid(subclass)
	));
template<typename subclass,typename superclass> 
  GiNaC::registered_class_info RegisteredSymbol<subclass,superclass>::reg_info =
	GiNaC::registered_class_info(GiNaC::registered_class_options(
		subclass::static_class_name(),
		className<superclass>,
		typeid(subclass)
	));





namespace internal {

/** Rotate bits of unsigned value by one bit to the left (from ginac/utils.h). */
inline unsigned rotate_left(unsigned n)
{
	return (n & 0x80000000U) ? (n << 1 | 0x00000001U) : (n << 1);
}

}





} /** @} */
#endif
