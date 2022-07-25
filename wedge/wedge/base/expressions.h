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

#ifndef EXPRESSIONS_H
#define EXPRESSIONS_H

#include "logging.h"
#include "wedgebase.h"
#include "wedgealgebraic.h"
#include <ginac/power.h>
/** @ingroup Base 
 * */

 /** @{ 
 * @file expressions.h
 * @brief Expressions as linear combinations of Vector objects. 
 */
 
namespace Wedge {
//TODO find out if this is really needed
class VectorBase {};

/** @brief The class representing a simple element of a generic vector space.
 *
 * To define a vector space, subclass from Vector.
 *
 */

class Vector : public Register<Vector,basic>::Algebraic, public VectorBase
{
public:
	enum {IsTypeVector};
	static const char* static_class_name() {return "Vector";}							\
};


/** @brief Common superclass for all template classes of the form Lambda<V>, where V derives from Vector
 *
 * @relates Lambda
 */

class LambdaVector : public Register<LambdaVector,ncmul>::Algebraic, public VectorBase
{
public:
	using RegClass::RegClass;
	enum {IsTypeVector};
	static const char* static_class_name() {return "LambdaVector";}							\
};


enum GetCoefficientsFlags {withRHS=0,withoutRHS=1};

/** @brief Return the coefficients appearing in a linear expression
 * @param e A linear combination of objects of type T
 * @param flags Flag specifying whether equations should be returned or expressions 
 * @param container (out) A container where the coefficients are to be stored
 * @returns A reference to container
 * 
 * - If T is derived from either Vector or LambdaVector, the simpler overloaded template function GetCoefficients<typename Container>
 * is equivalent to this function.
 * - if T is Poly<V>, the parameter e is allowed to be a polynomial in objects of type v
 * - if T is Poly1<V>, the parameter e is allowed to be a degree one polynomial in objects of type v
 * - If T is Lambda<V>, the parameter e is allowed to be a degree one polynomial combination of objects of type V, Lambda1<V> and Lambda<V>. 
 * 
 * @note In %Wedge, functions whose name starts with Get adhere to this convention: 
 * - The first parameter is a reference to a container
 * - The container is updated by inserting new elements
 * - A reference to the container is returned
  */ 
template<typename T, typename Container> Container& GetCoefficients(Container&, ex , GetCoefficientsFlags flags=withoutRHS);
template<typename T, typename Container> Container& GetCoefficients(Container&, const list<ex>&, GetCoefficientsFlags flags=withoutRHS);
template<typename T, typename Container> Container& GetCoefficients(Container&, const exvector&, GetCoefficientsFlags flags=withoutRHS);
template<typename T, typename Container> Container& GetCoefficientsComplex(Container&, ex , GetCoefficientsFlags flags=withoutRHS);
template<typename T, typename Container> Container& GetCoefficientsComplex(Container&, const list<ex>&, GetCoefficientsFlags flags=withoutRHS);
template<typename T, typename Container> Container& GetCoefficientsComplex(Container&, const exvector&, GetCoefficientsFlags flags=withoutRHS);

template<typename T,typename Container, typename InputIterator> Container& GetSimple(Container& c, InputIterator from, InputIterator to);
template<typename T,typename Container> Container& GetSimple(Container& c, const exvector&);
template<typename T,typename Container> Container& GetSimple(Container& c, ex);
template<typename T, typename Container> Container& GetSymbols(Container& c,ex e);
template<typename T, typename Container,typename InputIterator> Container& GetSymbols(Container& c,InputIterator from, InputIterator to);


/* @brief Return the coefficients appearing in a linear expression
 * @param e A linear combination or degree one polynomial in objects whose type is derived from Vector or LambdaVector
 * @param flags Flag specifying whether equations should be returned or expressions 
 * @param container (out) A container where the coefficients are to be stored
 * @returns A reference to container
 * 
 * This function can usually be used as a shorthand for GetCoefficients<T>. Explicit specification of the template
 * parameter is still necessary when T does not derive from Vector, e.g. T=VectorSpace::Coordinate
 */
template<typename Container> Container& GetCoefficients(Container& container, ex e, GetCoefficientsFlags flags=withoutRHS)
	{return GetCoefficients<VectorBase>(container,e,flags);}
template<typename Container> Container& GetCoefficients(Container& container, const list<ex>& e, GetCoefficientsFlags flags=withoutRHS)
	{return GetCoefficients<VectorBase>(container,e,flags);}
template<typename Container> Container& GetCoefficients(Container& container,  const exvector& e, GetCoefficientsFlags flags=withoutRHS)
	{return GetCoefficients<VectorBase>(container,e,flags);}
template<typename Container> Container& GetCoefficientsComplex(Container& container, ex e, GetCoefficientsFlags flags=withoutRHS)
	{return GetCoefficientsComplex<VectorBase>(container,e,flags);}
template<typename Container> Container& GetCoefficientsComplex(Container& container, const list<ex>& e, GetCoefficientsFlags flags=withoutRHS)
	{return GetCoefficientsComplex<VectorBase>(container,e,flags);}
template<typename Container> Container& GetCoefficientsComplex(Container& container,  const exvector& e, GetCoefficientsFlags flags=withoutRHS)
	{return GetCoefficientsComplex<VectorBase>(container,e,flags);}

/** @brief Dummy class describing the type of a polynomial in variables of type T 
 *
 * To be used as a template argument of GetCoefficients
 */
template<typename T> class Poly {};
/** @brief Dummy class describing the type of a degree one polynomial in variables of type T
 *
 * To be used as a template argument of GetCoefficients
 */
template<typename T> class Poly1 {};


namespace internal {
//needed by all template visitors that implement void visit(const T&), for compatibility with realsymbol and possymbol
template<typename T> struct VisitorType {
	typedef T Type;
};

template<> struct VisitorType<realsymbol> {
	typedef symbol Type;
};

template<> struct VisitorType<possymbol> {
	typedef symbol Type;
};


//needed in NormalFormHelper<VectorBase> and NormalFormHelper<Lambda<V> >
template<typename T> struct LogicalType {
	
	static inline bool Is(ex x) {return is_a<T>(x);}
};
template<typename T> struct LogicalType<Poly<T> > {
	
	static inline bool Is(ex x) {
		if (is_a<power>(x)) return LogicalType<T>::Is(x.op(0));
		else return is_a<T>(x);
	}
};
template<> struct LogicalType<Vector> {
	static inline bool Is(ex x) {return is_a<VectorBase>(x);}
};
//assumes expand has been called
//expects linear combinations (or polynomials of degree one)
template<typename T> struct FindSimpleHelper : public visitor,public T::visitor, public add::visitor, public mul::visitor, public basic::visitor
{
	set<ex, ex_is_less> occurrences;
	void visit(const T& x) {
		occurrences.insert(x);
	}
	void visit(const add& x)
	{
		for (int i=0;i<x.nops();i++)
			x.op(i).accept(*this);	
	}
	void visit(const mul& x)
	{
		int i=0;
		bool has_a_T=false;
		while (i<x.nops())
		{
			ex opi=x.op(i);
			if (LogicalType<T>::Is(opi)) {
				if (has_a_T) {
					LOG_ERROR(x);
					throw WedgeException<std::runtime_error>("Linear combination of elements of type T expected",__FILE__,__LINE__);
				}
				else {
					occurrences.insert(opi);
					has_a_T=true;
				}
			}			
			++i;
		}
		//has_a_T==false is not an error, because we allow degree-one polynomials. 
	}
	void visit (const basic& x)
	{
		//Not an error, because we allow degree-one polynomials. 
	}
};

//expects a linear combination
//assumes expand has been called
template<typename T> struct NormalFormHelper : public visitor,public T::visitor, public add::visitor, public mul::visitor, public basic::visitor
{
	exmap coeffs;
	void visit(const typename VisitorType<T>::Type& x) {
		++coeffs[x];
	}
	void visit(const mul& x)
	{
		int i=0;
		exvector ops(x.nops()-1);
		exvector::iterator k=ops.begin();
		ex theT;
		while (i<x.nops())
		{
			ex opi=x.op(i);
			if (LogicalType<T>::Is(opi)) {
				if (!theT.is_zero()) {
					LOG_ERROR(x);
					LOG_ERROR(ops);
					throw WedgeException<std::runtime_error>("Linear combination of elements of type T expected",__FILE__,__LINE__);
				}
				theT=opi;
			}
			else if (k==ops.end())	{
				LOG_ERROR(x);
				throw WedgeException<std::runtime_error>("Expression contains a scalar component, but linear combination of vectors was expected",__FILE__,__LINE__);
			}
			else *k++=opi;
			++i;
		}
		assert(k==ops.end());		
		coeffs[theT]+=mul(ops);
	}
	
	void visit(const add& x)
	{
		for (int i=0;i<x.nops();++i)
			x.op(i).accept(*this);	
	}
	void visit (const basic& x)
	{
		if (x!=0) {
			LOG_ERROR(x);
			throw WedgeException<std::runtime_error>("Expression contains a scalar component, but linear combination of vectors was expected",__FILE__,__LINE__);
		}
	}
};

//like NormalFormHelper<T>, but also allows polynomials in T
//assumes expand has been called
template<typename T> struct NormalFormHelper<Poly<T> > : public visitor,public T::visitor, public power::visitor, public add::visitor, public mul::visitor, public basic::visitor
{
	exmap coeffs;
	void visit(const typename VisitorType<T>::Type& x) {
		++coeffs[x];
	}
	void visit(const power& x) {
		if (LogicalType<Poly<T> >::Is(x)) ++coeffs[x];
		else coeffs[1]+=x;
	}
	void visit(const mul& x)
	{

		ex coeffpart=1;
		ex Tpart=1;
		for (int i=0; i<x.nops();++i)
		{
			ex opi=x.op(i);
			if (LogicalType<Poly<T> >::Is(opi)) Tpart*=opi;
			else coeffpart*=opi;
		}
		coeffs[Tpart]+=coeffpart;
	}
	
	void visit(const add& x)
	{
		for (int i=0;i<x.nops();++i)
			x.op(i).accept(*this);	
	}
	void visit (const basic& x)
	{
		coeffs[1]+=x;
	}
};

//like NormalFormHelper<T>, but also allows degree one polynomials in T
//assumes expand has been called
template<typename T> struct NormalFormHelper<Poly1<T> > : public visitor,public T::visitor, public add::visitor, public mul::visitor, public basic::visitor
{
	exmap coeffs;
	void visit(const typename VisitorType<T>::Type& x) {
		++coeffs[x];
	}
	void visit(const mul& x)
	{
		int i=0;
		exvector ops(x.nops()-1);
		exvector::iterator k=ops.begin();
		ex theT;
		while (i<x.nops())
		{
			ex opi=x.op(i);
			if (LogicalType<T>::Is(opi)) {
				if (!theT.is_zero()) {
					LOG_ERROR(x);
					LOG_ERROR(ops);
					throw WedgeException<std::runtime_error>("A degree one polynomial in elements of type T was expected",__FILE__,__LINE__);
				}
				theT=opi;
			}
			else if (k==ops.end())	{
				ops.push_back(opi);
				coeffs[1]+=mul(ops);
				return;
			}
			else *k++=opi;
			++i;
		}
		assert(k==ops.end());		
		coeffs[theT]+=mul(ops);
	}
	
	void visit(const add& x)
	{
		for (int i=0;i<x.nops();++i)
			x.op(i).accept(*this);	
	}
	void visit (const basic& x)
	{
		if (x!=0) coeffs[1]+=x;
	}
};

//Like NormalFormHelper, except it works with linear combinations of either Vector or LambdaVector-derived classes
//Won't work with VectorSpace::Coordinate.
//Notice that the vectors are converted to type Vector/LambdaVector, hence NormalForm<VectorBase> does not give an equivalent expression.
template<> struct NormalFormHelper<VectorBase>  : public NormalFormHelper<Poly1<Vector> >, public LambdaVector::visitor
{
	void visit(const LambdaVector& x) {
		++coeffs[x];
	}
};

template<typename T> struct SetOf
{
	typedef set<T,basic_is_less> type;
};

template<> struct SetOf<ex>
{
	typedef set<ex,ex_is_less> type;
};

template<typename T,typename StoreAs> struct FindSymbolsHelper : public visitor,public T::visitor, public basic::visitor
{
	typename SetOf<StoreAs>::type occurrences;
	void visit(const typename VisitorType<T>::Type& x) {
		occurrences.insert(x);
	}
	void visit(const basic& x)
	{
		for (int i=0;i<x.nops();i++)
			x.op(i).accept(*this);	
	}
};
}


/** @brief Cast an expression into "normal form"
 * @param e A polynomial in objects of type T
 * @returns An expression equivalent to e, with the form \f$\sum_i a_i v_i\f$ where the \f$v_i\f$ are monomials in T.
 * @exception WedgeException<std::runtime_error> Thrown if expression does not have the required type   
 */
template<typename T> ex NormalForm(ex e)
{
	e=e.expand();
	ex result;
	internal::NormalFormHelper<Poly<T> > v;
	e.accept(v);

	exmap inverse; 	//collect similar coefficients
	for (exmap::const_iterator i=v.coeffs.begin();i!=v.coeffs.end();++i)
	{
		ex normal=i->second.normal(); 
		if (inverse.find(-normal)!=inverse.end()) 
			inverse[-normal]-=i->first;
		else
			inverse[normal]+=i->first;
	}

	for (exmap::const_iterator i=inverse.begin();i!=inverse.end();++i)
		result+=i->first*i->second;
	return result;
}

/** @overload
 */
template<typename T> exvector NormalForm(exvector e)
{
	for (exvector::iterator i=e.begin();i!=e.end();i++)
		*i=NormalForm<T>(*i);
	return e;
}


/** @brief Produce a list of the simple components appearing in a linear combination of object of type T
 * @param c (out) A container where the components are to be stored
 * @param e A degree-one polynomial in objects of type T
 * @return e A reference to c
 * @exception WedgeException<std::runtime_error> Thrown if expression does not have the required type
 *  
 * - Degree-one polynomial means a linear combination plus, possibly, an extra `scalar' term, not involving the type T.
 * - If T is Lambda<V>, the parameter e should be a linear combination of objects of type V, Lambda1<V> or Lambda<V>, plus a scalar term.
 *  
 * @sa GetCoefficients
 */ 
template<typename T, typename Container> Container& GetSimple(Container& c, ex e)
{
	e=e.expand();
	return GetSimple<T>(c,&e, &e+1);	
}

/** @overload  
 */
template<typename T,typename Container, typename InputIterator> Container& GetSimple(Container& c, InputIterator from, InputIterator to)
{
	internal::FindSimpleHelper<T> v;
	while (from!=to)
	{
		ex e=from++->expand();
		assert (!is_a<relational>(e));		
		e.accept(v);
	}
	return Insert(c,v.occurrences.begin(),v.occurrences.end());	
}

/** @overload  
 */
 template<typename T, typename Container> Container& GetSimple(Container& c, const exvector& e)
{
	return GetSimple<T>(c,e.begin(),e.end());
}

/** @brief Produce a list of the elements of type T appearing in an expression
 * @param c (out) A container where the symbols are to be stored
 * @param e An expression (not necessarily linear)
 * @return A reference to c
 * 
 * - The value-type of the container may be either T or ex
 * - If symbols of type derived from T appear, using a container with value-type equal to T has the undesirable effect of converting objects to T
 * - If T is Lambda<V>, the value-type of the container must be ex; the matching types are V, Lambda1<V> and Lambda<V>.
 * @sa GetCoefficients
 */ 
template<typename T, typename Container> Container& GetSymbols(Container& c,ex e)
{
	return GetSymbols<T>(c,&e,&e+1);
}

/** @overload  
 */
template<typename T, typename Container,typename InputIterator> Container& GetSymbols(Container& c,InputIterator from, InputIterator to)
{
	internal::FindSymbolsHelper<T,typename Container::value_type> v;
	while (from!=to)
	{
		ex e=from++->expand();		
		e.accept(v);
	}
	return Insert(c,v.occurrences.begin(),v.occurrences.end());
}


template<typename T, typename Container> Container& GetCoefficients(Container& container, ex e,GetCoefficientsFlags noRHS)
{
	e=e.expand();	
	internal::NormalFormHelper<T> v;
	e.accept(v);
	
	for (exmap::const_iterator i=v.coeffs.begin();i!=v.coeffs.end();i++)
	{		
		if (!i->second.is_zero())
			Insert(container,noRHS ? i->second :ex(i->second==0));
	}
	return container;
}

/** @overload 
 */ 
template<typename T, typename Container> Container& GetCoefficients(Container& container,const list<ex>& l, GetCoefficientsFlags noRHS)
{	
	for (list<GiNaC::ex>::const_iterator i=l.begin();i!=l.end();i++)	
		GetCoefficients<T>(container,*i,noRHS);
	return container;
}

/** @overload 
 */ 
template<typename T, typename Container> Container& GetCoefficients(Container& container,const exvector& l, GetCoefficientsFlags noRHS)
{
	for (vector<GiNaC::ex>::const_iterator i=l.begin();i!=l.end();i++)
		GetCoefficients<T>(container,*i,noRHS);
	return container;
}


/** @brief Return the coefficients appearing in a linear expression
 * @param e A linear combination of objects of type T with complex coefficients
 * @param noRHS Flag specifying whether equations should be returned or expressions 
 * @param container (out) A container where the coefficients are to be stored
 * @returns A reference to container
 *
 * This function sets to zero the real and imaginary part separately.
 */
template<typename T, typename Container> Container& GetCoefficientsComplex(Container& container,ex e,GetCoefficientsFlags noRHS)
{
	ex real,imag;
	SplitIntoRealAndImaginary(e.expand(),real,imag);
	GetCoefficients<T>(container,real,noRHS);
	return GetCoefficients<T>(container,imag,noRHS);
}


/** @overload 
 */ 
template<typename T, typename Container> Container& GetCoefficientsComplex(Container& container,const list<ex>& l, GetCoefficientsFlags noRHS)
{
	for (list<GiNaC::ex>::const_iterator i=l.begin();i!=l.end();i++)	
		GetCoefficientsComplex<T>(container,*i,noRHS);
	return container;
}

/** @overload 
 */ 
 template<typename T, typename Container> Container& GetCoefficientsComplex(Container& container,const exvector& l, GetCoefficientsFlags noRHS)
{
	for (vector<GiNaC::ex>::const_iterator i=l.begin();i!=l.end();i++)	
		GetCoefficientsComplex<T>(container,*i,noRHS);
	return container;
}

} /** @} */
#endif /*EXPRESSIONS_H*/
