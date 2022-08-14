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
  #ifndef LAMBDA_H_
#define LAMBDA_H_
/** @ingroup LinearAlgebra */
/** @{ 
 * @file lambda.h
 * @brief Exterior algebra over a vector space
 * 
 * @warning The functionality provided in this file requires patching GiNaC. The simplest possible patch consists
 * in replacing every instance of is_exactly_a<ncmul> in ncmul.cpp with is_a_<ncmul>. 
 */
 
#include "wedge/base/wedgealgebraic.h"
#include "wedge/base/expressions.h"
#include "wedge/linearalgebra/bilinear.h"

namespace Wedge {
using namespace GiNaC;

/** @brief Degree of an element in the exterior algebra over T, or a polynomial in T
 * @exception InhomogeneousExpression thrown if element is not homogeneous
 *
 * If T has the form Poly<V>, the the expression is allowed to be nonhomogeneous, and degree returns the degree of the highest monomial, i.e. the polynomial's degree in the usual sense
 */
template<typename T> inline int Degree(ex);
/** @brief Parity of an element in the exterior algebra over T, or a polynomial in T
 * @exception InhomogeneousExpression thrown if element is not homogeneous
 */
template<typename T> inline bool IsOdd(ex);

// ****************************************************************************
// *							Lambda1						  *
// ****************************************************************************
/** @brief The type of generators of the exterior algebra over V
 *  @param V The type of simple elements in the vector space V
 * 
 * This class is a wrapper around V which allows to produce elements of the exterior algebra
 * using ordinary moltiplication.
 * 
 * Example:
 *   
 * Lambda1<V> a,b;
 * a*b-b*a; //simplifies to zero
 */

template<typename V> class Lambda1 : public Register<Lambda1<V>, V>::Algebraic {
public:
	Lambda1() {}
	Lambda1(const Name& name) : Lambda1::RegClass(name) {}
	/* Functions used internally by the GiNaC framework: */
	static const char* static_class_name() {return "Lambda1";}
	unsigned precedence() const {return 40;}
	int degree() {return 1;}
	int compare_same_type(const GiNaC::basic & other) const {
		return V::compare_same_type(other);
	}
protected:
	ex eval_ncmul(const exvector & v) const; ///< Implements skew-commutativity
	unsigned return_type() const {return return_types::noncommutative;}
};

// ****************************************************************************
// *							Lambda						  *
// ****************************************************************************
/** @brief Decomposable, or simple, elements of the exterior algebra over V
 *  @param V The type of simple elements in the vector space
 * 
 * An instance of Lambda<V> represents a decomposable element in the exterior algebra over V, i.e.
 * a product \f$ v_1\wedge\dots\wedge v_k\f$ of simple elements \f$v_i\f$
 * 
 * @internal @note Every template that depends on Lambda<V> must be specialized so that it accepts both objects of type V
 * and objects of type Lambda<V>, since an object of type V is also in the exterior algebra Lambda<V> 
  */

template<typename V> class Lambda : public Register<Lambda<V>, LambdaVector>::Algebraic {
public:
	Lambda() {}	///< Trivial constructor. Do not use directly.
	
/** @brief Construct a (decomposable) element as a product 
 *  @param v A sequence of elements of the exterior algebra
 * 
 *  Do not call this constructor directly; to create an element of the exterior algebra, use operator*.
 *  
 *  The constructor will then be invoked by the Wedge framework, with parameters of type Lambda1<V>,
 *  thus giving a decomposable element by definition.
 *  
 */
	Lambda(const exvector& v) : Lambda::RegClass (v) {}

	int degree() const {return this->nops();} ///< Degree of this element of the exterior algebra
	bool IsOdd() const {return degree%2!=0;}	///< Parity of this element of the exterior algebra

// Implementation stuff needed by GiNaC	
	static const char* static_class_name() {return "Lambda";} 
	int compare_same_type(const GiNaC::basic & other) const {
		return ncmul::compare_same_type(other); 
	}

private:
	static void  PrintVariableAndIndices(string varname, const vector<string>& indices,const print_context &c)
	{		
		assert(!varname.empty());
		c.s<<varname<<"^";
		int length=0;
		bool comma_needed=false;
		vector<string>::const_iterator i=indices.begin();
		while (i!=indices.end())
		{
			length+=i->size();
			if (i++->size()>1) comma_needed=true;
		}
		if (length>1) c.s<<"{";
		i=indices.begin();
		c.s<<*i++;
		while (i!=indices.end())
		{
			if (comma_needed) c.s<<",";
			c.s<<*i++;
		}
		if (length>1) c.s<<"}";
	}
	
	/** @brief Overloaded print function for nice \f$\text{\TeX}\f$ output
	 */
	void print(const print_context &c, unsigned level = 0) const
	{
		if (dynamic_cast<const print_latex*>(&c)==NULL)
		{ 
			ncmul::print(c,level); 
			return;
		}
		string varname;
		vector<string> indices;
		bool leavespace=false;
		for (int i=0;i<this->nops();i++)
		{
			string newvarname=ToString(this->op(i));
			const string alphabet="qwertyuiopasdfghjklzxcvbnmQWERTYUIOPASDFGHJKLZXCVBNM\\";
			int first_non_letter=newvarname.find_first_not_of(alphabet);
			if (first_non_letter==string::npos || first_non_letter==0)	//not in symbol+index form
			{
				if (leavespace) c.s<<"\\,";
				if (!varname.empty())
				{
					PrintVariableAndIndices(varname, indices,c);
					varname.clear();
					indices.clear();					
					c.s<<"\\,";
				}
				c.s<<newvarname;
				leavespace=true;
			}
			else {
				string alphapart=newvarname.substr(0,first_non_letter);
				if (alphapart!=varname)
				{
					if (!varname.empty())
					{					
						if (leavespace) c.s<<"\\,";
						PrintVariableAndIndices(varname, indices,c);
						leavespace=true;
					}
					varname=alphapart;
					indices.clear();
					indices.push_back(newvarname.substr(first_non_letter));					
				}
				else					
					indices.push_back(newvarname.substr(first_non_letter));				
			}
		}
		if (!varname.empty())
		{
			if (leavespace) c.s<<"\\,";
			PrintVariableAndIndices(varname, indices,c);
		}		
	}
}; 

// ****************************************************************************
// *							Template Specializations for wedgebase.h						  *
// ****************************************************************************
namespace internal {
	
template<typename T> struct SetOf<Lambda1<T> >
{
	typedef set<ex,ex_is_less> type;
};

//this definition causes a compile-time error if StoreAs differs from ex; this is intended because T and Lambda<T> are unrelated types.
template<typename T,typename StoreAs> struct FindSymbolsHelper<Lambda<T>, StoreAs > : public FindSymbolsHelper <Lambda1<T>,ex >, public T::visitor, public Lambda<T>::visitor
{
	void visit(const typename VisitorType<T>::Type& x) {
		this->occurrences.insert(x);
	}
	void visit(const Lambda<T>& x) {
		this->occurrences.insert(x);
	}
};

template<typename T> struct FindSimpleHelper<Lambda<T> > : public FindSimpleHelper <Lambda1<T> >, public T::visitor, public Lambda<T>::visitor
{
	void visit(const typename VisitorType<T>::Type& x) {
		this->occurrences.insert(x);
	}
	void visit(const Lambda<T>& x) {
		this->occurrences.insert(x);
	}
};

template<typename V> struct LogicalType<Lambda1<V> > {
	static inline bool Is(ex x) {return is_a<V>(x) || is_a<Lambda<V> >(x);}
};

template<typename T> struct NormalFormHelper<Lambda<T> > : public NormalFormHelper<Poly1<Lambda1<T> > >, public T::visitor, public Lambda<T>::visitor
{
	void visit(const typename VisitorType<T>::Type& x) {
		++this->coeffs[x];
	}
	void visit(const Lambda<T>& x) {
		++this->coeffs[x];
	}
};


}


// ****************************************************************************
// *							Template Specializations for operators.h						  *
// ****************************************************************************
/** @internal
@brief Template specialization \see Lambda

@note There is no need to specialize LinearOperator, since it derives from AdditiveOperator
*/ 
template<typename V> class AdditiveOperator<Lambda<V> > : public RecursiveVisitor<ex>, public basic::visitor ,public  add::visitor, public V::visitor, public Lambda<V>::visitor {
public:
	enum {ExpectsExpandedInput=0};
	typedef V OperatesOn;
	virtual ~AdditiveOperator() {}
	void visit(const add& alpha) {
		Result()=0;
		for (unsigned i=0;i<alpha.nops();i++)
			Result()+=this->RecursiveVisit(alpha.op(i));
	}
	void visit(const basic&) {Result()=0;}
	ex RecursiveVisit(const GiNaC::ex& e)	//overloaded for efficiency 
	{	
		ex oldresult=GetResult();
		if (is_exactly_a<add>(e))
			visit(ex_to<add>(e)); 
		else if (is_exactly_a<V>(e) || is_exactly_a<Lambda1<V> >(e))
			static_cast<typename V::visitor*>(this)->visit(ex_to<V>(e));			
		else if (is_exactly_a<Lambda<V> >(e))
			static_cast<typename Lambda<V>::visitor*>(this)->visit(ex_to<Lambda<V> >(e));
		else
			e.accept(*this);		
		ex newresult=GetResult();
		Result()=oldresult;
		return newresult;		
	}
};

namespace internal {
/** @internal
@brief Template specialization \see Lambda
*/	
template<typename LeftType,typename W,typename Bilinear> 
	class MyRightOperator<LeftType,AdditiveOperator<Lambda<W> >,Bilinear> :
 		public AdditiveOperator<Lambda<W> >	
{
	const Bilinear* bil;
	const LeftType& pv;
	typedef AdditiveOperator<Lambda<W> > RightOperator;
public:
	MyRightOperator(const LeftType& v,const Bilinear* bil) 	: pv(v)
	{
		this->bil=bil;
	} 
	void visit(const W& w) {
		this->Result()=bil->Apply(pv,w);
	}
	void visit(const Lambda<W>& w) {
		this->Result()=bil->Apply(pv,w);
	}
};

/** @internal
@brief Template specialization \see Lambda
*/
template<typename LeftType,typename W,typename Bilinear> 
	class MyRightOperator<LeftType,LinearOperator<Lambda<W> >,Bilinear > :
 		public LinearOperator<Lambda<W> >	
{
	const Bilinear* bil;
	const LeftType& pv;
	typedef LinearOperator<Lambda<W> > RightOperator;
public:
	MyRightOperator(const LeftType& v,const Bilinear* bil) : pv(v)
	{
		this->bil=bil;
	} 
	void visit(const W& w) {
		this->Result()=bil->Apply(pv,w);
	}
	void visit(const Lambda<W>& w) {
		this->Result()=bil->Apply(pv,w);
	}
};

/** @internal
@brief Template specialization \see Lambda
*/
template<typename V, typename RightOperator, typename Bilinear> 
	class MyLeftOperator<AdditiveOperator<Lambda<V> >,RightOperator, Bilinear > :
		public AdditiveOperator<Lambda<V> > 
{
	typedef AdditiveOperator<Lambda<V> > LeftOperator;	
	ex w;
	const Bilinear* bil;
public:
	MyLeftOperator(ex w,const Bilinear* bil) {this->w=w; this->bil=bil;}
	void visit(const V& v) {
		MyRightOperator<V,RightOperator,Bilinear> oper(v,bil);			
		w.accept(oper);
		this->Result()=oper.GetResult();
	}
	void visit(const Lambda<V>& v) {
		MyRightOperator<Lambda<V>,RightOperator,Bilinear> oper(v,bil);			
		w.accept(oper);
		this->Result()=oper.GetResult();
	}	
};

/** @internal
@brief Template specialization \see Lambda
*/
template<typename V, typename RightOperator, typename Bilinear> 
	class MyLeftOperator<LinearOperator<Lambda<V> >,RightOperator, Bilinear > :
		public LinearOperator<Lambda<V> > 
{
	typedef LinearOperator<Lambda<V> > LeftOperator;	
	ex w;
	const Bilinear* bil;
public:
	MyLeftOperator(ex w,const Bilinear* bil) {this->w=w; this->bil=bil;}
	void visit(const V& v) {
		MyRightOperator<V,RightOperator,Bilinear> oper(v,bil);			
		w.accept(oper);
		this->Result()=oper.GetResult();
	}
	void visit(const Lambda<V>& v) {
		MyRightOperator<Lambda<V>,RightOperator,Bilinear> oper(v,bil);			
		w.accept(oper);
		this->Result()=oper.GetResult();
	}	
};

}

/** @internal @brief Template specialization \see Lambda
*/	

template<typename V>  class TrivialPairingOperator<Lambda<V> > : public IBilinearOperator<
	LinearOperator<Lambda<V> >,
	LinearOperator<Lambda<V> >
  >
{
public:
/**  
 * @remark v and w may reference objects with different types.
 */ 
	ex Apply(const V& v, const V& w) const
	{	
		return v==w? 1 : 0;			
	}
	ex Apply(const Lambda<V>& v, const Lambda<V>& w) const
	{
		return v==w? 1 : 0;
	}
	ex Apply(const V& , const Lambda<V>& ) const {return 0;}
	ex Apply(const Lambda<V>&,const V&) const {return 0;}		
};

/** @brief Exception thrown by degree() when the element is not homogeneous
 */
class InhomogeneousExpression : public WedgeException<std::invalid_argument>
{
public:
	InhomogeneousExpression(const char* in_file, int at_line) :  WedgeException<std::invalid_argument>("Inhomogeneous expression",in_file,at_line) {}
};


namespace internal {
class IntMod2
{
	bool x;
public:
	IntMod2& operator=(unsigned n) {x=(n%2)!=0; return *this;}
	IntMod2& operator+=(IntMod2 n) {x=x xor n.x; return *this;}
	IntMod2& operator*=(unsigned n) {x=x && (n%2)!=0; return *this;}
	operator bool() const {return x;}
};


/** @brief Visitor class to compute the degree of a polynomial in V or a linear combination of lambda<V>'s
 */
template<typename T, typename DegreeType> class ComputeDegree: public RecursiveVisitor<DegreeType>,public power::visitor,public T::visitor,public mul::visitor,public add::visitor,public basic::visitor
{
	void visit(const Lambda<T> & x) {
		this->Result()=x.degree();
	}
	void visit(const power& x) {
		this->Result()=this->RecursiveVisit(x.op(0));
		if (this->Result()==0) return;
		if (is_a<numeric>(x.op(1))) {
			const numeric& exponent=ex_to<numeric>(x.op(1));
			if (exponent.is_integer()) {
				this->Result()*=exponent.to_int();
				return;
			}
		}
		LOG_MSG("ERROR: non-integer exponent");
		LOG_ERROR(x);
		throw InvalidArgument(__FILE__,__LINE__,x);
	}
	void visit(const mul & x) {
		this->Result()=0;
		for (unsigned i=0;i<x.nops();i++)
			this->Result()+=this->RecursiveVisit(x.op(i));
	}
	void visit(const typename  VisitorType<T>::Type&) {
		this->Result()=1;
	}

	void visit(const add & x) {
		unsigned i=0;
		this->Result()=this->RecursiveVisit(x.op(i++));
		while (i<x.nops())
			if (this->RecursiveVisit(x.op(i++))!=this->Result()) throw InhomogeneousExpression(__FILE__,__LINE__);
	}
	void visit(const basic&) {
		this->Result()=0;
	}
};

template<typename T> class ComputeDegree<Poly<T>, int> : public ComputeDegree<T,int>
{
	void visit(const add & x) {
		int maxdegree=-1;
		for (unsigned i=0;i<x.nops();++i)
			maxdegree=max(maxdegree,this->RecursiveVisit(x.op(i)));
		assert(maxdegree>=0);
		this->Result()=maxdegree;
	}	
};

//specialization for Lambda<V> case
template<typename T, typename DegreeType> class ComputeDegree<Lambda<T>,DegreeType > : public RecursiveVisitor<DegreeType>,public Lambda<T>::visitor,public T::visitor,public mul::visitor,public add::visitor,public basic::visitor
{
	void visit(const Lambda<T> & x) {
		this->Result()=x.degree();
	}
	void visit(const mul & x) {
		this->Result()=0;
		for (unsigned i=0;i<x.nops();i++)
			this->Result()+=this->RecursiveVisit(x.op(i));
	}
	void visit(const typename VisitorType<T>::Type&) {
		this->Result()=1;
	}

	void visit(const add & x) {
		unsigned i=0;
		this->Result()=this->RecursiveVisit(x.op(i++));
		while (i<x.nops())
			if (this->RecursiveVisit(x.op(i++))!=this->Result()) throw InhomogeneousExpression(__FILE__,__LINE__);
	}
	void visit(const basic&) {
		this->Result()=0;
	}
};


}


template<typename T> int Degree(ex e)
{
	internal::ComputeDegree<T,int> v;
	e.expand().accept(v);
	return v.GetResult();
}

template<typename T> bool IsOdd(ex e)
{
	internal::ComputeDegree<T,internal::IntMod2> v;
	e.expand().accept(v);
	return v.GetResult();
}


template<typename V> ex Lambda1<V>::eval_ncmul(const exvector & v) const
{
	exvector::const_iterator i=v.begin();
	while (i!=v.end() && is_a<Lambda1<V> >(*i))
	 	i++;
	if (i!=v.end()) return hold_ncmul(v).expand();	//may contain sums
 
	exvector s=v;
	bool change_sign=false;
	
	exvector::iterator end=s.end();
	while (s.begin()!=end)
	{			
		exvector::iterator j = s.begin() + 1;		
		while (j!=s.end())
		{
			int cmp=(j-1)->compare(*j);
			if (cmp==0)	return 0;
			else if (cmp>0) {
				j->swap(*(j-1));
				change_sign=not change_sign;
				}
			j++;
		}
		end--;
	}
	
	ex result=ex(
		(new Lambda<V>(s))->setflag(status_flags::dynallocated |status_flags::evaluated)
	);
	if (change_sign) result=-result;	
	return result;
}


} /** @} */

#include "wedge/linearalgebra/tensor.h"
#include "wedge/linearalgebra/tensorlambda.h"

#endif /*LAMBDA_H_*/
