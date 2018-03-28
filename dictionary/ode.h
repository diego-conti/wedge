/***************************************************************************
 *   Copyright (C) 2009 by Diego Conti, diego.conti@unimib.it              *
 *                                                                         *
 *   This file is part of Wedge.	                                   *
 *                                                                         *
 *   Wedge is free software; you can redistribute it and/or modify    *
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
#ifndef ODE_H_
#define ODE_H_

#include <wedge/manifold.h>
#include <wedge/liegroup.h>
#include <wedge/function.h>
#include <wedge/logging.h>
#include <wedge/leibniz.h>
#include <wedge/lambda.h>
#include <wedge/latex.h>
#include <wedge/polybasis.h>
//#include <wedge/cocoapolyalg.h>
#include <iterator>
#include <cctype>
#include <wedge/simplifier.h>

namespace Wedge {



/** @brief Helper class that applies the de l'Hopital rule to a fraction f/g
 * 
 * @warning This code has not been tested sufficiently and it has given rise to some inconsistent behaviour
*/
struct Limit {
	ex fx, x, f0, x0;
	const Simplifier& simplifier_;
	ex at0(ex f) {
		return simplifier_.Simplify(f.subs(f0,subs_options::algebraic).subs(x==x0).expand().subs(f0.subs(x==x0),subs_options::algebraic));
	}
	operator ex() {
		ex n=simplifier_.Simplify(fx.numer());
		ex d=simplifier_.Simplify(fx.denom());
		LOG_INFO(n);
		LOG_INFO(d);
		while (at0(d).is_zero() && at0(n).is_zero())
		{			
			n=diff(n,ex_to<symbol>(x)).expand();
			d=diff(d,ex_to<symbol>(x)).expand();
			if (d.is_zero()) break;
			ex x=simplifier_.Simplify(n/d).normal();
			n=x.numer();
			d=x.denom();
			n=simplifier_.Simplify(n);
			d=simplifier_.Simplify(d);
			LOG_INFO(n);
			LOG_INFO(d);
		}
		LOG_INFO(at0(n));
		LOG_INFO(at0(d));
		return (at0(n)/at0(d)).normal();
	}
/** @brief Construct a limit object
 * @param _fx A function object, e.g. sin(r)
 * @param _x The variable and the limit, e.g r==0
 * @param _f0 Relations satisfied by the function at the limit point, e.g. D[sin](r)=1
 */
	Limit(ex _fx, ex _x,  ex _f0, const Simplifier& simplifier=default_simplifier) : fx(_fx), f0(_f0), simplifier_(simplifier) {
		if (!_x.info(info_flags::relation_equal)) {
			LOG_ERROR(_x);
			throw InvalidArgument(__FILE__,__LINE__,_x);		
		}
		x=_x.lhs();
		x0=_x.rhs();
		assert(is_a<symbol>(x));
		LOG_INFO(fx);
		fx=simplifier.Simplify(fx).normal();
		LOG_INFO(fx);
	}
};



template<typename ScalarSimplifier> struct ScalarSimplifierClass {
	typedef ScalarSimplifier Simplifier;
};

//isolates "function" objects.
template<typename Simplifier> struct FunctionSimplifier : public ScalarSimplifierClass<Simplifier>::Simplifier {
	typedef typename ScalarSimplifierClass<Simplifier>::Simplifier ScalarSimplifier;
	ex Simplify(ex e) const {
		list<ex> functions;
		LOG_INFO(e);
		GetSymbols<function>(functions,e);
		lst subs;
		for (list<ex>::const_iterator i=functions.begin();i!=functions.end();++i) {
			for (int n=10;n>=1;--n) {
				ex lhs=pow(*i,wild(0));
				ex rhs=wild(0);
				for (int k=1;k<=n;++k) {
					lhs*=pow(*i,wild(k));
					rhs+=wild(k);
				}
				subs.append(lhs==pow(*i,rhs));
			}			
			//subs.append(pow(*i,wild(0))* *i==pow(*i,wild(0)+1));
		}
		e=ScalarSimplifier::Simplify(e);
		e=e.expand().subs(pow(pow(wild(0),wild(1)),wild(2))==pow(wild(0),wild(1)*wild(2)));
		e=ScalarSimplifier::Simplify(e);
		e=e.normal().subs(subs,subs_options::algebraic).expand().normal();
		for (list<ex>::const_iterator i=functions.begin();i!=functions.end();++i) 
			e=e.normal().subs(pow(*i,wild(0))* *i==pow(*i,wild(0)+1));
		e=e.expand().subs(subs,subs_options::algebraic);
		e=e.expand().subs(subs,subs_options::algebraic);
		return e;
	}
};

//ensures that FunctionSimplifier<FunctionSimplifier<S>> is equivalent to FunctionSimplifier<S>
template<typename ScalarSimplifier> struct ScalarSimplifierClass<FunctionSimplifier<ScalarSimplifier> > {
	typedef ScalarSimplifier Simplifier;
};

template<typename ParameterType, typename Iterator>  ex GenericPolynomial(Iterator begin, Iterator end, int degree, char parname='A') {
	ex result;
	vector<int> degrees(static_cast<int>(distance(begin,end)));	
	ex monomial=1;
	int n=0;	//degree of the current monomial
	int parno=0;
	while (true)
	{
		result+=ParameterType(parname+ToString(++parno))*monomial;
		Iterator i=begin;
		vector<int>::iterator j=degrees.begin(); //parallel iterators
		//increment
		while (n==degree && i!=end) {
			n-=*j;
		monomial=(monomial/pow(*i,*j)).expand();
		*j=0;
			++i; 
			++j;
		}
		if (i==end) break;
		//update monomial
		++*j;
		monomial*=*i;
		++n;
	}
	LOG_DEBUG(result);
	return result;	
}

template<typename ParameterType, typename Iterator>  ex GenericRational(Iterator begin, Iterator end, int degnum, int degden)
{
	return GenericPolynomial<ParameterType>(begin,end,degnum,'A')/GenericPolynomial<ParameterType>(begin,end,degden,'B');
}

WEDGE_DECLARE_NAMED_ALGEBRAIC(AnsatzParameter,symbol)

class Ansatz : public HasParameters<AnsatzParameter> {
	friend ostream& operator<<(ostream& , const Ansatz& );
	//a bunch of equations
	ExVector unknowns;	//unknowns, which can be of type function, symbol or whatever
	ExVector values;	//the values assigned to the unknowns, depending on constant parameters of type AnsatzParameter
	list<ex> equations;	//algebraic equations in the AnsatzParameter, no RHS
public:
	template<typename Iterator> Ansatz(Iterator unknowns_begin, Iterator unknowns_end, ex ansatz) : unknowns(unknowns_begin,unknowns_end) {
		LOG_INFO(ansatz);
		exvector parameters;
		GetSymbols<AnsatzParameter>(parameters, ansatz);
		values.resize(unknowns.size());
		
		for (int j=0;j<unknowns.size();++j)
		{
			ex f=ansatz;			
			for (int i=0;i<parameters.size();++i)
			{
				ex a=AnsatzParameter(ToString(parameters[i])+ToString(j+1));
				f=f.subs(parameters[i]==a);
			}
			values[j]=f;
		}
	}
/** 

 * If you use this constructor, you can set each function individually by using, say, DeclareZero(a-AnsatzParameter()) with a unknown
*/
	template<typename Iterator> Ansatz(Iterator unknowns_begin, Iterator unknowns_end) : unknowns(unknowns_begin,unknowns_end), values(unknowns) {}
	ex Value (ex unknown) const {
		lst subs;
		assert(unknowns.size()==values.size());
		for (int i=0;i<unknowns.size();++i)
			subs.append(unknowns[i]==values[i]);
		return unknown.subs(subs).expand();
	}
	void DeclareConditions(const lst& list_of_equations) {
		for (ExVector::iterator i=values.begin();i!=values.end();++i)
			*i=i->subs(list_of_equations);	
		for (list<ex>::iterator i=equations.begin();i!=equations.end();++i)
			*i=i->subs(list_of_equations);
	}
	void AddEquation(ex eqn) {
		equations.push_back(eqn);
	}
	template<typename Iterator> void AddEquations(Iterator begin, Iterator end) {
		equations.insert(equations.end(),begin,end);
	}
	ex Equations() const {return lst(equations);}
	void Simplify(const Simplifier& simplifier) {
		for (list<ex>::iterator i=equations.begin(); i!=equations.end();++i)
			*i=simplifier.ScalarEquationSimplify(*i);
	}
	template<typename Iterator> static ex PolynomialAnsatz(Iterator begin, Iterator end, int degree) {
		char parname='A';
		ex result;
		vector<int> degrees(static_cast<int>(distance(begin,end)));	
		ex monomial=1;
		int n=0;	//degree of the current monomial
		while (true)
		{
			result+=AnsatzParameter(ToString(parname++))*monomial;
			if (parname=='Z'+1) parname='a';
			if (parname=='z'+1) throw WedgeException<std::runtime_error>("Too many parameters", __FILE__,__LINE__);
			Iterator i=begin;
			vector<int>::iterator j=degrees.begin(); //parallel iterators
			//increment
			while (n==degree && i!=end) {
				n-=*j;
				monomial=(monomial/pow(*i,*j)).expand();
				*j=0;
				++i; 
				++j;
			}
			if (i==end) break;
			//update monomial
			++*j;
			monomial*=*i;
			++n;
		}
		LOG_DEBUG(result);
		return result;	
	}
/** @brief Alter the ansatz object by applying a given substitution to both equations and unknowns
 * @param substitution A substitution of the form x==2, or an lst of such
 * @warning When dealing with ODE's, use Ode::SetAnsatz instead, which also takes care of replacing the derivative.
 */
	void subs(ex substitution) {
		LOG_INFO(substitution);
		if (is_a<lst>(substitution)) {
			lst newsubs;
			const lst& s=ex_to<lst>(substitution);
			for (lst::const_iterator i=s.begin();i!=s.end();++i)
			{
				if (!is_a<relational>(*i)) {
					LOG_ERROR(*i);
					throw InvalidArgument(__FILE__,__LINE__,substitution);
				}
				newsubs.append(i->lhs()==Value(i->rhs()));			
			}
			substitution=newsubs;
		}
		else if (is_a<relational>(substitution))
		{
			substitution=(substitution.lhs()==Value(substitution.rhs()));
		}
		else throw InvalidArgument(__FILE__,__LINE__,substitution);
		LOG_INFO(substitution);
		for (exvector::iterator i=values.begin();i!=values.end();++i)
			*i=i->subs(substitution).expand();
		for (list<ex>::iterator i=equations.begin();i!=equations.end();++i)
			*i=i->subs(substitution).expand();
	} 
	
	bool isPopulated() const {return !equations.empty();}
};

ostream& operator<<(ostream& os, const Ansatz& ansatz)
{
	if (IsStreamTeX(os)) {
		for (int i=0;i<ansatz.unknowns.size();++i)
			os<<Equation(ToString(ansatz.unknowns[i]),ansatz.values[i]);
		for (list<ex>::const_iterator i=ansatz.equations.begin();i!=ansatz.equations.end();++i)
			os<<Equation(*i);
	}
	else {
		for (int i=0;i<ansatz.unknowns.size();++i)
			os<<ansatz.unknowns[i]<<"="<<ansatz.values[i]<<endl;
		os<<ansatz.equations<<endl;
	}
	return os;
}


/** @brief Base class for quasilinear ODE's
*/
class ODE  {
protected:
	ExVector unknowns;	//unknowns, which are functions of r; have type function (or derived class)
	ex r;			//variable; has type symbol (or derived class)
public:
/** @brief Construct a ODE container where the variable is the given variable */ 
	ODE(ex variable) : r(variable) {assert(is_a<symbol>(r));}
/** @brief Substitute a specific ansatz in the equations 
 * @param ansatz An ansatz
 * @return The equations (no RHS), with the ansatz replaced in
*/
	virtual ex Equations(const Ansatz& ansatz) const =0;
/** @brief Return the lst of equations
 * @return The equations (no RHS)
*/
	virtual ex Equations() const =0;
/** @brief Alter an Ansatz object by declaring an unknown to be of a certain form
 * @param ansatz An ansatz object
 * @param substitution A substitution of the form f(x)==cos(x)
 */
	Ansatz& SetAnsatz(Ansatz& ansatz, ex substitution) const {
		if (!is_a<relational>(substitution)) throw InvalidArgument(__FILE__,__LINE__,substitution);
		exvector::const_iterator i=unknowns.begin();
		while (i!=unknowns.end() && *i!=substitution.lhs()) ++i;
		if (i==unknowns.end()) throw InvalidArgument(__FILE__,__LINE__,substitution);
		ansatz.subs(substitution);
		ansatz.subs(i->diff(ex_to<symbol>(r))==substitution.rhs().diff(ex_to<symbol>(r)));		
		return ansatz;	
	}
	
/** @brief Populate a Ansatz object with the relevant equations
 * @param ansatz An Ansatz object whose unknowns correspond to this ODE's unknowns
 *
 * This function substitutes the ansatz in the unknowns.
 */
	template<typename Simplifier> Ansatz& GetEquations(Ansatz& ansatz) const {
		FunctionSimplifier<Simplifier> simplifier;
		ex neweqns=Equations(ansatz);
		LOG_INFO(neweqns);		
		for (int i=0;i<neweqns.nops();++i)
			ansatz.AddEquation(simplifier.Simplify(neweqns.op(i).numer()));
		return ansatz;
	}
/** @brief Create a (populated) Ansatz object corresponding to the corresponding ansatz
 * @param ansatz A function of r, depending on (constant) parameters of type AnsatzParameter
 *
 * This function substitutes the ansatz in each unknown, using different parameters for each function.
 */	
	template<typename Simplifier> Ansatz MakeAnsatz(ex ansatz) const {
		Ansatz result(unknowns.begin(),unknowns.end(),ansatz);
		return GetEquations<Simplifier>(result);
	}
	
	ex variable() const {return r;}

	virtual ex derivative(ex e) {throw NotImplemented(__FILE__,__LINE__);}
/** @brief Compute the derivatives of the functions at a point, given the initial condition
 * @param r0 The value of r at which the derivatives must be computed
 * @param IC A lst of substitutions of the form (f(r)==1) etc.
 * @param n The order up to which the derivatives must be computed
 * @return A lst of substitutions of the form (f(r)==1, f'(r)==0) etc.
 */	
	ex ComputeDerivativesUpTo(ex r0, ex IC, int n=6) const {
		ex subs;
		if (is_a<lst>(IC)) subs=IC;
		else subs=lst(IC);
		for (int i=1;i<=n;++i)
			subs=ComputeDerivativesAt(r0,subs,i);
		return subs;
	}


/** @brief Compute n-th derivatives of the functions at a point, given the derivatives up to order n-1
 * @param r0 The value of r at which the derivatives must be computed
 * @param IC A lst of substitutions of the form (f(r)==1, f'(r)==0) etc.
 * @param n The order of the derivatives to be computed
 * @return A lst of substitutions of the form (f(r)==1, f'(r)==0, f''(r)==0) etc.
 */	
	ex ComputeDerivativesAt(ex r0, ex IC, int n=1) const {
		lst Dsubs=ex_to<lst>(IC);	//a list of substitutions of the form f(0)==0, f'(0)==1, etc.
		lst symbols;
		lst subs_symbols_for_derivatives;
		for (int i=0;i<unknowns.size();++i) {
			ex a=symbol("D^"+ToString(n)+ToString(unknowns[i]));
			symbols.append(a);
			subs_symbols_for_derivatives.append(unknowns[i].diff(ex_to<symbol>(r),n).subs(r==r0)==a);
		}
		ex eqns=Equations();
		LOG_INFO(eqns);
		lst neweqns;

		for (int i=0;i<eqns.nops();++i)
		{
			list<ex> derivatives;
			int m=n-1;
			while (true)
			{
				ex eqn=diff(eqns.op(i).numer(),ex_to<symbol>(r),m);
				eqn=eqn.subs(Dsubs).subs(r==r0);
				eqn=eqn.subs(subs_symbols_for_derivatives);

				GetSymbols<fderivative>(derivatives,eqn);				
				//if the equation involves other (higher) derivatives than the ones we want to compute, ignore it and stop taking derivatives
				if (!derivatives.empty()) break;
				neweqns.append(eqn==0);
				++m;
			}
		}
		LOG_INFO(neweqns);
		ex sol=lsolve(neweqns,symbols);	//these equations may not be linear...
		if (sol==lst()) {
			LOG_ERROR(neweqns);
			LOG_ERROR(symbols);
			throw WedgeException<std::runtime_error>("Cannot compute the derivatives", __FILE__,__LINE__);
		}
		for (int i=0;i<unknowns.size();++i) 
			Dsubs.append(unknowns[i].diff(ex_to<symbol>(r),n)==symbols[i].subs(sol));		
		LOG_INFO(Dsubs);		
		return Dsubs;
	}
/** @brief Change the variable
 * @param x A new independent variable
 * @param eq An equation of the form r==2*x+1 etc.
 *
 * Replaces the variable x with the variable r
 */
	void ChangeVariable(ex x, ex eq) {
		assert(is_a<symbol>(x));
		assert(eq.lhs()==r);
		ex Dr=eq.rhs().diff(ex_to<symbol>(x));	// = dr/dx
		lst subs;
		lst fsubs;
		ex neweqns=Equations();		//this is a virtual function, so we call it before changing anything in *this
		for (exvector::iterator i=unknowns.begin();i!=unknowns.end();++i)
		{
			subs.append(i->diff(ex_to<symbol>(r))==i->subs(r==x).diff(ex_to<symbol>(x))/Dr); // df/dr (r) -> (df/dx)(x) /(dr/dx) (x)
			fsubs.append(*i==i->subs(r==x));
			*i=i->subs(r==x);	//f(r) -> f(x)
		}
		neweqns=neweqns.subs(fsubs).subs(subs).subs(eq);
		assert(is_a<lst>(neweqns));
		r=x;
		SetEquations(neweqns);
	}

/** @brief Change the variable
 * @param x A new independent variable
 * @param eq An equation of the form cos(r)==x
 *
 * Replaces the variable x with the variable r
 */
	void ChangeVariableInverse(ex x, ex eq, const Simplifier& simplifier=default_simplifier) {
		assert(is_a<symbol>(x));
		assert(eq.rhs()==x);
		ex Dx=eq.lhs().diff(ex_to<symbol>(r));	// = dx/dr
		lst subs(eq);
		lst fsubs;
		ex neweqns=Equations();		//this is a virtual function, so we call it before changing anything in *this
		for (exvector::iterator i=unknowns.begin();i!=unknowns.end();++i)
		{
			subs.append(i->diff(ex_to<symbol>(r))==i->subs(r==x).diff(ex_to<symbol>(x))*Dx); // df/dr (r) -> (df/dx)(x) *(dx/dr) (r)
			fsubs.append(*i==i->subs(r==x));
			*i=i->subs(r==x);	//f(r) -> f(x)
		}
		LOG_INFO(subs);
		LOG_INFO(neweqns);
		neweqns=neweqns.subs(subs,subs_options::algebraic).normal();
		LOG_INFO(neweqns);
		neweqns=simplifier.Simplify(neweqns);
		LOG_INFO(neweqns);
		neweqns=neweqns.subs(eq,subs_options::algebraic).normal();		
		LOG_INFO(neweqns);
		neweqns=neweqns.subs(fsubs);	//f(r) -> f(x)
		//@todo in principle one might obtain equations of the form f(r)g(x), so one can reduce to g(x) assuming f(r) is not identically zero
		ex to_eliminate=r;
		r=x;	//update object before invoking virtual functions
		SetEquations(neweqns);
		LOG_INFO(Equations());
		if (Equations().has(to_eliminate)) {
			LOG_ERROR(x);
			LOG_ERROR(eq);
			LOG_ERROR(neweqns);
			throw WedgeException<std::runtime_error>("Cannot eliminate r", __FILE__,__LINE__);
		}

	}
/** @brief Change an unknown
 * @param f A new unknown
 * @param eq An equation of the form g(r)==sin(f(r))+f etc.
 *
 * Replaces the unknown g with f
 */
	void ChangeUnknown(ex f, ex eq) {
		ex equations=Equations();	//invoke virtual function before altering object...
		assert(is_a<relational>(eq));
		exvector::iterator i=unknowns.begin();
		while (i!=unknowns.end() && *i!=eq.lhs()) ++i;
		if (i==unknowns.end()) throw InvalidArgument(__FILE__,__LINE__,eq);
		else unknowns.erase(i); 

		lst subs(eq);
		subs.append(diff(eq.lhs(),ex_to<symbol>(r))==diff(eq.rhs(),ex_to<symbol>(r)));
		unknowns.push_back(f);
		LOG_INFO(subs);		
		SetEquations(equations.subs(subs));
	}

	typedef exvector::const_iterator const_iterator;
	const_iterator unknowns_begin() const {return unknowns.begin();}
	const_iterator unknowns_end() const {return unknowns.end();}
protected:
/** @brief Change the equations
 * @param eqns A lst of equations (no RHS)
*/
	virtual void SetEquations(ex eqns)=0;
};


/** @brief Container holding a list of first order quasilinear ODE's.
*/
class QuasilinearODE : public ODE {
protected:
	lst eqns;		//equations
public:
/** @brief Construct a ODE container where the variable is the given variable */ 
	QuasilinearODE(ex variable) : ODE(variable) {}
	void AddEquation(ex eqn) {
		eqns.append(eqn);
	}
	ex Equations(const Ansatz& ansatz) const {
		lst subs;
		for (exvector::const_iterator i=unknowns.begin();i!=unknowns.end();++i)
		{
			subs.append(*i==ansatz.Value(*i));
			subs.append(i->diff(ex_to<symbol>(r))==ansatz.Value(*i).diff(ex_to<symbol>(r)));
		}
		return ex(eqns).subs(subs).expand();
	}
	ex Equations() const {
		return eqns;
	}
	void AddUnknown(ex f) {unknowns.push_back(f);}
protected:
	void SetEquations(ex neweqns) {
		assert(is_a<lst>(neweqns));
		eqns=ex_to<lst>(neweqns);
	}
};




/** @brief Container holding a list of first order linear ODE's.
*/
class FirstOrderODE : public ODE {
protected:
	ExVector derivatives;	//first derivatives of the functions.
	lst subs;		//substitutions corresponding to the equations
public:
/** @brief Construct a ODE container where the variable is the given variable */ 
	FirstOrderODE(ex variable) : ODE(variable) {}
/** @brief Add the equation corresponding to D(f)=DF 
 * @param f An unknown, such as f(r), where f is of type function
 * @param Df An expression depending on r, f(r)...
*/
	void Declare_D(ex f, ex Df) {
		if (!is_a<function>(f) || f.nops()!=1 || f.op(0)!=r)
		{
			LOG_ERROR(f);
			throw InvalidArgument(__FILE__,__LINE__,f);
		}
		if (find(unknowns.begin(),unknowns.end(),f)!=unknowns.end())
		{
			LOG_ERROR(subs);
			LOG_ERROR(f);
			LOG_ERROR(Df);
			throw InconsistentDeclaration(__FILE__,__LINE__,"same function appears twice on left hand side");
		}
		unknowns.push_back(f);
		derivatives.push_back(Df); 
		subs.append(diff(f,ex_to<symbol>(r))==Df);
	}

/** @brief Compute the derivative of some quantity depending on the unknowns and r
*/
	ex derivative(ex e) const	
	{
		LOG_DEBUG(e);
		ex de=e.diff(ex_to<symbol>(r));
		LOG_DEBUG(de);
		de=de.subs(subs);
		LOG_DEBUG(de);
		return de.expand().normal();
	}

	ex linearize(ex IC) const
	{
		lst result;
		for (int i=0;i<unknowns.size();++i)
		{
			ex Xi;
			for (int j=0;j<unknowns.size();++j)
			{
				symbol s("s");
				ex d=derivatives[i].subs(unknowns[j]==s).diff(s).subs(s==unknowns[j]);
				Xi+=d.subs(IC)*unknowns[j];
			}
			result.append(Xi);			
		}
		return result;	
	}

/** @brief Find monomial quantities in the unknowns whose derivatives has a simple expression
 *
 * In principle one could take a polynomial depending on symbols, then obtain a linear map, and look for eigenvalues
*/	
	void Reduce(int maxexp=5, int minexp=0) {
		if (minexp>maxexp) throw InvalidArgument(__FILE__,__LINE__,minexp);
		set<ex,ex_is_less> f(unknowns.begin(),unknowns.end());
//this is all the functions, so if I have sin(r), cos(r) and such in the expression they also appear
		for (int i=0;i<derivatives.size();++i)
			GetSymbols<function> (f,derivatives[i]);
		ExVector functions(f.begin(),f.end());
/*		int fsize=functions.size();
		for (int i=0;i<fsize;++i)
			for (int j=i+1;j<fsize;++j)
				functions.push_back(functions[i]-functions[j]);
*/		
		LOG_INFO(functions);
		vector<int> exponents(functions.size());
		for (int i=0;i<functions.size();++i)
			exponents[i]=minexp;
		
		ex m,dm;
		lst symbols;
		while (true) {
			ex monomial=1;
			for (int i=0;i<functions.size();++i)
				monomial*=pow(functions[i],exponents[i]);
			LOG_INFO(monomial);
			ex s=symbol();
			m+=s*monomial;
			ex der=derivative(monomial);
			dm+=(s*der).expand();
			{
				list<ex> coeffs;
				GetCoefficients<Poly<function> > (coeffs, der);
				list<ex>::const_iterator i=coeffs.begin();
				if (i!=coeffs.end()) {
					bool sign=*i>0;
					while (++i!=coeffs.end() && (sign xor static_cast<bool>(*i<0))) ;						
					if (i==coeffs.end()) {
						LOG_INFO(coeffs);
						LOG_INFO(der);
						cout<<Equation("m",monomial);
						cout<<Equation("m'",der);
					}
				}
			
			}
			symbols.append(s);			
			//if (is_a<numeric>(dm.denom()))
//				cout<<Equation("d"+ToString(monomial),dm);

			int r=0;
			++exponents[0];
			while (exponents[r]==maxexp) {
				if (r==functions.size()-1) break;
				exponents[r]=minexp;
				++exponents[++r];
			}
			if (exponents[r]==maxexp) break;
		}
		dm=dm.expand();
		LOG_INFO(m);
		LOG_INFO(dm);
		lst eqns;
		GetCoefficients<Poly<function> > (eqns,dm,withRHS);
		LOG_INFO(eqns);
		ex sol=lsolve(eqns,symbols);
		cout<<Equation(m.subs(sol));
	}

/** @brief Compute derivatives of the functions at a point, given the initial condition
 * @param r0 The value of r at which the derivatives must be computed
 * @param IC A lst of substitutions of the form (f(r)==1) etc.
 
	template<typename Simplifier> ex ComputeDerivativesAt(ex r0, ex IC, int upto=6) const {
		lst Dsubs=ex_to<lst>(IC);	//a list of substitutions of the form f(0)==0, f'(0)==1, etc.
		lst symbols;
		for (int i=0;i<unknowns.size();++i) {
			ex a=symbol("D^nf"+ToString(i));
			symbols.append(a);
		}
		for (int n=1;n<=upto;++n) {
			lst eqns;
			for (int i=0;i<unknowns.size();++i)
				eqns.append(symbols[i]==Limit<Simplifier>(diff(derivatives[i],ex_to<symbol>(r),n-1),r==r0,Dsubs));
			ex eqns2=eqns;
			for (int i=0;i<unknowns.size();++i) 
				eqns2=eqns2.subs(unknowns[i].diff(ex_to<symbol>(r),n).subs(r==r0)==symbols[i]);

			ex sol=lsolve(eqns2,symbols);
			LOG_INFO(eqns2);
			LOG_INFO(symbols); 
			LOG_INFO(sol);
			
			if (sol==lst()) {
				LOG_ERROR(eqns2);
				LOG_ERROR(symbols);
				throw WedgeException<std::runtime_error>("Cannot compute the derivatives", __FILE__,__LINE__);
			}
			for (int i=0;i<unknowns.size();++i) 
				Dsubs.append(unknowns[i].diff(ex_to<symbol>(r),n)==symbols[i].subs(sol));		
			LOG_INFO(Dsubs);
		}
		return Dsubs;
	}
	*/
	ex Equations(const Ansatz& ansatz) const {
		lst subs;
		for (exvector::const_iterator i=unknowns.begin();i!=unknowns.end();++i)
		{
			subs.append(*i==ansatz.Value(*i));
			subs.append(i->diff(ex_to<symbol>(r))==ansatz.Value(*i).diff(ex_to<symbol>(r)));
		}
		return Equations().subs(subs).expand();
	}
/** @brief Return the lst of equations
 * @return The equations (no RHS)
*/
	ex Equations() const {
		lst eqns;
		for (int i=0;i<unknowns.size();++i)
			eqns.append(derivatives[i]-unknowns[i].diff(ex_to<symbol>(r)));
		return eqns;	
	}
	
/** @brief Modify this ODE by an ansatz
* @param ansatz A partially specified ansatz
*/
	void Set(const Ansatz& ansatz) {
		ExVector new_unknowns, new_derivatives;
		for (int i=0;i<unknowns.size();++i)
		{
			ex x=ansatz.Value(unknowns[i]);
			if (unknowns[i]==x) {
				new_unknowns.push_back(unknowns[i]);
				new_derivatives.push_back(ansatz.Value(derivatives[i]));
			}
		}
		unknowns=new_unknowns;
		derivatives=new_derivatives;
		Update();
	}

/** @brief Compute the Jacobian at a point
 * @param point A substitution of the form y==1, etc.
 *
 * If the system has the form y'(t)=F(y,t), this returns the Jacobian of F(., t) at the point */

	matrix J(ex point, int order=0) const {
		matrix result(unknowns.size(),unknowns.size());
		for (int i=0;i<unknowns.size();++i)
		for (int j=0;j<unknowns.size();++j)
		{
			symbol s;
			ex derivative=	derivatives[i].subs(unknowns[j]==s).diff(s).subs(s==unknowns[j]);
			LOG_INFO(derivative);
			result(i,j)=derivative.subs(point).expand();
			if (order==1) {
				ex d;
				for (int k=0;k<unknowns.size();++k)
					d+=derivative.subs(unknowns[k]==s).diff(s).subs(s==unknowns[k]).subs(point).expand()*(unknowns[k]-unknowns[k].subs(point));
				result(i,j)+=d;
			}
			else if (order!=0) throw NotImplemented(__FILE__,__LINE__);
			LOG_INFO(result(i,j));
		}
		return result;
	}
	
protected: //this is used internally to change the equations
	void SetEquations(ex eqns) {
		lst unk;
		lst subs;
		for (int i=0;i<unknowns.size();++i)
		{
			ex s=symbol();
			unk.append(s);
			subs.append(unknowns[i].diff(ex_to<symbol>(r))==s);
		}
			
		assert(is_a<lst>(eqns));
		for (int i=0;i<eqns.nops();++i)
		{
			ex eq=eqns.op(i).subs(subs)==0;
			eqns.let_op(i)=eq;
		}
		ex sol=lsolve(eqns,unk);
		if (sol==lst()) {
			LOG_ERROR(eqns);
			LOG_ERROR(unk);
			throw InvalidArgument(__FILE__,__LINE__,eqns);
		}
		assert(derivatives.size()==unk.nops());
		this->subs.remove_all();
		for (int i=0;i<unknowns.size();++i) {
			derivatives[i]=unk.op(i).subs(sol);
			this->subs.append(unknowns[i].diff(ex_to<symbol>(r))==derivatives[i]);
		}
		LOG_INFO(eqns);
		LOG_INFO(sol);
		LOG_INFO(unknowns);
		LOG_INFO(derivatives);
		LOG_INFO(this->subs);
	}
	//call this after altering unknowns and derivatives
	void Update() {
		subs.remove_all();
		for (int i=0;i<unknowns.size();++i) {
			subs.append(unknowns[i].diff(ex_to<symbol>(r))==derivatives[i]);
		}	
	}

};


ostream& operator<<(ostream& os, const ODE& ode)
{
	ex eqns=ode.Equations();
	if (IsStreamTeX(os)) 
		for (int i=0;i<eqns.nops();++i)
			os<<Equation(eqns.op(i));
	else os<<eqns<<endl;
	return os;
}


}
#endif

