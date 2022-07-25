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
 
#ifndef COCOAPOLYALG_H_
#define COCOAPOLYALG_H_

/** @ingroup ExternalAlgorithms */ 

/** @{ 
 * @file cocoapolyalg.h
 * @brief Standard implementation of polynomial algebra algorithms, using %CoCoA
 * 
 */
#include "CoCoA/library.H"
#include "wedge/convenience/latex.h"
#include <sstream>
namespace Wedge {

namespace internal {



/** @brief Helper class to convert a %GiNaC polynomial into a %CoCoA polynomial */
template<typename Coordinate> class PolynomialVisitor : public visitor, public add::visitor, public Coordinate::visitor, public mul::visitor, public power::visitor, public numeric::visitor {
	CoCoA::ring Qx;
	CoCoA::RingElem result;
	const exvector& symbols;
	const vector<CoCoA::RingElem>& indets;
	 
public:	
	PolynomialVisitor(CoCoA::ring _Qx, const exvector& _symbols,const vector<CoCoA::RingElem>& _indets ) : Qx(_Qx), result(Qx),symbols(_symbols),indets(_indets)
	{
	} 
	
	CoCoA::RingElem RecursiveVisit(ex x)
	{
		CoCoA::RingElem old_result=result;
		x.accept(*this);
		CoCoA::RingElem new_result=result;
		result=old_result;
		return new_result;
	}
	
	void visit(const basic& o)
	{
		 throw WedgeException<std::runtime_error>("Invalid type in polynomial",__FILE__,__LINE__);
	}
	
	void visit(const add& o)
	{
		result=RecursiveVisit(o.op(0));
		for (int i=1;i<o.nops();i++)
			result+=RecursiveVisit(o.op(i));			
	}
	void visit(const mul& o)
	{
		result=RecursiveVisit(o.op(0));
		for (int i=1;i<o.nops();i++)
			result*=RecursiveVisit(o.op(i));			
	}
	void visit(const numeric& o)
	{
		if (o.is_real() && o.is_rational())	// i is considered a rational number
		{
			result=CoCoA::RingElem(Qx,o.numer().to_long())/CoCoA::RingElem(Qx,o.denom().to_long());
		}			
	}
	void visit(const power& o)
	{
		if (!is_a<numeric>(o.op(1))) {
			LOG_ERROR(o);
			throw WedgeException<std::runtime_error>("Non-numeric exponent in polynomial",__FILE__,__LINE__);						
		}
		const numeric& a=ex_to<numeric>(o.op(1));
		if (!a.is_integer()) {
			LOG_ERROR(o);
			throw WedgeException<std::runtime_error>("Non-integer exponent in polynomial",__FILE__,__LINE__);						
		}
		const int as_int=a.to_int();
		if (a<0) {
			LOG_ERROR(o);
			throw WedgeException<std::runtime_error>("Negative exponent in polynomial",__FILE__,__LINE__);						
		}
		result=CoCoA::power(RecursiveVisit(o.op(0)),as_int);
	}
	void visit(const typename internal::VisitorClass<Coordinate>::type& o)
	{			
		exvector::const_iterator i=find(symbols.begin(),symbols.end(),o);			
		if (i==symbols.end()) {
			LOG_ERROR(o);
			throw WedgeException<std::runtime_error>("Unknown symbol in polynomial",__FILE__,__LINE__);
		}
		assert(i-symbols.begin()>=0);
		assert(i-symbols.begin()<indets.size());
		result=indets[i-symbols.begin()];			
	}		
};	
}
/** @brief Default implementation of polynomial algebra algorithms (uses %CoCoA)
 * 
 * @remark This implementation is valid over the field \f$\mathbb{Q}\f$.
 */
struct CocoaPolyAlgorithms{
	struct Message {
		Message() {LOG_MSG("Initializing CoCoA");}
	};	
	struct Initializer {
		Initializer() {static Message message; static CoCoA::GlobalManager cocoaFoundations;}
	};
	
	static CoCoA::SparsePolyRing PolynomialRingOverQ(int indets) {
		vector<CoCoA::symbol> symbols;
		for (int i=0;i<indets;++i) symbols.emplace_back("x",i);
		CoCoA::FractionField Q=CoCoA::RingQQ();
		CoCoA::SparsePolyRing Qx = CoCoA::NewPolyRing(Q,symbols);
		return Qx;
	}

	template<typename Variable, typename Iterator> static exvector IdealReduce(const exvector& symbols, Iterator pol_begin, Iterator pol_end)
	{		
		if (symbols.empty())
		{ 		
			exvector result;
			if (IdealIsOne<Variable>(symbols,pol_begin,pol_end)) result.push_back(1);
			return result;   	
		}
		else {
			CoCoA::SparsePolyRing Qx = PolynomialRingOverQ(symbols.size());
			CoCoA::ideal I(Qx,Ginac2Cocoa<Variable>(Qx,symbols,pol_begin,pol_end));
			vector<CoCoA::RingElem> gens=CoCoA::TidyGens(I);
			return Cocoa2Ginac(symbols, gens.begin(),gens.end());	
		}		
	} 
	
	template<typename Variable, typename Iterator> static bool IdealContains(const exvector& symbols, Iterator pol_begin, Iterator pol_end, ex p)
	{	
		set<ex,ex_is_less> symbols_in_p;			
		GetSymbols<Variable>(symbols_in_p,p);
		
		exvector allsymbols=symbols;
		for (typename set<ex,ex_is_less>::const_iterator i=symbols_in_p.begin();i!=symbols_in_p.end();i++)		
			if (find(symbols.begin(),symbols.end(),*i)==symbols.end())
				allsymbols.push_back(*i);
			
		if (allsymbols.empty())		 		
			 return p.is_zero() || IdealIsOne<Variable>(allsymbols,pol_begin,pol_end);
		else {
			CoCoA::FractionField Q=CoCoA::RingQQ();
			CoCoA::SparsePolyRing Qx = CoCoA::NewPolyRing(Q,CoCoA::NewSymbols(allsymbols.size()));
			CoCoA::ideal I(Qx,Ginac2Cocoa<Variable>(Qx,allsymbols,pol_begin,pol_end));
			CoCoA::RingElem f=Ginac2Cocoa<Variable>(Qx,allsymbols,&p,(&p)+1)[0];
			return CoCoA::IsElem(f,I);	
		}		
	}

	
	template<typename Variable, typename Iterator> static bool RadicalContains(const exvector& symbols, Iterator pol_begin, Iterator pol_end, ex p)
	{	
		ex newpoly=p*Variable()-1;

		set<ex,ex_is_less> allsymbols(symbols.begin(),symbols.end());
		GetSymbols<Variable>(allsymbols,newpoly);
		exvector poly(pol_begin,pol_end);
		poly.push_back(newpoly);

		return IdealIsOne<Variable>(exvector(allsymbols.begin(),allsymbols.end()),poly.begin(),poly.end());

	//	alternatively, we might invoke CoCoA::RadicalMembership, which as of version 0.9933 gives incorrect results.
/*		set<ex,ex_is_less> symbols_in_p;			
		GetSymbols<Variable>(symbols_in_p,p);
		
		exvector allsymbols=symbols;
		for (typename set<ex,ex_is_less>::const_iterator i=symbols_in_p.begin();i!=symbols_in_p.end();i++)		
			if (find(symbols.begin(),symbols.end(),*i)==symbols.end())
				allsymbols.push_back(*i);
			
		if (allsymbols.empty())		 		
			 return p.is_zero() || IdealIsOne<Variable>(allsymbols,pol_begin,pol_end);
		else {
			CoCoA::FractionField Q=CoCoA::RingQQ();
			CoCoA::SparsePolyRing Qx = CoCoA::NewPolyRing(Q,allsymbols.size());
			vector<CoCoA::RingElem> I=Ginac2Cocoa<Variable>(Qx,allsymbols,pol_begin,pol_end);
			CoCoA::RingElem f=Ginac2Cocoa<Variable>(Qx,allsymbols,&p,(&p)+1)[0];			
			return CoCoA::RadicalMembership(I,f);	
		}
*/	
	}


/**  @brief
*/
	template<typename Variable, typename Iterator> static ex ElementModuloIdeal(const exvector& symbols, Iterator pol_begin, Iterator pol_end, ex p)
	{
					
		set<ex,ex_is_less> symbols_in_p;			
		GetSymbols<Variable>(symbols_in_p,p);
		
		exvector allsymbols=symbols;
		for (typename set<ex,ex_is_less>::const_iterator i=symbols_in_p.begin();i!=symbols_in_p.end();i++)		
			if (find(symbols.begin(),symbols.end(),*i)==symbols.end())
				allsymbols.push_back(*i);
			
		if (symbols_in_p.empty())		 		
			 return (p.is_zero() || IdealIsOne<Variable>(allsymbols,pol_begin,pol_end))? 0 : p;
		else {
			CoCoA::SparsePolyRing Qx = PolynomialRingOverQ(allsymbols.size());
			CoCoA::ideal I(Qx,Ginac2Cocoa<Variable>(Qx,allsymbols,pol_begin,pol_end));
			CoCoA::RingElem f=Ginac2Cocoa<Variable>(Qx,allsymbols,&p,(&p)+1)[0];
			CoCoA::RingElem r=f % I;							
			return Cocoa2Ginac(allsymbols,&r, (&r)+1)[0];
		}
	}
	
	
	template<typename Variable, typename Iterator> static bool IdealIsOne(const exvector& symbols, Iterator pol_begin, Iterator pol_end)
	{
		if (symbols.empty())
		{
			 while (pol_begin!=pol_end)
			 	if (!(pol_begin++->is_zero())) return true;
			 return false;
		}
		else {
			CoCoA::SparsePolyRing Qx = PolynomialRingOverQ(symbols.size());
			CoCoA::ideal I(Qx,Ginac2Cocoa<Variable>(Qx,symbols,pol_begin,pol_end));
			return CoCoA::IsOne(I);						
		}
	}

	template<typename Variable, typename Iterator1, typename Iterator2> static ExVector IdealIntersection(Iterator1 I_begin, Iterator1 I_end, Iterator2 J_begin, Iterator2 J_end)
	{
		set<ex,ex_is_less> symbols;
		GetSymbols<Variable>(symbols,I_begin,I_end);
		GetSymbols<Variable>(symbols,J_begin,J_end);
		exvector s(symbols.begin(),symbols.end());

		Iterator1 i;
		for (i=I_begin;i!=I_end && i->is_zero();++i) ;
		Iterator2 j;
		for (j=J_begin;j!=J_end && j->is_zero();++j) ;
		if (i==I_end || j==J_end)	//one of I, J is the zero ideal 
			return ExVector();
		else if (s.empty())	//neither is zero but no indet appears -> both are the full ring
		{
			ExVector v;
			v.push_back(1);
			return v;
		}
		else {				//neither is zero and at least one is not the full ring 
			CoCoA::SparsePolyRing Qx = PolynomialRingOverQ(s.size());
			CoCoA::ideal I(Qx,Ginac2Cocoa<Variable>(Qx,s,I_begin,I_end));
			CoCoA::ideal J(Qx,Ginac2Cocoa<Variable>(Qx,s,J_begin,J_end));
			CoCoA::ideal K=CoCoA::intersect(I,J);
			LOG_DEBUG(I);
			LOG_DEBUG(J);
			LOG_DEBUG(K);
			//vector<CoCoA::RingElem> gens=CoCoA::TidyGens(K); 			//triggers the bug 
			vector<CoCoA::RingElem> gens=CoCoA::gens(K);
			return Cocoa2Ginac(s, gens.begin(),gens.end());	
/*			CoCoA::PolyList I=Ginac2Cocoa<Variable>(Qx,s,I_begin,I_end);
			CoCoA::PolyList J=Ginac2Cocoa<Variable>(Qx,s,J_begin,J_end);
			CoCoA::PolyList IcapJ;
			CoCoA::ComputeIntersection(IcapJ,I,J);
			return Cocoa2Ginac(s,IcapJ.begin(),IcapJ.end());*/
		}
	}

//	template<typename Variable, typename Iterator> static bool IdealIsZero(const vector<Variable>& symbols, Iterator pol_begin, Iterator pol_end)
//	{
//		if (symbols.empty())
//		{
//			 while (pol_begin!=pol_end)
//			 	if (!pol_begin++->is_zero()) return false;
//			 return true;
//		}
//		else {
//			CoCoA::FractionField Q=CoCoA::RingQQ();
//			CoCoA::ring Qx = CoCoA::NewPolyRing(Q,symbols.size());
//			CoCoA::ideal I(Qx,Ginac2Cocoa(Qx,symbols,pol_begin,pol_end));
//			return CoCoA::IsZero(I);			
//		}
//	}

private:
	template<typename Variable, typename Iterator> static vector<CoCoA::RingElem> Ginac2Cocoa(const CoCoA::SparsePolyRing& Qx, const exvector& symbols,Iterator pol_begin, Iterator pol_end)
	{
		try {			
			vector<CoCoA::RingElem> x=CoCoA::indets(Qx);
			assert(x.size()==symbols.size());			
			internal::PolynomialVisitor<Variable> v(Qx,symbols,x);		
			vector<CoCoA::RingElem> polynomials; 
			polynomials.reserve(pol_end-pol_begin);
			while (pol_begin!=pol_end)
			{
				CoCoA::RingElem p=v.RecursiveVisit(*pol_begin++);
				if (p!=0) polynomials.push_back(p);	//CoCoA hangs horribly if you pass zero polynomials...
				assert(CoCoA::owner(p)==Qx);	
			}		
			return polynomials;
		}
		catch (CoCoA::ErrorInfo& exception) {
			stringstream ss;
			ss<<"CoCoA error in function Ginac2Cocoa: ";
			ss<<exception;
			throw WedgeException<std::runtime_error>(ss.str(),__FILE__,__LINE__);
		}
	}

	template<typename Iterator> static ExVector Cocoa2Ginac(const exvector& symbols,Iterator pol_begin, Iterator pol_end)
	{
		ExVector result;
		while (pol_begin!=pol_end)
		{	
			stringstream str;		
			str<<*pol_begin++;
			LOG_DEBUG(str.str());
			LOG_DEBUG(symbols);	
			result.push_back(ParseCocoaExpression(symbols,str.str().c_str()));
		}
		LOG_DEBUG(result);
		return result;
	}
};

/** @brief Alternative (experimental) implementation of polynomial algebra algorithms (uses %CoCoA)
 * 
 * @remark This implementation is valid over the extension of \f$\mathbb{Q}\f$ that contains \f$x_i^{1/n_i}\f$ for rational numbers x_i>0 and coprime integers n_i. 
 * @todo Extend to \f$\mathbb{R}\f$ by using ex::to_rational 
 */
struct CocoaPolyAlgorithms_R : public CocoaPolyAlgorithms
{
private:
/** @brief Convert ideals of polynomials over \f$\mathbb{R}\f$ to polynomials over \f$\mathbb{Q}\f$ by adding extra variables and polynomials
 * 
 * Example: \f$(x^2-\sqrt2)\f$ is transformed into \f$(x^2-y, y^2-2)\f$.
 */
	template<class Variable> class PolynomialsOverQ  {
		/** @brief used internally during initialization */
		struct Roots : public visitor, public power::visitor { 
			struct IntegralRoot {			//represents the roots n^{1/B} and n^{-1/B}; need n>0.
				numeric n;			//an integer
				ex positiveroot, negativeroot;	//really a Variable object
				IntegralRoot(const numeric& n) 
				{
					positiveroot=Variable();
					negativeroot=Variable();
					this->n=n;
				}
			};
			numeric B;	//an integer; the common exponent is 1/B.
			list<IntegralRoot> integralroots;		//this is a list of all the integral roots objects
			map<power,ex,basic_is_less> substitutions;	//map a power object to a product of Variables, each representing an integral root

			Roots() {B=1;}

			ex AddIntegralRoot(numeric n,bool positive) {
				if (n==1) return 1;
				for (typename list<IntegralRoot>::const_iterator i=integralroots.begin();i!=integralroots.end();++i)
					if (i->n==n)
						return positive ? i->positiveroot : i->negativeroot;
				IntegralRoot r(n);
				integralroots.push_back(r);
				return positive ? r.positiveroot : r.negativeroot;
			}

			void AddCommonFactors() {
				//rewrite all powers as x^{a/b}=x^{1/B}^{aB/b}; x^{1/B} will be added as a variable		
				for (typename map<power,ex,basic_is_less>::const_iterator i=substitutions.begin();i!=substitutions.end();++i)
				{				
					LOG_DEBUG(i->first);
					assert(is_a<numeric>(i->first.op(1)));
					B=lcm(B,ex_to<numeric>(i->first.op(1)).denom());
				}
				for (typename map<power,ex,basic_is_less>::iterator i=substitutions.begin();i!=substitutions.end();++i)
				{
					assert(is_a<numeric>(i->first.op(0)));
					const numeric& x=ex_to<numeric>(i->first.op(0));
					ex xnum=AddIntegralRoot(x.numer(),true);
					ex xden=AddIntegralRoot(x.denom(),false);
					assert(is_a<numeric>(i->first.op(1)*B));
					assert(ex_to<numeric>(i->first.op(1)*B).is_integer());
					i->second=pow(xnum,i->first.op(1)*B)*pow(xden,i->first.op(1)*B);
				}
				//if n^{1/B} and m^{1/B} appear as roots with (n,m)=d, then we add the root  d^{1/B}
				for (typename list<IntegralRoot>::iterator i=integralroots.begin();i!=integralroots.end();++i)
				{
					typename list<IntegralRoot>::iterator j=i;
					while (++j!=integralroots.end()) {
						numeric d=gcd(i->n,j->n);
						if (d!=1) {
							i->n/=d;
							j->n/=d;
							lst subs;
							subs.append(i->positiveroot==i->positiveroot*AddIntegralRoot(d,true));
							subs.append(i->negativeroot==i->negativeroot*AddIntegralRoot(d,false));
							subs.append(j->positiveroot==j->positiveroot*AddIntegralRoot(d,true));
							subs.append(j->negativeroot==j->negativeroot*AddIntegralRoot(d,false));
							for (map<power,ex,basic_is_less>::iterator i=substitutions.begin();i!=substitutions.end();++i)
								i->second=i->second.subs(subs);
						}
					}
				}
				//debugging stuff
/*
				lst subs;		
				for (typename list<IntegralRoot>::iterator i=integralroots.begin();i!=integralroots.end();++i)			
				{
					subs.append(i->positiveroot==symbol(ToString(i->n)+"^(1/"+ToString(B)+")"));
					subs.append(i->negativeroot==symbol(ToString(i->n)+"^(1/"+ToString(B)+")"));
					LOG_INFO(i->n);
				}
				for (map<power,ex,basic_is_less>::iterator i=substitutions.begin();i!=substitutions.end();++i)
				{
					LOG_INFO(i->first);
					LOG_INFO(i->second.subs(subs));
				}
*/
			}

			void visit(const power& p)
			{
				if (is_a<numeric>(p.op(0)) && is_a<numeric>(p.op(1)))
					substitutions[p]=0;
			}			

		};

		lst back_subs;	//substitutions to give back ordinary polynomials with no extra variables
	public:
		exvector symbols; ///<Variables appearing in the polynomials, including the "extra" variables introduced to convert 
		exvector polynomials;	///<Polynomials over \f$\mathbb{Q}\f$
		
		PolynomialsOverQ(const exvector& _symbols ) : symbols(_symbols) 
		{										
		}
	
		/** @brief Add the polynomials
		 * @param [pol_begin,pol_end) A range of polynomials over \f$\mathbb{R}\f$
		 * 
		 * The polynomial corresponding to the required extension of \f$\mathbb{Q}\f$ are added first, then the specified range
		 */
		template<typename Iterator> void AddPolynomials(Iterator pol_begin, Iterator pol_end)
		{		
			Roots roots;	
			LOG_DEBUG(polynomials);
			//rewrite the polynomials already in the container in terms of roots and reform the "roots" helper object.
			for (exvector::iterator i=polynomials.begin();i!=polynomials.end();++i)	
			{			
				*i=i->subs(back_subs);
				i->traverse(roots);
			}
			for (Iterator i=pol_begin; i!=pol_end;++i)
			{
				LOG_DEBUG(*i);
				i->traverse(roots);
			}
			roots.AddCommonFactors();
	
			back_subs.remove_all();
			//add the extra symbols and polynomials
			lst remove_trivial_subs;	//a list of substitutions to remove roots of the form 1^{1/B}
			for (typename list<typename Roots::IntegralRoot>::iterator i=roots.integralroots.begin();i!=roots.integralroots.end();++i)
			{
				if (i->n==1) {
					remove_trivial_subs.append(i->positiveroot==1);
					remove_trivial_subs.append(i->negativeroot==1);
				}
				else {
					symbols.push_back(i->positiveroot);
					symbols.push_back(i->negativeroot);
					polynomials.push_back(pow(i->positiveroot,roots.B)-i->n);
					polynomials.push_back(i->positiveroot*i->negativeroot-1);
					//polynomials.push_back(pow(i->negativeroot,roots.B)-ex(1)/i->n);
					typename list<typename Roots::IntegralRoot>::iterator j=i;
					while (++j!=roots.integralroots.end())
						assert (i->n!=j->n);
					back_subs.append(i->positiveroot==pow(i->n,ex(1)/roots.B));
					back_subs.append(i->negativeroot==pow(i->n,-ex(1)/roots.B));
				}
			}
			//replace the roots with the extra symbols
			lst subs;
			for (typename map<power,ex,basic_is_less>::const_iterator i=roots.substitutions.begin();i!=roots.substitutions.end();++i)
			{
				ex as_poly=i->second.subs(remove_trivial_subs);
				subs.append(i->first==as_poly);
			}
			//rewrite  the polynomials already in the container in terms of the new extra symbols
			for (exvector::iterator i=polynomials.begin();i!=polynomials.end();++i)	
				*i=i->subs(subs);
			//add the new polynomials
			for (Iterator i=pol_begin; i!=pol_end;++i)
				polynomials.push_back(i->subs(subs));
			LOG_DEBUG(polynomials);
		}
		
		/** @brief Remove the extra variables, and convert back to ordinary polynomials over \f$\mathbb{R}\f$
		 */
		exvector Convert(const exvector& polynomials) const
		{
			exvector result;
			for (exvector::const_iterator i=polynomials.begin(); i!=polynomials.end();++i)
				result.push_back(i->subs(back_subs));
			return result;			
		}
		
		/** @overload
		 */
		ex Convert(ex pol) const
		{
			return pol.subs(back_subs);
		}
	};
public:
	template<typename Variable, typename Iterator> static exvector IdealReduce(const exvector& symbols, Iterator pol_begin, Iterator pol_end)
	{
		PolynomialsOverQ<Variable> r(symbols);
		r.AddPolynomials(pol_begin,pol_end);		
		exvector overQ=CocoaPolyAlgorithms::IdealReduce<Variable>(r.symbols,r.polynomials.begin(),r.polynomials.end());
		return r.Convert(overQ);		
	}

	template<typename Variable, typename Iterator> static bool IdealContains(const exvector& symbols, Iterator pol_begin, Iterator pol_end, ex p)
	{
		PolynomialsOverQ<Variable> r(symbols);
		r.AddPolynomials(pol_begin,pol_end);
		r.AddPolynomials(&p,(&p)+1);
		return CocoaPolyAlgorithms::IdealContains<Variable>(r.symbols,r.polynomials.begin(),r.polynomials.end()-1,r.polynomials.back());
	}

	template<typename Variable, typename Iterator> static bool RadicalContains(const exvector& symbols, Iterator pol_begin, Iterator pol_end, ex p)
	{
		PolynomialsOverQ<Variable> r(symbols);
		r.AddPolynomials(pol_begin,pol_end);
		r.AddPolynomials(&p,(&p)+1);
		return CocoaPolyAlgorithms::RadicalContains<Variable>(r.symbols,r.polynomials.begin(),r.polynomials.end()-1,r.polynomials.back());
	}
	
	template<typename Variable, typename Iterator> static ex ElementModuloIdeal(const exvector& symbols, Iterator pol_begin, Iterator pol_end, ex p)
	{
		PolynomialsOverQ<Variable> r(symbols);
		r.AddPolynomials(pol_begin,pol_end);
		r.AddPolynomials(&p,(&p)+1);	
		ex mod=CocoaPolyAlgorithms::ElementModuloIdeal<Variable>(r.symbols,r.polynomials.begin(),r.polynomials.end()-1,r.polynomials.back());
		return r.Convert(mod);
	}			
	
	template<typename Variable, typename Iterator> static bool IdealIsOne(const exvector& symbols, Iterator pol_begin, Iterator pol_end)
	{
		PolynomialsOverQ<Variable> r(symbols);
		r.AddPolynomials(pol_begin,pol_end);
		return CocoaPolyAlgorithms::IdealIsOne<Variable>(r.symbols,r.polynomials.begin(),r.polynomials.end());		
	}
};

typedef CocoaPolyAlgorithms DefaultPolyAlgorithms;

} /** @} */
#endif /*COCOAPOLYALG_H_*/


