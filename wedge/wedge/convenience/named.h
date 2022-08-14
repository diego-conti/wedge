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
#ifndef NAMED_H
#define NAMED_H
#include <boost/lexical_cast.hpp>
#include <string>
#include <vector>

using namespace std;
using boost::lexical_cast;
namespace Wedge {

class NameAndIndex; 

/** This class provides a method to define the name of an algebraic object on construction.

Whenever V is a named class, V can be instantiated like V(Name("name")), or V(Name("name","texname")). Indices can be added using round brackets, as in V(Name("A")(i,j)). Superscripts (plus, minus, plusminus, prime, second) can be added by invoking the corresponding member function, as in V(Name("A").prime()); the same applies to the modifiers tilde and hat.

Single-character names are accessible as members of the global object N (@sa PredefinedNames).

@remark Objects of type Name are only meant to be used as arguments of a constructor.

@examples
- V(N.x)	-> x,x
- V(N.psi)	-> psi, \psi
- V(N.x(1))	-> x1, x_1
- V(N.x(10))	-> x10, x_{10}
- V(N.x(1,2))	-> x{1,2}, x_1^2
- V(N.x(1).tilde()) -> x~1, \tilde{x}_1
- V(N.x(1,2).plus()) -> x+{1,2}, x^+_{1,2}
- V(N.x(1,2,3)) -> x{1,2,3}, x_{1,2,3}
- V(N.x(1)(2)(3)(4)) -> x{1,2,3,4}, x_{1,2,3,4} 
*/



class Name {
	Name operator=(const Name&);
	Name(const Name&);	//disabled because it could convert a NameAndIndex to a Name
protected:
	string plain_, _tex;
	bool _hasSuperscript;
	Name(const string& plain, const string& tex, bool hasSuperscript) : plain_(plain), _tex(tex) {_hasSuperscript=hasSuperscript;}
public:
	Name(const string& plain, const string& tex) : plain_(plain), _tex(tex) {_hasSuperscript=false;}
	Name(const string& plain) : plain_(plain), _tex(plain) {_hasSuperscript=false;}

	NameAndIndex operator()  (int n) const;
	NameAndIndex operator()  (int n, int m) const;
	NameAndIndex operator()  (int n, int m, int k) const;
	Name Prime() const {
		return Name(plain_+"'",_tex+"'",true);
	}
	Name Second() const {
		return Name(plain_+"''",_tex+"''",true);
	}
	Name Plus() const {
		return Name(plain_+"+",_tex+"^+",true);
	}
	Name PlusMinus() const {
		return Name(plain_+"+-",_tex+"^\\pm",true);
	}
	Name Minus() const {
		return Name(plain_+"-",_tex+"^-",true);
	}

	Name Tilde() const {
		return Name(plain_+"~","\\tilde{"+_tex+"}",_hasSuperscript);
	}
	Name Hat() const {
		return Name(plain_+"^","\\hat{"+_tex+"}",_hasSuperscript);
	}

	virtual NameAndIndex PlaceHolder() const;	//return the same object as a NameAndIndex with no additional indices, to which an index can be added.
	virtual string plain() const {return plain_;}
	virtual string tex() const {return _tex;}
	virtual ~Name() {}
};

/** Subclass of Name representing a name with no indices. In addition to Name, it has a copy constructor. */
class NameWithNoIndices : public Name {
public:
	NameWithNoIndices(const NameWithNoIndices& o) : Name(o.plain_,o._tex,o._hasSuperscript) {}
	NameWithNoIndices(const string& plain, const string& tex) : Name(plain,tex) {}
	NameWithNoIndices(const string& plain) : Name(plain) {}
};

/** Subclass of Name representing a name and zero or more indices. NameAndIndex'es with zero indices can be used as placeholders.
*/
class NameAndIndex : public Name {
	friend class Name;
	vector<int> indices_;
	static string toString(int n) {
		string s=lexical_cast<string>(n);
		if (s.size()<2) return s; else return "{"+s+"}";
	}
	NameAndIndex(const string& plain, const string& tex, bool hasSuperscript) : Name(plain,tex,hasSuperscript) {}
	NameAndIndex(const string& plain, const string& tex, bool hasSuperscript, const vector<int>& indices) : Name(plain,tex,hasSuperscript),indices_(indices) {}
public:
	NameAndIndex();
	NameAndIndex(const string& plain, const string& tex, const vector<int>& indices) : Name(plain,tex),indices_(indices) {}
	NameAndIndex(const NameAndIndex& o) : Name(o.plain_,o._tex,o._hasSuperscript), indices_(o.indices_) {}
	NameAndIndex& operator=(const NameAndIndex& o) {plain_=o.plain_; _tex=o._tex; _hasSuperscript=o._hasSuperscript; indices_=o.indices_; return *this;}

	NameAndIndex operator() (int n) const
	{
		NameAndIndex result(plain_,_tex,_hasSuperscript,indices_);
		result.indices_.push_back(n);
		return result;
	}
	NameAndIndex operator()  (int n, int m) const
	{
		NameAndIndex result(plain_,_tex,_hasSuperscript,indices_);
		result.indices_.push_back(n);
		result.indices_.push_back(m);
		return result;
	}
	NameAndIndex operator()  (int n, int m, int k) const
	{
		NameAndIndex result(plain_,_tex,_hasSuperscript,indices_);
		result.indices_.push_back(n);
		result.indices_.push_back(m);
		result.indices_.push_back(k);
		return result;
	}
	NameAndIndex Plus() const {
		return NameAndIndex(plain_+"+",_tex+"^+",true,indices_);
	}
	NameAndIndex PlusMinus() const {
		return NameAndIndex(plain_+"+-",_tex+"^\\pm",true,indices_);
	}
	NameAndIndex Minus() const {
		return NameAndIndex(plain_+"-",_tex+"^-",true,indices_);
	}
	NameAndIndex Prime() const {
		return NameAndIndex(plain_+"'",_tex+"'",true,indices_);
	}
	NameAndIndex Second() const {
		return NameAndIndex(plain_+"''",_tex+"''",true,indices_);
	}
	NameAndIndex Tilde() const {
		return NameAndIndex(plain_+"~","\\tilde{"+_tex+"}",_hasSuperscript,indices_);
	}
	NameAndIndex Hat() const {
		return NameAndIndex(plain_+"^","\\hat{"+_tex+"}",_hasSuperscript,indices_);
	}
	NameAndIndex PlaceHolder() const  {return *this;}

	string plain() const {
		switch (indices_.size()) {
			case 0:
				return plain_;
			case 1:
				return plain_+lexical_cast<string>(indices_[0]);
			default:
				string result=plain_+lexical_cast<string>(indices_[0]);
				for (int k=1;k<indices_.size();++k) result+="$"+lexical_cast<string>(indices_[k]);
				return result;
		}
	}
	string tex() const {
		switch (indices_.size()) {
			case 0:
				return _tex;
			case 1:
				return _tex+"_"+toString(indices_[0]);
			case 2:
				if (!_hasSuperscript)
					return _tex+"_"+toString(indices_[0])+"^"+toString(indices_[1]);
				[[fallthrough]];
			default:
				string result=_tex+"_{"+lexical_cast<string>(indices_[0]);
				for (int k=1;k<indices_.size();++k) result+=","+lexical_cast<string>(indices_[k]);
				return result+"}";
		}
	}

};

/** Represents a range of variable names, e.g. a_1, ...,a_20. */
struct NameRange {
	const Name& name_;
	int begin_,end_;
	NameRange(const Name& name, int begin,int end) : name_(name) {begin_=begin, end_=end;}

	class const_iterator {
		const Name& name_;
		int i_;
	public:
		const_iterator(const Name& name, int i) : name_(name), i_(i) {}
		NameAndIndex operator*() const {return name_(i_);}
		bool operator==(const const_iterator& o) {return name_.plain()==o.name_.plain() && i_==o.i_;}
		bool operator!=(const const_iterator& o) {return name_.plain()!=o.name_.plain() || i_!=o.i_;}
		const_iterator& operator++() {++i_; return *this;}
	};
	const_iterator begin() const {return const_iterator(name_,begin_);}
	const_iterator end() const {return const_iterator(name_,end_);}
};


/** A single-instance class that contains every single letter that can be used as a name (English alphabet and Greek letters with a corresponding LaTeX command)
 */
extern struct PredefinedNames {
	Name alpha,beta,gamma,Gamma, delta,Delta,epsilon,zeta,eta,theta,Theta,kappa,lambda,Lambda,mu,nu,xi,Xi,rho,pi,Pi,sigma,Sigma,tau,upsilon,Upsilon,phi,Phi,chi,psi,Psi,omega,Omega;
	Name a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z;
	PredefinedNames();
} N;

}


#endif
