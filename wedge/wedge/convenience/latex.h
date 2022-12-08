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
#ifndef LATEX_H_
#define LATEX_H_

/** @ingroup WedgeBase */ 

/** @{ 
 * @file latex.h
 * @brief Latex output utilities
 */

#include "named.h"

using namespace std;
using namespace GiNaC;

namespace Wedge {

/** @brief Detect whether a stream expects \f$\text{\TeX}\f$ output
 * @param out An output stream
 * @return True if the stream expects \f$\text{\TeX}\f$ output
 */
template<class charT, class traits> bool IsStreamTeX(std::basic_ostream<charT,traits>& os)
{
	print_context* p=get_print_context(os);
	return dynamic_cast<print_latex*>(p)!=NULL;
}



namespace internal {
/** @brief Used internally to represent ranges in a STL container (for output)
 */
template<typename Iterator> struct Range
{
	Iterator begin,end;
};

}


class Label {
public:
	const NameAndIndex name_;
	Label(const Name& name) : name_(name.PlaceHolder()) {}
	static const Label* FromStream(std::ios_base& ios);
	static void Unlink(std::ios_base& ios);
	void Link(std::ios_base& ios) const;
private:
	const static int ios_index;
};


ostream& operator<< (ostream& os, const Label& label);

/** @brief Convert an STL range into an object that can be output to a stream
 * @param [b,e) a range
 *
 * Example: cout<<Range(begin,end); LOG_DEBUG(Range(begin,end));
 */
template<typename Iterator> typename internal::Range<Iterator> Range(Iterator b,Iterator e) {
	internal::Range<Iterator> result={b,e};
	return result;
};

/** @brief Overloaded output operator
 */
template<typename Iterator> std::ostream& operator<<(std::ostream& out, internal::Range<Iterator> range)
{
	const Label* label=Label::FromStream(out);

	if (range.begin==range.end) out<<"No elements"<<endl;
	else if (IsStreamTeX(out))
	{
		out<<"\\begin{gather*}"<<endl;
		int i=0;
		if (label) out<<label->name_(++i).tex()<<"=";
		out<<*range.begin;
		while (++range.begin!=range.end) {
			out<<"\\\\"<<endl;
			if (label) out<<label->name_(++i).tex()<<"=";
			out<<*range.begin;
		}
		out<<endl<<"\\end{gather*}"<<endl;
		if (label) {delete label; Label::Unlink(out);}
	}
	else {
		out<<"{{"<<endl;
		int i=0;
		while (range.begin!=range.end) {
			if (label) out<<label->name_(++i)<<"=";
			out<<*range.begin++<<endl;
		}
		out<<"}}"<<endl;
	}
	return out;
}

/** @brief Overloaded output operator
 */
template<typename T> std::ostream& operator<<(std::ostream& out, const std::vector<T>& V)
{
	return out<<Range(V.begin(),V.end());
}

std::ostream& operator<<(std::ostream& out, const ExVector& V);

/** @brief Overloaded output operator
 */
template<typename T> std::ostream& operator<<(std::ostream& out, const list<T>& V)
{
	return out<<Range(V.begin(),V.end());
}

/** @brief Overloaded output operator
 */
template<class Key, class Compare, class Alloc> std::ostream& operator<<(std::ostream& out, const set<Key,Compare,Alloc>& V)
{
	return out<<Range(V.begin(),V.end());
}

/** @brief Overloaded output operator
 */
template<typename Key, typename Data, typename Compare> std::ostream& operator<<(std::ostream& out, const map<Key,Data,Compare>& V)
{
	if (V.empty()) return out<<"No elements"<<endl;
	else if (IsStreamTeX(out))
	{
		out<<"\\begin{gather*}"<<endl;
		for (typename map<Key,Data,Compare>::const_iterator i=V.begin();i!=V.end();i++)
			out<<i->first<<" \\to "<<i->second<<"\\\\"<<endl;
		out<<"\\end{gather*}"<<endl;
	}
	else {
		out<<"map:"<<endl;
		for (typename map<Key,Data,Compare>::const_iterator i=V.begin();i!=V.end();i++)
			out<<i->first<<" -> "<<i->second<<endl;
		out<<endl;
	}
	return out;
}




/** @brief Helper class to output equations in latex
* 
* In latex output, equations are output as ordinary expressions except that they are enclosed in a multline environment (if needed) and line breaks are inserted
* @ex cout<<Equation(2*x);
*/
struct Equation {
	ex e;
	string text;
	Equation(ex rhs) : e(rhs) {}
	Equation(const string& txt, ex rhs) : e(rhs), text(txt) {}
	Equation(const Name& name, ex rhs) : e(rhs), text(name.tex()) {}
};

ostream& gather(ostream& os);



/** @brief Abstract base class to output lists of equations in latex
*/
class LatexContainer {
protected:
	mutable stringstream s;
	virtual void e() const=0;
public:
	LatexContainer() {s<<latex;}
	LatexContainer(const LatexContainer& c) {
		s<<latex;
		s.str(c.s.str());
		c.s.str("");
	}
	template<typename Iterator> void range(Iterator begin, Iterator end)  {		
		while (begin!=end) item(*begin++);
	}
	virtual void item(ex x)=0;
	virtual void item(exmap::value_type x)=0;
	void flush(ostream& os) const {
		if (!s.str().empty()) {
			e();
			os<<s.str();
			s.str("");
		}
		else os<<"{}"<<endl;
	}
	virtual ~LatexContainer();
};

ostream& operator<<(ostream& os, const LatexContainer& c);


/** @brief Helper class to output lists of equations in latex
 *
 * @ex cout<<Gather(lst(1,2));
 * @ex Gather x; x.item(1); x.item(2); x.flush(cout);
*/

class Gather : public LatexContainer {
	void b() {
		s<<"\\begin{gather*}"<<endl;
	}
	void e() const {
		s<<"\\end{gather*}"<<endl;
	}
public:
	template<typename Iterator> Gather(Iterator begin, Iterator end) {
		range(begin,end);
	}
	Gather(const lst& l) {
		range(l.begin(),l.end());
	} 
	Gather(const exvector& l) {
		range(l.begin(),l.end());
	} 
	Gather() {};
	void item(ex x) {
		if (s.str().empty()) b();
		s<<x<<"\\\\"<<endl;
	}
	void item(exmap::value_type x) {
		if (s.str().empty()) b();
		s<<x.first<<"\\to"<<x.second<<"\\\\"<<endl;
	}	
};


/** @brief Helper class to output lists of equations in latex
 *
 * @sa Gather
*/
class Enumerate : public LatexContainer {
	void b() {
		s<<"\\begin{enumerate*}"<<endl;
	}
	void e() {
		s<<"\\end{enumerate*}"<<endl;
	}
public:
	template<typename Iterator> Enumerate(Iterator begin, Iterator end) {
		range(begin,end);
	}
	Enumerate(lst l) {
		range(l.begin(),l.end());
	} 
	Enumerate(exvector l) {
		range(l.begin(),l.end());
	} 
	Enumerate() {};
	void item(ex x) {
		if (s.str().empty()) b();
		s<<"\\item"<<x<<"\\\\"<<endl;
	}
	void item(exmap::value_type x) {
		if (s.str().empty()) b();
		s<<x.first<<"\\to"<<x.second<<"\\\\"<<endl;
	}	
};

/** @brief Helper class to output lists of equations in latex
 *
 * @sa Gather
*/
class Itemize : public LatexContainer {
	void b() {
		s<<"\\begin{itemize*}"<<endl;
	}
	void e() {
		s<<"\\end{itemize*}"<<endl;
	}
public:
	template<typename Iterator> Itemize(Iterator begin, Iterator end) {
		range(begin,end);
	}
	Itemize(lst l) {
		range(l.begin(),l.end());
	} 
	Itemize(exvector l) {
		range(l.begin(),l.end());
	} 
	Itemize() {};
	void item(ex x) {
		if (s.str().empty()) b();
		s<<"\\item"<<x<<"\\\\"<<endl;
	}
	void item(exmap::value_type x) {
		if (s.str().empty()) b();
		s<<x.first<<"\\to"<<x.second<<"\\\\"<<endl;
	}	
};

/** @brief Helper class to output tables in latex
*/
class Array : public LatexContainer {
	void e() {
		hline();
		s<<"\\end{array}\\]"<<endl;
	}
	int cols, curcol;
public:
/** @brief Create a table with a specified format description
 * @param descr A string of the form {lrc}
 */
	Array(string descr) {
		if (descr.empty() || *descr.begin()!='{' || *descr.rbegin()!='}') throw InvalidArgument(__FILE__,__LINE__);
		s<<"\\[\\begin{array}"<<descr<<endl;
		hline();
		string columnindicator("lcr");
		cols=0;
		for (string::const_iterator i=descr.begin();i!=descr.end();++i)
			if (columnindicator.find(*i)!=string::npos) ++cols;
		curcol=0;
	}
/** @brief Create a table with a specified number of columns
 * @param cols Number of columns
 */
	Array(int cols) {
		this->cols=cols;
		s<<"\\[\\begin{array}{|";
		while (cols-->0) s<<"l|";
		s<<"}"<<endl;
		hline();
		curcol=0;
	}
/** @brief Output a horizontal line
*/
	void hline() {s<<"\\hline"<<endl; curcol=0;}
/** @brief Overloaded shift operator for using Array as a manipulator
 *
 * This function takes care of inserting line breaks and ampersands.
*/
	template<typename T> Array& operator<<(T t) {
		s<<t;
		if (++curcol==cols) {curcol=0; s<<"\\\\"<<endl;}
		else s<<"&";
		return *this;
	}
/** @brief Add a row to the table
*/
	void item(ex x1) {s<<x1<<"\\\\"<<endl;}
	void item(ex x1, ex x2) {s<<x1<<"&"<<x2<<"\\\\"<<endl;} ///< @overload
	void item(ex x1, ex x2, ex x3) {s<<x1<<"&"<<x2<<"&"<<x3<<"\\\\"<<endl;} ///< @overload
	void item(ex x1, ex x2, ex x3, ex x4) {s<<x1<<"&"<<x2<<"&"<<x3<<"&"<<x4<<"\\\\"<<endl;} ///< @overload
	void item(ex x1, ex x2, ex x3, ex x4, ex x5) {s<<x1<<"&"<<x2<<"&"<<x3<<"&"<<x4<<"&"<<x5<<"\\\\"<<endl;} ///< @overload
	void item(ex x1, ex x2, ex x3, ex x4, ex x5, ex x6) {s<<x1<<"&"<<x2<<"&"<<x3<<"&"<<x4<<"&"<<x5<<"&"<<x6<<"\\\\"<<endl;} ///< @overload
	void item(exmap::value_type x) {item(x.first,x.second);}
};


ostream& operator<<(ostream& os, Equation x);

/** @brief Manipulator for use at the beginning of a latex session
*/
ostream& latex_begin(ostream& os);

/** @brief Manipulator for use at the end of a latex session
*/
ostream& latex_end(ostream& os);

}

#endif
