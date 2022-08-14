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
#include "wedge/base/wedgebase.h"
#include "wedge/base/logging.h"
#include "latex.h"
#include "printcontext.h"

namespace Wedge {

std::ostream& operator<<(std::ostream& out, const ExVector& V)
{
	return out<<Range(V.begin(),V.end());
}


ostream& operator<<(ostream& os, const LatexContainer& c)
{
	c.flush(os);
	return os;
}

LatexContainer::~LatexContainer() {
	string k=s.str();
	if (!k.empty())
		LOG_WARN("The following LatexContainer object was not written to output: "<<k);
}

ostream& operator<<(ostream& os, Equation x)
{
	print_context* os_print_context=get_print_context(os);
	if (dynamic_cast<const GiNaC::print_latex*>(os_print_context)==NULL) {
		os<<latex;
		os_print_context=get_print_context(os);
	}
	ex e=x.e;
	if (!is_a<add>(e))
	{
		stringstream latex_stream;
		set_print_context(latex_stream, *os_print_context);
		latex_stream<<"\\[";
		if (!x.text.empty()) latex_stream<<x.text<<" = ";
		latex_stream<<e<<"\\]"<<endl;
		return os<<latex_stream.str();
	}
	int size=0;	//total size of lines already processed
	string output;
	bool multline=false;
	stringstream endline;
	endline<<"\\\\"<<endl;
	for (int i=0;i<e.nops();++i)
	{
		stringstream s;
		set_print_context(s, *os_print_context);
		s<<e.op(i)<<std::flush;
		assert (s.str().size()>0);
		size+=s.str().size();
		if (size>100) {
			output+=endline.str();
			size=s.str().size();
			multline=true;
		}
		if (i>0 && s.str()[0]!='-') output+="+";
		output+=s.str();
	}
	//handle multiple lines
	if (multline) {
		os<<"\\begin{multline*}"<<endl;
		if (!x.text.empty()) os<<x.text<<" = ";
		os<<output<<endl<<"\\end{multline*}"<<endl;
	}
	else {
		os<<"\\[";
		if (!x.text.empty()) os<<x.text<<" = ";
		os<<output<<"\\]"<<endl;
	}
	return os;
}


void event_callback(std::ios::event ev, ios_base& ios, int index)
{
	const Label* label=Label::FromStream(ios);
	switch (ev) {
		case ios_base::copyfmt_event:
		if (label)
			label->Link(ios);
		break;
		case ios_base::erase_event:
	    	if (label) {delete label; Label::Unlink(ios);}
		default:
			break;
	}
}

const Label* Label::FromStream(std::ios_base& ios)
{
	return static_cast<const Label*>(ios.pword(ios_index));
}
void Label::Unlink(std::ios_base& ios)
{
	ios.pword(ios_index)=0;
}
void Label::Link(std::ios_base& ios) const
{
	ios.pword(ios_index)=const_cast<Label*>(this);
	ios.register_callback(event_callback,0);
}

const int Label::ios_index=ios_base::xalloc();

ostream& operator<< (ostream& os, const Label& label)
{
	(new Label(label))->Link(os);
	return os;
}

/** @brief Manipulator for use at the beginning of a latex session
*/
ostream& latex_begin(ostream& os)
{
	os<<latex<<"\\documentclass{article}"<<endl;
	os<<"\\usepackage{amsmath}"<<endl;
	return os<<"\\begin{document}"<<endl;
}

/** @brief Manipulator for use at the end of a latex session
*/
ostream& latex_end(ostream& os)
{
	return os<<"\\end{document}"<<endl;
}

}
