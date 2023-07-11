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
#ifndef EXTRAS_OMITFUNCTIONARGUMENT_H_
#define EXTRAS_OMITFUNCTIONARGUMENT_H_


#include <ginac/ginac.h>

using namespace GiNaC;

namespace Wedge {

struct OmitArgument {
	ex argument;
	OmitArgument(ex arg) {argument=arg;}
};

class print_omit_function_argument : public print_latex
{
     GINAC_DECLARE_PRINT_CONTEXT(print_omit_function_argument, print_latex)
	struct Init {
		Init();
	};
	static Init init;
public:

	ex argument;	//only omit this argument.

	print_omit_function_argument (std::ostream & os, OmitArgument omit, unsigned opt = 0): print_latex(os, opt) {argument=omit.argument;}
	print_omit_function_argument (std::ostream & os, unsigned opt = 0): print_latex(os, opt) {}

	virtual void initialize_from(const print_context& o) {
		argument=static_cast<const print_omit_function_argument&>(o).argument;
	}

	bool ShouldOmit(const function& f) const
	{
		return f.nops()==1 && f.op(0)==argument && !is_order_function(f);
	}
	bool ShouldOmit(ex f) const
	{
		return is_a<function>(f) && ShouldOmit(ex_to<function>(f));
	}

};


std::ostream& operator<<(std::ostream& os, OmitArgument x);

}



#endif /* EXTRAS_OMITFUNCTIONARGUMENT_H_ */
