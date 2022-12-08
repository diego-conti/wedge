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

#include "wedge/representations/linearaction.h"

namespace Wedge {

ex LinearActionBase::eval_ncmul(const exvector & v) const
{

	exvector::const_reverse_iterator i=v.rbegin();
	assert(is_a<LinearActionBase>(*i));
	exvector product(ex_to<LinearActionBase>(*i).linear_subs.begin(),ex_to<LinearActionBase>(*i).linear_subs.end());
	while (++i!=v.rend())
		for (exvector::iterator k=product.begin();k!=product.end();++k)
		{
			ex rhs=k->rhs();
			rhs=rhs.subs(ex_to<LinearActionBase>(*i).linear_subs);
			k->let_op(1)=rhs;
		}
	return 	LinearActionBase (lst(product.begin(),product.end()));
}


ex GroupAction(ex g, ex w)
{
	internal::LinearActionVisitor v(w);
	return v.RecursiveVisit(g);
}


}
