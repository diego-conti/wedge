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
 #include "bilinearform.h"

namespace Wedge {

const ex AsDifferentialForm(ex x)
{
	if (is_a<DifferentialForm>(x)) return x;
	else if (is_a<DifferentialOneForm>(x)) return lst{x};
	else throw InvalidArgument(__FILE__,__LINE__,x);
}

ex BilinearForm::OnForms(ex v, ex w) const {
	v=v.expand(); w=w.expand();
	internal::NormalFormHelper<DifferentialForm> visit_v;
	v.accept(visit_v);
	internal::NormalFormHelper<DifferentialForm> visit_w;
	w.accept(visit_w);
	ex result;
	LOG_INFO(v);
	LOG_INFO(w);
	for (exmap::const_iterator i=visit_v.coeffs.begin();i!=visit_v.coeffs.end();++i)
	for (exmap::const_iterator j=visit_w.coeffs.begin();j!=visit_w.coeffs.end();++j)
	{
		result+=i->second * j->second * OnSimpleForms(AsDifferentialForm(i->first),AsDifferentialForm(j->first));
	}
	LOG_INFO(result);
	return result.expand();
}

ex BilinearForm::OnSimpleForms(ex v, ex w) const {
	if (v.nops()!=w.nops()) return 0;
	else {
		matrix m(w.nops(),w.nops());
		for (int i=0;i<w.nops();i++)
			for (int j=0;j<w.nops();j++)
				m(i,j)=	OnOneForms(v.op(i),w.op(j));
		return m.determinant();
	}
}


ex BilinearForm::Interior(ex v, ex w) const {
	v=Sharp(v).expand(); w=w.expand();
	internal::NormalFormHelper<DifferentialForm> visit_v;
	v.accept(visit_v);
	internal::NormalFormHelper<DifferentialForm> visit_w;
	w.accept(visit_w);
	ex result;
	for (exmap::const_iterator i=visit_v.coeffs.begin();i!=visit_v.coeffs.end();++i)
	for (exmap::const_iterator j=visit_w.coeffs.begin();j!=visit_w.coeffs.end();++j)
	{
		result+=i->second * j->second * Hook(i->first,j->first);
	}
	return result;
}

}
