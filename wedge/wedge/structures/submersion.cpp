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
#include "../structures/submersion.h"

namespace Wedge {

template<> ex Submersion::V<VectorField>(ex X) const
{
	return X.subs(vertical);
}

template<> ex Submersion::H<VectorField>(ex X) const
{
	return X.subs(horizontal);
}


template<> ex Submersion::V<DifferentialForm>(ex X) const
{
	return X.subs(verticalf);
}

template<> ex Submersion::H<DifferentialForm>(ex X) const
{
	return X.subs(horizontalf);
}


template<> bool Submersion::IsBasic<DifferentialForm>(ex X) const
{
	if (!V<DifferentialForm>(X).is_zero()) return false;
	for (exvector::const_iterator i=e().dual().begin()+coframe.size();i!=e().dual().end();++i)
		if (!LieDerivative(*i,X).expand().is_zero()) return false;
	return true;
}

template<> bool Submersion::IsBasic<VectorField>(ex X) const
{
	if (!V<VectorField>(X).is_zero()) return false;
	for (exvector::const_iterator i=e().dual().begin()+coframe.size();i!=e().dual().end();++i)
		if (!LieBracket(*i,X).expand().is_zero()) return false;
	return true;
}


}

