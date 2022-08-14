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
#include "../structures/pseudoriemannianstructure.h"

namespace Wedge {

pair<ex,matrix> StandardPseudoRiemannianStructure::DecomposeRicci(matrix ricci) const  {
	int dim=M()->Dimension();
	assert(ricci.rows()==dim);
	assert(ricci.cols()==dim);
	ex s;
	for (int i=0;i<p_;++i)
		s+=ricci(i,i);
	for (int i=p_;i<dim;++i)
		s-=ricci(i,i);
	ex normalized_s=s/dim;
	for (int i=0;i<p_;++i)
		ricci(i,i)-=normalized_s;
	for (int i=p_;i<dim;++i)
		ricci(i,i)+=normalized_s;
	return {normalized_s,ricci};
}

}
