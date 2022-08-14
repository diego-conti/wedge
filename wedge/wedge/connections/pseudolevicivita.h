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
#ifndef PSEUDOLEVICIVITA_H_
#define PSEUDOLEVICIVITA_H_

/** @ingroup RiemannianGeometry */ 

/** @{ 
 * @file pseudolevicivita.h
 * @brief Levi-Civita connection on a pseudoriemannian manifold
 */
 
#include "wedge/structures/pseudoriemannianstructure.h"
#include "wedge/connections/torsionfreeconnection.h"

namespace Wedge {
	
/** @brief The Levi-Civita connection  on a pseudo riemannian manifold that (already) knows how to 
 * take d of its forms.*/  
class PseudoLeviCivitaConnection : 
	public virtual Connection, public TorsionFreeConnection<true>
{
public:
/** @brief Construct the Levi-Civita connection
 *  @param manifold The manifold on which the connection is defined
 *  @param structure A pseudoriemannian structure on manifold
 *  @param christoffel The symbol to use for the connection parameters
 * 
 * @warning The pointer argument is stored internally. Caller is responsible for making sure 
 * that it remains valid.
 * @exception WedgeException<std::logic_error> Thrown if manifold does not know how to take d of its forms
*/
	PseudoLeviCivitaConnection(const Manifold* manifold, const PseudoRiemannianStructure& structure, const Name& christoffel=N.Gamma) : 
		Connection(manifold,structure.e(),true), 
		TorsionFreeConnection<true>(manifold,structure.e(),true,christoffel)
	{
		const int dimension=e().size();
		components.reserve(dimension);
		for (int i=0;i<dimension;i++)
		components.push_back(exvector(dimension));

		ExVector e_flat(dimension);	//i-th element represents ((e_i)^\flat))
		ExVector frame=structure.e().dual();
		for (int i=1;i<=dimension;++i)
		for (int j=1;j<=dimension;++j)
			e_flat(i)+=structure.ScalarProduct().OnVectors(frame(i),frame(j))*structure.e(j); 
		LOG_INFO(e_flat);
		exvector de;	//i-th element represents d((e_i)^\flat))
		de.reserve(dimension);
		for (int i=1;i<=dimension;++i)
			de.push_back(manifold->d(e_flat(i)).normal());
		for (int i=0;i<dimension;++i)
		for (int j=0;j<dimension;++j)
		for (int k=0;k<dimension;++k)
		{

		ex X=frame[i];
		ex Y=frame[j]; 
		ex Z=frame[k];
		ex XYZ=TrivialPairing<DifferentialForm>(X*Y,-de[k])+
				TrivialPairing<DifferentialForm>(Z*X,-de[j])+
				TrivialPairing<DifferentialForm>(Z*Y,-de[i]);
		for (int h=0;h<dimension;++h)
			(*this)(h,j)+=(e()[i]*XYZ/2*structure.ScalarProduct().OnOneForms(e()[h],e()[k])).expand();

		LOG_DEBUG((*this)(k,j));
		}
	}
};

}

#endif /*PSEUDOLEVICIVITA_H_*/

