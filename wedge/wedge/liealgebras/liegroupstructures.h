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
//hack to make doxygen document this file
#ifdef DOXYGEN_RUNNING
#define LIEGROUP_H
#define STRUCTURES_H_
#endif 

#ifdef LIEGROUP_H
#ifdef STRUCTURES_H_
#ifndef LIEGROUPSTRUCTURES_H_
#define LIEGROUPSTRUCTURES_H_

/** @ingroup RiemannianGeometry
 *  @{
 * 
 *  @file liegroupstructures.h
 *  @brief Structures on  Lie groups
*/


namespace Wedge {
/** @brief Template class for a fixed Lie group with a fixed G-structure
 * @param Structure A class derived from GStructure
 * 
 * The Lie group is defined in terms of structure constants.
 */
 
template<class Structure> class LieGroupWith : public AbstractLieGroup<> {
public:
	const GStructureHasParameters<Structure,false> P;	///< The G-structure associated to this group
	/** @brief Constructor  
	 *  @param structure_constants A string expressing the structure constants in Salamon's notation
	 *  @param frame A string expressing the frame in Salamon's notation
	 * 
	 * @sa ParseDifferentialForms(const exvector& frame, const char* to_parse)
	 */
	LieGroupWith(const char* structure_constants, const char* frame) : 
		AbstractLieGroup(structure_constants),
		P(this,ParseDifferentialForms(e(),frame))
		{
			assert(Structure::dimension==Dimension());
		}
};

} /** @} */
#endif /*LIEGROUPSTRUCTURES_H_*/

#endif
#endif
