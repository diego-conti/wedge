/***************************************************************************
 *   Copyright (C) 2007, 2008, 2009 by Diego Conti, diego.conti@unimib.it  *
 *                                                                         *
 *   This file is part of Wedge.                                           *
 *                                                                         *
 *   Wedge is free software; you can redistribute it and/or modify         *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 3 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   Wedge is distributed in the hope that it will be useful,              *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with Wedge; if not, write to the                                *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#ifndef COORDINATES_H_
#define COORDINATES_H_

/** @ingroup Manifolds 
 *  @{ */
/**  @file coordinates.h
 * @brief Manifold with coordinates
 */
#include "wedge/manifolds/concretemanifold.h"
namespace Wedge {

/** @brief A manifold represented by local coordinates
 * 
 * To define a manifold represented by both coordinates and a frame (e.g. the product of a Lie
 * group and a given "manifold with coordinates"), derive from  ManifoldWithCoordinates and
 * call the constructor ManifoldWithCoordinates::ManifoldWithCoordinates(exvector, int )    
 */
class ManifoldWithCoordinates : public ConcreteManifold, public virtual Has_dTable {
	ExVector coordinates;
	ExVector CreateFrame(exvector oneforms, int no_coordinates)
	{
		ExVector frame;
		frame.reserve(frame.size()+no_coordinates);
		frame=oneforms;
		for (int i=1;i<=no_coordinates;i++)
		{
			frame.push_back(DifferentialOneForm(Name("dx")(i)));
		}
		return frame;
	}
public:
/** @brief Construct a manifold, to be represented by a coordinate patch.  
 * @param dimension The dimension of the manifold
 */
	ManifoldWithCoordinates(int dimension) : ConcreteManifold(CreateFrame(exvector(),dimension)) {
		coordinates.reserve(dimension);
		for (int i=1;i<=dimension;i++)
		{
			Function x_i=Function(N.x(i));
			coordinates.push_back(x_i);
			Declare_d(x_i,e(i));
			Declare_d(e(i),0);
		}
	};

/** @brief Construct a manifold represented by a local frame, some components of which arise from coordinates   
 * @param oneforms A set of independent one-forms
 * @param no_coordinates The number of coordinate functions that are defined
 * 
 * The dimension of the manifold will be given by the size of oneforms incremented by no_coordinates
 */	
	ManifoldWithCoordinates(exvector oneforms, int no_coordinates) : ConcreteManifold(CreateFrame(oneforms, no_coordinates)) {
		coordinates.reserve(no_coordinates);
		for (int i=1;i<=no_coordinates;i++)
		{
			ex x_i=Function(N.x(i));
			coordinates.push_back(x_i);
			Declare_d(x_i,e(i+oneforms.size()));
			Declare_d(e(i+oneforms.size()),0);
		}
	};

	const ExVector& x() const {return coordinates;} ///< Return a vector of coordinates
	ex x(OneBased i) const {return coordinates.at(i-1);} ///< Return the i-th coordinate
};
} 
/** @} */
#endif /*COORDINATES_H_*/
