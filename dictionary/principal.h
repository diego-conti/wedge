/***************************************************************************
 *   Copyright (C) 2009 by Diego Conti, diego.conti@unimib.it              *
 *                                                                         *
 *   This file is part of Dictionary.                                      *
 *                                                                         *
 *   Dictionary is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 3 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   Dictionary is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with Dictionary; if not, write to the                           *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#ifndef PRINCIPAL_H_
#define PRINCIPAL_H_
/** @file principal.h
 * @brief Represents a principal bundle of the form \f$G\times V\to G\times_H V\f$, where \f$H\f$ is a subgroup of \f$G\f$
 */

#include <wedge/manifold.h>
#include <wedge/liegroup.h>
#include <wedge/function.h>
#include <wedge/logging.h>
#include "invariantforms.h"
#include "vvaluedforms.h"


using namespace Wedge;

template<class Representation, int dimOfStructureGroup> class Principal : public ConcreteManifold, public Has_dTable {
public:
	typedef ::NamedVValuedForm<Representation> NamedVValuedForm;	///< The type of "letters", i.e. elements of \f$\mathcal{L}\f$

	Principal(const LieGroup* _G) : ConcreteManifold(CreateFrame(_G)),  connectionForm(dimOfStructureGroup),  r(N.r), t(N.t), coordinates(N.a)
	{
		G=_G;
		for (int i=0;i<G->Dimension();++i)
		{
			Declare_d(e()[i],G->d(e()[i]));	
		}
		for (int i=0;i<Representation::FibreSize();i++)
		{
			ex a_i=Function(N.a(i+1));
			coordinates[i]=a_i;
			ex da_i=e()[i+G->Dimension()];
			Declare_d(a_i,da_i);
			Declare_d(da_i,0);
			dt+=2*a_i*da_i;
		}
	//not needed by Dictionary::d, but needed to construct a Submersion object.
		Declare_d(t,dt);	
		Declare_d(r,1/(2*r)*dt);	//can be overridden
	} 
	//not needed by Dictionary::d, but needed to construct a Submersion object.
	void Declare_dr(ex dr)
	{
		Declare_d(r,dr);
	}

/** @brief Returns a symbolic function  \f$t=aa\f$ that represents the radius squared
 */	
	ex RadiusSquared()  const	
	{
		return t;
	}	

/** @brief Returns a differential form representing   \f$dt=d(aa)\f$
 */	
	ex dRadiusSquared()  const	
	{
		return dt;
	}	

/** @brief Returns a radial function which may or may not coincide with the radius \f$r=\sqrt{aa}\f$
 */	
	ex RadialCoordinate()  const
	{
		return r;
	}

/** @brief Exterior covariant derivative
 * @param alpha A form in \f$\Omega^p(G\times V,V)\f$
 * @return The exterior covariant derivative \f$D\alpha\in \Omega^{p+1}(G\times V,V)\f$
 */ 
	NamedVValuedForm D(const NamedVValuedForm& alpha) const
	{		
		NamedVValuedForm result=Representation::LieAlgebraAction(connectionForm,alpha);
		result.set_name("D"+alpha.name(),"D"+alpha.texname());
		for (int i=0;i<Representation::FibreSize();i++)
		{
			result[i]+=d(alpha[i]);
		}
		return result;			
	}

	const NamedVValuedForm& a() const {return coordinates;} ///< Canonical element of \f$\mathcal{L}\f$

	ExVector connectionForm;			///< The connection form on the group G (initialized during construction by Dictionary-derived class)
	VectorSpace<DifferentialForm> horizontalForms;	///< The space of horizontal one-forms on the group G (initialized during construction by Dictionary-derived class)

private:
	static ExVector CreateFrame(const LieGroup* G)
	{
		ExVector frame(G->e());
		for (int i=0;i<Representation::FibreSize();++i)
		{
			ex da_i=DifferentialOneForm(Name("da")(i+1));
			frame.push_back(da_i);
		}
		return frame;
	}		
	ex dt;
	Function r,t;
	NamedVValuedForm coordinates;			///< Coordinates on the fibre
	const LieGroup* G;
};

#endif

