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

#ifndef LIEGROUPEXTENSION_H_
#define LIEGROUPEXTENSION_H_

#include "liegroup.h"
#include "polybasis.h"
namespace Wedge {


/** @ingroup Manifolds
 *  @{
 * 
 *  @file liegroupextension.h
 *  @brief Lie group extensions
*/

/** @brief Generic extension of a given Lie group
 * 
 * An extension of a Lie group \f$G\f$ is a Lie group that contains \f$G\f$. This class represents the generic extension of
 * a given dimension, and so its structure constants depend on parameters. 
 */
class LieGroupExtension :  public ConcreteManifold, public LieGroupHasParameters<true>, public virtual Has_dTable {
	ExVector CreateFrame(const LieGroup& G, int comprehensive_dimension)
	{
		ExVector e=G.e();
		if (comprehensive_dimension<G.Dimension()) throw InvalidArgument(__FILE__,__LINE__,comprehensive_dimension);
		for (int i=1;i<=comprehensive_dimension-G.Dimension();i++)
			e.push_back(DifferentialOneForm(N.E(i)));
		return e;
	}	
public:	
/** @brief Construct a Lie group extension \f$H\f$
 * @param G The Lie group to extend
 * @param comprehensive_dimension The dimension of \f$H\supset G\f$
 */
	LieGroupExtension(const LieGroup& G, int comprehensive_dimension) : ConcreteManifold(CreateFrame(G,comprehensive_dimension))
	{
		SubBasis<DifferentialOneForm> basis(G.e().begin(),G.e().end(),e().begin()+G.Dimension(),e().end());
		Subspace<DifferentialForm> V=TwoForms(basis);
		
		for (int i=1;i<=G.Dimension();++i) {
			ex de=G.d(e(i));
			int j=1;
			for (Frame::const_iterator k=V.complement_begin();k!=V.complement_end();++j,++k)
				de+=StructureConstant(N.a(i,j))* *k;
			Has_dTable::Declare_d(e(i),de);
		}		
		
		for (int i=G.Dimension()+1;i<=comprehensive_dimension;++i)
		{
			ex de;
			for (int j=1;j<=V.Dimension();++j)
				de+=StructureConstant(N.a(i,j))*V.e(j);
			Has_dTable::Declare_d(e(i),de);
		}		
	}
	
};


} /** @} */
#endif /*LIEGROUPEXTENSION_H_*/
