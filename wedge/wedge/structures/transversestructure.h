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
#ifndef TRANSVERSE_H
#define TRANSVERSE_H
/** @ingroup RiemannianGeometry 
 *  @{ */

/** @file transverse.h 
 *  @brief Tranverse Riemannian structures
 */

#include "../structures/pseudoriemannianstructure.h"
#include "wedge/connections/riemannianconnection.h"
#include "wedge/connections/torsionfreeconnection.h"
#include "wedge/linearalgebra/basis.h"
#include "wedge/convenience/simplifier.h"

namespace Wedge {

/** @brief A class for O(p)O(q) structures, i.e. Riemannian structure where the tangent space has an orthogonal splitting \f$TM=H\oplus V\f$

Can be used for Riemannian foliations, subriemannian structures...
*/
class TransverseRiemannianStructure : public StandardPseudoRiemannianStructure {
public:
	TransverseRiemannianStructure(const Manifold* M, const Frame& frame, int k);

	TransverseRiemannianStructure(const Manifold* M, const SubBasis<DifferentialOneForm>& frame) : TransverseRiemannianStructure {M,Frame{frame.full_begin(),frame.full_end()},frame.size()} {}

	template<typename Iterator, typename Iterator2> TransverseRiemannianStructure(const Manifold* M, Iterator hforms_begin, Iterator hforms_end, Iterator2 vforms_begin, Iterator2 vforms_end) : 
 		TransverseRiemannianStructure (M,CreateFrame(hforms_begin,hforms_end,vforms_begin,vforms_end),distance(hforms_begin,hforms_end)) {}

	void Reduce(const Simplifier& simplifier);
//	template<typename Iterator> TransverseRiemannianStructure(const Submersion* totalspace, Iterator hforms_begin, Iterator hforms_end) : RiemannianStructure((totalspace), coframe(hforms_begin,hforms_end,hforms_end,hforms_end)
 		
/** @brief Project a form or vector field on the vertical part
*/
	template<typename Section> ex V(ex X) const;

/** @brief Project a form or vector field on the horizontal part
*/
	template<typename Section> ex H(ex X) const;

	int BaseDimension() const {return k;}

/** @brief Throw if the frame does not define a Riemannian foliation (with leaves tangent to the horizontal distribution)
*/
	void AssertIntegrable() const 
	{
		lst result;
		GetEquationsIntegrable(result, true);
	}
/** @brief Return a list of equations that hold iff a Riemannian foliation is defined
 *
 * The conditions imposed are that the horizontal distribution is integrable and the transverse metric tensor is basic
 */
	template<typename Container> Container& GetEquationsIntegrable(Container& container, bool shouldthrow=false) const {
		//verify the distribution is integrable
		ex vol=e()[0];
		for (int i=1;i<k;++i)
			vol*=e()[i];
		vol=vol.expand();
		LOG_INFO(vol);
		for (exvector::const_iterator i=e().begin();i!=e().begin()+k;++i)
		{
			ex x=(M()->d(*i) * vol).expand();
			GetCoefficients<DifferentialForm>(container,x);
			if (!x.is_zero()) {
				LOG_WARN(e());
				LOG_WARN(k);
				LOG_WARN(x);
				if (shouldthrow)
					throw WedgeException<std::invalid_argument>("Horizontal distribution is not integrable",__FILE__,__LINE__);
			}
		}
		ExVector frame=e().dual();
		//verify the horizontal part of metric tensor is basic
		for (int i=1;i<=k;++i)
		for (int j=i+1;j<=k;++j) {
			ex x=Wedge::Hook(frame(j),M()->d(e(i)))+Wedge::Hook(frame(i),M()->d(e(j)));
			x=(x*vol).expand();
			GetCoefficients<DifferentialForm>(container,x);
			if (!x.is_zero()) {
				LOG_WARN(e());
				LOG_WARN(vol);
				LOG_WARN(k);
				LOG_WARN(e(i));
				LOG_WARN(e(j));
				LOG_WARN(x);
				if (shouldthrow)
					throw WedgeException<std::invalid_argument>("Horizontal metric tensor is not basic",__FILE__,__LINE__);
			}
		}
		return container;
	}

	ex HodgeStar(ex form) const {exvector h(e().begin(),e().begin()+k); return ScalarProduct().Interior(form,ncmul(h));}

private:
	const int k;
	lst horizontal,vertical;	//substitutions that map a vector field to its horizontal/vertical part
	lst horizontalf, verticalf;	//substitutions that map a forn to its horizontal/vertical part

	template<typename Iterator, typename Iterator2> exvector CreateFrame(Iterator hforms_begin, Iterator hforms_end, Iterator2 vforms_begin, Iterator2 vforms_end) {
		exvector f(hforms_begin, hforms_end);
		f.insert(f.end(),vforms_begin,vforms_end);
		return f;
	}
};

template<> ex TransverseRiemannianStructure::V<VectorField>(ex X) const;
template<> ex TransverseRiemannianStructure::H<VectorField>(ex X) const;
template<> ex TransverseRiemannianStructure::V<DifferentialForm>(ex X) const;
template<> ex TransverseRiemannianStructure::H<DifferentialForm>(ex X) const;




}

#endif

