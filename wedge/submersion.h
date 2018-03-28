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


#ifndef SUBMERSION_H
#define SUBMERSION_H
/** @ingroup RiemannianGeometry 
 *  @{ */

/** @file submersion.h 
 *  @brief The Levi Civita connection on the base of a Riemannian submersion
 */

#include "riemannianconnection.h"
#include "torsionfreeconnection.h"
#include "basis.h"

namespace Wedge {

typedef SubBasis<DifferentialOneForm,DefaultLinAlgAlgorithms> SubFrame;

/** @brief A class representing the base space of a submersion, or the space of leaves of a foliation
*
* This class is used to work with non-parallelizable manifolds
*/
class Submersion : public virtual Manifold , public Has_dTable {
	SubFrame coframe;	//one-forms
	lst horizontal,vertical;	//substitutions that map a vector field to its horizontal/vertical part
	lst horizontalf, verticalf;	//substitutions that map a form to its horizontal/vertical part
public:
/** @brief Create a Submersion object from a list of horizontal one-forms
 * @param hforms A global basis of horizontal one-forms
 *
 * Notice that the choice of the complement has no effect on the geometry
 */
	Submersion(const Has_dTable& totalspace, const exvector& hforms) :  Has_dTable(totalspace), coframe(hforms.begin(),hforms.end(),hforms.end(),hforms.end()) 
	{
		coframe.AddToComplement(totalspace.e().begin(),totalspace.e().end());
		assert(totalspace.Dimension()==coframe.DimensionOfContainingSpace());
		Initialize();
	}
/** @overload
*/
	template<typename Iterator> Submersion(const Has_dTable& totalspace, Iterator hforms_begin, Iterator hforms_end) : Has_dTable(totalspace), coframe(hforms_begin,hforms_end,hforms_end,hforms_end)
	{
		coframe.AddToComplement(totalspace.e().begin(),totalspace.e().end());
		assert(totalspace.Dimension()==coframe.DimensionOfContainingSpace());
		Initialize();
	}

/** @brief Create a Submersion object from a list of horizontal and vertical forms
 * @param hforms A global basis of horizontal one-forms
 * @param vforms A global basis of vertical one-forms
 */
	Submersion(const Has_dTable& totalspace, const exvector& hforms, const exvector& vforms) : Has_dTable(totalspace),  coframe(hforms.begin(),hforms.end(),vforms.begin(),vforms.end())
	{
		assert(totalspace.Dimension()==coframe.DimensionOfContainingSpace());
		Initialize();
	}

/** @overload
*/
	template<typename Iterator1, typename Iterator2> Submersion(const Has_dTable& totalspace, Iterator1 hforms_begin, Iterator1 hforms_end,
		Iterator2 vforms_begin, Iterator2 vforms_end) : Has_dTable(totalspace), coframe(hforms_begin,hforms_end,vforms_begin,vforms_end)
	{
		assert(totalspace.Dimension()==coframe.DimensionOfContainingSpace());
		Initialize();
	}
	
/** @brief Project a form or vector field on the vertical part
*/
	template<typename Section> ex V(ex X) const;

/** @brief Project a form or vector field on the horizontal part
*/
	template<typename Section> ex H(ex X) const;

/** @brief Determine whether a differential form or vector field is basic
 * @param X A differential form or vector field
 * @return true if X is basic
 */
	template<typename Section> bool IsBasic(ex X) const;

	const SubFrame& e() const {return coframe;}
	ex e(OneBased k) const {return Manifold::e(k);}

	int Dimension() const {return e().DimensionOfContainingSpace();}

	ex LieBracket(ex X, ex Y) const {
		ex XY;
		for (int k=1;k<=coframe.DimensionOfContainingSpace();k++)
		{
			ex XYk=LieDerivative(X,TrivialPairing<DifferentialForm>(e(k),Y))-LieDerivative(Y,TrivialPairing<DifferentialForm>(e(k),X));
			XYk-=Hook(Y,Hook(X,d(e(k)))); 
			XY+=XYk*e().dual()(k);
		}
		return XY;
	}
private:
	void Initialize()
	{
		LOG_INFO(coframe);
		exvector frame=coframe.dual();
		LOG_INFO(frame);
		{
			exvector images(frame.begin(),frame.begin()+coframe.size());
			images.resize(coframe.DimensionOfContainingSpace());	//pad with zeroes
			horizontal=LinearMapToSubstitutions<VectorField>(frame.begin(),frame.end(),images.begin(),images.end(),solve_algo::gauss);
		}
		{
			exvector images(coframe.size());	//initialize to zero vector
			images.insert(images.end(),frame.begin()+coframe.size(),frame.end());
			vertical=LinearMapToSubstitutions<VectorField>(frame.begin(),frame.end(),images.begin(),images.end(),solve_algo::gauss);
		}
		{
			exvector images(coframe.begin(),coframe.end());
			images.resize(coframe.DimensionOfContainingSpace());
			horizontalf=LinearMapToSubstitutions<DifferentialOneForm>(
				coframe.full_begin(),coframe.full_end(),images.begin(),images.end(),solve_algo::gauss);
		}
		{
			exvector images(coframe.size());
			images.insert(images.end(),coframe.complement_begin(),coframe.complement_end());
			verticalf=LinearMapToSubstitutions<DifferentialOneForm>(coframe.full_begin(),coframe.full_end(),images.begin(),images.end(),solve_algo::gauss);
		}
		LOG_INFO(horizontal);
		LOG_INFO(horizontalf);
		LOG_INFO(vertical);
		LOG_INFO(verticalf);
		
	}
};

template<> ex Submersion::V<VectorField>(ex X) const;
template<> ex Submersion::H<VectorField>(ex X) const;

template<> ex Submersion::V<DifferentialForm>(ex X) const;
template<> ex Submersion::H<DifferentialForm>(ex X) const;

template<> bool Submersion::IsBasic<DifferentialForm>(ex X) const;
template<> bool Submersion::IsBasic<VectorField>(ex X) const;


} /** @} */
#endif

