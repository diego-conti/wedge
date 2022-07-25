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


#ifndef SUBMERSIONWITH_H
#define SUBMERSIONWITH_H
/** @ingroup RiemannianGeometry 
 *  @{ */

/** @file submersionwith.h 
 *  @brief The Levi Civita connection on the base of a Riemannian submersion
 */

#include "../structures/submersion.h"
#include "wedge/connections/riemannianconnection.h"
#include "wedge/connections/torsionfreeconnection.h"
#include "wedge/linearalgebra/basis.h"

namespace Wedge {

/** @brief A GStructure on the space of leaves of a foliation, or base space of a submersion
 *
 * This class implements the O'Neill formulas for a Riemannian submersion, where \f$X\f$ is a parallelizable manifold. This enables
 * one to compute the curvature of the non-parallelizable manifold \f$M\f$. Also, one can compute the covariant derivative of a basic vector field or
 * differential form
 *
 * @todo Introduce support for spinors.
 * Notice that a submersion really represents a foliation where one looks at the space of leaves.
 * @todo In principle, Submersion should play the same role as Manifold. Classes GStructure and Connection should also take Submersion objects and behave accordingly.
 * The best way of doing this is probably making GStructure depend on a template parameter, e.g. the type of Manifold. Then the transverse G-structure template
 * specialization can also handle the problem of spinors. The class Connection should probably also depend on a parameter.
 * This would effectively introduce support for non-parallelizable manifolds in Wedge. For instance, it should be possible to compute the intrinsic torsion of a G-structure on a homogeneous space.
 * @notice This class is in no way analogous to ManifoldWith 
 */

template<typename Structure> class SubmersionWith : public Submersion, public Structure  {
	LeviCivitaConnection<true> omega;
public:
/** @brief Create a SubmersionWith object from a list of horizontal one-forms
 * @param hforms A global basis of horizontal one-forms
 *
 * Notice that the choice of the complement has no effect on the geometry
 */
	SubmersionWith(const Has_dTable& totalspace, const exvector& hforms) :  Submersion(totalspace, hforms), Structure(this, Submersion::e()), omega(this,*this)
	{
		LOG_INFO(omega);
	}
/** @overload
*/
	template<typename Iterator> SubmersionWith(const Has_dTable& totalspace, Iterator hforms_begin, Iterator hforms_end) :  Submersion(totalspace, hforms_begin,hforms_end), Structure(this, Submersion::e()), omega(this,*this)
	{
		LOG_INFO(omega);
	}

/** @brief Create a SubmersionWith object from a list of horizontal and vertical forms
 * @param hforms A global basis of horizontal one-forms
 * @param vforms A global basis of vertical one-forms
 */
	SubmersionWith(const Has_dTable& totalspace, const exvector& hforms, const exvector& vforms):  Submersion(totalspace, hforms, vforms), Structure(this, Submersion::e()), omega(this,*this)
	{
		LOG_INFO(omega);
	}

/** @overload
*/
	template<typename Iterator1, typename Iterator2> SubmersionWith(const Has_dTable& totalspace, Iterator1 hforms_begin, Iterator1 hforms_end,
		Iterator2 vforms_begin, Iterator2 vforms_end) :  Submersion(totalspace, hforms_begin,hforms_end, vforms_begin, vforms_end), Structure(this, Submersion::e()), omega(this,*this)
	{
		LOG_INFO(omega);
	}


/** @brief Compute the tensor T, \f$T(X,Y)=(\nabla_{X_V} Y_H)_V+(\nabla_{X_V} Y_V)_H\f$
 * @param X,Y %Vector fields on the total space 
*/
	ex T(ex X, ex Y) const;
/** @brief Compute the tensor A, \f$A(X,Y)=(\nabla_{X_H} Y_H)_V+(\nabla_{X_H} Y_V)_H\f$
 * @param X,Y %Vector fields on the total space 
*/
	ex A(ex X, ex Y) const;

/** @brief Covariant derivative of a basic form or tensor
 * @param X A horizontal vector field
 * @param alpha A basic form or tensor
 * @return The covariant derivative \f$\nabla_X\alpha\f$
 * 
 * @remark The return value is horizontal, but not necessarily basic, unless \f$X\f$ is also basic
*/
	template<typename Section> ex Nabla(ex X, ex alpha) const
	{
		assert(IsBasic<Section>(alpha));
		return omega.Nabla<Section>(X,alpha);
	}

/** @brief Compute the Ricci tensor of the Levi Civita connection on \f$M\f$
 * @return The Ricci tensor as a matrix
 * 
 * The result represents the Ricci tensor as a matrix. We use the formula of [Besse, Einstein manifolds]
**/
	matrix RicciAsMatrix() const;

/** @brief Compute the curvature 2-form of this connection
 *  @return The curvature 2-form as a GiNaC::matrix (whose indices are zero-based)
 * @exception WedgeException<std::logic_error> Thrown if the manifold does not know how to take d of its forms
 * 
 * Computes the curvature form  \f$\Omega^M\f$ in terms of the curvature of the total space \f$\Omega^M\f$ by the O'Neill formula
 * \f$2\Omega^X_{ij}(X,Y)=2\Omega^M_{ij}(X,Y)-2\langle A(X,Y),A_{ij}\rangle +\langle A(Y,e_i),A(X,e_j)\rangle -\langle A(X,e_i),A(Y,e_j)\rangle\f$. The
 * curvature \b tensor \f$R\f$ can be obtained by \f$\langle R(X,Y)e_i,e_j\rangle =2\Omega_{ji}(X,Y)\f$ 
 * 
 * @sa [O'Neill, The fundamental equation of a submersion], but notice O'Neill defines the curvature with the opposite sign
 */
	matrix CurvatureForm() const;
	
	const SubFrame& e() const {return Submersion::e();}
	ex e(OneBased i) const {return Manifold::e(i);}

//only handles vector fields. a bad hack because RiemannianStructure will not work correctly on the vertical component
	ex ScalarProduct(ex left, ex right) const {
		ex result;
		for (SubFrame::const_iterator i=e().full_begin();i!=e().full_end();++i)
			result+=Wedge::Hook(left,*i)*Wedge::Hook(right,*i);
		return result;
	}
};

template<typename Structure> ex SubmersionWith<Structure>::T(ex X, ex Y) const
{
	ex XV=V<VectorField>(X);
	return V<VectorField>(omega.Nabla<VectorField>(XV,H<VectorField>(Y)))+H<VectorField>(omega.Nabla<VectorField>(XV,V<VectorField>(Y)));
}

template<typename Structure> ex SubmersionWith<Structure>::A(ex X, ex Y) const
{
//	return V<VectorField>(LieBracket(X,Y))/2;	//this is only ok for horizontal vector fields
	ex XH=H<VectorField>(X);
	ex AXY=H<VectorField>(omega.Nabla<VectorField>(XH,V<VectorField>(Y)))+V<VectorField>(omega.Nabla<VectorField>(XH,H<VectorField>(Y)));
	LOG_INFO(X);
	LOG_INFO(Y);
	LOG_INFO(AXY);
	return AXY;
}

template<typename Structure> matrix SubmersionWith<Structure>::RicciAsMatrix() const
{
	matrix R=omega.RicciAsMatrix();
	LOG_INFO(R);
	matrix Ric(Dimension(),Dimension());
	ex N;
	SubBasis<VectorField> frame(e().dual().begin(),e().dual().begin()+Dimension(),e().dual().begin()+Dimension(),e().dual().end());
	for (exvector::const_iterator i=frame.complement_begin();i!=frame.complement_end();++i)
		N+=T(*i,*i);
	LOG_INFO(N);	
	for (int h=1;h<=frame.size();++h)
	for (int k=1;k<=frame.size();++k)
	{			
		ex X=frame(h),Y=frame(k);
		ex AXAY,TXTY;			
		for (exvector::const_iterator i=frame.begin();i!=frame.end();++i)
			AXAY+=this->ScalarProduct(A(X,*i),A(Y,*i));
		for (exvector::const_iterator i=frame.complement_begin();i!=frame.complement_end();++i)
			TXTY+=this->ScalarProduct(T(*i,X),T(*i,Y));
		Ric(h-1,k-1)=R(h-1,k-1)+2*AXAY+TXTY-
			(this->ScalarProduct(omega.Nabla<VectorField>(X,N),Y)+
					this->ScalarProduct(omega.Nabla<VectorField>(Y,N),X))/2;
	}
	return ex_to<matrix>(Ric.expand().normal());
}


template<typename Structure> matrix SubmersionWith<Structure>::CurvatureForm() const
{
	matrix R1=omega.CurvatureForm();
	LOG_DEBUG(R1);
	matrix R(Dimension(),Dimension());
	for (int i=0;i<Dimension();++i)
	for (int j=0;j<Dimension();++j)
		R(i,j)=H<DifferentialForm>(R1(i,j));
	LOG_DEBUG(R);
	ExVector frame=e().dual();
	for (int h=1;h<=Dimension();++h)
	for (int k=h+1;k<=Dimension();++k)
	{
		ex X=frame(h), Y=frame(k);
		for (int i=1;i<=Dimension();++i)
		for (int j=1;j<=Dimension();++j)
		{
			ex deltaXY=-2*this->ScalarProduct(A(X,Y), A(frame(i),frame(j)))
			+this->ScalarProduct(A(Y,frame(i)), A(X,frame(j)))
			-this->ScalarProduct(A(X,frame(i)), A(Y,frame(j)));
			LOG_DEBUG(deltaXY);
			R(i-1,j-1)-=deltaXY*e(h)*e(k);
		}
	}
	LOG_DEBUG(R);
	return R;
}

}
#endif 

