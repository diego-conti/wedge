/*******************************************************************************
 *  Copyright (C) 2007-2023 by Diego Conti, diego.conti@unipi.it 
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
#ifndef TORSIONFREECONNECTION_H_
#define TORSIONFREECONNECTION_H_

/** @ingroup RiemannianGeometry */ 

/** @{ 
 * @file torsionfreeconnection.h
 * @brief Torsion-free connections
 */
 
#include "wedge/connections/riemannianconnection.h"
#include "wedge/manifolds/concretemanifold.h"
namespace Wedge {


/** @brief Class template with explicit specializations.
 */
template<bool Use_d> class TorsionFreeConnection {};
	
/** @brief Template class for torsion-free connections on a manifold that (already) knows how to 
 * take d of its forms.*/  
template<> class TorsionFreeConnection<true> : 
	public virtual Connection
{
public:
/** @brief Construct a torsion-free connection
 *  @param manifold The manifold on which the connection is defined
 *  @param frame The frame used to represent the connection as a matrix
 *  @param christoffel The symbol to use for the connection parameters
 * 
 * @warning The pointer argument is stored internally. Caller is responsible for making sure 
 * that it remains valid.
 * @exception WedgeException<std::logic_error> Thrown if manifold does not know how to take d of its forms
*/
	TorsionFreeConnection(const Manifold* manifold, const Frame& frame, const Name& christoffel=N.Gamma) : Connection(manifold,frame, christoffel)
	{
		if (!manifold->KnowsHowToCompute_d())
		{
			LOG_MSG("Trying to construct a TorsionFreeConnection<true> object from a Manifold object that does not know how to take d of its forms. Was TorsionFreeConnection<false> intended?");
			throw WedgeException<std::logic_error>("KnowsHowToCompute_d returned false in TorsionFreeConnection<true>::TorsionFreeConnection()",__FILE__,__LINE__);
		}
		const ExVector torsion=Torsion();
		DeclareZero(torsion.begin(),torsion.end());
	}
protected:
//don't impose torsion zero
	TorsionFreeConnection(const Manifold* manifold, const Frame& frame, bool, const Name& christoffel=N.Gamma) : Connection(manifold,frame, christoffel)
	{
		if (!manifold->KnowsHowToCompute_d())
		{
			LOG_MSG("Trying to construct a TorsionFreeConnection<true> object from a Manifold object that does not know how to take d of its forms. Was TorsionFreeConnection<false> intended?");
			throw WedgeException<std::logic_error>("KnowsHowToCompute_d returned false in TorsionFreeConnection<true>::TorsionFreeConnection()",__FILE__,__LINE__);
		}
	}
};



/** @brief Specialization of TorsionFreeConnection for torsion-free connections on a manifold that doesn't know how to 
 * take d of its forms (yet).
 *
 * @note One cannot assume that manifold->KnowsTowToCompute_d() returns false, as it might be the case that the manifold uses this connection object
 * to compute d. @sa ManifoldWith
*/  
template<> class TorsionFreeConnection<false> : 
	public virtual Connection
{
public:
/** @brief Construct a torsion-free connection
 *  @param manifold The manifold on which the connection is defined
 *  @param frame The frame used to represent the connection as a matrix
 *  @param christoffel The symbol to use for the connection parameters
 * 
 * @warning The pointer argument is stored internally. Caller is responsible for making sure 
 * that it remains valid.
*/
	TorsionFreeConnection(const Manifold* manifold, const Frame& frame, const Name& christoffel=N.Gamma) : Connection(manifold,frame,christoffel) {}

/** @brief Impose conditions on the connection
 * @param alpha A differential form
 * @param beta A differential form
 * 
 * Impose conditions on the coefficients of the connection form in such a way that \f$d\alpha=\beta\f$.
 **/
	void Declare_d(ex alpha, ex beta);

/**
* @brief The exterior derivative, or d operator
* @param alpha A differential form on manifold
* @returns The exterior derivative of alpha
* 
* Use the fact that the connection is torsion-free to compute d from the covariant derivative
*/
	ex d(ex alpha) const;

/** @brief Curvature 2-form
 *  @return The curvature 2-form as a matrix
 * 
 *  This function differs from Connection::CurvatureForm() in that it is computed using solely
 *  the Christoffel symbols, rather than invoking Manifold::d 
 */
	matrix CurvatureForm() const;

/** @brief Compute the equations in the parameters corresponding to \f$d^2=0\f$
 * @param container (out) A reference to a container object where the equations are to be stored.
 * @return A reference to the container
 * 
 * @remark The simplest situation is when one knows the action of d on a basis of 1-forms, in which case
 * this function is useless. This class, however, allows for more general situations, e.g.
 * one might know the action of d on a basis of two-forms. Then, imposing \f$d^2=0\f$ gives extra
 * conditions, that are computed through this function.
 */
 	template<class Container> Container& GetEquations_ddZero(Container& container) const
 	{ 		
 		for (ExVector::const_iterator i=e().begin();i!=e().end();i++)
		{
			GetCoefficients<DifferentialForm>(container,d(d(*i)));			
		}		
		return container;
 	}
private:
/** Torsion() is a virtual function, so it has to be reimplemented to return zero.
 * However, it is defined as private because it makes no sense to invoke
 * Torsion() on a Connection object that is known to have type TorsionFreeConnection!
*/
	ExVector Torsion() const {return ExVector(manifold->Dimension());}

	TorsionFreeConnection(const Has_dTable* manifold, const Frame& frame);	///<Disabled. Use TorsionFreeConnection<true> instead.
};


/** @brief Torsion-free Riemannian connections. 
 * 
 * @remark The Levi-Civita connection is determined uniquely by the metric.
 */

//notice the order: RiemannianConnection's constructor must be called _before_ TorsionFreeConnection's
template<bool Use_d> class LeviCivitaConnection : public RiemannianConnection, public TorsionFreeConnection<Use_d>
{
public:
/** @brief Construct the Levi-Civita connection
 *  @param manifold The manifold on which the connection is defined
 *  @param g A riemannian structure on manifold
 *  @param christoffel The symbol to use for the connection parameters
 * 
 * @warning The pointer argument is stored internally. Caller is responsible for making sure 
 * that it remains valid.
*/
	LeviCivitaConnection(const Manifold* manifold, const RiemannianStructure& g,const Name& christoffel=N.Gamma);	
};

template<> LeviCivitaConnection<false>::LeviCivitaConnection(const Manifold* manifold, const RiemannianStructure& g,const Name& christoffel);
template<> LeviCivitaConnection<true>::LeviCivitaConnection(const Manifold* manifold, const RiemannianStructure& g,const Name& christoffel) ;
	
	
} /** @} */

#endif /*TORSIONFREECONNECTION_H_*/
