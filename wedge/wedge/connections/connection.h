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
#ifndef WEDGECONNECTION_H
#define WEDGECONNECTION_H

/** @ingroup RiemannianGeometry */ 

/** @{ 
 * @file connection.h
 * @brief Connections, covariant derivatives
 */


#include "wedge/structures/riemannianstructure.h"
#include "wedge/structures/spinor.h"
#include "wedge/manifolds/differentialform.h"
#include "wedge/manifolds/manifold.h"
#include "wedge/linearalgebra/linear.h"
#include "wedge/linearalgebra/derivation.h"
#include "wedge/base/parameters.h"

namespace Wedge {

class ConnectionParameter : public Register<ConnectionParameter, Function>::Algebraic
{
	friend class HasParameterList<ConnectionParameter>;
	ConnectionParameter(const Name& name) : Register<ConnectionParameter, Function>::Algebraic(name) {}
public:
	static const char* static_class_name() {return "ConnectionParameter";}	
};

template<class Structure> class ManifoldWith;

class Connection;

namespace internal {

template<typename Section> class CovariantDerivative;

template<> class CovariantDerivative<DifferentialForm> :
	public IBilinearOperator<LinearOperator<VectorField>, DerivationOver<DifferentialForm,Function,false> >
{
	const Connection& connection;
public:
	CovariantDerivative(const Connection& c) : connection(c) {}
	ex Apply(const VectorField& vfield, const VectorField& Alpha) const;
	ex Apply(const VectorField& vfield, const Function& Alpha) const;
};

template<> class CovariantDerivative<VectorField> :
	public IBilinearOperator<LinearOperator<VectorField>, Leibniz<VectorField,Function> >
{
	const Connection& connection;
public:
	CovariantDerivative(const Connection& c) : connection(c) {}
	ex Apply(const VectorField& vfield, const VectorField& Alpha) const; 
	ex Apply(const VectorField& vfield, const Function& Alpha) const;	
};

}


/** @brief A connection on the tangent bundle of a manifold
 * 
 * A connection is represented as a connection form \f$\omega_{jk}\f$ with respect to a given frame \f$e_i\f$, i.e. 
 * \f$\nabla e_j=\sum_i\omega_{ij}\otimes e_i\f$. The manifold and frame
 * are part of the Connection object, and the frame need not be the standard frame for the manifold.
 * 
 * In applications, the connection form of a connection is typically not completely known: this
 * is because the "adapted frame" may be uniquely determined only up to a group action. This is reflected in Connection
 * by the presence of the superclass HasParameters. In practice, a connection is initially created as a matrix depending 
 * on parameters of type ConnectionParameter, which are then determined in part or in full by imposing conditions on the covariant derivative.
 *
 * @note The components of \f$\omega_{jk}\f$ are represented as linear combinations of 1-forms on the manifold, i.e.
 * \f$\omega_{ij}=\Gamma_{kji}E^k\f$, where \f$E^k\f$ is the standard frame of the manifold. If  the standard frame of the manifold coincides with the 
 * adapted frame of the connection, i.e. \f$E^k=e^k\f$, the \f$\Gamma_{kji}\f$ are the usual Christoffel symbols.
*/

class Connection : 
	public HasParameterList<ConnectionParameter>
{
	template<typename Section> friend class internal::CovariantDerivative;
public:
/** @brief Construct a generic connection
 *  @param manifold The manifold on which the connection is defined
 *  @param frame The frame used to represent the connection as a matrix
 *  @param christoffel The symbol to use for the connection parameter
 *
 * @warning The pointer argument is stored internally. Caller is responsible for making sure 
 * that it remains valid.
*/
	Connection(const Manifold* manifold, const Frame& frame, const Name& christoffel=N.Gamma);

	virtual ~Connection() {}
/** @brief Return the connection form as a matrix (indices are zero-based)
  * @return A matrix object whose entries are the components of the connection form relative to the frame specified upon construction
 */
	matrix AsMatrix() const;

/** @brief More efficient alternative to AsMatrix()(i,j)
 * @param i,j Zero-based indices
 * @returns The (i,j) component of the connection form relative to the connection's frame
 */
	ex operator()(ZeroBased i,ZeroBased j) const {return components[i][j];}

/** @brief Let the caller modify the connection form directly
 * @param i,j zero-based indices
 * @returns The (i,j) component of the connection form relative to the connection's frame, as an lvalue
 */
	ex& operator()(ZeroBased i,ZeroBased j) {return components[i][j];}

/** @brief Impose conditions on the Christoffel symbols
 * @param X A vector field
 * @param Section The C++ type of a section of a vector bundle on which the covariant derivative is defined, i.e. DifferentialForm or VectorField
 * @param alpha A section of the bundle associated to Section
 * @param beta A section of the bundle associated to Section
 * @exception InconsistentDeclaration Thrown if the conditions are incompatible with the preexisting conditions
 * 
 * Impose conditions on the coefficients of the connection form in such a way that \f$\nabla_X\alpha=\beta\f$
**/
	template<typename Section> void DeclareNabla(ex X, ex alpha, ex beta)
	{
		alpha=alpha.expand().normal(); 
		beta=beta.expand().normal();
		if (alpha.is_zero()) {
			if (beta.is_zero()) return;
			else throw InconsistentDeclaration(__FILE__,__LINE__,"Christoffel symbols");
		}
		ex GenericNablaAlpha=Nabla<Section>(X,alpha);
		DeclareZero(GenericNablaAlpha-beta);
	}

/** @brief Compute covariant derivative
 * @param X A vector field
 * @param Section The C++ type of a section of a vector bundle on which the covariant derivative is defined, i.e. DifferentialForm or VectorField
 * @param alpha A section of the bundle associated to Section
 * @return The covariant derivative \f$\nabla_X\alpha\f$
 *
 * Implements the formulae, \f$\nabla e_j=\sum_i\omega_{ij}\otimes e_i\f$, \f$\nabla e^i=-\sum_i\omega_{ij}\otimes e^j\f$.
**/
	template<typename Section> ex Nabla(ex X, ex alpha) const
	{
		internal::CovariantDerivative<Section> d(*this);
		return internal::CovariantDerivative<Section>::BilinearOperator(X,alpha,&d).expand();
	}

/** @brief Compute the curvature 2-form of this connection
 *  @return The curvature 2-form as a GiNaC::matrix (whose indices are zero-based)
 * @exception WedgeException<std::logic_error> Thrown if the manifold does not know how to take d of its forms
 * 
 * Computes the curvature \b form  \f$\Omega=d\omega+\frac12[\omega,\omega]\f$, i.e. 
 * \f$\Omega_{ij}=d\omega_{ij}+\sum_k{\omega_{ik}\wedge\omega_{kj}}\f$. The
 * curvature \b tensor \f$R\f$ can be obtained by \f$R(e_i,e_j)=2\Omega(e_i,e_j)\f$ 
 * 
 * @sa [Kobayashi-Nomizu: Foundations of Differential Geometry, Wiley, 1963-69]
 */
	virtual matrix CurvatureForm() const;

/** @brief Compute the Ricci tensor of this connection, \f$\operatorname{Ric}=\sum_{i,h,k,h} 2\Omega_{ih}(e_i,e_k) e^k\otimes e^h\f$
 * @return The Ricci tensor as a section of \f$T^*M\otimes T^*M\f$
 * 
 * The result is a linear combination of Tensor<DifferentialOneForm,DifferentialOneForm>
**/
	ex Ricci() const;

/** @brief Compute the Ricci tensor of this connection, \f$\operatorname{Ric}_{kh}=\sum_{i,h,k,h} 2\Omega_{ih}(e_i,e_k)\f$
 * @return The Ricci tensor as a matrix
 * 
 * The result represents the Ricci tensor as a matrix
**/
	matrix RicciAsMatrix() const;
	
/** @brief Compute the torsion of this connection
 * @return The torsion as a vector of two-forms
 * @exception WedgeException<std::logic_error> Thrown if the manifold does not know how to take d of its forms
 *
 * Computes the torsion form \f$\Theta_i=de^i+\omega_{ij}\wedge e^j\f$.  The
 * torsion \b tensor \f$T\f$ can be obtained by \f$R(e_i,e_j)=2\Theta(e_i,e_j)\f$ 
**/
	virtual ExVector Torsion() const;
	
protected:
	const Frame& e() const {return frame;}	///< @brief The frame associated to this connection
/** @brief Returns the k-th element of the frame associated to this connection
 * @param k An integer in the interval [1,dimension]
 * @return An ex object representing a differential 1-form
*/
	ex e(OneBased k) const {
		return e()(k);
	};

	Connection(const Manifold* manifold, const Frame& frame,bool DoNotInitialize); ///< Overloaded constructor for subclasses that initialize the connection matrix themselves
	
	void DeclareConditions(const lst& eqns);

	const Manifold* manifold; 	///< Pointer to the manifold object to which this connection refers
	
	vector<exvector> components; ///< The element components[i][j] represents omega(i,j), where i,j are zero-based indices
private:
	template<class Structure> friend class ManifoldWith;
	Frame frame;	///< The "adapted" frame, i.e. the frame associated to this connection
};

/** @brief Overloaded output operator */
inline ostream& operator<<(ostream& os,const Connection& connection)
{
	return os<<connection.AsMatrix();
}

} /** @} */
#endif
