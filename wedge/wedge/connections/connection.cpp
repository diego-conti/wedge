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

#include "connection.h"

#include "wedge/structures/riemannianstructure.h"
#include "wedge/structures/submersion.h"
#include "wedge/connections/riemannianconnection.h"
#include "wedge/connections/torsionfreeconnection.h"
#include "wedge/manifolds/manifold.h"
#include "wedge/linearalgebra/tensor.h"

namespace Wedge {
 
////////////////////////////////////////////////////////////
//			Connection						 
////////////////////////////////////////////////////////////

Connection::Connection(const Manifold* manifold, const Frame& _frame,bool)  : frame(_frame)
{
	this->manifold=manifold;
}

Connection::Connection(const Manifold* manifold, const Frame& _frame, const Name& christoffel)  : frame(_frame)
{
	this->manifold=manifold;
	const int dimension=frame.size();
	components.reserve(dimension);
	for (int i=0;i<dimension;i++)
	{
		components.push_back(exvector(dimension));
		for (int j=0;j<dimension;j++)
		{
			ex omegaij;
			//here we would want to use dimension for a submersion, but manifold->Dimension() for a connection on a vector bundle
			for (int k=0;k<dimension;k++)
				omegaij+=Parameter(christoffel(k+1,j+1,i+1))*manifold->e()[k];
			components[i][j]=omegaij;
		}
	}
}

matrix Connection::AsMatrix() const
{
	const int dimension=e().size();
	matrix m(dimension,dimension);
	for (int i=0;i<dimension;i++)
		for (int j=0;j<dimension;j++)
			m(i,j)=components[i][j].normal();
	return m;
}

ex Connection::Ricci() const
{
	matrix R=CurvatureForm();
	ex ric;	
	for (int j=0;j<e().size();++j) {
		ex sum;
		for (int i=0;i<e().size();++i)
			sum+=Hook(e().dual()[i],R(i,j));
		ric+=TensorProduct<DifferentialOneForm,DifferentialOneForm>(sum.expand(),e()[j]);
	}
	return ric;
}

matrix Connection::RicciAsMatrix() const
{
	LOG_DEBUG(e());
	LOG_DEBUG(e().dual());	
	matrix R=CurvatureForm();
	matrix ric(e().size(),e().size());	
	for (int j=0;j<e().size();++j) {
		ex sum;
		for (int i=0;i<e().size();++i)
			sum+=Hook(e().dual()[i],R(i,j));
		for (int k=0;k<e().size();++k)
			ric(k,j)=Hook(e().dual()[k],sum);
	}
	return ric;
}

ExVector Connection::Torsion() const
{
	ExVector DTheta(frame.size());
	try {
		for (int i=0;i<frame.size();i++)
		{
			DTheta[i]=manifold->d(frame[i]);
			for (int j=0;j<frame.size();j++)
				DTheta[i]+=operator()(i,j)*frame[j];
		}
	}
	catch (const Manifold::dException&)
	{
		throw WedgeException<logic_error>("Manifold does not know how to take d of its forms",__FILE__,__LINE__);	
	}
	return DTheta;
}

matrix Connection::CurvatureForm() const
{
	const int dimension=e().size();
	matrix m(dimension,dimension);
	try {	
		for (int i=0;i<dimension;i++)
			for (int j=0;j<dimension;j++) {
			ex e=manifold->d(components[i][j]);
			LOG_DEBUG(e);
			for (int k=0;k<dimension;k++)
				e+=components[i][k]*components[k][j];
			LOG_DEBUG(e);
			m(i,j)=e;
			}
	}
	catch (const Manifold::dException&)
	{
		throw WedgeException<logic_error>("Manifold does not know how to take d of its forms",__FILE__,__LINE__);
	}
	return m;
}

void Connection::DeclareConditions(const lst& list_of_equations)
{
	for (int i=0;i<components.size();++i)
		for (int j=0;j<components.size();++j)
			components[i][j]=components[i][j].subs(list_of_equations).normal();	
}

//////////////////////////////////////////////////////////////////////
// 		         CovariantDerivative
//////////////////////////////////////////////////////////////////////

//covariant derivative
ex internal::CovariantDerivative<DifferentialForm>::Apply(const VectorField& X, const VectorField& alpha) const
{
	const Frame& frame=connection.e();
	exvector alpha_j=frame.Components(alpha);
	LOG_DEBUG(alpha);
	LOG_DEBUG(alpha_j);
	ex res;
	for (int j=0;j<alpha_j.size();++j)
		if (!alpha_j[j].is_zero())
		{
			for (int i=0;i<frame.size();++i)
				res-=alpha_j[j]*TrivialPairing<VectorField>(X,connection(j,i))*frame[i];
			res+=connection.manifold->LieDerivative(X,alpha_j[j])*frame[j];
		}
	LOG_DEBUG(res);
	return res;
}

ex internal::CovariantDerivative<DifferentialForm>::Apply(const VectorField& vfield, const Function& f) const
{
	return connection.manifold->LieDerivative(vfield,f);	
}

//covariant derivative
ex internal::CovariantDerivative<VectorField>::Apply(const VectorField& X, const VectorField& alpha) const
{
	const Frame& frame=connection.e();
	//compute the components \alpha_j by evaluating on the dual basis, i.e. \alpha_j=e^j(\alpha)
	exvector alpha_j, dalpha_j; alpha_j.reserve(frame.size());
	for (int i=0;i<frame.size();++i)
		alpha_j.push_back(TrivialPairing<VectorField>(frame[i],alpha));

	LOG_DEBUG(alpha);
	LOG_DEBUG(alpha_j);
	ex res=0;
	for (int j=0;j<alpha_j.size();++j)
		if (!alpha_j[j].is_zero())
		{
			for (int i=0;i<frame.size();++i)
				res+=alpha_j[j]*TrivialPairing<VectorField>(X,connection(i,j))*frame.dual()[i];
			res+=Hook(X,connection.manifold->d(alpha_j[j]))*frame.dual()[j];
		}
	LOG_DEBUG(res);
	return res;
}

ex internal::CovariantDerivative<VectorField>::Apply(const VectorField& X, const Function& f) const
{
	return connection.manifold->LieDerivative(X,f);	
}

ex internal::CovariantDerivative<Spinor>::Apply(const VectorField& X, const Spinor& psi) const
{
	ex res; 
	for (int i=0;i<connection.e().size();++i)
		for (int j=i+1;j<connection.e().size();++j)
			res+=connection.g.CliffordDot(TrivialPairing<VectorField>(X,connection(i,j))*connection.e()[i]*connection.e()[j],psi);
	return -res/2;
}

ex internal::CovariantDerivative<Spinor>::Apply(const VectorField& vfield, const Function& f) const
{
	return connection.manifold->LieDerivative(vfield,f);	
}

////////////////////////////////////////////////////////////
//		RiemannianConnection						  
////////////////////////////////////////////////////////////

RiemannianConnection::RiemannianConnection(const Manifold* manifold,const RiemannianStructure& _g, bool donotinitialize) : Connection(manifold,_g.e(),true), g(manifold,_g.e()) 
{
	const int dimension=e().size();
	components.reserve(dimension);
	for (int i=0;i<dimension;i++)
		components.push_back(exvector(dimension));
}

RiemannianConnection::RiemannianConnection(const Manifold* manifold,const RiemannianStructure& _g, const Name& christoffel) : Connection(manifold,_g.e(),true), g(manifold,_g.e())
{
	const int dimension=e().size();
	components.reserve(dimension);
	for (int i=0;i<dimension;i++)
	{
		components.push_back(exvector(dimension));
		for (int j=0;j<i;j++)
			components[i][j]=-components[j][i];	
		components[i][i]=0;
		for (int j=i+1;j<dimension;j++)
		{
			ex omegaij=0;
			//here we would want to use dimension for a submersion, but manifold->Dimension() for a connection on a vector bundle
			for (int k=0;k<dimension;k++)
				omegaij+=Parameter(christoffel(k+1,j+1,i+1))*manifold->e()[k];
			components[i][j]=omegaij;
		}
	}
}


////////////////////////////////////////////////////////////
//		TorsionFreeConnection						  
////////////////////////////////////////////////////////////
template<> LeviCivitaConnection<false>::LeviCivitaConnection(const Manifold* manifold, const RiemannianStructure& g,const Name& christoffel) :
	 Connection(manifold,g.e(),true), 
	 RiemannianConnection(manifold,g,christoffel), 
	 TorsionFreeConnection<false>(manifold,g.e(),christoffel)
{
}

template<> LeviCivitaConnection<true>::LeviCivitaConnection(const Manifold* manifold, const RiemannianStructure& g,const Name& christoffel) :
	 Connection(manifold,g.e(),true), 
	 RiemannianConnection(manifold,g,true), 
	 TorsionFreeConnection<true>(manifold,g.e(),true,christoffel)
{
	const int dimension=e().size();
		 for (int i=0;i<dimension;++i)
			 	for (int j=0;j<dimension;++j)
				 		(*this)(i,j)=0;
	exvector de;
	de.reserve(dimension);
	for (int i=0;i<dimension;++i)
		de.push_back(manifold->d(e()[i]));
	for (int i=0;i<dimension;++i)
	for (int j=0;j<dimension;++j)
	for (int k=j+1;k<dimension;++k)
	{
		ex X=e().dual()[i];
		ex Y=e().dual()[j];
		ex Z=e().dual()[k];
		ex XYZ=Hook(X*Y,-de[k])+
			Hook(Z*X,-de[j])+
			Hook(Z*Y,-de[i]);
		(*this)(k,j)+=(e()[i]*XYZ/2).expand();
		LOG_DEBUG((*this)(k,j));
	}
	for (int j=0;j<dimension;++j)
	for (int k=j+1;k<dimension;++k)
		(*this)(j,k)=-(*this)(k,j);
}


void TorsionFreeConnection<false>::Declare_d(ex x, ex dx)
{
	x=x.expand();	
	if (x.is_zero()) return;
	ex generic_d_x=d(x);
	
	ex conds[2]={generic_d_x-dx,d(dx)};
	DeclareZero(conds,conds+2);
	assert((d(x)-dx).expand().normal().is_zero());
}

//use the fact that connection is torsion free to compute d from Nabla
ex TorsionFreeConnection<false>::d(ex alpha) const
{
	ex result;
	for (int i=0;i<e().size();i++)
		result+=e()[i]*Nabla<DifferentialForm>(e().dual()[i],alpha);
	return result.expand();
}

matrix TorsionFreeConnection<false>::CurvatureForm() const
{
	const int dimension=e().size();
	assert(dimension>0);  
	matrix m(dimension,dimension);
	for (int i=0;i<dimension;i++)
		for (int j=0;j<dimension;j++) {
		ex e=d(components[i][j]);		
		for (int k=0;k<dimension;k++)
			e+=components[i][k]*components[k][j];
		m(i,j)=e.expand();
		}
	return m;
}

}
