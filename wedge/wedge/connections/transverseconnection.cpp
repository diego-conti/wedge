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
#include "wedge/connections/transverseconnection.h"

namespace Wedge 
{

ex TransverseLeviCivitaConnection::T(ex X, ex Y) const
{
	ex XV=V<VectorField>(X);
	return V<VectorField>(Nabla<VectorField>(XV,H<VectorField>(Y)))+H<VectorField>(Nabla<VectorField>(XV,V<VectorField>(Y)));
}

ex TransverseLeviCivitaConnection::A(ex X, ex Y) const
{
//	return V<VectorField>(LieBracket(X,Y))/2;	//this is only ok for horizontal vector fields
	ex XH=H<VectorField>(X);
	ex AXY=H<VectorField>(Nabla<VectorField>(XH,V<VectorField>(Y)))+V<VectorField>(Nabla<VectorField>(XH,H<VectorField>(Y)));
	return AXY;
}



matrix TransverseLeviCivitaConnection::BaseRicciAsMatrix(bool useRicciFormula,const Simplifier& simplifier) const
{
	matrix Ric(g.BaseDimension(),g.BaseDimension());
	if (useRicciFormula) {
		matrix R=RicciAsMatrix();
		LOG_INFO(R);
		ex N;
		SubBasis<VectorField> frame(e().dual().begin(),e().dual().begin()+g.BaseDimension(),e().dual().begin()+g.BaseDimension(),e().dual().end());
		LOG_INFO(frame);
		for (exvector::const_iterator i=frame.complement_begin();i!=frame.complement_end();++i)
			N+=T(*i,*i);
		LOG_INFO(N);	
		matrix Aij(frame.size(),frame.size());
		for (int h=1;h<=frame.size();++h)
		for (int k=1;k<=frame.size();++k)
			Aij(h-1,k-1)=A(frame(h),frame(k));
		Aij=simplifier.Simplify(Aij);
		for (int h=1;h<=frame.size();++h)
		for (int k=1;k<=frame.size();++k)
		{			
			ex X=frame(h),Y=frame(k);
			ex AXAY,TXTY;			
			for (int i=1;i<=frame.size();++i)
				AXAY+=g.ScalarProduct().OnVectors(Aij(h-1,i-1),Aij(k-1,i-1));
			for (exvector::const_iterator i=frame.complement_begin();i!=frame.complement_end();++i)
				TXTY+=g.ScalarProduct().OnVectors(T(*i,X),T(*i,Y));
			Ric(h-1,k-1)=R(h-1,k-1)+2*AXAY+TXTY-
				(g.ScalarProduct().OnVectors(Nabla<VectorField>(X,N),Y)+
						g.ScalarProduct().OnVectors(Nabla<VectorField>(Y,N),X))/2;
		}
		LOG_INFO(Ric);
		return ex_to<matrix>(Ric.expand());
	}
	else {
		matrix R=BaseCurvatureForm(simplifier);
		for (int j=0;j<g.BaseDimension();++j) {
			ex sum;
			for (int i=0;i<g.BaseDimension();++i)
				sum+=Wedge::Hook(e().dual()[i],R(i,j));
			for (int k=0;k<g.BaseDimension();++k)
				Ric(k,j)=Wedge::Hook(e().dual()[k],sum);
		}
	}
	return simplifier.Simplify(Ric);
}


matrix TransverseLeviCivitaConnection::BaseCurvatureForm(const Simplifier& simplifier) const
{
	matrix A(g.BaseDimension(),g.BaseDimension());
	ExVector frame=e().dual();
	for (int h=1;h<=g.BaseDimension();++h)
	for (int k=1;k<=g.BaseDimension();++k)
		A(h-1,k-1)=simplifier.Simplify(this->A(frame(h),frame(k)));	//symmetric?

	matrix R1=CurvatureForm();
	LOG_DEBUG(R1);
	matrix R(g.BaseDimension(),g.BaseDimension());
	for (int i=0;i<g.BaseDimension();++i)
	for (int j=0;j<g.BaseDimension();++j)
		R(i,j)=simplifier.Simplify(H<DifferentialForm>(R1(i,j)));
	LOG_DEBUG(R);
	for (int h=1;h<=g.BaseDimension();++h)
	for (int k=h+1;k<=g.BaseDimension();++k)
	{
		ex X=frame(h), Y=frame(k);
		for (int i=1;i<=g.BaseDimension();++i)
		for (int j=1;j<=g.BaseDimension();++j)
		{
			ex deltaXY=-2*g.ScalarProduct().OnVectors(A(h-1,k-1), A(i-1,j-1))
			+g.ScalarProduct().OnVectors(A(k-1,i-1), A(h-1,j-1))
			-g.ScalarProduct().OnVectors(A(h-1,i-1), A(k-1,j-1));
			LOG_DEBUG(deltaXY);
			R(i-1,j-1)-=deltaXY*e(h)*e(k);
		}
	}
	LOG_DEBUG(R);
	return R;
}


}

