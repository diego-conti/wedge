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
 *******************************************************************************///hack to make doxygen document this file
#ifdef DOXYGEN_RUNNING
#define LAMBDA_H_
#define TENSOR_H_
#endif 

#ifdef LAMBDA_H_
#ifdef TENSOR_H_
#ifndef TENSORLAMBDA_H_
#define TENSORLAMBDA_H_

/** @ingroup LinearAlgebra */ 

/** @{ 
 * @file tensorlambda.h
 * @brief Template specializations for tensor products of exterior algebras
 */
 
namespace Wedge {

namespace internal {

template<typename V,typename W> class TensorProductOperator<Lambda<V>,W> : public IBilinearOperator<LinearOperator<Lambda<V> >,LinearOperator<W> > {
public:
	ex Apply(const V& v, const W& w) const {
//		return Tensor<Lambda<V>,W>(exvector(&v,(&v)+1),w);
		return Tensor<Lambda<V>,W>(v,w);		
	}
	ex Apply(const Lambda<V>& v, const W& w) const {
		return Tensor<Lambda<V>,W>(v,w);		
	}
};

template<typename V,typename W> class TensorProductOperator<V,Lambda<W> > : public IBilinearOperator<LinearOperator<V>,LinearOperator<Lambda<W> > >  {
public:
	ex Apply(const V& v, const W& w) const {
//		return Tensor<V,Lambda<W> >(v,exvector(&w,(&w)+1));
		return Tensor<V,Lambda<W> >(v,w);		
	}
	ex Apply(const V& v,const Lambda<W>& w) const {
		return Tensor<V,Lambda<W> >(v,w);		
	}
};

template<typename V,typename W> class TensorProductOperator<Lambda<V>,Lambda<W> > : public IBilinearOperator<LinearOperator<Lambda<V> >,LinearOperator<Lambda<W> > > {
public:
	ex Apply(const V& v, const W& w) const {
//		return Tensor<Lambda<V>,Lambda<W> >(exvector(&v,(&v)+1),exvector(&w,(&w)+1));
		return Tensor<Lambda<V>,Lambda<W> >(v,w);		
	}
	ex Apply(const V& v,const Lambda<W>& w) const {
//		return Tensor<Lambda<V>,Lambda<W> >(exvector(&v,(&v)+1),w);
		return Tensor<Lambda<V>,Lambda<W> >(v,w);		
	}
	ex Apply(const Lambda<V>& v, const W& w) const {
//		return Tensor<Lambda<V>,Lambda<W> >(v,exvector(&w,(&w)+1));
		return Tensor<Lambda<V>,Lambda<W> >(v,w);		
	}
	ex Apply(const Lambda<V>& v, const Lambda<W>& w) const {
		return Tensor<Lambda<V>,Lambda<W> >(v,w);		

	}
};

}

}

/** @} */

#endif /*TENSORLAMBDA_H_*/
#endif /*TENSOR_H_*/
#endif /*LAMBDA_H_*/
