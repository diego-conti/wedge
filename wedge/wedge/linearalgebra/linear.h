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
#ifndef LINEAR_H
#define LINEAR_H
/** @defgroup LinearAlgebra Linear and polynomial algebra: vector spaces, (multi)linear operators, polynomial ideals */
/** @{ 
 * @file linear.h
 * @brief Linear operators on vector spaces
 */
#include "wedge/base/utilities.h"

namespace Wedge {
 using namespace  GiNaC;

// ****************************************************************************
// *							AdditiveOperator						  *
// ****************************************************************************
/** @brief Abstract base class for an additive operator on a space T
 * @param T The type of simple elements in T
 * 
 * visit(const basic&) must be implemented in subclass.
 * 
 * @warning An additive operator, by definition, satisfies \f[\phi(-a)=-\phi(a),\quad \phi(n\,a)=n\phi(a).\f] This must be implemented in the subclass (but see LinearOperator)
 */ 
template<typename T> class AdditiveOperator: public RecursiveVisitor<ex>, public basic::visitor ,public  add::visitor, public T::visitor {
protected:
/** @brief The ``basic'' type operated on by this operator
 * 
 * T can be a ``composite'' type, e.g. Lambda<V>. In this case, OperatesOn equals V. This 
 * is motivated by the fact that Lambda<V> contains V anyway, and the action on Lambda<V> may
 * be determined by the action on V in certain cases.
 */ 
	typedef T OperatesOn;

public:
	enum {ExpectsExpandedInput=0};

// Overloaded for efficiency 
	ex RecursiveVisit(const GiNaC::ex& e)	 
	{	
		ex oldresult=GetResult();
		if (is_exactly_a<basic>(e))
			visit(ex_to<basic>(e));		
		else if (is_exactly_a<add>(e))
			visit(ex_to<add>(e)); 
		else if (is_exactly_a<T>(e))
			static_cast<typename T::visitor*>(this)->visit(ex_to<T>(e));
		else
			e.accept(*this);		
		ex newresult=GetResult();
		Result()=oldresult;
		return newresult;		
	}

private: 

/** @brief Implements additivity
 */
	void visit(const add& alpha)
	{
		Result()=0;
		for (unsigned i=0;i<alpha.nops();i++)
			Result()+=this->RecursiveVisit(alpha.op(i));
	}	

/** @brief Called internally; applies operator to an object of type basic.
 * 
 * Strictly speaking, this should only occur in the case of an element of the base field in an algebra with 1.
 * However, this function is invoked internally with the coefficients as arguments.
 */
	void visit(const basic&) {Result()=0;}
};

 
// ****************************************************************************
// *							LinearOperator						  *
// ****************************************************************************
/** @brief Abstract base class for a linear operator on a vector space or algebra T
 *  @param T The type of simple elements in T
 */
template<typename T> class LinearOperator: public AdditiveOperator<T>, public mul::visitor
{
public:
	ex RecursiveVisit(const GiNaC::ex& e) 
	{	
		if (is_exactly_a<mul>(e))
		{
			ex oldresult=this->GetResult();
			visit(ex_to<mul>(e));
			ex newresult=this->GetResult();
			this->Result()=oldresult;
			return newresult;
		}
		else return AdditiveOperator<T>::RecursiveVisit(e);
	}

private:

/** @brief Implements linearity
 */
	void visit(const mul& alpha)
	{		
		exvector v;
		v.reserve(alpha.nops());
		bool vectorFound=false;
		for (unsigned i=0;i<alpha.nops();i++)
		{
			const ex& factor=alpha.op(i);
			ex image=RecursiveVisit(factor);
			if (image!=0) {
				if (vectorFound)
				{
					LOG_ERROR(alpha);
					LOG_ERROR(image);
					throw WedgeException<std::logic_error>("Linear expression expected instead of polynomial",__FILE__,__LINE__);
				}
				v.push_back(image);
				vectorFound=true;	
			}
			else
				v.push_back(factor);
		}
		if (vectorFound) 
			this->Result()=mul(v);
		else
			this->Result()=0;
	}
};

}  /** @} */
#endif
