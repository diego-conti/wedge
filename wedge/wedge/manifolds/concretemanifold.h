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
 
#ifndef CONCRETEMANIFOLD_H_
#define CONCRETEMANIFOLD_H_
#include "wedge/manifolds/manifold.h"

/** @ingroup Manifolds
 *  @{
 * 
 *  @file concretemanifold.h
 *  @brief Manifolds that know how to take d of their forms 
*/
namespace Wedge {

/** @brief Abstract base class implementing a dTable, i.e. a table that computes the action of on the standard basis of one-forms.
 * 
 * A Has_dTable contains a \em dTable, i.e. a map that associates to each simple form on the manifold
 * its exterior derivative. 
 *
 * @note This is not the only class to reimplement d(). @sa ManifoldWith.
 */
class Has_dTable : public virtual Manifold {
public:
	Has_dTable(const Has_dTable& o); 		///< Copy constructor
	Has_dTable();
	const exmap& dTable() const {return table;}	///< Return a read-only version of the dTable
	ex d(ex alpha) const;
	bool KnowsHowToCompute_d() const {return true;} 
	ex Apply(const VectorField& X,const Function& f) const;
protected:
/** @brief Declare that \f$d\alpha=\beta\f$
 *  @param alpha A simple differential form or function on this manifold
 *  @param beta A differential form on this manifold
 * 
 *  Updates the table defining the operator d so that \f$d\alpha=\beta\f$
 * 
 *  @note If the standard frame of the manifold does not consist of simple elements, the action of d on simple elements can be
 *  recovered using LinearMapToSubstitutions()
*/
	void Declare_d(ex alpha, ex beta) {assert(is_a<DifferentialOneForm>(alpha) || is_a<Function>(alpha)); table[alpha]=beta;}
private:
	unique_ptr<DerivationOver<DifferentialForm,Function> > the_d_operator; ///< Pointer to visitor class used to compute the action of d  on forms on this manifold. The type is really dOperator, which is however only defined in manifold.cpp	
	exmap table; ///< Table describing the action of d. Thus, table[e(i)] represents d(e(i))
	class dOperator; ///< Helper visitor class for the operator d
};


/**  @brief A specialization of Manifold that automatically declares a basis of 1-forms
 * 
 * @remark ConcreteManifold and its derived classes, unlike the base class Manifold, assume that the frame consists of simple elements.  
*/
class ConcreteManifold : public virtual Manifold {
public:
/**	@brief Create a ConcreteManifold
 * 	@param dimension The dimension of the manifold
 *  To construct an actual concrete manifold, subclass from ConcreteManifold and initialize it in the constructor
 *  by invoking Declare_d.
*/
	ConcreteManifold(int dimension);

/** @brief Create a ConcreteManifold with a given frame
 * 
 * For use by subclasses that define the frame themselves
*/	
 	ConcreteManifold(exvector frame);
	
/**
  * @brief Returns a frame for this manifold
  * @return A Frame object, consisting of a basis of one-forms
  * 
  * The elements of the frame are simple.
*/
	const Frame& e() const {return frame;}

	ex e(OneBased k) const {return frame(k);} // Redefined since masked by ConcreteManifold::e().
  
private:
	Frame frame;
	Frame CreateFrame(int dimension);	///<used internally in ctor
};

} /** @} */
#endif /*CONCRETEMANIFOLD_H_*/
