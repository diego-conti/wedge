/*******************************************************************************
 *  Copyright (C) 2007-2022 by Diego Conti, diego.conti@unipi.it 
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
 
#ifndef LIESUBGROUP_H_
#define LIESUBGROUP_H_
#include "wedge/liealgebras/liegroup.h"
#include "wedge/manifolds/concretemanifold.h"

/** @ingroup Manifolds
 *  @{
 * 
 *  @file liesubgroup.h  
 *  @brief Lie subgroups 
*/

namespace Wedge {

/** @brief Instances of LieSubgroup represent connected subgroups of a LieGroup.
 * 
 * @remark A connected Lie subgroup of a given group is determined by its Lie algebra.
 * @warning This class represents the Lie algebra of the subgroup as a subspace of the containing group. To make this work,
 * LieSubgroup has to compute the orthogonal complement of the subalgebra in the Lie group, which is computationally expensive
 * for high dimensions. @sa AbstractLieSubgroup    
*/

template<bool WithParams> class LieSubgroup : public LieGroupHasParameters<WithParams>, public virtual Has_dTable {	
	Subspace<DifferentialOneForm> lieAlgebra;	///< Lie algebra of this subgroup, as a subalgebra of the Lie algebra of containingGroup
public:	
/** @brief Create an object representing the unique connected subgroup of a Lie group with the given Lie algebra.
 *  @param G A Lie group
 *  @param subalgebra A set of generators for a subalgebra \f$\mathfrak{h}\f$ of the Lie algebra \f$\mathfrak{g}\f$ of G.
 * 
 * @note The elements of subalgebra are interpreted as vector fields, when defining the subalgebra. However, they also determine
 * the standard basis of one-forms of the group. So, for instance, if G has a basis e1,e2,e3 and subalgebra is {e1+e2,e2+e3}, then
 * \f[\mathfrak{h}=\langle e_1+e_2,e_2+e_3\rangle,\quad \mathfrak{h}^*=\langle e^1+e^2,e^2+e^3\rangle,\f]
 * which is consistent because the complement of \f$\mathfrak{h}\f$ in \f$\mathfrak{g}\f$ is spanned by \f$e_1-e_2+e_3\f$. 
*/ 
	LieSubgroup(const LieGroupHasParameters<WithParams>& G, const exvector& subalgebra) : LieGroupHasParameters<WithParams>(G),
		lieAlgebra(subalgebra, G.e(), TrivialPairingOperator<DifferentialOneForm>())
	{			
			Subspace<DifferentialForm> Lambda2(Wedge::TwoForms(lieAlgebra.e()));	// Space of two-forms on the Lie algebra of this subgroup
			exvector simple;	//the simple one-forms appearing in e()			
			GetSymbols<DifferentialOneForm>(simple, e().begin(),e().end());
			for (exvector::const_iterator i=simple.begin();i!=simple.end();i++)				
				Has_dTable::Declare_d(*i,Lambda2.Project(G.d(*i)));
			Check_ddZeroIfPossible();			
	} 
	const Frame& e() const {return lieAlgebra.e();}
	ex e(int k) const {return e()(k);}
	
private:
	void Check_ddZeroIfPossible();
};


/** @brief Instances of AbstractLieSubgroup represent connected subgroups of a LieGroup.
 * 
 * This class is analogous to LieSubgroup, except that a new frame is defined, consisting of simple elements, rather than
 * using one-forms on the containing group. This is considerably more efficient for high dimensions.
 * 
 * @remark This class may also be used to associate a new basis of left-invariant forms to a given Lie group.   
*/

template<bool WithParams> class AbstractLieSubgroup : public LieGroupHasParameters<WithParams>, public ConcreteManifold , public virtual Has_dTable{	
public:	
/** @brief Create an object representing the unique connected subgroup of a Lie group with the given Lie algebra.
 *  @param G A Lie group
 *  @param subalgebra A set of generators for a subalgebra \f$\mathfrak{h}\f$ of the Lie algebra \f$\mathfrak{g}\f$ of G.
 *  @exception WedgeException<std::runtime_error> Thrown if \f$\mathfrak{h}\f$ is not a subalgebra. 
 * 
 * This constructor computes the structure constants of the given subalgebra, and uses them to define a new LieGroup object
 * @note If HasParameters=true, then \f$\mathfrak{h}\f$ is required to be a subalgebra for every choice of the parameters
*/ 
	AbstractLieSubgroup(const LieGroupHasParameters<WithParams>& G, const exvector& subalgebra) :
		ConcreteManifold(subalgebra.size())
	{
		Initialize(G,subalgebra);
	}
	/** @brief Create an object representing the unique connected subgroup of a Lie group with the given Lie algebra.
	 *  @param G A Lie group
	 *  @param subalgebra A set of generators for a subalgebra \f$\mathfrak{h}\f$ of the Lie algebra \f$\mathfrak{g}\f$ of G.
	 *  @param oneforms A vector of
	 *  @exception WedgeException<std::runtime_error> Thrown if \f$\mathfrak{h}\f$ is not a subalgebra.
	 *
	 * This constructor computes the structure constants of the given subalgebra, and uses them to define a new LieGroup object
	 * @note If HasParameters=true, then \f$\mathfrak{h}\f$ is required to be a subalgebra for every choice of the parameters
	*/

	AbstractLieSubgroup(const LieGroupHasParameters<WithParams>& G, const exvector& subalgebra, const ExVector oneforms) :
		ConcreteManifold(oneforms)
	{
		if (oneforms.size()!=subalgebra.size()) throw WedgeException<std::invalid_argument>("Subalgebra should have same number of elements as frame in AbstractLieSubgroup",__FILE__,__LINE__);
		Initialize(G,subalgebra);
	}


protected:
	void Initialize(const LieGroupHasParameters<WithParams>& G, const exvector& subalgebra)
	{		
		Frame frame(subalgebra);
		assert(Dimension()==frame.size());
		exvector dei(Dimension());
		for (int i=0;i<Dimension();i++)
		for (int j=i+1;j<Dimension();j++)
		{
			ex eij=G.LieBracket(subalgebra[i],subalgebra[j]);			
			try {
				exvector comps=frame.Components(eij);
				assert(Dimension()==comps.size());
				for (int k=0;k<Dimension();k++)
					dei[k]-=comps[k]*e()[i]*e()[j];
			}
			catch (const InvalidArgument&)
			{
				LOG_ERROR(i<<" "<<j<<" "<<eij);
				throw WedgeException<std::runtime_error>("Not a Lie subalgebra",__FILE__,__LINE__);	
			}
		}
		for (int i=0;i<Dimension();i++)
			Has_dTable::Declare_d(e()[i],dei[i]); 
	} 

};




} /** @} */
#endif /*LIESUBGROUP_H_*/
