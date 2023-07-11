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
 #ifndef STRUCTURES_H_
#define STRUCTURES_H_

/** @ingroup RiemannianGeometry */ 

/** @{ 
 * @file structures.h 
 * @brief Concrete examples of G-structures
 * 
 * @note Classes derived from GStructure define the dimension of the corresponding manifold as an enum
 */

#include "wedge/structures/riemannianstructure.h"
#include "wedge/connections/torsionfreeconnection.h"
#include "wedge/manifolds/concretemanifold.h"

namespace Wedge {

/** @brief An \f$\rm SU(3)\f$-structure on a manifold of real dimension 6 (e.g. a Calabi-Yau 3-fold)
 * 
 * @sa [Chiossi-Salamon, The intrinsic torsion of \f$\rm SU(3)\f$ and \f$G\sb 2\f$ structures.  Differential geometry, Valencia, 2001,  115--133, World Sci. Publ., River Edge, NJ, 2002]
*/
class SU3Structure : public RiemannianStructure {
public:
	static IntrinsicTorsionClass W1Plus, W1Minus,W2Plus, W2Minus, W3,W4,W5;	// The Gray-Hervella torsion classes
	enum {dimension=6};	///< The dimension of the manifold
/** @brief Define an \f$\rm SU(3)\f$-structure on a 6-manifold, represented by the choice of a frame
 *  @param manifold The manifold on which the \f$\rm SU(3)\f$-structure is defined.
 *  @param frame The adapted frame defining the structure.
 * 
 *  @warning Caller must ensure that the pointer remains valid. 
 */
	SU3Structure(const Manifold* manifold, const Frame& frame)	: RiemannianStructure(manifold, frame) {}
	ex omega() const;	///< The almost-symplectic two-form
	ex psiplus() const; ///< The real part of the complex volume form
	ex psiminus() const;///< The imaginary part of the complex volume form
};

template<> IntrinsicTorsion GStructureHasParameters<SU3Structure,false>::GetIntrinsicTorsion() const;

/** @brief A \f$G_2\f$-structure on a manifold of real dimension 7
 * 
 * @sa [Fernández-Gray, Riemannian manifolds with structure group \f$G\sb{2}\f$ Ann. Mat. Pura Appl. (4) 132 (1982), 19--45 (1983)]  
*/

class G2Structure : public RiemannianStructure {
public:
	enum {dimension=7};	///<The dimension of the manifold
/** @brief Define a \f$G_2\f$-structure on a 7-manifold, represented by the choice of a frame
 *  @param manifold The manifold on which the \f$\rm SU(3)\f$-structure is defined.
 *  @param frame The adapted frame defining the structure.
 * 
 *  @warning Caller must ensure that the pointer remains valid. 
 */
	G2Structure(const Manifold* manifold, const Frame& frame)
		: RiemannianStructure(manifold, frame) {}
	ex phi() const;	///< The defining 3-form
	ex starphi() const; ///< The Hodge dual of phi
};

/** @brief An \f$\rm SU(2)\f$-structure on a 5-manifold
 * 
 * @sa [Conti-Salamon, Generalized Killing spinors in dimension 5, Trans. Amer. Math. Soc.  359 (2007), N. 11, 5319-5343] 
 */
class SU2Structure : public RiemannianStructure {
public:
	enum {dimension=5};	///< The dimension of the manifold
/** @brief Define a \f$\rm SU(2)\f$-structure on a 5-manifold, represented by the choice of a frame
 *  @param manifold The manifold on which the \f$\rm SU(2)\f$-structure is defined.
 *  @param frame The adapted frame defining the structure.
 * 
 *  @warning Caller must ensure that the pointer remains valid.
 */ 		
	SU2Structure(const Manifold* manifold, const Frame& frame)
		: RiemannianStructure(manifold, frame) {}

	ex omega1() const;	///< The fundamental 2-form
	ex omega2() const;	///< The real part of the complex volume on the almost-contact distribution
	ex omega3() const;	///< The imaginary part of the complex volume on the almost-contact distribution
	ex alpha() const;	///< The (almost) contact form
	ex psi() const; 	///< The spinor corresponding to the reduction \f$\mathrm{SO}(5)\supset \mathrm{SU}(2)\f$
};

/** @brief An \f$\rm SU(3)\f$-structure on a manifold of real dimension 7
 *  
 * @sa [Conti-Fino, Calabi-Yau cones from contact reduction, arXiv:0710.4441]
*/

class SU3StructureDim7 : public RiemannianStructure {
public:
	enum {dimension=7};///< The dimension of the manifold
		
/** @brief Define an \f$\rm SU(3)\f$-structure on a 7-manifold, represented by the choice of a frame
 *  @param manifold The manifold on which the \f$\rm SU(3)\f$-structure is defined.
 *  @param frame The adapted frame defining the structure.
 * 
 *  @warning Caller must ensure that the pointer remains valid.
 */ 		
	SU3StructureDim7(const Manifold* manifold, const Frame& frame)
: RiemannianStructure(manifold, frame) {}

	ex alpha() const;	///< The (almost)contact form
	ex F() const;	///< The fundamental 2-form
	ex OmegaPlus() const;	///< The real part of the complex volume on the almost-contact distribution
	ex OmegaMinus() const;///< The imaginary part of the complex volume on the almost-contact distribution
	ex Omega() const {return OmegaPlus()+I*OmegaMinus();}	///< The complex volume on the almost-contact distribution
	ex psi() const; 	///<The spinor corresponding to the reduction \f$\mathrm{SO}(7)\supset \mathrm{SU}(3)\f$
};

/** @brief A \f$\mathrm{PSU}(3)\f$-structure on a manifold of real dimension 8
 * 
 * \f$\mathrm{PSU}(3)\f$ stands for \f$\mathrm{SU}(3)/\mathbb{Z}_3\f$, acting on \f$\mathbb{R}^8=\mathfrak{su}(3)\f$ by the 
 * adjoint representation.
 * @sa [N. Hitchin. Stable forms and special metrics. In Global Differential Geometry: The Mathematical Legacy of Alfred Gray, volume 288 of Contemp. Math., pages 70-89. American Math. Soc., 2001]
*/

class PSU3Structure : public RiemannianStructure {
public:
	static IntrinsicTorsionClass Xi8Plus, Xi8Minus, Xi20, Xi27Plus, Xi27Minus, Xi14;
	enum LambdaComponentType {
		Lambda_2_8,		///< The 8-dimensional component in the \f$\mathrm{PSU}(3)\f$-module of 2-forms
		Lambda_2_20,	///< The 20-dimensional component in the \f$\mathrm{PSU}(3)\f$-module of 2-forms
		Lambda_3_8,		///< The 3-dimensional component in the \f$\mathrm{PSU}(3)\f$-module of 6-forms
		Lambda_3_20,	///< The 20-dimensional component in the \f$\mathrm{PSU}(3)\f$-module of 3-forms
	 	Lambda_4_8Plus,	///< The 8-dimensional component in the \f$\mathrm{PSU}(3)\f$-module of self-dual 4-forms
	 	Lambda_4_8Minus,///< The 8-dimensional component in the \f$\mathrm{PSU}(3)\f$-module of anti-self-dual 4-forms
	 	Lambda_4_27Plus,	///< The 27-dimensional component in the \f$\mathrm{PSU}(3)\f$-module of self-dual 4-forms
	 	Lambda_4_27Minus,	///< The 27-dimensional component in the \f$\mathrm{PSU}(3)\f$-module of anti-self-dual 4-forms
		Lambda_6_20,	///< The 20-dimensional component in the \f$\mathrm{PSU}(3)\f$-module of 6-forms
		Lambda_6_8,		///< The 8-dimensional component in the \f$\mathrm{PSU}(3)\f$-module of 6-forms			
	 	nLambdaComponents
	};
	enum {dimension=8};///< The dimension of the manifold
	
/** @brief Define a \f$\mathrm{PSU}(3)\f$-structure on an 8-manifold, represented by the choice of a frame
 *  @param manifold The manifold on which the \f$\mathrm{PSU}(3)\f$-structure is defined.
 *  @param frame The adapted frame defining the structure.
 * 
 *  @warning Caller must ensure that the pointer remains valid.
 */ 	
	PSU3Structure(const Manifold* manifold, const Frame& frame)
		: RiemannianStructure(manifold, frame) {}

	ex phi() const;			///< The defining 3-form
	ex starphi() const;		///< The Hodge-star of the defining 3-form
/** @brief Return one of the irreducible \f$\mathrm{PSU(3)}\f$-modules in the algebra of differential forms
 * @param type Identifies the component to return
 * @return A VectorSpace of differential forms on the manifold, corresponding to the requested \f$\mathrm{PSU(3)}\f$-module 
 * 
 * @remark Since the \f$\mathrm{PSU}(3)\f$-structure admits a global frame, the decomposition of \f$\Lambda(\mathbb{R}^8)^*\f$ into
 * irreducible \f$\mathrm{PSU}(3)\f$-modules carries over to a decomposition of the vector bundle \f$\Lambda(M)\f$
 */ 
	VectorSpace<DifferentialForm> LambdaComponent(LambdaComponentType type) const;

/** @brief Compute the metric tensor from a form with stabilizer \f$\mathrm{PSU}(3)\f$
 *
 * This static member computes a representative metric in the conformal class determined by a 3-form with stabilizer \f$\mathrm{PSU}(3)\f$
 *  
 * @remark G-structures are represented in %Wedge by the choice of a frame. However, a \f$\mathrm{PSU}(3)\f$-
 * structure is determined by the choice of a three-form with orbit type \f$\left(\frac{GL(8,\mathbb{R})}{\mathrm{PSU}(3)}\right)\f$.
*/ 
	static matrix Metric(const Manifold& M,ex phi);
protected:		
	bool IsGamma30Zero() const;
public:
	ex PiPlus(ex fourform) const; ///< The epimorphism  \f$\pi_+: \Lambda^4 \to \Lambda^1,\quad \pi_+(\alpha)=*([\alpha]_+\wedge\phi )\f$
	ex PiMinus(ex fourform) const;///< The epimorphism  \f$\pi_-: \Lambda^4 \to \Lambda^1,\quad \pi_-(\alpha)=*([\alpha]_-\wedge\phi )\f$
	ex PiZero(ex sixform) const;	///< The epimorphism \f$\pi_0: \Lambda^6 \to \Lambda^1,\quad \pi_0(\alpha)=*\phi\lrcorner\alpha\f$
	ex SigmaPlus(ex oneform) const;   ///< The monomorphism \f$\sigma_+: \Lambda^1 \to \Lambda^4_{8+},\quad \sigma_+(\alpha)=-\frac45[\alpha\wedge\phi]_+\f$
	ex SigmaMinus(ex oneform) const;  ///< The monomorphism \f$\sigma_-: \Lambda^1 \to \Lambda^4_{8-},\quad \sigma_-(\alpha)=\frac45[\alpha\wedge\phi]_-\f$
	ex SigmaZero(ex oneform) const;	///< The isomorphism \f$\sigma_0: \Lambda^1 \to \Lambda^6_{8},\quad \sigma_0(\alpha)=\frac23 *\phi\wedge\alpha\f$
};

template<> IntrinsicTorsion GStructureHasParameters<PSU3Structure,false>::GetIntrinsicTorsion() const;

/** @brief Output operator
 */ 
template<class charT, class traits> std::basic_ostream<charT,traits>& operator<<(std::basic_ostream<charT,traits>& os, const SU3Structure& gStructure)
{
	return os<<"omega="<<gStructure.omega()<<";"<<endl<<"psiplus="<<gStructure.psiplus()<<endl;
}

/** @brief Output operator
 */ 
template<class charT, class traits> std::basic_ostream<charT,traits>& operator<<(std::basic_ostream<charT,traits>& os, const G2Structure& gStructure)
{
	return os<<"phi="<<gStructure.phi()<<";"<<endl<<"*phi="<<gStructure.starphi()<<endl;
}

/** @brief Output operator
 */ 
template<class charT, class traits> std::basic_ostream<charT,traits>& operator<<(std::basic_ostream<charT,traits>& os, const SU2Structure& gStructure)
{
	os<<"alpha="<<gStructure.alpha()<<";"<<endl<<"omega1="<<gStructure.omega1()<<endl;
	return os<<"omega2="<<gStructure.omega2()<<endl;os<<"omega3="<<gStructure.omega3()<<endl;
}

/** @brief Output operator
 */ 
template<class charT, class traits> std::basic_ostream<charT,traits>& operator<<(std::basic_ostream<charT,traits>& out, const SU3StructureDim7& gStructure)
{
	return out<<"alpha="<<gStructure.alpha()<<";"<<endl<<"F="<<gStructure.F()<<";"<<endl<<"Omega="<<gStructure.OmegaPlus()+I*gStructure.OmegaMinus()<<endl;
}



////////////////////////////////////////////////////////////////////////////////////////////
//						Implementation of GetEquationsForTorsionIn
////////////////////////////////////////////////////////////////////////////////////////////

template<> template<typename Container> Container& GStructureHasParameters<SU3Structure,true>::GetEquationsForTorsionIn(Container& V, ex TorsionClass) const
{
	set<IntrinsicTorsionClass,basic_is_less> classes;	
	
	ex domega=M()->d(omega());
	ex dpsiplus=M()->d(psiplus());
	ex dpsiminus=M()->d(psiminus());
	
	ex zero_components;
	if (IntrinsicTorsion::AreTypesComplementary(TorsionClass,W1Plus+W1Minus+W3+W4))
	{
		GetCoefficients<DifferentialForm>(V,domega);
		zero_components+=W1Plus+W1Minus+W3+W4;
	}
	else if (IntrinsicTorsion::AreTypesComplementary(TorsionClass,W4))
	{
		GetCoefficients<DifferentialForm>(V,domega*omega());
		zero_components+=W4;
	}
	
	if (IntrinsicTorsion::AreTypesComplementary(TorsionClass,W1Plus+W2Plus+W5))
	{		
		GetCoefficients<DifferentialForm>(V,dpsiplus);
		zero_components+=W1Plus+W2Plus+W5;
	}
	if (IntrinsicTorsion::AreTypesComplementary(TorsionClass,W1Minus+W2Minus+W5))
	{
		GetCoefficients<DifferentialForm>(V,dpsiminus);
		zero_components+=W1Minus+W2Minus+W5;
	}
	
	if (IntrinsicTorsion::IsNonNegative(TorsionClass+zero_components-W1Plus+W1Minus+W2Plus+W2Minus+W3+W4+W5))
		return V;
	else throw NotImplemented(__FILE__,__LINE__);	
}

template<> template<typename Container> Container& GStructureHasParameters<PSU3Structure,true>::GetEquationsForTorsionIn(Container& V, ex TorsionClass) const
{
	if (IntrinsicTorsion::AreTypesComplementary(TorsionClass,Xi14))
	throw NotImplemented(__FILE__,__LINE__);	
	
	set<IntrinsicTorsionClass,basic_is_less> classes;	
	
	ex dphi=M()->d(phi());
	ex dstarphi=M()->d(starphi());
	VectorSpace<DifferentialForm> L_2_20=LambdaComponent(Lambda_2_20);
	if (IntrinsicTorsion::AreTypesComplementary(TorsionClass,Xi20))		
		for (IBasis<DifferentialForm>::const_iterator i=L_2_20.e_begin();i!=L_2_20.e_end();++i)
			GetCoefficients<DifferentialForm>(V,dstarphi * *i);
			
	
	if (IntrinsicTorsion::AreTypesComplementary(TorsionClass,Xi8Plus+Xi8Minus+Xi27Plus+Xi27Minus))
	{
		GetCoefficients<DifferentialForm>(V,dphi);		
	}
	else {
		exvector dphi_orth;	//basis of 4-forms orthogonal to dphi
		if (IntrinsicTorsion::AreTypesComplementary(TorsionClass,Xi8Plus))
		{	
			VectorSpace<DifferentialForm> V=LambdaComponent(Lambda_4_8Plus);
			dphi_orth.insert(dphi_orth.end(),V.e_begin(),V.e_end());
		}
		if (IntrinsicTorsion::AreTypesComplementary(TorsionClass,Xi8Minus))		
		{	
			VectorSpace<DifferentialForm> V=LambdaComponent(Lambda_4_8Minus);
			dphi_orth.insert(dphi_orth.end(),V.e_begin(),V.e_end());
		}
		if (IntrinsicTorsion::AreTypesComplementary(TorsionClass,Xi27Plus))		
		{	
			VectorSpace<DifferentialForm> V=LambdaComponent(Lambda_4_27Plus);
			dphi_orth.insert(dphi_orth.end(),V.e_begin(),V.e_end());
		}
		if (IntrinsicTorsion::AreTypesComplementary(TorsionClass,Xi27Minus))		
		{	
			VectorSpace<DifferentialForm> V=LambdaComponent(Lambda_4_27Minus);
			dphi_orth.insert(dphi_orth.end(),V.e_begin(),V.e_end());
		}
		for (ExVector::const_iterator i=dphi_orth.begin();i!=dphi_orth.end();++i)
			GetCoefficients<DifferentialForm>(V,dphi* *i);
	}
	return V;
}

} /** @} */

#include "wedge/liealgebras/liegroupstructures.h"
#endif /*STRUCTURES_H_*/
