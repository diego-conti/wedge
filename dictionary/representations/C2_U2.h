/*
 * C2_U2.h
 *
 *  Created on: Mar 27, 2015
 *      Author: diego
 */

#ifndef DICTIONARY_REPRESENTATIONS_C2_U2_H_
#define DICTIONARY_REPRESENTATIONS_C2_U2_H_

#include "../vvaluedforms.h"


/** @brief \f$\mathbb{C}^2\f$ as a representation of \f$U(2)\f$
*/
struct C2_U2
{
	static int FibreSize() {return 4;}	///<Real Dimension of \f$V=\mathbb{C}^2\f$

	/** @copydoc C2::ContractSym
	 */
	static NamedForm ContractSym(const NamedVValuedForm<C2_U2>& a, const NamedVValuedForm<C2_U2>& b)
	{
		ex result;
		for (int i=0;i<a.size();i++)
			result+=a[i]*b[i];
		return NamedForm(a.name()+b.name(),a.texname()+b.texname(),result.expand());
	}

	/** @copydoc C2::Contract
	 */
	static NamedForm Contract(const NamedVValuedForm<C2_U2>& a, const NamedVValuedForm<C2_U2>& b)
	{
		return 	ContractSym(a,b);
	}

	/** @copydoc C2::ContractSkew1
	 */
	static NamedForm ContractSkew(const NamedVValuedForm<C2_U2>& a, const NamedVValuedForm<C2_U2>& b)
	{
		ex result=a(1)*b(4)-a(4)*b(1)-a(2)*b(3)+a(3)*b(2);
		return NamedForm("("+a.name()+","+b.name()+")","\\sigma("+a.texname()+","+b.texname()+")",result.expand());
	}

	/** @copydoc R3::CreateRepresentationInfo
	 */
	static CohomogeneityOneInvariantForms::RepresentationInfo CreateRepresentationInfo(const NamedVValuedForm<C2_U2>& coordinates, ex r, int radius=1)
	{
		CohomogeneityOneInvariantForms::RepresentationInfo result;
		result.specialPoint= lst(coordinates[0]==0,coordinates[1]==0,coordinates[2]==0,coordinates[3]==0);
		result.principalPoint=lst(coordinates[0]==radius,coordinates[1]==0,coordinates[2]==0,coordinates[3]==0);
		result.genericPoint=lst(coordinates[0]==r,coordinates[1]==0,coordinates[2]==0,coordinates[3]==0);
		return result;
	}

	/** @copydoc R3::AllContractions
	 */
	template<typename InputIterator> static list<NamedForm> AllContractions(InputIterator from, InputIterator to)
	{
		list<NamedForm> simpleElements;
		for (InputIterator i=from;i!=to;i++)
			for (InputIterator j=i;j!=to;j++)
			{
				addElement(simpleElements,ContractSym(*i,*j));
				addElement(simpleElements,ContractSkew(*i,*j));
			}
		return simpleElements;
	}

	/** @copydoc R3::LieAlgebraAction
	 *
	 * @remark In quaternionic terms, the action of \f$\mathfrak{u(2)}\f$ on \f$\mathbb{C}^2\f$ is implemented as: \f$(L_i,L_j,L_k,R_k)\f$
	 */
	static NamedVValuedForm<C2_U2> LieAlgebraAction(const ExVector& omega, const NamedVValuedForm<C2_U2>& alpha)
	{
		//notice the use of [] to allow for zero-based indices on quaternions, as customary.
		NamedVValuedForm<C2_U2> result;
		//left multiplication
		result[0]=-omega(1)*alpha[1]-omega(2)*alpha[2]-(omega(3)+omega(4))*alpha[3];
		result[1]=omega(1)*alpha[0]+omega(2)*alpha[3]-(omega(3)-omega(4))*alpha[2];
		result[2]=omega(2)*alpha[0]+(omega(3)-omega(4))*alpha[1]-omega(1)*alpha[3];
		result[3]=(omega(3)+omega(4))*alpha[0]+omega(1)*alpha[2]-omega(2)*alpha[1];
		return result;
	}

	/** @copydoc R3::Contract
	 */
	ex Contract(string contractionname,const list<NamedVValuedForm<C2_U2> >& letters) const
	{
		if (letters.size()!=2)
			throw InvalidArgument(__FILE__,__LINE__);
		else if (contractionname=="") return ContractSym(letters.front(),letters.back());
		else if (contractionname=="sigma") return ContractSkew(letters.front(),letters.back());
		else throw InvalidArgument(__FILE__,__LINE__);
	}
private:
	static void addElement(list<NamedForm>& list, NamedForm elem)
	{
		if (!elem.expand().is_zero()) list.push_back(elem);
	}
};




#endif /* DICTIONARY_REPRESENTATIONS_C2_U2_H_ */
