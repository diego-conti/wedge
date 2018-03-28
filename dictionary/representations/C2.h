/*
 * C2.h
 *
 *  Created on: Mar 27, 2015
 *      Author: diego
 */

#ifndef DICTIONARY_REPRESENTATIONS_C2_H_
#define DICTIONARY_REPRESENTATIONS_C2_H_


#include "../vvaluedforms.h"


/** @brief \f$\mathbb{C}^2\f$ as a representation of \f$SU(2)\f$
*/
struct C2
{
	static int FibreSize() {return 4;}	///<Real Dimension of \f$V=\mathbb{C}^2\f$

	/** @brief The contraction corresponding to the scalar product
	 * @param a,b Elements of \f$\mathcal{L}\f$
	 * @return The contraction \f$a\cdot b\f$
	 */
	static NamedForm ContractSym(const NamedVValuedForm<C2>& a, const NamedVValuedForm<C2>& b)
	{
		ex result;
		for (int i=0;i<a.size();i++)
			result+=a[i]*b[i];
		return NamedForm(a.name()+b.name(),a.texname()+b.texname(),result.expand());
	}

	/** @brief The contraction corresponding to the two-form \f$e^{12}-e^{34}\f$
	 * @param a,b Elements of \f$\mathcal{L}\f$
	 * @return The contraction of \f$a\f$ and \f$b\f$
	 */
	static NamedForm ContractSkew1(const NamedVValuedForm<C2>& a, const NamedVValuedForm<C2>& b)
	{
		ex result=a(1)*b(2)-a(2)*b(1)-a(3)*b(4)+a(4)*b(3);
		return NamedForm("_1("+a.name()+","+b.name()+")","\\sigma_1("+a.texname()+","+b.texname()+")",result.expand());
	}

	/** @brief The contraction corresponding to the two-form \f$e^{13}-e^{42}\f$
	 * @param a,b Elements of \f$\mathcal{L}\f$
	 * @return The contraction of \f$a\f$ and \f$b\f$
	 */
	static NamedForm ContractSkew2(const NamedVValuedForm<C2>& a, const NamedVValuedForm<C2>& b)
	{
		ex result=a(1)*b(3)-a(3)*b(1)-a(4)*b(2)+a(2)*b(4);
		return NamedForm("_2("+a.name()+","+b.name()+")","\\sigma_2("+a.texname()+","+b.texname()+")",result.expand());
	}

	/** @brief The contraction corresponding to the two-form \f$e^{14}-e^{23}\f$
	 * @param a,b Elements of \f$\mathcal{L}\f$
	 * @return The contraction of \f$a\f$ and \f$b\f$
	 */
	static NamedForm ContractSkew3(const NamedVValuedForm<C2>& a, const NamedVValuedForm<C2>& b)
	{
		ex result=a(1)*b(4)-a(4)*b(1)-a(2)*b(3)+a(3)*b(2);
		return NamedForm("_3("+a.name()+","+b.name()+")","\\sigma_3("+a.texname()+","+b.texname()+")",result.expand());
	}


	/** @copydoc R3::CreateRepresentationInfo
	 */
	static CohomogeneityOneInvariantForms::RepresentationInfo CreateRepresentationInfo(const NamedVValuedForm<C2>& coordinates, ex r,int radius=1)
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
				addElement(simpleElements,ContractSkew1(*i,*j));
				addElement(simpleElements,ContractSkew2(*i,*j));
				addElement(simpleElements,ContractSkew3(*i,*j));
			}
		return simpleElements;
	}


	/** @copydoc R3::LieAlgebraAction
	 */
	static NamedVValuedForm<C2> LieAlgebraAction(const ExVector& omega, const NamedVValuedForm<C2>& alpha)
	{
		//notice the use of [] to allow for zero-based indices on quaternions, as customary.
		NamedVValuedForm<C2> result;
		//left multiplication
		result[0]=-omega(1)*alpha[1]-omega(2)*alpha[2]-omega(3)*alpha[3];
		result[1]=omega(1)*alpha[0]+omega(2)*alpha[3]-omega(3)*alpha[2];
		result[2]=omega(2)*alpha[0]+omega(3)*alpha[1]-omega(1)*alpha[3];
		result[3]=omega(3)*alpha[0]+omega(1)*alpha[2]-omega(2)*alpha[1];
		for (int i=0;i<4;i++) result[i]/=2;
		return result;
	}

	/** @copydoc R3::Contract
	 */
	ex Contract(string contractionname,const list<NamedVValuedForm<C2> >& letters) const
	{
		if (letters.size()!=2)
			throw InvalidArgument(__FILE__,__LINE__);
		else if (contractionname=="") return ContractSym(letters.front(),letters.back());
		else if (contractionname=="sigma1") return ContractSkew1(letters.front(),letters.back());
		else if (contractionname=="sigma2") return ContractSkew1(letters.front(),letters.back());
		else if (contractionname=="sigma3") return ContractSkew1(letters.front(),letters.back());
		else throw InvalidArgument(__FILE__,__LINE__);
	}
private:
	static void addElement(list<NamedForm>& list, NamedForm elem)
	{
		if (!elem.expand().is_zero()) list.push_back(elem);
	}
};


#endif /* DICTIONARY_REPRESENTATIONS_C2_H_ */
