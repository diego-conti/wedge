/*
 * R3.h
 *
 *  Created on: Mar 27, 2015
 *      Author: diego
 */

#ifndef DICTIONARY_REPRESENTATIONS_R3_H_
#define DICTIONARY_REPRESENTATIONS_R3_H_

#include "../vvaluedforms.h"

/** @brief \f$\mathbb{R}^3\f$ as a representation of \f$SO(3)\f$
*/
struct R3
{
	static int FibreSize() {return 3;} ///<Dimension of \f$V=\mathbb{R}^3\f$

	/** @brief The contraction corresponding to the determinant
	 * @param a,b,c Elements of \f$\mathcal{L}\f$
	 * @return The contraction of a, b and c.
	 */
	static NamedForm Contract(const NamedVValuedForm<R3>& a, const NamedVValuedForm<R3>& b,const NamedVValuedForm<R3>& c)
	{
		ex form=a(1)*b(2)*c(3)+a(2)*b(3)*c(1)+a(3)*b(1)*c(2)-
			a(1)*b(3)*c(2)-a(2)*b(1)*c(3)-a(3)*b(2)*c(1);
		return NamedForm(
			a.name()+b.name()+c.name(),a.texname()+b.texname()+c.texname(),form.expand());
	}

	/** @brief The contraction corresponding to the scalar product
	 * @param a,b Elements of \f$\mathcal{L}\f$
	 * @return The contraction \f$a\cdot b\f$
	 */
	static NamedForm Contract(const NamedVValuedForm<R3>& a, const NamedVValuedForm<R3>& b)
	{
		ex result;
		for (int i=0;i<a.size();i++)
			result+=a[i]*b[i];
		return NamedForm(a.name()+b.name(),a.texname()+b.texname(),result.expand());
	}

	/** @brief Construct a RepresentationInfo object relative to this representation
	 * @param coordinates The coordinates, obtained by Dictionary::a()
	 * @param radius The radius of the principal point to consider
	* @param r The radius function
	 */
	static CohomogeneityOneInvariantForms::RepresentationInfo CreateRepresentationInfo(const NamedVValuedForm<R3>& coordinates, ex r,int radius=1)
	{
		CohomogeneityOneInvariantForms::RepresentationInfo result;
		result.specialPoint= lst(coordinates[0]==0,coordinates[1]==0,coordinates[2]==0);
		result.principalPoint=lst(coordinates[0]==radius,coordinates[1]==0,coordinates[2]==0);
		result.genericPoint=lst(coordinates[0]==r,coordinates[1]==0,coordinates[2]==0);
		return result;
	}

	/** @brief Compute all the possible contractions obtained from  a given subset of \f$\mathcal{L}\f$
	 * @param [from,to) A range of elements of \f$\mathcal{L}\f$
	 * @return A list of all the contractions that can be obtained from the given range
	 */
	template<typename InputIterator> static list<NamedForm> AllContractions(InputIterator from, InputIterator to)
	{
		list<NamedForm> simpleElements;
		for (InputIterator i=from;i!=to;i++)
			for (InputIterator j=i;j!=to;j++)
			{
				simpleElements.push_back(Contract(*i,*j));
				for (InputIterator k=j;k!=to;k++)
					simpleElements.push_back(Contract(*i,*j,*k));
			}
		return simpleElements;
	}

	/** @brief Compute the contraction of a connection form with a \f$V\f$-valued pseudotensorial form
	 * @param omega The connection form
	 * @param alpha The V-valued pseudotensorial form
	*	@return The contraction \f$\omega\cdot\alpha\f$ induced by the infinitesimal action \f$\mathfrak{h}\otimes V\to V\f$
	 */
	static NamedVValuedForm<R3> LieAlgebraAction(const ExVector& omega, const NamedVValuedForm<R3>& alpha)
	{
		NamedVValuedForm<R3> result;
		result(1)=omega(2)*alpha(3)-omega(3)*alpha(2);
		result(2)=omega(3)*alpha(1)-omega(1)*alpha(3);
		result(3)=omega(1)*alpha(2)-omega(2)*alpha(1);
		return result;
	}

	/** @brief Compute a contraction indicated by name
	 *  @param contractionname The name of the contraction, if needed
	 * @param letters Elements of \f$\mathcal{L}\f$
	 * @return The result of contracting letters with contractionname
	 *
	 * In general, if there is only one contraction that takes a certain number of letters, the name is not used
	 */

	ex Contract(string contractionname,const list<NamedVValuedForm<R3> >& letters) const
	{
		if (letters.size()==3)
			return Contract(letters.front(),*++letters.begin(),letters.back());
		else if (letters.size()==2)
			return Contract(letters.front(),letters.back());
		else
			throw InvalidArgument(__FILE__,__LINE__);
	}

	/** @brief Return a vertical frame with a certain invariance property */
	static ExVector Frame(const NamedVValuedForm<R3>& a, const NamedVValuedForm<R3>& b, ex r)
	{
		ExVector x;
		x.push_back(Contract(a,b)/r);
/*		ex s=pow(a(1),2)+pow(a(2),2);
		x.push_back((-a(2)*b(1)+a(1)*b(2))/sqrt(s));
		x.push_back((a(1)*a(3)*b(1)+a(2)*a(3)*b(2)-s*b(3))/(r*sqrt(s)));
*/
//we are only going to take two derivatives, so we put in second order Taylor expansion at point (r,0,0)
		x.push_back((1-pow(a(2),2)/(2*pow(a(1),2)))*b(2)-a(2)/a(1)*b(1));
		x.push_back(a(3)/a(1)*b(1) + a(2)*a(3)/pow(a(1),2)*b(2)-(1-pow(a(3),2)/(2*pow(a(1),2)))*b(3));
		return x;
	}

};




#endif /* DICTIONARY_REPRESENTATIONS_R3_H_ */
