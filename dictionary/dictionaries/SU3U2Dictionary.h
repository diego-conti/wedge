/*
 * SU3U2Dictionary.h
 *
 *  Created on: Mar 27, 2015
 *      Author: diego
 */

#ifndef DICTIONARY_SU3U2DICTIONARY_H_
#define DICTIONARY_SU3U2DICTIONARY_H_

#include "dictionary/representations/C2_U2.h"
#include "dictionary/groups/SU3.h"
#include "dictionary/dictionary.h"

using namespace Wedge;
/** @brief The dictionary of \f$SU(3)\f$-invariant forms on \f$SU(3)\times_{U(2)}\mathbb{C}^2\f$
*/
class SU3U2Dictionary : public Dictionary<C2_U2,4> {
public:
	typedef NamedVValuedForm NamedC2ValuedForm;

/** @brief Construct a SU3U2Dictionary object, without starting the actual calculations
 * @param restrictToSphere If true, all computations are carried out restricting to the sphere bundle, i.e. a principal orbit
 * @param upToDegree If specified, forms of degree higher than upToDegree are ignored
 */
	SU3U2Dictionary(bool restrictToSphere=false, int upToDegree=1000) : Dictionary<C2_U2,4>(&Groups::SU3,restrictToSphere,upToDegree) {
		ConnectionForm(1)=Groups::SU3.e(1);
		ConnectionForm(2)=-Groups::SU3.e(6);
		ConnectionForm(3)=Groups::SU3.e(7);
		ConnectionForm(4)=Groups::SU3.e(8)*sqrt(ex(3));
		exvector horizontal=ParseDifferentialForms(Groups::SU3.e(),"2,3,4,5");
		addHorizontalForms(horizontal.begin(),horizontal.end());

		LOG_INFO(b());
		LOG_INFO(beta());
		LOG_INFO(starbeta());
		LOG_INFO(epsilon());
		LOG_INFO(gamma());
	}
/** @copydoc Dictionary::CreateAlgebra */
	void CreateAlgebra() {
		list<NamedC2ValuedForm> letters;
		letters.push_back(a());
		letters.push_back(b());
		letters.push_back(beta());
		letters.push_back(starbeta());
		letters.push_back(epsilon());
		letters.push_back(gamma());
		LOG_INFO(b());
		LOG_INFO(beta());
		LOG_INFO(starbeta());
		LOG_INFO(epsilon());
		LOG_INFO(gamma());

		list<NamedForm> words=AllContractions(letters.begin(),letters.end());
		Dictionary<C2_U2,4>::CreateAlgebra(words.begin(),words.end());

	}

	NamedC2ValuedForm beta() const {
		return NamedC2ValuedForm("beta","\\beta",ParseDifferentialForms(Groups::SU3.e(),"2,3,4,5"));
	}

	NamedC2ValuedForm starbeta() const {
		return NamedC2ValuedForm("*beta","\\tilde\\beta",ParseDifferentialForms(Groups::SU3.e(),"345,-245,235,-234"));
	}

	NamedC2ValuedForm epsilon() const {
		ExVector e=Groups::SU3.e();
		NamedC2ValuedForm epsilon("eps","\\epsilon");
		for (int i=1;i<5;i++)
			for (int j=1;j<5;j++)
				epsilon(i)+=a(j)*e(j+1)*e(i+1);
		return epsilon;
	}

	NamedC2ValuedForm gamma() const {
		return D(epsilon()).set_name("gamma","\\gamma");
	}

	list<NamedC2ValuedForm> Letters() const
	{
		list<NamedC2ValuedForm > result;
		result.push_back(beta());
		result.push_back(starbeta());
		result.push_back(epsilon());
		result.push_back(gamma());
		result.push_back(a());
		result.push_back(b());
		return result;
	}
};





#endif /* DICTIONARY_SU3U2DICTIONARY_H_ */
