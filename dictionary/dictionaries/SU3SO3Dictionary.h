/*
 * SU3SO3Dictionary.h
 *
 *  Created on: Mar 27, 2015
 *      Author: diego
 */

#ifndef DICTIONARY_SU3SO3DICTIONARY_H_
#define DICTIONARY_SU3SO3DICTIONARY_H_

#include "dictionary/representations/R3.h"
#include "dictionary/groups/SU3.h"
#include "dictionary/dictionary.h"

/** @brief The dictionary of \f$SU(3)\f$-invariant forms on \f$SU(3)\times_{SO(3)}\mathbb{R}^3\f$
*/
class SU3SO3Dictionary : public Dictionary<R3,3> {
public:
	typedef NamedVValuedForm NamedR3ValuedForm;
	bool witheps;

/** @brief Construct a SU3SO3Dictionary object, without starting the actual calculations
 * @param restrictToSphere If true, all computations are carried out restricting to the sphere bundle, i.e. a principal orbit
 * @param upToDegree If specified, forms of degree higher than upToDegree are ignored
 * @param witheps If true, use the letter \f$\epsilon\f$
 */
	SU3SO3Dictionary(bool restrictToSphere=false, int upToDegree=1000, bool _witheps=true) : Dictionary<R3,3>(&Groups::SU3,restrictToSphere,upToDegree) {
		witheps=_witheps;
		ConnectionForm(1)=Groups::SU3.e(1);
		ConnectionForm(2)=Groups::SU3.e(2);
		ConnectionForm(3)=Groups::SU3.e(3);
		exvector horizontal=ParseDifferentialForms(Groups::SU3.e(),"4,5,6,7,8");
		addHorizontalForms(horizontal.begin(),horizontal.end());

		LOG_INFO(a());
		LOG_INFO(b());
		LOG_INFO(Beta());
		LOG_INFO(StarBeta());
		LOG_INFO(epsilon());
		LOG_INFO(gamma());
	}
/** @copydoc Dictionary::CreateAlgebra */
	void CreateAlgebra()
	{
		list<NamedR3ValuedForm> letters;
		letters.push_back(a());
		letters.push_back(b());
		letters.push_back(Beta());
		letters.push_back(StarBeta());
		letters.push_back(gamma());
		if (witheps) letters.push_back(epsilon());


		list<NamedForm> words=AllContractions(letters.begin(),letters.end());
		words.push_back(Nu());
		Dictionary<R3,3>::CreateAlgebra(words.begin(),words.end());
	}
/** @brief The volume form on the base
*/
	NamedForm Nu() const
	{
		return NamedForm(N.nu,ParseDifferentialForm(Groups::SU3.e(),"45678"));
	}

	NamedR3ValuedForm Beta() const{
		return NamedR3ValuedForm(N.beta,ParseDifferentialForms(Groups::SU3.e(),"-2*76-45,[sqrt(3)]*85+75+46,74-[sqrt(3)]*84+65"));
	}

	NamedR3ValuedForm gamma() const {
		return D(epsilon()).set_name(N.gamma);
	}
	NamedR3ValuedForm epsilon() const {
		NamedR3ValuedForm epsilon(N.epsilon);
		ex sqrt3=sqrt(ex(3));
		ExVector e=Groups::SU3.e();
		epsilon(1)=-4*a(1)*e(8)-2*sqrt3*a(2)*e(4)-2*sqrt3*a(3)*e(5);
		epsilon(2)=-2*sqrt3*a(1)*e(4)+a(2)*(2*e(8)-2*sqrt3*e(7))-2*sqrt3*a(3)*e(6);
		epsilon(3)=-2*sqrt3*a(1)*e(5)-2*sqrt3*a(2)*e(6) +a(3)*(2*sqrt3*e(7)+2*e(8));
		return epsilon;
	}
	NamedR3ValuedForm StarBeta() const {
		return NamedR3ValuedForm("beta","\\tilde\\beta",ParseDifferentialForms(Groups::SU3.e(),"2*458-678,[-sqrt(3)]*467+468-578,-568-[sqrt(3)]*567-478"));
	}

	list<NamedForm> FormsOnBase() const
	{
		list<NamedForm> result;
		result.push_back(Nu());
		return result;
	}

	list<NamedR3ValuedForm> Letters() const
	{
		list<NamedR3ValuedForm > result;
		result.push_back(Beta());
		result.push_back(StarBeta());
		result.push_back(epsilon());
		result.push_back(gamma());
		result.push_back(a());
		result.push_back(b());
		return result;
	}
};



#endif /* DICTIONARY_SU3SO3DICTIONARY_H_ */
