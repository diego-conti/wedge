/*
 * SU3SU2Dictionary.h
 *
 *  Created on: Mar 27, 2015
 *      Author: diego
 */

#ifndef DICTIONARY_SU3SU2DICTIONARY_H_
#define DICTIONARY_SU3SU2DICTIONARY_H_

#include "dictionary/representations/R3.h"
#include "dictionary/groups/SU3.h"
#include "dictionary/dictionary.h"

/** @brief The dictionary of \f$SU(3)\f$-invariant forms on \f$SU(3)\times_{SU(2)}\mathbb{R}^3\f$
*/
class SU3SU2Dictionary : public Dictionary<R3,3> {
public:
	typedef NamedVValuedForm NamedR3ValuedForm;

/** @brief Construct a SU3SU2Dictionary object, without starting the actual calculations
 * @param restrictToSphere If true, all computations are carried out restricting to the sphere bundle, i.e. a principal orbit
 * @param upToDegree If specified, forms of degree higher than upToDegree are ignored
 */
	SU3SU2Dictionary(bool restrictToSphere=false, int upToDegree=1000) : Dictionary<R3,3>(&Groups::SU3,restrictToSphere,upToDegree) {
		ConnectionForm(1)=Groups::SU3.e(1)*2;
		ConnectionForm(2)=Groups::SU3.e(7)*2;
		ConnectionForm(3)=Groups::SU3.e(6)*2;
		exvector horizontal=ParseDifferentialForms(Groups::SU3.e(),"2,3,4,5,8");
		addHorizontalForms(horizontal.begin(),horizontal.end());
	}
/** @copydoc Dictionary::CreateAlgebra */
	void CreateAlgebra() {
		list<NamedR3ValuedForm> letters;
		letters.push_back(a());
		letters.push_back(b());
		letters.push_back(beta());
		list<NamedForm> words=AllContractions(letters.begin(),letters.end());

		NamedR3ValuedForm Dsigma=D(beta());
		assert(Dsigma(1).expand().is_zero());
		assert(Dsigma(2).expand().is_zero());
		assert(Dsigma(3).expand().is_zero());

		words.push_back(alpha());
		words.push_back(omega1());
		words.push_back(omega2());
		words.push_back(omega3());
		Dictionary<R3,3>::CreateAlgebra(words.begin(),words.end());
	}

//multiply everything by 24 for consistence with [Bryant-Salamon, On the construction of some complete metrics with exceptional holonomy.  Duke Math. J.  58  (1989),  no. 3, 829--850.]
	NamedForm alpha() const {ExVector e=Groups::SU3.e(); return NamedForm("eta","\\eta",e(8)*2/sqrt(ex(3)));}
	NamedForm omega1() const {ExVector e=Groups::SU3.e(); return NamedForm("omega_1","\\omega_1",(e(2)*e(5)+e(4)*e(3)));}
	NamedForm omega2() const {ExVector e=Groups::SU3.e(); return NamedForm("omega_2","\\omega_2",(e(2)*e(4)+e(3)*e(5)));}
	NamedForm omega3() const {ExVector e=Groups::SU3.e(); return NamedForm("omega_3","\\omega_3",(e(2)*e(3)+e(5)*e(4)));}

	NamedR3ValuedForm beta() const {
		ExVector beta=ParseDifferentialForms(Groups::SU3.e(),"23-54,25-43,-24+35");
		return NamedR3ValuedForm("beta","\\beta",beta);
	}

	list<NamedForm> FormsOnBase() const
	{
		list<NamedForm> result;
		result.push_back(alpha());
		result.push_back(omega1());
		result.push_back(omega2());
		result.push_back(omega3());
		return result;
	}

	list<NamedR3ValuedForm> Letters() const
	{
		list<NamedR3ValuedForm > result;
		result.push_back(beta());
		result.push_back(a());
		result.push_back(b());
		return result;
	}
};



#endif /* DICTIONARY_SU3SU2DICTIONARY_H_ */
