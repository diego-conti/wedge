/*
 * SU3U2DictionaryWithWeight.h
 *
 *  Created on: Mar 27, 2015
 *      Author: diego
 */

#ifndef DICTIONARY_DICTIONARIES_SU3U2DICTIONARYWITHWEIGHT_H_
#define DICTIONARY_DICTIONARIES_SU3U2DICTIONARYWITHWEIGHT_H_

#include "dictionary/representations/C2_U2.h"
#include "dictionary/groups/SU3.h"
#include "dictionary/dictionary.h"

using namespace Wedge;




/** @brief The dictionary of \f$SU(3)\f$-invariant forms on \f$SU(3)\times_{U(2)}\mathbb{C}^2\f$, where the action of U(2) on C^2 is (g,v)to det(g) (gv).
*/
class SU3U2DictionaryWithWeight : public Dictionary<C2_U2,4> {
public:
	ex dr;	//d of the radial coordinate.
	typedef NamedVValuedForm NamedC2ValuedForm;

/** @brief Construct a SU3U2Dictionary object, without starting the actual calculations
 * @param restrictToSphere If true, all computations are carried out restricting to the sphere bundle, i.e. a principal orbit
 * @param upToDegree If specified, forms of degree higher than upToDegree are ignored
 */
	SU3U2DictionaryWithWeight(bool restrictToSphere=false, int upToDegree=1000) : Dictionary<C2_U2,4>(&Groups::SU3,restrictToSphere,upToDegree) {
		ConnectionForm(1)=Groups::SU3.e(1);
		ConnectionForm(2)=-Groups::SU3.e(6);
		ConnectionForm(3)=Groups::SU3.e(7);
		ConnectionForm(4)=Groups::SU3.e(8)/sqrt(ex(3));

		ex t=RadialCoordinate();
		//change the default principal point to obtain correct stabilizer, e_1+ sqrt3*e_8
		representationInfo.specialPoint= lst(a(1)==0,a(2)==0,a(3)==0,a(4)==0);
		representationInfo.principalPoint= lst(a(1)==1,a(2)==0,a(3)==-1,a(4)==0);
		representationInfo.genericPoint= lst(a(1)==t,a(2)==0,a(3)==-t,a(4)==0);

		Groups::SU3.Declare_d(RadialCoordinate(),dRadialCoordinate());
		exvector horizontal=ParseDifferentialForms(Groups::SU3.e(),"2,3,4,5");
		addHorizontalForms(horizontal.begin(),horizontal.end());
		setRadiusSquared(2*t*t);
	}
/** @copydoc Dictionary::CreateAlgebra */
	void CreateAlgebra() {
		list<NamedC2ValuedForm> letters;
		letters.push_back(a());
		letters.push_back(b());
		letters.push_back(c());
		letters.push_back(epsilon());
		letters.push_back(gamma());
		letters.push_back(p());
		letters.push_back(q());
		letters.push_back(h());

		LOG_INFO(a());
		LOG_INFO(b());
		LOG_INFO(p());
		LOG_INFO(q());
		LOG_INFO(c());
		LOG_INFO(epsilon());
		LOG_INFO(gamma());

		list<NamedForm> words=AllContractions(letters.begin(),letters.end());
		words.push_back(kappa());
		Dictionary<C2_U2,4>::CreateAlgebra(words.begin(),words.end());

	}

	NamedForm kappa() const {
		ExVector T=ParseDifferentialForms(Groups::SU3.e(),"2,3,4,5");
		return NamedForm("k","\\kappa",T(1)*T(4)-T(2)*T(3));
	}

	NamedC2ValuedForm c() const {
		NamedC2ValuedForm c("c","c");
		ExVector T=ParseDifferentialForms(Groups::SU3.e(),"2,3,4,5");
		ex twoform=T(1)*T(4)-T(2)*T(3);
		for (int i=1;i<=4;++i)
			for (int j=1;j<=4;++j)
				c(i)+=a(j)*(T(i)*T(j)+::Hook(T(i),twoform)*::Hook(T(j),twoform));
		return c;
	}


	NamedC2ValuedForm epsilon() const {
		NamedC2ValuedForm epsilon("eps","\\epsilon");
		NamedC2ValuedForm c=this->c();
		ExVector T=ParseDifferentialForms(Groups::SU3.e(),"2,3,4,5");
		lst subs;
		subs=T(1)==b(1), T(2)==b(2), T(3)==b(3),T(4)==b(4);
		for (int i=1;i<=4;++i)
			epsilon(i)=c(i).subs(subs).expand();
		return epsilon;
	}
	NamedC2ValuedForm p() const {
		NamedC2ValuedForm p("p","p");
		ExVector T=ParseDifferentialForms(Groups::SU3.e(),"2,3,4,5");
		ex M=b(1)*T(1)+b(2)*T(2)+b(3)*T(3)+b(4)*T(4);
		ex N=-b(1)*T(4)+b(4)*T(1)+b(2)*T(3)-b(3)*T(2);
		p(1)=-a(2)*M -a(3)*N;
		p(2)=a(1)*M -a(4)*N;
		p(3)=a(4)*M+a(1)*N;
		p(4)=-a(3)*M+a(2)*N;
		return p;
	}
	NamedC2ValuedForm q() const {
		NamedC2ValuedForm q("q","q");
		ExVector T=ParseDifferentialForms(Groups::SU3.e(),"2,3,4,5");
		ex M=b(1)*T(1)+b(2)*T(2)+b(3)*T(3)+b(4)*T(4);
		ex N=-b(1)*T(4)+b(4)*T(1)+b(2)*T(3)-b(3)*T(2);
		q(1)=-a(2)*N +a(3)*M;
		q(2)=a(1)*N +a(4)*M;
		q(3)=a(4)*N-a(1)*M;
		q(4)=-a(3)*N-a(2)*M;
		return q;
	}

	NamedC2ValuedForm gamma() const {
		return D(epsilon()).set_name("gamma","\\gamma");
	}
	NamedC2ValuedForm h() const {
		return D(c()).set_name("h","h");
	}

	list<NamedC2ValuedForm> Letters() const
	{
		list<NamedC2ValuedForm > result;
		result.push_back(epsilon());
		result.push_back(gamma());
		result.push_back(a());
		result.push_back(b());
		result.push_back(c());
		result.push_back(p());
		result.push_back(q());
		return result;
	}
	list<NamedForm> FormsOnBase() const
	{
		list<NamedForm> result;
		result.push_back(kappa());
		return result;
	}

};





#endif /* DICTIONARY_DICTIONARIES_SU3U2DICTIONARYWITHWEIGHT_H_ */
