/*
 * SU3.h
 *
 *  Created on: Mar 27, 2015
 *      Author: diego
 */

#ifndef DICTIONARY_GROUPS_SU3_H_
#define DICTIONARY_GROUPS_SU3_H_

#include <wedge/liegroup.h>

using namespace Wedge;

namespace Groups {

class AbstractSU3 : public AbstractLieGroup<> {
public:
	AbstractSU3() :AbstractLieGroup<> ("-23-45+2*67,13+46-57-[sqrt(3)]*58,-12-47+[sqrt(3)]*48-56,"
			"15-26+37-[sqrt(3)]*38,"
			"-14+27+36+[sqrt(3)]*28,-2*17+24-35,2*16-25-34,-[sqrt(3)]*25+[sqrt(3)]*34") {}

	void Declare_d(ex x, ex dx)
	{
		AbstractLieGroup<>::Declare_d(x,dx);
	}

};

extern AbstractSU3 SU3;

}
#endif /* DICTIONARY_GROUPS_SU3_H_ */
