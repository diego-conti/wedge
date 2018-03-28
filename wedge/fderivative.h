/*
 * fderivative.h
 *
 *  Created on: Apr 27, 2015
 *      Author: diego
 */

#ifndef WEDGE_FDERIVATIVE_H_
#define WEDGE_FDERIVATIVE_H_

#include <ginac/fderivative.h>
#include <vector>
#include <cassert>

namespace Wedge {
using namespace GiNaC;
using namespace std;

/** Return the multiindex associated to a fderivative object.
 *
 * @param f an fderivative object
 * @result a vector with a number of elements equal to the number of arguments of f whose values
 * are the order of derivation
 *
 * e.g. D[0,0]f(x,y,z) is mapped to the vector [2,0,0]
 */
vector<int> FunctionDerivativeToMultiIndex(const fderivative& f);

/** Create an fderivative object from a multiindex, a serial number, and a list of arguments */
ex FunctionDerivativeFromMultiIndex(unsigned serial, const vector<int>& orders, const exvector& args);


}
#endif /* WEDGE_FDERIVATIVE_H_ */
