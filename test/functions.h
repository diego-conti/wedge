/***************************************************************************
 *   Copyright (C) 2015	by Diego Conti			  							*
 *   diego.conti@unimib.it                                                 *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 3 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#ifndef TEST_FUNCTIONS_H_
#define TEST_FUNCTIONS_H_
#include <list>
#include <vector>
#include <cxxtest/TestSuite.h>
#include <ginac/ginac.h>
#include "test.h"
#include "wedge/vectorspace.h"
#include "wedge/utilities.h"
#include "wedge/differentialform.h"
#include "wedge/tensor.h"
#include "wedge/lambda.h"
#include "wedge/pforms.h"
#include "wedge/wedgealgebraic.h"
#include "wedge/affinebasis.h"
#include "wedge/polybasis.h"
#include "wedge/fderivative.h"

using namespace GiNaC;
using namespace Wedge;

DECLARE_FUNCTION_1P(f)
REGISTER_FUNCTION(f, dummy())
DECLARE_FUNCTION_2P(g)
REGISTER_FUNCTION(g, dummy())
class FunctionsTestSuite : public CxxTest::TestSuite
{
public:
	void testFunctionDerivativeFromMultiIndex() {

		realsymbol r("r"), x("x");
		vector<int> indices;
		indices.push_back(1);
		indices.push_back(0);
		exvector args;
		args.push_back(r); args.push_back(x);
		unsigned serial = g(0,0).get_serial();
		ex y=FunctionDerivativeFromMultiIndex(serial, indices,args);
		TS_ASSERT_EQUALS(y,g(r,x).diff(r));
		indices[0]=2;
		y=FunctionDerivativeFromMultiIndex(serial, indices,args);
		TS_ASSERT_EQUALS(y,g(r,x).diff(r,2));
		indices[1]=3;
		y=FunctionDerivativeFromMultiIndex(serial, indices,args);
		TS_ASSERT_EQUALS(y,g(r,x).diff(r,2).diff(x,3));
		indices[0]=0;
		y=FunctionDerivativeFromMultiIndex(serial, indices,args);
		TS_ASSERT_EQUALS(y,g(r,x).diff(x,3));
	}
	void testFunctionDerivativeToMultiIndex() {
		realsymbol r("r"), x("x");
		ex y=f(r);
		ex dy=y.diff(r,2);
		vector<int> indices=FunctionDerivativeToMultiIndex(ex_to<fderivative>(dy));
		TS_ASSERT_EQUALS(indices.size(),1);
		TS_ASSERT_EQUALS(indices[0],2);
		dy=dy.diff(r);
		indices=FunctionDerivativeToMultiIndex(ex_to<fderivative>(dy));
		TS_ASSERT_EQUALS(indices.size(),1);
		TS_ASSERT_EQUALS(indices[0],3);
		dy=dy.diff(r,2);
		indices=FunctionDerivativeToMultiIndex(ex_to<fderivative>(dy));
		TS_ASSERT_EQUALS(indices.size(),1);
		TS_ASSERT_EQUALS(indices[0],5);

		dy=g(r,x).diff(r,2);
		indices=FunctionDerivativeToMultiIndex(ex_to<fderivative>(dy));
		TS_ASSERT_EQUALS(indices.size(),2);
		TS_ASSERT_EQUALS(indices[0],2);
		TS_ASSERT_EQUALS(indices[1],0);
		dy=dy.diff(x);
		indices=FunctionDerivativeToMultiIndex(ex_to<fderivative>(dy));
		TS_ASSERT_EQUALS(indices.size(),2);
		TS_ASSERT_EQUALS(indices[0],2);
		TS_ASSERT_EQUALS(indices[1],1);


	}
};


#endif /* TEST_FUNCTIONS_H_ */
