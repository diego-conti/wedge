/*******************************************************************************
 *  Copyright (C) 2007-2023 by Diego Conti, diego.conti@unipi.it 
 *  This file is part of Wedge.                                           
 *  Wedge is free software; you can redistribute it and/or modify         
 *  it under the terms of the GNU General Public License as published by  
 *  the Free Software Foundation; either version 3 of the License, or     
 *  (at your option) any later version.                                   
 *                                                                          
 *  Wedge is distributed in the hope that it will be useful,              
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of        
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         
 *  GNU General Public License for more details.                          
 *                                                                           
 *  You should have received a copy of the GNU General Public License     
 *  along with Wedge; if not, write to the                                
 *   Free Software Foundation, Inc.,                                       
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             
 *  
 *******************************************************************************/

#include <list>
#include <vector>
#include <cxxtest/TestSuite.h>
#include <ginac/ginac.h>
#include "wedge/base/expressions.h"

//To debug, run cxxtestgen with the flag --no-eh
#ifndef _CXXTEST_HAVE_EH
#undef TS_ASSERT_THROWS
#define TS_ASSERT_THROWS(p1,p2) 
#endif 
//#define LONG_TEST

using namespace GiNaC;
using namespace Wedge;

#define BASIC_DERIVED_EQUALS(a,b) ((a-b).is_zero())
#define BASIC_DERIVED_DIFFERS(a,b) (!(a-b).is_zero())
 
#ifdef CXXTEST_RUNNING
#include <cxxtest/ValueTraits.h>
#include <stdio.h>
#include <sstream>

namespace CxxTest
{
	CXXTEST_TEMPLATE_INSTANTIATION
	class ValueTraits<ex>
	{
		string s;
		public:
			ValueTraits( const ex &m ) {
				stringstream stream;
				stream<< m;
				s=stream.str();
			}
			const char *asString() const { return s.c_str(); }
	};
	CXXTEST_TEMPLATE_INSTANTIATION
	class ValueTraits<matrix>
	{
		string s;
		public:
			ValueTraits( const matrix &m ) {
				stringstream stream;
				stream<< m;
				s=stream.str();
			}
			const char *asString() const { return s.c_str(); }
	};
	CXXTEST_TEMPLATE_INSTANTIATION
	class ValueTraits<exvector>
	{
		string s;
		public:
			ValueTraits( const exvector &m ) {
				stringstream stream;
				stream<< m;
				s=stream.str();
			}
			const char *asString() const { return s.c_str(); }
	};
	CXXTEST_TEMPLATE_INSTANTIATION
	class ValueTraits<ExVector>
	{
		string s;
		public:
			ValueTraits( const ExVector &m ) {
				stringstream stream;
				stream<< m;
				s=stream.str();
			}
			const char *asString() const { return s.c_str(); }
	};
	CXXTEST_TEMPLATE_INSTANTIATION
	class ValueTraits<basic>
	{
		string s;
		public:
			ValueTraits( const basic &m ) {
				stringstream stream;
				stream<< ex(m);
				s=stream.str();
			}
			const char *asString() const { return s.c_str(); }
	};
}


class DummyTestSuite : public  CxxTest::TestSuite
{
public:
	void testNothing() {}
};
#endif // CXXTEST_RUNNING



