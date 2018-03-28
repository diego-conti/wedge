/***************************************************************************
 *   Copyright (C) 2007, 2008, 2009 by Diego Conti			   *
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

#ifndef WEDGETEST_H
#define WEDGETEST_H
#include <list>
#include <vector>
#include <cxxtest/TestSuite.h>
#include <ginac/ginac.h>
#include "wedge/expressions.h"
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

//To debug, run cxxtestgen with the flag --no-eh and compile with the symbol DEBUG defined:
#ifdef DEBUG
#warning "DEBUG MODE: not checking that exceptions are thrown as expected"
#undef TS_ASSERT_THROWS
#define TS_ASSERT_THROWS(p1,p2) 
#endif 
//#define LONG_TEST

using namespace GiNaC;
using namespace Wedge;

WEDGE_DECLARE_NAMED_ALGEBRAIC(V,Vector)

WEDGE_DECLARE_NAMED_ALGEBRAIC(W,Vector)


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

class Fixture : public CxxTest::GlobalFixture {
public:
	 bool setUpWorld() {LOG_INFO("Test started"<<MaxLength(10000)); return true; }
   bool tearDownWorld() {LOG_INFO("Test finished"); return true; }
};

static Fixture fixture;

class DummyTestSuite : public  CxxTest::TestSuite
{
public:
	void testNothing() {}
};
#endif // CXXTEST_RUNNING


#endif
