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

#include <ginac/ginac.h>
#include <ginac/version.h>
#include <iostream>
#include "wedge/wedgealgebraic.h"
#include <time.h>
#include "wedge/aragelilinalg.h"
#include "wedge/ginaclinalg.h"

//this program compares the speed of GinacLinAlg::IndependenceMatrix and ArageliLinAlg::IndependenceMatrix
//using GiNaC 1.3.9, Wedge 0.1.0 and Arageli 2.2.1.191, GinacLinAlg is much faster.
using namespace std;
using namespace Wedge;

struct Clock {
	clock_t start_clock;
	const char* testname;
	Clock(const char* _testname) : start_clock(clock()) {testname=_testname;}
	~Clock() {
		cout<<testname;
		cout<<": time elapsed: "<<(1000* (clock()-start_clock))/CLOCKS_PER_SEC <<"ms"<<endl;
	}
};

int main() {	
	for (int size=10;size<30;size+=5)
	{
		GinacLinAlgAlgorithms::IndependenceMatrix m(size,size);
		ArageliLinAlgAlgorithms::IndependenceMatrix n(size,size);
		for (int i=0;i<size;i++)
			for (int j=0;j<size;j++)
			{
				m.M(i,j)=rand();				
				n.M(i,j)=m.M(i,j);
			}
		cout<<"Matrix size = "<<size<<endl;
		{
			Clock c("GinacLinAlgAlgorithms::IndependenceMatrix");			
			m.ChooseLinearlyIndependentRows();
		}
		{
			Clock c("ArageliLinAlgAlgorithms::IndependenceMatrix");
			n.ChooseLinearlyIndependentRows();
		}
		ArageliLinAlgAlgorithms::IndependenceMatrix::const_iterator j=n.IndependentRowsBegin();
		for (GinacLinAlgAlgorithms::IndependenceMatrix::const_iterator i=m.IndependentRowsBegin();
			i!=m.IndependentRowsEnd() && j!=n.IndependentRowsEnd();i++,j++)
			assert(*i==*j);
	}
	return 0;
}

