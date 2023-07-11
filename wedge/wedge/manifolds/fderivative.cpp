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
#include "fderivative.h"
#include "wedge/base/wedgebase.h"

using namespace GiNaC;
using namespace std;

namespace Wedge {

// Current implementation uses fderivative::archive; in principle one could also parse the output of fderivative::print

vector<int> FunctionDerivativeToMultiIndex(const fderivative& f)
{
	vector<int> orders(f.nops());
	archive ar;
	archive_node node{ar};
	f.archive(node);
	unsigned i = 0;
	while (true) {
		unsigned u;
		if (node.find_unsigned("param", u, i)) {
			assert(u<f.nops());
			++orders[u];
		}
		else break;
		++i;
	}
/*
	if (f.nops()!=1)
		throw NotImplemented(__FILE__,__LINE__,"partial derivatives not implemented");
	ex arg=f.op(0);
	if (!is_a<symbol>(arg)) throw NotImplemented(__FILE__,__LINE__,"argument of derivative function should be a symbol");
	const symbol& parameter=ex_to<symbol>(arg);
	ex primitive=function(f.get_serial(),arg);
	for (int i=1;i<10;++i)
		if (f==primitive.diff(parameter,i)) {
			orders.push_back(i); return orders;
		}
	throw InvalidArgument(__FILE__,__LINE__,f);
	*/
	return orders;
}


ex FunctionDerivativeFromMultiIndex(unsigned serial, const vector<int>& orders, const exvector& args)
{
	paramset parameters;
	assert(orders.size()==args.size());
	for (int i=0;i<orders.size();++i)
		for (int k=0;k<orders[i];++k)
			parameters.insert(i);
	return fderivative(serial,parameters, args);
}

}


