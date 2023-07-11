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
#include "../structures/spinor.h"

#include "../structures/riemannianstructure.h"

namespace Wedge {

Spinor::Spinor(vector<bool> index)
{
	this->a=index;
}

int Spinor::compare_same_type(const basic &other) const
{
	const Spinor &o = static_cast<const Spinor &>(other);
//it is assumed that spinors appearing in the same expression belong to the same manifold. Here, we only check that the dimensions are consistent.
	assert(a.size()==o.a.size());
	for (int i=a.size()-1;i>=0;i--)
		if (a[i] && !o.a[i]) return 1;
		else if (o.a[i] && !a[i]) return -1;
	return 0;
}

void Spinor::print(const print_context &c, unsigned level) const {
	int index=0;
	for (auto i=a.rbegin();i!=a.rend();i++) index=*i? index*2+1:index*2;
	if (dynamic_cast<const GiNaC::print_latex*>(&c))
		c.s<<"u_{"<<index<<"}";
	else
		c.s<<"u"<<index;
}

Spinor Spinor::from_epsilons(const vector<int>& signs) {
	vector<bool> as_bools;
	transform(signs.begin(),signs.end(), back_inserter(as_bools), [] (int x) {return x<0;});
	return Spinor{as_bools};
}

/**
   @brief Returns the n-th element of a global basis of complex spinors
   @param n An index in the range [0,2^[dimension/2])
   @param dimension The dimension of the manifold
   @return The spinor \f$ u(\epsilon_m,\dots,\epsilon_1)\f$ where m=[dimension/2] and \epsilon_i=1 if the i-th least significant digit in base 2 of n is 0 and -1 otherwise
*/
Spinor Spinor::from_index_and_dimension(ZeroBased n, int dimension) {
	auto oldn=n;
	int m=dimension/2;
	if (n<0) throw OutOfRange(__FILE__,__LINE__,n);
	vector<int> signs;
	signs.reserve(m);
	for (int i=0;i<m;++i, n/=2)
		signs.push_back((n%2)? -1 : 1);
	if (n!=0) throw OutOfRange(__FILE__,__LINE__,oldn);
	return Spinor::from_epsilons(signs);
}

}
