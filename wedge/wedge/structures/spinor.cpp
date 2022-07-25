/***************************************************************************
 *   Copyright (C) 2007, 2008, 2009 by Diego Conti, diego.conti@unimib.it  *
 *                                                                         *
 *   This file is part of Wedge.                                           *
 *                                                                         *
 *   Wedge is free software; you can redistribute it and/or modify         *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 3 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   Wedge is distributed in the hope that it will be useful,              *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with Wedge; if not, write to the                                *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include "../structures/spinor.h"

#include "../structures/riemannianstructure.h"

namespace Wedge {

Spinor::Spinor(unsigned index, int dimension)
{
	a.reserve(dimension);
	//convert index to binary representation; least significant digit goes first.
	for (int i=0;i<dimension;i++) {
		a.push_back((index%2)!=0);
		index/=2;
	}
	if (index!=0) throw OutOfRange(__FILE__,__LINE__);
}

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
	for (vector<bool>::const_reverse_iterator i=a.rbegin();i!=a.rend();i++) index=*i? index*2+1:index*2;
	if (dynamic_cast<const GiNaC::print_latex*>(&c))
		c.s<<"u_{"<<index<<"}";
	else
		c.s<<"u"<<index;
}

}
