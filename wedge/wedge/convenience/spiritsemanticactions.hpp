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

#ifndef SEMANTIC_ACTIONS_H
#define SEMANTIC_ACTIONS_H
#include <boost/spirit/home/x3.hpp>
#include <boost/config/warning_disable.hpp>
#include <boost/fusion/include/at.hpp>



namespace x3 = boost::spirit::x3;
namespace ascii = boost::spirit::x3::ascii;
namespace fusion = boost::fusion;
namespace Wedge {


namespace SemanticActions {
	static auto read_zero=[](auto& ctx) {x3::_val(ctx)=0;};
	static auto read_opposite=[](auto& ctx) {x3::_val(ctx)=-x3::_attr(ctx);};

	static auto read_product=[](auto& ctx) {
		ex op1 = fusion::at_c<0>(x3::_attr(ctx));
		ex op2 = fusion::at_c<1>(x3::_attr(ctx));
		x3::_val(ctx)=op1*op2;
	};
	static auto read_mul=[](auto& ctx) {
		static_assert(is_convertible<decltype(x3::_attr(ctx)),exvector>::value,"synthesized attribute should be exvector");
		auto& factors = x3::_attr(ctx);
		x3::_val(ctx)= mul(factors);
	};
	static auto read_wedge_product=[](auto& ctx) {
		static_assert(is_convertible<decltype(x3::_attr(ctx)),exvector>::value,"synthesized attribute should be exvector");
		//DisplayType<decltype(x3::_attr(ctx))> type_display;
		exvector factors (x3::_attr(ctx));
		x3::_val(ctx)=ncmul(factors);
	};
	static auto read_quotient =[](auto& ctx) {
		ex numer = fusion::at_c<0>(x3::_attr(ctx));
		ex denom = fusion::at_c<1>(x3::_attr(ctx));
		x3::_val(ctx)=numer/denom;
	};
	static auto read_sum=[](auto& ctx) {
		static_assert(is_convertible<decltype(x3::_attr(ctx)),exvector>::value,"synthesized attribute should be exvector");
		exvector summands (x3::_attr(ctx));
		x3::_val(ctx)=add(summands);
		//DisplayType<decltype(x3::_attr(ctx))> type_display;
	};

	static auto read_power = [](auto& ctx) {
		auto op1= fusion::at_c<0> (x3::_attr(ctx));
		auto op2= fusion::at_c<1> (x3::_attr(ctx));
		x3::_val(ctx)= pow(op1,op2);
	};
	static auto read_equation = [] (auto& ctx) {
		ex lhs= fusion::at_c<0> (x3::_attr(ctx));
		ex rhs= fusion::at_c<1> (x3::_attr(ctx));
		x3::_val(ctx)= (lhs==rhs);
	};

}
}
#endif
