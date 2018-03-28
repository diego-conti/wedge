/*
 * linearaction.cpp
 *
 *  Created on: Jan 8, 2015
 *      Author: diego
 */




#include "linearaction.h"

namespace Wedge {

ex LinearActionBase::eval_ncmul(const exvector & v) const
{

	exvector::const_reverse_iterator i=v.rbegin();
	assert(is_a<LinearActionBase>(*i));
	exvector product(ex_to<LinearActionBase>(*i).linear_subs.begin(),ex_to<LinearActionBase>(*i).linear_subs.end());
	while (++i!=v.rend())
		for (exvector::iterator k=product.begin();k!=product.end();++k)
		{
			ex rhs=k->rhs();
			rhs=rhs.subs(ex_to<LinearActionBase>(*i).linear_subs);
			k->let_op(1)=rhs;
		}
	return 	LinearActionBase (lst(product.begin(),product.end()));
}


ex GroupAction(ex g, ex w)
{
	internal::LinearActionVisitor v(w);
	return v.RecursiveVisit(g);
}


}
