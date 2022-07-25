#include "wedge/polynomialalgebra/polybasis.h"

namespace Wedge {

namespace internal {

ex SimplifyPowerEqn(ex x)
{
	return is_a<power>(x)? x.op(0) : x;
}

ex SimplifyPolyEqn(ex x,lst variables)
{
	x=sqrfree(x.expand(),variables);	//bug in ginac requires expand before sqrfree??
	if (is_a<mul>(x)) {
		exvector ops;
		for (ex y : x) ops.push_back(SimplifyPowerEqn(y));
		return mul(ops);
	}
	else return SimplifyPowerEqn(x);
}

}
}
