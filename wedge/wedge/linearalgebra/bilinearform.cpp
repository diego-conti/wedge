#include "bilinearform.h"

namespace Wedge {

const ex AsDifferentialForm(ex x)
{
	if (is_a<DifferentialForm>(x)) return x;
	else if (is_a<DifferentialOneForm>(x)) return lst{x};
	else throw InvalidArgument(__FILE__,__LINE__,x);
}

ex BilinearForm::OnForms(ex v, ex w) const {
	v=v.expand(); w=w.expand();
	internal::NormalFormHelper<DifferentialForm> visit_v;
	v.accept(visit_v);
	internal::NormalFormHelper<DifferentialForm> visit_w;
	w.accept(visit_w);
	ex result;
	LOG_INFO(v);
	LOG_INFO(w);
	for (exmap::const_iterator i=visit_v.coeffs.begin();i!=visit_v.coeffs.end();++i)
	for (exmap::const_iterator j=visit_w.coeffs.begin();j!=visit_w.coeffs.end();++j)
	{
		result+=i->second * j->second * OnSimpleForms(AsDifferentialForm(i->first),AsDifferentialForm(j->first));
	}
	LOG_INFO(result);
	return result.expand();
}

ex BilinearForm::OnSimpleForms(ex v, ex w) const {
	if (v.nops()!=w.nops()) return 0;
	else {
		matrix m(w.nops(),w.nops());
		for (int i=0;i<w.nops();i++)
			for (int j=0;j<w.nops();j++)
				m(i,j)=	OnOneForms(v.op(i),w.op(j));
		return m.determinant();
	}
}


ex BilinearForm::Interior(ex v, ex w) const {
	v=Sharp(v).expand(); w=w.expand();
	internal::NormalFormHelper<DifferentialForm> visit_v;
	v.accept(visit_v);
	internal::NormalFormHelper<DifferentialForm> visit_w;
	w.accept(visit_w);
	ex result;
	for (exmap::const_iterator i=visit_v.coeffs.begin();i!=visit_v.coeffs.end();++i)
	for (exmap::const_iterator j=visit_w.coeffs.begin();j!=visit_w.coeffs.end();++j)
	{
		result+=i->second * j->second * Hook(i->first,j->first);
	}
	return result;
}

}
