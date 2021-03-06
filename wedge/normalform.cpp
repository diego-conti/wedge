#include "normalform.h"
#include "ginac/ginac.h"
#include "expressions.h"

namespace Wedge {
using namespace std;

template<bool ALLOW_SCALARS>
struct NormalFormHelper : public visitor,public Vector::visitor, public LambdaVector::visitor,public add::visitor, public mul::visitor, public basic::visitor {
	ex scalar_part;
	exmap& vector_part;
	exmap& lambda_vector_part;

	NormalFormHelper(exmap& vector_coeffs,exmap& lambda_vector_coeffs) : vector_part{vector_coeffs}, lambda_vector_part{lambda_vector_coeffs} {}
	void visit(const Vector& x) {
		++vector_part[x];
	}
	void visit(const LambdaVector& x) {
		++lambda_vector_part[x];
	}

	void SetScalarPart(ex x) {
		if (ALLOW_SCALARS) scalar_part+=x;
		else {
			LOG_ERROR(x);
			throw WedgeException<std::runtime_error>("Expression contains a scalar component, but linear combination of vectors was expected",__FILE__,__LINE__);
		}		
	}

	void visit(const mul& x)
	{
		exvector ops;
		ex the_vector;
		for (int i=0;i<x.nops();++i)
		{
			ex opi=x.op(i);
			if (is_a<VectorBase>(opi)) {
				if (!the_vector.is_zero()) {
					LOG_ERROR(x);
					LOG_ERROR(ops);
					throw WedgeException<std::runtime_error>("Linear combination of Vector's or LambdaVector's expected",__FILE__,__LINE__);
				}
				the_vector=opi;
			}
			else ops.push_back(opi); 
		}
		if (the_vector.is_zero()) SetScalarPart(x);
		else if (is_a<LambdaVector>(the_vector)) lambda_vector_part[the_vector]+=mul(ops);
		else vector_part[the_vector]+=mul(ops);
	}
	
	void visit(const add& x)
	{
		for (int i=0;i<x.nops();++i)
			x.op(i).accept(*this);	
	}
	void visit (const basic& x)
	{
		if (x!=0) SetScalarPart(x);
	}
};

VectorNormalForm::VectorNormalForm(ex linear_combination) {
	NormalFormHelper<false> visitor{coefficients, coefficients};
	linear_combination.expand().accept(visitor);
}

LambdaVectorNormalForm::LambdaVectorNormalForm(ex linear_combination) {
	NormalFormHelper<true> visitor{vector_part.coefficients, lambda_vector_part.coefficients};
	linear_combination.expand().accept(visitor);
	scalar_part=visitor.scalar_part;
}


ex VectorNormalForm::CollectCoefficients() const {
	exmap inverse;
	for (const auto& pair : coefficients) {
		ex coeff=pair.second.normal(); 
		if (inverse.find(-coeff)!=inverse.end()) 
			inverse[-coeff]-=pair.first;
		else
			inverse[coeff]+=pair.first;;
	}
	ex result;
	for (auto& pair : inverse)
		result+=pair.first*pair.second;
	return result;
}

}
