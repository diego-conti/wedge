#include "pseudolevicivita.h"

namespace Wedge {


class CovariantDerivativeSpinor: public IBilinearOperator<LinearOperator<VectorField>, Leibniz<Spinor,Function> > {
    const PseudoRiemannianStructureByOrthonormalFrame& g;
	const PseudoLeviCivitaConnection& connection;
    ex double_clifford(ZeroBased i, ZeroBased j, ex psi) const {
        ex e_i=g.e()[i];
        ex e_jsharp=g.ScalarProduct().Sharp(g.e().dual()[j]);
        return g.CliffordDot(e_jsharp, g.CliffordDot(e_i,psi));
    }
public:
	CovariantDerivativeSpinor(const PseudoRiemannianStructureByOrthonormalFrame& g, const PseudoLeviCivitaConnection& c) : g{g}, connection(c) {}
	ex Apply(const VectorField& X, const Spinor& psi) const {
        ex res; 
        for (int i=0;i<g.e().size();++i)
            for (int j=i+1;j<g.e().size();++j) 
                res+=TrivialPairing<VectorField>(X,connection(i,j))*double_clifford(i,j,psi);
        return res/2;
    }
	ex Apply(const VectorField& vfield, const Function& f) const {
    	return connection.Nabla<DifferentialForm>(vfield,f);	
    }

    static CovariantDerivativeSpinor* Create(const PseudoRiemannianStructure& g, const PseudoLeviCivitaConnection& c) {
        auto s = dynamic_cast<const PseudoRiemannianStructureByOrthonormalFrame*>(&g); 
        return s? new CovariantDerivativeSpinor(*s, c) : nullptr;
    }
};

void PseudoLeviCivitaConnection::Deleter::operator() (CovariantDerivativeSpinor* p) {
    delete p;
}

PseudoLeviCivitaConnection::PseudoLeviCivitaConnection(const Manifold* manifold, const PseudoRiemannianStructure& structure, const Name& christoffel) : 
		Connection(manifold,structure.e(),true), 
		TorsionFreeConnection<true>(manifold,structure.e(),true,christoffel),
        covariant_derivative_spinor{CovariantDerivativeSpinor::Create(structure,*this)}
         {
    const int dimension=e().size();
    components.reserve(dimension);
    for (int i=0;i<dimension;i++)
    components.push_back(exvector(dimension));

    ExVector e_flat(dimension);	//i-th element represents ((e_i)^\flat))
    ExVector frame=structure.e().dual();
    for (int i=1;i<=dimension;++i)
    for (int j=1;j<=dimension;++j)
        e_flat(i)+=structure.ScalarProduct().OnVectors(frame(i),frame(j))*structure.e(j); 
    exvector de;	//i-th element represents d((e_i)^\flat))
    de.reserve(dimension);
    for (int i=1;i<=dimension;++i)
        de.push_back(manifold->d(e_flat(i)).normal());
    for (int i=0;i<dimension;++i)
    for (int j=0;j<dimension;++j)
    for (int k=0;k<dimension;++k)
    {

    ex X=frame[i];
    ex Y=frame[j]; 
    ex Z=frame[k];
    ex XYZ=TrivialPairing<DifferentialForm>(X*Y,-de[k])+
            TrivialPairing<DifferentialForm>(Z*X,-de[j])+
            TrivialPairing<DifferentialForm>(Z*Y,-de[i]);
    for (int h=0;h<dimension;++h)
        (*this)(h,j)+=(e()[i]*XYZ/2*structure.ScalarProduct().OnOneForms(e()[h],e()[k])).expand();
    }
}

template<>
ex PseudoLeviCivitaConnection::Nabla<Spinor>(ex X, ex psi) const {
    return CovariantDerivativeSpinor::BilinearOperator(X,psi,covariant_derivative_spinor.get());
}

}
