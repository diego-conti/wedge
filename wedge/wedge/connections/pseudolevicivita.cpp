#include "pseudolevicivita.h"

namespace Wedge {


class CovariantDerivativeSpinor: public IBilinearOperator<LinearOperator<VectorField>, Leibniz<Spinor,Function> > {
    const PseudoRiemannianStructureByOrthonormalFrame& g;
	const PseudoLeviCivitaConnection& connection;
    ex double_clifford(ZeroBased i, ZeroBased j, ex psi) const {
        ex e_i=g.e().dual()[i];
        ex e_jsharp=g.ScalarProduct().Sharp(g.e()[j]);
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

    //uses the formula \nabla_{e_i}e_j = \sum_k \langle \nabla_{e_i}e_j, (e^k)^\sharp\rangle e_k
    //\omega_kj = \langle \nabla_{e_i}e_j, (e^k)^\sharp\rangle e^i

    ExVector e_flat(dimension);	//i-th element represents ((e_i)^\flat))
    ExVector frame=structure.e().dual();
    for (int i=1;i<=dimension;++i)
        e_flat(i)=structure.ScalarProduct().Flat(frame(i));
    exvector de;	//i-th element represents d((e_i)^\flat))
    de.reserve(dimension);
    for (int i=1;i<=dimension;++i)
        de.push_back(manifold->d(e_flat(i)).normal());
    for (int k=0;k<dimension;++k)  {
        ex ek=structure.e()[k];
        ex Z=structure.ScalarProduct().Sharp(ek);
        for (int i=0;i<dimension;++i)
        for (int j=0;j<dimension;++j) {
            ex X=frame[i];
            ex Y=frame[j];         
            ex XYZ=TrivialPairing<DifferentialForm>(X*Y,-manifold->d(ek))+
                TrivialPairing<DifferentialForm>(Z*X,-de[j])+
                TrivialPairing<DifferentialForm>(Z*Y,-de[i]); 
            (*this)(k,j)+=structure.e()[i]*XYZ/2;
        }
    }        
}

template<>
ex PseudoLeviCivitaConnection::Nabla<Spinor>(ex X, ex psi) const {
    return CovariantDerivativeSpinor::BilinearOperator(X,psi,covariant_derivative_spinor.get());
}

}
