/***************************************************************************
 *   Copyright (C) 2009  by Diego Conti					   *
 *   diego.conti@unimib.it                                                 *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 3 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#ifndef COHERENTSPLITTINGS_H
#define COHERENTSPLITTINGS_H 
#include "wedge/wedge.h"

namespace Wedge {

//Compute the intersection of two vector spaces
template<typename T> VectorSpace<T> Intersect(const VSpace<T>& a, const VSpace<T>& b)
{
	SubBasis<T> basis(a.e_begin(),a.e_end(),a.e_end(),a.e_end());
	basis.AddToComplement(b.e_begin(),b.e_end());
	ExVector comps=basis.AllComponents(b.GenericElement());
	return b.SubspaceFromEquations(comps.begin()+basis.size(),comps.end());
}

/** @brief A class that represents a coherent splitting \f$\mathfrak{g}=V_1\oplus V_2\f$ viewed as a product structure
 *
 * The "coherent" condition means that in terms of the doubly graded vector space \f$\Lambda^{p,q}=\Lambda^pV_1\otimes\Lambda^qV_2\f$, the exterior derivative d has only two bigraded components, e.g. (2,-1) + (1,1), effectively making it into a double complex
 * @note This class does not allow either the structure or the d operator to depend on parameters
*/
class pqStructure : public GStructure {
	int K; 			//dimension of V1
	int delta1p, delta2p;	//delta1 is the (delta1p, 1-delta1p) component of d
	class dOperator : public Derivation<DifferentialForm> {
	public:
		lst subs;
		void visit(const VectorField& alpha) {
			Result()=ex(alpha).subs(subs);
		}
	};
	mutable dOperator _delta1,_delta2;

public:
/** @brief Construct a pqStructure
 * @param manifold The manifold on which the structure is defined. Should be a LieGroupWithNoParameters.
 * @param adaptedframe A frame adapted to the splitting, meaning that the first k elements span \f$V_1\f$ and the rest span \f$V_2\f$
 * @param k The dimension of \f$V_1\f$
 * @param d1,d2 Integers between -1 and 2 that determine the bidegree of \f$\delta_1\f$ and \f$\delta_2\f$
 *
 * Assumes the following decomposition holds: \f$d=\delta_1+\delta_2\f$, where $\delta_i\Lambda^{p,q}\to \Lambda^{p+d_i,q+1-d_i}\f$.
 * @throw InvalidArgument if the decomposition does not hold
*/
	pqStructure(const Manifold* manifold, const Frame& adaptedframe, int k, int d1=1, int d2=2) : 
		GStructure(manifold,adaptedframe), hnk(manifold->Dimension()+1) 
	{
		if (k<1 || k>manifold->Dimension())  throw OutOfRange(__FILE__,__LINE__,k); 
		if (d1<-1 || d1>2) throw OutOfRange(__FILE__,__LINE__,d1);
		if (d2<-1 || d2>2) throw OutOfRange(__FILE__,__LINE__,d2);

		K=k;
		delta1p=d1;
		delta2p=d2;

		//\Lambda^1,0 goes to \Lambda^{1+d,1-d}
		//\Lambda^0,1 goes to \Lambda^{d,2-d}

		vector<Subspace<DifferentialForm> > L(3);
		L[0]=pqFormsAsSubspace(0,2);
		L[1]=pqFormsAsSubspace(1,1);
		L[2]=pqFormsAsSubspace(2,0);

		ExVector delta1e, delta2e;	
		for (Frame::const_iterator i=adaptedframe.begin();i!=adaptedframe.begin()+k;++i)
		{
			ex dei=manifold->d(*i);
			ex delta1ei,delta2ei;
			if (1+delta1p<=2) delta1ei=L[1+delta1p].Project(dei);			
			if (1+delta2p<=2) delta2ei=L[1+delta2p].Project(dei);
			delta1e.push_back(delta1ei);
			delta2e.push_back(delta2ei);
			if (dei!=delta1ei+delta2ei)
			{
				LOG_ERROR(adaptedframe);
				LOG_ERROR(k);
				LOG_ERROR(dei);
				LOG_ERROR(*i);
				LOG_ERROR(delta1ei);
				LOG_ERROR(delta2ei);
				throw InvalidArgument(__FILE__,__LINE__,*i);
			}
		}
		for (Frame::const_iterator i=adaptedframe.begin()+k;i!=adaptedframe.end();++i)
		{
			ex dei=manifold->d(*i);
			ex delta1ei,delta2ei;
			if (delta1p>=0) delta1ei=L[delta1p].Project(dei);			
			if (delta2p>=0) delta2ei=L[delta2p].Project(dei);
			delta1e.push_back(delta1ei);
			delta2e.push_back(delta2ei);
			if (dei!=delta1ei+delta2ei)
			{
				LOG_ERROR(adaptedframe);
				LOG_ERROR(k);
				LOG_ERROR(dei);
				LOG_ERROR(*i);
				LOG_ERROR(delta1ei);
				LOG_ERROR(delta2ei);
				throw InvalidArgument(__FILE__,__LINE__,*i);
			}
		}
		
		_delta1.subs=LinearMapToSubstitutions<DifferentialForm>(adaptedframe,delta1e);
		_delta2.subs=LinearMapToSubstitutions<DifferentialForm>(adaptedframe,delta2e);

		LOG_INFO(_delta1.subs);
		LOG_INFO(_delta2.subs);
		for (Frame::const_iterator i=adaptedframe.begin();i!=adaptedframe.end();++i)
		{
			LOG_INFO(M()->d(*i));
			LOG_INFO(delta1(*i));
			LOG_INFO(delta2(*i));
			LOG_INFO(delta1(delta1(*i)));
			LOG_INFO(delta2(delta2(*i)));
		}
/*
		for (int i=1;i<=M()->Dimension();++i)
		{
			VectorSpace<DifferentialForm> V=M()->pForms(i);
			for (exvector::const_iterator i=V.e_begin();i!=V.e_end();++i) {
				cout<<*i<<" -> "<<delta1(*i)<<" -- " << delta2(*i)<<endl;
				cout<<"      -> "<<delta1(delta1(*i))<<" -- " << delta2(delta2(*i))<<endl;
			}
		}
*/
		for (int i=0;i<=K;++i) {
			ex x=pqForms(i,K-i).GenericElement();
			if (M()->d(x).expand()!=(delta1(x)+delta2(x)).expand())
			{
				LOG_ERROR(adaptedframe);
				LOG_ERROR(k);
				LOG_ERROR(x);
				LOG_ERROR(M()->d(x));
				LOG_ERROR(delta1(x));
				LOG_ERROR(delta2(x));
				throw InvalidArgument(__FILE__,__LINE__,x);
			}
			else if (!delta1(delta1(x)).expand().is_zero())
			{
				LOG_ERROR(adaptedframe);
				LOG_ERROR(k);
				LOG_ERROR(x);
				LOG_ERROR(M()->d(x));
				LOG_ERROR(delta1(x).expand());
				LOG_ERROR(delta1(delta1(x)).expand());
				throw InvalidArgument(__FILE__,__LINE__,x);
			}
			else if (!delta2(delta2(x)).expand().is_zero())
			{
				LOG_ERROR(adaptedframe);
				LOG_ERROR(k);
				LOG_ERROR(x);
				LOG_ERROR(M()->d(x));
				LOG_ERROR(delta2(x).expand());
				LOG_ERROR(delta2(delta2(x)).expand());
				throw InvalidArgument(__FILE__,__LINE__,x);
			}
		}
	}

/** @brief Compute \f$\delta_1(\alpha)\f$
*/
	ex delta1(ex alpha) const	
	{
		return _delta1.RecursiveVisit(alpha.expand());
	}

/** @brief Compute \f$\delta_2(\alpha)\f$
*/
	ex delta2(ex alpha) const
	{
		return _delta2.RecursiveVisit(alpha.expand());
	}

/** @brief Compute the space of \delta_1-closed (p,q)-forms.
*/
	VectorSpace<DifferentialForm> d1closedpqforms(int p,int q) const
	{
		VectorSpace<DifferentialForm> forms=pqForms(p,q);
		lst eqns;
		GetCoefficients<DifferentialForm>(eqns,delta1(forms.GenericElement()));
		return forms.SubspaceFromEquations(eqns.begin(),eqns.end());
	}

/** @brief Compute the space of \delta_1-exact (p,q)-forms.
*/
	VectorSpace<DifferentialForm> d1exactpqforms(int p,int q) const
	{
		VectorSpace<DifferentialForm> result;
		VectorSpace<DifferentialForm> forms=pqForms(p-delta1p,q-(1-delta1p));
		for (exvector::const_iterator i=forms.e_begin();i!=forms.e_end();++i)
			result.AddGenerator(delta1(*i));
		return result;
	}
/** @brief Compute the dimension of the double cohomology space H^{p,q}
*/
	int h(int p, int q) const
	{
		if (p<0 || q<0 || p+q>M()->Dimension()) return 0;
		if (hnk[p+q].size()==0) Computehnk(p+q);
		int h1=hnk[p+q][p],h2;
		if (q>0) h2=hnk[p+q][p+1];
		else h2=0;
		if (h1<h2) {
			LOG_ERROR(h1);
			LOG_ERROR(h2);
		}
		return h1-h2;
	}

/** @brief Compute the dimension of the first complex in the spectral sequence.
*/
	int E1(int p, int q) const
	{
		return d1closedpqforms(p,q).Dimension()-d1exactpqforms(p,q).Dimension();
	}

/** @brief Compute the space of the forms of type (p,q), or zero if p=q=0
 *
 * In principle, \f$\Lambda^{00}=\mathbb{R}\f$, but setting it to zero has no effect on cohomology
*/
	VectorSpace<DifferentialForm> pqForms(int p,int q) const
	{
		VectorSpace<DifferentialForm> lambda0q;
		if (p<0 || q<0 || p>K || q>M()->Dimension()-K) return lambda0q;
		Frame V1,V2;
		for (int i=K+1;i<=M()->Dimension();++i)
			V2.push_back(e(i));
		for (int i=1;i<=K;++i)
			V1.push_back(e(i));
		if (q==0) {
			if (p==0) return lambda0q;
			else return Wedge::pForms(V1, p);
		}
		else {
			lambda0q=Wedge::pForms(V2, q);
			if (p==0) return lambda0q;
			VectorSpace<DifferentialForm> lambdap0=Wedge::pForms(V1, p);

			VectorSpace<DifferentialForm> result;
			for (exvector::const_iterator i=lambda0q.e_begin();i!=lambda0q.e_end();++i)
			for (exvector::const_iterator j=lambdap0.e_begin();j!=lambdap0.e_end();++j)
				result.AddGenerator(*j * *i);
			return result;
		}
	}
	

private:	
	//compute the space of the forms of type (p,q), or zero if p=q=0
	Subspace<DifferentialForm> pqFormsAsSubspace(int p,int q) const
	{
		VectorSpace<DifferentialForm> V=pqForms	(p,q);
		return M()->pForms(p+q).Subspace(V.e_begin(),V.e_end());
	}

	mutable vector<vector<int> > hnk;	// h[n][k] is the number of independent non-exact closed forms in lambda^{k,n-k}+lambda^{k+1,n-k-1}+\dots + \lambda^{n,0}
	void Computehnk(int n) const
	{
		if (n<=0 || n>M()->Dimension()) throw OutOfRange(__FILE__,__LINE__,n);
		hnk[n]=vector<int>(n+1);	
		const LieGroupHasParameters<false>* G = dynamic_cast<const LieGroupHasParameters<false>*> (M());
		if (G==NULL) throw NotImplemented(__FILE__,__LINE__);

		VectorSpace<DifferentialForm> Fpq;	//compute lambda^{p,q}+lambda^{p+1,q-1}+\dots + \lambda^{p+q,0}
		for (int p=n;p>=0;--p)
		{
			VectorSpace<DifferentialForm> lpq=pqForms(p,n-p);
			Fpq.AddGenerators(lpq.e_begin(),lpq.e_end());
			if (Fpq.Dimension()==0) {
				hnk[n][p]=0;
				continue;
			}

			list<ex> eqns;
			GetCoefficients<DifferentialForm>(eqns,M()->d(Fpq.GenericElement()));
			int z=Fpq.SubspaceFromEquations(eqns.begin(),eqns.end()).Dimension();
			
			VectorSpace<DifferentialForm> B=Intersect(Fpq,G->ExactForms(n));
			int b=B.Dimension();
			LOG_INFO(z<<" "<<b<<" "<<z-b);
			hnk[n][p]=z-b;
		}

	}

};


WEDGE_DECLARE_NAMED_ALGEBRAIC(SplittingParameter,realsymbol)

/** @brief A class that represents a coherent splitting \f$\mathfrak{g}=V_1\oplus V_2\f$
 *
 * A coherent splitting on a Lie algebra is represented by a decomposable k-form \f$\phi\Lambda^{k,0}\f$.
 *
 * Both the algebra and the splitting may have parameters.
*/
class CoherentSplitting : public HasParameters<SplittingParameter> {
	const LieGroup* M;
	ex definingForm;
	int k;
public:
/** @brief Construct the generic splitting on a given Lie algebra */
	CoherentSplitting(const LieGroupHasParameters<false>* algebra, int k)
	{
		M=algebra;
		this->k=k;
		if (k<=0 || k>=M->Dimension()) throw OutOfRange(__FILE__,__LINE__,k);
		VectorSpace<DifferentialForm> kforms=M->pForms(k);

		list<ex> eqns;
		for (int i=1;i<=M->Dimension();++i)
			GetCoefficients<DifferentialForm>(eqns, M->d(M->e(i))*kforms.GenericElement());
		GetCoefficients<DifferentialForm>(eqns, M->d(kforms.GenericElement()));

		Subspace<DifferentialForm> V=kforms.SubspaceFromEquations(eqns.begin(),eqns.end());
		for (int i=1;i<=V.Dimension();++i)
			definingForm+=SplittingParameter("a"+ToString(i))*V.e(i);
	}
/** @brief Construct the generic splitting on a given Lie algebra */
	CoherentSplitting(const LieGroupHasParameters<true>* algebra, int k)
	{
		M=algebra;
		this->k=k;
		if (k<=0 || k>=M->Dimension()) throw OutOfRange(__FILE__,__LINE__,k);
		VectorSpace<DifferentialForm> kforms=M->pForms(k);
		for (int i=1;i<=kforms.Dimension();++i)
			definingForm+=SplittingParameter("a"+ToString(i))*kforms.e(i);
	}

/** @brief Construct a fixed splitting on a given Lie algebra */
	CoherentSplitting(const LieGroup* algebra, ex definingForm) 
	{
		M=algebra;
		this->definingForm=definingForm;
		k=Degree<DifferentialForm>(definingForm);
	}


	template<typename Container> Container& GetDeclareCoherent(Container& container) 
	{
		bool has_simplified;
		set<ex,ex_is_less> eqns;
		do {
			has_simplified=false;
			eqns.clear();
			GetCoherentEquations(eqns);
			eqns.erase(0);
			set<ex,ex_is_less>::iterator i=eqns.begin();
			while (i!=eqns.end()) {
				exmap occ;
				ex w=wild(), y=wild();
				set<ex,ex_is_less>::iterator j=i;
				++i;
				if (j->match(pow(w,2),occ) ||
					j->match(-pow(w,2),occ) ||
						(j->match(y*pow(w,2),occ) && is_a<numeric>(occ[y]))) 
				{
					LOG_INFO(*j);
					LOG_INFO(occ);
					LOG_INFO(occ[w]);
					eqns.insert(occ[w]);					
					eqns.erase(j);
					has_simplified=true;
				}
			}
			list<ex> linear;
			for (set<ex,ex_is_less>::const_iterator i=eqns.begin();i!=eqns.end();++i)
				if (Degree<SplittingParameter>(*i)<=1)
					linear.push_back(*i);
			try {
				DeclareZero(linear.begin(),linear.end());
			}
			catch (InconsistentDeclaration&)
			{
				Insert(container,ex(1));
				return container;
			}
			if (!linear.empty()) has_simplified=true;
		} while (has_simplified);
		Insert(container,eqns.begin(),eqns.end());
		return container;
	}

/** @brief Compute the equations in the splitting parameters that make the splitting coherent
 *
 * Impose the conditions: \f$d\phi=0\f$, \f$d(\alpha)\wedge\phi=0\f$, \f$L_X\phi\in\mathbb{R}\phi\f$.
 * Does not check that the form is decomposable unless k=2
 */
	template<typename Container> Container& GetCoherentEquations(Container& container) const
	{
		//these will be trivial if the group does not have parameters 
		for (int i=1;i<=M->Dimension();++i)
			GetCoefficients<DifferentialForm>(container, M->d(M->e(i))*definingForm);
		GetCoefficients<DifferentialForm>(container, M->d(definingForm));

		//impose that L_x\phi be a multiple of \phi for all X
		VectorSpace<DifferentialForm> kforms=M->pForms(k);
		ExVector components=kforms.Components(definingForm);
		for (int i=1;i<=M->Dimension();++i)
		{
			ExVector components2=kforms.Components(M->LieDerivative(M->e(i),definingForm));
			for (int i=1;i<=components.size();++i)
			for (int j=i+1;j<=components.size();++j)
				Insert(container, (components(i)*components2(j)-components(j)*components2(i)).expand());
		}
		
		//impose that the form is decomposable 
		if (k==2) 
			GetCoefficients<DifferentialForm>(container,(definingForm*definingForm).expand());
		else if (k==3)	//TODO I suspect this is not a sufficient condition
			for (int i=1;i<=M->Dimension();++i)
			{
				ex a=Hook(M->e(i),definingForm);
				GetCoefficients<DifferentialForm>(container,(a*a).expand());
			}
		else LOG_WARN(k);
		return container;
	}

/** @brief Returns true if the splitting is coherent (for all values of the parameters)
*/
	bool IsCoherent() const {
		set<ex,ex_is_less> eqns;
		GetCoherentEquations(eqns);
		eqns.erase(0);
		if (!eqns.empty()) {
			LOG_WARN(eqns);
			return false;
		}
		return true;
	}
/** @brief Impose conditions on the splitting that make \f$H^{0,q}=0\f$
 *
 * @throw NotImplemented if the group depends on parameters
*/
	void DeclareH0qZero(int q)
	{
		const LieGroupHasParameters<false>* N=dynamic_cast<const LieGroupHasParameters<false>* >(M);
		if (N==NULL) throw NotImplemented(__FILE__,__LINE__,"DeclareH0qZero for groups depending on parameters");
		Subspace<DifferentialForm> closedforms=	N->ClosedForms(q);
		list<ex> eqns;
		for (exvector::const_iterator i=closedforms.e_begin();i!=closedforms.e_end();++i)
			GetCoefficients<DifferentialForm>(eqns,(definingForm * *i).expand());
		DeclareZero(eqns.begin(),eqns.end());
	}

	ex DefiningForm() const {return definingForm;}

/** @brief Return a pqstructure object associated to this coherent splitting
 *
 * @note This will not work if the splitting depends on parameters
*/
	pqStructure Structure() const {
		VectorSpace<DifferentialForm> TM=M->pForms(1);
		list<ex> eqns;
		GetCoefficients<DifferentialForm>(eqns,TM.GenericElement() * definingForm);
		Subspace<DifferentialForm> V=TM.SubspaceFromEquations(eqns.begin(),eqns.end());
		Frame adaptedframe(V.full_begin(),V.full_end());
		return pqStructure(M, adaptedframe, k,1, 2);
	}

/** @brief Return an "inverse" pqstructure object associated to this coherent splitting
 *
 * This is the pqStructure object whose double complex is obtained from the natural pqStructure's complex by swapping the rows and the columns
 *
 * @note This will not work if the splitting depends on parameters
*/
	pqStructure InverseStructure() const {
		VectorSpace<DifferentialForm> TM=M->pForms(1);
		list<ex> eqns;
		GetCoefficients<DifferentialForm>(eqns,TM.GenericElement() * definingForm);
		Subspace<DifferentialForm> V=TM.SubspaceFromEquations(eqns.begin(),eqns.end());
		Frame adaptedframe(V.complement_begin(),V.complement_end());
		adaptedframe.insert(adaptedframe.end(),V.e_begin(),V.e_end());
		return pqStructure(M, adaptedframe, M->Dimension()-k,-1, 0);
	}

	//construct a matrix whose columns are the vectors
	template<typename Iterator1, typename Iterator2> matrix MatrixFromVectors(Iterator1 basis_begin, Iterator1 basis_end, Iterator2 v_begin, Iterator2 v_end) const
	{
		exvector basis(basis_begin,basis_end);
		exvector v;
		for (Iterator2 i=v_begin;i!=v_end;++i)
			v.push_back(i->expand());
		if (basis.empty() || v.empty()) return matrix();
		matrix M(basis.size(),v.size());
		for (int i=0;i<basis.size();++i)
		for (int j=0;j<v.size();++j)
			M(i,j)=v[j].coeff(basis[i]);
		return M;
	}

/** @brief Print out the conditions that ensure that the cohomology group \f$H^{0,q}\f$ is zero
*/
	void PrintEquationsH0qZero(int q) const
	{
		VectorSpace<DifferentialForm> qforms=M->pForms(q);
		set<ex,ex_is_less> eqns;
		GetCoefficients<DifferentialForm>(eqns,qforms.GenericElement()*definingForm);
		set<ex,ex_is_less> eqns2;
		GetCoefficients<DifferentialForm>(eqns2,M->d(qforms.GenericElement()));
		cout<<"$H^{0,"<<q<<"}=0$ iff"<<endl;
		list<ex> coordinates;
		GetSymbols<VectorSpace<DifferentialForm>::Coordinate>(coordinates, eqns.begin(),eqns.end());
		GetSymbols<VectorSpace<DifferentialForm>::Coordinate>(coordinates, eqns2.begin(),eqns2.end());
		matrix L=MatrixFromVectors(coordinates.begin(),coordinates.end(), eqns.begin(),eqns.end());
		matrix D=MatrixFromVectors(coordinates.begin(),coordinates.end(), eqns2.begin(),eqns2.end());
		vector<bool> eliminaterow(D.rows());
		for (int i=0;i<eliminaterow.size();++i) eliminaterow[i]=false;

		for (int i=0;i<D.cols();++i) {
			int j=0; 
			while (j<D.rows() && (D(j,i).is_zero())) ++j;
			int k=j+1; 
			while (k<D.rows() && (D(k,i).is_zero())) ++k;
			if (j<D.rows() && is_a<numeric>(D(j,i)) && k==D.rows())
				eliminaterow[j]=true;
		}
		int newrows=0;
		for (int i=0;i<eliminaterow.size();++i) if (!eliminaterow[i]) ++newrows;
		cout<<newrows<<endl;
		matrix L2(newrows,L.cols());
		matrix D2(newrows,D.cols());
		int row=0;
		for (int i=0;i<eliminaterow.size();++i) 
			if (!eliminaterow[i]) {
				for (int j=0;j<L.cols();++j) {
					L2(row,j)=L(i,j);
				for (int j=0;j<D.cols();++j)
					D2(row,j)=D(i,j);
				}
				++row;
			}
		cout<<"\\["<<L2<<"\\text{ depends on }"<< D2<<"\\]"<<endl;
	}
private:
	void DeclareConditions(const lst& list_of_equations) {
		definingForm=definingForm.subs(list_of_equations);
	}
};

}

#endif
