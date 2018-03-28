/***************************************************************************
 *   Copyright (C) 2008, 2009 by Diego Conti				   *
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

#include <ginac/ginac.h>
#include <wedge/concretemanifold.h>
#include <wedge/affinebasis.h>
#include <wedge/liesubgroup.h>
#include <wedge/gl.h>
#include <wedge/su.h>
#include <sstream>
#include <stdlib.h>
#include <math.h>

using namespace Wedge;
using namespace  GiNaC;

/** @brief A class representing the bundle of frames over a manifold of a given dimension.
 */
 
class FrameBundle : public ConcreteManifold, public Has_dTable {	
public:
	int baseDimension;
	FrameBundle(int dimension) : ConcreteManifold(CreateFrame(baseDimension=dimension))
	{
		for (int i=1; i<=dimension;++i)
		{
			ex e;
			for (int j=1; j<=dimension;++j)
				e+=Theta(j)*Omega(i,j);
			Declare_d(Theta(i),e);
		}
	}
/** @brief The n-th component of the tautological form
*/
	ex Theta(OneBased n) const
	{
		if (n<1 || n>baseDimension) throw OutOfRange("FrameBundle::Theta",n);
		return e(n);
	}


/** @brief The (i,j) component of the flat connection form defined by the parallelism
*/
	ex Omega(OneBased i,OneBased j) const
	{
		if (i<1 || i>baseDimension) throw OutOfRange("FrameBundle::Omega",i);
		if (j<1 || j>baseDimension) throw OutOfRange("FrameBundle::Omega",j);
		return e(i*baseDimension+j);
	}

private:
	ExVector CreateFrame(int dimension)
	{
		ExVector frame(dimension*(dimension+1));
		for (int i=1; i<=dimension;++i)
		{	
			frame(i)=DifferentialOneForm("theta"+ToString(i));
			for (int j=1; j<=dimension;j++)
				frame(i*baseDimension+j)=DifferentialOneForm("omega"+ToString(i,j));
		}
		return frame;
	}
};


class EDSOnFrameBundle : public FrameBundle, public IterateOverPermutations
{
public:
	EDSOnFrameBundle(int dimensionOfBase) : FrameBundle(dimensionOfBase) {
		for (int i=1;i<=dimensionOfBase;++i)
			moduloIC.append(FrameBundle::e(i)==0);
		codim=-1; maxC=0;
		f=FrameBundle::e(); f.resize(baseDimension);
	}

	void ChooseFrame() {
		f.clear(); f.resize(baseDimension);
		cout<<"Choose a frame to determine a regular flag"<<endl;
		for (int i=1;i<=baseDimension;++i)
		{
			int k=0;
			while (k<1 || k>baseDimension)
			{
				cout<<"e_"<<i<<"="<<"? ";
				cin>>k;
			}
			f(i)=e(k);
		}
	}

/** @brief The n-th component of the tautological form

 Overridden, so as to allow for non-standard flags.
*/
	ex Theta(OneBased n) const
	{
		if (n<1 || n>baseDimension) throw OutOfRange("FrameBundle::Theta",n);
		return f(n);
	}
	
	ExVector Theta() const {return f;}

	void FindOrdinaryFlag(int startfrompermutationno=0) {
		if (startfrompermutationno==-1)	// randomize
			while (Iterate(baseDimension,random() % factorial(numeric(baseDimension)).to_int())) ;
		else
			Iterate(baseDimension,startfrompermutationno);
	}

	bool Apply(const std::vector<int>& permutation) {
		f.clear();
		for (std::vector<int>::const_iterator i=permutation.begin();i!=permutation.end();++i)
			f.push_back(e(*i+1));
		assert(f.size()==baseDimension);
		if (CartansTest(0)) {
			cout<<"Cartans test satisfied with tautological form "<<Theta()<<endl;
			CartansTest(1);
			return false;
		}
		return true;
	}

	virtual lst Ideal() const {return lst();}
	virtual string Geometry() const {return "";}

	bool CartansTest(int verbositylevel=1) const {
		if (verbositylevel>0) cout<<"Exterior differential system for "<<Geometry()<<" in dimension "<<baseDimension<<endl;
		
		lst I=this->Ideal();
		if (codim<0) {
			Basis<symbol> A;
			//a faster alternative than passing A to GetEquationsForVn
			set<ex,ex_is_less> v;
			GetEquationsForVn(v,I.begin(),I.end());
			LOG_INFO(v);
	  		A.insert(A.begin(),v.begin(),v.end());
			codim=A.size();
			LOG_INFO(codim);
			if (verbositylevel>1) cout<<"Variety of integral elements V_n:"<<endl<<A<<endl;	
		}
		if (verbositylevel>0) cout<<"Variety of integral elements V_n has codimension "<<codim<<endl;


		int C=0;
		for (int j=0;j<baseDimension;++j)
		{			
			Basis<DifferentialForm> V;
			for (lst::const_iterator i=I.begin();i!=I.end();++i)
				GetReducedPolarEquations(V,*i,j);
			if (verbositylevel>0) cout<<"The space of reduced polar equations for E_"<<j<<" has dimension "<<V.size()<<endl;
			if (verbositylevel>1) cout<<V<<endl;
			C+=V.size();
			
		}
		if (C>maxC) {
			maxC=C;
			LOG_INFO(C);
			LOG_INFO(Theta());
		}
		if (verbositylevel>0) cout<<"C="<<C<<"; Cartan's test "<< ((C==codim)?  "passed" :  "failed")<<endl;
		else LOG_INFO(C<<","<<codim);
		return C==codim;
	}

	template<typename Container, typename Iterator> Container& GetEquationsForVn(Container& container, Iterator begin, Iterator end) const
	{
		//coordinates on the grassmmannian: omega(i,j)=p(ijk)e(basis(k))
		lst substitutions;
		for (unsigned i=1; i<=baseDimension;++i)
		for (unsigned j=1; j<=baseDimension;++j)
		{
			ex ei;
			for (unsigned k=1;k<=baseDimension;++k)
				ei+=symbol("p_"+ToString(i,j,k))*e(k);
			substitutions.append(Omega(i,j)==ei);
		}
	 	for (Iterator i=begin; i!=end;++i)
			GetCoefficients(container,i->subs(substitutions).expand());
		return container;
	}

	template<typename Container> void GetReducedPolarEquations(Container& container, ex form, int j) const
	{
		if (Degree<DifferentialForm>(form)==1)
			container.push_back(form.subs(moduloIC));
		else if (j>0)
		{
			GetReducedPolarEquations(container, form,j-1);
			GetReducedPolarEquations(container, Hook(e(j),form),j-1);
		}
	}
	void AssumeStronglyAdmissible(int GroupDimension) const
	{
		if (codim<0) codim=baseDimension*((baseDimension*(baseDimension-1))/2-GroupDimension);
	}
	void AssumeCodimension(int c) {codim=c;}
private:
	//disallow access to e(), subclasses should use Theta instead.
	ex e(int n) const {return FrameBundle::e(n);}
	lst moduloIC;
	mutable int codim, maxC;
	ExVector f;	//components of the tautological form
};


//EDS for an SO(3)-structure on a 9-manifold
class SO3InSO9 : public EDSOnFrameBundle
{
	enum {dimension = 9};
public:
	string Geometry() {return "SO3 in SO9";}
	ex Gamma() const
	{
		ExVector e=Theta(); e.push_back(0);
		return ParseDifferentialForm(e,
"4679+2*3478-2789-3689+144 *3456+8*4569+144*5678-12*2468+2349-144*2457+72*2367-144*2358-144 *1467-128*1458-144*1368-4608*1357+144* 1278-8 *2589+8 *1249+144*1234-64512*1256-8 *1689"
		);

	}
	ex StarGamma() const
	{
		ExVector e=Theta(); e.push_back(0);
		return ParseDifferentialForm(e,
"2520 *35678+63 *34569+63 *56789-34789+2520*12378+24696*12569+35280*12358+63*12349-35280*12457+35280*15678-126*23679-63*23589-2520*13467-14*24689-35280*13456-63*24579-63*13689+63*12789-2520*23457-756*13579-63*14679+14*14589"
		);
	}

	lst Ideal() const {
		lst result;
		result.append(d(Gamma()));
		return result;
	}
	SO3InSO9() : EDSOnFrameBundle(dimension) {
	}
};



//EDS for an Einstein-Sasaki manifold of arbitrary dimension 2n+1
class EinsteinSasaki : public EDSOnFrameBundle
{
	int n;
public:
	string Geometry() const {return "Einstein-Sasaki";}
	
	EinsteinSasaki(int n) : EDSOnFrameBundle(2*n+1)
	{
		this->n=n;
	}

	lst Ideal() const {
		ex alpha=Theta(2*n+1);
		ex omega1;
		for (int i=1;i<=n;++i)
			omega1+=Theta(2*i-1)*Theta(2*i);
		ex Psi=1;
		ex omega2,omega3;
		for (int i=1;i<=n;++i)
			Psi=Psi*(Theta(2*i-1)+I*Theta(2*i));
		Psi=Psi.expand();
		SplitIntoRealAndImaginary(Psi,omega2,omega3);
	
		lst I;
		I.append(d(omega1));
		I.append(d(alpha*omega3));
		I.append(d(alpha*omega2));
		I.append(d(alpha));
//		I.append(d(alpha)+2*omega1);
//		I.append(d(omega2)-(n+1)*alpha*omega3);
//		I.append(d(omega3)+(n+1)*alpha*omega2);
		return I;
	}
};


//EDS for a parallel spinor in 7 dimensions i.e. a holonomy G2 manifold
class ParallelSeven : public EDSOnFrameBundle
{
public:
	ParallelSeven() : EDSOnFrameBundle(7) {}
	string Geometry() const {return "G_2";}
	lst Ideal() const {
		exvector phi=ParseDifferentialForms(Theta(),"567-512-534-613-642-714-723,1234-6712-6734-7513-7542-5614-5623");
		lst I;
		I=d(phi[0]), d(phi[1]);
		return I;
	}
};

class SU4 : public EDSOnFrameBundle
{
public:
	SU4() : EDSOnFrameBundle(15) {}
	string Geometry() const {return "SU_4";}
	lst Ideal() const {
		SU G(4);
		lst subs=LinearMapToSubstitutions<DifferentialForm>(G.e(),Theta());
		ex form=G.HodgeStar(G.ThreeForm()).subs(subs);
		lst I;
		I=d(form);
		return I;
	}
};

class SU3 : public EDSOnFrameBundle
{
public:
	SU3() : EDSOnFrameBundle(8) {}
	string Geometry() const {return "SU_3";}
	lst Ideal() const {
		SU G(3);
		lst subs=LinearMapToSubstitutions<DifferentialForm>(G.e(),Theta());
		ex form=G.HodgeStar(G.ThreeForm()).subs(subs);
		ex form2=(G.ThreeForm()).subs(subs);
		lst I;
		I=d(form), d(form2);
		return I;
	}
};

class Spn : public EDSOnFrameBundle
{
public:
	Spn(int n) : EDSOnFrameBundle(4*n) {}
	string Geometry() const {return "Sp_"+ToString(baseDimension/4);}
	lst Ideal() const {
		ex omega1,omega2,omega3;
		for (int i=1;i<baseDimension;i+=4)
		{
			omega1+=Theta(i)*Theta(i+1)+Theta(i+2)*Theta(i+3);
			omega2+=Theta(i)*Theta(i+2)+Theta(i+3)*Theta(i+1);
			omega3+=Theta(i)*Theta(i+3)+Theta(i+1)*Theta(i+2);
		}
		lst I;
		I=d(omega1), d(omega2),d(omega3);
		return I;
	}
};



//
// A more compact version of the same code follows
//
namespace Compact {

struct EDS :  ConcreteManifold, Has_dTable  {
	lst moduloIC;
	int dim;

	EDS(int dimension) : ConcreteManifold(dimension*(dimension+1)) {
		dim=dimension;
		for (int i=1;i<=dim;++i) {
			ex x;
			for (int j=1; j<=dim;++j)
				x+=e(j)*e(i*dim+j);
			Declare_d(e(i),x);
			moduloIC.append(e(i)==0);
		}
	}

	template<typename Container>
	void GetEquationsForVn(Container& container, lst I) {
		lst substitutions;
		for (int i=dim+1; i<=dim*(dim+1);++i) {
			ex ei;
			for (int j=1;j<=dim;++j)
				ei+=symbol("p"+ToString(i,j))*e(j);
			substitutions.append(e(i)==ei);
		}
	 	for (lst::const_iterator i=I.begin(); i!=I.end();++i)
			GetCoefficients(container,i->subs(substitutions));
	}

	template<typename Container>
	void GetReducedPolarEquations(Container& container, ex form, int j) {
		if (Degree<DifferentialForm>(form)==1)
			container.push_back(form.subs(moduloIC));
		else if (j>0) {
			GetReducedPolarEquations(container, form,j-1);
			GetReducedPolarEquations(container, Hook(e(j),form),j-1);
		}
	}
};


struct G2 : EDS {
	G2() : EDS (7) {
		lst I;
		I=d(ParseDifferentialForm(e(),"567-512-534-613-642-714-723")),
		  d(ParseDifferentialForm(e(),"1234-6712-6734-7513-7542-5614-5623"));
		Basis<symbol> A;
		GetEquationsForVn(A,I);

		for (int j=0;j<dim;++j)
		{				

			Basis<DifferentialForm> V;
			for (lst::const_iterator i=I.begin();i!=I.end();++i)
				GetReducedPolarEquations(V,*i,j);
			cout<<V.size()<<endl;
		}
		cout<<A.size()<<endl;
	}
};

}


int main()
{
	cout<<"Looking for an ordinary flag for the EDS of Sp(2)<SO(8)"<<endl;
	cout<<"(will not stop until one is found - probably never)"<<endl;
	while (true) {
		Spn X(2);
		X.ChooseFrame();
		X.CartansTest(1);
	}

	//Compact::G2 M;
	ParallelSeven B;
	B.FindOrdinaryFlag();

	for (int i=1;i<=3;++i)
	{		
		EinsteinSasaki M(i);
		cout<<"Looking for an ordinary flag for the Einstein-Sasaki structures in dimension "<<2*i+1<<endl;
		M.FindOrdinaryFlag();
	}


	cout<<"Looking for an ordinary flag for the EDS of SO(3)<SO(9)"<<endl;
	cout<<"(will not stop until one is found - probably never)"<<endl;
	SO3InSO9 X;	
	X.FindOrdinaryFlag(-1);

	SU4 A;
	A.AssumeCodimension(105);
	A.FindOrdinaryFlag();
	return 0;
}

