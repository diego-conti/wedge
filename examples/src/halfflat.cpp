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

//This program verifies some computations in [D. Conti, Half-flat nilmanifolds, arXiv:0903.1175]

#include <wedge/wedge.h>
#include "coherentsplittings.h"

using namespace Wedge;

class TestHalfFlat : public AbstractLieGroup<> {
public:	
	void Proposition3(const Frame& frame)
	{
		//check that we have a coherent splitting
		for (int i=1;i<=6;++i)
			if (!(d(e(i))*frame(1)*frame(2)).expand().is_zero()) return;
		if (!d(frame(1)).is_zero() || !d(frame(2)).is_zero()) return;
		//chech that b_6=1
		if (ExactForms(6).Dimension()>0) {
			cout<<"b_6=0"<<endl;
			return;
		}
		//compute the E_1
		cout<<"Computing h^{0,3} via Proposition 3"<<endl;
		cout<<"wrt frame: "<<frame<<endl;
		pqStructure P(this,frame,2);
		cout<<"\\dim E_1^{0,2}="<<P.E1(0,2)<<endl;
		cout<<"\\dim E_1^{1,2}="<<P.E1(1,2)<<endl;
		cout<<"\\dim E_1^{2,2}="<<P.E1(2,2)<<endl;

		VectorSpace<DifferentialForm> forms=P.pqForms(0,2);
		lst eqns;
		GetCoefficients<DifferentialForm>(eqns,P.delta1(forms.GenericElement()));
		cout<<"E_1^{0,2}="<<forms.SubspaceFromEquations(eqns.begin(),eqns.end()).e()<<endl;
		int n=forms.SubspaceFromEquations(eqns.begin(),eqns.end()).Dimension();
		assert (n==P.E1(0,2));
		vector<int> b=BettiNumbers();
		int m=b[3];
		assert (m==(1-b[1]+b[2])*2);
		assert (P.h(0,3)==m/2-n);
		cout << "h^{0,3}=="<<m/2<<" - "<<n<<" = "<< m/2-n<<endl;
	}

	//verify whether the algebra satisfies the conditions of Corollary 9.
	void Corollary9() {
		VectorSpace<DifferentialForm>  threeforms=ClosedForms(3),twoforms=ExactForms(2);		
		VectorSpace<DifferentialForm> t=pForms(2);
		bool ok;
		{
			list<ex> eqns;
			for (exvector::const_iterator j=threeforms.e_begin();j!=threeforms.e_end();++j)
				GetCoefficients<DifferentialForm>(eqns,t.GenericElement() **j);
			Subspace<DifferentialForm> S=t.SubspaceFromEquations(eqns.begin(),eqns.end()).e();
			ok=(S.Dimension()==0);
			if (ok)
				cout<<"Map \\Lambda^2\\to \\Hom(Z_3,\\Lambda^5) is injective"<<endl; 
			else
				cout<<"Map \\Lambda^2\\to \\Hom(Z_3,\\Lambda^5) has kernel "<<S.e()<<endl;
		}
		{
			list<ex> eqns;
			for (exvector::const_iterator j=twoforms.e_begin();j!=twoforms.e_end();++j)
				GetCoefficients<DifferentialForm>(eqns,t.GenericElement() **j);
			Subspace<DifferentialForm> S=t.SubspaceFromEquations(eqns.begin(),eqns.end()).e();
			if (S.Dimension()==0)
			{	
				ok=false; 
				cout<<"Map \\Lambda^2\\to \\Hom(B_2,\\Lambda^5) is injective"<<endl;
			}
			else 				
				cout<<"Map \\Lambda^2\\to \\Hom(B_2,\\Lambda^5) has kernel "<<S.e()<<endl;
		}
		if (ok) cout<<"Corollary 9 is satisfied"<<endl;
		else cout<<"Corollary 9 not satisfied"<<endl;
	}

	void ComputeCohomology()
	{

		int K=ClosedForms(1).Dimension();
		K=2;
		ex exacttwoform=ExactForms(2).GenericElement();
		for (int i=1;i<=K;++i)
		{
			assert(d(e(i))==0);	//we are assuming that ker d stays on front, i.e. lie{g}=(0,...,0, nonzero terms ).
			exacttwoform*=e(i);
		}
		if (!exacttwoform.expand().is_zero()) {
			cout<<"This algebra has reduced length greater than two, hence no coherent splitting"<<endl;
			return;
		}
		cout<<"Betti numbers relative to ""natural"" splitting ker d\\oplus \\ker d^\\perp and ""dual"" splitting ker d^\\perp\\oplus \\ker d:"<<endl;
		exvector dual(e().begin(),e().begin()+K);
		dual.insert(dual.begin(),e().begin()+K,e().end());
		assert(dual.size()==e().size());
		pqStructure P(this,e(),K), P2(this,dual,Dimension()-K,-1,0);

		for (int k=2;k<=5;++k)
		{
			int r=0;
			for (int p=0;p<=k;++p)
			{
				int h=P.h(p,k-p);
				r+=h;
				cout<<h;
				if (p<k) cout<<" + ";
			}
			if (BettiNumbers()[k]!=r) {
				LOG_ERROR("\\sum_{p+q=k} h^p,q}\\neq b_k");
				cout<<" != "<<BettiNumbers()[k]<<" !!!!!!!!!"<<endl;
			}
			else cout<<" = b_"<<k<<" = "<<r<<endl;
			r=0;
			for (int p=0;p<=k;++p)
			{
				int h=P2.h(p,k-p);
				r+=h;
				cout<<h;
				if (p<k) cout<<" + ";
			}
			if (BettiNumbers()[k]!=r) {
				LOG_ERROR("\\sum_{p+q=k} h^p,q}\\neq b_k");
				cout<<" != "<<BettiNumbers()[k]<<" !!!!!!!!!"<<endl;
			}
			else cout<<" = b_"<<k<<" = "<<r<<endl;
		}

	}

	void CoherentSplittings() {
		Frame kerd=exvector(ClosedForms(1).e());		

		VectorSpace<DifferentialForm>  twoforms=Wedge::pForms(kerd,2);
		list<ex> eqns;
		for (int i=1;i<=6;++i)
			GetCoefficients<DifferentialForm> (eqns,twoforms.GenericElement()*d(e(i)));
		Subspace<DifferentialForm> s=twoforms.SubspaceFromEquations(eqns.begin(),eqns.end());
		eqns.clear();
		
		if (s.Dimension()==0) {
			cout<<"No coherent splittings"<<endl;
			return;
		}
		cout<<"A two-form defines a coherent splitting iff it is simple and lies in  "<<s.e()<<endl;
/*
		for (int i=1;i<=s.Dimension();++i)
		{
			VectorSpace<DifferentialForm> k;
			for (exvector::const_iterator j=threeforms.e_begin();j!=threeforms.e_end();++j)
				k.AddGenerator(*j * s.e(i));		
			cout<<"V^1="<<s.e(i)<<" H03="<<k.Dimension();
			if ((ClosedForms(4).GenericElement() * s.e(i)).expand().is_zero())
				cout<<" H04=0"<<endl;
			else
				cout<<" H04=1"<<endl;
		}
*/
		VectorSpace<DifferentialForm>  threeforms=ClosedForms(3);
		for (exvector::const_iterator i=threeforms.e_begin();i!=threeforms.e_end();++i)
			GetCoefficients<DifferentialForm> (eqns,s.GenericElement()* *i);
		Subspace<DifferentialForm> h= s.SubspaceFromEquations(eqns.begin(),eqns.end());
		if (h.Dimension()!=0) cout<<"H_{03}=0 for V_1 spanned by a form in "<<h.e()<<endl;
		else  cout<<"H_{03}!=0 for all coherent splittings"<<endl;
	}

	TestHalfFlat (const char * str) : AbstractLieGroup<>(str) {
		cout<<"Testing algebra "<<str<<endl;
		CoherentSplittings();
		ComputeCohomology();
		Corollary9();
		Proposition3(e());
		Proposition3(ParseDifferentialForms(e(),"1,3,2,4,5,6"));
		cout<<endl;
	}
};


class SeekHalfFlat : public AbstractLieGroup<>, public IterateOverPermutations
{
public:
	int add;
	bool Apply (const vector<int>& permutation)
	{
		ExVector f;
		for (vector<int>::const_iterator i=permutation.begin();i!=permutation.end();++i)
			f.push_back(e()[*i]);
		for (int i=1;i<=add;++i)
			if (rand()%2==0)
				f(i)+=f((i+rand()%5)%6+1);
			else
				f(i)-=f((i+rand()%5)%6+1);
		ex omega=(f(1)*f(2)*f(3)*f(4)+f(1)*f(2)*f(5)*f(6)+f(3)*f(4)*f(5)*f(6)).expand();
		ex psiplus=(f(1)*f(3)*f(5)-f(2)*f(4)*f(5)-f(2)*f(3)*f(6)-f(1)*f(4)*f(6)).expand();
		if (d(omega).expand().is_zero() && d(psiplus).expand().is_zero())
		{
			cout <<" & "<<f<<"\\\\"<<endl;
			return false;
		}
		return true;
	}


//search for a "simple" half-flat structure
	SeekHalfFlat(const char* str) : AbstractLieGroup<>(str) 
	{		
		cout<<BettiNumbers()[2]<<"&";
		cout<<str;
		for (add=0;add<3; ++add) 
			if (!Iterate(6)) return;
		cout<<" no half-flat structure found\\\\"<<endl;
	}

//verify the given frame determines a half-flat structure
	SeekHalfFlat(const char* str, const char* frame) : AbstractLieGroup<>(str) 
	{
		if (Dimension()!=6) {
			LOG_ERROR(str);
			throw WedgeException<std::runtime_error>("Wrong dimension",__FILE__,__LINE__);
		}
		ExVector f=ParseDifferentialForms(e(),frame);
		if (f.size()!=6) {
			LOG_ERROR(frame);
			throw WedgeException<std::runtime_error>("Wrong size",__FILE__,__LINE__);
		}
		Frame basis(f);
		if (basis.size()!=6) {
			LOG_ERROR(frame);
			throw WedgeException<std::runtime_error>("Not a basis",__FILE__,__LINE__);
		} 
		ex omega=ParseDifferentialForm(f,"12+34+56");
		ex psiplus=ParseDifferentialForm(f,"135-146-245-236");
		cout<<BettiNumbers()[2]<<"&";
		cout<<str;
		if (d(omega*omega).expand().is_zero() && d(psiplus).expand().is_zero()) 
			cout <<" & "<<f<<"\\\\"<<endl;
		else
			cout<<" not half flat!\\\\"<<endl;

	}
};

class pqManifold : public AbstractLieGroup<> {
	int K;
public:
	pqManifold(const char* str) : AbstractLieGroup<>(str) {
		cout<<"Lie algebra:" <<StructureConstants()<<endl;
/*
		VectorSpace<DifferentialForm>  twoforms=pForms(2);
		VectorSpace<DifferentialForm>  closed3forms=ClosedForms(3);
		list<ex> eqns;
		for (exvector::const_iterator i=closed3forms.e_begin();i!=closed3forms.e_end();++i)
			GetCoefficients<DifferentialForm>(eqns,*i*twoforms.GenericElement());
		Subspace<DifferentialForm> s2=twoforms.SubspaceFromEquations(eqns.begin(),eqns.end());
		cout<<s2.e()<<endl;
		return; 
*/
		CS2();
		CoherentSplittings();
		ComputeCohomology();
	}

	void CS2() const {
		CoherentSplitting V1(this,2);
		set<ex,ex_is_less> eqns;
		//cout<<"Generic splitting: "<<V1.DefiningForm()<<endl;
		cout<<"Coherent equations: "<<endl<<V1.GetCoherentEquations(eqns)<<endl;
		eqns.clear();
		try {
			V1.DeclareH0qZero(2);
			cout<<"Coherent equations for H^{02}=0: "<<endl<<V1.GetCoherentEquations(eqns)<<endl;		
			eqns.clear();
			V1.DeclareH0qZero(3);
			cout<<"Coherent equations for H^{02}=0=H^{03}: "<<endl<<V1.GetCoherentEquations(eqns)<<endl;
		}
		catch (const InconsistentDeclaration&)
		{
			cout<<"No coherent splitting with H^{02}=0=H^{03}"<<endl;
		}
		cout<<V1.DefiningForm()<<endl<<endl;
	}
	//iterate through a list of coherent splittings (dim V_1=2) and compute cohomology for each
	void CoherentSplittings() {
//nilpotent case		Frame kerd=exvector(ClosedForms(1).e());		
		Frame kerd=e();	//TODO 
		VectorSpace<DifferentialForm>  twoforms=Wedge::pForms(kerd,2);
		list<ex> eqns;
		for (int i=1;i<=Dimension();++i)
			GetCoefficients<DifferentialForm> (eqns,twoforms.GenericElement()*d(e(i)));
		GetCoefficients<DifferentialForm>(eqns,d(twoforms.GenericElement()));
		Subspace<DifferentialForm> s=twoforms.SubspaceFromEquations(eqns.begin(),eqns.end());
		
		if (s.Dimension()==0) {
			cout<<"No coherent splittings"<<endl;
			return;
		}
		cout<<"A two-form defines a coherent splitting only if it is simple and lies in  "<<s.e()<<endl;
		for (int k=2;k<=Dimension()-2;++k) {
			eqns.clear();
			VectorSpace<DifferentialForm>  closed3forms=ClosedForms(k);
			for (exvector::const_iterator i=closed3forms.e_begin();i!=closed3forms.e_end();++i)
				GetCoefficients<DifferentialForm>(eqns,*i*s.GenericElement());
			LOG_INFO(eqns);
			Subspace<DifferentialForm> s2=s.SubspaceFromEquations(eqns.begin(),eqns.end());
			if (s2.Dimension()==0) cout<<" For all splittings h^{0,"<<k<<"} is non-zero"<<endl;
			else {
				cout<<"A two-form defines a coherent splitting with h^{0,"<<k<<"}=0 iff it is simple and lies in  "<<s2.e()<<endl;
				//cout<<closed3forms<<endl;
			}
			
			eqns.clear();
			for (exvector::const_iterator i=closed3forms.e_begin();i!=closed3forms.e_end();++i)
				GetCoefficients<DifferentialForm>(eqns,Hook(s.GenericElement(),*i));
			LOG_INFO(eqns);
			Subspace<DifferentialForm> s3=s.SubspaceFromEquations(eqns.begin(),eqns.end());
			if (s3.Dimension()==0) cout<<" For all splittings \\tilde h^{2,"<<k-2<<"} is non-zero"<<endl;
			else {
				cout<<"A two-form defines a coherent splitting with \\tilde h^{2,"<<k-2<<"}=0 iff it is simple and lies in  "<<s3.e()<<endl;
				//cout<<closed3forms<<endl;
			}

		} 
	return;
		//qui potrei fare qualcosa dipiu elaboarto come fissare e(i) e calcolare e(j)
		for (int i=1;i<=kerd.size();++i)
		for (int j=i+1;j<=kerd.size();++j)
			if (s.Contains(kerd(i)*kerd(j)))
			{
				cout<<"Cohomology for V^1 spanned by "<<kerd(i)<<","<<kerd(j)<<":"<<endl;
				SubBasis<DifferentialForm> f;
				f.push_back(kerd(i));
				f.push_back(kerd(j));
				f.AddToComplement(e().begin(),e().end());
				
				exvector dual(f.begin(),f.end());
				dual.insert(dual.begin(),f.complement_begin(),f.complement_end());
				assert(dual.size()==e().size());
				pqStructure P(this,exvector(f.full_begin(),f.full_end()),2), P2(this,dual,Dimension()-2,-1,0);
	
				for (int k=2;k<=Dimension()-1;++k)
				for (int p=0;p<=2;++p)
				{
					cout<<"h^{"<<p<<","<<k-p<<"}="<<P.h(p,k-p)<<endl;
					cout<<"\\tilde{h}^{"<<p<<","<<k-p<<"}="<<P.h(k-p,p)<<endl;
				}
			}			

	}

	void ComputeCohomology()
	{
		int K=ClosedForms(1).Dimension();
		ex exacttwoform=ExactForms(2).GenericElement();
		for (int i=1;i<=K;++i)
		{
			assert(d(e(i))==0);	//we are assuming that ker d stays on front, i.e. lie{g}=(0,...,0, nonzero terms ).
			exacttwoform*=e(i);
		}
		if (!exacttwoform.expand().is_zero()) {
			cout<<"This algebra has reduced length greater than two, hence no coherent splitting"<<endl;
			return;
		}
		cout<<"Betti numbers relative to ""natural"" splitting ker d\\oplus \\ker d^\\perp and ""dual"" splitting ker d^\\perp\\oplus \\ker d:"<<endl;
		exvector dual(e().begin(),e().begin()+K);
		dual.insert(dual.begin(),e().begin()+K,e().end());
		assert(dual.size()==e().size());
		pqStructure P(this,e(),K), P2(this,dual,Dimension()-K,-1,0);
		vector<int> bettiNumbers=BettiNumbers();

		for (int k=2;k<=3;++k)
		{
			int r=0;
			for (int p=0;p<=k;++p)
			{
				int h=P.h(p,k-p);
				r+=h;
				cout<<h;
				if (p<k) cout<<" + ";
			}
			if (bettiNumbers[k]!=r) {
				LOG_ERROR("\\sum_{p+q=k} h^p,q}\\neq b_k");
				cout<<" != "<<bettiNumbers[k]<<" !!!!!!!!!"<<endl;
			}
			else cout<<" = b_"<<k<<" = "<<r<<endl;
			r=0;
			for (int p=0;p<=k;++p)
			{
				int h=P2.h(p,k-p);
				r+=h;
				cout<<h;
				if (p<k) cout<<" + ";
			}
			if (bettiNumbers[k]!=r) {
				LOG_ERROR("\\sum_{p+q=k} h^p,q}\\neq b_k");
				cout<<" != "<<bettiNumbers[k]<<" !!!!!!!!!"<<endl;
			}
			else cout<<" = b_"<<k<<" = "<<r<<endl;
		}

	}
};


int main()
{	

	{

	AbstractLieGroup<> G("0, 0, 12, 0, 54, 0");
	VectorSpace<DifferentialForm> M=G.ClosedForms(3);
	VectorSpace<DifferentialForm> N=G.ClosedForms(4);
	ex m=M.GenericElement();
	ex n=N.GenericElement();
	list<ex> eqns;
	for (int j=1;j<=G.Dimension();++j)
	for (int i=1;i<=G.Dimension();++i)
		if ((Hook(G.e(i),m)*m*G.e(j)).expand().is_zero()) cout<<i<<j<<endl;
	cout<<(N.GenericElement()*G.e(1)*G.e(4)).expand()<<endl;
	cout<<(N.GenericElement()*G.e(2)*G.e(4)).expand()<<endl;
	}

	pqManifold("0, 0, 12, 13, 23, 14");
	pqManifold("0, 0, 12, 13, 23, 14 + 25");
	pqManifold("0, 0, 12, 13,23, 14 - 25");
	pqManifold("0,0,12,13,14+23,24+15");
	pqManifold("0, 0, 0, 12, 14, 15 + 23");
	pqManifold("0, 0, 0, 12, 14 - 23, 15 + 34");//no obstru
	pqManifold("0, 0, 0, 12, 14, 15");
	pqManifold("0, 0, 0, 12, 23, 14 + 35");	//no osbtru
	pqManifold("0, 0, 0, 12, 23, 14 - 35");//
	pqManifold("0, 0, 0, 12, 13, 14 + 35");
	pqManifold("0, 0, 0, 12, 13, 14 + 23");
	pqManifold("0, 0, 0, 12, 13, 24");    
	pqManifold("0, 0, 0, 12, 13, 23");    //
	pqManifold("0, 0, 0, 12, 14, 15 + 24");
	pqManifold("0, 0, 0, 12, 14, 15+ 23+ 24");
	pqManifold("0, 0, 0, 0, 12, 14 + 25");
	pqManifold("0, 0, 0, 0, 12, 15 + 34");//
	pqManifold("0, 0, 0, 0, 13 + 42, 14 + 23");//
	pqManifold("0, 0, 0, 0, 12, 14 + 23");//
	pqManifold("0, 0, 0, 0, 12, 13");//
	pqManifold("0, 0, 0, 0, 12, 34");   //
	pqManifold("0, 0, 0, 0, 0, 12 + 34");//
	pqManifold("0, 0, 0, 0, 0, 12");
	return 0;

/*
	pqManifold("0,0,12,13,14+23");

*/

/*
	SeekHalfFlat("0,12,-13,0,45,-46");
	pqManifold("0,12,-13,0,45,-46");
	return 0;
*/


	cout<<"The following algebras have a half-flat structure:"<<endl;
	cout<<"\\begin{array}{|l|l|l|}"<<endl;
	SeekHalfFlat("0, 0, 12, 13, 23, 14","1,5,2,4,3,6");
	SeekHalfFlat("0, 0, 12, 13, 23, 14 + 25","1-2,4,5,2,6,3");
	SeekHalfFlat("0, 0, 12, 13,23, 14 - 25","3,6,4,4-2,-5,1+5");
	SeekHalfFlat("0,0,12,13,14+23,24+15","-5,2,4,1,sqrt(2)*(3-5),sqrt(1/2)*6");
	SeekHalfFlat("0, 0, 0, 12, 14, 15 + 23","2+5,2+5+6,4,2,3,1");
	SeekHalfFlat("0, 0, 0, 12, 14 - 23, 15 + 34","2,4,3,1,6,5");
	SeekHalfFlat("0, 0, 0, 12, 14, 15","1,3,2,5,4,6");
	SeekHalfFlat("0, 0, 0, 12, 23, 14 + 35","1,3,4,5,6,2");
	SeekHalfFlat("0, 0, 0, 12, 23, 14 - 35","1,3,2,6,5,4");
	SeekHalfFlat("0, 0, 0, 12, 13, 14 + 35","2,6,-3,4,1,-2-5");
	SeekHalfFlat("0, 0, 0, 12, 13, 14 + 23","2-6,1+5,4,1,6,3");
	SeekHalfFlat("0, 0, 0, 12, 13, 24","1,6,2,3,4,5");    
	SeekHalfFlat("0, 0, 0, 12, 13, 23","1,4,2,5,3,6");    
	SeekHalfFlat("0, 0, 0, 12, 14, 15 + 24","1,3,2,4,3+5,-6");
	SeekHalfFlat("0, 0, 0, 12, 14, 15+ 23+ 24","1,3,2,4-2,3+5,-6");
	SeekHalfFlat("0, 0, 0, 0, 12, 14 + 25","1-6,4,5,2,6,3");
	SeekHalfFlat("0, 0, 0, 0, 12, 15 + 34","1,3,5,4,6,2");
	SeekHalfFlat("0, 0, 0, 0, 13 + 42, 14 + 23","1,2,3,4,5,6");
	SeekHalfFlat("0, 0, 0, 0, 12, 14 + 23","1,3,2,4,6,5");
	SeekHalfFlat("0, 0, 0, 0, 12, 13","1,4,2,3,5,6");
	SeekHalfFlat("0, 0, 0, 0, 12, 34","1+3,1,6,5,2,4");   
	SeekHalfFlat("0, 0, 0, 0, 0, 12 + 34","1,2,4,3,5,6");
	SeekHalfFlat("0, 0, 0, 0, 0, 12","1,3,2,4,5,6");
	cout<<"\\end{array}"<<endl;

	cout<<"The following algebras have a half-flat structure:"<<endl;
	TestHalfFlat("0, 0, 12, 13, 23, 14");
	TestHalfFlat("0, 0, 12, 13, 23, 14 + 25");
	TestHalfFlat("0, 0, 12, 13,23, 14 - 25");
	TestHalfFlat("0,0,12,13,14+23,24+15");
	TestHalfFlat("0, 0, 0, 12, 14, 15 + 23");
	TestHalfFlat("0, 0, 0, 12, 14 - 23, 15 + 34");
	TestHalfFlat("0, 0, 0, 12, 14, 15");
	TestHalfFlat("0, 0, 0, 12, 23, 14 + 35");
	TestHalfFlat("0, 0, 0, 12, 23, 14 - 35");
	TestHalfFlat("0, 0, 0, 12, 13, 14 + 35");
	TestHalfFlat("0, 0, 0, 12, 13, 14 + 23");
	TestHalfFlat("0, 0, 0, 12, 13, 24");    
	TestHalfFlat("0, 0, 0, 12, 13, 23");    
	TestHalfFlat("0, 0, 0, 12, 14, 15 + 24");
	TestHalfFlat("0, 0, 0, 12, 14, 15+ 23+ 24");
	TestHalfFlat("0, 0, 0, 0, 12, 14 + 25");
	TestHalfFlat("0, 0, 0, 0, 12, 15 + 34");
	TestHalfFlat("0, 0, 0, 0, 13 + 42, 14 + 23");
	TestHalfFlat("0, 0, 0, 0, 12, 14 + 23");
	TestHalfFlat("0, 0, 0, 0, 12, 13");
	TestHalfFlat("0, 0, 0, 0, 12, 34");   
	TestHalfFlat("0, 0, 0, 0, 0, 12 + 34");
	TestHalfFlat("0, 0, 0, 0, 0, 12");

	cout<<"The following algebras have no half-flat structure and no coherent splitting:"<<endl;
	TestHalfFlat("0,0,12,13,14+23,34+52");
	TestHalfFlat("0, 0, 12, 13, 14, 34 + 52");

	cout<<"The following algebras have a coherent splitting but no half-flat structure:"<<endl;
	TestHalfFlat("0, 0, 12, 13, 14, 15") ; 
	TestHalfFlat("0,0,12,13,14,23+15");
	TestHalfFlat("0,0,0,12,14,24");
	TestHalfFlat("0,0,0,12,13+42,14+23");
	TestHalfFlat("0,0,0,12,14,13+42");
	TestHalfFlat("0,0,0,12,13+14,24");
	TestHalfFlat("0,0,0,12,13,14");
	TestHalfFlat("0,0,0,0,12,15");

	cout<<"Now considering solvable Lie group (0,12,13,14,15,16)"<<endl;
	AbstractLieGroup<> G("0,12,13,14,15,16");
	pqStructure P(&G,G.e(),2);
	cout<<"h^{03}="<<P.h(0,3)<<endl;
	cout<<"h^{04}="<<P.h(0,4)<<endl;
	return 0;
}

