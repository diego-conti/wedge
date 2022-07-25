/***************************************************************************
 *   Copyright (C)2009 by Diego Conti					   *
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
#include <wedge/wedge.h>
#include <boost/date_time/posix_time/posix_time_types.hpp>
#include <climits>


using namespace Wedge;
using namespace GiNaC;
using namespace std;
using namespace boost::posix_time;

//#define GENERATE_THE_CODE

#ifdef GENERATE_THE_CODE
#define MEASURE(CLASS_NAME, PARAM, REFERENCETIME) \
cout<< "\tMEASURE("<<#CLASS_NAME<<", "<<PARAM<<", "<<Measure<CLASS_NAME>(PARAM)<<");"<<endl;
#else
#define MEASURE(CLASS_NAME, PARAM, REFERENCETIME) \
cout<< (Measure<CLASS_NAME >(PARAM, REFERENCETIME)).perf() <<"% - "<< #CLASS_NAME<<", parameter "<<PARAM<<endl;
#endif

struct Performance { 
#if LONG_MAX > 2147483647
	typedef long TimeSpan;
#else
	typedef double TimeSpan;
#endif
	TimeSpan elapsed;
	TimeSpan reference;
	Performance(TimeSpan s, TimeSpan ref) {elapsed=s; reference=ref;}
	int perf() const {LOG_INFO(elapsed); LOG_INFO(reference); LOG_INFO((reference*100)/elapsed); return static_cast<int>((reference*100)/elapsed);}
};

ostream& operator<< (ostream& os, Performance perf)
{
	return os<<perf.elapsed;
}

vector<Performance> performance;

template<typename T, typename V> Performance Measure(V v, Performance::TimeSpan referencetime=1)
{
	srand(0);	//we do not really want "random" elements.
	T t(v);
	int iterations=0;
	ptime cur_clock;
	ptime start_clock=microsec_clock::local_time();
	ptime end_clock=start_clock+seconds(2);
	do {	
		t.run();
		cur_clock=microsec_clock::local_time();
		++iterations;
	} while (cur_clock<end_clock);
	LOG_INFO((cur_clock-start_clock).total_nanoseconds());
	Performance perf((cur_clock-start_clock).total_nanoseconds()/iterations, referencetime);
	performance.push_back(perf);
	return perf;
}


//global variables

Frame bigGlobalFrame;
Frame nonSimpleGlobalFrame;
VectorSpace<DifferentialOneForm> nonSimpleGlobalVSpace;
SU globalLieGroup (10);
AbstractLieGroup<> globalRiemannianLieGroup ("0,0,12,13,14+23,24+15");
LeviCivitaConnection<true> globalRiemannianConnection (&globalRiemannianLieGroup,RiemannianStructure(&globalRiemannianLieGroup,globalRiemannianLieGroup.e()));

//Add some random elements to the container selected from the basis e.
template<typename Container> void AddRandomElements(Container& container, int space_size, int basis_size, int lengths)
{
	for (int i=0;i<basis_size;++i)
	{
		ex new_elem;
		for (int j=1; j<lengths;++j)
			new_elem+=(rand()%20)*bigGlobalFrame[rand()%space_size];
		container.push_back(new_elem);
	}
}


// the tests
struct EmptyTest {
	EmptyTest(int) {}
	void run() {}
};
struct DenseBasisTest {
	int N;
	DenseBasisTest(int n) {N=n;}
	void run() {
		Frame e;
		AddRandomElements(e,N,N/2+1,N/3+1);
	}
};
struct DenseBasisTest2 {
	exvector e;
	DenseBasisTest2(int N) {
		AddRandomElements(e,N,N/2+1,N/3+1);
	}
	void run() {
		Frame f;
		f.insert(f.begin(),e.begin(),e.end());
	}
};

struct SparseBasisTest {
	int N;
	SparseBasisTest(int n) {N=n;}
	void run() {
		Frame e;
		AddRandomElements(e,N,N/4+1,N/3+1);
	}
};
struct SparseBasisTest2 {
	exvector e;
	SparseBasisTest2(int N) {
		AddRandomElements(e,N,N/4+1,N/3+1);
	}
	void run() {
		Frame f;
		f.insert(f.begin(),e.begin(),e.end());
	}
};

struct DualTest {
	Frame e;
	DualTest(int N) : e(nonSimpleGlobalFrame.begin(),nonSimpleGlobalFrame.begin()+N) {}
	void run() {
		Frame (e).dual();
	}
};

struct ContainsTest {
	exvector e;
	ContainsTest(int N) {
		AddRandomElements(e,200,1,N);
	}
	void run() {
		nonSimpleGlobalVSpace.Contains(e[0]);
	}
};

struct GetSolutionsTest {
	VectorSpace<DifferentialOneForm> space;
	list<ex> eqns;
	GetSolutionsTest(int N) {
		Frame f(nonSimpleGlobalFrame.begin(),nonSimpleGlobalFrame.begin()+N);
		space.SetBasis(f);
		for (int i=0;i<N/3;++i)
		{
			ex eqn;
			for (int k=1;k<=5;++k)
				eqn+=(rand()%20-10)*space.coordinate(rand()%N+1);
			eqns.push_back(eqn);
		}
	}
	void run() {
		list<ex> sols;
		space.GetSolutions(sols,eqns.begin(),eqns.end());
	}
};

struct SubspaceFromEquationsTest {
	VectorSpace<DifferentialOneForm> space;
	list<ex> eqns;
	SubspaceFromEquationsTest(int N) {
		Frame f(nonSimpleGlobalFrame.begin(),nonSimpleGlobalFrame.begin()+N);
		space.SetBasis(f);
		for (int i=0;i<N/3;++i)
		{
			ex eqn;
			for (int k=1;k<=5;++k)
				eqn+=(rand()%20-10)*space.coordinate(rand()%N+1);
			eqns.push_back(eqn);
		}
	}
	void run() {
		space.SubspaceFromEquations(eqns.begin(),eqns.end());
	}
};

struct SubspaceFromEquationsOrthTest {
	VectorSpace<DifferentialOneForm> space;
	list<ex> eqns;
	SubspaceFromEquationsOrthTest(int N) {
		Frame f(nonSimpleGlobalFrame.begin(),nonSimpleGlobalFrame.begin()+N);
		space.SetBasis(f);
		for (int i=0;i<N/3;++i)
		{
			ex eqn;
			for (int k=1;k<=5;++k)
				eqn+=(rand()%20-10)*space.coordinate(rand()%N+1);
			eqns.push_back(eqn);
		}
	}
	void run() {
		space.SubspaceFromEquations(eqns.begin(),eqns.end(),TrivialPairingOperator<DifferentialOneForm>());
	}
};

struct SubspaceTest {
	VectorSpace<DifferentialOneForm> space;
	exvector vectors;
	SubspaceTest(int N) {
		Frame f(nonSimpleGlobalFrame.begin(),nonSimpleGlobalFrame.begin()+N);
		space.SetBasis(f);
		for (int i=0;i<N/3;++i)
		{
			ex v;
			for (int k=1;k<=5;++k)
				v+=(rand()%20-10)*space.e(rand()%N+1);
			vectors.push_back(v);
		}
	}
	void run() {
		space.Subspace(vectors.begin(),vectors.end());
	}
};

struct SubspaceOrthTest {
	VectorSpace<DifferentialOneForm> space;
	exvector vectors;
	SubspaceOrthTest(int N) {
		Frame f(nonSimpleGlobalFrame.begin(),nonSimpleGlobalFrame.begin()+N);
		space.SetBasis(f);
		for (int i=0;i<N/3;++i)
		{
			ex v;
			for (int k=1;k<=5;++k)
				v+=(rand()%20-10)*space.e(rand()%N+1);
			vectors.push_back(v);
		}
	}
	void run() {
		space.Subspace(vectors.begin(),vectors.end(),TrivialPairingOperator<DifferentialOneForm>());
	}
};

struct FormTest {
	ConcreteManifold M;
	int N;
	FormTest(int n) :  M(exvector(bigGlobalFrame.begin(),bigGlobalFrame.begin()+n)) {N=n;}
	void run() {
		ex sum;
		for (int k=0;k<=30;++k)
		{
			ex form=1;
			for (int i=0;i<N/2;++i)
				form*=(rand()%10-5)*M.e()[rand()%N];
			sum+=form;
		}
	}
};

struct FormTest2 {
	ConcreteManifold M;
	int N;
	FormTest2(int n) :  M(exvector(bigGlobalFrame.begin(),bigGlobalFrame.begin()+n)) {N=n;}
	void run() {
		ex product=1;
		for (int k=0;k<N/4;++k)
		{
			ex form=0;
			for (int i=0;i<N/2;++i)
				form+=(rand()%10-5)*M.e()[rand()%N]*M.e()[rand()%N];
			product*=form;
		}
	}
};

struct HookTest {
	ex left, right;
	HookTest(int N) {
		ConcreteManifold M(exvector(bigGlobalFrame.begin(),bigGlobalFrame.begin()+N));
		for (int i=0;i<N/2;++i)
			left+=(rand()%10-5)*M.e()[rand()%N]*M.e()[rand()%N];
		for (int i=0;i<N/2;++i)
			right+=(rand()%10-5)*M.e()[rand()%N]*M.e()[rand()%N]*M.e()[rand()%N]*M.e()[rand()%N]*M.e()[rand()%N];
	}
	void run() {
		Hook(left,right);
	}
};

struct dTest {
	ex form;
	dTest(int N)
	{
		for (int i=0;i<N/2;++i)
			form+=(rand()%10-5)*globalLieGroup.e()[rand()%N]*globalLieGroup.e()[rand()%N];
	}
	void run() {
		globalLieGroup.d(form);
	}
};

struct RiemannianTest {
	ex left,right;
	RiemannianStructure P;
	RiemannianTest(int N) : P(&globalLieGroup,globalLieGroup.e())
	{
		for (int i=0;i<N/2;++i)
			left+=(rand()%10-5)*P.e()[rand()%N]*P.e()[rand()%N];
		for (int i=0;i<N/2;++i)
			right+=(rand()%10-5)*P.e()[rand()%N]*P.e()[rand()%N];
	}
	void run() {
		P.ScalarProduct<DifferentialForm>(left,right);
	}
};

struct ConnectionTest {
	ConcreteManifold M;
	RiemannianStructure P;
	ConnectionTest(int N) : M(exvector(bigGlobalFrame.begin(),bigGlobalFrame.begin()+N)),  P(&M,M.e()) {}
	void run() {
		RiemannianConnection L(&M,P);
	}
};

struct ConnectionTest2 {
	ex left,right;
	ConnectionTest2(int N)
	{
		ex left;
		for (int i=0;i<N/2;++i)
			left+=(rand()%10-5)*globalRiemannianLieGroup.e()[rand()%globalRiemannianLieGroup.Dimension()];

		ex right;
		for (int i=0;i<N/2;++i)
			right+=(rand()%10-5)*globalRiemannianLieGroup.e()[rand()%globalRiemannianLieGroup.Dimension()]*globalRiemannianLieGroup.e()[rand()%globalRiemannianLieGroup.Dimension()];
	}
	void run() {
		for (int i=1;i<=10;++i)
			globalRiemannianConnection.Nabla<DifferentialForm>(left,right);
	}
};


struct LeviCivitaConnectionTest : public ConcreteManifold, public Has_dTable {
	RiemannianStructure g;
	LeviCivitaConnectionTest(int N) : ConcreteManifold(N) , g(this,e()) {
		for (int i=1;i<=N;++i) 
			Declare_d(e(i),e()[rand()%N]*e()[rand()%N]);
		ExVector f=e();
		for (int i=1;i<N;++i) {
			int k=i+1+rand()%(N-i);
			f(i)+=(rand()%10)*e(k);
		}
		g=RiemannianStructure(this,f);
	}	
	void run() {
		LeviCivitaConnection<true>(this,g);
	}
};

struct CurvatureTest : public ConcreteManifold, public Has_dTable {
	RiemannianStructure g;
	LeviCivitaConnection<true> omega;

	ExVector CreateFrame(int N) {
		for (int i=1;i<=N;++i)
			Declare_d(e(i),e()[rand()%N]*e()[rand()%N]);
		ExVector f=e();
		for (int i=1;i<N;++i) {
			int k=i+1+rand()%(N-i);
			f(i)+=(rand()%10)*e(k);
		}
		return f;
	}
	CurvatureTest(int N) : ConcreteManifold(N) , g(this,CreateFrame(N)), omega(this,g) {}
	void run() {
		matrix ricci=omega.RicciAsMatrix();
	}
};

struct ncmulTest {
	ex a;
	ncmulTest(int) {
		DifferentialOneForm u,v,w,s,t,r;
		symbol x,y,p;
		a=x*u*v+y*w*s+p*t*r;
	}
	void run() {
		for (int i=0;i<100;++i)
			a*a*a;      
	}
};

template<typename matrixType> struct matrixTest {
	matrixType m;
	matrixTest(int N) : m(4,4) {
		symbol x1("x"),x2("y"),x3("z"),x4("t");
		m=matrix{{	2/(1+pow(x2,2)+pow(x1,2)),0,0,0},
			{0,2/(1+pow(x2,2)+pow(x1,2)),0,0},
			{-2/(1+pow(x2,2)+pow(x1,2))*x4*x2,2/(1+pow(x2,2)+pow(x1,2))*x4*x1,1,0},
			{2/(1+pow(x2,2)+pow(x1,2))*x3*x2,-2/(1+pow(x2,2)+pow(x1,2))*x3*x1,0,1}};
		if (N==0) m=m.transpose();
	}
	void run() {
		m.inverse();
	}
};

struct ParseFormTest {
	ConcreteManifold M;
	ex x;
	ParseFormTest(int) : M(10) {}
	
	void run() {
		x=ParseDifferentialForm(M.e(),"2*4+3");
	}
};

int main()
{

	logging_level=LOGGING_LEVEL_ERROR;
	cout<<"Testing the performance of some Wedge functions (higher percentage means better performance):"<<endl;
	//initialization
	bigGlobalFrame=ConcreteManifold(500).e();
	exvector e;
	AddRandomElements(e,200,100,5);
	nonSimpleGlobalFrame.insert(nonSimpleGlobalFrame.begin(),e.begin(),e.end());
	nonSimpleGlobalVSpace.SetBasis(nonSimpleGlobalFrame);

	MEASURE(CurvatureTest, 6, 8322937);
	MEASURE(CurvatureTest, 12, 72291321);
	MEASURE(CurvatureTest, 24, 300474857);
	MEASURE(LeviCivitaConnectionTest, 6, 8265752);
	MEASURE(LeviCivitaConnectionTest, 12, 86354291);
	MEASURE(LeviCivitaConnectionTest, 24, 731307333);
	MEASURE(GetSolutionsTest, 10, 424618);
	MEASURE(GetSolutionsTest, 20, 1838460);
	MEASURE(GetSolutionsTest, 80, 58572142);
	MEASURE(SubspaceFromEquationsTest, 10, 4042961);
	MEASURE(SubspaceFromEquationsTest, 20, 27219283);
	MEASURE(SubspaceFromEquationsTest, 40, 175638583);
	MEASURE(SubspaceFromEquationsOrthTest, 10, 6871606);
	MEASURE(SubspaceFromEquationsOrthTest, 20, 36884727);
	MEASURE(SubspaceFromEquationsOrthTest, 40, 284950875);
	MEASURE(SubspaceTest, 10, 1971385);
	MEASURE(SubspaceTest, 20, 14715147);
	MEASURE(SubspaceTest, 40, 98584476);
	MEASURE(SubspaceOrthTest, 10, 4278452);
	MEASURE(SubspaceOrthTest, 20, 23975821);
	MEASURE(SubspaceOrthTest, 40, 187386090);
	MEASURE(ContainsTest, 2, 582722750);
	MEASURE(ContainsTest, 8, 583497000);
	MEASURE(ContainsTest, 20, 582291500);
	MEASURE(ContainsTest, 40, 582453750);
	MEASURE(SparseBasisTest, 20, 511434);
	MEASURE(SparseBasisTest, 40, 6976372);
	MEASURE(SparseBasisTest, 80, 104814000);
	MEASURE(SparseBasisTest2, 30, 638611);
	MEASURE(SparseBasisTest2, 50, 3092140);
	MEASURE(SparseBasisTest2, 100, 29565029);
	MEASURE(DenseBasisTest, 15, 907301);
	MEASURE(DenseBasisTest, 30, 15799960);
	MEASURE(DenseBasisTest, 60, 259607125);
	MEASURE(DenseBasisTest2, 30, 2806172);
	MEASURE(DenseBasisTest2, 50, 13964993);
	MEASURE(DenseBasisTest2, 100, 132707187);
	MEASURE(DualTest, 10, 17830469);
	MEASURE(DualTest, 20, 107793684);
	MEASURE(FormTest, 10, 260393);
	MEASURE(FormTest, 30, 593740);
	MEASURE(FormTest, 70, 1045252);
	MEASURE(FormTest2, 10, 25411);
	MEASURE(FormTest2, 30, 316325);
	MEASURE(FormTest2, 70, 2639846);
	MEASURE(dTest, 4, 103676);
	MEASURE(dTest, 8, 363578);
	MEASURE(dTest, 16, 933804);
	MEASURE(HookTest, 10, 43357);
	MEASURE(HookTest, 30, 335022);
	MEASURE(HookTest, 70, 2708648);
	MEASURE(RiemannianTest, 10, 1588821);
	MEASURE(RiemannianTest, 30, 14153626);
	MEASURE(RiemannianTest, 80, 123813235);
	MEASURE(ConnectionTest, 10, 805904);
	MEASURE(ConnectionTest, 30, 33329737);
	MEASURE(ConnectionTest, 50, 209037400);
	MEASURE(ConnectionTest2, 10, 2415);
	MEASURE(ConnectionTest2, 50, 2418);
	MEASURE(ConnectionTest2, 80, 2416);
	long perf=0;
	for (vector<Performance>::const_iterator i=performance.begin();i!=performance.end();++i)
		perf+=i->perf();
	cout<<"Average performance "<< perf/performance.size()<<"%"<<endl;
	cout<<"Empty test performance = "<<Measure<EmptyTest>(1,1038).perf()<<"%"<<endl;
	cout<<"Testing the performance of some GiNaC functions (higher percentage means better performance):"<<endl;
	MEASURE(ncmulTest, 0, 174521);
	MEASURE(matrixTest<matrix>, 0, 429320400);
	MEASURE(matrixTest<matrix>, 1, 1193580000);
	MEASURE(matrixTest<GinacLinAlgAlgorithms::InverseMatrix>, 0, 1027710);
	MEASURE(matrixTest<GinacLinAlgAlgorithms::InverseMatrix>, 1, 1210050000);
	return 0;
}

