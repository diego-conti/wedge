/***************************************************************************
 *   Copyright (C) 2008 by Diego Conti, diego.conti@unimib.it              *
 *                                                                         *
 *   This file is part of Dictionary.                                      *
 *                                                                         *
 *   Dictionary is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 3 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   Dictionary is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with Dictionary; if not, write to the                           *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include <dictionary/dictionaries/SU3SO3Dictionary.h>
#include <dictionary/dictionaries/SU3SU2Dictionary.h>
#include <dictionary/dictionaries/SU3U2Dictionary.h>
#include <dictionary/ode.h>
#include <dictionary/parsedictionary.h>
#include "latex.h"
#include "parambasis.h"
#include <wedge/polybasis.h>
#include <wedge/coordinates.h>
#include <wedge/submersion.h>
#include <wedge/structures.h>


DECLARE_FUNCTION_1P(f)
DECLARE_FUNCTION_1P(g)
DECLARE_FUNCTION_1P(h)
DECLARE_FUNCTION_1P(k)
DECLARE_FUNCTION_1P(m)
REGISTER_FUNCTION(f, dummy())
REGISTER_FUNCTION(g, dummy())
REGISTER_FUNCTION(h, dummy())
REGISTER_FUNCTION(k, dummy())
REGISTER_FUNCTION(m, dummy())

TrigonometricSimplifier trig;

void Wedge_ab()
{
	cout<<"**********"<<endl;
	SU3SO3Dictionary Y(false);	
	Y.CreateAlgebra();
	Y.PrintForms(cout);
	ex ab=ParseDictionaryElement(Y,"{a,b}");
	cout<<"Wedging forms on SU(3)\\times_SO(3)\\R^3 with ab;"<<endl;
	for (int deg=1;deg<=Y.Dimension();deg++)
	{
		cout<<"Degree "<<deg<<endl<<"\\begin{gather*}"<<endl;
		VectorSpace<CompositeElement> V=Y.p_forms(deg);
		for (int i=1;i<=V.Dimension();i++)
		{
			ex x=ab*V.e(i);
			cout<<x<<" = "<<Y.eval(x)<<"\\\\"<<endl;
		}
		cout<<"\\end{gather*}"<<endl;
	}
	SU3SU2Dictionary X(false);
	X.CreateAlgebra();
	ab=ParseDictionaryElement(X,"{a,b}");	
	cout<<"Wedging forms on SU(3)\\times_SU(2)\\R^3 with ab;"<<endl;
	for (int deg=1;deg<=X.Dimension();deg++)
	{
		cout<<"Degree = "<<deg<<endl;
		VectorSpace<CompositeElement> V=X.p_forms(deg);
		for (int i=1;i<=V.Dimension();i++)
		{
			ex x=ab*V.e(i);
			cout<<x<<" = "<<X.eval(x)<<endl;
		}
	}
}

void Translate()
{
	cout<<latex<<"\\section{Translating dictionaries}"<<endl;
	SU3SU2Dictionary X(false);	
	SU3SO3Dictionary Y(false);	
	X.CreateAlgebra();	
	cout<<"Forms on $SU(3)\\times_SU(2)\\R^3$:"<<endl;
	X.PrintForms(cout);
	Y.CreateAlgebra();
	ex Xab=ParseDictionaryElement(X,"{a,b}");
	ex Yab=X.TranslateTo(Y,Xab);
	cout<<"[ab]_X translates to ["<<Yab<<"]_Y"<<endl;
	Xab=Y.TranslateTo(X,Yab);
	cout<<"["<<Yab<<"]_Y translates to ["<<Xab<<"]_X"<<endl;
	cout<<"Forms on $SU(3)\\times_SU(2)\\R^3$ translate to forms on $SU(3)\\times_SO(3)\\R^3$ as follows:"<<endl;
	X.TranslateTo(Y,cout);
//	cout<<"Forms on SU(3)\\times_SO(3)\\R^3:"<<endl;
//	Y.PrintForms(cout);
	cout<<"Forms on $SU(3)\\times_SO(3)\\R^3$ translate to forms on $SU(3)\\times_SU(2)\\R^3$ as follows:"<<endl;
	Y.TranslateTo(X,cout);
}

//This proves that there are no weakly integrable generalized G2 structures on the sphere bundle inside TCP^2
//see [A. Fino and A. Tomassini - Generalized G_2 manifolds and SU(3) structures]
void ProveNoWeakly()
{
	cout<<latex<<"\\section{Generalized $G_2$ structure onsphere bundle inside $TCP^2$}"<<endl;
	cout<<"We seek an SU(3)-invariant weakly integrable generalized G2 structure on the sphere bundle inside $TCP^2$"<<endl;
	SU3U2Dictionary M(true);
	M.CreateAlgebra();

	VectorSpace<CompositeElement> V=M.p_forms(1);
	ex A=symbol("A"),B=symbol("B"), C=symbol("C");	
	ex alpha=A*V.e(1)+B*V.e(2)+C*V.e(3);
	cout<<Equation("\\alpha",alpha)<<endl;
	
	ex xi = (A*4)/5*V.e(1)+B*V.e(2)+C*V.e(3);
		
	ex omega=M.d(alpha);
	assert(M.Hook(xi,omega).is_zero());

	VectorSpace<CompositeElement> W=M.closed_p_forms(3);
	cout<<"The space of closed 3-forms is"<<endl<<W<<endl;
	
	ex H=W.GenericElement();
	ExVector eqns;
	GetCoefficients<CompositeElement>(eqns,M.eval(H*alpha*omega));
	cout<<"Taking a closed H such that $H\\wedge\\alpha\\wedge\\omega=0$ we get"<<endl;
	cout<<eqns;
	assert(W.Dimension()==3);
	assert(eqns.size()==3);

	matrix associatedmatrix(3,3);	//we have linear equations depending on parameters, we must show there is only one nonzero solution, so we compute the rank of associated matrix
	for (int i=0;i<3;++i)
		for (int j=0;j<W.Dimension();++j)
			associatedmatrix(i,j)=eqns[i].coeff(W.coordinate(j+1));
	cout<<associatedmatrix<<endl;
	cout<<"Associated matrix has characteristic polynomial"<<endl;
	ex lambda=symbol("lambda");
	ex charpoly=associatedmatrix.charpoly(lambda);
	cout<<Equation(charpoly);
	cout<<"divisible by $\\lambda^2$ only if"<<endl;
	cout<<Equation("0",factor(charpoly.coeff(lambda,1)));
	cout<<"and one can check directly that in that case it has rank 2"<<endl;
	assert(W.SubspaceFromEquations(eqns.begin(),eqns.end()).Dimension()==1);	//one solution for any choice of parameters
	H=A*W.SubspaceFromEquations(eqns.begin(),eqns.end()).e(1);
	cout<<"Hence we have one common solution"<<Equation("H",H);

	assert(M.eval(H*omega).is_zero());	//this condition tells us that d\psi_\pm\wedge\alpha=0
	assert(M.Hook(xi,H).is_zero());

	VectorSpace<CompositeElement> ThreeForms=M.p_forms(3);
	ex threeform=ThreeForms.GenericElement();
	
	ParamBasis<VectorSpace<CompositeElement>::Coordinate> Equations;
	GetCoefficients<CompositeElement>(Equations,M.eval(M.d(threeform)*alpha));
	GetCoefficients<CompositeElement>(Equations,M.Hook(xi,threeform)); 
	GetCoefficients<CompositeElement>(Equations,M.eval(omega*threeform));
	//these equations determine the space containing psi^pm

	eqns.clear();
	GetCoefficients<CompositeElement>(eqns,M.eval(threeform*H*alpha));
	assert(eqns.size()==1);
	cout<<Equations.Contains(eqns(1));
	//if true, it means that psi^pm\wedge H =0 necessarily
//	cout<<"Therefore there is no solution."<<endl;
}


void FindHypoContact7()
{
/*
	cout<<"**********"<<endl;
	SU3U2Dictionary M(true);
	M.CreateAlgebra();
	cout<<"We seek a generalized Killing spinor on the sphere bundle in TCP^2 compatible with a fixed contact structure"<<endl;

	VectorSpace<CohomogeneityOneInvariantForms::CompositeElement> V=M.p_forms(1);
	ex A=symbol("A"),B=symbol("B"), C=symbol("C");
	A=ex(1)/2;	
	ex alpha=A*V.e(1)+B*V.e(2)+C*V.e(3);
	cout<<"alpha="<<alpha<<endl;
	
	ex xi = (A*4)/5*V.e(1)+B*V.e(2)+C*V.e(3);
	cout<<"Characteristic vector field ="<<xi<<endl;
	xi = M.AtPrincipalPoint(M.CompositeElementToEx(xi));
		
	ex omega=M.AtPrincipalPoint(M.d(M.CompositeElementToEx(alpha)));
	alpha = M.AtPrincipalPoint(M.CompositeElementToEx(alpha));
	cout<<"Indeed Hook(xi,omega) ="<<Hook(xi,omega)<<endl;
	assert(Hook(xi,omega).is_zero());

	VectorSpace<CohomogeneityOneInvariantForms::CompositeElement> W=M.closed_p_forms(4);
	cout<<"Now the space of closed 4-forms is W="<<W<<endl;
	
	VectorSpace<CohomogeneityOneInvariantForms::CompositeElement> ThreeForms;
	for (int i=1;i<=W.Dimension();i++)
	{
		ex a=M.AtPrincipalPoint(M.CompositeElementToEx(W.e(i)));
		ThreeForms.AddGenerator(M.ExToCompositeElement(Hook(xi,a)));
	}
	cout<<"Since \\Omega\\wedge\\alpha lies in W, it follows that \\Omega is in "<<endl;
	cout<<ThreeForms;
	for (int i=1;i<=ThreeForms.Dimension();i++)
	{
		ex psi=M.AtPrincipalPoint(M.CompositeElementToEx(ThreeForms.e(i)));
		assert(Hook(xi,psi).is_zero());		
	}

	ex psi=M.AtPrincipalPoint(M.CompositeElementToEx(ThreeForms.GenericElement()));
	list<ex> eqns;
	GetCoefficients<DifferentialForm>(eqns,psi*omega);
	Subspace<CohomogeneityOneInvariantForms::CompositeElement> O=ThreeForms.SubspaceFromEquations(eqns.begin(),eqns.end());
	cout<<"Imposing further F\\wedge\\Omega=0, it follows that \\Omega is in "<<endl;
	cout<<O;
*/
}


void STCP2()
{
	cout<<latex<<"\\section{Hypo contact structures on the sphere bundle in $TCP^2$}"<<endl;
	SU3U2Dictionary M(true);	
	M.CreateAlgebra();
	cout<<"Dictionary of SU(3)-invariant forms on the sphere bundle in $TCP^2$"<<endl;
	M.PrintForms(cout);	//print the forms

	ex B=realsymbol("B"),C=realsymbol("C");
	//construct a specific frame 	
	ex beta2=M.beta()(2), beta3=M.beta()(3), b2=M.b(2), b3=M.b(3), delta=B*M.beta()(4)-C*M.beta()(1),gamma=2*M.beta()(4)-C*M.b(4);
	VectorSpace<CompositeElement> V=M.p_forms(1);
	ex alpha=V.e(1)+2*B*V.e(2)+2*C*V.e(3);

	cout<<"We fix the form"<<Equation("\\alpha",alpha)<<endl;
	ex omega=M.d(alpha);
	cout<<"It is a contact form because"<<Equation("\\alpha\\wedge(d\\alpha)^3",M.eval(omega*omega*omega*alpha));

	cout<<"Now consider the contact SU(3)-structure defined by"<<endl;
	//define a contact SU(3)-structure on the sphere bundle
	ex Psi =sqrt(ex(2))*(b2+B*beta3+C*beta2 +I*(b3+beta3* C-beta2* B))*(beta2-I*beta3)
		*(-sqrt(ex(2))*B/C*delta+(B*B+C*C)/C/sqrt(ex(2))*gamma+sqrt(ex(2))*(1+B*B+C*C)*I*delta);		
	
	ex psiplus=M.MakeGlobal(RealPart(Psi).normal());
	ex psiminus=M.MakeGlobal(RealPart(-I*Psi).normal());
	ex F =-M.d(alpha)/2;
	
	cout<<Equation("\\Omega^-",psiminus);
	cout<<Equation("\\Omega^+",psiplus);
	cout<<Equation("\\alpha",alpha);
	cout<<Equation("F",F);

	cout<<"Then"<<endl;
	cout<<Equation("\\Omega^+\\wedge\\alpha",M.eval(psiplus*alpha));
	cout<<Equation("\\Omega^-\\wedge\\alpha",M.eval(psiminus*alpha));
	cout<<Equation("\\Omega^+\\wedge F",M.eval(psiplus*F));
	cout<<Equation("\\Omega^-\\wedge F",M.eval(psiminus*F));	

	cout<<Equation("d\\Omega^+",M.d(psiplus));
	cout<<Equation("d\\Omega^-",M.d(psiminus));

	cout<<Equation("d(\\Omega^+\\wedge\\alpha)",M.d(psiplus*alpha));
	cout<<Equation("d(\\Omega^-\\wedge\\alpha)",M.d(psiminus*alpha));

	cout<<"Now consider the special case B=C=0"<<endl;
	{
		ex beta1=M.beta()(1),beta4=M.beta()(4);
		ex Psi =(b2 +I*b3)*(beta2-I*beta3)*(beta1+I*(beta4))/2;
		ex psiplus=M.MakeGlobal(RealPart(Psi).normal());
		ex psiminus=M.MakeGlobal(RealPart(-I*Psi).normal());
		alpha=alpha.subs(lst(B==0,C==0));
		cout<<Equation("\\Omega^-",psiminus);
		cout<<Equation("\\Omega^+",psiplus);
		cout<<Equation("\\alpha",alpha);
		cout<<Equation("\\Omega^+\\wedge\\alpha",M.eval(psiplus*alpha));
		cout<<Equation("\\Omega^-\\wedge\\alpha",M.eval(psiminus*alpha));
		cout<<Equation("d(\\Omega^+\\wedge\\alpha)",M.d(psiplus*alpha));
		cout<<Equation("d(\\Omega^-\\wedge\\alpha)",M.d(psiminus*alpha));
	}

}



void TCP2()
{
	cout<<latex<<"\\section{Hyperkaehler structure on $TCP^2$}"<<endl;
	SU3U2Dictionary M;	
	M.CreateAlgebra();
	cout<<"Dictionary of SU(3)-invariant forms on $TCP^2$"<<endl;
	M.PrintForms(cout);	//print the forms
	ex omega1=ParseDictionaryElement(M,"-1/4{sigma,b,b}+1/4*r^2{sigma,beta,beta}+1/2{sigma,a,eps}");
	ex omega2=ParseDictionaryElement(M,"-1/2*r{b,beta}-1/(2*r){a,b}{a,beta}");
	ex omega3=ParseDictionaryElement(M,"-1/2*r{sigma,b,beta}-1/(2*r){a,b}{sigma,a,beta}");
	cout<<"The Sp(2)-structure defined by the forms:"<<endl;	
	cout<<Equation("\\omega_1",omega1);
	cout<<Equation("\\omega_2",omega2);
	cout<<Equation("\\omega_3",omega3);
	cout<<"is hyperkaehler because"<<endl;
	cout<<Equation("d\\omega_1",M.d(omega1));
	cout<<Equation("d\\omega_2",M.d(omega2));
	cout<<Equation("d\\omega_3",M.d(omega3));
}



void su3form() {
	AbstractLieGroup<> X("-23-45+2*67,13+46-57-sqrt(3)*58,-12-47+sqrt(3)*48-56,"
					"15-26+37-sqrt(3)*38,"
					"-14+27+36+sqrt(3)*28,-2*17+24-35,2*16-25-34,-sqrt(3)*25+sqrt(3)*34");
	//produce the PSU(3)-invariant 3-form.
	ex gamma;
	for (int i=1;i<=8;++i)
		gamma+=(X.e(i)*X.d(X.e(i))/3).expand();
	ex stargamma=X.HodgeStar(gamma);
	
	cout<<Equation("\\gamma",gamma)<<endl;
	assert(X.d(gamma)==0);
	assert(X.d(stargamma)==0);

	//verify it's biinvariant
	for (int i=1;i<=8;++i)
		assert(X.LieDerivative(X.e(i),gamma)==0);

	SU3SO3Dictionary M;
	M.CreateAlgebra();
	//M.PrintForms(cout);
	ex l=M.RadialCoordinate();
	ExVector e=M.FrameAtPrincipalPoint();
	lst subs;
	list<ex> eqns;

	ex dr=M.d(M.RadialCoordinate());
	Function f(N.f),g(N.g),h(N.h),k(N.k);
	Function df(N.f.Prime()),dg(N.g.Prime()),dh(N.h.Prime()),dk(N.k.Prime());
	M.Declare_d(f,df*dr);
	M.Declare_d(g,dg*dr);
	M.Declare_d(h,dh*dr);
	M.Declare_d(k,dk*dr);
	subs=X.e(2)==f*X.e(2),X.e(3)==f*X.e(3),X.e(4)==f*X.e(4),X.e(5)==f*X.e(5),X.e(6)==h*X.e(6),X.e(7)==h*X.e(7),X.e(8)==k*X.e(8);
	gamma=gamma.subs(subs).subs(lst(f*f==g,f*f*f*f==g*g));	//everything depends on f^2, so set g=f^2
	stargamma=stargamma.subs(subs).subs(lst(f*f==g,f*f*f*f==g*g));

	subs.remove_all();
	subs=(X.e(1)==M.b(1)),(X.e(2)==-1/l*(cos(l)-1)*M.b(3)+1/l*sin(l)*M.b(2)),
		(X.e(3)==1/l*(cos(l)-1)*M.b(2)+1/l*sin(l)*M.b(3)),
		(X.e(4)==(cos(l)+1)*e(1)+sin(l)*e(2)),
		(X.e(5)==(cos(l)+1)*e(2)-sin(l)*e(1)),
		(X.e(6)==(cos(2*l)+1)*e(3)-sin(2*l)*e(4)),
		(X.e(7)==(cos(2*l)+1)*e(4)+sin(2*l)*e(3)),
		(X.e(8)==2*e(5));
	ex x=trig.Simplify(M.MakeGlobal(gamma.subs(subs)));

	//x=M.Rescale(x,Pi/2-1/M.RadialCoordinate());
	cout<<Equation("\\gamma",x);
	cout<<Equation("\\gamma",NormalForm<CompositeElement>(trig.Simplify(x)));
	ex dx=M.d(x).normal();
	cout<<Equation("d\\gamma",dx);
	cout<<Equation("d\\gamma",NormalForm<CompositeElement>(trig.Simplify(dx)));
	GetCoefficients<CompositeElement>(eqns, trig.Simplify(dx));

	ex starx=trig.Simplify(M.MakeGlobal(stargamma.subs(subs)));
	//starx=M.Rescale(starx,Pi/2-1/M.RadialCoordinate());
	cout<<Equation("*\\gamma",starx);
	cout<<Equation("*\\gamma",NormalForm<CompositeElement>(trig.Simplify(starx)));
	ex dstarx=M.d(starx);
	cout<<Equation("d*\\gamma",dstarx)<<endl;
	cout<<Equation("d*\\gamma",NormalForm<CompositeElement>(trig.Simplify(dstarx)));
	GetCoefficients<CompositeElement>(eqns, trig.Simplify(dstarx));

	exvector eqns2;
	for (list<ex>::const_iterator i=eqns.begin();i!=eqns.end();++i)
		eqns2.push_back(collect(*i,M.RadialCoordinate()).numer());
	eqns2[1]-=eqns2[0];	eqns2[3]-=eqns2[4]; eqns2[2]+=eqns2[4]; eqns2[3]+=eqns2[2]; eqns2[1]/=4; eqns2[0]+=eqns2[1]; eqns2[2]/=g; eqns2[3]/=4*h*cos(M.RadialCoordinate()); eqns2[4]+=eqns2[2]*(-g+2*h*h*cos(M.RadialCoordinate()));
	for (exvector::iterator i=eqns2.begin();i!=eqns2.end();++i)	
		*i=factor(i->expand(),factor_options::all);
	for (exvector::iterator i=eqns2.begin();i!=eqns2.end();++i)	
	{
		lst subs;
		subs=(f==symbol("f[r]")),(g==symbol("g[r]")),(h==symbol("h[r]")),(k==symbol("k[r]")),
			(df==symbol("f'[r]")),(dg==symbol("g'[r]")),(dh==symbol("h'[r]")),(dk==symbol("k'[r]"));
	}
	cout<<ExVector(eqns.begin(),eqns.end());
//	PolyBasis_impl<symbol,CocoaPolyAlgorithms> I(eqns2.begin(),eqns2.end());
//	I.Reduce();
	return;
//	x-=-4*sqrt(ex(3))/l*sin(l)*(b2*e(5)-b3*e(4))*e(8)+4*(cos(2*l)+1)*b1*e(6)*e(7)-2*(cos(l)+1)*b1*e(4)*e(5)+4*sin(l)*cos(l)/l*(b2*(e(4)*e(6)-e(5)*e(7))-b3*(e(4)*e(7)+e(5)*e(6)))-2/(l*l)*(1-cos(l))*b1*b2*b3;

	SU3SU2Dictionary N;
	N.CreateAlgebra();
	M.TranslateTo(N,cout);
	ex y=M.TranslateTo(N,x).normal();
	cout<<Equation("\\gamma",NormalForm<CompositeElement>(trig.Simplify(y)));
	cout<<Equation("d\\gamma",NormalForm<CompositeElement>(trig.Simplify(N.d(y))));
	y=M.TranslateTo(N,starx).normal();
	cout<<Equation("*\\gamma",NormalForm<CompositeElement>(trig.Simplify(y)));
	cout<<Equation("d*\\gamma",NormalForm<CompositeElement>(trig.Simplify(N.d(y))));
}

void su3ansatz() {
	AbstractLieGroup<> X("-23-45+2*67,13+46-57-sqrt(3)*58,-12-47+sqrt(3)*48-56,"
					"15-26+37-sqrt(3)*38,"
					"-14+27+36+sqrt(3)*28,-2*17+24-35,2*16-25-34,-sqrt(3)*25+sqrt(3)*34");
	//produce the PSU(3)-invariant 3-form.
	ex gamma;
	for (int i=1;i<=8;++i)
		gamma+=(X.e(i)*X.d(X.e(i))/3).expand();
	ex stargamma=X.HodgeStar(gamma);
	
	cout<<Equation("\\gamma",gamma)<<endl;
	assert(X.d(gamma)==0);
	assert(X.d(stargamma)==0);

	//verify it's biinvariant
	for (int i=1;i<=8;++i)
		assert(X.LieDerivative(X.e(i),gamma)==0);

	SU3SO3Dictionary M;
	M.CreateAlgebra();
	//M.PrintForms(cout);
	ex l=M.RadialCoordinate();
	ExVector e=M.FrameAtPrincipalPoint();
	lst subs;
	list<ex> eqns;

	ex dr=M.d(M.RadialCoordinate());
	Function f(N.f),g(N.g),h(N.h),k(N.k);
	Function df(N.f.Prime()),dg(N.g.Prime()),dh(N.h.Prime()),dk(N.k.Prime());
	M.Declare_d(f,df*dr);
	M.Declare_d(g,dg*dr);
	M.Declare_d(h,dh*dr);
	M.Declare_d(k,dk*dr);
	subs=X.e(2)==f*X.e(2),X.e(3)==f*X.e(3),X.e(4)==f*X.e(4),X.e(5)==f*X.e(5),X.e(6)==h*X.e(6),X.e(7)==h*X.e(7),X.e(8)==k*X.e(8);
	gamma=gamma.subs(subs).subs(lst(f*f==g,f*f*f*f==g*g));	//everything depends on f^2, so set g=f^2
	stargamma=stargamma.subs(subs).subs(lst(f*f==g,f*f*f*f==g*g));

	subs.remove_all();
	subs=(X.e(1)==M.b(1)),(X.e(2)==-1/l*(cos(l)-1)*M.b(3)+1/l*sin(l)*M.b(2)),
		(X.e(3)==1/l*(cos(l)-1)*M.b(2)+1/l*sin(l)*M.b(3)),
		(X.e(4)==(cos(l)+1)*e(1)+sin(l)*e(2)),
		(X.e(5)==(cos(l)+1)*e(2)-sin(l)*e(1)),
		(X.e(6)==(cos(2*l)+1)*e(3)-sin(2*l)*e(4)),
		(X.e(7)==(cos(2*l)+1)*e(4)+sin(2*l)*e(3)),
		(X.e(8)==2*e(5));
	ex x=trig.Simplify(M.MakeGlobal(gamma.subs(subs)));

	cout<<Equation("\\gamma",x);
	cout<<Equation("\\gamma",NormalForm<CompositeElement>(trig.Simplify(x)));
	ex dx=M.d(x).normal();
	cout<<Equation("d\\gamma",dx);
	cout<<Equation("d\\gamma",NormalForm<CompositeElement>(trig.Simplify(dx)));
	GetCoefficients<CompositeElement>(eqns, trig.Simplify(dx));

	ex starx=trig.Simplify(M.MakeGlobal(stargamma.subs(subs)));
	cout<<Equation("*\\gamma",starx);
	cout<<Equation("*\\gamma",NormalForm<CompositeElement>(trig.Simplify(starx)));
	ex dstarx=M.d(starx);
	cout<<Equation("d*\\gamma",dstarx)<<endl;
	cout<<Equation("d*\\gamma",NormalForm<CompositeElement>(trig.Simplify(dstarx)));
	GetCoefficients<CompositeElement>(eqns, trig.Simplify(dstarx));

	ExVector eqns2;
	for (list<ex>::const_iterator i=eqns.begin();i!=eqns.end();++i)
		eqns2.push_back(collect(*i,M.RadialCoordinate()).numer());
	eqns2[1]-=eqns2[0];	eqns2[3]-=eqns2[4]; eqns2[2]+=eqns2[4]; eqns2[3]+=eqns2[2]; eqns2[1]/=4; eqns2[0]+=eqns2[1]; eqns2[2]/=g; eqns2[3]/=4*h*cos(M.RadialCoordinate()); eqns2[4]+=eqns2[2]*(-g+2*h*h*cos(M.RadialCoordinate()));
	for (exvector::iterator i=eqns2.begin();i!=eqns2.end();++i)	
		*i=factor(i->expand(),factor_options::all);
	for (exvector::iterator i=eqns2.begin();i!=eqns2.end();++i)	
	{
		lst subs;
		subs=(f==symbol("f[r]")),(g==symbol("g[r]")),(h==symbol("h[r]")),(k==symbol("k[r]")),
			(df==symbol("f'[r]")),(dg==symbol("g'[r]")),(dh==symbol("h'[r]")),(dk==symbol("k'[r]"));
	}
	cout<<eqns2;
//	PolyBasis_impl<symbol,CocoaPolyAlgorithms> I(eqns2.begin(),eqns2.end());
//	I.Reduce();
	return;
}




/** @brief Perform some demo tasks (in no particular order)
 */
int main()
{
	Wedge_ab();
	SU3SO3Dictionary Y(false);	
	cout<<Y.gamma()<<endl;
	cout<<Y.epsilon()<<endl;
	Y.CreateAlgebra();
	Y.PrintForms(cout);
	TCP2();
	cout<<dflt<<"Switching to normal output"<<endl;
	Wedge_ab();
	STCP2();
	ProveNoWeakly();
	cout<<latex<<"Switching to latex output"<<endl;
	Translate();
	FindHypoContact7();
	return 0;
}
