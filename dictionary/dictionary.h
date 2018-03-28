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
#ifndef DICTIONARY_H_
#define DICTIONARY_H_
/** @file dictionary.h
 * @brief Base class for the actual dictionaries, where the Lie group and representation data appear as parameters
 */

#include "basedictionary.h"
#include "principal.h"
#include "productmanifold.h"
#include "wedge/transverseconnection.h"
#include "wedge/stabilizer.h"
#include "wedge/gl.h"
#include "wedge/repgl.h"
using namespace Wedge;


/** @brief The "dictionary" of invariant forms on an associated bundle \f$G\times_H V\f$
 * 
 * The template parameter Representation describes the fibre \f$V\f$; dimOfStructureGroup is the dimension of \f$H\f$.
 * @sa InvariantForms
*/
template<class Representation,int dimOfStructureGroup> class Dictionary : public BaseDictionary, public Representation {
	Principal<Representation,dimOfStructureGroup> P;
	ex dr;
public:
	ProductManifold GtimesR;
	ex dRadialCoordinate() {return dr;}
	typedef ::NamedVValuedForm<Representation> NamedVValuedForm;	///< The type of "letters", i.e. elements of \f$\mathcal{L}\f$

/** @brief Construct a Dictionary object, without starting the actual calculations
 * @param G A pointer to the Lie group G
 * @param restrictToSphere If true, all computations are carried out restricting to the sphere bundle, i.e. a principal orbit
 * @param upToDegree If specified, forms of degree higher than upToDegree are ignored
 *
 * @warning Caller must ensure that the pointer G remains valid
 */
	Dictionary(const LieGroup* G,  bool restrictToSphere=false, int upToDegree=1000) :
		P(G), GtimesR(*G,P.RadialCoordinate()), invariantForms(upToDegree)
	{
		ex r=P.RadialCoordinate();
		dr=GtimesR.d(r);
		setRadiusSquared(r*r);

		this->restrictToSphere=restrictToSphere;
		representationInfo=Representation::CreateRepresentationInfo(P.a(),r);
		if (restrictToSphere) {
			representationInfo.AppendCondition(AtPrincipalPoint(P.dRadiusSquared()/2)==0);
			representationInfo.genericPoint=representationInfo.principalPoint;
			LOG_INFO(representationInfo.genericPoint);
		}

	}

	/** @brief Initialize the dictionary to aa=f(r)
	 * @param radiussquared The radius squared aa as a function of the radial coordinate
	 * */
	void setRadiusSquared(ex radiussquared) {
		this->radiussquared=radiussquared;
		P.Declare_dr(P.dRadiusSquared()/RadialDerivative(radiussquared));
	}

	ex RadialDerivative(ex function)
	{
		const symbol& r=ex_to<symbol>(RadialCoordinate());
		return function.diff(r);
	}


/** @brief Compute a basis of generators for the "dictionary"
 * [from,to) A range of NamedForm's
*/
	template<typename InputIterator> void CreateAlgebra(InputIterator from, InputIterator to)
	{
		typedef list<typename iterator_traits<InputIterator>::value_type> List;
		//sort the generators by degree in the radius
		List sortedgens;
		while (from!=to)
		{
			typename List::iterator i=sortedgens.begin();
			while (i!=sortedgens.end() && CompareDegree(*i,*from)) ++i;
			sortedgens.insert(i,*from++);
		}
		invariantForms.CreateFromGenerators(sortedgens.begin(),sortedgens.end(), representationInfo);
	
		Declare_d(RadialCoordinate(),ExToCompositeElement(P.dRadiusSquared()/RadialDerivative(radiussquared)));
	}

/** @brief Print out a list of generators and the action of d on those generators
 * @param os The stream to be used for printing, e.g. cout
*/	
	void PrintForms(std::ostream& os)
	{	
		bool latex=IsStreamTeX(os); 
		os<<(latex? "\\[\\begin{array}{|l|l|}" : "{{")<<endl;
		for (int degree=1;degree<Dimension();degree++)
		{
			LOG_INFO(degree);
			if (latex) os<< "\\hline";
			os<<endl;
			VectorSpace<CompositeElement>	compositeGenerators=invariantForms.p_forms(degree);
			for (int i=0;i<compositeGenerators.Dimension();i++)
			{				
				os<<compositeGenerators.e()[i]<<(latex? " & " : " ->");
				os<<d(compositeGenerators.e()[i]);
				if (latex) os<<"\\\\";
				os<<endl;
			}
		}
		os<<(latex? "\\end{array}\\]" : "}}")<<endl<<endl;
	}

/** @brief Translate a single form of this Dictionary to another Dictionary
 * @param D A Dictionary object relative to the same Lie group G
 * @param form A CompositeElement representing a form on this dictionary
 * 
 * Assumes that the principal point is on \{aa=1\}.
 */
	ex TranslateTo(const Dictionary& D,ex form) const
	{
		form=CompositeElementToEx(form);

		VectorSpace<DifferentialForm> hor=P.horizontalForms;
		VectorSpace<DifferentialForm> Dhor=D.P.horizontalForms;
		ex r=RadialCoordinate();
		ex Dr=D.RadialCoordinate();
		ex ab=AtGenericPoint(P.dRadiusSquared()/2).subs(r==1/Dr,subs_options::algebraic);
		ex Dab=D.AtGenericPoint(D.P.dRadiusSquared()/2);
		
		for (int i=1;i<=Representation::FibreSize();++i)
		{
			hor.AddGenerator(AtGenericPoint(this->b(i)).subs(r==1/Dr,subs_options::algebraic));
			Dhor.AddGenerator(D.AtGenericPoint(D.b(i)));
		}
		assert(hor.Dimension()==Dimension());
		assert(Dhor.Dimension()==Dimension());
		
		LOG_INFO(hor);
		LOG_INFO(Dhor);
		
		exvector source, target;		
		for (int i=1;i<=hor.Dimension();++i)
		{			
			ex form=hor.e(i);
			exvector comps=P.e().Components(form);
			list<ex> eqns;
			exvector Dcomps=D.P.e().Components(Dhor.GenericElement());
			for (int i=0;i<P.Dimension()-Representation::FibreSize();i++)
				eqns.push_back(Dcomps[i]-comps[i]);
			//The formula is ab goes to -(aa)^{-2}ab; it is implemented by  
	
			eqns.push_back(
				TrivialPairing<DifferentialForm>(Dhor.GenericElement(),Dab)+
				TrivialPairing<DifferentialForm>(form,ab)
			);
			ex v;
			list<ex> W;
			Dhor.GetSolutions(W,eqns.begin(),eqns.end(),&v);
			assert(W.empty());
			
			source.push_back(form);
			target.push_back(v);
		}
		
		lst subs=LinearMapToSubstitutions<DifferentialForm>(source,target);
		LOG_DEBUG(source);
		LOG_DEBUG(target);
		LOG_DEBUG(subs);
		ex e=AtGenericPoint(form);
		e=e.subs(r==1/Dr,subs_options::algebraic);
		LOG_DEBUG(e);		
		e=e.subs(subs).expand();				
		LOG_DEBUG(e);
		return D.ExToCompositeElement(e);
	}

/** @brief Translate from this Dictionary to another Dictionary
 * @param os The stream to be used for printing, e.g. cout
 * @param D A Dictionary object relative to the same Lie group G
 * 
 * Assumes that the principal point is on \{aa=1\}.
*/	
	void TranslateTo(const Dictionary& D,std::ostream& os) const
	{
		VectorSpace<DifferentialForm> hor=P.horizontalForms;
		VectorSpace<DifferentialForm> Dhor=D.P.horizontalForms;
		ex r=RadialCoordinate();
		ex Dr=D.RadialCoordinate();
		ex ab=AtGenericPoint(P.dRadiusSquared()/2).subs(r==1/Dr,subs_options::algebraic);
		ex Dab=D.AtGenericPoint(D.P.dRadiusSquared()/2);
		
		for (int i=1;i<=Representation::FibreSize();i++)
		{
			hor.AddGenerator(AtGenericPoint(this->b(i)).subs(r==1/Dr,subs_options::algebraic));
			Dhor.AddGenerator(D.AtGenericPoint(D.b(i)));
		}
		assert(hor.Dimension()==Dimension());
		assert(Dhor.Dimension()==Dimension());
		
		LOG_INFO(hor);
		LOG_INFO(Dhor);
		
		exvector source, target;		
		for (int i=1;i<=hor.Dimension();i++)
		{			
			ex form=hor.e(i);
			exvector comps=P.e().Components(form);
			list<ex> eqns;

			exvector Dcomps=D.P.e().Components(Dhor.GenericElement());
			for (int i=0;i<P.Dimension()-Representation::FibreSize();i++)
				eqns.push_back(Dcomps[i]-comps[i]);

			//The formula is ab goes to -(aa)^{-2}ab; it is implemented by  
	
			eqns.push_back(
				TrivialPairing<DifferentialForm>(Dhor.GenericElement(),Dab)+
				TrivialPairing<DifferentialForm>(form,ab)
			);
			ex v;
			list<ex> W;
			Dhor.GetSolutions(W,eqns.begin(),eqns.end(),&v);
			assert(W.empty());
			
			source.push_back(form);
			target.push_back(v);
		}
		
		lst subs=LinearMapToSubstitutions<DifferentialForm>(source,target);
		
		LOG_INFO(source);
		LOG_INFO(target);		
		LOG_INFO(subs);
		
		for (int i=1;i<=P.e().size();i++)
			os<<P.e(i)<<" ->"<<P.e(i).subs(subs)<<endl;
		
		bool latex=IsStreamTeX(os); 
		os<<(latex? "\\[\\begin{array}{|l|l|}" : "{{")<<endl;

		for (int degree=1;degree<Dimension();degree++)
		{
			LOG_INFO(degree);
			if (latex) os<< "\\hline"; 
			os<<endl;
			VectorSpace<CompositeElement>	compositeGenerators=invariantForms.p_forms(degree);
			
			for (int i=0;i<compositeGenerators.Dimension();i++)
			{				
				os<<compositeGenerators.e()[i]<<(latex? " & " : "->");
				ex e=AtGenericPoint(invariantForms.CompositeElementToEx(compositeGenerators.e()[i]));
				e=e.subs(r==1/Dr,subs_options::algebraic);
				LOG_DEBUG(e);			
				e=e.subs(subs).expand();				
				LOG_DEBUG(e);
				try {
					os<<D.ExToCompositeElement(e);
				}
				catch (...)
				{
					os<< "?";
					LOG_WARN(e);	
				}
				if (latex) os<<"\\\\";
				os<<endl;
			}			
		}
		os<<(latex? "\\end{array}\\]" : "}}")<<endl<<endl;
	}
	
	const NamedVValuedForm& a() const {return P.a();} ///< Canonical element of \f$\mathcal{L}\f$
	ex a(OneBased i) const {return a()(i);} ///<The i-th component of \f$ a\colon G\times V\to V\f$
	NamedVValuedForm b() const {return D(a()).set_name(N.b);}  ///< Canonical element of \f$\mathcal{L}\f$
	ex b(OneBased i) const {return b()(i);}  ///< The i-th component of \f$ b\in\Omega^1(G\times V,V)\f$

/** @brief Invariant p-forms
  * @return The vector space of invariant p-forms as composite elements
 */
	VectorSpace<CompositeElement> p_forms(int degree) const
	{
		return invariantForms.p_forms(degree);
	}
	
/** @brief Compute the space of closed forms of a certain degree
 * 
 * @exception NotImplemented if restrictToSphere is false
 * @note This is only implemented if restrictToSphere is true, because Wedge cannot solve differential equations
 */
	VectorSpace<CompositeElement> closed_p_forms(int degree)
	{	
		if (!restrictToSphere) throw NotImplemented(__FILE__,__LINE__);
		if (degree<0 || degree>Dimension()) throw OutOfRange(__FILE__,__LINE__,degree);
		VectorSpace<CompositeElement> forms=p_forms(degree);		
		ex generic_p_form=forms.GenericElement();
		ex dalpha=d(generic_p_form);
		list<ex> eqns;
		GetCoefficients<CompositeElement>(eqns, dalpha);
		return forms.SubspaceFromEquations(eqns.begin(),eqns.end());
	}


/** @brief Compute the space of exact forms of a certain degree
 * 
 * @exception NotImplemented if restrictToSphere is false
 * @note This is only implemented if restrictToSphere is true, because Wedge cannot solve differential equations
 */
	VectorSpace<CompositeElement> exact_p_forms(int degree)
	{
		if (!restrictToSphere) throw NotImplemented(__FILE__,__LINE__);
		VectorSpace<CompositeElement> result;
		if (degree<0 || degree>Dimension()) throw OutOfRange(__FILE__,__LINE__,degree);
		VectorSpace<CompositeElement> forms=p_forms(degree-1);
		for (exvector::const_iterator i=forms.e_begin();i!=forms.e_end();++i)
		{
			result.AddGenerator(d(*i));
		}
		return result;
	}

/** @brief Turn a differential form defined on the slice into a global object
 *
 * @param alpha A DifferentialForm with coefficients functions of r, functions of a_1, or Function's
 * @return A CompositeElement with coefficients functions of r or Function's
*/
	ex MakeGlobal(ex alpha) const {
		return ExToCompositeElement(alpha);
	}

/** @brief Return a basis of one-forms of the associated bundle, computed at the principal point
 */
	Frame FrameAtPrincipalPoint() const
	{
/*		VectorSpace<DifferentialOneForm> T(P.e().begin(),P.e().end()-Representation::FibreSize());
		list<ex> eqns;
		T.GetOrthogonalEquations(eqns,P.connectionForm.begin(),P.connectionForm.end(),TrivialPairingOperator<DifferentialOneForm>());
		VectorSpace<DifferentialOneForm> TX=T.SubspaceFromEquations(eqns.begin(),eqns.end());
*/
		VectorSpace<DifferentialOneForm> TX(P.horizontalForms.e_begin(),P.horizontalForms.e_end());
		for (int i=1;i<=Representation::FibreSize();i++)
			TX.AddGenerator(AtPrincipalPoint(b()(i)));
		return TX.e();
	}

	ExVector PullBackFrame()
	{
		lst subs;
		for (int i=1;i<=GtimesR.Dimension()-1;i++)
			subs.append(P.e(i)==GtimesR.e(i));
		for (int i=1;i<=Representation::FibreSize();i++)
			subs.append(P.d(a(i))==RadialDerivative(AtGenericPoint(a(i)))*dr);
		Frame frame=FrameAtPrincipalPoint(RadialCoordinate());
		ExVector result;
		for (int i=1;i<=frame.size();++i) result.push_back(frame(i).subs(subs));
		return result;
	}
// pull back a CompositeElement to a form  on G\times \R
	ex PullBackCompositeElement (ex compositeElement)
	{
		lst subs;
		for (int i=1;i<=GtimesR.Dimension()-1;i++)
			subs.append(P.e(i)==GtimesR.e(i));
		for (int i=1;i<=Representation::FibreSize();i++)
			subs.append(P.d(a(i))==RadialDerivative(AtGenericPoint(a(i)))*dr);
		ex form= AtGenericPoint(CompositeElementToEx(compositeElement));
		return form.subs(subs);
	}

// take an element on G\times \R and transform it into a compositelement
	ex PushToCompositeElement(ex form)
	{
		lst subs;
		for (int i=1;i<=GtimesR.Dimension()-1;i++)
			subs.append(P.e(i)==GtimesR.e(i));
		for (int i=1;i<=Representation::FibreSize();i++)
			subs.append(P.d(a(i))==RadialDerivative(AtGenericPoint(a(i)))*dr);
		LOG_INFO(subs);
		Frame Pframe(P.horizontalForms.e_begin(),P.horizontalForms.e_end());
		for (int i=1;i<=Representation::FibreSize();i++)
			Pframe.push_back(AtGenericPoint(b()(i)));
		Frame newframe;
		for (int i=1;i<=Pframe.size();++i)
			newframe.push_back(Pframe(i).subs(subs));
		LOG_INFO(newframe);
		subs= LinearMapToSubstitutions<DifferentialForm>(newframe,Pframe);
		LOG_INFO(subs);
		return ExToCompositeElement(form.subs(subs));
	}

	Frame FrameAtPrincipalPoint(ex r) const
	{			
		VectorSpace<DifferentialOneForm> TX(P.horizontalForms.e_begin(),P.horizontalForms.e_end());
		for (int i=1;i<=Representation::FibreSize();i++)
			TX.AddGenerator(AtPrincipalPoint(b()(i),r));
		return TX.e();
	}
	
/** @brief Interior product, or "hook" operator
 * @param left A vector field or differential form
 * @param right A differential form
 * @return The interior product of left into right
 *
 * Here one uses the fact that \f$G\times V\f$ has a canonical basis consisting of the left-invariant basis of $\fG\f$ and the \f$da^i\f$
 */
	ex Hook (ex left, ex right) const
	{
		left=AtGenericPoint(CompositeElementToEx(left));
		right=AtGenericPoint(CompositeElementToEx(right));
		return ExToCompositeElement(Wedge::Hook(left,right));
	}

/** @brief Evaluate a form expressed in terms of CompositeElements
 * @param A CompositeElement with coefficients functions of r or Function's
 * @return A CompositeElement with coefficients functions of r or Function's
 */
	ex eval(ex alpha) const
	{
		return ExToCompositeElement(CompositeElementToEx(alpha));
	}

/** @brief The list of invariant forms on the homogeneous space G/H. 
 * 
 * Overload this function when appropriate.
 */
	list<NamedForm> FormsOnBase() const {return list<NamedForm>();}


/* @brief The dimension of the manifold \f$G\times_H V\f$
 */
	int Dimension() const {return P.Dimension() - dimOfStructureGroup;}

/** @brief Returns a symbolic function that represents the radial coordinate, may equal   \f$r=\sqrt{aa}\f$
 */	
	ex RadialCoordinate()  const {return P.RadialCoordinate();}

	ex d(ex alpha) const {return eval(BaseDictionary::d(alpha));}

/** @brief Rescales a form along the radial direction
 * @param alpha A CompositeElement
 * @param f A function of the radial coordinate, e.g. f=2*RadialCoordinate()
 * @return A CompositeElement obtained by pulling back by a global diffeomorphism obtained extending f invariantly.
 *
 * This function gives a different expression for the same form with respect to a dictionary which is obtained by a radius-rescaling diffeomorphism.
 * For instance if alpha is a function, the result is \f$f^*\alpha\f$.
 */	
	ex Rescale(ex alpha, ex f) const
	{
		//only works if generic point is (r,0,...,0).
		assert(AtGenericPoint(a(1))==RadialCoordinate());
		for (int i=2;i<=Representation::FibreSize();++i)
				assert(AtGenericPoint(a(i))==0);

		ex e=AtGenericPoint(CompositeElementToEx(alpha));
		LOG_INFO(e);
		lst subs;
		assert(is_a<symbol>(RadialCoordinate()));
		const symbol& r=ex_to<symbol>(RadialCoordinate());
		subs.append(P.d(a(1))==diff(f,r)*P.d(a(1)));
		for (int i=2;i<=Representation::FibreSize();++i)
			subs.append(P.d(a(i))==f/RadialCoordinate()*P.d(a(i)));
		
		e=e.subs(RadialCoordinate()==f).subs(subs);
		LOG_INFO(e);
		return ExToCompositeElement(e);
	}
/** @brief Return the i-th horizontal one-form on G
 * @param i A one-based index
 * This function should only be used in conjuction with InvariantStructure
*/
	ex hor(OneBased i) const
	{
		return P.horizontalForms.e(i);
	}
/** @brief Return the i-th vertical one-form
 * @param i A one-based index
 * This function should only be used in conjuction with InvariantStructure
*/
	ex vert(OneBased i) const
	{
		return Representation::Frame(a(),b(),RadialCoordinate())(i).expand();
	}
/** @brief Return a transverse riemannian structure on \f$G\times R\f$
 * @param Frame A basis of horizontal one-forms on \f$G\times R\f$
 *
 * The elements of the frame should be expressed in terms of PullBackFrame
*/
	TransverseRiemannianStructure InvariantStructure(const exvector& frame) const
	{
		VectorSpace<DifferentialForm> special_stabilizer(P.connectionForm.begin(),P.connectionForm.end());

		//the principal stabilizer is the subalgebra of the special stabilizer that fixes the principal point
		//so we compute its annihilator and take a complement
		ExVector omega_a=Representation::LieAlgebraAction(P.connectionForm,a());
		for (int i=1;i<=omega_a.size();++i)
			omega_a(i)=omega_a(i).subs(representationInfo.principalPoint);
		Subspace<DifferentialForm> annihilator=special_stabilizer.Subspace(omega_a.begin(),omega_a.end());
		Subspace<DifferentialForm> complement=annihilator.Complement();
		if (complement.Dimension()!=1) throw WedgeException<std::logic_error>("Principal stabilizer should be one-dimensional",__FILE__,__LINE__);

		return TransverseRiemannianStructure (&GtimesR,frame.begin(),frame.end(),complement.e_begin(),complement.e_end());
	}
	
	AbstractLieSubgroup<false> StabilizerAtPrincipalPoint(ex composite, ex r) {
		//can use the "standard" principal point even if r\neq1
		Frame frame=FrameAtPrincipalPoint();
		ConcreteManifold reference(frame.size());

		lst subs=LinearMapToSubstitutions<DifferentialForm>(frame,reference.e());
		ex form=AtPrincipalPoint(CompositeElementToEx(composite),r);
		cout<<form<<endl;
		LOG_INFO(form);
		GL gl(frame.size());
		GLRepresentation<VectorField> V(&gl,reference.e());
		return Stabilizer<DifferentialForm>(gl,V,form.subs(subs));
	}

	void testCompositeElement(ex form,lst phistarb) {
		form = CompositeElementToEx(form);
		ex atgeneric=AtGenericPoint(form);
		cout<<"at generic point "<<atgeneric<<endl;
		cout<<"dictionary element "<<MakeGlobal(form)<<endl;
		cout<<"d of form ="<<d(MakeGlobal(form))<<endl;
		cout<<"d of ex ="<<NormalForm<DifferentialForm>(AtGenericPoint(P.d(form)))<<endl;
		exvector basvector=b();
		lst subs=LinearMapToSubstitutions<DifferentialForm>(basvector.begin(),basvector.end(),phistarb.begin(),phistarb.end());
		cout<<"phi^*x = "<<NormalForm<DifferentialForm>(atgeneric.subs(subs))<<endl;
		cout<<"phi^* dx= "<<NormalForm<DifferentialForm>(AtGenericPoint(P.d(form)).subs(subs))<<endl;
	}
	void testEx(ex form) {
		ex atgeneric=AtGenericPoint(form);
		cout<<"at generic point "<<atgeneric<<endl;
		cout<<"dictionary element "<<MakeGlobal(form)<<endl;
		cout<<"d of form ="<<d(MakeGlobal(form))<<endl;
		cout<<"d of ex ="<<NormalForm<DifferentialForm>(AtGenericPoint(P.d(form)))<<endl;
	}

protected:
/** @brief Exterior covariant derivative
 * @param alpha A form in \f$\Omega^p(G\times V,V)\f$
 * @return The exterior covariant derivative \f$D\alpha\in \Omega^{p+1}(G\times V,V)\f$
 */ 
	NamedVValuedForm D(const NamedVValuedForm& alpha) const
	{
		return P.D(alpha);
	}
/** @brief Return a reference to the principal bundle object's connection form
 * @return a vector whose components are the connection forms
 * For use in initialization
*/
	ExVector& ConnectionForm() {return P.connectionForm;}
/** @brief Return a reference to a component of the connection form
 * For use in initialization
*/
	ex& ConnectionForm(OneBased i) {
		assert(i>0 && i<=dimOfStructureGroup);
		return P.connectionForm(i);
	}
/** @brief Add a range of one-forms to the set of horizontal forms
 * For use in initialization
*/
	template<typename Iterator> void addHorizontalForms(Iterator begin, Iterator end) {
		P.horizontalForms.AddGenerators(begin,end);
	}

	CohomogeneityOneInvariantForms invariantForms;	///< The space of invariant forms on the associated bundle (computed by CreateAlgebra)
/** @brief Info about the representation V of H (initialized during construction by derived class)
 * @remark In particular, the representation info contains specification of one principal point,
 * and also a "generic" point, which lies on a fixed half-line, and depends on a parameter. These will be referred to as
 * "the" principal and generic points.
 */ 
	CohomogeneityOneInvariantForms::RepresentationInfo representationInfo;
	bool restrictToSphere;	///<True if we are restricting to the sphere bundle

	ex radiussquared;	///< Radius squared as a function of the radial coordinate, i.e. aa=f(r)

/** @brief Conversion routine
  * @param A DifferentialForm with coefficients functions of r or Function's
  * @return composite A CompositeElement with coefficients functions of r or Function's
 */
	ex ExToCompositeElement(ex form) const
	{
		ex r=invariantForms.ExToCompositeElement(AtGenericPoint(form).expand());
		ex aa=invariantForms.InvariantFunction();

		if (restrictToSphere) r=r.subs(aa==1);
		else r=r.subs(aa==radiussquared);
		return r.normal().expand();
	}

/** @brief Conversion routine
  * @param composite A CompositeElement with coefficients functions of r or Function's
  * @return A DifferentialForm with coefficients functions of r or Function's
 */
	ex CompositeElementToEx(ex composite) const
	{
		ex aa=invariantForms.InvariantFunction();
		if (restrictToSphere) composite=composite.subs(aa==1);
		else composite=composite.subs(aa==radiussquared);
		return invariantForms.CompositeElementToEx(composite);
	}

	ex AtPrincipalPoint(ex form) const
	{
		return form.subs(representationInfo.principalPoint).expand();
	}

	ex AtPrincipalPoint(ex form, ex r) const
	{
		return form.subs(representationInfo.genericPoint).subs(RadialCoordinate()==r).expand();
	}

	ex AtGenericPoint(ex form) const
	{
		return form.subs(representationInfo.genericPoint).expand();
	}

	virtual const Manifold* TotalSpace() const {return &P;}

private:	
	ex d(const CompositeElement& alpha) const {
		ex form=invariantForms.CompositeElementToEx(alpha);
		ex dalpha=P.d(form);
		return ExToCompositeElement(dalpha);
	}

// Compare the degree of two forms as polynomials in the coordinates
	bool CompareDegree(ex x, ex y) const
	{
		x=x.subs(representationInfo.genericPoint);
		y=y.subs(representationInfo.genericPoint);
		return x.ldegree(RadialCoordinate())<=y.ldegree(RadialCoordinate());
	}
	ex approximate(ex x, ex r) const {
		LOG_INFO(x);
		ex res=AtPrincipalPoint(x,r);
		ExVector delta(Representation::FibreSize());
		for (int i=1;i<=Representation::FibreSize();++i)
		{
			delta(i)=a(i)-AtPrincipalPoint(a(i),r);
			res+=AtPrincipalPoint(x.diff(ex_to<symbol>(a(i))),r)*delta(i);
		}
		for (int i=1;i<=Representation::FibreSize();++i)
		for (int j=1;j<=Representation::FibreSize();++j)
			res+=AtPrincipalPoint(x.diff(ex_to<symbol>(a(i))).diff(ex_to<symbol>(a(j))),r)*delta(i)*delta(j)/2;
		res=res.expand();
		LOG_INFO(res);
		return res;
	}
};

#endif /*DICTIONARY_H_*/
