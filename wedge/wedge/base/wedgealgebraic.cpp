/*******************************************************************************
 *  Copyright (C) 2007-2022 by Diego Conti, diego.conti@unipi.it 
 *  This file is part of Wedge.                                           
 *  Wedge is free software; you can redistribute it and/or modify         
 *  it under the terms of the GNU General Public License as published by  
 *  the Free Software Foundation; either version 3 of the License, or     
 *  (at your option) any later version.                                   
 *                                                                          
 *  Wedge is distributed in the hope that it will be useful,              
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of        
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         
 *  GNU General Public License for more details.                          
 *                                                                           
 *  You should have received a copy of the GNU General Public License     
 *  along with Wedge; if not, write to the                                
 *   Free Software Foundation, Inc.,                                       
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             
 *  
 *******************************************************************************/
#include "wedgebase.h"
#include "wedgealgebraic.h"

namespace Wedge {
using namespace GiNaC;

#ifdef OLD_TINFO_METHOD
const unsigned TINFO_Vector =0x001ffff0U;
const unsigned TINFO_LambdaVector =TINFO_ncmul|WEDGE_TINFO_INCREMENT;
unsigned internal::TInfoHelper::TINFO_Last= WEDGE_TINFO_INCREMENT;
#endif

#if (GINACLIB_MAJOR_VERSION<=1) && (GINACLIB_MINOR_VERSION<=4)
void Vector::archive(GiNaC::archive_node &n) const {}
GiNaC::ex Vector::unarchive(const GiNaC::archive_node &n, GiNaC::lst &sym_lst)
{
	assert(false);
	return 0;
}
#endif

#if (GINACLIB_MAJOR_VERSION<=1) && (GINACLIB_MINOR_VERSION<=4)
void LambdaVector::archive(GiNaC::archive_node &n) const {}
GiNaC::ex LambdaVector::unarchive(const GiNaC::archive_node &n, GiNaC::lst &sym_lst)
{
	assert(false);
	return 0;
}
#endif



PredefinedNames::PredefinedNames() : alpha("alpha","\\alpha"), beta("beta","\\beta"), gamma ("gamma", "\\gamma"), Gamma ("Gamma", "\\Gamma"), delta ("delta", "\\delta"), 
		Delta("Delta", "\\Delta"), epsilon("epsilon", "\\epsilon"), zeta ("zeta", "\\zeta"),eta("eta", "\\eta"),theta ("theta", "\\theta"),Theta ("Theta", "\\Theta"),
		kappa ("kappa", "\\kappa"),lambda ("lambda", "\\lambda"),Lambda ("Lambda", "\\Lambda"),	mu("mu", "\\mu"),nu("nu", "\\nu"),xi ("xi", "\\xi"),
		Xi ("Xi", "\\Xi"),rho ("rho", "\\rho"),pi ("pi", "\\pi"),Pi ("Pi", "\\Pi"),sigma ("sigma", "\\sigma"),Sigma ("Sigma", "\\Sigma"),
		tau ("tau", "\\tau"),upsilon ("upsilon", "\\upsilon"),Upsilon ("Upsilon", "\\Upsilon"),phi ("phi", "\\phi"),Phi ("Phi", "\\Phi"),
		chi ("chi", "\\chi"),psi ("psi", "\\psi"),Psi ("Psi", "\\Psi"),omega ("omega", "\\omega"),Omega ("Omega", "\\Omega"),a("a"),b("b"),c("c"),d("d"),e("e"),
		f("f"),g("g"),h("h"),i("i"),j("j"),k("k"),l("l"),m("m"),n("n"),o("o"),p("p"),q("q"),r("r"),s("s"),t("t"),u("u"),v("v"),w("w"),x("x"),y("y"),z("z"),
		A("A"),B("B"),C("C"),D("D"),E("E"),F("F"),G("G"),H("H"),I("I"),J("J"),K("K"),L("L"),M("M"),N("N"),O("O"),P("P"),Q("Q"),R("R"),S("S"),T("T"),U("U"),V("V"),W("W"),X("X"),Y("Y"),Z("Z") {}

PredefinedNames N;


NameAndIndex Name::operator()  (int n )const
{
	return NameAndIndex(plain_, _tex,_hasSuperscript)(n); 		
}
NameAndIndex Name::operator()  (int n, int m) const
{
	return NameAndIndex(plain_, _tex,_hasSuperscript)(n,m); 
}
NameAndIndex Name::operator()  (int n, int m, int k) const
{
	return NameAndIndex(plain_, _tex,_hasSuperscript)(n,m, k); 
}
NameAndIndex Name::PlaceHolder() const
{
	return NameAndIndex(plain_,_tex,_hasSuperscript);
}

}
