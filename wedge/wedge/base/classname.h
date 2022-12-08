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


#ifndef WEDGE_CLASSNAME_H
#define WEDGE_CLASSNAME_H
/** @ingroup Base */
 
/** @{ 
 * @file classname.h
 * @brief Template function to retrieve a class name from a %GiNaC or %Wedge class
 */
 
namespace Wedge {

template<typename WedgeType,typename = decltype(WedgeType::static_class_name())> const char* wedgeClassName=WedgeType::static_class_name();
template<typename GiNaCType> const char* className=wedgeClassName<GiNaCType>;

template<> inline const char* className<print_csrc_double> = "print_csrc_double";
template<> inline const char* className<print_csrc_float> = "print_csrc_float";
template<> inline const char* className<print_csrc> = "print_csrc";
template<> inline const char* className<print_tree> = "print_tree";
template<> inline const char* className<print_python_repr> = "print_python_repr";
template<> inline const char* className<print_python> = "print_python";
template<> inline const char* className<print_latex> = "print_latex";
template<> inline const char* className<print_dflt> = "print_dflt";
template<> inline const char* className<print_context> = "print_context";
template<> inline const char* className<wildcard> = "wildcard";
template<> inline const char* className<tensepsilon> = "tensepsilon";
template<> inline const char* className<spinmetric> = "spinmetric";
template<> inline const char* className<minkmetric> = "minkmetric";
template<> inline const char* className<tensmetric> = "tensmetric";
template<> inline const char* className<tensdelta> = "tensdelta";
template<> inline const char* className<tensor> = "tensor";
template<> inline const char* className<symmetry> = "symmetry";
template<> inline const char* className<symbol> = "symbol";
template<> inline const char* className<realsymbol> = "symbol";
template<> inline const char* className<possymbol> = "symbol";
template<> inline const char* className<pseries> = "pseries";
template<> inline const char* className<relational> = "relational";
template<> inline const char* className<power> = "power";
template<> inline const char* className<numeric> = "numeric";
template<> inline const char* className<ncmul> = "ncmul";
template<> inline const char* className<mul> = "mul";
template<> inline const char* className<matrix> = "matrix";
template<> inline const char* className<lst> = "lst";
template<> inline const char* className<integral> = "integral";
template<> inline const char* className<user_defined_kernel> = "user_defined_kernel";
template<> inline const char* className<modular_form_kernel> = "modular_form_kernel";
template<> inline const char* className<Eisenstein_h_kernel> = "Eisenstein_h_kernel";
template<> inline const char* className<Eisenstein_kernel> = "Eisenstein_kernel";
template<> inline const char* className<Kronecker_dz_kernel> = "Kronecker_dz_kernel";
template<> inline const char* className<Kronecker_dtau_kernel> = "Kronecker_dtau_kernel";
template<> inline const char* className<Ebar_kernel> = "Ebar_kernel";
template<> inline const char* className<ELi_kernel> = "ELi_kernel";
template<> inline const char* className<multiple_polylog_kernel> = "multiple_polylog_kernel";
template<> inline const char* className<basic_log_kernel> = "basic_log_kernel";
template<> inline const char* className<integration_kernel> = "integration_kernel";
template<> inline const char* className<indexed> = "indexed";
template<> inline const char* className<spinidx> = "spinidx";
template<> inline const char* className<varidx> = "varidx";
template<> inline const char* className<idx> = "idx";
template<> inline const char* className<GiNaC::function> = "function";
template<> inline const char* className<fderivative> = "fderivative";
template<> inline const char* className<fail> = "fail";
template<> inline const char* className<exprseq> = "exprseq";
template<> inline const char* className<expairseq> = "expairseq";
template<> inline const char* className<constant> = "constant";
template<> inline const char* className<su3d> = "su3d";
template<> inline const char* className<su3f> = "su3f";
template<> inline const char* className<su3t> = "su3t";
template<> inline const char* className<su3one> = "su3one";
template<> inline const char* className<color> = "color";
template<> inline const char* className<diracgammaR> = "diracgammaR";
template<> inline const char* className<diracgammaL> = "diracgammaL";
template<> inline const char* className<diracgamma5> = "diracgamma5";
template<> inline const char* className<diracgamma> = "diracgamma";
template<> inline const char* className<cliffordunit> = "cliffordunit";
template<> inline const char* className<diracone> = "diracone";
template<> inline const char* className<clifford> = "clifford";
template<> inline const char* className<basic> = "basic";
template<> inline const char* className<add> = "add";


}
#endif
