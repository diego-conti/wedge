/*******************************************************************************
 *  Copyright (C) 2007-2023 by Diego Conti, diego.conti@unipi.it 
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
 
#include "liegrouptostring.h"
#include "wedge/convenience/horizontal.h"

namespace Wedge {

map<pair<int,int>,ex> coefficients_in_lie_group (ex de_k, const LieGroup& G) {
    map<pair<int,int>,ex> coeffs;    
    for (int i=1;i<=G.Dimension();++i)
        for (int j=i+1;j<=G.Dimension();++j) {
            ex coeff=Hook(G.e(i)*G.e(j),de_k).expand();
            if (!coeff.is_zero()) coeffs[make_pair(i,j)]	=coeff;
        }                
    return coeffs;
}

string coefficient_to_string(ex x) {
	stringstream s;
	s<<dflt;
	if (is_a<numeric>(x) && is_rational(ex_to<numeric>(x))) s<<x;	
	else s<<"["<<x<<"]";
	return s.str();	
}

void append_digit(ostream& os, int n) {
	if (n>=0 && n<=9) os<<n;
	else {
        char c='a'+static_cast<char>(n)-10;
        os<<c;
    }
}

string coefficient_and_form(string coefficient, int i, int j) {	
	stringstream s;
	if (coefficient=="-1") s<<"-";
	else if (coefficient!="1") s<<coefficient<<"*";
	append_digit(s,i);
	append_digit(s,j);
	return s.str();
}

void append_term(stringstream& os, const string& coefficient, int i, int j) {
	string term=coefficient_and_form(coefficient,i,j);
	if (os.str().empty() || term[0]=='-') os<<term;
	else os<<"+"<<term;
}

string coefficients_to_string(const map<pair<int,int>,ex>& coeffs) {
	stringstream s;
	for (auto indices_and_coeff : coeffs) {
		auto ij=indices_and_coeff.first;
		auto coeff=indices_and_coeff.second;		
		string coeff_string=coefficient_to_string(coeff);		
		append_term(s,coeff_string,ij.first, ij.second);			
	}
	if (s.str().empty()) s<<"0";
	return s.str();
}

string lie_group_to_string(const LieGroup& G) {
    vector<string> de;
    for (int k=1;k<=G.Dimension();++k) {
        ex de_k=G.d(G.e(k));
        auto coeffs=coefficients_in_lie_group(de_k,G);        
        de.push_back(coefficients_to_string(coeffs));
    }
    return horizontal(de);
}
}