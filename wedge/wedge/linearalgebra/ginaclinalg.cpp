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
 #include "ginaclinalg.h"
#include <ginac/matrix.h>
namespace Wedge {
using namespace GiNaC;

/** @brief Pivoting method for matrix elimination schemes. *  
 *  @param ro is the row to be inspected
 *  @param co is the column from where to begin
 *  @return true if there is a non-zero element in the range (r0,c0),...,(r0,n-1), false otherwise
 *  
 * If returns true, columns are swapped so that (r0,c0) is non-zero and the number of zeros below (r0,c0) is maximal
 *  
 * @note Code adapted from GiNaC::matrix
 */
bool GinacLinAlgAlgorithms::IndependenceMatrix::pivot(unsigned ro, unsigned co)
{
	unsigned k = co;
	int max_zeros_below=-1;
	int candidate_col=-1;
	// search first non-zero element in row ro beginning at column co
	while (true) {
		while ((k<n_cols) && (M(ro,k).is_zero()))
			++k;
		if (k==n_cols) break;
		int zeros_below=0;
		for (int i=ro+1;i<n_rows;i++)
			if (M(i,k).is_zero()) ++zeros_below;		
		if (zeros_below>max_zeros_below) {
			candidate_col=k;
			max_zeros_below=zeros_below;
		}
		++k;
	}
	if (candidate_col<0) // all elements in row ro after column co vanish 
		return false;
	if (candidate_col==co)
		// matrix needs no pivoting
		return true;
	// matrix needs pivoting, so swap columns
	SwapColumns(candidate_col,co);
	return true;
};


/** Use division free elimination to determine a maximal set of linearly independent rows.
 *
  */
void GinacLinAlgAlgorithms::IndependenceMatrix::ChooseLinearlyIndependentRows()
{	
	//prepare the matrix expanding everything. Necessary in order for is_zero() to return correct value
	for (int i=0;i<n_rows;i++)
		for (int j=0;j<n_cols;j++)
			M(i,j)=M(i,j).expand();
	ex last_pivot=1;
	unsigned c0=0;
	for (unsigned r0=0; r0<n_rows && c0<n_cols; ++r0) 	
		if (pivot(r0, c0)) {	//is there a non-zero element in the range (r0,c0),...,(r0,n-1)? If so, swap columns to put it in (r0,c0)
			LOG_DEBUG(*this);
			independentRows.push_back(r0);	
			for (unsigned r2=r0+1; r2<n_rows; ++r2) {
				for (unsigned c=c0+1; c<n_cols; ++c)
					M(r2,c)= ((M(r0,c0)*M(r2,c) - M(r2,c0)*M(r0,c)).expand()/last_pivot).expand();
				//M(r2,c0) = 0;
			}
			last_pivot=M(r0,c0);
			++c0;
		}
}
}
