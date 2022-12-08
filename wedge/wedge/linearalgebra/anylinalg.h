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
 #ifndef ANYLINALG_H_
#define ANYLINALG_H_

#include "wedge/base/wedgebase.h"
#include "wedge/convenience/latex.h"

/** @defgroup ExternalAlgorithms Template implementations of linear algebra and polynomial algorithms */ 

/** @{ 
 * @file anylinalg.h
 * @brief Output for types defined in the LinearAlgebraAlgorithms module
 */

namespace Wedge
{ 
namespace internal {
/** @brief Overloaded output operator for matrices
 *  @param os An output stream
 *  @param m A matrix
 *  @return A reference to os
 */
template<typename LinAlgAlgorithms>
std::ostream& Output(std::ostream& os,  const typename LinAlgAlgorithms::IndependenceMatrix& m)
{
	bool TeX=IsStreamTeX(os);
	assert (m.cols()>0 && m.rows()>0);
	if (TeX) os<<"\\begin{pmatrix}"<<std::endl;
	int r=0;
	while (true) {
		int c=0;
		os<<m.M(r,c++);
		while (c<m.cols())
		{
			if (TeX) os<< "&"; else os<< "\t";
			os<<m.M(r,c++);
		}
		if (++r==m.rows()) break;
		if (TeX) os<< "\\\\";
		os<<std::endl;
	}
	os<<std::endl;
	if (TeX) os<<"\\end{pmatrix}"<<std::endl;
	return os;
}

/** @brief Overloaded output operator for matrices
 *  @param os An output stream
 *  @param m A matrix
 *  @return A reference to os
 */

template<typename LinAlgAlgorithms>
std::ostream& Output(std::ostream& os,  const typename LinAlgAlgorithms::InverseMatrix& m)
{
	bool TeX=IsStreamTeX(os);
	assert (m.cols()>0 && m.rows()>0);
	if (TeX) os<<"\\begin{pmatrix}"<<std::endl;
	int r=0;
	while (true) {
		int c=0;
		os<<m(r,c++);
		while (c<m.cols())
		{
			if (TeX) os<< "&"; else os<< "\t";
			os<<m(r,c++);
		}
		if (++r==m.rows()) break;
		if (TeX) os<< "\\\\";
		os<<std::endl;
	}
	os<<std::endl;
	if (TeX) os<<"\\end{pmatrix}"<<std::endl;
	return os;
}


}
} /** @} */
#endif /*ANYLINALG_H_*/
