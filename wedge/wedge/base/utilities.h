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
#ifndef UTILITIES_H
#define UTILITIES_H

#include "wedge/base/wedgebase.h"

/** @ingroup Base */ 

/** @{ 
 * @file wedge/utilities.h
 * @brief Useful but not indispensable stuff 
 */

namespace Wedge {

using namespace  GiNaC;
using namespace std;
using GiNaC::function; //override std::function
// conversion of numeric types to strings

/** @brief Convert a value to a string
 */
template<class T> string ToString(const T& val)
{
    stringstream strm;
    strm << val;
	string result=strm.str();
	return result;
}

/** @brief Convert a pair of values to a \f$\text{\TeX}\f$-style string
 */
template<class T>
std::string ToString(const T& val1,const T& val2)
{
    	stringstream strm1,strm2;
    	strm1 << val1; strm2<<val2;
	string s1=strm1.str(),s2=strm2.str();
	return (s1.size()>1 || s2.size()>1) ?
		"{"+s1+","+s2+"}" : s1+s2;
}

/** @brief Convert a triple of values to a \f$\text{\TeX}\f$-style string
 */
template<class T>
std::string ToString(const T& val1,const T& val2, const T& val3)
{
    	stringstream strm1,strm2,strm3;
    	strm1 << val1; strm2<<val2; strm3<<val3;
	string s1=strm1.str(),s2=strm2.str(),s3=strm3.str();
	return (s1.size()>1 || s2.size()>1 || s3.size()>1) ?
		"{"+s1+","+s2+","+s3+"}" : s1+s2+s3;
}



/** @brief Extension of GiNaC::visitor allowing recursion
 */
template<typename ResultType> class RecursiveVisitor: public GiNaC::visitor
{
	ResultType res;
public:
	/** @brief The result of the last call to this visitor
	 * @return The result, returned by value
	 */
	inline ResultType GetResult() const {return res;}
	
	/** @brief Use this function instead of accept(*this)
	 * @param e The ex object on which to invoke this visitor
	 * @return The result of the recursive call
	 */
	inline ResultType RecursiveVisit(const GiNaC::ex& e) {
		ResultType oldresult=res;
		e.accept(*this);
		ResultType newresult=res;
		res=oldresult;
		return newresult;
	}
protected:
	/** @brief The result of the last call to this visitor
	 * @return The result as an lvalue	 
	*/
	inline ResultType& Result() {return res;}
};

/** @brief  Helper class to iterate over subsets of a finite set
*/
class IterateOverSubsets {
protected:
/** @brief Pure virtual function called by iterate
 * @param subset An increasing sequence of integers in the range [0,n-1), of length m
 * @return True if iteration should continue, false if iteration should stop
 * 
 * m and n refer to the parameters of iterate
 */
	virtual bool Apply(const vector<int>& subset)=0;
/** @brief Iterate over subsets of {0,..,n-1} of cardinality m
 *  @param m The cardinality of the subset
 *  @param n The size of the set
 
 * Invokes func once for each subset of {0,..,n-1} of cardinality m, until func returns false
*/
	void Iterate(int m, int n);

	virtual ~IterateOverSubsets() {}
};

/** @brief Helper class to iterate over permutations of a finite set
 */
class IterateOverPermutations {
protected:
/** @brief Pure virtual function called by iterate
 * @param permutation A sequence of n distinct integers in the range [0,n-1)
 * @return True if iteration should continue, false if iteration should stop
 * 
 * n refers to the parameter of iterate
 */
 	virtual bool Apply(const std::vector<int>& permutation)=0;
 	
/** @brief Iterate over permutations of {0,..,n-1} 
 *  @param n The size of the set
 *  @param [startFromPermutationNumber,endAtPermutationNumber) A range of permutations
 *  @return false if stopped early
 *  
 * Invokes func once for each permutation of {0,..,n-1} in the range  [startFromPermutationNumber,endAtPermutationNumber),
 * relative to a fixed ordering of the set of permutations, until func returns false.
*/
 	bool Iterate(int n, int startFromPermutationNumber=0, int endAtPermutationNumber=-1);
 		
 	virtual ~IterateOverPermutations() {}
};



/** @brief Simplify expression by rewriting roots and powers of numerics in a normal form.
 *
 * @deprecated This function is not optimized, and will hopefully become unnecessary with newer versions of %GiNaC
 *  
 * @todo pow(100,1/4) -> sqrt(10)
*/

ex NormalizeRoots(ex expression);

} /** @} */

#endif
