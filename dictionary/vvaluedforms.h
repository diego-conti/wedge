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
#ifndef VVALUEDFORMS_H_
#define VVALUEDFORMS_H_
/** @file vvaluedforms.h
 * @brief Vector bundle-valued differential forms; representations
*/

#include <wedge/liegroup.h>
#include <wedge/function.h>
#include <wedge/logging.h>
#include "invariantforms.h"

using namespace Wedge;

/** @brief A differential form taking values in the tautological fibre bundle, i.e. an element of \f$\mathcal{L}\f$.
 *
 * The type Representation describes the fibre V
*/
template<typename Representation> class NamedVValuedForm : public ExVector, public Named {
public:
	NamedVValuedForm(string name="") : ExVector(Representation::FibreSize()), Named(name) {}
	NamedVValuedForm(const Name& name) : ExVector(Representation::FibreSize()), Named(name) {}
	NamedVValuedForm(string name, string texname) : ExVector(Representation::FibreSize()), Named(name,texname) {}		
	NamedVValuedForm(const Name& name, const ExVector& v) : ExVector(v), Named(name) {}
	NamedVValuedForm(string name, string texname, const ExVector& v) : ExVector(v), Named(name,texname) {}

	NamedVValuedForm& set_name(const Name& name) {
		Named::set_name(name);
		return *this;
	}
	NamedVValuedForm& set_name(string name, string texname)
	{
		Named::set_name(name,texname);
		return *this;
	}
};



#endif /*VVALUEDFORMS_H_*/
