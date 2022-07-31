/***************************************************************************
 *   Copyright (C) 2007-2022        by Diego Conti, diego.conti@unimib.it  *
 *                                                                         *
 *   This file is part of Wedge.                                           *
 *                                                                         *
 *   Wedge is free software; you can redistribute it and/or modify         *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 3 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   Wedge is distributed in the hope that it will be useful,              *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with Wedge; if not, write to the                                *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#ifndef CANONICAL_PRINT_H
#define CANONICAL_PRINT_H
#include "wedge/base/wedgebase.h"

namespace Wedge {
void canonical_print(ostream& os, ex x);
string to_canonical_string(ex x);
string to_latex_canonical_string(ex x);

string to_string_using(const print_context* pc, ex x, int level=0);
string to_string_using(ostream& os, ex x, int level=0);

}
#endif
