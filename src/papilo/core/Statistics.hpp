/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*               This file is part of the program and library                */
/*    PaPILO --- Parallel Presolve for Integer and Linear Optimization       */
/*                                                                           */
/* Copyright (C) 2020  Konrad-Zuse-Zentrum                                   */
/*                     fuer Informationstechnik Berlin                       */
/*                                                                           */
/* This program is free software: you can redistribute it and/or modify      */
/* it under the terms of the GNU Lesser General Public License as published  */
/* by the Free Software Foundation, either version 3 of the License, or      */
/* (at your option) any later version.                                       */
/*                                                                           */
/* This program is distributed in the hope that it will be useful,           */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of            */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             */
/* GNU Lesser General Public License for more details.                       */
/*                                                                           */
/* You should have received a copy of the GNU Lesser General Public License  */
/* along with this program.  If not, see <https://www.gnu.org/licenses/>.    */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef _CORE_STATISTICS_HPP_
#define _CORE_STATISTICS_HPP_

namespace papilo
{

struct Statistics
{
   double presolvetime = 0;
   int ntsxapplied{0};
   int ntsxconflicts{0};
   int nboundchgs{0};
   int nsidechgs{0};
   int ncoefchgs{0};
   int nrounds{0};
   int ndeletedcols{0};
   int ndeletedrows{0};
};

inline Statistics
operator-( const Statistics& a, const Statistics& b )
{
   return Statistics{0.0,
                     a.ntsxapplied - b.ntsxapplied,
                     a.ntsxconflicts - b.ntsxconflicts,
                     a.nboundchgs - b.nboundchgs,
                     a.nsidechgs - b.nsidechgs,
                     a.ncoefchgs - b.ncoefchgs,
                     a.nrounds - b.nrounds,
                     a.ndeletedcols - b.ndeletedcols,
                     a.ndeletedrows - b.ndeletedrows};
}

} // namespace papilo

#endif