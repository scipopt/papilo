
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

#ifndef _PAPILO_MISC_KKT_CHECK_HELPER_HPP_
#define _PAPILO_MISC_KKT_CHECK_HELPER_HPP_

#include <iostream>

#include "papilo/core/Problem.hpp"
#include "papilo/io/Message.hpp"

namespace papilo
{

template <typename REAL>
struct KktRuleInfo
{
   int num_violated = 0;
   REAL max = -1;
   REAL sum = 0;

   std::vector<REAL> values;
};

template <typename REAL>
void
updateRuleInfo( KktRuleInfo<REAL>& info )
{
   if( info.values.size() == 0 )
      return;
   info.num_violated = info.values.size();
   info.max = info.values[0];
   info.sum = info.values[0];
   for( int i = 1; i < info.values.size(); i++ )
   {
      info.sum += info.values[i];
      if( info.max < info.values[i] )
         info.max = info.values[i];
   }
}

template <typename REAL>
struct KktInfo
{
   int num_col;
   int num_row;

   // Primal.
   KktRuleInfo<REAL> primal_col_bounds;
   KktRuleInfo<REAL> primal_row_bounds;

   // todo: dual.
};

template <typename REAL>
void
updateKktInfo( const int cols, const int rows, KktInfo<REAL>& info )
{
   info.num_col = cols;
   info.num_row = rows;

   // Primal.
   updateRuleInfo( info.primal_row_bounds );
   updateRuleInfo( info.primal_row_bounds );

   // todo: dual
}

template <typename REAL>
void
addRow( Problem<REAL>& problem, const int row, const int length,
        const REAL* values, const int* coeffs, const REAL lhs, const REAL rhs,
        const bool lb_inf, const bool ub_inf )
{
   // Assuming problem is expanded and row will be added at the preacllocated
   // space.

   // Modify lhs, rhs.
   // Modify rowFlags and colFlags.

   // Modify matrix.
}

// todo:
template <typename REAL>
void
addColToProblem()
{
}

} // namespace papilo

#endif