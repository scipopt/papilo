/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*               This file is part of the program and library                */
/*    PaPILO --- Parallel Presolve for Integer and Linear Optimization       */
/*                                                                           */
/* Copyright (C) 2020-2022 Konrad-Zuse-Zentrum                               */
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

#ifndef _PAPILO_CORE_SINGLE_BOUND_CHANGE_HPP_
#define _PAPILO_CORE_SINGLE_BOUND_CHANGE_HPP_

#include "papilo/core/Problem.hpp"
#include "papilo/misc/Num.hpp"

namespace papilo
{

template <typename REAL>
class SingleBoundChange
{

 public:

   int
   get_col() const
   {
      return col;
   }

   REAL
   get_new_bound_value() const
   {
      return new_bound_value;
   }

   int
   get_reason_row() const
   {
      return reason_row;
   }

   bool
   is_manually_triggered() const
   {
      return manually_triggered;
   }

   bool
   is_lower_bound() const
   {
      return lower_bound;
   }
   int
   get_depth_level() const
   {
      return depth_level;
   }

 private:
   int col;
   int reason_row; // row index or -1 if fixing
   REAL new_bound_value;
   bool manually_triggered;
   bool lower_bound; // if manually_triggered && is_lower_bound are false
                        // then it is an upper bound
   int depth_level;

 public:
   SingleBoundChange( int col_, int row_, REAL new_bound_value_,
                      bool manually_triggered_, bool lower_bound_,
                      int depth_level_ )
       : col( col_ ), reason_row(row_), new_bound_value( new_bound_value_ ),
         manually_triggered( manually_triggered_ ),
         lower_bound( lower_bound_ ), depth_level( depth_level_ )
   {
   }
};

} // namespace papilo
#endif
