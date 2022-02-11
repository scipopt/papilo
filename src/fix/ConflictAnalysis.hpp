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

#include "papilo/core/RowFlags.hpp"
#include "papilo/core/SingleBoundChange.hpp"
#include "papilo/io/Message.hpp"
#include "papilo/misc/Num.hpp"
#include "papilo/misc/Timer.hpp"
#include "fix/Constraint.hpp"
#include <cassert>
#include <cmath>
#include <fstream>

namespace papilo
{


template <typename REAL>
class ConflictAnalysis
{
   Message msg;
   Num<REAL> num;
   Timer timer;

 public:
   ConflictAnalysis( Message _msg, Num<REAL> _num, Timer timer_ ) : msg( _msg ), num( _num ), timer(timer_)
   {
   }

   bool
   perform_conflict_analysis( Vec<SingleBoundChange<REAL>>& bound_changes, Vec<std::pair<int,int>> infeasible_rows,
                              Vec<Constraint<REAL>>& constraints  )
   {
      // TODO: to be implemented

      // create an empty conflict set CS
      // if only fixings -> no-good-cut / or fix a variable to the other bound?
      // else BC = bound_changes.pop()
      // if FUIP -> CS.add(BC)
      // else -> resolve_bound_change(BC) (in our case maybe just delete?)

      // Maybe use generalized resolution and apply cardinality reduction?

      // should return a list of constraint to be added to the builder
      return true;
   }

   bool
   perform_conflict_analysis()
   {
      msg.info( "function call is dummy and waited to be implemented above\n" );
      return true;
   }

 private:
   bool
   resolve_bound_change()
   {
      return true;
   }

   bool
   get_conflict_set()
   {
      return true;
   }
};

} // namespace papilo
