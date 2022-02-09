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

#include "papilo/core/Problem.hpp"
#include "papilo/core/RowFlags.hpp"
#include "papilo/io/Message.hpp"
#include "papilo/misc/Num.hpp"
#include "papilo/misc/Timer.hpp"
#include <cassert>
#include <cmath>
#include <fstream>

namespace papilo
{

template <typename REAL>
class SingleBoundChange
{
 public:
   int col;
   REAL new_bound_value;
   int reason_row; // row index or -1 if fixing
   bool is_fixing;
   bool is_lower_bound; // if is_fixing && is_lower_bound are false
                        // then it is an upper bound
   int depth_level;
   SingleBoundChange( int _col, REAL _new_bound_value, int _reason_row,
                      bool _is_fixing, bool _is_lower_bound, int _depth_level )
       : col( _col ), new_bound_value( _new_bound_value ),
         reason_row( _reason_row ), is_fixing( _is_fixing ),
         is_lower_bound( _is_lower_bound ), depth_level( _depth_level )
   {
   }
};

template <typename REAL>
class ConflictAnalysis
{
   Message msg;
   Num<REAL> num;
   Timer timer;
   Problem<REAL> problem;

 public:
   ConflictAnalysis( Message _msg, Num<REAL> _num, Timer _timer,
                     Problem<REAL> _problem )
       : msg( _msg ), num( _num ), timer( _timer ), problem( _problem )
   {
   }
   bool
   simple_cut_from_fixings( Vec<SingleBoundChange<REAL>> bound_changes,
                            bool all_fixings_binary, Vec<int>& length,
                            Vec<int*>& indices, Vec<REAL*>& values,
                            Vec<RowFlags>& flags, Vec<REAL>& lhs,
                            Vec<REAL>& rhs )
   {
      length.push_back( 0 );
      lhs.push_back( 1 );
      rhs.push_back( 0 );
      RowFlags rf;
      rf.set( RowFlag::kRhsInf );
      flags.push_back( rf );
      // simplest (but also worse) cut is created from the fixings
      // e.g. Fixings x1 = x2 = x3 = 0 -> x1 + x2 + x3 >= 1
      // e.g. Fixings x1 = x2 = 0, x3 = 1 -> x1 + x2 + (1 - x3) >= 1
      // works only if all fixings are binary!
      if( all_fixings_binary )
      {
         std::vector<int> row_indices;
         std::vector<REAL> row_values;

         for( int i = 0; i < bound_changes.size(); i++ )
         {
            if( bound_changes[i].is_fixing )
            {
               length[0]++;
               row_indices.push_back( bound_changes[i].col );
               double coef =
                   bound_changes[i].new_bound_value > 0.5 ? -1.0 : 1.0;

               if( bound_changes[i].new_bound_value > 0.5 )
                  lhs[0]--;
               row_values.push_back( coef );
            }
         }
         indices.push_back( &row_indices[0] );
         values.push_back( &row_values[0] );
      }
      else
      {
         return false;
      }
      return true;
   }
   bool
   perform_conflict_analysis( Vec<SingleBoundChange<REAL>> bound_changes,
                              bool all_fixings_binary, Vec<int>& length,
                              Vec<int*>& indices, Vec<REAL*>& values,
                              Vec<RowFlags>& flags, Vec<REAL>& lhs,
                              Vec<REAL>& rhs )
   {

      // conf_con = conflict constraint (top of bound_changes)

      // Only fixings (works only for binaries)
      if( bound_changes.back().depth_level == bound_changes.size() - 1 )
      {
         simple_cut_from_fixings( bound_changes, all_fixings_binary, length,
                                  indices, values, flags, lhs, rhs );
         return true;
      }
      else
      {
         while( true )
         {
            SingleBoundChange<REAL> last_bound_change = bound_changes.back();
            bound_changes.pop_back();
         }
         // else BC = bound_changes.pop()
         // if FUIP -> CS.add(BC)
         // else -> resolve_bound_change(BC) (in our case maybe just delete?)
      }

      // should return a list of constraint to be added to the builder
      return true;
   }

   bool
   perform_conflict_analysis()
   {
      msg.info( "function call is dummy and waited to be implemented above" );
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
