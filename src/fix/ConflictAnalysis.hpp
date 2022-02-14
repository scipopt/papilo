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

#include "fix/Constraint.hpp"
#include "papilo/core/Problem.hpp"
#include "papilo/core/RowFlags.hpp"
#include "papilo/core/SingleBoundChange.hpp"
#include "papilo/core/SparseStorage.hpp"
#include "papilo/io/Message.hpp"
#include "papilo/misc/Num.hpp"
#include "papilo/misc/Timer.hpp"
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
   Problem<REAL>& problem;

 public:
   ConflictAnalysis( Message _msg, Num<REAL> _num, Timer _timer,
                     Problem<REAL>& _problem )
       : msg( _msg ), num( _num ), timer( _timer ), problem( _problem )
   {
   }
   bool
   perform_conflict_analysis( Vec<SingleBoundChange<REAL>>& bound_changes,
                              Vec<std::pair<int, int>>& infeasible_rows,
                              Vec<Constraint<REAL>>& constraints )
   {
      Vec<int> decision_levels;
      Vec<int> pos_in_bound_changes;
      Vec<int> conflict_set_candidates;
      Vec<bool> in_candidates;

      // bound change data as vectors for easier access
      decision_levels.resize( problem.getNCols(), 0 );
      pos_in_bound_changes.resize( problem.getNCols(), -1 );
      in_candidates.resize( problem.getNCols(), 0 );

      // row that led to infeasibility
      int conflict_row_index = infeasible_rows[0].second;
      // position in stack when infeasible
      int pos_at_infeasibility = infeasible_rows[0].first;

      int fixings = 0;
      // ToDo what about the integer case?
      for( int i = 0; i <= pos_at_infeasibility; i++ )
      {
         // consider only the latest bound change for each variable
         // (for binaries there is only one bound change possible)
         int col = bound_changes[i].get_col();
         // the decision level of the column as positive integers (compatible
         // for conflict analysis)
         int col_dec_level = abs( bound_changes[i].get_depth_level() ) - 1 - fixings;
         assert(col_dec_level > 0);
         if (bound_changes[i].is_manually_triggered()) fixings++;
         decision_levels[col] = col_dec_level;
         pos_in_bound_changes[col] = i;
      }
      int number_fixings = get_number_fixings( bound_changes );
      // Only fixings (works only for binaries)
      if( number_fixings == bound_changes.size() )
      {
         simple_cut_from_fixings( bound_changes, conflict_set_candidates,
                                  constraints );
         return true;
      }
      else
      {

         // Find subset of indices that explain the infeasibility
         // adds column indices in conflict_set_candidates
         explain_infeasibility( bound_changes, pos_in_bound_changes,
                                conflict_set_candidates, in_candidates,
                                conflict_row_index );
         // last depth level
         int last_depth_level = get_last_depth_level( bound_changes, conflict_set_candidates, decision_levels );

         int num_vars_last_depth_level = get_number_variables_depth_level(
             decision_levels, conflict_set_candidates, last_depth_level );

         if( num_vars_last_depth_level == 1 )
         {
            // Already at First UIP -> return conflict constraint
            // ToDo add constraint
            add_constraint( bound_changes, pos_in_bound_changes,
                            conflict_set_candidates, constraints );
            msg.info( "Only one variable at last depth level! \n" );
            return true;
         }
         // First-FUIP
         // Resolve as long as more than one bound changes at last depth level
         while( num_vars_last_depth_level > 1 )
         {
            //
            int col_index;
            int position_in_conflict_set;

            get_latest_col_index_in_depth_level(
                decision_levels, pos_in_bound_changes, conflict_set_candidates,
                last_depth_level, col_index, position_in_conflict_set );

            int antecedent_row_index =
                bound_changes[pos_in_bound_changes[col_index]].get_reason_row();
            if( antecedent_row_index == -1 )
            {
               // should not happen
               // ToDo add assert
               return false;
            }
            else
            {

               // remove latest_col_index
               in_candidates[col_index] = 0;
               conflict_set_candidates.erase( conflict_set_candidates.begin() +
                                              position_in_conflict_set );

               // resolve
               explain_infeasibility( bound_changes, pos_in_bound_changes,
                                      conflict_set_candidates, in_candidates,
                                      antecedent_row_index );
            }

            num_vars_last_depth_level = get_number_variables_depth_level(
                decision_levels, conflict_set_candidates, last_depth_level );
         }
      }

      add_constraint( bound_changes, pos_in_bound_changes,
                      conflict_set_candidates, constraints );
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
   is_rhs_reason( Vec<SingleBoundChange<REAL>>& bound_changes,
                  Vec<int>& pos_in_bound_changes, int row_idx )
   {
      // get conflict row
      SparseVectorView<REAL> conflict_row =
          problem.getConstraintMatrix().getRowCoefficients( row_idx );
      int row_length = conflict_row.getLength();
      const int* row_inds = conflict_row.getIndices();
      const REAL* row_vals = conflict_row.getValues();

      REAL rhs = problem.getConstraintMatrix().getRightHandSides()[row_idx];

      REAL min_activity = problem.getRowActivities()[row_idx].min;

      for( int i = 0; i < row_length; i++ )
      {
         int pos = pos_in_bound_changes[row_inds[i]];
         if( row_vals[i] > 0 && pos > 0 )
         {
            if( bound_changes[pos].get_new_bound_value() >
                problem.getLowerBounds()[row_vals[i]] )
            {
               min_activity +=
                   row_vals[i] * ( bound_changes[pos].get_new_bound_value() -
                                   problem.getLowerBounds()[row_vals[i]] );
            }
         }
         else if( row_vals[i] < 0 )
         {
            if( bound_changes[pos].get_new_bound_value() <
                problem.getUpperBounds()[row_vals[i]] )
            {
               min_activity -=
                   row_vals[i] * ( problem.getLowerBounds()[row_vals[i]] -
                                   bound_changes[pos].get_new_bound_value() );
            }
         }
         if( min_activity - rhs > 0 )
            return true;
      }
      return false;
   }

   bool
   is_lhs_reason( Vec<SingleBoundChange<REAL>>& bound_changes,
                  Vec<int>& pos_in_bound_changes, int row_idx )
   {
      // get conflict row
      SparseVectorView<REAL> conflict_row =
          problem.getConstraintMatrix().getRowCoefficients( row_idx );
      int row_length = conflict_row.getLength();
      const int* row_inds = conflict_row.getIndices();
      const REAL* row_vals = conflict_row.getValues();

      REAL lhs = problem.getConstraintMatrix().getLeftHandSides()[row_idx];

      REAL max_activity = problem.getRowActivities()[row_idx].max;
      for( int i = 0; i < row_length; i++ )
      {
         int pos = pos_in_bound_changes[row_inds[i]];
         if( row_vals[i] > 0 && pos > 0 )
         {
            if( bound_changes[pos].get_new_bound_value() <
                problem.getUpperBounds()[row_vals[i]] )
            {
               max_activity -=
                   row_vals[i] * ( problem.getUpperBounds()[row_vals[i]] -
                                   bound_changes[pos].get_new_bound_value() );
            }
         }
         else if( row_vals[i] < 0 )
         {
            if( bound_changes[pos].get_new_bound_value() >
                problem.getLowerBounds()[row_vals[i]] )
            {
               max_activity +=
                   row_vals[i] * ( bound_changes[pos].get_new_bound_value() -
                                   problem.getLowerBounds()[row_vals[i]] );
            }
         }
         if( max_activity - lhs < 0 )
            return true;
      }
      return false;
   }

   bool
   explain_infeasibility( Vec<SingleBoundChange<REAL>>& bound_changes,
                          Vec<int>& pos_in_bound_changes,
                          Vec<int>& conflict_set_candidates,
                          Vec<bool>& in_candidates, int row_idx )
   {
      RowFlags row_flag = problem.getConstraintMatrix().getRowFlags()[row_idx];

      // For GE constraints
      if( row_flag.test( RowFlag::kRhsInf ) )
      {
         explain_infeasibility_ge( bound_changes, pos_in_bound_changes,
                                   conflict_set_candidates, in_candidates,
                                   row_idx );
      }
      // For LE constraints
      else if( row_flag.test( RowFlag::kLhsInf ) )
      {
         explain_infeasibility_le( bound_changes, pos_in_bound_changes,
                                   conflict_set_candidates, in_candidates,
                                   row_idx );
      }
      // For equalities or ranged rows
      else
      {
         if( is_lhs_reason( bound_changes, pos_in_bound_changes, row_idx ) )
         {
            explain_infeasibility_ge( bound_changes, pos_in_bound_changes,
                                      conflict_set_candidates, in_candidates,
                                      row_idx );
         }
         else if( is_rhs_reason( bound_changes, pos_in_bound_changes,
                                 row_idx ) )
         {
            explain_infeasibility_le( bound_changes, pos_in_bound_changes,
                                      conflict_set_candidates, in_candidates,
                                      row_idx );
         }
      }

      return true;
   }
   bool
   explain_infeasibility_ge( Vec<SingleBoundChange<REAL>>& bound_changes,
                             Vec<int>& pos_in_bound_changes,
                             Vec<int>& conflict_set_candidates,
                             Vec<bool>& in_candidates, int row_idx )
   {
      // get conflict row
      SparseVectorView<REAL> conflict_row =
          problem.getConstraintMatrix().getRowCoefficients( row_idx );
      int row_length = conflict_row.getLength();
      const int* row_inds = conflict_row.getIndices();
      const REAL* row_vals = conflict_row.getValues();

      REAL lhs = problem.getConstraintMatrix().getLeftHandSides()[row_idx];

      REAL max_activity = problem.getRowActivities()[row_idx].max;
      for( int i = 0; i < row_length; i++ )
      {
         int pos = pos_in_bound_changes[row_inds[i]];
         if( row_vals[i] > 0 && pos > 0 )
         {
            if( bound_changes[pos].get_new_bound_value() <
                problem.getUpperBounds()[row_vals[i]] )
            {
               max_activity -=
                   row_vals[i] * ( problem.getUpperBounds()[row_vals[i]] -
                                   bound_changes[pos].get_new_bound_value() );
               if( !in_candidates[row_inds[i]] )
               {
                  in_candidates[row_inds[i]] = 1;
                  conflict_set_candidates.push_back( row_inds[i] );
               }
            }
         }
         else if( row_vals[i] < 0 )
         {
            if( bound_changes[pos].get_new_bound_value() >
                problem.getLowerBounds()[row_vals[i]] )
            {
               max_activity +=
                   row_vals[i] * ( bound_changes[pos].get_new_bound_value() -
                                   problem.getLowerBounds()[row_vals[i]] );
               if( !in_candidates[row_inds[i]] )
               {
                  in_candidates[row_inds[i]] = 1;
                  conflict_set_candidates.push_back( row_inds[i] );
               }
            }
         }
         if( max_activity - lhs < 0 )
            break;
      }
      return true;
   }

   bool
   explain_infeasibility_le( Vec<SingleBoundChange<REAL>>& bound_changes,
                             Vec<int>& pos_in_bound_changes,
                             Vec<int>& conflict_set_candidates,
                             Vec<bool>& in_candidates, int row_idx )
   {
      // get conflict row
      SparseVectorView<REAL> conflict_row =
          problem.getConstraintMatrix().getRowCoefficients( row_idx );
      int row_length = conflict_row.getLength();
      const int* row_inds = conflict_row.getIndices();
      const REAL* row_vals = conflict_row.getValues();

      REAL rhs = problem.getConstraintMatrix().getRightHandSides()[row_idx];

      REAL min_activity = problem.getRowActivities()[row_idx].min;
      for( int i = 0; i < row_length; i++ )
      {
         int pos = pos_in_bound_changes[row_inds[i]];
         if( row_vals[i] > 0 && pos > 0 )
         {
            if( bound_changes[pos].get_new_bound_value() >
                problem.getLowerBounds()[row_vals[i]] )
            {
               min_activity +=
                   row_vals[i] * ( bound_changes[pos].get_new_bound_value() -
                                   problem.getLowerBounds()[row_vals[i]] );
               if( !in_candidates[row_inds[i]] )
               {
                  in_candidates[row_inds[i]] = 1;
                  conflict_set_candidates.push_back( row_inds[i] );
               }
            }
         }
         else if( row_vals[i] < 0 )
         {
            if( bound_changes[pos].get_new_bound_value() <
                problem.getUpperBounds()[row_vals[i]] )
            {
               min_activity -=
                   row_vals[i] * ( problem.getLowerBounds()[row_vals[i]] -
                                   bound_changes[pos].get_new_bound_value() );
               if( !in_candidates[row_inds[i]] )
               {
                  in_candidates[row_inds[i]] = 1;
                  conflict_set_candidates.push_back( row_inds[i] );
               }
            }
         }
         if( min_activity - rhs > 0 )
            break;
      }
      return true;
   }

   int
   get_last_depth_level( Vec<SingleBoundChange<REAL>>& bound_changes,Vec<int>& conflict_set_candidates,Vec<int>& decision_levels )
   {
      int level = 0;
      for( int i = 0; i < conflict_set_candidates.size(); i++ )
      {
         level = decision_levels[conflict_set_candidates[i]] > level ? decision_levels[conflict_set_candidates[i]] : level; 
      }
      return level;
   }

   int
   get_number_variables_depth_level( Vec<int>& decision_levels,
                                     Vec<int>& conflict_set_candidates,
                                     int depth )
   {
      int number_variables_depth_level = 0;
      for( int i = 0; i < conflict_set_candidates.size(); i++ )
      {
         if( decision_levels[conflict_set_candidates[i]] == depth )
            number_variables_depth_level++;
      }
      return number_variables_depth_level;
   }
   void
   get_latest_col_index_in_depth_level( Vec<int>& decision_levels,
                                        Vec<int>& pos_in_bound_changes,
                                        Vec<int>& conflict_set_candidates,
                                        int depth, int& col,
                                        int& position_in_conflict_set )
   {

      int latest_position = -1;
      for( int i = 0; i < conflict_set_candidates.size(); i++ )
      {
         if( ( decision_levels[conflict_set_candidates[i]] == depth ) &&
             ( pos_in_bound_changes[conflict_set_candidates[i]] >
               latest_position ) )
         {
            position_in_conflict_set = i;
            col = conflict_set_candidates[i];
            latest_position = pos_in_bound_changes[conflict_set_candidates[i]];
         }
      }
      return;
   }

   int
   get_number_fixings( Vec<SingleBoundChange<REAL>>& bound_changes )
   {
      int number_fixings = 0;
      // in papilo the first depth is -2
      int curr_depth = 0;
      for( int i = 0; i < bound_changes.size(); i++ )
      {
         if( bound_changes[i].is_manually_triggered() &&
             ( curr_depth > bound_changes[i].get_depth_level() ) )
            number_fixings++;
         curr_depth = bound_changes[i].get_depth_level();
      }
      return number_fixings;
   }

   void
   add_constraint( Vec<SingleBoundChange<REAL>>& bound_changes,
                   Vec<int>& pos_in_bound_changes,
                   Vec<int>& conflict_set_candidates,
                   Vec<Constraint<REAL>>& constraints )
   {
      RowFlags row_flag;
      row_flag.set( RowFlag::kRhsInf );

      REAL added_con_vals[conflict_set_candidates.size()];
      int added_con_inds[conflict_set_candidates.size()];
      REAL lhs = 1.0;
      for( int i = 0; i < conflict_set_candidates.size(); i++ )
      {
         int col = conflict_set_candidates[i];
         int pos = pos_in_bound_changes[col];
         REAL coef =
             bound_changes[pos].get_new_bound_value() > 0.5 ? -1.0 : 1.0;
         if( coef < 0 )
            lhs--;
         added_con_inds[i] = col;
         added_con_vals[i] = coef;
      }

      SparseVectorView<REAL> row_data( &added_con_vals[0], &added_con_inds[0],
                                       (int)conflict_set_candidates.size() );

      Constraint<REAL> conf_con( row_data, row_flag, lhs, 0.0 );
      constraints.push_back( conf_con );

      return;
   }
   void
   simple_cut_from_fixings( Vec<SingleBoundChange<REAL>>& bound_changes,
                            Vec<int>& conflict_set_candidates,
                            Vec<Constraint<REAL>>& constraints )
   {

      // simplest (but also worse) cut is created from the fixings
      // e.g. Fixings x1 = x2 = x3 = 0 -> x1 + x2 + x3 >= 1
      // e.g. Fixings x1 = x2 = 0, x3 = 1 -> x1 + x2 + (1 - x3) >= 1
      RowFlags row_flag;
      row_flag.set( RowFlag::kRhsInf );

      REAL added_con_vals[conflict_set_candidates.size()];
      int added_con_inds[conflict_set_candidates.size()];

      REAL lhs = 1.0;

      for( int i = 0; i < bound_changes.size(); i++ )
      {
         if( bound_changes[i].is_manually_triggered() )
         {
            REAL coef =
                bound_changes[i].get_new_bound_value() > 0.5 ? -1.0 : 1.0;

            if( coef < 0 )
               lhs--;
            added_con_inds[i] = bound_changes[i].get_col();
            added_con_vals[i] = coef;
         }
      }
      SparseVectorView<REAL> row_data( &added_con_vals[0], &added_con_inds[0],
                                       conflict_set_candidates.size() );

      Constraint<REAL> conf_con( row_data, row_flag, lhs, 0.0 );
      constraints.push_back( conf_con );

      return;
   }
};

} // namespace papilo
