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
   void
   perform_conflict_analysis( Vec<SingleBoundChange<REAL>>& bound_changes,
                              Vec<std::pair<int, int>>& infeasible_rows,
                              Vec<Constraint<REAL>>& constraints )
   {
      Vec<int> decision_levels;
      Vec<int> pos_in_bound_changes;
      Vec<int> conflict_set_candidates;
      Vec<bool> was_in_candidates;

      // bound change data as vectors for easier access
      decision_levels.resize( problem.getNCols(), 0 );
      pos_in_bound_changes.resize( problem.getNCols(), -1 );
      was_in_candidates.resize( problem.getNCols(), 0 );

      assert( infeasible_rows.size() > 0 );
      // row that led to infeasibility
      int conflict_row_index = infeasible_rows[0].second;
      // position in stack when infeasible
      int pos_at_infeasibility = infeasible_rows[0].first;

      int last_col = -1;
      int curr_level = 0;

      // ToDo integer/continuous case:
      Vec<ColFlags> col_flags = problem.getColFlags();
      Vec<REAL> lbs = problem.getLowerBounds();
      Vec<REAL> ubs = problem.getUpperBounds();
      for( int i = 0; i < col_flags.size(); i++ )
      {
         if( !col_flags[i].test( ColFlag::kIntegral ) ||
             num.isGT( lbs[i], 0 ) || num.isLT( ubs[i], 1 ) )
            return;
      }
      for( int i = 0; i < pos_at_infeasibility; i++ )
      {
         int col = bound_changes[i].get_col();
         if( last_col == col )
            continue;
         if( bound_changes[i].is_manually_triggered() )
            curr_level++;
         decision_levels[col] = curr_level;
         pos_in_bound_changes[col] = i;
         last_col = col;
      }
      int number_fixings = get_number_fixings( bound_changes );
      // Only fixings (works only for binaries)
      if( number_fixings == bound_changes.size() )
      {
         simple_cut_from_fixings( bound_changes, conflict_set_candidates,
                                  constraints );
         return;
      }
      else
      {
         // Find subset of indices that explain the infeasibility
         // adds column indices in conflict_set_candidates
         explain_infeasibility( bound_changes, pos_in_bound_changes,
                                conflict_set_candidates, was_in_candidates,
                                conflict_row_index );
         // last decision level
         int last_decision_level = get_last_decision_level(
             bound_changes, conflict_set_candidates, decision_levels );

         int num_vars_last_decision_level = get_number_variables_decision_level(
             decision_levels, conflict_set_candidates, last_decision_level );

         if( num_vars_last_decision_level == 1 )
         {
            // return conflict constraint
            add_constraint( bound_changes, pos_in_bound_changes,
                            conflict_set_candidates, constraints );
            msg.info( "Only one variable at last decision level! \n" );
            return;
         }
         // First-FUIP
         // Resolve as long as more than one bound changes at last decision
         // level
         while( num_vars_last_decision_level > 1 )
         {
            int col_index;
            int position_in_conflict_set;

            get_latest_col_index_in_decision_level(
                decision_levels, pos_in_bound_changes, conflict_set_candidates,
                last_decision_level, col_index, position_in_conflict_set );

            int antecedent_row_index =
                bound_changes[pos_in_bound_changes[col_index]].get_reason_row();
            if( antecedent_row_index == -1 )
            {
               msg.info( "There should exist an antecedent row" );
               return;
            }
            else
            {

               // remove latest_col_index
               // was_in_candidates[col_index] = 0;
               conflict_set_candidates.erase( conflict_set_candidates.begin() +
                                              position_in_conflict_set );

               // resolve
               explain_infeasibility(
                   bound_changes, pos_in_bound_changes, conflict_set_candidates,
                   was_in_candidates, antecedent_row_index, col_index );
            }

            num_vars_last_decision_level = get_number_variables_decision_level(
                decision_levels, conflict_set_candidates, last_decision_level );
         }
      }

      add_constraint( bound_changes, pos_in_bound_changes,
                      conflict_set_candidates, constraints );
      return;
   }

   void
   perform_conflict_analysis()
   {
      msg.info( "function call is dummy and waited to be implemented above\n" );
      return;
   }

 private:
   bool
   is_rhs_reason( Vec<SingleBoundChange<REAL>>& bound_changes,
                  Vec<int>& pos_in_bound_changes, int row_idx, int col_idx )
   {
      // get conflict row
      SparseVectorView<REAL> conflict_row =
          problem.getConstraintMatrix().getRowCoefficients( row_idx );
      int row_length = conflict_row.getLength();
      const int* row_inds = conflict_row.getIndices();
      const REAL* row_vals = conflict_row.getValues();

      REAL rhs = problem.getConstraintMatrix().getRightHandSides()[row_idx];

      /* rhs is reason and coeff is positive -> lower bound */
      /* rhs is reason and coeff is negative -> upper bound */
      int pos;
      if( col_idx != -1 )
      {
         for( int i = 0; i < row_length; i++ )
         {
            if( row_inds[i] == col_idx )
            {

               pos = pos_in_bound_changes[row_inds[i]];
               assert( !num.isZero( row_vals[i] ) );
               assert( !( pos == -1 ) );
               return ( ( num.isGT( row_vals[i], 0 ) &&
                          !bound_changes[pos].is_lower_bound() ) ||
                        ( num.isLT( row_vals[i], 0 ) &&
                          bound_changes[pos].is_lower_bound() ) );
            }
         }
         msg.error( " is_rhs_reason should have terminated! \n" );
      }
      else
      {
         StableSum<REAL> min_activity{ 0 };
         min_activity.add( problem.getRowActivities()[row_idx].min );
         for( int i = 0; i < row_length; i++ )
         {
            assert( !num.isZero( row_vals[i] ) );
            pos = pos_in_bound_changes[row_inds[i]];
            if( num.isGT( row_vals[i], 0 ) && pos > 0 )
            {
               if( num.isGT( bound_changes[pos].get_new_bound_value(),
                             problem.getLowerBounds()[row_inds[i]] ) )
               {
                  min_activity.add( row_vals[i] *
                                    ( bound_changes[pos].get_new_bound_value() -
                                      problem.getLowerBounds()[row_inds[i]] ) );
               }
            }
            else if( num.isLT( row_vals[i], 0 ) && pos > 0 )
            {
               if( num.isLT( bound_changes[pos].get_new_bound_value(),
                             problem.getUpperBounds()[row_inds[i]] ) )
               {
                  min_activity.add(
                      -row_vals[i] *
                      ( problem.getUpperBounds()[row_inds[i]] -
                        bound_changes[pos].get_new_bound_value() ) );
               }
            }
            if( num.isGT( min_activity.get(), rhs ) )
               return true;
         }
      }
      return false;
   }

   bool
   is_lhs_reason( Vec<SingleBoundChange<REAL>>& bound_changes,
                  Vec<int>& pos_in_bound_changes, int row_idx, int col_idx )
   {
      // get conflict row
      SparseVectorView<REAL> conflict_row =
          problem.getConstraintMatrix().getRowCoefficients( row_idx );
      int row_length = conflict_row.getLength();
      const int* row_inds = conflict_row.getIndices();
      const REAL* row_vals = conflict_row.getValues();

      REAL lhs = problem.getConstraintMatrix().getLeftHandSides()[row_idx];

      /* lhs is reason and coeff is negative -> lower bound */
      /* lhs is reason and coeff is positive -> upper bound */
      int pos;
      if( col_idx != -1 )
      {
         for( int i = 0; i < row_length; i++ )
         {
            if( row_inds[i] == col_idx )
            {

               pos = pos_in_bound_changes[row_inds[i]];
               assert( !num.isZero( row_vals[i] ) );
               assert( !( pos == -1 ) );
               return ( ( num.isGT( row_vals[i], 0 ) &&
                          bound_changes[pos].is_lower_bound() ) ||
                        ( num.isLT( row_vals[i], 0 ) &&
                          !bound_changes[pos].is_lower_bound() ) );
            }
         }
         msg.error( " is_lhs_reason should have terminated! \n" );
      }
      else
      {
         StableSum<REAL> max_activity{ 0 };
         max_activity.add( problem.getRowActivities()[row_idx].max );
         for( int i = 0; i < row_length; i++ )
         {
            assert( !num.isZero( row_vals[i] ) );
            pos = pos_in_bound_changes[row_inds[i]];
            if( num.isGT( row_vals[i], 0 ) && pos > 0 )
            {
               if( num.isLT( bound_changes[pos].get_new_bound_value(),
                             problem.getUpperBounds()[row_inds[i]] ) )
               {
                  max_activity.add(
                      -row_vals[i] *
                      ( problem.getUpperBounds()[row_inds[i]] -
                        bound_changes[pos].get_new_bound_value() ) );
               }
            }
            else if( num.isLT( row_vals[i], 0 ) && pos > 0 )
            {
               if( num.isGT( bound_changes[pos].get_new_bound_value(),
                             problem.getLowerBounds()[row_inds[i]] ) )
               {
                  max_activity.add( row_vals[i] *
                                    ( bound_changes[pos].get_new_bound_value() -
                                      problem.getLowerBounds()[row_inds[i]] ) );
               }
            }
            if( num.isLT( max_activity.get(), lhs ) )
               return true;
         }
      }
      return false;
   }

   void
   explain_infeasibility( Vec<SingleBoundChange<REAL>>& bound_changes,
                          Vec<int>& pos_in_bound_changes,
                          Vec<int>& conflict_set_candidates,
                          Vec<bool>& was_in_candidates, int row_idx,
                          int col_idx = -1 )
   {
      RowFlags row_flag = problem.getConstraintMatrix().getRowFlags()[row_idx];

      // For GE constraints
      if( row_flag.test( RowFlag::kRhsInf ) )
      {
         explain_infeasibility_ge( bound_changes, pos_in_bound_changes,
                                   conflict_set_candidates, was_in_candidates,
                                   row_idx, col_idx );
      }
      // For LE constraints
      else if( row_flag.test( RowFlag::kLhsInf ) )
      {
         explain_infeasibility_le( bound_changes, pos_in_bound_changes,
                                   conflict_set_candidates, was_in_candidates,
                                   row_idx, col_idx );
      }
      // For equalities or ranged rows
      else
      {
         // rhs is reason and coeff is positive, or lhs is reason and coeff is
         // negative -> lower bound
         // lhs is reason and coeff is positive, or rhs is reason and coeff is
         // negative -> upper bound
         if( is_lhs_reason( bound_changes, pos_in_bound_changes, row_idx,
                            col_idx ) )
         {
            explain_infeasibility_ge( bound_changes, pos_in_bound_changes,
                                      conflict_set_candidates,
                                      was_in_candidates, row_idx, col_idx );
         }
         else if( is_rhs_reason( bound_changes, pos_in_bound_changes, row_idx,
                                 col_idx ) )
         {
            explain_infeasibility_le( bound_changes, pos_in_bound_changes,
                                      conflict_set_candidates,
                                      was_in_candidates, row_idx, col_idx );
         }
         else
         {
            msg.error( "Problem explaining infeasibility of row {} \n",
                       row_idx );
         }
      }

      return;
   }
   void
   explain_infeasibility_ge( Vec<SingleBoundChange<REAL>>& bound_changes,
                             Vec<int>& pos_in_bound_changes,
                             Vec<int>& conflict_set_candidates,
                             Vec<bool>& was_in_candidates, int row_idx,
                             int col_idx )
   {
      // get conflict row
      SparseVectorView<REAL> conflict_row =
          problem.getConstraintMatrix().getRowCoefficients( row_idx );
      int row_length = conflict_row.getLength();
      const int* row_inds = conflict_row.getIndices();
      const REAL* row_vals = conflict_row.getValues();

      REAL lhs = problem.getConstraintMatrix().getLeftHandSides()[row_idx];

      StableSum<REAL> max_activity{ 0 };
      max_activity.add( problem.getRowActivities()[row_idx].max );
      for( int i = 0; i < row_length; i++ )
      {
         if( col_idx == row_inds[i] )
            continue;

         assert( !num.isZero( row_vals[i] ) );
         int pos = pos_in_bound_changes[row_inds[i]];
         assert( !num.isZero( row_vals[i] ) );
         if( num.isGT( row_vals[i], 0 ) && pos > 0 )
         {
            if( num.isLT( bound_changes[pos].get_new_bound_value(),
                          problem.getUpperBounds()[row_inds[i]] ) )
            {
               max_activity.add( -row_vals[i] *
                                 ( problem.getUpperBounds()[row_inds[i]] -
                                   bound_changes[pos].get_new_bound_value() ) );
               if( !was_in_candidates[row_inds[i]] )
               {
                  was_in_candidates[row_inds[i]] = 1;
                  conflict_set_candidates.push_back( row_inds[i] );
               }
            }
         }
         else if( num.isLT( row_vals[i], 0 ) && pos > 0 )
         {
            if( num.isGT( bound_changes[pos].get_new_bound_value(),
                          problem.getLowerBounds()[row_inds[i]] ) )
            {
               max_activity.add( row_vals[i] *
                                 ( bound_changes[pos].get_new_bound_value() -
                                   problem.getLowerBounds()[row_inds[i]] ) );
               if( !was_in_candidates[row_inds[i]] )
               {
                  was_in_candidates[row_inds[i]] = 1;
                  conflict_set_candidates.push_back( row_inds[i] );
               }
            }
         }
         if( num.isLT( max_activity.get(), lhs ) )
            break;
      }
      return;
   }

   void
   explain_infeasibility_le( Vec<SingleBoundChange<REAL>>& bound_changes,
                             Vec<int>& pos_in_bound_changes,
                             Vec<int>& conflict_set_candidates,
                             Vec<bool>& was_in_candidates, int row_idx,
                             int col_idx )
   {
      // get conflict row
      SparseVectorView<REAL> conflict_row =
          problem.getConstraintMatrix().getRowCoefficients( row_idx );
      int row_length = conflict_row.getLength();
      const int* row_inds = conflict_row.getIndices();
      const REAL* row_vals = conflict_row.getValues();

      REAL rhs = problem.getConstraintMatrix().getRightHandSides()[row_idx];

      StableSum<REAL> min_activity{ 0 };
      min_activity.add( problem.getRowActivities()[row_idx].min );
      for( int i = 0; i < row_length; i++ )
      {
         if( col_idx == row_inds[i] )
            continue;
         assert( !num.isZero( row_vals[i] ) );
         int pos = pos_in_bound_changes[row_inds[i]];
         if( num.isGT( row_vals[i], 0 ) && pos > 0 )
         {
            if( num.isGT( bound_changes[pos].get_new_bound_value(),
                          problem.getLowerBounds()[row_inds[i]] ) )
            {
               min_activity.add( row_vals[i] *
                                 ( bound_changes[pos].get_new_bound_value() -
                                   problem.getLowerBounds()[row_inds[i]] ) );
               if( !was_in_candidates[row_inds[i]] )
               {
                  was_in_candidates[row_inds[i]] = 1;
                  conflict_set_candidates.push_back( row_inds[i] );
               }
            }
         }
         else if( num.isLT( row_vals[i], 0 ) && pos > 0 )
         {
            if( num.isLT( bound_changes[pos].get_new_bound_value(),
                          problem.getUpperBounds()[row_inds[i]] ) )
            {
               min_activity.add( -row_vals[i] *
                                 ( problem.getUpperBounds()[row_inds[i]] -
                                   bound_changes[pos].get_new_bound_value() ) );
               if( !was_in_candidates[row_inds[i]] )
               {
                  was_in_candidates[row_inds[i]] = 1;
                  conflict_set_candidates.push_back( row_inds[i] );
               }
            }
         }

         if( num.isGT( min_activity.get(), rhs ) )
            break;
      }
      return;
   }

   int
   get_last_decision_level( Vec<SingleBoundChange<REAL>>& bound_changes,
                            Vec<int>& conflict_set_candidates,
                            Vec<int>& decision_levels )
   {
      int level = 0;
      for( int i = 0; i < conflict_set_candidates.size(); i++ )
      {
         level = decision_levels[conflict_set_candidates[i]] > level
                     ? decision_levels[conflict_set_candidates[i]]
                     : level;
      }
      return level;
   }

   int
   get_number_variables_decision_level( Vec<int>& decision_levels,
                                        Vec<int>& conflict_set_candidates,
                                        int decision_level )
   {
      int number_variables_decision_level = 0;
      for( int i = 0; i < conflict_set_candidates.size(); i++ )
      {
         if( decision_levels[conflict_set_candidates[i]] == decision_level )
            number_variables_decision_level++;
      }
      return number_variables_decision_level;
   }
   void
   get_latest_col_index_in_decision_level( Vec<int>& decision_levels,
                                           Vec<int>& pos_in_bound_changes,
                                           Vec<int>& conflict_set_candidates,
                                           int decision_level, int& col,
                                           int& position_in_conflict_set )
   {

      int latest_position = -1;
      for( int i = 0; i < conflict_set_candidates.size(); i++ )
      {
         if( ( decision_levels[conflict_set_candidates[i]] ==
               decision_level ) &&
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
      int curr_decision_level = 0;
      for( int i = 0; i < bound_changes.size(); i++ )
      {
         if( bound_changes[i].is_manually_triggered() &&
             ( curr_decision_level > bound_changes[i].get_depth_level() ) )
            number_fixings++;
         curr_decision_level = bound_changes[i].get_depth_level();
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

      REAL* vals = new REAL[conflict_set_candidates.size()];
      int* inds = new int[conflict_set_candidates.size()];

      REAL lhs = 1.0;
      for( int i = 0; i < conflict_set_candidates.size(); i++ )
      {
         int col = conflict_set_candidates[i];
         int pos = pos_in_bound_changes[col];
         REAL coef =
             bound_changes[pos].get_new_bound_value() > 0.5 ? -1.0 : 1.0;
         if( coef < 0 )
            lhs--;
         inds[i] = col;
         vals[i] = coef;
      }

      SparseVectorView<REAL> row_data( &vals[0], &inds[0],
                                       (int)conflict_set_candidates.size() );

      Constraint<REAL> conf_con( row_data, row_flag, lhs, 0.0 );
      constraints.push_back( conf_con );
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

      REAL* vals = new REAL[conflict_set_candidates.size()];
      int* inds = new int[conflict_set_candidates.size()];

      REAL lhs = 1.0;

      for( int i = 0; i < bound_changes.size(); i++ )
      {
         assert( bound_changes[i].is_manually_triggered() );
         REAL coef = bound_changes[i].get_new_bound_value() > 0.5 ? -1.0 : 1.0;
         if( coef < 0 )
            lhs--;
         inds[i] = bound_changes[i].get_col();
         vals[i] = coef;
      }
      SparseVectorView<REAL> row_data( &vals[0], &inds[0],
                                       conflict_set_candidates.size() );

      Constraint<REAL> conf_con( row_data, row_flag, lhs, 0.0 );
      constraints.push_back( conf_con );

      return;
   }
};

} // namespace papilo
