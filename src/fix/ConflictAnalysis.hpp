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
#include <map>

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
                              Vec<Constraint<REAL>>& constraints,
                              REAL conflict_size_factor = 0.15 )
   {
      // col: (pos, is_lower)
      std::map<int, std::pair<int, bool>> current_conflict_set;
      // col: (pos_lower, pos_upper)
      // set pos_* to -1 if no bound change
      std::map<int, std::pair<int, int>> pos_in_bound_changes;
      // col: (decision_level_lower, decision_level_upper)
      // set decision_level_* to -1 if no bound change
      std::map<int, std::pair<int, int>> decision_levels;
      // if a general integer appears in the conflict we return no conflicts
      bool general_integers_in_conflict_set = false;

      assert( infeasible_rows.size() > 0 );
      // position in stack when infeasible
      int pos_at_infeasibility = infeasible_rows[0].first;
      int last_col = -1;
      int curr_level = 0;
      for( int i = 0; i < pos_at_infeasibility; i++ )
      {
         int col = bound_changes[i].get_col();
         // if no entry in positions
         if( pos_in_bound_changes.find( col ) == pos_in_bound_changes.end() )
         {
            pos_in_bound_changes.insert( { col, { -1, -1 } } );
            decision_levels.insert( { col, { -1, -1 } } );
         }
         if( bound_changes[i].is_manually_triggered() )
         {
            if( last_col != col )
               curr_level++;
         }
         if( bound_changes[i].is_lower_bound() )
         {
            pos_in_bound_changes[col].first = i;
            decision_levels[col].first = curr_level;
         }
         else if( !bound_changes[i].is_lower_bound() )
         {
            pos_in_bound_changes[col].second = i;
            decision_levels[col].second = curr_level;
         }
         else
            assert( false );
         last_col = col;
      }

      int index_infes = -1;
      int size_smallest_conflict_set = -1;
      // Find infeasible row with smallest initial conflict set
      for( int i = 0; i < infeasible_rows.size(); i++ )
      {
         int size_conflict_set;
         bool general_integers_found = false;
         if( infeasible_rows[i].first == pos_at_infeasibility )
         {
            size_conflict_set = conflict_set_size_of_row(
                bound_changes, pos_in_bound_changes, infeasible_rows[i].second,
                -1, general_integers_found );
            if( general_integers_found )
               continue;
            if( size_smallest_conflict_set == -1 )
            {
               size_smallest_conflict_set = size_conflict_set;
               index_infes = i;
            }
            else if( size_smallest_conflict_set > size_conflict_set )
            {
               size_smallest_conflict_set = size_conflict_set;
               index_infes = i;
            }
         }
         else
            break;
      }
      msg.detailed( "Index {}, Row {}, conflict set size {} \n", index_infes,
                    infeasible_rows[index_infes].second,
                    size_smallest_conflict_set );
      if( index_infes == -1 || size_smallest_conflict_set == -1 ||
          size_smallest_conflict_set == 0 )
      {
         msg.detailed(
             "\t\tConflict analysis returns 0 conflict constraints \n" );
         return;
      }
      // row that led to infeasibility
      int conflict_row_index = infeasible_rows[index_infes].second;
      assert( infeasible_rows[index_infes].first == pos_at_infeasibility );
      // Find subset of indices that explain the infeasibility
      // adds column indices in conflict_set_candidates
      // col_index is -1 since initially we do not resolve
      explain_infeasibility( bound_changes, pos_in_bound_changes,
                             current_conflict_set, conflict_row_index, -1,
                             general_integers_in_conflict_set );
      if( general_integers_in_conflict_set )
      {
         msg.detailed(
             "\t\tConflict analysis returns 0 conflict constraints \n" );
         return;
      }
      // last decision level
      int last_decision_level = get_last_decision_level(
          bound_changes, current_conflict_set, decision_levels );

      int num_vars_last_decision_level = get_number_variables_decision_level(
          decision_levels, current_conflict_set, last_decision_level );

      // col: (pos, is_lower)
      std::map<int, std::pair<int, bool>> all_fuips;
      // First-FUIP
      // Resolve as long as more than one bound changes at last decision level
      int col_index = -1;
      bool one_fuip_conflict = true;
      bool resolved_bounds = true;
      // maximal number of nonzeros in conflict constraint
      int max_size_conflict = conflict_size_factor * problem.getNumIntegralCols();
      while( last_decision_level > 0 )
      {
         while( num_vars_last_decision_level > 1 )
         {
            get_latest_col_index_in_decision_level(
                decision_levels, pos_in_bound_changes, current_conflict_set,
                last_decision_level, col_index );

            int antecedent_row_index =
                bound_changes[current_conflict_set.at( col_index ).first]
                    .get_reason_row();
            assert( antecedent_row_index != -1 );
            // resolve
            explain_infeasibility( bound_changes, pos_in_bound_changes,
                                   current_conflict_set, antecedent_row_index,
                                   col_index,
                                   general_integers_in_conflict_set );
            resolved_bounds = true;
            if( general_integers_in_conflict_set )
            {
               msg.detailed( "\t\tConflict analysis returns 0 conflict "
                             "constraints \n" );
               return;
            }
            num_vars_last_decision_level = get_number_variables_decision_level(
                decision_levels, current_conflict_set, last_decision_level );
         }
         // 1-FUIP conflict is always added!
         if( one_fuip_conflict )
         {
            if( current_conflict_set.size() <= max_size_conflict )
            {
               // add constraint and continue
               add_constraint( bound_changes, pos_in_bound_changes,
                               current_conflict_set, all_fuips, constraints );
               max_size_conflict = num.min( 2 * current_conflict_set.size(),
                                            max_size_conflict );
            }
            resolved_bounds = false;
            one_fuip_conflict = false;
            get_latest_col_index_in_decision_level(
                decision_levels, pos_in_bound_changes, current_conflict_set,
                last_decision_level, col_index );

            int pos_lower = pos_in_bound_changes.at( col_index ).first;
            int pos_upper = pos_in_bound_changes.at( col_index ).second;
            assert( !( pos_lower == -1 && pos_upper == -1 ) );
            if( pos_lower != -1 )
            {
               pos_in_bound_changes.at( col_index ).first = -1;
               all_fuips.insert( { col_index, { pos_lower, true } } );
            }
            else
            {
               pos_in_bound_changes.at( col_index ).second = -1;
               all_fuips.insert( { col_index, { pos_upper, false } } );
            }
            // remove fuip of last_decision_level from conflict set
            current_conflict_set.erase( col_index );
         }
         // All-FUIP
         else
         {
            get_latest_col_index_in_decision_level(
                decision_levels, pos_in_bound_changes, current_conflict_set,
                last_decision_level, col_index );
            int antecedent_row_index =
                bound_changes[current_conflict_set.at( col_index ).first]
                    .get_reason_row();
            // In case FUIP is a fixing
            if( antecedent_row_index == -1 )
            {
               int pos_lower = pos_in_bound_changes.at( col_index ).first;
               int pos_upper = pos_in_bound_changes.at( col_index ).second;
               assert( !( pos_lower == -1 && pos_upper == -1 ) );
               if( pos_lower != -1 )
               {
                  pos_in_bound_changes.at( col_index ).first = -1;
                  all_fuips.insert( { col_index, { pos_lower, true } } );
               }
               else
               {
                  pos_in_bound_changes.at( col_index ).second = -1;
                  all_fuips.insert( { col_index, { pos_upper, false } } );
               }
               // remove fuip of last_decision_level from conflict set
               current_conflict_set.erase( col_index );
            }
            else
            {
               int pos_lower = pos_in_bound_changes.at( col_index ).first;
               int pos_upper = pos_in_bound_changes.at( col_index ).second;
               assert( !( pos_lower == -1 && pos_upper == -1 ) );
               if( pos_lower != -1 )
               {
                  pos_in_bound_changes.at( col_index ).first = -1;
                  all_fuips.insert( { col_index, { pos_lower, true } } );
               }
               else
               {
                  pos_in_bound_changes.at( col_index ).second = -1;
                  all_fuips.insert( { col_index, { pos_upper, false } } );
               }
               // remove fuip of last_decision_level from conflict set
               current_conflict_set.erase( col_index );
            }
         }
         if( resolved_bounds &&
             current_conflict_set.size() + all_fuips.size() <=
                 max_size_conflict )
            add_constraint( bound_changes, pos_in_bound_changes,
                            current_conflict_set, all_fuips, constraints );
         resolved_bounds = false;
         // compute last decision level
         last_decision_level = get_last_decision_level(
             bound_changes, current_conflict_set, decision_levels );
         if( last_decision_level == 0 )
            break;
         // if not zero compute the number of elements in the last
         // decision level
         num_vars_last_decision_level = get_number_variables_decision_level(
             decision_levels, current_conflict_set, last_decision_level );
      }
      msg.detailed(
          "\t\tConflict analysis added {} cons (All-FUIP) of length: ",
          constraints.size() );
      for( int i = 0; i < constraints.size(); i++ )
      {
         msg.detailed( "{} ", constraints[i].get_data().getLength() );
      }
      msg.detailed( "\n" );
      return;
   }

 private:
   bool
   is_rhs_reason(
       const Vec<SingleBoundChange<REAL>>& bound_changes,
       const std::map<int, std::pair<int, int>>& pos_in_bound_changes,
       int row_idx, int col_idx, bool& general_integers_in_conflict_set )
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
      int pos_lower;
      int pos_upper;
      if( col_idx != -1 )
      {
         for( int i = 0; i < row_length; i++ )
         {
            if( row_inds[i] == col_idx )
            {
               assert( !( pos_in_bound_changes.find( col_idx ) ==
                          pos_in_bound_changes.end() ) );
               pos_lower = pos_in_bound_changes.at( col_idx ).first;
               pos_upper = pos_in_bound_changes.at( col_idx ).second;
               assert( !num.isZero( row_vals[i] ) );
               assert( !( pos_lower == -1 && pos_upper == -1 ) );
               return ( ( num.isGT( row_vals[i], 0 ) && pos_upper != -1 ) ||
                        ( num.isLT( row_vals[i], 0 ) && pos_lower != -1 ) );
            }
         }
         msg.error( "\t\tis_rhs_reason should have terminated! \n" );
      }
      else
      {
         StableSum<REAL> min_activity{ 0 };
         min_activity.add( problem.getRowActivities()[row_idx].min );
         for( int i = 0; i < row_length; i++ )
         {
            assert( !num.isZero( row_vals[i] ) );
            // check if there exists a bound change for row_inds[i]
            if( ( pos_in_bound_changes.find( row_inds[i] ) ==
                  pos_in_bound_changes.end() ) )
               continue;
            pos_lower = pos_in_bound_changes.at( row_inds[i] ).first;
            pos_upper = pos_in_bound_changes.at( row_inds[i] ).second;
            if( num.isGT( row_vals[i], 0 ) && pos_lower >= 0 )
            {
               assert( bound_changes[pos_lower].is_lower_bound() );
               if( !is_binary( row_inds[i] ) )
               {
                  general_integers_in_conflict_set = true;
                  return false;
               }
               assert( num.isGT( bound_changes[pos_lower].get_new_bound_value(),
                                 problem.getLowerBounds()[row_inds[i]] ) );
               min_activity.add(
                   row_vals[i] *
                   ( bound_changes[pos_lower].get_new_bound_value() -
                     problem.getLowerBounds()[row_inds[i]] ) );
            }
            else if( num.isLT( row_vals[i], 0 ) && pos_upper >= 0 )
            {
               assert( !bound_changes[pos_upper].is_lower_bound() );

               if( !is_binary( row_inds[i] ) )
               {
                  general_integers_in_conflict_set = true;
                  return false;
               }
               assert( num.isLT( bound_changes[pos_upper].get_new_bound_value(),
                                 problem.getUpperBounds()[row_inds[i]] ) );

               min_activity.add(
                   -row_vals[i] *
                   ( problem.getUpperBounds()[row_inds[i]] -
                     bound_changes[pos_upper].get_new_bound_value() ) );
            }
            if( num.isGT( min_activity.get(), rhs ) )
               return true;
         }
      }
      return false;
   }

   bool
   is_lhs_reason(
       const Vec<SingleBoundChange<REAL>>& bound_changes,
       const std::map<int, std::pair<int, int>>& pos_in_bound_changes,
       int row_idx, int col_idx, bool& general_integers_in_conflict_set )
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
      int pos_lower;
      int pos_upper;
      if( col_idx != -1 )
      {
         for( int i = 0; i < row_length; i++ )
         {
            if( row_inds[i] == col_idx )
            {
               assert( !( pos_in_bound_changes.find( col_idx ) ==
                          pos_in_bound_changes.end() ) );
               pos_lower = pos_in_bound_changes.at( col_idx ).first;
               pos_upper = pos_in_bound_changes.at( col_idx ).second;
               assert( !num.isZero( row_vals[i] ) );
               assert( !( pos_lower == -1 && pos_upper == -1 ) );
               return ( ( num.isGT( row_vals[i], 0 ) && ( pos_lower != -1 ) ) ||
                        ( num.isLT( row_vals[i], 0 ) && ( pos_upper != -1 ) ) );
            }
         }
         msg.error( "\t\tis_lhs_reason should have terminated! \n" );
      }
      else
      {
         StableSum<REAL> max_activity{ 0 };
         max_activity.add( problem.getRowActivities()[row_idx].max );
         for( int i = 0; i < row_length; i++ )
         {
            assert( !num.isZero( row_vals[i] ) );
            // check if there exists a bound change for row_inds[i]
            if( ( pos_in_bound_changes.find( row_inds[i] ) ==
                  pos_in_bound_changes.end() ) )
               continue;
            pos_lower = pos_in_bound_changes.at( row_inds[i] ).first;
            pos_upper = pos_in_bound_changes.at( row_inds[i] ).second;
            if( num.isGT( row_vals[i], 0 ) && pos_upper >= 0 )
            {
               assert( !bound_changes[pos_upper].is_lower_bound() );
               if( !is_binary( row_inds[i] ) )
               {
                  general_integers_in_conflict_set = true;
                  return false;
               }
               assert( num.isLT( bound_changes[pos_upper].get_new_bound_value(),
                                 problem.getUpperBounds()[row_inds[i]] ) );
               max_activity.add(
                   -row_vals[i] *
                   ( problem.getUpperBounds()[row_inds[i]] -
                     bound_changes[pos_upper].get_new_bound_value() ) );
            }
            else if( num.isLT( row_vals[i], 0 ) && pos_lower >= 0 )
            {
               assert( bound_changes[pos_lower].is_lower_bound() );
               if( !is_binary( row_inds[i] ) )
               {
                  general_integers_in_conflict_set = true;
                  return false;
               }
               assert( num.isGT( bound_changes[pos_lower].get_new_bound_value(),
                                 problem.getLowerBounds()[row_inds[i]] ) );
               max_activity.add(
                   row_vals[i] *
                   ( bound_changes[pos_lower].get_new_bound_value() -
                     problem.getLowerBounds()[row_inds[i]] ) );
            }
            if( num.isLT( max_activity.get(), lhs ) )
               return true;
         }
      }
      return false;
   }

   void
   explain_infeasibility(
       Vec<SingleBoundChange<REAL>>& bound_changes,
       std::map<int, std::pair<int, int>>& pos_in_bound_changes,
       std::map<int, std::pair<int, bool>>& current_conflict_set, int row_idx,
       int col_idx, bool& general_integers_in_conflict_set )
   {
      RowFlags row_flag = problem.getConstraintMatrix().getRowFlags()[row_idx];
      int pos_col = bound_changes.size();
      // For GE constraints
      if( row_flag.test( RowFlag::kRhsInf ) )
      {
         if( col_idx != -1 )
         {
            // the bound change will no longer be consider since it is
            // resolved
            assert( current_conflict_set.find( col_idx ) !=
                    current_conflict_set.end() );
            if( current_conflict_set.at( col_idx ).second == true )
            {
               pos_col = pos_in_bound_changes.at( col_idx ).first;
               pos_in_bound_changes.at( col_idx ).first = -1;
            }
            else
            {
               pos_col = pos_in_bound_changes.at( col_idx ).second;
               pos_in_bound_changes.at( col_idx ).second = -1;
            }
            current_conflict_set.erase( col_idx );
         }
         explain_infeasibility_ge( bound_changes, pos_in_bound_changes,
                                   current_conflict_set, row_idx, col_idx,
                                   general_integers_in_conflict_set, pos_col );
      }
      // For LE constraints
      else if( row_flag.test( RowFlag::kLhsInf ) )
      {
         if( col_idx != -1 )
         {
            // the bound change will no longer be consider since it is
            // resolved
            assert( current_conflict_set.find( col_idx ) !=
                    current_conflict_set.end() );
            if( current_conflict_set.at( col_idx ).second == true )
            {
               pos_col = pos_in_bound_changes.at( col_idx ).first;
               pos_in_bound_changes.at( col_idx ).first = -1;
            }
            else
            {
               pos_col = pos_in_bound_changes.at( col_idx ).second;
               pos_in_bound_changes.at( col_idx ).second = -1;
            }

            current_conflict_set.erase( col_idx );
         }
         explain_infeasibility_le( bound_changes, pos_in_bound_changes,
                                   current_conflict_set, row_idx, col_idx,
                                   general_integers_in_conflict_set, pos_col );
      }
      // For equalities or ranged rows
      else
      {
         // rhs is reason and coeff is positive, or lhs is reason and coeff is
         // negative -> lower bound
         // lhs is reason and coeff is positive, or rhs is reason and coeff is
         // negative -> upper bound
         if( is_lhs_reason( bound_changes, pos_in_bound_changes, row_idx,
                            col_idx, general_integers_in_conflict_set ) )
         {
            if( col_idx != -1 )
            {
               // the bound change will no longer be consider since it is
               // resolved
               assert( current_conflict_set.find( col_idx ) !=
                       current_conflict_set.end() );
               if( current_conflict_set.at( col_idx ).second == true )
               {
                  pos_col = pos_in_bound_changes.at( col_idx ).first;
                  pos_in_bound_changes.at( col_idx ).first = -1;
               }
               else
               {
                  pos_col = pos_in_bound_changes.at( col_idx ).second;
                  pos_in_bound_changes.at( col_idx ).second = -1;
               }

               current_conflict_set.erase( col_idx );
            }
            explain_infeasibility_ge(
                bound_changes, pos_in_bound_changes, current_conflict_set,
                row_idx, col_idx, general_integers_in_conflict_set, pos_col );
         }
         else if( is_rhs_reason( bound_changes, pos_in_bound_changes, row_idx,
                                 col_idx, general_integers_in_conflict_set ) )
         {
            if( col_idx != -1 )
            {
               // the bound change will no longer be consider since it is
               // resolved
               assert( current_conflict_set.find( col_idx ) !=
                       current_conflict_set.end() );
               if( current_conflict_set.at( col_idx ).second == true )
               {
                  pos_col = pos_in_bound_changes.at( col_idx ).first;
                  pos_in_bound_changes.at( col_idx ).first = -1;
               }
               else
               {
                  pos_col = pos_in_bound_changes.at( col_idx ).second;
                  pos_in_bound_changes.at( col_idx ).second = -1;
               }
               current_conflict_set.erase( col_idx );
            }
            explain_infeasibility_le(
                bound_changes, pos_in_bound_changes, current_conflict_set,
                row_idx, col_idx, general_integers_in_conflict_set, pos_col );
         }
         else
         {
            msg.error( "\t\tProblem explaining infeasibility of row {} \n",
                       row_idx );
         }
      }

      return;
   }
   int
   conflict_set_size_of_row(
       const Vec<SingleBoundChange<REAL>>& bound_changes,
       const std::map<int, std::pair<int, int>>& pos_in_bound_changes,
       int row_idx, int col_idx, bool& general_integers_in_conflict_set )
   {
      std::map<int, std::pair<int, bool>> current_conflict_set;
      RowFlags row_flag = problem.getConstraintMatrix().getRowFlags()[row_idx];
      int pos_col = bound_changes.size();
      // For GE constraints
      if( row_flag.test( RowFlag::kRhsInf ) )
      {
         explain_infeasibility_ge( bound_changes, pos_in_bound_changes,
                                   current_conflict_set, row_idx, col_idx,
                                   general_integers_in_conflict_set, pos_col );
      }
      // For LE constraints
      else if( row_flag.test( RowFlag::kLhsInf ) )
      {
         explain_infeasibility_le( bound_changes, pos_in_bound_changes,
                                   current_conflict_set, row_idx, col_idx,
                                   general_integers_in_conflict_set, pos_col );
      }
      // For equalities or ranged rows
      else
      {
         // rhs is reason and coeff is positive, or lhs is reason and coeff is
         // negative -> lower bound
         // lhs is reason and coeff is positive, or rhs is reason and coeff is
         // negative -> upper bound
         if( is_lhs_reason( bound_changes, pos_in_bound_changes, row_idx,
                            col_idx, general_integers_in_conflict_set ) )
         {
            explain_infeasibility_ge(
                bound_changes, pos_in_bound_changes, current_conflict_set,
                row_idx, col_idx, general_integers_in_conflict_set, pos_col );
         }
         else if( is_rhs_reason( bound_changes, pos_in_bound_changes, row_idx,
                                 col_idx, general_integers_in_conflict_set ) )
         {
            explain_infeasibility_le(
                bound_changes, pos_in_bound_changes, current_conflict_set,
                row_idx, col_idx, general_integers_in_conflict_set, pos_col );
         }
         else
         {
            msg.error( "\t\tProblem explaining infeasibility of row {} \n",
                       row_idx );
         }
      }

      return current_conflict_set.size();
   }

   void
   explain_infeasibility_ge(
       const Vec<SingleBoundChange<REAL>>& bound_changes,
       const std::map<int, std::pair<int, int>>& pos_in_bound_changes,
       std::map<int, std::pair<int, bool>>& current_conflict_set, int row_idx,
       int col_idx, bool& general_integers_in_conflict_set, int pos_col )
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
         // ToDo: Maybe use also the activity change from the current column to
         // get a shorter clause
         if( col_idx == row_inds[i] )
            continue;
         assert( !num.isZero( row_vals[i] ) );
         // check if there exists a bound change for row_inds[i]
         if( ( pos_in_bound_changes.find( row_inds[i] ) ==
               pos_in_bound_changes.end() ) )
            continue;
         int pos_lower = pos_in_bound_changes.at( row_inds[i] ).first;
         int pos_upper = pos_in_bound_changes.at( row_inds[i] ).second;
         assert( !num.isZero( row_vals[i] ) );
         if( num.isGT( row_vals[i], 0 ) && pos_upper >= 0 &&
             pos_upper < pos_col )
         {
            assert( !bound_changes[pos_upper].is_lower_bound() );

            if( !is_binary( row_inds[i] ) )
            {
               general_integers_in_conflict_set = true;
               return;
            }
            assert( num.isLT( bound_changes[pos_upper].get_new_bound_value(),
                              problem.getUpperBounds()[row_inds[i]] ) );
            max_activity.add(
                -row_vals[i] *
                ( problem.getUpperBounds()[row_inds[i]] -
                  bound_changes[pos_upper].get_new_bound_value() ) );
            current_conflict_set.insert(
                { row_inds[i], { pos_upper, false } } );
         }
         else if( num.isLT( row_vals[i], 0 ) && pos_lower >= 0 &&
                  pos_lower < pos_col )
         {
            assert( bound_changes[pos_lower].is_lower_bound() );
            if( !is_binary( row_inds[i] ) )
            {
               general_integers_in_conflict_set = true;
               return;
            }
            assert( num.isGT( bound_changes[pos_lower].get_new_bound_value(),
                              problem.getLowerBounds()[row_inds[i]] ) );
            max_activity.add( row_vals[i] *
                              ( bound_changes[pos_lower].get_new_bound_value() -
                                problem.getLowerBounds()[row_inds[i]] ) );
            current_conflict_set.insert( { row_inds[i], { pos_lower, true } } );
         }
         if( num.isLT( max_activity.get(), lhs ) )
            break;
      }
      return;
   }

   void
   explain_infeasibility_le(
       const Vec<SingleBoundChange<REAL>>& bound_changes,
       const std::map<int, std::pair<int, int>>& pos_in_bound_changes,
       std::map<int, std::pair<int, bool>>& current_conflict_set, int row_idx,
       int col_idx, bool& general_integers_in_conflict_set, int pos_col )
   { // get conflict row
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
         // ToDo: Maybe use also the activity change from the current column to
         // get a shorter clause
         if( col_idx == row_inds[i] )
            continue;
         assert( !num.isZero( row_vals[i] ) );
         // check if there exists a bound change for row_inds[i]
         if( ( pos_in_bound_changes.find( row_inds[i] ) ==
               pos_in_bound_changes.end() ) )
            continue;
         int pos_lower = pos_in_bound_changes.at( row_inds[i] ).first;
         int pos_upper = pos_in_bound_changes.at( row_inds[i] ).second;
         if( num.isGT( row_vals[i], 0 ) && pos_lower >= 0 &&
             pos_lower < pos_col )
         {
            assert( bound_changes[pos_lower].is_lower_bound() );
            if( !is_binary( row_inds[i] ) )
            {
               general_integers_in_conflict_set = true;
               return;
            }
            assert( num.isGT( bound_changes[pos_lower].get_new_bound_value(),
                              problem.getLowerBounds()[row_inds[i]] ) );

            min_activity.add( row_vals[i] *
                              ( bound_changes[pos_lower].get_new_bound_value() -
                                problem.getLowerBounds()[row_inds[i]] ) );
            current_conflict_set.insert( { row_inds[i], { pos_lower, true } } );
         }
         else if( num.isLT( row_vals[i], 0 ) && pos_upper >= 0 &&
                  pos_upper < pos_col )
         {
            assert( !bound_changes[pos_upper].is_lower_bound() );
            if( !is_binary( row_inds[i] ) )
            {
               general_integers_in_conflict_set = true;
               return;
            }
            assert( num.isLT( bound_changes[pos_upper].get_new_bound_value(),
                              problem.getUpperBounds()[row_inds[i]] ) );
            min_activity.add(
                -row_vals[i] *
                ( problem.getUpperBounds()[row_inds[i]] -
                  bound_changes[pos_upper].get_new_bound_value() ) );
            current_conflict_set.insert(
                { row_inds[i], { pos_upper, false } } );
         }

         if( num.isGT( min_activity.get(), rhs ) )
            break;
      }
      return;
   }

   int
   get_last_decision_level(
       Vec<SingleBoundChange<REAL>>& bound_changes,
       std::map<int, std::pair<int, bool>>& current_conflict_set,
       std::map<int, std::pair<int, int>>& decision_levels )
   {
      if( current_conflict_set.size() == 0 )
         return 0;
      int level = 0;
      int col;
      int pos;
      bool is_lower;
      int col_level;
      for( const std::pair<int, std::pair<int, bool>>& conflict :
           current_conflict_set )
      {
         col = conflict.first;
         pos = conflict.second.first;
         is_lower = conflict.second.second;
         assert( pos != -1 );
         col_level = is_lower ? decision_levels.at( col ).first
                              : decision_levels.at( col ).second;
         level = col_level > level ? col_level : level;
      }
      return level;
   }

   int
   get_number_variables_decision_level(
       std::map<int, std::pair<int, int>>& decision_levels,
       std::map<int, std::pair<int, bool>>& current_conflict_set,
       int decision_level )
   {
      if( current_conflict_set.size() == 0 )
         return 0;
      int number_variables_decision_level = 0;
      int col;
      for( const std::pair<int, std::pair<int, bool>>& conflict :
           current_conflict_set )
      {
         // ToDo for general case
         col = conflict.first;
         assert( decision_levels.at( col ).first != -1 ||
                 decision_levels.at( col ).second != -1 );
         if( decision_levels.at( col ).first == decision_level ||
             decision_levels.at( col ).second == decision_level )
            number_variables_decision_level++;
      }
      return number_variables_decision_level;
   }
   void
   get_latest_col_index_in_decision_level(
       std::map<int, std::pair<int, int>>& decision_levels,
       std::map<int, std::pair<int, int>>& pos_in_bound_changes,
       std::map<int, std::pair<int, bool>>& current_conflict_set,
       int decision_level, int& col )
   {
      int latest_position = -1;
      int curr_col;
      int pos;
      // loop over conflict candidates and check decision level and position
      for( const std::pair<int, std::pair<int, bool>>& conflict :
           current_conflict_set )
      {
         curr_col = conflict.first;
         pos = conflict.second.first;
         assert( pos != -1 );
         if( ( decision_levels.at( curr_col ).first == decision_level ) &&
             ( pos > latest_position ) )
         {
            col = curr_col;
            latest_position = pos;
         }
         if( ( decision_levels.at( curr_col ).second == decision_level ) &&
             ( pos > latest_position ) )
         {
            col = curr_col;
            latest_position = pos;
         }
      }
      return;
   }
   bool
   is_binary( int col )
   {
      if( problem.getColFlags()[col].test( ColFlag::kIntegral ) &&
          !problem.getColFlags()[col].test( ColFlag::kLbInf ) &&
          !problem.getColFlags()[col].test( ColFlag::kUbInf ) )
      {
         if( problem.getUpperBounds()[col] == 1 &&
             problem.getLowerBounds()[col] == 0 )
            return true;
         if( problem.getUpperBounds()[col] == 0 &&
             problem.getLowerBounds()[col] == -1 )
            return true;
      }

      return false;
   }
   // adds constraint based on current conflict set and fuips
   void
   add_constraint( Vec<SingleBoundChange<REAL>>& bound_changes,
                   std::map<int, std::pair<int, int>>& pos_in_bound_changes,
                   std::map<int, std::pair<int, bool>>& current_conflict_set,
                   std::map<int, std::pair<int, bool>>& all_fuips,
                   Vec<Constraint<REAL>>& constraints )
   {
      RowFlags row_flag;
      row_flag.set( RowFlag::kRhsInf );

      REAL* vals = new REAL[current_conflict_set.size() + all_fuips.size()];
      int* inds = new int[current_conflict_set.size() + all_fuips.size()];

      REAL lhs = 1.0;
      int i = 0;
      for( const std::pair<int, std::pair<int, bool>>& conflict :
           current_conflict_set )
      {
         int col = conflict.first;
         int pos = conflict.second.first;
         bool is_lower = conflict.second.second;
         REAL coef;
         assert( is_binary( col ) );
         coef = bound_changes[pos].get_new_bound_value() > 0.5 ? -1.0 : 1.0;
         if( coef < 0 )
            lhs--;
         inds[i] = col;
         vals[i] = coef;
         i++;
         msg.detailed( " old: {} ", bound_changes[pos].get_new_bound_value() );
         msg.detailed( " {}*{} ", col, coef );
      }
      for( const std::pair<int, std::pair<int, bool>>& conflict : all_fuips )
      {
         int col = conflict.first;
         int pos = conflict.second.first;
         bool is_lower = conflict.second.second;
         REAL coef;
         assert( is_binary( col ) );
         coef = bound_changes[pos].get_new_bound_value() > 0.5 ? -1.0 : 1.0;
         if( coef < 0 )
            lhs--;
         inds[i] = col;
         vals[i] = coef;
         i++;
         msg.detailed( " old: {} ", bound_changes[pos].get_new_bound_value() );
         msg.detailed( " {}*{} ", col, coef );
      }

      msg.detailed( ">= {}", lhs );
      msg.detailed( "\n" );
      int con_size = current_conflict_set.size() + all_fuips.size();
      SparseVectorView<REAL> row_data( &vals[0], &inds[0], con_size );

      Constraint<REAL> conf_con( row_data, row_flag, lhs, 0.0 );
      constraints.push_back( conf_con );
   }
};

} // namespace papilo
