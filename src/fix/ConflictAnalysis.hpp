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
   Problem<REAL> problem;

 public:
   ConflictAnalysis( Message _msg, Num<REAL> _num, Timer _timer,
                     Problem<REAL> _problem )
       : msg( _msg ), num( _num ), timer( _timer ), problem( _problem )
   {
   }
   bool
   perform_conflict_analysis( Vec<SingleBoundChange<REAL>>& bound_changes,
                              Vec<std::pair<int, int>> infeasible_rows,
                              Vec<Constraint<REAL>>& constraints )
   {
      // bound change data as vectors for easier access
      max_depths.resize( problem.getNCols(), 0 );
      is_fixing.resize( problem.getNCols(), 0 );
      is_lower_bound.resize( problem.getNCols(), 0 );
      reason_rows.resize( problem.getNCols(), 0 );
      pos_in_bound_changes.resize( problem.getNCols(), -1 );

      // ToDo what about the integer case?
      // define two more vectors lower_bounds, upper_bounds, multiple reasons?
      for( int i = 0; i < bound_changes.size(); i++ )
      {
         if( bound_changes[i].get_depth_level() >
             max_depths[bound_changes[i].get_col()] )
         {
            is_fixing[bound_changes[i].get_col()] =
                bound_changes[i].is_manually_triggered();
            is_lower_bound[bound_changes[i].get_col()] =
                bound_changes[i].is_lower_bound();
            max_depths[bound_changes[i].get_col()] =
                bound_changes[i].get_depth_level();
            reason_rows[bound_changes[i].get_col()] =
                bound_changes[i].get_reason_row();
            pos_in_bound_changes[bound_changes[i].get_col()] = i;
         }
      }
      int number_fixings = get_number_fixings( bound_changes );
      // Only fixings (works only for binaries)
      if( number_fixings == bound_changes.size() )
      {
         simple_cut_from_fixings( bound_changes, constraints );
         return true;
      }
      else
      {
         constraint_matrix = problem.getConstraintMatrix();
         row_flags = constraint_matrix.getRowFlags();

         // compute activities
         problem.recomputeAllActivities();
         // row that led to infeasibility
         int conflict_row_index = infeasible_rows[0].second;
         // position in stack when infeasible
         int pos_at_infeasibility = infeasible_rows[0].first;

         row_activities = problem.getRowActivities();

         // last depth level
         int last_depth_level = get_last_depth_level( bound_changes );

         // Find subset of indices that explain the infeasibility
         // adds column indices in conflict_set_candidates
         explain_infeasibility( bound_changes, conflict_row_index );

         int num_vars_last_depth_level =
             get_number_variables_depth_level( max_depths, last_depth_level );

         if( num_vars_last_depth_level == 1 )
         {
            // Already at First UIP -> return conflict constraint
            // ToDo add constraint
            add_constraint( bound_changes, constraints );
            msg.info( "Only one variable at last depth level!" );
            return true;
         }
         // First-FUIP
         // Resolve as long as more than one bound changes at last depth level
         while( num_vars_last_depth_level >= 1 )
         {
            //
            int col_index;
            int position_in_conflict_set;

            get_latest_col_index_in_depth_level(
                max_depths, pos_in_bound_changes, last_depth_level, col_index,
                position_in_conflict_set );

            int antecedent_row_index = reason_rows[col_index];
            if( antecedent_row_index == -1 )
            {
               // should not happen
               // ToDo add assert
               return false;
            }
            else
            {
               SparseVectorView<REAL> antecedent_row =
                   constraint_matrix.getRowCoefficients( antecedent_row_index );
               // const int* antecedent_row_inds = conflict_row.getIndices();
               // const REAL* antecedent_row_vals = conflict_row.getValues();

               // remove latest_col_index
               conflict_set_candidates.erase( conflict_set_candidates.begin() +
                                              position_in_conflict_set );

               // resolve
               explain_infeasibility( bound_changes, antecedent_row_index );

            }

            num_vars_last_depth_level = get_number_variables_depth_level(
                max_depths, last_depth_level );
         }
      }

      add_constraint(bound_changes, constraints);
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
   ConstraintMatrix<REAL> constraint_matrix;
   Vec<papilo::RowFlags> row_flags;
   Vec<papilo::RowActivity<REAL>> row_activities;

   Vec<int> max_depths;
   Vec<bool> is_fixing;
   Vec<bool> is_lower_bound;
   Vec<int> reason_rows;
   Vec<int> pos_in_bound_changes;

   // Make a set and not vec
   // std::set<int> conflict_set_candidates;
   Vec<int> conflict_set_candidates;

   bool
   explain_infeasibility( Vec<SingleBoundChange<REAL>>& bound_changes,
                          int row_idx )
   {
      // get conflict row
      SparseVectorView<REAL> conflict_row =
          constraint_matrix.getRowCoefficients( row_idx );
      int row_length = conflict_row.getLength();
      const int* row_inds = conflict_row.getIndices();
      const REAL* row_vals = conflict_row.getValues();

      REAL lhs = constraint_matrix.getLeftHandSides()[row_idx];
      REAL rhs = constraint_matrix.getRightHandSides()[row_idx];

      // ToDo for RowFlag::kEquation
      // ToDo for Ranged row go at the corresponding infeasible side

      // For GE constraints
      if( row_flags[row_idx].test( RowFlag::kRhsInf ) )
      {
         double activity = row_activities[row_idx].max;
         for( int i = 0; i < row_length; i++ )
         {
            if( activity - lhs < 0 )
               break;
            bool is_in_bound_changes = max_depths[row_inds[i]] > 0 ? 1 : 0;
            int pos = pos_in_bound_changes[row_inds[i]];
            if( row_vals[i] > 0 && is_in_bound_changes )
            {
               if( bound_changes[pos].get_new_bound_value() <
                   problem.getUpperBounds()[row_vals[i]] )
               {
                  activity -= row_vals[i] *
                              ( problem.getUpperBounds()[row_vals[i]] -
                                bound_changes[pos].get_new_bound_value() );
                  conflict_set_candidates.push_back( row_inds[i] );
               }
            }
            else if( row_vals[i] < 0 )
            {
               if( bound_changes[pos].get_new_bound_value() >
                   problem.getLowerBounds()[row_vals[i]] )
               {
                  activity +=
                      row_vals[i] * ( bound_changes[pos].get_new_bound_value() -
                                      problem.getLowerBounds()[row_vals[i]] );
                  conflict_set_candidates.push_back( row_inds[i] );
               }
            }
         }
      }

      // For LE constraints kRhsInf
      else if( row_flags[row_idx].test( RowFlag::kLhsInf ) )
      {
         double activity = row_activities[row_idx].min;
         for( int i = 0; i < row_length; i++ )
         {
            if( activity - rhs > 0 )
               break;
            bool is_in_bound_changes = max_depths[row_inds[i]] > 0 ? 1 : 0;
            int pos = pos_in_bound_changes[row_inds[i]];
            if( row_vals[i] > 0 && is_in_bound_changes )
            {
               if( bound_changes[pos].get_new_bound_value() >
                   problem.getLowerBounds()[row_vals[i]] )
               {
                  activity +=
                      row_vals[i] * ( bound_changes[pos].get_new_bound_value() -
                                      problem.getLowerBounds()[row_vals[i]] );
                  conflict_set_candidates.push_back( row_inds[i] );
               }
            }
            else if( row_vals[i] < 0 )
            {
               if( bound_changes[pos].get_new_bound_value() <
                   problem.getUpperBounds()[row_vals[i]] )
               {
                  activity -= row_vals[i] *
                              ( problem.getLowerBounds()[row_vals[i]] -
                                bound_changes[pos].get_new_bound_value() );
                  conflict_set_candidates.push_back( row_inds[i] );
               }
            }
         }
      }

      return true;
   }
   bool
   resolve_bound_change( int col, Vec<int> row_inds,
                         SparseVectorView<REAL> antecedent_row )
   {
      // ToDo
      return true;
   }
   int
   get_last_depth_level( Vec<SingleBoundChange<REAL>> bound_changes )
   {
      return bound_changes.back().get_depth_level();
   }

   int
   get_number_variables_depth_level( Vec<int> max_depths, int depth )
   {
      int number_variables_depth_level = 0;
      for( int i = 0; i < conflict_set_candidates.size(); i++ )
      {
         if( max_depths[conflict_set_candidates[i]] == depth )
            number_variables_depth_level++;
      }
      return number_variables_depth_level;
   }
   bool
   get_latest_col_index_in_depth_level( Vec<int> max_depths,
                                        Vec<int> pos_in_bound_changes,
                                        int depth, int& col,
                                        int& position_in_conflict_set )
   {

      int latest_position = -1;
      for( int i = 0; i < conflict_set_candidates.size(); i++ )
      {
         if( ( max_depths[conflict_set_candidates[i]] == depth ) &&
             ( pos_in_bound_changes[conflict_set_candidates[i]] >
               latest_position ) )
         {
            position_in_conflict_set = i;
            col = conflict_set_candidates[i];
            latest_position = pos_in_bound_changes[conflict_set_candidates[i]];
         }
      }
      return true;
   }

   int
   get_number_fixings( Vec<SingleBoundChange<REAL>> bound_changes )
   {
      int number_fixings = 0;
      for( int i = 0; i < bound_changes.size(); i++ )
      {
         if( bound_changes[i].is_manually_triggered() )
            number_fixings++;
      }
      return number_fixings;
   }

   // ToDo
   bool
   get_conflict_data( Vec<SingleBoundChange<REAL>> bound_changes )
   {

      return true;
   }
   bool
   add_constraint( Vec<SingleBoundChange<REAL>> bound_changes, Vec<Constraint<REAL>>& constraints )
   {
      RowFlags row_flag;
      row_flag.set( RowFlag::kRhsInf );

      Vec<REAL> vals;
      Vec<int> inds;

      REAL lhs = 1.0;
      for( int i = 0; i < conflict_set_candidates.size(); i++ )
      {
         int col = conflict_set_candidates[i];
         int pos = pos_in_bound_changes[col];
         double coef =
             bound_changes[pos].get_new_bound_value() > 0.5 ? -1.0 : 1.0;
         if( coef < 0 )
            lhs--;
         inds.push_back( col );
         vals.push_back( coef );
      }

      SparseVectorView<REAL> row_data( vals.data(), inds.data(), (int) inds.size() );

      Constraint<REAL> conf_con( row_data, row_flag, lhs, 0.0 );
      constraints.push_back( conf_con );

      return true;
   }
   bool
   simple_cut_from_fixings( Vec<SingleBoundChange<REAL>> bound_changes,
                            Vec<Constraint<REAL>>& constraints )
   {

      // simplest (but also worse) cut is created from the fixings
      // e.g. Fixings x1 = x2 = x3 = 0 -> x1 + x2 + x3 >= 1
      // e.g. Fixings x1 = x2 = 0, x3 = 1 -> x1 + x2 + (1 - x3) >= 1
      RowFlags row_flag;
      row_flag.set( RowFlag::kRhsInf );

      Vec<REAL> vals;
      Vec<int> inds;
      int len;

      REAL lhs = 1.0;

      for( int i = 0; i < bound_changes.size(); i++ )
      {
         if( bound_changes[i].is_manually_triggered() )
         {
            len++;
            double coef =
                bound_changes[i].get_new_bound_value() > 0.5 ? -1.0 : 1.0;

            if( coef < 0 )
               lhs--;
            inds.push_back( bound_changes[i].get_col() );
            vals.push_back( coef );
         }
      }
      SparseVectorView<REAL> row_data( vals.data(), inds.data(), len );

      Constraint<REAL> conf_con( row_data, row_flag, lhs, 0.0 );
      constraints.push_back( conf_con );

      return true;
   }
};

} // namespace papilo
