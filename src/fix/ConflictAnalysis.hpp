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
   simple_cut_from_fixings( Vec<SingleBoundChange<REAL>> bound_changes,
                            bool all_fixings_binary, Vec<int>& length,
                            Vec<int*>& indices, Vec<REAL*>& values,
                            Vec<RowFlags>& flags, Vec<REAL>& lhs,
                            Vec<REAL>& rhs )
   {
      length.push_back( 0 );
      lhs.push_back( 1 );
      rhs.push_back( 0 );
      RowFlags row_flag;
      row_flag.set( RowFlag::kRhsInf );
      flags.push_back( row_flag );
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
            if( bound_changes[i].is_manually_triggered() )
            {
               length[0]++;
               row_indices.push_back( bound_changes[i].get_col() );
               double coef =
                   bound_changes[i].get_new_bound_value() > 0.5 ? -1.0 : 1.0;

               if( bound_changes[i].get_new_bound_value() > 0.5 )
                  lhs[0]--;
               row_values.push_back( coef );
            }
         }
         indices.push_back( row_indices.data() );
         values.push_back( row_values.data() );
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

      // bound change data as vectors for easier access
      Vec<int> max_depths;
      max_depths.resize( problem.getNCols(), 0 );
      Vec<bool> is_fixing;
      is_fixing.resize( problem.getNCols(), 0 );
      Vec<bool> is_lower_bound;
      is_lower_bound.resize( problem.getNCols(), 0 );
      Vec<int> reason_rows;
      reason_rows.resize( problem.getNCols(), 0 );
      Vec<int> pos_in_bound_changes;
      pos_in_bound_changes.resize( problem.getNCols(), 0 );

      for( int i = 0; i < bound_changes.size() - 1; i++ )
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
         simple_cut_from_fixings( bound_changes, all_fixings_binary, length,
                                  indices, values, flags, lhs, rhs );
         return true;
      }
      else
      {
         // pointer to constraint matrix
         ConstraintMatrix<REAL> constraint_matrix =
             problem.getConstraintMatrix();

         // row that led to infeasibility
         int conflict_row_index = bound_changes.back().get_reason_row();
         SparseVectorView<REAL> conflict_row =
             constraint_matrix.getRowCoefficients( conflict_row_index );
         int row_length = conflict_row.getLength();
         const int* row_inds = conflict_row.getIndices();
         const REAL* row_vals = conflict_row.getValues();

         // last depth level
         int last_depth_level = get_last_depth_level( bound_changes );

         // initial conflict set (maybe just use indices for easy access?)
         Vec<int> conflict_set_candidates;
         for( int i = 0; i < row_length; i++ )
         {
            // ToDo add variables whose bound change affects the activities
            conflict_set_candidates.push_back( row_inds[i] );
         }

         int num_vars_last_depth_level = get_number_variables_depth_level(
             max_depths, conflict_set_candidates, last_depth_level );

         if (num_vars_last_depth_level == 1)
         {
            // Already at First UIP -> return conflict constraint
            // ToDo add constraint
            msg.info("Only one variable at last depth level!")
            return true
         }
         // First-FUIP
         // Resolve as long as more than one bound changes at last depth level
         while( num_vars_last_depth_level > 1 )
         {
            int col_index = get_latest_col_index_in_depth_level(
                max_depths, pos_in_bound_changes, conflict_set_candidates,
                last_depth_level );

            int antecedent_row_index = reason_rows[col_index];
            if (antecedent_row_index == -1) 
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

                  // resolve
                  resolve_bound_change( col_index, conflict_set_candidates,
                                     antecedent_row );

            }

            num_vars_last_depth_level = get_number_variables_depth_level(
                max_depths, conflict_set_candidates, last_depth_level );
         }
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
   get_number_variables_depth_level( Vec<int> max_depths, Vec<int> row_inds,
                                     int depth )
   {
      int number_variables_depth_level = 0;
      for( int i = 0; i < row_inds.size(); i++ )
      {
         if( max_depths[row_inds[i]] == depth )
            number_variables_depth_level++;
      }
      return number_variables_depth_level;
   }
   int
   get_latest_col_index_in_depth_level( Vec<int> max_depths,
                                        Vec<int> pos_in_bound_changes,
                                        Vec<int> row_inds, int depth )
   {

      int col = -1;
      int latest_position = -1;

      for( int i = 0; i < row_inds.size(); i++ )
      {
         if( ( max_depths[row_inds[i]] == depth ) &&
             ( pos_in_bound_changes[row_inds[i]] > latest_position ) )
         {
            col = row_inds[i];
            latest_position = pos_in_bound_changes[row_inds[i]];
         }
      }
      return col;
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
   get_conflict_set()
   {
      return true;
   }
};

} // namespace papilo
