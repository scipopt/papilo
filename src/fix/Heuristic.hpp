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

#ifndef FIX_FIX_AND_PROPAGATE_SERVICE_HPP
#define FIX_FIX_AND_PROPAGATE_SERVICE_HPP

#include "InfeasibleCopyStrategy.hpp"
#include "fix/ConflictAnalysis.hpp"
#include "fix/Constraint.hpp"
#include "fix/FixAndPropagate.hpp"
#include "fix/InfeasibleCopyStrategy.hpp"
#include "fix/strategy/FarkasRoundingStrategy.hpp"
#include "fix/strategy/FractionalRoundingStrategy.hpp"
#include "fix/strategy/RandomRoundingStrategy.hpp"
#include "papilo/core/Objective.hpp"
#include "papilo/core/Presolve.hpp"
#include "papilo/core/ProblemBuilder.hpp"

#include <cassert>
#include <fstream>
#include <string>

using namespace papilo;

template <typename REAL>
class Heuristic
{
   Message msg;
   Num<REAL> num;
   Timer timer;
   Vec<RoundingStrategy<REAL>*> strategies{};
   Vec<Vec<REAL>> int_solutions;
   Vec<ProbingView<REAL>> views;
   Vec<REAL> cols_sorted_by_obj;
   Vec<REAL> obj_value;
   Vec<bool> infeasible_arr;
   ConflictAnalysis<REAL> conflict_analysis;
   Vec<Vec<Constraint<REAL>>> constraints{};
   Vec<Constraint<REAL>> derived_conflicts{};
   PostsolveStorage<REAL>& postsolve_storage;
   bool calculate_original = false;
   FixAndPropagate<REAL> fixAndPropagate;

 public:
   Problem<REAL>& problem;

 public:
   Heuristic( Message msg_, Num<REAL> num_, RandomGenerator random,
              Timer& timer_, Problem<REAL>& problem_,
              PostsolveStorage<REAL>& postsolve_storage_,
              bool calculate_original_ = true )
       : msg( msg_ ), num( num_ ), timer( timer_ ), strategies( {} ),
         int_solutions( {} ), views( {} ), obj_value( {} ),
         infeasible_arr( {} ), cols_sorted_by_obj( {} ), problem( problem_ ),
         conflict_analysis( { msg, num, timer, problem_ } ),
         postsolve_storage( postsolve_storage_ ),
         calculate_original( calculate_original_ ),
         fixAndPropagate( { msg, num, random } )
   {
   }

   void
   setup( RandomGenerator random )
   {
#ifdef PAPILO_TBB
      auto s1 = new FarkasRoundingStrategy<REAL>{ 0, num, false };
      auto s2 = new FarkasRoundingStrategy<REAL>{ 0, num, true };
      auto s3 = new FractionalRoundingStrategy<REAL>{ num, problem };
      auto s4 = new RandomRoundingStrategy<REAL>{ random, num };
      strategies.push_back( s1 );
      strategies.push_back( s2 );
      strategies.push_back( s3 );
      strategies.push_back( s4 );

      Vec<REAL> int_solution{};
      int_solution.resize( problem.getNCols() );

      int_solutions.push_back( { int_solution } );
      int_solutions.push_back( { int_solution } );
      int_solutions.push_back( { int_solution } );
      int_solutions.push_back( { int_solution } );

      views.push_back( { problem, num } );
      views.push_back( { problem, num } );
      views.push_back( { problem, num } );
      views.push_back( { problem, num } );

      infeasible_arr.push_back( true );
      infeasible_arr.push_back( true );
      infeasible_arr.push_back( true );
      infeasible_arr.push_back( true );

      obj_value.push_back( 0 );
      obj_value.push_back( 0 );
      obj_value.push_back( 0 );
      obj_value.push_back( 0 );

      constraints.push_back( {} );
      constraints.push_back( {} );
      constraints.push_back( {} );
      constraints.push_back( {} );

      Vec<REAL>& objective = problem.getObjective().coefficients;
      cols_sorted_by_obj.reserve( objective.size() );
      for( int i = 0; i < objective.size(); i++ )
         cols_sorted_by_obj.push_back( i );
      pdqsort( cols_sorted_by_obj.begin(), cols_sorted_by_obj.end(),
               [&]( const int a, const int b )
               {
                  return objective[a] > objective[b] ||
                         ( objective[a] == objective[b] && a > b );
               } );
#else
      Vec<REAL> int_solution{};
      int_solution.resize( problem.getNCols() );
      auto s1 = new FarkasRoundingStrategy<REAL>{ 0, num, false };
      strategies.push_back( s1 );
      int_solutions.push_back( { int_solution } );
      views.push_back( { problem, num } );
      infeasible_arr.push_back( true );
      obj_value.push_back( 0 );
#endif
   }

   Message
   get_message(){ return msg; }

   bool
   find_initial_solution( REAL& current_objective,
                          Vec<REAL>& current_best_solution )
   {
#ifdef PAPILO_TBB
      tbb::parallel_for(
          tbb::blocked_range<int>( 0, 4 ),
          [&]( const tbb::blocked_range<int>& r )
          {
         for( int i = r.begin(); i != r.end(); ++i )
         {
#else
      int i = 0;
#endif
            infeasible_arr[i] = fixAndPropagate.find_initial_solution(
                i, views[i], int_solutions[i] );
            if( infeasible_arr[i] )
            {
               obj_value[i] = 0;
               msg.info( "\t\tInitial sol {} is infeasible!\n", i,
                         obj_value[i] );
#ifdef PAPILO_TBB
               break;
#else
         return false;
#endif
            }
            obj_value[i] = calculate_obj_value( int_solutions[i] );
            msg.info( "\t\tInitial sol {} found obj value {}!\n", i,
                      obj_value[i] );
#ifndef PAPILO_TBB
            current_best_solution = int_solutions[i];
            current_objective = obj_value[i];
            return true;
#else
   }
          } );
          return evaluate( current_objective, current_best_solution,
                           InfeasibleCopyStrategy::kNone );
#endif
         }

         bool perform_fix_and_propagate(
             const Vec<REAL>& primal_heur_sol, REAL& best_obj_val,
             Vec<REAL>& current_best_solution, bool perform_backtracking = true,
             bool perform_one_opt = true, bool stop_at_infeasible = true,
             InfeasibleCopyStrategy copy = InfeasibleCopyStrategy::kNone )
         {
#ifdef PAPILO_TBB
            tbb::parallel_for(
                tbb::blocked_range<int>( 0, 4 ),
                [&]( const tbb::blocked_range<int>& r )
                {
                   for( int i = r.begin(); i != r.end(); ++i )
                   {
                      int backtracks = 0;
                      infeasible_arr[i] = fixAndPropagate.fix_and_propagate(
                          primal_heur_sol, int_solutions[i], *( strategies[i] ),
                          views[i], backtracks, perform_backtracking,
                          stop_at_infeasible );
                      if( infeasible_arr[i] )
                      {
                         obj_value[i] = 0;
                         msg.info( "\t\tPropagating {} is infeasible! "
                                   "(backtracks {})\n",
                                   i, obj_value[i], backtracks );
                         break;
                      }
                      obj_value[i] = calculate_obj_value( int_solutions[i] );
                      msg.info( "\t\tPropagating {} found obj value {}! "
                                "(backtracks {})\n",
                                i, obj_value[i], backtracks );
                   }
                } );
#else
             int backtracks = 0;
             infeasible_arr[0] = fixAndPropagate.fix_and_propagate(
                 primal_heur_sol, int_solutions[0], *( strategies[0] ),
                 views[0], backtracks, perform_backtracking,
                 stop_at_infeasible );
             if( infeasible_arr[0] )
             {
                obj_value[0] = 0;
                return false;
             }
             obj_value[0] = calculate_obj_value( int_solutions[0] );
             msg.info(
                 "\t\tPropagating {} found obj value {}! (backtracks {})\n", 0,
                 obj_value[0], backtracks );
#endif
            one_opt( perform_one_opt, stop_at_infeasible );
            for( auto& c : constraints )
            {
               derived_conflicts.insert( derived_conflicts.end(), c.begin(),
                                         c.end() );
               c.clear();
            }
            return evaluate( best_obj_val, current_best_solution, copy );
         }

         void one_opt( bool perform_one_opt, bool perform_conflict_analysis )
         {
            if( !perform_one_opt && !perform_conflict_analysis )
               return;
            for( int i = 0; i < constraints.size(); i++ )
               constraints[i].clear();

            Vec<bool> infeas_copy{ infeasible_arr };
            Vec<REAL> coefficients = problem.getObjective().coefficients;
#ifdef PAPILO_TBB
            tbb::parallel_for(
                tbb::blocked_range<int>( 0, views.size() ),
                [&]( const tbb::blocked_range<int>& r )
                {
                   for( int i = r.begin(); i != r.end(); ++i )
#else
             for( int i = 0; i != views.size(); ++i )

#endif
                   {

                      if( infeas_copy[i] )
                      {
                         if( !perform_conflict_analysis )
                            continue;
                         assert( !views[i].get_infeasible_rows().empty() );
                         auto changes = views[i].get_changes();
                         auto infeasible_rows = views[i].get_infeasible_rows();
                         conflict_analysis.perform_conflict_analysis(
                             changes, infeasible_rows, constraints[i] );
                         assert( std::all_of(
                             constraints[i].begin(), constraints[i].end(),
                             []( Constraint<REAL>& c )
                             {
                                return c.get_row_flag().test(
                                           RowFlag::kEquation ) ||
                                       !c.get_row_flag().test(
                                           RowFlag::kLhsInf );
                             } ) );
                         assert( std::all_of(
                             constraints[i].begin(), constraints[i].end(),
                             [this, i]( Constraint<REAL>& c )
                             {
                                const int* indices = c.get_data().getIndices();
                                for( int j = 0; j < c.get_data().getLength();
                                     j++ )
                                {
                                   if( indices[j] < 0 ||
                                       indices[j] > views[i]
                                                        .getProbingUpperBounds()
                                                        .size() )
                                      return false;
                                }
                                return true;
                             } ) );
                      }
                      else if( perform_one_opt )
                      {
                         assert( !infeasible_arr[i] );
                         Vec<REAL> result = { int_solutions[i] };
                         for( int j = 0; j < cols_sorted_by_obj.size(); j++ )
                         {
                            views[i].reset();
                            if( num.isZero( coefficients[j] ) )
                               break;
                            if( !problem.getColFlags()[j].test(
                                    ColFlag::kIntegral ) ||
                                problem.getLowerBounds()[j] != 0 ||
                                problem.getUpperBounds()[j] != 1 )
                               continue;
                            REAL solution_value = int_solutions[i][j];
                            if( num.isGT( coefficients[j], 0 ) )
                            {
                               if( num.isZero( solution_value ) )
                                  continue;
                               bool infeasible = fixAndPropagate.one_opt(
                                   int_solutions[i], j, 0, views[i], result );
                               if( infeasible )
                               {
                                  msg.info(
                                      "\t\t{} - OneOpt flipping variable {}: "
                                      "infeasible\n",
                                      i, j );
                                  continue;
                               }
                               REAL value = calculate_obj_value( result );
                               if( num.isGE( value, obj_value[i] ) )
                                  msg.info(
                                      "\t\t{} - OneOpt flipping variable {}: "
                                      "unsuccessful -> worse obj {}: \n",
                                      i, j, value );
                               else if( num.isLT( value, obj_value[i] ) )
                               {
                                  msg.info(
                                      "\t\t{} - OneOpt flipping variable {}: "
                                      "successful -> better obj: {}\n",
                                      i, j, value );
                                  int_solutions[i] = result;
                                  obj_value[i] = value;
                               }
                            }
                            else
                            {
                               assert( num.isLT( coefficients[j], 0 ) );
                               if( num.isZero( solution_value ) )
                                  if( !num.isZero( solution_value ) )
                                     continue;
                               bool infeasible = fixAndPropagate.one_opt(
                                   int_solutions[i], j, 1, views[i], result );
                               if( infeasible )
                               {
                                  msg.info(
                                      "\t\t{} - OneOpt flipping variable {}: "
                                      "infeasible\n",
                                      i, j );
                                  continue;
                               }
                               REAL value = calculate_obj_value( result );
                               if( num.isGE( value, obj_value[i] ) )
                                  msg.info(
                                      "\t\t{} - OneOpt flipping variable {}: "
                                      "unsuccessful -> worse obj {}: \n",
                                      i, j, value );
                               else if( num.isLT( value, obj_value[i] ) )
                               {
                                  msg.info(
                                      "\t\t{} - OneOpt flipping variable {}: "
                                      "successful -> better obj: {}\n",
                                      i, j, value );
                                  int_solutions[i] = result;
                                  obj_value[i] = value;
                               }
                               break;
                            }
                         }
                      }
                   }
#ifdef PAPILO_TBB
                } );
#endif
         }

         Vec<Constraint<REAL>>& get_derived_conflicts()
         {
            return derived_conflicts;
         }

         Problem<REAL> copy_conflicts_to_problem(
             Problem<REAL>& old_problem, Vec<Constraint<REAL>> new_conflicts )
         {
            ProblemBuilder<REAL> builder;

            ConstraintMatrix<REAL>& matrix = old_problem.getConstraintMatrix();
            Vec<ColFlags>& colFlags = old_problem.getColFlags();
            Vec<RowFlags>& rowFlags = matrix.getRowFlags();
            Vec<int>& rowSizes = old_problem.getRowSizes();
            Vec<REAL>& coefficients = old_problem.getObjective().coefficients;
            Vec<REAL>& leftHandSides = matrix.getLeftHandSides();
            Vec<REAL>& rightHandSides = matrix.getRightHandSides();
            const Vec<RowActivity<REAL>>& activities =
                old_problem.getRowActivities();

            int nnz = matrix.getNnz();
            int ncols = old_problem.getNCols();
            int nrows = old_problem.getNRows();
            int new_nnz = 0;
            int new_rows = new_conflicts.size();
            for( const auto& c : new_conflicts )
               new_nnz += c.get_data().getLength();

            builder.reserve( nnz + new_nnz, nrows + new_rows, ncols );

            /* set up rows */
            builder.setNumRows( nrows + new_rows );
            for( int i = 0; i < nrows; ++i )
            {
               const SparseVectorView<REAL>& view =
                   matrix.getRowCoefficients( i );
               const int* rowcols = view.getIndices();
               const REAL* rowvals = view.getValues();
               int rowlen = view.getLength();
               REAL lhs = leftHandSides[i];
               REAL rhs = rightHandSides[i];

               builder.addRowEntries( i, rowlen, rowcols, rowvals );
               builder.setRowLhs( i, lhs );
               builder.setRowRhs( i, rhs );
               builder.setRowLhsInf( i, rowFlags[i].test( RowFlag::kLhsInf ) );
               builder.setRowRhsInf( i, rowFlags[i].test( RowFlag::kRhsInf ) );
               builder.setHardConstraint( i, rowFlags[i].test(RowFlag::kHardConstraint) );
            }

            for( int i = 0; i < new_rows; ++i )
            {
               const SparseVectorView<REAL>& view = new_conflicts[i].get_data();
               const int* rowcols = view.getIndices();
               const REAL* rowvals = view.getValues();
               int rowlen = view.getLength();
               REAL lhs = new_conflicts[i].get_lhs();
               REAL rhs = new_conflicts[i].get_rhs();

               builder.addRowEntries( i + nrows, rowlen, rowcols, rowvals );
               builder.setRowLhs( i + nrows, lhs );
               builder.setRowRhs( i + nrows, rhs );
               builder.setRowLhsInf(
                   i + nrows,
                   new_conflicts[i].get_row_flag().test( RowFlag::kLhsInf ) );
               assert(
                   !new_conflicts[i].get_row_flag().test( RowFlag::kLhsInf ) );
               builder.setRowRhsInf(
                   i + nrows,
                   new_conflicts[i].get_row_flag().test( RowFlag::kRhsInf ) );
               assert( !rowFlags[i].test(RowFlag::kHardConstraint) );
            }

            /* set up columns */
            builder.setNumCols( ncols );
            for( int i = 0; i < ncols; ++i )
            {
               builder.setColLb( i, old_problem.getLowerBounds()[i] );
               builder.setColUb( i, old_problem.getUpperBounds()[i] );

               auto flags = colFlags[i];
               builder.setColLbInf( i, flags.test( ColFlag::kLbInf ) );
               builder.setColUbInf( i, flags.test( ColFlag::kUbInf ) );
               builder.setColIntegral( i, flags.test( ColFlag::kIntegral ) );

               builder.setObj( i, coefficients[i] );
            }
            builder.setObjOffset( old_problem.getObjective().offset );
            return builder.build();
         }

       private:
         bool evaluate( REAL & best_obj_val, Vec<REAL> & current_best_solution,
                        InfeasibleCopyStrategy copy_infeasible_sol )
         {
            bool feasible =
                std::any_of( infeasible_arr.begin(), infeasible_arr.end(),
                             []( bool b ) { return !b; } );

            // TODO: copy the best solution;
            if( !feasible )
            {
               store_solution( copy_infeasible_sol, current_best_solution );
               msg.info( "\t\tFix and Propagate did not find a feasible "
                         "solution!\n" );
               return false;
            }

            int best_index = -1;
            for( int i = 0; i < obj_value.size(); i++ )
            {
               if( !infeasible_arr[i] &&
                   ( num.isLT( obj_value[i], best_obj_val ) ||
                     ( current_best_solution.empty() && best_index == -1 ) ) )
               {
                  best_index = i;
                  best_obj_val = obj_value[i];
               }
            }
            if( best_index == -1 )
            {
               msg.info( "\t\tFix and Propagate did not improve the current "
                         "solution!\n" );
               return false;
            }

            if( current_best_solution.empty() )
               msg.info(
                   "\t\tFix and Propagate found an initial solution: {}!\n",
                   best_obj_val );
            else
               msg.info( "\t\tFix and Propagate found a new solution: {}!\n",
                         best_obj_val );

            current_best_solution = int_solutions[best_index];
            assert( best_obj_val == obj_value[best_index] );
            return true;
         }

         void store_solution( const InfeasibleCopyStrategy& copy_infeasible_sol,
                              Vec<REAL>& current_best_solution )
         {
            int best_index = -1;
            REAL val{};
            switch( copy_infeasible_sol )
            {
            case InfeasibleCopyStrategy::kNone:
               return;
            case InfeasibleCopyStrategy::kBestObjective:
            {
               best_index = 0;
               val = calculate_obj_value( int_solutions[0], false );
               for( int i = 1; i < int_solutions.size(); i++ )
               {
                  REAL v = calculate_obj_value( int_solutions[i], false );
                  if( num.isLT( v, val ) )
                  {
                     val = v;
                     best_index = i;
                  }
               }
               break;
            }
            case InfeasibleCopyStrategy::kWorstObjective:
            {
               best_index = 0;
               val = calculate_obj_value( int_solutions[0], false );
               for( int i = 1; i < int_solutions.size(); i++ )
               {
                  REAL v = calculate_obj_value( int_solutions[i], false );
                  if( num.isGT( v, val ) )
                  {
                     val = v;
                     best_index = i;
                  }
               }
               break;
            }
            case InfeasibleCopyStrategy::kHighestDepthOfFirstConflict:
            {
               best_index = 0;
               assert( views[0].get_infeasible_rows().size() >= 1 );
               int aux = views[0].get_infeasible_rows()[0].second;
               assert( views[0].get_changes().size() >= aux );
               val = views[0].get_changes()[aux].get_depth_level();
               for( int i = 1; i < int_solutions.size(); i++ )
               {
                  assert( views[i].get_infeasible_rows().size() >= 1 );
                  aux = views[i].get_infeasible_rows()[0].second;
                  assert( views[i].get_changes().size() >= aux );
                  REAL v = views[i].get_changes()[aux].get_depth_level();
                  if( num.isGT( v, val ) )
                  {
                     val = v;
                     best_index = i;
                  }
               }
               break;
            }
            case InfeasibleCopyStrategy::kLowestDepthOfFirstConflict:
            {
               best_index = 0;
               assert( views[0].get_infeasible_rows().size() >= 1 );
               int aux = views[0].get_infeasible_rows()[0].second;
               assert( views[0].get_changes().size() >= aux );
               val = views[0].get_changes()[aux].get_depth_level();
               for( int i = 1; i < int_solutions.size(); i++ )
               {
                  assert( views[i].get_infeasible_rows().size() >= 1 );
                  aux = views[i].get_infeasible_rows()[0].second;
                  assert( views[i].get_changes().size() >= aux );
                  REAL v = views[i].get_changes()[aux].get_depth_level();
                  if( num.isLT( v, val ) )
                  {
                     val = v;
                     best_index = i;
                  }
               }
               break;
            }
            case InfeasibleCopyStrategy::kLeastInfeasibleRows:
            {
               best_index = 0;
               assert( views[0].get_infeasible_rows().size() >= 1 );
               val = views[0].get_infeasible_rows().size();
               for( int i = 1; i < int_solutions.size(); i++ )
               {
                  assert( views[i].get_changes().size() >= 1 );
                  REAL v = views[0].get_infeasible_rows().size();
                  if( num.isLT( v, val ) )
                  {
                     val = v;
                     best_index = i;
                  }
               }
               break;
            }
            case InfeasibleCopyStrategy::kMostInfeasibleRows:
            {
               best_index = 0;
               assert( views[0].get_infeasible_rows().size() >= 1 );
               val = views[0].get_infeasible_rows().size();
               for( int i = 1; i < int_solutions.size(); i++ )
               {
                  assert( views[i].get_changes().size() >= 1 );
                  REAL v = views[0].get_infeasible_rows().size();
                  if( num.isGT( v, val ) )
                  {
                     val = v;
                     best_index = i;
                  }
               }
               break;
            }
            default:
               assert( false );
            }
            assert( best_index != -1 );
            current_best_solution = int_solutions[best_index];
         }

         REAL calculate_obj_value( const Vec<REAL>& reduced,
                                   bool check_validation = true ) const
         {
            if( calculate_original )
            {
               Solution<REAL> original_solution{};
               Solution<REAL> reduced_solution{ reduced };
               Message quiet{};
               quiet.setVerbosityLevel( papilo::VerbosityLevel::kQuiet );
               Postsolve<REAL> postsolve{ quiet, num };
               auto status = postsolve.undo(
                   reduced_solution, original_solution, postsolve_storage );
               assert( !check_validation || status == PostsolveStatus::kOk );
               return postsolve_storage.getOriginalProblem()
                   .computeSolObjective( original_solution.primal );
            }
            else
            {
               StableSum<REAL> sum{};
               Vec<REAL>& coefficients = problem.getObjective().coefficients;
               for( int j = 0; j < reduced.size(); j++ )
                  sum.add( reduced[j] * coefficients[j] );
               return sum.get();
            }
         }
};

#endif
