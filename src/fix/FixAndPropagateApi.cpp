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
//#define FIX_DEBUG

#include "fix/FixAndPropagateApi.h"
#include "fix/Heuristic.hpp"

#include "papilo/io/MpsParser.hpp"
#include <string>


#ifdef FIX_DEBUG
#include "papilo/io/SolWriter.hpp"
#include "papilo/io/MpsWriter.hpp"
#endif

using namespace papilo;

Problem<double>
add_cutoff_objective( Problem<double>& problem, Num<double>& num );

void*
setup( const char* filename, int* result, int verbosity_level,
       double current_time_stamp, int add_cutoff_constraint )
{
   assert( verbosity_level <= 5 && verbosity_level >= 0 );
   assert( add_cutoff_constraint <= 1 && add_cutoff_constraint >= 0 );
   std::string filename_as_string( filename );
   boost::optional<Problem<double>> prob;
   {
      prob = MpsParser<double>::loadProblem( filename_as_string );
   }
   if( !prob )
   {
      fmt::print( "error loading problem {}\n", filename );
      *result = -1;
      return nullptr;
   }
   double tolerance = 1e-5;
   RandomGenerator random( 0 );
   Timer t{ current_time_stamp };
   Num<double> num{};
   num.setEpsilon( tolerance );
   num.setFeasTol( tolerance );
   Problem<double>* problem;
   if( add_cutoff_constraint )
   {
      problem = new Problem<double>(add_cutoff_objective(prob.get(), num));
   }
   else
      problem = new Problem<double>( prob.get() );

   problem->recomputeAllActivities();

   Message msg{};
   msg.setVerbosityLevel( static_cast<VerbosityLevel>(verbosity_level) );
   PostsolveStorage<double> storage{};

   assert( num.isEq( 80, 80 + tolerance / 10 ) );
   auto heuristic =
       new Heuristic<double>{ msg, num, random, t, *problem, storage, false };
   heuristic->setup( random );
   *result = 0;
   return heuristic;
}

Problem<double>
add_cutoff_objective( Problem<double>& problem, Num<double>& num )
{
   ProblemBuilder<double> builder;

   ConstraintMatrix<double>& matrix = problem.getConstraintMatrix();
   Vec<ColFlags>& colFlags = problem.getColFlags();
   Vec<RowFlags>& rowFlags = matrix.getRowFlags();
   Vec<int>& rowSizes = problem.getRowSizes();
   Vec<double>& coefficients = problem.getObjective().coefficients;
   Vec<double>& leftHandSides = matrix.getLeftHandSides();
   Vec<double>& rightHandSides = matrix.getRightHandSides();
   const Vec<RowActivity<double>>& activities = problem.getRowActivities();

   int nnz = matrix.getNnz();
   int ncols = problem.getNCols();
   int nrows = problem.getNRows();
   int new_nnz = 0;
   Vec<int> cut_off_indices{};
   Vec<double> cut_off_values{};
   for( int i = 0; i < ncols; i++ )
   {
      if( !num.isZero( coefficients[i] ) )
      {
         new_nnz++;
         cut_off_indices.push_back( i );
         cut_off_values.push_back( coefficients[i] );
      }
   }

   builder.reserve( nnz + new_nnz, nrows + 1, ncols );

   int rowcols_obj[new_nnz];
   double rowvals_obj[new_nnz];
   std::copy( cut_off_indices.begin(), cut_off_indices.end(), rowcols_obj );
   std::copy( cut_off_values.begin(), cut_off_values.end(), rowvals_obj );

   /* set up rows */
   builder.setNumRows( nrows + 1 );

   builder.addRowEntries( 0, new_nnz, rowcols_obj, rowvals_obj );
   builder.setRowLhs( 0, -1 );
   int rhsval = 2;
   builder.setRowRhs( 0, rhsval );
   bool infinite = true;
   builder.setRowLhsInf( 0, infinite );
   builder.setRowRhsInf( 0, infinite );
   builder.setCutoffConstraint( 0, true );
   builder.setConflictConstraint( 0, false );

   for( int i = 0; i < nrows; ++i )
   {
      const SparseVectorView<double>& view = matrix.getRowCoefficients( i );
      const int* rowcols = view.getIndices();
      const double* rowvals = view.getValues();
      int rowlen = view.getLength();
      double lhs = leftHandSides[i];
      double rhs = rightHandSides[i];

      builder.addRowEntries( i + 1, rowlen, rowcols, rowvals );
      builder.setRowLhs( i + 1, lhs );
      builder.setRowRhs( i + 1, rhs );
      builder.setRowLhsInf( i + 1, rowFlags[i].test( RowFlag::kLhsInf ) );
      builder.setRowRhsInf( i + 1, rowFlags[i].test( RowFlag::kRhsInf ) );
      builder.setConflictConstraint( i + 1, rowFlags[i].test(
                                                RowFlag::kConflictConstraint ) );
   }

   /* set up columns */
   builder.setNumCols( ncols );
   for( int i = 0; i < ncols; ++i )
   {
      builder.setColLb( i, problem.getLowerBounds()[i] );
      builder.setColUb( i, problem.getUpperBounds()[i] );

      auto flags = colFlags[i];
      builder.setColLbInf( i, flags.test( ColFlag::kLbInf ) );
      builder.setColUbInf( i, flags.test( ColFlag::kUbInf ) );
      builder.setColIntegral( i, flags.test( ColFlag::kIntegral ) );

      builder.setObj( i, coefficients[i] );
   }
   builder.setObjOffset( problem.getObjective().offset );
   return (builder.build());
}

void
delete_problem_instance( void* heuristic_void_ptr )
{
   auto heuristic = (Heuristic<double>*)( heuristic_void_ptr );
   delete heuristic;
}

int
call_algorithm( void* heuristic_void_ptr, double* cont_solution, double* result,
                int n_cols, double* current_obj_value,
                int infeasible_copy_strategy, int apply_conflicts,
                int size_of_constraints, int max_backtracks,
                int perform_one_opt, double remaining_time_in_sec )
{
   assert( infeasible_copy_strategy >= 0 && infeasible_copy_strategy <= 6 );
   assert( apply_conflicts >= 0 && apply_conflicts <= 1 );
   assert( max_backtracks >= 0 );
   assert( perform_one_opt >= 0 && perform_one_opt <= 2 );
   assert( size_of_constraints >= 0 );
   assert( n_cols ==
           ( (Heuristic<double>*)( heuristic_void_ptr ) )->problem.getNCols() );

#ifdef PAPILO_TBB
   tbb::task_arena arena( 7 );

   return arena.execute(
       [&]()
       {
#endif
          auto heuristic = (Heuristic<double>*)( heuristic_void_ptr );
          if( apply_conflicts == 1 &&
              size_of_constraints <= heuristic->get_derived_conflicts().size() )
          {
             heuristic->problem = heuristic->copy_conflicts_to_problem(
                 heuristic->problem, heuristic->get_derived_conflicts() );
             heuristic->get_message().info(
                 "added {} conflicts to the problem (rows: {})\n",
                 heuristic->get_derived_conflicts().size(),
                 heuristic->problem.getNRows() );
             heuristic->problem.recomputeAllActivities();
             heuristic->get_derived_conflicts().clear();
          }
          if( heuristic->problem.getRowFlags()[0].test(
                  RowFlag::kCutoffConstraint ) )
          {
             assert(
                 *current_obj_value <= heuristic->problem.getConstraintMatrix()
                                           .getRightHandSides()[0] ||
                 heuristic->problem.getRowFlags()[0].test( RowFlag::kRhsInf ) );
             assert(
                 heuristic->problem.getRowFlags()[0].test( RowFlag::kLhsInf ) );
             heuristic->get_message().info(
                 "update Cutoff constraint from {} ({}) to {}.\n",
                 heuristic->problem.getConstraintMatrix().getRightHandSides()[0],
                 heuristic->problem.getRowFlags()[0].test( RowFlag::kRhsInf ),
                 *current_obj_value
                 );
             heuristic->problem.getConstraintMatrix().getRightHandSides()[0] =
                 *current_obj_value;
             heuristic->problem.getRowFlags()[0].unset( RowFlag::kRhsInf );
          }
          Vec<double> sol( cont_solution, cont_solution + n_cols );
          Vec<double> res{};

#ifdef FIX_DEBUG
          SolWriter<double>::writePrimalSol(
              "lp_feasible.sol", sol,
              heuristic->problem.getObjective().coefficients, 0.0,
              heuristic->problem.getVariableNames() );
          Vec<int> row_mapping{};
          Vec<int> col_mapping{};
          for( int i = 0; i < heuristic->problem.getNRows(); i++ )
             row_mapping.push_back( i );
          for( int i = 0; i < heuristic->problem.getNCols(); i++ )
             col_mapping.push_back( i );
          MpsWriter<double>::writeProb( "test.mps", heuristic->problem,
                                        row_mapping, col_mapping );
#endif

          double i = heuristic->get_current_time();
          double time_limit = i + remaining_time_in_sec;
          double local_obj = *current_obj_value;
          bool better_solution_found = heuristic->perform_fix_and_propagate(
              sol, local_obj, res, time_limit, max_backtracks, perform_one_opt,
              false, (InfeasibleCopyStrategy)infeasible_copy_strategy );
          std::copy( res.begin(), res.end(), result );
          if( !better_solution_found )
             return false;
          assert( local_obj < *current_obj_value );
          *current_obj_value = local_obj;
          return true;
#ifdef PAPILO_TBB
       } );
#endif
}

int
call_simple_heuristic( void* heuristic_void_ptr, double* result,
                       double* current_obj_value )
{
#ifdef PAPILO_TBB
   tbb::task_arena arena( 7 );
   return arena.execute(
       [&]()
       {
#endif
          auto heuristic = (Heuristic<double>*)( heuristic_void_ptr );
          Vec<double> res{};
          heuristic->find_initial_solution( *current_obj_value, res );
          std::copy( res.begin(), res.end(), result );
          return !res.empty();
#ifdef PAPILO_TBB
       } );
#endif
}

void
perform_one_opt( void* heuristic_void_ptr, double* sol, int n_cols,
                 int perform_opt_one, double* current_obj_value,
                 double remaining_time_in_sec )
{
   double time = 0;
   Timer t = Timer(time);
   auto heuristic = (Heuristic<double>*)( heuristic_void_ptr );
   Vec<double> res{ sol, sol + n_cols };

   double i = heuristic->get_current_time();
   double time_limit = i + remaining_time_in_sec;
   ProbingView<double> view {heuristic-> problem, heuristic->get_num()};
   heuristic->perform_one_opt( perform_opt_one, res, view, *current_obj_value,
                               -1, time_limit );
   std::copy( res.begin(), res.end(), sol );
   heuristic->get_message().info("Spent {<3} in 1-opt call\n", t.getTime());
}
