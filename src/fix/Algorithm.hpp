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

#include "fix/ConflictAnalysis.hpp"
#include "fix/FixAndPropagate.hpp"
#include "fix/Heuristic.hpp"
#include "fix/VolumeAlgorithm.hpp"

#include "papilo/core/Presolve.hpp"
#include "papilo/core/Problem.hpp"
#include "papilo/core/ProblemBuilder.hpp"
#include "papilo/io/Message.hpp"
#include "papilo/misc/Num.hpp"
#include <cassert>
#include <cmath>
#include <fstream>

#include <cassert>
#include <fstream>

namespace papilo
{

// TODO: find a funny name for it
template <typename REAL>
class Algorithm
{
   Message msg;
   Num<REAL> num;
   Timer timer;
   double time_limit = 10 * 60;

 public:
   Algorithm( Message _msg, Num<REAL> _num, Timer t )
       : msg( _msg ), num( _num ), timer( t )
   {
   }

   void
   solve_problem( Problem<REAL>& problem,
                  VolumeAlgorithmParameter<REAL>& parameter )
   {
      // set up ProblemUpdate to trivialPresolve so that activities exist
      Presolve<REAL> presolve{};
      auto result = presolve.apply( problem, false );

      // TODO: add a check if presolve solved to optimality
      switch( result.status )
      {
      case papilo::PresolveStatus::kUnbounded:
      case papilo::PresolveStatus::kUnbndOrInfeas:
      case papilo::PresolveStatus::kInfeasible:
         fmt::print( "PaPILO detected infeasibility or unbounded-ness\n" );
         return;
      case papilo::PresolveStatus::kUnchanged:
      case papilo::PresolveStatus::kReduced:
         break;
      }

      if( problem.getNCols() == 0 )
      {
         msg.info( "Problem vanished during presolving\n" );
         Solution<REAL> original_solution{};
         Solution<REAL> reduced_solution{};
         Postsolve<REAL> postsolve{ msg, num };

         postsolve.undo( reduced_solution, original_solution,
                         result.postsolve );

         print_solution( original_solution.primal );

         return;
      }
      problem.recomputeAllActivities();

      Heuristic<REAL> service{msg, num, timer};
      VolumeAlgorithm<REAL> algorithm{ msg, num, timer, parameter };
      ConflictAnalysis<REAL> conflict_analysis{ msg, num, timer };

      service.setup( problem );
      REAL best_obj_value{};
      Vec<REAL> best_solution{};
      best_solution.reserve( problem.getNCols() );

      // setup data for the volume algorithm
      Vec<REAL> primal_heur_sol{};
      primal_heur_sol.reserve( problem.getNCols() );
      ProblemBuilder<REAL> builder = modify_problem( problem );
      Problem<REAL> reformulated = builder.build();
      Vec<REAL> pi;
      pi.reserve( reformulated.getNRows() );
      generate_initial_dual_solution( reformulated, pi );
      REAL min_val = calc_upper_bound_for_objective( problem );
      if( min_val == std::numeric_limits<double>::min() )
         return;

      while( true )
      {
         if(timer.getTime() >= parameter.time_limit )
            break;
         msg.info( "Starting volume algorithm\n" );
         primal_heur_sol = algorithm.volume_algorithm(
             reformulated.getObjective().coefficients,
             reformulated.getConstraintMatrix(),
             reformulated.getConstraintMatrix().getLeftHandSides(),
             reformulated.getVariableDomains(), pi, min_val );
         print_solution( primal_heur_sol );

         if(timer.getTime() >= parameter.time_limit )
            break;
         msg.info( "Starting fixing and propagating\n" );

         service.perform_fix_and_propagate( primal_heur_sol, best_obj_value, best_solution );

         if(timer.getTime() >= parameter.time_limit )
            break;

         msg.info( "Starting conflict analysis\n" );
         bool abort = conflict_analysis.perform_conflict_analysis();
         if( abort )
            return;
         // TODO: add constraint to builder and generate new problem
      }

      Solution<REAL> original_solution{};
      Solution<REAL> reduced_solution{ best_solution };
      Postsolve<REAL> postsolve{ msg, num };

      postsolve.undo( reduced_solution, original_solution, result.postsolve );

      print_solution( original_solution.primal );
      msg.info("Solving took {} seconds.\n", timer.getTime());
   }





   REAL
   calc_upper_bound_for_objective( const Problem<REAL>& problem ) const
   {
      /*
      StableSum<REAL> min_value{};
      for( int i = 0; i < problem.getNCols(); i++ )
      {
         if( num.isZero( problem.getObjective().coefficients[i] ) )
            continue;
         else if( num.isLT( problem.getObjective().coefficients[i], 0 ) )
         {
            if( problem.getColFlags()[i].test( ColFlag::kLbInf ) )
            {
               msg.error( "Could not calculate objective bound: variable {} "
                          "is unbounded",
                          i );
               return std::numeric_limits<double>::min();
            }
            min_value.add( problem.getObjective().coefficients[i] +
                           problem.getLowerBounds()[i] );
         }
         else
         {
            if( problem.getColFlags()[i].test( ColFlag::kUbInf ) )
            {
               msg.error( "Could not calculate objective bound: variable {} "
                          "is unbounded",
                          i );
               return std::numeric_limits<double>::min();
            }
            min_value.add( problem.getObjective().coefficients[i] +
                           problem.getUpperBounds()[i] );
         }
      }
      return min_value.get();
      */

      Vec<REAL> coefficients = problem.getObjective().coefficients;

      // TODO: not always a valid logic; change this!
      return *std::min_element( coefficients.begin(), coefficients.end() ) +
               problem.getObjective().offset;
   }

   void
   print_solution( const Vec<REAL>& sol )
   {
      msg.debug( "Primal solution:\n" );
      for( int i = 0; i < sol.size(); i++ )
         msg.debug( "   x[{}] = {}\n", i, sol[i] );
   }

   void
   generate_initial_dual_solution( const Problem<REAL>& problem,
                                   Vec<REAL>& dual_solution )
   {
      for( int i = 0; i < problem.getNRows(); i++ )
         dual_solution.push_back( 0 );
   }

   ProblemBuilder<REAL>
   modify_problem( Problem<REAL>& problem )
   {
      ProblemBuilder<REAL> builder;

      int nnz = 0;
      int ncols = problem.getNCols();
      int nrows = 0;
      ConstraintMatrix<REAL>& matrix = problem.getConstraintMatrix();
      Vec<ColFlags>& colFlags = problem.getColFlags();
      Vec<RowFlags>& rowFlags = matrix.getRowFlags();
      Vec<int>& rowSizes = problem.getRowSizes();
      Vec<REAL>& coefficients = problem.getObjective().coefficients;
      Vec<REAL>& leftHandSides = matrix.getLeftHandSides();
      Vec<REAL>& rightHandSides = matrix.getRightHandSides();
      const Vec<RowActivity<REAL>>& activities = problem.getRowActivities();

      int equations = 0;

      for( int i = 0; i < problem.getNRows(); i++ )
      {
         nrows++;
         int rowsize = rowSizes[i];
         nnz = nnz + rowsize;
         auto flags = rowFlags[i];
         if( flags.test( RowFlag::kEquation ) )
            equations++;
         if( flags.test( RowFlag::kEquation ) ||
             flags.test( RowFlag::kLhsInf ) || flags.test( RowFlag::kRhsInf ) )
            continue;
         nrows++;
         nnz = nnz + rowsize;
      }
      int slack_vars = nrows - equations;
      auto slack_var_upper_bounds = new double[slack_vars];

      builder.reserve( nnz, nrows, ncols + slack_vars );

      /* set up rows */
      builder.setNumRows( nrows );
      int counter = 0;
      int slack_var_counter = 0;
      for( int i = 0; i != problem.getNRows(); ++i )
      {
         const SparseVectorView<REAL>& view = matrix.getRowCoefficients( i );
         const int* rowcols = view.getIndices();
         const REAL* rowvals = view.getValues();
         int rowlen = view.getLength();
         auto flags = rowFlags[i];
         REAL lhs = leftHandSides[i];
         REAL rhs = rightHandSides[i];

         if( flags.test( RowFlag::kEquation ) )
         {
            builder.addRowEntries( counter, rowlen, rowcols, rowvals );
            builder.setRowLhs( counter, lhs );
            builder.setRowRhs( counter, rhs );
            builder.setRowLhsInf( counter, false );
            builder.setRowRhsInf( counter, false );
         }
         else if( flags.test( RowFlag::kLhsInf ) )
         {
            assert( !flags.test( RowFlag::kRhsInf ) );

            auto new_vals= new double[rowlen + 1];
            memcpy( new_vals, rowvals, rowlen * sizeof( double ) );
            new_vals[rowlen] = 1;
            auto new_indices= new int [rowlen + 1];
            memcpy( new_indices, rowcols, rowlen * sizeof( int ) );
            new_indices[rowlen] = ncols + slack_var_counter;

            builder.addRowEntries( counter, rowlen + 1, new_indices, new_vals );
            builder.setRowLhs( counter, rhs );
            builder.setRowRhs( counter, rhs );
            builder.setRowLhsInf( counter, false );
            builder.setRowRhsInf( counter, false );

            // TODO: what if activities[i].ninfmin != 0?
            slack_var_upper_bounds[slack_var_counter] = rhs - activities[i].min;

            slack_var_counter++;
         }
         else if( flags.test( RowFlag::kRhsInf ) )
         {
            assert( !flags.test( RowFlag::kLhsInf ) );

            auto new_vals = new double[rowlen + 1];
            memcpy( new_vals, rowvals, rowlen * sizeof( double ) );
            new_vals[rowlen] = -1;
            auto new_indices = new int [rowlen + 1];
            memcpy( new_indices, rowcols, rowlen * sizeof( int ) );
            new_indices[rowlen] = ncols + slack_var_counter;

            builder.addRowEntries( counter, rowlen + 1, new_indices, new_vals );
            builder.setRowLhs( counter, lhs );
            builder.setRowRhs( counter, lhs );
            builder.setRowLhsInf( counter, false );
            builder.setRowRhsInf( counter, false );

            // TODO: what if activities[i].ninfmax != 0?
            slack_var_upper_bounds[slack_var_counter] = activities[i].max - lhs;

            slack_var_counter++;
         }
         else
         {
            assert( !flags.test( RowFlag::kLhsInf ) );
            assert( !flags.test( RowFlag::kRhsInf ) );

            auto new_vals = new double[rowlen + 1];
            memcpy( new_vals, rowvals, rowlen * sizeof( double ) );
            new_vals[rowlen] = 1;
            auto new_indices = new int[rowlen + 1];
            memcpy( new_indices, rowcols, rowlen * sizeof( int ) );
            new_indices[rowlen] = ncols + slack_var_counter;

            builder.addRowEntries( counter, rowlen + 1, new_indices, new_vals );
            builder.setRowLhs( counter, rhs );
            builder.setRowRhs( counter, rhs );
            builder.setRowLhsInf( counter, false );
            builder.setRowRhsInf( counter, false );
            counter++;

            // TODO: what if activities[i].ninfmin != 0?
            slack_var_upper_bounds[slack_var_counter] = rhs - activities[i].min;

            slack_var_counter++;

            auto new_vals_2 = new double[rowlen + 1];
            memcpy( new_vals_2, rowvals, rowlen * sizeof( double ) );
            new_vals_2[rowlen] = -1;
            auto new_indices_2 = new int[rowlen + 1];
            memcpy( new_indices_2, rowcols, rowlen * sizeof( int ) );
            new_indices_2[rowlen] = ncols + slack_var_counter;

            builder.addRowEntries( counter, rowlen + 1, new_indices_2,
                  new_vals_2 );
            builder.setRowLhs( counter, lhs );
            builder.setRowRhs( counter, lhs );
            builder.setRowLhsInf( counter, false );
            builder.setRowRhsInf( counter, false );

            // TODO: what if activities[i].ninfmax != 0?
            slack_var_upper_bounds[slack_var_counter] = activities[i].max - lhs;

            slack_var_counter++;
         }
         counter++;
      }
      assert( slack_var_counter == slack_vars );

      /* set up columns */
      builder.setNumCols( ncols + slack_vars );
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

      // TODO: safe to assume slack_var ordering is consistent here and above
      //       while creating slack_var_upper_bounds?
      for( int i = ncols; i < ncols + slack_vars; ++i )
      {
         builder.setColLb( i, 0 );
         builder.setColUb( i, slack_var_upper_bounds[i] );
         builder.setColLbInf( i, false );
         builder.setColUbInf( i, false );
         builder.setColIntegral( i, false );
         builder.setObj( i, 0 );
      }

      return builder;
   }

   void
   invert( const REAL* pDouble, REAL* result, int length )
   {
      for( int i = 0; i < length; i++ )
         result[i] = pDouble[i] * -1;
   }
};

} // namespace papilo
