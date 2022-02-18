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

#include "fix/FixAndPropagate.hpp"
#include "fix/Heuristic.hpp"
#include "fix/Constraint.hpp"
#include "fix/VolumeAlgorithm.hpp"
#include "fix/strategy/FarkasRoundingStrategy.hpp"
#include "fix/strategy/FractionalRoundingStrategy.hpp"
#include "fix/strategy/RandomRoundingStrategy.hpp"
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

// TODO: statistics + evaluation

template <typename REAL>
class Algorithm
{
   Message msg;
   Num<REAL> num;
   Timer timer;
   AlgorithmParameter alg_parameter;
   PresolveOptions presolveOptions;

 public:
   Algorithm( Message _msg, Num<REAL> _num, Timer t )
       : msg( _msg ), num( _num ), timer( t )
   {
   }

   void
   solve_problem( Problem<REAL>& problem )
   {
#ifdef PAPILO_TBB
      tbb::task_arena arena( alg_parameter.threads );
#endif

#ifdef PAPILO_TBB
      return arena.execute(
          [this, &problem]()
          {
#endif
             msg.info( "Starting Algorithm\n Starting presolving:\n" );
             Presolve<REAL> presolve{};
             presolve.getPresolveOptions().threads = alg_parameter.threads;
             presolve.getPresolveOptions().tlim = alg_parameter.time_limit;
             presolve.addDefaultPresolvers();

             REAL box_upper_bound_volume = calc_box_upper_bound( problem );

             auto result = presolve.apply( problem, false );

             switch( result.status )
             {
             case papilo::PresolveStatus::kUnbounded:
             case papilo::PresolveStatus::kUnbndOrInfeas:
             case papilo::PresolveStatus::kInfeasible:
                msg.info( "PaPILO detected infeasibility or unbounded-ness.\n "
                          "Stopped after {}s\n",
                          timer.getTime() );
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

             RandomGenerator random(alg_parameter.seed);
             VolumeAlgorithm<REAL> algorithm{ msg, num, timer, alg_parameter,
                                              result.postsolve };

             // setup data for the volume algorithm
             int slack_vars = 0;
             ProblemBuilder<REAL> builder =
                 modify_problem( problem, slack_vars );
             Problem<REAL> reformulated = builder.build();

             Heuristic<REAL> service{ msg, num, random, timer, reformulated, result.postsolve };
             service.setup(random);
             reformulated.recomputeAllActivities();


             REAL best_obj_value = std::numeric_limits<REAL>::max();
             Vec<REAL> best_solution{};
             best_solution.reserve( problem.getNCols() );
             Vec<REAL> primal_heur_sol{};
             primal_heur_sol.reserve( reformulated.getNCols() );
             Vec<REAL> pi;
             pi.reserve( reformulated.getNRows() );
             generate_initial_dual_solution( reformulated, pi );
             Vec<REAL> pi_conflicts;
             Vec<Vec<SingleBoundChange<REAL>>> bound_changes;
             Vec<std::pair<int, int>> infeasible_rows;


             msg.info( "\tStarting primal heuristics - {:.3} s\n",
                       timer.getTime() );
             bool solution_found = service.find_initial_solution(best_obj_value,
                                                                 best_solution);

             if( solution_found )
                box_upper_bound_volume = best_obj_value;
             if( box_upper_bound_volume == std::numeric_limits<double>::min() )
                return;

             int round_counter = 0;
             int round_first_solution = -1;
             int round_best_solution = -1;
             while( true )
             {
                msg.info( "Starting round {} - {:.3}s\n", round_counter,
                          timer.getTime() );
                if( timer.getTime() >= alg_parameter.time_limit )
                   break;
                msg.info( "\tStarting volume algorithm - {:.3} s\n",
                          timer.getTime() );
                primal_heur_sol = algorithm.volume_algorithm(
                    reformulated.getObjective().coefficients,
                    reformulated.getConstraintMatrix(),
                    service.get_derived_conflicts(),
                    reformulated.getConstraintMatrix().getLeftHandSides(),
                    reformulated.getVariableDomains(),
                    reformulated.getNumIntegralCols(),
                    box_upper_bound_volume, solution_found,
                    pi, pi_conflicts );
                print_solution( primal_heur_sol );

                if( timer.getTime() >= alg_parameter.time_limit )
                   break;
                msg.info( "\tStarting fixing and propagating - {:.3} s\n",
                          timer.getTime() );

                //                Vec<REAL> sub( primal_heur_sol.begin(),
                //                               primal_heur_sol.end() -
                //                               slack_vars );
                //
                //                assert( sub.size() == problem.getNCols() );
                assert( problem.getNCols() == primal_heur_sol.size() );
                bool sol_updated = service.perform_fix_and_propagate(
                    primal_heur_sol, best_obj_value, best_solution );
                if( sol_updated )
                {
                   if( round_first_solution == -1 )
                      round_first_solution = round_counter;
                   round_best_solution = round_counter;
                }

                if( timer.getTime() >= alg_parameter.time_limit )
                   break;

                auto conflicts = (int) service.get_derived_conflicts().size();
                if( conflicts == 0 )
                   break;

                if( alg_parameter.copy_conflicts_to_problem &&
                    conflicts >
                        alg_parameter.size_of_conflicts_to_be_copied )
                {

                   reformulated =
                       service.copy_conflicts_to_problem( reformulated, service.get_derived_conflicts() );
                   msg.info(
                       "\tCopied {} conflicts to the (f&p) problem (constraints {}) - {:.3} s\n",
                       conflicts, reformulated.getNRows(), timer.getTime() );
                   pi_conflicts.clear();
                   reformulated.recomputeAllActivities();
                   service.get_derived_conflicts().clear();
                }
                else
                   msg.info( "\tFound {} conflicts (treated separately) - {:.3} s\n",
                             conflicts, timer.getTime() );

                round_counter++;

                pi_conflicts.resize( pi_conflicts.size() + conflicts, 0 );
             }

             Solution<REAL> original_solution{};
             Solution<REAL> reduced_solution{ best_solution };
             Postsolve<REAL> postsolve{ msg, num };

             auto status = postsolve.undo( reduced_solution, original_solution,
                                           result.postsolve );

             print_solution( original_solution.primal );
             msg.info( "Algorithm with objective value {:<3} finished after "
                       "{:.3} seconds.\n",
                       result.postsolve.problem.computeSolObjective(
                           original_solution.primal ),
                       timer.getTime() );
             msg.info( "First solution found in round {}.\n",
                       round_first_solution );
             msg.info( "Best solution found in round {}.\n",
                       round_best_solution );
#ifdef PAPILO_TBB
          } );
#endif
   }


   REAL
   calc_box_upper_bound( const Problem<REAL>& problem ) const
   {
      StableSum<REAL> min_value{};
      for( int i = 0; i < problem.getNCols(); i++ )
      {
         if( num.isZero( problem.getObjective().coefficients[i] ) )
            continue;
         else if( num.isLT( problem.getObjective().coefficients[i], 0 ) )
         {
            if( problem.getColFlags()[i].test( ColFlag::kLbInf ) )
            {
               // TODO: OK? - based on this UB's usage in volume algo code?
               continue;
            }
            min_value.add( problem.getObjective().coefficients[i] +
                           problem.getLowerBounds()[i] );
         }
         else
         {
            if( problem.getColFlags()[i].test( ColFlag::kUbInf ) )
            {
               continue;
            }
            min_value.add( problem.getObjective().coefficients[i] +
                           problem.getUpperBounds()[i] );
         }
      }
      return min_value.get();
   }

   void
   print_solution( const Vec<REAL>& sol )
   {
      if( msg.getVerbosityLevel() == VerbosityLevel::kDetailed )
      {
         msg.detailed( "Primal solution:\n" );
         for( int i = 0; i < sol.size(); i++ )
            msg.detailed( "   x[{}] = {}\n", i, sol[i] );
      }
   }

   void
   generate_initial_dual_solution( const Problem<REAL>& problem,
                                   Vec<REAL>& dual_solution )
   {
      for( int i = 0; i < problem.getNRows(); i++ )
         dual_solution.push_back( 0 );
   }

   ProblemBuilder<REAL>
   modify_problem( Problem<REAL>& problem, int& slack_vars )
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

      int hard_constraints = 0;
      for( int i = 0; i < problem.getNRows(); i++ )
      {
         int rowsize = rowSizes[i];

         if( !num.isEq( get_max_min_factor( matrix.getRowCoefficients( i ) ),
                        alg_parameter.threshold_hard_constraints ) )
         {
            hard_constraints++;
            rowFlags[i].set( RowFlag::kHardConstraint );
         }
         nnz = nnz + rowsize;
         nrows++;
         if( rowFlags[i].test( RowFlag::kEquation ) ||
             rowFlags[i].test( RowFlag::kLhsInf ) ||
             rowFlags[i].test( RowFlag::kRhsInf ) )
            continue;
         nrows++;
         nnz = nnz + rowsize;
      }

      msg.info( "\n{} of the {} rows were considered hard and were excluded.\n",
                hard_constraints, problem.getNRows() );

      slack_vars = 0;
      auto slack_var_upper_bounds = new double[slack_vars];

      builder.reserve( nnz, nrows, ncols + slack_vars );

      /* set up rows */
      builder.setNumRows( nrows );
      int counter = 0;
      int slack_var_counter = 0;
      for( int i = 0; i < problem.getNRows(); ++i )
      {
         auto flags = rowFlags[i];
         if( flags.test( RowFlag::kHardConstraint ) )
         {
            assert(
                !num.isEq( get_max_min_factor( matrix.getRowCoefficients( i ) ),
                           alg_parameter.threshold_hard_constraints ) );
         }
         const SparseVectorView<REAL>& view = matrix.getRowCoefficients( i );
         const int* rowcols = view.getIndices();
         const REAL* rowvals = view.getValues();
         int rowlen = view.getLength();
         REAL lhs = leftHandSides[i];
         REAL rhs = rightHandSides[i];

         if( flags.test( RowFlag::kEquation ) )
         {
            assert( num.isEq( lhs, rhs ) );

            builder.addRowEntries( counter, rowlen, rowcols, rowvals );
            builder.setRowLhs( counter, lhs );
            builder.setRowRhs( counter, rhs );
            builder.setRowLhsInf( counter, false );
            builder.setRowRhsInf( counter, false );
            builder.setHardConstraint( counter, flags.test( RowFlag::kHardConstraint ) );
         }
         else if( flags.test( RowFlag::kLhsInf ) )
         {
            assert( !flags.test( RowFlag::kRhsInf ) );

            auto neg_vals = new double[rowlen];
            invert( rowvals, neg_vals, rowlen );

            builder.addRowEntries( counter, rowlen, rowcols, neg_vals );
            builder.setRowLhs( counter, -1.0 * rhs );
            builder.setRowRhs( counter, 0 );
            builder.setRowLhsInf( counter, false );
            builder.setRowRhsInf( counter, true );
            builder.setHardConstraint( counter, flags.test( RowFlag::kHardConstraint ) );
         }
         else if( flags.test( RowFlag::kRhsInf ) )
         {
            assert( !flags.test( RowFlag::kLhsInf ) );

            builder.addRowEntries( counter, rowlen, rowcols, rowvals );
            builder.setRowLhs( counter, lhs );
            builder.setRowRhs( counter, 0 );
            builder.setRowLhsInf( counter, false );
            builder.setRowRhsInf( counter, true );
            builder.setHardConstraint( counter, flags.test( RowFlag::kHardConstraint ) );
         }
         else
         {
            assert( !flags.test( RowFlag::kLhsInf ) );
            assert( !flags.test( RowFlag::kRhsInf ) );

            auto neg_vals = new double[rowlen];
            invert( rowvals, neg_vals, rowlen );

            builder.addRowEntries( counter, rowlen, rowcols, neg_vals );
            builder.setRowLhs( counter, -1.0 * rhs );
            builder.setRowRhs( counter, 0 );
            builder.setRowLhsInf( counter, false );
            builder.setRowRhsInf( counter, true );
            builder.setHardConstraint( counter, flags.test( RowFlag::kHardConstraint ) );
            counter++;

            builder.addRowEntries( counter, rowlen, rowcols, rowvals );
            builder.setRowLhs( counter, lhs );
            builder.setRowRhs( counter, 0 );
            builder.setRowLhsInf( counter, false );
            builder.setRowRhsInf( counter, true );
            builder.setHardConstraint( counter, flags.test( RowFlag::kHardConstraint ) );
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

      for( int i = ncols; i < ncols + slack_vars; ++i )
      {
         builder.setColLb( i, 0 );
         builder.setColUb( i, 8 * slack_var_upper_bounds[i - ncols] );
         builder.setColLbInf( i, false );
         builder.setColUbInf( i, false );
         builder.setColIntegral( i, false );
         builder.setObj( i, 0 );
      }

      return builder;
   }

   REAL
   get_max_min_factor( const SparseVectorView<REAL>& row_data ) const
   {
      assert( row_data.getLength() > 0 );
      auto pair = row_data.getMinMaxAbsValue();
      return pair.second / pair.first;
   }

   void
   invert( const REAL* pDouble, REAL* result, int length )
   {
      for( int i = 0; i < length; i++ )
         result[i] = pDouble[i] * -1;
   }

   ParameterSet
   getParameters()
   {
      ParameterSet paramSet;
      msg.addParameters( paramSet );
      alg_parameter.addParameters( paramSet );
      presolveOptions.addParameters( paramSet );
      return paramSet;
   }
};

} // namespace papilo
