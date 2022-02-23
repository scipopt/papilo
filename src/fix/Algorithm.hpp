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
#include "fix/FixAndPropagate.hpp"
#include "fix/Heuristic.hpp"
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

             RandomGenerator random( alg_parameter.seed );
             VolumeAlgorithm<REAL> algorithm{ msg, num, timer, alg_parameter,
                                              result.postsolve };

             // setup data for the volume algorithm
             ProblemBuilder<REAL> builder = modify_problem( problem );
             Problem<REAL> reformulated = builder.build();

             REAL offset_for_cutoff = 1;
             if( alg_parameter.use_cutoff_constraint )
             {
                reformulated = add_cutoff_objective( reformulated );
                reformulated.recomputeAllActivities();
                offset_for_cutoff = calculate_cutoff_offset( reformulated );
             }
             assert( num.isGT( offset_for_cutoff, 0 ) &&
                     num.isLE( offset_for_cutoff, 1 ) );

             int n_rows_A_no_conflicts = reformulated.getNRows();
             int n_hard_constraints = 0;
             bool detect_hard_constraints =
                 alg_parameter.detect_hard_constraints;
             REAL threshold_hard_constraints =
                 alg_parameter.threshold_hard_constraints;
             REAL threshold_hard_constraints_incr_factor =
                 alg_parameter.threshold_hard_constraints_incr_factor;

             Heuristic<REAL> service{ msg,   num,          random,
                                      timer, reformulated, result.postsolve };
             service.setup( random );
             reformulated.recomputeAllActivities();

             assert( problem.getNCols() == reformulated.getNCols() );

             REAL best_obj_value = std::numeric_limits<REAL>::max();
             Vec<REAL> best_solution{};
             best_solution.reserve( problem.getNCols() );

             // volume algorithm
             Vec<REAL> primal_heur_sol{};
             primal_heur_sol.reserve( reformulated.getNCols() );
             Vec<REAL> pi;
             Vec<REAL> pi_conflicts;
             pi.reserve( n_rows_A_no_conflicts );
             generate_initial_dual_solution( reformulated, pi );

             // Fix and Propagate
             Vec<Vec<SingleBoundChange<REAL>>> bound_changes;
             Vec<std::pair<int, int>> infeasible_rows;

             msg.info( "\tStarting primal heuristics - {:.3} s\n",
                       timer.getTime() );
             bool solution_found =
                 service.find_initial_solution( best_obj_value, best_solution );

             if( solution_found )
             {
                box_upper_bound_volume = best_obj_value;
                if( alg_parameter.use_cutoff_constraint )
                {
                   REAL cutoff = calculate_objective_of_reduced_problem(
                                     reformulated, best_solution ) -
                                 offset_for_cutoff;
                   reformulated.getRowFlags()[0].unset( RowFlag::kRhsInf );
                   reformulated.getConstraintMatrix().getRightHandSides()[0] =
                       cutoff;
                   reformulated.recomputeAllActivities();
                }
             }

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
                    reformulated.getNumIntegralCols(), box_upper_bound_volume,
                    solution_found, threshold_hard_constraints,
                    n_rows_A_no_conflicts, detect_hard_constraints,
                    n_hard_constraints, pi, pi_conflicts );
                print_solution( primal_heur_sol );

                if( timer.getTime() >= alg_parameter.time_limit )
                   break;
                msg.info( "\tStarting fixing and propagating - {:.3} s\n",
                          timer.getTime() );

                assert( problem.getNCols() == primal_heur_sol.size() );
                auto old_conflicts =
                    (int)service.get_derived_conflicts().size();

                bool sol_updated = service.perform_fix_and_propagate(
                    primal_heur_sol, best_obj_value, best_solution );
                if( sol_updated )
                {
                   if( round_first_solution == -1 )
                      round_first_solution = round_counter;
                   round_best_solution = round_counter;
                   if( alg_parameter.use_cutoff_constraint )
                   {
                      REAL cutoff = calculate_objective_of_reduced_problem(
                                        problem, best_solution ) -
                                    offset_for_cutoff;
                      problem.getRowFlags()[0].unset( RowFlag::kRhsInf );
                      problem.getConstraintMatrix().getRightHandSides()[0] =
                          cutoff;
                      problem.recomputeAllActivities();
                   }
                }

                if( timer.getTime() >= alg_parameter.time_limit )
                   break;

                auto new_conflicts =
                    (int)service.get_derived_conflicts().size() - old_conflicts;
                if( new_conflicts == 0 )
                {
                   msg.info( "\tNo conflict could be generated - {:.3} s\n",
                             timer.getTime() );

                   // exit if no hard constraints can be added or the cutoff
                   // constraint was modified or cutoff was not updated
                   if( ( !alg_parameter.threshold_hard_constraints_vary ||
                       ( alg_parameter.threshold_hard_constraints_vary &&
                         n_hard_constraints == 0 ) ) &&
                       ( !alg_parameter.use_cutoff_constraint ||
                       ( alg_parameter.use_cutoff_constraint && !sol_updated ) )
                     )
                   {
                      msg.info( "\tNo new constraints could be added. Terminating now...\n" );
                      break;
                   }

                   // add hard constraints only if cutoff constraint was not modified
                   if( !alg_parameter.use_cutoff_constraint ||
                       ( alg_parameter.use_cutoff_constraint && !sol_updated ) )
                      if( alg_parameter.threshold_hard_constraints_vary )
                      {
                         assert( n_hard_constraints > 0 );
                         msg.info( "\tIncreasing the threshold for hard "
                                   "constraints from {} to {}.\n",
                                   threshold_hard_constraints,
                                   threshold_hard_constraints *
                                       threshold_hard_constraints_incr_factor );
                         detect_hard_constraints = true;
                         threshold_hard_constraints *=
                             threshold_hard_constraints_incr_factor;
                      }
                }
                else
                {
                   int cont_solution_not_feasible_for_conflicts = 0;
                   for( const Constraint<REAL> item :
                        service.get_derived_conflicts() )
                   {
                      auto data = item.get_data();
                      StableSum<REAL> sum{};
                      for( int i = 0; i < data.getLength(); i++ )
                         sum.add( data.getValues()[i] *
                                  primal_heur_sol[data.getIndices()[i]] );
                      REAL value = sum.get();
                      if( !item.get_row_flag().test( RowFlag::kRhsInf ) &&
                          num.isLE( value, item.get_rhs() ) )
                         continue;
                      if( !item.get_row_flag().test( RowFlag::kLhsInf ) &&
                          num.isGE( value, item.get_lhs() ) )
                         continue;
                      cont_solution_not_feasible_for_conflicts++;
                   }

                   if( alg_parameter.copy_conflicts_to_problem &&
                       ( new_conflicts + old_conflicts ) >
                           alg_parameter.size_of_conflicts_to_be_copied )
                   {
                      reformulated = service.copy_conflicts_to_problem(
                          reformulated, service.get_derived_conflicts() );
                      msg.info( "\tCopied {} conflicts ({} not feasible for "
                                "current solution)  to the (f&p) problem "
                                "(constraints {}) - {:.3} s\n",
                                new_conflicts + old_conflicts,
                                cont_solution_not_feasible_for_conflicts,
                                reformulated.getNRows(), timer.getTime() );
                      reformulated.recomputeAllActivities();
                      service.get_derived_conflicts().clear();
                      pi.insert( pi.end(), pi_conflicts.begin(),
                                 pi_conflicts.end() );
                      pi.resize( pi.size() + new_conflicts, 0 );
                      pi_conflicts.clear();
                   }
                   else
                   {
                      msg.info( "\tFound {} conflicts ({} not feasible for "
                                "current solution) (treated separately) - "
                                "{:.3} s\n",
                                new_conflicts,
                                cont_solution_not_feasible_for_conflicts,
                                timer.getTime() );
                      pi_conflicts.resize( pi_conflicts.size() + new_conflicts,
                                           0 );
                   }
                }

                round_counter++;
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
   calculate_objective_of_reduced_problem(
       const Problem<REAL>& problem, const Vec<REAL>& best_solution ) const
   {
      StableSum<REAL> obj{};
      for( int i = 0; i < problem.getNCols(); i++ )
         obj.add( problem.getObjective().coefficients[i] * best_solution[i] );
      return obj.get();
   }

   REAL
   calculate_cutoff_offset( const Problem<REAL>& problem ) const
   {
      for( int i = 0; i < problem.getNCols(); i++ )
      {
         if( !num.isZero( problem.getObjective().coefficients[i] ) &&
             !problem.getColFlags()[i].test( ColFlag::kIntegral ) &&
             !num.isIntegral( problem.getObjective().coefficients[i] ) )
            return 2 * num.getEpsilon();
      }
      return 1;
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

   Problem<REAL>
   add_cutoff_objective( Problem<REAL>& problem )
   {
      ProblemBuilder<REAL> builder;

      ConstraintMatrix<REAL>& matrix = problem.getConstraintMatrix();
      Vec<ColFlags>& colFlags = problem.getColFlags();
      Vec<RowFlags>& rowFlags = matrix.getRowFlags();
      Vec<int>& rowSizes = problem.getRowSizes();
      Vec<REAL>& coefficients = problem.getObjective().coefficients;
      Vec<REAL>& leftHandSides = matrix.getLeftHandSides();
      Vec<REAL>& rightHandSides = matrix.getRightHandSides();
      const Vec<RowActivity<REAL>>& activities = problem.getRowActivities();

      int nnz = matrix.getNnz();
      int ncols = problem.getNCols();
      int nrows = problem.getNRows();
      int new_nnz = 0;
      Vec<int> cut_off_indices{};
      Vec<int> cut_off_values{};
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
      REAL rowvals_obj[new_nnz];
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

      for( int i = 0; i < nrows; ++i )
      {
         const SparseVectorView<REAL>& view = matrix.getRowCoefficients( i );
         const int* rowcols = view.getIndices();
         const REAL* rowvals = view.getValues();
         int rowlen = view.getLength();
         REAL lhs = leftHandSides[i];
         REAL rhs = rightHandSides[i];

         builder.addRowEntries( i + 1, rowlen, rowcols, rowvals );
         builder.setRowLhs( i + 1, lhs );
         builder.setRowRhs( i + 1, rhs );
         builder.setRowLhsInf( i + 1, rowFlags[i].test( RowFlag::kLhsInf ) );
         builder.setRowRhsInf( i + 1, rowFlags[i].test( RowFlag::kRhsInf ) );
         builder.setHardConstraint( i + 1, rowFlags[i].test( RowFlag::kHardConstraint ) );
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
      return builder.build();
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

      for( int i = 0; i < problem.getNRows(); i++ )
      {
         int rowsize = rowSizes[i];

         nnz = nnz + rowsize;
         nrows++;
         if( rowFlags[i].test( RowFlag::kEquation ) ||
             rowFlags[i].test( RowFlag::kLhsInf ) ||
             rowFlags[i].test( RowFlag::kRhsInf ) )
            continue;
         nrows++;
         nnz = nnz + rowsize;
      }

      builder.reserve( nnz, nrows, ncols );

      /* set up rows */
      builder.setNumRows( nrows );
      int counter = 0;
      for( int i = 0; i < problem.getNRows(); ++i )
      {
         auto flags = rowFlags[i];
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
         }
         else if( flags.test( RowFlag::kRhsInf ) )
         {
            assert( !flags.test( RowFlag::kLhsInf ) );

            builder.addRowEntries( counter, rowlen, rowcols, rowvals );
            builder.setRowLhs( counter, lhs );
            builder.setRowRhs( counter, 0 );
            builder.setRowLhsInf( counter, false );
            builder.setRowRhsInf( counter, true );
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
            counter++;

            builder.addRowEntries( counter, rowlen, rowcols, rowvals );
            builder.setRowLhs( counter, lhs );
            builder.setRowRhs( counter, 0 );
            builder.setRowLhsInf( counter, false );
            builder.setRowRhsInf( counter, true );
         }
         counter++;
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

      return builder;
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
