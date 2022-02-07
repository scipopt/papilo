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
#ifdef PAPILO_TBB
      //#TODO:set to 0
      tbb::task_arena arena( 1 );
#endif

#ifdef PAPILO_TBB
      return arena.execute(
          [this, &problem, &parameter]()
          {
#endif
             // set up ProblemUpdate to trivialPresolve so that activities exist
             Presolve<REAL> presolve{};
             auto result = presolve.apply( problem, false );

             // TODO: add a check if presolve solved to optimality
             switch( result.status )
             {
             case papilo::PresolveStatus::kUnbounded:
             case papilo::PresolveStatus::kUnbndOrInfeas:
             case papilo::PresolveStatus::kInfeasible:
                fmt::print(
                    "PaPILO detected infeasibility or unbounded-ness\n" );
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

             Vec<REAL> primal_heur_sol{};
             primal_heur_sol.reserve( problem.getNCols() );

             REAL best_obj_value{};
             bool initialized = false;
             Vec<REAL> best_solution{};
             best_solution.reserve( problem.getNCols() );

             ProblemBuilder<REAL> builder = modify_problem( problem );
             Problem<REAL> reformulated = builder.build();

             // TODO: add same small heuristic
             Vec<REAL> pi;
             pi.reserve( reformulated.getNRows() );
             generate_initial_dual_solution( reformulated, pi );

             REAL min_val = calc_upper_bound_for_objective( problem );
             if( min_val == std::numeric_limits<double>::min() )
                return;

             VolumeAlgorithm<REAL> algorithm{ msg, num, timer, parameter };
             ConflictAnalysis<REAL> conflict_analysis{ msg, num, timer };

             problem.recomputeAllActivities();

             Vec<RoundingStrategy<REAL>*> strategies{};
             Vec<Vec<REAL>> int_solutions{};
             Vec<ProbingView<REAL>> views{};
             Vec<REAL> obj_value;
             Vec<bool> infeasible_arr;
             setup( problem, strategies, int_solutions, views, obj_value,
                    infeasible_arr );
             Vec<int> cols_sorted_by_obj{};
             sort( problem.getObjective().coefficients, cols_sorted_by_obj );
             while( true )
             {
                if( timer.getTime() >= parameter.time_limit )
                   break;
                msg.info( "Starting volume algorithm\n" );
                primal_heur_sol = algorithm.volume_algorithm(
                    reformulated.getObjective().coefficients,
                    reformulated.getConstraintMatrix(),
                    reformulated.getConstraintMatrix().getLeftHandSides(),
                    reformulated.getVariableDomains(), pi, min_val );
                print_solution( primal_heur_sol );

                msg.info( "Starting fixing and propagating\n" );

                if( timer.getTime() >= parameter.time_limit )
                   break;

                perform_fix_and_propagate( primal_heur_sol, strategies,
                                           int_solutions, views, obj_value,
                                           infeasible_arr );

                // TODO: consider non TBB version
                bool feasible = !infeasible_arr[0]
#ifdef PAPILO_TBB
                                || !infeasible_arr[1] || !infeasible_arr[2] ||
                                !infeasible_arr[3];
#else
             ;
#endif

                // TODO: copy the best solution;
                if( feasible )
                {
                   perform_one_opt( problem, int_solutions, views,
                                    infeasible_arr, cols_sorted_by_obj, obj_value );
                   int best_index = -1;
#ifdef PAPILO_TBB
                   for( int i = 0; i < 4; i++ )
#else
            for( int i = 0; i < 1; i++ )
#endif
                   {
                      if( !infeasible_arr[i] &&
                          ( num.isLT( obj_value[i], best_obj_value ) ||
                            !initialized ) )
                      {
                         initialized = true;
                         best_index = i;
                         best_obj_value = obj_value[i];
                      }
                   }
                   assert( best_index != -1 );
                   best_solution = int_solutions[best_index];
                   break;
                }

                if( timer.getTime() >= parameter.time_limit )
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

             postsolve.undo( reduced_solution, original_solution,
                             result.postsolve );

             print_solution( original_solution.primal );
             msg.info( "Solving took {} seconds.\n", timer.getTime() );
#ifdef PAPILO_TBB
          } );
#endif
   }

   void
   sort( const Vec<REAL>& objective, Vec<int>& cols_sorted_by_obj ) const
   {
      cols_sorted_by_obj.reserve( objective.size() );
      for( int i = 0; i < objective.size(); i++ )
         cols_sorted_by_obj.push_back( i );
      pdqsort( cols_sorted_by_obj.begin(), cols_sorted_by_obj.end(),
               [&]( const int a, const int b )
               {
                  return objective[a] > objective[b] ||
                         ( objective[a] == objective[b] && a > b );
               } );
   }

   void
   setup( const Problem<REAL>& problem,
          Vec<RoundingStrategy<REAL>*>& strategies,
          Vec<Vec<REAL>>& int_solutions, Vec<ProbingView<REAL>>& views,
          Vec<REAL>& obj_value, Vec<bool>& infeasible_arr )
   {
#ifdef PAPILO_TBB
      auto s1 = new FarkasRoundingStrategy<REAL>{ 0, num, false };
      auto s2 = new FarkasRoundingStrategy<REAL>{ 0, num, true };
      auto s3 = new FractionalRoundingStrategy<REAL>{ num };
      auto s4 = new RandomRoundingStrategy<REAL>{ 0, num };
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

   void
   perform_fix_and_propagate( const Vec<REAL>& primal_heur_sol,
                              Vec<RoundingStrategy<REAL>*>& strategies,
                              Vec<Vec<REAL>>& int_solutions,
                              Vec<ProbingView<REAL>>& views,
                              Vec<REAL>& obj_value,
                              Vec<bool>& infeasible_arr ) const
   {

#ifdef PAPILO_TBB
      for( auto view : views )
         view.reset();
      tbb::parallel_for(
          tbb::blocked_range<int>( 0, 4 ),
          [&]( const tbb::blocked_range<int>& r )
          {
             for( int i = r.begin(); i != r.end(); ++i )
             {
                // TODO: extract
                FixAndPropagate<REAL> fixAndPropagate{ msg, num, true };
                infeasible_arr[i] = fixAndPropagate.fix_and_propagate(
                    primal_heur_sol, int_solutions[i], *( strategies[i] ),
                    views[i] );
                if( infeasible_arr[i] )
                {
                   obj_value[i] = 0;
                   break;
                }
                obj_value[i] =
                    calculate_obj_value( int_solutions[i], views[i] );
                msg.info( "Diving {} found obj value {}!\n", i, obj_value[i] );
             }
          } );
#else
      for( auto view : views )
         view.reset();
      FixAndPropagate<REAL> fixAndPropagate{ msg, num, true };
      infeasible_arr[0] = fixAndPropagate.fix_and_propagate(
          primal_heur_sol, int_solutions[0], *( strategies[0] ), views[0] );
      if( infeasible_arr[0] )
      {
         obj_value[0] = 0;
         return;
      }
      StableSum<REAL> sum{};
      for( int j = 0; j < primal_heur_sol.size(); j++ )
         sum.add( int_solutions[0][j] * views[0].get_obj()[j] );
      obj_value[0] = sum.get();
      msg.info( "Diving {} found obj value {}!\n", 0, obj_value[0] );
#endif
   }
   REAL
   calculate_obj_value( const Vec<REAL>& int_solution,
                        const ProbingView<REAL>& view ) const
   {
      StableSum<REAL> sum{};
      for( int j = 0; j < int_solution.size(); j++ )
         sum.add( int_solution[j] * view.get_obj()[j] );
      return sum.get();
   }

   void
   perform_one_opt( const Problem<REAL>& problem, Vec<Vec<REAL>>& int_solutions,
                    Vec<ProbingView<REAL>>& views, Vec<bool>& infeasible_arr,
                    Vec<int>& cols_sorted_by_obj, Vec<REAL>& obj_value ) const
   {
      FixAndPropagate<REAL> fixAndPropagate{ msg, num, false };

#ifdef PAPILO_TBB
      for( auto view : views )
         view.reset();
      Vec<REAL> coefficients = problem.getObjective().coefficients;
      tbb::parallel_for(
          tbb::blocked_range<int>( 0, 4 ),
          [&]( const tbb::blocked_range<int>& r )
          {
             for( int i = r.begin(); i != r.end(); ++i )
             {
                Vec<REAL> result = { int_solutions[i] };
                for( int j = 0; j < cols_sorted_by_obj.size(); j++ )
                {
                   views[i].reset();
                   if( num.isZero( coefficients[j] ) )
                      break;
                   if( !problem.getColFlags()[j].test( ColFlag::kIntegral ) ||
                       problem.getLowerBounds()[j] != 0 ||
                       problem.getUpperBounds()[j] != 1 )
                      continue;
                   REAL solution_value = int_solutions[i][j];
                   if( num.isGT( coefficients[j], 0 ) )
                   {
                      if( num.isZero( solution_value ) )
                         continue;
                      // TODO: ignore the conflicts at this place?
                      bool infeasible = fixAndPropagate.one_opt(
                          int_solutions[i], j, 0, views[i], result );
                      REAL value = calculate_obj_value( result, views[i] );
                      msg.info( "OneOpt flipping variable {}: ", i);
                      if( infeasible )
                      {
                         msg.info( "infeasible: \n" );
                         continue;
                      }
                      if( num.isGE( value, obj_value[i] ) )
                         msg.info( "unsuccessful -> worse obj {}: \n", value );
                      else if( num.isGT( value, obj_value[i] ) )
                      {
                         msg.info( "successful -> better obj {}: \n", value );
                         int_solutions[i] = result;
                         obj_value[i] = value;
                      }
                   }
                   else
                   {
                      assert( num.isLT( coefficients[j], 0 ) );
                      if( !num.isZero( solution_value ) )
                         continue;
                      bool infeasible = fixAndPropagate.one_opt(
                          int_solutions[i], j, 1, views[i], result );
                      REAL value = calculate_obj_value( result, views[i] );
                      msg.info( "OneOpt flipping variable {}: ", i);
                      if( infeasible )
                      {
                         msg.info( "infeasible!\n" );
                         continue;
                      }
                      if( num.isGE( value, obj_value[i] ) )
                         msg.info( "unsuccessful -> worse obj {}!\n", value );
                      else if( num.isGT( value, obj_value[i] ) )
                      {
                         msg.info( "successful -> better obj {}!\n", value );
                         int_solutions[i] = result;
                         obj_value[i] = value;
                      }
                      break;
                   }
                }
             }
          } );

#endif
   }

   REAL
   calc_upper_bound_for_objective( const Problem<REAL>& problem ) const
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

      for( int i = 0; i < problem.getNRows(); i++ )
      {
         nrows++;
         int rowsize = rowSizes[i];
         nnz = nnz + rowsize;
         auto flags = rowFlags[i];
         if( flags.test( RowFlag::kEquation ) ||
             flags.test( RowFlag::kLhsInf ) || flags.test( RowFlag::kRhsInf ) )
            continue;
         nrows++;
         nnz = nnz + rowsize;
      }

      builder.reserve( nnz, nrows, ncols );

      /* set up columns */
      builder.setNumCols( ncols );
      for( int i = 0; i != ncols; ++i )
      {
         builder.setColLb( i, problem.getLowerBounds()[i] );
         builder.setColUb( i, problem.getUpperBounds()[i] );
         auto flags = colFlags[i];
         builder.setColLbInf( i, flags.test( ColFlag::kLbInf ) );
         builder.setColUbInf( i, flags.test( ColFlag::kUbInf ) );

         builder.setColIntegral( i, flags.test( ColFlag::kIntegral ) );
         builder.setObj( i, coefficients[i] );
      }

      /* set up rows */
      builder.setNumRows( nrows );
      int counter = 0;
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
            REAL* neg_rowvals = new REAL[rowlen];
            invert( rowvals, neg_rowvals, rowlen );
            builder.addRowEntries( counter, rowlen, rowcols, neg_rowvals );
            builder.setRowLhs( counter, -rhs );
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
            REAL* neg_rowvals = new REAL[rowlen];
            invert( rowvals, neg_rowvals, rowlen );

            builder.addRowEntries( counter, rowlen, rowcols, neg_rowvals );
            builder.setRowLhs( counter, -rhs );
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
