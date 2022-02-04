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
#include "fix/strategy/FarkasRoundingStrategy.hpp"
#include "fix/FixAndPropagate.hpp"
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

//TODO: find a funny name for it
template <typename REAL>
class Algorithm
{
   Message msg;
   Num<REAL> num;

 public:
   Algorithm( Message _msg, Num<REAL> _num ) : msg( _msg ), num( _num ) {}

   void
   solve_problem( Problem<REAL>& problem )
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

      Vec<REAL> primal_heur_sol{};
      primal_heur_sol.reserve( problem.getNCols() );

      Vec<REAL> int_solution{};
      int_solution.resize( problem.getNCols() );


      ProblemBuilder<REAL> builder = modify_problem( problem );
      Problem<REAL> reformulated = builder.build();

      // TODO: add same small heuristic
      Vec<REAL> pi;
      pi.reserve( reformulated.getNRows() );
      generate_initial_dual_solution( reformulated, pi );

      double min_val = calc_upper_bound_for_objective( problem );
      if( min_val == std::numeric_limits<double>::min() )
         return;

      // TODO: extract parameters
      VolumeAlgorithm<double> algorithm{ {},     {},    0.05, 0.1, 0.2,
                                         0.0005, 2,     2,    1.1, 0.66,
                                         0.01,   0.001, 0.02, 2,   20 };
      ConflictAnalysis<REAL> conflict_analysis{ msg, num };

      problem.recomputeAllActivities();

      while( true )
      {
         msg.info( "Starting volume algorithm\n" );
         primal_heur_sol = algorithm.volume_algorithm(
             reformulated.getObjective().coefficients,
             reformulated.getConstraintMatrix(),
             reformulated.getConstraintMatrix().getLeftHandSides(),
             reformulated.getVariableDomains(), pi, min_val );
         print_solution( primal_heur_sol );

         msg.info( "Starting fixing and propagating\n" );

         ProbingView<REAL> probing_view{ problem, num };
         FixAndPropagate<REAL> fixAndPropagate{ msg, num, probing_view, true };
         FarkasRoundingStrategy<REAL> strategy{ 0, {}, false };
         bool infeasible =
             fixAndPropagate.fix_and_propagate( primal_heur_sol, int_solution, strategy );
         if( !infeasible )
            break;

         msg.info( "Starting conflict analysis\n" );
         bool abort = conflict_analysis.perform_conflict_analysis();
         if( abort )
            return;
         //TODO: add constraint to builder and generate new problem
      }

      Solution<REAL> original_solution{};
      Solution<REAL> reduced_solution{ int_solution };
      Postsolve<REAL> postsolve{ msg, num };

      postsolve.undo( reduced_solution, original_solution, result.postsolve );

      print_solution( original_solution.primal );
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
   print_solution( const Vec<double>& sol )
   {
      msg.debug( "Primal solution:\n" );
      for( int i = 0; i < sol.size(); i++ )
         msg.debug( "   x[{}] = {}\n", i, sol[i] );
   }

   void
   generate_initial_dual_solution( const Problem<double>& problem,
                                   Vec<double>& dual_solution )
   {
      for( int i = 0; i < problem.getNRows(); i++ )
         dual_solution.push_back( 0 );
   }

   ProblemBuilder<double>
   modify_problem( Problem<double>& problem )
   {
      ProblemBuilder<double> builder;

      int nnz = 0;
      int ncols = problem.getNCols();
      int nrows = 0;
      ConstraintMatrix<double>& matrix = problem.getConstraintMatrix();
      Vec<ColFlags>& colFlags = problem.getColFlags();
      Vec<RowFlags>& rowFlags = matrix.getRowFlags();
      Vec<int>& rowSizes = problem.getRowSizes();
      Vec<double>& coefficients = problem.getObjective().coefficients;
      Vec<double>& leftHandSides = matrix.getLeftHandSides();
      Vec<double>& rightHandSides = matrix.getRightHandSides();

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
         const SparseVectorView<double>& view = matrix.getRowCoefficients( i );
         const int* rowcols = view.getIndices();
         const double* rowvals = view.getValues();
         int rowlen = view.getLength();
         auto flags = rowFlags[i];
         double lhs = leftHandSides[i];
         double rhs = rightHandSides[i];

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
            double* neg_rowvals = new double[rowlen];
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
            double* neg_rowvals = new double[rowlen];
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
   invert( const double* pDouble, double* result, int length )
   {
      for( int i = 0; i < length; i++ )
         result[i] = pDouble[i] * -1;
   }
};

} // namespace papilo
