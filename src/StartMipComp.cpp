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
#include "fix/VolumeAlgorithm.hpp"
#include "papilo/core/Problem.hpp"
#include "papilo/core/ProblemBuilder.hpp"
#include "papilo/io/MpsParser.hpp"
#include "papilo/misc/OptionsParser.hpp"
#include <boost/program_options.hpp>
#include <fstream>

using namespace papilo;

Problem<double>
modify_problem( Problem<double>& problem );

void
invert( const double* pDouble, double* result, int length );

int
main( int argc, char* argv[] )
{

   // get the options passed by the user
   OptionsInfo optionsInfo;
   try
   {
      optionsInfo = parseOptions( argc, argv );
   }
   catch( const boost::program_options::error& ex )
   {
      std::cerr << "Error while parsing the options.\n" << '\n';
      std::cerr << ex.what() << '\n';
      return 1;
   }

   if( !optionsInfo.is_complete )
      return 0;

   double readtime = 0;
   Problem<double> problem;
   Num<double> num{};
   Message msg{};
   boost::optional<Problem<double>> prob;

   {
      Timer t( readtime );
      prob = MpsParser<double>::loadProblem( optionsInfo.instance_file );
   }

   // Check whether reading was successful or not
   if( !prob )
   {
      fmt::print( "error loading problem {}\n", optionsInfo.instance_file );
      return 0;
   }
   problem = *prob;

   fmt::print( "reading took {:.3} seconds\n", readtime );

   // set up ProblemUpdate to trivialPresolve so that activities exist
   Presolve<double> presolve{};
   auto result = presolve.apply( problem, false );

   switch( result.status )
   {
   case papilo::PresolveStatus::kUnbounded:
   case papilo::PresolveStatus::kUnbndOrInfeas:
   case papilo::PresolveStatus::kInfeasible:
      fmt::print( "PaPILO detected infeasibility or unbounded-ness\n" );
      return 0;
   case papilo::PresolveStatus::kUnchanged:
   case papilo::PresolveStatus::kReduced:
      break;
   }


   Vec<double> primal_heur_sol{};
   primal_heur_sol.reserve( problem.getNCols() );
   {

      VolumeAlgorithm<double> algorithm{ {},  {},   0.5,  0.1,  1, 0.0005, 2,
                                         1.1, 0.66, 0.02, 0.01, 2, 20 };

      Problem<double> reformulated = modify_problem( problem );

      // TODO: add same small heuristic
      // generate pi
      Vec<double> pi{};
      for( int i = 0; i < reformulated.getNRows(); i++ )
         pi.push_back( 0 );

      // generate UB
      StableSum<double> min_value{};

      for( int i = 0; i < problem.getNCols(); i++ )
      {
         if( num.isZero( problem.getObjective().coefficients[i] ) )
            continue;
         else if( num.isLT( problem.getObjective().coefficients[i], 0 ) )
         {
            if( problem.getColFlags()[i].test( ColFlag::kLbInf ) )
            {
               fmt::print( "Could not calculate objective bound: variable {} is unbounded",
                           i );
               return 1;
            }
            min_value.add( problem.getObjective().coefficients[i] +
                           problem.getLowerBounds()[i] );
         }
         else
         {
            if( problem.getColFlags()[i].test( ColFlag::kUbInf ) )
            {
               fmt::print( "Could not calculate objective bound: variable {} is unbounded",
                           i );
               return 1;
            }
            min_value.add( problem.getObjective().coefficients[i] +
                           problem.getUpperBounds()[i] );
         }
      }

      primal_heur_sol = algorithm.volume_algorithm(
          reformulated.getObjective().coefficients,
          reformulated.getConstraintMatrix(),
          reformulated.getConstraintMatrix().getLeftHandSides(),
          reformulated.getVariableDomains(), pi, min_value.get() );
   }
   msg.info( "Primal heuristic solution:\n" );
   for( int i = 0; i < problem.getNCols(); i++ )
      msg.info( "   x[{}] = {}\n", i, primal_heur_sol[i] );

   // we continue with the presolved problem since the columns are the same
   problem.recomputeAllActivities();

   ProbingView<double> probing_view{ problem, num };
   FixAndPropagate<double> fixAndPropagate{ msg, num };
   Solution<double> sol = fixAndPropagate.fix_and_propagate(
       problem, probing_view, primal_heur_sol );

   Solution<double> original_solution{};
   Solution<double> reduced_solution{ primal_heur_sol };
   primal_heur_sol.reserve( problem.getNCols() );
   Postsolve<double> postsolve{ msg, num };
   msg.info( "\n" );

   postsolve.undo( sol, original_solution, result.postsolve );

   msg.info( "Solution\n" );
   for( int i = 0; i < original_solution.primal.size(); i++ )
      msg.info( "   x[{}] = {}\n", i, original_solution.primal[i] );

   return 0;
}

Problem<double>
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
      if( flags.test( RowFlag::kEquation ) || flags.test( RowFlag::kLhsInf ) ||
          flags.test( RowFlag::kRhsInf ) )
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
         double neg_rowvals[rowlen];
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
         double neg_rowvals[rowlen];
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
   return builder.build();
}

void
invert( const double* pDouble, double* result, int length )
{
   for( int i = 0; i < length; i++ )
      result[i] = pDouble[i] * -1;
}
