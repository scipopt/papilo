/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*               This file is part of the program and library                */
/*    PaPILO --- Parallel Presolve for Integer and Linear Optimization       */
/*                                                                           */
/* Copyright (C) 2020-2025 Zuse Institute Berlin (ZIB)                       */
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

#include "papilo/presolvers/SimpleProbing.hpp"
#include "papilo/external/catch/catch_amalgamated.hpp"
#include "papilo/core/PresolveMethod.hpp"
#include "papilo/core/Problem.hpp"
#include "papilo/core/ProblemBuilder.hpp"

using namespace papilo;

Problem<double>
setupProblemWithSimpleProbing( );

Problem<double>
setupProblemWithSimpleProbing( double coeff1, double coeff2, double coeff3,
                               int bin_coeff, double value );

TEST_CASE( "simple-probing-trivial-example", "[presolve]" )
{
   Message msg{};
   double time = 0.0;
   int cause = -1;
   Timer t{ time };
   Num<double> num{};
   Problem<double> problem = setupProblemWithSimpleProbing( );
   Statistics statistics{};
   PresolveOptions presolveOptions{};
   PostsolveStorage<double> postsolve =
       PostsolveStorage<double>( problem, num, presolveOptions );
   ProblemUpdate<double> problemUpdate( problem, postsolve, statistics,
                                        presolveOptions, num, msg );
   SimpleProbing<double> presolvingMethod{};
   Reductions<double> reductions{};
   problem.recomputeAllActivities();

   PresolveStatus presolveStatus =
       presolvingMethod.execute( problem, problemUpdate, num, reductions, t, cause );

   REQUIRE( presolveStatus == PresolveStatus::kReduced );
   REQUIRE( reductions.size() == 4 );

   // ub + x (lb - ub) = y => x = 1 -y
   REQUIRE( reductions.getReduction( 0 ).col == 1 );
   REQUIRE( reductions.getReduction( 0 ).row == ColReduction::REPLACE );
   REQUIRE( reductions.getReduction( 0 ).newval == -1 );

   REQUIRE( reductions.getReduction( 1 ).col == 0 );
   REQUIRE( reductions.getReduction( 1 ).row == papilo::ColReduction::NONE );
   REQUIRE( reductions.getReduction( 1 ).newval == 1 );

   // ub + x (lb - ub) = z => x = 1 -z
   REQUIRE( reductions.getReduction( 2 ).col == 2 );
   REQUIRE( reductions.getReduction( 2 ).row == papilo::ColReduction::REPLACE );
   REQUIRE( reductions.getReduction( 2 ).newval == -1 );

   REQUIRE( reductions.getReduction( 3 ).col == 0 );
   REQUIRE( reductions.getReduction( 3 ).row == papilo::ColReduction::NONE );
   REQUIRE( reductions.getReduction( 3 ).newval == 1 );
}


TEST_CASE( "simple-probing-negative-binary-coeff", "[presolve]" )
{
   double time = 0.0;
   int cause = -1;
   Timer t{ time };
   Message msg{};
   Num<double> num{};
   Problem<double> problem =
       setupProblemWithSimpleProbing( 1.0, 1.0, 1.0, -3.0, 0.0 );
   Statistics statistics{};
   PresolveOptions presolveOptions{};
   PostsolveStorage<double> postsolve =
       PostsolveStorage<double>( problem, num, presolveOptions );
   ProblemUpdate<double> problemUpdate( problem, postsolve, statistics,
                                        presolveOptions, num, msg );
   SimpleProbing<double> presolvingMethod{};
   Reductions<double> reductions{};
   problem.recomputeAllActivities();

   PresolveStatus presolveStatus =
       presolvingMethod.execute( problem, problemUpdate, num, reductions, t, cause );

   REQUIRE( presolveStatus == PresolveStatus::kReduced );
   REQUIRE( reductions.size() == 6 );

   for( int i = 0; i < 3; ++i )
   {
      REQUIRE( reductions.getReduction( 2 * i ).col == i );
      REQUIRE( reductions.getReduction( 2 * i ).row == ColReduction::REPLACE );
      REQUIRE( reductions.getReduction( 2 * i ).newval == 1 );

      REQUIRE( reductions.getReduction( 2 * i + 1 ).col == 3 );
      REQUIRE( reductions.getReduction( 2 * i + 1 ).row ==
               papilo::ColReduction::NONE );
      REQUIRE( reductions.getReduction( 2 * i + 1 ).newval == 0 );
   }
}

TEST_CASE( "simple-probing-positive-binary-coeff", "[presolve]" )
{
   Message msg{};
   double time = 0.0;
   int cause = -1;
   Timer t{ time };
   Num<double> num{};
   Problem<double> problem =
       setupProblemWithSimpleProbing( -1.0, -1.0, -1.0, 3, 0.0 );
   Statistics statistics{};
   PresolveOptions presolveOptions{};
   PostsolveStorage<double> postsolve =
       PostsolveStorage<double>( problem, num, presolveOptions );
   ProblemUpdate<double> problemUpdate( problem, postsolve, statistics,
                                        presolveOptions, num, msg );
   SimpleProbing<double> presolvingMethod{};
   Reductions<double> reductions{};
   problem.recomputeAllActivities();

   PresolveStatus presolveStatus =
       presolvingMethod.execute( problem, problemUpdate, num, reductions, t, cause );

   REQUIRE( presolveStatus == PresolveStatus::kReduced );
   REQUIRE( reductions.size() == 6 );

   for( int i = 0; i < 3; ++i )
   {
      REQUIRE( reductions.getReduction( 2 * i ).col == i );
      REQUIRE( reductions.getReduction( 2 * i ).row == ColReduction::REPLACE );
      REQUIRE( reductions.getReduction( 2 * i ).newval == 1 );

      REQUIRE( reductions.getReduction( 2 * i + 1 ).col == 3 );
      REQUIRE( reductions.getReduction( 2 * i + 1 ).row ==
               papilo::ColReduction::NONE );
      REQUIRE( reductions.getReduction( 2 * i + 1 ).newval == 0 );
   }
}

TEST_CASE( "simple-different-int-coeff", "[presolve]" )
{
   Message msg{};
   double time = 0.0;
   int cause = -1;
   Timer t{ time };
   Num<double> num{};
   Problem<double> problem =
       setupProblemWithSimpleProbing( -2.0, -2.0, 5, 1, 1.0 );
   Statistics statistics{};
   PresolveOptions presolveOptions{};
   PostsolveStorage<double> postsolve =
       PostsolveStorage<double>( problem, num, presolveOptions );
   ProblemUpdate<double> problemUpdate( problem, postsolve, statistics,
                                        presolveOptions, num, msg );
   SimpleProbing<double> presolvingMethod{};
   Reductions<double> reductions{};
   problem.recomputeAllActivities();

   PresolveStatus presolveStatus =
       presolvingMethod.execute( problem, problemUpdate, num, reductions, t, cause );

   REQUIRE( presolveStatus == PresolveStatus::kReduced );
   REQUIRE( reductions.size() == 6 );

   Vec<double> factors = { 1, 1, -1 };
   Vec<double> offset = { 0, 0, 1 };
   Vec<double> columns = { 0, 1, 3 };

   for( int i = 0; i < 3; ++i )
   {
      REQUIRE( reductions.getReduction( 2 * i ).col == columns[i] );
      REQUIRE( reductions.getReduction( 2 * i ).row == ColReduction::REPLACE );
      REQUIRE( reductions.getReduction( 2 * i ).newval == factors[i] );

      REQUIRE( reductions.getReduction( 2 * i + 1 ).col == 2 );
      REQUIRE( reductions.getReduction( 2 * i + 1 ).row ==
               papilo::ColReduction::NONE );
      REQUIRE( reductions.getReduction( 2 * i + 1 ).newval == offset[i] );
   }
}

TEST_CASE( "simple-different-int-coeff-pt2", "[presolve]" )
{
   Message msg{};
   double time = 0.0;
   int cause = -1;
   Timer t{ time };
   Num<double> num{};
   Problem<double> problem =
       setupProblemWithSimpleProbing( -2.0, -2.0, 1, 5, 1.0 );
   Statistics statistics{};
   PresolveOptions presolveOptions{};
   PostsolveStorage<double> postsolve =
       PostsolveStorage<double>( problem, num, presolveOptions );
   ProblemUpdate<double> problemUpdate( problem, postsolve, statistics,
                                        presolveOptions, num, msg );
   SimpleProbing<double> presolvingMethod{};
   Reductions<double> reductions{};
   problem.recomputeAllActivities();

   PresolveStatus presolveStatus =
       presolvingMethod.execute( problem, problemUpdate, num, reductions, t, cause );

   REQUIRE( presolveStatus == PresolveStatus::kReduced );
   REQUIRE( reductions.size() == 6 );

   Vec<double> factors = { 1, 1, -1 };
   Vec<double> offset = { 0, 0, 1 };
   Vec<double> columns = { 0, 1, 2 };

   for( int i = 0; i < 3; ++i )
   {
      REQUIRE( reductions.getReduction( 2 * i ).col == columns[i] );
      REQUIRE( reductions.getReduction( 2 * i ).row == ColReduction::REPLACE );
      REQUIRE( reductions.getReduction( 2 * i ).newval == factors[i] );

      REQUIRE( reductions.getReduction( 2 * i + 1 ).col == 3 );
      REQUIRE( reductions.getReduction( 2 * i + 1 ).row ==
               papilo::ColReduction::NONE );
      REQUIRE( reductions.getReduction( 2 * i + 1 ).newval == offset[i] );
   }
}

Problem<double>
setupProblemWithSimpleProbing ()
{
   // simple probing requires
   // - rhs = (sup -inf) / 2
   // futhermore for one column
   // - integral variables
   // - coeff = supp - rhs
   // i.e. 2x + y + z = 2 with (sup = 4 & x = binary)
   // -> ub + x (lb - ub) = y/z
   double coeff =  1.0 ;
   double bin_coeff = 2.0 ;
   Num<double> num{};
   Vec<double> coefficients{ 3.0, 1.0, 1.0, 1.0 };
   Vec<double> upperBounds{ 1.0, 1.0, 1.0, 1.0 };
   Vec<double> lowerBounds{ 0.0, 0.0, 0.0, 0.0 };
   Vec<uint8_t> isIntegral{ 1, 1, 1, 1 };

   Vec<double> rhs{ 2.0, 2.0 };
   Vec<double> lhs{ rhs[0], 3.0 };
   Vec<std::string> rowNames{ "A1", "A2" };
   Vec<std::string> columnNames{ "c1", "c2", "c3", "c4" };
   Vec<std::tuple<int, int, double>> entries{
       std::tuple<int, int, double>{ 0, 0, bin_coeff },
       std::tuple<int, int, double>{ 0, 1, coeff },
       std::tuple<int, int, double>{ 0, 2, coeff},
       std::tuple<int, int, double>{ 1, 1, 2.0 } };

   ProblemBuilder<double> pb;
   pb.reserve( (int)entries.size(), (int)rowNames.size(),
               (int)columnNames.size() );
   pb.setNumRows( (int)rowNames.size() );
   pb.setNumCols( (int)columnNames.size() );
   pb.setColUbAll( upperBounds );
   pb.setColLbAll( lowerBounds );
   pb.setObjAll( coefficients );
   pb.setObjOffset( 0.0 );
   pb.setColIntegralAll( isIntegral );
   pb.setRowRhsAll( rhs );
   pb.setRowLhsAll( lhs );
   pb.addEntryAll( entries );
   pb.setColNameAll( columnNames );
   pb.setProblemName( "matrix for testing simple probing" );
   Problem<double> problem = pb.build();
   problem.getConstraintMatrix().modifyLeftHandSide( 0, num, lhs[0] );
   return problem;
}

Problem<double>
setupProblemWithSimpleProbing( double coeff1, double coeff2,
                                                double coeff3, int bin_coeff, double value )
{
   Num<double> num{};
   Vec<double> coefficients{ 1.0, 1.0, 1.0, 1.0 };
   Vec<double> upperBounds{ 1.0, 1.0, 1.0, 1.0 };
   Vec<double> lowerBounds{ 0.0, 0.0, 0.0, 0.0 };
   Vec<uint8_t> isIntegral{ 1, 1, 1, 1 };

   Vec<double> rhs{ value };
   Vec<double> lhs{ rhs[0] };
   Vec<std::string> rowNames{ "A1" };
   Vec<std::string> columnNames{ "x1", "x2", "x3", "y" };
   Vec<std::tuple<int, int, double>> entries{
       std::tuple<int, int, double>{ 0, 0, coeff1 },
       std::tuple<int, int, double>{ 0, 1, coeff2 },
       std::tuple<int, int, double>{ 0, 2, coeff3 },
       std::tuple<int, int, double>{ 0, 3, bin_coeff } };

   ProblemBuilder<double> pb;
   pb.reserve( (int)entries.size(), (int)rowNames.size(),
               (int)columnNames.size() );
   pb.setNumRows( (int)rowNames.size() );
   pb.setNumCols( (int)columnNames.size() );
   pb.setColUbAll( upperBounds );
   pb.setColLbAll( lowerBounds );
   pb.setObjAll( coefficients );
   pb.setObjOffset( 0.0 );
   pb.setColIntegralAll( isIntegral );
   pb.setRowRhsAll( rhs );
   pb.setRowLhsAll( lhs );
   pb.addEntryAll( entries );
   pb.setColNameAll( columnNames );
   pb.setProblemName( "matrix for testing simple probing" );
   Problem<double> problem = pb.build();
   problem.getConstraintMatrix().modifyLeftHandSide( 0, num, lhs[0] );
   return problem;
}
