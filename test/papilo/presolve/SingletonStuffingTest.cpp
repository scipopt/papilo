/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*               This file is part of the program and library                */
/*    PaPILO --- Parallel Presolve for Integer and Linear Optimization       */
/*                                                                           */
/* Copyright (C) 2020-2026 Zuse Institute Berlin (ZIB)                       */
/*                                                                           */
/* Licensed under the Apache License, Version 2.0 (the "License");           */
/* you may not use this file except in compliance with the License.          */
/* You may obtain a copy of the License at                                   */
/*                                                                           */
/*     http://www.apache.org/licenses/LICENSE-2.0                            */
/*                                                                           */
/* Unless required by applicable law or agreed to in writing, software       */
/* distributed under the License is distributed on an "AS IS" BASIS,         */
/* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  */
/* See the License for the specific language governing permissions and       */
/* limitations under the License.                                            */
/*                                                                           */
/* You should have received a copy of the Apache-2.0 license                 */
/* along with PaPILO; see the file LICENSE. If not visit scipopt.org.        */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include "papilo/presolvers/SingletonStuffing.hpp"
#include "papilo/external/catch/catch_amalgamated.hpp"
#include "papilo/core/PresolveMethod.hpp"
#include "papilo/core/Problem.hpp"
#include "papilo/core/ProblemBuilder.hpp"

using namespace papilo;

Problem<double>
setupProblemWithSingletonStuffingColumn();

Problem<double>
setupProblemWithIntegralSingletonStuffingColumn( const Vec<double>& rowValues,
                                                 double side );

void
forceCalculationOfSingletonStuffingRows( Problem<double>& problem,
                                 ProblemUpdate<double>& problemUpdate )
{
   problem.recomputeLocks();
   problemUpdate.trivialColumnPresolve();
   problem.recomputeAllActivities();
}

TEST_CASE( "singleton-stuffing-make-sure-to-first-set-bounds-to-infinity", "[presolve]" )
{
   Num<double> num{};
   Message msg{};
   double time = 0.0;
   int cause = -1;
   Timer t{ time };
   Problem<double> problem = setupProblemWithSingletonStuffingColumn();
   Statistics statistics{};
   PresolveOptions presolveOptions{};
   auto postsolve = PostsolveStorage<double>( problem, num, presolveOptions );
   ProblemUpdate<double> problemUpdate( problem, postsolve, statistics,
                                        presolveOptions, num, msg );
   presolveOptions.dualreds = 0;
   forceCalculationOfSingletonStuffingRows( problem, problemUpdate );
   SingletonStuffing<double> presolvingMethod{};
   Reductions<double> reductions{};
   presolveOptions.dualreds = 2;


   PresolveStatus presolveStatus =
       presolvingMethod.execute( problem, problemUpdate, num, reductions, t, cause );

   REQUIRE( presolveStatus == PresolveStatus::kReduced );
   REQUIRE( reductions.size() == 7 );
   REQUIRE( reductions.getReduction( 0 ).col == 3 );
   REQUIRE( reductions.getReduction( 0 ).row == ColReduction::BOUNDS_LOCKED );
   REQUIRE( reductions.getReduction( 0 ).newval == 0 );

   REQUIRE( reductions.getReduction( 1 ).col == RowReduction::LOCKED );
   REQUIRE( reductions.getReduction( 1 ).row == 1 );
   REQUIRE( reductions.getReduction( 1 ).newval == 0 );

   REQUIRE( reductions.getReduction( 2 ).col == RowReduction::RHS );
   REQUIRE( reductions.getReduction( 2 ).row == 1 );
   REQUIRE( reductions.getReduction( 2 ).newval == 0 );

   REQUIRE( reductions.getReduction( 3 ).col == 3 );
   REQUIRE( reductions.getReduction( 3 ).row == ColReduction::SUBSTITUTE_OBJ );
   REQUIRE( reductions.getReduction( 3 ).newval == 1 );

   REQUIRE( reductions.getReduction( 4 ).col == 3 );
   REQUIRE( reductions.getReduction( 4 ).row == 1 );
   REQUIRE( reductions.getReduction( 4 ).newval == 0 );

   REQUIRE( reductions.getReduction( 5 ).col == RowReduction::RHS_INF );
   REQUIRE( reductions.getReduction( 5 ).row == 1 );
   REQUIRE( reductions.getReduction( 5 ).newval == 0 );

   REQUIRE( reductions.getReduction( 6 ).col == RowReduction::LHS );
   REQUIRE( reductions.getReduction( 6 ).row == 1 );
   REQUIRE( reductions.getReduction( 6 ).newval == 1.98 );
}

TEST_CASE( "singleton-stuffing-ensure-integrality", "[presolve]" )
{
   Num<double> num{};
   Message msg{};
   double time = 0.0;
   int cause = -1;
   Timer t{ time };
   Problem<double> problem = setupProblemWithIntegralSingletonStuffingColumn(
       { 2.0, 4.0, 8.0, 16.0 }, 1.0 );
   Statistics statistics{};
   PresolveOptions presolveOptions{};
   auto postsolve = PostsolveStorage<double>( problem, num, presolveOptions );
   ProblemUpdate<double> problemUpdate( problem, postsolve, statistics,
                                        presolveOptions, num, msg );
   presolveOptions.dualreds = 0;
   forceCalculationOfSingletonStuffingRows( problem, problemUpdate );
   SingletonStuffing<double> presolvingMethod{};
   Reductions<double> reductions{};
   presolveOptions.dualreds = 2;

   PresolveStatus presolveStatus =
       presolvingMethod.execute( problem, problemUpdate, num, reductions, t, cause );

   // 2x1 + 4x2 + 8x3 + 16x4 = 1 with integral x: the left-hand side is always
   // even, hence the equation cannot be satisfied
   REQUIRE( presolveStatus == PresolveStatus::kInfeasible );
}

TEST_CASE( "singleton-stuffing-ensure-integrality-not-divisible", "[presolve]" )
{
   Num<double> num{};
   Message msg{};
   double time = 0.0;
   int cause = -1;
   Timer t{ time };
   // 2x1 + 3x2 + 4x3 + 5x4 = 3 is feasible (x2 = 1), but no column can be
   // substituted without losing integrality since for every column some other
   // coefficient is not divisible by its coefficient
   Problem<double> problem = setupProblemWithIntegralSingletonStuffingColumn(
       { 2.0, 3.0, 4.0, 5.0 }, 3.0 );
   Statistics statistics{};
   PresolveOptions presolveOptions{};
   auto postsolve = PostsolveStorage<double>( problem, num, presolveOptions );
   ProblemUpdate<double> problemUpdate( problem, postsolve, statistics,
                                        presolveOptions, num, msg );
   presolveOptions.dualreds = 0;
   forceCalculationOfSingletonStuffingRows( problem, problemUpdate );
   SingletonStuffing<double> presolvingMethod{};
   Reductions<double> reductions{};
   presolveOptions.dualreds = 2;

   PresolveStatus presolveStatus =
       presolvingMethod.execute( problem, problemUpdate, num, reductions, t, cause );

   REQUIRE( presolveStatus == PresolveStatus::kUnchanged );
   REQUIRE( reductions.size() == 0 );
}

TEST_CASE( "singleton-stuffing-ensure-integrality-side-not-exactly-integral", "[presolve]" )
{
   Num<double> num{};
   Message msg{};
   double time = 0.0;
   int cause = -1;
   Timer t{ time };
   // the scaled side of c1 is integral w.r.t. the feasibility tolerance but not
   // w.r.t. epsilon: no infeasibility, but c1 must not be substituted either
   Problem<double> problem = setupProblemWithIntegralSingletonStuffingColumn(
       { 2.0, 4.0, 8.0, 16.0 }, 2.0 + 1e-7 );
   Statistics statistics{};
   PresolveOptions presolveOptions{};
   auto postsolve = PostsolveStorage<double>( problem, num, presolveOptions );
   ProblemUpdate<double> problemUpdate( problem, postsolve, statistics,
                                        presolveOptions, num, msg );
   presolveOptions.dualreds = 0;
   forceCalculationOfSingletonStuffingRows( problem, problemUpdate );
   SingletonStuffing<double> presolvingMethod{};
   Reductions<double> reductions{};
   presolveOptions.dualreds = 2;

   PresolveStatus presolveStatus =
       presolvingMethod.execute( problem, problemUpdate, num, reductions, t, cause );

   REQUIRE( presolveStatus == PresolveStatus::kUnchanged );
   REQUIRE( reductions.size() == 0 );
}

Problem<double>
setupProblemWithIntegralSingletonStuffingColumn( const Vec<double>& rowValues,
                                                 double side )
{
   assert( rowValues.size() == 4 );
   const Vec<double> coefficients{ 0.0, 0.0, 0.0, 0.0 };
   const Vec<double> upperBounds{ 1.0, 1.0, 1.0, 1.0 };
   const Vec<double> lowerBounds{ 0.0, 0.0, 0.0, 0.0 };
   const Vec<uint8_t> upper_bound_infinity{ 0, 0, 0, 0 };

   const Vec<uint8_t> isIntegral{ 1, 1, 1, 1 };

   const Vec<uint8_t> isLefthandsideInfinity{ 0 };
   const Vec<uint8_t> isRighthandsideInfinity{ 0 };
   const Vec<double> rhs{ side };
   const Vec<double> lhs{ side };
   const Vec<std::string> rowNames{ "A1" };
   const Vec<std::string> columnNames{ "c1", "c2", "c3", "c4" };
   const Vec<std::tuple<int, int, double>> entries{
      std::tuple<int, int, double>{ 0, 0, rowValues[0] },
      std::tuple<int, int, double>{ 0, 1, rowValues[1] },
      std::tuple<int, int, double>{ 0, 2, rowValues[2] },
      std::tuple<int, int, double>{ 0, 3, rowValues[3] },
   };

   ProblemBuilder<double> pb;
   pb.reserve( static_cast<int>(entries.size()), static_cast<int>(rowNames.size()), static_cast<int>(columnNames.size()) );
   pb.setNumRows( static_cast<int>(rowNames.size()) );
   pb.setNumCols( static_cast<int>(columnNames.size()) );
   pb.setColUbAll( upperBounds );
   pb.setColUbInfAll( upper_bound_infinity );
   pb.setColLbAll( lowerBounds );
   pb.setObjAll( coefficients );
   pb.setObjOffset( 0.0 );
   pb.setColIntegralAll( isIntegral );
   pb.setRowRhsAll( rhs );
   pb.setRowLhsInfAll( isLefthandsideInfinity );
   pb.setRowRhsInfAll( isRighthandsideInfinity );
   pb.addEntryAll( entries );
   pb.setColNameAll( columnNames );
   pb.setProblemName( "integral singleton column" );
   Problem<double> problem = pb.build();
   const Num<double> num{};
   problem.getConstraintMatrix().modifyLeftHandSide( 0, num, lhs[0] );
   return problem;
}

Problem<double>
setupProblemWithSingletonStuffingColumn()
{
   const Vec<double> coefficients{ 0.0, 0.0, 0.0, -9.0699679999999994 };
   const Vec<double> upperBounds{ 0.0, 0.0, 0.0, 0.0 };
   const Vec<double> lowerBounds{ 0.0, 0.0, 0.0, 1.98 };
   const Vec<uint8_t> upper_bound_infinity{ 1, 1, 1, 1 };

   const Vec<uint8_t> isIntegral{ 0, 0, 0, 0 };

   const Vec<uint8_t> isLefthandsideInfinity{ 0, 0 };
   const Vec<uint8_t> isRighthandsideInfinity{ 1, 1 };
   const Vec<double> rhs{ 1.0, 0.0 };
   const Vec<std::string> rowNames{ "A1", "row with SingletonRow" };
   const Vec<std::string> columnNames{ "c1", "c2", "c3", "c4" };
   const Vec<std::tuple<int, int, double>> entries{
       std::tuple<int, int, double>{ 0, 0, 1.0 },
       std::tuple<int, int, double>{ 0, 1, 1.0 },
       std::tuple<int, int, double>{ 0, 2, 2.0 },
       std::tuple<int, int, double>{ 1, 0, 1.0 },
       std::tuple<int, int, double>{ 1, 1, 1.0 },
       std::tuple<int, int, double>{ 1, 2, 1.0 },
       std::tuple<int, int, double>{ 1, 3, -1.0 },
   };

   ProblemBuilder<double> pb;
   pb.reserve( static_cast<int>(entries.size()), static_cast<int>(rowNames.size()), static_cast<int>(columnNames.size()) );
   pb.setNumRows( static_cast<int>(rowNames.size()) );
   pb.setNumCols( static_cast<int>(columnNames.size()) );
   pb.setColUbAll( upperBounds );
   pb.setColUbInfAll( upper_bound_infinity );
   pb.setColLbAll( lowerBounds );
   pb.setObjAll( coefficients );
   pb.setObjOffset( 0.0 );
   pb.setColIntegralAll( isIntegral );
   pb.setRowRhsAll( rhs );
   pb.setRowLhsInfAll( isLefthandsideInfinity );
   pb.setRowRhsInfAll( isRighthandsideInfinity );
   pb.addEntryAll( entries );
   pb.setColNameAll( columnNames );
   pb.setProblemName( "singleton column" );
   Problem<double> problem = pb.build();

   return problem;
}
