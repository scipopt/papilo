/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*               This file is part of the program and library                */
/*    PaPILO --- Parallel Presolve for Integer and Linear Optimization       */
/*                                                                           */
/* Copyright (C) 2020-2025 Zuse Institute Berlin (ZIB)                       */
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

#include "papilo/core/PresolveMethod.hpp"
#include "papilo/core/Problem.hpp"
#include "papilo/core/ProblemBuilder.hpp"
#include "papilo/presolvers/FixContinuous.hpp"
#include "papilo/external/catch/catch_amalgamated.hpp"

using namespace papilo;

Problem<double>
setupProblemForTest( double upperBoundForVar3 );

TEST_CASE( "happy-path-presolve-fix-continuous", "[presolve]" )
{
   double time = 0.0;
   int cause = -1;
   Timer t{time};
   Num<double> num{};
   Message msg{};
   Problem<double> problem = setupProblemForTest( num.getEpsilon() / 4 );
   Statistics statistics{};
   PresolveOptions presolveOptions{};
   PostsolveStorage<double> postsolve =
       PostsolveStorage<double>( problem, num, presolveOptions );
   ProblemUpdate<double> problemUpdate( problem, postsolve, statistics,
                                        presolveOptions, num, msg );
   FixContinuous<double> presolvingMethod{};
   Reductions<double> reductions{};

   PresolveStatus presolveStatus =
       presolvingMethod.execute( problem, problemUpdate, num, reductions, t, cause );

   REQUIRE( presolveStatus == PresolveStatus::kReduced );
   REQUIRE( reductions.size() == 2 );
   REQUIRE( reductions.getReduction( 0 ).col == 2 );
   REQUIRE( reductions.getReduction( 0 ).row == ColReduction::BOUNDS_LOCKED );
   REQUIRE( reductions.getReduction( 0 ).newval == 0 );
   REQUIRE( reductions.getReduction( 1 ).row == ColReduction::FIXED );
   REQUIRE( reductions.getReduction( 1 ).col == 2 );
   REQUIRE( reductions.getReduction( 1 ).newval == 0 );
}

TEST_CASE( "happy-path-no-presolve-fix-continuous", "[presolve]" )
{
   double time = 0.0;
   int cause = -1;
   Timer t{time};
   Num<double> num{};
   Message msg{};
   Problem<double> problem = setupProblemForTest( num.getEpsilon() );
   Statistics statistics{};
   PresolveOptions presolveOptions{};
   PostsolveStorage<double> postsolve =
       PostsolveStorage<double>( problem, num, presolveOptions );
   ProblemUpdate<double> problemUpdate( problem, postsolve, statistics,
                                        presolveOptions, num, msg );
   FixContinuous<double> presolvingMethod{};
   Reductions<double> reductions{};
   PresolveStatus presolveStatus =
       presolvingMethod.execute( problem, problemUpdate, num, reductions, t, cause );
   REQUIRE( presolveStatus == PresolveStatus::kUnchanged );
   REQUIRE( reductions.size() == 0 );
}

Problem<double>
setupProblemForTest( double upperBoundForVar3 )
{
   Vec<double> coefficients{ 1.0, 1.0, 1.0 };
   Vec<double> upperBounds{ 10.0, 10.0, upperBoundForVar3 };
   Vec<double> lowerBounds{ 0.0, 0.0, 0.0 };
   Vec<uint8_t> isIntegral{ 1, 1, 0 };

   Vec<double> rhs{ 1.0, 2.0 };
   Vec<std::string> rowNames{ "A1", "A2" };
   Vec<std::string> columnNames{ "c1", "c2", "c3" };
   Vec<std::tuple<int, int, double>> entries{
       std::tuple<int, int, double>{ 0, 0, 1.0 },
       std::tuple<int, int, double>{ 0, 1, 2.0 },
       std::tuple<int, int, double>{ 0, 2, 3.0 },
       std::tuple<int, int, double>{ 1, 1, 2.0 },
   };

   ProblemBuilder<double> pb;
   pb.reserve( 4, 2, 3 );
   pb.setNumRows( 2 );
   pb.setNumCols( 3 );
   pb.setColUbAll( upperBounds );
   pb.setColLbAll( lowerBounds );
   pb.setObjAll( coefficients );
   pb.setObjOffset( 0.0 );
   pb.setColIntegralAll( isIntegral );
   pb.setRowRhsAll( rhs );
   pb.addEntryAll( entries );
   pb.setColNameAll( columnNames );
   pb.setProblemName( "matrix with potential continuous fix" );
   Problem<double> problem = pb.build();
   return problem;
}
