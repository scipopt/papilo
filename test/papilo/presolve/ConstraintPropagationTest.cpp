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

#include "papilo/presolvers/ConstraintPropagation.hpp"
#include "papilo/external/catch/catch_amalgamated.hpp"
#include "papilo/core/PresolveMethod.hpp"
#include "papilo/core/Problem.hpp"
#include "papilo/core/ProblemBuilder.hpp"

using namespace papilo;

Problem<double>
setupProblemWithConstraintPropagation( bool integer_values );

TEST_CASE( "constraint-propagation-happy-path", "[presolve]" )
{
   double time = 0.0;
   int cause = -1;
   Timer t{time};
   Num<double> num{};
   Message msg{};
   Problem<double> problem = setupProblemWithConstraintPropagation( true );
   Statistics statistics{};
   PresolveOptions presolveOptions{};
   presolveOptions.dualreds = 0;
   PostsolveStorage<double> postsolve =
       PostsolveStorage<double>( problem, num, presolveOptions );
   ProblemUpdate<double> problemUpdate( problem, postsolve, statistics,
                                        presolveOptions, num, msg );
   ConstraintPropagation<double> presolvingMethod{};
   Reductions<double> reductions{};

   problem.recomputeAllActivities();
   problemUpdate.trivialPresolve();
   problemUpdate.clearChangeInfo();

   PresolveStatus presolveStatus =
       presolvingMethod.execute( problem, problemUpdate, num, reductions, t, cause );

   REQUIRE( presolveStatus == PresolveStatus::kReduced );
   REQUIRE( reductions.size() == 8 );

   REQUIRE( reductions.getReduction( 0 ).col == RowReduction::SAVE_ROW );
   REQUIRE( reductions.getReduction( 0 ).row == 0 );

   REQUIRE( reductions.getReduction( 1 ).col == 0 );
   REQUIRE( reductions.getReduction( 1 ).row == ColReduction::UPPER_BOUND );
   REQUIRE( reductions.getReduction( 1 ).newval == 1 );

   REQUIRE( reductions.getReduction( 2 ).col == RowReduction::SAVE_ROW );
   REQUIRE( reductions.getReduction( 2 ).row == 0 );

   REQUIRE( reductions.getReduction( 3 ).row == ColReduction::UPPER_BOUND );
   REQUIRE( reductions.getReduction( 3 ).col == 1 );
   REQUIRE( reductions.getReduction( 3 ).newval == 1 );

   REQUIRE( reductions.getReduction( 4 ).col == RowReduction::SAVE_ROW );
   REQUIRE( reductions.getReduction( 4 ).row == 1 );

   REQUIRE( reductions.getReduction( 5 ).col == 1 );
   REQUIRE( reductions.getReduction( 5 ).row == ColReduction::FIXED );
   REQUIRE( reductions.getReduction( 5 ).newval == 0 );

   REQUIRE( reductions.getReduction( 6 ).col == RowReduction::SAVE_ROW );
   REQUIRE( reductions.getReduction( 6 ).row == 1 );

   REQUIRE( reductions.getReduction( 7 ).col == 2 );
   REQUIRE( reductions.getReduction( 7 ).newval == 0.1 );
   REQUIRE( reductions.getReduction( 7 ).row == ColReduction::UPPER_BOUND );
}

TEST_CASE( "constraint-propagation-no-tightening-for-lp", "[presolve]" )
{
   double time = 0.0;
   int cause = -1;
   Timer t{ time };
   Num<double> num{};
   Message msg{};
   Problem<double> problem = setupProblemWithConstraintPropagation( false );
   Statistics statistics{};
   PresolveOptions presolveOptions{};
   presolveOptions.dualreds = 0;
   PostsolveStorage<double> postsolve =
       PostsolveStorage<double>( problem, num , presolveOptions);
   ProblemUpdate<double> problemUpdate( problem, postsolve, statistics,
                                        presolveOptions, num, msg );
   ConstraintPropagation<double> presolvingMethod{};
   Reductions<double> reductions{};

   problem.recomputeAllActivities();
   problemUpdate.trivialPresolve();
   problemUpdate.clearChangeInfo();

   PresolveStatus presolveStatus =
       presolvingMethod.execute( problem, problemUpdate, num, reductions, t, cause);

   REQUIRE( presolveStatus == PresolveStatus::kUnchanged );
}

Problem<double>
setupProblemWithConstraintPropagation( bool integer_values )
{
   Vec<double> coefficients{ 1.0, 1.0, 1.0 };
   Vec<double> upperBounds{ 10.0, 10.0, 10.0 };
   Vec<double> lowerBounds{ 0.0, 0.0, 0.0 };
   Vec<uint8_t> isIntegral{ 0, integer_values, 0 };

   Vec<double> rhs{ 1.0, 2.0 };
   Vec<std::string> rowNames{ "A1", "A2" };
   Vec<std::string> columnNames{ "c1", "c2", "c3" };
   Vec<std::tuple<int, int, double>> entries{
       std::tuple<int, int, double>{ 0, 0, 1.0 },
       std::tuple<int, int, double>{ 0, 1, 1.0 },
       std::tuple<int, int, double>{ 1, 1, 20.0 },
       std::tuple<int, int, double>{ 1, 2, 20.0 },
   };

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
   pb.addEntryAll( entries );
   pb.setColNameAll( columnNames );
   pb.setProblemName( "matrix for testing constraint propagation" );
   Problem<double> problem = pb.build();
   return problem;
}
