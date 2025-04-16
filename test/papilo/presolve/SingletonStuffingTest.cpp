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

#include "papilo/presolvers/SingletonStuffing.hpp"
#include "papilo/external/catch/catch_amalgamated.hpp"
#include "papilo/core/PresolveMethod.hpp"
#include "papilo/core/Problem.hpp"
#include "papilo/core/ProblemBuilder.hpp"
#include "papilo/core/RowFlags.hpp"

using namespace papilo;

Problem<double>
setupProblemWithSingletonStuffingColumn();

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
   PostsolveStorage<double> postsolve = PostsolveStorage<double>( problem, num, presolveOptions );
   ProblemUpdate<double> problemUpdate( problem, postsolve, statistics,
                                        presolveOptions, num , msg);
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

Problem<double>
setupProblemWithSingletonStuffingColumn()
{
   Vec<double> coefficients{ 0.0, 0.0, 0.0, -9.0699679999999994 };
   Vec<double> upperBounds{ 0.0, 0.0, 0.0, 0.0 };
   Vec<double> lowerBounds{ 0.0, 0.0, 0.0, 1.98 };
   Vec<uint8_t> upper_bound_infinity{ 1, 1, 1, 1 };

   Vec<uint8_t> isIntegral{ 0, 0, 0, 0 };

   Vec<uint8_t> isLefthandsideInfinity{ 0, 0 };
   Vec<uint8_t> isRighthandsideInfinity{ 1, 1 };
   Vec<double> rhs{ 1.0, 0 };
   Vec<double> lhs{ 1.0, -1.6239999999999999 };
   Vec<std::string> rowNames{ "A1", "row with SingletonRow" };
   Vec<std::string> columnNames{ "c1", "c2", "c3", "c4" };
   Vec<std::tuple<int, int, double>> entries{
       std::tuple<int, int, double>{ 0, 0, 1.0 },
       std::tuple<int, int, double>{ 0, 1, 1.0 },
       std::tuple<int, int, double>{ 0, 2, 2.0 },
       std::tuple<int, int, double>{ 1, 0, 1.0 },
       std::tuple<int, int, double>{ 1, 1, 1.0 },
       std::tuple<int, int, double>{ 1, 2, 1.0 },
       std::tuple<int, int, double>{ 1, 3, -1.0 },
   };

   ProblemBuilder<double> pb;
   pb.reserve( (int) entries.size(), (int) rowNames.size(), (int) columnNames.size() );
   pb.setNumRows( (int) rowNames.size() );
   pb.setNumCols( (int) columnNames.size() );
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
