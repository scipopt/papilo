/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*               This file is part of the program and library                */
/*    PaPILO --- Parallel Presolve for Integer and Linear Optimization       */
/*                                                                           */
/* Copyright (C) 2020  Konrad-Zuse-Zentrum                                   */
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


#include "fix/Heuristic.hpp"
#include "catch/catch.hpp"
#include "papilo/core/Problem.hpp"
#include "papilo/core/ProblemBuilder.hpp"

Problem<double>
setupProblemForConflictAnalysis_3();

TEST_CASE( "heuristics-all-false-best-objective", "[fix]" )
{
#ifndef PAPILO_TBB
    return;
#endif

   Problem<double> problem = setupProblemForConflictAnalysis_3();

   problem.recomputeAllActivities();

   Message msg {};
   msg.setVerbosityLevel(papilo::VerbosityLevel::kDetailed);
   double time = 5;
   Timer t{time};
   PostsolveStorage<double> storage{};
   RandomGenerator random( 0 );
   auto heuristic =
       new Heuristic<double>{ msg, {}, random, t, problem, storage, false };
   Vec<double> res{};
   Vec<double> sol = { 0.9, 0.9, 0.6, 0.3, 0.2 };

   double objVal = 50;
   heuristic->setup(random);
   bool infeasible =heuristic->perform_fix_and_propagate(
       sol, objVal, res, true, true, false,
       InfeasibleCopyStrategy::kBestObjective );
   REQUIRE( !infeasible );
   REQUIRE( res[0] == 0 );
   REQUIRE( res[1] == 0 );
   REQUIRE( res[2] == 1 );
   REQUIRE( res[3] == 0 );
   REQUIRE( res[4] == 0 );

}

TEST_CASE( "heuristics-all-false-worst-objective", "[fix]" )
{
#ifndef PAPILO_TBB
   return;
#endif
   Problem<double> problem = setupProblemForConflictAnalysis_3();

   problem.recomputeAllActivities();

   Message msg {};
   msg.setVerbosityLevel(papilo::VerbosityLevel::kDetailed);
   double time = 5;
   Timer t{time};
   PostsolveStorage<double> storage{};
   RandomGenerator random( 0 );
   auto heuristic =
       new Heuristic<double>{ msg, {}, random, t, problem, storage, false };
   Vec<double> res{};
   Vec<double> sol = { 0.9, 0.9, 0.6, 0.3, 0.2 };

   double objVal = 50;
   heuristic->setup(random);
   bool infeasible =heuristic->perform_fix_and_propagate(
       sol, objVal, res, true, true, false,
       InfeasibleCopyStrategy::kWorstObjective);
   REQUIRE( !infeasible );
   REQUIRE( res[0] == 1 );
   REQUIRE( res[1] == 1 );
   REQUIRE( res[2] == 0 );
   REQUIRE( res[3] == 1 );
   REQUIRE( res[4] == 1 );

}

TEST_CASE( "heuristics-all-false-highest-depth", "[fix]" )
{
#ifndef PAPILO_TBB
   return;
#endif
   Problem<double> problem = setupProblemForConflictAnalysis_3();

   problem.recomputeAllActivities();

   Message msg {};
   msg.setVerbosityLevel(papilo::VerbosityLevel::kDetailed);
   double time = 5;
   Timer t{time};
   PostsolveStorage<double> storage{};
   RandomGenerator random( 0 );
   auto heuristic =
       new Heuristic<double>{ msg, {}, random, t, problem, storage, false };
   Vec<double> res{};
   Vec<double> sol = { 0.9, 0.9, 0.6, 0.3, 0.2 };

   double objVal = 50;
   heuristic->setup(random);
   bool infeasible =heuristic->perform_fix_and_propagate(
       sol, objVal, res, true, true, false,
       InfeasibleCopyStrategy::kHighestDepthOfFirstConflict);
   REQUIRE( !infeasible );
   REQUIRE( res[0] == 0 );
   REQUIRE( res[1] == 0 );
   REQUIRE( res[2] == 1 );
   REQUIRE( res[3] == 0 );
   REQUIRE( res[4] == 0 );
}

TEST_CASE( "heuristics-all-false-lowest-depth", "[fix]" )
{
#ifndef PAPILO_TBB
   return;
#endif
   Problem<double> problem = setupProblemForConflictAnalysis_3();

   problem.recomputeAllActivities();

   Message msg {};
   msg.setVerbosityLevel(papilo::VerbosityLevel::kDetailed);
   double time = 5;
   Timer t{time};
   PostsolveStorage<double> storage{};
   RandomGenerator random( 0 );
   auto heuristic =
       new Heuristic<double>{ msg, {}, random, t, problem, storage, false };
   Vec<double> res{};
   Vec<double> sol = { 0.9, 0.9, 0.6, 0.3, 0.2 };

   double objVal = 50;
   heuristic->setup(random);
   bool infeasible =heuristic->perform_fix_and_propagate(
       sol, objVal, res, true, true, false,
       InfeasibleCopyStrategy::kLowestDepthOfFirstConflict);
   REQUIRE( !infeasible );
   REQUIRE( res[0] == 0 );
   REQUIRE( res[1] == 0 );
   REQUIRE( res[2] == 1 );
   REQUIRE( res[3] == 0 );
   REQUIRE( res[4] == 0 );
}

Problem<double>
setupProblemForConflictAnalysis_3()
{
   Vec<double> coefficients{ 1.0, 1.0, 1.0, 1.0, 1.0 };
   Vec<double> upperBounds{ 1.0, 1.0, 1.0, 1.0, 1.0 };
   Vec<double> lowerBounds{ 0.0, 0.0, 0.0, 0.0, 0.0 };
   Vec<uint8_t> isIntegral{ 1, 1, 1, 1, 1 };

   Vec<double> rhs{0.0, 2.0, 3.0, 2.0 };
   Vec<std::string> rowNames{ "A1", "A2", "A3", "A4" };
   Vec<std::string> columnNames{ "c1", "c2", "c3", "c4", "c5" };
   Vec<std::tuple<int, int, double>> entries{
       std::tuple<int, int, double>{ 0, 0, 1.0 },
       std::tuple<int, int, double>{ 0, 2, 1.0 },

       std::tuple<int, int, double>{ 1, 0, 1.0 },
       std::tuple<int, int, double>{ 1, 1, 1.0 },
       std::tuple<int, int, double>{ 1, 2, 1.0 },

       std::tuple<int, int, double>{ 2, 1, 1.0 },
       std::tuple<int, int, double>{ 2, 2, 1.0 },
       std::tuple<int, int, double>{ 2, 3, 1.0 },
       std::tuple<int, int, double>{ 2, 4, 1.0 },

       std::tuple<int, int, double>{ 3, 3, 1.0 },
       std::tuple<int, int, double>{ 3, 4, 1.0 },
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
   pb.setProblemName( "example for conflict analysis" );
   Problem<double> problem = pb.build();
   problem.getConstraintMatrix().modifyLeftHandSide( 0, {}, rhs[0] );
   problem.getConstraintMatrix().modifyLeftHandSide( 1, {}, rhs[1] );
   problem.getConstraintMatrix().modifyLeftHandSide( 2, {}, rhs[2] );
   problem.getConstraintMatrix().modifyLeftHandSide( 3, {}, rhs[3] );
   return problem;
}
