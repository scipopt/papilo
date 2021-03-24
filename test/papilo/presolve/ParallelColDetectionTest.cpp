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

#include "papilo/presolvers/ParallelColDetection.hpp"
#include "catch/catch.hpp"
#include "papilo/core/PresolveMethod.hpp"
#include "papilo/core/Problem.hpp"
#include "papilo/core/ProblemBuilder.hpp"

using namespace papilo;

Problem<double>
setupProblemWithParallelColumns();

Problem<double>
setupExample8ofChapter6Dot3InPresolveReductions();

TEST_CASE( "happy-path-parallel-column-detection", "[presolve]" )
{
   Num<double> num{};
   Message msg{};
   Problem<double> problem = setupProblemWithParallelColumns();
   Statistics statistics{};
   PresolveOptions presolveOptions{};
   Postsolve<double> postsolve =
       Postsolve<double>( problem, num );
   ProblemUpdate<double> problemUpdate( problem, postsolve, statistics,
                                                presolveOptions, num, msg );
   problemUpdate.checkChangedActivities();
   ParallelColDetection<double> presolvingMethod{};
   Reductions<double> reductions{};

   PresolveStatus presolveStatus =
       presolvingMethod.execute( problem, problemUpdate, num, reductions );

   REQUIRE( presolveStatus == PresolveStatus::kReduced );
   REQUIRE( reductions.size() == 3 );
   REQUIRE( reductions.getReduction( 0 ).row == ColReduction::LOCKED );
   REQUIRE( reductions.getReduction( 0 ).col == 2 );
   REQUIRE( reductions.getReduction( 0 ).newval == 0 );

   REQUIRE( reductions.getReduction( 1 ).col == 1 );
   REQUIRE( reductions.getReduction( 1 ).row == ColReduction::LOCKED );
   REQUIRE( reductions.getReduction( 1 ).newval == 0 );

   REQUIRE( reductions.getReduction( 2 ).row ==
            ColReduction::PARALLEL );
   REQUIRE( reductions.getReduction( 2 ).col == 2 );
   REQUIRE( reductions.getReduction( 2 ).newval == 1 );
}

TEST_CASE( "example-8-from-6.3-Presolve-Reductions-in-MIP", "[presolve]" )
{
   Num<double> num{};
   Message msg{};
   Problem<double> problem =
       setupExample8ofChapter6Dot3InPresolveReductions();
   Statistics statistics{};
   PresolveOptions presolveOptions{};
   Postsolve<double> postsolve =
       Postsolve<double>( problem, num );
   ProblemUpdate<double> problemUpdate( problem, postsolve, statistics,
                                                presolveOptions, num, msg );
   problemUpdate.checkChangedActivities();
   ParallelColDetection<double> presolvingMethod{};
   Reductions<double> reductions{};

   PresolveStatus presolveStatus =
       presolvingMethod.execute( problem, problemUpdate, num, reductions );

   REQUIRE( presolveStatus == PresolveStatus::kUnchanged );
}

Problem<double>
setupProblemWithParallelColumns()
{
   Vec<double> coefficients{ 1.0, 1.0, 2.0 };
   Vec<double> upperBounds{ 10.0, 10.0, 10.0 };
   Vec<double> lowerBounds{ 0.0, 0.0, 0.0 };
   Vec<uint8_t> isIntegral{ 1, 1, 1 };

   Vec<double> rhs{ 1.0, 2.0, 3.0 };
   Vec<std::string> rowNames{ "A1", "A2", "A3" };
   Vec<std::string> columnNames{ "c1", "c2", "c3" };
   Vec<std::tuple<int, int, double>> entries{
       std::tuple<int, int, double>{ 0, 0, 1.0 },
       std::tuple<int, int, double>{ 1, 1, 1.0 },
       std::tuple<int, int, double>{ 1, 2, 2.0 },
       std::tuple<int, int, double>{ 2, 1, -2.0 },
       std::tuple<int, int, double>{ 2, 2, -4.0 },
   };

   ProblemBuilder<double> pb;
   pb.reserve( entries.size(), rowNames.size(), columnNames.size() );
   pb.setNumRows( rowNames.size() );
   pb.setNumCols( columnNames.size() );
   pb.setColUbAll( upperBounds );
   pb.setColLbAll( lowerBounds );
   pb.setObjAll( coefficients );
   pb.setObjOffset( 0.0 );
   pb.setColIntegralAll( isIntegral );
   pb.setRowRhsAll( rhs );
   pb.addEntryAll( entries );
   pb.setColNameAll( columnNames );
   pb.setProblemName( "matrix with parallel columns (1 and 2)" );
   Problem<double> problem = pb.build();
   return problem;
}

Problem<double>
setupExample8ofChapter6Dot3InPresolveReductions()
{
   Vec<double> coefficients{ 2.0, 4.0, 1.0 };
   Vec<double> lowerBounds{ 0, 0, 0 };
   Vec<uint8_t> isIntegral{ 0, 0, 0 };

   Vec<double> rhs{ -10 };
   Vec<std::string> rowNames{ "A1" };
   Vec<std::string> columnNames{ "x", "y", "z" };
   Vec<std::tuple<int, int, double>> entries{
       std::tuple<int, int, double>{ 0, 0, -1 },
       std::tuple<int, int, double>{ 0, 1, -2 },
       std::tuple<int, int, double>{ 0, 2, -1 } };

   ProblemBuilder<double> pb;
   pb.reserve( entries.size(), rowNames.size(), columnNames.size() );
   pb.setNumRows( rowNames.size() );
   pb.setNumCols( columnNames.size() );
   pb.setColLbAll( lowerBounds );
   pb.setObjAll( coefficients );
   pb.setObjOffset( 0.0 );
   pb.setColIntegralAll( isIntegral );
   pb.setRowRhsAll( rhs );
   pb.addEntryAll( entries );
   pb.setColNameAll( columnNames );
   pb.setProblemName(
       "matrix Example 8 of chapter 6.3 in Presolve Reductions in MIP" );
   Problem<double> problem = pb.build();
   return problem;
}