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

papilo::Problem<double>
setupProblemWithParallelColumns();

papilo::Problem<double>
setupExample8ofChapter6Dot3InPresolveReductions();

TEST_CASE( "happy-path-parallel-column-detection", "[presolve]" )
{
   papilo::Num<double> num{};
   papilo::Problem<double> problem = setupProblemWithParallelColumns();
   papilo::Statistics statistics{};
   papilo::PresolveOptions presolveOptions{};
   papilo::Postsolve<double> postsolve =
       papilo::Postsolve<double>( problem, num );
   papilo::ProblemUpdate<double> problemUpdate( problem, postsolve, statistics,
                                                presolveOptions, num );
   problemUpdate.checkChangedActivities();
   papilo::ParallelColDetection<double> presolvingMethod{};
   papilo::Reductions<double> reductions{};

   papilo::PresolveStatus presolveStatus =
       presolvingMethod.execute( problem, problemUpdate, num, reductions );

   REQUIRE( presolveStatus == papilo::PresolveStatus::kReduced );
   REQUIRE( reductions.size() == 3 );
   REQUIRE( reductions.getReduction( 0 ).row == papilo::ColReduction::LOCKED );
   REQUIRE( reductions.getReduction( 0 ).col == 2 );
   REQUIRE( reductions.getReduction( 0 ).newval == 0 );

   REQUIRE( reductions.getReduction( 1 ).col == 1 );
   REQUIRE( reductions.getReduction( 1 ).row == papilo::ColReduction::LOCKED );
   REQUIRE( reductions.getReduction( 1 ).newval == 0 );

   REQUIRE( reductions.getReduction( 2 ).row ==
            papilo::ColReduction::PARALLEL );
   REQUIRE( reductions.getReduction( 2 ).col == 2 );
   REQUIRE( reductions.getReduction( 2 ).newval == 1 );
}

TEST_CASE( "example-8-from-6.3-Presolve-Reductions-in-MIP", "[presolve]" )
{
   papilo::Num<double> num{};
   papilo::Problem<double> problem =
       setupExample8ofChapter6Dot3InPresolveReductions();
   papilo::Statistics statistics{};
   papilo::PresolveOptions presolveOptions{};
   papilo::Postsolve<double> postsolve =
       papilo::Postsolve<double>( problem, num );
   papilo::ProblemUpdate<double> problemUpdate( problem, postsolve, statistics,
                                                presolveOptions, num );
   problemUpdate.checkChangedActivities();
   papilo::ParallelColDetection<double> presolvingMethod{};
   papilo::Reductions<double> reductions{};

   papilo::PresolveStatus presolveStatus =
       presolvingMethod.execute( problem, problemUpdate, num, reductions );

   REQUIRE( presolveStatus == papilo::PresolveStatus::kUnchanged );
}

papilo::Problem<double>
setupProblemWithParallelColumns()
{
   papilo::Vec<double> coefficients{ 1.0, 1.0, 2.0 };
   papilo::Vec<double> upperBounds{ 10.0, 10.0, 10.0 };
   papilo::Vec<double> lowerBounds{ 0.0, 0.0, 0.0 };
   papilo::Vec<uint8_t> isIntegral{ 1, 1, 1 };

   papilo::Vec<double> rhs{ 1.0, 2.0, 3.0 };
   papilo::Vec<std::string> rowNames{ "A1", "A2", "A3" };
   papilo::Vec<std::string> columnNames{ "c1", "c2", "c3" };
   papilo::Vec<std::tuple<int, int, double>> entries{
       std::tuple<int, int, double>{ 0, 0, 1.0 },
       std::tuple<int, int, double>{ 1, 1, 1.0 },
       std::tuple<int, int, double>{ 1, 2, 2.0 },
       std::tuple<int, int, double>{ 2, 1, -2.0 },
       std::tuple<int, int, double>{ 2, 2, -4.0 },
   };

   papilo::ProblemBuilder<double> pb;
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
   papilo::Problem<double> problem = pb.build();
   return problem;
}

papilo::Problem<double>
setupExample8ofChapter6Dot3InPresolveReductions()
{
   papilo::Vec<double> coefficients{ 2.0, 4.0, 1.0 };
   papilo::Vec<double> lowerBounds{ 0, 0, 0 };
   papilo::Vec<uint8_t> isIntegral{ 0, 0, 0 };

   papilo::Vec<double> rhs{ -10 };
   papilo::Vec<std::string> rowNames{ "A1" };
   papilo::Vec<std::string> columnNames{ "x", "y", "z" };
   papilo::Vec<std::tuple<int, int, double>> entries{
       std::tuple<int, int, double>{ 0, 0, -1 },
       std::tuple<int, int, double>{ 0, 1, -2 },
       std::tuple<int, int, double>{ 0, 2, -1 } };

   papilo::ProblemBuilder<double> pb;
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
   papilo::Problem<double> problem = pb.build();
   return problem;
}