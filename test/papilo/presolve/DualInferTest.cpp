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

#include "papilo/presolvers/DualInfer.hpp"
#include "catch/catch.hpp"
#include "papilo/core/PresolveMethod.hpp"
#include "papilo/core/Problem.hpp"
#include "papilo/core/ProblemBuilder.hpp"

using namespace papilo;

Problem<double>
setupExample10ofChapter7Dot5InPresolveReductions();


TEST_CASE( "example-10-from-7.5-Presolve-Reductions-in-MIP", "[presolve]" )
{
   Num<double> num{};
   Problem<double> problem = setupExample10ofChapter7Dot5InPresolveReductions();
   Statistics statistics{};
   PresolveOptions presolveOptions{};
   Postsolve<double> postsolve = Postsolve<double>( problem, num );
   ProblemUpdate<double> problemUpdate( problem, postsolve, statistics,
                                        presolveOptions, num );
   papilo::DualInfer<double> presolvingMethod{};

   problem.recomputeAllActivities();

   Reductions<double> reductions{};

   PresolveStatus presolveStatus =
       presolvingMethod.execute( problem, problemUpdate, num, reductions );

   REQUIRE( presolveStatus == PresolveStatus::kReduced );
   REQUIRE( reductions.getReductions().size() == 4 );
   REQUIRE( reductions.getReduction( 0 ).newval == 0 );
   REQUIRE( reductions.getReduction( 0 ).row == 0 );
   REQUIRE( reductions.getReduction( 0 ).col == papilo::RowReduction::LOCKED );

   REQUIRE( reductions.getReduction( 1 ).newval == 1 );
   REQUIRE( reductions.getReduction( 1 ).row == 0 );
   REQUIRE( reductions.getReduction( 1 ).col == papilo::RowReduction::RHS );

   REQUIRE( reductions.getReduction( 2 ).newval == 0 );
   REQUIRE( reductions.getReduction( 2 ).row ==
            papilo::ColReduction::BOUNDS_LOCKED );
   REQUIRE( reductions.getReduction( 2 ).col == 0 );

   REQUIRE( reductions.getReduction( 3 ).newval == 0 );
   REQUIRE( reductions.getReduction( 3 ).row ==
            papilo::ColReduction::FIXED );
   REQUIRE( reductions.getReduction( 3 ).col == 0 );
}

Problem<double>
setupExample10ofChapter7Dot5InPresolveReductions()
{
   // min 2.1 x + 5.9 y + 8.9 z - w
   // 1.05 y + 2.2 z - w >= 1
   // x + 2.5 y + 2w >= 3
   // y + z - 0.1 w >= 1
   Vec<double> coefficients{ 2.1, 5.9, 8.9 , -1 };
   Vec<double> lower_bounds{ 0, 0, 0, 0};
   Vec<uint8_t> isIntegral{ 0, 0, 0, 0 };
   Vec<uint8_t> lhsInfinity{ 0, 0, 0 };
   Vec<uint8_t> rhsInfinity{ 1, 1, 1 };
   Vec<uint8_t> upperBoundInfinity{ 1, 1, 1, 1 };
   Vec<uint8_t> lowerBoundInfinity{ 0, 0, 0, 0 };

   Vec<double> lhs{ 1, -3, 0.5 };
   Vec<std::string> rowNames{ "r1", "r2", "r3" };
   Vec<std::string> columnNames{ "x", "y", "z", "w" };
   Vec<std::tuple<int, int, double>> entries{
       std::tuple<int, int, double>{ 0, 1, 1.05 },
       std::tuple<int, int, double>{ 0, 2, 2.2 },
       std::tuple<int, int, double>{ 0, 3, -1.0 },
       std::tuple<int, int, double>{ 1, 0, 1.0 },
       std::tuple<int, int, double>{ 1, 1, 2.5 },
       std::tuple<int, int, double>{ 1, 3, 2.0 },
       std::tuple<int, int, double>{ 2, 1, 1.0 },
       std::tuple<int, int, double>{ 2, 2, 1.0 },
       std::tuple<int, int, double>{ 2, 2, -0.1 }
   };

   ProblemBuilder<double> pb;
   pb.reserve( entries.size(), rowNames.size(), columnNames.size() );
   pb.setNumRows( rowNames.size() );
   pb.setNumCols( columnNames.size() );
   pb.setColLbInfAll( lowerBoundInfinity );
   pb.setColLbAll( lower_bounds );
   pb.setColUbInfAll( upperBoundInfinity );
   pb.setObjAll( coefficients );
   pb.setObjOffset( 0.0 );
   pb.setColIntegralAll( isIntegral );
   pb.setRowLhsAll( lhs );
   pb.setRowLhsInfAll( lhsInfinity );
   pb.setRowRhsInfAll( rhsInfinity );
   pb.addEntryAll( entries );
   pb.setColNameAll( columnNames );
   pb.setProblemName(
       "matrix Example 10 of chapter 7.5 in Presolve Reductions in MIP" );
   Problem<double> problem = pb.build();
   return problem;
}

