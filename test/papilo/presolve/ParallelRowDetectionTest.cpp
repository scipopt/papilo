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

#include "papilo/presolvers/ParallelRowDetection.hpp"
#include "catch/catch.hpp"
#include "papilo/core/PresolveMethod.hpp"
#include "papilo/core/Problem.hpp"
#include "papilo/core/ProblemBuilder.hpp"

using namespace papilo;

Problem<double>
setupProblemWithParallelRow();

TEST_CASE( "happy-path-parallel-row-detection", "[presolve]" )
{
   Num<double> num{};
   Message msg{};
   Problem<double> problem = setupProblemWithParallelRow();
   Statistics statistics{};
   PresolveOptions presolveOptions{};
   Postsolve<double> postsolve =
       Postsolve<double>( problem, num );
   ProblemUpdate<double> problemUpdate( problem, postsolve, statistics,
                                                presolveOptions, num, msg );
   problemUpdate.checkChangedActivities();
   ParallelRowDetection<double> presolvingMethod{};
   Reductions<double> reductions{};

   PresolveStatus presolveStatus =
       presolvingMethod.execute( problem, problemUpdate, num, reductions );

   REQUIRE( presolveStatus == PresolveStatus::kReduced );
   REQUIRE( reductions.size() == 3 );
   REQUIRE( reductions.getReduction( 0 ).col == RowReduction::LOCKED );
   REQUIRE( reductions.getReduction( 0 ).row == 2 );
   REQUIRE( reductions.getReduction( 0 ).newval == 0 );

   REQUIRE( reductions.getReduction( 1 ).row == 0 );
   REQUIRE( reductions.getReduction( 1 ).col == RowReduction::LOCKED );
   REQUIRE( reductions.getReduction( 1 ).newval == 0 );

   REQUIRE( reductions.getReduction( 2 ).col ==
            RowReduction::REDUNDANT );
   REQUIRE( reductions.getReduction( 2 ).newval == 0 );
   REQUIRE( reductions.getReduction( 2 ).row == 0 );
}

Problem<double>
setupProblemWithParallelRow()
{
   Vec<double> coefficients{ 1.0, 1.0, 1.0 };
   Vec<double> upperBounds{ 10.0, 10.0, 10.0 };
   Vec<double> lowerBounds{ 0.0, 0.0, 0.0 };
   Vec<uint8_t> isIntegral{ 1, 1, 1 };

   Vec<double> rhs{ 1.0, 2.0, 3.0 };
   Vec<std::string> rowNames{ "A1", "A2", "A3" };
   Vec<std::string> columnNames{ "c1", "c2", "c3" };
   Vec<std::tuple<int, int, double>> entries{
       std::tuple<int, int, double>{ 0, 0, 1.0 },
       std::tuple<int, int, double>{ 0, 1, 1.0 },
       std::tuple<int, int, double>{ 0, 2, 2.0 },
       std::tuple<int, int, double>{ 1, 1, 2.0 },
       std::tuple<int, int, double>{ 1, 1, 2.0 },
       std::tuple<int, int, double>{ 2, 0, 3.0 },
       std::tuple<int, int, double>{ 2, 1, 3.0 },
       std::tuple<int, int, double>{ 2, 2, 6.0 },
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
   pb.setProblemName( "matrix with parallel rows (0 and 2)" );
   Problem<double> problem = pb.build();
   return problem;
}
