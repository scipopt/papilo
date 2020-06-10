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

#include "papilo/core/Problem.hpp"
#include "papilo/io/MpsParser.hpp"
#include "catch/catch.hpp"
#include "papilo/misc/NumericalStatistics.hpp"
#include "papilo/misc/fmt.hpp"

#include <string>
#include <cmath>

using namespace papilo;

TEST_CASE( "accurate-numerical-statistics",
           "[misc]" )
{
   // Load problems and check values
   const std::string pathInstDir = "instances/";

   const std::string instance1 = "bell5.mps";
   const std::string instance2 = "blend2.mps";

   const std::string filename1 = pathInstDir + instance1;
   Problem<double> problem1 = MpsParser<double>::loadProblem( filename1 );

   NumericalStatistics<double> nstats1(problem1);
   const Num_stats<double>& stats = nstats1.getNum_stats();
   nstats1.printStatistics();
   //bell5.mps
   REQUIRE( fmt::format("{:.0e}", double( stats.matrixMin ) ) == "8e-05" );
   REQUIRE( fmt::format("{:.0e}", stats.matrixMax ) == "1e+03" );
   REQUIRE( fmt::format("{:.0e}", stats.objMin ) == "2e-01" );
   REQUIRE( fmt::format("{:.0e}", stats.objMax ) == "6e+04" );
   REQUIRE( fmt::format("{:.0e}", stats.boundsMin ) == "1e+00" );
   REQUIRE( fmt::format("{:.0e}", stats.boundsMax ) == "1e+04" );
   REQUIRE( fmt::format("{:.0e}", stats.rhsMin ) == "1e+00" );
   REQUIRE( fmt::format("{:.0e}", stats.rhsMax ) == "7e+03" );
   // REQUIRE( round( stats.colDynamism ) ==  );
   // REQUIRE( round( stats.rowDynamism ) ==  );

   Problem<double> problem2 = MpsParser<double>::loadProblem( pathInstDir + instance2 );

   NumericalStatistics<double> nstats2(problem2);
   const Num_stats<double>& stats2 = nstats2.getNum_stats();
   // blend2.mps
   REQUIRE( fmt::format("{:.0e}", stats2.matrixMin ) == "1e+00" );
   REQUIRE( fmt::format("{:.0e}", stats2.matrixMax ) == "7e+03" );
   REQUIRE( fmt::format("{:.0e}", stats2.objMin ) == "5e-01" );
   REQUIRE( fmt::format("{:.0e}", stats2.objMax ) == "2e+01" );
   REQUIRE( fmt::format("{:.0e}", stats2.boundsMin ) == "1e+00" );
   REQUIRE( fmt::format("{:.0e}", stats2.boundsMax ) == "2e+04" );
   REQUIRE( fmt::format("{:.0e}", stats2.rhsMin ) == "9e+00" );
   REQUIRE( fmt::format("{:.0e}", stats2.rhsMax ) == "1e+03" );
   // REQUIRE( round( stats.colDynamism ) ==  );
   // REQUIRE( round( stats.rowDynamism ) ==  );
}