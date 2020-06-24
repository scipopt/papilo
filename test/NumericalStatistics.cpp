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
   const std::string pathInstDir = "../../test/instances/";

   const std::string instance1 = "bell5.mps";
   const std::string instance2 = "blend2.mps";

   const std::string filename1 = pathInstDir + instance1;
   Problem<double> problem1 = MpsParser<double>::loadProblem( filename1 );

     // enum declaration, only needed once
   enum class boundtype{ kLE, kEq, kGE };
   // Variable declaration
   Vec<Triplet<double>> entries;
   Vec<double> rowlhs;
   Vec<double> rowrhs;
   Vec<std::string> rownames;
   Vec<std::string> colnames;

   HashMap<std::string, int> rowname2idx;
   HashMap<std::string, int> colname2idx;
   Vec<double> lb4cols;
   Vec<double> ub4cols;
   Vec<boundtype> row_type;
   Vec<RowFlags> row_flags;
   Vec<ColFlags> col_flags;
   double objoffset = 0;

   Problem<double> problem;
   // Objective
   Vec<double> coeffobj{ 0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,43000.0,43000.0,43000.0,43000.0,43000.0,43000.0,43000.0,43000.0,43000.0,43000.0,43000.0,43000.0,43000.0,43000.0,58000.0,58000.0,58000.0,59000.0,59000.0,59000.0,59000.0,60000.0,59000.0,59000.0,59000.0,58000.0,58000.0,58000.0,10000.0,10000.0,10000.0,10000.0,10000.0,10000.0,10000.0,10000.0,10000.0,10000.0,10000.0,10000.0,10000.0,10000.0,24.5645,20.3962,14.1693,50.2605,58.0423,36.6095,39.201,48.034,29.4336,36.0182,18.7245,30.3169,5.3655,25.55,20.7977,1.825,2.45645,2.03962,1.41693,5.02605,5.80423,3.66095,3.9201,4.8034,2.94336,3.60182,1.87245,3.03169,0.53655,2.555,2.07977,0.1825,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,};
   problem.setObjective( coeffobj, 0.0 );

   const Objective<double> o1 = problem.getObjective();
   const Objective<double> o2 = problem1.getObjective();

   REQUIRE( o1.coefficients.size() == o2.coefficients.size() );
   REQUIRE( o1.coefficients.size() == o2.coefficients.size() );

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