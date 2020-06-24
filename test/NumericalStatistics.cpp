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
   int nCols = 104; int nRows = 91;
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

   // Constraint Matrix
   Vec<std::tuple<int, int, double>> entries{{0,0,-1.0},{0,1,1.0},{1,1,-1.0},{1,2,1.0},{2,2,-1.0},{2,3,1.0},{3,3,-1.0},{3,4,1.0},{4,4,-1.0},{4,5,1.0},{5,4,-1.0},{5,6,1.0},{6,6,-1.0},{6,7,1.0},{7,7,-1.0},{7,8,1.0},{8,3,-1.0},{8,9,1.0},{9,9,-1.0},{9,10,1.0},{10,10,-1.0},{10,11,1.0},{11,1,-1.0},{11,12,1.0},{12,12,-1.0},{12,13,1.0},{13,13,-1.0},{13,14,1.0},{14,0,-1.0},{14,15,1.0},{15,0,-20.0},{15,16,1.0},{15,30,1.0},{16,1,-20.0},{16,17,1.0},{16,31,1.0},{17,2,-20.0},{17,18,1.0},{17,32,1.0},{18,3,-20.0},{18,19,1.0},{18,33,1.0},{19,4,-20.0},{19,20,1.0},{19,34,1.0},{20,5,-20.0},{20,21,1.0},{20,35,1.0},{21,6,-20.0},{21,22,1.0},{21,36,1.0},{22,8,-20.0},{22,23,1.0},{22,37,1.0},{23,9,-20.0},{23,24,1.0},{23,38,1.0},{24,10,-20.0},{24,25,1.0},{24,39,1.0},{25,11,-20.0},{25,26,1.0},{25,40,1.0},{26,12,-20.0},{26,27,1.0},{26,41,1.0},{27,13,-20.0},{27,28,1.0},{27,42,1.0},{28,15,-20.0},{28,29,1.0},{28,43,1.0},{29,16,-672.0},{29,30,-1344.0},{29,74,-1.0},{29,75,1.0},{29,89,1.0},{29,90,1.0},{30,17,-672.0},{30,31,-1344.0},{30,75,-1.0},{30,76,1.0},{30,86,1.0},{30,91,1.0},{31,18,-672.0},{31,32,-1344.0},{31,76,-1.0},{31,77,1.0},{31,92,1.0},{32,19,-672.0},{32,33,-1344.0},{32,77,-1.0},{32,78,1.0},{32,83,1.0},{32,93,1.0},{33,20,-672.0},{33,34,-1344.0},{33,78,-1.0},{33,79,1.0},{33,80,1.0},{33,94,1.0},{34,21,-672.0},{34,35,-1344.0},{34,79,-1.0},{34,95,1.0},{35,22,-672.0},{35,36,-1344.0},{35,80,-1.0},{35,81,1.0},{35,96,1.0},{36,81,-1.0},{36,82,1.0},{37,23,-672.0},{37,37,-1344.0},{37,82,-1.0},{37,97,1.0},{38,24,-672.0},{38,38,-1344.0},{38,83,-1.0},{38,84,1.0},{38,98,1.0},{39,25,-672.0},{39,39,-1344.0},{39,84,-1.0},{39,85,1.0},{39,99,1.0},{40,26,-672.0},{40,40,-1344.0},{40,85,-1.0},{40,100,1.0},{41,27,-672.0},{41,41,-1344.0},{41,86,-1.0},{41,87,1.0},{41,101,1.0},{42,28,-672.0},{42,42,-1344.0},{42,87,-1.0},{42,88,1.0},{42,102,1.0},{43,88,-1.0},{44,29,-672.0},{44,43,-1344.0},{44,89,-1.0},{44,103,1.0},{45,44,-24.0},{45,90,1.0},{46,45,-24.0},{46,91,1.0},{47,46,-24.0},{47,92,1.0},{48,47,-24.0},{48,93,1.0},{49,48,-24.0},{49,94,1.0},{50,49,-24.0},{50,95,1.0},{51,50,-24.0},{51,96,1.0},{52,51,-24.0},{52,97,1.0},{53,52,-24.0},{53,98,1.0},{54,53,-24.0},{54,99,1.0},{55,54,-24.0},{55,100,1.0},{56,55,-24.0},{56,101,1.0},{57,56,-24.0},{57,102,1.0},{58,57,-24.0},{58,103,1.0},{59,58,-1.0},{59,59,1.0},{59,73,1.0},{59,90,-1.0},{60,59,-1.0},{60,60,1.0},{60,70,1.0},{60,91,-1.0},{61,60,-1.0},{61,61,1.0},{61,92,-1.0},{62,61,-1.0},{62,62,1.0},{62,67,1.0},{62,93,-1.0},{63,62,-1.0},{63,63,1.0},{63,64,1.0},{63,94,-1.0},{64,63,-1.0},{64,95,-1.0},{65,64,-1.0},{65,65,1.0},{65,96,-1.0},{66,65,-1.0},{66,66,1.0},{67,66,-1.0},{67,97,-1.0},{68,67,-1.0},{68,68,1.0},{68,98,-1.0},{69,68,-1.0},{69,69,1.0},{69,99,-1.0},{70,69,-1.0},{70,100,-1.0},{71,70,-1.0},{71,71,1.0},{71,101,-1.0},{72,71,-1.0},{72,72,1.0},{72,102,-1.0},{73,72,-1.0},{74,73,-1.0},{74,103,-1.0},{75,0,1.0},{75,58,0.000833},{75,74,8.3e-05},{76,1,1.0},{76,59,0.000833},{76,75,8.3e-05},{77,2,1.0},{77,60,0.000833},{77,76,8.3e-05},{78,3,1.0},{78,61,0.000833},{78,77,8.3e-05},{79,4,1.0},{79,62,0.000833},{79,78,8.3e-05},{80,5,1.0},{80,63,0.000833},{80,79,8.3e-05},{81,6,1.0},{81,64,0.000833},{81,80,8.3e-05},{82,7,1.0},{82,65,0.000833},{82,81,8.3e-05},{83,8,1.0},{83,66,0.000833},{83,82,8.3e-05},{84,9,1.0},{84,67,0.000833},{84,83,8.3e-05},{85,10,1.0},{85,68,0.000833},{85,84,8.3e-05},{86,11,1.0},{86,69,0.000833},{86,85,8.3e-05},{87,12,1.0},{87,70,0.000833},{87,86,8.3e-05},{88,13,1.0},{88,71,0.000833},{88,87,8.3e-05},{89,14,1.0},{89,72,0.000833},{89,88,8.3e-05},{90,15,1.0},{90,73,0.000833},{90,89,8.3e-05},};

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