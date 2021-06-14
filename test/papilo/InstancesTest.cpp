/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*               This file is part of the program and library                */
/*    PaPILO --- Parallel Presolve for Integer and Linear Optimization       */
/*                                                                           */
/* Copyright (C) 2020-2021  Konrad-Zuse-Zentrum                              */
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

#include "papilo/misc/NumericalStatistics.hpp"
#include "../instances/Instances.hpp"
#include "catch/catch.hpp"
#include "papilo/misc/fmt.hpp"

using namespace papilo;

TEST_CASE( "accurate-numerical-statistics", "[misc]" )
{
   /// Check if statistics are printed correctly

   // bell5
   NumericalStatistics<double> nstats( instances::bell5() );
   const Num_stats<double>& stats = nstats.getNum_stats();
   REQUIRE( fmt::format( "{:.0e}", double( stats.matrixMin ) ) == "8e-05" );
   REQUIRE( fmt::format( "{:.0e}", stats.matrixMax ) == "1e+03" );
   REQUIRE( fmt::format( "{:.0e}", stats.objMin ) == "2e-01" );
   REQUIRE( fmt::format( "{:.0e}", stats.objMax ) == "6e+04" );
   REQUIRE( fmt::format( "{:.0e}", stats.boundsMin ) == "1e+00" );
   REQUIRE( fmt::format( "{:.0e}", stats.boundsMax ) == "1e+04" );
   REQUIRE( fmt::format( "{:.0e}", stats.rhsMin ) == "1e+00" );
   REQUIRE( fmt::format( "{:.0e}", stats.rhsMax ) == "7e+03" );

   // blend2
   NumericalStatistics<double> nstats1( instances::blend2() );
   const Num_stats<double>& stats1 = nstats1.getNum_stats();
   REQUIRE( fmt::format( "{:.0e}", stats1.matrixMin ) == "1e+00" );
   REQUIRE( fmt::format( "{:.0e}", stats1.matrixMax ) == "7e+03" );
   REQUIRE( fmt::format( "{:.0e}", stats1.objMin ) == "5e-01" );
   REQUIRE( fmt::format( "{:.0e}", stats1.objMax ) == "2e+01" );
   REQUIRE( fmt::format( "{:.0e}", stats1.boundsMin ) == "1e+00" );
   REQUIRE( fmt::format( "{:.0e}", stats1.boundsMax ) == "2e+04" );
   REQUIRE( fmt::format( "{:.0e}", stats1.rhsMin ) == "9e+00" );
   REQUIRE( fmt::format( "{:.0e}", stats1.rhsMax ) == "1e+03" );

   // dcmulti
   NumericalStatistics<double> nstats2( instances::dcmulti() );
   const Num_stats<double>& stats2 = nstats2.getNum_stats();
   REQUIRE( fmt::format( "{:.0e}", stats2.matrixMin ) == "1e+00" );
   REQUIRE( fmt::format( "{:.0e}", stats2.matrixMax ) == "6e+02" );
   REQUIRE( fmt::format( "{:.0e}", stats2.objMin ) == "4e-01" );
   REQUIRE( fmt::format( "{:.0e}", stats2.objMax ) == "2e+03" );
   REQUIRE( fmt::format( "{:.0e}", stats2.boundsMin ) == "1e+00" );
   REQUIRE( fmt::format( "{:.0e}", stats2.boundsMax ) == "1e+00" );
   REQUIRE( fmt::format( "{:.0e}", stats2.rhsMin ) == "1e+00" );
   REQUIRE( fmt::format( "{:.0e}", stats2.rhsMax ) == "3e+02" );

   // egout
   NumericalStatistics<double> nstats3( instances::egout() );
   const Num_stats<double>& stats3 = nstats3.getNum_stats();
   REQUIRE( fmt::format( "{:.0e}", stats3.matrixMin ) == "1e+00" );
   REQUIRE( fmt::format( "{:.0e}", stats3.matrixMax ) == "1e+02" );
   REQUIRE( fmt::format( "{:.0e}", stats3.objMin ) == "1e-03" );
   REQUIRE( fmt::format( "{:.0e}", stats3.objMax ) == "4e+01" );
   REQUIRE( fmt::format( "{:.0e}", stats3.boundsMin ) == "7e-02" );
   REQUIRE( fmt::format( "{:.0e}", stats3.boundsMax ) == "2e+01" );
   REQUIRE( fmt::format( "{:.0e}", stats3.rhsMin ) == "0e+00" );
   REQUIRE( fmt::format( "{:.0e}", stats3.rhsMax ) == "0e+00" );

   // enigma
   Problem<double> prob4 = instances::enigma();
   NumericalStatistics<double> nstats4( prob4 );
   const Num_stats<double>& stats4 = nstats4.getNum_stats();
   REQUIRE( fmt::format( "{:.0e}", stats4.matrixMin ) == "1e+00" );
   REQUIRE( fmt::format( "{:.0e}", stats4.matrixMax ) == "9e+05" );
   REQUIRE( fmt::format( "{:.0e}", stats4.objMin ) == "1e+00" );
   REQUIRE( fmt::format( "{:.0e}", stats4.objMax ) == "9e+00" );
   REQUIRE( fmt::format( "{:.0e}", stats4.boundsMin ) == "1e+00" );
   REQUIRE( fmt::format( "{:.0e}", stats4.boundsMax ) == "1e+00" );
   REQUIRE( fmt::format( "{:.0e}", stats4.rhsMin ) == "1e+00" );
   REQUIRE( fmt::format( "{:.0e}", stats4.rhsMax ) == "1e+00" );

   // flugpl
   Problem<double> prob5 = instances::flugpl();
   NumericalStatistics<double> nstats5( prob5 );
   const Num_stats<double>& stats5 = nstats5.getNum_stats();
   REQUIRE( fmt::format( "{:.0e}", stats5.matrixMin ) == "9e-01" );
   REQUIRE( fmt::format( "{:.0e}", stats5.matrixMax ) == "2e+02" );
   REQUIRE( fmt::format( "{:.0e}", stats5.objMin ) == "3e+01" );
   REQUIRE( fmt::format( "{:.0e}", stats5.objMax ) == "3e+03" );
   REQUIRE( fmt::format( "{:.0e}", stats5.boundsMin ) == "2e+01" );
   REQUIRE( fmt::format( "{:.0e}", stats5.boundsMax ) == "8e+01" );
   REQUIRE( fmt::format( "{:.0e}", stats5.rhsMin ) == "6e+01" );
   REQUIRE( fmt::format( "{:.0e}", stats5.rhsMax ) == "1e+04" );
   REQUIRE( fmt::format( "{:.0e}", stats5.colDynamism ) == "2e+02" );
   REQUIRE( fmt::format( "{:.0e}", stats5.rowDynamism ) == "2e+02" );

   // gt2
   Problem<double> prob6 = instances::gt2();
   NumericalStatistics<double> nstats6( prob6 );
   const Num_stats<double>& stats6 = nstats6.getNum_stats();
   REQUIRE( fmt::format( "{:.0e}", stats6.matrixMin ) == "1e+00" );
   REQUIRE( fmt::format( "{:.0e}", stats6.matrixMax ) == "3e+03" );
   REQUIRE( fmt::format( "{:.0e}", stats6.objMin ) == "1e+03" );
   REQUIRE( fmt::format( "{:.0e}", stats6.objMax ) == "8e+03" );
   REQUIRE( fmt::format( "{:.0e}", stats6.boundsMin ) == "1e+00" );
   REQUIRE( fmt::format( "{:.0e}", stats6.boundsMax ) == "2e+01" );
   REQUIRE( fmt::format( "{:.0e}", stats6.rhsMin ) == "1e+00" );
   REQUIRE( fmt::format( "{:.0e}", stats6.rhsMax ) == "6e+03" );

   // lseu
   Problem<double> prob7 = instances::lseu();
   NumericalStatistics<double> nstats7( prob7 );
   const Num_stats<double>& stats7 = nstats7.getNum_stats();
   REQUIRE( fmt::format( "{:.0e}", stats7.matrixMin ) == "1e+00" );
   REQUIRE( fmt::format( "{:.0e}", stats7.matrixMax ) == "5e+02" );
   REQUIRE( fmt::format( "{:.0e}", stats7.objMin ) == "6e+00" );
   REQUIRE( fmt::format( "{:.0e}", stats7.objMax ) == "5e+02" );
   REQUIRE( fmt::format( "{:.0e}", stats7.boundsMin ) == "1e+00" );
   REQUIRE( fmt::format( "{:.0e}", stats7.boundsMax ) == "1e+00" );
   REQUIRE( fmt::format( "{:.0e}", stats7.rhsMin ) == "1e+00" );
   REQUIRE( fmt::format( "{:.0e}", stats7.rhsMax ) == "3e+03" );

   // misc03
   Problem<double> prob8 = instances::misc03();
   NumericalStatistics<double> nstats8( prob8 );
   const Num_stats<double>& stats8 = nstats8.getNum_stats();
   REQUIRE( fmt::format( "{:.0e}", stats8.matrixMin ) == "1e+00" );
   REQUIRE( fmt::format( "{:.0e}", stats8.matrixMax ) == "1e+03" );
   REQUIRE( fmt::format( "{:.0e}", stats8.objMin ) == "1e+00" );
   REQUIRE( fmt::format( "{:.0e}", stats8.objMax ) == "1e+00" );
   REQUIRE( fmt::format( "{:.0e}", stats8.boundsMin ) == "1e+00" );
   REQUIRE( fmt::format( "{:.0e}", stats8.boundsMax ) == "1e+00" );
   REQUIRE( fmt::format( "{:.0e}", stats8.rhsMin ) == "1e+00" );
   REQUIRE( fmt::format( "{:.0e}", stats8.rhsMax ) == "2e+02" );

   // p0548
   NumericalStatistics<double> nstats9( instances::p0548() );
   const Num_stats<double>& stats9 = nstats9.getNum_stats();
   REQUIRE( fmt::format( "{:.0e}", stats9.matrixMin ) == "1e+00" );
   REQUIRE( fmt::format( "{:.0e}", stats9.matrixMax ) == "1e+04" );
   REQUIRE( fmt::format( "{:.0e}", stats9.objMin ) == "5e+00" );
   REQUIRE( fmt::format( "{:.0e}", stats9.objMax ) == "1e+04" );
   REQUIRE( fmt::format( "{:.0e}", stats9.boundsMin ) == "1e+00" );
   REQUIRE( fmt::format( "{:.0e}", stats9.boundsMax ) == "1e+00" );
   REQUIRE( fmt::format( "{:.0e}", stats9.rhsMin ) == "1e+00" );
   REQUIRE( fmt::format( "{:.0e}", stats9.rhsMax ) == "1e+04" );

   // rgn
   NumericalStatistics<double> nstats10( instances::rgn() );
   const Num_stats<double>& stats10 = nstats10.getNum_stats();
   REQUIRE( fmt::format( "{:.0e}", stats10.matrixMin ) == "1e+00" );
   REQUIRE( fmt::format( "{:.0e}", stats10.matrixMax ) == "5e+00" );
   REQUIRE( fmt::format( "{:.0e}", stats10.objMin ) == "1e+00" );
   REQUIRE( fmt::format( "{:.0e}", stats10.objMax ) == "3e+00" );
   REQUIRE( fmt::format( "{:.0e}", stats10.boundsMin ) == "1e+00" );
   REQUIRE( fmt::format( "{:.0e}", stats10.boundsMax ) == "1e+02" );
   REQUIRE( fmt::format( "{:.0e}", stats10.rhsMin ) == "1e+00" );
   REQUIRE( fmt::format( "{:.0e}", stats10.rhsMax ) == "4e+00" );
}
