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
#include "papilo/core/ProblemBuilder.hpp"
#include "papilo/io/MpsParser.hpp"
#include "catch/catch.hpp"
#include "papilo/misc/NumericalStatistics.hpp"
#include "papilo/misc/fmt.hpp"
#include "../test/instances/Instances.hpp"

#include <string>
#include <cmath>

#include "papilo/misc/Flags.hpp"
#include "papilo/misc/String.hpp"
#include "papilo/core/ConstraintMatrix.hpp"
#include "papilo/core/Objective.hpp"
#include <memory>
#include <tuple>

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

   Problem<double> problem2 = instances::bell5();

   /// CHECK IF CONVMPS WORKS CORRECTLY

   // Various variables
   REQUIRE( problem1.getNCols() == problem2.getNCols() );
   REQUIRE( problem1.getNRows() == problem2.getNRows() );
   REQUIRE( problem1.getNumIntegralCols() == problem2.getNumIntegralCols() );
   REQUIRE( problem1.getNumContinuousCols() == problem2.getNumContinuousCols() );
   REQUIRE( problem1.getNRows() == problem2.getNRows() );
   // Objective
   const Objective<double> o1 = problem1.getObjective();
   const Objective<double> o2 = problem2.getObjective();

   REQUIRE( o1.coefficients.size() == o2.coefficients.size() );
   REQUIRE( o1.coefficients == o2.coefficients );

   const ConstraintMatrix<double>& cm1 = problem1.getConstraintMatrix();
   const ConstraintMatrix<double>& cm2 = problem2.getConstraintMatrix();
   // Rhs Lhs
   REQUIRE( cm2.getLeftHandSides() == cm1.getLeftHandSides() );
   REQUIRE( cm2.getRightHandSides() == cm1.getRightHandSides() );

   // Row Flags
   Vec<RowFlags> row_flags1 = cm1.getRowFlags();
   Vec<RowFlags> row_flags2 = cm2.getRowFlags();

   REQUIRE( row_flags2.size() == row_flags1.size() );
   for(int r = 0; r < row_flags1.size(); ++r)
   {
      REQUIRE( row_flags2[r].test( RowFlag::NONE ) == row_flags1[r].test( RowFlag::NONE ) );
      REQUIRE( row_flags2[r].test( RowFlag::kLhsInf ) == row_flags1[r].test( RowFlag::kLhsInf ) );
      REQUIRE( row_flags2[r].test( RowFlag::kRhsInf ) == row_flags1[r].test( RowFlag::kRhsInf ) );
      REQUIRE( row_flags2[r].test( RowFlag::kIntegral ) == row_flags1[r].test( RowFlag::kIntegral ) );
      REQUIRE( row_flags2[r].test( RowFlag::kRedundant ) == row_flags1[r].test( RowFlag::kRedundant ) );
   }

   // Variable Domains
   VariableDomains<double> vd1 = problem1.getVariableDomains();
   VariableDomains<double> vd2 = problem2.getVariableDomains();

   REQUIRE( vd1.lower_bounds == vd2.lower_bounds );
   REQUIRE( vd1.upper_bounds == vd2.upper_bounds );
   for( int c = 0; c < problem1.getNCols();  ++c )
   {
      REQUIRE( vd1.flags[c].test( ColFlag::kNone ) == vd2.flags[c].test( ColFlag::kNone ) );
      REQUIRE( vd1.flags[c].test( ColFlag::kLbInf ) == vd2.flags[c].test( ColFlag::kLbInf ) );
      REQUIRE( vd1.flags[c].test( ColFlag::kLbHuge ) == vd2.flags[c].test( ColFlag::kLbHuge ) );
      REQUIRE( vd1.flags[c].test( ColFlag::kUbInf ) == vd2.flags[c].test( ColFlag::kUbInf ) );
      REQUIRE( vd1.flags[c].test( ColFlag::kUbHuge ) == vd2.flags[c].test( ColFlag::kUbHuge ) );
      REQUIRE( vd1.flags[c].test( ColFlag::kIntegral ) == vd2.flags[c].test( ColFlag::kIntegral ) );
      REQUIRE( vd1.flags[c].test( ColFlag::kFixed ) == vd2.flags[c].test( ColFlag::kFixed ) );
      REQUIRE( vd1.flags[c].test( ColFlag::kSubstituted ) == vd2.flags[c].test( ColFlag::kSubstituted ) );
      REQUIRE( vd1.flags[c].test( ColFlag::kImplInt ) == vd2.flags[c].test( ColFlag::kImplInt ) );
   }

   // Coefficients
   const SparseStorage<double>& ss1 = cm1.getConstraintMatrix();
   const SparseStorage<double>& ss2 = cm2.getConstraintMatrix();

   REQUIRE( ss1.getValuesVec() == ss2.getValuesVec() );

   // Row/ColNames
   REQUIRE( problem1.getVariableNames() == problem2.getVariableNames() );
   REQUIRE( problem1.getConstraintNames() == problem2.getConstraintNames() );

   /// END CHECK

   /// Check if statistics are printed correctly

   //bell5
   Problem<double> prob = instances::bell5();
   NumericalStatistics<double> nstats(prob);
   const Num_stats<double>& stats = nstats.getNum_stats();
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

   // blend2
   Problem<double> prob1 = instances::blend2();
   NumericalStatistics<double> nstats1(prob1);
   const Num_stats<double>& stats1 = nstats1.getNum_stats();
   REQUIRE( fmt::format("{:.0e}", stats1.matrixMin ) == "1e+00" );
   REQUIRE( fmt::format("{:.0e}", stats1.matrixMax ) == "7e+03" );
   REQUIRE( fmt::format("{:.0e}", stats1.objMin ) == "5e-01" );
   REQUIRE( fmt::format("{:.0e}", stats1.objMax ) == "2e+01" );
   REQUIRE( fmt::format("{:.0e}", stats1.boundsMin ) == "1e+00" );
   REQUIRE( fmt::format("{:.0e}", stats1.boundsMax ) == "2e+04" );
   REQUIRE( fmt::format("{:.0e}", stats1.rhsMin ) == "9e+00" );
   REQUIRE( fmt::format("{:.0e}", stats1.rhsMax ) == "1e+03" );
   // REQUIRE( round( stats.colDynamism ) ==  );
   // REQUIRE( round( stats.rowDynamism ) ==  );

   // dcmulti
   Problem<double> prob2 = instances::dcmulti();
   NumericalStatistics<double> nstats2(prob2);
   const Num_stats<double>& stats2 = nstats2.getNum_stats();
   REQUIRE( fmt::format("{:.0e}", stats2.matrixMin ) == "1e+00" );
   REQUIRE( fmt::format("{:.0e}", stats2.matrixMax ) == "6e+02" );
   REQUIRE( fmt::format("{:.0e}", stats2.objMin ) == "4e-01" );
   REQUIRE( fmt::format("{:.0e}", stats2.objMax ) == "2e+03" );
   REQUIRE( fmt::format("{:.0e}", stats2.boundsMin ) == "1e+00" );
   REQUIRE( fmt::format("{:.0e}", stats2.boundsMax ) == "1e+00" );
   REQUIRE( fmt::format("{:.0e}", stats2.rhsMin ) == "1e+00" );
   REQUIRE( fmt::format("{:.0e}", stats2.rhsMax ) == "3e+02" );
   // REQUIRE( round( stats.colDynamism ) ==  );
   // REQUIRE( round( stats.rowDynamism ) ==  );

   // egout
   Problem<double> prob3 = instances::egout();
   NumericalStatistics<double> nstats3(prob3);
   const Num_stats<double>& stats3 = nstats3.getNum_stats();
   REQUIRE( fmt::format("{:.0e}", stats3.matrixMin ) == "1e+00" );
   REQUIRE( fmt::format("{:.0e}", stats3.matrixMax ) == "1e+02" );
   REQUIRE( fmt::format("{:.0e}", stats3.objMin ) == "1e-03" );
   REQUIRE( fmt::format("{:.0e}", stats3.objMax ) == "4e+01" );
   REQUIRE( fmt::format("{:.0e}", stats3.boundsMin ) == "7e-02" );
   REQUIRE( fmt::format("{:.0e}", stats3.boundsMax ) == "2e+01" );
   REQUIRE( fmt::format("{:.0e}", stats3.rhsMin ) == "0e+00" );
   REQUIRE( fmt::format("{:.0e}", stats3.rhsMax ) == "0e+00" );
   // REQUIRE( round( stats.colDynamism ) ==  );
   // REQUIRE( round( stats.rowDynamism ) ==  );

   // enigma
   Problem<double> prob4 = instances::enigma();
   NumericalStatistics<double> nstats4(prob4);
   const Num_stats<double>& stats4 = nstats4.getNum_stats();
   REQUIRE( fmt::format("{:.0e}", stats4.matrixMin ) == "1e+00" );
   REQUIRE( fmt::format("{:.0e}", stats4.matrixMax ) == "9e+05" );
   REQUIRE( fmt::format("{:.0e}", stats4.objMin ) == "1e+00" );
   REQUIRE( fmt::format("{:.0e}", stats4.objMax ) == "9e+00" );
   REQUIRE( fmt::format("{:.0e}", stats4.boundsMin ) == "1e+00" );
   REQUIRE( fmt::format("{:.0e}", stats4.boundsMax ) == "1e+00" );
   REQUIRE( fmt::format("{:.0e}", stats4.rhsMin ) == "1e+00" );
   REQUIRE( fmt::format("{:.0e}", stats4.rhsMax ) == "1e+00" );
   // REQUIRE( round( stats.colDynamism ) ==  );
   // REQUIRE( round( stats.rowDynamism ) ==  );

   // flugpl
   Problem<double> prob5 = instances::flugpl();
   NumericalStatistics<double> nstats5(prob5);
   const Num_stats<double>& stats5 = nstats5.getNum_stats();
   REQUIRE( fmt::format("{:.0e}", stats5.matrixMin ) == "9e-01" );
   REQUIRE( fmt::format("{:.0e}", stats5.matrixMax ) == "2e+02" );
   REQUIRE( fmt::format("{:.0e}", stats5.objMin ) == "3e+01" );
   REQUIRE( fmt::format("{:.0e}", stats5.objMax ) == "3e+03" );
   REQUIRE( fmt::format("{:.0e}", stats5.boundsMin ) == "2e+01" );
   REQUIRE( fmt::format("{:.0e}", stats5.boundsMax ) == "8e+01" );
   REQUIRE( fmt::format("{:.0e}", stats5.rhsMin ) == "6e+01" );
   REQUIRE( fmt::format("{:.0e}", stats5.rhsMax ) == "1e+04" );
   // REQUIRE( round( stats.colDynamism ) ==  );
   // REQUIRE( round( stats.rowDynamism ) ==  );

   // gt2
   Problem<double> prob6 = instances::gt2();
   NumericalStatistics<double> nstats6(prob6);
   const Num_stats<double>& stats6 = nstats6.getNum_stats();
   REQUIRE( fmt::format("{:.0e}", stats6.matrixMin ) == "1e+00" );
   REQUIRE( fmt::format("{:.0e}", stats6.matrixMax ) == "3e+03" );
   REQUIRE( fmt::format("{:.0e}", stats6.objMin ) == "1e+03" );
   REQUIRE( fmt::format("{:.0e}", stats6.objMax ) == "8e+03" );
   REQUIRE( fmt::format("{:.0e}", stats6.boundsMin ) == "1e+00" );
   REQUIRE( fmt::format("{:.0e}", stats6.boundsMax ) == "2e+01" );
   REQUIRE( fmt::format("{:.0e}", stats6.rhsMin ) == "1e+00" );
   REQUIRE( fmt::format("{:.0e}", stats6.rhsMax ) == "6e+03" );
   // REQUIRE( round( stats.colDynamism ) ==  );
   // REQUIRE( round( stats.rowDynamism ) ==  );

   // lseu
   Problem<double> prob7 = instances::lseu();
   NumericalStatistics<double> nstats7(prob7);
   const Num_stats<double>& stats7 = nstats7.getNum_stats();
   REQUIRE( fmt::format("{:.0e}", stats7.matrixMin ) == "1e+00" );
   REQUIRE( fmt::format("{:.0e}", stats7.matrixMax ) == "5e+02" );
   REQUIRE( fmt::format("{:.0e}", stats7.objMin ) == "6e+00" );
   REQUIRE( fmt::format("{:.0e}", stats7.objMax ) == "5e+02" );
   REQUIRE( fmt::format("{:.0e}", stats7.boundsMin ) == "1e+00" );
   REQUIRE( fmt::format("{:.0e}", stats7.boundsMax ) == "1e+00" );
   REQUIRE( fmt::format("{:.0e}", stats7.rhsMin ) == "1e+00" );
   REQUIRE( fmt::format("{:.0e}", stats7.rhsMax ) == "3e+03" );
   // REQUIRE( round( stats.colDynamism ) ==  );
   // REQUIRE( round( stats.rowDynamism ) ==  );

   // misc03
   Problem<double> prob8 = instances::misc03();
   NumericalStatistics<double> nstats8(prob8);
   const Num_stats<double>& stats8 = nstats8.getNum_stats();
   REQUIRE( fmt::format("{:.0e}", stats8.matrixMin ) == "1e+00" );
   REQUIRE( fmt::format("{:.0e}", stats8.matrixMax ) == "1e+03" );
   REQUIRE( fmt::format("{:.0e}", stats8.objMin ) == "1e+00" );
   REQUIRE( fmt::format("{:.0e}", stats8.objMax ) == "1e+00" );
   REQUIRE( fmt::format("{:.0e}", stats8.boundsMin ) == "1e+00" );
   REQUIRE( fmt::format("{:.0e}", stats8.boundsMax ) == "1e+00" );
   REQUIRE( fmt::format("{:.0e}", stats8.rhsMin ) == "1e+00" );
   REQUIRE( fmt::format("{:.0e}", stats8.rhsMax ) == "2e+02" );
   // REQUIRE( round( stats.colDynamism ) ==  );
   // REQUIRE( round( stats.rowDynamism ) ==  );

   // p0548
   Problem<double> prob9 = instances::p0548();
   NumericalStatistics<double> nstats9(prob9);
   const Num_stats<double>& stats9 = nstats9.getNum_stats();
   REQUIRE( fmt::format("{:.0e}", stats9.matrixMin ) == "1e+00" );
   REQUIRE( fmt::format("{:.0e}", stats9.matrixMax ) == "1e+04" );
   REQUIRE( fmt::format("{:.0e}", stats9.objMin ) == "5e+00" );
   REQUIRE( fmt::format("{:.0e}", stats9.objMax ) == "1e+04" );
   REQUIRE( fmt::format("{:.0e}", stats9.boundsMin ) == "1e+00" );
   REQUIRE( fmt::format("{:.0e}", stats9.boundsMax ) == "1e+00" );
   REQUIRE( fmt::format("{:.0e}", stats9.rhsMin ) == "1e+00" );
   REQUIRE( fmt::format("{:.0e}", stats9.rhsMax ) == "1e+04" );
   // REQUIRE( round( stats.colDynamism ) ==  );
   // REQUIRE( round( stats.rowDynamism ) ==  );

   // rgn
   Problem<double> prob10 = instances::rgn();
   NumericalStatistics<double> nstats10(prob10);
   const Num_stats<double>& stats10 = nstats10.getNum_stats();
   REQUIRE( fmt::format("{:.0e}", stats10.matrixMin ) == "1e+00" );
   REQUIRE( fmt::format("{:.0e}", stats10.matrixMax ) == "5e+00" );
   REQUIRE( fmt::format("{:.0e}", stats10.objMin ) == "1e+00" );
   REQUIRE( fmt::format("{:.0e}", stats10.objMax ) == "3e+00" );
   REQUIRE( fmt::format("{:.0e}", stats10.boundsMin ) == "1e+00" );
   REQUIRE( fmt::format("{:.0e}", stats10.boundsMax ) == "1e+02" );
   REQUIRE( fmt::format("{:.0e}", stats10.rhsMin ) == "1e+00" );
   REQUIRE( fmt::format("{:.0e}", stats10.rhsMax ) == "4e+00" );
   // REQUIRE( round( stats.colDynamism ) ==  );
   // REQUIRE( round( stats.rowDynamism ) ==  );
}