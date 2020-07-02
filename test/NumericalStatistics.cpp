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
      REQUIRE( row_flags2[r].test( RowFlag::kEquation ) == row_flags1[r].test( RowFlag::kEquation ) );
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

   //bell5.mps
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

   // blend2.mps
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
}