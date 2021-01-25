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

#include "catch/catch.hpp"
#include "papilo/core/PresolveMethod.hpp"
#include "papilo/core/Problem.hpp"
#include "papilo/core/ProblemBuilder.hpp"
#include "papilo/core/Reductions.hpp"
#include "papilo/presolvers/FreeVarSubstitution.hpp"

#include <tuple>

papilo::Problem<double>
setupProblemForFreeVariableSubstitution();

TEST_CASE( "happy-path-test-free-variable-detection", "[presolve]" )
{
   papilo::Problem<double> problem = setupProblemForFreeVariableSubstitution();

   papilo::Num<double> num{};
   papilo::Statistics statistics{};
   papilo::PresolveOptions presolveOptions{};
   papilo::Postsolve<double> postsolve =
       papilo::Postsolve<double>( problem, num );
   problem.recomputeAllActivities();
   auto& activities = problem.getRowActivities();

   papilo::ProblemUpdate<double> problemUpdate( problem, postsolve, statistics,
                                                presolveOptions, num );

   papilo::Substitution<double> presolvingMethod{};
   papilo::Reductions<double> reductions{};

   presolvingMethod.initialize( problem, presolveOptions );

   papilo::PresolveStatus presolveStatus =
       presolvingMethod.execute( problem, problemUpdate, num, reductions );

   REQUIRE( presolveStatus == papilo::PresolveStatus::kReduced );

   REQUIRE( reductions.size() == 3 );
   REQUIRE( reductions.getReduction( 0 ).col == papilo::RowReduction::LOCKED );
   REQUIRE( reductions.getReduction( 0 ).row == 2 );
   REQUIRE( reductions.getReduction( 0 ).newval == 0 );

   REQUIRE( reductions.getReduction( 1 ).col == 3 );
   REQUIRE( reductions.getReduction( 1 ).row ==
            papilo::ColReduction::BOUNDS_LOCKED );
   REQUIRE( reductions.getReduction( 1 ).newval == 0 );

   REQUIRE( reductions.getReduction( 2 ).col == 3 );
   REQUIRE( reductions.getReduction( 2 ).row ==
            papilo::ColReduction::SUBSTITUTE );
   REQUIRE( reductions.getReduction( 2 ).newval == 2 );
}

papilo::Problem<double>
setupProblemForFreeVariableSubstitution()
{
   // min x + y + z + v + w
   // 2x + y >= 1
   // x + 2z <= 2
   // x + v + w = 1
   // |x| <= 3; y <= 1 ; z >= 0
   double inf = 10000000;
   papilo::Vec<double> coefficients{ 1.0, 1.0, 1.0, 1.0, 1.0 };
   papilo::Vec<double> upperBounds{ 3.0, 1.0, inf, inf, inf };
   papilo::Vec<double> lowerBounds{ -3.0, -inf, 0.0, -inf, -inf };
   papilo::Vec<uint8_t> upperBoundsInfinity{ 0, 0, 1, 1, 1 };
   papilo::Vec<uint8_t> lowerBoundsInfinity{ 0, 1, 0, 1, 1 };
   papilo::Vec<uint8_t> isIntegral{ 0, 0, 0, 0, 0 };
   papilo::Vec<double> rhs{ inf, 2.0, 1.0 };
   papilo::Vec<double> lhs{ 1.0, -inf, rhs[2] };
   papilo::Vec<std::string> rowNames{ "A1", "A2", "A3" };
   papilo::Vec<std::string> columnNames{ "x", "y", "z", "v", "w" };
   papilo::Vec<std::tuple<int, int, double>> entries{
       std::tuple<int, int, double>{ 0, 0, 2.0 },
       std::tuple<int, int, double>{ 0, 1, 1.0 },
       std::tuple<int, int, double>{ 1, 0, 1.0 },
       std::tuple<int, int, double>{ 1, 2, 1.0 },
       std::tuple<int, int, double>{ 2, 0, 1.0 },
       std::tuple<int, int, double>{ 2, 3, 1.0 },
       std::tuple<int, int, double>{ 2, 4, -1.0 },
   };

   papilo::ProblemBuilder<double> pb;
   pb.reserve( entries.size(), rowNames.size(), columnNames.size() );
   pb.setNumRows( rowNames.size() );
   pb.setNumCols( columnNames.size() );
   pb.setColUbAll( upperBounds );
   pb.setColLbAll( lowerBounds );
   pb.setObjAll( coefficients );
   pb.setColLbInfAll( lowerBoundsInfinity );
   pb.setColUbInfAll( upperBoundsInfinity );
   pb.setObjOffset( 0.0 );
   pb.setColIntegralAll( isIntegral );
   pb.setRowRhsAll( rhs );
   pb.addEntryAll( entries );
   pb.setColNameAll( columnNames );
   pb.setProblemName( "matrix with free variables (3,4)" );
   papilo::Problem<double> problem = pb.build();
   problem.getConstraintMatrix().modifyLeftHandSide( 2, lhs[2] );
   return problem;
}