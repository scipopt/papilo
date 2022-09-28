/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*               This file is part of the program and library                */
/*    PaPILO --- Parallel Presolve for Integer and Linear Optimization       */
/*                                                                           */
/* Copyright (C) 2020-2022 Konrad-Zuse-Zentrum                               */
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

#include "papilo/io/PboParser.hpp"
#include "papilo/core/PresolveMethod.hpp"
#include "papilo/core/Problem.hpp"
#include "papilo/external/catch/catch.hpp"

using namespace papilo;

Problem<double>
setupProblemForCoefficientStrengthening();

TEST_CASE( "pbo-parser-loading-neg-objective-problem", "[io]" )
{
   boost::optional<Problem<double>> optional =
       PboParser<double>::loadProblem( "./resources/neg_objective.opb" );
   REQUIRE( optional.is_initialized() == true );
   Problem<double> problem = optional.get();
   const double three = 3.0;
   REQUIRE( problem.getObjective().offset == three );
   REQUIRE( problem.getObjective().coefficients.size() == 1 );
   REQUIRE( problem.getObjective().coefficients[0] == -three );
}

TEST_CASE( "pbo-parser-loading-pos-objective-problem", "[io]" )
{
   boost::optional<Problem<double>> optional =
       PboParser<double>::loadProblem( "./resources/pos_objective.obp" );
   REQUIRE( optional.is_initialized() == true );
   Problem<double> problem = optional.get();
   const double three = 3.0;
   const double zero = 0.0;
   REQUIRE( problem.getObjective().offset == zero );
   REQUIRE( problem.getObjective().coefficients.size() == 1 );
   REQUIRE( problem.getObjective().coefficients[0] == three );
}

TEST_CASE( "pbo-parser-loading-knapsack-problem", "[io]" )
{
   boost::optional<Problem<double>> optional =
       PboParser<double>::loadProblem( "./resources/knapsack.obp" );
   REQUIRE( optional.is_initialized() == true );
   Problem<double> problem = optional.get();
   REQUIRE( problem.getNCols() == 5 );
   REQUIRE( problem.getNRows() == 1 );
}
