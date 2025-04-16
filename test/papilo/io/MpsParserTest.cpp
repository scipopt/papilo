/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*               This file is part of the program and library                */
/*    PaPILO --- Parallel Presolve for Integer and Linear Optimization       */
/*                                                                           */
/* Copyright (C) 2020-2025 Zuse Institute Berlin (ZIB)                       */
/*                                                                           */
/* Licensed under the Apache License, Version 2.0 (the "License");           */
/* you may not use this file except in compliance with the License.          */
/* You may obtain a copy of the License at                                   */
/*                                                                           */
/*     http://www.apache.org/licenses/LICENSE-2.0                            */
/*                                                                           */
/* Unless required by applicable law or agreed to in writing, software       */
/* distributed under the License is distributed on an "AS IS" BASIS,         */
/* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  */
/* See the License for the specific language governing permissions and       */
/* limitations under the License.                                            */
/*                                                                           */
/* You should have received a copy of the Apache-2.0 license                 */
/* along with PaPILO; see the file LICENSE. If not visit scipopt.org.        */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <memory>
#include "papilo/io/MpsParser.hpp"
#include "papilo/external/catch/catch_amalgamated.hpp"
#include "papilo/core/PresolveMethod.hpp"
#include "papilo/core/Problem.hpp"

using namespace papilo;

Problem<double>
setupProblemForCoefficientStrengthening();

TEST_CASE( "mps-parser-loading-simple-problem", "[io]" )
{
   boost::optional<Problem<double>> optional = MpsParser<double>::loadProblem(
       "./resources/dual_fix_neg_inf.mps" );
   REQUIRE( optional.is_initialized() == true );
   Problem<double> problem = optional.get();
   Vec<int> expected_row_sizes{2,2,3};
   Vec<int> expected_col_sizes{3,2,2};
   REQUIRE(problem.getConstraintMatrix().getRowSizes() == expected_row_sizes);
   REQUIRE(problem.getConstraintMatrix().getColSizes() == expected_col_sizes);
}
