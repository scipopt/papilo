/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*               This file is part of the program and library                */
/*    PaPILO --- Parallel Presolve for Integer and Linear Optimization       */
/*                                                                           */
/* Copyright (C) 2020-22  Konrad-Zuse-Zentrum */
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

#include "fix/FixAndPropagateApi.cpp"
#include "catch/catch.hpp"
#include "papilo/io/SolParser.hpp"

TEST_CASE( "fix-and-propagate-api", "[fix]" )
{
   int result = 1;
   auto problem_ptr = setup( "./resources/api_test.mps", &result, 4 );
   assert( result == 0 );
   int n_cols = 3;
   auto primal_solution = new double[n_cols];
   auto sol = new double[n_cols];
   for( int i = 0; i < n_cols; i++ )
      primal_solution[i] = ( 1.0 + i ) / 10.0;
   double val = 50;
   bool success = call_algorithm( problem_ptr, primal_solution, sol, n_cols,
                                  &val, 0, 0 , 0, 1, 1);
   REQUIRE( success );
   REQUIRE( val == 9 );
   delete_problem_instance( problem_ptr );
}

TEST_CASE( "fix-and-propagate-api-simple-heuristic", "[fix]" )
{
   int result = 1;
   auto problem_ptr = setup( "./resources/api_test.mps", &result, 4 );
   assert( result == 0 );
   int n_cols = 3;
   auto sol = new double[n_cols];

   double val = 50;
   bool success = call_simple_heuristic( problem_ptr, sol, &val );
   delete_problem_instance( problem_ptr );
   REQUIRE( success );
   REQUIRE( val == 9 );
}

// TEST_CASE( "fix-and-propagate-api-specify-test", "[fix]" )
//{
//    int result = 1;
//    const char* path_to_file = "/home/alexander/git_repositories/mipcomp22/"
//                               "presolved/neos-1354092.mps.gz";
//    auto problem_ptr = setup( path_to_file,
//                              &result, 3 );
//    assert( result == 0 );
//    auto heuristic = (Heuristic<double>*)( problem_ptr );
//
//    int n_cols = heuristic->problem.getNCols();
//
//    Vec<int> mapping{};
//    for( int i = 0; i < n_cols; i++ )
//       mapping.push_back( i );
//    Solution<double> solution;
//    auto sol = new double[n_cols];
//    bool success = SolParser<double>::read(
//        "/home/alexander/Downloads/sol.txt",
//        mapping, heuristic->problem.getVariableNames(), solution );
//    double* a = &solution.primal[0];
//
//    assert(success);
//    double val = 50;
//    success = call_algorithm( problem_ptr, a, sol, n_cols, &val,  0, 0, 0, 1, 0 );
//    delete_problem_instance( problem_ptr );
//    assert(success);
// }
