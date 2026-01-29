/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*               This file is part of the program and library                */
/*    PaPILO --- Parallel Presolve for Integer and Linear Optimization       */
/*                                                                           */
/* Copyright (C) 2020-2026 Zuse Institute Berlin (ZIB)                       */
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

#include "papilo/core/postsolve/Postsolve.hpp"
#include "papilo/external/catch/catch_amalgamated.hpp"
#include "papilo/core/postsolve/PostsolveStatus.hpp"
#include <boost/archive/binary_iarchive.hpp>

using namespace papilo;

TEST_CASE( "finding-the-right-value-in-postsolve-for-a-column-fixed-neg-inf",
           "[core]" )
{

   const Num<double> num{};
   Message msg{};
   PostsolveStorage<double> postsolveStorage{};

   std::ifstream inArchiveFile( "./resources/dual_fix_neg_inf.postsolve",
                                std::ios_base::binary );
   boost::archive::binary_iarchive inputArchive( inArchiveFile );
   inputArchive >> postsolveStorage;
   inArchiveFile.close();
   Solution<double> reduced_solution{};
   Solution<double> original_solution{};
   Postsolve<double> postsolve{msg, num};

   REQUIRE( postsolve.undo( reduced_solution, original_solution, postsolveStorage) ==
            PostsolveStatus::kOk );
   papilo::Vec<double> values = original_solution.primal;
   papilo::Vec<double> expected_values{ -11, -5, -5 };
   REQUIRE( values == expected_values );
}

TEST_CASE( "finding-the-right-value-in-postsolve-for-a-column-fixed-pos-inf",
           "[core]" )
{

   const Num<double> num{};
   Message msg{};
   PostsolveStorage<double> postsolveStorage{};

   std::ifstream inArchiveFile( "./resources/dual_fix_pos_inf.postsolve",
                                std::ios_base::binary );
   boost::archive::binary_iarchive inputArchive( inArchiveFile );
   inputArchive >> postsolveStorage;
   inArchiveFile.close();
   Solution<double> reduced_solution{};
   Solution<double> original_solution{};
   Postsolve<double> postsolve{msg, num};

   REQUIRE( postsolve.undo( reduced_solution, original_solution, postsolveStorage ) ==
            PostsolveStatus::kOk );
   papilo::Vec<double> values = original_solution.primal;
   papilo::Vec<double> expected_values{ 13, 9, -5, -2.5 };
   REQUIRE( values == expected_values );

}
