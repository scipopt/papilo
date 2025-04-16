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

#include "papilo/presolvers/CliqueMerging.hpp"
#include "papilo/core/PresolveMethod.hpp"
#include "papilo/core/Problem.hpp"
#include "papilo/core/ProblemBuilder.hpp"
#include "papilo/external/catch/catch_amalgamated.hpp"

#include <set>
#include <unordered_map>

using namespace papilo;

Problem<double>
setupSmallMatrixForCliqueMerging();

Problem<double>
setupSmallerMatrixForCliqueMerging();

TEST_CASE( "clique-merging-basic", "[presolve]" )
{

   CliqueMerging<double> presolvingMethod{};

   double time = 0.0;
   Timer t{ time };
   Problem<double> problem = setupSmallMatrixForCliqueMerging();
   Statistics statistics{};
   PresolveOptions presolveOptions{};
   PostsolveStorage<double> postsolve =
       PostsolveStorage<double>( problem, {}, presolveOptions );
   ProblemUpdate<double> problemUpdate( problem, postsolve, statistics,
                                        presolveOptions, {}, {} );

   Reductions<double> reductions{};
   presolvingMethod.setParameters( 1000000, 100000, 100, 10000 );
   int cause = -1;
   PresolveStatus status = presolvingMethod.execute(
       problem, problemUpdate, { }, reductions, t, cause );
    /*
    for( int i = 0; (unsigned int)i < reductions.size(); ++i)
    {
        std::cout << "\n";
        std::cout << reductions.getReduction(i).row;
        std::cout << " ";
        std::cout << reductions.getReduction(i).col;
    }*/

   REQUIRE( status == PresolveStatus::kReduced );
#ifdef PAPILO_TBB
    //large number of reductions due to the grouping of transactions into O(nlog(n)) by the clique merging presolver
    REQUIRE( reductions.size() == 72 );
    
    REQUIRE( reductions.getReduction(0).row == ColReduction::LOCKED );
    REQUIRE( reductions.getReduction(0).col == 1 );

    REQUIRE( reductions.getReduction(1).row == ColReduction::BOUNDS_LOCKED );
    REQUIRE( reductions.getReduction(1).col == 1 );
    
    REQUIRE( reductions.getReduction(2).row == ColReduction::LOCKED );
    REQUIRE( reductions.getReduction(2).col == 2 );

    REQUIRE( reductions.getReduction(3).row == ColReduction::BOUNDS_LOCKED );
    REQUIRE( reductions.getReduction(3).col == 2 );
    
    REQUIRE( reductions.getReduction(4).row == ColReduction::LOCKED );
    REQUIRE( reductions.getReduction(4).col == 0 );

    REQUIRE( reductions.getReduction(5).row == ColReduction::BOUNDS_LOCKED );
    REQUIRE( reductions.getReduction(5).col == 0 );

    REQUIRE( reductions.getReduction(6).row == 0 );
    REQUIRE( reductions.getReduction(6).col == RowReduction::LOCKED );
    
    REQUIRE( reductions.getReduction(7).row == 1 );
    REQUIRE( reductions.getReduction(7).col == RowReduction::LOCKED );
    
    REQUIRE( reductions.getReduction(8).row == 2 );
    REQUIRE( reductions.getReduction(8).col == RowReduction::LOCKED );
        
    REQUIRE( reductions.getReduction(9).row == ColReduction::LOCKED );
    REQUIRE( reductions.getReduction(9).col == 0 );

    REQUIRE( reductions.getReduction(10).row == ColReduction::BOUNDS_LOCKED );
    REQUIRE( reductions.getReduction(10).col == 0 );
    
    REQUIRE( reductions.getReduction(11).row == ColReduction::LOCKED );
    REQUIRE( reductions.getReduction(11).col == 2 );

    REQUIRE( reductions.getReduction(12).row == ColReduction::BOUNDS_LOCKED );
    REQUIRE( reductions.getReduction(12).col == 2 );
    
    REQUIRE( reductions.getReduction(13).row == ColReduction::LOCKED );
    REQUIRE( reductions.getReduction(13).col == 1 );

    REQUIRE( reductions.getReduction(14).row == ColReduction::BOUNDS_LOCKED );
    REQUIRE( reductions.getReduction(14).col == 1 );
    
#else 
    REQUIRE( reductions.size() == 12 );
    
    REQUIRE( reductions.getReduction(0).row == ColReduction::LOCKED );
    REQUIRE( reductions.getReduction(0).col == 0 );

    REQUIRE( reductions.getReduction(1).row == ColReduction::BOUNDS_LOCKED );
    REQUIRE( reductions.getReduction(1).col == 0 );

    REQUIRE( reductions.getReduction(2).row == ColReduction::LOCKED );
    REQUIRE( reductions.getReduction(2).col == 1 );

    REQUIRE( reductions.getReduction(3).row == ColReduction::BOUNDS_LOCKED );
    REQUIRE( reductions.getReduction(3).col == 1 );

    REQUIRE( reductions.getReduction(4).row == ColReduction::LOCKED );
    REQUIRE( reductions.getReduction(4).col == 2 );

    REQUIRE( reductions.getReduction(5).row == ColReduction::BOUNDS_LOCKED );
    REQUIRE( reductions.getReduction(5).col == 2 );

    REQUIRE( reductions.getReduction(6).row == 1 );
    REQUIRE( reductions.getReduction(6).col == RowReduction::LOCKED );

    REQUIRE( reductions.getReduction(7).row == 2 );
    REQUIRE( reductions.getReduction(7).col == RowReduction::LOCKED );

    REQUIRE( reductions.getReduction(8).row == 0 );
    REQUIRE( reductions.getReduction(8).col == RowReduction::LOCKED );

    REQUIRE( reductions.getReduction(9).row == 0 );
    REQUIRE( reductions.getReduction(9).col == 2 );

    REQUIRE( reductions.getReduction(10).row == 1 );
    REQUIRE( reductions.getReduction(10).col == RowReduction::REDUNDANT );

    REQUIRE( reductions.getReduction(11).row == 2 );
    REQUIRE( reductions.getReduction(11).col == RowReduction::REDUNDANT );
#endif
}

TEST_CASE( "clique-merging-cover", "[presolve]" )
{

   CliqueMerging<double> presolvingMethod{};

   double time = 0.0;
   Timer t{ time };
   Problem<double> problem = setupSmallerMatrixForCliqueMerging();
   Statistics statistics{};
   PresolveOptions presolveOptions{};
   PostsolveStorage<double> postsolve =
       PostsolveStorage<double>( problem, {}, presolveOptions );
   ProblemUpdate<double> problemUpdate( problem, postsolve, statistics,
                                        presolveOptions, {}, {} );

   Reductions<double> reductions{};
   presolvingMethod.setParameters( 1000000, 100000, 100, 10000 );
   int cause = -1;
   PresolveStatus status = presolvingMethod.execute(
       problem, problemUpdate, { }, reductions, t, cause );
/*
    for( int i = 0; (unsigned int)i < reductions.size(); ++i)
    {
        std::cout << "\n";
        std::cout << reductions.getReduction(i).row;
        std::cout << " ";
        std::cout << reductions.getReduction(i).col;
    }*/

   REQUIRE( status == PresolveStatus::kReduced );
#ifdef PAPILO_TBB
    //large number of reductions due to the grouping of transactions into O(nlog(n)) by the clique merging presolver
    REQUIRE( reductions.size() == 38 );

    REQUIRE( reductions.getReduction(0).row == ColReduction::LOCKED );
    REQUIRE( reductions.getReduction(0).col == 0 );

    REQUIRE( reductions.getReduction(1).row == ColReduction::BOUNDS_LOCKED );
    REQUIRE( reductions.getReduction(1).col == 0 );

    REQUIRE( reductions.getReduction(2).row == ColReduction::LOCKED );
    REQUIRE( reductions.getReduction(2).col == 1 );

    REQUIRE( reductions.getReduction(3).row == ColReduction::BOUNDS_LOCKED );
    REQUIRE( reductions.getReduction(3).col == 1 );

    REQUIRE( reductions.getReduction(4).row == ColReduction::LOCKED );
    REQUIRE( reductions.getReduction(4).col == 2 );

    REQUIRE( reductions.getReduction(5).row == ColReduction::BOUNDS_LOCKED );
    REQUIRE( reductions.getReduction(5).col == 2 );

    REQUIRE( reductions.getReduction(6).row == 1 );
    REQUIRE( reductions.getReduction(6).col == RowReduction::LOCKED );

    REQUIRE( reductions.getReduction(7).row == 0 );
    REQUIRE( reductions.getReduction(7).col == RowReduction::LOCKED );

    REQUIRE( reductions.getReduction(8).row == ColReduction::LOCKED );
    REQUIRE( reductions.getReduction(8).col == 0 );

    REQUIRE( reductions.getReduction(9).row == ColReduction::BOUNDS_LOCKED );
    REQUIRE( reductions.getReduction(9).col == 0 );

    REQUIRE( reductions.getReduction(10).row == ColReduction::LOCKED );
    REQUIRE( reductions.getReduction(10).col == 1 );

    REQUIRE( reductions.getReduction(11).row == ColReduction::BOUNDS_LOCKED );
    REQUIRE( reductions.getReduction(11).col == 1 );

    REQUIRE( reductions.getReduction(12).row == ColReduction::LOCKED );
    REQUIRE( reductions.getReduction(12).col == 2 );

    REQUIRE( reductions.getReduction(13).row == ColReduction::BOUNDS_LOCKED );
    REQUIRE( reductions.getReduction(13).col == 2 );

#else
    REQUIRE( reductions.size() == 9 );
    
    REQUIRE( reductions.getReduction(0).row == ColReduction::LOCKED );
    REQUIRE( reductions.getReduction(0).col == 0 );

    REQUIRE( reductions.getReduction(1).row == ColReduction::BOUNDS_LOCKED );
    REQUIRE( reductions.getReduction(1).col == 0 );

    REQUIRE( reductions.getReduction(2).row == ColReduction::LOCKED );
    REQUIRE( reductions.getReduction(2).col == 1 );

    REQUIRE( reductions.getReduction(3).row == ColReduction::BOUNDS_LOCKED );
    REQUIRE( reductions.getReduction(3).col == 1 );

    REQUIRE( reductions.getReduction(4).row == ColReduction::LOCKED );
    REQUIRE( reductions.getReduction(4).col == 2 );

    REQUIRE( reductions.getReduction(5).row == ColReduction::BOUNDS_LOCKED );
    REQUIRE( reductions.getReduction(5).col == 2 );

    REQUIRE( reductions.getReduction(6).row == 1 );
    REQUIRE( reductions.getReduction(6).col == RowReduction::LOCKED );

    REQUIRE( reductions.getReduction(7).row == 0 );
    REQUIRE( reductions.getReduction(7).col == RowReduction::LOCKED );

    REQUIRE( reductions.getReduction(8).row == 1 );
    REQUIRE( reductions.getReduction(8).col == RowReduction::REDUNDANT );

#endif
}

Problem<double>
setupSmallMatrixForCliqueMerging()
{
   // Clique x, y, z
   // min -x -y -z
   // A: x + y <= 1
   // B: x + z <= 1
   // C: y + z <= 1

   Vec<std::string> columnNames{ "x", "y", "z" };

   Vec<double> coefficients{ -1.0, -1.0, -1.0 };
   Vec<double> upperBounds{ 1.0, 1.0, 1.0 };
   Vec<double> lowerBounds{ 0.0, 0.0, 0.0 };
   Vec<uint8_t> isIntegral{ 1, 1, 1 };

   Vec<double> rhs{ 1.0, 1.0, 1.0 };
   Vec<std::string> rowNames{ "A", "B", "C" };
   Vec<uint8_t> lhsInfinity{ 1, 1, 1 };
   Vec<uint8_t> rhsInfinity{ 0, 0, 0 };
   Vec<std::tuple<int, int, double>> entries{
       std::tuple<int, int, double>{ 0, 0, 1.0 },
       std::tuple<int, int, double>{ 0, 1, 1.0 },

       std::tuple<int, int, double>{ 1, 0, 1.0 },
       std::tuple<int, int, double>{ 1, 2, 1.0 },

       std::tuple<int, int, double>{ 2, 1, 1.0 },
       std::tuple<int, int, double>{ 2, 2, 1.0 }
   };

   ProblemBuilder<double> pb;
   pb.reserve( (int)entries.size(), (int)rowNames.size(),
               (int)columnNames.size() );
   pb.setNumRows( (int)rowNames.size() );
   pb.setNumCols( (int)columnNames.size() );
   pb.setColUbAll( upperBounds );
   pb.setColLbAll( lowerBounds );
   pb.setObjAll( coefficients );
   pb.setObjOffset( 0.0 );
   pb.setColIntegralAll( isIntegral );
   pb.setRowRhsAll( rhs );
   pb.setRowLhsInfAll( lhsInfinity );
   pb.setRowRhsInfAll( rhsInfinity );
   pb.addEntryAll( entries );
   pb.setColNameAll( columnNames );
   pb.setProblemName( "small matrix for testing Clique Merging" );
   Problem<double> problem = pb.build();
   return problem;
}

Problem<double>
setupSmallerMatrixForCliqueMerging()
{
   // Clique x, y, z
   // min -x -y -z
   // A: x + y <= 1
   // B: x + z <= 1
   // C: y + z <= 1

   Vec<std::string> columnNames{ "x", "y", "z" };

   Vec<double> coefficients{ -1.0, -1.0, -1.0 };
   Vec<double> upperBounds{ 1.0, 1.0, 1.0 };
   Vec<double> lowerBounds{ 0.0, 0.0, 0.0 };
   Vec<uint8_t> isIntegral{ 1, 1, 1 };

   Vec<double> rhs{ 1.0, 1.0 };
   Vec<std::string> rowNames{ "A", "B" };
   Vec<uint8_t> lhsInfinity{ 1, 1 };
   Vec<uint8_t> rhsInfinity{ 0, 0 };
   Vec<std::tuple<int, int, double>> entries{
       std::tuple<int, int, double>{ 0, 0, 1.0 },
       std::tuple<int, int, double>{ 0, 1, 1.0 },
       std::tuple<int, int, double>{ 0, 2, 1.0 },

       std::tuple<int, int, double>{ 1, 0, 1.0 },
       std::tuple<int, int, double>{ 1, 1, 1.0 },
   };

   ProblemBuilder<double> pb;
   pb.reserve( (int)entries.size(), (int)rowNames.size(),
               (int)columnNames.size() );
   pb.setNumRows( (int)rowNames.size() );
   pb.setNumCols( (int)columnNames.size() );
   pb.setColUbAll( upperBounds );
   pb.setColLbAll( lowerBounds );
   pb.setObjAll( coefficients );
   pb.setObjOffset( 0.0 );
   pb.setColIntegralAll( isIntegral );
   pb.setRowRhsAll( rhs );
   pb.setRowLhsInfAll( lhsInfinity );
   pb.setRowRhsInfAll( rhsInfinity );
   pb.addEntryAll( entries );
   pb.setColNameAll( columnNames );
   pb.setProblemName( "smaller matrix for testing Clique Merging" );
   Problem<double> problem = pb.build();
   return problem;
}
