/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*               This file is part of the program and library                */
/*    PaPILO --- Parallel Presolve for Integer and Linear Optimization       */
/*                                                                           */
/* Copyright (C) 2025 Zuse Institute Berlin (ZIB)                       */
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

#include "papilo/presolvers/CliqueMerging.hpp"
#include "papilo/external/catch/catch_amalgamated.hpp"
#include "papilo/core/PresolveMethod.hpp"
#include "papilo/core/Problem.hpp"
#include "papilo/core/ProblemBuilder.hpp"

#include <map>
#include <set>
#include <vector>

using namespace papilo;

Problem<double>
setupMatrixForCliqueMerging();

TEST_CASE( "clique-merging-basic", "[presolve]" )
{
   int cause = -1;
   double time = 0.0;
   Timer t{time};
   Num<double> num{};
   Message msg{};
   Problem<double> problem = setupMatrixForCliqueMerging();
   Statistics statistics{};
   PresolveOptions presolveOptions{};
   PostsolveStorage<double> postsolve =
       PostsolveStorage<double>( problem, num, presolveOptions );
   ProblemUpdate<double> problemUpdate( problem, postsolve, statistics,
                                        presolveOptions, num, msg );

   CliqueMerging<double> presolvingMethod{};
   Reductions<double> reductions{};

    Vec<int> newClique;
    Vec<int> noClique;
    Vec<int> Cliques;
    const auto& matrix = problem.getConstraintMatrix();
    const std::vector<RowFlags> rowFlags = matrix.getRowFlags();
   PresolveStatus presolveStatus =
       presolvingMethod.execute( problem, problemUpdate, num, reductions, t, cause);

   REQUIRE( presolveStatus == PresolveStatus::kReduced );
    bool correctness1 = false;
    bool correctness2 = false;
    for( unsigned int red = 0; red < reductions.size(); ++ red )
    {
        if( reductions.getReduction(red).row ==  1 && reductions.getReduction(red).col == RowReduction::REDUNDANT )
        {
            correctness1 = true;
        }
        if( reductions.getReduction(red).row ==  1 && reductions.getReduction(red).col == RowReduction::REDUNDANT )
        {
            correctness2 = true;
        }
    }
    REQUIRE( correctness1 );
    REQUIRE( correctness2 );
}


TEST_CASE( "clique-merging-functions", "[presolve]" )
{   
    double time = 0.0;
   Timer t{time};
   Num<double> num{};
   Message msg{};
   Problem<double> problem = setupMatrixForCliqueMerging();
   Statistics statistics{};
   PresolveOptions presolveOptions{};
   PostsolveStorage<double> postsolve =
       PostsolveStorage<double>( problem, num, presolveOptions );
   ProblemUpdate<double> problemUpdate( problem, postsolve, statistics,
                                        presolveOptions, num, msg );

   CliqueMerging<double> presolvingMethod{};

    Vec<int> newClique;
    Vec<int> noClique;
    Vec<int> Cliques;
    const auto& matrix = problem.getConstraintMatrix();
    const std::vector<RowFlags> rowFlags = matrix.getRowFlags();
   
    std::set<int> imaginaryclique1 = {0,1,2,3,4};
    std::set<int> imaginaryclique2 = {1,2};
    REQUIRE( presolvingMethod.isCovered( matrix, 0, imaginaryclique1 ) );
    REQUIRE( !presolvingMethod.isCovered( matrix, 0, imaginaryclique2 ) );

    std::set<std::pair<int, int>> edges = {{0,1},{1,0},{0,2},{2,0},{1,2},{2,1}};
    std::map<int, std::set<int>> Neighbourhoodlists;
    std::set<int> r1 = {1,2};
    std::set<int> r2 = {0,2};
    std::set<int> r3 = {0,1};
    Neighbourhoodlists.emplace(0, r1);
    Neighbourhoodlists.emplace(1, r2);
    Neighbourhoodlists.emplace(2, r3);
    auto result = presolvingMethod.greedyClique( matrix, edges,
    Neighbourhoodlists, 0 );
    std::set<int> r4 = {0,1,2};
    std::vector<int> r5 = {2};
    std::pair<std::set<int>,std::vector<int>> expectedresult = {r4,r5};
    REQUIRE( result == expectedresult );
    
   REQUIRE( problem.is_clique( problem.getConstraintMatrix(), 1,  num) );
   REQUIRE( problem.is_clique( problem.getConstraintMatrix(), 2,  num) );
   REQUIRE( problem.is_clique( problem.getConstraintMatrix(), 3,  num) );
   REQUIRE( problem.is_clique( problem.getConstraintMatrix(), 0,  num) );
   REQUIRE( problem.is_clique( problem.getConstraintMatrix(), 4,  num) );
   REQUIRE( problem.is_clique( problem.getConstraintMatrix(), 5,  num) );
}

Problem<double>
setupMatrixForCliqueMerging()
{
   // Clique x, y, z, no Clique a, b, c
   // min -x -y -z -a -b -c
   // A: x + y <= 1
   // B: x + z <= 1
   // C: y + z <= 1
   // D: a + b <= 1
   // E: a + c <= 1
   // F: a + z <= 1

   Vec<std::string> columnNames{ "x", "y", "z", "a", "b", "c" };

   Vec<double> coefficients{ -1.0, -1.0, -1.0, -1.0, -1.0, -1.0 };
   Vec<double> upperBounds{ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };
   Vec<double> lowerBounds{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
   Vec<uint8_t> isIntegral{ 1, 1, 1, 1, 1, 1 };

   Vec<double> rhs{ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };
   Vec<std::string> rowNames{ "A", "B", "C", "D", "E", "F" };
   Vec<uint8_t> lhsInfinity{ 1, 1, 1, 1, 1, 1 };
   Vec<uint8_t> rhsInfinity{ 0, 0, 0, 0, 0, 0 };
   Vec<std::tuple<int, int, double>> entries{
       std::tuple<int, int, double>{ 0, 0, 1.0 },
       std::tuple<int, int, double>{ 0, 1, 1.0 },

       std::tuple<int, int, double>{ 1, 0, 1.0 },
       std::tuple<int, int, double>{ 1, 2, 1.0 },

       std::tuple<int, int, double>{ 2, 1, 1.0 },
       std::tuple<int, int, double>{ 2, 2, 1.0 },

       std::tuple<int, int, double>{ 3, 3, 1.0 },
       std::tuple<int, int, double>{ 3, 4, 1.0 },

       std::tuple<int, int, double>{ 4, 3, 1.0 },
       std::tuple<int, int, double>{ 4, 5, 1.0 },

       std::tuple<int, int, double>{ 5, 3, 1.0 },
       std::tuple<int, int, double>{ 5, 2, 1.0 },
   };

   ProblemBuilder<double> pb;
   pb.reserve( (int) entries.size(), (int) rowNames.size(), (int) columnNames.size() );
   pb.setNumRows( (int) rowNames.size() );
   pb.setNumCols( (int) columnNames.size() );
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
   pb.setProblemName( "matrix x dom y" );
   Problem<double> problem = pb.build();
   return problem;
}
