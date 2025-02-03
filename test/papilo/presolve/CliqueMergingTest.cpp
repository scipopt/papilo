/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*               This file is part of the program and library                */
/*    PaPILO --- Parallel Presolve for Integer and Linear Optimization       */
/*                                                                           */
/* Copyright (C) 2020-2025 Zuse Institute Berlin (ZIB)                       */
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
#include "papilo/core/PresolveMethod.hpp"
#include "papilo/core/Problem.hpp"
#include "papilo/core/ProblemBuilder.hpp"
#include "papilo/external/catch/catch_amalgamated.hpp"

#include <set>
#include <unordered_map>

using namespace papilo;

Problem<double>
setupMatrixForCliqueMerging();

TEST_CASE( "clique-merging-basic", "[presolve]" )
{

   CliqueMerging<double> presolvingMethod{};

   double time = 0.0;
   Timer t{ time };
   Problem<double> problem = setupMatrixForCliqueMerging();
   Statistics statistics{};
   PresolveOptions presolveOptions{};
   PostsolveStorage<double> postsolve =
       PostsolveStorage<double>( problem, {}, presolveOptions );
   ProblemUpdate<double> problemUpdate( problem, postsolve, statistics,
                                        presolveOptions, {}, {} );

   Reductions<double> reductions{};

   //   Vec<int> newClique, noClique, Cliques;
//   const auto& matrix = problem.getConstraintMatrix();
   presolvingMethod.setParameters( 1000000, 100000, 100, 10000 );
   int cause = -1;
   PresolveStatus status = presolvingMethod.execute(
       problem, problemUpdate, { }, reductions, t, cause );

   REQUIRE( status == PresolveStatus::kReduced );
    REQUIRE( reductions.size() <= 300 );
#ifdef PAPILO_TBB
/*
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
    REQUIRE( reductions.getReduction(6).row == ColReduction::LOCKED );
    REQUIRE( reductions.getReduction(6).col == 3 );
    REQUIRE( reductions.getReduction(7).row == ColReduction::BOUNDS_LOCKED );
    REQUIRE( reductions.getReduction(7).col == 3 );*/
    REQUIRE( reductions.getReduction(0).row == ColReduction::LOCKED );
    REQUIRE( reductions.getReduction(0).col == 2 );
    REQUIRE( reductions.getReduction(1).row == ColReduction::BOUNDS_LOCKED );
    REQUIRE( reductions.getReduction(1).col == 2 );
    REQUIRE( reductions.getReduction(2).row == ColReduction::LOCKED );
    REQUIRE( reductions.getReduction(2).col == 3 );
    REQUIRE( reductions.getReduction(3).row == ColReduction::BOUNDS_LOCKED );
    REQUIRE( reductions.getReduction(3).col == 3 );
    REQUIRE( reductions.getReduction(4).row == ColReduction::LOCKED );
    REQUIRE( reductions.getReduction(4).col == 0 );
    REQUIRE( reductions.getReduction(5).row == ColReduction::BOUNDS_LOCKED );
    REQUIRE( reductions.getReduction(5).col == 0 );
    REQUIRE( reductions.getReduction(6).row == ColReduction::LOCKED );
    REQUIRE( reductions.getReduction(6).col == 1 );
    REQUIRE( reductions.getReduction(7).row == ColReduction::BOUNDS_LOCKED );
    REQUIRE( reductions.getReduction(7).col == 1 );
#else
    REQUIRE( reductions.getReduction(0).row == ColReduction::LOCKED );
    REQUIRE( reductions.getReduction(0).col == 2 );
    REQUIRE( reductions.getReduction(1).row == ColReduction::BOUNDS_LOCKED );
    REQUIRE( reductions.getReduction(1).col == 2 );
    REQUIRE( reductions.getReduction(2).row == ColReduction::LOCKED );
    REQUIRE( reductions.getReduction(2).col == 3 );
    REQUIRE( reductions.getReduction(3).row == ColReduction::BOUNDS_LOCKED );
    REQUIRE( reductions.getReduction(3).col == 3 );
    REQUIRE( reductions.getReduction(4).row == ColReduction::LOCKED );
    REQUIRE( reductions.getReduction(4).col == 0 );
    REQUIRE( reductions.getReduction(5).row == ColReduction::BOUNDS_LOCKED );
    REQUIRE( reductions.getReduction(5).col == 0 );
    REQUIRE( reductions.getReduction(6).row == ColReduction::LOCKED );
    REQUIRE( reductions.getReduction(6).col == 1 );
    REQUIRE( reductions.getReduction(7).row == ColReduction::BOUNDS_LOCKED );
    REQUIRE( reductions.getReduction(7).col == 1 );
#endif
    REQUIRE( reductions.getReduction(8).row == 1 );
    REQUIRE( reductions.getReduction(8).col == RowReduction::LOCKED );
    REQUIRE( reductions.getReduction(9).row == 2 );
    REQUIRE( reductions.getReduction(9).col == RowReduction::LOCKED );
    REQUIRE( reductions.getReduction(10).row == 5 );
    REQUIRE( reductions.getReduction(10).col == RowReduction::LOCKED );
    REQUIRE( reductions.getReduction(11).row == 6 );
    REQUIRE( reductions.getReduction(11).col == RowReduction::LOCKED );
    REQUIRE( reductions.getReduction(11).row == 0 );
    REQUIRE( reductions.getReduction(11).col == RowReduction::LOCKED );
    REQUIRE( reductions.getReduction(13).row == 1 );
    REQUIRE( reductions.getReduction(13).col == RowReduction::REDUNDANT );
    REQUIRE( reductions.getReduction(14).row == 2 );
    REQUIRE( reductions.getReduction(14).col == RowReduction::REDUNDANT );
    REQUIRE( reductions.getReduction(15).row == 5 );
    REQUIRE( reductions.getReduction(15).col == RowReduction::REDUNDANT );
    REQUIRE( reductions.getReduction(16).row == 6 );
    REQUIRE( reductions.getReduction(16).col == RowReduction::REDUNDANT );


//First Lock Col, Lock Colbounds of each index that already exists.
//Then Lock col, lock col bounds of all new vertices
// then lock row of all covered cliques
// then lock row of original clique
// then change row
// then mark rows redundant
//row 0, with indices 0, 1 is original clique, indices 2,3 get added.
//Then rows 1,2 5, 6 are covered
   ////TODO: please encode the reductions like this.
   //   REQUIRE( reductions.size() <= 3 );
   //   REQUIRE( reductions.getReduction(0).row ==  1 );
   //   REQUIRE( reductions.getReduction(1).col == RowReduction::REDUNDANT );
   //   REQUIRE( reductions.getReduction(1).row ==  1 );
   //   REQUIRE( reductions.getReduction(2).col == RowReduction::REDUNDANT );
}

TEST_CASE( "clique-merging-functions", "[presolve]" )
{
   double time = 0.0;
   Timer t{ time };
   Problem<double> problem = setupMatrixForCliqueMerging();
   Statistics statistics{};
   PresolveOptions presolveOptions{};
   PostsolveStorage<double> postsolve =
       PostsolveStorage<double>( problem, { }, presolveOptions );
   ProblemUpdate<double> problemUpdate( problem, postsolve, statistics,
                                        presolveOptions, { }, { } );

   CliqueMerging<double> presolvingMethod{};

   Vec<int> newClique;
   Vec<int> noClique;
   Vec<int> Cliques;
   const auto& matrix = problem.getConstraintMatrix();
   const std::vector<RowFlags> rowFlags = matrix.getRowFlags();

   std::set<int> imaginaryclique1 = { 0, 1, 2, 3, 4 };
   std::set<int> imaginaryclique2 = { 1, 2 };
   //TODO: this function should not be visible (only within the class) -> test on presolveMethod
//   REQUIRE( presolvingMethod.isCovered( matrix, 0, imaginaryclique1 ) );
//   REQUIRE( !presolvingMethod.isCovered( matrix, 0, imaginaryclique2 ) );

   std::set<std::pair<int, int>> edges = { { 0, 1 }, { 1, 0 }, { 0, 2 },
                                           { 2, 0 }, { 1, 2 }, { 2, 1 } };
   std::unordered_map<int, std::set<int>> Neighbourhoodlists;
   std::set<int> r1 = { 1, 2 };
   std::set<int> r2 = { 0, 2 };
   std::set<int> r3 = { 0, 1 };
   Neighbourhoodlists[0] = r1;
   Neighbourhoodlists[1] = r2;
   Neighbourhoodlists[2] = r3;
   //TODO: this function should not be visible (only within the class) -> test on presolveMethod
//   auto result =
//       presolvingMethod.greedyClique( matrix, edges, Neighbourhoodlists, 0 );
//   std::set<int> r4 = { 0, 1, 2 };
//   std::vector<int> r5 = { 2 };
//   REQUIRE( status == PresolveStatus::kUnchanged);
   //TODO: in contrast to general coding expectations should not be encoded in variables (bad practice). You can write down the expectations immedaitely.
//   std::pair<std::set<int>, std::vector<int>> expectedresult = { r4, r5 };
   //   REQUIRE( result == expectedresult );
   //
   //   REQUIRE( problem.is_clique( problem.getConstraintMatrix(), 1, num ) );
   //   REQUIRE( problem.is_clique( problem.getConstraintMatrix(), 2, num ) );
   //   REQUIRE( problem.is_clique( problem.getConstraintMatrix(), 3, num ) );
   //   REQUIRE( problem.is_clique( problem.getConstraintMatrix(), 0, num ) );
   //   REQUIRE( problem.is_clique( problem.getConstraintMatrix(), 4, num ) );
   //   REQUIRE( problem.is_clique( problem.getConstraintMatrix(), 5, num ) );
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
   // G: x + y + z + a <= 1

   Vec<std::string> columnNames{ "x", "y", "z", "a", "b", "c" };

   Vec<double> coefficients{ -1.0, -1.0, -1.0, -1.0, -1.0, -1.0 };
   Vec<double> upperBounds{ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };
   Vec<double> lowerBounds{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
   Vec<uint8_t> isIntegral{ 1, 1, 1, 1, 1, 1 };

   Vec<double> rhs{ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };
   Vec<std::string> rowNames{ "A", "B", "C", "D", "E", "F", "G" };
   Vec<uint8_t> lhsInfinity{ 1, 1, 1, 1, 1, 1, 1 };
   Vec<uint8_t> rhsInfinity{ 0, 0, 0, 0, 0, 0, 0 };
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

       std::tuple<int, int, double>{ 6, 0, 1.0 },
       std::tuple<int, int, double>{ 6, 1, 1.0 },
       std::tuple<int, int, double>{ 6, 2, 1.0 },
       std::tuple<int, int, double>{ 6, 3, 1.0 },

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
   pb.setProblemName( "matrix for testing Clique Merging" );
   Problem<double> problem = pb.build();
   return problem;
}
