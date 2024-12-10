/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*               This file is part of the program and library                */
/*    PaPILO --- Parallel Presolve for Integer and Linear Optimization       */
/*                                                                           */
/* Copyright (C) 2020-2024 Zuse Institute Berlin (ZIB)                       */
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

using namespace papilo;

Problem<double>
setupMatrixForCliqueMerging();

TEST_CASE( "clique-merging-basic", "[presolve]" )
{
    
   double time = 0.0;
   int cause = -1;
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
    const std::vector<RowFlags> rowFlags = matrix.getRowFlags();/*
    for( int i = 0; i < 6; ++i ){
        Cliques.push_back(i);
        REQUIRE( rowFlags[i].test(RowFlag::kClique) );
    }
    int col = 0;
    int neighbournumber = 1;
    int neigh = presolvingMethod.getNeighbour( matrix, col, neighbournumber );
    REQUIRE( neigh == 2 );
    
    int clique = 0;
    newClique = presolvingMethod.greedyClique( matrix, Cliques[clique] );
    REQUIRE( newClique.size() == 1 );
    REQUIRE( presolvingMethod.greedyClique( matrix, 3 ).size() == 0 );
    REQUIRE( presolvingMethod.greedyClique( matrix, 4 ).size() == 0 );
    REQUIRE( presolvingMethod.greedyClique( matrix, 5 ).size() == 0 );
    REQUIRE( presolvingMethod.greedyClique( matrix, 1 ).size() == 1 );
    REQUIRE( presolvingMethod.greedyClique( matrix, 2 ).size() == 1 );

    int neighnumber = presolvingMethod.getNeighbourhoodSize(matrix, col);
    REQUIRE( neighnumber == 2);
    REQUIRE( presolvingMethod.getNeighbourhoodSize(matrix, 3) == 3 );
    REQUIRE( presolvingMethod.getNeighbourhoodSize(matrix, 0) == 2 );
    REQUIRE( presolvingMethod.isNeighbour(matrix, 0, 2));
    REQUIRE( !presolvingMethod.isNeighbour(matrix, 0, 3));
    REQUIRE( presolvingMethod.isCovered( matrix, 1, newClique, 0));
    REQUIRE( !presolvingMethod.isCovered( matrix, 3, newClique, 0));

    REQUIRE( !presolvingMethod.isCovered( matrix, 1, noClique, 0));
    REQUIRE( !presolvingMethod.isCovered( matrix, 2, noClique, 0));
    REQUIRE( !presolvingMethod.isCovered( matrix, 4, noClique, 0));
    REQUIRE( !presolvingMethod.isCovered( matrix, 5, noClique, 0));
    REQUIRE( !presolvingMethod.isCovered( matrix, 3, noClique, 0));

    REQUIRE( !presolvingMethod.isCovered( matrix, 0, noClique, 1));
    REQUIRE( !presolvingMethod.isCovered( matrix, 2, noClique, 1));
    REQUIRE( !presolvingMethod.isCovered( matrix, 1, noClique, 2));
    REQUIRE( !presolvingMethod.isCovered( matrix, 4, noClique, 1));
    REQUIRE( !presolvingMethod.isCovered( matrix, 5, noClique, 1));
    REQUIRE( !presolvingMethod.isCovered( matrix, 3, noClique, 1));
    for( int i = 0; i < 6; ++i )
    {
        for( int j = 0; j < 6; ++j )
        {
            if( i == j )
            {
                //REQUIRE( presolvingMethod.isCovered( matrix, i, noClique, j));
            }
            else
            {
                //REQUIRE( !presolvingMethod.isCovered( matrix, i, noClique, j));
            }
        }
    }*/
    
   PresolveStatus presolveStatus =
       presolvingMethod.execute( problem, problemUpdate, num, reductions, t, cause);

   REQUIRE( presolveStatus == PresolveStatus::kReduced );
   //REQUIRE( reductions.size() == 9 );
   /*
    REQUIRE( reductions.getReduction( 0 ).row == 0 );
    REQUIRE( reductions.getReduction( 0 ).col == RowReduction::LOCKED );

    REQUIRE( reductions.getReduction( 1 ).row == 1 );
    REQUIRE( reductions.getReduction( 1 ).col == RowReduction::LOCKED );

    REQUIRE( reductions.getReduction( 2 ).row == 2 );
    REQUIRE( reductions.getReduction( 2 ).col == RowReduction::LOCKED );

    REQUIRE( reductions.getReduction( 3 ).row == ColReduction::BOUNDS_LOCKED );
    REQUIRE( reductions.getReduction( 3 ).col == 2 );

    REQUIRE( reductions.getReduction( 4 ).row == ColReduction::BOUNDS_LOCKED );
    REQUIRE( reductions.getReduction( 4 ).col == 0 );

    REQUIRE( reductions.getReduction( 5 ).row == ColReduction::BOUNDS_LOCKED );
    REQUIRE( reductions.getReduction( 5 ).col == 1 );

    REQUIRE( reductions.getReduction( 6 ).row == 1 );
    REQUIRE( reductions.getReduction( 6 ).col == RowReduction::REDUNDANT );

    REQUIRE( reductions.getReduction( 7 ).row == 2 );
    REQUIRE( reductions.getReduction( 7 ).col == RowReduction::REDUNDANT );

    REQUIRE( reductions.getReduction( 8 ).row == 0 );
    REQUIRE( reductions.getReduction( 8 ).col == 2 );
    REQUIRE( reductions.getReduction( 8 ).newval == 1.0 );*/
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
