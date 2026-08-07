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

#include "papilo/core/ConstraintMatrix.hpp"
#include "papilo/core/Problem.hpp"
#include "papilo/core/ProblemBuilder.hpp"
#include "papilo/external/catch/catch_amalgamated.hpp"

using namespace papilo;

static Vec<Triplet<double>>
setupTriplets()
{
   // 1  1  0  0  0
   // 0  0  1  1  1
   return { Triplet<double>{ 0, 0, 1.0 }, Triplet<double>{ 0, 1, 1.0 },
            Triplet<double>{ 1, 2, 1.0 }, Triplet<double>{ 1, 3, 1.0 },
            Triplet<double>{ 1, 4, 1.0 } };
}

// ProblemBuilder always installs the matrix with the default fill-in headroom,
// so the matrix is set a second time with an explicit spare ratio. With a spare
// ratio of one, the spare space of every row and column is minInterRowSpace.
static Problem<double>
setupProblem( double spareRatio, int minInterRowSpace )
{
   ProblemBuilder<double> problemBuilder;
   problemBuilder.setNumRows( 2 );
   problemBuilder.setNumCols( 5 );
   problemBuilder.setColLbAll( { 0.0, 0.0, 0.0, 0.0, 0.0 } );
   problemBuilder.setColUbAll( { 1.0, 1.0, 1.0, 1.0, 1.0 } );
   problemBuilder.setObjAll( { 0.0, 0.0, 0.0, 0.0, 0.0 } );
   problemBuilder.setRowLhsAll( { 0.0, 0.0 } );
   problemBuilder.setRowRhsAll( { 2.0, 3.0 } );
   problemBuilder.addEntryAll( setupTriplets() );

   Problem<double> problem = problemBuilder.build();

   SparseStorage<double> storage{ setupTriplets(), 2, 5, true, spareRatio,
                                  minInterRowSpace };
   problem.setConstraintMatrix( std::move( storage ), { 0.0, 0.0 },
                                { 2.0, 3.0 }, Vec<RowFlags>( 2 ) );
   problem.recomputeAllActivities();

   return problem;
}

static int
spareSpace( const SparseStorage<double>& storage, int index )
{
   const IndexRange* ranges = storage.getRowRanges();
   return ranges[index + 1].start - ranges[index].end;
}

static int
rowLength( const SparseStorage<double>& storage, int index )
{
   const IndexRange* ranges = storage.getRowRanges();
   return ranges[index].end - ranges[index].start;
}

// Adds every given column to a row in one call. This is how aggregate(),
// sparsify() and clique merging consume a row's fill-in space, and it is the
// only way a row can end up with no spare space at all: change_coefficient()
// inserts one entry at a time and never leaves a row completely full.
static void
growRow( SparseStorage<double>& storage, int row, const Vec<int>& columns )
{
   Vec<double> valbuffer;
   Vec<int> indbuffer;

   storage.changeRow(
       row, 0, static_cast<int>( columns.size() ),
       [&]( int k ) { return columns[k]; }, []( int ) { return 1.0; },
       []( const double&, const double& newval ) { return newval; },
       []( int, int, double, double ) {}, valbuffer, indbuffer );
}

static bool
addCoefficient( Problem<double>& problem, int row, int col )
{
   Num<double> num{};
   Vec<int> indbuffer;
   Vec<double> valbuffer;
   Vec<int> changedActivities;

   return problem.getConstraintMatrix().change_coefficient(
       num, row, col, 1.0, problem.getVariableDomains(), indbuffer, valbuffer,
       changedActivities, problem.getRowActivities(), 0 );
}

TEST_CASE( "change-coefficient-rejects-row-without-spare-space", "[core]" )
{
   // Two spare slots in every row and every column.
   Problem<double> problem = setupProblem( 1.0, 2 );
   ConstraintMatrix<double>& matrix = problem.getConstraintMatrix();
   SparseStorage<double>& storage = matrix.getConstraintMatrix();

   REQUIRE( spareSpace( storage, 0 ) == 2 );

   // Fill row 0 up to exactly its allocation.
   growRow( storage, 0, { 2, 3 } );
   REQUIRE( spareSpace( storage, 0 ) == 0 );
   REQUIRE( rowLength( storage, 0 ) == 4 );

   // Column 4 is not part of row 0, so this insertion has to grow the row, but
   // row 0 has no slot left to grow into. Column 4 itself still has spare
   // space, so the transposed storage cannot be what rejects the insertion.
   REQUIRE( spareSpace( matrix.getMatrixTranspose(), 4 ) == 2 );

   REQUIRE( addCoefficient( problem, 0, 4 ) == false );
   REQUIRE( rowLength( storage, 0 ) == 4 );
}

TEST_CASE( "change-coefficient-uses-the-last-spare-slot", "[core]" )
{
   // Exactly one spare slot in every row and every column.
   Problem<double> problem = setupProblem( 1.0, 1 );
   ConstraintMatrix<double>& matrix = problem.getConstraintMatrix();

   REQUIRE( spareSpace( matrix.getConstraintMatrix(), 0 ) == 1 );
   REQUIRE( spareSpace( matrix.getMatrixTranspose(), 4 ) == 1 );

   // One free slot is exactly enough room for one coefficient, in both the row
   // and the column, so this insertion has to be accepted.
   REQUIRE( addCoefficient( problem, 0, 4 ) == true );
   REQUIRE( rowLength( matrix.getConstraintMatrix(), 0 ) == 3 );
   REQUIRE( matrix.getRowSizes()[0] == 3 );
   REQUIRE( matrix.getColSizes()[4] == 2 );
}
