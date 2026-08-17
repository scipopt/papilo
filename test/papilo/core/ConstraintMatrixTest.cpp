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

namespace papilo
{
static Vec<Triplet<double>>
setupTriplets()
{
   // 1  1  0  0  0
   // 0  0  1  1  1
   return { Triplet<double>{ 0, 0, 1.0 }, Triplet<double>{ 0, 1, 1.0 },
            Triplet<double>{ 1, 2, 1.0 }, Triplet<double>{ 1, 3, 1.0 },
            Triplet<double>{ 1, 4, 1.0 } };
}

// ProblemBuilder is used for the bounds and sides while the constraint matrix is
// set separately with an explicit spare ratio. With a spare ratio of one, the
// spare space of every row and column is minInterRowSpace. For transposed
// enabled the same matrix is installed transposed, so a column of the original
// matrix becomes a row.
static Problem<double>
setupProblem( double spareRatio, int minInterRowSpace, bool transposed = false )
{
   ProblemBuilder<double> problemBuilder;
   problemBuilder.setNumRows( transposed ? 5 : 2 );
   problemBuilder.setNumCols( transposed ? 2 : 5 );
   problemBuilder.setColLbAll( Vec<double>( transposed ? 2 : 5, 0.0 ) );
   problemBuilder.setColUbAll( Vec<double>( transposed ? 2 : 5, 1.0 ) );
   problemBuilder.setObjAll( Vec<double>( transposed ? 2 : 5, 0.0 ) );
   problemBuilder.setRowLhsAll( Vec<double>( transposed ? 5 : 2, 0.0 ) );
   problemBuilder.setRowRhsAll( transposed ? Vec<double>( 5, 1.0 )
         : Vec<double>{ 2.0, 3.0 } );

   Problem<double> problem = problemBuilder.build();
   SparseStorage<double> storage{ setupTriplets(), 2, 5, true, spareRatio,
         minInterRowSpace };
   problem.setConstraintMatrix( std::move( storage ),
         Vec<double>( transposed ? 5 : 2, 0.0 ),
         transposed ? Vec<double>( 5, 1.0 ) : Vec<double>{ 2.0, 3.0 },
         Vec<RowFlags>( transposed ? 5 : 2 ), transposed );
   problem.recomputeAllActivities();

   return problem;
}

static int
getRowSpare( const SparseStorage<double>& storage, int index )
{
   const IndexRange* ranges = storage.getRowRanges();
   return ranges[index + 1].start - ranges[index].end;
}

static int
getRowLength( const SparseStorage<double>& storage, int index )
{
   const IndexRange* ranges = storage.getRowRanges();
   return ranges[index].end - ranges[index].start;
}

// Adds every given column to a row in one call, which is how aggregate(),
// sparsify(), and clique merging consume row spare space.
static void
growRow( SparseStorage<double>& storage, int row, const Vec<int>& columns )
{
   Vec<double> valbuffer;
   Vec<int> indbuffer;

   storage.changeRow( row, 0, static_cast<int>( columns.size() ),
         [&]( int k ) { return columns[k]; }, []( int ) { return 1.0; },
         []( const double&, const double& newval ) { return newval; },
         []( int, int, double, double ) {}, valbuffer, indbuffer );
}

// Adds every given row to a column in one call, which is how aggregate(),
// sparsify(), and clique merging consume column spare space.
static void
growCol( SparseStorage<double>& storage, int col, const Vec<int>& rows )
{
   Vec<double> valbuffer;
   Vec<int> indbuffer;

   storage.changeRow( col, 0, static_cast<int>( rows.size() ),
         [&]( int k ) { return rows[k]; }, []( int ) { return 1.0; },
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

   return problem.getConstraintMatrix().change_coefficient( num, row, col, 1.0,
         problem.getVariableDomains(), indbuffer, valbuffer, changedActivities,
         problem.getRowActivities(), 0 );
}

TEST_CASE( "change-coefficient-rejects-row-without-spare-space", "[core]" )
{
   // Two spare slots in every row and every column.
   Problem<double> problem = setupProblem( 1.0, 2 );
   ConstraintMatrix<double>& matrix = problem.getConstraintMatrix();
   SparseStorage<double>& storage = matrix.getConstraintMatrix();

   // Fill row 0 up to exactly its allocation.
   REQUIRE( getRowSpare( storage, 0 ) == 2 );
   growRow( storage, 0, { 2, 3 } );
   REQUIRE( getRowSpare( storage, 0 ) == 0 );
   REQUIRE( getRowLength( storage, 0 ) == 4 );

   // Column 4 is not part of row 0, so this insertion has to grow the row, but
   // row 0 has no slot left to grow into. Column 4 itself still has spare space,
   // so the transposed sparse storage cannot be what rejects the insertion.
   REQUIRE( getRowSpare( matrix.getMatrixTranspose(), 4 ) == 2 );
   const bool inserted = addCoefficient( problem, 0, 4 );
   REQUIRE( !inserted );
   REQUIRE( getRowLength( storage, 0 ) == 4 );
}

TEST_CASE( "change-coefficient-rejects-column-without-spare-space", "[core]" )
{
   // Two spare slots in every column and every row.
   Problem<double> problem = setupProblem( 1.0, 2, true );
   ConstraintMatrix<double>& matrix = problem.getConstraintMatrix();
   SparseStorage<double>& storage = matrix.getMatrixTranspose();

   // Fill column 0 up to exactly its allocation.
   REQUIRE( getRowSpare( storage, 0 ) == 2 );
   growCol( storage, 0, { 2, 3 } );
   REQUIRE( getRowSpare( storage, 0 ) == 0 );
   REQUIRE( getRowLength( storage, 0 ) == 4 );

   // Row 4 is not part of column 0, so this insertion has to grow the column, but
   // column 0 has no slot left to grow into. Row 4 itself still has spare space,
   // so the sparse storage cannot be what rejects the insertion.
   REQUIRE( getRowSpare( matrix.getConstraintMatrix(), 4 ) == 2 );
   const bool inserted = addCoefficient( problem, 4, 0 );
   REQUIRE( !inserted );
   REQUIRE( getRowLength( storage, 0 ) == 4 );
}

TEST_CASE( "change-coefficient-uses-the-last-spare-slot", "[core]" )
{
   // Exactly one spare slot in every row and every column.
   Problem<double> problem = setupProblem( 1.0, 1 );
   ConstraintMatrix<double>& matrix = problem.getConstraintMatrix();

   // One free slot is exactly enough room for one coefficient, in both the row
   // and the column, so this insertion has to be accepted.
   REQUIRE( getRowSpare( matrix.getConstraintMatrix(), 0 ) == 1 );
   REQUIRE( getRowSpare( matrix.getMatrixTranspose(), 4 ) == 1 );
   const bool inserted = addCoefficient( problem, 0, 4 );
   REQUIRE( inserted );
   REQUIRE( getRowLength( matrix.getConstraintMatrix(), 0 ) == 3 );
   REQUIRE( matrix.getRowSizes()[0] == 3 );
   REQUIRE( matrix.getColSizes()[4] == 2 );
}
} // namespace papilo
