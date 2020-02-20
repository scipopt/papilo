/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*               This file is part of the program and library                */
/*    PaPILO --- Parallel Presolve for Integer and Linear Optimization       */
/*                                                                           */
/* Copyright (C) 2020  Konrad-Zuse-Zentrum                                   */
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

#define CATCH_CONFIG_MAIN
#include "papilo/core/SparseStorage.hpp"
#include "catch/catch.hpp"
#include "papilo/misc/compress_vector.hpp"

TEST_CASE( "sparse storage can be created from triplets and compressed",
           "[core]" )
{
   // build the triplets for the following matrix:
   // 1  2  0  0  0  0  0  0  0
   // 0  3  4  5  6  7  0  0  0
   // 0  8  0  0  0  0  0  0  0
   // 0  0  0  0  0  0  0  0  0
   // 9  10 11 0  0  0  0  12 13
   int nrows = 5;
   int ncols = 9;
   Vec<Triplet<double>> triplets = {
       Triplet<double>{0, 0, 1.0},  Triplet<double>{0, 1, 2.0},
       Triplet<double>{1, 1, 3.0},  Triplet<double>{1, 2, 4.0},
       Triplet<double>{1, 3, 5.0},  Triplet<double>{1, 4, 6.0},
       Triplet<double>{1, 5, 7.0},  Triplet<double>{2, 1, 8.0},
       Triplet<double>{4, 0, 9.0},  Triplet<double>{4, 1, 10.0},
       Triplet<double>{4, 2, 11.0}, Triplet<double>{4, 7, 12.0},
       Triplet<double>{4, 8, 13.0}};

   Vec<int> rowsize = {2, 5, 1, -1, 5};
   Vec<int> colsize = {2, 4, 2, 1, 1, 1, -1, 1, 1};
   SparseStorage<double> matrix{triplets, nrows, ncols, true};
   SparseStorage<double> transpose = matrix.getTranspose();

   REQUIRE( matrix.getNnz() == triplets.size() );
   REQUIRE( matrix.getNRows() == nrows );
   REQUIRE( matrix.getNCols() == ncols );

   REQUIRE( transpose.getNnz() == triplets.size() );
   REQUIRE( transpose.getNRows() == ncols );
   REQUIRE( transpose.getNCols() == nrows );

   auto rowranges = matrix.getRowRanges();
   auto rowvalues = matrix.getValues();
   auto columns = matrix.getColumns();
   for( int i = 0; i < nrows; ++i )
   {
      REQUIRE( rowranges[i].end - rowranges[i].start == rowsize[i] );
      double* row = rowvalues + rowranges[i].start;
      int* rowcols = columns + rowranges[i].start;
      switch( i )
      {
      case 0:
         REQUIRE( row[0] == 1.0 );
         REQUIRE( row[1] == 2.0 );
         REQUIRE( rowcols[0] == 0 );
         REQUIRE( rowcols[1] == 1 );
         break;
      case 1:
         REQUIRE( row[0] == 3.0 );
         REQUIRE( row[1] == 4.0 );
         REQUIRE( row[2] == 5.0 );
         REQUIRE( row[3] == 6.0 );
         REQUIRE( row[4] == 7.0 );
         REQUIRE( rowcols[0] == 1 );
         REQUIRE( rowcols[1] == 2 );
         REQUIRE( rowcols[2] == 3 );
         REQUIRE( rowcols[3] == 4 );
         REQUIRE( rowcols[4] == 5 );
         break;
      case 2:
         REQUIRE( row[0] == 8.0 );
         REQUIRE( rowcols[0] == 1 );
         break;
      case 3:
         continue;
      case 4:
         REQUIRE( row[0] == 9.0 );
         REQUIRE( row[1] == 10.0 );
         REQUIRE( row[2] == 11.0 );
         REQUIRE( row[3] == 12.0 );
         REQUIRE( row[4] == 13.0 );
         REQUIRE( rowcols[0] == 0 );
         REQUIRE( rowcols[1] == 1 );
         REQUIRE( rowcols[2] == 2 );
         REQUIRE( rowcols[3] == 7 );
         REQUIRE( rowcols[4] == 8 );
      }

      // todo check the rowvalues and columns are correct
   }

   auto colranges = transpose.getRowRanges();
   auto colvalues = transpose.getValues();
   auto rows = transpose.getColumns();
   for( int i = 0; i < ncols; ++i )
   {
      REQUIRE( colranges[i].end - colranges[i].start == colsize[i] );
      double* col = colvalues + colranges[i].start;
      int* colrows = rows + colranges[i].start;
      switch( i )
      {
      case 0:
         REQUIRE( col[0] == 1.0 );
         REQUIRE( col[1] == 9.0 );
         REQUIRE( colrows[0] == 0 );
         REQUIRE( colrows[1] == 4 );
         break;
      case 1:
         REQUIRE( col[0] == 2.0 );
         REQUIRE( col[1] == 3.0 );
         REQUIRE( col[2] == 8.0 );
         REQUIRE( col[3] == 10.0 );
         REQUIRE( colrows[0] == 0 );
         REQUIRE( colrows[1] == 1 );
         REQUIRE( colrows[2] == 2 );
         REQUIRE( colrows[3] == 4 );
         break;
      case 2:
         REQUIRE( col[0] == 4.0 );
         REQUIRE( col[1] == 11.0 );
         REQUIRE( colrows[0] == 1 );
         REQUIRE( colrows[1] == 4 );
         break;
      case 3:
         REQUIRE( col[0] == 5.0 );
         REQUIRE( colrows[0] == 1 );
         break;
      case 4:
         REQUIRE( col[0] == 6.0 );
         REQUIRE( colrows[0] == 1 );
         break;
      case 5:
         REQUIRE( col[0] == 7.0 );
         REQUIRE( colrows[0] == 1 );
         break;
      case 6:
         continue;
      case 7:
         REQUIRE( col[0] == 12.0 );
         REQUIRE( colrows[0] == 4 );
         break;
      case 8:
         REQUIRE( col[0] == 13.0 );
         REQUIRE( colrows[0] == 4 );
         break;
      }
   }

   Vec<int> col_mapping = matrix.compress( colsize );
   Vec<int> row_mapping = transpose.compress( rowsize );

   REQUIRE( col_mapping.size() == ncols );
   REQUIRE( row_mapping.size() == nrows );

   /* check if empty row at index 3 was removed properly */
   for( int i = 0; i < nrows; ++i )
   {
      if( i < 3 )
         REQUIRE( row_mapping[i] == i );
      else if( i == 3 )
         REQUIRE( row_mapping[i] == -1 );
      else
         REQUIRE( row_mapping[i] == i - 1 );
   }

   /* check if empty col at index 6 was removed properly */
   for( int i = 0; i < ncols; ++i )
   {
      if( i < 6 )
         REQUIRE( col_mapping[i] == i );
      else if( i == 6 )
         REQUIRE( col_mapping[i] == -1 );
      else
         REQUIRE( col_mapping[i] == i - 1 );
   }

   --nrows;
   --ncols;

   compress_vector( col_mapping, colsize );
   compress_vector( row_mapping, rowsize );

   rowranges = matrix.getRowRanges();
   rowvalues = matrix.getValues();
   columns = matrix.getColumns();
   for( int i = 0; i < nrows; ++i )
   {
      REQUIRE( rowranges[i].end - rowranges[i].start == rowsize[i] );
      double* row = rowvalues + rowranges[i].start;
      int* rowcols = columns + rowranges[i].start;
      switch( i )
      {
      case 0:
         REQUIRE( row[0] == 1.0 );
         REQUIRE( row[1] == 2.0 );
         REQUIRE( rowcols[0] == 0 );
         REQUIRE( rowcols[1] == 1 );
         break;
      case 1:
         REQUIRE( row[0] == 3.0 );
         REQUIRE( row[1] == 4.0 );
         REQUIRE( row[2] == 5.0 );
         REQUIRE( row[3] == 6.0 );
         REQUIRE( row[4] == 7.0 );
         REQUIRE( rowcols[0] == 1 );
         REQUIRE( rowcols[1] == 2 );
         REQUIRE( rowcols[2] == 3 );
         REQUIRE( rowcols[3] == 4 );
         REQUIRE( rowcols[4] == 5 );
         break;
      case 2:
         REQUIRE( row[0] == 8.0 );
         REQUIRE( rowcols[0] == 1 );
         break;
      case 3:
         REQUIRE( row[0] == 9.0 );
         REQUIRE( row[1] == 10.0 );
         REQUIRE( row[2] == 11.0 );
         REQUIRE( row[3] == 12.0 );
         REQUIRE( row[4] == 13.0 );
         REQUIRE( rowcols[0] == 0 );
         REQUIRE( rowcols[1] == 1 );
         REQUIRE( rowcols[2] == 2 );
         REQUIRE( rowcols[3] == 6 );
         REQUIRE( rowcols[4] == 7 );
      }

      // todo check the rowvalues and columns are correct
   }

   colranges = transpose.getRowRanges();
   colvalues = transpose.getValues();
   rows = transpose.getColumns();
   for( int i = 0; i < ncols; ++i )
   {
      REQUIRE( colranges[i].end - colranges[i].start == colsize[i] );
      double* col = colvalues + colranges[i].start;
      int* colrows = rows + colranges[i].start;
      switch( i )
      {
      case 0:
         REQUIRE( col[0] == 1.0 );
         REQUIRE( col[1] == 9.0 );
         REQUIRE( colrows[0] == 0 );
         REQUIRE( colrows[1] == 3 );
         break;
      case 1:
         REQUIRE( col[0] == 2.0 );
         REQUIRE( col[1] == 3.0 );
         REQUIRE( col[2] == 8.0 );
         REQUIRE( col[3] == 10.0 );
         REQUIRE( colrows[0] == 0 );
         REQUIRE( colrows[1] == 1 );
         REQUIRE( colrows[2] == 2 );
         REQUIRE( colrows[3] == 3 );
         break;
      case 2:
         REQUIRE( col[0] == 4.0 );
         REQUIRE( col[1] == 11.0 );
         REQUIRE( colrows[0] == 1 );
         REQUIRE( colrows[1] == 3 );
         break;
      case 3:
         REQUIRE( col[0] == 5.0 );
         REQUIRE( colrows[0] == 1 );
         break;
      case 4:
         REQUIRE( col[0] == 6.0 );
         REQUIRE( colrows[0] == 1 );
         break;
      case 5:
         REQUIRE( col[0] == 7.0 );
         REQUIRE( colrows[0] == 1 );
         break;
      case 6:
         REQUIRE( col[0] == 12.0 );
         REQUIRE( colrows[0] == 3 );
         break;
      case 7:
         REQUIRE( col[0] == 13.0 );
         REQUIRE( colrows[0] == 3 );
         break;
      }
   }
}