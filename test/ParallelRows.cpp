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
#include "catch/catch.hpp"
#include "papilo/core/ParallelRow.hpp"
#include "papilo/misc/compress_vector.hpp"
#include "papilo/misc/fmt.hpp"

// TODO fix

static bool
operator==( const Reductions<double>::Reduction& lhs,
            const Reductions<double>::Reduction& rhs )
{
   double max = std::max( double{1.0}, std::max( std::fabs( lhs.newval ),
                                                 std::fabs( rhs.newval ) ) );
   bool equalityCheck = std::fabs( lhs.newval - rhs.newval ) <=
                        ( std::numeric_limits<double>::epsilon() * max );

   return equalityCheck && ( lhs.row == rhs.row ) && ( lhs.col == rhs.col );
}

TEST_CASE( "Parallel row detection", "[core]" )
{
   double coef1[] = {1, 1, 1, 1, 1};
   double coef2[] = {2, 2, 2, 2, 2};
   int id1[] = {1, 2, 3, 4, 5};
   int id2[] = {1, 2, 3, 4, 5};

   double coef3[] = {1.5, -1.5, 1.5, 1.5, -1.5};
   double coef4[] = {-1, 1, -1, -1, 1};
   int id3[] = {1, 2, 3, 4, 5};
   int id4[] = {1, 2, 3, 4, 5};

   double coef5[] = {1, 1, 1, 1, 1};
   int id5[] = {1, 2, 4, 6, 9};

   double coef6[] = {1, 1, 1, -1, -1};
   double coef7[] = {2, 2, 2, -2, -2};
   int id6[] = {1, 2, 4, 5, 6};
   int id7[] = {1, 2, 4, 5, 6};

   double coef8[] = {1, 1, -1, -1, -1};
   double coef9[] = {1.0000001, 1, -1, -1, -1};
   int id8[] = {1, 2, 4, 5, 6};
   int id9[] = {1, 2, 4, 5, 6};

   const double* coefs[] = {coef1, coef2, coef3, coef4, coef5,
                            coef6, coef7, coef8, coef9};
   const int* ids[] = {id1, id2, id3, id4, id5, id6, id7, id8, id9};
   const int lens[] = {5, 5, 5, 5, 5, 5, 5, 5, 5};

   double inf = infinity<double>();
   Vec<double> lhs_values{-1, -2, 0, -1, -1, -inf, -inf, -1, -3};
   Vec<double> rhs_values{4, 1, 0, 1, 1, 2, 1, 1, 3};

   PresolveResult result;
   Reductions<double> reductions;

   // test lambda used in the presolver
   auto handlerows = [&reductions, &result, &lhs_values,
                      &rhs_values]( int row1, int row2, double ratio ) {
      bool firstconsEquality = ( rhs_values[row1] == lhs_values[row1] );
      bool secondconsEquality = ( rhs_values[row2] == lhs_values[row2] );
      assert( ratio != 0.0 );
      // TODO use this
      double adjustedLHS = lhs_values[row2] * ratio;
      double adjustedRHS = rhs_values[row2] * ratio;

      if( firstconsEquality && secondconsEquality )
      {
         // TODO fp comparaison
         if( rhs_values[row1] == adjustedRHS )
         {
            TransactionGuard<double> guard{reductions};
            reductions.lockRow( row1 );
            reductions.lockRow( row2 );
            reductions.markRowRedundant( row2 );
         }
         else
            result = PresolveResult::INFEASIBLE;
      }
      else if( firstconsEquality && !secondconsEquality )
      {
         if( ratio > 0 )
         {
            if( rhs_values[row1] <= adjustedRHS &&
                rhs_values[row1] >= adjustedLHS )
            {
               TransactionGuard<double> guard{reductions};
               reductions.lockRow( row1 );
               reductions.lockRow( row2 );
               reductions.markRowRedundant( row2 );
            }
            else
               result = PresolveResult::INFEASIBLE;
         }
         else
         {
            if( rhs_values[row1] <= adjustedLHS &&
                rhs_values[row1] >= adjustedRHS )
            {
               TransactionGuard<double> guard{reductions};
               reductions.lockRow( row1 );
               reductions.lockRow( row2 );
               reductions.markRowRedundant( row2 );
            }
            else
               result = PresolveResult::INFEASIBLE;
         }
      }
      else if( !firstconsEquality && secondconsEquality )
      {
         if( ratio > 0 )
         {
            if( rhs_values[row1] >= adjustedRHS &&
                lhs_values[row1] <= adjustedRHS )
            {
               TransactionGuard<double> guard{reductions};
               reductions.lockRow( row1 );
               reductions.lockRow( row2 );
               reductions.markRowRedundant( row1 );
            }
            else
               result = PresolveResult::INFEASIBLE;
         }
         else
         {
            if( rhs_values[row1] >= adjustedRHS &&
                rhs_values[row1] <= adjustedRHS )
            {
               TransactionGuard<double> guard{reductions};
               reductions.lockRow( row1 );
               reductions.lockRow( row2 );
               reductions.markRowRedundant( row1 );
            }
            else
               result = PresolveResult::INFEASIBLE;
         }
      }
      else
      {
         if( ratio > 0 )
         {
            if( rhs_values[row1] < adjustedLHS ||
                adjustedRHS < lhs_values[row1] )
               result = PresolveResult::INFEASIBLE;
            else
            {
               TransactionGuard<double> guard{reductions};
               reductions.lockRow( row1 );
               reductions.lockRow( row2 );
               if( rhs_values[row1] < adjustedRHS )
               {
                  reductions.markRowRedundant( row2 );
                  if( lhs_values[row1] < adjustedLHS )
                     reductions.changeRowLHS( row1, adjustedLHS );
               }
               else
               {
                  reductions.markRowRedundant( row1 );
                  if( lhs_values[row1] > adjustedLHS )
                     reductions.changeRowLHS( row2, lhs_values[row1] );
               }
            }
         }
         else
         {
            if( rhs_values[row1] < adjustedLHS ||
                adjustedRHS < lhs_values[row1] )
               result = PresolveResult::INFEASIBLE;
            else
            {
               TransactionGuard<double> guard{reductions};
               reductions.lockRow( row1 );
               reductions.lockRow( row2 );
               if( rhs_values[row1] < adjustedLHS )
               {
                  reductions.markRowRedundant( row2 );
                  if( lhs_values[row1] < adjustedRHS )
                     reductions.changeRowLHS( row1, adjustedRHS );
               }
               else
               {
                  reductions.markRowRedundant( row1 );
                  if( lhs_values[row1] > adjustedRHS )
                     reductions.changeRowLHS( row2, lhs_values[row1] );
               }
            }
         }
      }
   };

   getParallelRows( coefs, ids, 9, lens, handlerows );

   REQUIRE( reductions.size() == 12 );
   const Reductions<double>::Reduction& reduction =
       reductions.getReduction( 0 );
   REQUIRE( reduction == Reductions<double>::Reduction( 0.0, 0, -4 ) );

   const Reductions<double>::Reduction& reduction1 =
       reductions.getReduction( 1 );
   REQUIRE( reduction1 == Reductions<double>::Reduction( 0.0, 1, -4 ) );

   const Reductions<double>::Reduction& reduction2 =
       reductions.getReduction( 2 );
   REQUIRE( reduction2 == Reductions<double>::Reduction( 0.0, 0, -3 ) );

   const Reductions<double>::Reduction& reduction3 =
       reductions.getReduction( 3 );
   REQUIRE( reduction3 == Reductions<double>::Reduction( 0.0, 2, -4 ) );

   const Reductions<double>::Reduction& reduction4 =
       reductions.getReduction( 4 );
   REQUIRE( reduction4 == Reductions<double>::Reduction( 0.0, 3, -4 ) );

   const Reductions<double>::Reduction& reduction5 =
       reductions.getReduction( 5 );
   REQUIRE( reduction5 == Reductions<double>::Reduction( 0.0, 3, -3 ) );

   const Reductions<double>::Reduction& reduction6 =
       reductions.getReduction( 6 );
   REQUIRE( reduction6 == Reductions<double>::Reduction( 0.0, 5, -4 ) );

   const Reductions<double>::Reduction& reduction7 =
       reductions.getReduction( 7 );
   REQUIRE( reduction7 == Reductions<double>::Reduction( 0.0, 6, -4 ) );

   const Reductions<double>::Reduction& reduction8 =
       reductions.getReduction( 8 );
   REQUIRE( reduction8 == Reductions<double>::Reduction( 0.0, 5, -3 ) );

   const Reductions<double>::Reduction& reduction9 =
       reductions.getReduction( 9 );
   REQUIRE( reduction9 == Reductions<double>::Reduction( 0.0, 7, -4 ) );

   const Reductions<double>::Reduction& reduction10 =
       reductions.getReduction( 10 );
   REQUIRE( reduction10 == Reductions<double>::Reduction( 0.0, 8, -4 ) );

   const Reductions<double>::Reduction& reduction11 =
       reductions.getReduction( 11 );
   REQUIRE( reduction11 == Reductions<double>::Reduction( 0.0, 8, -3 ) );
}
