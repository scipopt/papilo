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
#include "papilo/core/Common.hpp"
#include "catch/catch.hpp"

TEST_CASE( "test activity computation and constraint propagation", "[core]" )
{
   Vec<int> colinds{0, 1, 2, 3, 4, 5, 6, 7};
   Vec<double> rowvalues{1.0, -1.0, 2.0, -2.0, 3.0, -3.0, 4.0, -4.0};

   Vec<double> lbs{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
   Vec<double> ubs{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

   int colrows = 0;

   Vec<RowActivity<double>> activities;
   activities.emplace_back( compute_row_activity(
       rowvalues.data(), colinds.data(), rowvalues.size(), lbs, ubs ) );

   REQUIRE( activities[0].ninfmin == 0 );
   REQUIRE( activities[0].ninfmax == 0 );
   REQUIRE( activities[0].min == -activities[0].max );
   REQUIRE( activities[0].max == ( 1.0 + 2.0 + 3.0 + 4.0 ) );

   int ncalls = 0;
   // fix column 0 to upper bound
   auto activity_callback = [&]( ActivityChange actChange, int row,
                                 RowActivity<double>& activity ) {
      REQUIRE( row == 0 );
      REQUIRE( &activity == &activities[0] );

      if( ncalls == 0 )
      {
         REQUIRE( actChange == ActivityChange::kMin );
      }
      else if( ncalls == 1 )
      {
         REQUIRE( actChange == ActivityChange::kMax );
      }

      ++ncalls;
   };
   update_activities_after_boundchange( &rowvalues[0], &colrows, 1, lbs[0],
                                        ubs[0], ubs[0], ubs[0], activities,
                                        activity_callback );
   lbs[0] = ubs[0];

   REQUIRE( ncalls == 1 );
   REQUIRE( activities[0].ninfmin == 0 );
   REQUIRE( activities[0].ninfmax == 0 );
   REQUIRE( activities[0].max == ( 1.0 + 2.0 + 3.0 + 4.0 ) );
   REQUIRE( activities[0].min == -activities[0].max + 1.0 );

   // fix column 1 to upper bound
   update_activities_after_boundchange( &rowvalues[1], &colrows, 1, lbs[1],
                                        ubs[1], ubs[1], ubs[1], activities,
                                        activity_callback );
   lbs[1] = ubs[1];

   REQUIRE( activities[0].ninfmin == 0 );
   REQUIRE( activities[0].ninfmax == 0 );
   REQUIRE( activities[0].max == ( 1.0 - 1.0 + 2.0 + 3.0 + 4.0 ) );
   REQUIRE( activities[0].min == -activities[0].max );
   REQUIRE( ncalls == 2 );

   update_activities_after_coeffchange(
       lbs[6], ubs[6], rowvalues[6], 3.0, activities[0],
       [&]( ActivityChange actChange, RowActivity<double>& activity ) {
          REQUIRE( &activity == &activities[0] );
          REQUIRE( actChange == ActivityChange::kMax );
          ++ncalls;
       } );

   rowvalues[6] = 3.0;

   REQUIRE( activities[0].max == ( 1.0 - 1.0 + 2.0 + 3.0 + 4.0 - 1.0 ) );
   REQUIRE( activities[0].min == -( 1.0 - 1.0 + 2.0 + 3.0 + 4.0 ) );
   REQUIRE( activities[0].ninfmin == 0 );
   REQUIRE( activities[0].ninfmax == 0 );
   REQUIRE( ncalls == 3 );

   bool maxchanged = false;
   bool minchanged = false;
   update_activities_after_coeffchange(
       lbs[7], ubs[7], rowvalues[7], 3.0, activities[0],
       [&]( ActivityChange actChange, RowActivity<double>& activity ) {
          REQUIRE( &activity == &activities[0] );
          if( actChange == ActivityChange::kMax )
          {
             REQUIRE( maxchanged == false );
             maxchanged = true;
          }
          else if( actChange == ActivityChange::kMin )
          {
             REQUIRE( minchanged == false );
             minchanged = true;
          }
          ++ncalls;
       } );

   rowvalues[7] = 3.0;

   REQUIRE( ncalls == 5 );
   REQUIRE( minchanged );
   REQUIRE( maxchanged );
   REQUIRE( activities[0].ninfmin == 0 );
   REQUIRE( activities[0].ninfmax == 0 );
   REQUIRE( activities[0].max == ( 1.0 - 1.0 + 2.0 + 3.0 + 4.0 - 1.0 + 3.0 ) );
   REQUIRE( activities[0].min == -( 1.0 - 1.0 + 2.0 + 3.0 ) );

   RowActivity<double> recompute = compute_row_activity(
       rowvalues.data(), colinds.data(), rowvalues.size(), lbs, ubs );

   REQUIRE( recompute.min == activities[0].min );
   REQUIRE( recompute.max == activities[0].max );
   REQUIRE( recompute.ninfmin == activities[0].ninfmin );
   REQUIRE( recompute.ninfmax == activities[0].ninfmax );

   Vec<double> lbs_cpy( lbs );
   Vec<double> ubs_cpy( ubs );

   ncalls = 0;
   propagate_row( rowvalues.data(), colinds.data(), rowvalues.size(),
                  activities[0], -infinity<double>(), activities[0].min, lbs,
                  ubs, [&]( BoundChange bndChg, int colid, double newbnd ) {
                     if( bndChg == BoundChange::kLower )
                     {
                        lbs_cpy[colid] = newbnd;
                     }
                     else if( bndChg == BoundChange::kUpper )
                     {
                        ubs_cpy[colid] = newbnd;
                     }
                     ++ncalls;
                  } );
   REQUIRE( ncalls == 6 );

   // everything should be fixed to its bounds
   for( size_t i = 0; i < lbs_cpy.size(); ++i )
   {
      REQUIRE( lbs_cpy[i] == ubs_cpy[i] );
   }

   // since everything is fixed maximum and minimum activity should match and be
   // equal to activities[0].min since this has been passed as side
   recompute = compute_row_activity( rowvalues.data(), colinds.data(),
                                     rowvalues.size(), lbs_cpy, ubs_cpy );

   REQUIRE( recompute.min == activities[0].min );
   REQUIRE( recompute.min == recompute.max );

   lbs_cpy = lbs;
   ubs_cpy = ubs;

   ncalls = 0;
   propagate_row( rowvalues.data(), colinds.data(), rowvalues.size(),
                  activities[0], activities[0].max, infinity<double>(), lbs,
                  ubs, [&]( BoundChange bndChg, int colid, double newbnd ) {
                     if( bndChg == BoundChange::kLower )
                     {
                        lbs_cpy[colid] = newbnd;
                     }
                     else if( bndChg == BoundChange::kUpper )
                     {
                        ubs_cpy[colid] = newbnd;
                     }
                     ++ncalls;
                  } );
   REQUIRE( ncalls == 6 );

   // everything should be fixed to its bounds
   for( size_t i = 0; i < lbs_cpy.size(); ++i )
   {
      REQUIRE( lbs_cpy[i] == ubs_cpy[i] );
   }

   // since everything is fixed maximum and minimum activity should match and be
   // equal to activities[0].min since this has been passed as side
   recompute = compute_row_activity( rowvalues.data(), colinds.data(),
                                     rowvalues.size(), lbs_cpy, ubs_cpy );

   REQUIRE( recompute.min == activities[0].max );
   REQUIRE( recompute.min == recompute.max );
}