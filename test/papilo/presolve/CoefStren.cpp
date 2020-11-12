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
#include "papilo/core/PresolveMethod.hpp"
#include "papilo/core/Strengthening.hpp"
#include "papilo/misc/compress_vector.hpp"

#include <cmath>
#include <limits>

static bool
operator==( const Reductions<float>::Reduction& lhs,
            const Reductions<float>::Reduction& rhs )
{
   float max = std::max( float{1.0}, std::max( std::fabs( lhs.newval ),
                                               std::fabs( rhs.newval ) ) );
   bool equalityCheck = std::fabs( lhs.newval - rhs.newval ) <=
                        ( std::numeric_limits<float>::epsilon() * max );

   return equalityCheck && ( lhs.row == rhs.row ) && ( lhs.col == rhs.col );
}

TEST_CASE( "test coefficient tightening reduction", "[core]" )
{
   Reductions<float> reductions;

   auto constraintchange = [&reductions](
                               Action constchange, int i,
                               Vec<std::pair<int, float>> rowcoefchanges,
                               float newbound ) {
      assert( !rowcoefchanges.empty() );

      TransactionGuard<float> guard{reductions};
      reductions.lockRow( i );

      for( auto coefchange : rowcoefchanges )
      {
         reductions.changeMatrixEntry( i, coefchange.first, coefchange.second );
      }

      if( constchange == Action::CHANGE_COEF_UPPER_BOUND )
         reductions.changeRowRHS( i, newbound );
      else if( constchange == Action::CHANGE_COEF_LOWER_BOUND )
         reductions.changeRowLHS( i, newbound );
   };

   // test constraint x + 2y <= 5.1
   // 0 <= x <= 2
   // 0 <= y <= 2
   // x,y integer
   float coefficients[] = {1, 2};
   int indices[] = {0, 1};
   RowActivity<float> activity;
   activity.max = 6.0;
   activity.min = 0.0;
   float lhs = -infinity<float>();
   float rhs = 5.1;
   VariableDomains<float> domains;
   domains.lower_bounds = Vec<float>{0.0, 0.0};
   domains.upper_bounds = Vec<float>{2.0, 2.0};
   boost::dynamic_bitset<> integral( 2 );
   integral[0] = 1;
   integral[1] = 1;
   domains.is_integral = integral;

   strengthen_coefficients( coefficients, indices, 2, lhs, rhs, activity,
                            domains, 0, constraintchange );

   REQUIRE( reductions.size() == 4 );
   REQUIRE( reductions.getReduction( 0 ) ==
            Reductions<float>::Reduction( 0.0, 0, -4 ) );
   REQUIRE( reductions.getReduction( 1 ) ==
            Reductions<float>::Reduction( 0.9, 0, 0 ) );
   REQUIRE( reductions.getReduction( 2 ) ==
            Reductions<float>::Reduction( 0.9, 0, 1 ) );
   REQUIRE( reductions.getReduction( 3 ) ==
            Reductions<float>::Reduction( 2.7, 0, -1 ) );

   // test constraint  -x - 2y => -5.1
   // 0 <= x <= 2
   // 0 <= y <= 2
   // x,y integer
   coefficients[0] = -1.0;
   coefficients[1] = -2.0;
   indices[0] = 0;
   indices[1] = 1;
   activity.max = 0.0;
   activity.min = -6.0;
   lhs = -5.1;
   rhs = infinity<float>();
   domains.lower_bounds = Vec<float>{0.0, 0.0};
   domains.upper_bounds = Vec<float>{2.0, 2.0};
   integral[0] = 1;
   integral[1] = 1;
   domains.is_integral = integral;

   reductions.clear();
   strengthen_coefficients( coefficients, indices, 2, lhs, rhs, activity,
                            domains, 0, constraintchange );

   REQUIRE( reductions.size() == 4 );
   REQUIRE( reductions.getReduction( 0 ) ==
            Reductions<float>::Reduction( 0.0, 0, -4 ) );
   REQUIRE( reductions.getReduction( 1 ) ==
            Reductions<float>::Reduction( -0.9, 0, 0 ) );
   REQUIRE( reductions.getReduction( 2 ) ==
            Reductions<float>::Reduction( -0.9, 0, 1 ) );
   REQUIRE( reductions.getReduction( 3 ) ==
            Reductions<float>::Reduction( -2.7, 0, -2 ) );

   // test constraint  -x + 2y <= +5.1
   // -2 <= x <= 0
   // 0 <= y <= 2
   // x,y integer
   coefficients[0] = -1.0;
   coefficients[1] = 2.0;
   indices[0] = 0;
   indices[1] = 1;
   activity.max = 6.0;
   activity.min = 0.0;
   lhs = -infinity<float>();
   rhs = 5.1;
   domains.lower_bounds = Vec<float>{-2.0, 0.0};
   domains.upper_bounds = Vec<float>{0.0, 2.0};
   integral[0] = 1;
   integral[1] = 1;
   domains.is_integral = integral;

   reductions.clear();
   strengthen_coefficients( coefficients, indices, 2, lhs, rhs, activity,
                            domains, 0, constraintchange );

   REQUIRE( reductions.size() == 4 );
   REQUIRE( reductions.getReduction( 0 ) ==
            Reductions<float>::Reduction( 0.0, 0, -4 ) );
   REQUIRE( reductions.getReduction( 1 ) ==
            Reductions<float>::Reduction( -0.9, 0, 0 ) );
   REQUIRE( reductions.getReduction( 2 ) ==
            Reductions<float>::Reduction( 0.9, 0, 1 ) );
   REQUIRE( reductions.getReduction( 3 ) ==
            Reductions<float>::Reduction( 2.7, 0, -1 ) );

   // test constraint  x - 2y >= -5.1
   // -2 <= x <= 0
   // 0 <= y <= 2
   // x,y integer
   coefficients[0] = 1.0;
   coefficients[1] = -2.0;
   indices[0] = 0;
   indices[1] = 1;
   activity.max = 0.0;
   activity.min = -6.0;
   lhs = -5.1;
   rhs = infinity<float>();
   domains.lower_bounds = Vec<float>{-2.0, 0.0};
   domains.upper_bounds = Vec<float>{0.0, 2.0};
   integral[0] = 1;
   integral[1] = 1;
   domains.is_integral = integral;

   reductions.clear();
   strengthen_coefficients( coefficients, indices, 2, lhs, rhs, activity,
                            domains, 0, constraintchange );

   REQUIRE( reductions.size() == 4 );
   REQUIRE( reductions.getReduction( 0 ) ==
            Reductions<float>::Reduction( 0.0, 0, -4 ) );
   REQUIRE( reductions.getReduction( 1 ) ==
            Reductions<float>::Reduction( 0.9, 0, 0 ) );
   REQUIRE( reductions.getReduction( 2 ) ==
            Reductions<float>::Reduction( -0.9, 0, 1 ) );
   REQUIRE( reductions.getReduction( 3 ) ==
            Reductions<float>::Reduction( -2.7, 0, -2 ) );
}
