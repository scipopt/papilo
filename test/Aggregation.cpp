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
#include "papilo/core/Presolve.hpp"
#include "papilo/core/PresolveMethod.hpp"
#include "papilo/core/Problem.hpp"
#include "papilo/misc/compress_vector.hpp"
#include "papilo/presolvers/FreeVarSubstitution.hpp"

#include <tuple>

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

TEST_CASE( "test free variable detection ", "[core]" )
{
   // 2*x1 + x2               >= 1
   //   x1      + x3          <= 2
   //   x1           + x4 - x5 = 1
   // -3 <= x1 <= +3; x2 <= 1; 0 <= x3
   Problem<float> problem;

   Vec<float> lhs_values( 3 );
   Vec<float> rhs_values( 3 );

   lhs_values[0] = 1.0;
   rhs_values[0] = infinity<float>();

   lhs_values[1] = -infinity<float>();
   rhs_values[1] = 2.0;

   lhs_values[2] = 1.0;
   rhs_values[2] = 1.0;

   Vec<Triplet<float>> entries( 7 );
   entries[0] = std::make_tuple( 0, 0, 2.0 );
   entries[1] = std::make_tuple( 0, 1, 1.0 );

   entries[2] = std::make_tuple( 1, 0, 1.0 );
   entries[3] = std::make_tuple( 1, 2, 1.0 );

   entries[4] = std::make_tuple( 2, 0, 1.0 );
   entries[5] = std::make_tuple( 2, 3, 1.0 );
   entries[6] = std::make_tuple( 2, 4, -1.0 );

   problem.setConstraintMatrix( SparseStorage<float>( entries, 3, 5, 1.0 ),
                                lhs_values, rhs_values );

   Vec<float> lower_bounds( 5 );
   Vec<float> upper_bounds( 5 );
   boost::dynamic_bitset<> integral( 5 );

   lower_bounds[0] = -3;
   upper_bounds[0] = 3;

   lower_bounds[1] = -infinity<float>();
   upper_bounds[1] = 1.0;

   lower_bounds[2] = 0.0;
   upper_bounds[2] = infinity<float>();

   lower_bounds[3] = -infinity<float>();
   upper_bounds[3] = infinity<float>();

   lower_bounds[4] = -infinity<float>();
   upper_bounds[4] = infinity<float>();

   integral[0] = 0;
   integral[1] = 0;
   integral[2] = 0;
   integral[3] = 0;
   integral[4] = 0;

   problem.setVariableDomains( lower_bounds, upper_bounds, integral );

   auto& activities = problem.getRowActivities();
   activities = Vec<RowActivity<float>>( 3 );

   activities[0].min = -6.0;
   activities[0].max = 7.0;
   activities[0].ninfmin = 1;
   activities[0].ninfmax = 0;

   activities[1].min = -3.0;
   activities[1].max = 3.0;
   activities[1].ninfmin = 0;
   activities[1].ninfmax = 1;

   activities[2].min = 0.0;
   activities[2].max = 1.0;
   activities[2].ninfmin = 2;
   activities[2].ninfmax = 2;

   problem.setObjective( Vec<float>( 5, 1.0 ) );

   // now we test the presolver and then the core functions involved
   Reductions<float> reductions;
   Substitution<float> presolver;
   Presolve<float> presolve;
   auto& constMatrix = problem.getConstraintMatrix();

   // test the presolver; make sure it add a reduction
   PresolveResult result = presolver.execute( problem, reductions );
   REQUIRE( reductions.size() == 5 );

   // the equality row is locked
   const Reductions<float>::Reduction& reduction = reductions.getReduction( 0 );
   REQUIRE( reduction == Reductions<float>::Reduction( 0.0, 2, -4 ) );

   // the free column has locks on the bounds
   const Reductions<float>::Reduction& reduction1 =
       reductions.getReduction( 1 );
   REQUIRE( reduction1 == Reductions<float>::Reduction( 0.0, -8, 0 ) );

   // the origin row of the implied lower bounds is locked
   const Reductions<float>::Reduction& reduction2 =
       reductions.getReduction( 2 );
   REQUIRE( reduction2 == Reductions<float>::Reduction( 0.0, 0, -4 ) );

   // the origin row of the implied upper bounds is locked
   const Reductions<float>::Reduction& reduction3 =
       reductions.getReduction( 3 );
   REQUIRE( reduction3 == Reductions<float>::Reduction( 0.0, 1, -4 ) );

   // the substitution reduction is passed
   const Reductions<float>::Reduction& reduction4 =
       reductions.getReduction( 4 );
   REQUIRE( reduction4 ==
            Reductions<float>::Reduction( static_cast<float>( 2 ), -7, 0 ) );

   auto colsorted_changed_coefs = Vec<Triplet<float>>( 0 );
   auto changed_coefs =
       constMatrix.getSubstitutionChanges( 0, 2, colsorted_changed_coefs );

   // test the row sorted coefficient changes
   REQUIRE( changed_coefs.size() == 6 );
   REQUIRE( changed_coefs[0] == std::make_tuple( 0, 0, 0.0 ) );
   REQUIRE( changed_coefs[1] == std::make_tuple( 0, 3, -2.0 ) );
   REQUIRE( changed_coefs[2] == std::make_tuple( 0, 4, 2.0 ) );
   REQUIRE( changed_coefs[3] == std::make_tuple( 1, 0, 0.0 ) );
   REQUIRE( changed_coefs[4] == std::make_tuple( 1, 3, -1.0 ) );
   REQUIRE( changed_coefs[5] == std::make_tuple( 1, 4, 1.0 ) );

   // test the column sorted coefficient changes
   REQUIRE( colsorted_changed_coefs.size() == 6 );
   REQUIRE( colsorted_changed_coefs[0] == std::make_tuple( 0, 0, 0.0 ) );
   REQUIRE( colsorted_changed_coefs[1] == std::make_tuple( 1, 0, 0.0 ) );
   REQUIRE( colsorted_changed_coefs[2] == std::make_tuple( 0, 3, -2.0 ) );
   REQUIRE( colsorted_changed_coefs[3] == std::make_tuple( 1, 3, -1.0 ) );
   REQUIRE( colsorted_changed_coefs[4] == std::make_tuple( 0, 4, 2.0 ) );
   REQUIRE( colsorted_changed_coefs[5] == std::make_tuple( 1, 4, 1.0 ) );

   // test constraint bound changes
   auto& lhs = problem.getConstraintMatrix().getLeftHandSides();
   auto& rhs = problem.getConstraintMatrix().getRightHandSides();
   REQUIRE( lhs[0] == -1.0 );
   REQUIRE( rhs[1] == 1.0 );

   problem.substituteVarInObj( 0, 2 );

   // test objective changes
   auto& obj = problem.getObjective();
   REQUIRE( obj.coefficients.size() == 5 );
   REQUIRE( obj.coefficients == Vec<float>{0.0, 1.0, 1.0, 0.0, 2.0} );
   REQUIRE( obj.offset == 1.0 );
}
