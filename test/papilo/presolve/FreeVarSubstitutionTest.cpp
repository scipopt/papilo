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
#include "papilo/core/ProblemBuilder.hpp"
#include "papilo/core/Reductions.hpp"
#include "papilo/misc/compress_vector.hpp"
#include "papilo/presolvers/FreeVarSubstitution.hpp"

#include <tuple>

//static bool
//operator==( const papilo::Reductions<float>::Reduction& lhs,
//            const papilo::Reductions<float>::Reduction& rhs )
//{
//   float max = std::max( float{1.0}, std::max( std::fabs( lhs.newval ),
//                                               std::fabs( rhs.newval ) ) );
//   bool equalityCheck = std::fabs( lhs.newval - rhs.newval ) <=
//                        ( std::numeric_limits<float>::epsilon() * max );
//
//   return equalityCheck && ( lhs.row == rhs.row ) && ( lhs.col == rhs.col );
//}

papilo::Problem<double>
setupProblemForFreeVariableSubstitution()
{
   double inf = 10000000;
   papilo::Vec<double> coefficients{ 1.0, 1.0, 1.0, 1.0, 1.0 };
   papilo::Vec<double> upperBounds{ 3.0, 1.0, inf, inf, inf };
   papilo::Vec<double> lowerBounds{ -3.0, -inf, 0.0, -inf, -inf };
   papilo::Vec<uint8_t> isIntegral{ 0,0,0,0,0 }; //2x1 + x2 >= 1
   papilo::Vec<double> rhs{ inf, 2.0, 1.0, inf, inf };
   papilo::Vec<double> lhs{ 1.0, -inf, 1.0, -inf, -inf };
   papilo::Vec<std::string> rowNames{ "A1", "A2", "A3" };
   papilo::Vec<std::string> columnNames{ "c1", "c2", "c3", "c4", "c5" };
   papilo::Vec<std::tuple<int, int, double>> entries{
       std::tuple<int, int, double>{ 0, 0, 2.0 },
       std::tuple<int, int, double>{ 0, 1, 1.0 },
       std::tuple<int, int, double>{ 1, 0, 1.0 },
       std::tuple<int, int, double>{ 1, 2, 1.0 },
       std::tuple<int, int, double>{ 2, 0, 1.0 },
       std::tuple<int, int, double>{ 2, 3, 1.0 },
       std::tuple<int, int, double>{ 2, 4, -1.0 },
   };

   papilo::ProblemBuilder<double> pb;
   pb.reserve( entries.size(), rowNames.size(), columnNames.size() );
   pb.setNumRows( rowNames.size() );
   pb.setNumCols( columnNames.size() );
   pb.setColUbAll( upperBounds );
   pb.setColLbAll( lowerBounds );
   pb.setObjAll( coefficients );
   pb.setObjOffset( 0.0 );
   pb.setColIntegralAll( isIntegral );
   pb.setRowRhsAll( rhs );
   pb.addEntryAll( entries );
   pb.setColNameAll( columnNames );
   pb.setProblemName( "matrix with free variables (3,4)" );
   papilo::Problem<double> problem = pb.build();
   return problem;
}

TEST_CASE( "test free variable detection ", "[core]" )
{
   papilo::Problem<double> problem= setupProblemForFreeVariableSubstitution();

   papilo::Num<double> num{};
   papilo::Statistics statistics{};
   papilo::PresolveOptions presolveOptions{};
   papilo::Postsolve<double> postsolve =
       papilo::Postsolve<double>( problem, num );
   auto& activities = problem.getRowActivities();
//   activities = papilo::Vec<papilo::RowActivity<float>>( 3 );
//
//   activities[0].min = -6.0;
//   activities[0].max = 7.0;
//   activities[0].ninfmin = 1;
//   activities[0].ninfmax = 0;
//
//   activities[1].min = -3.0;
//   activities[1].max = 3.0;
//   activities[1].ninfmin = 0;
//   activities[1].ninfmax = 1;
//
//   activities[2].min = 0.0;
//   activities[2].max = 1.0;
//   activities[2].ninfmin = 2;
//   activities[2].ninfmax = 2;
   papilo::ProblemUpdate<double> problemUpdate( problem, postsolve, statistics,
                                                presolveOptions, num );
   //TODO: why this class can't be found?
   papilo::FreeVarSubstitution<double> presolvingMethod{};
   papilo::Reductions<double> reductions{};


//   papilo::PresolveStatus presolveStatus =
//       presolvingMethod.execute( problem, problemUpdate, num, reductions );

//   // test the presolver; make sure it add a reduction
//   papilo::PresolveResult result = presolver.execute( problem, reductions );
//   REQUIRE( reductions.size() == 5 );
//
//   // the equality row is locked
//   const papilo::Reductions<float>::Reduction& reduction = reductions.getReduction( 0 );
//   REQUIRE( reduction == papilo::Reductions<float>::Reduction( 0.0, 2, -4 ) );
//
//   // the free column has locks on the bounds
//   const papilo::Reductions<float>::Reduction& reduction1 =
//       reductions.getReduction( 1 );
//   REQUIRE( reduction1 == papilo::Reductions<float>::Reduction( 0.0, -8, 0 ) );
//
//   // the origin row of the implied lower bounds is locked
//   papilo::Reductions<float>::Reduction& reduction2 =
//       reductions.getReduction( 2 );
//   REQUIRE( reduction2 == Reductions<float>::Reduction( 0.0, 0, -4 ) );
//
//   // the origin row of the implied upper bounds is locked
//   const papilo::Reductions<float>::Reduction& reduction3 =
//       reductions.getReduction( 3 );
//   REQUIRE( reduction3 == Reductions<float>::Reduction( 0.0, 1, -4 ) );
//
//   // the substitution reduction is passed
//   const papilo::Reductions<float>::Reduction& reduction4 =
//       reductions.getReduction( 4 );
//   REQUIRE( reduction4 ==
//                papilo::Reductions<float>::Reduction( static_cast<float>( 2 ), -7, 0 ) );
//
//   auto colsorted_changed_coefs = papilo::Vec<papilo::Triplet<float>>( 0 );
//   auto changed_coefs =
//       constMatrix.getSubstitutionChanges( 0, 2, colsorted_changed_coefs );
//
//   // test the row sorted coefficient changes
//   REQUIRE( changed_coefs.size() == 6 );
//   REQUIRE( changed_coefs[0] == std::make_tuple( 0, 0, 0.0 ) );
//   REQUIRE( changed_coefs[1] == std::make_tuple( 0, 3, -2.0 ) );
//   REQUIRE( changed_coefs[2] == std::make_tuple( 0, 4, 2.0 ) );
//   REQUIRE( changed_coefs[3] == std::make_tuple( 1, 0, 0.0 ) );
//   REQUIRE( changed_coefs[4] == std::make_tuple( 1, 3, -1.0 ) );
//   REQUIRE( changed_coefs[5] == std::make_tuple( 1, 4, 1.0 ) );
//
//   // test the column sorted coefficient changes
//   REQUIRE( colsorted_changed_coefs.size() == 6 );
//   REQUIRE( colsorted_changed_coefs[0] == std::make_tuple( 0, 0, 0.0 ) );
//   REQUIRE( colsorted_changed_coefs[1] == std::make_tuple( 1, 0, 0.0 ) );
//   REQUIRE( colsorted_changed_coefs[2] == std::make_tuple( 0, 3, -2.0 ) );
//   REQUIRE( colsorted_changed_coefs[3] == std::make_tuple( 1, 3, -1.0 ) );
//   REQUIRE( colsorted_changed_coefs[4] == std::make_tuple( 0, 4, 2.0 ) );
//   REQUIRE( colsorted_changed_coefs[5] == std::make_tuple( 1, 4, 1.0 ) );
//
//   // test constraint bound changes
//   auto& lhs = problem.getConstraintMatrix().getLeftHandSides();
//   auto& rhs = problem.getConstraintMatrix().getRightHandSides();
//   REQUIRE( lhs[0] == -1.0 );
//   REQUIRE( rhs[1] == 1.0 );
//
////   problem.substituteVarInObj( 0, 2 );
//
//   // test objective changes
//   auto& obj = problem.getObjective();
//   REQUIRE( obj.coefficients.size() == 5 );
////   REQUIRE( obj.coefficients == papilo::Vec<float>{0.0, 1.0, 1.0, 0.0, 2.0} );
//   REQUIRE( obj.offset == 1.0 );
}
