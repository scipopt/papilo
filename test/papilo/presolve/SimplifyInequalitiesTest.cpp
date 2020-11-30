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

#include "papilo/presolvers/SimplifyInequalities.hpp"
#include "catch/catch.hpp"
#include "papilo/core/PresolveMethod.hpp"
#include "papilo/core/Problem.hpp"
#include "papilo/core/ProblemBuilder.hpp"

using namespace papilo;

Problem<double>
setupExample1ofChapter3Dot5InPresolveReductions();

Problem<double>
setupExample2ofChapter3Dot5InPresolveReductions();

Problem<double>
setupProblemForSimplifyingInequalities();

TEST_CASE( "happy-path-simplify-inequalities-only-greatest-divisor",
           "[presolve]" )
{
   Num<double> num{};
   Problem<double> problem = setupProblemForSimplifyingInequalities();
   Statistics statistics{};
   PresolveOptions presolveOptions{};
   Postsolve<double> postsolve = Postsolve<double>( problem, num );
   ProblemUpdate<double> problemUpdate( problem, postsolve, statistics,
                                        presolveOptions, num );
   SimplifyInequalities<double> presolvingMethod{};
   Reductions<double> reductions{};

   problem.recomputeAllActivities();

   PresolveStatus presolveStatus =
       presolvingMethod.execute( problem, problemUpdate, num, reductions );

   // x -3y + 6z = 9 -> x = 3 w; w in Z -> w - y + 2 = 6
   // potential clash with implied integer
   REQUIRE( presolveStatus == PresolveStatus::kUnchanged );
}

TEST_CASE( "example-1-from-3.5-Presolve-Reductions-in-MIP", "[presolve]" )
{
   Num<double> num{};
   Problem<double> problem = setupExample1ofChapter3Dot5InPresolveReductions();
   Statistics statistics{};
   PresolveOptions presolveOptions{};
   Postsolve<double> postsolve = Postsolve<double>( problem, num );
   ProblemUpdate<double> problemUpdate( problem, postsolve, statistics,
                                        presolveOptions, num );
   SimplifyInequalities<double> presolvingMethod{};
   Reductions<double> reductions{};

   PresolveStatus presolveStatus =
       presolvingMethod.execute( problem, problemUpdate, num, reductions );

   // TODO:
   REQUIRE( presolveStatus == PresolveStatus::kUnchanged );
}

/***
 * 1867x + 1913y = 3618894
 * 1867 & 1913 are primal and relatively primal (gcd(1867,1913)=1) and
 * 3618894/1913 in Z
 *
 * 1206 * 3618894 = 1009 (mod 1913) (modular multiplative inverse
 * x = 1913z + 1009 -> 3571571z +1913y = 1735091 -> 1867z + y = 907
 * -> y = 907; z= 0; x = 1009 (instead of 2000 b&b decisions))
 */
TEST_CASE( "example-2-from-3.5-Presolve-Reductions-in-MIP", "[presolve]" )
{
   Num<double> num{};
   Problem<double> problem = setupExample2ofChapter3Dot5InPresolveReductions();
   Statistics statistics{};
   PresolveOptions presolveOptions{};
   Postsolve<double> postsolve = Postsolve<double>( problem, num );
   ProblemUpdate<double> problemUpdate( problem, postsolve, statistics,
                                        presolveOptions, num );
   SimplifyInequalities<double> presolvingMethod{};
   Reductions<double> reductions{};

   PresolveStatus presolveStatus =
       presolvingMethod.execute( problem, problemUpdate, num, reductions );

   // TODO:
   REQUIRE( presolveStatus == PresolveStatus::kUnchanged );
}

Problem<double>
setupExample1ofChapter3Dot5InPresolveReductions()
{
   Vec<double> coefficients{ 1.0, 1.0, 1.0 };
   Vec<double> lowerBounds{ 0.0, 0.0, 0.0 };
   Vec<uint8_t> isIntegral{ 0, 1, 1 };

   Vec<double> rhs{ 9.0 };
   Vec<double> lhs{ rhs[0] };
   Vec<std::string> rowNames{
       "A1",
   };
   Vec<std::string> columnNames{ "c1", "c2", "c3" };
   Vec<std::tuple<int, int, double>> entries{
       std::tuple<int, int, double>{ 0, 0, 1.0 },
       std::tuple<int, int, double>{ 0, 1, 3.0 },
       std::tuple<int, int, double>{ 0, 2, 6.0 } };

   ProblemBuilder<double> pb;
   pb.reserve( entries.size(), rowNames.size(), columnNames.size() );
   pb.setNumRows( rowNames.size() );
   pb.setNumCols( columnNames.size() );
   pb.setColLbAll( lowerBounds );
   pb.setObjAll( coefficients );
   pb.setObjOffset( 0.0 );
   pb.setColIntegralAll( isIntegral );
   pb.setRowRhsAll( rhs );
   pb.setRowLhsAll( lhs );
   pb.addEntryAll( entries );
   pb.setColNameAll( columnNames );
   pb.setProblemName(
       "matrix Example 1 of chapter 3.5 in Presolve Reductions in MIP" );
   Problem<double> problem = pb.build();
   problem.getConstraintMatrix().modifyLeftHandSide( 0, lhs[0] );
   return problem;
}

Problem<double>
setupExample2ofChapter3Dot5InPresolveReductions()
{
   Vec<double> coefficients{ 1.0, 1.0 };
   Vec<double> lowerBounds{ 0, 0 };
   Vec<uint8_t> isIntegral{ 0, 0 };

   Vec<double> rhs{ 3618894 };
   Vec<double> lhs{ rhs[0] };
   Vec<std::string> rowNames{
       "A1",
   };
   Vec<std::string> columnNames{ "c1", "c2" };
   Vec<std::tuple<int, int, double>> entries{
       std::tuple<int, int, double>{ 0, 0, 1867.0 },
       std::tuple<int, int, double>{ 0, 1, 1913.0 } };

   ProblemBuilder<double> pb;
   pb.reserve( entries.size(), rowNames.size(), columnNames.size() );
   pb.setNumRows( rowNames.size() );
   pb.setNumCols( columnNames.size() );
   pb.setColLbAll( lowerBounds );
   pb.setObjAll( coefficients );
   pb.setObjOffset( 0.0 );
   pb.setColIntegralAll( isIntegral );
   pb.setRowRhsAll( rhs );
   pb.setRowLhsAll( lhs );
   pb.addEntryAll( entries );
   pb.setColNameAll( columnNames );
   pb.setProblemName(
       "matrix Example 2 of chapter 3.5 in Presolve Reductions in MIP" );
   Problem<double> problem = pb.build();
   problem.getConstraintMatrix().modifyLeftHandSide( 0, lhs[0] );
   return problem;
}

Problem<double>
setupProblemForSimplifyingInequalities()
{
   Vec<double> coefficients{ 1.0, 1.0 };
   Vec<double> lowerBounds{ 0.0, 0.0 };
   Vec<uint8_t> lhsInf{ 1 };
   Vec<uint8_t> isIntegral{ 0, 0 };

   Vec<double> rhs{ 6 };
   Vec<std::string> rowNames{"A1" };
   Vec<std::string> columnNames{ "c1", "c2" };
   Vec<std::tuple<int, int, double>> entries{
       std::tuple<int, int, double>{ 0, 0, 2.0 },
       std::tuple<int, int, double>{ 0, 1, 4.0 } };

   ProblemBuilder<double> pb;
   pb.reserve( entries.size(), rowNames.size(), columnNames.size() );
   pb.setNumRows( rowNames.size() );
   pb.setNumCols( columnNames.size() );
   pb.setColLbAll( lowerBounds );
   pb.setObjAll( coefficients );
   pb.setObjOffset( 0.0 );
   pb.setColIntegralAll( isIntegral );
   pb.setRowRhsAll( rhs );
   pb.setRowLhsInfAll( lhsInf );
   pb.setRowLhs(0, 1);
   pb.addEntryAll( entries );
   pb.setColNameAll( columnNames );
   pb.setProblemName( "matrix constraint 1 can be divided by 2" );
   Problem<double> problem = pb.build();
   return problem;
}
