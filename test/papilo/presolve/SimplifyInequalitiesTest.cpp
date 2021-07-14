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
#include "papilo/io/MpsParser.hpp"

using namespace papilo;

Problem<double>
setupExample1ofChapter3Dot5InPresolveReductions();

Problem<double>
setupExample2ofChapter3Dot5InPresolveReductions();

Problem<double>
setupProblemForSimplifyingInequalities();

Problem<double>
setup_simplify_ineq_reduce_rhs();

Problem<double>
setup_simple_problem_for_simplify_inequalities_2();

TEST_CASE( "happy-path-simplify-inequalities", "[presolve]" )
{
   // 15x1 +15x2 +7x3 +3x4 +y1 <= 26
   Num<double> num{};
   Message msg{};
   Problem<double> problem = setupProblemForSimplifyingInequalities();
   problem.recomputeAllActivities();
   Statistics statistics{};
   PresolveOptions presolveOptions{};
   PostsolveListener<double> postsolve =
       PostsolveListener<double>( problem, num );
   ProblemUpdate<double> problemUpdate( problem, postsolve, statistics,
                                        presolveOptions, num, msg );
   SimplifyInequalities<double> presolvingMethod{};
   Reductions<double> reductions{};


   PresolveStatus presolveStatus =
       presolvingMethod.execute( problem, problemUpdate, num, reductions );

   REQUIRE( presolveStatus == PresolveStatus::kReduced );

   REQUIRE( reductions.getReductions().size() == 5 );

   REQUIRE( reductions.getReduction( 0 ).newval == 0 );
   REQUIRE( reductions.getReduction( 0 ).row == 0 );
   REQUIRE( reductions.getReduction( 0 ).col == RowReduction::LOCKED );

   REQUIRE( reductions.getReduction( 1 ).newval == 0 );
   REQUIRE( reductions.getReduction( 1 ).row == 0 );
   REQUIRE( reductions.getReduction( 1 ).col == 2 );

   REQUIRE( reductions.getReduction( 2 ).newval == 0 );
   REQUIRE( reductions.getReduction( 2 ).row == 0 );
   REQUIRE( reductions.getReduction( 2 ).col == 3 );

   REQUIRE( reductions.getReduction( 3 ).newval == 0 );
   REQUIRE( reductions.getReduction( 3 ).row == 0 );
   REQUIRE( reductions.getReduction( 3 ).col == 4 );

   REQUIRE( reductions.getReduction( 4 ).newval == 15 );
   REQUIRE( reductions.getReduction( 4 ).row == 0 );
   REQUIRE( reductions.getReduction( 4 ).col == RowReduction::RHS );

}

TEST_CASE( "simplify_inequ_doesnt_lock_more_rows", "[presolve]" )
{
   Num<double> num{};
   Message msg{};
   Problem<double> problem = setup_simplify_ineq_reduce_rhs();
   Statistics statistics{};
   PresolveOptions presolveOptions{};
   PostsolveListener<double> postsolve =
       PostsolveListener<double>( problem, num );
   ProblemUpdate<double> problemUpdate( problem, postsolve, statistics,
                                        presolveOptions, num, msg );
   problem.recomputeAllActivities();
   SimplifyInequalities<double> presolvingMethod{};
   Reductions<double> reductions{};

   PresolveStatus presolveStatus =
       presolvingMethod.execute( problem, problemUpdate, num, reductions );

   REQUIRE( presolveStatus == PresolveStatus::kReduced );
   REQUIRE( reductions.getReduction( 1 ).newval == -275 );
   REQUIRE( reductions.getReduction( 1 ).row == 1 );
}

TEST_CASE( "simplify_inequ_doesnt_apply_lb_and_ub_on_one_row", "[presolve]" )
{
   Num<double> num{};
   Message msg{};
   Problem<double> problem = setup_simple_problem_for_simplify_inequalities_2();
   Statistics statistics{};
   PresolveOptions presolveOptions{};
   PostsolveListener<double> postsolve =
       PostsolveListener<double>( problem, num );
   ProblemUpdate<double> problemUpdate( problem, postsolve, statistics,
                                        presolveOptions, num, msg );
   problem.recomputeAllActivities();
   SimplifyInequalities<double> presolvingMethod{};
   Reductions<double> reductions{};

   PresolveStatus presolveStatus =
       presolvingMethod.execute( problem, problemUpdate, num, reductions );

   REQUIRE( presolveStatus == PresolveStatus::kUnchanged );
}

TEST_CASE( "example-1-from-3.5-Presolve-Reductions-in-MIP", "[presolve]" )
{
   Num<double> num{};
   Message msg{};
   Problem<double> problem = setupExample1ofChapter3Dot5InPresolveReductions();
   Statistics statistics{};
   PresolveOptions presolveOptions{};
   PostsolveListener<double> postsolve =
       PostsolveListener<double>( problem, num );
   ProblemUpdate<double> problemUpdate( problem, postsolve, statistics,
                                        presolveOptions, num, msg );
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
   Message msg{};
   Problem<double> problem = setupExample2ofChapter3Dot5InPresolveReductions();
   Statistics statistics{};
   PresolveOptions presolveOptions{};
   PostsolveListener<double> postsolve =
       PostsolveListener<double>( problem, num );
   ProblemUpdate<double> problemUpdate( problem, postsolve, statistics,
                                        presolveOptions, num, msg );
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
   Num<double> num{};
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
   problem.getConstraintMatrix().modifyLeftHandSide( 0,num, lhs[0] );
   return problem;
}

Problem<double>
setupExample2ofChapter3Dot5InPresolveReductions()
{
   Num<double> num{};
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
   problem.getConstraintMatrix().modifyLeftHandSide( 0,num, lhs[0] );
   return problem;
}

Problem<double>
setupProblemForSimplifyingInequalities()
{
   Vec<double> coefficients{ 1.0, 1.0, 1, 1, 1 };
   Vec<double> lowerBounds{ 0.0, 0.0, 0, 0, 0 };
   Vec<double> upperBounds{ 1.0, 1.0, 1, 1, 1 };
   Vec<uint8_t> lhsInf{ 1 };
   Vec<uint8_t> isIntegral{ 1, 1, 1, 1, 0 };

   Vec<double> rhs{ 26 };
   Vec<std::string> rowNames{ "A1" };
   Vec<std::string> columnNames{ "x", "y", "z", "w", "v" };
   Vec<std::tuple<int, int, double>> entries{
       std::tuple<int, int, double>{ 0, 0, 15.0 },
       std::tuple<int, int, double>{ 0, 1, 15.0 },
       std::tuple<int, int, double>{ 0, 2, 7.0 },
       std::tuple<int, int, double>{ 0, 3, 3.0 },
       std::tuple<int, int, double>{ 0, 4, 1.0 } };

   ProblemBuilder<double> pb;
   pb.reserve( entries.size(), rowNames.size(), columnNames.size() );
   pb.setNumRows( rowNames.size() );
   pb.setNumCols( columnNames.size() );
   pb.setColLbAll( lowerBounds );
   pb.setColUbAll( upperBounds );
   pb.setObjAll( coefficients );
   pb.setObjOffset( 0.0 );
   pb.setColIntegralAll( isIntegral );
   pb.setRowRhsAll( rhs );
   pb.setRowLhsInfAll( lhsInf );
   pb.setRowLhs( 0, 1 );
   pb.addEntryAll( entries );
   pb.setColNameAll( columnNames );
   pb.setProblemName( "variables v,w,z can be deletedd" );
   Problem<double> problem = pb.build();
   return problem;
}

Problem<double>
setup_simplify_ineq_reduce_rhs()
{
   Vec<double> coefficients{ 0.0, 0.0, 0.0, 0.0 };
   Vec<double> lowerBounds{ 0.0, 0.0, 0, 0 };
   Vec<double> upperBounds{ 1.0, 1.0, 1.0, 1.0 };
   Vec<uint8_t> lhsInf{ 1, 1 };
   Vec<uint8_t> isIntegral{ 1, 1, 1, 1 };

   Vec<double> rhs{ -270, -271 };
   Vec<std::string> rowNames{ "R128", "R_test" };
   Vec<std::string> columnNames{ "C151", "C163", "C188", "C189" };
   Vec<std::tuple<int, int, double>> entries{
       std::tuple<int, int, double>{ 0, 0, -300.0 },
       std::tuple<int, int, double>{ 0, 1, -285.0 },
       std::tuple<int, int, double>{ 0, 2, -200.0 },
       std::tuple<int, int, double>{ 0, 3, -400.0 },
       std::tuple<int, int, double>{ 1, 0, -300.0 },
       std::tuple<int, int, double>{ 1, 1, -285.0 },
       std::tuple<int, int, double>{ 1, 2, -200.0 },
       std::tuple<int, int, double>{ 1, 3, -400.0 } };

   ProblemBuilder<double> pb;
   pb.reserve( entries.size(), rowNames.size(), columnNames.size() );
   pb.setNumRows( rowNames.size() );
   pb.setNumCols( columnNames.size() );
   pb.setColLbAll( lowerBounds );
   pb.setColUbAll( upperBounds );
   pb.setObjAll( coefficients );
   pb.setObjOffset( 0.0 );
   pb.setColIntegralAll( isIntegral );
   pb.setRowRhsAll( rhs );
   pb.setRowLhsInfAll( lhsInf );
   pb.setRowLhs( 0, 1 );
   pb.addEntryAll( entries );
   pb.setColNameAll( columnNames );
   pb.setProblemName( "matrix constraint 1 can be divided by 2" );
   Problem<double> problem = pb.build();
   return problem;
}

Problem<double>
setup_simple_problem_for_simplify_inequalities_2()
{
   Vec<double> coefficients{ 0.0, 0.0, 0.0, 0.0 };
   Vec<double> lowerBounds{ 0.0, 0.0, 0, 0 };
   Vec<double> upperBounds{ 1.0, 1.0, 1.0, 1.0 };
   Vec<uint8_t> lhsInf{ 0, 0 };
   Vec<uint8_t> isIntegral{ 1, 1, 1, 1 };

   Vec<double> rhs{ -270, -271 };
   Vec<double> lhs{ -800, -804 };
   Vec<std::string> rowNames{ "R128", "R_test" };
   Vec<std::string> columnNames{ "C151", "C163", "C188", "C189" };
   Vec<std::tuple<int, int, double>> entries{
       std::tuple<int, int, double>{ 0, 0, -300.0 },
       std::tuple<int, int, double>{ 0, 1, -285.0 },
       std::tuple<int, int, double>{ 0, 2, -200.0 },
       std::tuple<int, int, double>{ 0, 3, -400.0 },
       std::tuple<int, int, double>{ 1, 0, -300.0 },
       std::tuple<int, int, double>{ 1, 1, -285.0 },
       std::tuple<int, int, double>{ 1, 2, -200.0 },
       std::tuple<int, int, double>{ 1, 3, -400.0 } };

   ProblemBuilder<double> pb;
   pb.reserve( entries.size(), rowNames.size(), columnNames.size() );
   pb.setNumRows( rowNames.size() );
   pb.setNumCols( columnNames.size() );
   pb.setColLbAll( lowerBounds );
   pb.setColUbAll( upperBounds );
   pb.setObjAll( coefficients );
   pb.setObjOffset( 0.0 );
   pb.setColIntegralAll( isIntegral );
   pb.setRowRhsAll( rhs );
   pb.setRowLhsAll( lhs );
   pb.setRowLhsInfAll( lhsInf );
   pb.addEntryAll( entries );
   pb.setColNameAll( columnNames );
   pb.setProblemName( "rhs and lhs of constraint2 can be tightened" );
   Problem<double> problem = pb.build();
   return problem;
}
