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

#include "papilo/presolvers/SimplifyInequalities.hpp"
#include "papilo/external/catch/catch_amalgamated.hpp"
#include "papilo/core/PresolveMethod.hpp"
#include "papilo/core/Problem.hpp"
#include "papilo/core/ProblemBuilder.hpp"

using namespace papilo;

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
   double time = 0.0;
   int cause = -1;
   Timer t{ time };
   Message msg{};
   Problem<double> problem = setupProblemForSimplifyingInequalities();
   problem.recomputeAllActivities();
   Statistics statistics{};
   PresolveOptions presolveOptions{};
   PostsolveStorage<double> postsolve =
       PostsolveStorage<double>( problem, num, presolveOptions );
   ProblemUpdate<double> problemUpdate( problem, postsolve, statistics,
                                        presolveOptions, num, msg );
   SimplifyInequalities<double> presolvingMethod{};
   Reductions<double> reductions{};

#ifndef PAPILO_TBB
   presolveOptions.threads = 1;
#endif

   PresolveStatus presolveStatus =
       presolvingMethod.execute( problem, problemUpdate, num, reductions, t, cause );

   REQUIRE( presolveStatus == PresolveStatus::kReduced );

   REQUIRE( reductions.getReductions().size() == 6 );

   REQUIRE( reductions.getReduction( 0 ).newval == 0 );
   REQUIRE( reductions.getReduction( 0 ).row == 0 );
   REQUIRE( reductions.getReduction( 0 ).col == RowReduction::LOCKED );

   REQUIRE( reductions.getReduction( 1 ).newval == 15 );
   REQUIRE( reductions.getReduction( 1 ).row == 0 );
   REQUIRE( reductions.getReduction( 1 ).col == RowReduction::CERTIFICATE_RHS_GCD );

   REQUIRE( reductions.getReduction( 2 ).newval == 0 );
   REQUIRE( reductions.getReduction( 2 ).row == 0 );
   REQUIRE( reductions.getReduction( 2 ).col == 2 );

   REQUIRE( reductions.getReduction( 3 ).newval == 0 );
   REQUIRE( reductions.getReduction( 3 ).row == 0 );
   REQUIRE( reductions.getReduction( 3 ).col == 3 );

   REQUIRE( reductions.getReduction( 4 ).newval == 0 );
   REQUIRE( reductions.getReduction( 4 ).row == 0 );
   REQUIRE( reductions.getReduction( 4 ).col == 4 );

   REQUIRE( reductions.getReduction( 5 ).newval == 15 );
   REQUIRE( reductions.getReduction( 5 ).row == 0 );
   REQUIRE( reductions.getReduction( 5 ).col == RowReduction::RHS );

}

TEST_CASE( "simplify_inequ_doesnt_lock_more_rows", "[presolve]" )
{
   Num<double> num{};
   double time = 0.0;
   int cause = -1;
   Timer t{ time };
   Message msg{};
   Problem<double> problem = setup_simplify_ineq_reduce_rhs();
   Statistics statistics{};
   PresolveOptions presolveOptions{};
   PostsolveStorage<double> postsolve =
       PostsolveStorage<double>( problem, num, presolveOptions );
   ProblemUpdate<double> problemUpdate( problem, postsolve, statistics,
                                        presolveOptions, num, msg );
   problem.recomputeAllActivities();
   SimplifyInequalities<double> presolvingMethod{};
   Reductions<double> reductions{};

#ifndef PAPILO_TBB
   presolveOptions.threads = 1;
#endif

   PresolveStatus presolveStatus =
       presolvingMethod.execute( problem, problemUpdate, num, reductions, t, cause );

   REQUIRE( presolveStatus == PresolveStatus::kReduced );
   REQUIRE( reductions.getReduction( 2 ).newval == -275 );
   REQUIRE( reductions.getReduction( 2 ).row == 1 );
}

TEST_CASE( "simplify_inequ_doesnt_apply_lb_and_ub_on_one_row", "[presolve]" )
{
   Num<double> num{};
   double time = 0.0;
   int cause = -1;
   Timer t{ time };
   Message msg{};
   Problem<double> problem = setup_simple_problem_for_simplify_inequalities_2();
   Statistics statistics{};
   PresolveOptions presolveOptions{};
   PostsolveStorage<double> postsolve =
       PostsolveStorage<double>( problem, num, presolveOptions );
   ProblemUpdate<double> problemUpdate( problem, postsolve, statistics,
                                        presolveOptions, num, msg );
   problem.recomputeAllActivities();
   SimplifyInequalities<double> presolvingMethod{};
   Reductions<double> reductions{};

   PresolveStatus presolveStatus =
       presolvingMethod.execute( problem, problemUpdate, num, reductions, t, cause );

#ifndef PAPILO_TBB
   presolveOptions.threads = 1;
#endif


   REQUIRE( presolveStatus == PresolveStatus::kUnchanged );
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
   pb.reserve( (int) entries.size(), (int) rowNames.size(), (int) columnNames.size() );
   pb.setNumRows( (int) rowNames.size() );
   pb.setNumCols( (int) columnNames.size() );
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
   pb.reserve( (int) entries.size(), (int) rowNames.size(), (int) columnNames.size() );
   pb.setNumRows( (int) rowNames.size() );
   pb.setNumCols( (int) columnNames.size() );
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
   pb.reserve( (int) entries.size(), (int) rowNames.size(), (int) columnNames.size() );
   pb.setNumRows( (int) rowNames.size() );
   pb.setNumCols( (int) columnNames.size() );
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
