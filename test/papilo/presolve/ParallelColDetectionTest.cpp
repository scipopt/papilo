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

#include "papilo/presolvers/ParallelColDetection.hpp"
#include "papilo/core/PresolveMethod.hpp"
#include "papilo/core/Problem.hpp"
#include "papilo/core/ProblemBuilder.hpp"
#include "papilo/core/postsolve/Postsolve.hpp"
#include "papilo/external/catch/catch_amalgamated.hpp"
#include "papilo/core/postsolve/Postsolve.hpp"

#include <cmath>

using namespace papilo;

Problem<double>
setupProblemWithParallelColumns( bool first_col_int, bool second_col_int,
                                 double factor, double ub_first_col,
                                 double ub_second_col, double lb_first_col,
                                 double lb_second_col, bool objectiveZero );

Problem<double>
setupProblemWithMultipleParallelColumns();

Problem<double>
setupProblemWithParallelSingletonColumns( bool integral, double factor,
                                          double obj_factor,
                                          int n_inactive_empty_cols );

int
countParallelReductions( Reductions<double>& reductions );

PresolveStatus
executeParallelColDetection( Problem<double>& problem,
                             Reductions<double>& reductions );

Vec<double>
postsolveParallelContIntColumns( double lb_cont_col, double ub_cont_col,
                                 double lb_int_col, double ub_int_col,
                                 double scale, double merged_value );

Vec<double>
postsolveParallelColumns( bool col1_integral, bool col2_integral,
                          bool col1_lb_inf, double lb_col1, bool col1_ub_inf,
                          double ub_col1, double lb_col2, double ub_col2,
                          double scale, double merged_value,
                          bool col2_lb_inf = false, bool col2_ub_inf = false );

Problem<double>
setupProblemWithContAndTwoIntParallelColumns( double lb, double rhs,
                                              bool equality );

Vec<double>
presolveAndPostsolveContAndTwoIntParallelColumns( Problem<double>& problem,
                                                  double lb,
                                                  double merged_value );

TEST_CASE( "parallel_col_detection_2_integer_columns", "[presolve]" )
{
   double time = 0.0;
   int cause = -1;
   Timer t{time};
   Num<double> num{};
   Message msg{};
   Problem<double> problem = setupProblemWithParallelColumns(
       true, true, 3.0, 10.0, 10.0, 0.0, 0.0, false );
   Statistics statistics{};
   PresolveOptions presolveOptions{};
   PostsolveStorage<double> postsolve =
       PostsolveStorage<double>( problem, num, presolveOptions );
   ProblemUpdate<double> problemUpdate( problem, postsolve, statistics,
                                        presolveOptions, num, msg );
   problemUpdate.checkChangedActivities();
   ParallelColDetection<double> presolvingMethod{};
   Reductions<double> reductions{};

   PresolveStatus presolveStatus =
       presolvingMethod.execute( problem, problemUpdate, num, reductions, t, cause );

   REQUIRE( presolveStatus == PresolveStatus::kReduced );
   REQUIRE( reductions.size() == 3 );
   REQUIRE( reductions.getReduction( 0 ).row == ColReduction::LOCKED );
   REQUIRE( reductions.getReduction( 0 ).col == 1 );
   REQUIRE( reductions.getReduction( 0 ).newval == 0 );

   REQUIRE( reductions.getReduction( 1 ).row == ColReduction::LOCKED );
   REQUIRE( reductions.getReduction( 1 ).col == 0 );
   REQUIRE( reductions.getReduction( 1 ).newval == 0 );

   REQUIRE( reductions.getReduction( 2 ).row == ColReduction::PARALLEL );
   REQUIRE( reductions.getReduction( 2 ).col == 1 );
   REQUIRE( reductions.getReduction( 2 ).newval == 0 );
}

TEST_CASE( "parallel_col_detection_objective_zero", "[presolve]" )
{
   double time = 0.0;
   int cause = -1;
   Timer t{time};
   Num<double> num{};
   Message msg{};
   Problem<double> problem = setupProblemWithParallelColumns(
       true, true, 0.5, 10.0, 10.0, 0.0, 0.0, true );
   Statistics statistics{};
   PresolveOptions presolveOptions{};
   PostsolveStorage<double> postsolve =
       PostsolveStorage<double>( problem, num, presolveOptions );
   ProblemUpdate<double> problemUpdate( problem, postsolve, statistics,
                                        presolveOptions, num, msg );
   problemUpdate.checkChangedActivities();
   ParallelColDetection<double> presolvingMethod{};
   Reductions<double> reductions{};

   PresolveStatus presolveStatus =
       presolvingMethod.execute( problem, problemUpdate, num, reductions, t, cause );

   REQUIRE( presolveStatus == PresolveStatus::kReduced );
   REQUIRE( reductions.size() == 3 );
   REQUIRE( reductions.getReduction( 0 ).row == ColReduction::LOCKED );
   REQUIRE( reductions.getReduction( 0 ).col == 0 );
   REQUIRE( reductions.getReduction( 0 ).newval == 0 );

   REQUIRE( reductions.getReduction( 1 ).row == ColReduction::LOCKED );
   REQUIRE( reductions.getReduction( 1 ).col == 1 );
   REQUIRE( reductions.getReduction( 1 ).newval == 0 );

   REQUIRE( reductions.getReduction( 2 ).row == ColReduction::PARALLEL );
   REQUIRE( reductions.getReduction( 2 ).col == 0 );
   REQUIRE( reductions.getReduction( 2 ).newval == 1 );
}


TEST_CASE( "parallel_col_detection_2_continuous_columns", "[presolve]" )
{
   double time = 0.0;
   int cause = -1;
   Timer t{time};
   Num<double> num{};
   Message msg{};
   Problem<double> problem = setupProblemWithParallelColumns(
       false, false, 2.0, 10.0, 10.0, 0.0, 0.0, false );
   Statistics statistics{};
   PresolveOptions presolveOptions{};
   PostsolveStorage<double> postsolve =
       PostsolveStorage<double>( problem, num, presolveOptions );
   ProblemUpdate<double> problemUpdate( problem, postsolve, statistics,
                                        presolveOptions, num, msg );
   problemUpdate.checkChangedActivities();
   ParallelColDetection<double> presolvingMethod{};
   Reductions<double> reductions{};

   PresolveStatus presolveStatus =
       presolvingMethod.execute( problem, problemUpdate, num, reductions, t, cause );

   REQUIRE( presolveStatus == PresolveStatus::kReduced );
   REQUIRE( reductions.size() == 3 );
   REQUIRE( reductions.getReduction( 0 ).row == ColReduction::LOCKED );
   REQUIRE( reductions.getReduction( 0 ).col == 1 );
   REQUIRE( reductions.getReduction( 0 ).newval == 0 );

   REQUIRE( reductions.getReduction( 1 ).row == ColReduction::LOCKED );
   REQUIRE( reductions.getReduction( 1 ).col == 0 );
   REQUIRE( reductions.getReduction( 1 ).newval == 0 );

   REQUIRE( reductions.getReduction( 2 ).row == ColReduction::PARALLEL );
   REQUIRE( reductions.getReduction( 2 ).col == 1 );
   REQUIRE( reductions.getReduction( 2 ).newval == 0 );
}

TEST_CASE( "parallel_col_detection_int_cont_merge_possible", "[presolve]" )
{
   double time = 0.0;
   int cause = -1;
   Timer t{time};
   Num<double> num{};
   Message msg{};
   Problem<double> problem = setupProblemWithParallelColumns(
       true, false, 2.0, 10.0, 10.0, 0.0, 0.0, false );
   Statistics statistics{};
   PresolveOptions presolveOptions{};
   PostsolveStorage<double> postsolve =
       PostsolveStorage<double>( problem, num, presolveOptions );
   ProblemUpdate<double> problemUpdate( problem, postsolve, statistics,
                                        presolveOptions, num, msg );
   problemUpdate.checkChangedActivities();
   ParallelColDetection<double> presolvingMethod{};
   Reductions<double> reductions{};

   PresolveStatus presolveStatus =
       presolvingMethod.execute( problem, problemUpdate, num, reductions, t, cause );

   REQUIRE( presolveStatus == PresolveStatus::kReduced );
   REQUIRE( reductions.size() == 4 );
   REQUIRE( reductions.getReduction( 0 ).row == ColReduction::LOCKED );
   REQUIRE( reductions.getReduction( 0 ).col == 1 );
   REQUIRE( reductions.getReduction( 0 ).newval == 0 );

   REQUIRE( reductions.getReduction( 1 ).row == ColReduction::LOCKED );
   REQUIRE( reductions.getReduction( 1 ).col == 0 );
   REQUIRE( reductions.getReduction( 1 ).newval == 0 );

   REQUIRE( reductions.getReduction( 2 ).row == ColReduction::BOUNDS_LOCKED );
   REQUIRE( reductions.getReduction( 2 ).col == 1 );
   REQUIRE( reductions.getReduction( 2 ).newval == 0 );

   REQUIRE( reductions.getReduction( 3 ).row == ColReduction::PARALLEL );
   REQUIRE( reductions.getReduction( 3 ).col == 1 );
   REQUIRE( reductions.getReduction( 3 ).newval == 0 );
}

TEST_CASE( "parallel_col_detection_cont_int_merge_possible", "[presolve]" )
{
   double time = 0.0;
   int cause = -1;
   Timer t{time};
   Num<double> num{};
   Message msg{};
   Problem<double> problem = setupProblemWithParallelColumns(
       false, true, 0.5, 10.0, 10.0, 0.0, 0.0, false );
   Statistics statistics{};
   PresolveOptions presolveOptions{};
   PostsolveStorage<double> postsolve =
       PostsolveStorage<double>( problem, num, presolveOptions );
   ProblemUpdate<double> problemUpdate( problem, postsolve, statistics,
                                        presolveOptions, num, msg );
   problemUpdate.checkChangedActivities();
   ParallelColDetection<double> presolvingMethod{};
   Reductions<double> reductions{};

   PresolveStatus presolveStatus =
       presolvingMethod.execute( problem, problemUpdate, num, reductions, t, cause );

   REQUIRE( presolveStatus == PresolveStatus::kReduced );
   REQUIRE( reductions.size() == 4 );
   REQUIRE( reductions.getReduction( 0 ).row == ColReduction::LOCKED );
   REQUIRE( reductions.getReduction( 0 ).col == 0 );
   REQUIRE( reductions.getReduction( 0 ).newval == 0 );

   REQUIRE( reductions.getReduction( 1 ).row == ColReduction::LOCKED );
   REQUIRE( reductions.getReduction( 1 ).col == 1 );
   REQUIRE( reductions.getReduction( 1 ).newval == 0 );

   REQUIRE( reductions.getReduction( 2 ).row == ColReduction::BOUNDS_LOCKED );
   REQUIRE( reductions.getReduction( 2 ).col == 0 );
   REQUIRE( reductions.getReduction( 2 ).newval == 0 );

   REQUIRE( reductions.getReduction( 3 ).row == ColReduction::PARALLEL );
   REQUIRE( reductions.getReduction( 3 ).col == 0 );
   REQUIRE( reductions.getReduction( 3 ).newval == 1 );
}

TEST_CASE( "parallel_col_detection_cont_int_merge_failed", "[presolve]" )
{
   double time = 0.0;
   int cause = -1;
   Timer t{time};
   Num<double> num{};
   Message msg{};
   Problem<double> problem = setupProblemWithParallelColumns(
       false, true, 1.0, 0.9, 10.0, 0.0, 0.0, false );
   Statistics statistics{};
   PresolveOptions presolveOptions{};
   PostsolveStorage<double> postsolve =
       PostsolveStorage<double>( problem, num, presolveOptions );
   ProblemUpdate<double> problemUpdate( problem, postsolve, statistics,
                                        presolveOptions, num, msg );
   problemUpdate.checkChangedActivities();
   ParallelColDetection<double> presolvingMethod{};
   Reductions<double> reductions{};

   PresolveStatus presolveStatus =
       presolvingMethod.execute( problem, problemUpdate, num, reductions, t, cause );

   REQUIRE( presolveStatus == PresolveStatus::kUnchanged );
}

TEST_CASE( "parallel_col_detection_int_cont_merge_failed", "[presolve]" )
{
   double time = 0.0;
   int cause = -1;
   Timer t{time};
   Num<double> num{};
   Message msg{};
   Problem<double> problem = setupProblemWithParallelColumns(
       true, false, 1.0, 10.0, 0.9, 0.0, 0.0, false );
   Statistics statistics{};
   PresolveOptions presolveOptions{};
   PostsolveStorage<double> postsolve =
       PostsolveStorage<double>( problem, num, presolveOptions );
   ProblemUpdate<double> problemUpdate( problem, postsolve, statistics,
                                        presolveOptions, num, msg );
   problemUpdate.checkChangedActivities();
   ParallelColDetection<double> presolvingMethod{};
   Reductions<double> reductions{};

   PresolveStatus presolveStatus =
       presolvingMethod.execute( problem, problemUpdate, num, reductions, t, cause );

   REQUIRE( presolveStatus == PresolveStatus::kUnchanged );
}

TEST_CASE( "parallel_col_detection_int_merge_failed_hole", "[presolve]" )
{
   double time = 0.0;
   int cause = -1;
   Timer t{time};
   Num<double> num{};
   Message msg{};
   Problem<double> problem = setupProblemWithParallelColumns(
       true, true, 1.33333, 5.0, 7, 3.0, 5.0, false );
   Statistics statistics{};
   PresolveOptions presolveOptions{};
   PostsolveStorage<double> postsolve =
       PostsolveStorage<double>( problem, num, presolveOptions );
   ProblemUpdate<double> problemUpdate( problem, postsolve, statistics,
                                        presolveOptions, num, msg );
   problemUpdate.checkChangedActivities();
   ParallelColDetection<double> presolvingMethod{};
   Reductions<double> reductions{};

   PresolveStatus presolveStatus =
       presolvingMethod.execute( problem, problemUpdate, num, reductions, t, cause );

   REQUIRE( presolveStatus == PresolveStatus::kUnchanged );
}

TEST_CASE( "parallel_col_detection_obj_not_parallel", "[presolve]" )
{
   double time = 0.0;
   int cause = -1;
   Timer t{time};
   Num<double> num{};
   Message msg{};
   Vec<double> obj = { 3, 2 };
   Problem<double> problem = setupProblemWithParallelColumns(
       true, true, 1.0, 10.0, 10.0, 0.0, 0.0, false );
   problem.setObjective( obj, 0 );
   Statistics statistics{};
   PresolveOptions presolveOptions{};
   PostsolveStorage<double> postsolve =
       PostsolveStorage<double>( problem, num, presolveOptions );
   ProblemUpdate<double> problemUpdate( problem, postsolve, statistics,
                                        presolveOptions, num, msg );
   problemUpdate.checkChangedActivities();
   ParallelColDetection<double> presolvingMethod{};
   Reductions<double> reductions{};

   PresolveStatus presolveStatus =
       presolvingMethod.execute( problem, problemUpdate, num, reductions, t, cause );

   REQUIRE( presolveStatus == PresolveStatus::kUnchanged );
}

TEST_CASE( "parallel_col_detection_multiple_parallel_columns", "[presolve]" )
{
   double time = 0.0;
   int cause = -1;
   Timer t{time};
   Num<double> num{};
   Message msg{};
   Problem<double> problem = setupProblemWithMultipleParallelColumns();
   Statistics statistics{};
   PresolveOptions presolveOptions{};
   PostsolveStorage<double> postsolve =
       PostsolveStorage<double>( problem, num, presolveOptions );
   ProblemUpdate<double> problemUpdate( problem, postsolve, statistics,
                                        presolveOptions, num, msg );
   problemUpdate.checkChangedActivities();
   ParallelColDetection<double> presolvingMethod{};
   Reductions<double> reductions{};

   PresolveStatus presolveStatus =
       presolvingMethod.execute( problem, problemUpdate, num, reductions, t, cause );

   REQUIRE( presolveStatus == PresolveStatus::kReduced );
   REQUIRE( reductions.size() == 13 );
   REQUIRE( reductions.getTransactions().size() == 2 );
   REQUIRE( reductions.getTransactions()[0].start == 0 );
   REQUIRE( reductions.getTransactions()[0].end == 10 );
   REQUIRE( reductions.getTransactions()[1].start == 10 );
   REQUIRE( reductions.getTransactions()[1].end == 13 );
   int smallest_cont_column = 4;
   int remaining_integer_col = 1;
   Vec<int> locked_rows = { smallest_cont_column, 6, 2, 0,
                            remaining_integer_col };
   for( int i = 0; i < (int) locked_rows.size(); i++ )
   {

      REQUIRE( reductions.getReduction( i ).row == ColReduction::LOCKED );
      REQUIRE( reductions.getReduction( i ).col == locked_rows[i] );
      REQUIRE( reductions.getReduction( i ).newval == 0 );
   }
   REQUIRE( reductions.getReduction( 5 ).row == ColReduction::BOUNDS_LOCKED );
   REQUIRE( reductions.getReduction( 5 ).col == smallest_cont_column );
   REQUIRE( reductions.getReduction( 5 ).newval == 0 );

   for( int i = 1; i < (int) locked_rows.size() - 1; i++ )
   {
      REQUIRE( reductions.getReduction( 5 + i ).row == ColReduction::PARALLEL );
      REQUIRE( reductions.getReduction( 5 + i ).col == locked_rows[i] );
      REQUIRE( reductions.getReduction( 5 + i ).newval ==
               smallest_cont_column );
   }

   REQUIRE( reductions.getReduction( 10 ).row == ColReduction::LOCKED );
   REQUIRE( reductions.getReduction( 10 ).col == 5 );
   REQUIRE( reductions.getReduction( 10 ).newval == 0 );

   REQUIRE( reductions.getReduction( 11 ).row == ColReduction::LOCKED );
   REQUIRE( reductions.getReduction( 11 ).col == remaining_integer_col );
   REQUIRE( reductions.getReduction( 11 ).newval == 0 );

   REQUIRE( reductions.getReduction( 12 ).row == ColReduction::PARALLEL );
   REQUIRE( reductions.getReduction( 12 ).col == 5 );
   REQUIRE( reductions.getReduction( 12 ).newval == remaining_integer_col );
}

TEST_CASE( "parallel_col_detection_cont_and_two_int_columns", "[presolve]" )
{
   // x_0 continuous, x_1 and x_2 integer, all in [0, 10] and parallel: x_0 is
   // merged into x_1 first, which turns x_1 continuous, before x_2 is merged
   // into x_1. The resulting column must be continuous (issue #93)
   Problem<double> problem =
       setupProblemWithContAndTwoIntParallelColumns( 0, 20, false );
   Vec<double> values =
       presolveAndPostsolveContAndTwoIntParallelColumns( problem, 0, 0 );
   REQUIRE( values[0] == 0 );
   REQUIRE( values[1] == 0 );
   REQUIRE( values[2] == 0 );
}

TEST_CASE( "parallel_col_detection_cont_and_two_int_columns_fractional",
           "[presolve]" )
{
   // x_0 + x_1 + x_2 = -0.2 is feasible, e.g., (0.8, 0, -1), hence the merged
   // column must not be integral (issue #93)
   Problem<double> problem =
       setupProblemWithContAndTwoIntParallelColumns( -10, -0.2, true );
   Vec<double> values =
       presolveAndPostsolveContAndTwoIntParallelColumns( problem, -10, -0.2 );
   REQUIRE( values[1] == std::round( values[1] ) );
   REQUIRE( values[2] == std::round( values[2] ) );
   REQUIRE( values[0] + values[1] + values[2] == Catch::Approx( -0.2 ) );
}

TEST_CASE( "parallel_col_detection_singleton_integer_columns", "[presolve]" )
{
   // knapsack row with two integer columns with the same coefficient and
   // objective (as in p2m2p1m1p0n100), only one entry per column (#115)
   Problem<double> problem =
       setupProblemWithParallelSingletonColumns( true, 1.0, 1.0, 0 );
   Reductions<double> reductions{};

   PresolveStatus presolveStatus =
       executeParallelColDetection( problem, reductions );

   REQUIRE( presolveStatus == PresolveStatus::kReduced );
   REQUIRE( countParallelReductions( reductions ) == 1 );
}

TEST_CASE( "parallel_col_detection_singleton_continuous_columns",
           "[presolve]" )
{
   // continuous singleton columns are parallel for any scale if the
   // objective is scaled by the same factor (#115)
   Problem<double> problem =
       setupProblemWithParallelSingletonColumns( false, 2.5, 2.5, 0 );
   Reductions<double> reductions{};

   PresolveStatus presolveStatus =
       executeParallelColDetection( problem, reductions );

   REQUIRE( presolveStatus == PresolveStatus::kReduced );
   REQUIRE( countParallelReductions( reductions ) == 1 );
}

TEST_CASE( "parallel_col_detection_singleton_obj_not_parallel", "[presolve]" )
{
   // same coefficient in the row but a different objective: the singleton
   // columns are not parallel
   Problem<double> problem =
       setupProblemWithParallelSingletonColumns( true, 1.0, 2.0, 0 );
   Reductions<double> reductions{};

   PresolveStatus presolveStatus =
       executeParallelColDetection( problem, reductions );

   REQUIRE( presolveStatus == PresolveStatus::kUnchanged );
   REQUIRE( countParallelReductions( reductions ) == 0 );
}

TEST_CASE( "parallel_col_detection_singleton_with_inactive_empty_columns",
           "[presolve]" )
{
   // inactive columns have no entries and end up in a common bucket; they
   // must neither be merged nor prevent detecting the parallel singletons
   Problem<double> problem =
       setupProblemWithParallelSingletonColumns( true, 1.0, 1.0, 2 );
   Reductions<double> reductions{};

   PresolveStatus presolveStatus =
       executeParallelColDetection( problem, reductions );

   REQUIRE( presolveStatus == PresolveStatus::kReduced );
   REQUIRE( countParallelReductions( reductions ) == 1 );
   for( int i = 0; i < (int) reductions.size(); ++i )
      REQUIRE( reductions.getReduction( i ).col < 2 );
}

PresolveStatus
executeParallelColDetection( Problem<double>& problem,
                             Reductions<double>& reductions )
{
   double time = 0.0;
   int cause = -1;
   Timer t{ time };
   Num<double> num{};
   Message msg{};
   Statistics statistics{};
   PresolveOptions presolveOptions{};
   PostsolveStorage<double> postsolve =
       PostsolveStorage<double>( problem, num, presolveOptions );
   ProblemUpdate<double> problemUpdate( problem, postsolve, statistics,
                                        presolveOptions, num, msg );
   problemUpdate.checkChangedActivities();
   ParallelColDetection<double> presolvingMethod{};
   return presolvingMethod.execute( problem, problemUpdate, num, reductions, t,
                                    cause );
}

int
countParallelReductions( Reductions<double>& reductions )
{
   int count = 0;
   for( int i = 0; i < (int) reductions.size(); ++i )
      if( reductions.getReduction( i ).row == ColReduction::PARALLEL )
         ++count;
   return count;
}

Problem<double>
setupProblemWithParallelSingletonColumns( bool integral, double factor,
                                          double obj_factor,
                                          int n_inactive_empty_cols )
{
   // min 3 c1 + 3 obj_factor c2 + 7 c3
   // s.t. 3 c1 + 3 factor c2 + 7 c3 >= 5
   //      c4 + c3 <= 1
   // the first two columns are singletons, the last n_inactive_empty_cols
   // columns are empty and marked as inactive
   const int ncols = 4 + n_inactive_empty_cols;
   Vec<double> obj( ncols, 0.0 );
   obj[0] = 3.0;
   obj[1] = 3.0 * obj_factor;
   obj[2] = 7.0;
   obj[3] = 1.0;
   Vec<double> lowerBounds( ncols, 0.0 );
   Vec<double> upperBounds( ncols, 1.0 );
   upperBounds[0] = 4.0;
   upperBounds[1] = 4.0;
   Vec<uint8_t> isIntegral( ncols, integral );
   Vec<std::string> columnNames;
   for( int i = 0; i < ncols; ++i )
      columnNames.emplace_back( "c" + std::to_string( i + 1 ) );

   Vec<std::tuple<int, int, double>> entries{
       std::tuple<int, int, double>{ 0, 0, 3.0 },
       std::tuple<int, int, double>{ 0, 1, 3.0 * factor },
       std::tuple<int, int, double>{ 0, 2, 7.0 },
       std::tuple<int, int, double>{ 1, 2, 1.0 },
       std::tuple<int, int, double>{ 1, 3, 1.0 },
   };

   ProblemBuilder<double> pb;
   pb.reserve( entries.size(), 2, ncols );
   pb.setNumRows( 2 );
   pb.setNumCols( ncols );
   pb.setColUbAll( upperBounds );
   pb.setColLbAll( lowerBounds );
   pb.setObjAll( obj );
   pb.setObjOffset( 0.0 );
   pb.setColIntegralAll( isIntegral );
   pb.setRowLhs( 0, 5.0 );
   pb.setRowRhsInf( 0, true );
   pb.setRowLhsInf( 1, true );
   pb.setRowRhs( 1, 1.0 );
   pb.addEntryAll( entries );
   pb.setColNameAll( columnNames );
   pb.setProblemName( "matrix with parallel singleton columns (1 and 2)" );
   Problem<double> problem = pb.build();
   for( int i = 4; i < ncols; ++i )
      problem.getColFlags()[i].set( ColFlag::kInactive );
   return problem;
}

TEST_CASE( "parallel_col_postsolve_cont_int_recompute_after_rounding",
           "[presolve]" )
{
   // x_int + 0.2 * x_cont = 5.8: x_cont is clamped to 25, x_int = 0.8 is
   // rounded up to 1, hence x_cont has to be recomputed as 24
   Vec<double> values = postsolveParallelContIntColumns( 0, 25, 0, 1, 0.2, 5.8 );
   REQUIRE( values[0] == Catch::Approx( 24 ) );
   REQUIRE( values[1] == 1 );

   // x_int + 0.25 * x_cont = 6.75: x_cont is clamped to 25, x_int = 0.5 is
   // rounded up to 1, hence x_cont has to be recomputed as 23
   values = postsolveParallelContIntColumns( 16.999999999999989, 25, 0, 1, 0.25,
                                             6.75 );
   REQUIRE( values[0] == Catch::Approx( 23 ) );
   REQUIRE( values[1] == 1 );
}

Vec<double>
postsolveParallelContIntColumns( double lb_cont_col, double ub_cont_col,
                                 double lb_int_col, double ub_int_col,
                                 double scale, double merged_value )
{
   ProblemBuilder<double> pb;
   pb.reserve( 0, 0, 2 );
   pb.setNumRows( 0 );
   pb.setNumCols( 2 );
   pb.setColLbAll( { lb_cont_col, lb_int_col } );
   pb.setColUbAll( { ub_cont_col, ub_int_col } );
   pb.setObjAll( { 0.0, 0.0 } );
   pb.setColIntegralAll( { 0, 1 } );
   pb.setColNameAll( Vec<std::string>{ "c1", "i2" } );
   Problem<double> problem = pb.build();

   Num<double> num{};
   Message msg{};
   PostsolveStorage<double> storage{ problem, num, PresolveOptions{} };
   storage.storeParallelCols( 0, false, false, lb_cont_col, false, ub_cont_col,
                              1, true, false, lb_int_col, false, ub_int_col,
                              scale );

   // the merged column takes over the index of the integer column
   Solution<double> reduced_solution{ Vec<double>{ 0, merged_value } };
   Solution<double> original_solution{};
   Postsolve<double> postsolve{ msg, num };
   auto status = postsolve.undo( reduced_solution, original_solution, storage );
   REQUIRE( ( status == PostsolveStatus::kOk) );
   return original_solution.primal;
}

TEST_CASE( "parallel_col_postsolve_int_int_not_on_bounds", "[presolve]" )
{
   // x_2 + scale * x_1 = merged_value with two integer columns: the postsolved
   // values are feasible but neither column needs to be at one of its bounds
   // (issue #95)
   struct Instance
   {
      double lb1, ub1, lb2, ub2, scale, merged_value;
   };
   const Vec<Instance> instances{
       { -8, 3, -5, 3, -1, 0 }, // fiball
       { 0, 2, -2, 3, -3, -3 }, // nu25-pr12
       { 0, 5, 0, 5, 2, 6 },    // neos-3024952-loue
       // only feasible integer solution is (2, 1), i.e., in the interior
       { 0, 4, 0, 2, 3, 7 } };

   for( const Instance& i : instances )
   {
      Vec<double> values =
          postsolveParallelColumns( true, true, false, i.lb1, false, i.ub1,
                                    i.lb2, i.ub2, i.scale, i.merged_value );
      REQUIRE( values[0] >= i.lb1 );
      REQUIRE( values[0] <= i.ub1 );
      REQUIRE( values[1] >= i.lb2 );
      REQUIRE( values[1] <= i.ub2 );
      REQUIRE( values[0] == std::round( values[0] ) );
      REQUIRE( values[1] == std::round( values[1] ) );
      REQUIRE( values[1] + i.scale * values[0] ==
               Catch::Approx( i.merged_value ) );
   }
}

TEST_CASE( "parallel_col_postsolve_cont_free_col1_col2_on_bounds",
           "[presolve]" )
{
   // x_2 + x_1 = 5 with x_1 free and x_2 in [0, 10]: x_2 is set to its lower
   // bound, which has to be detected using the bound flags of x_2 (issue #95)
   Vec<double> values = postsolveParallelColumns( false, false, true, 0, true,
                                                  0, 0, 10, 1, 5 );
   REQUIRE( values[0] == Catch::Approx( 5 ) );
   REQUIRE( values[1] == 0 );
}

TEST_CASE( "parallel_col_postsolve_int_into_cont", "[presolve]" )
{
   // x_2 + scale * x_1 = merged_value with x_1 integer and x_2 continuous:
   // x_1 has to be rounded towards the violated bound of x_2 since
   // merged_value is fractional (issue #93)
   struct Instance
   {
      double lb1, ub1, lb2, ub2, scale, merged_value, x1, x2;
   };
   const Vec<Instance> instances{
       // x_2 = 25.3 > 20: x_1 = ceil(5.3) = 6
       { 0, 10, 0, 20, 1, 25.3, 6, 19.3 },
       // x_2 = -5.3 < 0: x_1 = floor(-5.3) = -6
       { -10, 10, 0, 20, 1, -5.3, -6, 0.7 },
       // x_2 = -15.3 < 0: x_1 = ceil(7.65) = 8
       { 0, 10, 0, 20, -2, -15.3, 8, 0.7 } };

   for( const Instance& i : instances )
   {
      Vec<double> values =
          postsolveParallelColumns( true, false, false, i.lb1, false, i.ub1,
                                    i.lb2, i.ub2, i.scale, i.merged_value );
      REQUIRE( values[0] == i.x1 );
      REQUIRE( values[1] == Catch::Approx( i.x2 ) );
   }
}

TEST_CASE( "parallel_col_postsolve_int_into_cont_infinite_bounds",
           "[presolve]" )
{
   // x_2 + x_1 = merged_value with x_1 in [0, 10] integer and x_2 continuous
   // with infinite bounds: this happens if a free continuous column was merged
   // into x_2 earlier in the same round (issue #93)

   // x_2 in (-inf, 20]: x_2 = 25.3 > 20, hence x_1 = ceil(5.3) = 6
   Vec<double> values = postsolveParallelColumns(
       true, false, false, 0, false, 10, 0, 20, 1, 25.3, true, false );
   REQUIRE( values[0] == 6 );
   REQUIRE( values[1] == Catch::Approx( 19.3 ) );

   // x_2 in (-inf, 20]: the lower bound must be ignored
   values = postsolveParallelColumns( true, false, false, 0, false, 10, 0, 20,
                                      1, -5.3, true, false );
   REQUIRE( values[0] == 0 );
   REQUIRE( values[1] == Catch::Approx( -5.3 ) );

   // x_2 free: x_1 stays at 0
   values = postsolveParallelColumns( true, false, false, 0, false, 10, 0, 0, 1,
                                      -5.3, true, true );
   REQUIRE( values[0] == 0 );
   REQUIRE( values[1] == Catch::Approx( -5.3 ) );
}

Vec<double>
postsolveParallelColumns(bool col1_integral, bool col2_integral,
                          bool col1_lb_inf, double lb_col1, bool col1_ub_inf,
                          double ub_col1, double lb_col2, double ub_col2,
                          double scale, double merged_value, bool col2_lb_inf,
                          bool col2_ub_inf )
{
   ProblemBuilder<double> pb;
   pb.reserve( 0, 0, 2 );
   pb.setNumRows( 0 );
   pb.setNumCols( 2 );
   pb.setColLbAll( { lb_col1, lb_col2 } );
   pb.setColUbAll( { ub_col1, ub_col2 } );
   pb.setColLbInf( 0, col1_lb_inf );
   pb.setColUbInf( 0, col1_ub_inf );
   pb.setColLbInf( 1, col2_lb_inf );
   pb.setColUbInf( 1, col2_ub_inf );
   pb.setObjAll( { 0.0, 0.0 } );
   pb.setColIntegralAll( { col1_integral, col2_integral } );
   pb.setColNameAll( Vec<std::string>{ "x1", "x2" } );
   Problem<double> problem = pb.build();

   Num<double> num{};
   Message msg{};
   PostsolveStorage<double> storage{ problem, num, PresolveOptions{} };
   storage.storeParallelCols( 0, col1_integral, col1_lb_inf, lb_col1,
                              col1_ub_inf, ub_col1, 1, col2_integral,
                              col2_lb_inf, lb_col2, col2_ub_inf, ub_col2,
                              scale );

   // the merged column takes over the index of the second column
   Solution<double> reduced_solution{ Vec<double>{ 0, merged_value } };
   Solution<double> original_solution{};
   Postsolve<double> postsolve{ msg, num };
   auto status = postsolve.undo( reduced_solution, original_solution, storage );
   REQUIRE( (status == PostsolveStatus::kOk) );

   return original_solution.primal;
}

Problem<double>
setupProblemWithParallelColumns( bool first_col_int, bool second_col_int,
                                 double factor, double ub_first_col,
                                 double ub_second_col, double lb_first_col,
                                 double lb_second_col, bool objectiveZero )
{
   Vec<double> coefficients{ ! objectiveZero ? 1.0 : 0.0,
                             ! objectiveZero ? 1.0 * factor : 0.0 };
   Vec<double> lowerBounds{ lb_first_col, lb_second_col };
   Vec<double> upperBounds{ ub_first_col, ub_second_col };
   Vec<uint8_t> isIntegral{ first_col_int, second_col_int };

   Vec<double> rhs{ 1.0, 2.0 };
   Vec<std::string> rowNames{
       "A1",
       "A2",
   };
   Vec<std::string> columnNames{ "c1", "c2" };
   Vec<std::tuple<int, int, double>> entries{
       std::tuple<int, int, double>{ 0, 0, 1.0 },
       std::tuple<int, int, double>{ 0, 1, 1.0 * factor },
       std::tuple<int, int, double>{ 1, 0, 2.0 },
       std::tuple<int, int, double>{ 1, 1, 2.0 * factor },
   };

   ProblemBuilder<double> pb;
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
   pb.setProblemName( "matrix with parallel columns (1 and 2)" );
   Problem<double> problem = pb.build();
   return problem;
}

Problem<double>
setupProblemWithMultipleParallelColumns()
{
   Vec<double> coefficients{ 3.0, 2.0, 3.0, 3.0, 1.0, 4.0, 2.0 };
   Vec<double> lowerBounds{ 0, 1, 2, 3, 4, 5, 6 };
   Vec<double> upperBounds{ 2, 3, 4, 5, 6, 7, 8 };
   Vec<uint8_t> isIntegral{ 0, 1, 0, 1, 0, 1, 0 };

   Vec<double> rhs{ 1.0, 2.0 };
   Vec<std::string> rowNames{
       "A1",
       "A2",
   };
   Vec<std::string> columnNames{ "c1", "i2", "c3", "i4", "c5", "i6", "c7" };
   Vec<std::tuple<int, int, double>> entries{
       std::tuple<int, int, double>{ 0, 0, 3.0 },
       std::tuple<int, int, double>{ 1, 0, 6.0 },

       std::tuple<int, int, double>{ 0, 1, 2.0 },
       std::tuple<int, int, double>{ 1, 1, 4.0 },

       std::tuple<int, int, double>{ 0, 2, 3.0 },
       std::tuple<int, int, double>{ 1, 2, 6.0 },

       std::tuple<int, int, double>{ 0, 3, 3.0 },
       std::tuple<int, int, double>{ 1, 3, 6.0 },

       std::tuple<int, int, double>{ 0, 4, 1.0 },
       std::tuple<int, int, double>{ 1, 4, 2.0 },

       std::tuple<int, int, double>{ 0, 5, 4.0 },
       std::tuple<int, int, double>{ 1, 5, 8.0 },

       std::tuple<int, int, double>{ 0, 6, 2.0 },
       std::tuple<int, int, double>{ 1, 6, 4.0 },
   };

   ProblemBuilder<double> pb;
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
   pb.setProblemName( "matrix with multiple parallel columns" );
   Problem<double> problem = pb.build();
   return problem;
}

Problem<double>
setupProblemWithContAndTwoIntParallelColumns( double lb, double rhs,
                                              bool equality )
{
   ProblemBuilder<double> pb;
   pb.reserve( 6, 2, 3 );
   pb.setNumRows( 2 );
   pb.setNumCols( 3 );
   pb.setColLbAll( { lb, lb, lb } );
   pb.setColUbAll( { 10, 10, 10 } );
   pb.setColIntegralAll( { 0, 1, 1 } );
   pb.setObjAll( { 1, 1, 1 } );
   pb.setRowRhsAll( { rhs, 2 * rhs } );
   if( equality )
      pb.setRowLhsAll( { rhs, 2 * rhs } );
   pb.addEntryAll( { { 0, 0, 1 }, { 0, 1, 1 }, { 0, 2, 1 },
                     { 1, 0, 2 }, { 1, 1, 2 }, { 1, 2, 2 } } );
   pb.setColNameAll( Vec<std::string>{ "c0", "i1", "i2" } );
   pb.setProblemName( "continuous and two integer parallel columns" );
   return pb.build();
}

Vec<double>
presolveAndPostsolveContAndTwoIntParallelColumns( Problem<double>& problem,
                                                  double lb,
                                                  double merged_value )
{
   double time = 0.0;
   int cause = -1;
   Timer t{ time };
   Num<double> num{};
   Message msg{};
   Statistics statistics{};
   PresolveOptions presolveOptions{};
   PostsolveStorage<double> postsolve{ problem, num, presolveOptions };
   ProblemUpdate<double> problemUpdate( problem, postsolve, statistics,
                                        presolveOptions, num, msg );
   problemUpdate.checkChangedActivities();
   ParallelColDetection<double> presolvingMethod{};
   Reductions<double> reductions{};

   REQUIRE( presolvingMethod.execute( problem, problemUpdate, num, reductions,
                                      t, cause ) == PresolveStatus::kReduced );

   const auto& reds = reductions.getReductions();
   for( const auto& tx : reductions.getTransactions() )
      problemUpdate.applyTransaction( reds.data() + tx.start, reds.data() + tx.end,
                                      ArgumentType::kPrimal );

   // all three columns are merged into one continuous column
   int remaining_col = -1;
   for( int col = 0; col < problem.getNCols(); ++col )
   {
      if( problem.getColFlags()[col].test( ColFlag::kInactive ) )
         continue;
      REQUIRE( remaining_col == -1 );
      remaining_col = col;
   }
   REQUIRE( remaining_col != -1 );
   REQUIRE( !problem.getColFlags()[remaining_col].test( ColFlag::kIntegral ) );
   REQUIRE( problem.getLowerBounds()[remaining_col] == 3 * lb );
   REQUIRE( problem.getUpperBounds()[remaining_col] == 30 );

   Vec<double> reduced( problem.getNCols(), 0.0 );
   reduced[remaining_col] = merged_value;
   Solution<double> reduced_solution{ reduced };
   Solution<double> original_solution{};
   Postsolve<double> postsolve_engine{ msg, num };
   auto status = postsolve_engine.undo( reduced_solution, original_solution, postsolve );
   REQUIRE( ( status == PostsolveStatus::kOk) );
   return original_solution.primal;
}
