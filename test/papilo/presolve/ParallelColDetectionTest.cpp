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

Vec<double>
postsolveParallelContIntColumns( double lb_cont_col, double ub_cont_col,
                                 double lb_int_col, double ub_int_col,
                                 double scale, double merged_value );

Vec<double>
postsolveParallelColumns( bool integral, bool col1_lb_inf, double lb_col1,
                          bool col1_ub_inf, double ub_col1, double lb_col2,
                          double ub_col2, double scale, double merged_value );

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
   REQUIRE( (postsolve.undo( reduced_solution, original_solution, storage ) ==
            PostsolveStatus::kOk) );
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
          postsolveParallelColumns( true, false, i.lb1, false, i.ub1, i.lb2,
                                    i.ub2, i.scale, i.merged_value );
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
   Vec<double> values = postsolveParallelColumns( false, true, 0, true, 0, 0,
                                                  10, 1, 5 );
   REQUIRE( values[0] == Catch::Approx( 5 ) );
   REQUIRE( values[1] == 0 );
}

Vec<double>
postsolveParallelColumns( bool integral, bool col1_lb_inf, double lb_col1,
                          bool col1_ub_inf, double ub_col1, double lb_col2,
                          double ub_col2, double scale, double merged_value )
{
   ProblemBuilder<double> pb;
   pb.reserve( 0, 0, 2 );
   pb.setNumRows( 0 );
   pb.setNumCols( 2 );
   pb.setColLbAll( { lb_col1, lb_col2 } );
   pb.setColUbAll( { ub_col1, ub_col2 } );
   pb.setColLbInf( 0, col1_lb_inf );
   pb.setColUbInf( 0, col1_ub_inf );
   pb.setObjAll( { 0.0, 0.0 } );
   pb.setColIntegralAll( { integral, integral } );
   pb.setColNameAll( Vec<std::string>{ "x1", "x2" } );
   Problem<double> problem = pb.build();

   Num<double> num{};
   Message msg{};
   PostsolveStorage<double> storage{ problem, num, PresolveOptions{} };
   storage.storeParallelCols( 0, integral, col1_lb_inf, lb_col1, col1_ub_inf,
                              ub_col1, 1, integral, false, lb_col2, false,
                              ub_col2, scale );

   // the merged column takes over the index of the second column
   Solution<double> reduced_solution{ Vec<double>{ 0, merged_value } };
   Solution<double> original_solution{};
   Postsolve<double> postsolve{ msg, num };
   REQUIRE( (postsolve.undo( reduced_solution, original_solution, storage ) ==
            PostsolveStatus::kOk) );
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
