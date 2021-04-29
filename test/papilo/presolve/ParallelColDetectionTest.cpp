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

#include "papilo/presolvers/ParallelColDetection.hpp"
#include "catch/catch.hpp"
#include "papilo/core/PresolveMethod.hpp"
#include "papilo/core/Problem.hpp"
#include "papilo/core/ProblemBuilder.hpp"

using namespace papilo;

Problem<double>
setupProblemWithParallelColumns( bool first_col_int, bool second_col_int,
                                 double factor, double ub_first_col,
                                 double ub_second_col, double lb_first_col,
                                 double lb_second_col );


TEST_CASE( "parallel_col_detection_2_integer_columns", "[presolve]" )
{
   Num<double> num{};
   Message msg{};
   Problem<double> problem =
       setupProblemWithParallelColumns( true, true, 2.0, 10.0, 10.0, 0.0, 0.0 );
   Statistics statistics{};
   PresolveOptions presolveOptions{};
   Postsolve<double> postsolve = Postsolve<double>( problem, num );
   ProblemUpdate<double> problemUpdate( problem, postsolve, statistics,
                                        presolveOptions, num, msg );
   problemUpdate.checkChangedActivities();
   ParallelColDetection<double> presolvingMethod{};
   Reductions<double> reductions{};

   PresolveStatus presolveStatus =
       presolvingMethod.execute( problem, problemUpdate, num, reductions );

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

TEST_CASE( "parallel_col_detection_2_integer_columns-hole", "[presolve]" )
{
   Num<double> num{};
   Message msg{};
   Problem<double> problem =
       setupProblemWithParallelColumns( true, true, 3.0, 1.0, 2.0, 0.0, 1.0 );
   Statistics statistics{};
   PresolveOptions presolveOptions{};
   Postsolve<double> postsolve = Postsolve<double>( problem, num );
   ProblemUpdate<double> problemUpdate( problem, postsolve, statistics,
                                        presolveOptions, num, msg );
   problemUpdate.checkChangedActivities();
   ParallelColDetection<double> presolvingMethod{};
   Reductions<double> reductions{};

   PresolveStatus presolveStatus =
       presolvingMethod.execute( problem, problemUpdate, num, reductions );

   REQUIRE( presolveStatus == PresolveStatus::kUnchanged );
}

TEST_CASE( "parallel_col_detection_2_continuous_columns", "[presolve]" )
{
   Num<double> num{};
   Message msg{};
   Problem<double> problem = setupProblemWithParallelColumns(
       false, false, 2.0, 10.0, 10.0, 0.0, 0.0 );
   Statistics statistics{};
   PresolveOptions presolveOptions{};
   Postsolve<double> postsolve = Postsolve<double>( problem, num );
   ProblemUpdate<double> problemUpdate( problem, postsolve, statistics,
                                        presolveOptions, num, msg );
   problemUpdate.checkChangedActivities();
   ParallelColDetection<double> presolvingMethod{};
   Reductions<double> reductions{};

   PresolveStatus presolveStatus =
       presolvingMethod.execute( problem, problemUpdate, num, reductions );

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
   Num<double> num{};
   Message msg{};
   Problem<double> problem = setupProblemWithParallelColumns(
       true, false, 2.0, 10.0, 10.0, 0.0, 0.0 );
   Statistics statistics{};
   PresolveOptions presolveOptions{};
   Postsolve<double> postsolve = Postsolve<double>( problem, num );
   ProblemUpdate<double> problemUpdate( problem, postsolve, statistics,
                                        presolveOptions, num, msg );
   problemUpdate.checkChangedActivities();
   ParallelColDetection<double> presolvingMethod{};
   Reductions<double> reductions{};

   PresolveStatus presolveStatus =
       presolvingMethod.execute( problem, problemUpdate, num, reductions );

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
   Num<double> num{};
   Message msg{};
   Problem<double> problem = setupProblemWithParallelColumns(
       false, true, 2.0, 10.0, 10.0, 0.0, 0.0 );
   Statistics statistics{};
   PresolveOptions presolveOptions{};
   Postsolve<double> postsolve = Postsolve<double>( problem, num );
   ProblemUpdate<double> problemUpdate( problem, postsolve, statistics,
                                        presolveOptions, num, msg );
   problemUpdate.checkChangedActivities();
   ParallelColDetection<double> presolvingMethod{};
   Reductions<double> reductions{};

   PresolveStatus presolveStatus =
       presolvingMethod.execute( problem, problemUpdate, num, reductions );

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
   Num<double> num{};
   Message msg{};
   Problem<double> problem =
       setupProblemWithParallelColumns( false, true, 1.0, 0.9, 10.0, 0.0, 0.0 );
   Statistics statistics{};
   PresolveOptions presolveOptions{};
   Postsolve<double> postsolve = Postsolve<double>( problem, num );
   ProblemUpdate<double> problemUpdate( problem, postsolve, statistics,
                                        presolveOptions, num, msg );
   problemUpdate.checkChangedActivities();
   ParallelColDetection<double> presolvingMethod{};
   Reductions<double> reductions{};

   PresolveStatus presolveStatus =
       presolvingMethod.execute( problem, problemUpdate, num, reductions );

   REQUIRE( presolveStatus == PresolveStatus::kUnchanged );
}

TEST_CASE( "parallel_col_detection_int_cont_merge_failed", "[presolve]" )
{
   Num<double> num{};
   Message msg{};
   Problem<double> problem =
       setupProblemWithParallelColumns( true, false, 1.0, 10.0, 0.9, 0.0, 0.0 );
   Statistics statistics{};
   PresolveOptions presolveOptions{};
   Postsolve<double> postsolve = Postsolve<double>( problem, num );
   ProblemUpdate<double> problemUpdate( problem, postsolve, statistics,
                                        presolveOptions, num, msg );
   problemUpdate.checkChangedActivities();
   ParallelColDetection<double> presolvingMethod{};
   Reductions<double> reductions{};

   PresolveStatus presolveStatus =
       presolvingMethod.execute( problem, problemUpdate, num, reductions );

   REQUIRE( presolveStatus == PresolveStatus::kUnchanged );
}

TEST_CASE( "parallel_col_detection_obj_not_parallel", "[presolve]" )
{
   Num<double> num{};
   Message msg{};
   Vec<double> obj = {3,2};
   Problem<double> problem =
       setupProblemWithParallelColumns( true, true, 1.0, 10.0, 10.0, 0.0, 0.0 );
   problem.setObjective(obj, 0);
   Statistics statistics{};
   PresolveOptions presolveOptions{};
   Postsolve<double> postsolve = Postsolve<double>( problem, num );
   ProblemUpdate<double> problemUpdate( problem, postsolve, statistics,
                                        presolveOptions, num, msg );
   problemUpdate.checkChangedActivities();
   ParallelColDetection<double> presolvingMethod{};
   Reductions<double> reductions{};

   PresolveStatus presolveStatus =
       presolvingMethod.execute( problem, problemUpdate, num, reductions );

   REQUIRE( presolveStatus == PresolveStatus::kUnchanged );
}

Problem<double>
setupProblemWithParallelColumns( bool first_col_int, bool second_col_int,
                                 double factor, double ub_first_col,
                                 double ub_second_col, double lb_first_col,
                                 double lb_second_col )
{
   Vec<double> coefficients{ 1.0, 1.0 * factor };
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
       std::tuple<int, int, double>{ 1, 0, 1.0 * factor },
       std::tuple<int, int, double>{ 0, 1, 2.0 },
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
