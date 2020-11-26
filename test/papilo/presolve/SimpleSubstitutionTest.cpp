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

#include "papilo/presolvers/SimpleSubstitution.hpp"
#include "catch/catch.hpp"
#include "papilo/core/PresolveMethod.hpp"
#include "papilo/core/Problem.hpp"
#include "papilo/core/ProblemBuilder.hpp"

using namespace papilo;

Problem<double>
setupProblemWithSimpleSubstitution( uint8_t is_x_integer, uint8_t is_y_integer,
                                    double a_y );

TEST_CASE( "happy-path-simple-substitution-for-2-int", "[presolve]" )
{
   Num<double> num{};
   Problem<double> problem = setupProblemWithSimpleSubstitution( 1, 1, 1.0 );
   Statistics statistics{};
   PresolveOptions presolveOptions{};
   presolveOptions.dualreds = 0;
   Postsolve<double> postsolve = Postsolve<double>( problem, num );
   ProblemUpdate<double> problemUpdate( problem, postsolve, statistics,
                                        presolveOptions, num );
   SimpleSubstitution<double> presolvingMethod{};
   Reductions<double> reductions{};
   problem.recomputeAllActivities();

   PresolveStatus presolveStatus =
       presolvingMethod.execute( problem, problemUpdate, num, reductions );

   // Reduction => x = 2 - y/2 -> 0.5 (1 for int) <= x <= 2
   BOOST_ASSERT( presolveStatus == PresolveStatus::kReduced );
   BOOST_ASSERT( reductions.size() == 5 );

   BOOST_ASSERT( reductions.getReduction( 0 ).col == RowReduction::LOCKED );
   BOOST_ASSERT( reductions.getReduction( 0 ).row == 0 );
   BOOST_ASSERT( reductions.getReduction( 0 ).newval == 0 );

   BOOST_ASSERT( reductions.getReduction( 1 ).row ==
                 ColReduction::BOUNDS_LOCKED );
   BOOST_ASSERT( reductions.getReduction( 1 ).col == 1 );
   BOOST_ASSERT( reductions.getReduction( 1 ).newval == 0 );

   BOOST_ASSERT( reductions.getReduction( 2 ).col == 0 );
   BOOST_ASSERT( reductions.getReduction( 2 ).row ==
                 papilo::ColReduction::UPPER_BOUND );
   BOOST_ASSERT( reductions.getReduction( 2 ).newval == 2 );

   BOOST_ASSERT( reductions.getReduction( 3 ).col == 0 );
   BOOST_ASSERT( reductions.getReduction( 3 ).row ==
                 papilo::ColReduction::LOWER_BOUND );
   // TODO: this could be rounded down for an integer
   BOOST_ASSERT( reductions.getReduction( 3 ).newval == 0.5 );

   BOOST_ASSERT( reductions.getReduction( 4 ).col == 1 );
   BOOST_ASSERT( reductions.getReduction( 4 ).row ==
                 papilo::ColReduction::SUBSTITUTE );
   BOOST_ASSERT( reductions.getReduction( 4 ).newval == 0 );
}

TEST_CASE( "happy-path-simple-substitution-for-2-continuous", "[presolve]" )
{
   Num<double> num{};
   Problem<double> problem = setupProblemWithSimpleSubstitution( 0, 0, 1.0 );
   Statistics statistics{};
   PresolveOptions presolveOptions{};
   presolveOptions.dualreds = 0;
   Postsolve<double> postsolve = Postsolve<double>( problem, num );
   ProblemUpdate<double> problemUpdate( problem, postsolve, statistics,
                                        presolveOptions, num );
   SimpleSubstitution<double> presolvingMethod{};
   Reductions<double> reductions{};
   problem.recomputeAllActivities();

   PresolveStatus presolveStatus =
       presolvingMethod.execute( problem, problemUpdate, num, reductions );

   // Reduction => x = 4 - 2y -> 0 <= x <= 4 (no further bound relaxation)
   BOOST_ASSERT( presolveStatus == PresolveStatus::kReduced );
   BOOST_ASSERT( reductions.size() == 3 );

   BOOST_ASSERT( reductions.getReduction( 0 ).col == RowReduction::LOCKED );
   BOOST_ASSERT( reductions.getReduction( 0 ).row == 0 );
   BOOST_ASSERT( reductions.getReduction( 0 ).newval == 0 );

   BOOST_ASSERT( reductions.getReduction( 1 ).row ==
                 ColReduction::BOUNDS_LOCKED );
   BOOST_ASSERT( reductions.getReduction( 1 ).col == 0 );
   BOOST_ASSERT( reductions.getReduction( 1 ).newval == 0 );

   BOOST_ASSERT( reductions.getReduction( 2 ).col == 0 );
   BOOST_ASSERT( reductions.getReduction( 2 ).row ==
                 papilo::ColReduction::SUBSTITUTE );
   BOOST_ASSERT( reductions.getReduction( 2 ).newval == 0 );
}

TEST_CASE( "happy-path-simple-substitution-for-continuous-and-integer",
           "[presolve]" )
{
   Num<double> num{};
   Problem<double> problem = setupProblemWithSimpleSubstitution( 0, 1, 1.0 );
   Statistics statistics{};
   PresolveOptions presolveOptions{};
   presolveOptions.dualreds = 0;
   Postsolve<double> postsolve = Postsolve<double>( problem, num );
   ProblemUpdate<double> problemUpdate( problem, postsolve, statistics,
                                        presolveOptions, num );
   SimpleSubstitution<double> presolvingMethod{};
   Reductions<double> reductions{};
   problem.recomputeAllActivities();

   PresolveStatus presolveStatus =
       presolvingMethod.execute( problem, problemUpdate, num, reductions );

   // Reduction => x = 4 - 2y -> 0 <= x <= 4 (no further bound relaxation)
   BOOST_ASSERT( presolveStatus == PresolveStatus::kReduced );
   BOOST_ASSERT( reductions.size() == 3 );

   BOOST_ASSERT( reductions.getReduction( 0 ).col == RowReduction::LOCKED );
   BOOST_ASSERT( reductions.getReduction( 0 ).row == 0 );
   BOOST_ASSERT( reductions.getReduction( 0 ).newval == 0 );

   BOOST_ASSERT( reductions.getReduction( 1 ).row ==
                 ColReduction::BOUNDS_LOCKED );
   BOOST_ASSERT( reductions.getReduction( 1 ).col == 0 );
   BOOST_ASSERT( reductions.getReduction( 1 ).newval == 0 );

   BOOST_ASSERT( reductions.getReduction( 2 ).col == 0 );
   BOOST_ASSERT( reductions.getReduction( 2 ).row ==
                 papilo::ColReduction::SUBSTITUTE );
   BOOST_ASSERT( reductions.getReduction( 2 ).newval == 0 );
}

TEST_CASE( "failed-path-simple-substitution-for-2-int", "[presolve]" )
{
   Num<double> num{};
   Problem<double> problem = setupProblemWithSimpleSubstitution( 1, 1, 3.0 );
   Statistics statistics{};
   PresolveOptions presolveOptions{};
   presolveOptions.dualreds = 0;
   Postsolve<double> postsolve = Postsolve<double>( problem, num );
   ProblemUpdate<double> problemUpdate( problem, postsolve, statistics,
                                        presolveOptions, num );
   SimpleSubstitution<double> presolvingMethod{};
   Reductions<double> reductions{};
   problem.recomputeAllActivities();

   PresolveStatus presolveStatus =
       presolvingMethod.execute( problem, problemUpdate, num, reductions );

   BOOST_ASSERT( presolveStatus == PresolveStatus::kUnchanged );
}

Problem<double>
setupProblemWithSimpleSubstitution( uint8_t is_x_integer, uint8_t is_y_integer,
                                    double a_y )
{
   // 2x + y = 4
   // 0<= x,y y= 3
   Vec<double> coefficients{ 3.0, 1.0 };
   Vec<double> upperBounds{ 3.0, 3.0 };
   Vec<double> lowerBounds{ 0.0, 0.0 };
   Vec<uint8_t> isIntegral{ is_x_integer, is_y_integer };

   Vec<double> rhs{ 4.0 };
   Vec<std::string> rowNames{ "A1" };
   Vec<std::string> columnNames{ "c1", "c2" };
   Vec<std::tuple<int, int, double>> entries{
       std::tuple<int, int, double>{ 0, 0, 2.0 },
       std::tuple<int, int, double>{ 0, 1, a_y },
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
   pb.setProblemName( "matrix for testing simple probing" );
   Problem<double> problem = pb.build();
   problem.getConstraintMatrix().modifyLeftHandSide( 0, rhs[0] );
   return problem;
}
