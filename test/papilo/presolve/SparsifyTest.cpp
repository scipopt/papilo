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

#include "papilo/presolvers/Sparsify.hpp"
#include "catch/catch.hpp"
#include "papilo/core/PresolveMethod.hpp"
#include "papilo/core/Problem.hpp"
#include "papilo/core/ProblemBuilder.hpp"

using namespace papilo;

Problem<double>
setupProblemWithSparsify();

Problem<double>
setupProblemWithSparsifyMoreThanOneColumn();

TEST_CASE( "happy-path-sparsify", "[presolve]" )
{
   Num<double> num{};
   Message msg{};
   Problem<double> problem = setupProblemWithSparsify( );
   Statistics statistics{};
   PresolveOptions presolveOptions{};
   presolveOptions.dualreds = 0;
   Postsolve<double> postsolve = Postsolve<double>( problem, num );
   ProblemUpdate<double> problemUpdate( problem, postsolve, statistics,
                                        presolveOptions, num, msg );
   Sparsify<double> presolvingMethod{};
   Reductions<double> reductions{};
   problem.recomputeAllActivities();

   PresolveStatus presolveStatus =
       presolvingMethod.execute( problem, problemUpdate, num, reductions );
   REQUIRE( presolveStatus == PresolveStatus::kReduced );
   REQUIRE( reductions.size() == 3 );
   REQUIRE( reductions.getReduction( 0 ).row == 0 );
   REQUIRE( reductions.getReduction( 0 ).col == RowReduction::LOCKED );
   REQUIRE( reductions.getReduction( 0 ).newval == 0 );

   REQUIRE( reductions.getReduction( 1 ).row == 0 );
   REQUIRE( reductions.getReduction( 1 ).col == RowReduction::SPARSIFY );
   REQUIRE( reductions.getReduction( 1 ).newval == 1 );

   REQUIRE( reductions.getReduction( 2 ).row == 1 );
   REQUIRE( reductions.getReduction( 2 ).col == RowReduction::NONE );
   REQUIRE( reductions.getReduction( 2 ).newval == -1 );
}

TEST_CASE( "failed-path-sparsify", "[presolve]" )
{
   Num<double> num{};
   Message msg{};

   Problem<double> problem = setupProblemWithSparsifyMoreThanOneColumn( );
   Statistics statistics{};
   PresolveOptions presolveOptions{};
   presolveOptions.dualreds = 0;
   Postsolve<double> postsolve = Postsolve<double>( problem, num );
   ProblemUpdate<double> problemUpdate( problem, postsolve, statistics,
                                        presolveOptions, num, msg );
   Sparsify<double> presolvingMethod{};
   Reductions<double> reductions{};
   problem.recomputeAllActivities();

   PresolveStatus presolveStatus =
       presolvingMethod.execute( problem, problemUpdate, num, reductions );
   // TODO: REQUIRE( presolveStatus == PresolveStatus::kUnchanged );
}

Problem<double>
setupProblemWithSparsify()
{
   Num<double> num{};
   Vec<double> coefficients{ 3.0, 1.0, 1.0 };
   Vec<double> upperBounds{ 3.0, 3.0, 3.0 };
   Vec<double> lowerBounds{ 0.0, 0.0, 0.0 };

   Vec<double> rhs{ 4.0, 2.0 };
   Vec<std::string> rowNames{ "r1", "r2" };
   Vec<std::string> columnNames{ "c1", "c2", "c3" };
   Vec<std::tuple<int, int, double>> entries{
       std::tuple<int, int, double>{ 0, 0, 1.0 },
       std::tuple<int, int, double>{ 0, 1, 1.0 },
       std::tuple<int, int, double>{ 1, 0, 1.0 },
       std::tuple<int, int, double>{ 1, 1, 1.0 },
       std::tuple<int, int, double>{ 1, 2, 1.0 },
   };

   ProblemBuilder<double> pb;
   pb.reserve( entries.size(), rowNames.size(), columnNames.size() );
   pb.setNumRows( rowNames.size() );
   pb.setNumCols( columnNames.size() );
   pb.setColUbAll( upperBounds );
   pb.setColLbAll( lowerBounds );
   pb.setObjAll( coefficients );
   pb.setObjOffset( 0.0 );
   pb.setRowRhsAll( rhs );
   pb.addEntryAll( entries );
   pb.setColNameAll( columnNames );
   pb.setProblemName( "matrix for testing sparsify" );
   Problem<double> problem = pb.build();
   problem.getConstraintMatrix().modifyLeftHandSide( 0,num, rhs[0] );
   return problem;
}

Problem<double>
setupProblemWithSparsifyMoreThanOneColumn()
{
   // 2x + y = 4
   // 0<= x,y y= 3
   Num<double> num{};
   Vec<double> coefficients{ 3.0, 1.0, 1.0, 1.0 };
   Vec<double> upperBounds{ 3.0, 3.0, 3.0, 3.0 };
   Vec<double> lowerBounds{ 0.0, 0.0, 0.0, 0.0 };

   Vec<double> rhs{ 4.0, 2.0 };
   Vec<std::string> rowNames{ "r1", "r2" };
   Vec<std::string> columnNames{ "c1", "c2", "c3", "c4" };
   Vec<std::tuple<int, int, double>> entries{
       std::tuple<int, int, double>{ 0, 0, 1.0 },
       std::tuple<int, int, double>{ 0, 1, 1.0 },
       std::tuple<int, int, double>{ 1, 0, 1.0 },
       std::tuple<int, int, double>{ 1, 1, 1.0 },
       std::tuple<int, int, double>{ 1, 2, 1.0 },
       std::tuple<int, int, double>{ 1, 3, 1.0 },
   };

   ProblemBuilder<double> pb;
   pb.reserve( entries.size(), rowNames.size(), columnNames.size() );
   pb.setNumRows( rowNames.size() );
   pb.setNumCols( columnNames.size() );
   pb.setColUbAll( upperBounds );
   pb.setColLbAll( lowerBounds );
   pb.setObjAll( coefficients );
   pb.setObjOffset( 0.0 );
   pb.setRowRhsAll( rhs );
   pb.addEntryAll( entries );
   pb.setColNameAll( columnNames );
   pb.setProblemName(
       "matrix for testing sparsify -> will fail because l1 !=l2 +1" );
   Problem<double> problem = pb.build();
   problem.getConstraintMatrix().modifyLeftHandSide( 0, num, rhs[0] );
   return problem;
}
