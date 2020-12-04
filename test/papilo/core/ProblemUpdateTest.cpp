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

#include "catch/catch.hpp"
#include "papilo/core/Presolve.hpp"
#include "papilo/core/Problem.hpp"
#include "papilo/core/ProblemBuilder.hpp"
#include "papilo/core/Reductions.hpp"
#include "papilo/presolvers/ImplIntDetection.hpp"

papilo::Problem<double>
setupProblemPresolveSingletonRow();

papilo::Problem<double>
setupProblemPresolveSingletonRowFixed();

TEST_CASE( "happy-path-presolve-singleton-row", "[core]" )
{
   papilo::Num<double> num{};
   papilo::Problem<double> problem = setupProblemPresolveSingletonRow();
   papilo::Statistics statistics{};
   papilo::PresolveOptions presolveOptions{};
   papilo::Postsolve<double> postsolve =
       papilo::Postsolve<double>( problem, num );
   papilo::ProblemUpdate<double> problemUpdate( problem, postsolve, statistics,
                                                presolveOptions, num );
   papilo::PresolveStatus status = problemUpdate.trivialPresolve();
   //   TODO:why is the status not changed REQUIRE(status ==
   //   papilo::PresolveStatus::kReduced);
   REQUIRE( problem.getUpperBounds()[2] == 1 );
   REQUIRE( problem.getRowFlags()[1].test( papilo::RowFlag::kRedundant ) );
}

TEST_CASE( "happy-path-presolve-singleton-row-fixed", "[core]" )
{
   papilo::Num<double> num{};
   papilo::Problem<double> problem = setupProblemPresolveSingletonRowFixed();
   papilo::Statistics statistics{};
   papilo::PresolveOptions presolveOptions{};
   papilo::Postsolve<double> postsolve =
       papilo::Postsolve<double>( problem, num );
   papilo::ProblemUpdate<double> problemUpdate( problem, postsolve, statistics,
                                                presolveOptions, num );
   papilo::PresolveStatus status = problemUpdate.trivialPresolve();
   //   TODO:why is the status not changed REQUIRE(status ==
   //   papilo::PresolveStatus::kReduced);
   REQUIRE( problem.getUpperBounds()[2] == 1 );
   REQUIRE( problem.getLowerBounds()[2] == 1 );
   REQUIRE( problem.getRowFlags()[1].test( papilo::RowFlag::kRedundant ) );
   // TODO: should this variable be fixed
   REQUIRE( problem.getColFlags()[1].test( papilo::ColFlag::kFixed ) );
}

papilo::Problem<double>
setupProblemPresolveSingletonRow()
{

   const papilo::Vec<double> coefficients{ 3.0, 1.0, 1.0 };
   papilo::Vec<std::string> rowNames{ "A1", "A2" };
   papilo::Vec<std::string> columnNames{ "x", "y", "z" };
   const papilo::Vec<double> rhs{ 3.0, 1.0 };
   const papilo::Vec<double> upperBounds{ 3.0, 7.0, 7.0 };
   const papilo::Vec<double> lowerBounds{ 0.0, 0.0, 0.0 };
   papilo::Vec<uint8_t> integral = papilo::Vec<uint8_t>{ 1, 1, 1 };

   papilo::Vec<std::tuple<int, int, double>> entries{
       std::tuple<int, int, double>{ 0, 0, 2.0 },
       std::tuple<int, int, double>{ 0, 1, 1.0 },
       std::tuple<int, int, double>{ 0, 2, 1.0 },
       std::tuple<int, int, double>{ 1, 2, 1.0 } };

   papilo::ProblemBuilder<double> pb;
   pb.reserve( entries.size(), rowNames.size(), columnNames.size() );
   pb.setNumRows( rowNames.size() );
   pb.setNumCols( columnNames.size() );
   pb.setColUbAll( upperBounds );
   pb.setColLbAll( lowerBounds );
   pb.setObjAll( coefficients );
   pb.setObjOffset( 0.0 );
   pb.setColIntegralAll( integral );
   pb.setRowRhsAll( rhs );
   pb.addEntryAll( entries );
   pb.setColNameAll( columnNames );
   pb.setProblemName( "matrix for singleton row" );
   papilo::Problem<double> problem = pb.build();
   return problem;
}

papilo::Problem<double>
setupProblemPresolveSingletonRowFixed()
{

   const papilo::Vec<double> coefficients{ 3.0, 1.0, 1.0 };
   papilo::Vec<std::string> rowNames{ "A1", "A2" };
   papilo::Vec<std::string> columnNames{ "x", "y", "z" };
   const papilo::Vec<double> rhs{ 3.0, 1.0 };
   const papilo::Vec<double> upperBounds{ 3.0, 7.0, 7.0 };
   const papilo::Vec<double> lowerBounds{ 0.0, 0.0, 0.0 };
   papilo::Vec<uint8_t> integral = papilo::Vec<uint8_t>{ 1, 1, 1 };

   papilo::Vec<std::tuple<int, int, double>> entries{
       std::tuple<int, int, double>{ 0, 0, 2.0 },
       std::tuple<int, int, double>{ 0, 1, 1.0 },
       std::tuple<int, int, double>{ 0, 2, 1.0 },
       std::tuple<int, int, double>{ 1, 2, 1.0 } };

   papilo::ProblemBuilder<double> pb;
   pb.reserve( entries.size(), rowNames.size(), columnNames.size() );
   pb.setNumRows( rowNames.size() );
   pb.setNumCols( columnNames.size() );
   pb.setColUbAll( upperBounds );
   pb.setColLbAll( lowerBounds );
   pb.setObjAll( coefficients );
   pb.setObjOffset( 0.0 );
   pb.setColIntegralAll( integral );
   pb.setRowRhsAll( rhs );
   pb.addEntryAll( entries );
   pb.setColNameAll( columnNames );
   pb.setProblemName( "matrix for singleton row fixed" );
   papilo::Problem<double> problem = pb.build();
   problem.getConstraintMatrix().modifyLeftHandSide( 1, rhs[1] );
   return problem;
}
