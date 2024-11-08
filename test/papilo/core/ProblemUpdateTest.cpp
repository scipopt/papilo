/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*               This file is part of the program and library                */
/*    PaPILO --- Parallel Presolve for Integer and Linear Optimization       */
/*                                                                           */
/* Copyright (C) 2020-2024 Zuse Institute Berlin (ZIB)                       */
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

#include "papilo/external/catch/catch_amalgamated.hpp"
#include "papilo/core/Problem.hpp"
#include "papilo/core/ProblemBuilder.hpp"
#include "papilo/core/Reductions.hpp"
#include "papilo/external/catch/catch.hpp"
#include "papilo/presolvers/ImplIntDetection.hpp"

namespace papilo
{
Problem<double>
setupProblemPresolveSingletonRow();

Problem<double>
setupProblemPresolveSingletonRowFixed();

Problem<double>
setupProblemWIthCliques();

TEST_CASE( "trivial-presolve-singleton-row", "[core]" )
{
   Num<double> num{};
   Message msg{};
   Problem<double> problem = setupProblemPresolveSingletonRow();
   Statistics statistics{};
   PresolveOptions presolveOptions{};
   PostsolveStorage<double> postsolve =
       PostsolveStorage<double>( problem, num, presolveOptions );
   ProblemUpdate<double> problemUpdate( problem, postsolve, statistics,
                                        presolveOptions, num, msg );
   problemUpdate.trivialPresolve();
   REQUIRE( problem.getUpperBounds()[2] == 1 );
   REQUIRE( problem.getRowFlags()[1].test( RowFlag::kRedundant ) );
}

TEST_CASE( "trivial-presolve-singleton-row-pt-2", "[core]" )
{
   Num<double> num{};
   Problem<double> problem = setupProblemPresolveSingletonRowFixed();
   Message msg{};
   Statistics statistics{};
   PresolveOptions presolveOptions{};
   PostsolveStorage<double> postsolve =
       PostsolveStorage<double>( problem, num, presolveOptions );
   ProblemUpdate<double> problemUpdate( problem, postsolve, statistics,
                                        presolveOptions, num, msg );
   problemUpdate.trivialPresolve();
   REQUIRE( problem.getUpperBounds()[2] == 1 );
   REQUIRE( problem.getLowerBounds()[2] == 1 );
   REQUIRE( problem.getRowFlags()[1].test( RowFlag::kRedundant ) );
   REQUIRE( problemUpdate.getSingletonCols().size() == 2 );
}

TEST_CASE( "clique-row-flag-detection1", "[core]" )
{
   Num<double> num{};
   Problem<double> problem = setupProblemWIthCliques();
   REQUIRE( !problem.getRowFlags()[0].test( RowFlag::kClique ) );
   REQUIRE( problem.getRowFlags()[0].test( RowFlag::kSOS1 ) );
}

TEST_CASE( "clique-row-flag-detection2", "[core]" )
{
   Num<double> num{};
   Problem<double> problem = setupProblemWIthCliques();
   REQUIRE( !problem.getRowFlags()[1].test( RowFlag::kClique ) );
   REQUIRE( !problem.getRowFlags()[1].test( RowFlag::kSOS1 ) );
}

TEST_CASE( "clique-row-flag-detection3", "[core]" )
{
   Num<double> num{};
   Problem<double> problem = setupProblemWIthCliques();
   REQUIRE( !problem.getRowFlags()[2].test( RowFlag::kClique ) );
   REQUIRE( problem.getRowFlags()[2].test( RowFlag::kSOS1 ) );
}

TEST_CASE( "clique-row-flag-detection4", "[core]" )
{
   Num<double> num{};
   Problem<double> problem = setupProblemWIthCliques();
   REQUIRE( !problem.getRowFlags()[3].test( RowFlag::kClique ) );
   REQUIRE( problem.getRowFlags()[3].test( RowFlag::kSOS1 ) );
}

TEST_CASE( "clique-row-flag-detection5", "[core]" )
{
   Num<double> num{};
   Problem<double> problem = setupProblemWIthCliques();
   REQUIRE( !problem.getRowFlags()[4].test( RowFlag::kClique ) );
   REQUIRE( !problem.getRowFlags()[4].test( RowFlag::kSOS1 ) );
}

TEST_CASE( "clique-row-flag-detection6", "[core]" )
{
   Num<double> num{};
   Problem<double> problem = setupProblemWIthCliques();
   REQUIRE( !problem.getRowFlags()[5].test( RowFlag::kClique ) );
   REQUIRE( !problem.getRowFlags()[5].test( RowFlag::kSOS1 ) );
}

TEST_CASE( "clique-row-flag-detection7", "[core]" )
{
   Num<double> num{};
   Problem<double> problem = setupProblemWIthCliques();
   REQUIRE( !problem.getRowFlags()[6].test( RowFlag::kClique ) );
   REQUIRE( problem.getRowFlags()[6].test( RowFlag::kSOS1 ) );
}

TEST_CASE( "clique-row-flag-detection8", "[core]" )
{
   Num<double> num{};
   Problem<double> problem = setupProblemWIthCliques();
   REQUIRE( !problem.getRowFlags()[7].test( RowFlag::kClique ) );
   REQUIRE( !problem.getRowFlags()[7].test( RowFlag::kSOS1 ) );
}

TEST_CASE( "clique-row-flag-detection9", "[core]" )
{
   Num<double> num{};
   Problem<double> problem = setupProblemWIthCliques();
   // REQUIRE( problem.getRowFlags()[8].test( RowFlag::kClique ) );
   // REQUIRE( !problem.getRowFlags()[8].test( RowFlag::kSOS1 ) );
}

Problem<double>
setupProblemPresolveSingletonRow()
{

   const Vec<double> coefficients{ 3.0, 1.0, 1.0 };
   Vec<std::string> rowNames{ "A1", "A2" };
   Vec<std::string> columnNames{ "x", "y", "z" };
   const Vec<double> rhs{ 3.0, 1.0 };
   const Vec<double> upperBounds{ 3.0, 7.0, 7.0 };
   const Vec<double> lowerBounds{ 0.0, 0.0, 0.0 };
   Vec<uint8_t> integral = Vec<uint8_t>{ 1, 1, 1 };

   Vec<std::tuple<int, int, double>> entries{
       std::tuple<int, int, double>{ 0, 0, 2.0 },
       std::tuple<int, int, double>{ 0, 1, 1.0 },
       std::tuple<int, int, double>{ 0, 2, 1.0 },
       std::tuple<int, int, double>{ 1, 2, 1.0 } };

   ProblemBuilder<double> pb;
   pb.reserve( (int)entries.size(), (int)rowNames.size(),
               (int)columnNames.size() );
   pb.setNumRows( (int)rowNames.size() );
   pb.setNumCols( (int)columnNames.size() );
   pb.setColUbAll( upperBounds );
   pb.setColLbAll( lowerBounds );
   pb.setObjAll( coefficients );
   pb.setObjOffset( 0.0 );
   pb.setColIntegralAll( integral );
   pb.setRowRhsAll( rhs );
   pb.addEntryAll( entries );
   pb.setColNameAll( columnNames );
   pb.setProblemName( "matrix for singleton row" );
   Problem<double> problem = pb.build();
   return problem;
}

Problem<double>
setupProblemPresolveSingletonRowFixed()
{
   Num<double> num{};
   const Vec<double> coefficients{ 3.0, 1.0, 1.0 };
   Vec<std::string> rowNames{ "A1", "A2" };
   Vec<std::string> columnNames{ "x", "y", "z" };
   const Vec<double> rhs{ 3.0, 1.0 };
   const Vec<double> upperBounds{ 3.0, 7.0, 7.0 };
   const Vec<double> lowerBounds{ 0.0, 0.0, 0.0 };
   Vec<uint8_t> integral = Vec<uint8_t>{ 1, 1, 1 };

   Vec<std::tuple<int, int, double>> entries{
       std::tuple<int, int, double>{ 0, 0, 2.0 },
       std::tuple<int, int, double>{ 0, 1, 1.0 },
       std::tuple<int, int, double>{ 0, 2, 1.0 },
       std::tuple<int, int, double>{ 1, 2, 1.0 } };

   ProblemBuilder<double> pb;
   pb.reserve( (int)entries.size(), (int)rowNames.size(),
               (int)columnNames.size() );
   pb.setNumRows( (int)rowNames.size() );
   pb.setNumCols( (int)columnNames.size() );
   pb.setColUbAll( upperBounds );
   pb.setColLbAll( lowerBounds );
   pb.setObjAll( coefficients );
   pb.setObjOffset( 0.0 );
   pb.setColIntegralAll( integral );
   pb.setRowRhsAll( rhs );
   pb.addEntryAll( entries );
   pb.setColNameAll( columnNames );
   pb.setProblemName( "matrix for singleton row fixed" );
   Problem<double> problem = pb.build();
   problem.getConstraintMatrix().modifyLeftHandSide( 1, num, rhs[1] );
   return problem;
}

Problem<double>
setupProblemWIthCliques()
{
   Num<double> num{};
   const Vec<double> coefficients{ 1.0, 1.0, -1.0, 1.0, 1.0, -1.0, 1.0 };
   Vec<std::string> rowNames{ "A1", "A2", "A3", "A4", "A5",
                              "A6", "A7", "A8", "A9" };
   Vec<std::string> columnNames{ "x", "y", "z", "a", "b", "c", "d" };
   const Vec<double> rhs{ 1.0, 2.0, 2.0, 2.0, 1.0, 1.0, 1.0, 1.0, 1.0 };
   const Vec<double> lhs{ -10.0, -10.0, -10.0, -10.0, -10.0,
                          -10.0, -10.0, -10.0, -10.0 };
   const Vec<double> upperBounds{ 1.0, 1.0, 0.0, 5.0, 1.0, 0.0, 1.0 };
   const Vec<double> lowerBounds{ 0.0, 0.0, -1.0, 0.0, 0.0, -5.0, 0.0 };
   Vec<uint8_t> integral = Vec<uint8_t>{ 1, 1, 1, 1, 0, 1, 1 };

   Vec<std::tuple<int, int, double>> entries{
       std::tuple<int, int, double>{ 0, 0, 1.0 },
       std::tuple<int, int, double>{ 0, 1, 1.0 },
       std::tuple<int, int, double>{ 0, 2, -1.0 },
       // First Row is a SOS1, but not a Clique
       std::tuple<int, int, double>{ 1, 0, 1.0 },
       std::tuple<int, int, double>{ 1, 1, 1.0 },
       std::tuple<int, int, double>{ 1, 2, -1.0 },
       // Second Row is not a Clique, not SOS1
       std::tuple<int, int, double>{ 2, 0, 1.0 },
       std::tuple<int, int, double>{ 2, 1, 2.0 },
       std::tuple<int, int, double>{ 2, 2, -1.5 },
       // Third Row is a SOS1, but no Clique
       std::tuple<int, int, double>{ 3, 0, 1.5 },
       std::tuple<int, int, double>{ 3, 1, 2.0 },
       std::tuple<int, int, double>{ 3, 3, 1.0 },
       // Fourth row is a SOS1 but no clique
       std::tuple<int, int, double>{ 4, 0, 2.0 },
       std::tuple<int, int, double>{ 4, 1, 1.0 },
       std::tuple<int, int, double>{ 4, 2, -1.0 },
       // Fifth Row is not a Clique, no SOS1
       std::tuple<int, int, double>{ 5, 0, 1.0 },
       std::tuple<int, int, double>{ 5, 1, 0.5 },
       std::tuple<int, int, double>{ 5, 4, -1.0 },
       // Sixth Row is not a Clique, no SOS1
       std::tuple<int, int, double>{ 6, 0, -5.0 },
       std::tuple<int, int, double>{ 6, 2, 6.0 },
       std::tuple<int, int, double>{ 6, 5, 5.5 },
       // Seventh Row is not a Clique, but a SOS1
       std::tuple<int, int, double>{ 7, 0, -5.0 },
       std::tuple<int, int, double>{ 7, 2, 6.0 },
       std::tuple<int, int, double>{ 7, 5, 5.0 },
       // Eigth Row is not a Clique, no SOS1
       std::tuple<int, int, double>{ 7, 0, 1.0 },
       std::tuple<int, int, double>{ 7, 1, 1.0 },
       std::tuple<int, int, double>{ 7, 6, 1.0 },
       // Ninth Row is a clique, no SOS1
   };

   ProblemBuilder<double> pb;
   pb.reserve( (int)entries.size(), (int)rowNames.size(),
               (int)columnNames.size() );
   pb.setNumRows( (int)rowNames.size() );
   pb.setNumCols( (int)columnNames.size() );
   pb.setColUbAll( upperBounds );
   pb.setColLbAll( lowerBounds );
   pb.setObjAll( coefficients );
   pb.setObjOffset( 0.0 );
   pb.setColIntegralAll( integral );
   pb.setRowRhsAll( rhs );
   pb.setRowLhsAll( lhs );
   pb.addEntryAll( entries );
   pb.setColNameAll( columnNames );
   pb.setProblemName( "matrix for cliques" );
   Problem<double> problem = pb.build();
   return problem;
}

} // namespace papilo
