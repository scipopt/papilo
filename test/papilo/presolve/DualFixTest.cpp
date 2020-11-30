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

#include "papilo/presolvers/DualFix.hpp"
#include "catch/catch.hpp"
#include "papilo/core/PresolveMethod.hpp"
#include "papilo/core/Problem.hpp"
#include "papilo/core/ProblemBuilder.hpp"

using namespace papilo;

Problem<double>
setupExample4ofChapter4Dot4InPresolveReductions();

Problem<double>
setupExample5ofChapter4Dot4InPresolveReductions();

Problem<double>
setupMatrixForDualFixFirstColumnOnlyPositiveValues();

Problem<double>
setupMatrixForDualSubstitution();

Problem<double>
setupMatrixForDualSubstitutionEquation();

Problem<double>
setupMatrixForDualSubstitutionWithUnboundedVar();

Problem<double>
setupMatrixForDualSubstitutionIntegerRounding();

TEST_CASE( "trivial-column-presolve-does-dual-presolve-already", "[presolve]" )
{
   Num<double> num{};
   Problem<double> problem =
       setupMatrixForDualFixFirstColumnOnlyPositiveValues();
   Statistics statistics{};
   PresolveOptions presolveOptions{};
   Postsolve<double> postsolve = Postsolve<double>( problem, num );
   ProblemUpdate<double> problemUpdate( problem, postsolve, statistics,
                                        presolveOptions, num );

   //   problem.recomputeAllActivities(); // TODO: what does trivial presolve
   //   and if recalculates the dual fix why is it needed
   problemUpdate.trivialPresolve();

   REQUIRE( problem.getColFlags()[0].test( ColFlag::kFixed ) );
   REQUIRE( problem.getLowerBounds()[0] == 1 );
   REQUIRE( problem.getUpperBounds()[0] == 1 );
}

TEST_CASE( "happy-path-dual-fix", "[presolve]" )
{
   Num<double> num{};
   Problem<double> problem =
       setupMatrixForDualFixFirstColumnOnlyPositiveValues();
   Statistics statistics{};
   PresolveOptions presolveOptions{};
   Postsolve<double> postsolve = Postsolve<double>( problem, num );
   ProblemUpdate<double> problemUpdate( problem, postsolve, statistics,
                                        presolveOptions, num );

   problem.recomputeAllActivities();
   papilo::DualFix<double> presolvingMethod{};
   Reductions<double> reductions{};

   PresolveStatus presolveStatus =
       presolvingMethod.execute( problem, problemUpdate, num, reductions );

   REQUIRE( presolveStatus == PresolveStatus::kReduced );
   REQUIRE( reductions.getReduction( 0 ).col == 0 );
   REQUIRE( reductions.getReduction( 0 ).row == papilo::ColReduction::LOCKED );
   REQUIRE( reductions.getReduction( 0 ).newval == 0 );
   REQUIRE( reductions.getReduction( 1 ).col == 0 );
   REQUIRE( reductions.getReduction( 1 ).row == papilo::ColReduction::FIXED );
   REQUIRE( reductions.getReduction( 1 ).newval == 1 );
}

TEST_CASE( "happy_path_dual_substitution", "[presolve]" )
{
   Num<double> num{};
   Problem<double> problem = setupMatrixForDualSubstitution();
   Statistics statistics{};
   PresolveOptions presolveOptions{};
   Postsolve<double> postsolve = Postsolve<double>( problem, num );
   ProblemUpdate<double> problemUpdate( problem, postsolve, statistics,
                                        presolveOptions, num );

   problemUpdate.trivialPresolve();
   papilo::DualFix<double> presolvingMethod{};
   Reductions<double> reductions{};

   PresolveStatus presolveStatus =
       presolvingMethod.execute( problem, problemUpdate, num, reductions );

   REQUIRE( presolveStatus == PresolveStatus::kReduced );
   REQUIRE( reductions.getReduction( 0 ).col == 0 );
   REQUIRE( reductions.getReduction( 0 ).row == papilo::ColReduction::LOCKED );
   REQUIRE( reductions.getReduction( 0 ).newval == 0 );

   REQUIRE( reductions.getReduction( 1 ).col == 0 );
   REQUIRE( reductions.getReduction( 1 ).row ==
            papilo::ColReduction::BOUNDS_LOCKED );
   REQUIRE( reductions.getReduction( 1 ).newval == 0 );

   REQUIRE( reductions.getReduction( 2 ).col == 0 );
   REQUIRE( reductions.getReduction( 2 ).row ==
            papilo::ColReduction::UPPER_BOUND );
   REQUIRE( reductions.getReduction( 2 ).newval == 6 );
}

TEST_CASE( "happy_path_dual_substitution_rounding", "[presolve]" )
{
   Num<double> num{};
   Problem<double> problem = setupMatrixForDualSubstitutionIntegerRounding();
   Statistics statistics{};
   PresolveOptions presolveOptions{};
   Postsolve<double> postsolve = Postsolve<double>( problem, num );
   ProblemUpdate<double> problemUpdate( problem, postsolve, statistics,
                                        presolveOptions, num );

   problemUpdate.trivialPresolve();
   papilo::DualFix<double> presolvingMethod{};
   Reductions<double> reductions{};

   PresolveStatus presolveStatus =
       presolvingMethod.execute( problem, problemUpdate, num, reductions );

   REQUIRE( presolveStatus == PresolveStatus::kReduced );
   REQUIRE( reductions.getReduction( 0 ).col == 0 );
   REQUIRE( reductions.getReduction( 0 ).row == papilo::ColReduction::LOCKED );
   REQUIRE( reductions.getReduction( 0 ).newval == 0 );

   REQUIRE( reductions.getReduction( 1 ).col == 0 );
   REQUIRE( reductions.getReduction( 1 ).row ==
            papilo::ColReduction::BOUNDS_LOCKED );
   REQUIRE( reductions.getReduction( 1 ).newval == 0 );

   REQUIRE( reductions.getReduction( 2 ).col == 0 );
   REQUIRE( reductions.getReduction( 2 ).row ==
            papilo::ColReduction::UPPER_BOUND );
   REQUIRE( reductions.getReduction( 2 ).newval == 6 );
}

TEST_CASE( "happy_path_dual_substitution_unbounded_variables", "[presolve]" )
{
   Num<double> num{};
   Problem<double> problem = setupMatrixForDualSubstitutionWithUnboundedVar();
   Statistics statistics{};
   PresolveOptions presolveOptions{};
   Postsolve<double> postsolve = Postsolve<double>( problem, num );
   ProblemUpdate<double> problemUpdate( problem, postsolve, statistics,
                                        presolveOptions, num );

   problem.recomputeAllActivities();
   papilo::DualFix<double> presolvingMethod{};
   Reductions<double> reductions{};

   PresolveStatus presolveStatus =
       presolvingMethod.execute( problem, problemUpdate, num, reductions );

   REQUIRE( presolveStatus == PresolveStatus::kUnchanged );
}

TEST_CASE( "happy_path_dual_substitution_for_equations", "[presolve]" )
{
   Num<double> num{};
   Problem<double> problem = setupMatrixForDualSubstitutionEquation();
   Statistics statistics{};
   PresolveOptions presolveOptions{};
   Postsolve<double> postsolve = Postsolve<double>( problem, num );
   ProblemUpdate<double> problemUpdate( problem, postsolve, statistics,
                                        presolveOptions, num );

   problem.recomputeAllActivities();
   papilo::DualFix<double> presolvingMethod{};
   Reductions<double> reductions{};

   PresolveStatus presolveStatus =
       presolvingMethod.execute( problem, problemUpdate, num, reductions );

   // TODO:
   REQUIRE( presolveStatus == PresolveStatus::kReduced );
   REQUIRE( reductions.getReduction( 0 ).col == 1 );
   REQUIRE( reductions.getReduction( 0 ).row == papilo::ColReduction::LOCKED );
   REQUIRE( reductions.getReduction( 0 ).newval == 0 );

   REQUIRE( reductions.getReduction( 1 ).col == 1 );
   REQUIRE( reductions.getReduction( 1 ).row ==
            papilo::ColReduction::BOUNDS_LOCKED );
   REQUIRE( reductions.getReduction( 1 ).newval == 0 );

   REQUIRE( reductions.getReduction( 2 ).col == 1 );
   REQUIRE( reductions.getReduction( 2 ).row ==
            papilo::ColReduction::LOWER_BOUND );
   REQUIRE( reductions.getReduction( 2 ).newval == 2 );
}

TEST_CASE( "example-4-from-4.4-Presolve-Reductions-in-MIP", "[presolve]" )
{

   Num<double> num{};
   Problem<double> problem = setupExample4ofChapter4Dot4InPresolveReductions();
   Statistics statistics{};
   PresolveOptions presolveOptions{};
   Postsolve<double> postsolve = Postsolve<double>( problem, num );
   ProblemUpdate<double> problemUpdate( problem, postsolve, statistics,
                                        presolveOptions, num );
   problemUpdate.trivialPresolve();
   papilo::DualFix<double> presolvingMethod{};
   Reductions<double> reductions{};

   PresolveStatus presolveStatus =
       presolvingMethod.execute( problem, problemUpdate, num, reductions );

   REQUIRE( presolveStatus == PresolveStatus::kUnchanged );
}

TEST_CASE( "example-5-from-4.4-Presolve-Reductions-in-MIP", "[presolve]" )
{
   Num<double> num{};
   Problem<double> problem = setupExample5ofChapter4Dot4InPresolveReductions();
   Statistics statistics{};
   PresolveOptions presolveOptions{};
   Postsolve<double> postsolve = Postsolve<double>( problem, num );
   ProblemUpdate<double> problemUpdate( problem, postsolve, statistics,
                                        presolveOptions, num );
   papilo::DualFix<double> presolvingMethod{};

   problemUpdate.trivialPresolve();

   Reductions<double> reductions{};

   PresolveStatus presolveStatus =
       presolvingMethod.execute( problem, problemUpdate, num, reductions );

   REQUIRE( presolveStatus == PresolveStatus::kReduced );
   REQUIRE( reductions.getReduction( 0 ).col == 0 );
   REQUIRE( reductions.getReduction( 0 ).row == papilo::ColReduction::LOCKED );
   REQUIRE( reductions.getReduction( 0 ).newval == 0 );

   REQUIRE( reductions.getReduction( 1 ).col == 0 );
   REQUIRE( reductions.getReduction( 1 ).row ==
            papilo::ColReduction::BOUNDS_LOCKED );
   REQUIRE( reductions.getReduction( 1 ).newval == 0 );

   REQUIRE( reductions.getReduction( 2 ).col == 0 );
   REQUIRE( reductions.getReduction( 2 ).row ==
            papilo::ColReduction::UPPER_BOUND );
   REQUIRE( reductions.getReduction( 2 ).newval == 2 );
}

Problem<double>
setupExample4ofChapter4Dot4InPresolveReductions()
{
   // min x + y + z
   // 2x + 4 y -3 z <= 8
   // - y - z <= -4
   // - x - y +8 z <= 0
   // 0<= x,y <= 4 z in {0,1}
   Vec<double> coefficients{ 1.0, 1.0, 1.0 };
   Vec<double> upperBounds{ 4.0, 4.0, 1.0 };
   Vec<double> lowerBounds{ 0.0, 0.0, 0.0 };
   Vec<uint8_t> isIntegral{ 0, 0, 1 };
   Vec<uint8_t> lhsInfinity{ 1, 1, 1 };
   Vec<uint8_t> rhsInfinity{ 0, 0, 0 };

   Vec<double> rhs{ 8.0, -4.0, 0.0 };
   Vec<std::string> rowNames{ "A1", "A2", "A3" };
   Vec<std::string> columnNames{ "x", "y", "z" };
   Vec<std::tuple<int, int, double>> entries{
       std::tuple<int, int, double>{ 0, 0, 2.0 },
       std::tuple<int, int, double>{ 0, 1, 4.0 },
       std::tuple<int, int, double>{ 0, 2, -3.0 },
       std::tuple<int, int, double>{ 1, 1, -1.0 },
       std::tuple<int, int, double>{ 1, 2, -1.0 },
       std::tuple<int, int, double>{ 2, 0, -1.0 },
       std::tuple<int, int, double>{ 2, 1, -1.0 },
       std::tuple<int, int, double>{ 2, 2, 8.0 } };

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
   pb.setRowLhsInfAll( lhsInfinity );
   pb.setRowRhsInfAll( rhsInfinity );
   pb.addEntryAll( entries );
   pb.setColNameAll( columnNames );
   pb.setProblemName(
       "matrix Example 4 of chapter 4.4 in Presolve Reductions in MIP" );
   Problem<double> problem = pb.build();
   return problem;
}

Problem<double>
setupExample5ofChapter4Dot4InPresolveReductions()
{
   // min x + y + z
   // 2x + 4 y -3 z <= 8
   // - y - z <= -4
   // - 2x - 2y + z <= 6
   // 0<= x,y,z <= 4
   Vec<double> coefficients{ 1.0, 1.0, 1.0 };
   Vec<double> upperBounds{ 10.0, 10.0, 10.0 };
   Vec<double> lowerBounds{ 0.0, 0.0, 0.0 };
   Vec<uint8_t> isIntegral{ 0, 0, 0 };
   Vec<uint8_t> lhsInfinity{ 1, 1, 1 };
   Vec<uint8_t> rhsInfinity{ 0, 0, 0 };

   Vec<double> rhs{ 8.0, -4.0, 6.0 };
   Vec<std::string> rowNames{ "A1", "A2", "A3" };
   Vec<std::string> columnNames{ "x", "y", "z" };
   Vec<std::tuple<int, int, double>> entries{
       std::tuple<int, int, double>{ 0, 0, 2.0 },
       std::tuple<int, int, double>{ 0, 1, 4.0 },
       std::tuple<int, int, double>{ 0, 2, -3.0 },
       std::tuple<int, int, double>{ 1, 1, -1.0 },
       std::tuple<int, int, double>{ 1, 2, -1.0 },
       std::tuple<int, int, double>{ 2, 0, -2.0 },
       std::tuple<int, int, double>{ 2, 1, -2.0 },
       std::tuple<int, int, double>{ 2, 2, 1.0 } };

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
   pb.setRowLhsInfAll( lhsInfinity );
   pb.setRowRhsInfAll( rhsInfinity );
   pb.addEntryAll( entries );
   pb.setColNameAll( columnNames );
   pb.setProblemName(
       "matrix Example 5 of chapter 4.4 in Presolve Reductions in MIP" );
   Problem<double> problem = pb.build();
   return problem;
}

Problem<double>
setupMatrixForDualFixFirstColumnOnlyPositiveValues()
{
   // min x + y + z
   // 2x + 4 y -3 z <= 8
   // - y - z <= -4
   // 2x - 2y + z <= 6
   // 0<= x,y,z <= 4
   Vec<double> coefficients{ 1.0, 1.0, 1.0 };
   Vec<double> upperBounds{ 10.0, 10.0, 10.0 };
   Vec<double> lowerBounds{ 1.0, 1.0, 1.0 };
   Vec<uint8_t> isIntegral{ 0, 0, 0 };
   Vec<uint8_t> lhsInfinity{ 1, 1, 1 };
   Vec<uint8_t> rhsInfinity{ 0, 0, 0 };

   Vec<double> rhs{ 8.0, -4.0, 6.0 };
   Vec<std::string> rowNames{ "A1", "A2", "A3" };
   Vec<std::string> columnNames{ "x", "y", "z" };
   Vec<std::tuple<int, int, double>> entries{
       std::tuple<int, int, double>{ 0, 0, 2.0 },
       std::tuple<int, int, double>{ 0, 1, 4.0 },
       std::tuple<int, int, double>{ 0, 2, -3.0 },
       std::tuple<int, int, double>{ 1, 1, -1.0 },
       std::tuple<int, int, double>{ 1, 2, -1.0 },
       std::tuple<int, int, double>{ 2, 0, 2.0 },
       std::tuple<int, int, double>{ 2, 1, -2.0 },
       std::tuple<int, int, double>{ 2, 2, 1.0 } };

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
   pb.setRowLhsInfAll( lhsInfinity );
   pb.setRowRhsInfAll( rhsInfinity );
   pb.addEntryAll( entries );
   pb.setColNameAll( columnNames );
   pb.setProblemName( "matrix for dual fix 1st row positive" );
   Problem<double> problem = pb.build();
   return problem;
}

Problem<double>
setupMatrixForDualSubstitution()
{
   // min y - z
   // 4 y -3 z <= 6
   // - y - z <= -3
   //- y + z <= 4
   // 0<= x,y,z <= 10
   Vec<double> coefficients{ 1.0, -1.0 };
   Vec<double> upperBounds{ 10.0, 10.0 };
   Vec<double> lowerBounds{ 0.0, 0.0 };
   Vec<uint8_t> isIntegral{ 0, 0 };
   Vec<uint8_t> lhsInfinity{ 1, 1, 1 };
   Vec<uint8_t> rhsInfinity{ 0, 0, 0 };

   Vec<double> rhs{ 6.0, -3.0, 4.0 };
   Vec<std::string> rowNames{ "A1", "A2", "A3" };
   Vec<std::string> columnNames{ "y", "z" };
   Vec<std::tuple<int, int, double>> entries{
       std::tuple<int, int, double>{ 0, 0, 4.0 },
       std::tuple<int, int, double>{ 0, 1, -3.0 },
       std::tuple<int, int, double>{ 1, 0, -1.0 },
       std::tuple<int, int, double>{ 1, 1, -1.0 },
       std::tuple<int, int, double>{ 2, 0, -1.0 },
       std::tuple<int, int, double>{ 2, 1, 1.0 } };

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
   pb.setRowLhsInfAll( lhsInfinity );
   pb.setRowRhsInfAll( rhsInfinity );
   pb.addEntryAll( entries );
   pb.setColNameAll( columnNames );
   pb.setProblemName( "matrix for dual substitution" );
   Problem<double> problem = pb.build();
   return problem;
}

Problem<double>
setupMatrixForDualSubstitutionIntegerRounding()
{
   // min y - z
   // 4 y -3 z <= 6
   //- y + z <= 4.2
   // 0<= x,y,z <= 10
   Vec<double> coefficients{ 1.0, -1.0 };
   Vec<double> upperBounds{ 10.0, 10.0 };
   Vec<double> lowerBounds{ 0.0, 0.0 };
   Vec<uint8_t> isIntegral{ 1, 1 };
   Vec<uint8_t> lhsInfinity{ 1, 1 };
   Vec<uint8_t> rhsInfinity{ 0, 0 };

   Vec<double> rhs{ 6.0, 4.2 };
   Vec<std::string> rowNames{ "A1", "A2" };
   Vec<std::string> columnNames{ "y", "z" };
   Vec<std::tuple<int, int, double>> entries{
       std::tuple<int, int, double>{ 0, 0, 4.0 },
       std::tuple<int, int, double>{ 0, 1, -3.0 },
       std::tuple<int, int, double>{ 1, 0, -1.0 },
       std::tuple<int, int, double>{ 1, 1, 1.0 } };

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
   pb.setRowLhsInfAll( lhsInfinity );
   pb.setRowRhsInfAll( rhsInfinity );
   pb.addEntryAll( entries );
   pb.setColNameAll( columnNames );
   pb.setProblemName(
       "matrix for dual substitution (new non integer bound on integer)" );
   Problem<double> problem = pb.build();
   return problem;
}

Problem<double>
setupMatrixForDualSubstitutionWithUnboundedVar()
{
   // min y - z
   // 4 y -3 z <= 6
   // - y - z <= -3
   //- y + z <= 4
   Vec<double> coefficients{ 1.0, -1.0 };
   Vec<uint8_t> isIntegral{ 0, 0 };
   Vec<uint8_t> lhsInfinity{ 1, 1, 1 };
   Vec<uint8_t> infinity{ 1, 1 };
   Vec<uint8_t> rhsInfinity{ 0, 0, 0 };

   Vec<double> rhs{ 6.0, -3.0, 4.0 };
   Vec<std::string> rowNames{ "A1", "A2", "A3" };
   Vec<std::string> columnNames{ "y", "z" };
   Vec<std::tuple<int, int, double>> entries{
       std::tuple<int, int, double>{ 0, 0, 4.0 },
       std::tuple<int, int, double>{ 0, 1, -3.0 },
       std::tuple<int, int, double>{ 1, 0, -1.0 },
       std::tuple<int, int, double>{ 1, 1, -1.0 },
       std::tuple<int, int, double>{ 2, 0, -1.0 },
       std::tuple<int, int, double>{ 2, 1, 1.0 } };

   ProblemBuilder<double> pb;
   pb.reserve( entries.size(), rowNames.size(), columnNames.size() );
   pb.setNumRows( rowNames.size() );
   pb.setNumCols( columnNames.size() );
   pb.setColLbInfAll( infinity );
   pb.setColUbInfAll( infinity );
   pb.setObjAll( coefficients );
   pb.setObjOffset( 0.0 );
   pb.setColIntegralAll( isIntegral );
   pb.setRowRhsAll( rhs );
   pb.setRowLhsInfAll( lhsInfinity );
   pb.setRowRhsInfAll( rhsInfinity );
   pb.addEntryAll( entries );
   pb.setColNameAll( columnNames );
   pb.setProblemName( "matrix for dual substitution" );
   Problem<double> problem = pb.build();
   return problem;
}

Problem<double>
setupMatrixForDualSubstitutionEquation()
{
   Vec<double> coefficients{ 0.0, 0.0 };
   Vec<double> upperBounds{ 1.0, 3.0 };
   Vec<double> lowerBounds{ 0.0, 0.0 };
   Vec<uint8_t> isIntegral{ 1, 1 };
   Vec<uint8_t> lhsInfinity{ 0, 0 };
   Vec<uint8_t> rhsInfinity{ 0, 0 };

   Vec<double> rhs{ 3.0, 4.0 };
   Vec<std::string> rowNames{ "A1", "a" };
   Vec<std::string> columnNames{ "y", "z" };
   Vec<std::tuple<int, int, double>> entries{
       std::tuple<int, int, double>{ 0, 0, 1.0 },
       std::tuple<int, int, double>{ 0, 1, 1.0 },
       std::tuple<int, int, double>{ 1, 0, 2.0 },
       std::tuple<int, int, double>{ 1, 1, 1.0 },
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
   pb.setRowLhsInfAll( lhsInfinity );
   pb.setRowRhsInfAll( rhsInfinity );
   pb.addEntryAll( entries );
   pb.setColNameAll( columnNames );
   pb.setProblemName( "matrix for dual substitution only equations" );
   Problem<double> problem = pb.build();
   problem.getConstraintMatrix().modifyLeftHandSide( 0, rhs[0] );
   problem.getConstraintMatrix().modifyLeftHandSide( 1, rhs[1] );
   return problem;
}
