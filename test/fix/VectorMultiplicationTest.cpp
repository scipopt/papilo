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
#include "fix/VectorMultiplication.hpp"
#include "papilo/core/Problem.hpp"
#include "papilo/core/ProblemBuilder.hpp"

using namespace papilo;

Problem<double>
setupProblemForVectorMultiplication();

TEST_CASE( "vector-calc_b_minus_Ax", "[fix]" )
{
   VectorMultiplication<double> multiplication{};
   Problem<double> problem = setupProblemForVectorMultiplication();

   Vec<double> scalar{};
   scalar.push_back(2);
   scalar.push_back(3);
   scalar.push_back(3);

   Vec<double> subtract{};
   subtract.push_back(1);
   subtract.push_back(2);


   Vec<double> res = multiplication.calc_b_minus_Ax(
       problem.getConstraintMatrix(), scalar, subtract );

   REQUIRE( res.size() == problem.getNRows() );
   REQUIRE( res[0] == -7 );
   REQUIRE( res[1] == -19 );

}

TEST_CASE( "vector-calc_b_minus_xA", "[fix]" )
{
   VectorMultiplication<double> multiplication{};
   Problem<double> problem = setupProblemForVectorMultiplication();

   Vec<double> scalar{};
   scalar.push_back(2);
   scalar.push_back(3);

   Vec<double> subtract{};
   subtract.push_back(1);
   subtract.push_back(2);
   subtract.push_back(3);


   Vec<double> res = multiplication.calc_b_minus_xA(
       problem.getConstraintMatrix(), scalar, subtract );

   REQUIRE( res.size() == problem.getNCols() );
   REQUIRE( res[0] == -1 );
   REQUIRE( res[1] == -11 );
   REQUIRE( res[2] == -9 );
}

TEST_CASE( "vector-l2-norm", "[fix]" )
{
   VectorMultiplication<double> multiplication{};

   Vec<double> subtract{};
   subtract.push_back(1);
   subtract.push_back(2);
   subtract.push_back(2);

   REQUIRE( multiplication.l2_norm( subtract ) == 3 );
}

TEST_CASE( "vector-addition", "[fix]" )
{
   VectorMultiplication<double> multiplication{};
   Problem<double> problem = setupProblemForVectorMultiplication();

   Vec<double> b{};
   b.push_back(2);
   b.push_back(3);
   b.push_back(4);

   Vec<double> x{};
   x.push_back(1);
   x.push_back(2);
   x.push_back(3);


   Vec<double> res = multiplication.calc_b_plus_sx(
       b, 2, x );

   REQUIRE( res.size() == problem.getNCols() );
   REQUIRE( res[0] == 4 );
   REQUIRE( res[1] == 7 );
   REQUIRE( res[2] == 10 );
}

TEST_CASE( "vector-addition-2", "[fix]" )
{
   VectorMultiplication<double> multiplication{};
   Problem<double> problem = setupProblemForVectorMultiplication();

   Vec<double> b{};
   b.push_back(2);
   b.push_back(3);
   b.push_back(4);

   Vec<double> x{};
   x.push_back(1);
   x.push_back(2);
   x.push_back(3);


   Vec<double> res = multiplication.calc_qb_plus_sx(3,
       b, 2, x );

   REQUIRE( res.size() == problem.getNCols() );
   REQUIRE( res[0] == 8 );
   REQUIRE( res[1] == 13 );
   REQUIRE( res[2] == 18 );
}

Problem<double>
setupProblemForVectorMultiplication()
{
//   1x + 2y <= 2
   Vec<double> coefficients{ 1.0, 1.0, 1.0 };
   Vec<double> upperBounds{ 1.0, 1.0, 1.0 };
   Vec<double> lowerBounds{ 0.0, 0.0, 0.0 };
   Vec<uint8_t> isIntegral{ 1, 1, 1 };

   Vec<double> rhs{ 2.0, 3.0};
   Vec<std::string> rowNames{ "A1", "A2" };
   Vec<std::string> columnNames{ "c1", "c2", "c3" };
   Vec<std::tuple<int, int, double>> entries{
       std::tuple<int, int, double>{ 0, 0, 1.0 },
       std::tuple<int, int, double>{ 0, 1, 2.0 },
       std::tuple<int, int, double>{ 1, 1, 3.0 },
       std::tuple<int, int, double>{ 1, 2, 4.0 },
   };

   ProblemBuilder<double> pb;
   pb.reserve( (int) entries.size(), (int) rowNames.size(), (int) columnNames.size() );
   pb.setNumRows( (int) rowNames.size() );
   pb.setNumCols( (int) columnNames.size() );
   pb.setColUbAll( upperBounds );
   pb.setColLbAll( lowerBounds );
   pb.setObjAll( coefficients );
   pb.setObjOffset( 0.0 );
   pb.setColIntegralAll( isIntegral );
   pb.setRowRhsAll( rhs );
   pb.addEntryAll( entries );
   pb.setColNameAll( columnNames );
   pb.setProblemName( "coefficient strengthening matrix" );
   Problem<double> problem = pb.build();
   return problem;
}
