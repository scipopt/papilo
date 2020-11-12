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
#include "papilo/core/PresolveMethod.hpp"
#include "papilo/core/Problem.hpp"
#include "papilo/core/ProblemBuilder.hpp"
#include "papilo/core/RowFlags.hpp"
#include "papilo/presolvers/SingletonCols.hpp"

using namespace papilo;

Problem<double>
setupProblemWithOnlyOneEntryIn1stRowAndColumn();

Problem<double>
setupProblemWithSingletonColumn();

void
forceCalculationOfSingletonRows( Problem<double>& problem,
                                 ProblemUpdate<double>& problemUpdate )
{
   problem.recomputeLocks();
   problemUpdate.trivialColumnPresolve();
}

TEST_CASE( "happy path - singleton column", "[presolve]" )
{
   Num<double> num{};
   Problem<double> problem = setupProblemWithSingletonColumn();
   Statistics statistics{};
   PresolveOptions presolveOptions{};
   Postsolve<double> postsolve = Postsolve<double>( problem, num );
   ProblemUpdate<double> problemUpdate( problem, postsolve, statistics,
                                        presolveOptions, num );
   forceCalculationOfSingletonRows( problem, problemUpdate );
   problem.recomputeAllActivities();
   SingletonCols<double> presolvingMethod{};
   Reductions<double> reductions{};

   PresolveStatus presolveStatus =
       presolvingMethod.execute( problem, problemUpdate, num, reductions );

   BOOST_ASSERT( presolveStatus == PresolveStatus::kReduced );
   BOOST_ASSERT( reductions.size() == 5 );
   BOOST_ASSERT( reductions.getReduction( 0 ).col ==0);
   BOOST_ASSERT( reductions.getReduction( 0 ).row == ColReduction::BOUNDS_LOCKED);
   BOOST_ASSERT( reductions.getReduction( 0 ).newval == 0 );

   BOOST_ASSERT( reductions.getReduction( 1 ).row == 0 );
   BOOST_ASSERT( reductions.getReduction( 1 ).col ==
                 papilo::RowReduction::LOCKED );
   BOOST_ASSERT( reductions.getReduction( 1 ).newval == 0 );

   BOOST_ASSERT( reductions.getReduction( 2 ).col == 0 );
   BOOST_ASSERT( reductions.getReduction( 2 ).newval == 0 );
   BOOST_ASSERT( reductions.getReduction( 2 ).row == ColReduction::SUBSTITUTE_OBJ);

   //in matrix entry (0,0) new value 0
   BOOST_ASSERT( reductions.getReduction( 3 ).col == 0 );
   BOOST_ASSERT( reductions.getReduction( 3 ).newval == 0 );
   BOOST_ASSERT( reductions.getReduction( 3 ).row == 0);

   BOOST_ASSERT( reductions.getReduction( 4 ).col == RowReduction::LHS_INF );
   BOOST_ASSERT( reductions.getReduction( 4 ).newval == 0 );
   BOOST_ASSERT( reductions.getReduction( 4 ).row == 0);
}

TEST_CASE( "failed path - singleton column & row", "[presolve]" )
{
   Num<double> num{};
   Problem<double> problem = setupProblemWithOnlyOneEntryIn1stRowAndColumn();
   Statistics statistics{};
   PresolveOptions presolveOptions{};
   Postsolve<double> postsolve = Postsolve<double>( problem, num );
   ProblemUpdate<double> problemUpdate( problem, postsolve, statistics,
                                        presolveOptions, num );
   forceCalculationOfSingletonRows( problem, problemUpdate );
   SingletonCols<double> presolvingMethod{};
   Reductions<double> reductions{};

   PresolveStatus presolveStatus =
       presolvingMethod.execute( problem, problemUpdate, num, reductions );

   BOOST_ASSERT( presolveStatus == PresolveStatus::kUnchanged );
}

Problem<double>
setupProblemWithOnlyOneEntryIn1stRowAndColumn()
{
   Vec<double> coefficients{ 1.0, 1.0, 1.0 };
   Vec<double> upperBounds{ 10.0, 10.0, 10.0 };
   Vec<double> lowerBounds{ 0.0, 0.0, 0.0 };
   Vec<uint8_t> isIntegral{ 1, 1, 1 };

   Vec<uint8_t> isLefthandsideInfinity{ 0, 1, 1 };
   Vec<uint8_t> isRighthandsideInfinity{ 0, 0, 0 };
   Vec<double> rhs{ 1.0, 2.0, 3.0 };
   Vec<std::string> rowNames{ "A1", "A2", "A3" };
   Vec<std::string> columnNames{ "c1", "c2", "c3" };
   Vec<std::tuple<int, int, double>> entries{
       std::tuple<int, int, double>{ 0, 0, 1.0 },
       std::tuple<int, int, double>{ 1, 1, 2.0 },
       std::tuple<int, int, double>{ 2, 2, 3.0 },
       std::tuple<int, int, double>{ 1, 2, 3.0 },
       std::tuple<int, int, double>{ 2, 1, 4.0 },
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
   pb.setRowLhsInfAll( isLefthandsideInfinity );
   pb.setRowRhsInfAll( isRighthandsideInfinity );
   pb.addEntryAll( entries );
   pb.setColNameAll( columnNames );
   pb.setProblemName( "singleton column & row matrix with equation" );
   Problem<double> problem = pb.build();
   // TODO: if constructing the problem like this, why is there noch check on Equation?
   problem.getConstraintMatrix().modifyLeftHandSide( 0, 1 );
   return problem;
}

Problem<double>
setupProblemWithSingletonColumn()
{
   Vec<double> coefficients{ 1.0, 1.0, 1.0 };
   Vec<double> upperBounds{ 10.0, 10.0, 10.0 };
   Vec<double> lowerBounds{ 0.0, 0.0, 0.0 };
   Vec<uint8_t> isIntegral{ 1, 1, 1 };

   Vec<uint8_t> isLefthandsideInfinity{ 0, 1, 1 };
   Vec<uint8_t> isRighthandsideInfinity{ 0, 0, 0 };
   Vec<double> rhs{ 1.0, 2.0, 3.0 };
   Vec<std::string> rowNames{ "A1", "A2", "A3" };
   Vec<std::string> columnNames{ "c1", "c2", "c3" };
   Vec<std::tuple<int, int, double>> entries{
       std::tuple<int, int, double>{ 0, 0, 1.0 },
       std::tuple<int, int, double>{ 0, 1, 1.0 },
       std::tuple<int, int, double>{ 1, 1, 2.0 },
       std::tuple<int, int, double>{ 2, 2, 3.0 },
       std::tuple<int, int, double>{ 1, 2, 3.0 },
       std::tuple<int, int, double>{ 2, 1, 4.0 },
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
   pb.setRowLhsInfAll( isLefthandsideInfinity );
   pb.setRowRhsInfAll( isRighthandsideInfinity );
   pb.addEntryAll( entries );
   pb.setColNameAll( columnNames );
   pb.setProblemName( "singleton column matrix with equation" );
   Problem<double> problem = pb.build();
   // TODO: is there no check if the l
   problem.getConstraintMatrix().modifyLeftHandSide( 0, 1 );
   return problem;
}

