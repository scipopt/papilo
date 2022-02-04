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

#include "fix/FixAndPropagate.hpp"

#include "fix/strategy/FarkasRoundingStrategy.hpp"
#include "catch/catch.hpp"
#include "papilo/core/ProbingView.hpp"
#include "papilo/core/Problem.hpp"
#include "papilo/core/ProblemBuilder.hpp"
#include "fix/strategy/FractionalRoundingStrategy.hpp"
#include "fix/strategy/RandomRoundingStrategy.hpp"

Problem<double>
setupProblemForFixAndPropagation();


TEST_CASE( "fix-and-propagate-integer-variable", "[fix]" )
{
   Problem<double> problem = setupProblemForFixAndPropagation();

   problem.getUpperBounds()[0] = 3;
   problem.getConstraintMatrix().getRightHandSides()[0] = 4;
   problem.getConstraintMatrix().getLeftHandSides()[0] = 4;

   problem.recomputeAllActivities();

   Vec<double> primal_solution = { 0.1, 0.2, 0.3, 0.4 };

   primal_solution[0] = 2.6;
   Vec<double> res{ primal_solution };

   FixAndPropagate<double> fixAndPropagate{ {}, {}, problem, {problem, {}}, true  };
   FractionalRoundingStrategy<double> strategy{ {} };

   bool infeasible =
       fixAndPropagate.fix_and_propagate( primal_solution, res, strategy );

   assert( !infeasible );
   assert( res[0] == 3 );
   assert( res[1] == 1 );
   assert( res[2] == 0 );
   assert( res[3] == 0 );
}

TEST_CASE( "fix-and-propagate-cont-variable-stays-cont", "[fix]" )
{
   Problem<double> problem = setupProblemForFixAndPropagation();

   problem.getColFlags()[0].unset(ColFlag::kIntegral);
   problem.getColFlags()[1].unset(ColFlag::kIntegral);
   problem.getColFlags()[2].unset(ColFlag::kIntegral);
   problem.recomputeAllActivities();

   Vec<double> primal_solution = { 0.1, 0.2, 0.3, 0.4 };

   Vec<double> res{ primal_solution };

   FixAndPropagate<double> fixAndPropagate{ {}, {}, problem, {problem, {}}, true  };
   FractionalRoundingStrategy<double> strategy{ {} };

   bool infeasible =
       fixAndPropagate.fix_and_propagate( primal_solution, res, strategy );

   assert( !infeasible );
   assert( res[0] == 0.1 );
   assert( res[1] == 0.8999999999999999 );
   assert( res[2] == 1 );
   assert( res[3] == 0 );
}



TEST_CASE( "fix-and-propagate-all-integer-solutions", "[fix]" )
{
   Problem<double> problem = setupProblemForFixAndPropagation();

   problem.recomputeAllActivities();

   Vec<double> primal_solution = { 0, 0, 0, 0 };

   Vec<double> res{ primal_solution };

   FixAndPropagate<double> fixAndPropagate{ {}, {}, problem, {problem, {}}, true  };
   FractionalRoundingStrategy<double> strategy{ {} };

   bool infeasible =
       fixAndPropagate.fix_and_propagate( primal_solution, res, strategy );

   assert( !infeasible );
   assert( res[0] == 0 );
   assert( res[1] == 0 );
   assert( res[2] == 1 );
   assert( res[3] == 1 );
}

TEST_CASE( "fix-and-propagate-integer-can-not-be-fixed", "[fix]" )
{
   Problem<double> problem = setupProblemForFixAndPropagation();

   problem.recomputeAllActivities();

   Vec<double> primal_solution = { 0, 0, 0, 0 };

   Vec<double> res{ primal_solution };

   FixAndPropagate<double> fixAndPropagate{ {}, {}, problem, {problem, {}}, true  };
   FractionalRoundingStrategy<double> strategy{ {} };

   bool infeasible =
       fixAndPropagate.fix_and_propagate( primal_solution, res, strategy );

   assert( !infeasible );
   assert( res[0] == 0 );
   assert( res[1] == 0 );
   assert( res[2] == 1 );
   assert( res[3] == 1 );
}

TEST_CASE( "fix-and-propagate-frac-backtrack", "[fix]" )
{
   Problem<double> problem = setupProblemForFixAndPropagation();

   problem.recomputeAllActivities();

   Vec<double> primal_solution = { 0.1, 0.2, 0.3, 0.4 };

   Vec<double> res{ primal_solution };

   ProbingView<double> view{ problem, {} };
   FixAndPropagate<double> fixAndPropagate{ {}, {}, view, true };
   FractionalRoundingStrategy<double> strategy{ {} };

   bool infeasible =
       fixAndPropagate.fix_and_propagate( primal_solution, res, strategy );

   assert( !infeasible );
   assert( res[0] == 1 );
   assert( res[1] == 1 );
   assert( res[2] == 0 );
   assert( res[3] == 0 );
}

TEST_CASE( "fix-and-propagate-frac-check-within-bounds", "[fix]" )
{
   Problem<double> problem = setupProblemForFixAndPropagation();
   problem.getUpperBounds()[0] = 3;
   problem.getUpperBounds()[1] = 3;
   problem.getUpperBounds()[2] = 3;
   problem.getUpperBounds()[3] = 3;
   problem.getConstraintMatrix().getRightHandSides()[0] = 4;
   problem.getConstraintMatrix().getLeftHandSides()[0] = 4;
   problem.recomputeAllActivities();

   Vec<double> primal_solution = { 1.8, 1.8, 2.8, 2.8 };

   Vec<double> res{ primal_solution };

   FixAndPropagate<double> fixAndPropagate{ {}, {}, problem, {problem, {}}, true  };
   FractionalRoundingStrategy<double> strategy{ {} };

   bool infeasible =
       fixAndPropagate.fix_and_propagate( primal_solution, res, strategy );

   assert( !infeasible );
   assert( res[0] == 1 );
   assert( res[1] == 0 );
   assert( res[2] == 3 );
   assert( res[3] == 0 );
}

TEST_CASE( "fix-and-propagate-random-backtrack", "[fix]" )
{
   Problem<double> problem = setupProblemForFixAndPropagation();

   problem.recomputeAllActivities();

   Vec<double> primal_solution = { 0.1, 0.2, 0.3, 0.4 };

   Vec<double> res{ primal_solution };

   ProbingView<double> view{ problem, {} };
   FixAndPropagate<double> fixAndPropagate{ {}, {}, problem, view, true };
   RandomRoundingStrategy<double> strategy{ 0, {} };

   bool infeasible =
       fixAndPropagate.fix_and_propagate( primal_solution, res, strategy );

   assert( !infeasible );
   assert( res[0] == 0 );
   assert( res[1] == 0 );
   assert( res[2] == 1 );
   assert( res[3] == 1 );
}


TEST_CASE( "fix-and-propagate-random", "[fix]" )
{
   Problem<double> problem = setupProblemForFixAndPropagation();

   problem.recomputeAllActivities();

   Vec<double> primal_solution = { 0.1, 0.2, 0.3, 0.4 };

   Vec<double> res{ primal_solution };

   ProbingView<double> view{ problem, {} };
   FixAndPropagate<double> fixAndPropagate{ {}, {}, problem, view, false };
   RandomRoundingStrategy<double> strategy{ 0, {} };

   bool infeasible =
       fixAndPropagate.fix_and_propagate( primal_solution, res, strategy );

   assert( infeasible );
   assert( res[0] == 1 );
   assert( res[1] == 1 );
   assert( res[2] == 1 );
   assert( res[3] == 1 );
}

TEST_CASE( "fix-and-propagate-random-check-within-bounds", "[fix]" )
{
   Problem<double> problem = setupProblemForFixAndPropagation();
   problem.getUpperBounds()[0] = 3;
   problem.getUpperBounds()[1] = 3;
   problem.getUpperBounds()[2] = 3;
   problem.getUpperBounds()[3] = 3;
   problem.getConstraintMatrix().getRightHandSides()[0] = 4;
   problem.getConstraintMatrix().getLeftHandSides()[0] = 4;
   problem.recomputeAllActivities();

   Vec<double> primal_solution = { 1.8, 1.8, 2.8, 2.8 };

   Vec<double> res{ primal_solution };

   FixAndPropagate<double> fixAndPropagate{ {}, {}, problem, {problem, {}}, true  };
   RandomRoundingStrategy<double> strategy{ 0, {} };

   bool infeasible =
       fixAndPropagate.fix_and_propagate( primal_solution, res, strategy );

   assert( !infeasible );
   assert( res[0] == 1 );
   assert( res[1] == 0 );
   assert( res[2] == 3 );
   assert( res[3] == 0 );
}

TEST_CASE( "fix-and-propagate-farkas-backtrack", "[fix]" )
{
   Problem<double> problem = setupProblemForFixAndPropagation();

   problem.recomputeAllActivities();

   Vec<double> primal_solution = { 0.1, 0.2, 0.3, 0.4 };

   Vec<double> res{ primal_solution };

   ProbingView<double> view{ problem, {} };
   FixAndPropagate<double> fixAndPropagate{ {}, {}, problem, view, true };
   FarkasRoundingStrategy<double> strategy{ 0, {}, false };

   bool infeasible =
       fixAndPropagate.fix_and_propagate( primal_solution, res, strategy );

   assert( !infeasible );
   assert( res[0] == 1 );
   assert( res[1] == 1 );
   assert( res[2] == 0 );
   assert( res[3] == 0 );
}

TEST_CASE( "fix-and-propagate-farkas-check-within-bounds", "[fix]" )
{
   Problem<double> problem = setupProblemForFixAndPropagation();
   problem.getUpperBounds()[0] = 1;
   problem.getUpperBounds()[1] = 3;
   problem.getUpperBounds()[2] = 3;
   problem.getUpperBounds()[3] = 3;
   problem.getConstraintMatrix().getRightHandSides()[0] = 4;
   problem.getConstraintMatrix().getLeftHandSides()[0] = 4;
   problem.recomputeAllActivities();

   Vec<double> primal_solution = { 0.8, 1.8, 2.3, 2.3 };

   Vec<double> res{ primal_solution };

   FixAndPropagate<double> fixAndPropagate{ {}, {}, problem, {problem, {}}, true  };
   FarkasRoundingStrategy<double> strategy{ 0, {}, true };

   bool infeasible =
       fixAndPropagate.fix_and_propagate( primal_solution, res, strategy );

   assert( !infeasible );
   assert( res[0] == 0 );
   assert( res[1] == 1 );
   assert( res[2] == 1 );
   assert( res[3] == 2 );
}

Problem<double>
setupProblemForFixAndPropagation()
{
   Vec<double> coefficients{ 1.0, 2.0, 3.0, 4.0 };
   Vec<double> upperBounds{ 1.0, 1.0, 1.0, 1.0 };
   Vec<double> lowerBounds{ 0.0, 0.0, 0.0, 0.0 };
   Vec<uint8_t> isIntegral{ 1, 1, 1, 1 };

   Vec<double> rhs{ 2.0 };
   Vec<std::string> rowNames{ "A1" };
   Vec<std::string> columnNames{ "c1", "c2", "c3", "c4" };
   Vec<std::tuple<int, int, double>> entries{
       std::tuple<int, int, double>{ 0, 0, 1.0 },
       std::tuple<int, int, double>{ 0, 1, 1.0 },
       std::tuple<int, int, double>{ 0, 2, 1.0 },
       std::tuple<int, int, double>{ 0, 3, 1.0 },
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
   pb.setColIntegralAll( isIntegral );
   pb.setRowRhsAll( rhs );
   pb.addEntryAll( entries );
   pb.setColNameAll( columnNames );
   pb.setProblemName( "coefficient strengthening matrix" );
   Problem<double> problem = pb.build();
   problem.getConstraintMatrix().modifyLeftHandSide( 0, {}, rhs[0] );
   return problem;
}
