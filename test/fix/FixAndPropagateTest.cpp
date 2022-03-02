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
#include "fix/strategy/MostFractionalRoundingStrategy.hpp"
#include "fix/strategy/LeastFractionalRoundingStrategy.hpp"

Problem<double>
setupProblemForFixAndPropagation();

Problem<double>
setupProblemForConflictAnalysis_2();

// TODO: add tests for conflict diving

TEST_CASE( "fix-and-propagate-most-frac-backtrack", "[fix]" )
{
   Problem<double> problem = setupProblemForFixAndPropagation();

   problem.recomputeAllActivities();

   Vec<double> primal_solution = { 0.1, 0.2, 0.3, 0.6 };

   Vec<double> res{ primal_solution };

   Message msg{};
   msg.setVerbosityLevel(papilo::VerbosityLevel::kDetailed);
   RandomGenerator random( 0 );
   FixAndPropagate<double> fixAndPropagate{ msg, {}, random };
   MostFractionalRoundingStrategy<double> strategy{ {} };

   ProbingView<double> view {problem, {}};
   int backtracks = 0;
   bool infeasible =
       fixAndPropagate.fix_and_propagate( primal_solution, res, strategy, view,
                                          backtracks, true, false );

   REQUIRE( !infeasible );
   REQUIRE( res[0] == 1 );
   REQUIRE( res[1] == 0 );
   REQUIRE( res[2] == 0 );
   REQUIRE( res[3] == 1 );
}

TEST_CASE( "fix-and-propagate-most-frac-only-0.5-backtrack", "[fix]" )
{
   Problem<double> problem = setupProblemForFixAndPropagation();

   problem.recomputeAllActivities();

   Vec<double> primal_solution = { 0.5, 0.5, 0.5, 0.5 };

   Vec<double> res{ primal_solution };

   Message msg{};
   msg.setVerbosityLevel(papilo::VerbosityLevel::kDetailed);
   RandomGenerator random( 0 );
   FixAndPropagate<double> fixAndPropagate{ msg, {}, random };
   MostFractionalRoundingStrategy<double> strategy{ {} };

   ProbingView<double> view {problem, {}};
   int backtracks = 0;
   bool infeasible =
       fixAndPropagate.fix_and_propagate( primal_solution, res, strategy, view,
                                          backtracks, true, false );

   REQUIRE( !infeasible );
   REQUIRE( res[0] == 0 );
   REQUIRE( res[1] == 0 );
   REQUIRE( res[2] == 0 );
   REQUIRE( res[3] == 1 );
}

TEST_CASE( "fix-and-propagate-least-frac-backtrack", "[fix]" )
{
   Problem<double> problem = setupProblemForFixAndPropagation();

   problem.recomputeAllActivities();

   Vec<double> primal_solution = { 0.1, 0.8, 0.3, 0.5 };

   Vec<double> res{ primal_solution };

   Message msg{};
   msg.setVerbosityLevel(papilo::VerbosityLevel::kDetailed);
   RandomGenerator random( 0 );
   FixAndPropagate<double> fixAndPropagate{ msg, {}, random };
   LeastFractionalRoundingStrategy<double> strategy{ {} };

   ProbingView<double> view {problem, {}};
   int backtracks = 0;
   bool infeasible =
       fixAndPropagate.fix_and_propagate( primal_solution, res, strategy, view,
                                          backtracks, true, false );

   REQUIRE( !infeasible );
   REQUIRE( res[0] == 0 );
   REQUIRE( res[1] == 1 );
   REQUIRE( res[2] == 0 );
   REQUIRE( res[3] == 1 );
}

TEST_CASE( "fix-and-propagate-least-frac-only-0.5-backtrack", "[fix]" )
{
   Problem<double> problem = setupProblemForFixAndPropagation();

   problem.recomputeAllActivities();

   Vec<double> primal_solution = { 0.5, 0.5, 0.5, 0.5 };

   Vec<double> res{ primal_solution };

   Message msg{};
   msg.setVerbosityLevel(papilo::VerbosityLevel::kDetailed);
   RandomGenerator random( 0 );
   FixAndPropagate<double> fixAndPropagate{ msg, {}, random };
   LeastFractionalRoundingStrategy<double> strategy{ {} };

   ProbingView<double> view {problem, {}};
   int backtracks = 0;
   bool infeasible =
       fixAndPropagate.fix_and_propagate( primal_solution, res, strategy, view,
                                          backtracks, true, false );

   REQUIRE( !infeasible );
   REQUIRE( res[0] == 0 );
   REQUIRE( res[1] == 0 );
   REQUIRE( res[2] == 1 );
   REQUIRE( res[3] == 1 );
}

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

   Message msg{};
   msg.setVerbosityLevel(papilo::VerbosityLevel::kDetailed);
   RandomGenerator random( 0 );
   FixAndPropagate<double> fixAndPropagate{ msg, {}, random };   FractionalRoundingStrategy<double> strategy{ {}, problem };

   ProbingView<double> view {problem, {}};
   int backtracks = 0;
   bool infeasible =
       fixAndPropagate.fix_and_propagate( primal_solution, res, strategy, view, backtracks, true, false  );

   REQUIRE( !infeasible );
   REQUIRE( res[0] == 3 );
   REQUIRE( res[1] == 1 );
   REQUIRE( res[2] == 0 );
   REQUIRE( res[3] == 0 );
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

   Message msg{};
   msg.setVerbosityLevel(papilo::VerbosityLevel::kDetailed);
   RandomGenerator random( 0 );
   FixAndPropagate<double> fixAndPropagate{ msg, {}, random };
   FractionalRoundingStrategy<double> strategy{ {}, problem };

   ProbingView<double> view {problem, {}};
   int backtracks = 0;

   bool infeasible =
       fixAndPropagate.fix_and_propagate( primal_solution, res, strategy, view, backtracks, true, false  );

   REQUIRE( !infeasible );
   REQUIRE( res[0] == 0.1 );
   REQUIRE( res[1] == 0.8999999999999999 );
   REQUIRE( res[2] == 1 );
   REQUIRE( res[3] == 0 );
}



TEST_CASE( "fix-and-propagate-all-integer-solutions", "[fix]" )
{
   Problem<double> problem = setupProblemForFixAndPropagation();

   problem.recomputeAllActivities();

   Vec<double> primal_solution = { 0, 0, 0, 0 };

   Vec<double> res{ primal_solution };

   Message msg{};
   msg.setVerbosityLevel(papilo::VerbosityLevel::kDetailed);
   RandomGenerator random( 0 );
   FixAndPropagate<double> fixAndPropagate{ msg, {}, random };   FractionalRoundingStrategy<double> strategy{ {}, problem };

   ProbingView<double> view {problem, {}};
   int backtracks = 0;
   bool infeasible =
       fixAndPropagate.fix_and_propagate( primal_solution, res, strategy, view, backtracks, true, false  );


   REQUIRE( !infeasible );
   REQUIRE( res[0] == 0 );
   REQUIRE( res[1] == 0 );
   REQUIRE( res[2] == 1 );
   REQUIRE( res[3] == 1 );
}

TEST_CASE( "fix-and-propagate-integer-can-not-be-fixed", "[fix]" )
{
   Problem<double> problem = setupProblemForFixAndPropagation();

   problem.recomputeAllActivities();

   Vec<double> primal_solution = { 0, 0, 0, 0 };

   Vec<double> res{ primal_solution };

   Message msg{};
   msg.setVerbosityLevel(papilo::VerbosityLevel::kDetailed);
   RandomGenerator random( 0 );
   FixAndPropagate<double> fixAndPropagate{ msg, {}, random };
   FractionalRoundingStrategy<double> strategy{ {}, problem };

   ProbingView<double> view {problem, {}};
   int backtracks = 0;
   bool infeasible =
       fixAndPropagate.fix_and_propagate( primal_solution, res, strategy, view, backtracks, true, false  );


   REQUIRE( !infeasible );
   REQUIRE( res[0] == 0 );
   REQUIRE( res[1] == 0 );
   REQUIRE( res[2] == 1 );
   REQUIRE( res[3] == 1 );
}

TEST_CASE( "fix-and-propagate-frac-backtrack", "[fix]" )
{
   Problem<double> problem = setupProblemForFixAndPropagation();

   problem.recomputeAllActivities();

   Vec<double> primal_solution = { 0.1, 0.2, 0.3, 0.4 };

   Vec<double> res{ primal_solution };

   Message msg{};
   msg.setVerbosityLevel(papilo::VerbosityLevel::kDetailed);
   RandomGenerator random( 0 );
   FixAndPropagate<double> fixAndPropagate{ msg, {}, random };
   FractionalRoundingStrategy<double> strategy{ {}, problem };

   ProbingView<double> view {problem, {}};
   int backtracks = 0;
   bool infeasible =
       fixAndPropagate.fix_and_propagate( primal_solution, res, strategy, view,backtracks, true, false  );


   REQUIRE( !infeasible );
   REQUIRE( res[0] == 1 );
   REQUIRE( res[1] == 1 );
   REQUIRE( res[2] == 0 );
   REQUIRE( res[3] == 0 );
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

   Message msg{};
   msg.setVerbosityLevel(papilo::VerbosityLevel::kDetailed);
   RandomGenerator random( 0 );
   FixAndPropagate<double> fixAndPropagate{ msg, {}, random };
   FractionalRoundingStrategy<double> strategy{ {}, problem };

   ProbingView<double> view {problem, {}};
   int backtracks = 0;
   bool infeasible =
       fixAndPropagate.fix_and_propagate( primal_solution, res, strategy, view, backtracks,true, false  );

   REQUIRE( !infeasible );
   REQUIRE( res[0] == 1 );
   REQUIRE( res[1] == 0 );
   REQUIRE( res[2] == 3 );
   REQUIRE( res[3] == 0 );
}

TEST_CASE( "fix-and-propagate-random-backtrack", "[fix]" )
{
   Problem<double> problem = setupProblemForFixAndPropagation();

   problem.recomputeAllActivities();

   Vec<double> primal_solution = { 0.1, 0.2, 0.3, 0.4 };

   Vec<double> res{ primal_solution };

   Message msg{};
   msg.setVerbosityLevel(papilo::VerbosityLevel::kDetailed);
   RandomGenerator random( 0 );
   FixAndPropagate<double> fixAndPropagate{ msg, {}, random };
   RandomRoundingStrategy<double> strategy{ 0, {} };

   ProbingView<double> view {problem, {}};
   int backtracks = 0;
   bool infeasible =
       fixAndPropagate.fix_and_propagate( primal_solution, res, strategy, view,backtracks, true, false  );

   REQUIRE( !infeasible );
   REQUIRE( res[0] == 0 );
   REQUIRE( res[1] == 0 );
   REQUIRE( res[2] == 1 );
   REQUIRE( res[3] == 1 );
}


TEST_CASE( "fix-and-propagate-random", "[fix]" )
{
   Problem<double> problem = setupProblemForFixAndPropagation();

   problem.recomputeAllActivities();

   Vec<double> primal_solution = { 0.1, 0.2, 0.3, 0.4 };

   Vec<double> res{ primal_solution };

   Message msg{};
   msg.setVerbosityLevel(papilo::VerbosityLevel::kDetailed);
   RandomGenerator random( 0 );
   FixAndPropagate<double> fixAndPropagate{ msg, {}, random };
   RandomRoundingStrategy<double> strategy{ 0, {} };

   ProbingView<double> view {problem, {}};
   int backtracks = 0;
   bool infeasible =
       fixAndPropagate.fix_and_propagate( primal_solution, res, strategy, view,backtracks, false, false  );

   REQUIRE( !infeasible );
   REQUIRE( res[0] == 0 );
   REQUIRE( res[1] == 0 );
   REQUIRE( res[2] == 1 );
   REQUIRE( res[3] == 1 );
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

   Message msg{};
   msg.setVerbosityLevel(papilo::VerbosityLevel::kDetailed);
   RandomGenerator random( 0 );
   FixAndPropagate<double> fixAndPropagate{ msg, {}, random };
   RandomRoundingStrategy<double> strategy{ 0, {} };

   ProbingView<double> view {problem, {}};
   int backtracks = 0;
   bool infeasible =
       fixAndPropagate.fix_and_propagate( primal_solution, res, strategy, view,backtracks, true, false  );

   REQUIRE( !infeasible );
   REQUIRE( res[0] == 1 );
   REQUIRE( res[1] == 0 );
   REQUIRE( res[2] == 3 );
   REQUIRE( res[3] == 0 );
}

TEST_CASE( "fix-and-propagate-farkas-backtrack", "[fix]" )
{
   Problem<double> problem = setupProblemForFixAndPropagation();

   problem.recomputeAllActivities();

   Vec<double> primal_solution = { 0.1, 0.2, 0.3, 0.4 };

   Vec<double> res{ primal_solution };

   Message msg{};
   msg.setVerbosityLevel(papilo::VerbosityLevel::kDetailed);
   RandomGenerator random( 0 );
   FixAndPropagate<double> fixAndPropagate{ msg, {}, random };
   FarkasRoundingStrategy<double> strategy{ 0, {}, false };

   ProbingView<double> view {problem, {}};
   int backtracks = 0;
   bool infeasible =
       fixAndPropagate.fix_and_propagate( primal_solution, res, strategy, view, backtracks,true, false  );

   REQUIRE( !infeasible );
   REQUIRE( res[0] == 1 );
   REQUIRE( res[1] == 1 );
   REQUIRE( res[2] == 0 );
   REQUIRE( res[3] == 0 );
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

   Message msg{};
   msg.setVerbosityLevel(papilo::VerbosityLevel::kDetailed);
   RandomGenerator random( 0 );
   FixAndPropagate<double> fixAndPropagate{ msg, {}, random };
   FarkasRoundingStrategy<double> strategy{ 0, {}, true };

   ProbingView<double> view {problem, {}};
   int backtracks = 0;
   bool infeasible =
       fixAndPropagate.fix_and_propagate( primal_solution, res, strategy, view, backtracks, true, false  );

   REQUIRE( !infeasible );
   REQUIRE( res[0] == 0 );
   REQUIRE( res[1] == 1 );
   REQUIRE( res[2] == 1 );
   REQUIRE( res[3] == 2 );
}

TEST_CASE( "fix-and-propagate-stop-at-infeas-false", "[fix]" )
{
   Problem<double> problem = setupProblemForConflictAnalysis_2();
   problem.recomputeAllActivities();
   // Binary problem with constraints
   // A1: x1 + x3 = 1
   // A2: x1 + x2 + x3 = 2
   // A3: x2 + x3 + x4 + x5 = 3
   // A4: x4 + x5 = 2

   Vec<double> primal_solution = { 0.9, 0.9, 0.6, 0.3, 0.2 };

   Vec<double> res{ primal_solution };

   Message msg{};
   msg.setVerbosityLevel(papilo::VerbosityLevel::kDetailed);
   RandomGenerator random( 0 );
   FixAndPropagate<double> fixAndPropagate{ msg, {}, random };
   FractionalRoundingStrategy<double> strategy{ {}, problem };

   ProbingView<double> view {problem, {}};
   int backtracks = 0;
   bool infeasible =
       fixAndPropagate.fix_and_propagate( primal_solution, res, strategy, view, backtracks, true, true  );

   REQUIRE( infeasible );
   // if stop at infeasible is true the res should not be modified if infeasible
   for( int i = 0; i < primal_solution.size(); i++ )
      REQUIRE( primal_solution[i] == res[i] );
}


TEST_CASE( "fix-and-propagate-initial-solution-fix-to-zero", "[fix]" )
{
   Problem<double> problem = setupProblemForFixAndPropagation();
   problem.recomputeAllActivities();

   Message msg{};
   msg.setVerbosityLevel(papilo::VerbosityLevel::kDetailed);
   RandomGenerator random( 0 );
   FixAndPropagate<double> fixAndPropagate{ msg, {}, random };
   FractionalRoundingStrategy<double> strategy{ {}, problem };
   Vec<double> res = { 0.9, 0.9, 0.6, 0.3, 0.2 };

   ProbingView<double> view {problem, {}};
   int backtracks = 0;
   bool infeasible =
       fixAndPropagate.find_initial_solution( 0, view, res  );

   REQUIRE( !infeasible );
   REQUIRE( res[0] == 0 );
   REQUIRE( res[1] == 0 );
   REQUIRE( res[2] == 1 );
   REQUIRE( res[3] == 1 );
}

TEST_CASE( "fix-and-propagate-initial-solution-fix-to-lowerbound", "[fix]" )
{
   Problem<double> problem = setupProblemForFixAndPropagation();
   problem.recomputeAllActivities();

   Message msg{};
   msg.setVerbosityLevel(papilo::VerbosityLevel::kDetailed);
   RandomGenerator random( 0 );
   FixAndPropagate<double> fixAndPropagate{ msg, {}, random };
   FractionalRoundingStrategy<double> strategy{ {}, problem };
   Vec<double> res = { 0.9, 0.9, 0.6, 0.3, 0.2 };

   ProbingView<double> view {problem, {}};
   int backtracks = 0;
   bool infeasible =
       fixAndPropagate.find_initial_solution( 1, view, res  );

   REQUIRE( !infeasible );
   REQUIRE( res[0] == 0 );
   REQUIRE( res[1] == 0 );
   REQUIRE( res[2] == 1 );
   REQUIRE( res[3] == 1 );
}

TEST_CASE( "fix-and-propagate-initial-solution-fix-to-upperbound", "[fix]" )
{
   Problem<double> problem = setupProblemForFixAndPropagation();
   problem.recomputeAllActivities();

   Message msg{};
   msg.setVerbosityLevel(papilo::VerbosityLevel::kDetailed);
   RandomGenerator random( 0 );
   FixAndPropagate<double> fixAndPropagate{ msg, {}, random };
   FractionalRoundingStrategy<double> strategy{ {}, problem };
   Vec<double> res = { 0.9, 0.9, 0.6, 0.3, 0.2 };

   ProbingView<double> view {problem, {}};
   int backtracks = 0;
   bool infeasible =
       fixAndPropagate.find_initial_solution( 2, view, res  );

   REQUIRE( !infeasible );
   REQUIRE( res[0] == 1 );
   REQUIRE( res[1] == 1 );
   REQUIRE( res[2] == 0 );
   REQUIRE( res[3] == 0 );
}


TEST_CASE( "fix-and-propagate-check-conflict-analysis-data", "[fix]" )
{
   Problem<double> problem = setupProblemForConflictAnalysis_2();
   problem.recomputeAllActivities();
   // Binary problem with constraints
   // A1: x1 + x3 = 1
   // A2: x1 + x2 + x3 = 2
   // A3: x2 + x3 + x4 + x5 = 3
   // A4: x4 + x5 = 2

   Vec<double> primal_solution = { 0.9, 0.9, 0.6, 0.3, 0.2 };

   Vec<double> res{ primal_solution };

   Message msg{};
   msg.setVerbosityLevel(papilo::VerbosityLevel::kDetailed);
   RandomGenerator random( 0 );
   FixAndPropagate<double> fixAndPropagate{ msg, {}, random };
   FractionalRoundingStrategy<double> strategy{ {}, problem };

   ProbingView<double> view {problem, {}};
   int backtracks = 0;
   bool infeasible =
       fixAndPropagate.fix_and_propagate( primal_solution, res, strategy, view, backtracks, true, false  );

   auto changes = view.get_changes();
   auto conflicts = view.get_infeasible_rows();

   REQUIRE( infeasible );
   REQUIRE( changes.size() == 5);
   REQUIRE( conflicts.size() == 1);

   REQUIRE( conflicts[0].first == 5);
   REQUIRE( conflicts[0].second == 3);

   REQUIRE( changes[0].get_new_bound_value() == 1);
   REQUIRE( changes[0].get_reason_row() == -1);
   REQUIRE( changes[0].is_lower_bound() );
   REQUIRE( changes[0].is_manually_triggered() );
   REQUIRE( changes[0].get_depth_level() == -2 );

   REQUIRE( changes[1].get_col() == 0);
   REQUIRE( changes[1].get_new_bound_value() == 0);
   REQUIRE( changes[1].get_reason_row() == 0);
   REQUIRE( !changes[1].is_lower_bound() );
   REQUIRE( !changes[1].is_manually_triggered() );
   REQUIRE( changes[1].get_depth_level() == -3 );

   REQUIRE( changes[2].get_col() == 1);
   REQUIRE( changes[2].get_new_bound_value() == 1);
   REQUIRE( changes[2].get_reason_row() == 1);
   REQUIRE( changes[2].is_lower_bound() );
   REQUIRE( !changes[2].is_manually_triggered() );
   REQUIRE( changes[2].get_depth_level() == -3 );

   REQUIRE( changes[3].get_col() == 3);
   REQUIRE( changes[3].get_new_bound_value() == 1);
   REQUIRE( changes[3].get_reason_row() == -1);
   REQUIRE( changes[3].is_lower_bound() );
   REQUIRE( changes[3].is_manually_triggered() );
   REQUIRE( changes[3].get_depth_level() == -4 );

   REQUIRE( changes[4].get_col() == 4);
   REQUIRE( changes[4].get_new_bound_value() == 0);
   REQUIRE( changes[4].get_reason_row() == 2);
   REQUIRE( !changes[4].is_lower_bound() );
   REQUIRE( !changes[4].is_manually_triggered() );
   REQUIRE( changes[4].get_depth_level() == -5 );

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

Problem<double>
setupProblemForConflictAnalysis_2()
{
   Vec<double> coefficients{ 1.0, 1.0, 1.0, 1.0, 1.0 };
   Vec<double> upperBounds{ 1.0, 1.0, 1.0, 1.0, 1.0 };
   Vec<double> lowerBounds{ 0.0, 0.0, 0.0, 0.0, 0.0 };
   Vec<uint8_t> isIntegral{ 1, 1, 1, 1, 1 };

   Vec<double> rhs{ 1.0, 2.0, 3.0, 2.0 };
   Vec<std::string> rowNames{ "A1", "A2", "A3", "A4" };
   Vec<std::string> columnNames{ "c1", "c2", "c3", "c4", "c5" };
   Vec<std::tuple<int, int, double>> entries{
       std::tuple<int, int, double>{ 0, 0, 1.0 },
       std::tuple<int, int, double>{ 0, 2, 1.0 },
       std::tuple<int, int, double>{ 1, 0, 1.0 },
       std::tuple<int, int, double>{ 1, 1, 1.0 },
       std::tuple<int, int, double>{ 1, 2, 1.0 },
       std::tuple<int, int, double>{ 2, 1, 1.0 },
       std::tuple<int, int, double>{ 2, 2, 1.0 },
       std::tuple<int, int, double>{ 2, 3, 1.0 },
       std::tuple<int, int, double>{ 2, 4, 1.0 },
       std::tuple<int, int, double>{ 3, 3, 1.0 },
       std::tuple<int, int, double>{ 3, 4, 1.0 },
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
   pb.setProblemName( "example for conflict analysis" );
   Problem<double> problem = pb.build();
   problem.getConstraintMatrix().modifyLeftHandSide( 0, {}, rhs[0] );
   problem.getConstraintMatrix().modifyLeftHandSide( 1, {}, rhs[1] );
   problem.getConstraintMatrix().modifyLeftHandSide( 2, {}, rhs[2] );
   problem.getConstraintMatrix().modifyLeftHandSide( 3, {}, rhs[3] );
   return problem;
}
