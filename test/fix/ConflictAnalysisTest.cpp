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
#include "fix/ConflictAnalysis.hpp"
#include "fix/FixAndPropagate.hpp"

#include "catch/catch.hpp"
#include "papilo/core/ProbingView.hpp"
#include "papilo/core/Problem.hpp"
#include "papilo/core/ProblemBuilder.hpp"

Problem<double>
setupProblemForConflictAnalysis();

Problem<double>
setupProblemForConflictAnalysisNr2();

TEST_CASE( "conflict-analysis-binary-no-resolution", "[conflict]" )
{
   Problem<double> problem = setupProblemForConflictAnalysis();
   problem.recomputeAllActivities();
   // Binary problem with constraints
   // A1: x1 + x3 <= 1
   // A2: x1 + x2 + x3 >= 2
   // A3: x2 + x3 + x4 + x5 <= 3
   // A4: x4 + x5 = 2

   // Assume that fix_and_propagate does the following:
   // Fix: x3 = 1 (decision level 1)
   // propagate A1: x1 = 0 (reason row 0, decision level 1)
   // propagate A2: x2 = 1 (reason row 1, decision level 1)
   // Fix: x4 = 1 (decision level 2)
   // propagate A3: x5 = 0 (reason row 2, decision level 2)
   // propagating A4 leads to conflict -> apply conflict analysis

   Vec<SingleBoundChange<double>> bound_changes;
   // bound changes for x1, x2, x3 (decision level 1)
   SingleBoundChange<double> bound_change_1( 2, -1, 1.0, true, true, -2 );
   bound_changes.push_back( bound_change_1 );
   SingleBoundChange<double> bound_change_2( 0, 0, 0.0, false, false, -3 );
   bound_changes.push_back( bound_change_2 );
   SingleBoundChange<double> bound_change_3( 1, 1, 1.0, false, true, -3 );
   bound_changes.push_back( bound_change_3 );
   // bound changes for x4, x5 (decision level 2)
   SingleBoundChange<double> bound_change_4( 3, -1, 1.0, true, true, -4 );
   bound_changes.push_back( bound_change_4 );
   SingleBoundChange<double> bound_change_5( 4, 2, 0.0, false, false, -5 );
   bound_changes.push_back( bound_change_5 );

   double t = 0;
   Timer timer{ t };

   // initialize conflict analysis
   ConflictAnalysis<double> conflictAnalysis( {}, {}, timer, problem );

   papilo::Vec<std::pair<int, int>> infeasible_rows;
   std::pair<int, int> infeas_pair;
   infeas_pair.first = 5;
   infeas_pair.second = 3;
   infeasible_rows.push_back( infeas_pair );
   Vec<Constraint<double>> conflict_constraints;
   conflictAnalysis.perform_conflict_analysis( bound_changes, infeasible_rows,
                                               conflict_constraints );

   REQUIRE( conflict_constraints.size() == 2 );

   // conflict analysis returns the constraint x5 >= 1 (1-FUIP)
   SparseVectorView<double> data_c1 = conflict_constraints[0].get_data();
   const double* vals_c1 = data_c1.getValues();
   const int* inds_c1 = data_c1.getIndices();
   int len_c1 = data_c1.getLength();

   REQUIRE( len_c1 == 1 );
   REQUIRE( vals_c1[0] == 1.0 );
   REQUIRE( inds_c1[0] == 4 );

   // conflict analysis returns the constraint x5 >= 1 (All-FUIP)
   SparseVectorView<double> data_c2 = conflict_constraints[1].get_data();
   const double* vals_c2 = data_c2.getValues();
   const int* inds_c2 = data_c2.getIndices();
   int len_c2 = data_c2.getLength();

   REQUIRE( len_c2 == 1 );
   REQUIRE( vals_c2[0] == 1.0 );
   REQUIRE( inds_c2[0] == 4 );
}

TEST_CASE( "conflict-analysis-binary-with-resolution", "[conflict]" )
{
   Problem<double> problem = setupProblemForConflictAnalysisNr2();
   problem.recomputeAllActivities();
   // Binary problem with constraints
   // A1: x1 + x2 +x3 <= 1
   // A2: x1 + x3 + x4 + x5 <= 2
   // A3: 1 <= x3 + x5 + x6
   // A4: -x2 -x5 + x6 <= 0

   // Assume that fix_and_propagate does the following:
   // Fix: x1 = 1
   // Prop A1: x2 = 0, x3 = 0
   // Fix: x4 = 1
   // Prop A2: x5 = 0
   // Prop A3: x6 = 1
   // Conflict at A4

   Vec<SingleBoundChange<double>> bound_changes;
   // bound changes for x1, x2, x3 (decision level 1)
   SingleBoundChange<double> bound_change_1( 0, -1, 1.0, true, true, -2 );
   bound_changes.push_back( bound_change_1 );
   SingleBoundChange<double> bound_change_2( 1, 0, 0.0, false, false, -3 );
   bound_changes.push_back( bound_change_2 );
   SingleBoundChange<double> bound_change_3( 2, 0, 0.0, false, false, -3 );
   bound_changes.push_back( bound_change_3 );
   // bound changes for x4, x5, x6 (decision level 2)
   SingleBoundChange<double> bound_change_4( 3, -1, 1.0, true, true, -4 );
   bound_changes.push_back( bound_change_4 );
   SingleBoundChange<double> bound_change_5( 4, 1, 0.0, false, false, -5 );
   bound_changes.push_back( bound_change_5 );
   SingleBoundChange<double> bound_change_6( 5, 2, 1.0, false, true, -5 );
   bound_changes.push_back( bound_change_6 );

   double t = 0;
   Timer timer{ t };

   // initialize conflict analysis
   ConflictAnalysis<double> conflictAnalysis( {}, {}, timer, problem );

   papilo::Vec<std::pair<int, int>> infeasible_rows;
   std::pair<int, int> infeas_pair;
   infeas_pair.first = 6;
   infeas_pair.second = 3;
   infeasible_rows.push_back( infeas_pair );
   Vec<Constraint<double>> conflict_constraints;
   conflictAnalysis.perform_conflict_analysis( bound_changes, infeasible_rows,
                                               conflict_constraints );

   REQUIRE( conflict_constraints.size() == 2 );

   // conflict analysis returns the constraint x2 + x3 + x5 >= 1 (1-FUIP)
   SparseVectorView<double> data_c1 = conflict_constraints[0].get_data();
   const double* vals_c1 = data_c1.getValues();
   const int* inds_c1 = data_c1.getIndices();
   int len_c1 = data_c1.getLength();

   REQUIRE( len_c1 == 3 );
   REQUIRE( vals_c1[0] == 1.0 );
   REQUIRE( vals_c1[1] == 1.0 );
   REQUIRE( vals_c1[2] == 1.0 );
   REQUIRE( inds_c1[0] == 1 );
   REQUIRE( inds_c1[1] == 2 );
   REQUIRE( inds_c1[2] == 4 );

   // ToDo All-FUIP
}

Problem<double>
setupProblemForConflictAnalysis()
{
   Vec<double> coefficients{ 1.0, 1.0, 1.0, 1.0, 1.0 };
   Vec<double> upperBounds{ 1.0, 1.0, 1.0, 1.0, 1.0 };
   Vec<double> lowerBounds{ 0.0, 0.0, 0.0, 0.0, 0.0 };
   Vec<uint8_t> isIntegral{ 1, 1, 1, 1, 1 };

   Vec<uint8_t> lhsInf{ 1, 0, 1, 0 };
   Vec<uint8_t> rhsInf{ 0, 1, 0, 0 };
   Vec<double> lhs{ 0.0, 2.0, 0.0, 2.0 };
   Vec<double> rhs{ 1.0, 0.0, 3.0, 2.0 };
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
   pb.setRowLhsInfAll( lhsInf );
   pb.setRowRhsInfAll( rhsInf );
   pb.setRowLhsAll( lhs );
   pb.setRowRhsAll( rhs );
   pb.addEntryAll( entries );
   pb.setColNameAll( columnNames );
   pb.setProblemName( "example for conflict analysis" );
   Problem<double> problem = pb.build();
   return problem;
}

Problem<double>
setupProblemForConflictAnalysisNr2()
{
   Vec<double> coefficients{ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };
   Vec<double> upperBounds{ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };
   Vec<double> lowerBounds{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
   Vec<uint8_t> isIntegral{ 1, 1, 1, 1, 1, 1 };

   Vec<uint8_t> lhsInf{ 1, 1, 0, 1 };
   Vec<uint8_t> rhsInf{ 0, 0, 1, 0 };
   Vec<double> lhs{ 0.0, 0.0, 1.0, 0.0 };
   Vec<double> rhs{ 1.0, 2.0, 0.0, 0.0 };
   Vec<std::string> rowNames{ "A1", "A2", "A3", "A4" };
   Vec<std::string> columnNames{ "c1", "c2", "c3", "c4", "c5", "c6" };
   Vec<std::tuple<int, int, double>> entries{
       std::tuple<int, int, double>{ 0, 0, 1.0 },
       std::tuple<int, int, double>{ 0, 1, 1.0 },
       std::tuple<int, int, double>{ 0, 2, 1.0 },

       std::tuple<int, int, double>{ 1, 0, 1.0 },
       std::tuple<int, int, double>{ 1, 2, 1.0 },
       std::tuple<int, int, double>{ 1, 3, 1.0 },
       std::tuple<int, int, double>{ 1, 4, 1.0 },

       std::tuple<int, int, double>{ 2, 2, 1.0 },
       std::tuple<int, int, double>{ 2, 4, 1.0 },
       std::tuple<int, int, double>{ 2, 5, 1.0 },

       std::tuple<int, int, double>{ 3, 1, -1.0 },
       std::tuple<int, int, double>{ 3, 4, -1.0 },
       std::tuple<int, int, double>{ 3, 5, 1.0 },
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
   pb.setRowLhsInfAll( lhsInf );
   pb.setRowRhsInfAll( rhsInf );
   pb.setRowLhsAll( lhs );
   pb.setRowRhsAll( rhs );
   pb.addEntryAll( entries );
   pb.setColNameAll( columnNames );
   pb.setProblemName( "example 2 for conflict analysis" );
   Problem<double> problem = pb.build();
   return problem;
}

Problem<double>
setupSimpleProblem()
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
