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
#include "catch/catch.hpp"
#include "papilo/core/Problem.hpp"
#include "papilo/core/ProblemBuilder.hpp"
#include "papilo/core/ProbingView.hpp"

Problem<double>
setupProblemForFixAndPropagation();

TEST_CASE( "fix-and-propagate", "[fix]" )
{
   Problem<double> problem = setupProblemForFixAndPropagation();

   problem.recomputeAllActivities();

   Vec<double> primal_solution;
   for( int i = 0; i < problem.getNCols(); i++ )
   {
      double random_number = ( 1.0 + i ) / 10.0;
      primal_solution.push_back( random_number );
   }
   Vec<double> res{ primal_solution };

   ProbingView<double> view{ problem, {}};
   FixAndPropagate<double> fixAndPropagate{ {}, {}, problem, view };
   bool success = fixAndPropagate.fix_and_propagate( primal_solution, res );

   assert( success );
   assert( res[0] == 1 );
   assert( res[1] == 1 );
   assert( res[2] == 0 );
}

Problem<double>
setupProblemForFixAndPropagation()
{
   Vec<double> coefficients{ 1.0, 1.0, 1.0 };
   Vec<double> upperBounds{ 1.0, 1.0, 1.0 };
   Vec<double> lowerBounds{ 0.0, 0.0, 0.0 };
   Vec<uint8_t> isIntegral{ 1, 1, 1 };

   Vec<double> rhs{ 2.0 };
   Vec<std::string> rowNames{ "A1" };
   Vec<std::string> columnNames{ "c1", "c2", "c3" };
   Vec<std::tuple<int, int, double>> entries{
       std::tuple<int, int, double>{ 0, 0, 1.0 },
       std::tuple<int, int, double>{ 0, 1, 1.0 },
       std::tuple<int, int, double>{ 0, 2, 1.0 },
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
