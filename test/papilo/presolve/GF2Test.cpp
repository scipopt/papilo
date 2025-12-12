/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*               This file is part of the program and library                */
/*    PaPILO --- Parallel Presolve for Integer and Linear Optimization       */
/*                                                                           */
/* Copyright (C) 2020-2025 Zuse Institute Berlin (ZIB)                       */
/*                                                                           */
/* Licensed under the Apache License, Version 2.0 (the "License");           */
/* you may not use this file except in compliance with the License.          */
/* You may obtain a copy of the License at                                   */
/*                                                                           */
/*     http://www.apache.org/licenses/LICENSE-2.0                            */
/*                                                                           */
/* Unless required by applicable law or agreed to in writing, software       */
/* distributed under the License is distributed on an "AS IS" BASIS,         */
/* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  */
/* See the License for the specific language governing permissions and       */
/* limitations under the License.                                            */
/*                                                                           */
/* You should have received a copy of the Apache-2.0 license                 */
/* along with PaPILO; see the file LICENSE. If not visit scipopt.org.        */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include "papilo/core/PresolveMethod.hpp"
#include "papilo/core/Problem.hpp"
#include "papilo/core/ProblemBuilder.hpp"
#include "papilo/presolvers/GF2.hpp"
#include "papilo/external/catch/catch_amalgamated.hpp"

using namespace papilo;

Problem<double>
setupProblemForGF2();

TEST_CASE( "happy-path-gf2", "[presolve]" )
{
   double time = 0.0;
   int cause = -1;
   Timer t{time};
   Num<double> num{};
   Message msg{};
   Problem<double> problem = setupProblemForGF2();
   Statistics statistics{};
   PresolveOptions presolveOptions{};
   PostsolveStorage<double> postsolve =
       PostsolveStorage<double>( problem, num, presolveOptions );
   ProblemUpdate<double> problemUpdate( problem, postsolve, statistics,
                                        presolveOptions, num, msg);
   GF2<double> presolvingMethod{};
   Reductions<double> reductions{};

   // PresolveStatus presolveStatus =
       presolvingMethod.execute( problem, problemUpdate, num, reductions, t, cause);

   //TODO:
}

Problem<double>
setupProblemForGF2()
{
   // x + y + z <= 1
   // x + y     <= 1
   // ⇒ z = 0 (GF(2) elimination)
   Vec<double> obj{ 0.0, 0.0, 0.0 };
   Vec<double> ub{ 1.0, 1.0, 1.0 };
   Vec<double> lb{ 0.0, 0.0, 0.0 };
   Vec<uint8_t> isIntegral{ 1, 1, 1 };

   Vec<double> rhs{ 1.0, 1.0 };
   Vec<double> lhs{ 0.0, 0.0 };
   Vec<std::string> rowNames{ "A1", "A2" };
   Vec<std::string> colNames{ "x", "y", "z" };

   Vec<std::tuple<int, int, double>> entries{
          { 0, 0, 1.0 }, // x
          { 0, 1, 1.0 }, // y
          { 0, 2, 1.0 }, // z

          { 1, 0, 1.0 }, // x
          { 1, 1, 1.0 }  // y
   };

   ProblemBuilder<double> pb;
   pb.reserve( entries.size(), rowNames.size(), colNames.size() );
   pb.setNumRows( rowNames.size() );
   pb.setNumCols( colNames.size() );
   pb.setColUbAll( ub );
   pb.setColLbAll( lb );
   pb.setObjAll( obj );
   pb.setObjOffset( 0.0 );
   pb.setColIntegralAll( isIntegral );
   pb.setRowRhsAll( rhs );
   pb.setRowLhsAll( lhs );
   pb.addEntryAll( entries );
   pb.setColNameAll( colNames );
   pb.setProblemName( "gf2 presolve test" );

   return pb.build();
}
