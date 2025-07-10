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

#include "papilo/presolvers/Probing.hpp"
#include "papilo/external/catch/catch_amalgamated.hpp"
#include "papilo/core/PresolveMethod.hpp"
#include "papilo/core/Problem.hpp"
#include "papilo/core/ProblemBuilder.hpp"

using namespace papilo;

Problem<double>
setupProblemWithProbing();

Problem<double>
setupProblemWithProbingWithNoBinary();

Problem<double>
setupProblemWithCliqueProbing();

Problem<double>
setupProblemWithCliqueProbingSubstitution();

TEST_CASE( "happy-path-probing", "[presolve]" )
{
   Num<double> num{};
   double time = 0.0;
   int cause = -1;
   Timer t{ time };
   Message msg{};
   Problem<double> problem = setupProblemWithProbing();
   Statistics statistics{};
   PresolveOptions presolveOptions{};
   presolveOptions.dualreds = 0;
   PostsolveStorage<double> postsolve =
       PostsolveStorage<double>( problem, num, presolveOptions );
   ProblemUpdate<double> problemUpdate( problem, postsolve, statistics,
                                        presolveOptions, num, msg );
   Probing<double> presolvingMethod{};
   Reductions<double> reductions{};
   problem.recomputeAllActivities();

   PresolveStatus presolveStatus =
       presolvingMethod.execute( problem, problemUpdate, num, reductions, t, cause );

   REQUIRE( presolveStatus == PresolveStatus::kReduced );
   REQUIRE( reductions.size() == 1 );
   REQUIRE( reductions.getReduction( 0 ).col == 1 );
   REQUIRE( reductions.getReduction( 0 ).row ==
            papilo::ColReduction::UPPER_BOUND );
   REQUIRE( reductions.getReduction( 0 ).newval == 0 );
}

TEST_CASE( "failed-path-probing-on-not-binary-variables", "[presolve]" )
{
   double time = 0.0;
   int cause = -1;
   Timer t{ time };
   Num<double> num{};
   Message msg{};
   Problem<double> problem = setupProblemWithProbingWithNoBinary();
   Statistics statistics{};
   PresolveOptions presolveOptions{};
   presolveOptions.dualreds = 0;
   PostsolveStorage<double> postsolve =
       PostsolveStorage<double>( problem, num, presolveOptions );
   ProblemUpdate<double> problemUpdate( problem, postsolve, statistics,
                                        presolveOptions, num, msg );
   Probing<double> presolvingMethod{};
   Reductions<double> reductions{};
   problem.recomputeAllActivities();

   PresolveStatus presolveStatus =
       presolvingMethod.execute( problem, problemUpdate, num, reductions, t, cause );

   REQUIRE( presolveStatus == PresolveStatus::kUnchanged );
}

TEST_CASE( "clique-probing-1", "[presolve]" )
{
   Num<double> num{};
   double time = 0.0;
   int cause = -1;
   Timer t{ time };
   Message msg{};
   Problem<double> problem = setupProblemWithCliqueProbing();
   Statistics statistics{};
   PresolveOptions presolveOptions{};
   presolveOptions.dualreds = 0;
   PostsolveStorage<double> postsolve =
       PostsolveStorage<double>( problem, num, presolveOptions );
   ProblemUpdate<double> problemUpdate( problem, postsolve, statistics,
                                        presolveOptions, num, msg );
   Probing<double> presolvingMethod{};
   Reductions<double> reductions{};
   problem.recomputeAllActivities();

   PresolveStatus presolveStatus =
       presolvingMethod.execute( problem, problemUpdate, num, reductions, t, cause );
    std::cout<<"\nTESTRESULTS:\n";
    for( int i = 0; i < static_cast<int>(reductions.size()); ++i )
    {
        std::cout<<"\n";
        std::cout<<reductions.getReduction(i).col;
        std::cout<<"\n";
        std::cout<<reductions.getReduction(i).row;
        std::cout<<"\n";
        std::cout<<reductions.getReduction(i).newval;
    }
    
   REQUIRE( presolveStatus == PresolveStatus::kReduced );/*
#ifdef PAPLIO_TBB
   REQUIRE( reductions.size() == 8 );

   REQUIRE( reductions.getReduction( 0 ).col == 2 );
   REQUIRE( reductions.getReduction( 0 ).row == -3 );
   REQUIRE( reductions.getReduction( 0 ).newval == 1 );

   REQUIRE( reductions.getReduction( 1 ).col == 0 );
   REQUIRE( reductions.getReduction( 1 ).row == -4 );
   REQUIRE( reductions.getReduction( 1 ).newval == 0 );

   REQUIRE( reductions.getReduction( 2 ).col == 1 );
   REQUIRE( reductions.getReduction( 2 ).row == -4 );
   REQUIRE( reductions.getReduction( 2 ).newval == 0 );

   REQUIRE( reductions.getReduction( 3 ).col == 3 );
   REQUIRE( reductions.getReduction( 3 ).row == -4 );
   REQUIRE( reductions.getReduction( 3 ).newval == 1 );

   REQUIRE( reductions.getReduction( 4 ).col == 4 );
   REQUIRE( reductions.getReduction( 4 ).row == -4 );
   REQUIRE( reductions.getReduction( 4 ).newval == 0 );

   REQUIRE( reductions.getReduction( 5 ).col == 3 );
   REQUIRE( reductions.getReduction( 5 ).row == -4 );
   REQUIRE( reductions.getReduction( 5 ).newval == 1 );

   REQUIRE( reductions.getReduction( 6 ).col == 4 );
   REQUIRE( reductions.getReduction( 6 ).row == -10 );
   REQUIRE( reductions.getReduction( 6 ).newval == 1 );

   REQUIRE( reductions.getReduction( 7 ).col == 1 );
   REQUIRE( reductions.getReduction( 7 ).row == -1 );
   REQUIRE( reductions.getReduction( 7 ).newval == 0 );

#else

REQUIRE( reductions.size() == 8 );

REQUIRE( reductions.getReduction( 0 ).col == 2 );
REQUIRE( reductions.getReduction( 0 ).row == -3 );
REQUIRE( reductions.getReduction( 0 ).newval == 1 );

REQUIRE( reductions.getReduction( 1 ).col == 0 );
REQUIRE( reductions.getReduction( 1 ).row == -4 );
REQUIRE( reductions.getReduction( 1 ).newval == 0 );

REQUIRE( reductions.getReduction( 2 ).col == 1 );
REQUIRE( reductions.getReduction( 2 ).row == -4 );
REQUIRE( reductions.getReduction( 2 ).newval == 0 );

REQUIRE( reductions.getReduction( 3 ).col == 3 );
REQUIRE( reductions.getReduction( 3 ).row == -4 );
REQUIRE( reductions.getReduction( 3 ).newval == 1 );

REQUIRE( reductions.getReduction( 4 ).col == 4 );
REQUIRE( reductions.getReduction( 4 ).row == -4 );
REQUIRE( reductions.getReduction( 4 ).newval == 0 );

REQUIRE( reductions.getReduction( 5 ).col == 3 );
REQUIRE( reductions.getReduction( 5 ).row == -4 );
REQUIRE( reductions.getReduction( 5 ).newval == 1 );

REQUIRE( reductions.getReduction( 6 ).col == 4 );
REQUIRE( reductions.getReduction( 6 ).row == -10 );
REQUIRE( reductions.getReduction( 6 ).newval == 1 );

REQUIRE( reductions.getReduction( 7 ).col == 1 );
REQUIRE( reductions.getReduction( 7 ).row == -1 );
REQUIRE( reductions.getReduction( 7 ).newval == 0 );*/
/*
   REQUIRE( reductions.size() == 3 );

   REQUIRE( reductions.getReduction( 0 ).col == 3 );
   REQUIRE( reductions.getReduction( 0 ).row == papilo::ColReduction::UPPER_BOUND );
   REQUIRE( reductions.getReduction( 0 ).newval == 1 );

   REQUIRE( reductions.getReduction( 1 ).col == 4 );
   REQUIRE( reductions.getReduction( 1 ).row == papilo::ColReduction::REPLACE );
   REQUIRE( reductions.getReduction( 1 ).newval == 1 );

   REQUIRE( reductions.getReduction( 2 ).col == 1 );
   REQUIRE( reductions.getReduction( 2 ).row == papilo::ColReduction::NONE );
   REQUIRE( reductions.getReduction( 2 ).newval == 0 );*/

#endif
}


TEST_CASE( "clique-probing-2", "[presolve]" )
{
   Num<double> num{};
   double time = 0.0;
   int cause = -1;
   Timer t{ time };
   Message msg{};
   Problem<double> problem = setupProblemWithCliqueProbingSubstitution();
   Statistics statistics{};
   PresolveOptions presolveOptions{};
   presolveOptions.dualreds = 0;
   PostsolveStorage<double> postsolve =
       PostsolveStorage<double>( problem, num, presolveOptions );
   ProblemUpdate<double> problemUpdate( problem, postsolve, statistics,
                                        presolveOptions, num, msg );
   Probing<double> presolvingMethod{};
   Reductions<double> reductions{};
   problem.recomputeAllActivities();

   PresolveStatus presolveStatus =
       presolvingMethod.execute( problem, problemUpdate, num, reductions, t, cause );
    std::cout<<"\nTESTRESULTS:\n";
    for( int i = 0; i < static_cast<int>(reductions.size()); ++i )
    {
        std::cout<<"\n";
        std::cout<<reductions.getReduction(i).col;
        std::cout<<"\n";
        std::cout<<reductions.getReduction(i).row;
        std::cout<<"\n";
        std::cout<<reductions.getReduction(i).newval;
    }

   REQUIRE( presolveStatus == PresolveStatus::kReduced );
   REQUIRE( reductions.size() == 1 );

   REQUIRE( reductions.getReduction( 0 ).col == 9 );
   REQUIRE( reductions.getReduction( 0 ).row == papilo::ColReduction::REPLACE );
   REQUIRE( reductions.getReduction( 0 ).newval == 1 );

   REQUIRE( reductions.getReduction( 1 ).col == 2 );
   REQUIRE( reductions.getReduction( 1 ).row == papilo::ColReduction::NONE );
   REQUIRE( reductions.getReduction( 1 ).newval == 0 );
   
   REQUIRE( reductions.getReduction( 2 ).col == 7 );
   REQUIRE( reductions.getReduction( 2 ).row == papilo::ColReduction::REPLACE );
   REQUIRE( reductions.getReduction( 2 ).newval == 1 );

   REQUIRE( reductions.getReduction( 3 ).col == 3 );
   REQUIRE( reductions.getReduction( 3 ).row == papilo::ColReduction::NONE );
   REQUIRE( reductions.getReduction( 3 ).newval == 0 );
   
   REQUIRE( reductions.getReduction( 4 ).col == 5 );
   REQUIRE( reductions.getReduction( 4 ).row == papilo::ColReduction::REPLACE );
   REQUIRE( reductions.getReduction( 4 ).newval == 1 );

   REQUIRE( reductions.getReduction( 5 ).col == 0 );
   REQUIRE( reductions.getReduction( 5 ).row == papilo::ColReduction::NONE );
   REQUIRE( reductions.getReduction( 5 ).newval == 0 );

}

Problem<double>
setupProblemWithProbing()
{
   // x + 2y<=1
   // y + z + w <=1
   Vec<double> coefficients{ 1.0, 1.0, 1.0, 1.0 };
   Vec<double> upperBounds{ 1.0, 1.0, 1.0, 1.0 };
   Vec<double> lowerBounds{ 0.0, 0.0, 0.0, 0.0 };
   Vec<uint8_t> isIntegral{ 1, 1, 1, 1 };

   Vec<double> rhs{ 1.0 };
   Vec<std::string> rowNames{ "A1" };
   Vec<std::string> columnNames{ "c1", "c2", "c3", "c4" };
   Vec<std::tuple<int, int, double>> entries{
       std::tuple<int, int, double>{ 0, 0, 1.0 },
       std::tuple<int, int, double>{ 0, 1, 2.0 },
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
   pb.addEntryAll( entries );
   pb.setColNameAll( columnNames );
   pb.setProblemName( "matrix for testing probing" );
   Problem<double> problem = pb.build();
   return problem;
}

Problem<double>
setupProblemWithProbingWithNoBinary()
{
   // 2x + y<=1
   // y + z + w <=1
   Vec<double> coefficients{ 1.0, 1.0, 1.0, 1.0 };
   Vec<double> upperBounds{ 4.0, 3.0, 2.0, 3.0 };
   Vec<double> lowerBounds{ 0.0, 0.0, 0.0, 0.0 };
   Vec<uint8_t> isIntegral{ 1, 1, 1, 1 };

   Vec<double> rhs{ 3.0, 3.0 };
   Vec<std::string> rowNames{ "A1", "A2" };
   Vec<std::string> columnNames{ "c1", "c2", "c3", "c4" };
   Vec<std::tuple<int, int, double>> entries{
       std::tuple<int, int, double>{ 0, 0, 1.0 },
       std::tuple<int, int, double>{ 0, 1, 1.0 },
       std::tuple<int, int, double>{ 1, 1, 1.0 },
       std::tuple<int, int, double>{ 1, 2, 1.0 } };

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
   pb.addEntryAll( entries );
   pb.setColNameAll( columnNames );
   pb.setProblemName( "matrix for testing probing no binaries" );
   Problem<double> problem = pb.build();
   return problem;
}

Problem<double>
setupProblemWithCliqueProbing()
{
   // x + y + z <= 1
   // x + w <= 2
   // -x + w <= 1
   // y - v = 0
   Vec<double> coefficients{ 1.0, 1.0, 1.0, 1.0, 1.0 };
   Vec<double> upperBounds{ 1.0, 1.0, 1.0, 2.0, 1.0 };
   Vec<double> lowerBounds{ 0.0, 0.0, 0.0, 0.0, 0.0 };
   Vec<uint8_t> isIntegral{ 1, 1, 1, 1, 1 };

   Vec<double> rhs{ 1.0, 2.0, 1.0, 0.0 };
   Vec<double> lhs{ 0.0, 0.0, 0.0, 0.0 };
   Vec<std::string> rowNames{ "A1", "A2", "A3", "A4" };
   Vec<std::string> columnNames{ "c1", "c2", "c3", "c4", "c5" };
   Vec<std::tuple<int, int, double>> entries{
       std::tuple<int, int, double>{ 0, 0, 1.0 },
       std::tuple<int, int, double>{ 0, 1, 1.0 },
       std::tuple<int, int, double>{ 0, 2, 1.0 },

       std::tuple<int, int, double>{ 1, 0, 1.0 },
       std::tuple<int, int, double>{ 1, 3, 1.0 },

       std::tuple<int, int, double>{ 2, 0, -1.0 },
       std::tuple<int, int, double>{ 2, 3, 1.0 },

       std::tuple<int, int, double>{ 3, 1, 1.0 },
       std::tuple<int, int, double>{ 3, 4, -1.0 },
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
   pb.setRowLhsAll( lhs );
   pb.addEntryAll( entries );
   pb.setColNameAll( columnNames );
   pb.setProblemName( "matrix for testing clique probing" );
   Problem<double> problem = pb.build();
   return problem;
}

Problem<double>
setupProblemWithCliqueProbingSubstitution()
{
   // x1 + x2 + x3 + x4 + x5 <= 1
   // x1 - x6 = 0
   // x2 + x7 = 1
   // x3 - x10 = 0
   // x4 - x8 = 0
   // x5 + x9 = 1
   // x1 - x2 + x3 - x4 + x5 - x11 == 0
   // x1 - x2 + x3 - x4 + x5 - x12 == 0
   // x1 - x2 + x3 - x4 + x5 - x13 == 0
   // x1 - x2 + x3 - x4 + x5 - x14 == 0
   // x1 - x2 + x3 - x4 + x5 - x15 == 0
   // x1 - x2 + x3 - x4 + x5 - x16 == 0
   
   Vec<double> coefficients{ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };
   Vec<double> upperBounds{ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };
   Vec<double> lowerBounds{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
   Vec<uint8_t> isIntegral{ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };

   Vec<double> rhs{ 1.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
   Vec<double> lhs{ 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
   Vec<std::string> rowNames{ "A1", "A2", "A3", "A4", "A5", "A6", "A7", "A8", "A9", "A10", "A11", "A12" };
   Vec<std::string> columnNames{ "c1", "c2", "c3", "c4", "c5", "c6", "c7", "c8", "c9", "c10", "c11", "c12", "c13", "c14", "c15", "c16" };
   Vec<std::tuple<int, int, double>> entries{
       std::tuple<int, int, double>{ 0, 0, 1.0 },
       std::tuple<int, int, double>{ 0, 1, 1.0 },
       std::tuple<int, int, double>{ 0, 2, 1.0 },
       std::tuple<int, int, double>{ 0, 3, 1.0 },
       std::tuple<int, int, double>{ 0, 4, 1.0 },

       std::tuple<int, int, double>{ 1, 0, 1.0 },
       std::tuple<int, int, double>{ 1, 5, -1.0 },

       std::tuple<int, int, double>{ 2, 1, 1.0 },
       std::tuple<int, int, double>{ 2, 6, 1.0 },

       std::tuple<int, int, double>{ 3, 2, 1.0 },
       std::tuple<int, int, double>{ 3, 9, -1.0 },

       std::tuple<int, int, double>{ 4, 3, 1.0 },
       std::tuple<int, int, double>{ 4, 7, -1.0 },

       std::tuple<int, int, double>{ 5, 4, 1.0 },
       std::tuple<int, int, double>{ 5, 8, 1.0 },

       std::tuple<int, int, double>{ 6, 0, 1.0 },
       std::tuple<int, int, double>{ 6, 1, -1.0 },
       std::tuple<int, int, double>{ 6, 2, 1.0 },
       std::tuple<int, int, double>{ 6, 3, -1.0 },
       std::tuple<int, int, double>{ 6, 4, 1.0 },
       std::tuple<int, int, double>{ 6, 10, -1.0 },
       
       std::tuple<int, int, double>{ 7, 0, 1.0 },
       std::tuple<int, int, double>{ 7, 1, -1.0 },
       std::tuple<int, int, double>{ 7, 2, 1.0 },
       std::tuple<int, int, double>{ 7, 3, -1.0 },
       std::tuple<int, int, double>{ 7, 4, 1.0 },
       std::tuple<int, int, double>{ 7, 11, -1.0 },
       
       std::tuple<int, int, double>{ 8, 0, 1.0 },
       std::tuple<int, int, double>{ 8, 1, -1.0 },
       std::tuple<int, int, double>{ 8, 2, 1.0 },
       std::tuple<int, int, double>{ 8, 3, -1.0 },
       std::tuple<int, int, double>{ 8, 4, 1.0 },
       std::tuple<int, int, double>{ 8, 12, -1.0 },
       
       std::tuple<int, int, double>{ 9, 0, 1.0 },
       std::tuple<int, int, double>{ 9, 1, -1.0 },
       std::tuple<int, int, double>{ 9, 2, 1.0 },
       std::tuple<int, int, double>{ 9, 3, -1.0 },
       std::tuple<int, int, double>{ 9, 4, 1.0 },
       std::tuple<int, int, double>{ 9, 13, -1.0 },
       
       std::tuple<int, int, double>{ 10, 0, 1.0 },
       std::tuple<int, int, double>{ 10, 1, -1.0 },
       std::tuple<int, int, double>{ 10, 2, 1.0 },
       std::tuple<int, int, double>{ 10, 3, -1.0 },
       std::tuple<int, int, double>{ 10, 4, 1.0 },
       std::tuple<int, int, double>{ 10, 14, -1.0 },
       
       std::tuple<int, int, double>{ 11, 0, 1.0 },
       std::tuple<int, int, double>{ 11, 1, -1.0 },
       std::tuple<int, int, double>{ 11, 2, 1.0 },
       std::tuple<int, int, double>{ 11, 3, -1.0 },
       std::tuple<int, int, double>{ 11, 4, 1.0 },
       std::tuple<int, int, double>{ 11, 15, -1.0 },

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
   pb.setRowLhsAll( lhs );
   pb.addEntryAll( entries );
   pb.setColNameAll( columnNames );
   pb.setProblemName( "matrix for testing clique probing substitutions" );
   Problem<double> problem = pb.build();
   return problem;
}