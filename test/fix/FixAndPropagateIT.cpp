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
#include "papilo/io/SolParser.hpp"

#include "catch/catch.hpp"
#include "fix/strategy/FarkasRoundingStrategy.hpp"
#include "fix/strategy/FractionalRoundingStrategy.hpp"
#include "fix/strategy/RandomRoundingStrategy.hpp"
#include "papilo/core/ProbingView.hpp"
#include "papilo/core/Problem.hpp"

Problem<double>
read_prob( const std::string& filename );

Solution<double>
parse_solution( const std::string& path_to_sol, const Problem<double>& problem );

TEST_CASE( "fix-and-propagate-it-solution-is-feasible", "[fix]" )
{
   std::string path_to_mps = "/home/alexander/git_repositories/mipcomp22/MIP "
                             "2022 Open Instances/rococoC10-001000.mps.gz";
   std::string path_to_sol = "/home/alexander/Downloads/solutions/"
                             "rococoC10-001000/1/rococoC10-001000.sol.gz";
   auto problem = read_prob( path_to_mps );
   problem.recomputeAllActivities();
   Solution<double> solution = parse_solution( path_to_sol, problem );
   FixAndPropagate<double> fixAndPropagate{
       {}, {}, problem, { problem, {} }, true };
   FarkasRoundingStrategy<double> strategy{ 0, {} };
   Vec<double> res{ solution.primal };

   bool infeasible =
       fixAndPropagate.fix_and_propagate( solution.primal, res, strategy );

   assert( !infeasible );
   for( int i = 0; i < problem.getNCols(); i++ )
      assert( solution.primal[i] == res[i] );
}

TEST_CASE( "fix-and-propagate-it-modified-solution-is-feasible", "[fix]" )
{
   std::string path_to_mps = "/home/alexander/git_repositories/mipcomp22/MIP "
                             "2022 Open Instances/rococoC10-001000.mps.gz";
   std::string path_to_sol = "/home/alexander/Downloads/solutions/"
                             "rococoC10-001000/1/rococoC10-001000.sol.gz";
   auto problem = read_prob( path_to_mps );
   problem.recomputeAllActivities();
   Solution<double> solution = parse_solution( path_to_sol, problem );
   FixAndPropagate<double> fixAndPropagate{
       {}, {}, problem, { problem, {} }, true };
   FarkasRoundingStrategy<double> strategy{ 0, {} };
   Vec<double> modified_primal{ solution.primal };
   modified_primal[1] = 0.6;

   Vec<double> res{ modified_primal };
   bool infeasible =
       fixAndPropagate.fix_and_propagate( modified_primal, res, strategy );

   assert( !infeasible );
   for( int i = 0; i < problem.getNCols(); i++ )
      assert( solution.primal[i] == res[i] );
}

Solution<double>
parse_solution( const std::string& path_to_sol, const Problem<double>& problem )
{
   Solution<double> solution;

   Vec<int> mapping{};
   for( int i = 0; i < problem.getNCols(); i++ )
      mapping.push_back( i );

   bool success = SolParser<double>::read(
       path_to_sol, mapping, problem.getVariableNames(), solution );
   assert( success );
   return solution;
}

Problem<double>
read_prob( const std::string& filename )
{
   boost::optional<Problem<double>> prob;
   {
      prob = MpsParser<double>::loadProblem( filename );
   }
   if( !prob )
   {
      fmt::print( "error loading problem {}\n", filename );
      assert( false );
   }
   return prob.get();
}