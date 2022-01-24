/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*               This file is part of the program and library                */
/*    PaPILO --- Parallel Presolve for Integer and Linear Optimization       */
/*                                                                           */
/* Copyright (C) 2020-2022 Konrad-Zuse-Zentrum                               */
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
#include "fix/VolumeAlgorithm.hpp"
#include "papilo/core/Problem.hpp"
#include "papilo/io/MpsParser.hpp"
#include "papilo/misc/OptionsParser.hpp"
#include <boost/program_options.hpp>
#include <fstream>

using namespace papilo;

Solution<double>
generate_random_solution( const Problem<double>& problem );

int
main( int argc, char* argv[] )
{

   // get the options passed by the user
   OptionsInfo optionsInfo;
   try
   {
      optionsInfo = parseOptions( argc, argv );
   }
   catch( const boost::program_options::error& ex )
   {
      std::cerr << "Error while parsing the options.\n" << '\n';
      std::cerr << ex.what() << '\n';
      return 1;
   }

   if( !optionsInfo.is_complete )
      return 0;

   double readtime = 0;
   Problem<double> problem;
   boost::optional<Problem<double>> prob;

   {
      Timer t( readtime );
      prob = MpsParser<double>::loadProblem( optionsInfo.instance_file );
   }

   // Check whether reading was successful or not
   if( !prob )
   {
      fmt::print( "error loading problem {}\n", optionsInfo.instance_file );
      return 0;
   }
   problem = *prob;

   fmt::print( "reading took {:.3} seconds\n", readtime );

   // set up ProblemUpdate to trivialPresolve so that activities exist
   Presolve<double> presolve{};
   presolve.apply( problem, false );

   VolumeAlgorithm<double> algorithm{ {}, {}, 0.5, 0.1, 1, 0.0005, 2, 1.1, 0.66,
      0.02, 0.01, 2, 20 };

   // generate pi
   Vec<double> pi{};
   for( int i = 0; i < problem.getNRows(); i++ )
   {
      pi.push_back( 0 );
   }

   // generate UB
   StableSum<double> min_value{};
   Num<double> num{};
   for( int i = 0; i < problem.getNCols(); i++ )
   {
      if( num.isZero( problem.getObjective().coefficients[i] ) )
         continue;
      else if( num.isLT( problem.getObjective().coefficients[i], 0 ) )
      {
         if( problem.getColFlags()[i].test( ColFlag::kLbInf ) )
         {
            fmt::print(
                "Could not calculate objective bound: variable {} is unbounded",
                i );
            return 1;
         }
         min_value.add( problem.getObjective().coefficients[i] +
                        problem.getLowerBounds()[i] );
      }
      else
      {
         if( problem.getColFlags()[i].test( ColFlag::kUbInf ) )
         {
            fmt::print(
                "Could not calculate objective bound: variable {} is unbounded",
                i );
            return 1;
         }
         min_value.add( problem.getObjective().coefficients[i] +
                        problem.getUpperBounds()[i] );
      }
   }
   // TODO: Suresh we need to discuss how we treat the inequalities
   // currently all constraints are considered as equalities
   algorithm.volume_algorithm(
       problem.getObjective().coefficients, problem.getConstraintMatrix(),
       problem.getConstraintMatrix().getRightHandSides(),
       problem.getVariableDomains(), pi, min_value.get() );

   //   Solution<double> random_solution = generate_random_solution( problem );
   //
   //   ProbingView<double> probing_view{ problem, num };
   //   FixAndPropagate<double> fixAndPropagate{ msg, num };
   //   fixAndPropagate.fix_and_propagate( probUpdate.getProblem(),
   //                                      probing_view, random_solution );

   return 0;
}

Solution<double>
generate_random_solution( const Problem<double>& problem )
{
   //   std::random_device dev;
   //   std::mt19937 rng( dev() );

   Vec<double> solution;
   for( int i = 0; i < problem.getNCols(); i++ )
   {
      double random_number = ( 1.0 + i ) / 10.0;
      solution.push_back( random_number );
   }
   return { SolutionType::kPrimal, solution };
}
