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

#include "fix/Algorithm.hpp"
#include "papilo/core/Problem.hpp"
#include "papilo/io/MpsParser.hpp"
#include "papilo/misc/OptionsParser.hpp"
#include <boost/program_options.hpp>
#include <fstream>

using namespace papilo;

/***
 * Missing parts:
 * Overall:
 * - Statistics
 *  - first solution found, optimal solution found
 *  - optimal value of solution/ optimum
 *  - test on MipComp instances
 *  (- profiling)
 *
 * Volume algorithm:
 * - restarting possibility (use vector pi should work already)
 * - add constraints of conflict analysis (extend pi)
 * - fractional of integer variables
 *
 * Conflict Analysis:
 * - TBD
 *
 * Fix and Propagate:
 * - optimize parallel scheme of OneOp
 */


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

//   if( !optionsInfo.param_settings_file.empty() )
//   {
//      ParameterSet paramSet = presolve.getParameters();
//
//      if( !optionsInfo.param_settings_file.empty() && !opts.print_params )
//      {
//         std::ifstream input( opts.param_settings_file );
//         if( input )
//         {
//            String theoptionstr;
//            String thevaluestr;
//            for( String line; getline( input, line ); )
//            {
//               std::size_t pos = line.find_first_of( '#' );
//               if( pos != String::npos )
//                  line = line.substr( 0, pos );
//
//               pos = line.find_first_of( '=' );
//
//               if( pos == String::npos )
//                  continue;
//
//               theoptionstr = line.substr( 0, pos - 1 );
//               thevaluestr = line.substr( pos + 1 );
//
//               boost::algorithm::trim( theoptionstr );
//               boost::algorithm::trim( thevaluestr );
//
//               try
//               {
//                  paramSet.parseParameter( theoptionstr.c_str(),
//                                           thevaluestr.c_str() );
//                  fmt::print( "set {} = {}\n", theoptionstr, thevaluestr );
//               }
//               catch( const std::exception& e )
//               {
//                  fmt::print( "parameter '{}' could not be set: {}\n", line,
//                              e.what() );
//               }
//            }
//         }
//         else
//         {
//            fmt::print( "could not read parameter file '{}'\n",
//                        opts.param_settings_file );
//         }
//      }
//   }


   double readtime = 0;
   Timer t( readtime );

   Problem<double> problem;
   Num<double> num{};
   Message msg{};
   boost::optional<Problem<double>> prob;

   {
      prob = MpsParser<double>::loadProblem( optionsInfo.instance_file );
   }

   // Check whether reading was successful or not
   if( !prob )
   {
      fmt::print( "error loading problem {}\n", optionsInfo.instance_file );
      return 0;
   }
   problem = *prob;

   fmt::print( "reading took {:.3} seconds\n", t.getTime() );

   double time_limit = 10 * 60;
   VolumeAlgorithmParameter<double> para{ 0.05, 0.1, 0.2,  0.0005, 2,
                                          2,    1.1, 0.66, 0.02,   0.001,
                                          0.02, 2,   20,  time_limit, 20, 0.7 };

#ifndef PAPILO_TBB
   msg.error("Please build with TBB to use the parallel feature!!!");
#endif

   Algorithm<double> alg{ msg, num, t, time_limit, 1};
   alg.solve_problem( problem, para );
   return 0;
}
