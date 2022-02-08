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
                                          2,    1.1, 0.66, 0.01,   0.001,
                                          0.02, 2,   20,   time_limit };

#ifndef PAPILO_TBB
   msg.error("Please build with TBB to use the parallel feature!!!");
#endif

   Algorithm<double> alg{ msg, num, t, time_limit};
   alg.solve_problem( problem, para );
   return 0;
}
