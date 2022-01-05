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

   FixAndPropagate<double> fixAndPropagate{};

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
   Num<double> num{};
   Statistics stats{};
   Message msg{};
   PresolveOptions presolve_options{};
   PostsolveStorage<double> postsolve_storage;
   ProblemUpdate<double> probUpdate( problem, postsolve_storage, stats,
                                   presolve_options, num, msg );
   probUpdate.trivialPresolve();

   ProbingView<double> probing_view{ problem, num };
   fixAndPropagate.fix_and_propagate( probUpdate.getProblem(), num,
                                      probing_view );

   return 0;
}
