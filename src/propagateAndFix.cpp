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


#include "papilo/core/Objective.hpp"
#include "papilo/core/Presolve.hpp"
#include "papilo/misc/MultiPrecision.hpp"
#include "papilo/misc/OptionsParser.hpp"
#include "papilo/misc/VersionLogger.hpp"
#include "papilo/misc/tbb.hpp"

#include <boost/program_options.hpp>
#include <cassert>
#include <fstream>

void
let_go(papilo::OptionsInfo& opts);

void
fix_and_propagate( papilo::Problem<double> _problem,
                   papilo::Num<double>& _num );

std::pair<int, double>
select_diving_variable( papilo::Problem<double> _problem,
                        papilo::Vec<int> _fixed_variables );
int
main( int argc, char* argv[] )
{

   // get the options passed by the user
   papilo::OptionsInfo optionsInfo;
   try
   {
      optionsInfo = papilo::parseOptions( argc, argv );
   }
   catch( const boost::program_options::error& ex )
   {
      std::cerr << "Error while parsing the options.\n" << '\n';
      std::cerr << ex.what() << '\n';
      return 1;
   }

   if( !optionsInfo.is_complete )
      return 0;

   let_go( optionsInfo );

   return 0;
}

void
let_go(papilo::OptionsInfo& opts)
{
      double readtime = 0;
      papilo::Problem<double> problem;
      boost::optional< papilo::Problem<double>> prob;

      {
         papilo::Timer t( readtime );
         prob =  papilo::MpsParser<double>::loadProblem( opts.instance_file );
      }

      // Check whether reading was successful or not
      if( !prob )
      {
         fmt::print( "error loading problem {}\n", opts.instance_file );
         return;
      }
      problem = *prob;

      fmt::print( "reading took {:.3} seconds\n", readtime );

      //do trivial presolve so that the activities
      papilo::Num<double> num{};
      papilo::Statistics stats{};
      papilo::Message msg{};
      papilo::PresolveOptions presolve_options{};
      papilo::PostsolveStorage<double> postsolve_storage;
      papilo::ProblemUpdate<double> probUpdate( problem, postsolve_storage, stats,
                                                presolve_options, num, msg );
      probUpdate.trivialPresolve();
      fix_and_propagate( probUpdate.getProblem(), num );

}

void
fix_and_propagate( papilo::Problem<double> _problem, papilo::Num<double>& _num )
{
   papilo::Vec<int> fixed_variables;
   papilo::Vec<double> fixed_values;
   papilo::ProbingView<double> probing_view{_problem, _num };
   while(fixed_variables.size() != _problem.getNCols())
   {
      std::pair<int, double> value =
          select_diving_variable( _problem, fixed_variables );
      fixed_variables.push_back(value.first);
      fixed_values.push_back(value.second);
      fmt::print("{} {}\n", value.first, value.second);

      //TODO: does only work with binary currently
      probing_view.setProbingColumn(value.first, value.second == 1);
      probing_view.propagateDomains();
      probing_view.storeImplications();
      if( probing_view.isInfeasible() )
      {
         fmt::print("infeasible\n");
         // TODO: implement reverting
      }
   }
}


std::pair<int, double>
select_diving_variable( papilo::Problem<double> _problem,
                        papilo::Vec<int> _fixed_variables )
{

   for(int i=0; i< _problem.getNCols(); i++)
      if(std::find(_fixed_variables.begin(), _fixed_variables.end(), i) == _fixed_variables.end())
         return {i, 0};

   assert(false);
}
