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

#include <boost/program_options.hpp>
#include <cassert>
#include <fstream>

using namespace papilo;

void
let_go(const papilo::OptionsInfo& opts);

void
fix_and_propagate( const papilo::Problem<double> _problem,
                   const papilo::Num<double>& _num );

std::pair<int, double>
select_diving_variable( const papilo::Problem<double>& _problem,
                        const papilo::ProbingView<double>& _fixed_variables );

papilo::ProbingView<double>
propagate_to_leaf_or_infeasibility( const papilo::Problem<double>& _problem, const papilo::Num<double>& _num );

papilo::Solution<double>
create_solution( const papilo::ProbingView<double>& _view );
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
let_go(const papilo::OptionsInfo& opts)
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
fix_and_propagate( const papilo::Problem<double> _problem, const papilo::Num<double>& _num )
{
   while(true)
   {
      papilo::ProbingView<double> probing_view =
          propagate_to_leaf_or_infeasibility( _problem, _num );
      if( probing_view.isInfeasible())
      {
         //TODO: backtrack
         //TODO: pass conflict
      }
      else
      {
         //TODO: store objective value and solution
         papilo::Solution<double> solution = create_solution(probing_view);
         fmt::print("found solution {}", _problem.computeSolObjective(solution.primal));
      }

      break;
   }
}

papilo::Solution<double>
create_solution( const papilo::ProbingView<double>& _view )
{
   papilo::Vec<double> values{};
   for(int i=0; i < _view.getProbingUpperBounds().size(); i++)
   {
      assert( _view.getProbingUpperBounds()[i] ==
              _view.getProbingLowerBounds()[i]);
      values.push_back( _view.getProbingUpperBounds()[i] );
   }
   papilo::Solution<double> solution {papilo::SolutionType::kPrimal, values};
   return solution;
}

papilo::ProbingView<double>
propagate_to_leaf_or_infeasibility( const papilo::Problem<double>& _problem, const papilo::Num<double>& _num )
{
   papilo::Vec<int> fixed_variables;
   papilo::Vec<double> fixed_values;
   papilo::ProbingView<double> probing_view{_problem, _num };
   //TODO: don't know if order in which the variables are applied can be extracted
   while(true)
   {
      std::pair<int, double> value =
          select_diving_variable( _problem, probing_view );
      if(value.first == -1)
         return probing_view;
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
         return probing_view;
      }
   }
}

std::pair<int, double>
select_diving_variable( const papilo::Problem<double>& _problem,
                        const papilo::ProbingView<double>& _probing_view )
{

   for( int i = 0; i < _problem.getNCols(); i++ )
   {
      _probing_view.getProbingUpperBounds();
      if( _probing_view.getProbingUpperBounds()[i]!=
          _probing_view.getProbingLowerBounds()[i] )
         return { i, 1 };
   }
   return {-1, -1};
}
