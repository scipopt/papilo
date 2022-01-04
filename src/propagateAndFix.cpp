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
let_go( const OptionsInfo& opts );

void
fix_and_propagate( const Problem<double>& _problem, const Num<double>& _num,
                   ProbingView<double>& probing_view );

Fixing<double>
select_diving_variable( const Problem<double>& _problem,
                        const ProbingView<double>& _probing_view );

void
propagate_to_leaf_or_infeasibility( const Problem<double>& _problem,
                                    const Num<double>& _num,
                                    ProbingView<double>& view );

Solution<double>
create_solution( const ProbingView<double>& _view );

double
modify_value_due_to_backtrack( double value );

// TODO: Probing does only work with boolean values
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

   let_go( optionsInfo );

   return 0;
}

void
let_go( const OptionsInfo& opts )
{
   double readtime = 0;
   Problem<double> problem;
   boost::optional<Problem<double>> prob;

   {
      Timer t( readtime );
      prob = MpsParser<double>::loadProblem( opts.instance_file );
   }

   // Check whether reading was successful or not
   if( !prob )
   {
      fmt::print( "error loading problem {}\n", opts.instance_file );
      return;
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
   fix_and_propagate( probUpdate.getProblem(), num, probing_view );
}

void
fix_and_propagate( const Problem<double>& _problem, const Num<double>& _num,
                   ProbingView<double>& probing_view )
{
   papilo::Vec<Fixing<double>> fixings{};
   while( true )
   {
      probing_view.reset();
      if( !fixings.empty() )
      {
         for(auto & fixing : fixings)
            probing_view.setProbingColumn( fixing.get_column_index(),
                                           fixing.get_value() );
      }

      propagate_to_leaf_or_infeasibility( _problem, _num, probing_view );

      if( probing_view.isInfeasible() )
      {
         // TODO: store fixings since they code the conflict
         fixings = probing_view.get_fixings();
         assert( !fixings.empty() );
         Fixing<double> infeasible_fixing = fixings[fixings.size() - 1];
         fixings[fixings.size() - 1] =
             ( Fixing<double> ){ infeasible_fixing.get_column_index(),
             modify_value_due_to_backtrack( infeasible_fixing.get_value() ) };
      }
      else
      {
         // TODO: store objective value and solution
         Solution<double> solution = create_solution( probing_view );
         fmt::print( "found solution {}",
                     _problem.computeSolObjective( solution.primal ) );
         break;
      }
   }
}

double
modify_value_due_to_backtrack( double value )
{
   return value == 1 ? 0 : 1;
}

Solution<double>
create_solution( const ProbingView<double>& _view )
{
   Vec<double> values{};
   for( int i = 0; i < _view.getProbingUpperBounds().size(); i++ )
   {
      assert( _view.getProbingUpperBounds()[i] ==
              _view.getProbingLowerBounds()[i] );
      values.push_back( _view.getProbingUpperBounds()[i] );
   }
   Solution<double> solution{ SolutionType::kPrimal, values };
   return solution;
}

void
propagate_to_leaf_or_infeasibility( const Problem<double>& _problem,
                                    const Num<double>& _num,
                                    ProbingView<double>& probing_view )
{
   while( true )
   {
      Fixing<double> fixing = select_diving_variable( _problem, probing_view );
      if( fixing.is_invalid() )
         return;

      fmt::print( "{} {}\n", fixing.get_column_index(), fixing.get_value() );

      probing_view.setProbingColumn( fixing.get_column_index(),
                                     fixing.get_value() == 1 );
      probing_view.propagateDomains();
      probing_view.storeImplications();
      if( probing_view.isInfeasible() )
      {
         fmt::print( "infeasible\n" );
         return;
      }
   }
}

Fixing<double>
select_diving_variable( const Problem<double>& _problem,
                        const ProbingView<double>& _probing_view )
{

   // TODO: currently a draft
   for( int i = 0; i < _problem.getNCols(); i++ )
   {
      _probing_view.getProbingUpperBounds();
      if( _probing_view.getProbingUpperBounds()[i] !=
          _probing_view.getProbingLowerBounds()[i] )
      {
         return { i, 0 };
      }
   }
   return { -1, -1 };
}
