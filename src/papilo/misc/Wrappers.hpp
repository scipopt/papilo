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

#ifndef _PAPILO_MISC_WRAPPERS_HPP_
#define _PAPILO_MISC_WRAPPERS_HPP_

#include "papilo/core/Postsolve.hpp"
#include "papilo/core/Presolve.hpp"
#include "papilo/io/MpsParser.hpp"
#include "papilo/io/MpsWriter.hpp"
#include "papilo/io/SolParser.hpp"
#include "papilo/io/SolWriter.hpp"
#include "papilo/misc/OptionsParser.hpp"
#include "papilo/misc/tbb.hpp"

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/lexical_cast.hpp>
#include <fstream>
#include <string>
#include <utility>

namespace papilo
{

enum class ResultStatus
{
   kOk = 0,
   kUnbndOrInfeas,
   kError
};

template <typename REAL>
ResultStatus
presolve_and_solve(
    const OptionsInfo& opts,
    std::unique_ptr<SolverFactory<REAL>> lpSolverFactory = nullptr,
    std::unique_ptr<SolverFactory<REAL>> mipSolverFactory = nullptr )
{

   double readtime = 0;
   Problem<REAL> problem;

   {
      Timer t( readtime );
      problem = MpsParser<REAL>::loadProblem( opts.instance_file );
   }

   fmt::print( "reading took {:.3} seconds\n", readtime );

   Presolve<REAL> presolve;
   presolve.addDefaultPresolvers();
   presolve.getPresolveOptions().threads = std::max( 0, opts.nthreads );

   if( !opts.param_settings_file.empty() || !opts.unparsed_options.empty() ||
       opts.print_params )
   {
      ParameterSet paramSet = presolve.getParameters();

      if( !opts.param_settings_file.empty() && !opts.print_params )
      {
         std::ifstream input( opts.param_settings_file );
         if( input )
         {
            String theoptionstr;
            String thevaluestr;
            for( String line; getline( input, line ); )
            {
               std::size_t pos = line.find_first_of( '#' );
               if( pos != String::npos )
                  line = line.substr( 0, pos );

               pos = line.find_first_of( '=' );

               if( pos == String::npos )
                  continue;

               theoptionstr = line.substr( 0, pos - 1 );
               thevaluestr = line.substr( pos + 1 );

               boost::algorithm::trim( theoptionstr );
               boost::algorithm::trim( thevaluestr );

               try
               {
                  paramSet.parseParameter( theoptionstr.c_str(),
                                           thevaluestr.c_str() );
                  fmt::print( "set {} = {}\n", theoptionstr, thevaluestr );
               }
               catch( const std::exception& e )
               {
                  fmt::print( "parameter '{}' could not be set: {}\n", line,
                              e.what() );
               }
            }
         }
         else
         {
            fmt::print( "could not read parameter file '{}'\n",
                        opts.param_settings_file );
         }
      }

      if( !opts.unparsed_options.empty() )
      {
         String theoptionstr;
         String thevaluestr;

         for( const auto& option : opts.unparsed_options )
         {
            std::size_t pos = option.find_first_of( '=' );
            if( pos != String::npos && pos > 2 )
            {
               theoptionstr = option.substr( 2, pos - 2 );
               thevaluestr = option.substr( pos + 1 );
               try
               {
                  paramSet.parseParameter( theoptionstr.c_str(),
                                           thevaluestr.c_str() );
                  fmt::print( "set {} = {}\n", theoptionstr, thevaluestr );
               }
               catch( const std::exception& e )
               {
                  fmt::print( "parameter '{}' could not be set: {}\n", option,
                              e.what() );
               }
            }
            else
            {
               fmt::print( "parameter '{}' could not be set: value expected\n",
                           option );
            }
         }
      }

      if( opts.print_params )
      {
         if( !opts.param_settings_file.empty() )
         {
            std::ofstream outfile( opts.param_settings_file );

            if( outfile )
            {
               std::ostream_iterator<char> out_it( outfile );
               paramSet.printParams( out_it );
            }
            else
            {
               fmt::print( "could not write to parameter file '{}'\n",
                           opts.param_settings_file );
            }
         }
         else
         {
            String paramDesc;
            paramSet.printParams( std::back_inserter( paramDesc ) );
            puts( paramDesc.c_str() );
         }
      }
   }

   presolve.setLPSolverFactory( std::move( lpSolverFactory ) );
   presolve.setMIPSolverFactory( std::move( mipSolverFactory ) );

   presolve.getPresolveOptions().tlim = opts.tlim;

   auto result = presolve.apply( problem );

   switch( result.status )
   {
   case PresolveStatus::kInfeasible:
      fmt::print( "presolve detected infeasible problem\n" );
      return ResultStatus::kUnbndOrInfeas;
   case PresolveStatus::kUnbndOrInfeas:
      fmt::print( "presolve detected unbounded or infeasible problem\n" );
      return ResultStatus::kUnbndOrInfeas;
   case PresolveStatus::kUnbounded:
      fmt::print( "presolve detected unbounded problem\n" );
      return ResultStatus::kUnbndOrInfeas;
   case PresolveStatus::kUnchanged:
   case PresolveStatus::kReduced:
      break;
   }

   fmt::print( "\npresolving finished after {:.3f} seconds\n\n",
               presolve.getStatistics().presolvetime );

   double writetime = 0;
   if( !opts.reduced_problem_file.empty() )
   {
      Timer t( writetime );
      const auto t0 = tbb::tick_count::now();

      MpsWriter<REAL>::writeProb( opts.reduced_problem_file, problem,
                                  result.postsolve.origrow_mapping,
                                  result.postsolve.origcol_mapping );

      const auto t1 = tbb::tick_count::now();
      fmt::print( "reduced problem written to {} in {:.3f} seconds\n\n",
                  opts.reduced_problem_file, t.getTime() );
   }

   if( !opts.postsolve_archive_file.empty() )
   {

      Timer t( writetime );
      std::ofstream ofs( opts.postsolve_archive_file, std::ios_base::binary );
      boost::archive::binary_oarchive oa( ofs );

      // write class instance to archive
      oa << result.postsolve;
      fmt::print( "postsolve archive written to {} in {:.3f} seconds\n\n",
                  opts.postsolve_archive_file, t.getTime() );
   }

   if( opts.command == Command::kPresolve )
      return ResultStatus::kOk;

   double solvetime = 0;
   {
      Timer t( solvetime );
      std::unique_ptr<SolverInterface<REAL>> solver;

      if( result.postsolve.getOriginalProblem().getNumIntegralCols() == 0 &&
          presolve.getLPSolverFactory() )
         solver = presolve.getLPSolverFactory()->newSolver(
             presolve.getVerbosityLevel() );
      else if( presolve.getMIPSolverFactory() )
         solver = presolve.getMIPSolverFactory()->newSolver(
             presolve.getVerbosityLevel() );
      else
      {
         fmt::print( "no solver available for solving\n" );
         return ResultStatus::kError;
      }

      solver->setUp( problem, result.postsolve.origrow_mapping,
                     result.postsolve.origcol_mapping );

      if( opts.tlim != std::numeric_limits<double>::max() )
      {
         double tlim =
             opts.tlim - presolve.getStatistics().presolvetime - writetime;
         if( tlim <= 0 )
         {
            fmt::print( "time limit reached in presolving\n" );
            return ResultStatus::kOk;
         }
         solver->setTimeLimit( tlim );
      }

      solver->solve();

      SolverStatus status = solver->getStatus();

      if( opts.print_stats )
         solver->printDetails();

      Solution<REAL> solution;
      solution.type = SolutionType::kPrimal;

      if( result.postsolve.getOriginalProblem().getNumIntegralCols() == 0 )
         solution.type = SolutionType::kPrimalDual;

      if( ( status == SolverStatus::kOptimal ||
            status == SolverStatus::kInterrupted ) &&
          solver->getSolution( solution ) )
         postsolve( result.postsolve, solution, opts.objective_reference,
                    opts.orig_solution_file );
   }

   fmt::print( "\nsolving finished after {:.3f} seconds\n",
               presolve.getStatistics().presolvetime + solvetime + writetime );

   return ResultStatus::kOk;
}

template <typename REAL>
void
postsolve( Postsolve<REAL>& postsolve, const Solution<REAL>& reduced_sol,
           const std::string& objective_reference = "",
           const std::string& solution_output = "" )
{
   Solution<REAL> original_sol;

   auto t0 = tbb::tick_count::now();
   PostsolveStatus status = postsolve.undo( reduced_sol, original_sol );
   auto t1 = tbb::tick_count::now();

   fmt::print( "\npostsolve finished after {:.3f} seconds\n",
               ( t1 - t0 ).seconds() );

   const Problem<REAL>& origprob = postsolve.getOriginalProblem();
   REAL origobj = origprob.computeSolObjective( original_sol.primal );

   REAL boundviol;
   REAL intviol;
   REAL rowviol;
   bool origfeas = origprob.computeSolViolations(
       postsolve.getNum(), original_sol.primal, boundviol, rowviol, intviol );

   fmt::print( "feasible: {}\nobjective value: {:.15}\n", origfeas,
               double( origobj ) );

   fmt::print( "\nviolations:\n" );
   fmt::print( "  bounds:      {:.15}\n", double( boundviol ) );
   fmt::print( "  constraints: {:.15}\n", double( rowviol ) );
   fmt::print( "  integrality: {:.15}\n\n", double( intviol ) );

   if( !solution_output.empty() )
   {
      auto t0 = tbb::tick_count::now();
      SolWriter<REAL>::writeSol( solution_output, original_sol.primal,
                                 origprob.getObjective().coefficients, origobj,
                                 origprob.getVariableNames() );
      auto t1 = tbb::tick_count::now();

      fmt::print( "solution written to file {} in {:.3} seconds\n",
                  solution_output, ( t1 - t0 ).seconds() );
   }

   if( !objective_reference.empty() )
   {
      if( origfeas &&
          postsolve.num.isFeasEq(
              boost::lexical_cast<double>( objective_reference ), origobj ) )
         fmt::print( "validation: SUCCESS\n" );
      else
         fmt::print( "validation: FAILURE\n" );
   }
}

template <typename REAL>
void
postsolve( const OptionsInfo& opts )
{
   Postsolve<REAL> ps;
   std::ifstream inArchiveFile( opts.postsolve_archive_file,
                                std::ios_base::binary );
   boost::archive::binary_iarchive inputArchive( inArchiveFile );
   inputArchive >> ps;
   inArchiveFile.close();

   SolParser<REAL> parser;
   std::ifstream solFile( opts.reduced_solution_file );
   Solution<REAL> reduced_solution =
       parser.read( solFile, ps.origcol_mapping,
                    ps.getOriginalProblem().getVariableNames() );
   solFile.close();

   postsolve( ps, reduced_solution, opts.objective_reference,
              opts.orig_solution_file );
}

} // namespace papilo

#endif
