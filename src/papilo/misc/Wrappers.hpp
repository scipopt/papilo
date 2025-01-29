/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*               This file is part of the program and library                */
/*    PaPILO --- Parallel Presolve for Integer and Linear Optimization       */
/*                                                                           */
/* Copyright (C) 2020-2025 Zuse Institute Berlin (ZIB)                       */
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

#include "papilo/Config.hpp"
#include "papilo/core/Presolve.hpp"
#include "papilo/core/postsolve/Postsolve.hpp"
#include "papilo/io/MpsWriter.hpp"
#include "papilo/io/OpbWriter.hpp"
#include "papilo/io/Parser.hpp"
#include "papilo/io/SolParser.hpp"
#include "papilo/io/SolWriter.hpp"
#include "papilo/misc/NumericalStatistics.hpp"
#include "papilo/misc/OptionsParser.hpp"
#include "papilo/misc/Validation.hpp"
#ifdef PAPILO_TBB
#include "papilo/misc/tbb.hpp"
#else
#include <chrono>
#endif

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
    std::unique_ptr<SolverFactory<REAL>> mipSolverFactory = nullptr,
    std::unique_ptr<SolverFactory<REAL>> satSolverFactory = nullptr )
{
   try
   {
      double readtime = 0;
      Problem<REAL> problem;
      boost::optional<Problem<REAL>> prob;

      {
         Timer t( readtime );
         prob = Parser<REAL>::loadProblem( opts.instance_file );
      }

      // Check whether reading was successful or not
      if( !prob )
      {
         fmt::print( "error loading problem {}\n", opts.instance_file );
         return ResultStatus::kError;
      }
      problem = *prob;

      fmt::print( "reading took {:.3} seconds\n", readtime );

      NumericalStatistics<REAL> nstats( problem );
      nstats.printStatistics();

      Presolve<REAL> presolve;
      presolve.addDefaultPresolvers();
      presolve.getPresolveOptions().threads = std::max( 0, opts.nthreads );

      if( !opts.param_settings_file.empty() || !opts.unparsed_options.empty() ||
          opts.print_params )
      {
         ParameterSet paramSet = presolve.getParameters();
         if(lpSolverFactory != nullptr)
            lpSolverFactory->add_parameters(paramSet);
         if(mipSolverFactory != nullptr)
            mipSolverFactory->add_parameters(paramSet);
         if(satSolverFactory != nullptr)
            satSolverFactory->add_parameters(paramSet);

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
                     fmt::print( "parameter '{}' could not be set: {}\n",
                                 option, e.what() );
                  }
               }
               else
               {
                  fmt::print(
                      "parameter '{}' could not be set: value expected\n",
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
      presolve.setSATSolverFactory( std::move( satSolverFactory ) );

      presolve.getPresolveOptions().tlim =
          std::min( opts.tlim, presolve.getPresolveOptions().tlim );

      bool store_dual = false;
      std::unique_ptr<SolverInterface<REAL>> solver;
      if(opts.command == Command::kSolve)
      {
         if( problem.getNumIntegralCols() == 0 &&
             presolve.getLPSolverFactory() )
         {
            solver = presolve.getLPSolverFactory()->newSolver(
                presolve.getVerbosityLevel() );
            store_dual = solver->is_dual_solution_available();
         }
         else if( presolve.getSATSolverFactory() &&
                  problem.test_problem_type(ProblemFlag::kBinary) )
         {
            solver = presolve.getSATSolverFactory()->newSolver(
                presolve.getVerbosityLevel() );
         }
         else if( presolve.getMIPSolverFactory() )
            solver = presolve.getMIPSolverFactory()->newSolver(
                presolve.getVerbosityLevel() );
         else
         {
            fmt::print( "no solver available for solving; aborting\n" );
            return ResultStatus::kError;
         }
      }
      else
      {
         store_dual = true;
      }

      auto result = presolve.apply( problem, store_dual );

      if( !opts.optimal_solution_file.empty() )
      {
         Validation<REAL>::validateProblem( problem, result.postsolve,
                                            opts.optimal_solution_file,
                                            result.status );
      }

      switch( result.status )
      {
      case PresolveStatus::kInfeasible:
         fmt::print(
             "presolving detected infeasible problem after {:.3f} seconds\n",
             presolve.getStatistics().presolvetime );
         return ResultStatus::kUnbndOrInfeas;
      case PresolveStatus::kUnbndOrInfeas:
         fmt::print(
             "presolving detected unbounded or infeasible problem after "
             "{:.3f} seconds\n",
             presolve.getStatistics().presolvetime );
         return ResultStatus::kUnbndOrInfeas;
      case PresolveStatus::kUnbounded:
         fmt::print(
             "presolving detected unbounded problem after {:.3f} seconds\n",
             presolve.getStatistics().presolvetime );
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

         if(boost::algorithm::ends_with( opts.reduced_problem_file, ".opb" )
             || boost::algorithm::ends_with( opts.reduced_problem_file, ".opb.bz2" )
             || boost::algorithm::ends_with( opts.reduced_problem_file, ".opb.gz" ))
         {
            Num<REAL> num{};
            num.setFeasTol( REAL{ presolve.getPresolveOptions().feastol } );
            num.setEpsilon( REAL{ presolve.getPresolveOptions().epsilon } );
            num.setHugeVal( REAL{ presolve.getPresolveOptions().hugeval } );
            num.setUseAbsFeas( presolve.getPresolveOptions().useabsfeas );
            bool success = OpbWriter<REAL>::writeProb( opts.reduced_problem_file, problem,
                                        result.postsolve.origcol_mapping, presolve.getRowScalingFactors(), num );
            //TODO: change name
            if(!success)
               MpsWriter<REAL>::writeProb( opts.reduced_problem_file, problem,
                                           result.postsolve.origrow_mapping,
                                           result.postsolve.origcol_mapping );
         }
         else
            MpsWriter<REAL>::writeProb( opts.reduced_problem_file, problem,
                                        result.postsolve.origrow_mapping,
                                        result.postsolve.origcol_mapping );

         fmt::print( "reduced problem written to {} in {:.3f} seconds\n\n",
                     opts.reduced_problem_file, t.getTime() );
      }

      if( !opts.postsolve_archive_file.empty() )
      {

#ifdef PAPILO_SERIALIZATION_AVAILABLE
         Timer t( writetime );
         std::ofstream ofs( opts.postsolve_archive_file,
                            std::ios_base::binary );
         boost::archive::binary_oarchive oa( ofs );

         // write class instance to archive
         oa << result.postsolve;
         fmt::print( "postsolve archive written to {} in {:.3f} seconds\n\n",
                     opts.postsolve_archive_file, t.getTime() );
#else
         fmt::print( "writing to potsolve archive is not possible without Boost Serialization!\n");
#endif
      }

      if( opts.command == Command::kPresolve || problem.getNCols() == 0 )
         return ResultStatus::kOk;

      double solvetime = 0;
      {
         Timer t( solvetime );

         if(presolve.getPresolveOptions().verification_with_VeriPB)
            solver->setRowScalingFactor(presolve.getRowScalingFactors());
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

         if( result.postsolve.getOriginalProblem().getNumIntegralCols() == 0 && store_dual )
            solution.type = SolutionType::kPrimalDual;

         if( ( status == SolverStatus::kOptimal ||
               status == SolverStatus::kInterrupted ) &&
             solver->getSolution( solution, result.postsolve ) )
            postsolve( result.postsolve, solution, opts.objective_reference,
                       opts.orig_solution_file, opts.orig_dual_solution_file,
                       opts.orig_reduced_costs_file, opts.orig_basis_file );
         solvetime = t.getTime();
         double time =
             presolve.getStatistics().presolvetime + solvetime + writetime;
         fmt::print("Solving time {:.3f} seconds\n", time);
         switch( status )
         {
      case SolverStatus::kInfeasible:
            fmt::print(
                "\nsolving detected infeasible problem after {:.3f} seconds\n",
                time );
            break;
         case SolverStatus::kUnbounded:
            fmt::print(
                "\nsolving detected unbounded problem after {:.3f} seconds\n",
                time );
            break;
         case SolverStatus::kUnbndOrInfeas:
            fmt::print(
                "\nsolving detected unbounded or infeasible problem after "
                "{:.3f} seconds\n",
                time );
            break;
         case SolverStatus::kInterrupted:
         case SolverStatus::kError:
            fmt::print( "\nsolving interrupted after {:.3f} seconds\n", time );
            break;
         default:
            fmt::print( "\nsolving finished after {:.3f} seconds\n", time );
            break;
         }
      }
   }
   catch( std::bad_alloc& ex )
   {
      fmt::print( "Memory out exception occurred! Please assign more memory\n" );
      return ResultStatus::kError;
   }

   return ResultStatus::kOk;
}

template <typename REAL>
void
postsolve( PostsolveStorage<REAL>& postsolveStorage,
           const Solution<REAL>& reduced_sol,
           const std::string& objective_reference = "",
           const std::string& primal_solution_output = "",
           const std::string& dual_solution_output = "",
           const std::string& reduced_solution_output = "",
           const std::string& basis_output = "" )
{
   Solution<REAL> original_sol;

#ifdef PAPILO_TBB
   auto t0 = tbb::tick_count::now();
#else
   auto t0 = std::chrono::steady_clock::now();
#endif
   const Message msg{};
   Postsolve<REAL> postsolve{ msg, postsolveStorage.getNum() };
   PostsolveStatus status =
       postsolve.undo( reduced_sol, original_sol, postsolveStorage );
#ifdef PAPILO_TBB
   auto t1 = tbb::tick_count::now();
   double sec1 = ( t1 - t0 ).seconds();
#else
   auto t1 = std::chrono::steady_clock::now();
   double sec1 =
       std::chrono::duration_cast<std::chrono::milliseconds>( t1 - t0 )
           .count() /
       1000;
#endif
   fmt::print( "\npostsolve finished after {:.3f} seconds\n", sec1 );

   const Problem<REAL>& origprob = postsolveStorage.getOriginalProblem();
   REAL origobj = origprob.computeSolObjective( original_sol.primal );

   REAL boundviol = 0;
   REAL intviol = 0;
   REAL rowviol = 0;
   bool origfeas = origprob.computeSolViolations( postsolveStorage.getNum(),
                                                  original_sol.primal,
                                                  boundviol, rowviol, intviol );

   fmt::print( "feasible: {}\nobjective value: {:.15}\n", origfeas,
               double( origobj ) );

   fmt::print( "\nviolations:\n" );
   fmt::print( "  bounds:      {:.15}\n", double( boundviol ) );
   fmt::print( "  constraints: {:.15}\n", double( rowviol ) );
   fmt::print( "  integrality: {:.15}\n\n", double( intviol ) );

   if( !primal_solution_output.empty() )
   {
#ifdef PAPILO_TBB
      auto t2 = tbb::tick_count::now();
#else
      auto t2 = std::chrono::steady_clock::now();
#endif
      SolWriter<REAL>::writePrimalSol( primal_solution_output,
                                       original_sol.primal,
                                       origprob.getObjective().coefficients,
                                       origobj, origprob.getVariableNames() );
#ifdef PAPILO_TBB
      auto t3 = tbb::tick_count::now();
      double sec3 = ( t3 - t2 ).seconds();
#else
      auto t3 = std::chrono::steady_clock::now();
      double sec3 =
          std::chrono::duration_cast<std::chrono::milliseconds>( t3 - t2 )
              .count() /
          1000;
#endif

      fmt::print( "solution written to file {} in {:.3} seconds\n",
                  primal_solution_output, sec3 );
   }

   if( ( !reduced_solution_output.empty() || !dual_solution_output.empty() ||
         !basis_output.empty() ) &&
       original_sol.type == SolutionType::kPrimal )
   {
      fmt::print(
          "dual solution does not exist and therefore can not be saved.\n" );
   }

   if( !dual_solution_output.empty() &&
       original_sol.type == SolutionType::kPrimalDual )
   {
#ifdef PAPILO_TBB
      auto t2 = tbb::tick_count::now();
#else
      auto t2 = std::chrono::steady_clock::now();
#endif
      const ConstraintMatrix<REAL>& constraintMatrix =
          origprob.getConstraintMatrix();
      SolWriter<REAL>::writeDualSol( dual_solution_output, original_sol.dual,
                                     constraintMatrix.getRightHandSides(),
                                     constraintMatrix.getLeftHandSides(),
                                     origobj, origprob.getConstraintNames() );
#ifdef PAPILO_TBB
      auto t3 = tbb::tick_count::now();
      double sec3 = ( t3 - t2 ).seconds();
#else
      auto t3 = std::chrono::steady_clock::now();
      double sec3 =
          std::chrono::duration_cast<std::chrono::milliseconds>( t3 - t2 )
              .count() /
          1000;
#endif

      fmt::print( "dual solution written to file {} in {:.3} seconds\n",
                  dual_solution_output, sec3 );
   }

   if( !reduced_solution_output.empty() &&
       original_sol.type == SolutionType::kPrimalDual )
   {
#ifdef PAPILO_TBB
      auto t2 = tbb::tick_count::now();
#else
      auto t2 = std::chrono::steady_clock::now();
#endif
      SolWriter<REAL>::writeReducedCostsSol(
          reduced_solution_output, original_sol.reducedCosts,
          origprob.getUpperBounds(), origprob.getLowerBounds(), origobj,
          origprob.getVariableNames() );
#ifdef PAPILO_TBB
      auto t3 = tbb::tick_count::now();
      double sec3 = ( t3 - t2 ).seconds();
#else
      auto t3 = std::chrono::steady_clock::now();
      double sec3 =
          std::chrono::duration_cast<std::chrono::milliseconds>( t3 - t2 )
              .count() /
          1000;
#endif

      fmt::print( "reduced solution written to file {} in {:.3} seconds\n",
                  reduced_solution_output, sec3 );
   }

   if( !basis_output.empty() && original_sol.type == SolutionType::kPrimalDual )
   {
#ifdef PAPILO_TBB
      auto t2 = tbb::tick_count::now();
#else
      auto t2 = std::chrono::steady_clock::now();
#endif

      SolWriter<REAL>::writeBasis( basis_output, original_sol.varBasisStatus,
                                   original_sol.rowBasisStatus,
                                   origprob.getVariableNames(),
                                   origprob.getConstraintNames() );
#ifdef PAPILO_TBB
      auto t3 = tbb::tick_count::now();
      double sec3 = ( t3 - t2 ).seconds();
#else
      auto t3 = std::chrono::steady_clock::now();
      double sec3 =
          std::chrono::duration_cast<std::chrono::milliseconds>( t3 - t2 )
              .count() /
          1000;
#endif
      fmt::print( "basis written to file {} in {:.3} seconds\n", basis_output,
                  sec3 );
   }

   if( !objective_reference.empty() )
   {
      if( origfeas && status == PostsolveStatus::kOk &&
          postsolveStorage.num.isFeasEq(
              boost::lexical_cast<double>( objective_reference ),
              double( origobj ) ) )
         fmt::print( "validation: SUCCESS\n" );
      else
         fmt::print( "validation: FAILURE\n" );
   }
}

template <typename REAL>
void
postsolve( const OptionsInfo& opts )
{
#ifdef PAPILO_SERIALIZATION_AVAILABLE
   PostsolveStorage<REAL> ps;
   std::ifstream inArchiveFile( opts.postsolve_archive_file,
                                std::ios_base::binary );
   boost::archive::binary_iarchive inputArchive( inArchiveFile );
   inputArchive >> ps;
   inArchiveFile.close();

   SolParser<REAL> parser;
   Vec<REAL> primal_solution;
   bool success = parser.read( opts.reduced_solution_file, ps.origcol_mapping,
                               ps.getOriginalProblem().getVariableNames(),
                               primal_solution );
   Solution<REAL> reduced_solution{primal_solution};
   if( success )
   {
      if( !opts.reduced_dual_solution_file.empty() !=
          !opts.reduced_reduced_costs_solution_file.empty() )
      {
         fmt::print( "Dual solution and reduced costs must be provided. Only "
                     "original primal solution reconstructed.\n" );
      }
      else if( !opts.reduced_dual_solution_file.empty() &&
          !opts.reduced_reduced_costs_solution_file.empty() )
      {
         if(ps.postsolveType == PostsolveType::kFull)
         {
            bool dual_success = parser.read(
                opts.reduced_dual_solution_file, ps.origrow_mapping,
                ps.getOriginalProblem().getConstraintNames(),
                reduced_solution.dual );

            if( !dual_success )
               return;

            bool reduced_costs_success = parser.read(
                opts.reduced_reduced_costs_solution_file, ps.origcol_mapping,
                ps.getOriginalProblem().getVariableNames(),
                reduced_solution.reducedCosts );
            if( !reduced_costs_success )
               return;
            reduced_solution.type = SolutionType::kPrimalDual;

            if( !opts.reduced_basis_solution_file.empty() && ps.presolveOptions.calculate_basis_for_dual )
            {
                  reduced_solution.basisAvailabe = true;
                  parser.read_basis( opts.reduced_basis_solution_file, ps,
                                     reduced_solution.varBasisStatus,
                                     reduced_solution.rowBasisStatus );
            }
            else
            {
               if(!opts.reduced_basis_solution_file.empty())
                  fmt::print("Postsolve storage not suitable for basis calculation.\n");
               reduced_solution.basisAvailabe = false;
               reduced_solution.varBasisStatus =
                   Vec<VarBasisStatus>{ reduced_solution.primal.size() };
               reduced_solution.rowBasisStatus =
                   Vec<VarBasisStatus>( reduced_solution.dual.size() );
            }
         }
         else
         {
            fmt::print(
                "Postsolve does not contain information about dual solution. "
                "Only original primal solution reconstructed.\n" );
         }
      }
      postsolve( ps, reduced_solution, opts.objective_reference,
                 opts.orig_solution_file, opts.orig_dual_solution_file,
                 opts.orig_reduced_costs_file, opts.orig_basis_file );


   }
#endif
   }

} // namespace papilo

#endif
