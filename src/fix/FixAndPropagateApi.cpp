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
//#define FIX_DEBUG

#include "fix/FixAndPropagateApi.h"
#include "fix/Heuristic.hpp"

#include "papilo/io/MpsParser.hpp"
#include <string>

#ifdef FIX_DEBUG
#include "papilo/io/SolWriter.hpp"
#include "papilo/io/MpsWriter.hpp"
#endif

using namespace papilo;

void*
setup( const char* filename, int* result, int verbosity_level, double current_time_stamp )
{

   std::string filename_as_string( filename );
   boost::optional<Problem<double>> prob;
   {
      prob = MpsParser<double>::loadProblem( filename_as_string );
   }
   if( !prob )
   {
      fmt::print( "error loading problem {}\n", filename );
      *result = -1;
      return nullptr;
   }
   double tolerance = 1e-5;
   RandomGenerator random( 0 );
   Timer t{ current_time_stamp };
   auto problem = new Problem<double>( prob.get() );
   problem->recomputeAllActivities();
   Message msg{};
   msg.setVerbosityLevel( static_cast<VerbosityLevel>(verbosity_level) );
   PostsolveStorage<double> storage{};
   Num<double> num{};
   num.setEpsilon( tolerance );
   num.setFeasTol( tolerance );
   assert( num.isEq( 80, 80 + tolerance / 10 ) );
   auto heuristic =
       new Heuristic<double>{ msg, num, random, t, *problem, storage, false };
   heuristic->setup( random );
   *result = 0;
   return heuristic;
}

void
delete_problem_instance( void* heuristic_void_ptr )
{
   auto heuristic = (Heuristic<double>*)( heuristic_void_ptr );
   delete heuristic;
}

int
call_algorithm( void* heuristic_void_ptr, double* cont_solution, double* result,
                int n_cols, double* current_obj_value,
                int infeasible_copy_strategy, int apply_conflicts,
                int size_of_constraints, int max_backtracks,
                int perform_one_opt, double remaining_time_in_sec )
{
   assert( infeasible_copy_strategy >= 0 && infeasible_copy_strategy <= 6 );
   assert( apply_conflicts >= 0 && apply_conflicts <= 1 );
   assert( max_backtracks >= 0 );
   assert( perform_one_opt >= 0 && perform_one_opt <= 2 );
   assert( size_of_constraints >= 0 );
   assert( n_cols ==
           ( (Heuristic<double>*)( heuristic_void_ptr ) )->problem.getNCols() );

#ifdef PAPILO_TBB
   tbb::task_arena arena( 7 );

   return arena.execute(
       [&]()
       {
#endif
          auto heuristic = (Heuristic<double>*)( heuristic_void_ptr );
          if( apply_conflicts == 1 &&
              size_of_constraints <= heuristic->get_derived_conflicts().size() )
          {
             heuristic->problem = heuristic->copy_conflicts_to_problem(
                 heuristic->problem, heuristic->get_derived_conflicts() );
             heuristic->get_message().info(
                 "added {} conflicts to the problem (rows: {})\n",
                 heuristic->get_derived_conflicts().size(),
                 heuristic->problem.getNRows() );
             heuristic->problem.recomputeAllActivities();
             heuristic->get_derived_conflicts().clear();
          }
          Vec<double> sol( cont_solution, cont_solution + n_cols );
          Vec<double> res{};

#ifdef FIX_DEBUG
          SolWriter<double>::writePrimalSol(
              "lp_feasible.sol", sol,
              heuristic->problem.getObjective().coefficients, 0.0,
              heuristic->problem.getVariableNames() );
          Vec<int> row_mapping{};
          Vec<int> col_mapping{};
          for( int i = 0; i < heuristic->problem.getNRows(); i++ )
             row_mapping.push_back( i );
          for( int i = 0; i < heuristic->problem.getNCols(); i++ )
             col_mapping.push_back( i );
          MpsWriter<double>::writeProb( "test.mps", heuristic->problem,
                                        row_mapping, col_mapping );
#endif

          double i = heuristic->get_current_time();
          double time_limit = i + remaining_time_in_sec;
          double local_obj = *current_obj_value;
          bool better_solution_found = heuristic->perform_fix_and_propagate(
              sol, local_obj, res, time_limit, max_backtracks, perform_one_opt,
              false, (InfeasibleCopyStrategy)infeasible_copy_strategy );
          std::copy( res.begin(), res.end(), result );
          if( !better_solution_found )
             return false;
          assert( local_obj < *current_obj_value );
          *current_obj_value = local_obj;
          return true;
#ifdef PAPILO_TBB
       } );
#endif
}

int
call_simple_heuristic( void* heuristic_void_ptr, double* result,
                       double* current_obj_value )
{
#ifdef PAPILO_TBB
   tbb::task_arena arena( 7 );
   return arena.execute(
       [&]()
       {
#endif
          auto heuristic = (Heuristic<double>*)( heuristic_void_ptr );
          Vec<double> res{};
          heuristic->find_initial_solution( *current_obj_value, res );
          std::copy( res.begin(), res.end(), result );
          return !res.empty();
#ifdef PAPILO_TBB
       } );
#endif
}

void
perform_one_opt( void* heuristic_void_ptr, double* sol, int n_cols,
                 int perform_opt_one, double* current_obj_value,
                 double remaining_time_in_sec )
{
   auto heuristic = (Heuristic<double>*)( heuristic_void_ptr );
   Vec<double> res{ sol, sol + n_cols };

   double i = heuristic->get_current_time();
   double time_limit = i + remaining_time_in_sec;
   ProbingView<double> view {heuristic-> problem, heuristic->get_num()};
   heuristic->perform_one_opt( perform_opt_one, res, view, *current_obj_value,
                               -1, time_limit );
   std::copy( res.begin(), res.end(), sol );
}
