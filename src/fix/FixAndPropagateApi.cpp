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
#endif

using namespace papilo;

void*
setup( const char* filename, int* result, int verbosity_level )
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
   double time = 0;
   Timer t{ time };
   auto problem = new Problem<double>( prob.get() );
   problem->recomputeAllActivities();
   Message msg{};
   switch( verbosity_level )
   {
   case 0:
      msg.setVerbosityLevel( papilo::VerbosityLevel::kQuiet );
      break;
   case 1:
      msg.setVerbosityLevel( papilo::VerbosityLevel::kError );
      break;
   case 2:
      msg.setVerbosityLevel( papilo::VerbosityLevel::kWarning );
      break;
   case 3:
      msg.setVerbosityLevel( papilo::VerbosityLevel::kInfo );
      break;
   case 4:
      msg.setVerbosityLevel( papilo::VerbosityLevel::kDetailed );
      break;
   default:
      assert( false );
   }
   PostsolveStorage<double> storage{};
   auto heuristic =
       new Heuristic<double>{ msg, {}, t, *problem, storage, false };
   heuristic->setup();
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
                int n_cols, double* current_obj_value )
{
#ifdef PAPILO_TBB
   tbb::task_arena arena( 8 );

   return arena.execute(
       [&]()
       {
#endif
          auto heuristic = (Heuristic<double>*)( heuristic_void_ptr );
          Vec<double> sol( cont_solution, cont_solution + n_cols );
          Vec<double> res{};

#ifdef FIX_DEBUG
          SolWriter<double>::writePrimalSol(
              "test.mps", sol, heuristic->problem.getObjective().coefficients, 0.0,
              heuristic->problem.getVariableNames() );
#endif

          double local_obj = *current_obj_value;
          heuristic->perform_fix_and_propagate( sol, local_obj, res, true, true,
                                                false, true );

          if( local_obj < *current_obj_value )
             *current_obj_value = local_obj;
          std::copy( res.begin(), res.end(), result );
          return !res.empty();
#ifdef PAPILO_TBB
       } );
#endif
}

int
call_simple_heuristic( void* heuristic_void_ptr, double* result,
                       double* current_obj_value )
{
#ifdef PAPILO_TBB
   tbb::task_arena arena( 8 );
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
