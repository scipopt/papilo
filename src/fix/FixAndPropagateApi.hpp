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

#include "papilo/io/MpsParser.hpp"
#include <cassert>
#include <fstream>
#include <string>

using namespace papilo;

void*
setup( const char* filename, int* result )
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
   *result = 0;
   auto problem = new Problem<double>( prob.get() );
   problem->recomputeAllActivities();
   return problem;
   //   problem.recomputeAllActivities();
   //   Num<double> num{};
   //   Message msg{};
   //   auto fixAndPropagate = new FixAndPropagate<double>{ msg, num, problem };

   //   if( !fixAndPropagate->init( filename ) )
   //   {
   //      result = -1;
   //      return nullptr;
   //   }

   //   result = 0;
   //   return fixAndPropagate;
}

void
delete_problem_instance( void* problem_ptr )
{
   auto problem = (Problem<double>*)( problem_ptr );
   delete problem;
}

bool
call_algorithm( void* problem_ptr, double* cont_solution, double* result,
                int n_cols )
{
   //   auto f = (FixAndPropagate<double>*)( fix_ );
   auto problem = (Problem<double>*)( problem_ptr );
   ProbingView<double> view{ *problem, {} };
   FixAndPropagate<double> f{ {}, {}, *problem, view };
   Vec<double> sol( cont_solution, cont_solution + n_cols );
   Vec<double> res( result, result + n_cols );

   bool is_infeasible = f.fix_and_propagate( sol, res );
   result = &res[0];
   return is_infeasible;
}
