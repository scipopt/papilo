/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*               This file is part of the program and library                */
/*    PaPILO --- Parallel Presolve for Integer and Linear Optimization       */
/*                                                                           */
/* Copyright (C) 2020-22  Konrad-Zuse-Zentrum */
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

#include "fix/FixAndPropagateApi.h"
#include <assert.h>
#include <stdlib.h>

int
main( void )
{
   int result = 1;
   void* heuristic = setup( "./../../resources/api_test.mps", &result, 4, 0, 0 );
   assert( result == 0 );
   int n_cols = 3;
   double* primal_solution = malloc( n_cols * sizeof( int ) );
   double* sol = malloc( n_cols * sizeof( double ) );
   for( int i = 0; i < n_cols; i++ )
      primal_solution[i] = ( 1.0 + i ) / 10.0;
   double current_solution = 50;

   int generated_conflicts = 0;
   int success = call_algorithm( heuristic, primal_solution, sol, n_cols,
                                 &current_solution,0,  0, 0, 0, 1, 1, 10000, &generated_conflicts );

   // should find a better solution
   assert( success );
   assert(generated_conflicts == 0);
   assert( sol[0] == 0 );
   assert( sol[1] == 0 );
   assert( sol[2] == 1 );
   assert( current_solution == 9 );

   double* sol2 = malloc( n_cols * sizeof( double ) );

   success = call_algorithm( heuristic, primal_solution, sol2, n_cols,
                                 &current_solution, 1, 0, 0, 0, 1, 1, 10000, &generated_conflicts );
   assert( !success );
   assert( current_solution == 9 );
   assert(generated_conflicts == 0);
   assert( sol2[0] == 0 );
   assert( sol2[1] == 0 );
   assert( sol2[2] == 1 );

   sol2[0] = 1;
   current_solution = 10;
   perform_one_opt(heuristic, sol2, n_cols,
                    2, &current_solution, 1000 );

   assert(current_solution == 9);
   assert( sol2[0] == 0 );
   assert( sol2[1] == 0 );
   assert( sol2[2] == 1 );

   delete_problem_instance( heuristic );

}
