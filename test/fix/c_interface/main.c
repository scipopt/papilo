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

int
main( void )
{
   int result = 1;
   void* problem_ptr = setup( "./resources/api_test.mps", &result );
   assert( result == 0 );
   int n_cols = 3;
   double primal_solution[n_cols];
   double sol[n_cols];
   for( int i = 0; i < n_cols; i++ )
      primal_solution[i] = ( 1.0 + i ) / 10.0;
   int success = call_algorithm( problem_ptr, primal_solution, sol, n_cols );
   delete_problem_instance( problem_ptr );
   assert( success );
}