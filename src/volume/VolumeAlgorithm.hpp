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

#include "papilolib.h"
#include <cassert>
#include <fstream>
#include <cmath>

using namespace papilo;

template <typename REAL>
class VolumeAlgorithm
{
   Message msg;
   Num<REAL> num;

 public:
   VolumeAlgorithm( Message _msg, Num<REAL> _num ) : msg( _msg ), num( _num ) {}

   Vec<REAL>
   volume_algorithm( const Problem<REAL>& orig_prob,
                      Solution<REAL> dual_sol )
   {
      //TODO: create pi_bar from dual_sol
      Vec<double> pi_bar;

      //TODO: all constraints in equality form?
      //TODO: construct formulation (6) from given orig_prob and dual_sol
      //TODO: add correct nnz s
      PAPILO_PROBLEM* prob =
         papilo_problem_create( 1e30, "volume_algo", 0, 0, 2 );

      // solving (6) to obtain x_bar and z_bar
      PAPILO_SOLVER* solver = papilo_solver_create();
      papilo_solver_load_problem( solver, prob );
      PAPILO_SOLVING_INFO* result = papilo_solver_start( solver );

      assert( result->solve_result == PAPILO_SOLVE_RESULT_OPTIMAL );
      assert( result->bestsol != nullptr );
      //TODO: Vec<double> correct type of bestsol?
      Vec<double> x_bar = result->bestsol;
      double z_bar = result->bestsol_obj;

      VectorMultiplication<double> vec_mult;
      Vec<double> minus_v_t;
      double f;
      double UB;
      double v_norm;
      double s;
      Vec<double> pi_t;
      Vec<double> x_t;
      double z_t;
      double alpha;

      while ( true )
      {
         // computing -v_t
         //TODO: obtain A and b
         minus_v_t = vec_mult.multiplication( A, x_bar, b );

         // computing -pi_t 
         f = 1.0;
         //TODO: right value of UB?
         UB = 1e10;
         v_norm = l2_norm(minus_v_t);
         s = f * ( UB - z_bar ) / v_norm;
         pi_t = addition( pi_bar, 1.0, minus_v_t, -1.0 * s );

         //TODO: update prob's objective function using papilo_problem_change_col_obj
         //TODO: a different way to allow warm starting?
         papilo_solver_load_problem( solver, prob );
         PAPILO_SOLVING_INFO* result = papilo_solver_start( solver );

         assert( result->solve_result == PAPILO_SOLVE_RESULT_OPTIMAL );
         assert( result->bestsol != nullptr );
         //TODO: Vec<double> correct type of bestsol?
         x_t = result->bestsol;
         z_t = result->bestsol_obj;

         alpha = 0.5;
         x_bar = addition( x_t, alpha, x_bar, (1.0 - alpha) );

         //TODO: right function for this check?
         if ( z_t > z_bar )
         {
            pi_bar = pi_t;
            z_bar = z_t;
         }
      }

      papilo_problem_free( prob );

      return x_t;
   }

 private:
   double
   l2_norm( const Vec<REAL>& vec )
   {
      double squared_norm = 0.0;

      for( int i = 0; i < vec.size(); i++ )
      {
         squared_norm += vec[i] * vec[i];
      }

      return sqrt(squared_norm);
   }

   Vec<REAL>
   addition( const Vec<REAL>& vec1, const double mult1,  const Vec<REAL>& vec2,
         const double mult2 )
   {
      assert( vec1.size() == vec2.size() );
      Vec<REAL> result( vec1 );
#ifdef PAPILO_TBB
      tbb::parallel_for( tbb::blocked_range<int>( 0, vec1.size() ),
            [&]( const tbb::blocked_range<int>& r )
            {
               for( int i = r.begin(); i < r.end(); ++i )
#else
               for( int i = 0; i < vec1.size(); ++i )
#endif
               {
                  result[i] = mult1 * vec1[i] + mult2 * vec2[i];
               }
#ifdef PAPILO_TBB
            } );
#endif
      return result;
   }

};
