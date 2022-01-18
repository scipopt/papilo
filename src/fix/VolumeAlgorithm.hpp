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

#include "fix/VectorMultiplication.hpp"
#include "papilo/core/Presolve.hpp"
#include "papilo/core/Problem.hpp"
#include "papilo/io/Message.hpp"
#include "papilo/misc/Num.hpp"
#include <cassert>
#include <cmath>
#include <fstream>

#include <cassert>
#include <fstream>

namespace papilo
{

template <typename REAL>
class VolumeAlgorithm
{
   Message msg;
   Num<REAL> num;
   VectorMultiplication<REAL> op;

   REAL alpha = 0.5;
   REAL f = 1;
   REAL con_threshold = 0.01;
   REAL obj_threshold = 0.02;

 public:
   VolumeAlgorithm( Message _msg, Num<REAL> _num, REAL _alpha, REAL _f )
       : msg( _msg ), num( _num ), alpha( _alpha ), f( _f )
   {
      op = {};
   }

   // TODO: define data structure
   /***
    * minimize cx s.t. Ax = b, Dx = e, x ≥ 0.
    * @param c
    * @param A
    * @param b
    * @param problem min 0 s.t. Dx = e, x ≥ 0.
    * @param pi
    * @return
    */
   Vec<REAL>
   volume_algorithm( const Vec<REAL> c, const ConstraintMatrix<REAL>& A,
                     const Vec<REAL>& b, Problem<REAL> problem,
                     const Vec<REAL> pi )
   {

      // TODO: define/determine UB
      REAL best_bound_on_obj = 0;
      REAL n_rows_A = A.getNRows();

      // Step 0
      // We start with a vector π̄ and solve (6) to obtain x̄ and z̄.
      std::pair<Vec<REAL>, REAL> sol =
          create_problem_6_and_solve_it( c, A, b, problem, pi );

      // Set x_0 = x_bar, z_0 = z_bar̄, t = 1
      int counter = 1;
      int non_improvement_iter_counter = 0;
      Vec<REAL> v_t( b );
      Vec<REAL> pi_t( pi );
      Vec<REAL> pi_bar( pi );
      Vec<REAL> residual_t( b );
      REAL best_objective = sol.second;

      while( stopping_criteria( v_t, n_rows_A, c, sol.first, best_objective ) )
      {
         // STEP 1:
         // Compute v_t = b − A x_bar and π_t = pi_bar + sv_t for a step size s given by (7).
         op.calc_b_minus_Ax( A, sol.first, b, v_t );
         REAL step_size = f * ( best_bound_on_obj - best_objective ) /
                          pow( op.l2_norm( v_t ), 2.0 );
         op.calc_b_plus_sx( pi_bar, step_size, v_t, pi_t );
         // Solve (6) with π_t , let x_t and z_t be the solutions obtained.
         std::pair<Vec<REAL>, REAL> sol_t =
             create_problem_6_and_solve_it( c, A, b, problem, pi_t );
         // x_bar ← αx_t + (1 − α)x_bar,
         op.calc_qb_plus_sx( alpha, sol_t.first, 1 - alpha, sol.first,
                             sol.first );

         // Step 2:
         // If z_t > z_bar update π_bar and z_bar̄ as
         if( num.isGT( sol_t.second, best_objective ) )
         {
            // π̄ ← π t , z̄ ← z t .
            best_objective = sol_t.second;
            pi_bar = pi_t;

            // If d (= v_t . (b - A x_t)) >= 0, then f = 1.1 * f
            op.calc_b_minus_Ax( A, sol_t.first, b, residual_t )
            if( num.isGE( op.multi( v_t, residual_t ), 0.0 ) )
               f = 1.1 * f;
            // TODO: need to verify if f <= 2?
         }
         else
         {
            if( ++non_improvement_iter_counter >= 20 )
               f = 0.66 * f;
         }
         // Let t ← t + 1 and go to Step 1.
         counter = counter + 1;
      }
      // TODO: return the list of x_t
      return sol.first;
   }

 private:
   bool
   stopping_criteria( const Vec<REAL>& v, const REAL n_rows_A,
                      const Vec<REAL>& c, const Vec<REAL>& x_bar,
                      const REAL z_bar )
   {
      return num.isGE( op.l1_norm( v ), n_rows_A * con_threshold ) ||
          num.isGE( ( op.multi( c, x_bar ) - z_bar ), z_bar * obj_threshold );
   }

   std::pair<Vec<REAL>, REAL>
   create_problem_6_and_solve_it( const Vec<REAL>& c,
                                  const ConstraintMatrix<REAL>& A,
                                  const Vec<REAL>& b, Problem<REAL>& problem,
                                  Vec<REAL> pi )
   {
      // TODO: z = (c − π̄ A)x + π̄b. π̄ is a transposed vector?
      //      problem.getObjective().coefficients = op.calc_b_minus_Ax(c, A,
      //      pi);
      Vec<REAL> vector( c );
      op.calc_b_minus_xA( A, pi, c, vector );
      problem.getObjective().coefficients = vector;
      problem.getObjective().offset = op.multi( b, pi );
      // TODO: the problem has now a new objective function and could be
      // presolved
      // TODO: missing
      return { pi, 0 };
   }
};

} // namespace papilo
