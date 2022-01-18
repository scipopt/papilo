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
      // Step 0
      std::pair<Vec<REAL>, REAL> sol =
          create_problem_6_and_solve_it( c, A, b, problem, pi );
      Vec<REAL> pi_bar { pi };
      Vec<REAL> x_bar { sol.first };
      REAL z_bar = sol.second;
      int t = 1;
      Vec<REAL> v_t( b );
      REAL n_rows_A = A.getNRows();

      while( stopping_criteria(v_t, n_rows_A, c, x_bar, z_bar) )
      {
         // STEP 1:
         v_t = op.calc_b_minus_Ax( A, sol.first, b );
         REAL step_size =
             f * ( best_bound_on_obj - sol.second ) / op.l2_norm( v_t );
         Vec<REAL> pi_t = op.calc_b_plus_sx( pi_bar, step_size, v_t );

         std::pair<Vec<REAL>, REAL> sol_t =
             create_problem_6_and_solve_it( c, A, b, problem, pi_t );
         x_bar = op.calc_qb_plus_sx( alpha, sol_t.first, 1 - alpha, x_bar );

         // Step 2:
         if( num.isGT( sol_t.second, z_bar ) )
         {
            //TODO: make quick check if the value gets overwritten
            sol = sol_t;
            z_bar = sol_t.second;
         }
         t = t + 1;
      }
      return x_bar;
   }

 private:
   bool
   stopping_criteria( const Vec<REAL>& v, const REAL n_rows_A,
         const Vec<REAL>& c, const Vec<REAL>& x_bar, const Vec<REAL>& z_bar )
   {
      if( num.isLT(op.l1_norm(v), n_rows_A * con_threshold) &&
          num.isLT((op.multi(c, x_bar) - z_bar), z_bar * obj_threshold) )
         return false;

      return true;
   }

   std::pair<Vec<REAL>, REAL>
   create_problem_6_and_solve_it( const Vec<REAL> c,
                                  const ConstraintMatrix<REAL>& A,
                                  const Vec<REAL>& b, Problem<REAL> problem,
                                  Vec<REAL> pi )
   {
      // TODO: z = (c − π̄ A)x + π̄b. π̄ is a transposed vector?
      //      problem.getObjective().coefficients = op.calc_b_minus_Ax(c, A,
      //      pi);
      const Vec<REAL>& vector = op.calc_b_minus_xA( A, pi, c );
      problem.getObjective().coefficients = vector;
      problem.getObjective().offset = op.multi( b, pi );
      // TODO: the problem has now a new objective function and could be
      // presolved
      // TODO: missing
      return { pi, 0 };
   }
};

} // namespace papilo
