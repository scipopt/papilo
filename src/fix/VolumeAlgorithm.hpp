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

   REAL alpha;
   REAL alpha_max;
   REAL f;
   REAL f_min;
   REAL f_max;
   REAL f_incr_factor;
   REAL f_decr_factor;
   REAL obj_threshold;
   REAL con_threshold;
   int weak_improvement_iter_limit;
   int non_improvement_iter_limit;

 public:
   VolumeAlgorithm( Message _msg, Num<REAL> _num,
                    REAL _alpha, REAL _alpha_max,
                    REAL _f, REAL _f_min, REAL _f_max,
                    REAL _f_incr_factor, REAL _f_decr_factor,
                    REAL _obj_threshold, REAL _con_threshold,
                    int _weak_improvement_iter_limit,
                    int _non_improvement_iter_limit )
       : msg( _msg ), num( _num ),
         alpha( _alpha ), alpha_max( _alpha_max ),
         f( _f ), f_min( _f_min ), f_max( _f_max ),
         f_incr_factor( _f_incr_factor ), f_decr_factor( _f_decr_factor ),
         obj_threshold( _obj_threshold ), con_threshold( _con_threshold ),
         weak_improvement_iter_limit( _weak_improvement_iter_limit ),
         non_improvement_iter_limit( _non_improvement_iter_limit ), op( {} )
   {
   }

   /**
    * minimize cx s.t. Ax = b, Dx = e (D = empty), x ≥ 0.
    * @param c objective function
    * @param A equation or at least one finte bound for every constraint
    * @param b not needed
    * @param domains variables domains (lb/ub/flags)
    * @param pi initial dual multiplier
    * @param best_bound_on_obj max bound of c^T x
    * @return
    */
   Vec<REAL>
   volume_algorithm( const Vec<REAL> c, const ConstraintMatrix<REAL>& A,
                     const Vec<REAL>& b, const VariableDomains<REAL>& domains,
                     const Vec<REAL> pi, REAL best_bound_on_obj )
   {
      int n_rows_A = A.getNRows();

      // Step 0
      // Set x_0 = x_bar, z_0 = z_bar, t = 1
      int counter = 1;
      bool improvement_indicator = false;
      int weak_improvement_iter_counter = 0;
      int non_improvement_iter_counter = 0;
      Vec<REAL> v_t( b );
      Vec<REAL> x_t( c );
      Vec<REAL> pi_t( pi );
      Vec<REAL> pi_bar( pi );
      modify_pi( n_rows_A, A, pi_bar );
      Vec<REAL> residual_t( b );

      // We start with a vector π̄ and solve (6) to obtain x̄ and z̄.
      REAL z_bar = create_problem_6_and_solve_it( c, A, b, domains, pi, x_t );
      Vec<REAL> x_bar( x_t );

      do
      {
         msg.info( "Round of volume algorithm: {}\n", counter );
         // STEP 1:
         // Compute v_t = b − A x_bar and π_t = pi_bar + sv_t for a step size s
         // given by (7).
         op.calc_b_minus_Ax( A, x_bar, b, v_t );
         update_best_bound_on_obj( z_bar, best_bound_on_obj );
         REAL step_size =
             f * ( best_bound_on_obj - z_bar ) / pow( op.l2_norm( v_t ), 2.0 );
         msg.info( "   Step size: {}\n", step_size );
         op.calc_b_plus_sx( pi_bar, step_size, v_t, pi_t );
         modify_pi( n_rows_A, A, pi_t );

         // Solve (6) with π_t , let x_t and z_t be the solutions obtained.
         REAL z_t =
             create_problem_6_and_solve_it( c, A, b, domains, pi_t, x_t );

         // Update alpha
         op.calc_b_minus_Ax( A, x_t, b, residual_t );
         update_alpha( residual_t, v_t );

         // x_bar ← αx_t + (1 − α)x_bar
         op.calc_qb_plus_sx( alpha, x_t, 1 - alpha, x_bar, x_bar );

         // Step 2:
         // If z_t > z_bar update π_bar and z_bar
         if( num.isGT( z_t, z_bar ) )
         {
            improvement_indicator = true;

            // π̄ ← π t , z̄ ← z t .
            z_bar = z_t;
            pi_bar = pi_t;
         }
         else
            improvement_indicator = false;

         // Update f
         update_f( improvement_indicator, v_t, residual_t,
                   weak_improvement_iter_counter, non_improvement_iter_counter
                   );

         // Let t ← t + 1 and go to Step 1.
         counter = counter + 1;
      }
      while( stopping_criteria( v_t, n_rows_A, c, x_bar, z_bar ) );

      return x_bar;
   }

 private:
   // Assumptions:
   // 1. Minimization objective sense
   // 2. Variable lower bounds: x >= 0
   // 3. A constraint is either an = or >= type.
   // TODO: Simplify this function further upon finalzing assumptions
   void
   modify_pi( const int n_rows_A, const ConstraintMatrix<REAL>& A, Vec<REAL>& pi )
   {
      for( int i = 0; i < n_rows_A; i++ )
      {
         if( A.getRowFlags()[i].test( RowFlag::kRhsInf ) )
            pi[i] = num.max( pi[i], REAL{ 0.0 } );
      }
   }

   bool
   stopping_criteria( const Vec<REAL>& v, const int n_rows_A,
                      const Vec<REAL>& c, const Vec<REAL>& x_bar,
                      const REAL z_bar )
   {
      msg.info( "   sc_1: {}\n", op.l1_norm( v ) / n_rows_A );
      msg.info( "   sc_2: {}\n", abs( op.multi( c, x_bar ) - z_bar ) / z_bar );
      return num.isGE( op.l1_norm( v ), n_rows_A * con_threshold ) ||
             num.isGE( abs( op.multi( c, x_bar ) - z_bar ),
                       z_bar * obj_threshold );
   }

   REAL
   create_problem_6_and_solve_it( const Vec<REAL>& c,
                                  const ConstraintMatrix<REAL>& A,
                                  const Vec<REAL>& b,
                                  const VariableDomains<REAL>& domains,
                                  const Vec<REAL>& pi, Vec<REAL>& solution )
   {
      Vec<REAL> updated_objective( c );
      op.calc_b_minus_xA( A, pi, c, updated_objective );
      StableSum<REAL> obj_value {};
      obj_value.add(op.multi( b, pi ));

      for( int i = 0; i < updated_objective.size(); i++ )
      {
         if( num.isZero( updated_objective[i] ) )
         {
            solution[i] = domains.lower_bounds[i];
            continue;
         }
         else if( num.isGT( updated_objective[i], REAL{ 0.0 } ) )
         {
            if( domains.flags[i].test( ColFlag::kLbInf ) )
               return std::numeric_limits<REAL>::min();
            solution[i] = domains.lower_bounds[i];
         }
         else
         {
            if( domains.flags[i].test( ColFlag::kUbInf ) )
               return std::numeric_limits<REAL>::min();
            solution[i] = domains.upper_bounds[i];
         }
         obj_value.add( updated_objective[i] * solution[i]);
      }

      msg.info( "   opt_val: {}\n", obj_value.get() );
      return obj_value.get();
   }

   void
   update_best_bound_on_obj( const REAL z_bar, REAL& best_bound_on_obj )
   {
      if( num.isGE( z_bar, 0.95 * best_bound_on_obj ) )
      {
         best_bound_on_obj = 1.05 * z_bar;
         msg.info( "   increased best bound: {}\n", best_bound_on_obj );
      }
   }

   void
   update_alpha( const Vec<REAL>& residual_t, const Vec<REAL>& residual_bar )
   {
      //TODO: introduce some logic for varying alpha_max
      // alpha_opt = minimizer of || alpha * residual_t + ( 1 - alpha ) *
      //                               residual_bar ||
      REAL t_t_prod = op.multi( residual_t, residual_t );

      REAL t_bar_prod = op.multi( residual_t, residual_bar );

      REAL bar_bar_prod = op.multi( residual_bar, residual_bar );

      REAL alpha_opt = ( bar_bar_prod - t_bar_prod ) /
                       ( t_t_prod + bar_bar_prod - 2.0 * t_bar_prod );
      alpha = num.isLT( alpha_opt, REAL{ 0.0 } )
                  ? alpha_max / 10.0
                  : num.min( alpha_opt, alpha_max );
      msg.info( "   alpha: {}\n", alpha );
   }

   void
   update_f( const bool improvement_indicator, const Vec<REAL>& v_t,
             const Vec<REAL>& residual_t, int& weak_improvement_iter_counter,
             int& non_improvement_iter_counter )
   {
      // +1 for increase in f, 0 for no change in f, -1 for decrease in f
      int change_f = 0;

      if( improvement_indicator )
      {
         //  If d (= v_t . (b - A x_t)) >= 0, then increase f
         if( num.isGE( op.multi( v_t, residual_t ), REAL{ 0.0 } ) )
            change_f = 1;
         else
         {
            ++( weak_improvement_iter_counter );
            if( weak_improvement_iter_counter >= weak_improvement_iter_limit )
            {
               weak_improvement_iter_counter = 0;
               change_f = 1;
            }
         }
      }
      else
      {
         ++( non_improvement_iter_counter );
         if( non_improvement_iter_counter >= non_improvement_iter_limit )
         {
            non_improvement_iter_counter = 0;
            change_f = -1;
         }
      }

      if( change_f >= 1 )
      {
         f = num.min( f_incr_factor * f, f_max );
         msg.info( "   increased f: {}\n", f );
      }
      else if( change_f <= -1 && num.isGE( f_decr_factor * f, f_min ) )
      {
         f = f_decr_factor * f;
         msg.info( "   decreased f: {}\n", f );
      }
   }
};

} // namespace papilo
