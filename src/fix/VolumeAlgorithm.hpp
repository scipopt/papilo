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
#include "fix/VolumeAlgorithmParameter.hpp"
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

 private:
   Message msg;
   Num<REAL> num;
   VectorMultiplication<REAL> op;
   Timer timer;
   VolumeAlgorithmParameter<REAL>& parameter;
   REAL alpha;
   REAL alpha_max;
   REAL f;

 public:
   VolumeAlgorithm( Message _msg, Num<REAL> _num, Timer t,
                    VolumeAlgorithmParameter<REAL>& parameter_ )
       : msg( _msg ), num( _num ), timer( t ), parameter( parameter_ ), op( {} )
   {
      alpha = parameter.alpha;
      alpha_max = parameter.alpha_max;
      f = parameter.f;
   }

   /**
    * minimize cx s.t. Ax = b, Dx = e (D = empty), x ≥ 0.
    * @param c objective function
    * @param A equation or at least one finte bound for every constraint
    * @param b not needed
    * @param domains variables domains (lb/ub/flags)
    * @param pi initial dual multiplier
    * @param box_upper_bound max box bound of c^T x
    * @return
    */
   Vec<REAL>
   volume_algorithm( const Vec<REAL> c, const ConstraintMatrix<REAL>& A,
                     const Vec<REAL>& b, const VariableDomains<REAL>& domains,
                     const Vec<REAL>& pi, REAL box_upper_bound )
   {
      std::all_of( A.getRowFlags().begin(), A.getRowFlags().end(),
                   []( RowFlags row_flag ) { return row_flag.test(RowFlag::kEquation); } );

      int n_rows_A = A.getNRows();

//      assert_pi( n_rows_A, A, pi );

      // Step 0
      // Set x_0 = x_bar, z_0 = z_bar, t = 1
      int counter = 1;
      bool improvement_indicator = false;
      int weak_improvement_iter_counter = 0;
      int non_improvement_iter_counter = 0;
      Vec<REAL> v_t( b );
      Vec<REAL> viol_t( b );
      Vec<REAL> x_t( c );
      Vec<REAL> pi_t( pi );
      Vec<REAL> pi_bar( pi );
      Vec<REAL> residual_t( b );

      // We start with a vector π̄ and solve (6) to obtain x̄ and z̄.
      REAL z_bar = create_problem_6_and_solve_it( c, A, b, domains, pi, x_t );
      Vec<REAL> x_bar( x_t );
      REAL z_bar_old = z_bar;
      // TODO: ok?
      REAL upper_bound_reset_val = num.isGE( box_upper_bound, REAL{ 1.0 } ) ?
                                   1.0 : box_upper_bound;
      // TODO: how to set -inf (or a large negative value)?
      REAL upper_bound = -1e30;

      op.calc_b_minus_Ax( A, x_bar, b, v_t );
      calc_violations( n_rows_A, A, pi_t, v_t, viol_t );

      while( stopping_criteria( viol_t, n_rows_A, c, x_bar, z_bar ) )
      {
         msg.detailed( "Round of volume algorithm: {}\n", counter );
         // STEP 1:
         // Compute v_t = b − A x_bar and π_t = pi_bar + sv_t for a step size s
         // given by (7).
         update_upper_bound( z_bar, upper_bound_reset_val, upper_bound );
         REAL step_size = f * ( upper_bound - z_bar ) /
                          pow( op.l2_norm( viol_t ), 2.0 );
         msg.debug( "   Step size: {}\n", step_size );
         op.calc_b_plus_sx( pi_bar, step_size, viol_t, pi_t );
//         update_pi( n_rows_A, A, pi_t );

         // Solve (6) with π_t , let x_t and z_t be the solutions obtained.
         REAL z_t =
             create_problem_6_and_solve_it( c, A, b, domains, pi_t, x_t );

         // Update alpha
         op.calc_b_minus_Ax( A, x_t, b, residual_t );
         calc_alpha( residual_t, v_t );

         // x_bar ← αx_t + (1 − α)x_bar
         op.calc_qb_plus_sx( alpha, x_t, 1 - alpha, x_bar,
                             x_bar );


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
                   weak_improvement_iter_counter,
                   non_improvement_iter_counter );

         // Update z_bar_old if needed
         if( counter % 100 == 0 )
         {
            update_alpha_max( z_bar, z_bar_old );
            z_bar_old = z_bar;
         }

         op.calc_b_minus_Ax( A, x_bar, b, v_t );
         calc_violations( n_rows_A, A, pi_t, v_t, viol_t );

         // Let t ← t + 1 and go to Step 1.
         counter = counter + 1;
      };
      // TODO: ahoen@suresh -> overwrite pi with current pi to be able to warm
      // restart the algorithm?
      //TODO: print some more useful data @Suresh
      msg.info("\t\tVol alg performed {} rounds.\n", counter);
      return x_bar;
   }

 private:
   // Assumptions:
   // 1. Each pi_i is either free or >= 0.
   void
   assert_pi( const int n_rows_A, const ConstraintMatrix<REAL>& A,
              const Vec<REAL>& pi )
   {
      for( int i = 0; i < n_rows_A; i++ )
      {
         if( A.getRowFlags()[i].test( RowFlag::kRhsInf ) )
         {
            assert( !A.getRowFlags()[i].test( RowFlag::kLhsInf ) );
            // Note: add another assert for LB if assumption 1 is invalid.
         }
      }
   }

   // Assumptions:
   // 1. Minimization objective sense
   // 2. Variable lower bounds: x >= 0
   // 3. A constraint is either an = or >= type.
   // 4. All non-free dual variables pi are >= 0 (i.e., no general bounds
   //    such as lb_i <= pi_i <= ub_i).
   // TODO: simplify this function further upon finalzing assumptions
   void
   update_pi( const int n_rows_A, const ConstraintMatrix<REAL>& A,
              Vec<REAL>& pi )
   {
      for( int i = 0; i < n_rows_A; i++ )
      {
         if( A.getRowFlags()[i].test( RowFlag::kRhsInf ) )
         {
            // Note: change following max if assumption 4 is invalid.
            pi[i] = num.max( pi[i], REAL{ 0.0 } );
         }
      }
   }

   bool
   stopping_criteria( const Vec<REAL>& v, const int n_rows_A,
                      const Vec<REAL>& c, const Vec<REAL>& x_bar,
                      const REAL z_bar )
   {
      msg.detailed( "   sc_1: {}\n", op.l1_norm( v ) / n_rows_A );
      msg.detailed( "   sc_2: {}\n",
                num.isZero( z_bar )
                    ? abs( op.multi( c, x_bar ) )
                    : abs( op.multi( c, x_bar ) - z_bar ) / z_bar );
      return num.isGE( op.l1_norm( v ), n_rows_A * parameter.con_abstol ) ||
             ( num.isZero( z_bar )
                   ? num.isGE( abs( op.multi( c, x_bar ) ),
                               parameter.obj_abstol )
                   : num.isGE( abs( op.multi( c, x_bar ) - z_bar ),
                               z_bar * parameter.obj_reltol ) );
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
      StableSum<REAL> obj_value{};
      obj_value.add( op.multi( b, pi ) );

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
         obj_value.add( updated_objective[i] * solution[i] );
      }

      msg.debug( "   opt_val: {}\n", obj_value.get() );
      return obj_value.get();
   }

   void
   calc_violations( const int n_rows_A, const ConstraintMatrix<REAL>& A,
                    const Vec<REAL>& pi, const Vec<REAL>& residual,
                    Vec<REAL>& viol_residual )
   {
      viol_residual = residual;
      /*
      for( int i = 0; i < n_rows_A; i++ )
      {
         // Note: isZero check would be different in case of non-zero LB on pi
         if( A.getRowFlags()[i].test( RowFlag::kRhsInf ) &&
             ( num.isLT( residual[i], REAL{ 0.0 } ) && num.isZero( pi[i] ) ) )
            viol_residual[i] = 0;
      }
      */
   }

   void
   update_upper_bound( const REAL z_bar, const REAL upper_bound_reset_val,
                       REAL& upper_bound )
   {
      // TODO: shall we make 0.05, 0.03, 0.06, 1.0 global params same as f_min?
      if( num.isGE( z_bar,
                    upper_bound - abs( upper_bound ) * 0.05 ) )
      {
         // TODO: replace 1.0 with some logical value
         upper_bound = num.isZero( z_bar ) ? upper_bound_reset_val :
             num.max( upper_bound + abs( upper_bound ) * 0.03,
                      z_bar + abs( z_bar ) * 0.06 );
         msg.debug( "   increased best bound: {}\n", upper_bound );
      }
   }

   void
   calc_alpha( const Vec<REAL>& residual_t, const Vec<REAL>& residual_bar )
   {
      // alpha_opt = minimizer of || alpha * residual_t + ( 1 - alpha ) *
      //                               residual_bar ||
      REAL t_t_prod = op.multi( residual_t, residual_t );

      REAL t_bar_prod = op.multi( residual_t, residual_bar );

      REAL bar_bar_prod = op.multi( residual_bar, residual_bar );

      REAL alpha_opt = alpha_max;
      if( num.isGT( t_t_prod + bar_bar_prod - 2.0 * t_bar_prod, REAL{ 0.0 } ) )
         alpha_opt = ( bar_bar_prod - t_bar_prod ) /
                     ( t_t_prod + bar_bar_prod - 2.0 * t_bar_prod );

      if( num.isLT( alpha_opt, alpha_max / 10.0 ) )
         alpha = alpha_max / 10.0;
      else if( num.isGT( alpha_opt, alpha_max ) )
         alpha = alpha_max;
      else
         alpha = alpha_opt;
      /*
      alpha = num.isLT( alpha_opt, REAL{ 0.0 } )
                  ? alpha_max / 10.0
                  : num.min( alpha_opt, alpha_max );
      */
      msg.detailed( "   alpha_opt: {},\t alpha_max: {},\t alpha: {}\n",
                    alpha_opt, alpha_max, alpha );
   }

   void
   update_f( const bool improvement_indicator, const Vec<REAL>& v_t,
             const Vec<REAL>& residual_t, int& weak_improvement_iter_counter,
             int& non_improvement_iter_counter )
   {
      // +2 for strong increase in f
      // +1 for weak increase in f
      // 0 for no change in f
      // -1 for decrease in f
      int change_f = 0;

      if( improvement_indicator )
      {
         // If d (= v_t . (b - A x_t)) >= 0, then increase f
         // TODO: should this be different for ineq cons?
         if( num.isGE( op.multi( v_t, residual_t ), REAL{ 0.0 } ) )
            change_f = 2;
         else
         {
            ++( weak_improvement_iter_counter );
            if( weak_improvement_iter_counter >=
                parameter.weak_improvement_iter_limit )
            {
               weak_improvement_iter_counter = 0;
               change_f = 1;
            }
         }
      }
      else
      {
         ++( non_improvement_iter_counter );
         if( non_improvement_iter_counter >=
             parameter.non_improvement_iter_limit )
         {
            non_improvement_iter_counter = 0;
            change_f = -1;
         }
      }

      if( change_f == 2 )
      {
         f = num.min( parameter.f_strong_incr_factor * f,
                                parameter.f_max );
         msg.debug( "   increased f: {}\n", f );
      }
      else if( change_f == 1 )
      {
         f = num.min( parameter.f_weak_incr_factor * f,
                                parameter.f_max );
         msg.debug( "   increased f: {}\n", f );
      }
      else if( change_f <= -1 &&
               num.isGE( parameter.f_decr_factor * f,
                         parameter.f_min ) )
      {
         f = parameter.f_decr_factor * f;
         msg.debug( "   decreased f: {}\n", f );
      }
   }

   void
   update_alpha_max( const REAL z_bar, const REAL z_bar_old )
   {
      // TODO: change 0.01, 1e-5, and 2.0 as global params?
      if( num.isLT( z_bar, z_bar_old + 0.01 * abs( z_bar_old ) ) &&
          num.isGE( alpha_max, REAL{ 1e-5 } ) )
         alpha_max = alpha_max / 2.0;
   }
};

} // namespace papilo