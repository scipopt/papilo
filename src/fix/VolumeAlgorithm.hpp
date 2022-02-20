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
#include "fix/Constraint.hpp"
#include "fix/AlgorithmParameter.hpp"
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
   AlgorithmParameter& parameter;
   PostsolveStorage<REAL>& postsolve_storage;
   REAL alpha;
   REAL alpha_max;
   REAL f;

 public:
   VolumeAlgorithm( Message _msg, Num<REAL> _num, Timer t,
                    AlgorithmParameter& parameter_,
                    PostsolveStorage<REAL>& postsolve_storage_ )
       : msg( _msg ), num( _num ), timer( t ), parameter( parameter_ ),
         postsolve_storage( postsolve_storage_ ), op( {} )
   {
      alpha = parameter.alpha;
      alpha_max = parameter.alpha_max;
      f = parameter.f;
   }

   /**
    * minimize cx s.t. Ax = b, Dx = e (D = empty), x ≥ 0.
    * @param c objective function
    * @param A equation or at least one finte bound for every constraint
    * @param derived_conflicts
    * @param b rhs of constraints
    * @param domains variables domains (lb/ub/flags)
    * @param pi initial dual multiplier
    * @param init_upper_bound value >= optimal objective value (c^T x)
    * @param init_primal_sol whether an initial primal solution was found
    * @return
    */
   Vec<REAL>
   volume_algorithm( const Vec<REAL> c, const ConstraintMatrix<REAL>& A,
                     const Vec<Constraint<REAL>>& derived_conflicts,
                     const Vec<REAL>& b, const VariableDomains<REAL>& domains,
                     const int num_int_vars, const REAL init_upper_bound,
                     const bool init_primal_sol, const int n_hard_constraints,
                     Vec<REAL>& pi, Vec<REAL>& pi_conflicts )
   {
      REAL st = timer.getTime();
      int n_rows_A = A.getNRows();
      int n_conflicts = derived_conflicts.size();

      if( !n_rows_A )
      {
         msg.error( "\t\tZero constraints detected! Volume algorithm "
                   "cannot handle this case currently.\n" );
         assert( false );
      }

      assert( pi.size() == n_rows_A );
      assert( pi_conflicts.size() == n_conflicts );

      assert( constraints_E_or_GE( n_rows_A, A, n_conflicts, derived_conflicts ) );

      Vec<REAL> b_conflicts( n_conflicts );
      build_conflict_data( n_conflicts, derived_conflicts, b_conflicts );

      // Step 0
      // Set x_0 = x_bar, z_0 = z_bar, t = 1
      int counter = 0;
      bool improvement_indicator = false;
      int weak_improvement_iter_counter = 0;
      int non_improvement_iter_counter = 0;
      Vec<REAL> v_t( b );
      Vec<REAL> viol_t( b );
      Vec<REAL> v_t_conflicts( n_conflicts );
      Vec<REAL> viol_t_conflicts( n_conflicts );
      Vec<REAL> x_t( c );
      Vec<REAL> pi_t( pi );
      Vec<REAL> pi_t_conflicts( pi_conflicts );
      Vec<REAL> pi_bar_conflicts( pi_conflicts );
      update_pi( n_rows_A, A, n_conflicts, derived_conflicts, pi_t,
                 pi_t_conflicts );
      Vec<REAL> pi_bar( pi_t );
      Vec<REAL> residual_t( b );
      Vec<REAL> residual_t_conflicts( n_conflicts );

      // We start with a vector π̄ and solve (6) to obtain x̄ and z̄.
      REAL z_bar = create_problem_6_and_solve_it( c, A, b, domains,
            derived_conflicts, b_conflicts, pi_t, pi_conflicts, x_t );
      Vec<REAL> x_bar( x_t );
      REAL z_bar_old = z_bar;
      // TODO: move away from box bound if possible
      REAL upper_bound_reset_val = num.isZero( init_upper_bound ) ? 0.1 :
         ( init_primal_sol || num.isLT( init_upper_bound, REAL{ 1.0 } ) ) ?
         init_upper_bound : 1.0;
      REAL upper_bound;
      bool finite_upper_bound = false;

      Vec<REAL> x_bar_last_iter( x_bar );
      Vec<int> fixed_int_vars_count( x_bar.size() );
      Vec<int> num_fixed_int_vars( parameter.max_iterations );
      int fixed_int_var_check_counter = 0;
      init_fixed_int_count( x_bar, domains, fixed_int_vars_count );

      op.calc_b_minus_Ax( A, x_bar, b, derived_conflicts, v_t, v_t_conflicts );
      calc_violations( n_rows_A, A, pi_bar, v_t, n_conflicts,
                       derived_conflicts, pi_bar_conflicts, v_t_conflicts,
                       viol_t, viol_t_conflicts );

      while( stopping_criteria( viol_t, n_rows_A, n_hard_constraints,
                                viol_t_conflicts, n_conflicts,
                                c, x_bar, z_bar, num_int_vars,
                                num_fixed_int_vars, counter ) )
      {
         msg.detailed( "Round of volume algorithm: {}\n", counter + 1 );
         // STEP 1:
         // Compute v_t = b − A x_bar and π_t = pi_bar + sv_t for a step size s
         // given by (7).
         update_upper_bound( z_bar, upper_bound_reset_val, init_primal_sol,
                             upper_bound, finite_upper_bound );
         assert( num.isGT( upper_bound, z_bar ) );
         // TODO: what is norm( v_t ) = 0?
         REAL step_size = f * ( upper_bound - z_bar ) /
                          ( pow( op.l2_norm( v_t ), 2.0 ) +
                            pow( op.l2_norm( v_t_conflicts ), 2.0 ) );
         msg.debug( "   Step size: {}\n", step_size );
         op.calc_b_plus_sx( pi_bar, step_size, v_t, pi_t );
         op.calc_b_plus_sx( pi_bar_conflicts, step_size, v_t_conflicts,
                            pi_t_conflicts );
         update_pi( n_rows_A, A, n_conflicts, derived_conflicts, pi_t,
                    pi_t_conflicts );

         // Solve (6) with π_t , let x_t and z_t be the solutions obtained.
         REAL z_t = create_problem_6_and_solve_it( c, A, b, domains,
                                    derived_conflicts, b_conflicts, pi,
                                    pi_conflicts, x_t );

         // Update alpha
         op.calc_b_minus_Ax( A, x_t, b, derived_conflicts, residual_t,
                             residual_t_conflicts );
         calc_alpha( residual_t, v_t, residual_t_conflicts, v_t_conflicts );

         x_bar_last_iter = x_bar;
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
            pi_bar_conflicts = pi_t_conflicts;
         }
         else
            improvement_indicator = false;

         update_fixed_int_count( x_bar, x_bar_last_iter, domains,
                                 counter,
                                 fixed_int_vars_count,
                                 num_fixed_int_vars );

         op.calc_b_minus_Ax( A, x_bar, b, derived_conflicts, v_t, v_t_conflicts );
         calc_violations( n_rows_A, A, pi_bar, v_t, n_conflicts,
                          derived_conflicts, pi_bar_conflicts, v_t_conflicts,
                          viol_t, viol_t_conflicts );

         // Update f
         update_f( improvement_indicator, v_t, residual_t,
                   v_t_conflicts, residual_t_conflicts,
                   weak_improvement_iter_counter,
                   non_improvement_iter_counter );

         // Update z_bar_old if needed
         if( counter && ( counter % 100 == 0 ) )
         {
            update_alpha_max( z_bar, z_bar_old );
            z_bar_old = z_bar;
         }

         // check integrality of x_bar
         integrality_check( x_bar, x_bar_last_iter, domains );

         // Let t ← t + 1 and go to Step 1.
         counter = counter + 1;
      };

      if( msg.getVerbosityLevel() >= VerbosityLevel::kInfo )
      {
         msg.info( "\t\tVol. alg. iterations: {} ( {} )\n", counter,
               parameter.max_iterations );
         msg.info( "\t\tAvg. (easy) constraint violation: {} ( {} )\n",
               ( op.l1_norm( viol_t ) + op.l1_norm( viol_t_conflicts ) ) /
                 ( n_rows_A + n_conflicts ), parameter.con_abstol );
         msg.info( "\t\tPrimal absolute objective value: {} ( {} )\n",
               calculate_orig_obj_value( x_bar ), parameter.obj_abstol );
         msg.info( "\t\tDuality gap: {} ( {} )\n", abs( op.multi( c, x_bar ) -
                  z_bar ) / abs( z_bar ), parameter.obj_reltol );
         int num_iters_check = parameter.num_iters_fixed_int_vars;
         /*
         msg.info( "\t\tNumber of fixed int vars: {} ( {} )\n",
               std::count_if( fixed_int_vars_count.begin(),
                  fixed_int_vars_count.end(),
                  [num_iters_check] ( int val ) { return val > num_iters_check;
                  } ), num_int_vars * parameter.fixed_int_var_threshold );
         */
         msg.info( "\t\tTotal time: {} ( {} )\n", timer.getTime() - st,
               parameter.time_limit );
      }

      pi = pi_bar;
      pi_conflicts = pi_bar_conflicts;

      return x_bar;
   }

 private:
   bool
   constraints_E_or_GE( const int n_rows_A, const ConstraintMatrix<REAL>& A,
              const int n_conflicts,
              const Vec<Constraint<REAL>>& derived_conflicts )
   {
      for( int i = 0; i < n_rows_A; i++ )
      {
         if( A.getRowFlags()[i].test( RowFlag::kLhsInf ) )
            return false;
      }

      for( int i = 0; i < n_conflicts; i++ )
      {
         if( derived_conflicts[i].get_row_flag().test( RowFlag::kLhsInf ) )
            return false;
      }

      return true;
   }

   void
   build_conflict_data( const int n_conflicts,
                        const Vec<Constraint<REAL>>& derived_conflicts,
                        Vec<REAL>& b_conflicts )
   {
      assert( b_conflicts.size() == n_conflicts );

#ifdef PAPILO_TBB
      tbb::parallel_for( tbb::blocked_range<int>( 0, derived_conflicts.size() ),
                         [&]( const tbb::blocked_range<int>& r )
                         {
                            for( int i = r.begin(); i < r.end(); ++i )
#else
      for( int i = 0; i < derived_conflicts.size(); ++i )
#endif
                            {
                               b_conflicts[i] = derived_conflicts[i].get_lhs();
                            }
#ifdef PAPILO_TBB
                         } );
#endif
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
              const int n_conflicts,
              const Vec<Constraint<REAL>>& derived_conflicts,
              Vec<REAL>& pi, Vec<REAL>& pi_conflicts )
   {
      for( int i = 0; i < n_rows_A; i++ )
      {
         if( A.getRowFlags()[i].test( RowFlag::kHardConstraint ) )
         {
            pi[i] = 0;
         }
         else if( A.getRowFlags()[i].test( RowFlag::kRhsInf ) )
         {
            // Note: change following max if assumption 4 is invalid.
            pi[i] = num.max( pi[i], REAL{ 0.0 } );
         }
      }

      for( int i = 0; i < n_conflicts; i++ )
      {
         if( derived_conflicts[i].get_row_flag().test( RowFlag::kRhsInf ) )
         {
            // Note: change following max if assumption 4 is invalid.
            pi_conflicts[i] = num.max( pi_conflicts[i], REAL{ 0.0 } );
         }
      }
   }

   bool
   stopping_criteria( const Vec<REAL>& v, const int n_rows_A,
                      const int n_hard_constraints,
                      const Vec<REAL>& v_conflicts, const int n_conflicts,
                      const Vec<REAL>& c, const Vec<REAL>& x_bar,
                      const REAL z_bar, const int num_int_vars,
                      const Vec<int>& num_fixed_int_vars,
                      const int iter_counter )
   {
      bool primal_feas_term = num.isLT( op.l1_norm( v ) +
                                        op.l1_norm( v_conflicts ),
                                        ( n_rows_A - n_hard_constraints +
                                          n_conflicts ) *
                                        parameter.con_abstol );

      bool duality_gap_abs_term = num.isLT( abs( op.multi( c, x_bar ) ),
                                           parameter.obj_abstol );
      bool duality_gap_rel_term = num.isLT( abs( op.multi( c, x_bar ) - z_bar ),
                                            abs( z_bar ) * parameter.obj_reltol
                                          );
      bool duality_gap_term = num.isZero( z_bar ) ? duality_gap_abs_term :
                                                    duality_gap_rel_term;

      // fixed integer variables with integer values
      int num_iters_fixed = parameter.num_iters_fixed_int_vars;
      int num_iters_percent = parameter.num_iters_fixed_int_vars_percent;
      bool fixed_int_var_term = false;
      int max_fixed = -1;
      int min_fixed = -1;
      if( iter_counter >= std::max( num_iters_fixed, num_iters_percent ) )
      {
         max_fixed = *std::max_element(
               std::next( num_fixed_int_vars.begin(), iter_counter -
                  num_iters_percent + 1 ),
               std::next( num_fixed_int_vars.begin(), iter_counter + 1 ) );
         min_fixed = *std::min_element(
               std::next( num_fixed_int_vars.begin(), iter_counter -
                  num_iters_percent + 1 ),
               std::next( num_fixed_int_vars.begin(), iter_counter + 1 ) );

         // TODO: OK even if max_fixed = 0 = min_fixed?
         fixed_int_var_term = ( (max_fixed - min_fixed ) <= num_int_vars *
                                         parameter.fixed_int_var_threshold );
      }

      msg.detailed( "   cons: {}\n", op.l1_norm( v ) + op.l1_norm( v_conflicts )
                                     / ( n_rows_A - n_hard_constraints +
                                         n_conflicts ) );
      msg.detailed( "   zbar: {}\n", z_bar );
      msg.detailed( "   objA: {}\n", abs( op.multi( c, x_bar ) ) );
      msg.detailed( "   objR: {}\n", abs( op.multi( c, x_bar ) - z_bar ) /
                                        abs( z_bar ) );
      msg.detailed( "   max_fixed: {}\t min_fixed: {}\t threshold: {}\n",
                        max_fixed, min_fixed,
                        num_int_vars * parameter.fixed_int_var_threshold );

      bool time_limit_term = ( num.isGE( timer.getTime(),
                                         parameter.time_limit ) );

      bool iter_limit_term = ( iter_counter >= parameter.max_iterations );

      return !( (primal_feas_term && duality_gap_term ) ||
                fixed_int_var_term ||
                time_limit_term ||
                iter_limit_term );
   }

   REAL
   create_problem_6_and_solve_it( const Vec<REAL>& c,
                                  const ConstraintMatrix<REAL>& A,
                                  const Vec<REAL>& b,
                                  const VariableDomains<REAL>& domains,
                                  const Vec<Constraint<REAL>>& derived_conflicts,
                                  const Vec<REAL>& b_conflicts,
                                  const Vec<REAL>& pi,
                                  const Vec<REAL>& pi_conflicts,
                                  Vec<REAL>& solution )
   {
      Vec<REAL> updated_objective( c );
      op.calc_b_minus_xA( A, pi, c, updated_objective );
      op.calc_b_minus_xA( derived_conflicts, pi_conflicts, updated_objective,
                          updated_objective );
      StableSum<REAL> obj_value{};
      obj_value.add( op.multi( b, pi ) );
      obj_value.add( op.multi( b_conflicts, pi_conflicts )  );

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
               // TODO: how to deal with this?
               return std::numeric_limits<REAL>::min();
            solution[i] = domains.lower_bounds[i];
         }
         else
         {
            if( domains.flags[i].test( ColFlag::kUbInf ) )
               // TODO: how to deal with this?
               return std::numeric_limits<REAL>::min();
            solution[i] = domains.upper_bounds[i];
         }
         obj_value.add( updated_objective[i] * solution[i] );
      }

      msg.debug( "   opt_val: {}\n", obj_value.get() );
      return obj_value.get();
   }

   void
   init_fixed_int_count( const Vec<REAL>& x_bar,
                         const VariableDomains<REAL>& domains,
                         Vec<int>& fixed_int_vars_count )
   {
      int x_bar_size = x_bar.size();
      assert( fixed_int_vars_count.size() == x_bar_size );

      for( int i = 0; i < x_bar_size; i++ )
      {
         if( domains.flags[i].test( ColFlag::kIntegral ) &&
               num.isIntegral( x_bar[i] ) )
            fixed_int_vars_count[i] = 1;
      }
   }

   // TODO: create an array for int var indices only once in Algorithm.hpp and
   //       pass it to the volume algo call?
   void
   update_fixed_int_count( const Vec<REAL>& x_bar,
                           const Vec<REAL> x_bar_last_iter,
                           const VariableDomains<REAL>& domains,
                           const int iter_counter,
                           Vec<int>& fixed_int_vars_count,
                           Vec<int>& num_fixed_int_vars )
   {
      int x_bar_size = x_bar.size();
      assert( fixed_int_vars_count.size() == x_bar_size);

      for( int i = 0; i < x_bar_size; i++ )
      {
         if( domains.flags[i].test( ColFlag::kIntegral ) &&
               num.isIntegral( x_bar[i] ) &&
               num.isEq( x_bar[i], x_bar_last_iter[i] ) )
            fixed_int_vars_count[i]++;
         else
            fixed_int_vars_count[i] = 0;
      }

      // number of integer variables with fixed integer values
      int num_iters_fixed_check = parameter.num_iters_fixed_int_vars;
      int n_fixed_int_vars = std::count_if( fixed_int_vars_count.begin(),
                                          fixed_int_vars_count.end(),
                                          [num_iters_fixed_check]( int val )
                                          { return val >= num_iters_fixed_check;
                                          } );
      num_fixed_int_vars[iter_counter] = n_fixed_int_vars;
   }

   void
   calc_violations( const int n_rows_A, const ConstraintMatrix<REAL>& A,
                    const Vec<REAL>& pi, const Vec<REAL>& residual,
                    const int n_conflicts,
                    const Vec<Constraint<REAL>>& derived_conflicts,
                    const Vec<REAL>& pi_conflicts,
                    const Vec<REAL>& residual_conflicts,
                    Vec<REAL>& viol_residual,
                    Vec<REAL>& viol_residual_conflicts )
   {
      viol_residual = residual;
      for( int i = 0; i < n_rows_A; i++ )
      {
         // Note: isZero check would be different in case of non-zero LB on pi
         if( A.getRowFlags()[i].test( RowFlag::kHardConstraint ) ||
             ( A.getRowFlags()[i].test( RowFlag::kRhsInf ) &&
               ( num.isLT( residual[i], REAL{ 0.0 } ) && num.isZero( pi[i] ) )
             ) )
            viol_residual[i] = 0;
      }

      viol_residual_conflicts = residual_conflicts;
      for( int i = 0; i < n_conflicts; i++ )
      {
         // Note: isZero check would be different in case of non-zero LB on pi
         if( derived_conflicts[i].get_row_flag().test( RowFlag::kRhsInf ) &&
             ( num.isLT( residual_conflicts[i], REAL{ 0.0 } ) &&
               num.isZero( pi_conflicts[i] ) ) )
            viol_residual_conflicts[i] = 0;
      }
   }

   void
   update_upper_bound( const REAL z_bar, const REAL upper_bound_reset_val,
                       const bool init_primal_sol,
                       REAL& upper_bound, bool& finite_upper_bound )
   {
      if( finite_upper_bound )
      {
         // TODO: shall we make 0.05, 0.03, 0.06, 1.0 global params same as f_min?
         if( num.isGE( z_bar,
                  upper_bound - abs( upper_bound ) * 0.05 ) )
         {
            upper_bound = num.isZero( z_bar ) ? upper_bound_reset_val :
                          num.max( upper_bound + abs( upper_bound ) * 0.03,
                          z_bar + abs( z_bar ) * 0.06 );
            msg.debug( "   updated best bound: {}\n", upper_bound );
         }
      }
      else
      {
         upper_bound = ( init_primal_sol || num.isZero( z_bar ) ) ?
                       upper_bound_reset_val : z_bar + abs( z_bar ) * 0.06;
         finite_upper_bound = true;
         msg.debug( "   updated best bound: {}\n", upper_bound );
      }
      assert( finite_upper_bound );
   }

   void
   calc_alpha( const Vec<REAL>& residual_t, const Vec<REAL>& residual_bar,
               const Vec<REAL>& residual_t_conflicts,
               const Vec<REAL>& residual_bar_conflicts )
   {
      // alpha_opt = minimizer of || alpha * residual_t + ( 1 - alpha ) *
      //                               residual_bar ||
      REAL t_t_prod = op.multi( residual_t, residual_t ) +
                      op.multi( residual_t_conflicts, residual_t_conflicts );

      REAL t_bar_prod = op.multi( residual_t, residual_bar ) +
                        op.multi( residual_t_conflicts, residual_bar_conflicts );

      REAL bar_bar_prod = op.multi( residual_bar, residual_bar ) +
                          op.multi( residual_bar_conflicts, residual_bar_conflicts );

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

      msg.detailed( "   alpha_opt: {},\t alpha_max: {},\t alpha: {}\n",
                    alpha_opt, alpha_max, alpha );
   }

   void
   update_f( const bool improvement_indicator, const Vec<REAL>& v_t,
             const Vec<REAL>& residual_t, const Vec<REAL>& v_t_conflicts,
             const Vec<REAL>& residual_t_conflicts,
             int& weak_improvement_iter_counter,
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
         if( num.isGE( op.multi( v_t, residual_t ) +
                       op.multi( v_t_conflicts, residual_t_conflicts ),
                       REAL{ 0.0 } ) )
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
   integrality_check( const Vec<REAL>& x_bar,
                      const Vec<REAL>& x_bar_last_iter,
                      const VariableDomains<REAL>& domains )
   {
      if( msg.getVerbosityLevel() == VerbosityLevel::kDetailed )
      {
         int num_integral = 0;
         int num_fixed_int = 0;
         for( int i = 0; i < x_bar.size(); i++ )
         {
            if( domains.flags[i].test( ColFlag::kIntegral ) &&
                  num.isIntegral( x_bar[i] ) )
               {
                  num_integral++;
                  if( num.isEq( x_bar[i], x_bar_last_iter[i] ) )
                     num_fixed_int++;
               }
         }
         msg.detailed( "   numInt: {}\t numFixedInt: {}\n", num_integral,
               num_fixed_int );
      }
   }

   void
   update_alpha_max( const REAL z_bar, const REAL z_bar_old )
   {
      // TODO: change 0.01, 1e-5, and 2.0 as global params?
      if( ( ( num.isZero( z_bar ) && num.isZero( z_bar_old ) ) ||
              num.isLT( z_bar, z_bar_old + 0.01 * abs( z_bar_old ) ) ) &&
          num.isGE( alpha_max / 2.0, REAL{ 1e-4 } ) )
         alpha_max = alpha_max / 2.0;
   }

   REAL
   calculate_orig_obj_value( const Vec<REAL>& reduced_sol ) const
   {
      Solution<REAL> original_solution{};
      Solution<REAL> reduced_solution{ reduced_sol };
      Message quiet{};
      quiet.setVerbosityLevel( papilo::VerbosityLevel::kQuiet );
      Postsolve<REAL> postsolve{ quiet, num };
      auto status = postsolve.undo( reduced_solution, original_solution,
            postsolve_storage );
      return postsolve_storage.getOriginalProblem().computeSolObjective(
            original_solution.primal );
   }
};

} // namespace papilo
