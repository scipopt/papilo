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

#include "papilo/io/Message.hpp"
#include "papilo/misc/Num.hpp"

namespace papilo
{

struct AlgorithmParameter
{
 public:
   // overall parameters
   double time_limit = 10 * 60;
   int threads = 8;
   int seed = 0;
   bool copy_conflicts_to_problem = true;
   int size_of_conflicts_to_be_copied = 0;

   // vol algorithm parameters
   bool threshold_hard_constraints_vary = true;
   double threshold_hard_constraints_incr_factor = 2;
   double threshold_hard_constraints = 1;
   double alpha = 0.5;
   double alpha_max = 0.1;
   double f = 0.2;
   double f_min = 0.0005;
   double f_max = 2;
   double f_strong_incr_factor = 2;
   double f_weak_incr_factor = 1.1;
   double f_decr_factor = 0.66;
   double obj_reltol = 0.01;
   double obj_abstol = 0.01;
   double con_abstol = 0.02;
   double fixed_int_var_threshold = 0.05;
   int num_iters_fixed_int_vars = 20;
   int num_iters_fixed_int_vars_percent = 10001;
   int weak_improvement_iter_limit = 2;
   int non_improvement_iter_limit = 20;
   int max_iterations = 10000;

   // fix and propagate parameters

   // conflict analysis parameters

   void
   addParameters( ParameterSet& paramSet )
   {
      paramSet.addParameter( "time_limit", "", time_limit, 0.0 );
      paramSet.addParameter( "seed", "", seed, 0.0 );
      paramSet.addParameter( "threads",
                             "maximal number of threads to use (0: automatic)",
                             threads, 0 );

      paramSet.addParameter(
          "vol.threshold_hard_constraints_vary",
          "whether to vary the threshold for hard constraints or not",
          threshold_hard_constraints_vary );
      paramSet.addParameter(
          "vol.threshold_hard_constraints_incr_factor",
          "multiplier for increasing the threshold for hard constraints",
          threshold_hard_constraints_incr_factor, 1.0, 10.0 );
      paramSet.addParameter(
          "vol.threshold_hard_constraints",
          "constraint for which "
          "max(abs(coeff))/min(abs(coeff)) != x are excluded",
          threshold_hard_constraints, 1.0, 10.0 );
      paramSet.addParameter( "vol.alpha",
                             "multiplier for the convex "
                             "combination of primal solutions",
                             alpha, 0, 1.0 );
      paramSet.addParameter( "vol.alpha_max",
                             "upper bound for the parameter "
                             "alpha",
                             alpha_max, 0, 1.0 );
      paramSet.addParameter( "vol.f",
                             "multiplier for evaluating the step "
                             "size",
                             f, 0.0, 2.0 );
      paramSet.addParameter( "vol.f_min", "lower bound for the parameter f",
                             f_min, 0.0, 1.0 );
      paramSet.addParameter( "vol.f_max", "upper bound for the parameter f",
                             f_max, 0.0, 2.0 );
      paramSet.addParameter( "vol.f_strong_incr_factor",
                             "multiplier for "
                             "varying the parameter f in green iterations",
                             f_strong_incr_factor, 1.0, 3.0 );
      paramSet.addParameter( "vol.f_weak_incr_factor",
                             "multiplier for "
                             "varying the parameter f in yellow iterations",
                             f_weak_incr_factor, 1.0, 2.0 );
      paramSet.addParameter( "vol.f_decr_factor", "multiplier for varying the "
                             "parameter f in red iterations", f_decr_factor,
                             0.0, 1.0 );
      paramSet.addParameter( "vol.obj_reltol", "relative tolerance for "
                             "duality gap", obj_reltol, 0.0, 1.0 );
      paramSet.addParameter( "vol.obj_abstol", "absolute tolerance for "
                             "duality gap", obj_abstol, 0.0, 1.0 );
      paramSet.addParameter( "vol.con_abstol", "absolute tolerance for average "
                             "primal feasibility", con_abstol, 0.0, 1.0 );
      paramSet.addParameter( "vol.fixed_int_var_threshold", "threshold for "
                             "change in percentage of integer variables with "
                             "fixed integer values", fixed_int_var_threshold,
                             0.0, 1.0 );
      paramSet.addParameter( "vol.num_iters_fixed_int_vars", "number of"
                             "iterations the integer variable values need to"
                             "be fixed",
                             num_iters_fixed_int_vars, 1, 22222 );
      paramSet.addParameter( "vol.num_iters_fixed_int_vars_percent", "number of"
                             "iterations over which the change in percentage "
                             "of integer variables with fixed integer values "
                             "should be less than fixed_int_var_threshold",
                             num_iters_fixed_int_vars_percent, 1, 22222 );
      paramSet.addParameter( "vol.weak_improvement_iter_limit", "number of "
                             "yellow iterations after which the parameter f "
                             "is updated",
                             weak_improvement_iter_limit, 1, 22 );
      paramSet.addParameter( "vol.non_improvement_iter_limit",
                             "number of "
                             "red iterations after which the parameter f is "
                             "updated",
                             non_improvement_iter_limit, 2, 44 );
      paramSet.addParameter( "vol.max_iterations",
                             "maximum number of volume algorithm iterations in"
                             "one round",
                             max_iterations, 0, 44444 );
      paramSet.addParameter(
          "copy_conflicts_to_problem",
          "should the conflicts be copied to the problem of fix & propagate",
          copy_conflicts_to_problem );
      paramSet.addParameter(
          "size_of_conflicts_to_be_copied",
          "if the conflict have a size of x should they copied to the current problem?",
          size_of_conflicts_to_be_copied );
   }
};

} // namespace papilo
