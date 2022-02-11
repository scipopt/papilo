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

   // vol algorithm parameters
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
   int weak_improvement_iter_limit = 2;
   int non_improvement_iter_limit = 20;

   // fix and propagate parameters

   // conflict analysis parameters

   void
   addParameters( ParameterSet& paramSet )
   {
      paramSet.addParameter( "vol.alpha", "@Suresh", alpha, 0, 2 );
      paramSet.addParameter( "vol.alpha_max", "@Suresh", alpha_max, 0, 2 );
      paramSet.addParameter( "vol.f", "@Suresh", f, 0.0, 2.0 );
      paramSet.addParameter( "vol.f_min", "@Suresh", f_min, 0.0, 2.0 );
      paramSet.addParameter( "vol.f_max", "@Suresh", f_max, 0.0, 2.0 );
      paramSet.addParameter( "vol.f_strong_incr_factor", "@Suresh",
                             f_strong_incr_factor, 0.0, 1.0 );
      paramSet.addParameter( "vol.f_weak_incr_factor", "@Suresh",
                             f_weak_incr_factor, 0.0, 1.0 );
      paramSet.addParameter( "vol.f_decr_factor", "@Suresh", f_decr_factor, 0.0,
                             1.0 );
      paramSet.addParameter( "vol.obj_reltol", "@Suresh", obj_reltol, 0.0,
                             1.0 );
      paramSet.addParameter( "vol.obj_abstol", "@Suresh", obj_abstol, 0.0,
                             1.0 );
      paramSet.addParameter( "vol.con_abstol", "@Suresh", con_abstol, 0.0,
                             1.0 );
      paramSet.addParameter( "vol.weak_improvement_iter_limit", "@Suresh",
                             weak_improvement_iter_limit, 0.0, 1.0 );
      paramSet.addParameter( "vol.non_improvement_iter_limit", "@Suresh",
                             non_improvement_iter_limit, 0.0, 1.0 );
      paramSet.addParameter( "vol.threshold_hard_constraints",
                             "constraint for which "
                             "max(abs(coeff))/max(abs(coeff)) > x are excluded",
                             threshold_hard_constraints, 0.0, 10.0 );
      paramSet.addParameter( "time_limit", "", time_limit, 0.0 );
      paramSet.addParameter( "threads",
                             "maximal number of threads to use (0: automatic)",
                             time_limit, 0.0 );
   }
};

} // namespace papilo
