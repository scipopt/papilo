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

template <typename REAL>
struct VolumeAlgorithmParameter
{

   REAL alpha;
   REAL alpha_max;
   REAL f;
   REAL f_min;
   REAL f_max;
   REAL f_strong_incr_factor;
   REAL f_weak_incr_factor;
   REAL f_decr_factor;
   REAL obj_reltol;
   REAL obj_abstol;
   REAL con_abstol;
   int weak_improvement_iter_limit;
   int non_improvement_iter_limit;

   VolumeAlgorithmParameter( REAL _alpha, REAL _alpha_max, REAL _f, REAL _f_min,
                             REAL _f_max, REAL _f_strong_incr_factor,
                             REAL _f_weak_incr_factor, REAL _f_decr_factor,
                             REAL _obj_reltol, REAL _obj_abstol,
                             REAL _con_abstol, int _weak_improvement_iter_limit,
                             int _non_improvement_iter_limit )
       : alpha( _alpha ), alpha_max( _alpha_max ), f( _f ), f_min( _f_min ),
         f_max( _f_max ), f_strong_incr_factor( _f_strong_incr_factor ),
         f_weak_incr_factor( _f_weak_incr_factor ),
         f_decr_factor( _f_decr_factor ), obj_reltol( _obj_reltol ),
         obj_abstol( _obj_abstol ), con_abstol( _con_abstol ),
         weak_improvement_iter_limit( _weak_improvement_iter_limit ),
         non_improvement_iter_limit( _non_improvement_iter_limit )
   {
   }
};

} // namespace papilo
