/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*               This file is part of the program and library                */
/*    PaPILO --- Parallel Presolve for Integer and Linear Optimization       */
/*                                                                           */
/* Copyright (C) 2020-2025 Zuse Institute Berlin (ZIB)                       */
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

#ifndef _PAPILO_MISC_STABLE_SUM_HPP_
#define _PAPILO_MISC_STABLE_SUM_HPP_

#include "papilo/misc/Num.hpp"

namespace papilo
{

template <typename REAL, bool isfp = num_traits<REAL>::is_floating_point>
class StableSum;

template <typename REAL>
class StableSum<REAL, true>
{
   REAL sum = 0;
   REAL c = 0;

 public:
   StableSum() = default;

   explicit StableSum( const REAL& init ) : sum( init ), c( 0 ) {}

   void
   add( const REAL& input )
   {
      REAL t = sum + input;
      REAL z = t - sum;
      REAL y = ( sum - ( t - z ) ) + ( input - z );
      c += y;

      sum = t;
   }

   REAL
   get() const
   {
      return sum + c;
   }
};

template <typename REAL>
class StableSum<REAL, false>
{
   REAL sum = 0;

 public:
   StableSum() = default;

   explicit StableSum( const REAL& init ) : sum( init ) {}

   void
   add( const REAL& input )
   {
      sum += input;
   }

   REAL
   get() const
   {
      return sum;
   }
};

} // namespace papilo

#endif
