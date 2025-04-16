/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*               This file is part of the program and library                */
/*    PaPILO --- Parallel Presolve for Integer and Linear Optimization       */
/*                                                                           */
/* Copyright (C) 2020-2025 Zuse Institute Berlin (ZIB)                       */
/*                                                                           */
/* Licensed under the Apache License, Version 2.0 (the "License");           */
/* you may not use this file except in compliance with the License.          */
/* You may obtain a copy of the License at                                   */
/*                                                                           */
/*     http://www.apache.org/licenses/LICENSE-2.0                            */
/*                                                                           */
/* Unless required by applicable law or agreed to in writing, software       */
/* distributed under the License is distributed on an "AS IS" BASIS,         */
/* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  */
/* See the License for the specific language governing permissions and       */
/* limitations under the License.                                            */
/*                                                                           */
/* You should have received a copy of the Apache-2.0 license                 */
/* along with PaPILO; see the file LICENSE. If not visit scipopt.org.        */
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
