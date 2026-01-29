/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*               This file is part of the program and library                */
/*    PaPILO --- Parallel Presolve for Integer and Linear Optimization       */
/*                                                                           */
/* Copyright (C) 2020-2026 Zuse Institute Berlin (ZIB)                       */
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

#ifndef _PAPILO_MISC_SIGNATURE_HPP_
#define _PAPILO_MISC_SIGNATURE_HPP_

#include "papilo/misc/Hash.hpp"
#include <cmath>
#include <cstdint>
#include <type_traits>

namespace papilo
{

template <typename T>
class Signature
{
 public:
   Signature() : state( 0 ) {}

   template <typename U>
   void
   add( U elem )
   {
      state |=
          1 << ( ( uint32_t( elem ) *
                   HashHelpers<uint32_t>::fibonacci_muliplier() ) >>
                 ( 32 - static_cast<int>( std::log2( 8 * sizeof( T ) ) ) ) );
   }

   // template <typename U>
   // void
   // add( U elem )
   //{
   //   state |= 1 << ( uint32_t( elem ) % ( sizeof( T ) * 8 ) );
   //}

   bool
   isSubset( Signature other )
   {
      return ( state & ~other.state ) == 0;
   }

   bool
   isSuperset( Signature other )
   {
      return ( other.state & ~state ) == 0;
   }

   bool
   isEqual( Signature other )
   {
      return state == other.state;
   }

 private:
   T state;
};

using Signature32 = Signature<uint32_t>;
using Signature64 = Signature<uint64_t>;

} // namespace papilo

#endif
