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

#ifndef _PAPILO_MISC_COMPRESS_VECTOR_HPP_
#define _PAPILO_MISC_COMPRESS_VECTOR_HPP_

#include "papilo/misc/Vec.hpp"
#include <cassert>

namespace papilo
{

/// helper function to compress a vector-like container using the given mapping
template <typename VEC>
void
compress_vector( const Vec<int>& mapping, VEC& vec )
{
   assert( vec.size() == mapping.size() );

   int newSize = 0;
   for( int i = 0; i != static_cast<int>( vec.size() ); ++i )
   {
      assert( mapping[i] <= i );

      if( mapping[i] != -1 )
      {
         vec[mapping[i]] = vec[i];
         newSize++;
      }
   }
   vec.resize( newSize );
}

/// helper function to compress a vector-like container of indicies using the
/// given mapping
template <typename VEC>
void
compress_index_vector( const Vec<int>& mapping, VEC& vec )
{
   int offset = 0;
   for( std::size_t i = 0; i < vec.size(); ++i )
   {
      int newindex = mapping[vec[i]];
      if( newindex != -1 )
         vec[i - offset] = newindex;
      else
         ++offset;
   }

   vec.resize( vec.size() - offset );
}

} // namespace papilo

#endif
