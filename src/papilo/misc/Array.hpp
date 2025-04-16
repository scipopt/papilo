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

#ifndef _PAPILO_MISC_ARRAY_HPP_
#define _PAPILO_MISC_ARRAY_HPP_

#include "papilo/misc/Alloc.hpp"
#include <cstdint>
#include <memory>

namespace papilo
{

template <typename T>
struct ArrayDeleter
{
   std::size_t size;

   ArrayDeleter( std::size_t _size ) : size( _size ) {}

   void
   operator()( T* p )
   {
      Allocator<T>().deallocate( p, size );
   }
};

template <typename T>
class Array
{
 public:
   Array( std::size_t n )
       : ptr( Allocator<T>().allocate( n ), ArrayDeleter<T>( n ) )
   {
   }

   std::size_t
   getSize() const
   {
      return ptr.get_deleter().size;
   }

   T&
   operator[]( int i )
   {
      return ptr[i];
   }

   const T&
   operator[]( int i ) const
   {
      return ptr[i];
   }

 private:
   std::unique_ptr<T[], ArrayDeleter<T>> ptr;
};

} // namespace papilo

#endif
