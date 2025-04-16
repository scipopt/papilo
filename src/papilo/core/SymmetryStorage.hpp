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

#ifndef _PAPILO_CORE_SYMMETRIE_HPP_
#define _PAPILO_CORE_SYMMETRIE_HPP_

#include "papilo/core/SparseStorage.hpp"
#include "papilo/misc/Vec.hpp"
#include <cassert>
#include <cstdint>

namespace papilo
{

/// possible types of post solving
enum class SymmetryType : int
{
   kXgeY = 0,
   kXplusYge1 = 1,
};

struct Symmetry
{
//   REAL val;
   int dominating_col;
   int dominated_col;
   SymmetryType symmetryType;

   int
   getDominatingCol() const
   {
      return dominating_col;
   }

   int
   getDominatedCol() const
   {
      return dominated_col;
   }

   SymmetryType
   getSymmetryType() const
   {
      return symmetryType;
   }

   Symmetry() = default;

   Symmetry( int _dominating_col, int _dominated_col, SymmetryType _type )
       : dominating_col( _dominating_col ), dominated_col( _dominated_col ), symmetryType(_type)
   {
   }


   template <typename Archive>
   void
   serialize( Archive& ar, const unsigned int version )
   {
      ar& dominating_col;
      ar& dominated_col;
      ar& symmetryType;
   }
};

struct SymmetryStorage
{
   void
   addSymmetry( int dominating_dominated_col, int dominated_col, SymmetryType symmetryType )
   {
      symmetries.emplace_back( dominating_dominated_col, dominated_col, symmetryType );
   }

   Symmetry
   find_symmetry( int dominating_col, int dominated_col ) const
   {
      for(auto symmetry : symmetries)
      {
         if( ( symmetry.dominated_col == dominating_col &&
               dominated_col == symmetry.dominating_col ) ||
             ( symmetry.dominated_col == dominated_col &&
               dominating_col == symmetry.dominating_col ) )
            return symmetry;
      }
      return { -1, -1, SymmetryType::kXgeY };
   }

   bool
   contains_symmetry( int dominating_col, int dominated_col ) const
   {
      for(auto symmetry : symmetries)
      {
         if( ( symmetry.dominated_col == dominating_col &&
               dominated_col == symmetry.dominating_col ) ||
             ( symmetry.dominated_col == dominated_col &&
               dominating_col == symmetry.dominating_col ) )
            return true;
      }
      return false;
   }

   void
   reserve( int size )
   {
      symmetries.reserve( size );
   }

   SymmetryStorage()
   {
      symmetries = {};
   }

   void
   compress( const Vec<int>& colmapping, bool full = false )
   {
      int newSize = 0;
      for( int i = 0; i < static_cast<int>( symmetries.size() ); ++i )
      {
         if( colmapping[symmetries[i].getDominatingCol()] == -1 ||
             colmapping[symmetries[i].getDominatedCol()] == -1 )
            continue;
         symmetries[newSize] = { colmapping[symmetries[i].getDominatingCol()],
                                 colmapping[symmetries[i].getDominatedCol()],
                                 symmetries[i].getSymmetryType()};
         newSize++;
      }
      symmetries.resize( newSize );
      if( full )
         symmetries.shrink_to_fit();
   }


   template <typename Archive>
   void
   serialize( Archive& ar, const unsigned int version )
   {
      ar& symmetries;
   }

   Vec<Symmetry> symmetries;
};




} // namespace papilo

#endif
