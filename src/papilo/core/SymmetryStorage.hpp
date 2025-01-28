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
