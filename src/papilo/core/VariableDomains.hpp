/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*               This file is part of the program and library                */
/*    PaPILO --- Parallel Presolve for Integer and Linear Optimization       */
/*                                                                           */
/* Copyright (C) 2020  Konrad-Zuse-Zentrum                                   */
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

#ifndef _PAPILO_CORE_VARIABLE_DOMAINS_HPP_
#define _PAPILO_CORE_VARIABLE_DOMAINS_HPP_

#include "papilo/misc/Flags.hpp"
#include "papilo/misc/MultiPrecision.hpp"
#include "papilo/misc/Vec.hpp"
#include "papilo/misc/compress_vector.hpp"
#include "tbb/parallel_invoke.h"

namespace papilo
{

enum class ColFlag : uint8_t
{
   NONE = 0,
   LB_INF = 1 << 0,
   LB_HUGE = 1 << 1,
   UB_INF = 1 << 2,
   UB_HUGE = 1 << 3,
   INTEGRAL = 1 << 4,
   FIXED = 1 << 5,
   SUBSTITUTED = 1 << 6,
   IMPL_INT = 1 << 7,
   UNBOUNDED = static_cast<uint8_t>( ColFlag::LB_INF ) |
               static_cast<uint8_t>( ColFlag::UB_INF ),
   INACTIVE = static_cast<uint8_t>( ColFlag::FIXED ) |
              static_cast<uint8_t>( ColFlag::SUBSTITUTED ),
   LB_USELESS = static_cast<uint8_t>( ColFlag::LB_INF ) |
                static_cast<uint8_t>( ColFlag::LB_HUGE ),
   UB_USELESS = static_cast<uint8_t>( ColFlag::UB_INF ) |
                static_cast<uint8_t>( ColFlag::UB_HUGE ),
};

using ColFlags = Flags<ColFlag>;

/// Type to store the domains for the variables in a problem.
/// This includes the lower and upper bounds, and whether the
/// variable is constraint to integral values.
template <typename REAL>
struct VariableDomains
{
   Vec<REAL> lower_bounds;
   Vec<REAL> upper_bounds;
   Vec<ColFlags> flags;

   void
   compress( const Vec<int>& colmapping, bool full = false );

   bool
   isBinary( int col ) const
   {
      return flags[col].test( ColFlag::INTEGRAL ) &&
             !flags[col].test( ColFlag::LB_USELESS, ColFlag::UB_USELESS,
                               ColFlag::INACTIVE ) &&
             lower_bounds[col] == 0 && upper_bounds[col] == 1;
   }

   template <typename Archive>
   void
   serialize( Archive& ar, const unsigned int version )
   {
      ar& lower_bounds;
      ar& upper_bounds;
      ar& flags;
   }
};

#ifdef PAPILO_USE_EXTERN_TEMPLATES
extern template struct VariableDomains<double>;
extern template struct VariableDomains<Quad>;
extern template struct VariableDomains<Rational>;
#endif

template <typename REAL>
void
VariableDomains<REAL>::compress( const Vec<int>& colmapping, bool full )
{
   tbb::parallel_invoke(
       [this, &colmapping, full]() {
          compress_vector( colmapping, lower_bounds );
          if( full )
             lower_bounds.shrink_to_fit();
       },
       [this, &colmapping, full]() {
          compress_vector( colmapping, upper_bounds );
          if( full )
             upper_bounds.shrink_to_fit();
       },
       [this, &colmapping, full]() {
          compress_vector( colmapping, flags );
          if( full )
             flags.shrink_to_fit();
       } );
}
} // namespace papilo

#endif
