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

#ifndef FIX_RANDOM_ROUNDING_STRATEGY_HPP
#define FIX_RANDOM_ROUNDING_STRATEGY_HPP

#include "fix/RandomGenerator.hpp"
#include "fix/strategy/RoundingStrategy.hpp"

template <typename REAL>
class RandomRoundingStrategy : public RoundingStrategy<REAL>
{

   const Num<REAL> num;
   RandomGenerator random;

 public:
   RandomRoundingStrategy( RandomGenerator random_, Num<REAL> num_ )
       : random( random_ ), num( num_ )
   {
   }

   void
   recompute_locks()
   {
   }

   Fixing<REAL>
   select_rounding_variable( const Vec<REAL>& cont_solution,
                             const ProbingView<REAL>& view ) override
   {
      Vec<int> remaining_unfixed_cols{};
      for( int i = 0; i < cont_solution.size(); i++ )
      {
         if( num.isIntegral( cont_solution[i] ) ||
             num.isEq( view.getProbingUpperBounds()[i],
                       view.getProbingLowerBounds()[i] ) ||
             !view.is_integer_variable( i ) ||
             !view.is_within_bounds( i, cont_solution[i] ) )
            continue;
         remaining_unfixed_cols.push_back( i );
      }
      if( remaining_unfixed_cols.empty() )
         return { -1, -1 };

      std::uniform_int_distribution<uint32_t> dist_variable(
          0, remaining_unfixed_cols.size() - 1 );
      std::uniform_int_distribution<uint32_t> dist_rounding( 0, 1 );
      int variable =
          remaining_unfixed_cols[random.get_random_int( dist_variable )];
      REAL value = -1;
      if( random.get_random_int( dist_rounding ) )
         value = num.epsCeil( cont_solution[variable] );
      else
         value = num.epsFloor( cont_solution[variable] );

      return { variable, value };
   }
};

#endif
