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

#include "fix/strategy/RoundingStrategy.hpp"

template <typename REAL>
class RandomRoundingStrategy : public RoundingStrategy<REAL>
{

   const Num<REAL> num;

   typedef std::mt19937 MyRNG;
   uint32_t seed;

   MyRNG random_generator;

 public:
   RandomRoundingStrategy( uint32_t seed_, Num<REAL> num_ )
       : seed( seed_ ), num( num_ )
   {
      random_generator.seed( seed );
   }

   Fixing<REAL>
   select_rounding_variable( const Vec<REAL>& cont_solution,
                             const ProbingView<REAL>& view ) override
   {
      // TODO: this does not work since fixed variable could be obtained

      Vec<int> remaining_unfixed_cols{};
      for( int i = 0; i < cont_solution.size(); i++ )
      {
         if( num.isEq( view.getProbingUpperBounds()[i],
                       view.getProbingLowerBounds()[i] ) ||
             !view.is_integer( i ) )
            continue;
         remaining_unfixed_cols.push_back( i );
      }
      if( remaining_unfixed_cols.empty() )
         return { -1, -1 };

      std::uniform_int_distribution<uint32_t> dist_variable(
          0, remaining_unfixed_cols.size() - 1 );
      std::uniform_int_distribution<uint32_t> dist_rounding( 0, 1 );
      int variable = remaining_unfixed_cols[dist_variable( random_generator )];
      REAL value = -1;
      if( dist_rounding( random_generator ) )
         value = round( cont_solution[variable] + 0.5 );
      else
         value = round( cont_solution[variable] - 0.5 );

      return { variable, value };
   }
};
#include "RoundingStrategy.hpp"
#endif // define FIX_RANDOM_ROUNDING_STRATEGY_HPP
