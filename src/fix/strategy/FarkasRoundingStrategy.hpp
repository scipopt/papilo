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

#ifndef FIX_FARKAS_ROUNDING_STRATEGY_HPP
#define FIX_FARKAS_ROUNDING_STRATEGY_HPP

#include "papilo/core/ProbingView.hpp"

#include <cassert>
#include <random>

using namespace papilo;

template <typename REAL>
class FarkasRoundingStrategy : public RoundingStrategy<REAL>
{

   const Num<REAL> num;

   typedef std::mt19937 MyRNG;
   uint32_t seed;

   MyRNG random_generator;

 public:
   FarkasRoundingStrategy( uint32_t seed_, Num<REAL> num_ )
       : seed( seed_ ), num( num_ )
   {
      random_generator.seed( seed );
   }

   Fixing<REAL>
   select_rounding_variable( const Vec<REAL>& cont_solution,
                             const ProbingView<REAL>& view ) override
   {
      // this is currently fractional diving
      REAL value = -1;
      int variable = -1;
      REAL score = -1;
      std::uniform_int_distribution<int> dist_rounding( 0, 1e5 );
      auto obj = view.get_obj();

      for( int i = 0; i < cont_solution.size(); i++ )
      {
         if( num.isIntegral( cont_solution[i] ) ||
             num.isEq( view.getProbingUpperBounds()[i],
                       view.getProbingLowerBounds()[i] ) ||
             !view.is_within_bounds( i, cont_solution[i] ) ||
             !view.is_integer_variable( i ) )
            continue;

         REAL current_score = abs( obj[i] ) + dist_rounding( random_generator );

         /* prefer decisions on binary variables */
         if( view.getProbingUpperBounds()[i] != 1 )
            current_score = -1 / current_score;

         if( current_score > score || variable == -1 )
         {
            if( num.isLT( obj[i], 0 ) )
            {
               variable = i;
               value = num.epsCeil( cont_solution[i] );
            }
            else if( num.isGT( obj[i], 0 ) )
            {
               variable = i;
               value = num.epsFloor( cont_solution[i] );
            }
            else
            {
               REAL frac = cont_solution[i] - num.epsFloor( cont_solution[i] );
               REAL proposed_value;
               if( frac > 0.5 )
                  proposed_value = num.epsCeil( cont_solution[i] );
               else
                  proposed_value = num.epsFloor( cont_solution[i] );
               variable = i;
               value = proposed_value;
            }
         }
      }
      return { variable, value };
   }
};

#endif
