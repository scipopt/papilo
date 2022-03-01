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

#ifndef FIX_MOST_FRACTIONAL_ROUNDING_STRATEGY_HPP
#define FIX_MOST_FRACTIONAL_ROUNDING_STRATEGY_HPP

#include "fix/strategy/RoundingStrategy.hpp"

template <typename REAL>
class MostFractionalRoundingStrategy : public RoundingStrategy<REAL>
{

   const Num<REAL> num;

 public:
   MostFractionalRoundingStrategy( Num<REAL> num_ )
       : num( num_ )
   {
   }

   void
   update_data_structure_before_dive() override
   {
   }

   Fixing<REAL>
   select_rounding_variable( const Vec<REAL>& cont_solution,
                             const ProbingView<REAL>& view ) override
   {
      REAL value = -1;
      int variable = -1;
      REAL max_frac = -1;

      for( int i = 0; i < cont_solution.size(); i++ )
      {
         if( num.isIntegral( cont_solution[i] ) ||
             num.isEq( view.getProbingUpperBounds()[i],
                       view.getProbingLowerBounds()[i] ) ||
             !view.is_within_bounds( i, cont_solution[i] ) ||
             !view.is_integer_variable( i ) )
            continue;

         REAL frac = REAL{ 0.5 } -
                     abs( cont_solution[i] - num.epsFloor( cont_solution[i] ) -
                          REAL{ 0.5 } );
         if( num.isGT( frac, max_frac ) )
         {
            max_frac = frac;
            variable = i;
            if( num.isGT( cont_solution[i] - num.epsFloor( cont_solution[i] ),
                          REAL{ 0.5 } ) )
               value = num.epsCeil( cont_solution[i] );
            else
               value = num.epsFloor( cont_solution[i] );
         }
      }
      return { variable, value };
   }
};

#endif
