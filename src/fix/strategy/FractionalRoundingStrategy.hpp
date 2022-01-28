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


#ifndef FIX_FRACTIONAL_ROUNDING_STRATEGY_HPP
#define FIX_FRACTIONAL_ROUNDING_STRATEGY_HPP

#include "fix/strategy/RoundingStrategy.hpp"


template <typename REAL>
class FractionalRoundingStrategy : public RoundingStrategy<REAL>
{

   const Num<REAL> num;

 public:
   FractionalRoundingStrategy( Num<REAL> num_ ) : num( num_ ) {}

   Fixing<REAL>
   select_diving_variable( const Vec<REAL>& cont_solution,
                           const ProbingView<REAL>& view ) override
   {

      // this is currently fractional diving
      REAL value = -1;
      int variable = -1;
      REAL score = -1;

      for( int i = 0; i < cont_solution.size(); i++ )
      {
         REAL frac = cont_solution[i] - floor( cont_solution[i] );
         if( frac == 0 || num.isEq( view.getProbingUpperBounds()[i],
                                    view.getProbingLowerBounds()[i] ) )
            continue;
         else if( frac > 0.5 )
         {
            if( variable == -1 || 1 - frac > score )
            {
               score = 1 - frac;
               variable = i;
               value = ceil( cont_solution[i] );
            }
         }
         else
         {
            if( variable == -1 || frac > score )
            {
               score = frac;
               variable = i;
               value = floor( cont_solution[i] );
            }
         }
      }
      return { variable, value };
   }
};

#endif // PAPILO_FRACTIONALROUNDINGSTRATEGY_HPP
