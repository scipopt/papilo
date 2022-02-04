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
   select_rounding_variable( const Vec<REAL>& cont_solution,
                             const ProbingView<REAL>& view ) override
   {
      REAL value = -1;
      int variable = -1;
      REAL score = -1;
      auto obj = view.get_obj();

      REAL norm = *max_element( std::begin( obj ), std::end( obj ) );

      for( int i = 0; i < cont_solution.size(); i++ )
      {
         if( num.isIntegral( cont_solution[i] ) ||
             num.isEq( view.getProbingUpperBounds()[i],
                       view.getProbingLowerBounds()[i] ) ||
             !view.is_within_bounds( i, cont_solution[i] ) ||
             !view.is_integer_variable( i ) )
            continue;

         // TODO: implement it more efficient and the norm
         std::pair<bool, bool> has_locks = view.has_locks( i );
         // may round down if now down locks
         bool may_round_down = has_locks.second;
         bool may_round_up = has_locks.first;
         REAL frac = cont_solution[i] - num.epsFloor( cont_solution[i] );

         /* divide by objective norm to normalize obj into [-1,1] */
         REAL new_val;
         REAL gain;
         if( may_round_down && !may_round_up )
         {
            new_val = num.epsCeil( cont_solution[i] );
            gain = obj[i] / norm * ( 1 - frac );
            frac = 1 - frac;
         }
         else if( !may_round_down && may_round_up )
         {
            new_val = num.epsFloor( cont_solution[i] );
            gain = -obj[i] / norm * frac;
         }
         else if( frac > 0.5 )
         {
            assert( may_round_up == may_round_down );
            new_val = num.epsCeil( cont_solution[i] );
            gain = obj[i] / norm * ( 1 - frac );
            frac = 1 - frac;
         }
         else
         {
            assert( may_round_up == may_round_down );
            new_val = num.epsFloor( cont_solution[i] );
            gain = -obj[i] / norm * frac;
         }

         /* penalize too small fractions */
         if( frac < 0.01 )
            frac += 10.0;

         /* prefer decisions on binary variables */
         if( view.getProbingUpperBounds()[i] != 1 )
            frac = frac * 1000;

         /* prefer variables which cannot be rounded by scoring their
          * fractionality */
         if( !( may_round_down || may_round_up ) )
         {
            if( frac > score || variable == -1 )
            {
               value = new_val;
               variable = i;
               score = frac;
            }
         }
         else
         {
            if( -2.0 - gain > score || variable == -1 )
            {
               value = new_val;
               variable = i;
               score = frac;
            }
         }

         //         assert( !num.isZero( frac ) );
         //         if( frac > 0.5 )
         //         {
         //            REAL current_score = ( 1 - frac ) * obj[i];
         //            REAL prosposed_value = num.epsCeil( cont_solution[i] );
         //            if( ( variable == -1 || current_score > score ) &&
         //                view.is_within_bounds( i, prosposed_value ) )
         //            {
         //               score = current_score;
         //               variable = i;
         //               value = prosposed_value;
         //            }
         //         }
         //         else
         //         {
         //            REAL current_score = frac * obj[i];
         //            REAL prosposed_value = num.epsFloor( cont_solution[i] );
         //            if( ( variable == -1 || current_score > score ) &&
         //                view.is_within_bounds( i, prosposed_value ) )
         //            {
         //               score = current_score;
         //               variable = i;
         //               value = prosposed_value;
         //            }
         //         }
      }
      return { variable, value };
   }
};

#endif
