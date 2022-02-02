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
   bool scale_score;

   typedef std::mt19937 MyRNG;
   uint32_t seed;

   MyRNG random_generator;

 public:
   FarkasRoundingStrategy( uint32_t seed_, Num<REAL> num_, bool scale_score_ )
       : seed( seed_ ), num( num_ ), scale_score( scale_score_ )
   {
      random_generator.seed( seed );
   }

   Fixing<REAL>
   select_rounding_variable( const Vec<REAL>& cont_solution,
                             const ProbingView<REAL>& view ) override
   {
      REAL value = -1;
      int variable = -1;
      REAL score = -1;
      bool stored_var_binary = false;
      std::uniform_real_distribution<double> small_number_generator( 0.0000001,
                                                                     0.000001 );
      auto obj = view.get_obj();

      for( int i = 0; i < cont_solution.size(); i++ )
      {
         if( num.isIntegral( cont_solution[i] ) ||
             num.isEq( view.getProbingUpperBounds()[i],
                       view.getProbingLowerBounds()[i] ) ||
             !view.is_within_bounds( i, cont_solution[i] ) ||
             !view.is_integer_variable( i ) )
            continue;

         /* prefer decisions on binary variables */
         bool current_var_binary = view.getProbingUpperBounds()[i] == 1 &&
                                   view.getProbingLowerBounds()[i] == 0 &&
                                   !num.isZero( obj[i] );
         if( !current_var_binary && stored_var_binary )
            continue;

         REAL current_score =
             abs( obj[i] ) + small_number_generator( random_generator );

         REAL frac = cont_solution[i] - num.epsFloor( cont_solution[i] );
         bool round_up =
             num.isLT( obj[i], 0 ) || ( num.isZero( obj[i] ) && frac > 0.5 );
         if( scale_score )
         {
            if( round_up )
               current_score = current_score * ( 1 - frac );
            else
               current_score = current_score * frac;
         }
         if( ( current_var_binary && !stored_var_binary ) || variable == -1 ||
             current_score > score )
         {
            assert( ( current_var_binary && !stored_var_binary ) ||
                    current_var_binary == stored_var_binary );
            score = current_score;
            variable = i;
            stored_var_binary = current_var_binary;
            if( round_up )
               value = num.epsCeil( cont_solution[i] );
            else
               value = num.epsFloor( cont_solution[i] );
         }
      }
      return { variable, value };
   }
};

#endif
