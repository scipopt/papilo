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
   Vec<bool> no_down_locks;
   Vec<bool> no_up_locks;
   REAL norm;


 public:
   FractionalRoundingStrategy( Num<REAL> num_, Problem<REAL> problem_ )
       : num( num_ ), no_down_locks( {} ), no_up_locks()
   {
      auto obj = problem_.getObjective().coefficients;
      norm = *max_element( std::begin( obj ), std::end( obj ) );
      auto rflags = problem_.getRowFlags();
      const ConstraintMatrix<REAL>& matrix = problem_.getConstraintMatrix();
      int n_cols = problem_.getNCols();
      no_up_locks.resize( n_cols );
      no_down_locks.resize( n_cols );

#ifdef PAPILO_TBB
      tbb::parallel_for(
          tbb::blocked_range<int>( 0, n_cols ),
          [this]( const tbb::blocked_range<int>& c ) {
             for( int col = c.begin(); col != c.end(); ++col )
#else
      for( int col = 0; col < n_cols; ++col )
#endif
             {
                int n_up_locks = 0;
                int n_down_locks = 0;

                auto colvec = matrix.getColumnCoefficients( col );
                const REAL* values = colvec.getValues();
                const int* rowinds = colvec.getIndices();

                bool skip = false;
                for( int j = 0; j < colvec.getLength(); ++j )
                {
                   count_locks( values[j], rflags[rowinds[j]], n_down_locks,
                                n_up_locks );
                   if( n_up_locks != 0 && n_down_locks != 0 )
                   {
                      no_down_locks[col] = false;
                      no_up_locks[col] = false;
                      skip = true;
                      break;
                   }
                }
                if(skip)
                   continue;
                assert( n_down_locks != 0 || n_up_locks != 0 );
                no_down_locks[col] = ( n_down_locks == 0 );
                no_up_locks[col] = ( n_up_locks == 0 );
             }
#ifdef PAPILO_TBB
          } );
#endif
      assert( no_down_locks.size() == no_up_locks.size() );
      assert( no_down_locks.size() == problem_.getNCols() );
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
      REAL score = -1;
      auto obj = view.get_obj();


      for( int i = 0; i < cont_solution.size(); i++ )
      {
         if( num.isIntegral( cont_solution[i] ) ||
             num.isEq( view.getProbingUpperBounds()[i],
                       view.getProbingLowerBounds()[i] ) ||
             !view.is_within_bounds( i, cont_solution[i] ) ||
             !view.is_integer_variable( i ) )
            continue;

         // may round down if now down locks
         bool may_round_down = no_up_locks[i];
         bool may_round_up = no_down_locks[i];
         REAL frac = cont_solution[i] - num.epsFloor( cont_solution[i] );

         /* divide by objective norm to normalize obj into [-1,1] */
         REAL new_val;
         REAL gain;
         assert(!(may_round_up && may_round_down));
         if( !may_round_down && may_round_up )
         {
            new_val = num.epsCeil( cont_solution[i] );
            gain = obj[i] / norm * ( 1 - frac );
            frac = 1 - frac;
         }
         else if( may_round_down && !may_round_up )
         {
            new_val = num.epsFloor( cont_solution[i] );
            gain = -obj[i] / norm * frac;
         }
         else if( frac > 0.5 )
         {
            assert( !may_round_up );
            assert( !may_round_down );
            new_val = num.epsCeil( cont_solution[i] );
            gain = obj[i] / norm * ( 1 - frac );
            frac = 1 - frac;
         }
         else
         {
            assert( !may_round_up );
            assert( !may_round_down );
            new_val = num.epsFloor( cont_solution[i] );
            gain = -obj[i] / norm * frac;
         }

         /* penalize too small fractions */
         if( frac < 0.01 )
            frac += 10.0;

         assert(view.is_integer_variable(i));
         /* prefer decisions on binary variables */
         if( view.getProbingUpperBounds()[i] != 1 || view.getProbingLowerBounds()[i] != 0 )
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
      }
      return { variable, value };
   }
};

#endif
