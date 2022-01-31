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

#ifndef FIX_FIX_AND_PROPAGATE_HPP
#define FIX_FIX_AND_PROPAGATE_HPP

#include "papilo/core/Objective.hpp"
#include "papilo/core/Presolve.hpp"
#include "papilo/core/ProbingView.hpp"

#include "fix/strategy/RoundingStrategy.hpp"
#include "papilo/io/MpsParser.hpp"
#include <cassert>
#include <fstream>
#include <string>

using namespace papilo;

// TODO: Probing does only work with boolean values

template <typename REAL>
class FixAndPropagate
{
   Message msg;
   Num<REAL> num;
   ProbingView<REAL> probing_view;

 public:
   FixAndPropagate( Message msg_, Num<REAL> num_, Problem<REAL>& problem_,
                    ProbingView<REAL> view_ )
       : msg( msg_ ), num( num_ ), probing_view( view_ )
   {
   }

   bool
   fix_and_propagate( const Vec<REAL>& cont_solution, Vec<REAL>& result,
                      RoundingStrategy<REAL>& strategy )
   {
      while( true )
      {
         propagate_to_leaf_or_infeasibility( cont_solution, strategy );

         if( probing_view.isInfeasible() )
         {
            // TODO: store fixings since they code the conflict

            Vec<Fixing<REAL>> fixings = probing_view.get_fixings();
            assert( !fixings.empty() );
            Fixing<REAL> last_fix = fixings[fixings.size() - 1];

            probing_view.reset();
            // TODO: maybe there is an more efficient implementation
            for( int i = 0; i < fixings.size() - 1; i++ )
               probing_view.setProbingColumn( fixings[i].get_column_index(),
                                              fixings[i].get_value() );
            probing_view.setProbingColumn(
                last_fix.get_column_index(),
                modify_value_due_to_backtrack(
                    last_fix.get_value(),
                    cont_solution[last_fix.get_column_index()] ) );
            bool infeasible = perform_probing_step();
            if( infeasible )
               return false;
         }
         else
         {
            if( !fix_remaining_integer_solutions( cont_solution ) )
               return false;
            // TODO: store objective value and solution
            create_solution( result );
            return true;
         }
      }
      return false;
   }

 private:
   void
   propagate_to_leaf_or_infeasibility( Vec<REAL> cont_solution,
                                       RoundingStrategy<REAL>& strategy )
   {
      while( true )
      {
         Fixing<REAL> fixing =
             strategy.select_rounding_variable( cont_solution, probing_view );
         // dive until all vars are fixed (and returned fixing is invalid)
         if( fixing.is_invalid() )
            return;

         msg.info( "Fix var {} to {}\n", fixing.get_column_index(),
                   fixing.get_value() );

         //TODO: check if rounded variable is valid

         probing_view.setProbingColumn( fixing.get_column_index(),
                                        fixing.get_value() );
         bool infeasibility_detected = perform_probing_step();
//         if( infeasibility_detected )
//            return;
      }
   }

   bool
   perform_probing_step()
   {
      if( probing_view.isInfeasible() )
      {
         msg.info( "changing bound of variable is infeasible row: {} col {} \n",
                   probing_view.get_row_causing_infeasibility(),
                   probing_view.get_col_causing_infeasibility() );
//         return true;
      }
      probing_view.propagateDomains();
      //      probing_view.storeImplications();
      if( probing_view.isInfeasible() )
      {
         msg.info( "propagation is infeasible row: {} col {} \n",
                   probing_view.get_row_causing_infeasibility(),
                   probing_view.get_col_causing_infeasibility() );
         return true;
      }
      return false;
   }

   // TODO fix this
   REAL
   modify_value_due_to_backtrack( REAL value, REAL solution_value )
   {
      if( num.isGE( value, solution_value ) )
      {
         assert( num.feasFloor( solution_value ) == value - 1 );
         return value - 1;
      }
      else
      {
         assert( num.feasCeil( solution_value ) == value + 1 );
         return value + 1;
      }
   }

   bool
   fix_remaining_integer_solutions( const Vec<REAL>& cont_solution )
   {
      for( int i = 0; i < cont_solution.size(); i++ )
      {
         auto lowerBounds = probing_view.getProbingLowerBounds();
         auto upperBounds = probing_view.getProbingUpperBounds();
         if( !num.isEq( upperBounds[i], lowerBounds[i] ) )
         {
            // only variable with integer value in the solution should be
            // non fixed
            assert(
                num.isEq( cont_solution[i], num.round( cont_solution[i] ) ) );

            REAL value;
            bool ge_lb = num.isGE( cont_solution[i], lowerBounds[i] );
            bool le_ub = num.isLE( cont_solution[i], upperBounds[i] );
            if( ge_lb && le_ub )
               value = cont_solution[i];
            else if( ge_lb )
               value = upperBounds[i];
            else
            {
               assert( le_ub );
               value = lowerBounds[i];
            }
            probing_view.setProbingColumn( i, value );
            msg.info( "Fix integer var {} to {}\n", i, value );

            bool infeasibility_detected = perform_probing_step();
            if( infeasibility_detected )
               return false;
         }
      }
      return true;
   }

   void
   create_solution( Vec<REAL>& result )
   {
      for( int i = 0; i < probing_view.getProbingUpperBounds().size(); i++ )
      {
         assert( probing_view.getProbingUpperBounds()[i] ==
                 probing_view.getProbingLowerBounds()[i] );
         result[i] = probing_view.getProbingUpperBounds()[i];
      }
   }
};

#endif
