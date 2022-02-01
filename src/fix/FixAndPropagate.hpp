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

/***
 * This class performs a fix-and-propagate algorithm:
 *
 * V = all integer variables with non integer solution and proposed value is within bounds
 * propagate domains does not propagate violated rows
 *
 * while V is not empty
 *  max, var_max, val_max = max_{ v in V} score  (defined by strategy)
 *  fix var_max to value val_max
 *  propagate domains
 *  if perform_backtrack:
 *      [if propagation or fixing is infeasible to backtrack by fixing var_max to val_max +/-1]
 *      [if this is still infeasible then perform no more backtracks]
 *
 * for all non fixed variables v
 *  if lb_v < sol(v) < ub_v
 *      fix v to sol(v)
 *  else lb_v > sol(v)
 *      fix v to lb_v
 *  else
 *      fix v to ub_v
 *  propagate domains*
 *

 * @tparam REAL the arithmetic parameter
 */
template <typename REAL>
class FixAndPropagate
{
   bool perform_backtracking = true;
   Message msg;
   Num<REAL> num;
   ProbingView<REAL> probing_view;

 public:
   FixAndPropagate( Message msg_, Num<REAL> num_, Problem<REAL>& problem_,
                    ProbingView<REAL> view_, bool perform_backtracking_ )
       : msg( msg_ ), num( num_ ), probing_view( view_ ),
         perform_backtracking( perform_backtracking_ )
   {
   }

   bool
   fix_and_propagate( const Vec<REAL>& cont_solution, Vec<REAL>& result,
                      RoundingStrategy<REAL>& strategy )
   {
      // if no backtrack just "dive" to the node whether it is infeasible or not
      if( !perform_backtracking )
      {
         propagate_to_leaf_or_infeasibility( cont_solution, strategy, false );
         fix_remaining_integer_solutions( cont_solution );
         create_solution( result );
         return probing_view.isInfeasible();
      }

      while( true )
      {
         propagate_to_leaf_or_infeasibility( cont_solution, strategy, true );

         if( probing_view.isInfeasible() )
         {
            if( perform_backtracking )
            {
               Vec<Fixing<REAL>> fixings = probing_view.get_fixings();
               assert( !fixings.empty() );
               Fixing<REAL> last_fix = fixings[fixings.size() - 1];

               probing_view.reset();
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
               {
                  propagate_to_leaf_or_infeasibility( cont_solution, strategy,
                                                      false );
                  fix_remaining_integer_solutions( cont_solution );
                  create_solution( result );
                  return probing_view.isInfeasible();
               }
            }
         }
         else
         {
            fix_remaining_integer_solutions( cont_solution );
            create_solution( result );
            return probing_view.isInfeasible();
         }
      }
   }

 private:
   void
   propagate_to_leaf_or_infeasibility( Vec<REAL> cont_solution,
                                       RoundingStrategy<REAL>& strategy,
                                       bool stop_at_infeasibility )
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

         probing_view.setProbingColumn( fixing.get_column_index(),
                                        fixing.get_value() );
         bool infeasibility_detected = perform_probing_step();
         if( stop_at_infeasibility && infeasibility_detected )
            return;
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
         bool ge_lb = num.isGE( cont_solution[i], lowerBounds[i] );
         bool le_ub = num.isLE( cont_solution[i], upperBounds[i] );
         if( !num.isEq( upperBounds[i], lowerBounds[i] ) )
         {
            REAL value;
            if( !probing_view.is_integer_variable( i ) ) {
               if( ge_lb && le_ub )
                  value = cont_solution[i];
               else if( ge_lb )
                  value = upperBounds[i];
               else
               {
                  assert( le_ub );
                  value = lowerBounds[i];
               }
            }
            else
            {
               // only variable with integer value in the solution should be
               // non fixed
               //TODO: this assertion is not correct
               assert( num.isEq( cont_solution[i],
                                 num.round( cont_solution[i] ) ) );
               if( ge_lb && le_ub )
                  value = cont_solution[i];
               else if( ge_lb )
                  value = upperBounds[i];
               else
               {
                  assert( le_ub );
                  value = lowerBounds[i];
               }
            }
            probing_view.setProbingColumn( i, value );
            msg.info( "Fix integer var {} to {}\n", i, value );

            perform_probing_step();
         }
   }
   return probing_view.isInfeasible();
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
}
;

#endif
