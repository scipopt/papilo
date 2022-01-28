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
      Vec<Fixing<REAL>> fixings{};
      fixings.reserve( cont_solution.size() );
      int infeasible_var = -1;

      // fix first variables with integer value in the solution
      for( int i = 0; i < cont_solution.size(); i++ )
      {
         if( num.isEq( cont_solution[i], num.round( cont_solution[i] ) ) )
         {
            Fixing<REAL> fixing = { i, cont_solution[i] };
            probing_view.setProbingColumn( fixing.get_column_index(),
                                           fixing.get_value() == 1 );
            msg.info( "Fix integer var {} to {}\n", fixing.get_column_index(),
                      fixing.get_value() );
            bool infeasibility_detected = perform_probing_step();
            if( infeasibility_detected )
            {
               probing_view.reset();
               for( auto& f : fixings )
                  probing_view.setProbingColumn( f.get_column_index(),
                                                 f.get_value() );
               fixing = { i, determine_value_for_failed_integer_solution(
                                 cont_solution[i] ) };
               msg.info( "Fix integer var {} to complementary {}\n", fixing.get_column_index(),
                         fixing.get_value() );
               probing_view.setProbingColumn( i, fixing.get_value() );
               bool infeasibility_detected_again = perform_probing_step();
               if( infeasibility_detected_again )
                  return true;
               fixings.push_back( fixing );
            }
            else
            {
               fixings.push_back( fixing );
            }
         }
      }

      while( true )
      {
         probing_view.reset();
         // TODO: maybe there is an more efficient implementation
         if( !fixings.empty() )
         {
            for( auto& fixing : fixings )
               probing_view.setProbingColumn( fixing.get_column_index(),
                                              fixing.get_value() );
         }

         propagate_to_leaf_or_infeasibility( cont_solution, strategy );

         if( probing_view.isInfeasible() )
         {
            // TODO: store fixings since they code the conflict

            fixings = probing_view.get_fixings();
            assert( !fixings.empty() );
            Fixing<REAL> infeasible_fixing = fixings[fixings.size() - 1];
            if( infeasible_var == infeasible_fixing.get_column_index() )
               return false;
            const Fixing<REAL>& fixing{ infeasible_fixing.get_column_index(),
                                        modify_value_due_to_backtrack(
                                            infeasible_fixing.get_value() ) };

            infeasible_var = fixing.get_column_index();
            fixings[fixings.size() - 1] = fixing;
         }
         else
         {
            // TODO: store objective value and solution
            create_solution( result );
            return true;
         }
      }
      return false;
   }

   bool
   perform_probing_step()
   {
      if( probing_view.isInfeasible() )
      {
         msg.info( "changing bound of variable is infeasible row: {} col {} \n",
                   probing_view.get_row_causing_infeasibility(),
                   probing_view.get_col_causing_infeasibility() );
         return true;
      }
      probing_view.propagateDomains();
      probing_view.storeImplications();
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
   determine_value_for_failed_integer_solution( REAL value )
   {
      if( num.isZero( value ) )
         return REAL{ 1 };
      return REAL{ 0 };
   }

 private:
   double
   modify_value_due_to_backtrack( double value )
   {
      return value == 1 ? 0 : 1;
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

   void
   propagate_to_leaf_or_infeasibility( Vec<REAL> cont_solution,
                                       RoundingStrategy<REAL>& strategy )
   {
      while( true )
      {
         Fixing<REAL> fixing =
             strategy.select_diving_variable( cont_solution, probing_view );
         // dive until all vars are fixed (and returned fixing is invalid)
         if( fixing.is_invalid() )
            return;

         msg.info( "Fix var {} to {}\n", fixing.get_column_index(),
                   fixing.get_value() );

         probing_view.setProbingColumn( fixing.get_column_index(),
                                        fixing.get_value() == 1 );
         bool infeasibility_detected = perform_probing_step();
         if( infeasibility_detected )
            return;
      }
   }
};

#endif
