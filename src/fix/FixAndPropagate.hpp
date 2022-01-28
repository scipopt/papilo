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
      int infeasible_var = -1;
      while( true )
      {
         probing_view.reset();
         // TODO: this means recalculating all propagations. Makes that sense or
         // is reverting the propagation better?
         if( !fixings.empty() )
         {
            for( auto& fixing : fixings )
               probing_view.setProbingColumn( fixing.get_column_index(),
                                              fixing.get_value() );
            // TODO: it may be better necessary to update the activities
            // immediately
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
         Fixing<REAL> fixing = strategy.select_diving_variable( cont_solution, probing_view );
         // dive until all vars are fixed (and returned fixing is invalid)
         if( fixing.is_invalid() )
            return;

         msg.info( "Fix var {} to {}\n", fixing.get_column_index(),
                   fixing.get_value() );

         probing_view.setProbingColumn( fixing.get_column_index(),
                                        fixing.get_value() == 1 );
         if( probing_view.isInfeasible() )
         {
            msg.info(
                "changing bound of variable is infeasible row: {} col {} \n",
                probing_view.get_row_causing_infeasibility(),
                probing_view.get_col_causing_infeasibility() );
            return;
         }
         probing_view.propagateDomains();
         probing_view.storeImplications();
         if( probing_view.isInfeasible() )
         {
            msg.info( "propagation is infeasible row: {} col {} \n",
                      probing_view.get_row_causing_infeasibility(),
                      probing_view.get_col_causing_infeasibility() );
            return;
         }
      }
   }
};

#endif
