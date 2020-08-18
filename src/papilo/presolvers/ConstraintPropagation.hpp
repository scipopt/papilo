/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*               This file is part of the program and library                */
/*    PaPILO --- Parallel Presolve for Integer and Linear Optimization       */
/*                                                                           */
/* Copyright (C) 2020  Konrad-Zuse-Zentrum                                   */
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

#ifndef _PAPILO_PRESOLVERS_CONSTRAINT_PROPAGATION_HPP_
#define _PAPILO_PRESOLVERS_CONSTRAINT_PROPAGATION_HPP_

#include "papilo/core/PresolveMethod.hpp"
#include "papilo/core/Problem.hpp"
#include "papilo/core/ProblemUpdate.hpp"
#include "papilo/core/SingleRow.hpp"

namespace papilo
{

template <typename REAL>
class ConstraintPropagation : public PresolveMethod<REAL>
{
 public:
   ConstraintPropagation() : PresolveMethod<REAL>()
   {
      this->setName( "propagation" );
      this->setTiming( PresolverTiming::kFast );
   }

   /// todo how to communicate about postsolve information
   virtual PresolveStatus
   execute( const Problem<REAL>& problem,
            const ProblemUpdate<REAL>& problemUpdate, const Num<REAL>& num,
            Reductions<REAL>& reductions ) override;
};

#ifdef PAPILO_USE_EXTERN_TEMPLATES
extern template class ConstraintPropagation<double>;
extern template class ConstraintPropagation<Quad>;
extern template class ConstraintPropagation<Rational>;
#endif

template <typename REAL>
PresolveStatus
ConstraintPropagation<REAL>::execute( const Problem<REAL>& problem,
                                      const ProblemUpdate<REAL>& problemUpdate,
                                      const Num<REAL>& num,
                                      Reductions<REAL>& reductions )
{
   const auto& domains = problem.getVariableDomains();
   const auto& activities = problem.getRowActivities();
   const auto& changedactivities = problemUpdate.getChangedActivities();
   const auto& consMatrix = problem.getConstraintMatrix();
   const auto& lhsValues = consMatrix.getLeftHandSides();
   const auto& rhsValues = consMatrix.getRightHandSides();
   const auto& rflags = consMatrix.getRowFlags();

   PresolveStatus result = PresolveStatus::kUnchanged;

   // for LP constraint propagation we might want to weaken the bounds by some
   // small amount above the feasibility tolerance
   REAL weakenbounds =
       problem.getNumIntegralCols() == 0
           ? REAL{ problemUpdate.getPresolveOptions().weakenlpvarbounds *
                   num.getFeasTol() }
           : REAL{ 0 };

   auto add_boundchange = [&]( BoundChange boundChange, int col, REAL val ) {
      // do not accept huge values as bounds
      if( num.isHugeVal( val ) )
         return;

      if( boundChange == BoundChange::kLower )
      {
         assert( domains.flags[col].test( ColFlag::kLbInf ) ||
                 val > domains.lower_bounds[col] );

         if( domains.flags[col].test( ColFlag::kIntegral, ColFlag::kImplInt ) )
            val = num.feasCeil( val );

         if( !domains.flags[col].test( ColFlag::kUbInf ) )
         {
            // compute distance of new lower bound to the upper bound
            REAL bnddist = domains.upper_bounds[col] - val;

            // bound exceeded by more then feastol means infeasible
            if( bnddist < -num.getFeasTol() )
            {
               result = PresolveStatus::kInfeasible;
               return;
            }

            // if the upper bound is reached, or reached within tolerances and
            // the change of feasibility is also within tolerances fix to the
            // upper bound
            if( bnddist <= 0 || ( bnddist <= num.getFeasTol() &&
                                  consMatrix.getMaxFeasChange( col, bnddist ) <=
                                      num.getFeasTol() ) )
            {
               // todo reductions.forcingRowToUpper(currentrow, col);
               reductions.fixCol( col, domains.upper_bounds[col] );
               result = PresolveStatus::kReduced;
               return;
            }
         }

         val -= weakenbounds;
         if( domains.flags[col].test( ColFlag::kLbInf ) ||
             val - domains.lower_bounds[col] > +1000 * num.getFeasTol() )
         {
            reductions.changeColLB( col, val );
            result = PresolveStatus::kReduced;
         }
      }
      else
      {
         assert( boundChange == BoundChange::kUpper );
         assert( domains.flags[col].test( ColFlag::kUbInf ) ||
                 val < domains.upper_bounds[col] );

         if( domains.flags[col].test( ColFlag::kIntegral, ColFlag::kImplInt ) )
            val = num.feasFloor( val );

         if( !domains.flags[col].test( ColFlag::kLbInf ) )
         {
            // compute distance of new upper bound to the lower bound
            REAL bnddist = val - domains.lower_bounds[col];

            // bound exceeded by more then feastol means infeasible
            if( bnddist < -num.getFeasTol() )
            {
               result = PresolveStatus::kInfeasible;
               return;
            }

            // if the lower bound is reached, or reached within tolerances and
            // the change of feasibility is also within tolerances fix to the
            // lower bound
            if( bnddist <= 0 || ( bnddist <= num.getFeasTol() &&
                                  consMatrix.getMaxFeasChange( col, bnddist ) <=
                                      num.getFeasTol() ) )
            {
               // todo reductions.forcingRowToLower(currentrow, col);
               reductions.fixCol( col, domains.lower_bounds[col] );
               result = PresolveStatus::kReduced;
               return;
            }
         }

         val += weakenbounds;
         if( domains.flags[col].test( ColFlag::kUbInf ) ||
             val - domains.upper_bounds[col] < -1000 * num.getFeasTol() )
         {
            reductions.changeColUB( col, val );
            result = PresolveStatus::kReduced;
         }
      }
   };

   for( int row : changedactivities )
   {
      auto rowvec = consMatrix.getRowCoefficients( row );

      assert( !consMatrix.isRowRedundant( row ) );

      switch( rowvec.getLength() )
      {
      case 0:
         if( ( !rflags[row].test( RowFlag::kLhsInf ) &&
               num.isFeasGT( lhsValues[row], 0 ) ) ||
             ( !rflags[row].test( RowFlag::kRhsInf ) &&
               num.isFeasLT( rhsValues[row], 0 ) ) )
            result = PresolveStatus::kInfeasible;
         else
            reductions.markRowRedundant( row );
         break;
      case 1:
         // do nothing, singleton row presolver handles this bound change
         break;
      default:
         propagate_row( rowvec.getValues(), rowvec.getIndices(),
                        rowvec.getLength(), activities[row], lhsValues[row],
                        rhsValues[row], rflags[row], domains.lower_bounds,
                        domains.upper_bounds, domains.flags, add_boundchange );
      }

      if( result == PresolveStatus::kInfeasible )
         break;
   }

   return result;
}

} // namespace papilo

#endif
