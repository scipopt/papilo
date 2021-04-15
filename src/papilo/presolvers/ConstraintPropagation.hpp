/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*               This file is part of the program and library                */
/*    PaPILO --- Parallel Presolve for Integer and Linear Optimization       */
/*                                                                           */
/* Copyright (C) 2020-2021 Konrad-Zuse-Zentrum                               */
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

   PresolveStatus
   execute( const Problem<REAL>& problem,
            const ProblemUpdate<REAL>& problemUpdate, const Num<REAL>& num,
            Reductions<REAL>& reductions ) override;

 private:
   PresolveStatus
   perform_propagation_step( const Num<REAL>& num, Reductions<REAL>& reductions,
                             const VariableDomains<REAL>& domains,
                             const Vec<RowActivity<REAL>>& activities,
                             const ConstraintMatrix<REAL>& consMatrix,
                             const Vec<REAL>& lhsValues,
                             const Vec<REAL>& rhsValues,
                             const Vec<RowFlags>& rflags, REAL weakenbounds,
                             int row ) const;
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

   if( problemUpdate.getPresolveOptions().runs_sequentiell() )
   {
      for( int row : changedactivities )
      {
         PresolveStatus local_result = perform_propagation_step(
             num, reductions, domains, activities, consMatrix, lhsValues,
             rhsValues, rflags, weakenbounds, row );
         if( local_result == PresolveStatus::kInfeasible )
            return PresolveStatus::kInfeasible;
         if( local_result == PresolveStatus::kReduced )
            result = PresolveStatus::kReduced;
      }
   }
   else
   {
      bool infeasible = false;
      Vec<Reductions<REAL>> stored_reductions( changedactivities.size() );
      tbb::parallel_for(
          tbb::blocked_range<int>( 0, changedactivities.size() ),
          [&]( const tbb::blocked_range<int>& r ) {
             for( int j = r.begin(); j < r.end(); ++j )
             {
                PresolveStatus local_result = perform_propagation_step(
                    num, stored_reductions[j], domains, activities, consMatrix, lhsValues,
                    rhsValues, rflags, weakenbounds, changedactivities[j] );
                assert( local_result == PresolveStatus::kInfeasible ||
                        local_result == PresolveStatus::kReduced ||
                        local_result == PresolveStatus::kUnchanged );
                if( local_result == PresolveStatus::kInfeasible )
                   infeasible = true;
                else if( local_result == PresolveStatus::kReduced )
                   result = PresolveStatus::kReduced;
             }
          } );
      if( infeasible )
         return PresolveStatus::kInfeasible;
      else if( result == PresolveStatus::kUnchanged )
         return PresolveStatus::kUnchanged;

      for( int i = 0; i < stored_reductions.size(); ++i )
      {
         Reductions<REAL> reds = stored_reductions[i];
         if( reds.size() > 0 )
         {
            for( const auto& reduction : reds.getReductions() )
            {
               reductions.add_reduction( reduction.row, reduction.col,
                                         reduction.newval );
            }
         }
      }
   }

   return result;
}
template <typename REAL>
PresolveStatus
ConstraintPropagation<REAL>::perform_propagation_step(
    const Num<REAL>& num, Reductions<REAL>& reductions,
    const VariableDomains<REAL>& domains,
    const Vec<RowActivity<REAL>>& activities,
    const ConstraintMatrix<REAL>& consMatrix, const Vec<REAL>& lhsValues,
    const Vec<REAL>& rhsValues, const Vec<RowFlags>& rflags, REAL weakenbounds,
    int row ) const
{
   PresolveStatus result = PresolveStatus::kUnchanged;
   auto rowvec = consMatrix.getRowCoefficients( row );

   assert( !consMatrix.isRowRedundant( row ) );

   switch( rowvec.getLength() )
   {
   case 0:
      if( ( !rflags[row].test( RowFlag::kLhsInf ) &&
            num.isFeasGT( lhsValues[row], 0 ) ) ||
          ( !rflags[row].test( RowFlag::kRhsInf ) &&
            num.isFeasLT( rhsValues[row], 0 ) ) )
         return PresolveStatus::kInfeasible;
      else
         reductions.markRowRedundant( row );
      break;
   case 1:
      // do nothing, singleton row presolver handles this bound change
      break;
   default:
      auto add_boundchange = [&]( BoundChange boundChange, int col, REAL val ) {
         // do not accept huge values as bounds
         if( num.isHugeVal( val ) )
            return;

         if( boundChange == BoundChange::kLower )
         {
            assert( domains.flags[col].test( ColFlag::kLbInf ) ||
                    val > domains.lower_bounds[col] );

            if( domains.flags[col].test( ColFlag::kIntegral,
                                         ColFlag::kImplInt ) )
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

               // if the upper bound is reached, or reached within tolerances
               // and the change of feasibility is also within tolerances fix to
               // the upper bound
               if( bnddist <= 0 || ( bnddist <= num.getFeasTol() &&
                                     consMatrix.getMaxFeasChange(
                                         col, bnddist ) <= num.getFeasTol() ) )
               {
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

            if( domains.flags[col].test( ColFlag::kIntegral,
                                         ColFlag::kImplInt ) )
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

               // if the lower bound is reached, or reached within tolerances
               // and the change of feasibility is also within tolerances fix to
               // the lower bound
               if( bnddist <= 0 || ( bnddist <= num.getFeasTol() &&
                                     consMatrix.getMaxFeasChange(
                                         col, bnddist ) <= num.getFeasTol() ) )
               {
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
      propagate_row( rowvec.getValues(), rowvec.getIndices(),
                     rowvec.getLength(), activities[row], lhsValues[row],
                     rhsValues[row], rflags[row], domains.lower_bounds,
                     domains.upper_bounds, domains.flags, add_boundchange );
   }
   return result;
}

} // namespace papilo

#endif
