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

#ifndef _PAPILO_PRESOLVERS_DUAL_FIX_HPP_
#define _PAPILO_PRESOLVERS_DUAL_FIX_HPP_

#include "papilo/core/PresolveMethod.hpp"
#include "papilo/core/Problem.hpp"
#include "papilo/core/ProblemUpdate.hpp"

namespace papilo
{

/// dual-fixing presolve method which looks at the coefficients of the objective
/// and the column entries and performs a dual fixing if possible
/// If fixing is not possible, it tries to strengthen the bounds.
template <typename REAL>
class DualFix : public PresolveMethod<REAL>
{
 public:
   DualFix() : PresolveMethod<REAL>()
   {
      this->setName( "dualfix" );
      this->setTiming( PresolverTiming::kMedium );
   }

   bool
   initialize( const Problem<REAL>& problem,
               const PresolveOptions& presolveOptions ) override
   {
      if( presolveOptions.dualreds == 0 )
         this->setEnabled( false );
      return false;
   }

   virtual PresolveStatus
   execute( const Problem<REAL>& problem,
            const ProblemUpdate<REAL>& problemUpdate, const Num<REAL>& num,
            Reductions<REAL>& reductions ) override;

   std::pair<bool, REAL>
   calc_upper_bound_with_dual_substitution(
       const Num<REAL>& num, const ConstraintMatrix<REAL>& consMatrix,
       const Vec<RowActivity<REAL>>& activities, const Vec<REAL>& lbs,
       const Vec<REAL>& ubs, const Vec<REAL>& rhs, const Vec<REAL>& lhs, int i,
       int collen, const REAL* values, const int* rowinds,
       const Vec<ColFlags>& cflags, const Vec<RowFlags>& rflags );

   std::pair<bool, REAL>
   calc_lower_bound_with_dual_substitution(
       const Num<REAL>& num, const ConstraintMatrix<REAL>& consMatrix,
       const Vec<RowActivity<REAL>>& activities, const Vec<REAL>& lbs,
       const Vec<REAL>& ubs, const Vec<REAL>& lhs, const Vec<REAL>& rhs, int i,
       int collen, const REAL* values, const int* rowinds,
       const Vec<ColFlags>& cflags, const Vec<RowFlags>& rflags );
};

#ifdef PAPILO_USE_EXTERN_TEMPLATES
extern template class DualFix<double>;
extern template class DualFix<Quad>;
extern template class DualFix<Rational>;
#endif

template <typename REAL>
PresolveStatus
DualFix<REAL>::execute( const Problem<REAL>& problem,
                        const ProblemUpdate<REAL>& problemUpdate,
                        const Num<REAL>& num, Reductions<REAL>& reductions )
{
   const auto& consMatrix = problem.getConstraintMatrix();
   const Vec<RowActivity<REAL>>& activities = problem.getRowActivities();
   const Vec<ColFlags>& cflags = problem.getColFlags();
   const Vec<REAL>& objective = problem.getObjective().coefficients;
   const Vec<REAL>& lbs = problem.getLowerBounds();
   const Vec<REAL>& ubs = problem.getUpperBounds();
   const int ncols = consMatrix.getNCols();
   const Vec<RowFlags>& rflags = consMatrix.getRowFlags();
   const Vec<int>& colsize = consMatrix.getColSizes();
   const Vec<REAL>& lhs = consMatrix.getLeftHandSides();
   const Vec<REAL>& rhs = consMatrix.getRightHandSides();

   PresolveStatus result = PresolveStatus::kUnchanged;

   for( int i = 0; i < ncols; ++i )
   {
      // skip inactive columns
      if( cflags[i].test( ColFlag::kInactive ) )
         continue;

      // if strong dual reductions are not allowed, we cannot dual fix variables
      // with zero objective
      if( problemUpdate.getPresolveOptions().dualreds < 2 && objective[i] == 0 )
         continue;

      // TODO: this can be replaced with    Vec<Locks>& locks =
      // problem.getLocks(); or nah? probably because it isn't updated?
      auto colvec = consMatrix.getColumnCoefficients( i );
      int collen = colvec.getLength();
      const REAL* values = colvec.getValues();
      const int* rowinds = colvec.getIndices();

      int nuplocks = 0;
      int ndownlocks = 0;

      // count "lock" for objective function
      if( objective[i] < 0 )
         ++ndownlocks;
      else if( objective[i] > 0 )
         ++nuplocks;

      for( int j = 0; j < collen; ++j )
      {
         count_locks( values[j], rflags[rowinds[j]], ndownlocks, nuplocks );

         if( nuplocks != 0 && ndownlocks != 0 )
            break;
      }

      // TODO add presolve result to return unbounded or infeasible when
      // column would be fixed to an infinite value

      // fix column to lower or upper bound
      if( ndownlocks == 0 )
      {
         assert( cflags[i].test( ColFlag::kUnbounded ) || ubs[i] != lbs[i] );

         // use a transaction and lock the column to protect it from changes
         // of other presolvers
         if( !cflags[i].test( ColFlag::kLbInf ) )
         {
            TransactionGuard<REAL> guard{ reductions };

            reductions.lockCol( i );
            reductions.fixCol( i, lbs[i] );

            result = PresolveStatus::kReduced;
         }
         else if( objective[i] != 0 )
         {
            return PresolveStatus::kUnbndOrInfeas;
         }
         // TODO else case
      }
      else if( nuplocks == 0 )
      {
         assert( cflags[i].test( ColFlag::kUnbounded ) || ubs[i] != lbs[i] );

         // fmt::print("col {} with bounds [{},{}] can be fixed to upper
         // bound\n", i, lbs[i], ubs[i]);

         // use a transaction and lock the column to protect it from changes
         // of other presolvers
         if( !cflags[i].test( ColFlag::kUbInf ) )
         {
            TransactionGuard<REAL> guard{ reductions };

            reductions.lockCol( i );
            reductions.fixCol( i, ubs[i] );

            result = PresolveStatus::kReduced;
         }
         else if( objective[i] != 0 )
         {
            return PresolveStatus::kUnbndOrInfeas;
         }
         // TODO else case
         // else
         // {
         //    remove variable and all constraints with a_ij != 0
         //    in postprocessing step we can always find a finit value
         // }
      }
      // apply dual substitution
      else
      {
         // If c_i >= 0, we might derive a tighter upper bound.
         // We consider only rows of
         // M := { (a_ji < 0 and rhs != inf) or (a_ji > 0 and lhs != inf)}.
         // If all constraints in M get redundant for x_i = new_UB, the upper
         // bound can be set to new_UB.
         if( objective[i] >= 0 )
         {
            const std::pair<bool, REAL>& new_upper_bound =
                calc_upper_bound_with_dual_substitution(
                    num, consMatrix, activities, lbs, ubs, rhs, lhs, i, collen,
                    values, rowinds, cflags, rflags );

            if( new_upper_bound.first )
            {
               REAL new_UB = new_upper_bound.second;

               // set new upper bound
               if( !num.isHugeVal( new_UB ) )
               {
                  assert( cflags[i].test( ColFlag::kUbInf ) ||
                          new_UB < ubs[i] );

                  // cannot detect infeasibility with this method, so at most
                  // tighten the bound to the lower bound
                  if( !cflags[i].test( ColFlag::kLbInf ) )
                     new_UB = num.max( lbs[i], new_UB );

               // A transaction is only needed to group several reductions that
               // belong together
               // TODO are the locks to strict?
               TransactionGuard<REAL> guard{ reductions };

               reductions.lockCol( i );
               reductions.lockColBounds( i );
               reductions.changeColUB( i, new_UB );
               Message::debug( this, "tightened upper bound of col {} to {}\n",
                               i, double( new_UB ) );

               result = PresolveStatus::kReduced;

               // If new upper bound is set, we continue with the next column.
               // Although, If c=0, we can try to derive an additional lower
               // bound it will conflict with the locks of this reduction and
               // hence will never be applied.
               continue;
            }
            }
         }

         // If c_i <= 0, we might derive a tighter lower bound.
         // We consider only rows of
         // M := { (a_ji > 0 and rhs != inf) or (a_ji < 0 and lhs != inf)}.
         // If all constraints in M get redundant for x_i = new_LB, the lower
         // bound can be set to new_LB.
         if( objective[i] <= 0 )
         {
            const std::pair<bool, REAL>& new_lower_bound =
                calc_lower_bound_with_dual_substitution(
                    num, consMatrix, activities, lbs, ubs, lhs, rhs, i, collen,
                    values, rowinds, cflags, rflags );

            if( !new_lower_bound.first )
               continue;
            REAL new_LB = new_lower_bound.second;
            // set new lower bound
            if( !num.isHugeVal( new_LB ) )
            {
               assert( cflags[i].test( ColFlag::kLbInf ) || new_LB > lbs[i] );

               // cannot detect infeasibility with this method, so at most
               // tighten the bound to the upper bound
               if( !cflags[i].test( ColFlag::kUbInf ) )
                  new_LB = num.min( ubs[i], new_LB );

               // A transaction is only needed to group several reductions that
               // belong together
               // TODO are any locks needed?
               TransactionGuard<REAL> guard{ reductions };

               reductions.lockCol( i );
               reductions.lockColBounds( i );
               reductions.changeColLB( i, new_LB );

               Message::debug( this, "tightened lower bound of col {} to {}\n",
                               i, double( new_LB ) );

               result = PresolveStatus::kReduced;
            }
         }
      }
   }

   return result;
}

template <typename REAL>
std::pair<bool, REAL>
DualFix<REAL>::calc_upper_bound_with_dual_substitution(
    const Num<REAL>& num, const ConstraintMatrix<REAL>& consMatrix,
    const Vec<RowActivity<REAL>>& activities, const Vec<REAL>& lbs,
    const Vec<REAL>& ubs, const Vec<REAL>& rhs, const Vec<REAL>& lhs, int i,
    int collen, const REAL* values, const int* rowinds,
    const Vec<ColFlags>& cflags, const Vec<RowFlags>& rflags )
{

   bool new_UB_init = false;
   REAL new_UB;

   for( int j = 0; j != collen; ++j )
   {
      int row = rowinds[j];

      // candidate for new upper bound
      REAL cand_bound;

      if( consMatrix.isRowRedundant( row ) )
         continue;

      if( values[j] < 0.0 && !rflags[row].test( RowFlag::kRhsInf ) &&
          activities[row].ninfmax == 0 )
      {
         cand_bound =
             ( rhs[row] - ( activities[row].max - values[j] * lbs[i] ) ) /
             values[j];
      }
      else if( values[j] > 0.0 && !rflags[row].test( RowFlag::kLhsInf ) &&
               activities[row].ninfmin == 0 )
      {
         cand_bound =
             ( lhs[row] - ( activities[row].min - values[j] * ubs[i] ) ) /
             values[j];
      }
      else if( ( values[j] < 0.0 && !rflags[row].test( RowFlag::kRhsInf ) &&
                 ( activities[row].ninfmax > 1 ||
                   ( activities[row].ninfmax == 1 &&
                     !cflags[i].test( ColFlag::kLbInf ) ) ) ) ||
               values[j] >= 0.0 && !rflags[row].test( RowFlag::kLhsInf ) &&
                   ( activities[row].ninfmin > 1 ||
                     ( activities[row].ninfmin == 1 &&
                       !cflags[i].test( ColFlag::kLbInf ) ) ) )
         return { false, 0 };
      else
         continue;

      // Only if variable is greater than or equal to new_UB, all rows
      // in M are redundant.
      // I. e. we round up for integer variables.
      if( cflags[i].test( ColFlag::kIntegral ) )
         cand_bound = num.epsCeil( cand_bound );

      if( !new_UB_init || cand_bound > new_UB )
      {
         new_UB = cand_bound;
         new_UB_init = true;

         // check if bound is already equal or worse than current bound
         // and abort in that case
         if( ( !cflags[i].test( ColFlag::kUbInf ) &&
               num.isGE( new_UB, ubs[i] ) ) ||
             new_UB >= num.getHugeVal() )
            return { false, 0 };
      }
   }
   return { true, new_UB };
}

template <typename REAL>
std::pair<bool, REAL>
DualFix<REAL>::calc_lower_bound_with_dual_substitution(
    const Num<REAL>& num, const ConstraintMatrix<REAL>& consMatrix,
    const Vec<RowActivity<REAL>>& activities, const Vec<REAL>& lbs,
    const Vec<REAL>& ubs, const Vec<REAL>& lhs, const Vec<REAL>& rhs, int i,
    int collen, const REAL* values, const int* rowinds,
    const Vec<ColFlags>& cflags, const Vec<RowFlags>& rflags )
{
   bool new_LB_init = false;
   REAL new_LB;

   for( int j = 0; j != collen; ++j )
   {
      int row = rowinds[j];
      // candidate for new lower bound
      REAL cand_bound;

      if( consMatrix.isRowRedundant( row ) )
         continue;

      if( values[j] < 0.0 && !rflags[row].test( RowFlag::kLhsInf ) &&
          activities[row].ninfmax == 0 )
      {
         cand_bound =
             ( lhs[row] - ( activities[row].min - values[j] * lbs[i] ) ) /
             values[j];
      }
      else if( values[j] > 0.0 && !rflags[row].test( RowFlag::kRhsInf ) &&
               activities[row].ninfmin == 0 )
      {
         cand_bound =
             ( rhs[row] - ( activities[row].max - values[j] * ubs[i] ) ) /
             values[j];
      }
      else if( ( values[j] < 0.0 && !rflags[row].test( RowFlag::kLhsInf ) &&
                     ( activities[row].ninfmax == 1 &&
                       !cflags[i].test( ColFlag::kUbInf ) ) ||
                 activities[row].ninfmax > 1 ) ||
               ( values[j] > 0.0 && !rflags[row].test( RowFlag::kRhsInf ) &&
                     ( activities[row].ninfmin == 1 &&
                       !cflags[i].test( ColFlag::kUbInf ) ) ||
                 activities[row].ninfmin > 1 ) )
         return { false, 0 };
      else
         continue;

      // Only if variable is less than or equal to new_LB, all rows in
      // M are redundant. I. e. we round down for integer variables.
      if( cflags[i].test( ColFlag::kIntegral ) )
         cand_bound = num.epsFloor( cand_bound );

      if( !new_LB_init || cand_bound < new_LB )
      {
         new_LB = cand_bound;
         new_LB_init = true;

         // check if bound is already equal or worse than current bound
         // and abort in that case
         if( ( !cflags[i].test( ColFlag::kLbInf ) &&
               num.isLE( new_LB, lbs[i] ) ) ||
             new_LB <= -num.getHugeVal() )
         {
            return { false, 0 };
         }
      }
   }
   return { true, new_LB };
}

} // namespace papilo

#endif
