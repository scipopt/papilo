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

   /// todo how to communicate about postsolve information
   virtual PresolveStatus
   execute( const Problem<REAL>& problem,
            const ProblemUpdate<REAL>& problemUpdate, const Num<REAL>& num,
            Reductions<REAL>& reductions ) override;
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

   for( int i = 0; i != ncols; ++i )
   {
      // skip inactive columns
      if( cflags[i].test( ColFlag::kInactive ) )
         continue;

      // if strong dual reductions are not allowed, we cannot dual fix variables
      // with zero objective
      if( problemUpdate.getPresolveOptions().dualreds < 2 && objective[i] == 0 )
         continue;

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

      for( int j = 0; j != collen; ++j )
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
      // try to derive tighter bounds
      else
      {
         // Function checks if considered row allows dual bound strengthening
         // and calculates tightest bound for this row.
         auto check_row = []( int ninf, REAL activity, const REAL& side,
                              const REAL& coeff, const REAL& boundval,
                              bool boundinf, bool& skip, REAL& cand_bound ) {
            switch( ninf )
            {
            case 0:
               assert( !boundinf );
               // calculate residual activity
               activity -= boundval * coeff;
               break;
            case 1:
               if( boundinf )
                  break;
            default:
               // If one of the other variables with non-zero entry is
               // unbounded, dual bound strengthening is not possible for this
               // column; skip column.
               skip = true;
               return;
            }

            // calculate candidate for new bound
            cand_bound = ( side - activity ) / coeff;
         };

         // If c_i >= 0, we might derive a tighter upper bound.
         // We consider only rows of
         // M := { (a_ji < 0 and rhs != inf) or (a_ji > 0 and lhs != inf)}.
         // If all constraints in M get redundant for x_i = new_UB, the upper
         // bound can be set to new_UB.
         if( objective[i] >= 0 )
         {
            bool skip = false;
            bool new_UB_init = false;
            REAL new_UB;

            // go through all rows with non-zero entry
            for( int j = 0; j != collen; ++j )
            {
               int row = rowinds[j];
               // candidate for new upper bound
               REAL cand_bound;

               if( consMatrix.isRowRedundant( row ) )
                  continue;

               // if row is in M, calculate candidate for new upper bound
               if( values[j] < 0.0 )
               {
                  if( !rflags[row].test( RowFlag::kRhsInf ) )
                  {
                     check_row( activities[row].ninfmax, activities[row].max,
                                rhs[row], values[j], lbs[i],
                                cflags[i].test( ColFlag::kLbInf ), skip,
                                cand_bound );
                  }
                  else
                     // row is not in M
                     continue;
               }
               else
               {
                  if( !rflags[row].test( RowFlag::kLhsInf ) )
                  {
                     check_row( activities[row].ninfmin, activities[row].min,
                                lhs[row], values[j], lbs[i],
                                cflags[i].test( ColFlag::kLbInf ), skip,
                                cand_bound );
                  }
                  else
                     // row is not in M
                     continue;
               }

               if( skip == true )
                  break;

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
                  {
                     skip = true;
                     break;
                  }
               }
            }

            // set new upper bound
            if( !skip && new_UB_init && !num.isHugeVal( new_UB ) )
            {
               assert( cflags[i].test( ColFlag::kUbInf ) || new_UB < ubs[i] );

               // cannot detect infeasibility with this method, so at most
               // tighten the bound to the lower bound
               if( !cflags[i].test( ColFlag::kLbInf ) )
                  new_UB = num.max( lbs[i], new_UB );

               // A transaction is only needed to group several reductions that
               // belong together
               // TODO are any locks needed?
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

         // If c_i <= 0, we might derive a tighter lower bound.
         // We consider only rows of
         // M := { (a_ji > 0 and rhs != inf) or (a_ji < 0 and lhs != inf)}.
         // If all constraints in M get redundant for x_i = new_LB, the lower
         // bound can be set to new_LB.
         if( objective[i] <= 0 )
         {
            bool skip = false;
            bool new_LB_init = false;
            REAL new_LB;

            // go through all rows with non-zero entry
            for( int j = 0; j != collen; ++j )
            {
               int row = rowinds[j];
               // candidate for new lower bound
               REAL cand_bound;

               if( consMatrix.isRowRedundant( row ) )
                  continue;

               // if row is in M, calculate candidate for new lower bound
               if( values[j] > 0.0 )
               {
                  if( !rflags[row].test( RowFlag::kRhsInf ) )
                  {
                     check_row( activities[row].ninfmax, activities[row].max,
                                rhs[row], values[j], ubs[i],
                                cflags[i].test( ColFlag::kUbInf ), skip,
                                cand_bound );
                  }
                  else
                     // row is not in M
                     continue;
               }
               else
               {
                  if( !rflags[row].test( RowFlag::kLhsInf ) )
                  {
                     check_row( activities[row].ninfmin, activities[row].min,
                                lhs[row], values[j], ubs[i],
                                cflags[i].test( ColFlag::kUbInf ), skip,
                                cand_bound );
                  }
                  else
                     // row is not in M
                     continue;
               }

               if( skip == true )
                  break;

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
                     skip = true;
                     break;
                  }
               }
            }

            // set new lower bound
            if( !skip && new_LB_init && !num.isHugeVal( new_LB ) )
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

} // namespace papilo

#endif
