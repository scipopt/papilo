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

#ifndef _PAPILO_PRESOLVERS_DUAL_INFER_HPP_
#define _PAPILO_PRESOLVERS_DUAL_INFER_HPP_

#include "papilo/core/PresolveMethod.hpp"
#include "papilo/core/Problem.hpp"
#include "papilo/core/ProblemUpdate.hpp"
#include "papilo/core/SingleRow.hpp"

namespace papilo
{

template <typename REAL>
class DualInfer : public PresolveMethod<REAL>
{
 public:
   DualInfer() : PresolveMethod<REAL>()
   {
      this->setName( "dualinfer" );
      this->setTiming( PresolverTiming::EXHAUSTIVE );
      this->setType( PresolverType::CONTINUOUS_COLS );
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
extern template class DualInfer<double>;
extern template class DualInfer<Quad>;
extern template class DualInfer<Rational>;
#endif

template <typename REAL>
PresolveStatus
DualInfer<REAL>::execute( const Problem<REAL>& problem,
                          const ProblemUpdate<REAL>& problemUpdate,
                          const Num<REAL>& num, Reductions<REAL>& reductions )
{
   assert( problem.getNumContinuousCols() != 0 );

   const auto& obj = problem.getObjective().coefficients;
   const auto& consMatrix = problem.getConstraintMatrix();
   const auto& lbValues = problem.getLowerBounds();
   const auto& ubValues = problem.getUpperBounds();
   const auto& lhsValues = consMatrix.getLeftHandSides();
   const auto& rhsValues = consMatrix.getRightHandSides();
   const auto& rflags = consMatrix.getRowFlags();
   const auto& cflags = problem.getColFlags();
   const auto& rowsize = consMatrix.getRowSizes();
   const auto& colsize = consMatrix.getColSizes();
   const int ncols = problem.getNCols();
   const int nrows = problem.getNRows();
   const auto& activities = problem.getRowActivities();

   PresolveStatus result = PresolveStatus::UNCHANGED;

   // get called less and less over time regardless of success since the
   // presolver can be expensive otherwise
   this->skipRounds( this->getNCalls() );

   Vec<RowActivity<REAL>> activitiesCopy( activities );
   auto checkNonImpliedBounds = [&]( int col, bool& lbinf, bool& ubinf ) {
      const REAL& lb = lbValues[col];
      const REAL& ub = ubValues[col];
      ColFlags colf = cflags[col];

      auto colvec = consMatrix.getColumnCoefficients( col );
      const int len = colvec.getLength();
      const int* inds = colvec.getIndices();
      const REAL* vals = colvec.getValues();

      int i = 0;
      while(
          ( !colf.test( ColFlag::LB_INF ) || !colf.test( ColFlag::UB_INF ) ) &&
          i != len )
      {
         int row = inds[i];

         assert( !consMatrix.isRowRedundant( row ) );

         bool lhsinf = rflags[row].test( RowFlag::LHS_INF );
         bool rhsinf = rflags[row].test( RowFlag::RHS_INF );
         const REAL& lhs = lhsValues[row];
         const REAL& rhs = rhsValues[row];

         if( !colf.test( ColFlag::UB_INF ) )
         {
            if( row_implies_UB( num, lhs, rhs, rflags[row], activitiesCopy[row],
                                vals[i], lb, ub, colf ) )
            {
               colf.set( ColFlag::UB_INF );

               if( !colf.test( ColFlag::UB_HUGE ) )
                  update_activities_remove_finite_bound(
                      inds, vals, len, BoundChange::UPPER, ub, activitiesCopy );
            }
         }

         if( !colf.test( ColFlag::LB_INF ) )
         {
            if( row_implies_LB( num, lhs, rhs, rflags[row], activitiesCopy[row],
                                vals[i], lb, ub, colf ) )
            {
               colf.set( ColFlag::LB_INF );

               if( !colf.test( ColFlag::LB_HUGE ) )
                  update_activities_remove_finite_bound(
                      inds, vals, len, BoundChange::LOWER, lb, activitiesCopy );
            }
         }

         ++i;
      }

      lbinf = colf.test( ColFlag::LB_INF );
      ubinf = colf.test( ColFlag::UB_INF );
   };

   // initialize dual variable domains
   Vec<REAL> dualLB;
   Vec<REAL> dualUB;
   Vec<ColFlags> dualColFlags( nrows );

   dualLB.resize( nrows, REAL{0} );
   dualUB.resize( nrows, REAL{0} );

   for( int i = 0; i != nrows; ++i )
   {
      if( consMatrix.isRowRedundant( i ) )
         continue;

      assert( rowsize[i] > 0 );

      if( !rflags[i].test( RowFlag::RHS_INF ) )
         dualColFlags[i].set( ColFlag::LB_INF );
      if( !rflags[i].test( RowFlag::LHS_INF ) )
         dualColFlags[i].set( ColFlag::UB_INF );
   }

   // initialize dual constraint sides, mark constraints of integral columns
   // redundant to skip them for propagation
   const Vec<REAL>& dualLHS = obj;
   const Vec<REAL>& dualRHS = obj;
   Vec<RowFlags> dualRowFlags( ncols );

   Vec<std::pair<int, int>> checkRedundantBounds;
   checkRedundantBounds.reserve( ncols );
   const Vec<int>& colperm = problemUpdate.getRandomColPerm();
   for( int c = 0; c != ncols; ++c )
   {
      assert( !dualRowFlags[c].test( RowFlag::LHS_INF ) );
      assert( !dualRowFlags[c].test( RowFlag::RHS_INF ) );

      if( cflags[c].test( ColFlag::INTEGRAL, ColFlag::INACTIVE ) ||
          colsize[c] == 0 )
      {
         dualRowFlags[c].set( RowFlag::RHS_INF );
         dualRowFlags[c].set( RowFlag::LHS_INF );
         dualRowFlags[c].set( RowFlag::REDUNDANT );
         continue;
      }

      int64_t colweight = int64_t( !cflags[c].test( ColFlag::LB_INF ) +
                                   !cflags[c].test( ColFlag::UB_INF ) ) *
                              colsize[c] * ncols +
                          colperm[c];

      checkRedundantBounds.emplace_back(
          int( std::min( int64_t( std::numeric_limits<int>::max() ),
                         colweight ) ),
          c );
   }

   pdqsort( checkRedundantBounds.begin(), checkRedundantBounds.end(),
            []( std::pair<int, int> c1, std::pair<int, int> c2 ) {
               return c1.first < c2.first;
            } );

   for( std::pair<int, int> pair : checkRedundantBounds )
   {
      int c = pair.second;
      assert( colsize[c] > 0 );

      bool lbinf;
      bool ubinf;
      checkNonImpliedBounds( c, lbinf, ubinf );

      if( !lbinf && !ubinf )
      {
         dualRowFlags[c].set( RowFlag::RHS_INF );
         dualRowFlags[c].set( RowFlag::LHS_INF );
         dualRowFlags[c].set( RowFlag::REDUNDANT );
         continue;
      }

      if( !lbinf )
         dualRowFlags[c].set( RowFlag::LHS_INF );

      if( !ubinf )
         dualRowFlags[c].set( RowFlag::RHS_INF );
   }

   // compute initial activities, skip redundant rows
   Vec<RowActivity<REAL>> dualActivities( ncols );

   Vec<int> changedActivity;
   Vec<int> nextChangedActivity;

   auto checkRedundancy = [&]( int dualRow ) {
      if( ( dualRowFlags[dualRow].test( RowFlag::LHS_INF ) ||
            ( dualActivities[dualRow].ninfmin == 0 &&
              num.isFeasGE( dualActivities[dualRow].min,
                            dualLHS[dualRow] ) ) ) &&
          ( dualRowFlags[dualRow].test( RowFlag::RHS_INF ) ||
            ( dualActivities[dualRow].ninfmax == 0 &&
              num.isFeasLE( dualActivities[dualRow].max,
                            dualRHS[dualRow] ) ) ) )
         return true;
      return false;
   };

   for( int i = 0; i != ncols; ++i )
   {
      if( dualRowFlags[i].test( RowFlag::REDUNDANT ) )
      {
         dualActivities[i].min = 0;
         dualActivities[i].max = 0;
         dualActivities[i].ninfmin = colsize[i];
         dualActivities[i].ninfmax = colsize[i];
         continue;
      }

      auto colvec = consMatrix.getColumnCoefficients( i );

      dualActivities[i] = compute_row_activity(
          colvec.getValues(), colvec.getIndices(), colvec.getLength(), dualLB,
          dualUB, dualColFlags );

      if( checkRedundancy( i ) )
      {
         dualRowFlags[i].set( RowFlag::REDUNDANT );
         continue;
      }

      if( ( dualActivities[i].ninfmax <= 1 &&
            !dualRowFlags[i].test( RowFlag::LHS_INF ) ) ||
          ( dualActivities[i].ninfmin <= 1 &&
            !dualRowFlags[i].test( RowFlag::RHS_INF ) ) )
         changedActivity.push_back( i );
   }

   // do domain propagation
   int nrounds = 0;
   auto boundChanged = [&]( BoundChange boundChg, int dualCol, REAL newbound ) {
      auto rowvec = consMatrix.getRowCoefficients( dualCol );

      bool oldboundinf;
      REAL oldbound;

      if( boundChg == BoundChange::LOWER )
      {
         oldboundinf = dualColFlags[dualCol].test( ColFlag::LB_INF );
         oldbound = dualLB[dualCol];

         // check against other bound for infeasibility
         if( !dualColFlags[dualCol].test( ColFlag::UB_INF ) )
         {
            REAL bnddist = dualUB[dualCol] - newbound;

            // bound exceeded by more then feastol means infeasible
            if( bnddist < -num.getFeasTol() )
            {
               result = PresolveStatus::UNBND_OR_INFEAS;
               return;
            }

            // bound is equal, or almost equal
            if( bnddist <= num.getFeasTol() )
            {
               bool forcingrow = true;
               if( bnddist > 0 )
               {
                  // bound is almost equal, check if fixing the column is
                  // numerically negligible
                  auto dualColVec = consMatrix.getRowCoefficients( dualCol );
                  REAL maxabsval = dualColVec.getMaxAbsValue();
                  if( bnddist * maxabsval > num.getFeasTol() )
                     forcingrow = false;
               }

               // fix column
               if( forcingrow )
                  newbound = dualUB[dualCol];
            }
         }

         // reject too small bound changes
         if( !oldboundinf && newbound - oldbound <= +1000 * num.getFeasTol() )
            return;

         dualColFlags[dualCol].unset( ColFlag::LB_INF );
         dualLB[dualCol] = newbound;

         assert( dualColFlags[dualCol].test( ColFlag::UB_INF ) ||
                 dualLB[dualCol] <= dualUB[dualCol] );
      }
      else
      {
         oldboundinf = dualColFlags[dualCol].test( ColFlag::UB_INF );
         oldbound = dualUB[dualCol];

         // check against other bound for infeasibility
         if( !dualColFlags[dualCol].test( ColFlag::LB_INF ) )
         {
            REAL bnddist = newbound - dualLB[dualCol];

            // bound exceeded by more then feastol means infeasible
            if( bnddist < -num.getFeasTol() )
            {
               result = PresolveStatus::UNBND_OR_INFEAS;
               return;
            }

            if( bnddist <= num.getFeasTol() )
            {
               bool forcingrow = true;
               if( bnddist > 0 )
               {
                  auto dualColVec = consMatrix.getRowCoefficients( dualCol );
                  REAL maxabsval = dualColVec.getMaxAbsValue();
                  if( bnddist * maxabsval > num.getFeasTol() )
                     forcingrow = false;
               }

               if( forcingrow )
                  newbound = dualLB[dualCol];
            }
         }

         // reject too small bound changes
         if( !oldboundinf && newbound - oldbound >= -1000 * num.getFeasTol() )
            return;

         dualColFlags[dualCol].unset( ColFlag::UB_INF );
         dualUB[dualCol] = newbound;

         assert( dualColFlags[dualCol].test( ColFlag::LB_INF ) ||
                 dualUB[dualCol] >= dualLB[dualCol] );
      }

      const REAL* vals = rowvec.getValues();
      const int* inds = rowvec.getIndices();
      const int len = rowvec.getLength();

      for( int i = 0; i != len; ++i )
      {
         int dualRow = inds[i];

         if( dualRowFlags[dualRow].test( RowFlag::REDUNDANT ) )
            continue;

         RowActivity<REAL>& activity = dualActivities[dualRow];

         ActivityChange actChange = update_activity_after_boundchange(
             vals[i], boundChg, oldbound, newbound, oldboundinf, activity );

         if( checkRedundancy( dualRow ) )
         {
            dualRowFlags[dualRow].set( RowFlag::REDUNDANT );
            continue;
         }

         if( activity.lastchange != nrounds )
         {
            if( actChange == ActivityChange::MIN &&
                !dualRowFlags[dualRow].test( RowFlag::RHS_INF ) &&
                activity.ninfmin <= 1 )
            {
               activity.lastchange = nrounds;
               nextChangedActivity.push_back( dualRow );
            }
            else if( actChange == ActivityChange::MAX &&
                     !dualRowFlags[dualRow].test( RowFlag::LHS_INF ) &&
                     activity.ninfmax <= 1 )
            {
               activity.lastchange = nrounds;
               nextChangedActivity.push_back( dualRow );
            }
         }
      }
   };

   using std::swap;
   while( !changedActivity.empty() )
   {
      Message::debug( this, "dual progation round {} on {} dual rows\n",
                      nrounds, changedActivity.size() );
      for( int dualRow : changedActivity )
      {
         if( dualRowFlags[dualRow].test( RowFlag::REDUNDANT ) )
            continue;

         auto colvec = consMatrix.getColumnCoefficients( dualRow );

         propagate_row( colvec.getValues(), colvec.getIndices(),
                        colvec.getLength(), dualActivities[dualRow],
                        dualLHS[dualRow], dualRHS[dualRow],
                        dualRowFlags[dualRow], dualLB, dualUB, dualColFlags,
                        boundChanged );

         if( result == PresolveStatus::UNBND_OR_INFEAS )
            return result;
      }

      swap( changedActivity, nextChangedActivity );
      nextChangedActivity.clear();
      ++nrounds;
   }

   // analyze dual domains:
   int impliedeqs = 0;
   int fixedints = 0;
   int fixedconts = 0;

   // change inequality rows to equations if their dual value cannot be zero
   for( int i = 0; i < nrows; ++i )
   {
      if( consMatrix.isRowRedundant( i ) )
         continue;

      if( !rflags[i].test( RowFlag::LHS_INF ) &&
          !rflags[i].test( RowFlag::RHS_INF ) && lhsValues[i] == rhsValues[i] )
         continue;

      if( ( !dualColFlags[i].test( ColFlag::LB_INF ) &&
            num.isFeasGT( dualLB[i], 0 ) ) )
      {
         assert( !rflags[i].test( RowFlag::LHS_INF ) );

         if( activities[i].ninfmax != 0 ||
             num.isFeasLT( lhsValues[i], activities[i].max ) )
         {
            TransactionGuard<REAL> tg{reductions};
            reductions.lockRow( i );
            reductions.changeRowRHS( i, lhsValues[i] );
            ++impliedeqs;
         }
      }
      else if( !dualColFlags[i].test( ColFlag::UB_INF ) &&
               num.isFeasLT( dualUB[i], 0 ) )
      {
         assert( !rflags[i].test( RowFlag::RHS_INF ) );

         if( activities[i].ninfmin != 0 ||
             num.isFeasGT( rhsValues[i], activities[i].min ) )
         {
            TransactionGuard<REAL> tg{reductions};
            reductions.lockRow( i );
            reductions.changeRowLHS( i, rhsValues[i] );
            ++impliedeqs;
         }
      }
   }

   // use dualActivities to compute reduced cost bounds and use them to fix
   // variables
   for( int i = 0; i < ncols; ++i )
   {
      if( colsize[i] <= 0 )
         continue;

      if( dualRowFlags[i].test( RowFlag::REDUNDANT ) )
      {
         auto colvec = consMatrix.getColumnCoefficients( i );

         dualActivities[i] = compute_row_activity(
             colvec.getValues(), colvec.getIndices(), colvec.getLength(),
             dualLB, dualUB, dualColFlags );
      }

      if( dualActivities[i].ninfmax == 0 &&
          num.isFeasLT( dualActivities[i].max, obj[i] ) &&
          num.isSafeLT( dualActivities[i].max, obj[i] ) )
      {
         assert( !cflags[i].test( ColFlag::LB_INF ) );
         TransactionGuard<REAL> tg{reductions};
         reductions.lockColBounds( i );
         reductions.fixCol( i, lbValues[i] );

         if( cflags[i].test( ColFlag::INTEGRAL ) )
            ++fixedints;
         else
            ++fixedconts;
      }
      else if( dualActivities[i].ninfmin == 0 &&
               num.isFeasGT( dualActivities[i].min, obj[i] ) &&
               num.isSafeGT( dualActivities[i].min, obj[i] ) )
      {
         assert( !cflags[i].test( ColFlag::UB_INF ) );
         TransactionGuard<REAL> tg{reductions};
         reductions.lockColBounds( i );
         reductions.fixCol( i, ubValues[i] );

         if( cflags[i].test( ColFlag::INTEGRAL ) )
            ++fixedints;
         else
            ++fixedconts;
      }
   }

   // set result if reductions where found
   if( impliedeqs > 0 || fixedints > 0 || fixedconts > 0 )
      result = PresolveStatus::REDUCED;

   // if( nrounds != 0 && result == PresolveStatus::REDUCED )
   //    fmt::print(
   //        "   Dual infer: {} propagation rounds, {} implied equations, "
   //        "{} integers fixed, {} conts fixed\n",
   //        nrounds, impliedeqs, fixedints, fixedconts );

   return result;
}

} // namespace papilo

#endif
