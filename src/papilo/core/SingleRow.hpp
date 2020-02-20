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

#ifndef _CORE_SINGLE_ROW_HPP_
#define _CORE_SINGLE_ROW_HPP_

#include "papilo/core/RowFlags.hpp"
#include "papilo/core/VariableDomains.hpp"
#include "papilo/misc/Flags.hpp"
#include "papilo/misc/Num.hpp"
#include "papilo/misc/Vec.hpp"
#include <tuple>

namespace papilo
{

enum class BoundChange
{
   LOWER,
   UPPER
};

enum class ActivityChange
{
   MIN,
   MAX
};

enum class RowStatus
{
   INFEASIBLE,
   REDUNDANT,
   REDUNDANT_LHS,
   REDUNDANT_RHS,
   UNKNOWN,
};

template <typename REAL>
struct RowActivity
{
   /// minimal activity of the row
   REAL min;

   /// maximal activity of the row
   REAL max;

   /// number of variables that contribute with an infinite obund to the minimal
   /// activity of this row
   int ninfmin;

   /// number of variables that contribute with an infinite obund to the maximal
   /// activity of this row
   int ninfmax;

   /// last presolving round where this activity changed
   int lastchange;

   bool
   repropagate( ActivityChange actChange, RowFlags rflags )
   {
      if( actChange == ActivityChange::MIN &&
          !rflags.test( RowFlag::RHS_INF ) && ninfmin <= 1 )
         return true;

      if( actChange == ActivityChange::MAX &&
          !rflags.test( RowFlag::LHS_INF ) && ninfmax <= 1 )
         return true;

      return false;
   }

   RowStatus
   checkStatus( const Num<REAL>& num, RowFlags rflags, const REAL& lhs,
                const REAL& rhs ) const
   {
      RowStatus status = RowStatus::REDUNDANT;

      if( !rflags.test( RowFlag::LHS_INF ) )
      {
         if( ninfmax == 0 && num.isFeasLT( max, lhs ) &&
             num.isSafeLT( max, lhs ) )
            return RowStatus::INFEASIBLE;

         if( ninfmin == 0 && num.isFeasGE( min, lhs ) )
            status = RowStatus::REDUNDANT_LHS;
         else
            status = RowStatus::UNKNOWN;
      }

      if( !rflags.test( RowFlag::RHS_INF ) )
      {
         if( ninfmin == 0 && num.isFeasGT( min, rhs ) &&
             num.isSafeGT( min, rhs ) )
            return RowStatus::INFEASIBLE;

         if( ninfmax == 0 && num.isFeasLE( max, rhs ) )
         {
            if( status == RowStatus::UNKNOWN )
               status = RowStatus::REDUNDANT_RHS;
            else
               status = RowStatus::REDUNDANT;
         }
         else if( status == RowStatus::REDUNDANT )
            status = RowStatus::UNKNOWN;
      }
      else if( status == RowStatus::REDUNDANT_LHS )
         status = RowStatus::REDUNDANT;

      return status;
   }

   template <typename Archive>
   void
   serialize( Archive& ar, const unsigned int version )
   {
      ar& min;
      ar& max;
      ar& ninfmin;

      ar& ninfmax;
      ar& lastchange;
   }

   RowActivity() {}
};

/// counts the locks for the given row entry
template <typename REAL>
void
count_locks( const REAL& val, RowFlags rflags, int& ndownlocks, int& nuplocks )
{
   assert( val != 0 );

   if( val < 0 )
   {
      if( !rflags.test( RowFlag::LHS_INF ) )
         ++nuplocks;

      if( !rflags.test( RowFlag::RHS_INF ) )
         ++ndownlocks;
   }
   else
   {
      if( !rflags.test( RowFlag::LHS_INF ) )
         ++ndownlocks;

      if( !rflags.test( RowFlag::RHS_INF ) )
         ++nuplocks;
   }
}

/// computes activity of single row from scratch
template <typename REAL>
RowActivity<REAL>
compute_row_activity( const REAL* rowvals, const int* colindices, int rowlen,
                      const Vec<REAL>& lower_bounds,
                      const Vec<REAL>& upper_bounds, const Vec<ColFlags>& flags,
                      int presolveround = -1 )
{
   RowActivity<REAL> activity;

   activity.min = 0.0;
   activity.max = 0.0;
   activity.ninfmin = 0;
   activity.ninfmax = 0;
   activity.lastchange = presolveround;

   for( int j = 0; j < rowlen; ++j )
   {
      int col = colindices[j];
      if( !flags[col].test( ColFlag::UB_USELESS ) )
      {
         if( rowvals[j] < 0 )
            activity.min += rowvals[j] * upper_bounds[col];
         else
            activity.max += rowvals[j] * upper_bounds[col];
      }
      else
      {
         assert( flags[col].test( ColFlag::UB_USELESS ) );
         if( rowvals[j] < 0 )
            ++activity.ninfmin;
         else
            ++activity.ninfmax;
      }

      if( !flags[col].test( ColFlag::LB_USELESS ) )
      {
         if( rowvals[j] < 0 )
            activity.max += rowvals[j] * lower_bounds[col];
         else
            activity.min += rowvals[j] * lower_bounds[col];
      }
      else
      {
         assert( flags[col].test( ColFlag::LB_USELESS ) );
         if( rowvals[j] < 0 )
            ++activity.ninfmax;
         else
            ++activity.ninfmin;
      }
   }

   return activity;
}

/// update the vector of row activities after lower or upper bounds of a column
/// changed. The last argument must be callable with arguments (ActivityChange,
/// rowid, RowActivity) and is called to inform about row activities that
/// changed
template <typename REAL>
ActivityChange
update_activity_after_boundchange( const REAL& colval, BoundChange type,
                                   const REAL& oldbound, const REAL& newbound,
                                   bool oldbound_inf,
                                   RowActivity<REAL>& activity )
{
   assert( oldbound_inf ||
           ( type == BoundChange::LOWER && newbound != oldbound ) ||
           ( type == BoundChange::UPPER && newbound != oldbound ) );

   if( type == BoundChange::LOWER )
   {
      if( colval < REAL{0.0} )
      {
         if( oldbound_inf )
         {
            assert( activity.ninfmax > 0 );
            --activity.ninfmax;

            activity.max += newbound * colval;
         }
         else
         {
            activity.max += ( newbound - oldbound ) * colval;
         }

         return ActivityChange::MAX;
      }
      else
      {
         if( oldbound_inf )
         {
            assert( activity.ninfmin > 0 );
            --activity.ninfmin;

            activity.min += newbound * colval;
         }
         else
         {
            activity.min += ( newbound - oldbound ) * colval;
         }

         return ActivityChange::MIN;
      }
   }
   else
   {
      if( colval < REAL{0.0} )
      {
         if( oldbound_inf )
         {
            assert( activity.ninfmin > 0 );
            --activity.ninfmin;

            activity.min += newbound * colval;
         }
         else
         {
            activity.min += ( newbound - oldbound ) * colval;
         }

         return ActivityChange::MIN;
      }
      else
      {
         if( oldbound_inf )
         {
            assert( activity.ninfmax > 0 );
            --activity.ninfmax;

            activity.max += newbound * colval;
         }
         else
         {
            activity.max += ( newbound - oldbound ) * colval;
         }

         return ActivityChange::MAX;
      }
   }
}

/// update the vector of row activities after removing a finite lower or upper
/// bound of a column
template <typename REAL>
void
update_activities_remove_finite_bound( const int* colinds, const REAL* colvals,
                                       int collen, BoundChange type,
                                       const REAL& oldbound,
                                       Vec<RowActivity<REAL>>& activities )
{
   if( type == BoundChange::LOWER )
   {
      for( int i = 0; i != collen; ++i )
      {
         const REAL& colval = colvals[i];
         RowActivity<REAL>& activity = activities[colinds[i]];

         if( colval < REAL{0.0} )
         {
            activity.max -= oldbound * colval;
            ++activity.ninfmax;
         }
         else
         {
            activity.min -= oldbound * colval;
            ++activity.ninfmin;
         }
      }
   }
   else
   {
      for( int i = 0; i != collen; ++i )
      {
         const REAL& colval = colvals[i];
         RowActivity<REAL>& activity = activities[colinds[i]];

         if( colval < REAL{0.0} )
         {
            activity.min -= oldbound * colval;
            ++activity.ninfmin;
         }
         else
         {
            activity.max -= oldbound * colval;
            ++activity.ninfmax;
         }
      }
   }
}

/// update the vector of row activities after lower or upper bounds of a column
/// changed. The last argument must be callable with arguments (ActivityChange,
/// rowid, RowActivity) and is called to inform about row activities that
/// changed
template <typename REAL, typename ACTIVITYCHANGE>
void
update_activities_after_boundchange( const REAL* colvals, const int* colrows,
                                     int collen, BoundChange type,
                                     REAL oldbound, REAL newbound,
                                     bool oldbound_inf,
                                     Vec<RowActivity<REAL>>& activities,
                                     ACTIVITYCHANGE&& activityChange,
                                     bool watchInfiniteActivities = false )
{
   assert( oldbound_inf ||
           ( type == BoundChange::LOWER && newbound != oldbound ) ||
           ( type == BoundChange::UPPER && newbound != oldbound ) );

   for( int i = 0; i < collen; ++i )
   {
      RowActivity<REAL>& activity = activities[colrows[i]];

      ActivityChange actChange = update_activity_after_boundchange(
          colvals[i], type, oldbound, newbound, oldbound_inf, activity );

      if( actChange == ActivityChange::MIN &&
          ( activity.ninfmin == 0 || watchInfiniteActivities ) )
         activityChange( ActivityChange::MIN, colrows[i], activity );

      if( actChange == ActivityChange::MAX &&
          ( activity.ninfmax == 0 || watchInfiniteActivities ) )
         activityChange( ActivityChange::MAX, colrows[i], activity );
   }
}

/// update the vector of row activities after the coefficient of a column within
/// a row changed. The last argument must be callable with arguments
/// (ActivityChange, RowActivity) and is called to inform about row activities
/// that changed
template <typename REAL, typename ACTIVITYCHANGE>
void
update_activities_after_coeffchange( REAL collb, REAL colub, ColFlags cflags,
                                     REAL oldcolcoef, REAL newcolcoef,
                                     RowActivity<REAL>& activity,
                                     ACTIVITYCHANGE&& activityChange )
{
   assert( oldcolcoef != newcolcoef );

   if( oldcolcoef * newcolcoef <= 0.0 )
   { // the sign of the coefficient flipped, so the column bounds now contribute
     // to the opposite activity bound

      // remember old activity
      RowActivity<REAL> oldactivity = activity;

      if( oldcolcoef != 0.0 )
      { // if the old coefficient was not 0.0 we remove its contributions to the
        // minimum and maximum activity
         // remove old contributions of the lower bound
         if( cflags.test( ColFlag::LB_USELESS ) )
         {
            if( oldcolcoef < 0.0 )
               --activity.ninfmax;
            else
               --activity.ninfmin;
         }
         else
         {
            if( oldcolcoef < 0.0 )
               activity.max -= oldcolcoef * collb;
            else
               activity.min -= oldcolcoef * collb;
         }

         // remove old contributions of the upper bound
         if( cflags.test( ColFlag::UB_USELESS ) )
         {
            if( oldcolcoef < 0.0 )
               --activity.ninfmin;
            else
               --activity.ninfmax;
         }
         else
         {
            if( oldcolcoef < 0.0 )
               activity.min -= oldcolcoef * colub;
            else
               activity.max -= oldcolcoef * colub;
         }
      }

      if( newcolcoef != 0.0 )
      { // if the new coefficient is not 0.0 we add its contributions to the
        // minimum and maximum activity
         // add new contributions of the lower bound
         if( cflags.test( ColFlag::LB_USELESS ) )
         {
            if( newcolcoef < 0.0 )
               ++activity.ninfmax;
            else
               ++activity.ninfmin;
         }
         else
         {
            if( newcolcoef < 0.0 )
               activity.max += newcolcoef * collb;
            else
               activity.min += newcolcoef * collb;
         }

         // addnewold contributions of the upper bound
         if( cflags.test( ColFlag::UB_USELESS ) )
         {
            if( newcolcoef < 0.0 )
               ++activity.ninfmin;
            else
               ++activity.ninfmax;
         }
         else
         {
            if( newcolcoef < 0.0 )
               activity.min += newcolcoef * colub;
            else
               activity.max += newcolcoef * colub;
         }
      }

      if( ( oldactivity.ninfmin != 0 && activity.ninfmin == 0 ) ||
          ( oldactivity.ninfmin == 0 && activity.ninfmin == 0 &&
            oldactivity.min != activity.min ) )
         activityChange( ActivityChange::MIN, activity );

      if( ( oldactivity.ninfmax != 0 && activity.ninfmax == 0 ) ||
          ( oldactivity.ninfmax == 0 && activity.ninfmax == 0 &&
            oldactivity.max != activity.max ) )
         activityChange( ActivityChange::MAX, activity );
   }
   else
   { // the sign of the coefficient did not flip, so the column bounds still
     // contribute to the same activity bound
      if( !cflags.test( ColFlag::LB_USELESS ) && collb != 0.0 )
      {
         if( newcolcoef < REAL{0.0} )
         {
            activity.max += collb * ( newcolcoef - oldcolcoef );
            if( activity.ninfmax == 0 )
               activityChange( ActivityChange::MAX, activity );
         }
         else
         {
            activity.min += collb * ( newcolcoef - oldcolcoef );
            if( activity.ninfmin == 0 )
               activityChange( ActivityChange::MIN, activity );
         }
      }

      if( !cflags.test( ColFlag::UB_USELESS ) && colub != 0.0 )
      {
         if( newcolcoef < REAL{0.0} )
         {
            activity.min += colub * ( newcolcoef - oldcolcoef );
            if( activity.ninfmin == 0 )
               activityChange( ActivityChange::MIN, activity );
         }
         else
         {
            activity.max += colub * ( newcolcoef - oldcolcoef );
            if( activity.ninfmax == 0 )
               activityChange( ActivityChange::MAX, activity );
         }
      }
   }
}

/// propagate domains of variables using the given a row and its activity. The
/// last argument must be callable with arguments (BoundChange, colid, newbound)
/// and is called to inform about column bounds that changed.
template <typename REAL, typename BOUNDCHANGE>
void
propagate_row( const REAL* rowvals, const int* colindices, int rowlen,
               const RowActivity<REAL>& activity, REAL lhs, REAL rhs,
               RowFlags rflags, const Vec<REAL>& lower_bounds,
               const Vec<REAL>& upper_bounds, const Vec<ColFlags>& domainFlags,
               BOUNDCHANGE&& boundchange )
{
   if( activity.ninfmin == 1 && activity.ninfmax == 0 &&
       rflags.test( RowFlag::RHS_INF ) )
      rhs = activity.max;

   if( activity.ninfmax == 1 && activity.ninfmin == 0 &&
       rflags.test( RowFlag::LHS_INF ) )
      lhs = activity.min;

   if( !rflags.test( RowFlag::RHS_INF ) && activity.ninfmin <= 1 )
   {
      for( int j = 0; j < rowlen; ++j )
      {
         int col = colindices[j];
         REAL lb = lower_bounds[col];
         REAL ub = upper_bounds[col];
         REAL minresact = activity.min;
         REAL val = rowvals[j];

         if( val < REAL{0.0} )
         {
            if( activity.ninfmin == 1 )
            {
               if( !domainFlags[col].test( ColFlag::UB_USELESS ) )
                  continue;

               j = rowlen;
            }
            else
            {
               assert( !domainFlags[col].test( ColFlag::UB_USELESS ) );
               minresact -= val * ub;
            }

            REAL newlb = ( rhs - minresact ) / val;
            if( domainFlags[col].test( ColFlag::LB_INF ) || newlb > lb )
               boundchange( BoundChange::LOWER, col, newlb );
         }
         else
         {
            if( activity.ninfmin == 1 )
            {
               if( !domainFlags[col].test( ColFlag::LB_USELESS ) )
                  continue;

               j = rowlen;
            }
            else
            {
               assert( !domainFlags[col].test( ColFlag::LB_USELESS ) );
               minresact -= val * lb;
            }

            REAL newub = ( rhs - minresact ) / val;
            if( domainFlags[col].test( ColFlag::UB_INF ) || newub < ub )
               boundchange( BoundChange::UPPER, col, newub );
         }
      }
   }

   if( !rflags.test( RowFlag::LHS_INF ) && activity.ninfmax <= 1 )
   {
      for( int j = 0; j < rowlen; ++j )
      {
         int col = colindices[j];
         REAL lb = lower_bounds[col];
         REAL ub = upper_bounds[col];
         REAL maxresact = activity.max;
         REAL val = rowvals[j];

         if( val < REAL{0.0} )
         {
            if( activity.ninfmax == 1 )
            {
               if( !domainFlags[col].test( ColFlag::LB_USELESS ) )
                  continue;

               j = rowlen;
            }
            else
            {
               assert( !domainFlags[col].test( ColFlag::LB_USELESS ) );
               maxresact -= val * lb;
            }

            REAL newub = ( lhs - maxresact ) / val;
            if( domainFlags[col].test( ColFlag::UB_INF ) || newub < ub )
               boundchange( BoundChange::UPPER, col, newub );
         }
         else
         {
            if( activity.ninfmax == 1 )
            {
               if( !domainFlags[col].test( ColFlag::UB_USELESS ) )
                  continue;

               j = rowlen;
            }
            else
            {
               assert( !domainFlags[col].test( ColFlag::UB_USELESS ) );
               maxresact -= val * ub;
            }

            REAL newlb = ( lhs - maxresact ) / val;
            if( domainFlags[col].test( ColFlag::LB_INF ) || newlb > lb )
               boundchange( BoundChange::LOWER, col, newlb );
         }
      }
   }
}

template <typename REAL>
bool
row_implies_LB( const Num<REAL>& num, REAL lhs, REAL rhs, RowFlags rflags,
                const RowActivity<REAL>& activity, REAL colcoef, REAL collb,
                REAL colub, ColFlags cflags )

{
   if( cflags.test( ColFlag::LB_INF ) )
      return true;

   REAL resact;
   REAL side;

   if( colcoef > 0.0 && !rflags.test( RowFlag::LHS_INF ) )
   {
      if( activity.ninfmax == 0 )
      {
         assert( !cflags.test( ColFlag::UB_USELESS ) );
         resact = activity.max - colub * colcoef;
      }
      else if( activity.ninfmax == 1 && cflags.test( ColFlag::UB_USELESS ) )
         resact = activity.max;
      else
         return false;

      side = lhs;
   }
   else if( colcoef < 0.0 && !rflags.test( RowFlag::RHS_INF ) )
   {
      if( activity.ninfmin == 0 )
      {
         assert( !cflags.test( ColFlag::UB_USELESS ) );
         resact = activity.min - colub * colcoef;
      }
      else if( activity.ninfmin == 1 && cflags.test( ColFlag::UB_USELESS ) )
         resact = activity.min;
      else
         return false;

      side = rhs;
   }
   else
      return false;

   return num.isFeasGE( ( side - resact ) / colcoef, collb );
}

template <typename REAL>
bool
row_implies_UB( const Num<REAL>& num, REAL lhs, REAL rhs, RowFlags rflags,
                const RowActivity<REAL>& activity, REAL colcoef, REAL collb,
                REAL colub, ColFlags cflags )
{
   if( cflags.test( ColFlag::UB_INF ) )
      return true;

   REAL resact;
   REAL side;

   if( colcoef > 0.0 && !rflags.test( RowFlag::RHS_INF ) )
   {
      if( activity.ninfmin == 0 )
      {
         assert( !cflags.test( ColFlag::LB_USELESS ) );
         resact = activity.min - collb * colcoef;
      }
      else if( activity.ninfmin == 1 && cflags.test( ColFlag::LB_USELESS ) )
         resact = activity.min;
      else
         return false;

      side = rhs;
   }
   else if( colcoef < 0.0 && !rflags.test( RowFlag::LHS_INF ) )
   {
      if( activity.ninfmax == 0 )
      {
         assert( !cflags.test( ColFlag::LB_USELESS ) );
         resact = activity.max - collb * colcoef;
      }
      else if( activity.ninfmax == 1 && cflags.test( ColFlag::LB_USELESS ) )
         resact = activity.max;
      else
         return false;

      side = lhs;
   }
   else
      return false;

   return num.isFeasLE( ( side - resact ) / colcoef, colub );
}

} // namespace papilo

#endif
