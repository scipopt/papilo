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

#ifndef _PAPILO_PRESOLVERS_GCD_REDUCTIONS_HPP_
#define _PAPILO_PRESOLVERS_GCD_REDUCTIONS_HPP_

#include "papilo/core/PresolveMethod.hpp"
#include "papilo/core/Problem.hpp"
#include "papilo/core/ProblemUpdate.hpp"
#include "pdqsort/pdqsort.h"
#include <boost/integer/common_factor.hpp>

namespace papilo
{

//TODO: maybe rename to Euclidian reducation for unbounded inequalities?
//since test don't show expected results, functionality of the presolver is questionable
// or the solver behaves in another manner than the testwriter expected

template <typename REAL>
class SimplifyInequalities : public PresolveMethod<REAL>
{
   REAL
   computeGreatestCommonDivisor( REAL val1, REAL val2, const Num<REAL>& num );

   void
   simplify( const REAL* values, const int* colinds, int rowLength,
             const RowActivity<REAL>& activity, const RowFlags& rflag,
             const Vec<ColFlags>& cflags, const REAL& rhs, const REAL& lhs,
             const Vec<REAL>& lbs, const Vec<REAL>& ubs, Vec<int>& colOrder,
             Vec<int>& coeffDelete, REAL& gcd, bool& change,
             const Num<REAL>& num );

 public:
   SimplifyInequalities() : PresolveMethod<REAL>()
   {
      this->setName( "simplifyineq" );
      this->setTiming( PresolverTiming::kMedium );
      this->setType( PresolverType::kIntegralCols );
   }

   virtual PresolveStatus
   execute( const Problem<REAL>& problem,
            const ProblemUpdate<REAL>& problemUpdate, const Num<REAL>& num,
            Reductions<REAL>& reductions ) override;

   bool
   isUnbounded( int row, const Vec<RowFlags>& rowFlags ) const;

   bool
   isRedundant( int row, const Vec<RowFlags>& rflags ) const;

   bool
   isInfiniteActivity( const Vec<RowActivity<REAL>>& activities,
                       int row ) const;
};

#ifdef PAPILO_USE_EXTERN_TEMPLATES
extern template class SimplifyInequalities<double>;
extern template class SimplifyInequalities<Quad>;
extern template class SimplifyInequalities<Rational>;
#endif

/***
 * to calculate the Greatest Common Divisor heuristics are used according to
 * "Presolve Reductions in Mixed Integer Programming" from T. Achterberg et. al.
 *
 * - Euclidian algorithm for integral values (numerical issues for flaoting point)
 *
 * 1. Divide all coefficients by a_min = min{|a_ij| j in supp(A_i)}. If this
 * leads to integer values for all coefficients return
 *  d= a_min * gcd(a_i1/a_min,..., a_in /a_min)
 *
 * 2. Use a_min = 1/600 (multiply by 600), because it is a multiple of many
 * small integer values that arise as denominators in real-world problems
 *
 * @tparam REAL
 * @param val1
 * @param val2
 * @param num
 * @return gcd (with heuristics for floating points)
 */
 //TODO: why are only two numbers compared? the parameter suggests that there should/can be more?
template <typename REAL>
REAL
SimplifyInequalities<REAL>::computeGreatestCommonDivisor( REAL val1, REAL val2,
                                        const Num<REAL>& num )
{
   //TODO: I removed the 2nd parameter because it is not necessary?
   auto isIntegral = [&num]( REAL val ) {
      if( val > std::numeric_limits<int64_t>::max() ||
          val < std::numeric_limits<int64_t>::min() )
         return false;
      if( !num.isEq( num.round( val ), val ) )
         return false;
      return true;
   };

   if( num.isZero( val1 ) || num.isZero( val2 ) )
      return 0;

   // gcd for integer values
   if( isIntegral( val1 ) && isIntegral( val2 ) )
   {
      return boost::gcd( static_cast<int64_t>( val1 ), static_cast<int64_t>( val2 ) );
   }

   // heuristic for fractional values
   // if max(abs(val1), abs(val2)) divided by d:=min(abs(val1), abs(val2)) is
   // integral, return d
   if( abs( val2 ) < abs( val1 ) )
   {
      if( isIntegral( val1 / val2) )
         return abs( val2 );
   }
   else
   {
      if( isIntegral( val2 / val1 ) )
         return abs( val1 );
   }

   double multiplier = 600;
   if( isIntegral( multiplier * val1 ) &&isIntegral( multiplier * val2 ) )
      return boost::gcd( static_cast<int64_t>( val1 * multiplier ) , static_cast<int64_t>( val2 * multiplier )  )
             / REAL{ multiplier };

   // applied heuristics didn't find an greatest common divisor
   return 0;
}

template <typename REAL>
void
SimplifyInequalities<REAL>::simplify(
    const REAL* values, const int* colinds, int rowLength,
    const RowActivity<REAL>& activity, const RowFlags& rflag,
    const Vec<ColFlags>& cflags, const REAL& rhs, const REAL& lhs,
    const Vec<REAL>& lbs, const Vec<REAL>& ubs, Vec<int>& colOrder,
    Vec<int>& coeffDelete, REAL& gcd, bool& change, const Num<REAL>& num )
{
   auto maxActivity = activity.max;
   auto minActivity = activity.min;

   // 'colOrder' contains indices of 'values'; colOrder[0] is index of biggest
   // absolut coefficient in 'values' (of integer variables)

   // order variables
   //TODO: hm I don't think that colOrder contains the described information.
   // I see no ordering, just inserting the indices of the rows?
   // the values of the test (first column 2; 2nd 4) are always in the order [0,1]
   for( int i = 0; i < rowLength; ++i )
   {
      colOrder.push_back( i );
   }
   // continuous variables to the end
   Vec<int>::iterator start_cont;
   start_cont = partition(
       colOrder.begin(), colOrder.end(), [&colinds, &cflags]( int const& a ) {
          return cflags[colinds[a]].test( ColFlag::kIntegral );
       } );
   // integer variables after non-increasing absolute value of the
   // coefficients
   pdqsort( colOrder.begin(), start_cont,
            [&values]( int const& a, int const& b ) {
               return abs( values[a] ) > abs( values[b] );
            } );

   // check if continuous variables or variables with small absolut value
   // always fit into the constraint
   REAL resmaxact = maxActivity;
   REAL resminact = minActivity;
   assert( num.isGE( resmaxact, resminact ) );

   // start value important for first variable
   gcd = values[colOrder[0]];
   assert( gcd != 0 );
   REAL siderest;
   bool redundant = false;
   // i is index of last non-redundant variable
   int i = 0;

   // iterate over ordered non-zero entries
   for( ; i != rowLength; ++i )
   {
      // index of variable in rowvec
      int v = colOrder[i];

      // break if variable not integral
      if( !cflags[colinds[v]].test( ColFlag::kIntegral ) )
         break;

      // update gcd
      gcd = computeGreatestCommonDivisor( gcd, values[v], num );
      if( num.isLE( gcd, 1 ) )
         break;

      assert( !cflags[colinds[v]].test( ColFlag::kLbInf, ColFlag::kUbInf ) );

      // update residual activities
      // attention: the calculation inaccuracy can be greater than epsilon
      if( values[v] > 0 )
      {
         resmaxact -= values[v] * ubs[colinds[v]];
         resminact -= values[v] * lbs[colinds[v]];
      }
      else
      {
         resmaxact -= values[v] * lbs[colinds[v]];
         resminact -= values[v] * ubs[colinds[v]];
      }

      // calculate siderest
      if( !rflag.test( RowFlag::kRhsInf ) )
      {
         siderest = rhs - num.epsFloor( rhs / gcd ) * gcd;
      }
      else
      {
         siderest = lhs - num.epsFloor( lhs / gcd ) * gcd;
         if( num.isZero( siderest ) )
            siderest = gcd;
      }

      // check if the ordered variables on the right of i are redundant
      if( ( !rflag.test( RowFlag::kRhsInf ) && resmaxact <= siderest &&
            num.isFeasLT( siderest - gcd, resminact ) ) ||
          ( !rflag.test( RowFlag::kLhsInf ) && resminact >= siderest - gcd &&
            num.isFeasGT( siderest, resmaxact ) ) )
      {
         redundant = true;
         break;
      }
   }

   if( redundant )
   {
      change = true;
      // safe indices of redundant variables
      for( int w = i + 1; w < rowLength; ++w )
      {
         coeffDelete.push_back( colOrder[w] );
      }
   }
}

template <typename REAL>
PresolveStatus
SimplifyInequalities<REAL>::execute( const Problem<REAL>& problem,
                                     const ProblemUpdate<REAL>& problemUpdate,
                                     const Num<REAL>& num,
                                     Reductions<REAL>& reductions )
{
   const auto& consMatrix = problem.getConstraintMatrix();
   const Vec<RowActivity<REAL>>& activities = problem.getRowActivities();
   const Vec<RowFlags>& rflags = consMatrix.getRowFlags();
   const Vec<ColFlags>& cflags = problem.getColFlags();
   const Vec<REAL>& lhs = consMatrix.getLeftHandSides();
   const Vec<REAL>& rhs = consMatrix.getRightHandSides();
   const int nrows = consMatrix.getNRows();
   const Vec<REAL>& lbs = problem.getLowerBounds();
   const Vec<REAL>& ubs = problem.getUpperBounds();

   PresolveStatus result = PresolveStatus::kUnchanged;

   // allocate only once
   Vec<int> colOrder;
   Vec<int> coefficientsThatCanBeDeleted;

   // iterate over all constraints and try to simplify them
   for( int row = 0; row != nrows; ++row )
   {
      auto rowvec = consMatrix.getRowCoefficients( row );
      int rowLength = rowvec.getLength();

      if( isRedundant( row, rflags ) || isUnbounded( row, rflags )  ||
          isInfiniteActivity(activities, row) ||
          // ignore empty or bound-constraints
          rowLength < 2  )
         continue;

      const int* colinds = rowvec.getIndices();

      REAL greatestCommonDivisor = 0;
      bool isSimplificationPossible = false;

      colOrder.clear();
      coefficientsThatCanBeDeleted.clear();

      // if variables always fit into the constraint, delete them
      // e.g. x are binary and y is continuous with 0 <= y <= 1
      // 15x1 +15x2 +7x3 +3x4 +y1 <= 26
      // <=> 15x1 +15x2 <= 26  # delete variables
      // <=> x1 +x2 <=1  # divide by greatestCommonDivisor and round right side down
      //
      // if no variables can be deleted, but the greatestCommonDivisor of all coefficients is
      // greater than 1, round side to multiple of greatestCommonDivisor
      // e.g. x are binary
      // 15x1 +15x2 +10x3 +5x4 <= 18
      // 15x1 +15x2 +10x3 +5x4 <= 15  # round right side down
      simplify( rowvec.getValues(), colinds, rowLength, activities[row], rflags[row], cflags,
                rhs[row], lhs[row], lbs, ubs, colOrder,
                coefficientsThatCanBeDeleted,
                greatestCommonDivisor, isSimplificationPossible, num );

      // simplification is possible
      if( isSimplificationPossible )
      {
         assert( greatestCommonDivisor >= 1 );

         TransactionGuard<REAL> guard{ reductions };
         reductions.lockRow( row );
         // TODO: are the locks sufficient?

         // remove redundant variables
         for( int col : coefficientsThatCanBeDeleted )
         {
            reductions.changeMatrixEntry( row, colinds[col], 0 );

            Message::debug( this, "removed variable {} in row {}\n", col, row );

            result = PresolveStatus::kReduced;
         }

         // round side to multiple of greatestCommonDivisor; don't divide row by greatestCommonDivisor
         if( !rflags[row].test( RowFlag::kRhsInf ) && rhs[row] != 0 )
         {
            REAL newrhs = num.feasFloor( rhs[row] / greatestCommonDivisor ) *
                          greatestCommonDivisor;
            // side is really changed
            if( newrhs != rhs[row] )
            {
               assert( rhs[row] != 0 );
               reductions.changeRowRHS( row, newrhs );

               Message::debug( this, "changed rhs of row {}\n", row );

               result = PresolveStatus::kReduced;
            }
         }
         else if( !rflags[row].test( RowFlag::kLhsInf ) && lhs[row] != 0 )
         {
            REAL newlhs = num.feasCeil( lhs[row] / greatestCommonDivisor ) *
                          greatestCommonDivisor;
            // side is really changed
            if( newlhs != lhs[row] )
            {
               assert( lhs[row] != 0 );
               reductions.changeRowLHS( row, newlhs );

               Message::debug( this, "changed lhs of row {}\n", row );

               result = PresolveStatus::kReduced;
            }
         }
      }
   }

   return result;
}

template <typename REAL>
bool
SimplifyInequalities<REAL>::isInfiniteActivity(
    const Vec<RowActivity<REAL>>& activities, int row ) const
{
   return activities[row].ninfmax != 0 || activities[row].ninfmin != 0;
}

template <typename REAL>
bool
SimplifyInequalities<REAL>::isRedundant( int row, const Vec<RowFlags>& rflags ) const
{
   return rflags[row].test( RowFlag::kRedundant );
}

template <typename REAL>
bool
SimplifyInequalities<REAL>::isUnbounded( int row, const Vec<RowFlags>& rowFlags ) const
{
   return !rowFlags[row].test( RowFlag::kRhsInf, RowFlag::kLhsInf );
}

} // namespace papilo

#endif