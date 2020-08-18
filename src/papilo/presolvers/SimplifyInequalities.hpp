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

template <typename REAL>
class SimplifyInequalities : public PresolveMethod<REAL>
{
   REAL
   computeGCD( REAL gcd, REAL val, const Num<REAL>& num );

   void
   simplify( const REAL* values, const int* colinds, int rowlen,
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

   /// todo how to communicate about postsolve information
   virtual PresolveStatus
   execute( const Problem<REAL>& problem,
            const ProblemUpdate<REAL>& problemUpdate, const Num<REAL>& num,
            Reductions<REAL>& reductions ) override;
};

#ifdef PAPILO_USE_EXTERN_TEMPLATES
extern template class SimplifyInequalities<double>;
extern template class SimplifyInequalities<Quad>;
extern template class SimplifyInequalities<Rational>;
#endif

template <typename REAL>
REAL
SimplifyInequalities<REAL>::computeGCD( REAL val1, REAL val2,
                                        const Num<REAL>& num )
{
   // check if value is integral
   auto isIntegral = [&num]( REAL val, int64_t& intval ) {
      if( val > std::numeric_limits<int64_t>::max() ||
          val < std::numeric_limits<int64_t>::min() )
         return false;
      REAL intval_real = num.round( val );
      if( !num.isEq( intval_real, val ) )
         return false;
      intval = static_cast<int64_t>( intval_real );
      return true;
   };

   if( num.isZero( val1 ) || num.isZero( val2 ) )
      return 0;

   int64_t intval1;
   int64_t intval2;

   // gcd for integer values
   if( isIntegral( val1, intval1 ) && isIntegral( val2, intval2 ) )
   {
      return boost::gcd( intval1, intval2 );
   }

   // heuristic for fractional values
   // if max(abs(val1), abs(val2)) divided by d:=min(abs(val1), abs(val2)) is
   // integral, return d
   if( abs( val2 ) < abs( val1 ) )
   {
      int64_t intval3;
      if( isIntegral( val1 / val2, intval3 ) )
         return abs( val2 );
   }
   else
   {
      int64_t intval3;
      if( isIntegral( val2 / val1, intval3 ) )
         return abs( val1 );
   }
   // multiply with 600; if values are integral, return gcd
   int64_t intval4;
   int64_t intval5;
   if( isIntegral( 600 * val1, intval4 ) && isIntegral( 600 * val2, intval5 ) )
      return boost::gcd( intval4, intval5 ) / REAL{ 600 };

   // gcd not defined
   return 0;
}

template <typename REAL>
void
SimplifyInequalities<REAL>::simplify(
    const REAL* values, const int* colinds, int rowlen,
    const RowActivity<REAL>& activity, const RowFlags& rflag,
    const Vec<ColFlags>& cflags, const REAL& rhs, const REAL& lhs,
    const Vec<REAL>& lbs, const Vec<REAL>& ubs, Vec<int>& colOrder,
    Vec<int>& coeffDelete, REAL& gcd, bool& change, const Num<REAL>& num )
{
   auto maxact = activity.max;
   auto minact = activity.min;

   // 'colOrder' contains indices of 'values'; colOrder[0] is index of biggest
   // absolut coefficient in 'values' (of integer variables)

   // order variables
   for( int i = 0; i != rowlen; ++i )
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
   REAL resmaxact = maxact;
   REAL resminact = minact;
   assert( num.isGE( resmaxact, resminact ) );

   // start value important for first variable
   gcd = values[colOrder[0]];
   assert( gcd != 0 );
   REAL siderest = 0;
   bool redundant = false;
   // i is index of last non-redundant variable
   int i = 0;

   // iterate over ordered non-zero entries
   for( ; i != rowlen; ++i )
   {
      // index of variable in rowvec
      int v = colOrder[i];

      // break if variable not integral
      if( !cflags[colinds[v]].test( ColFlag::kIntegral ) )
         break;

      // update gcd
      gcd = computeGCD( gcd, values[v], num );
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
      for( int w = i + 1; w < rowlen; ++w )
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
   Vec<int> coeffDelete;

   // iterate over all constraints and try to simplify it
   for( int row = 0; row != nrows; ++row )
   {
      auto rowvec = consMatrix.getRowCoefficients( row );
      int rowlen = rowvec.getLength();
      const REAL* values = rowvec.getValues();
      const int* colinds = rowvec.getIndices();

      if( rflags[row].test( RowFlag::kRedundant ) )
         continue;
      // don't check empty or bound-constraints
      if( rowlen < 2 )
         continue;
      // cannot work with infinite activities
      if( activities[row].ninfmax != 0 || activities[row].ninfmin != 0 )
         continue;
      // consider only inequalities
      if( !rflags[row].test( RowFlag::kRhsInf, RowFlag::kLhsInf ) )
         continue;

      REAL gcd = 0;
      bool change = false;

      colOrder.clear();
      coeffDelete.clear();

      // if variables always fit into the constraint, delete them
      // e.g. x are binary and y is continuous with 0 <= y <= 1
      // 15x1 +15x2 +7x3 +3x4 +y1 <= 26
      // <=> 15x1 +15x2 <= 26  # delete variables
      // <=> x1 +x2 <=1  # divide by gcd and round right side down
      //
      // if no variables can be deleted, but the gcd of all coefficients is
      // greater than 1, round side to multiple of gcd
      // e.g. x are binary
      // 15x1 +15x2 +10x3 +5x4 <= 18
      // 15x1 +15x2 +10x3 +5x4 <= 15  # round right side down
      simplify( values, colinds, rowlen, activities[row], rflags[row], cflags,
                rhs[row], lhs[row], lbs, ubs, colOrder, coeffDelete, gcd,
                change, num );

      // simplification is possible
      if( change )
      {
         assert( gcd >= 1 );

         TransactionGuard<REAL> guard{ reductions };
         reductions.lockRow( row );
         // TODO other locks needed?

         // remove redundant variables
         for( int col : coeffDelete )
         {
            reductions.changeMatrixEntry( row, colinds[col], 0 );

            Message::debug( this, "removed variable {} in row {}\n", col, row );

            result = PresolveStatus::kReduced;
         }

         // round side to multiple of gcd; don't divide row by gcd
         if( !rflags[row].test( RowFlag::kRhsInf ) && rhs[row] != 0 )
         {
            REAL newrhs = num.feasFloor( rhs[row] / gcd ) * gcd;
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
            REAL newlhs = num.feasCeil( lhs[row] / gcd ) * gcd;
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

} // namespace papilo

#endif