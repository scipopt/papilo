/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*               This file is part of the program and library                */
/*    PaPILO --- Parallel Presolve for Integer and Linear Optimization       */
/*                                                                           */
/* Copyright (C) 2020-2025 Zuse Institute Berlin (ZIB)                       */
/*                                                                           */
/* Licensed under the Apache License, Version 2.0 (the "License");           */
/* you may not use this file except in compliance with the License.          */
/* You may obtain a copy of the License at                                   */
/*                                                                           */
/*     http://www.apache.org/licenses/LICENSE-2.0                            */
/*                                                                           */
/* Unless required by applicable law or agreed to in writing, software       */
/* distributed under the License is distributed on an "AS IS" BASIS,         */
/* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  */
/* See the License for the specific language governing permissions and       */
/* limitations under the License.                                            */
/*                                                                           */
/* You should have received a copy of the Apache-2.0 license                 */
/* along with PaPILO; see the file LICENSE. If not visit scipopt.org.        */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef _PAPILO_PRESOLVERS_COEFFICIENT_STRENGTHENING_HPP_
#define _PAPILO_PRESOLVERS_COEFFICIENT_STRENGTHENING_HPP_

#include "papilo/Config.hpp"
#include "papilo/core/PresolveMethod.hpp"
#include "papilo/core/Problem.hpp"
#include "papilo/core/ProblemUpdate.hpp"

namespace papilo
{

template <typename REAL>
class CoefficientStrengthening : public PresolveMethod<REAL>
{
 public:
   CoefficientStrengthening() : PresolveMethod<REAL>()
   {
      this->setName( "coefftightening" );
      this->setType( PresolverType::kIntegralCols );
      this->setTiming( PresolverTiming::kFast );
      this->setArgument( ArgumentType::kSaturation );
   }

   PresolveStatus
   execute( const Problem<REAL>& problem,
            const ProblemUpdate<REAL>& problemUpdate,
            const Num<REAL>& num, Reductions<REAL>& reductions,
            const Timer& timer, int& reason_of_infeasibility) override;

 private:
   PresolveStatus
   perform_coefficient_tightening(
       const Num<REAL>& num, const VariableDomains<REAL>& domains,
       const Vec<RowActivity<REAL>>& activities, int changed_activity,
       const ConstraintMatrix<REAL>& constMatrix, const Vec<REAL>& lhs_values,
       const Vec<REAL>& rhs_values, const Vec<RowFlags>& rflags,
       const Vec<ColFlags>& cflags, Reductions<REAL>& reductions,
       Vec<std::pair<REAL, int>>& integerCoefficients ) const;
};

#ifdef PAPILO_USE_EXTERN_TEMPLATES
extern template class CoefficientStrengthening<double>;
extern template class CoefficientStrengthening<Quad>;
extern template class CoefficientStrengthening<Rational>;
#endif

template <typename REAL>
PresolveStatus
CoefficientStrengthening<REAL>::execute( const Problem<REAL>& problem,
                                         const ProblemUpdate<REAL>& problemUpdate,
                                         const Num<REAL>& num, Reductions<REAL>& reductions,
                                         const Timer& timer, int& reason_of_infeasibility){
   assert( problem.getNumIntegralCols() != 0 );

   const auto& domains = problem.getVariableDomains();
   const auto& cflags = domains.flags;
   const auto& activities = problem.getRowActivities();
   const auto& changedActivities = problemUpdate.getChangedActivities();

   const auto& constMatrix = problem.getConstraintMatrix();
   const auto& lhs_values = constMatrix.getLeftHandSides();
   const auto& rhs_values = constMatrix.getRightHandSides();
   const auto& rflags = constMatrix.getRowFlags();

   PresolveStatus result = PresolveStatus::kUnchanged;

#ifndef PAPILO_TBB
   assert( problemUpdate.getPresolveOptions().runs_sequential() );
#endif

   if( problemUpdate.getPresolveOptions().runs_sequential() ||
      !problemUpdate.getPresolveOptions().coefficient_strengthening_parallel )
   {
      Vec<std::pair<REAL, int>> integerCoefficients;
      for( int i : changedActivities )
      {
         if( reductions.size() >= problemUpdate.getPresolveOptions().max_reduction )
            break;
         if( perform_coefficient_tightening(
                 num, domains, activities, i, constMatrix, lhs_values,
                 rhs_values, rflags, cflags, reductions,
                 integerCoefficients ) == PresolveStatus::kReduced )
            result = PresolveStatus::kReduced;
      }
   }
#ifdef PAPILO_TBB
   else
   {
      Vec<Reductions<REAL>> stored_reductions( changedActivities.size() );
      tbb::parallel_for(
          tbb::blocked_range<unsigned int>( 0, changedActivities.size() ),
          [&]( const tbb::blocked_range<unsigned int>& r ) {
            Vec<std::pair<REAL, int>> integerCoefficients;
             for( unsigned int j = r.begin(); j < r.end(); ++j )
             {
                if( perform_coefficient_tightening(
                        num, domains, activities, changedActivities[j],
                        constMatrix, lhs_values, rhs_values, rflags, cflags,
                        stored_reductions[j], integerCoefficients ) == PresolveStatus::kReduced )
                   result = PresolveStatus::kReduced;
             }
          } );

      if( result == PresolveStatus::kUnchanged )
         return PresolveStatus::kUnchanged;

      for( int i = 0; i < (int) stored_reductions.size(); ++i )
      {
         Reductions<REAL> reds = stored_reductions[i];
         if( reds.size() > 0 )
         {
            for( const auto& transaction : reds.getTransactions() )
            {
               int start = transaction.start;
               int end = transaction.end;
               TransactionGuard<REAL> guard{ reductions };
               for( int c = start; c < end; c++ )
               {
                  Reduction<REAL>& reduction = reds.getReduction( c );
                  reductions.add_reduction( reduction.row,
                                            reduction.col, reduction.newval );
               }
            }
         }
      }
   }
#endif
   return result;
}

template <typename REAL>
PresolveStatus
CoefficientStrengthening<REAL>::perform_coefficient_tightening(
    const Num<REAL>& num, const VariableDomains<REAL>& domains,
    const Vec<RowActivity<REAL>>& activities, int changed_activity,
    const ConstraintMatrix<REAL>& constMatrix, const Vec<REAL>& lhs_values,
    const Vec<REAL>& rhs_values, const Vec<RowFlags>& rflags,
    const Vec<ColFlags>& cflags, Reductions<REAL>& reductions,
    Vec<std::pair<REAL, int>>& integerCoefficients ) const
{
   auto rowCoefficients = constMatrix.getRowCoefficients( changed_activity );
   const REAL* coefficients = rowCoefficients.getValues();
   const int len = rowCoefficients.getLength();
   const int* coefindices = rowCoefficients.getIndices();

   auto& rowFlag = rflags[changed_activity];
   if( ( !rowFlag.test( RowFlag::kLhsInf ) &&
         !rowFlag.test( RowFlag::kRhsInf ) ) ||
       len <= 1 )
      return PresolveStatus::kUnchanged;

   REAL rhs;
   REAL maxact;
   int scale;

   // normalize constraint to a * x <= b constraint, remember if it
   // was scaled by -1
   if( !rowFlag.test( RowFlag::kLhsInf ) )
   {
      assert( rowFlag.test( RowFlag::kRhsInf ) );

      if( activities[changed_activity].ninfmin == 0 )
         maxact = -activities[changed_activity].min;
      else
         return PresolveStatus::kUnchanged;

      rhs = -lhs_values[changed_activity];
      scale = -1;
   }
   else
   {
      assert( !rowFlag.test( RowFlag::kRhsInf ) );
      assert( rowFlag.test( RowFlag::kLhsInf ) );

      if( activities[changed_activity].ninfmax == 0 )
         maxact = activities[changed_activity].max;
      else
         return PresolveStatus::kUnchanged;

      rhs = rhs_values[changed_activity];
      scale = 1;
   }

   // remember the integer coefficients and the maximum absolute value
   // of an integer variable
   integerCoefficients.clear();
   REAL newabscoef = maxact - rhs;

   if( num.isZero( newabscoef ) )
      newabscoef = 0;
   else if( num.isFeasEq( newabscoef, ceil( newabscoef ) ) )
      newabscoef = ceil( newabscoef );

   for( int k = 0; k != len; ++k )
   {
      if( !cflags[coefindices[k]].test( ColFlag::kIntegral,
                                        ColFlag::kImplInt ) ||
          cflags[coefindices[k]].test( ColFlag::kFixed ) )
         continue;

      if( num.isFeasGE( newabscoef, abs( coefficients[k] ) ) )
         continue;

      integerCoefficients.emplace_back( coefficients[k] * scale,
                                        coefindices[k] );
   }

   if( integerCoefficients.empty() )
      return PresolveStatus::kUnchanged;

   assert( num.isFeasGT( maxact, rhs ) );

   // adjust side and qualified coefficients
   for( std::pair<REAL, int>& intCoef : integerCoefficients )
   {
      assert( domains.flags[intCoef.second].test( ColFlag::kLbInf ) ||
              domains.flags[intCoef.second].test( ColFlag::kUbInf ) ||
              domains.lower_bounds[intCoef.second] !=
                  domains.upper_bounds[intCoef.second] );
      assert( domains.flags[intCoef.second].test( ColFlag::kLbInf ) ||
              domains.flags[intCoef.second].test( ColFlag::kUbInf ) ||
              domains.upper_bounds[intCoef.second] -
                      domains.lower_bounds[intCoef.second] >=
                  1.0 );
      assert( newabscoef < abs( intCoef.first ) );

      if( intCoef.first < REAL{ 0 } )
      {
         assert( !domains.flags[intCoef.second].test( ColFlag::kLbInf ) );
         rhs -= ( newabscoef + intCoef.first ) *
                domains.lower_bounds[intCoef.second];
         intCoef.first = -newabscoef;
      }
      else
      {
         assert( !domains.flags[intCoef.second].test( ColFlag::kUbInf ) );
         rhs += ( newabscoef - intCoef.first ) *
                domains.upper_bounds[intCoef.second];
         intCoef.first = newabscoef;
      }
   }

   // add reduction to change side and coefficients
   TransactionGuard<REAL> guard{ reductions };
   reductions.lockRow( changed_activity );

   if( scale == -1 )
   {
      for( const std::pair<REAL, int>& intCoef : integerCoefficients )
         reductions.changeMatrixEntry( changed_activity, intCoef.second,
                                       -intCoef.first );
      assert( rowFlag.test( RowFlag::kRhsInf ) );
      assert( !rowFlag.test( RowFlag::kLhsInf ) );

      if( lhs_values[changed_activity] != -rhs )
         reductions.changeRowLHS( changed_activity, -rhs );
   }
   else
   {
      assert( scale == 1 );
      for( const std::pair<REAL, int>& intCoef : integerCoefficients )
         reductions.changeMatrixEntry( changed_activity, intCoef.second,
                                       intCoef.first );
      assert( rowFlag.test( RowFlag::kLhsInf ) );
      assert( !rowFlag.test( RowFlag::kRhsInf ) );

      if( rhs_values[changed_activity] != rhs )
         reductions.changeRowRHS( changed_activity, rhs );
   }
   //   stored_reductions[j] = reduction;
   return PresolveStatus::kReduced;
}

} // namespace papilo

#endif
