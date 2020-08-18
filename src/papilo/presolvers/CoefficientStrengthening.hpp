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

#ifndef _PAPILO_PRESOLVERS_COEFFICIENT_STRENGTHENING_HPP_
#define _PAPILO_PRESOLVERS_COEFFICIENT_STRENGTHENING_HPP_

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
   }

   /// todo how to communicate about postsolve information
   virtual PresolveStatus
   execute( const Problem<REAL>& problem,
            const ProblemUpdate<REAL>& problemUpdate, const Num<REAL>& num,
            Reductions<REAL>& reductions ) override;
};

#ifdef PAPILO_USE_EXTERN_TEMPLATES
extern template class CoefficientStrengthening<double>;
extern template class CoefficientStrengthening<Quad>;
extern template class CoefficientStrengthening<Rational>;
#endif

template <typename REAL>
PresolveStatus
CoefficientStrengthening<REAL>::execute(
    const Problem<REAL>& problem, const ProblemUpdate<REAL>& problemUpdate,
    const Num<REAL>& num, Reductions<REAL>& reductions )
{
   assert( problem.getNumIntegralCols() != 0 );

   const auto& domains = problem.getVariableDomains();
   const auto& cflags = domains.flags;
   const auto& activities = problem.getRowActivities();
   const auto& changedactivities = problemUpdate.getChangedActivities();

   const auto& constMatrix = problem.getConstraintMatrix();
   const auto& lhs_values = constMatrix.getLeftHandSides();
   const auto& rhs_values = constMatrix.getRightHandSides();
   const auto& rflags = constMatrix.getRowFlags();

   PresolveStatus result = PresolveStatus::kUnchanged;

   Vec<std::pair<REAL, int>> integerCoefficients;

   for( int i : changedactivities )
   {
      auto rowcoefficients = constMatrix.getRowCoefficients( i );
      const REAL* coefficients = rowcoefficients.getValues();
      const int len = rowcoefficients.getLength();
      const int* coefindices = rowcoefficients.getIndices();

      if( ( !rflags[i].test( RowFlag::kLhsInf ) &&
            !rflags[i].test( RowFlag::kRhsInf ) ) ||
          len <= 1 )
         continue;

      REAL rhs;
      REAL maxact;
      int scale;

      // normalize constraint to a * x <= b constraint, remember if it was
      // scaled by -1
      if( !rflags[i].test( RowFlag::kLhsInf ) )
      {
         assert( rflags[i].test( RowFlag::kRhsInf ) );

         if( activities[i].ninfmin == 0 )
            maxact = -activities[i].min;
         else
            continue;

         rhs = -lhs_values[i];
         scale = -1;
      }
      else
      {
         assert( !rflags[i].test( RowFlag::kRhsInf ) );
         assert( rflags[i].test( RowFlag::kLhsInf ) );

         if( activities[i].ninfmax == 0 )
            maxact = activities[i].max;
         else
            continue;

         rhs = rhs_values[i];
         scale = 1;
      }

      // remember the integer coefficients and the maximum absolute value of
      // an integer variable
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
         continue;

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
      reductions.lockRow( i );

      if( scale == -1 )
      {
         for( const std::pair<REAL, int>& intCoef : integerCoefficients )
            reductions.changeMatrixEntry( i, intCoef.second, -intCoef.first );
         assert( rflags[i].test( RowFlag::kRhsInf ) );
         assert( !rflags[i].test( RowFlag::kLhsInf ) );

         if( lhs_values[i] != -rhs )
            reductions.changeRowLHS( i, -rhs );
      }
      else
      {
         assert( scale == 1 );
         for( const std::pair<REAL, int>& intCoef : integerCoefficients )
            reductions.changeMatrixEntry( i, intCoef.second, intCoef.first );
         assert( rflags[i].test( RowFlag::kLhsInf ) );
         assert( !rflags[i].test( RowFlag::kRhsInf ) );

         if( rhs_values[i] != rhs )
            reductions.changeRowRHS( i, rhs );
      }

      result = PresolveStatus::kReduced;
   }
   return result;
}

} // namespace papilo

#endif
