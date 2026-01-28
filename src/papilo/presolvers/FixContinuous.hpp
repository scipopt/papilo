/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*               This file is part of the program and library                */
/*    PaPILO --- Parallel Presolve for Integer and Linear Optimization       */
/*                                                                           */
/* Copyright (C) 2020-2026 Zuse Institute Berlin (ZIB)                       */
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

#ifndef _PAPILO_PRESOLVERS_FIX_CONTINUOUS_HPP_
#define _PAPILO_PRESOLVERS_FIX_CONTINUOUS_HPP_

#include "papilo/Config.hpp"
#include "papilo/core/PresolveMethod.hpp"
#include "papilo/core/Problem.hpp"
#include "papilo/core/ProblemUpdate.hpp"

namespace papilo
{

/// presolver to fix continuous variables whose bounds are very close
template <typename REAL>
class FixContinuous : public PresolveMethod<REAL>
{
 public:
   FixContinuous() : PresolveMethod<REAL>()
   {
      this->setName( "fixcontinuous" );
      this->setTiming( PresolverTiming::kMedium );
      this->setType( PresolverType::kContinuousCols );
   }

   virtual PresolveStatus
   execute( const Problem<REAL>& problem,
            const ProblemUpdate<REAL>& problemUpdate,
            const Num<REAL>& num, Reductions<REAL>& reductions,
            const Timer& timer, int& reason_of_infeasibility)
       override;
};

#ifdef PAPILO_USE_EXTERN_TEMPLATES
extern template class FixContinuous<double>;
extern template class FixContinuous<Quad>;
extern template class FixContinuous<Rational>;
#endif

template <typename REAL>
PresolveStatus
FixContinuous<REAL>::execute( const Problem<REAL>& problem,
                              const ProblemUpdate<REAL>& problemUpdate,
                              const Num<REAL>& num, Reductions<REAL>& reductions,
                              const Timer& timer, int& reason_of_infeasibility){
   assert( problem.getNumContinuousCols() != 0 );

   const auto& consMatrix = problem.getConstraintMatrix();
   const auto& domains = problem.getVariableDomains();
   const auto& cflags = domains.flags;
   const auto& objective = problem.getObjective();
   const auto& lbs = problem.getLowerBounds();
   const auto& ubs = problem.getUpperBounds();

   const int ncols = consMatrix.getNCols();

   PresolveStatus result = PresolveStatus::kUnchanged;

   if( num.getEpsilon() == 0 )
      return result;

   for( int i = 0; i < ncols; ++i )
   {
      if( reductions.size() >= problemUpdate.getPresolveOptions().max_reduction_seq )
         break;
      // do not fix columns which are inactive, unbounded, integral,
      // or have feasibly distinct bounds
      if( cflags[i].test( ColFlag::kInactive, ColFlag::kUnbounded, ColFlag::kIntegral )
          || num.isFeasLT( lbs[i], ubs[i] ) )
         continue;

      assert( consMatrix.getColSizes()[i] >= 0 );
      assert( lbs[i] < ubs[i] );

      // if the change in activity due to fixing this column is at most
      // epsilon in every row we can fix it
      if( num.isLE( (ubs[i] - lbs[i]) * num.max( abs( objective.coefficients[i] ),
         consMatrix.getColumnCoefficients( i ).getMaxAbsValue() ), 0 ) )
      {
         REAL fixval;

         // if one bound is an integral value, use that one
         if( floor( ubs[i] ) == lbs[i] )
            fixval = lbs[i];
         else if( ceil( lbs[i] ) == ubs[i] )
            fixval = ubs[i];
         else // otherwise take the midpoint
            fixval = REAL{ 0.5 } * ( ubs[i] + lbs[i] );

         TransactionGuard<REAL> tg{ reductions };
         reductions.lockColBounds( i );
         reductions.fixCol( i, fixval );
         result = PresolveStatus::kReduced;
      }
   }

   return result;
}

} // namespace papilo

#endif
