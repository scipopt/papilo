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

#ifndef _PAPILO_PRESOLVERS_SIMPLE_FREEVAR_HPP_
#define _PAPILO_PRESOLVERS_SIMPLE_FREEVAR_HPP_

#include "papilo/core/PresolveMethod.hpp"
#include "papilo/core/Problem.hpp"
#include "papilo/core/ProblemUpdate.hpp"
#include "papilo/misc/Num.hpp"
#include "papilo/misc/fmt.hpp"

namespace papilo
{

template <typename REAL>
class SimpleSubstitution : public PresolveMethod<REAL>
{
 public:
   SimpleSubstitution() : PresolveMethod<REAL>()
   {
      this->setName( "doubletoneq" );
      this->setTiming( PresolverTiming::kMedium );
   }

   virtual PresolveStatus
   execute( const Problem<REAL>& problem,
            const ProblemUpdate<REAL>& problemUpdate, const Num<REAL>& num,
            Reductions<REAL>& reductions ) override;
};

#ifdef PAPILO_USE_EXTERN_TEMPLATES
extern template class SimpleSubstitution<double>;
extern template class SimpleSubstitution<Quad>;
extern template class SimpleSubstitution<Rational>;
#endif

template <typename REAL>
PresolveStatus
SimpleSubstitution<REAL>::execute( const Problem<REAL>& problem,
                                   const ProblemUpdate<REAL>& problemUpdate,
                                   const Num<REAL>& num,
                                   Reductions<REAL>& reductions )
{
   // go over the rows and get the equalities, extract the columns that
   // verify the conditions add them to a hash map, loop over the hash map
   // and compute the implied bounds and finally look for implied free
   // variables and add reductions
   const auto& domains = problem.getVariableDomains();
   const auto& lower_bounds = domains.lower_bounds;
   const auto& upper_bounds = domains.upper_bounds;
   const auto& cflags = domains.flags;
   const auto& domainFlags = domains.flags;
   const auto& obj = problem.getObjective().coefficients;

   const auto& activities = problem.getRowActivities();

   const auto& constMatrix = problem.getConstraintMatrix();
   const auto& lhs_values = constMatrix.getLeftHandSides();
   const auto& rhs_values = constMatrix.getRightHandSides();
   const auto& rflags = constMatrix.getRowFlags();
   const auto& nrows = constMatrix.getNRows();
   const auto& ncols = constMatrix.getNCols();
   const auto& rowperm = problemUpdate.getRandomRowPerm();

   PresolveStatus result = PresolveStatus::kUnchanged;

   for( int k = 0; k < nrows; ++k )
   {
      int i = rowperm[k];
      // check that equality flag is correct or row is redundant
      assert( rflags[i].test( RowFlag::kRedundant ) ||
              ( !rflags[i].test( RowFlag::kEquation ) &&
                ( lhs_values[i] != rhs_values[i] ||
                  rflags[i].test( RowFlag::kLhsInf, RowFlag::kRhsInf ) ) ) ||
              ( rflags[i].test( RowFlag::kEquation ) &&
                !rflags[i].test( RowFlag::kLhsInf, RowFlag::kRhsInf ) &&
                lhs_values[i] == rhs_values[i] ) );

      if( rflags[i].test( RowFlag::kRedundant ) ||
          !rflags[i].test( RowFlag::kEquation ) ||
          constMatrix.getRowSizes()[i] != 2 )
         continue;

      auto rowvec = constMatrix.getRowCoefficients( i );
      assert( rowvec.getLength() == 2 );
      const REAL* vals = rowvec.getValues();
      const int* inds = rowvec.getIndices();
      REAL rhs = rhs_values[i];

      int subst;
      int stay;

      if( cflags[inds[0]].test( ColFlag::kIntegral ) !=
          cflags[inds[1]].test( ColFlag::kIntegral ) )
      {
         if( cflags[inds[0]].test( ColFlag::kIntegral ) )
            subst = 1;
         else
            subst = 0;

         stay = 1 - subst;
      }
      else if( cflags[inds[0]].test( ColFlag::kIntegral ) )
      {
         assert( cflags[inds[1]].test( ColFlag::kIntegral ) );
         if( abs( vals[0] ) < abs( vals[1] ) ||
             ( abs( vals[0] ) == abs( vals[1] ) &&
               problemUpdate.isColBetterForSubstitution( inds[0], inds[1] ) ) )
            subst = 0;
         else
            subst = 1;

         stay = 1 - subst;

         if( !num.isIntegral( vals[stay] / vals[subst] ) )
            continue;

         if( !num.isFeasIntegral( rhs / vals[subst] ) )
            return PresolveStatus::kInfeasible;
      }
      else
      {
         REAL absval0 = abs( vals[0] );
         REAL absval1 = abs( vals[1] );
         if( absval0 * problemUpdate.getPresolveOptions().markowitz_tolerance >
             absval1 )
            subst = 0;
         else if( absval1 *
                      problemUpdate.getPresolveOptions().markowitz_tolerance >
                  absval0 )
            subst = 1;
         else if( problemUpdate.isColBetterForSubstitution( inds[0], inds[1] ) )
            subst = 0;
         else
            subst = 1;

         stay = 1 - subst;
      }

      result = PresolveStatus::kReduced;

      TransactionGuard<REAL> guard{ reductions };

      reductions.lockRow( i );

      reductions.lockColBounds( inds[subst] );

      REAL s = vals[subst] * vals[stay];
      if( !cflags[inds[subst]].test( ColFlag::kLbInf ) )
      {
         REAL staybound =
             ( rhs - vals[subst] * domains.lower_bounds[inds[subst]] ) /
             vals[stay];
         if( s < 0 &&
             ( cflags[inds[stay]].test( ColFlag::kLbInf ) ||
               num.isGT( staybound, domains.lower_bounds[inds[stay]] ) ) )
            reductions.changeColLB( inds[stay], staybound );

         if( s > 0 &&
             ( cflags[inds[stay]].test( ColFlag::kUbInf ) ||
               num.isLT( staybound, domains.upper_bounds[inds[stay]] ) ) )
            reductions.changeColUB( inds[stay], staybound );
      }

      if( !cflags[inds[subst]].test( ColFlag::kUbInf ) )
      {
         REAL staybound =
             ( rhs - vals[subst] * domains.upper_bounds[inds[subst]] ) /
             vals[stay];
         if( s > 0 && ( cflags[inds[stay]].test( ColFlag::kLbInf ) ||
                        staybound > domains.lower_bounds[inds[stay]] ) )
            reductions.changeColLB( inds[stay], staybound );

         if( s < 0 && ( cflags[inds[stay]].test( ColFlag::kUbInf ) ||
                        staybound < domains.upper_bounds[inds[stay]] ) )
            reductions.changeColUB( inds[stay], staybound );
      }

      reductions.aggregateFreeCol( inds[subst], i );
   }

   return result;
}

} // namespace papilo

#endif
