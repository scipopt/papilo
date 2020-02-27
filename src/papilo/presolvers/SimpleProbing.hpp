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

#ifndef _PAPILO_PRESOLVERS_SIMPLE_PROBING_HPP_
#define _PAPILO_PRESOLVERS_SIMPLE_PROBING_HPP_

#include "papilo/core/PresolveMethod.hpp"
#include "papilo/core/Problem.hpp"
#include "papilo/core/ProblemUpdate.hpp"
#include "papilo/io/Message.hpp"

namespace papilo
{

template <typename REAL>
class SimpleProbing : public PresolveMethod<REAL>
{
 public:
   SimpleProbing() : PresolveMethod<REAL>()
   {
      this->setName( "simpleprobing" );
      this->setType( PresolverType::kIntegralCols );
      this->setTiming( PresolverTiming::kMedium );
   }

   virtual PresolveStatus
   execute( const Problem<REAL>& problem,
            const ProblemUpdate<REAL>& problemUpdate, const Num<REAL>& num,
            Reductions<REAL>& reductions ) override;
};

#ifdef PAPILO_USE_EXTERN_TEMPLATES
extern template class SimpleProbing<double>;
extern template class SimpleProbing<Quad>;
extern template class SimpleProbing<Rational>;
#endif

template <typename REAL>
PresolveStatus
SimpleProbing<REAL>::execute( const Problem<REAL>& problem,
                              const ProblemUpdate<REAL>& problemUpdate,
                              const Num<REAL>& num,
                              Reductions<REAL>& reductions )
{
   assert( problem.getNumIntegralCols() != 0 );

   const auto& domains = problem.getVariableDomains();
   const auto& cflags = domains.flags;
   const auto& activities = problem.getRowActivities();
   const auto& changedactivities = problemUpdate.getChangedActivities();

   const auto& rowsize = problem.getRowSizes();

   const auto& constMatrix = problem.getConstraintMatrix();
   const auto& lhs_values = constMatrix.getLeftHandSides();
   const auto& rhs_values = constMatrix.getRightHandSides();
   const auto& rflags = constMatrix.getRowFlags();

   PresolveStatus result = PresolveStatus::kUnchanged;
   int nrows = problem.getNRows();

   for( int i = 0; i != nrows; ++i )
   {
      if( !rflags[i].test( RowFlag::kEquation ) || rowsize[i] <= 2 ||
          activities[i].ninfmin != 0 || activities[i].ninfmax != 0 ||
          !num.isEq( activities[i].min + activities[i].max,
                     2 * rhs_values[i] ) )
         continue;

      assert( rflags[i].test( RowFlag::kEquation ) );
      assert( activities[i].ninfmin == 0 && activities[i].ninfmax == 0 );
      assert( num.isEq( activities[i].min + activities[i].max,
                        2 * rhs_values[i] ) );

      auto rowvec = constMatrix.getRowCoefficients( i );
      const REAL* rowvals = rowvec.getValues();
      const int* rowcols = rowvec.getIndices();
      const int rowlen = rowvec.getLength();

      REAL bincoef = activities[i].max - rhs_values[i];
      int bincol = -1;

      for( int k = 0; k != rowlen; ++k )
      {
         int col = rowcols[k];
         assert( !cflags[col].test( ColFlag::kUnbounded ) );
         if( !cflags[col].test( ColFlag::kIntegral ) ||
             domains.lower_bounds[col] != 0 || domains.upper_bounds[col] != 1 ||
             !num.isEq( abs( rowvals[k] ), bincoef ) )
            continue;

         bincol = col;
         // could be negative
         bincoef = rowvals[k];
         break;
      }

      if( bincol != -1 )
      {
         assert(
             num.isEq( abs( bincoef ), activities[i].max - rhs_values[i] ) );
         assert( domains.lower_bounds[bincol] == 0 );
         assert( domains.upper_bounds[bincol] == 1 );
         assert( cflags[bincol].test( ColFlag::kIntegral ) );

         Message::debug(
             this, "probing on simple equation detected {} subsitutions\n",
             rowlen - 1 );

         result = PresolveStatus::kReduced;
         for( int k = 0; k != rowlen; ++k )
         {
            int col = rowcols[k];
            if( col == bincol )
               continue;

            REAL factor;
            REAL offset;
            if( bincoef * rowvals[k] > 0 )
            {
               factor =
                   -( domains.upper_bounds[col] - domains.lower_bounds[col] );
               offset = domains.upper_bounds[col];
            }
            else
            {
               factor = domains.upper_bounds[col] - domains.lower_bounds[col];
               offset = domains.lower_bounds[col];
            }

            reductions.replaceCol( col, bincol, factor, offset );
         }
      }
   }
   return result;
}

} // namespace papilo

#endif
