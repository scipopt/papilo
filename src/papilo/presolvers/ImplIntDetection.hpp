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

#ifndef _PAPILO_PRESOLVERS_IMPL_INT_DETECTION_HPP_
#define _PAPILO_PRESOLVERS_IMPL_INT_DETECTION_HPP_

#include "papilo/core/PresolveMethod.hpp"
#include "papilo/core/Problem.hpp"
#include "papilo/core/ProblemUpdate.hpp"
#include "papilo/misc/Num.hpp"
#include "papilo/misc/fmt.hpp"

namespace papilo
{

template <typename REAL>
class ImplIntDetection : public PresolveMethod<REAL>
{
 public:
   ImplIntDetection() : PresolveMethod<REAL>()
   {
      this->setName( "implint" );
      this->setTiming( PresolverTiming::EXHAUSTIVE );
      this->setType( PresolverType::MIXED_COLS );
   }

   virtual PresolveStatus
   execute( const Problem<REAL>& problem,
            const ProblemUpdate<REAL>& problemUpdate, const Num<REAL>& num,
            Reductions<REAL>& reductions ) override;
};

#ifdef PAPILO_USE_EXTERN_TEMPLATES
extern template class ImplIntDetection<double>;
extern template class ImplIntDetection<Quad>;
extern template class ImplIntDetection<Rational>;
#endif

template <typename REAL>
PresolveStatus
ImplIntDetection<REAL>::execute( const Problem<REAL>& problem,
                                 const ProblemUpdate<REAL>& problemUpdate,
                                 const Num<REAL>& num,
                                 Reductions<REAL>& reductions )
{
   const auto& domains = problem.getVariableDomains();
   const auto& lower_bounds = domains.lower_bounds;
   const auto& upper_bounds = domains.upper_bounds;
   const auto& cflags = domains.flags;
   const auto& domainFlags = domains.flags;

   const auto& activities = problem.getRowActivities();

   const auto& consmatrix = problem.getConstraintMatrix();
   const auto& lhs_values = consmatrix.getLeftHandSides();
   const auto& rhs_values = consmatrix.getRightHandSides();
   const auto& rflags = consmatrix.getRowFlags();
   const auto& nrows = consmatrix.getNRows();
   const auto& ncols = consmatrix.getNCols();

   PresolveStatus result = PresolveStatus::UNCHANGED;

   for( int col = 0; col != ncols; ++col )
   {
      if( cflags[col].test( ColFlag::INTEGRAL, ColFlag::IMPL_INT,
                            ColFlag::INACTIVE ) )
         continue;

      bool testinequalities =
          problemUpdate.getPresolveOptions().dualreds == 2 ? true : false;
      bool impliedint = false;

      auto colvec = consmatrix.getColumnCoefficients( col );
      int collen = colvec.getLength();
      const int* colrows = colvec.getIndices();
      const REAL* colvals = colvec.getValues();

      for( int i = 0; i != collen; ++i )
      {
         int row = colrows[i];

         if( rflags[row].test( RowFlag::REDUNDANT ) ||
             !rflags[row].test( RowFlag::EQUALITY ) )
            continue;

         testinequalities = false;
         REAL scale = 1 / colvals[i];
         if( !num.isIntegral( scale * rhs_values[row] ) )
            continue;

         auto rowvec = consmatrix.getRowCoefficients( row );
         int rowlen = rowvec.getLength();
         const int* rowcols = rowvec.getIndices();
         const REAL* rowvals = rowvec.getValues();

         impliedint = true;

         for( int j = 0; j != rowlen; ++j )
         {
            int rowcol = rowcols[j];

            if( rowcol == col )
               continue;

            if( !cflags[rowcol].test( ColFlag::INTEGRAL, ColFlag::IMPL_INT ) ||
                !num.isIntegral( scale * rowvals[j] ) )
            {
               impliedint = false;
               break;
            }
         }

         if( impliedint )
            break;
      }

      if( impliedint )
      {
         // add reduction
         reductions.impliedInteger( col );
         result = PresolveStatus::REDUCED;
         continue;
      }

      if( !testinequalities )
         continue;

      if( !cflags[col].test( ColFlag::LB_INF ) &&
          !num.isIntegral( lower_bounds[col] ) )
         continue;

      if( !cflags[col].test( ColFlag::UB_INF ) &&
          !num.isIntegral( upper_bounds[col] ) )
         continue;

      impliedint = true;

      for( int i = 0; i != collen; ++i )
      {
         int row = colrows[i];

         if( rflags[row].test( RowFlag::REDUNDANT ) )
            continue;

         REAL scale = 1 / colvals[i];

         if( !rflags[row].test( RowFlag::RHS_INF ) &&
             !num.isIntegral( scale * rhs_values[row] ) )
         {
            impliedint = false;
            break;
         }

         if( !rflags[row].test( RowFlag::LHS_INF ) &&
             !num.isIntegral( scale * lhs_values[row] ) )
         {
            impliedint = false;
            break;
         }

         auto rowvec = consmatrix.getRowCoefficients( row );
         int rowlen = rowvec.getLength();
         const int* rowcols = rowvec.getIndices();
         const REAL* rowvals = rowvec.getValues();

         for( int j = 0; j != rowlen; ++j )
         {
            int rowcol = rowcols[j];

            if( rowcol == col )
               continue;

            if( !cflags[rowcol].test( ColFlag::INTEGRAL, ColFlag::IMPL_INT ) ||
                !num.isIntegral( scale * rowvals[j] ) )
            {
               impliedint = false;
               break;
            }
         }

         if( !impliedint )
            break;
      }

      if( impliedint )
      {
         // add reduction
         reductions.impliedInteger( col );
         result = PresolveStatus::REDUCED;
      }
   }

   return result;
}

} // namespace papilo

#endif
