/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*               This file is part of the program and library                */
/*    PaPILO --- Parallel Presolve for Integer and Linear Optimization       */
/*                                                                           */
/* Copyright (C) 2020-2022 Konrad-Zuse-Zentrum                               */
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

#ifndef _FIX_VECTOR_MULTIPLICATION_HPP_
#define _FIX_VECTOR_MULTIPLICATION_HPP_

#include "fix/Constraint.hpp"
#include "papilo/core/Objective.hpp"
#include "papilo/core/Presolve.hpp"

#include <cassert>
#include <fstream>

using namespace papilo;

template <typename REAL>
class VectorMultiplication
{

 public:
   VectorMultiplication() = default;

   void
   calc_b_minus_Ax( const ConstraintMatrix<REAL>& A, const Vec<REAL>& x,
                    const Vec<REAL>& b,
                    const Vec<Constraint<REAL>>& constraints,
                    Vec<REAL>& result1, Vec<REAL>& result2 )
   {
      int n_rows_A = A.getNRows();
      assert( n_rows_A == b.size() );
      assert( A.getNCols() == x.size() );
      assert( constraints.size() == result2.size() );
      // TODO: add another assertion for x.size()=n_cols_constraints

      const Vec<RowFlags>& rowFlags = A.getRowFlags();
#ifdef PAPILO_TBB
      tbb::parallel_for( tbb::blocked_range<int>( 0, n_rows_A +
                                                     constraints.size() ),
                         [&]( const tbb::blocked_range<int>& r )
                         {
                            for( int i = r.begin(); i < r.end(); ++i )
#else
      for( int i = 0; i < n_rows_A + constraints.size(); ++i )
#endif
                            if( i < n_rows_A )
                            {
                               if( !rowFlags[i].test( RowFlag::kHardConstraint ) )
                               {
                                  auto coeff = A.getRowCoefficients( i );
                                  StableSum<REAL> aux( b[i] );
                                  for( int j = 0; j < coeff.getLength(); j++ )
                                     aux.add( -coeff.getValues()[j] *
                                              x[coeff.getIndices()[j]] );
                                  result1[i] = aux.get();
                               }
                               else
                                  result1[i] = 0;
                            }
                            else
                            {
                               auto coeff = constraints[i - n_rows_A].get_data();
                               StableSum<REAL> aux( constraints[i - n_rows_A].
                                                    get_lhs() );
                               for( int j = 0; j < coeff.getLength(); j++ )
                                  aux.add( -coeff.getValues()[j] *
                                           x[coeff.getIndices()[j]] );
                               result2[i - n_rows_A] = aux.get();
                            }
#ifdef PAPILO_TBB
                         } );
#endif
   }

   void
   calc_b_minus_xA( const ConstraintMatrix<REAL>& A, const Vec<REAL>& x,
                    const Vec<REAL>& b, Vec<REAL>& result )
   {
      assert( A.getNCols() == b.size() );
      assert( A.getNRows() == x.size() );

      const Vec<RowFlags>& rowFlags = A.getRowFlags();
#ifdef PAPILO_TBB
      tbb::parallel_for( tbb::blocked_range<int>( 0, A.getNCols() ),
                         [&]( const tbb::blocked_range<int>& r )
                         {
                            for( int i = r.begin(); i < r.end(); ++i )
#else
      for( int i = 0; i < A.getNCols(); ++i )
#endif
                            {
                               auto coeff = A.getColumnCoefficients( i );
                               StableSum<REAL> aux( b[i] );
                               for( int j = 0; j < coeff.getLength(); j++ )
                               {
                                  if( !rowFlags[coeff.getIndices()[j]].
                                       test( RowFlag::kHardConstraint ) )
                                     aux.add( -coeff.getValues()[j] *
                                              x[coeff.getIndices()[j]] );
                               }
                               result[i] = aux.get();
                            }
#ifdef PAPILO_TBB
                         } );
#endif
   }

   void
   calc_b_minus_xA( const Vec<Constraint<REAL>>& constraints,
                    const Vec<REAL>& x, const Vec<REAL>& b, Vec<REAL>& result )
   {
      // TODO: add another assertion for result.size()=n_cols_A=b.size()
      assert( constraints.size() == x.size() );

      result = b;
      for( int i = 0; i < constraints.size(); ++i )
      {
         auto coeff = constraints[i].get_data();
         for( int j = 0; j < coeff.getLength(); j++ )
            result[coeff.getIndices()[j]] -= coeff.getValues()[j] * x[i];
      }
   }

   void
   calc_qb_plus_sx( const REAL q, const Vec<REAL>& b, const REAL s,
                    const Vec<REAL>& x, Vec<REAL>& result )
   {
      assert( b.size() == x.size() );
#ifdef PAPILO_TBB
      tbb::parallel_for( tbb::blocked_range<int>( 0, b.size() ),
                         [&]( const tbb::blocked_range<int>& r )
                         {
                            for( int i = r.begin(); i < r.end(); ++i )
#else
      for( int i = 0; i < b.size(); ++i )
#endif
                            {
                               result[i] = q * b[i] + s * x[i];
                            }
#ifdef PAPILO_TBB
                         } );
#endif
   }

   void
   calc_b_plus_sx( const Vec<REAL>& b, const REAL s, const Vec<REAL>& x,
                   Vec<REAL>& result )
   {
      assert( b.size() == x.size() );
#ifdef PAPILO_TBB
      tbb::parallel_for( tbb::blocked_range<int>( 0, b.size() ),
                         [&]( const tbb::blocked_range<int>& r )
                         {
                            for( int i = r.begin(); i < r.end(); ++i )
#else
      for( int i = 0; i < b.size(); ++i )
#endif
                            {
                               result[i] = b[i] + s * x[i];
                            }
#ifdef PAPILO_TBB
                         } );
#endif
   }

   REAL
   multi( const Vec<REAL>& b, const Vec<REAL>& x )
   {
      assert( b.size() == x.size() );
      StableSum<REAL> result;

      for( int i = 0; i < b.size(); ++i )
         result.add( b[i] * x[i] );
      return result.get();
   }

   REAL
   l2_norm( const Vec<REAL>& vec )
   {
      StableSum<REAL> squared_norm;
      for( int i = 0; i < vec.size(); ++i )
         squared_norm.add( vec[i] * vec[i] );
      return sqrt( squared_norm.get() );
   }

   REAL
   l1_norm( const Vec<REAL>& vec )
   {
      StableSum<REAL> result;
      for( int i = 0; i < vec.size(); ++i )
         result.add( abs( vec[i] ) );
      return result.get();
   }


};

#endif
