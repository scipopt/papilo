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


   Vec<REAL>
   calc_b_minus_Ax( const ConstraintMatrix<REAL>& A, const Vec<REAL>& x,
                    const Vec<REAL>& b )
   {
      assert( A.getNRows() == b.size() );
      assert( A.getNCols() == x.size() );
      Vec<REAL> result( b );
#ifdef PAPILO_TBB
      tbb::parallel_for( tbb::blocked_range<int>( 0, A.getNRows() ),
                         [&]( const tbb::blocked_range<int>& r )
                         {
                            for( int i = r.begin(); i < r.end(); ++i )
#else
      for( int i = 0; i < A.getNRows(); ++i )
#endif
                            {
                               auto coeff = A.getRowCoefficients( i );
                               StableSum<REAL> aux( b[i] );
                               for( int j = 0; j < coeff.getLength(); j++ )
                                  aux.add( -coeff.getValues()[j] *
                                           x[coeff.getIndices()[j]] );
                               result[i] = aux.get();
                            }
#ifdef PAPILO_TBB
                         } );
#endif
      return result;
   }

   Vec<REAL>
   calc_b_minus_xA( const ConstraintMatrix<REAL>& A, const Vec<REAL>& x,
                    const Vec<REAL>& b )
   {
      assert( A.getNCols() == b.size() );
      assert( A.getNRows() == x.size() );
      Vec<REAL> result( b );
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
                                  aux.add( -coeff.getValues()[j] *
                                           x[coeff.getIndices()[j]] );
                               result[i] = aux.get();
                            }
#ifdef PAPILO_TBB
                         } );
#endif
      return result;
   }

   REAL
   l2_norm( const Vec<REAL>& vec )
   {
      StableSum<REAL> squared_norm;

//#ifdef PAPILO_TBB
//         tbb::parallel_for( tbb::blocked_range<int>( 0, vec.size() ),
//                            [&]( const tbb::blocked_range<int>& r )
//                            {
//                               for( int i = r.begin(); i < r.end(); ++i )
//#else
         for( int i = 0; i < vec.size(); ++i )
//#endif
                               {
                                  squared_norm.add( vec[i] * vec[i] );
                               }

//#ifdef PAPILO_TBB
//                            } );
//#endif
      return sqrt( squared_norm.get() );
   }

   REAL
   l1_norm( const Vec<REAL>& vec )
   {
      StableSum<REAL> result;

//#ifdef PAPILO_TBB
//         tbb::parallel_for( tbb::blocked_range<int>( 0, vec.size() ),
//                            [&]( const tbb::blocked_range<int>& r )
//                            {
//                               for( int i = r.begin(); i < r.end(); ++i )
//#else
         for( int i = 0; i < vec.size(); ++i )
//#endif
                               {
                                  result.add( abs(vec[i]) );
                               }

//#ifdef PAPILO_TBB
//                            } );
//#endif
      return result.get();
   }

   Vec<REAL>
   calc_qb_plus_sx( const REAL q, const Vec<REAL>& b, const REAL s, const Vec<REAL>& x )
   {
      assert( b.size() == x.size() );
      Vec<REAL> result( b );
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
      return result;
   }

   Vec<REAL>
   calc_b_plus_sx( const Vec<REAL>& b, const REAL s, const Vec<REAL>& x )
   {
      assert( b.size() == x.size() );
      Vec<REAL> result( b );
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
      return result;
   }

   REAL
   multi( const Vec<REAL>& b, const Vec<REAL>& x )
   {
      assert( b.size() == x.size() );
      StableSum<REAL> result;
#ifdef PAPILO_TBB
      tbb::parallel_for( tbb::blocked_range<int>( 0, b.size() ),
                         [&]( const tbb::blocked_range<int>& r )
                         {
                            for( int i = r.begin(); i < r.end(); ++i )
#else
      for( int i = 0; i < b.size(); ++i )
#endif
                            {
                               result.add( b[i] * x[i] );
                            }
#ifdef PAPILO_TBB
                         } );
#endif
      return result.get();
   }




};
