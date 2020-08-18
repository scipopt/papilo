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

#ifndef _PAPILO_PARALLEL_COL_DETECTION_HPP_
#define _PAPILO_PARALLEL_COL_DETECTION_HPP_

#include "papilo/core/PresolveMethod.hpp"
#include "papilo/core/Problem.hpp"
#include "papilo/core/ProblemUpdate.hpp"
#include "papilo/misc/Hash.hpp"
#include "papilo/misc/tbb.hpp"
#include "pdqsort/pdqsort.h"

namespace papilo
{

template <typename REAL>
class ParallelColDetection : public PresolveMethod<REAL>
{
   struct SupportHashCompare
   {
      SupportHashCompare() {}

      static size_t
      hash( const std::pair<int, const int*>& row )
      {
         Hasher<size_t> hasher( row.first );

         const int* support = row.second;

         for( int i = 0; i != row.first; ++i )
         {
            hasher.addValue( support[i] );
         }

         return hasher.getHash();
      }

      static bool
      equal( const std::pair<int, const int*>& row1,
             const std::pair<int, const int*>& row2 )
      {
         int length = row1.first;

         if( length != row2.first )
            return false;

         return memcmp( static_cast<const void*>( row1.second ),
                        static_cast<const void*>( row2.second ),
                        length * sizeof( int ) ) == 0;
      }
   };

   struct SupportHash
   {
      std::size_t
      operator()( const std::pair<int, const int*>& row ) const
      {
         return SupportHashCompare::hash( row );
      }
   };

   struct SupportEqual
   {
      bool
      operator()( const std::pair<int, const int*>& row1,
                  const std::pair<int, const int*>& row2 ) const
      {
         return SupportHashCompare::equal( row1, row2 );
      }
   };

   void
   findParallelCols( const Num<REAL>& num, const int* bucket, int bucketsize,
                     const ConstraintMatrix<REAL>& constMatrix,
                     const Vec<REAL>& obj, const VariableDomains<REAL>& domains,
                     Vec<std::pair<int, int>>& parallelCols );

   void
   computeColHashes( const ConstraintMatrix<REAL>& constMatrix,
                     const Vec<REAL>& obj, unsigned int* colhashes );

   void
   computeSupportId( const ConstraintMatrix<REAL>& constMatrix,
                     unsigned int* supporthashes );

 public:
   ParallelColDetection() : PresolveMethod<REAL>()
   {
      this->setName( "parallelcols" );
      this->setTiming( PresolverTiming::kMedium );
   }

   bool
   initialize( const Problem<REAL>& problem,
               const PresolveOptions& presolveOptions ) override
   {
      if( presolveOptions.dualreds < 2 )
         this->setEnabled( false );
      return false;
   }

   /// todo how to communicate about postsolve information
   virtual PresolveStatus
   execute( const Problem<REAL>& problem,
            const ProblemUpdate<REAL>& problemUpdate, const Num<REAL>& num,
            Reductions<REAL>& reductions ) override;
};

#ifdef PAPILO_USE_EXTERN_TEMPLATES
extern template class ParallelColDetection<double>;
extern template class ParallelColDetection<Quad>;
extern template class ParallelColDetection<Rational>;
#endif

template <typename REAL>
void
ParallelColDetection<REAL>::findParallelCols(
    const Num<REAL>& num, const int* bucket, int bucketsize,
    const ConstraintMatrix<REAL>& constMatrix, const Vec<REAL>& obj,
    const VariableDomains<REAL>& domains,
    Vec<std::pair<int, int>>& parallelCols )
{
   // TODO if bucketsize too large do gurobi trick
   const Vec<ColFlags>& cflags = domains.flags;
   const Vec<REAL>& lbs = domains.lower_bounds;
   const Vec<REAL>& ubs = domains.upper_bounds;

   auto checkDomainsForHoles = [&]( int col1, int col2, const REAL& scale2 ) {
      // test whether we need to check that the domain of the merged
      // column has no holes
      assert( cflags[col1].test( ColFlag::kIntegral ) );
      assert( cflags[col2].test( ColFlag::kIntegral ) );
      assert( !cflags[col1].test( ColFlag::kLbInf, ColFlag::kUbInf ) );
      assert( !cflags[col2].test( ColFlag::kLbInf, ColFlag::kUbInf ) );

      // compute the domains of the merged column
      REAL mergeval = lbs[col2];
      REAL mergeub = ubs[col2];

      if( scale2 < 0 )
      {
         mergeval += scale2 * ubs[col1];
         mergeub += scale2 * lbs[col1];
      }
      else
      {
         mergeval += scale2 * lbs[col1];
         mergeub += scale2 * ubs[col1];
      }

      // scan domain of new variable for holes:
      // test every value in domain of column 1 if it allows to
      // find a value in the domain of column 2 to constitute the
      // value the all values in the domain of the merged column.
      // If no such value is found the columns cannot be merged
      bool foundhole = false;
      while( num.isLE( mergeval, mergeub ) )
      {
         // initialize col1val with the lower bound of column 1
         REAL col1val = lbs[col1];

         // test if col1val + scale2 * col2val = mergeval implies a
         // value for col2val that is within the domain of column 2.
         // If that is the case we can stop and increase mergeval by
         // 1, otherwise we found a hole If that check failed for
         // the current col1val we increase it by 1 the upper bound
         // of column 1 is reached
         foundhole = true;
         while( num.isLE( col1val, ubs[col1] ) )
         {
            REAL col2val = mergeval - col1val * scale2;

            if( num.isIntegral( col2val ) && num.isGE( col2val, lbs[col2] ) &&
                num.isLE( col2val, ubs[col2] ) )
            {
               foundhole = false;
               break;
            }

            col1val += 1;
         }

         // if a hole was found we can stop
         if( foundhole )
            break;

         // test next value in domain
         mergeval += 1;
      }

      return foundhole;
   };

   for( int i = 0; i < bucketsize; ++i )
   {
      int col1 = bucket[i];
      auto col1vec = constMatrix.getColumnCoefficients( col1 );

      const int length = col1vec.getLength();
      const REAL* coefs1 = col1vec.getValues();
      bool col1integral = cflags[col1].test( ColFlag::kIntegral );

      if( length < 2 )
         return;

      for( int j = i + 1; j < bucketsize; ++j )
      {
         int col2 = bucket[j];
         auto col2vec = constMatrix.getColumnCoefficients( col2 );

         if( ( obj[col1] == 0 ) != ( obj[col2] == 0 ) )
            continue;

         // support should already be checked
         assert( length == col2vec.getLength() );
         assert( std::memcmp( static_cast<const void*>( col1vec.getIndices() ),
                              static_cast<const void*>( col2vec.getIndices() ),
                              length * sizeof( int ) ) == 0 );

         const REAL* coefs2 = col2vec.getValues();
         bool col2integral = cflags[col2].test( ColFlag::kIntegral );

         bool checkdomains = false;
         bool scalecol2;

         if( col1integral != col2integral )
         {
            scalecol2 = col2integral;
            if( scalecol2 )
            {
               assert( !col1integral );
               if( !cflags[col1].test( ColFlag::kLbInf, ColFlag::kUbInf ) &&
                   num.isLT(
                       abs( ( ubs[col1] - lbs[col1] ) * coefs1[0] / coefs2[0] ),
                       1 ) )
                  continue;
            }
            else
            {
               assert( !col2integral );

               if( !cflags[col2].test( ColFlag::kLbInf, ColFlag::kUbInf ) &&
                   num.isLT(
                       abs( ( ubs[col2] - lbs[col2] ) * coefs2[0] / coefs1[0] ),
                       1 ) )
                  continue;
            }
         }
         else
         {
            scalecol2 = num.isGE( abs( coefs1[0] ), abs( coefs2[0] ) );

            if( col1integral &&
                !num.isEq( abs( coefs1[0] ), abs( coefs2[0] ) ) )
            {
               assert( col2integral );

               if( cflags[col1].test( ColFlag::kLbInf, ColFlag::kUbInf ) ||
                   cflags[col2].test( ColFlag::kLbInf, ColFlag::kUbInf ) )
                  continue;

               if( scalecol2 && !num.isIntegral( coefs1[0] / coefs2[0] ) )
                  continue;

               if( !scalecol2 && !num.isIntegral( coefs2[0] / coefs1[0] ) )
                  continue;

               checkdomains = true;
            }
         }

         bool parallel = true;

         if( scalecol2 )
         {
            REAL scale2 = coefs1[0] / coefs2[0];

            if( !num.isEq( obj[col1], scale2 * obj[col2] ) )
               continue;

            for( int k = 1; k < length; ++k )
            {
               if( !num.isEq( coefs1[k], scale2 * coefs2[k] ) )
               {
                  parallel = false;
                  break;
               }
            }

            if( parallel )
            {
               if( checkdomains && checkDomainsForHoles( col1, col2, scale2 ) )
                  continue;

               parallelCols.push_back( std::make_pair( col1, col2 ) );
            }
         }
         else
         {
            REAL scale1 = coefs2[0] / coefs1[0];

            if( !num.isEq( scale1 * obj[col1], obj[col2] ) )
               continue;

            for( int k = 1; k < length; ++k )
            {
               if( !num.isEq( scale1 * coefs1[k], coefs2[k] ) )
               {
                  parallel = false;
                  break;
               }
            }

            if( parallel )
            {
               if( checkdomains && checkDomainsForHoles( col2, col1, scale1 ) )
                  continue;

               parallelCols.push_back( std::make_pair( col2, col1 ) );
            }
         }
      }
   }
}

template <typename REAL>
void
ParallelColDetection<REAL>::computeColHashes(
    const ConstraintMatrix<REAL>& constMatrix, const Vec<REAL>& obj,
    unsigned int* colhashes )
{
   tbb::parallel_for(
       tbb::blocked_range<int>( 0, constMatrix.getNCols() ),
       [&]( const tbb::blocked_range<int>& r ) {
          for( int i = r.begin(); i != r.end(); ++i )
          {
             // compute hash-value for coefficients

             auto colcoefs = constMatrix.getColumnCoefficients( i );
             const REAL* vals = colcoefs.getValues();
             const int len = colcoefs.getLength();

             Hasher<unsigned int> hasher( len );

             if( len > 1 )
             {
                // compute scale such that the first coefficient is
                // positive 1/golden ratio. The choice of
                // the constant is arbitrary and is used to make cases
                // where two coefficients that are equal
                // within epsilon get different values are
                // more unlikely by choosign some irrational number
                REAL scale = REAL( 2.0 / ( 1.0 + sqrt( 5.0 ) ) ) / vals[0];

                // add scaled coefficients of other row
                // entries to compute the hash
                for( int j = 1; j != len; ++j )
                   hasher.addValue( Num<REAL>::hashCode( vals[j] * scale ) );

                if( obj[i] != 0 )
                   hasher.addValue( Num<REAL>::hashCode( obj[i] * scale ) );
             }

             colhashes[i] = hasher.getHash();
          }
       } );
}

template <typename REAL>
void
ParallelColDetection<REAL>::computeSupportId(
    const ConstraintMatrix<REAL>& constMatrix, unsigned int* supporthashes )
{
   using SupportMap =
       HashMap<std::pair<int, const int*>, int, SupportHash, SupportEqual>;

   SupportMap supportMap(
       static_cast<std::size_t>( constMatrix.getNCols() * 1.1 ) );

   for( int i = 0; i < constMatrix.getNCols(); ++i )
   {
      auto col = constMatrix.getColumnCoefficients( i );
      int length = col.getLength();
      const int* support = col.getIndices();

      auto insResult =
          supportMap.emplace( std::make_pair( length, support ), i );

      if( insResult.second )
         supporthashes[i] = i;
      else // support already exists, use the previous support id
         supporthashes[i] = insResult.first->second;
   }
}

/// todo how to communicate about postsolve information
template <typename REAL>
PresolveStatus
ParallelColDetection<REAL>::execute( const Problem<REAL>& problem,
                                     const ProblemUpdate<REAL>& problemUpdate,
                                     const Num<REAL>& num,
                                     Reductions<REAL>& reductions )
{
   const auto& constMatrix = problem.getConstraintMatrix();
   const auto& obj = problem.getObjective().coefficients;
   const auto& lhs_values = constMatrix.getLeftHandSides();
   const auto& rhs_values = constMatrix.getRightHandSides();
   const auto& rflags = constMatrix.getRowFlags();
   const auto& cflags = problem.getColFlags();
   const int ncols = constMatrix.getNCols();
   const Vec<int>& colperm = problemUpdate.getRandomColPerm();

   PresolveStatus result = PresolveStatus::kUnchanged;

   // get called less and less over time regardless of success since the
   // presolver can be expensive otherwise
   this->skipRounds( this->getNCalls() );

   assert( ncols > 0 );

   std::unique_ptr<unsigned int[]> supportid{ new unsigned int[ncols] };
   std::unique_ptr<unsigned int[]> coefhash{ new unsigned int[ncols] };
   std::unique_ptr<int[]> col{ new int[ncols] };

   tbb::parallel_invoke(
       [ncols, &col]() {
          for( int i = 0; i < ncols; ++i )
             col[i] = i;
       },
       [&constMatrix, &coefhash, &obj, this]() {
          computeColHashes( constMatrix, obj, coefhash.get() );
       },
       [&constMatrix, &supportid, this]() {
          computeSupportId( constMatrix, supportid.get() );
       } );

   pdqsort( col.get(), col.get() + ncols, [&]( int a, int b ) {
      return supportid[a] < supportid[b] ||
             ( supportid[a] == supportid[b] && coefhash[a] < coefhash[b] ) ||
             ( supportid[a] == supportid[b] && coefhash[a] == coefhash[b] &&
               colperm[a] < colperm[b] );
   } );

   Vec<std::pair<int, int>> parallelCols;

   for( int i = 0; i < ncols; )
   {
      // determine size of bucket
      int j;
      for( j = i + 1; j < ncols; ++j )
      {
         if( coefhash[col[i]] != coefhash[col[j]] ||
             supportid[col[i]] != supportid[col[j]] )
            break;
      }
      int len = j - i;

      // if more  than one col is in the bucket try to find parallel
      // cols
      if( len > 1 )
      {
         // fmt::print( "bucket of length {} starting at {}\n", len, i
         // );
         findParallelCols( num, col.get() + i, len, constMatrix, obj,
                           problem.getVariableDomains(), parallelCols );
      }

      assert( j > i );
      // set i to start of next bucket
      i = j;
   }

   if( !parallelCols.empty() )
   {
      pdqsort( parallelCols.begin(), parallelCols.end(),
               [&colperm]( const std::pair<int, int>& a,
                           const std::pair<int, int>& b ) {
                  return std::make_pair( colperm[a.first], colperm[a.second] ) <
                         std::make_pair( colperm[b.first], colperm[b.second] );
               } );
      result = PresolveStatus::kReduced;
      // fmt::print( "found {} parallel columns\n", parallelCols.size() );
      for( const std::pair<int, int>& parallelCol : parallelCols )
      {
         int col1 = parallelCol.first;
         int col2 = parallelCol.second;

         TransactionGuard<REAL> tg{ reductions };
         reductions.lockCol( col1 );
         reductions.lockCol( col2 );

         if( cflags[col2].test( ColFlag::kIntegral ) !=
             cflags[col1].test( ColFlag::kIntegral ) )
            reductions.lockColBounds( col1 );

         reductions.parallelCols( col1, col2 );
      }
   }

   return result;
}

} // namespace papilo

#endif
