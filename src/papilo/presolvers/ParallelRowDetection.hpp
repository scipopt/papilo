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

#ifndef _PAPILO_PARALLEL_ROW_DETECTION_HPP_
#define _PAPILO_PARALLEL_ROW_DETECTION_HPP_

#include "papilo/core/PresolveMethod.hpp"
#include "papilo/core/Problem.hpp"
#include "papilo/core/ProblemUpdate.hpp"
#include "papilo/misc/Hash.hpp"
#include "papilo/misc/tbb.hpp"
#include "pdqsort/pdqsort.h"

namespace papilo
{

template <typename REAL>
class ParallelRowDetection : public PresolveMethod<REAL>
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
   findParallelRows( const Num<REAL>& num, const int* bucket, int bucketsize,
                     const ConstraintMatrix<REAL>& constMatrix,
                     Vec<std::tuple<int, int, REAL>>& parallelRows );

   void
   computeRowHashes( const ConstraintMatrix<REAL>& constMatrix,
                     unsigned int* rowhashes );

   void
   computeSupportId( const ConstraintMatrix<REAL>& constMatrix,
                     unsigned int* supporthashes );

   void
   computeSupportIdParallel( const ConstraintMatrix<REAL>& constMatrix,
                             unsigned int* supportid );

 public:
   ParallelRowDetection() : PresolveMethod<REAL>()
   {
      this->setName( "parallelrows" );
      this->setTiming( PresolverTiming::kMedium );
   }

   /// todo how to communicate about postsolve information
   virtual PresolveStatus
   execute( const Problem<REAL>& problem,
            const ProblemUpdate<REAL>& problemUpdate, const Num<REAL>& num,
            Reductions<REAL>& reductions ) override;
};

#ifdef PAPILO_USE_EXTERN_TEMPLATES
extern template class ParallelRowDetection<double>;
extern template class ParallelRowDetection<Quad>;
extern template class ParallelRowDetection<Rational>;
#endif

template <typename REAL>
void
ParallelRowDetection<REAL>::findParallelRows(
    const Num<REAL>& num, const int* bucket, int bucketsize,
    const ConstraintMatrix<REAL>& constMatrix,
    Vec<std::tuple<int, int, REAL>>& parallelRows )
{
   // TODO if bucketsize too large do gurobi trick

   for( int i = 0; i < bucketsize; ++i )
   {
      auto row1 = constMatrix.getRowCoefficients( bucket[i] );

      const int length = row1.getLength();
      const REAL* coefs1 = row1.getValues();

      if( length < 2 )
         return;

      // TODO handle case of multiple parallel rows with one transaction:

      for( int j = i + 1; j < bucketsize; ++j )
      {
         auto row2 = constMatrix.getRowCoefficients( bucket[j] );

         // support should already be checked
         assert( length == row2.getLength() );
         assert( std::memcmp( static_cast<const void*>( row1.getIndices() ),
                              static_cast<const void*>( row2.getIndices() ),
                              length * sizeof( int ) ) == 0 );

         const REAL* coefs2 = row2.getValues();

         bool parallel = true;

         if( num.isGE( abs( coefs1[0] ), abs( coefs2[0] ) ) )
         {
            REAL scale2 = coefs1[0] / coefs2[0];

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
               parallelRows.push_back(
                   std::make_tuple( bucket[i], bucket[j], scale2 ) );
            }
         }
         else
         {
            REAL scale1 = coefs2[0] / coefs1[0];

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
               parallelRows.push_back(
                   std::make_tuple( bucket[j], bucket[i], scale1 ) );
            }
         }
      }
   }
}

template <typename REAL>
void
ParallelRowDetection<REAL>::computeRowHashes(
    const ConstraintMatrix<REAL>& constMatrix, unsigned int* rowhashes )
{
   tbb::parallel_for(
       tbb::blocked_range<int>( 0, constMatrix.getNRows() ),
       [&]( const tbb::blocked_range<int>& r ) {
          for( int i = r.begin(); i != r.end(); ++i )
          {
             // compute hash-value for coefficients

             auto rowcoefs = constMatrix.getRowCoefficients( i );
             const REAL* rowvals = rowcoefs.getValues();
             const int len = rowcoefs.getLength();

             Hasher<unsigned int> hasher( len );
             // only makes sense for non-singleton rows
             // (should not occur after redundant rows are
             // already deleted)
             if( len > 1 )
             {
                // compute scale such that the first coefficient is
                // positive 1/golden ratio. The choice of
                // the constant is arbitrary and is used to make cases
                // where two coefficients that are equal
                // within epsilon get different values are
                // more unlikely by choosing some irrational number
                REAL scale = REAL( 2.0 / ( 1.0 + sqrt( 5.0 ) ) ) / rowvals[0];

                // add scaled coefficients of other row
                // entries to compute the hash
                for( int j = 1; j != len; ++j )
                {
                   hasher.addValue( Num<REAL>::hashCode( rowvals[j] * scale ) );
                }
             }

             rowhashes[i] = hasher.getHash();
          }
       } );
}

template <typename REAL>
void
ParallelRowDetection<REAL>::computeSupportId(
    const ConstraintMatrix<REAL>& constMatrix, unsigned int* supporthashes )
{
   using SupportMap =
       HashMap<std::pair<int, const int*>, int, SupportHash, SupportEqual>;

   SupportMap supportMap(
       static_cast<std::size_t>( constMatrix.getNRows() * 1.1 ) );

   for( int i = 0; i < constMatrix.getNRows(); ++i )
   {
      auto row = constMatrix.getRowCoefficients( i );
      int length = row.getLength();
      const int* support = row.getIndices();

      auto insResult =
          supportMap.emplace( std::make_pair( length, support ), i );

      if( insResult.second )
         supporthashes[i] = i;
      else // support already exists, use the previous support id
         supporthashes[i] = insResult.first->second;
   }
}

template <typename REAL>
void
ParallelRowDetection<REAL>::computeSupportIdParallel(
    const ConstraintMatrix<REAL>& constMatrix, unsigned int* supportid )
{
   using SupportMap =
       tbb::concurrent_hash_map<std::pair<int, const int*>, unsigned int,
                                SupportHashCompare>;

   SupportMap supportMap( constMatrix.getNRows() * 2 );

   tbb::parallel_for(
       tbb::blocked_range<int>( 0, constMatrix.getNRows() ),
       [&]( const tbb::blocked_range<int>& r ) {
          for( int i = r.begin(); i != r.end(); ++i )
          {
             unsigned int thissupportid;
             auto row = constMatrix.getRowCoefficients( i );
             int length = row.getLength();
             const int* support = row.getIndices();

             {
                typename SupportMap::const_accessor a;
                if( supportMap.insert(
                        a, std::make_pair( std::make_pair( length, support ),
                                           i ) ) )
                {
                   thissupportid = i;
                }
                else
                   thissupportid = a->second;
             }

             supportid[i] = thissupportid;
          }
       } );
}

template <typename REAL>
PresolveStatus
ParallelRowDetection<REAL>::execute( const Problem<REAL>& problem,
                                     const ProblemUpdate<REAL>& problemUpdate,
                                     const Num<REAL>& num,
                                     Reductions<REAL>& reductions )
{
   const auto& constMatrix = problem.getConstraintMatrix();
   const auto& lhs_values = constMatrix.getLeftHandSides();
   const auto& rhs_values = constMatrix.getRightHandSides();
   const auto& rflags = constMatrix.getRowFlags();
   const int nrows = constMatrix.getNRows();
   const Vec<int>& rowperm = problemUpdate.getRandomRowPerm();

   PresolveStatus result = PresolveStatus::kUnchanged;

   // get called less and less over time regardless of success since the
   // presolver can be expensive otherwise
   this->skipRounds( this->getNCalls() );

   // lambda to handle the parallel rows
   auto handlerows = [&reductions, &result, &lhs_values, &rhs_values, &rflags,
                      &num]( int row1, int row2, REAL ratio ) {
      bool firstconsEquality = rflags[row1].test( RowFlag::kEquation );
      bool secondconsEquality = rflags[row2].test( RowFlag::kEquation );

      assert( !firstconsEquality || ( !rflags[row1].test( RowFlag::kRhsInf ) &&
                                      !rflags[row1].test( RowFlag::kLhsInf ) &&
                                      rhs_values[row1] == lhs_values[row1] ) );
      assert( !secondconsEquality || ( !rflags[row2].test( RowFlag::kRhsInf ) &&
                                       !rflags[row2].test( RowFlag::kLhsInf ) &&
                                       rhs_values[row2] == lhs_values[row2] ) );
      assert( ratio != 0.0 );

      REAL adjustedLHS = lhs_values[row2] * ratio;
      REAL adjustedRHS = rhs_values[row2] * ratio;
      bool adjustedLHSInf = rflags[row2].test( RowFlag::kLhsInf );
      bool adjustedRHSInf = rflags[row2].test( RowFlag::kRhsInf );

      using std::swap;
      if( ratio < REAL{ 0.0 } )
      {
         swap( adjustedLHS, adjustedRHS );
         swap( adjustedLHSInf, adjustedRHSInf );
      }

      // l1 <= A1 <= r1, first row
      // adjusted second row: s = A1 / A2
      // s*l2 <= s*A2 <= s*r2, if s > 0
      // s*r2 <= s*A2 <= s*l2, if s < 0
      if( firstconsEquality && secondconsEquality )
      {
         // l1 = r1, l2 = r2
         // if r1 != r2 infeasible
         if( num.isFeasEq( rhs_values[row1], adjustedRHS ) )
         {
            TransactionGuard<REAL> guard{ reductions };
            reductions.lockRow( row1 );
            reductions.lockRow( row2 );
            reductions.markRowRedundant( row2 );
         }
         else
         {
            result = PresolveStatus::kInfeasible;
         }
      }
      else if( firstconsEquality && !secondconsEquality )
      {
         // l1 = r1, s > 0
         // if r1 not in [s*l2, s*r2], infeasible
         // else inequality redundant
         if( ( adjustedRHSInf ||
               num.isFeasLE( rhs_values[row1], adjustedRHS ) ) &&
             ( adjustedLHSInf ||
               num.isFeasGE( rhs_values[row1], adjustedLHS ) ) )
         {
            TransactionGuard<REAL> guard{ reductions };
            reductions.lockRow( row1 );
            reductions.lockRow( row2 );
            reductions.markRowRedundant( row2 );
         }
         else
         {
            result = PresolveStatus::kInfeasible;
         }
      }
      else if( !firstconsEquality && secondconsEquality )
      {
         // same as previous case, but row1 and row2 are flipped
         if( ( rflags[row1].test( RowFlag::kRhsInf ) ||
               num.isFeasGE( rhs_values[row1], adjustedRHS ) ) &&
             ( rflags[row1].test( RowFlag::kLhsInf ) ||
               num.isFeasLE( lhs_values[row1], adjustedRHS ) ) )
         {
            TransactionGuard<REAL> guard{ reductions };
            reductions.lockRow( row1 );
            reductions.lockRow( row2 );
            reductions.markRowRedundant( row1 );
         }
         else
         {
            result = PresolveStatus::kInfeasible;
         }
      }
      else
      {
         // l1 != r1, l2 != r2, s > 0
         // if [l1, r1] inter [s*l2, s*r2] emplty, infeasible
         // else we keep one row and give it the thightest bounds
         if( ( !rflags[row1].test( RowFlag::kRhsInf ) && !adjustedLHSInf &&
               num.isFeasLT( rhs_values[row1], adjustedLHS ) ) ||
             ( !rflags[row1].test( RowFlag::kLhsInf ) && !adjustedRHSInf &&
               num.isFeasLT( adjustedRHS, lhs_values[row1] ) ) )
         {
            result = PresolveStatus::kInfeasible;
         }
         else
         {
            TransactionGuard<REAL> guard{ reductions };
            reductions.lockRow( row1 );
            reductions.lockRow( row2 );

            if( !adjustedRHSInf && ( rflags[row1].test( RowFlag::kRhsInf ) ||
                                     adjustedRHS < rhs_values[row1] ) )
               reductions.changeRowRHS( row1, adjustedRHS );

            if( !adjustedLHSInf && ( rflags[row1].test( RowFlag::kLhsInf ) ||
                                     adjustedLHS > lhs_values[row1] ) )
               reductions.changeRowLHS( row1, adjustedLHS );

            reductions.markRowRedundant( row2 );
         }
      }
   };

   assert( nrows > 0 );

   std::unique_ptr<unsigned int[]> supportid{ new unsigned int[nrows] };
   std::unique_ptr<unsigned int[]> coefhash{ new unsigned int[nrows] };
   std::unique_ptr<int[]> row{ new int[nrows] };

#if 0
      tbb::parallel_invoke(
          [nrows, &row]() {
             for( int i = 0; i < nrows; ++i )
                row[i] = i;
          },
          [&constMatrix, &coefhash, this]() {
             computeRowHashes( constMatrix, coefhash.get() );
          },
          [&constMatrix, &supportid, this]() {
             computeSupportIdParallel( constMatrix, supportid.get() );
          } );
#else

   tbb::parallel_invoke(
       [nrows, &row]() {
          for( int i = 0; i < nrows; ++i )
             row[i] = i;
       },
       [&constMatrix, &coefhash, this]() {
          computeRowHashes( constMatrix, coefhash.get() );
       },
       [&constMatrix, &supportid, this]() {
          computeSupportId( constMatrix, supportid.get() );
       } );
#endif

   pdqsort( row.get(), row.get() + nrows, [&]( int a, int b ) {
      return supportid[a] < supportid[b] ||
             ( supportid[a] == supportid[b] && coefhash[a] < coefhash[b] ) ||
             ( supportid[a] == supportid[b] && coefhash[a] == coefhash[b] &&
               rowperm[a] < rowperm[b] );
   } );

   Vec<std::tuple<int, int, REAL>> parallelRows;

   for( int i = 0; i < nrows; )
   {
      // determine size of bucket
      int j;
      for( j = i + 1; j < nrows; ++j )
      {
         if( coefhash[row[i]] != coefhash[row[j]] ||
             supportid[row[i]] != supportid[row[j]] )
            break;
      }
      int len = j - i;

      // if more  than one row is in the bucket try to find parallel rows
      if( len > 1 )
      {
         // fmt::print( "bucket of length {} starting at {}\n", len, i );
         findParallelRows( num, row.get() + i, len, constMatrix, parallelRows );
      }

      assert( j > i );
      // set i to start of next bucket
      i = j;
   }

   if( !parallelRows.empty() )
   {
      pdqsort( parallelRows.begin(), parallelRows.end(),
               [&rowperm]( const std::tuple<int, int, REAL>& a,
                           const std::tuple<int, int, REAL>& b ) {
                  return std::make_pair( rowperm[std::get<0>( a )],
                                         rowperm[std::get<1>( a )] ) <
                         std::make_pair( rowperm[std::get<0>( b )],
                                         rowperm[std::get<1>( b )] );
               } );

      result = PresolveStatus::kReduced;

      for( const std::tuple<int, int, REAL>& parallelRow : parallelRows )
      {
         int row1;
         int row2;
         REAL ratio;

         std::tie( row1, row2, ratio ) = parallelRow;

         handlerows( row1, row2, ratio );

         if( result == PresolveStatus::kInfeasible )
            break;
      }
   }

   return result;
}

} // namespace papilo

#endif
