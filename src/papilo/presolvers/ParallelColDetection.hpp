/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*               This file is part of the program and library                */
/*    PaPILO --- Parallel Presolve for Integer and Linear Optimization       */
/*                                                                           */
/* Copyright (C) 2020-2025 Zuse Institute Berlin (ZIB)                       */
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

#include "papilo/Config.hpp"
#include "papilo/core/PresolveMethod.hpp"
#include "papilo/core/Problem.hpp"
#include "papilo/core/ProblemUpdate.hpp"
#include "papilo/external/pdqsort/pdqsort.h"
#include "papilo/misc/Hash.hpp"
#ifdef PAPILO_TBB
#include "papilo/misc/tbb.hpp"
#endif

namespace papilo
{

template <typename REAL>
class ParallelColDetection : public PresolveMethod<REAL>
{
   struct SupportHashCompare
   {
      SupportHashCompare() = default;

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
   findParallelCols( const Num<REAL>& num, const int* bucket, int bucketSize,
                     const ConstraintMatrix<REAL>& constMatrix,
                     const Vec<REAL>& obj, const VariableDomains<REAL>& domains,
                     const SymmetryStorage& symmetries,
                     bool is_binary, Reductions<REAL>& reductions );

   void
   computeColHashes( const ConstraintMatrix<REAL>& constMatrix,
                     const Vec<REAL>& obj, unsigned int* columnHashes );

   void
   computeSupportId( const ConstraintMatrix<REAL>& constMatrix,
                     unsigned int* supportHashes );

   void
   addPresolverParams( ParameterSet& paramSet ) override
   {
      paramSet.addParameter( "parallelcols.symmetries_enabled",
          "should Parallel Cols search for symmetries at the end (for example "
          "for binary problems where merging columns does not work)",
          symmetries );
   }

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
      if( presolveOptions.dualreds < 2 ||
          problem.test_problem_type( ProblemFlag::kBinary ) )
         this->setEnabled( false );
      return false;
   }

   PresolveStatus
   execute( const Problem<REAL>& problem,
            const ProblemUpdate<REAL>& problemUpdate,
            const Num<REAL>& num, Reductions<REAL>& reductions,
            const Timer& timer, int& reason_of_infeasibility)
       override;

   PresolveStatus
   execute_symmetries( const Problem<REAL>& problem,
                       const ProblemUpdate<REAL>& problemUpdate,
                       const Num<REAL>& num, Reductions<REAL>& reductions,
                       const Timer& timer ) override;

 private:
   int
   determineBucketSize( int nColumns,
                        std::unique_ptr<unsigned int[]>& supportid,
                        std::unique_ptr<unsigned int[]>& coefficentHashes,
                        std::unique_ptr<int[]>& column, int i );

   bool
   check_parallelity( const Num<REAL>& num, const Vec<REAL>& obj, int col1,
                      int length, const REAL* coefs1, int col2,
                      const REAL* coefs2 ) const;

   bool
   can_be_merged( const Num<REAL>& num, const Vec<REAL>& lbs,
                  const Vec<REAL>& ubs, int col1, const REAL* coefs1,
                  const REAL* coefs2, const Vec<ColFlags>& cflags ) const;

   bool
   determineOderingForZeroObj( REAL val1, REAL val2, int colpermCol1,
            int colpermCol2 ) const;

   bool symmetries = false;
};

#ifdef PAPILO_USE_EXTERN_TEMPLATES
extern template class ParallelColDetection<double>;
extern template class ParallelColDetection<Quad>;
extern template class ParallelColDetection<Rational>;
#endif

template <typename REAL>
void
ParallelColDetection<REAL>::findParallelCols(
    const Num<REAL>& num, const int* bucket, int bucketSize,
    const ConstraintMatrix<REAL>& constMatrix, const Vec<REAL>& obj,
    const VariableDomains<REAL>& domains, const SymmetryStorage& symmetries,
    bool is_binary, Reductions<REAL>& reductions )
{
   // TODO if bucketSize too large do gurobi trick
   const Vec<ColFlags>& cflags = domains.flags;
   const Vec<REAL>& lbs = domains.lower_bounds;
   const Vec<REAL>& ubs = domains.upper_bounds;

   auto checkDomainForHoles = [&]( int col, const REAL& scale )
   {
      assert( cflags[col].test( ColFlag::kIntegral ) );
      assert( !cflags[col].test( ColFlag::kLbInf, ColFlag::kUbInf ) );
      assert( num.isIntegral(scale) );

      // there are holes iff the range of the unit column does not cover the scale
      return num.isLT(ubs[col] - lbs[col], abs(scale) - 1);
   };

   if( constMatrix.getColumnCoefficients( bucket[0] ).getLength() <= 1 )
      return;

   assert( !cflags[bucket[bucketSize - 1]].test( ColFlag::kInactive ) );
   // only continuous columns
   if( !cflags[bucket[bucketSize - 1]].test( ColFlag::kIntegral ) )
   {
      assert(!is_binary);
      int col1 = bucket[0];
      auto col1vec = constMatrix.getColumnCoefficients( col1 );
      const REAL* coefs1 = col1vec.getValues();
      const int length = col1vec.getLength();

      for( int j = 1; j < bucketSize; ++j )
      {
         int col2 = bucket[j];
         assert( !cflags[col2].test( ColFlag::kIntegral ) );
         const REAL* coefs2 = constMatrix.getColumnCoefficients( col2 ).getValues();
         assert( num.isLE( abs( coefs1[0] ), abs( coefs2[0] ) ) );

         if( check_parallelity( num, obj, col2, length, coefs2, col1, coefs1 ) )
         {
            TransactionGuard<REAL> tg{ reductions };
            reductions.lockCol( col2 );
            reductions.lockCol( col1 );
            reductions.mark_parallel_cols( col2, col1 );
         }
      }
      return;
   }
   // only integer columns
   else if( cflags[bucket[0]].test( ColFlag::kIntegral ) )
   {
      // to avoid adding redundant tsxs in case of multiple parallel cols
      // all tsxs only refer to the first col in the bucket IF it is possible
      // to construct a tsxs between the first and the remaining ones
      bool abort_after_first_loop = !is_binary;
      int first_col_with_finite_bounds = -1;

      for( int i = 0; i < bucketSize; ++i )
      {
         int col1 = bucket[i];
         assert( cflags[col1].test( ColFlag::kIntegral ) );

         if( cflags[col1].test( ColFlag::kLbInf, ColFlag::kUbInf ) )
            continue;

         if( first_col_with_finite_bounds == -1 )
            first_col_with_finite_bounds = i;

         if( i == first_col_with_finite_bounds + 1 && abort_after_first_loop )
            break;

         auto col1vec = constMatrix.getColumnCoefficients( col1 );
         const REAL* coefs1 = col1vec.getValues();
         const int length = col1vec.getLength();

         for( int j = i + 1; j < bucketSize; ++j )
         {
            int col2 = bucket[j];
            assert( cflags[col2].test( ColFlag::kIntegral ) );

            if( cflags[col2].test( ColFlag::kLbInf, ColFlag::kUbInf ) )
               continue;

            auto col2vec = constMatrix.getColumnCoefficients( col2 );
            assert( col2vec.getLength() == length );
            const REAL* coefs2 = col2vec.getValues();
            assert( num.isLE( abs( coefs1[0] ), abs( coefs2[0] ) ) );
            REAL scale = coefs2[0] / coefs1[0];

            if( num.isIntegral(scale)
                  && check_parallelity( num, obj, col2, length, coefs2, col1, coefs1 )
                  && !checkDomainForHoles( col1, scale )
//                  && !symmetries.contains_symmetry( col2, col1 )
                  )
            {
               TransactionGuard<REAL> tg{ reductions };
               reductions.lockCol( col2 );
               reductions.lockCol( col1 );
               reductions.mark_parallel_cols( col2, col1 );
            }
            else
               abort_after_first_loop = false;
         }
      }
      return;
   }
   else
   {
      assert(!is_binary);
      int col1 = bucket[0];
      auto col1vec = constMatrix.getColumnCoefficients( col1 );
      const REAL* coefs1 = col1vec.getValues();
      const int length = col1vec.getLength();
      int ncontinuous_cols;

      assert( !cflags[col1].test( ColFlag::kIntegral ) );
      {
         Vec<std::pair<int, int>> continuous_parallel_cols;

         for( ncontinuous_cols = 1; ncontinuous_cols < bucketSize; ++ncontinuous_cols )
         {
            int col2 = bucket[ncontinuous_cols];

            if( cflags[col2].test( ColFlag::kIntegral ) )
               break;

            const REAL* coefs2 = constMatrix.getColumnCoefficients( col2 ).getValues();
            assert( num.isLE( abs( coefs1[0] ), abs( coefs2[0] ) ) );

            if( check_parallelity( num, obj, col2, length, coefs2, col1, coefs1 ) )
               continuous_parallel_cols.emplace_back( col2, col1 );
         }

         // TODO: use the updated bounds of the cont parallel columns
         int col2 = bucket[ncontinuous_cols];
         assert( cflags[col2].test( ColFlag::kIntegral ) );
         const REAL* coefs2 = constMatrix.getColumnCoefficients( col2 ).getValues();

         if( can_be_merged( num, lbs, ubs, col1, coefs1, coefs2, cflags )
               && check_parallelity( num, obj, col1, length, coefs1, col2, coefs2 ) )
         {
            TransactionGuard<REAL> tg{ reductions };
            reductions.lockCol( col1 );
            for( std::pair<int, int> pair : continuous_parallel_cols )
               reductions.lockCol( pair.first );
            reductions.lockCol( col2 );
            reductions.lockColBounds( col1 );
            for( std::pair<int, int> pair : continuous_parallel_cols )
               reductions.mark_parallel_cols( pair.first, pair.second );
            reductions.mark_parallel_cols( col1, col2 );
         }
         else
         {
            for( std::pair<int, int> pair : continuous_parallel_cols )
            {
               TransactionGuard<REAL> tg{ reductions };
               reductions.lockCol( pair.first );
               reductions.lockCol( pair.second );
               reductions.mark_parallel_cols( pair.first, pair.second );
            }
         }
      }

      for( int i = ncontinuous_cols; i < bucketSize; ++i )
      {
         col1 = bucket[i];
         assert( cflags[col1].test( ColFlag::kIntegral ) );

         if( cflags[col1].test( ColFlag::kLbInf, ColFlag::kUbInf ) )
            continue;

         col1vec = constMatrix.getColumnCoefficients( col1 );
         assert( col1vec.getLength() == length );
         coefs1 = col1vec.getValues();

         for( int j = i + 1; j < bucketSize; ++j )
         {
            int col2 = bucket[j];
            assert( cflags[col2].test( ColFlag::kIntegral ) );

            if( cflags[col2].test( ColFlag::kLbInf, ColFlag::kUbInf ) )
               continue;

            auto col2vec = constMatrix.getColumnCoefficients( col2 );
            assert( col2vec.getLength() == length );
            const REAL* coefs2 = col2vec.getValues();
            assert( num.isLE( abs( coefs1[0] ), abs( coefs2[0] ) ) );
            REAL scale = coefs2[0] / coefs1[0];

            if( num.isIntegral(scale)
                  && check_parallelity( num, obj, col2, length, coefs2, col1, coefs1 )
                  && !checkDomainForHoles( col1, scale )
//                  && !symmetries.contains_symmetry( col2, col1 )
                  )
            {
               TransactionGuard<REAL> tg{ reductions };
               reductions.lockCol( col2 );
               reductions.lockCol( col1 );
               reductions.mark_parallel_cols( col2, col1 );
            }
         }
      }
      return;
   }
}
template <typename REAL>
bool
ParallelColDetection<REAL>::can_be_merged(
    const Num<REAL>& num, const Vec<REAL>& lbs, const Vec<REAL>& ubs, int col1,
    const REAL* coefs1, const REAL* coefs2, const Vec<ColFlags>& cflags ) const
{
   return cflags[col1].test( ColFlag::kLbInf, ColFlag::kUbInf ) ||
          !num.isLT(
              abs( ( ubs[col1] - lbs[col1] ) * coefs1[0] / coefs2[0] ), 1 );
}

template <typename REAL>
bool
ParallelColDetection<REAL>::check_parallelity( const Num<REAL>& num,
                                               const Vec<REAL>& obj, int col1,
                                               int length, const REAL* coefs1,
                                               int col2,
                                               const REAL* coefs2 ) const
{
   REAL scale = coefs1[0] / coefs2[0];
   if( !num.isEq( obj[col1], obj[col2] * scale ) )
      return false;
   for( int k = 1; k < length; ++k )
      if( !num.isEq( coefs1[k], coefs2[k] * scale ) )
         return false;
   return true;
}

template <typename REAL>
void
ParallelColDetection<REAL>::computeColHashes(
    const ConstraintMatrix<REAL>& constMatrix, const Vec<REAL>& obj,
    unsigned int* columnHashes )
{
#ifdef PAPILO_TBB
   tbb::parallel_for(
       tbb::blocked_range<int>( 0, constMatrix.getNCols() ),
       [&]( const tbb::blocked_range<int>& r ) {
          for( int i = r.begin(); i < r.end(); ++i )
#else
   for( int i = 0; i< constMatrix.getNCols(); i++)
#endif
          {
             // compute hash-value for coefficients
             auto columnCoefficients = constMatrix.getColumnCoefficients( i );
             const REAL* values = columnCoefficients.getValues();
             const int len = columnCoefficients.getLength();

             Hasher<unsigned int> hasher( len );

             if( len > 1 )
             {
                // compute scale such that the first coefficient is
                // positive 1/golden ratio. The choice of
                // the constant is arbitrary and is used to make cases
                // where two coefficients that are equal
                // within epsilon get different values are
                // more unlikely by choose some irrational number
                REAL scale = REAL( 2.0 / ( 1.0 + sqrt( 5.0 ) ) ) / values[0];

                // add scaled coefficients of other row
                // entries to compute the hash
                for( int j = 1; j != len; ++j )
                {
                   hasher.addValue( Num<REAL>::hashCode( values[j] * scale ) );
                }
                if( obj[i] != 0 )
                   hasher.addValue( Num<REAL>::hashCode( obj[i] * scale ) );
             }

             columnHashes[i] = hasher.getHash();
          }
#ifdef PAPILO_TBB
       } );
#endif
}

template <typename REAL>
void
ParallelColDetection<REAL>::computeSupportId(
    const ConstraintMatrix<REAL>& constMatrix, unsigned int* supportHashes )
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
         supportHashes[i] = i;
      else // support already exists, use the previous support id
         supportHashes[i] = insResult.first->second;
   }
}

template <typename REAL>
PresolveStatus
ParallelColDetection<REAL>::execute( const Problem<REAL>& problem,
                                     const ProblemUpdate<REAL>& problemUpdate,
                                     const Num<REAL>& num, Reductions<REAL>& reductions,
                                     const Timer& timer, int& reason_of_infeasibility){
   const auto& constMatrix = problem.getConstraintMatrix();
   const auto& obj = problem.getObjective().coefficients;
   const auto& cflags = problem.getColFlags();
   const auto& symmetries = problem.getSymmetries();
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

#ifdef PAPILO_TBB
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
#else
   for( int i = 0; i < ncols; ++i )
      col[i] = i;
   computeColHashes( constMatrix, obj, coefhash.get() );
   computeSupportId( constMatrix, supportid.get() );
#endif

   pdqsort(
       col.get(), col.get() + ncols,
       [&]( int a, int b )
       {
          if( cflags[a].test( ColFlag::kInactive ) && cflags[b].test( ColFlag::kInactive ) )
             return a < b;
          if( cflags[a].test( ColFlag::kInactive ) )
             return true;
          if( cflags[b].test( ColFlag::kInactive ) )
             return false;
          assert(constMatrix.getColumnCoefficients( a ).getLength() > 0);
          assert(constMatrix.getColumnCoefficients( b ).getLength() > 0);

          if( supportid[a] < supportid[b] ||
              ( supportid[a] == supportid[b] && coefhash[a] < coefhash[b] ) )
             return true;
          else if( !( supportid[a] == supportid[b] &&
                        coefhash[a] == coefhash[b] ) )
             return false;
          assert( supportid[a] == supportid[b] && coefhash[a] == coefhash[b] );

          bool flag_a_integer = cflags[a].test( ColFlag::kIntegral );
          bool flag_b_integer = cflags[b].test( ColFlag::kIntegral );
          if( flag_a_integer != flag_b_integer )
             return !flag_a_integer;

          return // sort by scale factor
              abs( obj[a] ) < abs( obj[b] ) ||
              // sort by scale factor if obj is zero
              ( abs( obj[a] ) == abs( obj[b] ) && obj[a] == 0 &&
                determineOderingForZeroObj(
                    constMatrix.getColumnCoefficients( a ).getValues()[0],
                    constMatrix.getColumnCoefficients( b ).getValues()[0],
                    colperm[a], colperm[b] ) ) ||
              // sort by permutation
              ( abs( obj[a] ) == abs( obj[b] ) && obj[a] != 0 &&
                colperm[a] < colperm[b] );
       } );

   const bool is_binary = problem.test_problem_type( ProblemFlag::kBinary );
   for( int i = 0; i < ncols; )
   {
      int bucketSize =
          determineBucketSize( ncols, supportid, coefhash, col, i );

      // if more than one col is in the bucket find parallel cols
      if( bucketSize > 1 )
         findParallelCols( num, col.get() + i, bucketSize, constMatrix, obj,
                           problem.getVariableDomains(), symmetries, is_binary, reductions );
      i = i + bucketSize;
   }
   if( reductions.getTransactions().size() > 0 )
      result = PresolveStatus::kReduced;
   return result;
}

template <typename REAL>
PresolveStatus
ParallelColDetection<REAL>::execute_symmetries( const Problem<REAL>& problem,
                    const ProblemUpdate<REAL>& problemUpdate,
                    const Num<REAL>& num, Reductions<REAL>& reductions,
                    const Timer& timer )
{
   if( !symmetries )
      return PresolveStatus::kUnchanged;
   if( this->isEnabled() )
   {
      fmt::print( "For Symmetries parallel columns need to be "
                  "disabled!\n" );
      return PresolveStatus::kUnchanged;
   }

   int cause = -1;
   return execute( problem, problemUpdate, num, reductions, timer,  cause);
}

template <typename REAL>
int
ParallelColDetection<REAL>::determineBucketSize(
    int nColumns, std::unique_ptr<unsigned int[]>& supportid,
    std::unique_ptr<unsigned int[]>& coefficentHashes,
    std::unique_ptr<int[]>& column, int i )
{
   int j;
   for( j = i + 1; j < nColumns; ++j )
   {
      if( coefficentHashes[column[i]] != coefficentHashes[column[j]] ||
          supportid[column[i]] != supportid[column[j]] )
      {
         break;
      }
   }
   assert( j > i );
   return j - i;
}

template <typename REAL>
bool
ParallelColDetection<REAL>::determineOderingForZeroObj( REAL val1, REAL val2,
                                     int colpermCol1, int colpermCol2 ) const
{
   if(val1 == val2)
      return colpermCol1 < colpermCol2;
   return abs(val1) < abs(val2);
}

} // namespace papilo

#endif
