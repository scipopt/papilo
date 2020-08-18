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

#ifndef _PAPILO_PRESOLVERS_SPARSIFY_HPP_
#define _PAPILO_PRESOLVERS_SPARSIFY_HPP_

#include "papilo/core/PresolveMethod.hpp"
#include "papilo/core/Problem.hpp"
#include "papilo/core/ProblemUpdate.hpp"
#include "papilo/core/SingleRow.hpp"
#include "papilo/misc/Hash.hpp"
#include "papilo/misc/tbb.hpp"

namespace papilo
{

template <typename REAL>
class Sparsify : public PresolveMethod<REAL>
{
   double maxscale = 1000;

   using HitCount = uint16_t;

   struct SparsifyData
   {
      Vec<HitCount> candrowhits;
      Vec<int> candrows;
      Vec<std::pair<int, REAL>> sparsify;
      Vec<std::tuple<int, int, int>> reductionBuffer;

      SparsifyData( int nrows ) : candrowhits( nrows )
      {
         candrows.reserve( nrows );
      }
   };

 public:
   Sparsify() : PresolveMethod<REAL>()
   {
      this->setName( "sparsify" );
      this->setTiming( PresolverTiming::kExhaustive );
      this->setDelayed( true );
   }

   void
   addPresolverParams( ParameterSet& paramSet ) override
   {
      paramSet.addParameter(
          "sparsify.maxscale",
          "maximum absolute scale to use for cancelling nonzeros",
          this->maxscale, 1.0 );
   }

   virtual PresolveStatus
   execute( const Problem<REAL>& problem,
            const ProblemUpdate<REAL>& problemUpdate, const Num<REAL>& num,
            Reductions<REAL>& reductions ) override;
};

#ifdef PAPILO_USE_EXTERN_TEMPLATES
extern template class Sparsify<double>;
extern template class Sparsify<Quad>;
extern template class Sparsify<Rational>;
#endif

template <typename REAL>
PresolveStatus
Sparsify<REAL>::execute( const Problem<REAL>& problem,
                         const ProblemUpdate<REAL>& problemUpdate,
                         const Num<REAL>& num, Reductions<REAL>& reductions )
{
   // go over the rows and get the equalities, extract the columns that
   // verify the conditions add them to a hash map, loop over the hash map
   // and compute the implied bounds and finally look for implied free
   // variables and add reductions
   const auto& domains = problem.getVariableDomains();
   const auto& lower_bounds = domains.lower_bounds;
   const auto& upper_bounds = domains.upper_bounds;
   const auto& cflags = domains.flags;

   const auto& activities = problem.getRowActivities();

   const ConstraintMatrix<REAL>& consmatrix = problem.getConstraintMatrix();

   const auto& lhs_values = consmatrix.getLeftHandSides();
   const auto& rhs_values = consmatrix.getRightHandSides();
   const auto& rflags = consmatrix.getRowFlags();
   const auto& rowsize = consmatrix.getRowSizes();
   const auto& colsize = consmatrix.getColSizes();
   const auto& nrows = consmatrix.getNRows();
   const auto& ncols = consmatrix.getNCols();

   auto isBinaryCol = [&]( int col ) {
      return cflags[col].test( ColFlag::kIntegral ) &&
             !cflags[col].test( ColFlag::kUnbounded ) &&
             lower_bounds[col] == 0 && upper_bounds[col] == 1;
   };

   PresolveStatus result = PresolveStatus::kUnchanged;

   // after each call skip more rounds to not call sparsify too often
   this->skipRounds( this->getNCalls() );

   // TODO: identify the equality rows, and sort them by length
   // loop over the sorted rows, and loop over the columns in reverse order
   // check numerical condition, implied freeness (check if row implies a
   // tighter bound, don't check it if it does later) and we substitute the
   // first col that we can and break (use a vector instead of a hashmap) and
   // add reductions in place
   // add a map that tells if the implied bounds of a col was already
   // computed and failed (using dynamicbitset)

   Vec<int> equalities;
   equalities.reserve( nrows );

   for( int i = 0; i < nrows; ++i )
   {
      if( rflags[i].test( RowFlag::kRedundant ) ||
          !rflags[i].test( RowFlag::kEquation ) || rowsize[i] <= 1 ||
          rowsize[i] > std::numeric_limits<HitCount>::max() )
         continue;

      assert( !rflags[i].test( RowFlag::kLhsInf, RowFlag::kRhsInf ) &&
              lhs_values[i] == rhs_values[i] );

      equalities.emplace_back( i );
   }

   tbb::combinable<SparsifyData> sparsifyData(
       [nrows]() { return SparsifyData( nrows ); } );

   tbb::parallel_for(
       tbb::blocked_range<int>( 0, static_cast<int>( equalities.size() ) ),
       [&]( const tbb::blocked_range<int>& r ) {
          std::size_t sparsifyStart;

          SparsifyData& localData = sparsifyData.local();

          auto& candrowhits = localData.candrowhits;
          auto& candrows = localData.candrows;
          auto& sparsify = localData.sparsify;
          auto& reductionBuffer = localData.reductionBuffer;

          for( int i = r.begin(); i != r.end(); ++i )
          {
             int eqrow = equalities[i];

             auto rowvec = consmatrix.getRowCoefficients( eqrow );

             int eqlen = rowvec.getLength();
             const int* eqcols = rowvec.getIndices();
             bool cancelint = true;
             int minhits = eqlen - 1;
             int nint = 0;
             Message::debug(
                 this,
                 "trying sparsification with equality row {} of length {}\n",
                 eqrow, eqlen );

             if( problem.getNumIntegralCols() != 0 )
             {
                int ncont = 0;
                int nbin = 0;

                for( int i = 0; i != eqlen; ++i )
                {
                   int col = eqcols[i];

                   if( !cflags[col].test( ColFlag::kIntegral ) )
                      ++ncont;
                   else if( isBinaryCol( col ) )
                      ++nbin;
                   else
                   {
                      ++nint;
                      continue;
                   }

                   auto colvec = consmatrix.getColumnCoefficients( col );
                   const int* colrows = colvec.getIndices();
                   int collen = colvec.getLength();

                   for( int j = 0; j != collen; ++j )
                   {
                      int row = colrows[j];

                      if( row == eqrow )
                         continue;

                      if( candrowhits[row] == 0 )
                      {
                         if( nbin + ncont > 2 )
                            continue;

                         candrows.push_back( row );
                      }

                      ++candrowhits[row];
                   }
                }

                if( nbin + nint == 0 )
                {
                   auto it = std::remove_if( candrows.begin(), candrows.end(),
                                             [&]( int r ) {
                                                if( candrowhits[r] < ncont - 1 )
                                                {
                                                   candrowhits[r] = 0;
                                                   return true;
                                                }
                                                return false;
                                             } );

                   cancelint = false;

                   candrows.erase( it, candrows.end() );
                }
                else
                {
                   auto it = std::remove_if(
                       candrows.begin(), candrows.end(), [&]( int r ) {
                          if( candrowhits[r] < nbin + ncont - 1 )
                          {
                             candrowhits[r] = 0;
                             return true;
                          }
                          if( cancelint )
                             candrowhits[r] = nbin + ncont;
                          return false;
                       } );

                   candrows.erase( it, candrows.end() );

                   minhits = eqlen;
                }
             }

             if( problem.getNumIntegralCols() == 0 || nint != 0 )
             {
                for( int i = 0; i != eqlen; ++i )
                {
                   int col = eqcols[i];

                   if( problem.getNumIntegralCols() != 0 &&
                       ( !cflags[col].test( ColFlag::kIntegral ) ||
                         isBinaryCol( col ) ) )
                      continue;

                   auto colvec = consmatrix.getColumnCoefficients( col );
                   const int* colrows = colvec.getIndices();
                   int collen = colvec.getLength();

                   for( int j = 0; j != collen; ++j )
                   {
                      int row = colrows[j];

                      if( row == eqrow )
                         continue;

                      if( candrowhits[row] == 0 )
                      {
                         if( i > eqlen - minhits )
                            continue;

                         candrows.push_back( row );
                      }

                      ++candrowhits[row];
                   }
                }

                auto it = std::remove_if( candrows.begin(), candrows.end(),
                                          [&]( int r ) {
                                             if( candrowhits[r] < minhits )
                                             {
                                                candrowhits[r] = 0;
                                                return true;
                                             }
                                             return false;
                                          } );

                candrows.erase( it, candrows.end() );
             }

             if( !candrows.empty() )
             {
                Vec<REAL> scales( eqlen );
                const REAL* eqvals = rowvec.getValues();

                sparsifyStart = sparsify.size();
                sparsify.reserve( sparsifyStart + candrows.size() );

                for( int candrow : candrows )
                {
                   auto candrowvec = consmatrix.getRowCoefficients( candrow );
                   const int* candcols = candrowvec.getIndices();
                   const REAL* candvals = candrowvec.getValues();
                   int candlen = candrowvec.getLength();

                   if( !cancelint && candrowhits[candrow] != eqlen )
                   {
                      bool skip = false;
                      for( int j = 0; j != candlen; ++j )
                      {
                         if( cflags[candcols[j]].test( ColFlag::kIntegral ) )
                         {
                            skip = true;
                            break;
                         }
                      }

                      if( skip )
                         continue;
                   }

                   int i = 0;
                   int j = 0;

                   int currcancel = 0;

                   while( i != eqlen && j != candlen )
                   {
                      if( eqcols[i] == candcols[j] )
                      {
                         scales[i] = -candvals[j] / eqvals[i];

                         ++i;
                         ++j;
                      }
                      else if( eqcols[i] < candcols[j] )
                      {
                         --currcancel;
                         scales[i] = 0;
                         ++i;
                      }
                      else
                      {
                         ++j;
                      }
                   }

                   while( i != eqlen )
                   {
                      --currcancel;
                      scales[i] = 0;
                      ++i;
                   }

                   pdqsort( scales.begin(), scales.end() );

                   int bestcancel = 0;
                   REAL bestscale = 0;

                   for( int k = 0; k != eqlen - 1; ++k )
                   {
                      if( scales[k] == 0 || abs( scales[k] ) > maxscale )
                         continue;

                      int ncancel = currcancel;

                      for( int j = k + 1; j != eqlen; ++j )
                      {
                         if( num.isEq( scales[k], scales[j] ) )
                            ++ncancel;
                         else
                            break;
                      }

                      if( ncancel > bestcancel )
                      {
                         bestcancel = ncancel;
                         bestscale = scales[k];
                      }
                   }

                   if( bestcancel > 0 )
                   {
                      Message::debug(
                          this,
                          "equation row{} cancels {} nonzeros on row{} "
                          "with scale {}\n",
                          eqrow, bestcancel, candrow, bestscale );

                      sparsify.emplace_back( candrow, bestscale );
                   }
                }

                for( int r : candrows )
                   candrowhits[r] = 0;
                candrows.clear();

                if( sparsify.size() != sparsifyStart )
                   reductionBuffer.emplace_back( eqrow, int( sparsifyStart ),
                                                 int( sparsify.size() ) );
             }
          }
       } );

   int nreductions = 0;
   sparsifyData.combine_each( [&]( const SparsifyData& localData ) {
      nreductions += localData.reductionBuffer.size();
   } );

   if( nreductions != 0 )
   {
      result = PresolveStatus::kReduced;

      Vec<std::tuple<int, int, std::pair<int, REAL>*>> reductionData;

      reductionData.reserve( nreductions );

      sparsifyData.combine_each( [&]( SparsifyData& localData ) {
         for( const std::tuple<int, int, int>& reductionTuple :
              localData.reductionBuffer )
         {
            int eqrow = std::get<0>( reductionTuple );
            int start = std::get<1>( reductionTuple );
            int end = std::get<2>( reductionTuple );
            reductionData.emplace_back( eqrow, end - start,
                                        &localData.sparsify[start] );
         }
      } );

      const Vec<int>& rowperm = problemUpdate.getRandomRowPerm();

      pdqsort( reductionData.begin(), reductionData.end(),
               [&]( const std::tuple<int, int, std::pair<int, REAL>*>& a,
                    const std::tuple<int, int, std::pair<int, REAL>*>& b ) {
                  int eqrowA = std::get<0>( a );
                  int eqrowB = std::get<0>( b );

                  int numCandsA = std::get<1>( a );
                  int numCandsB = std::get<1>( b );

                  return std::make_tuple( rowsize[eqrowA], -numCandsA,
                                          rowperm[eqrowA] ) <
                         std::make_tuple( rowsize[eqrowB], -numCandsB,
                                          rowperm[eqrowB] );
               } );

      for( const std::tuple<int, int, std::pair<int, REAL>*>& reductionTuple :
           reductionData )
      {
         int eqrow = std::get<0>( reductionTuple );
         int num = std::get<1>( reductionTuple );
         const std::pair<int, REAL>* sparsify = std::get<2>( reductionTuple );

         TransactionGuard<REAL> tg{ reductions };
         reductions.lockRow( eqrow );
         reductions.sparsify( eqrow, num, sparsify );
      }
   }

   Message::debug( this, "sparsify finished\n" );

   return result;
}

} // namespace papilo

#endif
