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

#ifndef _PAPILO_PRESOLVERS_DOMINATED_COLS_HPP_
#define _PAPILO_PRESOLVERS_DOMINATED_COLS_HPP_

#include "papilo/Config.hpp"
#include "papilo/core/PresolveMethod.hpp"
#include "papilo/core/Problem.hpp"
#include "papilo/core/ProblemUpdate.hpp"
#include "papilo/core/SingleRow.hpp"
#include "papilo/misc/Hash.hpp"
#include "papilo/misc/Signature.hpp"
#ifdef PAPILO_TBB
#include "papilo/misc/tbb.hpp"
#endif

namespace papilo
{

template <typename REAL>
class DominatedCols : public PresolveMethod<REAL>
{
 public:
   DominatedCols() : PresolveMethod<REAL>()
   {
      this->setName( "domcol" );
      /// argument case can be primal or dual depending on the sign of the domination
      this->setArgument( ArgumentType::kPrimal );
      this->setTiming( PresolverTiming::kExhaustive );
   }

   bool
   initialize( const Problem<REAL>& problem,
               const PresolveOptions& presolveOptions ) override
   {
      if( presolveOptions.dualreds < 2 )
         this->setEnabled( false );
      return false;
   }

   /// stores implied bound information and signatures for a column
   struct ColInfo
   {
      Signature32 pos;
      Signature32 neg;
      int lbfree = 0;
      int ubfree = 0;

      Signature32
      getNegSignature( int scale ) const
      {
         assert( scale == 1 || scale == -1 );
         return scale == 1 ? neg : pos;
      }

      Signature32
      getPosSignature( int scale ) const
      {
         assert( scale == 1 || scale == -1 );
         return scale == 1 ? pos : neg;
      }

      bool
      allowsDomination( int scale, const ColInfo& other, int otherscale ) const
      {
         return getNegSignature( scale ).isSuperset(
                    other.getNegSignature( otherscale ) ) &&
                getPosSignature( scale ).isSubset(
                    other.getPosSignature( otherscale ) );
      }
   };

   struct DomcolReduction
   {
      int col1;
      int col2;
      int implrowlock;
      BoundChange boundchg;
   };

   PresolveStatus
   execute( const Problem<REAL>& problem,
            const ProblemUpdate<REAL>& problemUpdate, const Num<REAL>& num,
            Reductions<REAL>& reductions, const Timer& timer,
            int& reason_of_infeasibility ) override;
};

#ifdef PAPILO_USE_EXTERN_TEMPLATES
extern template class DominatedCols<double>;
extern template class DominatedCols<Quad>;
extern template class DominatedCols<Rational>;
#endif

template <typename REAL>
PresolveStatus
DominatedCols<REAL>::execute( const Problem<REAL>& problem,
                              const ProblemUpdate<REAL>& problemUpdate,
                              const Num<REAL>& num, Reductions<REAL>& reductions,
                              const Timer& timer, int& reason_of_infeasibility){
   const size_t ncols = problem.getNCols();

   //TODO: Reduce skips by one to rerun initially
   // do not call dominated column presolver too often, since it can be
   // expensive
   this->skipRounds( this->getNCalls() );

   if( ncols <= 1 )
      return PresolveStatus::kUnchanged;

   const auto& obj = problem.getObjective().coefficients;
   const auto& consMatrix = problem.getConstraintMatrix();
   const auto& lbValues = problem.getLowerBounds();
   const auto& ubValues = problem.getUpperBounds();
   const auto& lhsValues = consMatrix.getLeftHandSides();
   const auto& rhsValues = consMatrix.getRightHandSides();
   const auto& rflags = consMatrix.getRowFlags();
   const auto& cflags = problem.getColFlags();
   const auto& activities = problem.getRowActivities();
   const auto& rowsize = consMatrix.getRowSizes();
   const size_t nrows = problem.getNRows();
   Vec<ColInfo> colinfo( ncols );

#ifdef PAPILO_TBB
   std::atomic<int> nunboundedcols { };
#else
   int nunboundedcols = 0;
#endif

   // compute signatures and implied bound information of all columns in
   // parallel
#ifdef PAPILO_TBB
   tbb::parallel_for(
       tbb::blocked_range<unsigned int>( 0, ncols ),
       [&]( const tbb::blocked_range<unsigned int>& r ) {
          for( unsigned int col = r.begin(); col < r.end(); ++col )
#else
   for( int col = 0; col < (int)ncols; ++col )
#endif
          {
             auto colvec = consMatrix.getColumnCoefficients( col );
             int collen = colvec.getLength();
             const int* colrows = colvec.getIndices();
             const REAL* colvals = colvec.getValues();

             if( cflags[col].test( ColFlag::kLbInf ) )
                colinfo[col].lbfree = -1;
             if( cflags[col].test( ColFlag::kUbInf ) )
                colinfo[col].ubfree = -1;

             for( int j = 0; j < collen; ++j )
             {
                int row = colrows[j];
                if( colinfo[col].ubfree == 0 &&
                    row_implies_UB( num, lhsValues[row], rhsValues[row],
                                    rflags[row], activities[row], colvals[j],
                                    lbValues[col], ubValues[col],
                                    cflags[col] ) )
                   colinfo[col].ubfree = j + 1;

                if( colinfo[col].lbfree == 0 &&
                    row_implies_LB( num, lhsValues[row], rhsValues[row],
                                    rflags[row], activities[row], colvals[j],
                                    lbValues[col], ubValues[col],
                                    cflags[col] ) )
                   colinfo[col].lbfree = j + 1;

                if( !rflags[row].test( RowFlag::kLhsInf, RowFlag::kRhsInf ) )
                {
                   // ranged row or equality, add to positive and negative
                   // signature
                   colinfo[col].pos.add( row );
                   colinfo[col].neg.add( row );
                }
                else if( rflags[row].test( RowFlag::kLhsInf ) )
                {
                   // <= constraint, add to positive signature for positive
                   // coefficient and negative signature otherwise
                   if( colvals[j] < 0 )
                      colinfo[col].neg.add( row );
                   else
                      colinfo[col].pos.add( row );
                }
                else
                {
                   // >= constraint, add to positive signature for negative
                   // coefficient and negative signature otherwise
                   assert( rflags[row].test( RowFlag::kRhsInf ) );
                   if( colvals[j] < 0 )
                      colinfo[col].pos.add( row );
                   else
                      colinfo[col].neg.add( row );
                }
             }

             if( colinfo[col].lbfree != 0 || colinfo[col].ubfree != 0 )
                ++nunboundedcols;
          }
#ifdef PAPILO_TBB
       } );
#endif

   Vec<int> unboundedcols(nunboundedcols);

   nunboundedcols = 0;

   for( int col = 0; nunboundedcols < (int)unboundedcols.size(); ++col )
   {
      if( colinfo[col].lbfree != 0 || colinfo[col].ubfree != 0 )
      {
         unboundedcols[nunboundedcols] = col;
         ++nunboundedcols;
      }
   }

   auto checkDominance = [&]( int col1, int col2, int scal1, int scal2 ) {
      assert( !cflags[col1].test( ColFlag::kIntegral ) ||
              cflags[col2].test( ColFlag::kIntegral ) );

      // first check if the signatures rule out domination
      if( !colinfo[col1].allowsDomination( scal1, colinfo[col2], scal2 ) )
         return false;

      auto col1vec = consMatrix.getColumnCoefficients( col1 );
      int col1len = col1vec.getLength();
      const int* col1rows = col1vec.getIndices();
      const REAL* col1vals = col1vec.getValues();

      auto col2vec = consMatrix.getColumnCoefficients( col2 );
      int col2len = col2vec.getLength();
      const int* col2rows = col2vec.getIndices();
      const REAL* col2vals = col2vec.getValues();

      int i = 0;
      int j = 0;

      while( i != col1len && j != col2len )
      {
         REAL val1;
         REAL val2;
         RowFlags rowf;

         if( col1rows[i] == col2rows[j] )
         {
            val1 = col1vals[i] * scal1;
            val2 = col2vals[j] * scal2;
            rowf = rflags[col1rows[i]];

            ++i;
            ++j;
         }
         else if( col1rows[i] < col2rows[j] )
         {
            val1 = col1vals[i] * scal1;
            val2 = 0;
            rowf = rflags[col1rows[i]];
            ++i;
         }
         else
         {
            assert( col1rows[i] > col2rows[j] );
            val1 = 0;
            val2 = col2vals[j] * scal2;
            rowf = rflags[col2rows[j]];
            ++j;
         }

         if( !rowf.test( RowFlag::kLhsInf, RowFlag::kRhsInf ) )
         {
            // ranged row or equality, values must be equal
            if( !num.isEq( val1, val2 ) )
               return false;
         }
         else if( rowf.test( RowFlag::kLhsInf ) )
         {
            // <= constraint, col1 must have a smaller or equal coefficient
            if( num.isGT( val1, val2 ) )
               return false;
         }
         else
         {
            // >= constraint, col1 must have a larger or equal coefficient
            assert( rowf.test( RowFlag::kRhsInf ) );
            if( num.isLT( val1, val2 ) )
               return false;
         }
      }

      while( i != col1len )
      {
         REAL val1 = col1vals[i] * scal1;
         RowFlags rowf = rflags[col1rows[i]];
         ++i;

         if( !rowf.test( RowFlag::kLhsInf, RowFlag::kRhsInf ) )
         {
            // ranged row or equality, values must be equal
            return false;
         }
         else if( rowf.test( RowFlag::kLhsInf ) )
         {
            // <= constraint, col1 must have a smaller or equal coefficient
            if( num.isGT( val1, 0 ) )
               return false;
         }
         else
         {
            // >= constraint, col1 must have a larger or equal coefficient
            assert( rowf.test( RowFlag::kRhsInf ) );
            if( num.isLT( val1, 0 ) )
               return false;
         }
      }

      while( j != col2len )
      {
         REAL val2 = col2vals[j] * scal2;
         RowFlags rowf = rflags[col2rows[j]];
         ++j;

         if( !rowf.test( RowFlag::kLhsInf, RowFlag::kRhsInf ) )
         {
            // ranged row or equality, values must be equal
            return false;
         }
         else if( rowf.test( RowFlag::kLhsInf ) )
         {
            // <= constraint, col1 must have a smaller or equal coefficient
            if( num.isGT( 0, val2 ) )
               return false;
         }
         else
         {
            // >= constraint, col1 must have a larger or equal coefficient
            assert( rowf.test( RowFlag::kRhsInf ) );
            if( num.isLT( 0, val2 ) )
               return false;
         }
      }

      if(problemUpdate.getPresolveOptions().dualreds <= 1 && num.isEq( obj[col1], obj[col2] ) )
         return false;
      return true;
   };

   Vec<int> domcol(ncols, -1);
   Vec<int> nchildren(ncols, 0);
   Vec<int> leaves(0);
   size_t ndomcolsbound = ncols;

   for( int col = 0; col < (int)ncols; ++col )
   {
      if( cflags[col].test( ColFlag::kLbInf ) && cflags[col].test( ColFlag::kUbInf ) )
      {
         domcol[col] = -2;
         --ndomcolsbound;
      }
   }

   assert(ncols > 1);
   assert(nrows >= 1);

   if( ndomcolsbound >= ncols )
      ndomcolsbound = ncols - 1;

   Vec<DomcolReduction> domcols(0);
   Vec<Vec<DomcolReduction>> domcolsbuffers(0);
   size_t ndomcols = 0;
   size_t start = 0;

   // repeat finding and filtering dominations to bound memory demand
   while( ndomcols < ndomcolsbound && start < unboundedcols.size() )
   {
      int ndomcolsbuffers = 0;
      unsigned int base = start;
      unsigned int stopp;

      // find dominations until number of columns is reached
      while( ndomcols < ncols && start < unboundedcols.size() )
      {
         stopp = std::min(start + nrows, unboundedcols.size());
         ndomcolsbuffers = stopp - base;

         if( (int)domcolsbuffers.size() < ndomcolsbuffers )
            domcolsbuffers.resize(ndomcolsbuffers);

#ifdef PAPILO_TBB
   // scan unbounded columns if they dominate other columns
   tbb::parallel_for(
       tbb::blocked_range<unsigned int>( start, stopp ),
       [&]( const tbb::blocked_range<unsigned int>& r ) {
          for( unsigned int k = r.begin(); k < r.end(); ++k )
#else
   for( int k = start; k < (int)stopp; ++k )
#endif
          {
             int unbounded_col = unboundedcols[k];
             int lbfree = colinfo[unbounded_col].lbfree;
             int ubfree = colinfo[unbounded_col].ubfree;
             assert( lbfree != 0 || ubfree != 0 );
             auto colvec = consMatrix.getColumnCoefficients( unbounded_col );
             int collen = colvec.getLength();
             const int* colrows = colvec.getIndices();
             const REAL* colvals = colvec.getValues();
             int lbrowlock = lbfree > 0 ? colrows[lbfree - 1] : -1;
             int ubrowlock = ubfree > 0 ? colrows[ubfree - 1] : -1;
             int bestrowsize = std::numeric_limits<int>::max();
             int bestrowlock = std::numeric_limits<int>::max();
             int bestrow = -1;
             int bestscale = 0;
             int row;

             for( int i = 0; i < collen; ++i )
             {
                row = colrows[i];

                if( bestrowsize < rowsize[row] )
                   continue;

                // determine the scale of the dominating column depending on
                // whether the upper or lower bound is free, and remember which
                // row needs to be locked to protect the implied bound (if any)
                int rowlock = std::numeric_limits<int>::max();
                int scale = 0;

                if( ubfree != 0 && rowlock > ubrowlock
                      && !rflags[row].test( colvals[i] < 0 ? RowFlag::kLhsInf : RowFlag::kRhsInf ) )
                {
                   rowlock = ubrowlock;
                   scale = 1;
                }

                if( lbfree != 0 && rowlock > lbrowlock
                      && !rflags[row].test( colvals[i] < 0 ? RowFlag::kRhsInf : RowFlag::kLhsInf ) )
                {
                   rowlock = lbrowlock;
                   scale = -1;
                }

                assert(rowsize[row] >= 1);

                if( scale == 0 || ( bestrowsize == rowsize[row] && bestrowlock <= rowlock ) )
                   continue;

                if( rowsize[row] == 1 )
                {
                   bestrow = -1;
                   break;
                }
                else
                   bestrow = i;

                bestrowsize = rowsize[row];
                bestrowlock = rowlock;
                bestscale = scale;
             }

             if( bestrow == -1 )
                continue;

             row = colrows[bestrow];

             REAL scaled_val = colvals[bestrow] * bestscale;
             REAL scaled_obj = obj[unbounded_col] * bestscale;
             auto rowvec = consMatrix.getRowCoefficients( row );
             const int* rowcols = rowvec.getIndices();
             const REAL* rowvals = rowvec.getValues();

             for( int j = bestrowsize - 1; j >= 0; --j )
             {
                int col = rowcols[j];

                if( domcol[col] != -1 || col == unbounded_col
                      || ( cflags[unbounded_col].test( ColFlag::kIntegral )
                      && !cflags[col].test( ColFlag::kIntegral ) ) )
                   continue;

                bool to_lb = false;
                bool to_ub = false;

                if( !rflags[row].test( RowFlag::kLhsInf,
                                           RowFlag::kRhsInf ) )
                {
                   if( !cflags[col].test( ColFlag::kLbInf ) &&
                       num.isEq( scaled_val, rowvals[j] ) &&
                       num.isLE( scaled_obj, obj[col] ) &&
                       checkDominance( unbounded_col, col, bestscale, 1 ) )
                      to_lb = true;
                   else if( !cflags[col].test( ColFlag::kUbInf ) &&
                       num.isEq( scaled_val, -rowvals[j] ) &&
                       num.isLE( scaled_obj, -obj[col] ) &&
                       checkDominance( unbounded_col, col, bestscale, -1 ) )
                      to_ub = true;
                }
                else if( rflags[row].test( RowFlag::kLhsInf ) )
                {
                   assert( scaled_val > 0 &&
                           !rflags[row].test( RowFlag::kRhsInf ) );
                   if( !cflags[col].test( ColFlag::kLbInf ) &&
                       num.isLE( scaled_val, rowvals[j] ) &&
                       num.isLE( scaled_obj, obj[col] ) &&
                       checkDominance( unbounded_col, col, bestscale, 1 ) )
                      to_lb = true;
                   else if( !cflags[col].test( ColFlag::kUbInf ) &&
                       num.isLE( scaled_val, -rowvals[j] ) &&
                       num.isLE( scaled_obj, -obj[col] ) &&
                       checkDominance( unbounded_col, col, bestscale, -1 ) )
                      to_ub = true;
                }
                else
                {
                   assert( scaled_val < 0 &&
                           rflags[row].test( RowFlag::kRhsInf ) );
                   if( !cflags[col].test( ColFlag::kLbInf ) &&
                       num.isGE( scaled_val, rowvals[j] ) &&
                       num.isLE( scaled_obj, obj[col] ) &&
                       checkDominance( unbounded_col, col, bestscale, 1 ) )
                      to_lb = true;
                   else if( !cflags[col].test( ColFlag::kUbInf ) &&
                       num.isGE( scaled_val, -rowvals[j] ) &&
                       num.isLE( scaled_obj, -obj[col] ) &&
                       checkDominance( unbounded_col, col, bestscale, -1 ) )
                      to_ub = true;
                }

                if( to_lb || to_ub )
                {
                   domcolsbuffers[k - base].emplace_back( DomcolReduction{ unbounded_col, col,
                         bestrowlock, to_lb ? BoundChange::kUpper : BoundChange::kLower } );
                }
             }
          }
#ifdef PAPILO_TBB
       } );
#endif

         for( int i = start - base; i < ndomcolsbuffers; ++i )
            ndomcols += domcolsbuffers[i].size();

         start = stopp;
      }

      bool lock = false;

      // filter dominations avoiding cyclic conflicts
      do
      {
         for( int i = 0; i < ndomcolsbuffers; ++i )
         {
            if( domcolsbuffers[i].empty() || ( !lock && domcolsbuffers[i][0].implrowlock != -1 ) )
               continue;
            for( int j = 0; j < (int)domcolsbuffers[i].size(); ++j )
            {
               int source = domcolsbuffers[i][j].col2;
               if( domcol[source] != -1 )
                  continue;
               int sink = domcolsbuffers[i][j].col1;
               int node = sink;
               while( domcol[node] >= 0 )
                  node = domcols[domcol[node]].col1;
               if( node == source )
                  continue;
               domcol[source] = domcols.size();
               domcols.emplace_back(std::move(domcolsbuffers[i][j]));
               if( nchildren[sink] >= 1 || domcol[sink] <= -1 )
               {
                  if( nchildren[source] == 0 )
                  {
                     nchildren[source] = -leaves.size();
                     leaves.push_back(source);
                  }
                  ++nchildren[sink];
               }
               else
               {
                  if( nchildren[source] == 0 )
                  {
                     nchildren[source] = nchildren[sink];
                     leaves[-nchildren[source]] = source;
                  }
                  nchildren[sink] = 1;
               }
            }
            domcolsbuffers[i].clear();
         }
         lock = !lock;
      }
      while(lock);

      ndomcols = domcols.size();
   }

   domcolsbuffers.clear();
   domcolsbuffers.shrink_to_fit();
   domcols.shrink_to_fit();
   leaves.shrink_to_fit();

   // add reductions in topological order
   for( auto& node : leaves )
   {
      while( nchildren[node] <= 0 && domcol[node] >= 0 )
      {
         const DomcolReduction& reduction = domcols[domcol[node]];
         domcol[node] = -1;
         node = reduction.col1;
         --nchildren[node];
         TransactionGuard<REAL> tg{ reductions };
         reductions.lockCol( reduction.col1 );
         reductions.lockColBounds( reduction.col1 );
         reductions.lockCol( reduction.col2 );
         reductions.lockColBounds( reduction.col2 );
         if( reduction.implrowlock >= 0 )
            reductions.lockRow( reduction.implrowlock );
         // upper bound is changed to lower bound
         if( reduction.boundchg == BoundChange::kUpper )
         {
            reductions.dominance(reduction.col2, reduction.col1);
            reductions.fixCol( reduction.col2, lbValues[reduction.col2], reduction.implrowlock );
         }
         // lower bound is changed to upper bound
         else
         {
            reductions.dominance(reduction.col1, reduction.col2);
            reductions.fixCol( reduction.col2, ubValues[reduction.col2], reduction.implrowlock );
         }
      }
   }

   return domcols.empty() ? PresolveStatus::kUnchanged : PresolveStatus::kReduced;
}


} // namespace papilo

#endif
