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

#ifndef _PAPILO_PRESOLVERS_PROBING_HPP_
#define _PAPILO_PRESOLVERS_PROBING_HPP_

#include "papilo/core/PresolveMethod.hpp"
#include "papilo/core/ProbingView.hpp"
#include "papilo/core/Problem.hpp"
#include "papilo/core/ProblemUpdate.hpp"
#include "papilo/core/SingleRow.hpp"
#include "papilo/misc/Array.hpp"
#include "papilo/misc/Hash.hpp"
#include "papilo/misc/Vec.hpp"
#include "papilo/misc/compress_vector.hpp"
#include "papilo/misc/fmt.hpp"
#include "papilo/misc/tbb.hpp"
#include <atomic>
#include <boost/functional/hash.hpp>

namespace papilo
{

template <typename REAL>
class Probing : public PresolveMethod<REAL>
{
   Vec<int> nprobed;
   int maxinitialbadgesize = 1000;
   int minbadgesize = 10;
   double mincontdomred = 0.3;

 public:
   Probing() : PresolveMethod<REAL>()
   {
      this->setName( "probing" );
      this->setTiming( PresolverTiming::kExhaustive );
      this->setType( PresolverType::kIntegralCols );
   }

   void
   compress( const Vec<int>& rowmap, const Vec<int>& colmap ) override
   {
      assert( colmap.size() == nprobed.size() );
      compress_vector( colmap, nprobed );
      Message::debug( this,
                      "compress was called, compressed nprobed vector from "
                      "size {} to size {}\n",
                      colmap.size(), nprobed.size() );
   }

   bool
   initialize( const Problem<REAL>& problem,
               const PresolveOptions& presolveOptions ) override
   {
      nprobed.clear();
      nprobed.resize( problem.getNCols(), 0 );

      Message::debug( this, "initialized nprobed vector to size {}\n",
                      nprobed.size() );

      return true;
   }

   void
   addPresolverParams( ParameterSet& paramSet ) override
   {
      paramSet.addParameter( "probing.maxinitialbadgesize",
                             "maximum number of probing candidates probed in "
                             "the first badge of candidates",
                             maxinitialbadgesize, 1 );

      paramSet.addParameter( "probing.minbadgesize",
                             "minimum number of probing candidates probed in "
                             "a single badge of candidates",
                             minbadgesize, 1 );

      paramSet.addParameter(
          "probing.mincontdomred",
          "minimum fraction of domain that needs to be reduced for continuous "
          "variables to accept a bound change in probing",
          mincontdomred, 0.0, 1.0 );
   }

   virtual PresolveStatus
   execute( const Problem<REAL>& problem,
            const ProblemUpdate<REAL>& problemUpdate, const Num<REAL>& num,
            Reductions<REAL>& reductions ) override;
};

#ifdef PAPILO_USE_EXTERN_TEMPLATES
extern template class Probing<double>;
extern template class Probing<Quad>;
extern template class Probing<Rational>;
#endif

template <typename REAL>
PresolveStatus
Probing<REAL>::execute( const Problem<REAL>& problem,
                        const ProblemUpdate<REAL>& problemUpdate,
                        const Num<REAL>& num, Reductions<REAL>& reductions )
{
   if( problem.getNumIntegralCols() == 0 )
      return PresolveStatus::kUnchanged;

   const auto& domains = problem.getVariableDomains();
   const auto& lower_bounds = domains.lower_bounds;
   const auto& upper_bounds = domains.upper_bounds;
   const auto& cflags = domains.flags;

   const auto& consMatrix = problem.getConstraintMatrix();
   const auto& lhs = consMatrix.getLeftHandSides();
   const auto& rhs = consMatrix.getRightHandSides();
   const auto& rflags = consMatrix.getRowFlags();
   const auto& activities = problem.getRowActivities();
   const int ncols = problem.getNCols();
   const auto& colsize = consMatrix.getColSizes();
   const auto& colperm = problemUpdate.getRandomColPerm();

   Vec<int> probing_cands;
   probing_cands.reserve( ncols );

   for( int i = 0; i != ncols; ++i )
   {
      if( !cflags[i].test( ColFlag::kUnbounded ) &&
          cflags[i].test( ColFlag::kIntegral ) && colsize[i] > 0 &&
          lower_bounds[i] == 0 && upper_bounds[i] == 1 )
         probing_cands.push_back( i );
   }

   if( probing_cands.size() == 0 )
      return PresolveStatus::kUnchanged;

   Array<std::atomic_int> probing_scores( ncols );

   for( int i = 0; i != ncols; ++i )
      probing_scores[i].store( 0, std::memory_order_relaxed );

   if( nprobed.size() == 0 )
   {
      nprobed.resize( size_t( ncols ), 0 );

      assert( static_cast<int>( nprobed.size() ) == ncols );
      assert( std::all_of( nprobed.begin(), nprobed.end(),
                           []( int n ) { return n == 0; } ) );
   }

   tbb::parallel_for(
       tbb::blocked_range<int>( 0, problem.getNRows() ),
       [&]( const tbb::blocked_range<int>& r ) {
          Vec<std::pair<REAL, int>> binvarsRow;
          for( int row = r.begin(); row != r.end(); ++row )
          {
             if( consMatrix.isRowRedundant( row ) )
                continue;

             if( ( activities[row].ninfmin != 0 ||
                   rflags[row].test( RowFlag::kRhsInf ) ) &&
                 ( activities[row].ninfmax != 0 ||
                   rflags[row].test( RowFlag::kLhsInf ) ) )
                continue;

             auto rowvec = consMatrix.getRowCoefficients( row );
             const int* colinds = rowvec.getIndices();
             const REAL* rowvals = rowvec.getValues();
             const int rowlen = rowvec.getLength();

             binvarsRow.reserve( rowlen );

             for( int i = 0; i != rowlen; ++i )
             {
                if( cflags[colinds[i]].test( ColFlag::kIntegral ) &&
                    lower_bounds[colinds[i]] == 0.0 &&
                    upper_bounds[colinds[i]] == 1.0 )
                   binvarsRow.emplace_back( rowvals[i], colinds[i] );
             }

             const int nbinvarsrow = static_cast<int>( binvarsRow.size() );

             if( nbinvarsrow == 0 )
                continue;

             pdqsort( binvarsRow.begin(), binvarsRow.end(),
                      []( const std::pair<REAL, int>& a,
                          const std::pair<REAL, int>& b ) {
                         return abs( a.first ) > abs( b.first );
                      } );

             for( int i = 0; i != nbinvarsrow; ++i )
             {
                int col = binvarsRow[i].second;
                REAL abscoef = abs( binvarsRow[i].first );
                REAL minimplcoef = abscoef;

                if( activities[row].ninfmin == 0 &&
                    !rflags[row].test( RowFlag::kRhsInf ) )
                   minimplcoef = std::min(
                       minimplcoef,
                       REAL( rhs[row] - activities[row].min - abscoef ) );

                if( activities[row].ninfmax == 0 &&
                    !rflags[row].test( RowFlag::kLhsInf ) )
                   minimplcoef =
                       std::min( minimplcoef, REAL( activities[row].max -
                                                    abscoef - lhs[row] ) );

                if( num.isFeasLE( abscoef, minimplcoef ) )
                   break;

                int nimplbins = 0;
                for( int j = i + 1; j != nbinvarsrow; ++j )
                {
                   if( num.isFeasGT( abs( binvarsRow[j].first ), minimplcoef ) )
                      ++nimplbins;
                   else
                      break;
                }

                if( nimplbins != 0 )
                   probing_scores[col].fetch_add( nimplbins,
                                                  std::memory_order_relaxed );
                else
                   break;
             }

             binvarsRow.clear();
          }
       } );

   pdqsort( probing_cands.begin(), probing_cands.end(),
            [this, &probing_scores, &colsize, &colperm]( int col1, int col2 ) {
               std::pair<double, double> s1;
               std::pair<double, double> s2;
               if( nprobed[col2] == 0 && probing_scores[col2] > 0 )
                  s2.first = probing_scores[col2] /
                             static_cast<double>( colsize[col2] );
               else
                  s2.first = 0;
               if( nprobed[col1] == 0 && probing_scores[col1] > 0 )
                  s1.first = probing_scores[col1] /
                             static_cast<double>( colsize[col1] );
               else
                  s1.first = 0;

               s1.second =
                   ( probing_scores[col1].load( std::memory_order_relaxed ) /
                     static_cast<double>( 1 + nprobed[col1] * colsize[col1] ) );
               s2.second =
                   ( probing_scores[col2].load( std::memory_order_relaxed ) /
                     static_cast<double>( 1 + nprobed[col2] * colsize[col2] ) );
               return s1 > s2 || ( s1 == s2 && colperm[col1] < colperm[col2] );
            } );

   const auto& rowsize = consMatrix.getRowSizes();

   HashMap<std::pair<int, int>, int, boost::hash<std::pair<int, int>>>
       substitutionsPos;
   Vec<ProbingSubstitution<REAL>> substitutions;
   Vec<int> boundPos( size_t( 2 * ncols ), 0 );
   Vec<ProbingBoundChg<REAL>> boundChanges;
   boundChanges.reserve( ncols );

   std::atomic_bool infeasible{ false };

   int currentbadgestart = 0;

   int64_t workinglimit = consMatrix.getNnz() * 2;

   const int nprobingcands = static_cast<int>( probing_cands.size() );
   int badgesize = 0;
   int initialbadgelim = 0.1 * workinglimit;
   for( int i : probing_cands )
   {
      ++badgesize;

      if( badgesize == maxinitialbadgesize )
         break;

      initialbadgelim -= colsize[i];
      if( initialbadgelim <= 0 )
         break;

      auto colvec = consMatrix.getColumnCoefficients( i );
      const int* rowinds = colvec.getIndices();
      for( int k = 0; k != colvec.getLength(); ++k )
      {
         initialbadgelim -= ( rowsize[rowinds[k]] - 1 );

         if( initialbadgelim <= 0 )
            break;
      }

      if( initialbadgelim <= 0 )
         break;
   }

   badgesize = std::max( std::min( nprobingcands, minbadgesize ), badgesize );

   int currentbadgeend = currentbadgestart + badgesize;
   int nuseless = 0;
   bool abort = false;

   // use tbb combinable so that each thread will copy the activities and
   // bounds at most once
   tbb::combinable<ProbingView<REAL>> probing_views( [this, &problem, &num]() {
      ProbingView<REAL> probingView( problem, num );
      probingView.setMinContDomRed( mincontdomred );
      return probingView;
   } );

   do
   {
      Message::debug( this, "probing candidates {} to {}\n", currentbadgestart,
                      currentbadgeend );

      tbb::parallel_for(
          tbb::blocked_range<int>( currentbadgestart, currentbadgeend ),
          [&]( const tbb::blocked_range<int>& r ) {
             ProbingView<REAL>& probingView = probing_views.local();

             for( int i = r.begin(); i != r.end(); ++i )
             {
                const int col = probing_cands[i];

                assert( cflags[col].test( ColFlag::kIntegral ) &&
                            lower_bounds[col] == 0 ||
                        upper_bounds[col] == 1 );

                if( infeasible.load( std::memory_order_relaxed ) )
                   break;

                assert( !probingView.isInfeasible() );
                probingView.setProbingColumn( col, true );
                probingView.propagateDomains();
                probingView.storeImplications();
                probingView.reset();

                if( infeasible.load( std::memory_order_relaxed ) )
                   break;

                assert( !probingView.isInfeasible() );
                probingView.setProbingColumn( col, false );
                probingView.propagateDomains();

                bool globalInfeasible = probingView.analyzeImplications();
                probingView.reset();

                ++nprobed[col];

                if( globalInfeasible )
                {
                   infeasible.store( true, std::memory_order_relaxed );
                   break;
                }
             }
          } );

      if( infeasible.load( std::memory_order_relaxed ) )
         return PresolveStatus::kInfeasible;

      int64_t amountofwork = 0;
      int nfixings = 0;
      int nboundchgs = 0;
      int nsubstitutions = -substitutions.size();

      probing_views.combine_each( [&]( ProbingView<REAL>& probingView ) {
         const auto& probingBoundChgs = probingView.getProbingBoundChanges();
         const auto& probingSubstitutions =
             probingView.getProbingSubstitutions();

         amountofwork += probingView.getAmountOfWork();

         for( const ProbingSubstitution<REAL>& subst : probingSubstitutions )
         {
            auto insres = substitutionsPos.emplace(
                std::make_pair( subst.col1, subst.col2 ),
                substitutions.size() );

            if( insres.second )
               substitutions.push_back( subst );
         }

         for( const ProbingBoundChg<REAL>& boundChg : probingBoundChgs )
         {
            if( boundPos[2 * boundChg.col + boundChg.upper] == 0 )
            {
               // found new bound change
               boundChanges.emplace_back( boundChg );
               boundPos[2 * boundChg.col + boundChg.upper] =
                   boundChanges.size();

               // check if column is now fixed
               if( ( boundChg.upper &&
                     boundChg.bound == lower_bounds[boundChg.col] ) ||
                   ( !boundChg.upper &&
                     boundChg.bound == upper_bounds[boundChg.col] ) )
                  ++nfixings;
               else
                  ++nboundchgs;
            }
            else
            {
               // already changed that bound
               ProbingBoundChg<REAL>& otherBoundChg =
                   boundChanges[boundPos[2 * boundChg.col + boundChg.upper] -
                                1];

               if( boundChg.upper && boundChg.bound < otherBoundChg.bound )
               {
                  // new upper bound change is tighter
                  otherBoundChg.bound = boundChg.bound;

                  // check if column is now fixed
                  if( boundChg.bound == lower_bounds[boundChg.col] )
                     ++nfixings;
               }
               else if( !boundChg.upper &&
                        boundChg.bound > otherBoundChg.bound )
               {
                  // new lower bound change is tighter
                  otherBoundChg.bound = boundChg.bound;

                  // check if column is now fixed
                  if( boundChg.bound == upper_bounds[boundChg.col] )
                     ++nfixings;
               }

               // do only count fixings in this case for two reasons:
               // 1) the number of bound changes depends on the order and
               // would make probing non deterministic 2) the boundchange was
               // already counted in previous rounds and will only be added
               // once
            }
         }

         probingView.clearResults();
      } );

      nsubstitutions += substitutions.size();
      currentbadgestart = currentbadgeend;

      if( nfixings == 0 && nboundchgs == 0 && nsubstitutions == 0 )
         nuseless += amountofwork;
      else
         nuseless = 0;

      Message::debug(
          this,
          "probing found: {} fixings, {} substitutions, {} bound changes\n",
          nfixings, nsubstitutions, nboundchgs );

      int64_t extrawork =
          ( ( 0.1 * ( nfixings + nsubstitutions ) + 0.01 * nboundchgs ) *
            consMatrix.getNnz() );

      workinglimit -= amountofwork;
      workinglimit += extrawork;

      badgesize = static_cast<int>(
          ceil( badgesize * static_cast<double>( workinglimit + extrawork ) /
                amountofwork ) );
      badgesize = std::max( badgesize, minbadgesize );
      badgesize = std::min( nprobingcands - currentbadgestart, badgesize );
      currentbadgeend = currentbadgestart + badgesize;

      abort = nuseless >= consMatrix.getNnz() * 2 || workinglimit < 0 ||
              currentbadgestart == currentbadgeend;
   } while( !abort );

   PresolveStatus result = PresolveStatus::kUnchanged;

   if( !boundChanges.empty() )
   {
      pdqsort(
          boundChanges.begin(), boundChanges.end(),
          []( const ProbingBoundChg<REAL>& a, const ProbingBoundChg<REAL>& b ) {
             return ( a.col << 1 | a.upper ) < ( b.col << 1 | b.upper );
          } );

      for( const ProbingBoundChg<REAL>& boundChg : boundChanges )
      {
         if( boundChg.upper )
            reductions.changeColUB( boundChg.col, boundChg.bound );
         else
            reductions.changeColLB( boundChg.col, boundChg.bound );
      }

      result = PresolveStatus::kReduced;
   }

   if( !substitutions.empty() )
   {
      pdqsort( substitutions.begin(), substitutions.end(),
               []( const ProbingSubstitution<REAL>& a,
                   const ProbingSubstitution<REAL>& b ) {
                  return std::make_pair( a.col1, a.col2 ) >
                         std::make_pair( b.col1, b.col2 );
               } );

      int lastsubstcol = -1;

      for( const ProbingSubstitution<REAL>& subst : substitutions )
      {
         if( subst.col1 == lastsubstcol )
            continue;

         lastsubstcol = subst.col1;

         reductions.replaceCol( subst.col1, subst.col2, subst.col2scale,
                                subst.col2const );
      }

      result = PresolveStatus::kReduced;
   }

   return result;
}

} // namespace papilo

#endif
