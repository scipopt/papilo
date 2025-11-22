/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*               This file is part of the program and library                */
/*    PaPILO --- Parallel Presolve for Integer and Linear Optimization       */
/*                                                                           */
/* Copyright (C) 2020-2025 Zuse Institute Berlin (ZIB)                       */
/*                                                                           */
/* Licensed under the Apache License, Version 2.0 (the "License");           */
/* you may not use this file except in compliance with the License.          */
/* You may obtain a copy of the License at                                   */
/*                                                                           */
/*     http://www.apache.org/licenses/LICENSE-2.0                            */
/*                                                                           */
/* Unless required by applicable law or agreed to in writing, software       */
/* distributed under the License is distributed on an "AS IS" BASIS,         */
/* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  */
/* See the License for the specific language governing permissions and       */
/* limitations under the License.                                            */
/*                                                                           */
/* You should have received a copy of the Apache-2.0 license                 */
/* along with PaPILO; see the file LICENSE. If not visit scipopt.org.        */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef _PAPILO_PRESOLVERS_PROBING_HPP_
#define _PAPILO_PRESOLVERS_PROBING_HPP_

#include "papilo/Config.hpp"
#include "papilo/core/CliqueProbingView.hpp"
#include "papilo/core/PresolveMethod.hpp"
#include "papilo/core/ProbingView.hpp"
#include "papilo/core/Problem.hpp"
#include "papilo/core/ProblemUpdate.hpp"
#include "papilo/misc/Array.hpp"
#include "papilo/misc/Hash.hpp"
#include "papilo/misc/Vec.hpp"
#include "papilo/misc/compress_vector.hpp"
#include <algorithm>
#include <atomic>
#include <boost/functional/hash.hpp>
#include <set>
#include <tuple>
#ifdef PAPILO_TBB
#include "papilo/misc/tbb.hpp"
#endif

namespace papilo
{

const static int DEFAULT_MAX_BADGE_SIZE = -1;

template <typename REAL>
class Probing : public PresolveMethod<REAL>
{
   Vec<int> nprobed;
   Vec<int> nprobedcliques;
   int unsuccessfulcliqueprobing = 0;
   int maxinitialbadgesize = 1000;
   int minbadgesize = 10;
   int max_badge_size = DEFAULT_MAX_BADGE_SIZE;
   double mincontdomred = 0.3;
   int maxCliqueLength = 150;
   int max_probed_clique_vars = 3000;
   int ratiocoveredcliquevars = 2;
   int initialbatchsize = 2;
   int cliquereductionfactor = 3;
   int minabortedvariables = 0;
   int numcliquefails = 0;

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
      nprobedcliques.clear();
      unsuccessfulcliqueprobing = 0;
      nprobed.resize( problem.getNCols(), 0 );
      nprobedcliques.resize( problem.getNRows(), 0 );

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

      paramSet.addParameter( "probing.maxbadgesize",
                             "maximal number of probing candidates probed in "
                             "a single badge of candidates (-1, 0: unlimited)",
                             max_badge_size, DEFAULT_MAX_BADGE_SIZE );

      paramSet.addParameter(
          "probing.mincontdomred",
          "minimum fraction of domain that needs to be reduced for continuous "
          "variables to accept a bound change in probing",
          mincontdomred, 0.0, 1.0 );

      paramSet.addParameter( "probing.maxCliqueLength",
                             "maximal size of cliques that are probed on ",
                             maxCliqueLength, 150 );

      paramSet.addParameter( "probing.max_probed_clique_vars",
                             "maximal number of variables in cliques that are probed on ",
                             max_probed_clique_vars, 3000 );

      paramSet.addParameter( "probing.ratiocoveredcliquevars",
                             "maximal ratio of variables in clique that are also covered by other probed cliques ",
                             ratiocoveredcliquevars, 2 );

      paramSet.addParameter( "probing.initialbatchsize",
                             "initial number of probed cliques ",
                             initialbatchsize, 2 );

      paramSet.addParameter( "probing.cliquereductionfactor",
                             "number of reductions per variable in probed clique to deem clique probing successful ",
                             cliquereductionfactor, 3 );

      paramSet.addParameter( "probing.minabortedvariables",
                             "minimum number of variables left in clique before we consider abort probing on the clique",
                             minabortedvariables, 0 );

      paramSet.addParameter( "probing.numcliquefails",
                             "number of times clique probing may fail before being disabled",
                             numcliquefails, -1 );
   }

   PresolveStatus
   execute( const Problem<REAL>& problem,
            const ProblemUpdate<REAL>& problemUpdate, const Num<REAL>& num,
            Reductions<REAL>& reductions, const Timer& timer,
            int& reason_of_infeasibility ) override;

   bool
   isBinaryVariable( REAL upper_bound, REAL lower_bound, int column_size,
                     const Flags<ColFlag>& colFlag ) const;

   void
   set_max_badge_size( int val );
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
                        const Num<REAL>& num, Reductions<REAL>& reductions,
                        const Timer& timer, int& reason_of_infeasibility )
{
   //auto initstarttime = timer.getTime();
   if( problem.getNumIntegralCols() == 0 )
      return PresolveStatus::kUnchanged;

   const auto& domains = problem.getVariableDomains();
   const Vec<REAL>& lower_bounds = domains.lower_bounds;
   const Vec<REAL>& upper_bounds = domains.upper_bounds;
   const Vec<ColFlags>& cflags = domains.flags;

   const auto& consMatrix = problem.getConstraintMatrix();
   const auto& lhs = consMatrix.getLeftHandSides();
   const auto& rhs = consMatrix.getRightHandSides();
   const Vec<RowFlags>& rowFlags = consMatrix.getRowFlags();
   const auto& activities = problem.getRowActivities();
   const int ncols = problem.getNCols();
   const Vec<int>& colsize = consMatrix.getColSizes();
   const auto& colperm = problemUpdate.getRandomColPerm();

   const int nrows = problem.getNRows();
   Vec<std::pair<int,std::pair<int,bool>>> cliques;
   cliques.reserve( nrows );
   Vec<int> probing_cands;
   probing_cands.reserve( ncols );
   std::cout<<"\nNumber of unsuccessful clique probing attempts:" << unsuccessfulcliqueprobing;

   if( unsuccessfulcliqueprobing <= numcliquefails )
   {
      //auto cliquefindstarttime = timer.getTime();
      for( int row = 0; row != nrows; ++row )
      {
         assert( row >= 0 && row < nrows );
         auto rowvec = consMatrix.getRowCoefficients( row );
         auto cliquecheck = problem.is_clique_and_equation( consMatrix, row, num );
         if( cliquecheck.first && rowvec.getLength() < maxCliqueLength )
         {
            cliques.emplace_back( row, std::make_pair(0, cliquecheck.second) );
         }
      }
   }

   for( int i = 0; i != ncols; ++i )
   {
      if( isBinaryVariable( upper_bounds[i], lower_bounds[i], colsize[i],
                            cflags[i] ) )
         probing_cands.push_back( i );
   }


   if( probing_cands.empty() )
      return PresolveStatus::kUnchanged;

   Array<std::atomic_int> probing_scores( ncols );

   for( int i = 0; i != ncols; ++i )
      probing_scores[i].store( 0, std::memory_order_relaxed );

   if( nprobed.empty() )
   {
      nprobed.resize( size_t( ncols ), 0 );

      assert( static_cast<int>( nprobed.size() ) == ncols );
      assert( std::all_of( nprobed.begin(), nprobed.end(),
                           []( int n ) { return n == 0; } ) );
   }

   if( unsuccessfulcliqueprobing <= numcliquefails )
   {
      if( nprobedcliques.empty() )
      {
         nprobedcliques.resize( size_t( nrows ), 0 );
      }
   }

#ifdef PAPILO_TBB
   tbb::parallel_for(
       tbb::blocked_range<int>( 0, problem.getNRows() ),
       [&]( const tbb::blocked_range<int>& r )
       {
          Vec<std::pair<REAL, int>> binary_variables_in_row;
          for( int row = r.begin(); row != r.end(); ++row )
#else
   Vec<std::pair<REAL, int>> binary_variables_in_row;
   for( int row = 0; row != problem.getNRows(); ++row )
#endif
          {
             if( consMatrix.isRowRedundant( row ) )
                continue;

             if( ( activities[row].ninfmin != 0 ||
                   rowFlags[row].test( RowFlag::kRhsInf ) ) &&
                 ( activities[row].ninfmax != 0 ||
                   rowFlags[row].test( RowFlag::kLhsInf ) ) )
                continue;
             assert( row >= 0 && row < nrows );
             auto rowvec = consMatrix.getRowCoefficients( row );
             const int* colinds = rowvec.getIndices();
             const REAL* rowvals = rowvec.getValues();
             const int rowlen = rowvec.getLength();

             binary_variables_in_row.reserve( rowlen );

             for( int i = 0; i != rowlen; ++i )
             {
                if( isBinaryVariable( upper_bounds[i], lower_bounds[i],
                                      colsize[i], cflags[i] ) )
                   binary_variables_in_row.emplace_back( rowvals[i],
                                                         colinds[i] );
             }

             const int nbinvarsrow =
                 static_cast<int>( binary_variables_in_row.size() );

             if( nbinvarsrow == 0 )
                continue;

             pdqsort( binary_variables_in_row.begin(),
                      binary_variables_in_row.end(),
                      []( const std::pair<REAL, int>& a,
                          const std::pair<REAL, int>& b )
                      { return abs( a.first ) > abs( b.first ); } );

             for( int i = 0; i < nbinvarsrow; ++i )
             {
                // TODO: wouldn't be simpler: calculate minimplcoef, if greater
                // equals than abs, then discard
                int col = binary_variables_in_row[i].second;
                REAL abscoef = abs( binary_variables_in_row[i].first );
                REAL minimplcoef = abscoef;

                if( activities[row].ninfmin == 0 &&
                    !rowFlags[row].test( RowFlag::kRhsInf ) )
                   minimplcoef = std::min(
                       minimplcoef,
                       REAL( rhs[row] - activities[row].min - abscoef ) );

                if( activities[row].ninfmax == 0 &&
                    !rowFlags[row].test( RowFlag::kLhsInf ) )
                   minimplcoef =
                       std::min( minimplcoef, REAL( activities[row].max -
                                                    abscoef - lhs[row] ) );

                if( num.isFeasLE( abscoef, minimplcoef ) )
                   break;

                int nimplbins = 0;
                for( int j = i + 1; j != nbinvarsrow; ++j )
                {
                   if( num.isFeasGT( abs( binary_variables_in_row[j].first ),
                                     minimplcoef ) )
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

             binary_variables_in_row.clear();
          }

#ifdef PAPILO_TBB
       } );
#endif

   Vec<std::pair<int,bool>> probingCliques;

   if( unsuccessfulcliqueprobing <= numcliquefails )
   {
#ifdef PAPILO_TBB
      tbb::parallel_for(
      tbb::blocked_range<int>( 0, cliques.end() - cliques.begin() ),
      [&]( const tbb::blocked_range<int>& r ) {
      for( int clique = r.begin(); clique != r.end(); ++clique )
#else
      for( int clique = 0; clique < cliques.end() - cliques.begin(); ++clique )
#endif
      {
         assert( cliques[clique].first >= 0 && cliques[clique].first < nrows );
         auto rowvec = consMatrix.getRowCoefficients( cliques[clique].first );
         auto rowinds = rowvec.getIndices();
         for( int ind = 0; ind < rowvec.getLength(); ++ind )
         {
            cliques[clique].second.first += probing_scores[rowinds[ind]];
            assert( isBinaryVariable( problem.getUpperBounds()[rowinds[ind]], problem.getLowerBounds()[rowinds[ind]],
            problem.getColSizes()[rowinds[ind]], problem.getColFlags()[rowinds[ind]]  ));
         }
         cliques[clique].second.first = cliques[clique].second.first
         / ( rowvec.getLength() * ( nprobedcliques[cliques[clique].first] + 1 ) ) ;
      }
#ifdef PAPILO_TBB
      } );
#endif

      pdqsort( cliques.begin(), cliques.end(),
      []( const std::pair<int,std::pair<int,bool>>& clique1, const std::pair<int,std::pair<int,bool>>& clique2 )
      {
         return clique1.second.first > clique2.second.first ;
      } );

      int cliquevars = 0;
      Vec<bool> probedCliqueVars(ncols, false);
      probingCliques.reserve( cliques.end() - cliques.begin() );
      for( int clique = 0; clique < static_cast<int>(cliques.end() - cliques.begin()); ++clique )
      {
         assert( cliques[clique].first >= 0 && cliques[clique].first < nrows );
         auto rowvec = consMatrix.getRowCoefficients( cliques[clique].first );
         auto rowinds = rowvec.getIndices();
         int covered = 0;
         for( int ind = 0; ind < rowvec.getLength() && ratiocoveredcliquevars * covered <= rowvec.getLength(); ++ind )
         {
            if( probedCliqueVars[rowinds[ind]] )
               covered += 1;
         }
         if( ratiocoveredcliquevars * covered <= rowvec.getLength() )
         {
            assert( cliques[clique].first >= 0 && cliques[clique].first < nrows );
            probingCliques.emplace_back( cliques[clique].first, cliques[clique].second.second );
            for( int ind = 0; ind < static_cast<int>(rowvec.getLength()); ++ind )
            {
               if( !probedCliqueVars[rowinds[ind]] )
               {
                  probedCliqueVars[rowinds[ind]] = true;
                  cliquevars += 1;
               }
            }
         }
         if( cliquevars > max_probed_clique_vars)
         {
            break;
         }
      }
   }

   std::set<int> probedvars;
   int nprobedvars = 0;
   std::atomic_bool infeasible{ false };
   std::atomic_int infeasible_variable{ -1 };
   int batchend = 0;
   int totalnumpropagations = 0;
   auto cliqueprobingtime = timer.getTime();
   Vec<CliqueProbingSubstitution<REAL>> cliquesubstitutions;
   Vec<CliqueProbingBoundChg<REAL>> cliqueBoundChanges;
   int ncliquefixings = 0;
   int ncliqueboundchgs = 0;
   int ncliquesubstitutions = 0;
   int64_t amountofwork = 0;
#ifdef PAPILO_TBB
   Vec<int> change_to_equation_comb;
#endif

      HashMap<std::pair<int, int>, int, boost::hash<std::pair<int, int>>>
      cliqueSubstitutionsPos;
      Vec<int> cliqueBoundPos( size_t( 2 * ncols ), 0 );
      cliqueBoundChanges.reserve( ncols );

#ifdef PAPILO_TBB

      tbb::combinable<Vec<CliqueProbingBoundChg<REAL>>> clique_probing_bound_changes;
      tbb::combinable<Vec<CliqueProbingSubstitution<REAL>>> clique_probing_subs;
      tbb::combinable<Vec<int>> change_to_equation;
      tbb::combinable<int> amounts_of_work;
#else
      Vec<CliqueProbingBoundChg<REAL>> clique_probing_bound_changes;
      Vec<CliqueProbingSubstitution<REAL>> clique_probing_subs;
      Vec<int> change_to_equation;
#endif
   Vec<int> finalinds( static_cast<int>(probingCliques.end() - probingCliques.begin()), 0 );

   bool earlycliqueabort = false;

   if( unsuccessfulcliqueprobing <= numcliquefails )
   {
      auto propagate_cliques = [&]( int cliquestart, int cliqueend )
      {
#ifdef PAPILO_TBB
         tbb::enumerable_thread_specific<int> numpropagations(0);
         tbb::parallel_for( tbb::blocked_range<int>( cliquestart, cliqueend ),
         [&]( const tbb::blocked_range<int>& r )
         {
            for( int i = r.begin(); i < r.end(); ++i )
#else
            int numpropagations = 0;
            for( int i = cliquestart; i < cliqueend; ++i )
#endif
            {
               int clique = probingCliques[i].first;
               assert( clique >= 0 && clique < nrows );
               auto cliquevec = consMatrix.getRowCoefficients( clique );
               auto cliqueind = cliquevec.getIndices();
               auto cliquelen = cliquevec.getLength();
               for( int j = 0; j < cliquelen; ++j )
               {
                  nprobed[cliqueind[j]] +=1;
               }
               nprobedcliques[clique] += 1;
               CliqueProbingView<REAL> local_clique_probing_view( problem, num );
               local_clique_probing_view.setMinContDomRed( mincontdomred );


               if( infeasible.load( std::memory_order_relaxed ) )
                  break;

               std::tuple<bool,bool,int> cliqueProbingResult = local_clique_probing_view.probeClique(clique, cliqueind, cliquelen,
                  probing_cands, probingCliques[i].second, probing_scores, colsize, colperm, nprobed, cliquereductionfactor,
                  minabortedvariables );


               bool globalInfeasible = std::get<0>(cliqueProbingResult);
               finalinds[i] = std::get<2>(cliqueProbingResult);
               if( std::get<1>(cliqueProbingResult) && !probingCliques[i].second )
               {
#ifdef PAPILO_TBB
                  change_to_equation.local().emplace_back( clique );
#else
                  change_to_equation.emplace_back( clique );
#endif
               }
               if( !globalInfeasible )
                  globalInfeasible = local_clique_probing_view.analyzeImplications();
               if( globalInfeasible )
                     {
                        local_clique_probing_view.resetClique();
                        infeasible.store( true, std::memory_order_relaxed );
                        infeasible_variable.store( cliqueind[0] );
                        break;
                     }
#ifdef PAPILO_TBB
               numpropagations.local() += local_clique_probing_view.getNumPropagations();
               amounts_of_work.local() += local_clique_probing_view.getAmountOfWork();
#else
               numpropagations += local_clique_probing_view.getNumPropagations();
               amountofwork += local_clique_probing_view.getAmountOfWork();
#endif
               local_clique_probing_view.clearAmountOfWork();
               local_clique_probing_view.resetClique();
               Vec<CliqueProbingBoundChg<REAL>> local_bound_chgs = local_clique_probing_view.getProbingBoundChanges();
               Vec<CliqueProbingSubstitution<REAL>> local_subs = local_clique_probing_view.getProbingSubstitutions();
#ifdef PAPILO_TBB
               clique_probing_bound_changes.local().insert(clique_probing_bound_changes.local().end(),
                  local_bound_chgs.begin(), local_bound_chgs.end() );
               clique_probing_subs.local().insert(clique_probing_subs.local().end(),
                  local_subs.begin(), local_subs.end() );
#else
               clique_probing_bound_changes.insert(clique_probing_bound_changes.end(),
                  local_bound_chgs.begin(), local_bound_chgs.end() );
               clique_probing_subs.insert(clique_probing_subs.end(),
                  local_subs.begin(), local_subs.end() );
#endif
            }
#ifdef PAPILO_TBB
         } );
      totalnumpropagations += std::accumulate(numpropagations.begin(), numpropagations.end(), 0);
#else
      totalnumpropagations += numpropagations;
#endif
      };
      auto cliqueprobinstarttime = timer.getTime();

      int batchsize = initialbatchsize;
      int batchstart = 0;
      batchend = std::min(batchstart + batchsize, static_cast<int>(probingCliques.end() - probingCliques.begin()));
      bool successlasttime = true;

      while( batchstart < static_cast<int>(probingCliques.end() - probingCliques.begin()) )
      {
         propagate_cliques( batchstart, std::min(batchend, static_cast<int>(probingCliques.end() - probingCliques.begin())) );

#ifdef PAPILO_TBB
         int numcliquereductions = 0;
         int numcliquebc = 0;
         int numcliquesubs = 0;
         clique_probing_bound_changes.combine_each([&numcliquebc](
            Vec<CliqueProbingBoundChg<REAL>>& clique_probing_bound_changes_local) {
            numcliquebc += static_cast<int>(clique_probing_bound_changes_local.size());
         });
         clique_probing_subs.combine_each([&numcliquesubs](
            Vec<CliqueProbingSubstitution<REAL>>& clique_probing_substitutions_local) {
            numcliquesubs += static_cast<int>(clique_probing_substitutions_local.size());
         });
         amounts_of_work.combine_each([&amountofwork]( int work )
         {
            amountofwork += work;
         });
         amounts_of_work.clear();
         numcliquereductions = numcliquebc + numcliquesubs;
         if( infeasible || numcliquereductions
            <= totalnumpropagations * cliquereductionfactor )
#else
         if( infeasible || (static_cast<int>(clique_probing_bound_changes.size()) + static_cast<int>(clique_probing_subs.size()))
            <= totalnumpropagations * cliquereductionfactor )
#endif
         {
            if( !successlasttime )
               break;
            else
            {
               successlasttime = false;
               batchstart = batchend;
               batchsize /= 2;
               batchend = batchstart + batchsize;
#ifdef PAPILO_TBB
               if( totalnumpropagations * static_cast<double>( consMatrix.getNnz() * 2 + ( ( 0.1 * 
                  ( numcliquebc + numcliquesubs ) 
                  + 0.01 * numcliquebc ) *
                  consMatrix.getNnz() ) ) / amountofwork < 0.1 * totalnumpropagations )
#else 
               if( totalnumpropagations * static_cast<double>( consMatrix.getNnz() * 2 + ( ( 0.1 * 
                  ( static_cast<int>(clique_probing_bound_changes.size()) + static_cast<int>(clique_probing_subs.size()) ) 
                  + 0.01 * static_cast<int>(clique_probing_bound_changes.size()) ) *
                  consMatrix.getNnz() ) ) / amountofwork < 0.1 * totalnumpropagations )
#endif
               {
                  earlycliqueabort = true;
                  break;
               }   
            }
         }
         else
         {
            successlasttime = true;
            batchstart = batchend;
            batchsize *= 2;
            batchend = batchstart + batchsize;
         }
      }
      if( infeasible )
      {
         reason_of_infeasibility =
             infeasible_variable.load( std::memory_order_relaxed );

         return PresolveStatus::kInfeasible;
      }
      cliqueprobingtime = timer.getTime() - cliqueprobinstarttime;
      for( int clique = 0; clique < std::min(batchend, static_cast<int>(probingCliques.end() - probingCliques.begin())); ++clique )
      {
         auto cliquevec = consMatrix.getRowCoefficients( probingCliques[clique].first );
         auto cliqueind = cliquevec.getIndices();
         for( int ind = 0; ind < finalinds[clique]; ++ind )
         {
            probing_scores[cliqueind[ind]] = -100000;
            probedvars.emplace(cliqueind[ind]);
         }
         nprobedvars += finalinds[clique];
      }

      ncliquesubstitutions = -cliquesubstitutions.size();

#ifdef PAPILO_TBB
      change_to_equation_comb = change_to_equation.combine(
         [](const Vec<int>& a, const Vec<int>& b) {
         Vec<int> result = a;
         result.insert(result.end(), b.begin(), b.end() );
         return result;
      } );
      pdqsort( change_to_equation_comb.begin(), change_to_equation_comb.end() );
      for( int i = 0; i < static_cast<int>(change_to_equation_comb.size()); ++i )
      {
         assert( change_to_equation_comb[i] >= 0 && change_to_equation_comb[i] < nrows );
         auto cliquevec = consMatrix.getRowCoefficients( change_to_equation_comb[i] );
         auto cliquelen = cliquevec.getLength();
         auto vals = cliquevec.getValues();
         auto maxcoeff = vals[0];
         auto mincoeff = vals[0];
         for( int j = 1; j < cliquelen; ++j )
         {
            if( num.isGT(vals[j], maxcoeff) )
               maxcoeff = vals[j];
            else if( num.isLT(vals[j], mincoeff) )
               mincoeff = vals[j];
         }
         if( problem.is_rhs_clique( consMatrix, change_to_equation_comb[i], num )
         && num.isLT(consMatrix.getLeftHandSides()[change_to_equation_comb[i]], mincoeff) )
         {
            assert( mincoeff != consMatrix.getLeftHandSides()[change_to_equation_comb[i]]);
            reductions.changeRowLHS( change_to_equation_comb[i],  mincoeff );
         }
         else if( num.isGT(consMatrix.getRightHandSides()[change_to_equation_comb[i]], maxcoeff) )
         {
            assert( maxcoeff != consMatrix.getRightHandSides()[change_to_equation_comb[i]]);
            reductions.changeRowRHS( change_to_equation_comb[i],  maxcoeff );
         }
      }
#else
      pdqsort( change_to_equation.begin(), change_to_equation.end() );
      for( int i = 0; i < static_cast<int>(change_to_equation.size()); ++i )
      {
         assert( change_to_equation[i] >= 0 && change_to_equation[i] < nrows );
         auto cliquevec = consMatrix.getRowCoefficients( change_to_equation[i] );
         auto cliquelen = cliquevec.getLength();
         auto vals = cliquevec.getValues();
         auto maxcoeff = vals[0];
         auto mincoeff = vals[0];
         for( int j = 1; j < cliquelen; ++j )
         {
            if( num.isGT(vals[j], maxcoeff) )
               maxcoeff = vals[j];
            else if( num.isLT(vals[j], mincoeff) )
               mincoeff = vals[j];
         }
         if( problem.is_rhs_clique( consMatrix, change_to_equation[i], num )
         && num.isLT(consMatrix.getLeftHandSides()[change_to_equation[i]], mincoeff) )
         {
            assert( mincoeff != consMatrix.getLeftHandSides()[change_to_equation[i]]);
            reductions.changeRowLHS( change_to_equation[i],  mincoeff );
         }
         else if( num.isGT(consMatrix.getRightHandSides()[change_to_equation[i]], maxcoeff) )
         {
            assert( maxcoeff != consMatrix.getRightHandSides()[change_to_equation[i]]);
            reductions.changeRowRHS( change_to_equation[i],  maxcoeff );
         }
   }
#endif


#ifdef PAPILO_TBB
      Vec<CliqueProbingBoundChg<REAL>> cliqueProbingBoundChgs;
      clique_probing_bound_changes.combine_each(
      [&cliqueProbingBoundChgs]( Vec<CliqueProbingBoundChg<REAL>>& cliqueProbingBoundChgsLocal )
      {
         cliqueProbingBoundChgs.insert(cliqueProbingBoundChgs.end(),
            cliqueProbingBoundChgsLocal.begin(), cliqueProbingBoundChgsLocal.end());
      });

      Vec<CliqueProbingSubstitution<REAL>> cliqueProbingSubstitutions;
      clique_probing_subs.combine_each(
      [&cliqueProbingSubstitutions]( Vec<CliqueProbingSubstitution<REAL>>& cliqueProbingSubsLocal )
      {
         cliqueProbingSubstitutions.insert(cliqueProbingSubstitutions.end(),
            cliqueProbingSubsLocal.begin(), cliqueProbingSubsLocal.end());
      });
#else
      auto& cliqueProbingBoundChgs = clique_probing_bound_changes;
      auto& cliqueProbingSubstitutions = clique_probing_subs;
#endif

               for( const CliqueProbingSubstitution<REAL>& subst :
                     cliqueProbingSubstitutions )
               {
                  auto insres = cliqueSubstitutionsPos.emplace(
                     std::make_pair( subst.col1, subst.col2 ),
                     cliquesubstitutions.size() );

                  if( insres.second )
                     cliquesubstitutions.push_back( subst );
               }

               for( const CliqueProbingBoundChg<REAL>& boundChg : cliqueProbingBoundChgs )
               {
                  if( cliqueBoundPos[2 * boundChg.col + boundChg.upper] == 0 )
                  {
                     // found new bound change
                     cliqueBoundChanges.emplace_back( boundChg );
                     cliqueBoundPos[2 * boundChg.col + boundChg.upper] =
                        cliqueBoundChanges.size();

                     // check if column is now fixed
                     if( ( boundChg.upper &&
                           boundChg.bound == lower_bounds[boundChg.col] ) ||
                        ( !boundChg.upper &&
                           boundChg.bound == upper_bounds[boundChg.col] ) )
                        ++ncliquefixings;
                     else
                        ++ncliqueboundchgs;
                  }
                  else
                  {
                     // already changed that bound
                     CliqueProbingBoundChg<REAL>& cliqueOtherBoundChg = cliqueBoundChanges
                        [cliqueBoundPos[2 * boundChg.col + boundChg.upper] - 1];

                     if( boundChg.upper && boundChg.bound < cliqueOtherBoundChg.bound )
                     {
                        // new upper bound change is tighter
                        cliqueOtherBoundChg.bound = boundChg.bound;

                        // check if column is now fixed
                        if( boundChg.bound == lower_bounds[boundChg.col] )
                           ++ncliquefixings;

                        if( problemUpdate.getPresolveOptions()
                                 .verification_with_VeriPB )
                        {
                           if( boundChg.probing_col == -1 )
                              cliqueOtherBoundChg.probing_col = -1;
                        }
                     }
                     else if( !boundChg.upper &&
                              boundChg.bound > cliqueOtherBoundChg.bound )
                     {
                        // new lower bound change is tighter
                        cliqueOtherBoundChg.bound = boundChg.bound;

                        // check if column is now fixed
                        if( boundChg.bound == upper_bounds[boundChg.col] )
                           ++ncliquefixings;
                     }

                     // do only count fixings in this case for two reasons:
                     // 1) the number of bound changes depends on the order and
                     // would make probing non deterministic 2) the boundchange
                     // was already counted in previous rounds and will only be
                     // added once
                  }
               }

      ncliquesubstitutions += cliquesubstitutions.size();
   }

   PresolveStatus result = PresolveStatus::kUnchanged;

   if( unsuccessfulcliqueprobing <= numcliquefails )
   {
      std::cout<<"\n\nClique Probing on ";
      std::cout<<std::min(static_cast<int>(probingCliques.size()), batchend);
      std::cout<<" Cliques with ";
      std::cout<< totalnumpropagations;
      std::cout<<" propagations on ";
      std::cout<<nprobedvars<<" duplicate variables and " << probedvars.size() << " individual variables led to ";
      std::cout<<ncliquefixings;
      std::cout<<" fixings, ";
#ifdef PAPILO_TBB
      std::cout<<static_cast<int>(change_to_equation_comb.size());
#else
      std::cout<<static_cast<int>(change_to_equation.size());
#endif
      std::cout<<" changed lhs/rhs, ";
      std::cout<<static_cast<int>(cliqueBoundChanges.size());
      std::cout<<" ";
      std::cout<<ncliqueboundchgs;
      std::cout<<" Bound Changes and ";
      std::cout<<static_cast<int>(cliquesubstitutions.size());
      std::cout<<" ";
      std::cout<<ncliquesubstitutions;
      std::cout<<" Substitutions in \n";
      std::cout<<cliqueprobingtime;
      std::cout<<" seconds.";
      std::cout<<"\n\n\nPerformance Ratio: ";
      std::cout<< static_cast<float>(ncliquefixings + ncliqueboundchgs + ncliquesubstitutions)
         / static_cast<float>(totalnumpropagations);
      std::cout<<"\n\n\n";
      if( ncliquefixings + ncliqueboundchgs + ncliquesubstitutions == 0 || earlycliqueabort )
            unsuccessfulcliqueprobing += 1;
      else
            unsuccessfulcliqueprobing = 0;

      if( !cliqueBoundChanges.empty() )
      {
         pdqsort( cliqueBoundChanges.begin(), cliqueBoundChanges.end(),
            []( const CliqueProbingBoundChg<REAL>& a, const CliqueProbingBoundChg<REAL>& b )
            {
               return std::make_pair( a.col, a.bound ) > std::make_pair( b.col, b.bound );
            } );
         for( const CliqueProbingBoundChg<REAL>& boundChg : cliqueBoundChanges )
         {
            bool binary = problem.getVariableDomains().isBinary( boundChg.col );
            if( boundChg.upper )
            {
               if( problemUpdate.getPresolveOptions().verification_with_VeriPB &&
                  boundChg.probing_col != -1 )
                  reductions.reason_probing_upper_bound_change(
                     boundChg.probing_col, boundChg.col );
               reductions.changeColUB( boundChg.col, boundChg.bound );
               if( binary )
               {
                  probing_scores[boundChg.col] = -100000;
               }
            }
            else
            {
               if( problemUpdate.getPresolveOptions().verification_with_VeriPB &&
                  boundChg.probing_col != -1 )
                  reductions.reason_probing_lower_bound_change(
                     boundChg.probing_col, boundChg.col );
               reductions.changeColLB( boundChg.col, boundChg.bound );
               if( binary )
               {
                  probing_scores[boundChg.col] = -100000;
               }
            }
         }

      result = PresolveStatus::kReduced;
      }

      if( !cliquesubstitutions.empty() )
      {
         pdqsort( cliquesubstitutions.begin(), cliquesubstitutions.end(),
                  []( const CliqueProbingSubstitution<REAL>& a,
                     const CliqueProbingSubstitution<REAL>& b )
                  {
                     return std::make_pair( a.col1, a.col2 ) >
                           std::make_pair( b.col1, b.col2 );
                  } );

         int lastsubstcol = -1;

         for( const CliqueProbingSubstitution<REAL>& subst : cliquesubstitutions )
         {
            if( subst.col1 == lastsubstcol )
               continue;

            lastsubstcol = subst.col1;
            if(false)
               reductions.replaceCol( subst.col1, subst.col2, subst.col2scale,
                                    subst.col2const );
         }

         result = PresolveStatus::kReduced;
      }
   }

   pdqsort( probing_cands.begin(), probing_cands.end(),
            [this, &probing_scores, &colsize, &colperm]( int col1, int col2 )
            {
               std::pair<double, double> s1;
               std::pair<double, double> s2;
               if( nprobed[col2] == 0 && probing_scores[col2] != 0 )
                  s2.first = probing_scores[col2] /
                             static_cast<double>( colsize[col2] );
               else
                  s2.first = 0;
               if( nprobed[col1] == 0 && probing_scores[col1] != 0 )
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


   int clique_cutoff_ub = 0;

   if( unsuccessfulcliqueprobing <= numcliquefails )
   {
      clique_cutoff_ub = static_cast<int>(probing_cands.size())-1;
      int clique_cutoff_lb = 0;

      assert( clique_cutoff_ub < static_cast<int>(probing_cands.size()));
      if( clique_cutoff_ub != -1 && probing_scores[probing_cands[clique_cutoff_ub]] < 0 )
      {
         while (clique_cutoff_ub - clique_cutoff_lb > 1 )
         {
            if( probing_scores[probing_cands[ ( clique_cutoff_ub + clique_cutoff_lb ) / 2 ]] >= 0 )
            {
               clique_cutoff_lb = ( clique_cutoff_ub + clique_cutoff_lb ) / 2;
            }
            else if( probing_scores[probing_cands[ ( clique_cutoff_ub + clique_cutoff_lb ) / 2 ]] < 0 )
            {
               clique_cutoff_ub = ( clique_cutoff_ub + clique_cutoff_lb ) / 2;
            }
         }
      }
      if( probing_scores[probing_cands[0]] < 0 )
      {
         probing_cands.clear();
         return result;
      } 

      assert( clique_cutoff_ub + 1 == static_cast<int>(probing_cands.size()) 
              || probing_scores[ probing_cands[ clique_cutoff_ub + 1 ]] < 0 );
      probing_cands.resize(clique_cutoff_ub+1);
      assert( probing_scores[ probing_cands[ static_cast<int>(probing_cands.size()) - 1 ] ] >= 0 );
   }

   const Vec<int>& rowsize = consMatrix.getRowSizes();

   int current_badge_start = 0;
   amountofwork = 0;

   int64_t working_limit = consMatrix.getNnz() * 2;
   int initial_badge_limit = 0.1 * working_limit;

   const int nprobingcands = static_cast<int>( probing_cands.size() );
   int badge_size = 0;
   for( int i : probing_cands )
   {
      ++badge_size;

      if( badge_size == maxinitialbadgesize )
         break;

      initial_badge_limit -= colsize[i];
      if( initial_badge_limit <= 0 )
         break;

      auto colvec = consMatrix.getColumnCoefficients( i );
      const int* rowinds = colvec.getIndices();
      for( int k = 0; k != colvec.getLength(); ++k )
      {
         initial_badge_limit -= ( rowsize[rowinds[k]] - 1 );

         if( initial_badge_limit <= 0 )
            break;
      }

      if( initial_badge_limit <= 0 )
         break;
   }


   badge_size = std::max( std::min( nprobingcands, minbadgesize ), badge_size );

   int current_badge_end = std::max( 0, std::min(static_cast<int>(probing_cands.size()), current_badge_start + badge_size));
   assert(current_badge_end <= static_cast<int>(probing_cands.size()));
   int n_useless = 0;
   bool abort = false;

   HashMap<std::pair<int, int>, int, boost::hash<std::pair<int, int>>>
   substitutionsPos;
   Vec<ProbingSubstitution<REAL>> substitutions;
   Vec<int> boundPos( size_t( 2 * ncols ), 0 );
   Vec<ProbingBoundChg<REAL>> boundChanges;
   boundChanges.reserve( ncols );


   // use tbb combinable so that each thread will copy the activities and
   // bounds at most once
#ifdef PAPILO_TBB
   tbb::combinable<ProbingView<REAL>> probing_views(
       [this, &problem, &num, &cliqueBoundChanges]()
       {
          ProbingView<REAL> probingView( problem, num, cliqueBoundChanges );
          probingView.setMinContDomRed( mincontdomred );
          return probingView;
       } );
#else
   ProbingView<REAL> probingView( problem, num, cliqueBoundChanges );
   probingView.setMinContDomRed( mincontdomred );
#endif

   do
   {
      Message::debug( this, "probing candidates {} to {}\n",
                      current_badge_start, current_badge_end );

      auto propagate_variables = [&]( int start, int end )
      {
#ifdef PAPILO_TBB
         tbb::parallel_for(
             tbb::blocked_range<int>( start, end ),
             [&]( const tbb::blocked_range<int>& r )
             {
                ProbingView<REAL>& probingView = probing_views.local();

                for( int i = r.begin(); i != r.end(); ++i )
#else
         for( int i = start; i < end; i++ )
#endif
                {
                   if( PresolveMethod<REAL>::is_interrupted(
                           timer, problemUpdate.getPresolveOptions().tlim, problemUpdate.getPresolveOptions().early_exit_callback ) )
                      break;
                  assert(i >= 0);
                  assert(i < static_cast<int>(probing_cands.size()));
                   const int col = probing_cands[i];

                   assert(
                       cflags[col].test( ColFlag::kIntegral ) &&
                       ( lower_bounds[col] == 0 || upper_bounds[col] == 1 ) );

                   if(  probingView.origin_upper_bounds[col] == probingView.origin_lower_bounds[col] )
                   {
                      assert( probing_scores[col] == -100000 );
                      continue;
                   }

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
                      infeasible_variable.store( col );
                      break;
                   }
                }
#ifdef PAPILO_TBB
             } );
#endif
      };

      int nfixings = 0;
      int nboundchgs = 0;
      int nsubstitutions = -substitutions.size();

      assert(current_badge_end <= static_cast<int>(probing_cands.size()));
      assert(current_badge_end >= 0);
      assert(current_badge_start >= 0 );
      assert(current_badge_start <= current_badge_end );
      auto probingstarttime = timer.getTime();
      propagate_variables( current_badge_start, current_badge_end);
      auto probingtime = timer.getTime() - probingstarttime;

      if( PresolveMethod<REAL>::is_time_exceeded(
              timer, problemUpdate.getPresolveOptions().tlim ) )
      {
         std::cout<<"\nTime limit of probing exceeded!";
         return PresolveStatus::kUnchanged;
      }

      if( infeasible.load( std::memory_order_relaxed ) )
      {
         reason_of_infeasibility =
             infeasible_variable.load( std::memory_order_relaxed );
         return PresolveStatus::kInfeasible;
      }

#ifdef PAPILO_TBB
      probing_views.combine_each(
          [&]( ProbingView<REAL>& probingView )
          {
#endif
             const auto& probingBoundChgs =
                 probingView.getProbingBoundChanges();
             const auto& probingSubstitutions =
                 probingView.getProbingSubstitutions();

             amountofwork += probingView.getAmountOfWork();

             for( const ProbingSubstitution<REAL>& subst :
                  probingSubstitutions )
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
                   ProbingBoundChg<REAL>& otherBoundChg = boundChanges
                       [boundPos[2 * boundChg.col + boundChg.upper] - 1];

                   if( boundChg.upper && boundChg.bound < otherBoundChg.bound )
                   {
                      // new upper bound change is tighter
                      otherBoundChg.bound = boundChg.bound;

                      // check if column is now fixed
                      if( boundChg.bound == lower_bounds[boundChg.col] )
                         ++nfixings;

                      if( problemUpdate.getPresolveOptions()
                              .verification_with_VeriPB )
                      {
                         if( boundChg.probing_col == -1 )
                            otherBoundChg.probing_col = -1;
                      }
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
                   // would make probing non deterministic 2) the boundchange
                   // was already counted in previous rounds and will only be
                   // added once
                }
             }

             probingView.clearResults();
#ifdef PAPILO_TBB
          } );
#endif
      nsubstitutions += substitutions.size();

      std::cout<<"\nNormal probing on ";
      std::cout<<badge_size;
      std::cout<<" variables found ";
      std::cout<< nfixings;
      std::cout<<" fixings, ";
      std::cout<< nsubstitutions;
      std::cout<<" Substitutions and ";
      std::cout<< nboundchgs;
      std::cout<<" Boundchanges in ";
      std::cout<< probingtime;
      std::cout<<" seconds.\n";
      std::cout<<"\n\n\nPerfomance ratio: ";
      std::cout<< static_cast<float>(nsubstitutions + nboundchgs + nfixings) / static_cast<float>( 2 * ( current_badge_end - current_badge_start ) );
      std::cout<<"\n\n";

      current_badge_start = current_badge_end;

      if( nfixings == 0 && nboundchgs == 0 && nsubstitutions == 0 )
         n_useless += amountofwork;
      else
         n_useless = 0;

      Message::debug(
          this,
          "probing found: {} fixings, {} substitutions, {} bound changes\n",
          nfixings, nsubstitutions, nboundchgs );

      int64_t extrawork =
          ( ( 0.1 * ( nfixings + nsubstitutions ) + 0.01 * nboundchgs ) *
            consMatrix.getNnz() );

      working_limit -= amountofwork;
      working_limit += extrawork;

      if (amountofwork != 0)
      badge_size = static_cast<int>(
          ceil( badge_size * static_cast<double>( working_limit + extrawork ) /
                (double)amountofwork ) );
      else
         badge_size = nprobingcands - current_badge_start;
      badge_size = std::min( nprobingcands - current_badge_start, badge_size );
      if( max_badge_size > 0 )
         badge_size = std::min( max_badge_size, badge_size );
      current_badge_end = std::max(current_badge_start + badge_size, 0);

      abort = n_useless >= consMatrix.getNnz() * 2 || working_limit < 0 ||
              current_badge_start == current_badge_end ||
              PresolveMethod<REAL>::is_time_exceeded(
                  timer, problemUpdate.getPresolveOptions().tlim );
   } while( !abort );

   if( !boundChanges.empty() )
   {
      pdqsort( boundChanges.begin(), boundChanges.end(),
            []( const ProbingBoundChg<REAL>& a, const ProbingBoundChg<REAL>& b )
            {
               return std::make_pair( a.col, a.bound ) > std::make_pair( b.col, b.bound );
            } );
      for( const ProbingBoundChg<REAL>& boundChg : boundChanges )
      {
         bool binary = problem.getVariableDomains().isBinary( boundChg.col );
         if( boundChg.upper )
         {
            if( problemUpdate.getPresolveOptions().verification_with_VeriPB &&
                boundChg.probing_col != -1 )
               reductions.reason_probing_upper_bound_change(
                   boundChg.probing_col, boundChg.col );
            reductions.changeColUB( boundChg.col, boundChg.bound );
            if(binary)
            {
               probing_scores[boundChg.col] = -100000;
               nprobed[boundChg.col] = -100000;
            }
         }
         else
         {
            if( problemUpdate.getPresolveOptions().verification_with_VeriPB &&
                boundChg.probing_col != -1 )
               reductions.reason_probing_lower_bound_change(
                   boundChg.probing_col, boundChg.col );
            reductions.changeColLB( boundChg.col, boundChg.bound );
            if(binary)
            {
               probing_scores[boundChg.col] = -100000;
               nprobed[boundChg.col] = -100000;
            }
         }
      }

      result = PresolveStatus::kReduced;
   }

   if( !substitutions.empty() )
   {
      pdqsort( substitutions.begin(), substitutions.end(),
               []( const ProbingSubstitution<REAL>& a,
                   const ProbingSubstitution<REAL>& b )
               {
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

   if( unsuccessfulcliqueprobing <= numcliquefails )
   {
      assert( ncliquefixings + ncliqueboundchgs + ncliquesubstitutions == 0
         || result == PresolveStatus::kInfeasible || result == PresolveStatus::kReduced);
   }

   return result;
}

template <typename REAL>
bool
Probing<REAL>::isBinaryVariable( REAL upper_bound, REAL lower_bound,
                                 int column_size,
                                 const Flags<ColFlag>& colFlag ) const
{
   return !colFlag.test( ColFlag::kUnbounded ) &&
          colFlag.test( ColFlag::kIntegral ) && column_size > 0 &&
          lower_bound == 0 && upper_bound == 1;
}

template <typename REAL>
void
Probing<REAL>::set_max_badge_size( int val )
{
   max_badge_size = val;
}

} // namespace papilo

#endif
