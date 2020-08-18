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

#ifndef _PAPILO_CORE_PRESOLVE_HPP_
#define _PAPILO_CORE_PRESOLVE_HPP_

#include <algorithm>
#include <cctype>
#include <fstream>
#include <initializer_list>
#include <memory>
#include <utility>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include "papilo/core/PresolveMethod.hpp"
#include "papilo/core/PresolveOptions.hpp"
#include "papilo/core/Problem.hpp"
#include "papilo/core/ProblemUpdate.hpp"
#include "papilo/core/Statistics.hpp"
#include "papilo/interfaces/SolverInterface.hpp"
#include "papilo/io/Message.hpp"
#include "papilo/misc/DependentRows.hpp"
#include "papilo/misc/ParameterSet.hpp"
#include "papilo/misc/Timer.hpp"
#include "papilo/misc/Vec.hpp"
#include "papilo/misc/tbb.hpp"
#include "papilo/presolvers/CoefficientStrengthening.hpp"
#include "papilo/presolvers/ConstraintPropagation.hpp"
#include "papilo/presolvers/DominatedCols.hpp"
#include "papilo/presolvers/DualFix.hpp"
#include "papilo/presolvers/DualInfer.hpp"
#include "papilo/presolvers/FixContinuous.hpp"
#include "papilo/presolvers/FreeVarSubstitution.hpp"
#include "papilo/presolvers/ImplIntDetection.hpp"
#include "papilo/presolvers/ParallelColDetection.hpp"
#include "papilo/presolvers/ParallelRowDetection.hpp"
#include "papilo/presolvers/Probing.hpp"
#include "papilo/presolvers/SimpleProbing.hpp"
#include "papilo/presolvers/SimpleSubstitution.hpp"
#include "papilo/presolvers/SimplifyInequalities.hpp"
#include "papilo/presolvers/SingletonCols.hpp"
#include "papilo/presolvers/SingletonStuffing.hpp"
#include "papilo/presolvers/Sparsify.hpp"

namespace papilo
{

template <typename REAL>
struct PresolveResult
{
   Postsolve<REAL> postsolve;
   PresolveStatus status;
};

template <typename REAL>
class Presolve
{
 public:
   void
   addDefaultPresolvers()
   {
      using uptr = std::unique_ptr<PresolveMethod<REAL>>;

      addPresolveMethod( uptr( new SingletonCols<REAL>() ) );
      addPresolveMethod( uptr( new CoefficientStrengthening<REAL>() ) );
      addPresolveMethod( uptr( new SimpleProbing<REAL>() ) );
      addPresolveMethod( uptr( new ConstraintPropagation<REAL>() ) );
      addPresolveMethod( uptr( new SingletonStuffing<REAL>() ) );
      addPresolveMethod( uptr( new DualFix<REAL>() ) );
      addPresolveMethod( uptr( new ImplIntDetection<REAL>() ) );
      addPresolveMethod( uptr( new FixContinuous<REAL>() ) );
      addPresolveMethod( uptr( new ParallelRowDetection<REAL>() ) );
      addPresolveMethod( uptr( new ParallelColDetection<REAL>() ) );
      addPresolveMethod( uptr( new SimpleSubstitution<REAL>() ) );
      addPresolveMethod( uptr( new DualInfer<REAL> ) );
      addPresolveMethod( uptr( new Substitution<REAL>() ) );
      addPresolveMethod( uptr( new Probing<REAL>() ) );
      addPresolveMethod( uptr( new DominatedCols<REAL>() ) );
      addPresolveMethod( uptr( new Sparsify<REAL>() ) );
      addPresolveMethod( uptr( new SimplifyInequalities<REAL>() ) );
   }

   ParameterSet
   getParameters()
   {
      ParameterSet paramSet;
      msg.addParameters( paramSet );
      presolveOptions.addParameters( paramSet );

      for( const std::unique_ptr<PresolveMethod<REAL>>& presolver : presolvers )
         presolver->addParameters( paramSet );

      return paramSet;
   }

   /// apply presolving to problem
   PresolveResult<REAL>
   apply( Problem<REAL>& problem );

   /// add presolve method to presolving
   void
   addPresolveMethod( std::unique_ptr<PresolveMethod<REAL>> presolveMethod )
   {
      presolvers.emplace_back( std::move( presolveMethod ) );
   }

   void
   setLPSolverFactory( std::unique_ptr<SolverFactory<REAL>> lpSolverFactory )
   {
      this->lpSolverFactory = std::move( lpSolverFactory );
   }

   void
   setMIPSolverFactory( std::unique_ptr<SolverFactory<REAL>> mipSolverFactory )
   {
      this->mipSolverFactory = std::move( mipSolverFactory );
   }

   const std::unique_ptr<SolverFactory<REAL>>&
   getLPSolverFactory() const
   {
      return this->lpSolverFactory;
   }

   const std::unique_ptr<SolverFactory<REAL>>&
   getMIPSolverFactory() const
   {
      return this->mipSolverFactory;
   }

   void
   setPresolverOptions( const PresolveOptions& presolveOptions )
   {
      this->presolveOptions = presolveOptions;
   }

   const PresolveOptions&
   getPresolveOptions() const
   {
      return this->presolveOptions;
   }

   PresolveOptions&
   getPresolveOptions()
   {
      return this->presolveOptions;
   }

   /// get epsilon value for numerical comparisons
   const REAL&
   getEpsilon() const
   {
      return num.getEpsilon();
   }

   /// get feasibility tolerance value
   const REAL&
   getFeasTol() const
   {
      return num.getFeasTol();
   }

   /// set the verbosity level
   void
   setVerbosityLevel( VerbosityLevel verbosity )
   {
      msg.setVerbosityLevel( verbosity );
   }

   /// get the verbosity level
   VerbosityLevel
   getVerbosityLevel() const
   {
      return msg.getVerbosityLevel();
   }

   const Message&
   message() const
   {
      return msg;
   }

   Message&
   message()
   {
      return msg;
   }

   /// access statistics of presolving
   const Statistics&
   getStatistics() const
   {
      return stats;
   }

 private:
   /// evaluate result array of each presolver, return the largest result value
   PresolveStatus
   evaluateResults();

   std::pair<int, int>
   applyReductions( int p, const Reductions<REAL>& reductions,
                    ProblemUpdate<REAL>& probUpdate );

   void
   finishRound( ProblemUpdate<REAL>& probUpdate );

   void
   applyPostponed( ProblemUpdate<REAL>& probUpdate );

   bool
   updateRoundCounter( Problem<REAL>& problem, ProblemUpdate<REAL>& probUpdate,
                       const Statistics& roundStats, const Timer& presolvetimer,
                       bool unchanched = false );

   bool
   applyPresolversReductions( ProblemUpdate<REAL>& probUpdate );

   void
   printRoundStats( bool unchanged = false );

   void
   printPresolversStats();

   // data to perform presolving
   Vec<PresolveStatus> results;
   Vec<std::unique_ptr<PresolveMethod<REAL>>> presolvers;
   Vec<Reductions<REAL>> reductions;
   int roundCounter;

   Vec<std::pair<const Reduction<REAL>*, const Reduction<REAL>*>>
       postponedReductions;
   Vec<int> postponedReductionToPresolver;

   // settings for presolve behavior
   Num<REAL> num;
   Message msg;
   PresolveOptions presolveOptions;
   // statistics
   Statistics stats;

   std::unique_ptr<SolverFactory<REAL>> lpSolverFactory;
   std::unique_ptr<SolverFactory<REAL>> mipSolverFactory;

   Vec<std::pair<int, int>> presolverStats;
   bool lastRoundReduced;
   int nunsuccessful;
   bool rundelayed;
};

#ifdef PAPILO_USE_EXTERN_TEMPLATES
extern template class Presolve<double>;
extern template class Presolve<Quad>;
extern template class Presolve<Rational>;
#endif

/// evaluate result array of each presolver, return the largest result value
template <typename REAL>
PresolveStatus
Presolve<REAL>::evaluateResults()
{
   int result = static_cast<int>( PresolveStatus::kUnchanged );

   for( std::size_t i = 0; i < results.size(); ++i )
      result = std::max( result, static_cast<int>( results[i] ) );

   return static_cast<PresolveStatus>( result );
}

template <typename REAL>
std::pair<int, int>
Presolve<REAL>::applyReductions( int p, const Reductions<REAL>& reductions,
                                 ProblemUpdate<REAL>& probUpdate )
{
   int k = 0;
   ApplyResult result;
   int nbtsxAppliedStart = stats.ntsxapplied;
   int nbtsxTotal = 0;

   const auto& reds = reductions.getReductions();
   const auto& tsx = reductions.getTransactions();

   for( const auto& transaction : reductions.getTransactions() )
   {
      int start = transaction.start;
      int end = transaction.end;

      for( ; k != start; ++k )
      {
         result = probUpdate.applyTransaction( &reds[k], &reds[k + 1] );
         if( result == ApplyResult::kApplied )
            ++stats.ntsxapplied;
         else if( result == ApplyResult::kRejected )
            ++stats.ntsxconflicts;
         else if( result == ApplyResult::kInfeasible )
            return std::make_pair( -1, -1 );
         else if( result == ApplyResult::kPostponed )
            postponedReductions.emplace_back( &reds[k], &reds[k + 1] );

         ++nbtsxTotal;
      }

      result = probUpdate.applyTransaction( &reds[start], &reds[end] );
      if( result == ApplyResult::kApplied )
         ++stats.ntsxapplied;
      else if( result == ApplyResult::kRejected )
         ++stats.ntsxconflicts;
      else if( result == ApplyResult::kInfeasible )
         return std::make_pair( -1, -1 );
      else if( result == ApplyResult::kPostponed )
         postponedReductions.emplace_back( &reds[start], &reds[end] );

      k = end;
      ++nbtsxTotal;
   }

   for( ; k != static_cast<int>( reds.size() ); ++k )
   {
      result = probUpdate.applyTransaction( &reds[k], &reds[k + 1] );
      if( result == ApplyResult::kApplied )
         ++stats.ntsxapplied;
      else if( result == ApplyResult::kRejected )
         ++stats.ntsxconflicts;
      else if( result == ApplyResult::kInfeasible )
         return std::make_pair( -1, -1 );
      else if( result == ApplyResult::kPostponed )
         postponedReductions.emplace_back( &reds[k], &reds[k + 1] );

      ++nbtsxTotal;
   }

   return std::pair<int, int>( nbtsxTotal,
                               ( stats.ntsxapplied - nbtsxAppliedStart ) );
}

template <typename REAL>
void
Presolve<REAL>::finishRound( ProblemUpdate<REAL>& probUpdate )
{
   probUpdate.clearStates();

   // clear reductions
   for( auto& reduction : reductions )
      reduction.clear();

   std::fill( results.begin(), results.end(), PresolveStatus::kUnchanged );

   // TODO compress if problem size decreased by some factor
}

template <typename REAL>
void
Presolve<REAL>::applyPostponed( ProblemUpdate<REAL>& probUpdate )
{
   probUpdate.setPostponeSubstitutions( false );

   // apply all postponed reductions
   for( int presolver = 0; presolver != presolvers.size(); ++presolver )
   {
      int first = postponedReductionToPresolver[presolver];
      int last = postponedReductionToPresolver[presolver + 1];
      for( int i = first; i != last; ++i )
      {
         const auto& ptrpair = postponedReductions[i];

         ApplyResult r =
             probUpdate.applyTransaction( ptrpair.first, ptrpair.second );
         if( r == ApplyResult::kApplied )
         {
            ++stats.ntsxapplied;
            ++presolverStats[presolver].second;
         }
         else if( r == ApplyResult::kRejected )
            ++stats.ntsxconflicts;
      }
   }

   postponedReductions.clear();
   postponedReductionToPresolver.clear();
}

template <typename REAL>
bool
Presolve<REAL>::updateRoundCounter( Problem<REAL>& problem,
                                    ProblemUpdate<REAL>& probUpdate,
                                    const Statistics& roundStats,
                                    const Timer& presolvetimer,
                                    bool unchanched )
{
   if( presolveOptions.tlim != std::numeric_limits<double>::max() )
   {
      if( presolvetimer.getTime() >= presolveOptions.tlim )
         return true;
   }

   if( !unchanched )
   {
      double abortfac = problem.getNumIntegralCols() == 0
                            ? presolveOptions.lpabortfac
                            : presolveOptions.abortfac;
      // update statistics
      bool increment =
          ( 0.1 * roundStats.nboundchgs + roundStats.ndeletedcols ) <=
          abortfac * probUpdate.getNActiveCols();
      increment =
          increment && ( roundStats.nsidechgs + roundStats.ndeletedrows ) <=
                           abortfac * probUpdate.getNActiveRows();
      increment =
          increment && ( roundStats.ncoefchgs <=
                         abortfac * problem.getConstraintMatrix().getNnz() );

      if( increment )
      {
         lastRoundReduced =
             lastRoundReduced || roundStats.nsidechgs > 0 ||
             roundStats.nboundchgs > 0 || roundStats.ndeletedcols > 0 ||
             roundStats.ndeletedrows > 0 || roundStats.ncoefchgs > 0;
         ++roundCounter;
      }
      else
      {
         printRoundStats();
         lastRoundReduced = true;
         roundCounter = 0;
         nunsuccessful = 0;
      }
   }
   else
      ++roundCounter;

   bool abort = false;

   if( roundCounter == 3 )
   {
      ++nunsuccessful;

      abort = rundelayed && ( !lastRoundReduced || nunsuccessful == 2 );

      if( !abort )
      {
         roundCounter = 2;
         printRoundStats( !lastRoundReduced );

         if( !rundelayed )
         {
            msg.info( "activating delayed presolvers\n" );
            for( auto& p : presolvers )
               p->setDelayed( false );
            rundelayed = true;
         }

         roundCounter = 0;
      }
      else
         printRoundStats( !lastRoundReduced );
   }

   if( roundCounter == 0 )
      ++stats.nrounds;

   assert( roundCounter != 3 || abort );

   return abort;
}

template <typename REAL>
bool
Presolve<REAL>::applyPresolversReductions( ProblemUpdate<REAL>& probUpdate )
{
   probUpdate.setPostponeSubstitutions( true );

   postponedReductionToPresolver.push_back( 0 );

   for( std::size_t i = 0; i < presolvers.size(); ++i )
   {
      if( results[i] == PresolveStatus::kReduced )
      {
         Message::debug( this, "applying reductions of presolver {}\n",
                         presolvers[i]->getName() );

         auto stats = applyReductions( i, reductions[i], probUpdate );

         if( stats.first < 0 || stats.second < 0 )
            return false;

         presolverStats[i].first += stats.first;
         presolverStats[i].second += stats.second;
         results[i] = PresolveStatus::kUnchanged;
      }

      postponedReductionToPresolver.push_back( postponedReductions.size() );
   }

   probUpdate.flushChangedCoeffs();

   applyPostponed( probUpdate );

   if( probUpdate.flush() == PresolveStatus::kInfeasible )
      return false;

   return true;
}

template <typename REAL>
void
Presolve<REAL>::printRoundStats( bool unchanged )
{
   std::string rndtype;
   switch( roundCounter )
   {
   case 0:
      rndtype = "Fast";
      break;
   case 1:
      rndtype = "Medium";
      break;
   case 2:
      rndtype = "Exhaustive";
      break;
   case 3:
      rndtype = "Final";
      break;
   case 4:
      rndtype = "Trivial";
   }

   if( unchanged )
   {
      msg.info( "round {:<3} ({:^10}): Unchanged\n", stats.nrounds, rndtype );
      return;
   }

   msg.info( "round {:<3} ({:^10}): {:>4} del cols, {:>4} del rows, "
             "{:>4} chg bounds, {:>4} chg sides, {:>4} chg coeffs, "
             "{:>4} tsx applied, {:>4} tsx conflicts\n",
             stats.nrounds, rndtype, stats.ndeletedcols, stats.ndeletedrows,
             stats.nboundchgs, stats.nsidechgs, stats.ncoefchgs,
             stats.ntsxapplied, stats.ntsxconflicts );
}

template <typename REAL>
void
Presolve<REAL>::printPresolversStats()
{
   msg.info( "presolved {} rounds: {:>4} del cols, {:>4} del rows, "
             "{:>4} chg bounds, {:>4} chg sides, {:>4} chg coeffs, "
             "{:>4} tsx applied, {:>4} tsx conflicts\n",
             stats.nrounds, stats.ndeletedcols, stats.ndeletedrows,
             stats.nboundchgs, stats.nsidechgs, stats.ncoefchgs,
             stats.ntsxapplied, stats.ntsxconflicts );
   msg.info( "\n {:>18} {:>12} {:>18} {:>18} {:>18} {:>18} \n", "presolver",
             "nb calls", "success calls(%)", "nb transactions",
             "tsx applied(%)", "execution time(s)" );
   for( std::size_t i = 0; i < presolvers.size(); ++i )
   {
      presolvers[i]->printStats( msg, presolverStats[i] );
   }

   msg.info( "\n" );
}

/// apply presolving to problem
template <typename REAL>
PresolveResult<REAL>
Presolve<REAL>::apply( Problem<REAL>& problem )
{
   tbb::task_arena arena( presolveOptions.threads == 0
                              ? tbb::task_arena::automatic
                              : presolveOptions.threads );

   return arena.execute( [this, &problem]() {
      stats = Statistics();
      num.setFeasTol( REAL{ presolveOptions.feastol } );
      num.setEpsilon( REAL{ presolveOptions.epsilon } );
      num.setHugeVal( REAL{ presolveOptions.hugeval } );

      Timer timer( stats.presolvetime );

      ConstraintMatrix<REAL>& constraintMatrix = problem.getConstraintMatrix();
      VariableDomains<REAL>& variableDomains = problem.getVariableDomains();
      Vec<REAL>& lhsVals = constraintMatrix.getLeftHandSides();
      Vec<REAL>& rhsVals = constraintMatrix.getRightHandSides();
      Vec<RowFlags>& rflags = constraintMatrix.getRowFlags();
      const Vec<int>& rowsize = constraintMatrix.getRowSizes();
      Vec<RowActivity<REAL>>& rowActivities = problem.getRowActivities();

      msg.info( "\nstarting presolve of problem {}:\n", problem.getName() );
      msg.info( "  rows:     {}\n", problem.getNRows() );
      msg.info( "  columns:  {} ({} int., {} cont.)\n", problem.getNCols(),
                problem.getNumIntegralCols(), problem.getNumContinuousCols() );
      msg.info( "  nonzeros: {}\n\n", problem.getConstraintMatrix().getNnz() );

      PresolveResult<REAL> result;

      result.postsolve = Postsolve<REAL>( problem, num );
      result.postsolve.getChecker().setOriginalProblem( problem );

      result.status = PresolveStatus::kUnchanged;

      std::stable_sort( presolvers.begin(), presolvers.end(),
                        []( const std::unique_ptr<PresolveMethod<REAL>>& a,
                            const std::unique_ptr<PresolveMethod<REAL>>& b ) {
                           return static_cast<int>( a->getTiming() ) <
                                  static_cast<int>( b->getTiming() );
                        } );

      std::pair<int, int> fastPresolvers;
      std::pair<int, int> mediumPresolvers;
      std::pair<int, int> exhaustivePresolvers;

      int npresolvers = static_cast<int>( presolvers.size() );

      fastPresolvers.first = fastPresolvers.second = 0;
      while( fastPresolvers.second < npresolvers &&
             presolvers[fastPresolvers.second]->getTiming() ==
                 PresolverTiming::kFast )
         ++fastPresolvers.second;

      mediumPresolvers.first = mediumPresolvers.second = fastPresolvers.second;
      while( mediumPresolvers.second < npresolvers &&
             presolvers[mediumPresolvers.second]->getTiming() ==
                 PresolverTiming::kMedium )
         ++mediumPresolvers.second;

      exhaustivePresolvers.first = exhaustivePresolvers.second =
          mediumPresolvers.second;
      while( exhaustivePresolvers.second < npresolvers &&
             presolvers[exhaustivePresolvers.second]->getTiming() ==
                 PresolverTiming::kExhaustive )
         ++exhaustivePresolvers.second;

      reductions.resize( presolvers.size() );
      results.resize( presolvers.size() );

      roundCounter = 0;

      presolverStats.resize( presolvers.size(), std::pair<int, int>( 0, 0 ) );

      // todo move trivial round 0 to initiliaze function

      ProblemUpdate<REAL> probUpdate( problem, result.postsolve, stats,
                                      presolveOptions, num );

      for( int i = 0; i != npresolvers; ++i )
      {
         if( presolvers[i]->isEnabled() )
         {
            if( presolvers[i]->initialize( problem, presolveOptions ) )
               probUpdate.observeCompress( presolvers[i].get() );
         }
      }

      result.status = probUpdate.trivialPresolve();

      if( result.status == PresolveStatus::kInfeasible ||
          result.status == PresolveStatus::kUnbndOrInfeas ||
          result.status == PresolveStatus::kUnbounded )
         return result;

      roundCounter = 4;
      printRoundStats();
      roundCounter = 0;

      finishRound( probUpdate );
      ++stats.nrounds;

// #define PARALLEL_FAST_PRESOLVERS
#define PARALLEL_MEDIUM_PRESOLVERS
#define PARALLEL_EXHAUSTIVE_PRESOLVERS

      nunsuccessful = 0;
      rundelayed = true;
      for( int i = 0; i < npresolvers; ++i )
      {
         if( presolvers[i]->isEnabled() && presolvers[i]->isDelayed() )
         {
            rundelayed = false;
            break;
         }
      }
      bool abort = false;
      do
      {
         // if problem is trivial abort here
         if( probUpdate.getNActiveCols() == 0 ||
             probUpdate.getNActiveRows() == 0 )
            break;

         // call presolvers
         switch( roundCounter )
         {
         case 0:
#ifdef PARALLEL_FAST_PRESOLVERS
            tbb::parallel_for(
                tbb::blocked_range<int>( fastPresolvers.first,
                                         fastPresolvers.second ),
                [&]( const tbb::blocked_range<int>& r ) {
                   for( int i = r.begin(); i != r.end(); ++i )
                   {
                      assert( presolvers[i]->runInRound( roundCounter ) );
                      results[i] = presolvers[i]->run( problem, probUpdate, num,
                                                       reductions[i] );
                   }
                },
                tbb::simple_partitioner() );
#else
            for( int i = fastPresolvers.first; i != fastPresolvers.second; ++i )
            {
               assert( presolvers[i]->runInRound( roundCounter ) );
               results[i] = presolvers[i]->run( problem, probUpdate, num,
                                                reductions[i] );
            }
#endif
            break;
         case 1:
#ifdef PARALLEL_MEDIUM_PRESOLVERS
            tbb::parallel_for(
                tbb::blocked_range<int>( mediumPresolvers.first,
                                         mediumPresolvers.second ),
                [&]( const tbb::blocked_range<int>& r ) {
                   for( int i = r.begin(); i != r.end(); ++i )
                   {
                      assert( presolvers[i]->runInRound( roundCounter ) );
                      results[i] = presolvers[i]->run( problem, probUpdate, num,
                                                       reductions[i] );
                   }
                },
                tbb::simple_partitioner() );
#else
            for( int i = mediumPresolvers.first; i != mediumPresolvers.second;
                 ++i )
            {
               assert( presolvers[i]->runInRound( roundCounter ) );
               results[i] = presolvers[i]->run( problem, probUpdate, num,
                                                reductions[i] );
            }
#endif
            break;
         case 2:
#ifdef PARALLEL_EXHAUSTIVE_PRESOLVERS
            tbb::parallel_for(
                tbb::blocked_range<int>( exhaustivePresolvers.first,
                                         exhaustivePresolvers.second ),
                [&]( const tbb::blocked_range<int>& r ) {
                   for( int i = r.begin(); i != r.end(); ++i )
                   {
                      assert( presolvers[i]->runInRound( roundCounter ) );
                      results[i] = presolvers[i]->run( problem, probUpdate, num,
                                                       reductions[i] );
                   }
                },
                tbb::simple_partitioner() );
#else
            for( int i = exhaustivePresolvers.first;
                 i != exhaustivePresolvers.second; ++i )
            {
               assert( presolvers[i]->runInRound( roundCounter ) );
               results[i] = presolvers[i]->run( problem, probUpdate, num,
                                                reductions[i] );
            }
#endif
            break;
         case 3:
            assert( false );
         }

         if( roundCounter == 0 )
         {
            probUpdate.clearChangeInfo();
            lastRoundReduced = false;
         }

         Statistics oldstats = stats;

         // evaluate results
         result.status = evaluateResults();
         switch( result.status )
         {
         case PresolveStatus::kUnbndOrInfeas:
            // in case of unbounded or infeasible results we return immediately
            printPresolversStats();
            Message::debug(
                this,
                "[{}:{}] presolvers detected infeasibility or unboundedness\n",
                __FILE__, __LINE__ );
            return result;
         case PresolveStatus::kUnbounded:
            // in case of unbounded or infeasible results we return immediately
            printPresolversStats();
            Message::debug( this,
                            "[{}:{}] presolvers detected unbounded problem\n",
                            __FILE__, __LINE__ );
            return result;
         case PresolveStatus::kInfeasible:
            // in case of unbounded or infeasible results we return immediately
            printPresolversStats();
            Message::debug( this, "[{}:{}] presolvers detected infeasibility\n",
                            __FILE__, __LINE__ );
            return result;
         case PresolveStatus::kUnchanged:
            // printRoundStats( true );
            //++roundCounter;
            abort = updateRoundCounter( problem, probUpdate,
                                        ( stats - oldstats ), timer, true );
            break;
         case PresolveStatus::kReduced:
            // problem reductions where found by at least one presolver
            if( !applyPresolversReductions( probUpdate ) )
            {
               result.status = PresolveStatus::kInfeasible;
               return result;
            }

            abort = updateRoundCounter( problem, probUpdate,
                                        ( stats - oldstats ), timer );

            // end round
            finishRound( probUpdate );

#if 0
         assertCorrectness( problem, num );
#endif
         }

      } while( !abort );

      if( stats.ntsxapplied > 0 || stats.nboundchgs > 0 ||
          stats.ncoefchgs > 0 || stats.ndeletedcols > 0 ||
          stats.ndeletedrows > 0 || stats.nsidechgs > 0 )
      {
         result.status = probUpdate.trivialPresolve();

         if( result.status == PresolveStatus::kInfeasible ||
             result.status == PresolveStatus::kUnbndOrInfeas ||
             result.status == PresolveStatus::kUnbounded )
            return result;

         probUpdate.clearStates();
      }

      printPresolversStats();

      if( DependentRows<REAL>::Enabled &&
          ( presolveOptions.detectlindep == 2 ||
            ( problem.getNumIntegralCols() == 0 &&
              presolveOptions.detectlindep == 1 ) ) )
      {
         ConstraintMatrix<REAL>& consMatrix = problem.getConstraintMatrix();
         Vec<int> equations;

         equations.reserve( problem.getNRows() );
         size_t eqnnz = 0;

         for( int i = 0; i != problem.getNRows(); ++i )
         {
            if( rflags[i].test( RowFlag::kRedundant ) ||
                !rflags[i].test( RowFlag::kEquation ) )
               continue;

            equations.push_back( i );
            eqnnz += rowsize[i] + 1;
         }

         if( !equations.empty() )
         {
            DependentRows<REAL> depRows( equations.size(), problem.getNCols(),
                                         eqnnz );

            for( size_t i = 0; i != equations.size(); ++i )
               depRows.addRow( i, consMatrix.getRowCoefficients( equations[i] ),
                               REAL( rhsVals[equations[i]] ) );

            Vec<int> dependentEqs;
            double factorTime = 0.0;
            msg.info( "found {} equations, checking for linear dependency\n",
                      equations.size() );
            {
               Timer t{ factorTime };
               dependentEqs = depRows.getDependentRows( msg, num );
            }
            msg.info(
                "{} equations are redundant, factorization took {} seconds\n",
                dependentEqs.size(), factorTime );

            if( !dependentEqs.empty() )
            {
               for( int dependentEq : dependentEqs )
               {
                  probUpdate.markRowRedundant( equations[dependentEq] );
               }
               probUpdate.flush();
            }
         }

         if( presolveOptions.dualreds == 2 )
         {
            Vec<int> freeCols;
            freeCols.reserve( problem.getNCols() );
            size_t freeColNnz = 0;

            const Vec<ColFlags>& cflags = problem.getColFlags();
            const Vec<int>& colsize = problem.getColSizes();
            const Vec<REAL>& obj = problem.getObjective().coefficients;
            const Vec<REAL>& lbs = problem.getLowerBounds();
            const Vec<REAL>& ubs = problem.getUpperBounds();

            for( int col = 0; col != problem.getNCols(); ++col )
            {
               if( cflags[col].test( ColFlag::kInactive, ColFlag::kIntegral ) ||
                   !cflags[col].test( ColFlag::kLbInf ) ||
                   !cflags[col].test( ColFlag::kUbInf ) )
                  continue;

               freeCols.push_back( col );
               freeColNnz += colsize[col] + 1;
            }

            if( !freeCols.empty() )
            {
               DependentRows<REAL> depRows( freeCols.size(), problem.getNRows(),
                                            freeColNnz );

               for( size_t i = 0; i != freeCols.size(); ++i )
                  depRows.addRow(
                      i, consMatrix.getColumnCoefficients( freeCols[i] ),
                      obj[freeCols[i]] );

               Vec<int> dependentFreeCols;
               double factorTime = 0.0;
               msg.info(
                   "found {} free columns, checking for linear dependency\n",
                   freeCols.size(), freeColNnz );

               {
                  Timer t{ factorTime };
                  dependentFreeCols = depRows.getDependentRows( msg, num );
               }

               msg.info( "{} free columns are redundant, factorization took {} "
                         "seconds\n",
                         dependentFreeCols.size(), factorTime );

               if( !dependentFreeCols.empty() )
               {
                  for( int dependentFreeCol : dependentFreeCols )
                     probUpdate.fixCol( freeCols[dependentFreeCol], 0 );

                  probUpdate.flush();
               }
            }
         }
      }

      // finally compress problem fully and release excess storage even if
      // problem was not reduced
      probUpdate.compress( true );

      // check whether problem was reduced
      if( stats.ntsxapplied > 0 || stats.nboundchgs > 0 ||
          stats.ncoefchgs > 0 || stats.ndeletedcols > 0 ||
          stats.ndeletedrows > 0 || stats.nsidechgs > 0 )
      {
         if( presolveOptions.boundrelax && problem.getNumIntegralCols() == 0 )
         {
            int nremoved;
            int nnewfreevars;
            // todo check if lp solver is simplex solver / add options
            std::tie( nremoved, nnewfreevars ) =
                probUpdate.removeRedundantBounds();
            if( nremoved != 0 )
               msg.info( "removed {} redundant column bounds, got {} new free "
                         "variables\n",
                         nremoved, nnewfreevars );
         }

         bool detectComponents = presolveOptions.componentsmaxint != -1;

         if( !lpSolverFactory && problem.getNumContinuousCols() != 0 )
            detectComponents = false;

         if( !mipSolverFactory && problem.getNumIntegralCols() != 0 )
            detectComponents = false;

         if( problem.getNCols() == 0 )
            detectComponents = false;

         if( detectComponents )
         {
            assert( problem.getNCols() != 0 && problem.getNRows() != 0 );
            Components components;

            int ncomponents = components.detectComponents( problem );

            if( ncomponents > 1 )
            {
               const Vec<ComponentInfo>& compInfo =
                   components.getComponentInfo();

               msg.info( "found {} disconnected components\n", ncomponents );
               msg.info(
                   "largest component has {} cols ({} int., {} cont.) and "
                   "{} nonzeros\n",
                   compInfo[ncomponents - 1].nintegral +
                       compInfo[ncomponents - 1].ncontinuous,
                   compInfo[ncomponents - 1].nintegral,
                   compInfo[ncomponents - 1].ncontinuous,
                   compInfo[ncomponents - 1].nnonz );

               Solution<REAL> solution;
               solution.primal.resize( problem.getNCols() );
               Vec<uint8_t> componentSolved( ncomponents );

               tbb::parallel_for(
                   tbb::blocked_range<int>( 0, ncomponents - 1 ),
                   [this, &components, &solution, &problem, &result, &compInfo,
                    &componentSolved,
                    &timer]( const tbb::blocked_range<int>& r ) {
                      for( int i = r.begin(); i != r.end(); ++i )
                      {
                         if( compInfo[i].nintegral == 0 )
                         {
                            std::unique_ptr<SolverInterface<REAL>> solver =
                                lpSolverFactory->newSolver(
                                    VerbosityLevel::kQuiet );

                            solver->setUp( problem,
                                           result.postsolve.origrow_mapping,
                                           result.postsolve.origcol_mapping,
                                           components, compInfo[i] );

                            if( presolveOptions.tlim !=
                                std::numeric_limits<double>::max() )
                            {
                               double tlim =
                                   presolveOptions.tlim - timer.getTime();
                               if( tlim <= 0 )
                                  break;
                               solver->setTimeLimit( tlim );
                            }

                            solver->solve();

                            SolverStatus status = solver->getStatus();

                            if( status == SolverStatus::kOptimal )
                            {
                               if( solver->getSolution( components,
                                                        compInfo[i].componentid,
                                                        solution ) )
                                  componentSolved[compInfo[i].componentid] =
                                      true;
                            }
                         }
                         else if( compInfo[i].nintegral <=
                                  presolveOptions.componentsmaxint )
                         {
                            std::unique_ptr<SolverInterface<REAL>> solver =
                                mipSolverFactory->newSolver(
                                    VerbosityLevel::kQuiet );

                            solver->setGapLimit( 0 );
                            solver->setNodeLimit(
                                problem.getConstraintMatrix().getNnz() /
                                std::max( compInfo[i].nnonz, 1 ) );

                            solver->setUp( problem,
                                           result.postsolve.origrow_mapping,
                                           result.postsolve.origcol_mapping,
                                           components, compInfo[i] );

                            if( presolveOptions.tlim !=
                                std::numeric_limits<double>::max() )
                            {
                               double tlim =
                                   presolveOptions.tlim - timer.getTime();
                               if( tlim <= 0 )
                                  break;
                               solver->setTimeLimit( tlim );
                            }

                            solver->solve();

                            SolverStatus status = solver->getStatus();

                            if( status == SolverStatus::kOptimal )
                            {
                               if( solver->getSolution( components,
                                                        compInfo[i].componentid,
                                                        solution ) )
                                  componentSolved[compInfo[i].componentid] =
                                      true;
                            }
                         }
                      }
                   },
                   tbb::simple_partitioner() );

               int nsolved = 0;

               int oldndelcols = stats.ndeletedcols;
               int oldndelrows = stats.ndeletedrows;

               auto& lbs = problem.getLowerBounds();
               auto& ubs = problem.getUpperBounds();
               for( int i = 0; i != ncomponents; ++i )
               {
                  if( componentSolved[i] )
                  {
                     ++nsolved;

                     const int* compcols = components.getComponentsCols( i );
                     int numcompcols = components.getComponentsNumCols( i );

                     for( int j = 0; j != numcompcols; ++j )
                     {
                        lbs[compcols[j]] = solution.primal[compcols[j]];
                        ubs[compcols[j]] = solution.primal[compcols[j]];
                        probUpdate.markColFixed( compcols[j] );
                     }

                     const int* comprows = components.getComponentsRows( i );
                     int numcomprows = components.getComponentsNumRows( i );

                     for( int j = 0; j != numcomprows; ++j )
                     {
                        probUpdate.markRowRedundant( comprows[j] );
                     }
                  }
               }

               if( nsolved != 0 )
               {
                  if( probUpdate.flush() == PresolveStatus::kInfeasible )
                     assert( false );

                  probUpdate.compress();

                  msg.info(
                      "solved {} components: {} cols fixed, {} rows deleted\n",
                      nsolved, stats.ndeletedcols - oldndelcols,
                      stats.ndeletedrows - oldndelrows );
               }
            }
         }

         msg.info( "presolved problem:\n" );
         msg.info( "  rows:     {}\n", problem.getNRows() );
         msg.info( "  columns:  {} ({} int., {} cont.)\n", problem.getNCols(),
                   problem.getNumIntegralCols(),
                   problem.getNumContinuousCols() );
         msg.info( "  nonzeros: {}\n", problem.getConstraintMatrix().getNnz() );

         result.status = PresolveStatus::kReduced;
         result.postsolve.getChecker().setReducedProblem( problem );
         return result;
      }

      msg.info( "presolved problem:\n" );
      msg.info( "  rows:     {}\n", problem.getNRows() );
      msg.info( "  columns:  {} ({} int., {} cont.)\n", problem.getNCols(),
                problem.getNumIntegralCols(), problem.getNumContinuousCols() );
      msg.info( "  nonzeros: {}\n", problem.getConstraintMatrix().getNnz() );

      // problem was not changed
      result.status = PresolveStatus::kUnchanged;
      return result;
   } );
}

} // namespace papilo

#endif
