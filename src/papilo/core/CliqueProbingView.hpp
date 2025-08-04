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

#ifndef _PAPILO_CORE_CLIQUE_PROBING_VIEW_HPP_
#define _PAPILO_CORE_CLIQUE_PROBING_VIEW_HPP_

#include "papilo/Config.hpp"
#include "papilo/core/Problem.hpp"
#include "papilo/core/SingleRow.hpp"
#include "papilo/io/Message.hpp"
#include "papilo/misc/Array.hpp"
#include "papilo/misc/MultiPrecision.hpp"
#include "papilo/misc/Vec.hpp"
#include <memory>
#include <list>

namespace papilo
{

template <typename REAL>
struct CliqueProbingBoundChg
{
   REAL bound;
   unsigned int col : 31;
   unsigned int upper : 1;
   int probing_col : 32;

   CliqueProbingBoundChg( bool upper_, int col_, REAL bound_, int probing_col_ )
   {
      this->upper = upper_ ? 1 : 0;
      this->col = static_cast<unsigned int>( col_ );
      this->bound = bound_;
      this->probing_col = ( probing_col_ );
   }
};


template <typename REAL>
struct CliqueProbingSubstitution
{
   REAL col2scale;
   REAL col2const;
   int col1;
   int col2;

   CliqueProbingSubstitution( int col1_, REAL col2scale_, int col2_, REAL col2const_ )
       : col2scale( col2scale_ ), col2const( col2const_ ), col1( col1_ ),
         col2( col2_ )
   {
   }
};

template <typename REAL>
class CliqueProbingView
{
 public:
   CliqueProbingView( const Problem<REAL>& problem, const Num<REAL>& num );

   void
   setMinIntDomRed( const REAL& value )
   {
      this->minintdomred = value;
   }

   void
   setMinContDomRed( const REAL& value )
   {
      this->mincontdomred = value;
   }

   void
   reset();

#ifdef PAPILO_TBB
   void
   parallelProbe( const tbb::blocked_range<int>& r, bool& initbounds, std::list<std::pair<int, REAL>>& changed_clique_lbs_inds_vals_thread_local, 
      std::list<std::pair<int,REAL>>& changed_clique_ubs_inds_vals_thread_local, Vec<std::pair<int,int>>& lb_implications_thread_local,
      Vec<std::pair<int,int>>& ub_implications_thread_local, Vec<int>& fix_to_zero_thread_local, bool& cliqueEquation, const Vec<int>& indices,
      const int& clique, const Vec<int>& binary_indices, const int& len )
   {
      assert( indices.size() > 0 );
      assert( r.begin() >= -static_cast<int>(!cliqueEquation) );
      assert( r.end() <= static_cast<int>(indices.size()) );
      probingClique = clique;
      cliqueind = indices;
      cliquelen = len;
      binary_inds = binary_indices;
      assert( cliqueind.size() > 0 );
      assert( cliquelen == static_cast<int>(cliqueind.size()) );
      //////std::cout<<"\nBinary inds at start of parallel clique probing: " << static_cast<int>(binary_inds.end() - binary_inds.begin());
      assert( ub_implications_thread_local.size() == binary_inds.size() );
      assert( lb_implications_thread_local.size() == binary_inds.size() );
      for( int i = r.begin(); i < r.end(); ++i )
      {
         //////////std::cout<<"\nTest1\n";
         if( i == -1 )
         {
            assert( cliqueEquation == false && indices.size() > 0 );
            setProbingColumn(-1);
            propagateDomains();
            if( isInfeasible() )
            {
               cliqueEquation = true;
               ////std::cout<<"\n\nTurned row " << clique << " into equation due to infeasibility.";
               ////std::cout.flush();
            }
            else
            {
               cliqueEquation = false;
               ////std::cout<<"\n\nTurned row " << clique << " into nonequation due to feasibility.";
               ////std::cout.flush();
               for( int var = 0; var != static_cast<int>(probing_lower_bounds.size()); ++var )
               {
                  //////////std::cout<<"\nTest2\n";
                  if( num.isGT(probing_lower_bounds[var], problem.getLowerBounds()[var]) )
                  {
                     changed_clique_lbs_inds_vals_thread_local.emplace_back(var, probing_lower_bounds[var]);
                     ////////std::cout<<"\nProbing ";
                     ////////std::cout<<" all on zero ";
                     ////////std::cout<<"yields stronger lower bounds ";
                     ////////std::cout<<probing_lower_bounds[var];
                     ////////std::cout<<" ";
                     ////////std::cout<<problem.getLowerBounds()[var];
                     ////////std::cout<<" for variable ";
                     ////////std::cout<<var;
                  }
                  if( num.isLT(probing_upper_bounds[var], problem.getUpperBounds()[var]) )
                  {
                     changed_clique_ubs_inds_vals_thread_local.emplace_back(var, probing_upper_bounds[var]);
                     ////////std::cout<<"\nProbing ";
                     ////////std::cout<<" all on zero ";
                     ////////std::cout<<"yields stronger upper bounds ";
                     ////////std::cout<<probing_upper_bounds[var];
                     ////////std::cout<<" ";
                     ////////std::cout<<problem.getUpperBounds()[var];
                     ////////std::cout<<" for variable ";
                     ////////std::cout<<var;
                  }
               }
               initbounds = true;
            }
            for( unsigned int ind = 0; ind !=  binary_inds.size() ; ++ind )
            {
               //////////std::cout<<"\nTest3\n";
               assert( ind < binary_inds.size() );
               assert( binary_inds[ind] < static_cast<int>(probing_lower_bounds.size()) );
               if( num.isEq( 1.0, probing_lower_bounds[binary_inds[ind]] ) )
               {
                  assert( ind < lb_implications_thread_local.size() );
                  lb_implications_thread_local[ind].first += 1;
                  lb_implications_thread_local[ind].second = -1;
               }
               assert( ind < binary_inds.size());
               assert( binary_inds[ind] < static_cast<int>(probing_upper_bounds.size()) );
               if( num.isEq( 0.0, probing_upper_bounds[binary_inds[ind]] ) )
               {
                  assert( ind < ub_implications_thread_local.size() );
                  ub_implications_thread_local[ind].first += 1;
                  ub_implications_thread_local[ind].second = -1;
               }
            }
            reset();
         }
         else
         {
            assert( indices.size() > 0 && i >= 0 && i < static_cast<int>(indices.size()) );
            setProbingColumn(i);
            propagateDomains();
            if( isInfeasible() )
            {
               fix_to_zero_thread_local.emplace_back( indices[i] );
               reset();
               continue;
            }
            for( unsigned int ind = 0; ind !=  binary_inds.size() ; ++ind )
            {
               //////////std::cout<<"\nTest4\n";
               assert( ind < binary_inds.size() );
               assert( binary_inds[ind] < static_cast<int>(probing_lower_bounds.size()) );
               if( num.isEq( 1.0, probing_lower_bounds[binary_inds[ind]] ) )
               {
                  assert( ind < lb_implications_thread_local.size() );
                  lb_implications_thread_local[ind].first += 1;
                  lb_implications_thread_local[ind].second = probingCol;
               }
               assert( ind < binary_inds.size());
               assert( binary_inds[ind] < static_cast<int>(probing_upper_bounds.size()) );
               if( num.isEq( 0.0, probing_upper_bounds[binary_inds[ind]] ) )
               {
                  assert( ind < ub_implications_thread_local.size() );
                  ub_implications_thread_local[ind].first += 1;
                  ub_implications_thread_local[ind].second = probingCol;
               }
            }
            //found new global bounds
            if( !initbounds )
            {
               for( unsigned int var = 0; var != probing_lower_bounds.size(); ++var )
               {
                  //////////std::cout<<"\nTest5\n";
                  if( num.isGT( probing_lower_bounds[var], problem.getLowerBounds()[var] ) )
                  {
                     changed_clique_lbs_inds_vals_thread_local.emplace_back(std::pair<int,REAL> {var, probing_lower_bounds[var] } );
                     ////////std::cout<<"\nProbing ";
                     ////////std::cout<<cliqueind[i];
                     ////////std::cout<<" on one and rest on zero ";
                     ////////std::cout<<"yields stronger lower bounds ";
                     ////////std::cout<<probing_lower_bounds[var];
                     ////////std::cout<<" ";
                     ////////std::cout<<problem.getLowerBounds()[var];
                     ////////std::cout<<" for variable ";
                     ////////std::cout<<var;
                  }
                  if( num.isLT( probing_upper_bounds[var], problem.getUpperBounds()[var] ) )
                  {
                     changed_clique_ubs_inds_vals_thread_local.emplace_back(std::pair<int,REAL> {var, probing_upper_bounds[var] } );
                     ////////std::cout<<"\nProbing ";
                     ////////std::cout<<cliqueind[i];
                     ////////std::cout<<" on one and rest on zero ";
                     ////////std::cout<<"yields stronger upper bounds ";
                     ////////std::cout<<probing_upper_bounds[var];
                     ////////std::cout<<" ";
                     ////////std::cout<<problem.getUpperBounds()[var];
                     ////////std::cout<<" for variable ";
                     ////////std::cout<<var;
                  }
               }
               initbounds = true;
               reset();
               continue;
            }
            
            typename std::list<std::pair<int,REAL>>::iterator ind = changed_clique_lbs_inds_vals_thread_local.begin();
            while( ind != changed_clique_lbs_inds_vals_thread_local.end() )
            {
               //////////std::cout<<"\nTest6\n";
               if( num.isLT( probing_lower_bounds[(*ind).first], (*ind).second ) )
               {
                  ////////std::cout<<"\nProbing ";
                  ////////std::cout<<cliqueind[i];
                  ////////std::cout<<" on one and rest on zero ";
                  ////////std::cout<<"yields weaker lower bounds ";
                  ////////std::cout<<probing_lower_bounds[(*ind).first];
                  ////////std::cout<<" ";
                  ////////std::cout<<(*ind).second;
                  ////////std::cout<<" for variable ";
                  ////////std::cout<<(*ind).first;
                  (*ind).second = probing_lower_bounds[(*ind).first];
               }
               if( num.isLE(probing_lower_bounds[(*ind).first], problem.getLowerBounds()[(*ind).first] ) )
               {
                  ind = changed_clique_lbs_inds_vals_thread_local.erase(ind);
               }
               else
                  std::advance(ind, 1);
            }
            ind = changed_clique_ubs_inds_vals_thread_local.begin();
            while( ind != changed_clique_ubs_inds_vals_thread_local.end() )
            {
               //////////std::cout<<"\nTest7\n";
               if( num.isGT( probing_upper_bounds[(*ind).first], (*ind).second ) )
               {
                  ////////std::cout<<"\nProbing ";
                  ////////std::cout<<cliqueind[i];
                  ////////std::cout<<" on one and rest on zero ";
                  ////////std::cout<<"yields weaker upper bounds ";
                  ////////std::cout<<probing_upper_bounds[(*ind).first];
                  ////////std::cout<<" ";
                  ////////std::cout<<(*ind).second;
                  ////////std::cout<<" for variable ";
                  ////////std::cout<<(*ind).first;
                  (*ind).second = probing_upper_bounds[(*ind).first];
               }
               if( num.isGE(probing_upper_bounds[(*ind).first], problem.getUpperBounds()[(*ind).first] ) )
               {
                  ind = changed_clique_ubs_inds_vals_thread_local.erase(ind);
               }
               else
                  std::advance(ind, 1);
            }
            reset();
         }
      }
      assert( cliqueind.size() > 0 );
   }
#endif

   std::pair<bool,bool>
   probeClique( const int clique, const int*& indices, const int len, const Vec<int>& binary_indices, 
      bool equation, Array<std::atomic_int>& probing_scores, const Vec<int>& colsize, const Vec<int>& colperm,
      Vec<int>& nprobed )
   {
      /*std::cout<<"\nStarting Probing of Clique.";
      std::cout.flush();*/
      if( alreadyinitialized == true )
      {
         std::cout<<"\nERROR WITH RESET: " << probingClique << " " << clique;
         std::cout.flush();
      }
      assert( alreadyinitialized = false );
      alreadyinitialized = true;
      binary_inds = binary_indices;
      fewreductions = false;
      probingClique = clique;
      for( int ind = 0; ind < len; ++ind )
      {
         cliqueind.emplace_back( indices[ind] );
      }
      const Vec<int> checkInd = cliqueind;
      const int checkLen = len;
      ////std::cout<<"\nLen and cliqueind size: " << len <<" " <<static_cast<int>(cliqueind.size());
      ////std::cout.flush();
      assert( cliqueind.size() > 0 );
      assert(len == static_cast<int>(cliqueind.size()));

      /*pdqsort( cliqueind.begin(), cliqueind.end(),
            [&probing_scores, &colsize, &colperm, &nprobed]( int col1, int col2 )
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
               return !(s1 > s2 || ( s1 == s2 && colperm[col1] < colperm[col2] ));
            } );*/

      cliquelen = len;
      ////std::cout<<"\nLen and cliqueind size after sort: " << len <<" " <<static_cast<int>(cliqueind.size());
      ////std::cout.flush();
      assert(len == static_cast<int>(cliqueind.size()));
      lb_implications.reserve( static_cast<int>(binary_inds.size()) );
      ub_implications.reserve( static_cast<int>( binary_inds.size() ) );
      for( unsigned int ind = 0; ind != binary_inds.size(); ++ind )
      {
         lb_implications.emplace_back( 0, -1 );
         ub_implications.emplace_back( 0, -1 );
      }
      assert(changed_clique_ubs_inds_vals.empty());
      assert(changed_clique_lbs_inds_vals.empty());
      equationBefore = equation;
      cliqueEquation = equation;
      //if( equation )
      //{//std::cout<<"\n\nRow " << clique << " is equation due to presets: " << problem.getConstraintMatrix().getLeftHandSides()[clique] 
      //<< " " << problem.getConstraintMatrix().getRightHandSides()[clique];}
      //std::cout.flush();
      bool initbounds = false;
#ifdef PAPILO_TBB
      tbb::combinable<Vec<std::pair<int,int>>> lb_implications_thread;
      tbb::combinable<Vec<std::pair<int,int>>> ub_implications_thread;
      tbb::combinable<std::pair<std::list<std::pair<int,REAL>>,bool>> changed_clique_ubs_inds_vals_initbounds_thread;
      tbb::combinable<std::pair<std::list<std::pair<int,REAL>>,bool>> changed_clique_lbs_inds_vals_initbounds_thread;
      tbb::combinable<Vec<int>> fix_to_zero_thread;
      tbb::combinable<int> numprobings;
      Vec<int> fix_to_zero_combined;
      Vec<std::pair<int,int>> lb_implications_combined;
      Vec<std::pair<int,int>> ub_implications_combined;
      std::list<std::pair<int,REAL>> changed_clique_lbs_inds_vals_combined;
      std::list<std::pair<int,REAL>> changed_clique_ubs_inds_vals_combined;
      for( unsigned int ind = 0; ind != binary_inds.size(); ++ind )
      {
         lb_implications_combined.emplace_back( 0, -1 );
         ub_implications_combined.emplace_back( 0, -1 );
      }
      
      assert( ub_implications_combined.size() == binary_inds.size() );
      assert( lb_implications_combined.size() == binary_inds.size() );
      assert( cliqueind.size() > 0 );

      assert( checkInd.size() == cliqueind.size() );
      assert( checkLen == cliquelen );

      int batchstart = -(!equation);
      int batchend = std::min( batchstart + 24, len );
      while( batchstart != len )
      {
         
         assert( checkInd.size() == cliqueind.size() );
         assert( checkLen == cliquelen );

         if( ( static_cast<int>(changed_clique_lbs_inds_vals_combined.size())
             + static_cast<int>(changed_clique_ubs_inds_vals_combined.size()) - cliquelen + static_cast<int>(batchstart) ) 
             < cliquelen * cliquereductionfactor && initbounds )
         {     
            //fewreductions = true;
            //return { false, cliqueEquation && !equationBefore } ;
         }
         assert( cliqueind.size() > 0 );
         tbb::parallel_for( tbb::blocked_range<int>( batchstart, batchend ),
            [&]( const tbb::blocked_range<int>& r )
            {

               
               assert( checkInd.size() == cliqueind.size() );
               assert( checkLen == cliquelen );

               assert( cliqueind.size() > 0 );
               Vec<int> localcliqueind = cliqueind;
               ////////std::cout << "Thread ID: " << std::this_thread::get_id() << " for range [" << r.begin() << ", " << r.end() << ")\n";
               if( ub_implications_thread.local().size() != binary_inds.size() 
                || lb_implications_thread.local().size() != binary_inds.size()  )
               {
                  if( lb_implications_thread.local().size() != 0 || ub_implications_thread.local().size() != 0 )
                  {
                     ////////std::cout<<"\n";
                     ////////std::cout<<"\n";
                     ////////std::cout<<lb_implications_thread.local().size();
                     ////////std::cout<<"\n";
                     ////////std::cout<<ub_implications_thread.local().size();
                     ////////std::cout<<"\n";
                     ////////std::cout<<binary_inds.size();
                     ////////std::cout<<"\n";
                     ////////std::cout<<"\n";
                  }
                  assert( lb_implications_thread.local().size() == 0 );
                  assert( ub_implications_thread.local().size() == 0 );
                  lb_implications_thread.local() = lb_implications_combined;
                  ub_implications_thread.local() = ub_implications_combined;
               }
               bool initbounds_thread_local = initbounds;
               if( !changed_clique_lbs_inds_vals_initbounds_thread.local().second )
               {
                  changed_clique_lbs_inds_vals_initbounds_thread.local().first = changed_clique_lbs_inds_vals_combined;
                  changed_clique_ubs_inds_vals_initbounds_thread.local().first = changed_clique_ubs_inds_vals_combined;
               }
               else
                  initbounds_thread_local = true;
               //assert( cliqueind.size() > 0 );
               CliqueProbingView<REAL> local_clique_probing( problem, num );
               //assert( cliqueind.size() > 0 );
               local_clique_probing.setMinContDomRed( mincontdomred );
               //assert( cliqueind.size() > 0 );

               assert( ub_implications_thread.local().size() == binary_inds.size() );
               assert( lb_implications_thread.local().size() == binary_inds.size() );
               
               assert( ub_implications_combined.size() == binary_inds.size() );
               assert( lb_implications_combined.size() == binary_inds.size() );
               //assert( cliqueind.size() > 0 );

               assert( localcliqueind.size() > 0 );
               
               assert( checkInd.size() == cliqueind.size() );
               assert( checkLen == cliquelen );

               local_clique_probing.parallelProbe( r, initbounds_thread_local, changed_clique_lbs_inds_vals_initbounds_thread.local().first, 
               changed_clique_ubs_inds_vals_initbounds_thread.local().first, lb_implications_thread.local(),
                  ub_implications_thread.local(), fix_to_zero_thread.local(), cliqueEquation, localcliqueind,
                  clique, binary_inds, cliquelen );
               numprobings.local() += 1;

               if( checkInd.size() != cliqueind.size() )
               {
                  std::cout<<"\nERROR in clique " << clique;
                  std::cout.flush();
               }
               assert( checkInd.size() == cliqueind.size() );
               assert( checkLen == cliquelen );
               assert( localcliqueind.size() > 0 );
               changed_clique_lbs_inds_vals_initbounds_thread.local().second = initbounds_thread_local;
               changed_clique_ubs_inds_vals_initbounds_thread.local().second = initbounds_thread_local;
               ////////std::cout << "\nLocal thread lb list contents: ";
               /*for (const auto& pair : changed_clique_lbs_inds_vals_initbounds_thread.local().first) {
                  ////////std::cout << "(" << pair.first << "," << pair.second << ") ";
               }
               ////////std::cout << "\n";
               ////////std::cout << "\nLocal thread ub list contents: ";
               for (const auto& pair : changed_clique_ubs_inds_vals_initbounds_thread.local().first) {
                  ////////std::cout << "(" << pair.first << "," << pair.second << ") ";
               }*/
               ////////std::cout << "\n";
               //assert( cliqueind.size() > 0 );
            }
         );

         ////////std::cout << "=== AFTER PARALLEL_FOR ===\n";
         //int post_count = 0;
         /*changed_clique_lbs_inds_vals_initbounds_thread.combine_each([&](const auto& data) {
            ////////std::cout << "Post thread data #" << post_count++ << ", initbounds: " << data.second 
            //         << ", list size: " << data.first.size() << "\n";
         });

         //int count = 0;
         changed_clique_lbs_inds_vals_initbounds_thread.combine_each([&](const auto& data) {
            ////////std::cout << "combine_each call #" << count++ << ", list size: " << data.first.size() << "\n";
         });*/

         numpropagations += static_cast<int>(batchend) - static_cast<int>(batchstart);

         
      assert( checkInd.size() == cliqueind.size() );
      assert( checkLen == cliquelen );

         fix_to_zero_thread.combine_each([&](const std::vector<int>& fix_to_zero_local ) {
            fix_to_zero_combined.insert(fix_to_zero_combined.end(), fix_to_zero_local.begin(), fix_to_zero_local.end());
         });
         fix_to_zero_thread.clear();

         bool initlowerbounds = initbounds;

         
      assert( checkInd.size() == cliqueind.size() );
      assert( checkLen == cliquelen );
         //static int combine_call_count = 0;

         changed_clique_lbs_inds_vals_initbounds_thread.combine_each([&]( std::pair<std::list<std::pair<int,REAL>>,bool> changed_clique_lbs_inds_vals_initbounds_local ) 
         {
            /*++combine_call_count;
            ////////std::cout << "\n=== COMBINE CALL #" << combine_call_count << " ===\n";
            ////////std::cout << "initlowerbounds: " << initlowerbounds << ", local.second: " << changed_clique_lbs_inds_vals_initbounds_local.second << "\n";
            
            // Debug: Print contents of local list
            ////////std::cout << "Local list contents: ";
            for (const auto& pair : changed_clique_lbs_inds_vals_initbounds_local.first) {
               ////////std::cout << "(" << pair.first << "," << pair.second << ") ";
            }
            ////////std::cout << "\n";
            
            // Debug: Print contents of combined list BEFORE processing
            ////////std::cout << "Combined list BEFORE: ";
            for (const auto& pair : changed_clique_lbs_inds_vals_combined) {
               ////////std::cout << "(" << pair.first << "," << pair.second << ") ";
            }
            ////////std::cout << "\n";
            
            // Debug: Check if lists are sorted
            auto check_sorted = [](const std::list<std::pair<int,REAL>>& lst, const std::string& name) {
               bool sorted = true;
               int prev = -1;
               for (const auto& pair : lst) {
                  if (pair.first <= prev) {
                     ////////std::cout << "ERROR: " << name << " is NOT sorted! Found " << pair.first << " after " << prev << "\n";
                     sorted = false;
                  }
                  prev = pair.first;
               }
               if (sorted) ////////std::cout << name << " is properly sorted\n";
               return sorted;
            };
            
            check_sorted(changed_clique_lbs_inds_vals_initbounds_local.first, "Local list");
            check_sorted(changed_clique_lbs_inds_vals_combined, "Combined list");*/
            
            if( !initlowerbounds && changed_clique_lbs_inds_vals_initbounds_local.second )
            {
               ////////std::cout << "INITIALIZING combined list with local list\n";
               changed_clique_lbs_inds_vals_combined = changed_clique_lbs_inds_vals_initbounds_local.first;
               initlowerbounds = true;
            }
            else if( initlowerbounds && changed_clique_lbs_inds_vals_initbounds_local.second )
            {
               ////////std::cout << "INTERSECTING lists\n";
               typename std::list<std::pair<int,REAL>>::iterator ind_local = changed_clique_lbs_inds_vals_initbounds_local.first.begin();
               typename std::list<std::pair<int,REAL>>::iterator ind_combined = changed_clique_lbs_inds_vals_combined.begin();
               
               while( ind_combined != changed_clique_lbs_inds_vals_combined.end() )
               {
                  ////////std::cout << "Comparing: combined(" << (*ind_combined).first << "," << (*ind_combined).second << ")";
                  if (ind_local != changed_clique_lbs_inds_vals_initbounds_local.first.end()) {
                     ////////std::cout << " vs local(" << (*ind_local).first << "," << (*ind_local).second << ")";
                  } else {
                     ////////std::cout << " vs local(END)";
                  }
                  ////////std::cout << "\n";
                  
                  if( ind_local == changed_clique_lbs_inds_vals_initbounds_local.first.end() )
                  {
                     ////////std::cout << "Local list exhausted, removing remaining combined elements\n";
                     while( ind_combined != changed_clique_lbs_inds_vals_combined.end() )
                     {
                        ////////std::cout << "Throwing out lower bound reduction: " << (*ind_combined).first << " " << (*ind_combined).second << "\n";
                        ind_combined = changed_clique_lbs_inds_vals_combined.erase(ind_combined);
                     }
                  }
                  else if( (*ind_combined).first < (*ind_local).first )
                  {
                     ////////std::cout << "Combined element not in local, removing: " << (*ind_combined).first << " " << (*ind_combined).second << "\n";
                     ind_combined = changed_clique_lbs_inds_vals_combined.erase(ind_combined);
                  }
                  else if( (*ind_combined).first > (*ind_local).first )
                  {
                     ////////std::cout << "Local element not in combined, advancing local\n";
                     std::advance(ind_local, 1);
                  }
                  else if( num.isGT((*ind_combined).second, (*ind_local).second) )
                  {
                     ////////std::cout << "Weakening lower bound reduction: " << (*ind_combined).first << " " << (*ind_combined).second << " -> " << (*ind_local).second << "\n";
                     (*ind_combined).second = (*ind_local).second;
                     ind_local = changed_clique_lbs_inds_vals_initbounds_local.first.erase(ind_local);
                     std::advance(ind_combined, 1);
                  }
                  else
                  {
                     ////////std::cout << "Combined bound is better or equal, keeping it\n";
                     ind_local = changed_clique_lbs_inds_vals_initbounds_local.first.erase(ind_local);
                     std::advance(ind_combined, 1);
                  }
               }
            }
            else
            {
               ////////std::cout << "Skipping this thread's results (initbounds_local.second = false)\n";
            }
            
            // Debug: Print contents of combined list AFTER processing
            ////////std::cout << "Combined list AFTER: ";
            /*for (const auto& pair : changed_clique_lbs_inds_vals_combined) {
               ////////std::cout << "(" << pair.first << "," << pair.second << ") ";
            }*/
            ////////std::cout << "\n=== END COMBINE CALL #" << combine_call_count << " ===\n";
         });

         // Final verification
         ////////std::cout << "\n=== FINAL RESULT ===\n";
         ////////std::cout << "Final combined list: ";
         /*for (const auto& pair : changed_clique_lbs_inds_vals_combined) {
            ////////std::cout << "(" << pair.first << "," << pair.second << ") ";
         }*/
         ////////std::cout << "\n";

         bool initupperbounds = initbounds;
         //static int ub_combine_call_count = 0;
         
      assert( checkInd.size() == cliqueind.size() );
      assert( checkLen == cliquelen );

         changed_clique_ubs_inds_vals_initbounds_thread.combine_each([&]( std::pair<std::list<std::pair<int,REAL>>,bool> changed_clique_ubs_inds_vals_initbounds_local ) 
         {
            /*++ub_combine_call_count;
            ////////std::cout << "\n=== UB COMBINE CALL #" << ub_combine_call_count << " ===\n";
            ////////std::cout << "initupperbounds: " << initupperbounds << ", local.second: " << changed_clique_ubs_inds_vals_initbounds_local.second << "\n";
            
            // Debug: Print contents of local list
            ////////std::cout << "Local UB list contents: ";
            for (const auto& pair : changed_clique_ubs_inds_vals_initbounds_local.first) {
               ////////std::cout << "(" << pair.first << "," << pair.second << ") ";
            }
            ////////std::cout << "\n";
            
            // Debug: Print contents of combined list BEFORE processing
            ////////std::cout << "Combined UB list BEFORE: ";
            for (const auto& pair : changed_clique_ubs_inds_vals_combined) {
               ////////std::cout << "(" << pair.first << "," << pair.second << ") ";
            }
            ////////std::cout << "\n";
            
            // Debug: Check if lists are sorted
            auto check_sorted = [](const std::list<std::pair<int,REAL>>& lst, const std::string& name) {
               bool sorted = true;
               int prev = -1;
               for (const auto& pair : lst) {
                  if (pair.first <= prev) {
                     ////////std::cout << "ERROR: " << name << " is NOT sorted! Found " << pair.first << " after " << prev << "\n";
                     sorted = false;
                  }
                  prev = pair.first;
               }
               if (sorted) ////////std::cout << name << " is properly sorted\n";
               return sorted;
            };
            
            check_sorted(changed_clique_ubs_inds_vals_initbounds_local.first, "Local UB list");
            check_sorted(changed_clique_ubs_inds_vals_combined, "Combined UB list");*/
            
            if( !initupperbounds && changed_clique_ubs_inds_vals_initbounds_local.second )
            {
               ////////std::cout << "INITIALIZING combined UB list with local list\n";
               changed_clique_ubs_inds_vals_combined = changed_clique_ubs_inds_vals_initbounds_local.first;
               initupperbounds = true;
            }
            else if( initupperbounds && changed_clique_ubs_inds_vals_initbounds_local.second )
            {
               ////////std::cout << "INTERSECTING UB lists\n";
               typename std::list<std::pair<int,REAL>>::iterator ind_local = changed_clique_ubs_inds_vals_initbounds_local.first.begin();
               typename std::list<std::pair<int,REAL>>::iterator ind_combined = changed_clique_ubs_inds_vals_combined.begin();
               
               while( ind_combined != changed_clique_ubs_inds_vals_combined.end() )
               {
                  ////////std::cout << "Comparing UB: combined(" << (*ind_combined).first << "," << (*ind_combined).second << ")";
                  if (ind_local != changed_clique_ubs_inds_vals_initbounds_local.first.end()) {
                     ////////std::cout << " vs local(" << (*ind_local).first << "," << (*ind_local).second << ")";
                  } else {
                     ////////std::cout << " vs local(END)";
                  }
                  ////////std::cout << "\n";
                  
                  if( ind_local == changed_clique_ubs_inds_vals_initbounds_local.first.end() )
                  {
                     ////////std::cout << "Local UB list exhausted, removing remaining combined elements\n";
                     while( ind_combined != changed_clique_ubs_inds_vals_combined.end() )
                     {
                        ////////std::cout << "Throwing out upper bound reduction: " << (*ind_combined).first << " " << (*ind_combined).second << "\n";
                        ind_combined = changed_clique_ubs_inds_vals_combined.erase(ind_combined);   
                     }
                  }
                  else if( (*ind_combined).first < (*ind_local).first )
                  {
                     ////////std::cout << "Combined UB element not in local, removing: " << (*ind_combined).first << " " << (*ind_combined).second << "\n";
                     ind_combined = changed_clique_ubs_inds_vals_combined.erase(ind_combined);
                  }
                  else if( (*ind_combined).first > (*ind_local).first )
                  {
                     ////////std::cout << "Local UB element not in combined, advancing local\n";
                     std::advance(ind_local, 1);
                  }
                  else if( num.isLT((*ind_combined).second, (*ind_local).second) )
                  {
                     ////////std::cout << "Weakening upper bound reduction: " << (*ind_combined).first << " " << (*ind_combined).second << " -> " << (*ind_local).second << "\n";
                     (*ind_combined).second = (*ind_local).second;
                     ind_local = changed_clique_ubs_inds_vals_initbounds_local.first.erase(ind_local);
                     std::advance(ind_combined, 1);
                  }
                  else
                  {
                     ////////std::cout << "Combined UB bound is better or equal, keeping it\n";
                     ind_local = changed_clique_ubs_inds_vals_initbounds_local.first.erase(ind_local);
                     std::advance(ind_combined, 1);
                  }
               }
            }
            else
            {
               ////////std::cout << "Skipping this thread's UB results (initbounds_local.second = false)\n";
            }
            
            // Debug: Print contents of combined list AFTER processing
            ////////std::cout << "Combined UB list AFTER: ";
            /*for (const auto& pair : changed_clique_ubs_inds_vals_combined) {
               ////////std::cout << "(" << pair.first << "," << pair.second << ") ";
            }*/
            ////////std::cout << "\n=== END UB COMBINE CALL #" << ub_combine_call_count << " ===\n";
         });

         // Final verification
         ////////std::cout << "\n=== FINAL UB RESULT ===\n";
         ////////std::cout << "Final combined UB list: ";
         /*for (const auto& pair : changed_clique_ubs_inds_vals_combined) {
            ////////std::cout << "(" << pair.first << "," << pair.second << ") ";
         }*/
         ////////std::cout << "\n";

         initbounds = initupperbounds;
         /*
         ////////std::cout<<"\n";
         ////////std::cout<<batchstart;
         ////////std::cout<<"\n";
         ////////std::cout<<batchend;
         ////////std::cout<<"\n";
         ////////std::cout<<tbb::this_task_arena::max_concurrency();
         ////////std::cout<<"\n";
         */
         changed_clique_lbs_inds_vals_initbounds_thread.clear();
         changed_clique_ubs_inds_vals_initbounds_thread.clear();

         batchstart = batchend;
         batchend = std::min( batchstart + 24, len );
      }

      
      assert( checkInd.size() == cliqueind.size() );
      assert( checkLen == cliquelen );

      lb_implications_thread.combine_each([&](const std::vector<std::pair<int,int>>& lb_implications_local )
      {
         for( unsigned int ind = 0; ind != binary_inds.size(); ++ind )
         {
            //////////std::cout<<"\nTest13\n";
            lb_implications_combined[ind].first += lb_implications_local[ind].first;
            if( lb_implications_local[ind].second > lb_implications_combined[ind].second )
               lb_implications_combined[ind].second = lb_implications_local[ind].second;
         }
      });
      lb_implications_thread.clear();

      
      assert( checkInd.size() == cliqueind.size() );
      assert( checkLen == cliquelen );

      ub_implications_thread.combine_each([&](const std::vector<std::pair<int,int>>& ub_implications_local )
      {
         for( unsigned int ind = 0; ind != binary_inds.size(); ++ind )
         {
            //////////std::cout<<"\nTest14\n";
            ub_implications_combined[ind].first += ub_implications_local[ind].first;
            if( ub_implications_local[ind].second > ub_implications_combined[ind].second )
               ub_implications_combined[ind].second = ub_implications_local[ind].second;
         }
      });
      
      assert( checkInd.size() == cliqueind.size() );
      assert( checkLen == cliquelen );

      int totalnumprobings = 0;
      numprobings.combine_each([&](const int& numprobingslocal )
      {
         totalnumprobings+=numprobingslocal;
      });
      if( totalnumprobings != cliquelen + 1 - static_cast<int>(equationBefore) )
      {
         std::cout<< "\n"<< totalnumprobings<< " "<< cliquelen << " "<< cliqueind.size()<< " "<< equationBefore
         << " " << checkInd.size() << " " << checkLen << "\n";
         std::cout.flush();
      }
      assert( checkInd.size() == cliqueind.size() );
      assert( checkLen == cliquelen );
      assert( totalnumprobings == cliquelen + 1 - static_cast<int>(equationBefore) );

      ub_implications_thread.clear();

      
      changed_clique_lbs_inds_vals = changed_clique_lbs_inds_vals_combined;
      changed_clique_ubs_inds_vals = changed_clique_ubs_inds_vals_combined;
      ub_implications = ub_implications_combined;
      lb_implications = lb_implications_combined;
      fix_to_zero = fix_to_zero_combined;
      /*bool feas = fix_to_zero.end() - fix_to_zero.begin() == cliquelen && cliqueEquation;
      std::cout<<"\nFinished initial probing, Infeasibility: " << feas;
      std::cout.flush();*/
      return { fix_to_zero.end() - fix_to_zero.begin() == cliquelen && cliqueEquation, cliqueEquation && !equationBefore } ;
#else
      if(!equation)
      {
         reset();
         setProbingColumn(-1);

         cliqueEquation = false;

         propagateDomains();

         numpropagations += 1;
         if( isInfeasible() )
         {
            cliqueEquation = true;
         }
         else
         {
            for( int var = 0; var != static_cast<int>(probing_lower_bounds.size()); ++var )
            {
               if( num.isGT(probing_lower_bounds[var], problem.getLowerBounds()[var]) )
               {
                  changed_clique_lbs_inds_vals.emplace_back(var, probing_lower_bounds[var]);
               }
               if( num.isLT(probing_upper_bounds[var], problem.getUpperBounds()[var]) )
               {
                  changed_clique_ubs_inds_vals.emplace_back(var, probing_upper_bounds[var]);
               }
            }
            initbounds = true;
         }
         for( unsigned int ind = 0; ind !=  binary_inds.size() ; ++ind )
         {
            assert( ind < binary_inds.size() );
            assert( binary_inds[ind] < static_cast<int>(probing_lower_bounds.size()) );
            if( num.isEq( 1.0, probing_lower_bounds[binary_inds[ind]] ) )
            {
               assert( ind < lb_implications.size() );
               lb_implications[ind].first += 1;
               lb_implications[ind].second = -1;
            }
            assert( ind < binary_inds.size());
            assert( binary_inds[ind] < static_cast<int>(probing_upper_bounds.size()) );
            if( num.isEq( 0.0, probing_upper_bounds[binary_inds[ind]] ) )
            {
               assert( ind < ub_implications.size() );
               ub_implications[ind].first += 1;
               ub_implications[ind].second = -1;
            }
         }
         reset();
      }

      for( int i = 0; i < cliquelen; ++i )
      {
         if( ( static_cast<int>(changed_clique_lbs_inds_vals.size())
             + static_cast<int>(changed_clique_ubs_inds_vals.size()) - cliquelen + i ) < cliquelen * cliquereductionfactor
             && initbounds )
         {     
            fewreductions = true;
            return { false, cliqueEquation && !equationBefore } ;
         }
         reset();
         assert( probing_upper_bounds[cliqueind[i]] == 1.0 );
         assert( probing_lower_bounds[cliqueind[i]] == 0.0 );
         setProbingColumn( i );

         propagateDomains();

         numpropagations += 1;
         if( isInfeasible() )
         {
            fix_to_zero.emplace_back( probingCol );
            reset();
            continue;
         }
         for( unsigned int ind = 0; ind !=  binary_inds.size() ; ++ind )
         {
            assert( ind < binary_inds.size() );
            assert( binary_inds[ind] < static_cast<int>(probing_lower_bounds.size()) );
            if( num.isEq( 1.0, probing_lower_bounds[binary_inds[ind]] ) )
            {
               assert( ind < lb_implications.size() );
               lb_implications[ind].first += 1;
               lb_implications[ind].second = probingCol;
            }
            assert( ind < binary_inds.size());
            assert( binary_inds[ind] < static_cast<int>(probing_upper_bounds.size()) );
            if( num.isEq( 0.0, probing_upper_bounds[binary_inds[ind]] ) )
            {
               assert( ind < ub_implications.size() );
               ub_implications[ind].first += 1;
               ub_implications[ind].second = probingCol;
            }
         }
         //found new global bounds
         if( !initbounds )
         {
            for( unsigned int var = 0; var != probing_lower_bounds.size(); ++var )
            {
               if( num.isGT( probing_lower_bounds[var], problem.getLowerBounds()[var] ) )
               {
                  changed_clique_lbs_inds_vals.emplace_back(std::pair<int,REAL> {var, probing_lower_bounds[var] } );
               }
               if( num.isLT( probing_upper_bounds[var], problem.getUpperBounds()[var] ) )
               {
                  changed_clique_ubs_inds_vals.emplace_back(std::pair<int,REAL> {var, probing_upper_bounds[var] } );
               }
            }
            initbounds = true;
            reset();
            continue;
         }
         //check for substitutions
         typename std::list<std::pair<int,REAL>>::iterator ind = changed_clique_lbs_inds_vals.begin();
         while( ind != changed_clique_lbs_inds_vals.end() )
         {
            if( num.isLT( probing_lower_bounds[(*ind).first], (*ind).second ) )
            {
               (*ind).second = probing_lower_bounds[(*ind).first];
            }
            if( num.isLE(probing_lower_bounds[(*ind).first], problem.getLowerBounds()[(*ind).first] ) )
            {
               ind = changed_clique_lbs_inds_vals.erase(ind);
            }
            else
               std::advance(ind, 1);
         }
         ind = changed_clique_ubs_inds_vals.begin();
         while( ind != changed_clique_ubs_inds_vals.end() )
         {
            if( num.isGT( probing_upper_bounds[(*ind).first], (*ind).second ) )
            {
               (*ind).second = probing_upper_bounds[(*ind).first];
            }
            if( num.isGE(probing_upper_bounds[(*ind).first], problem.getUpperBounds()[(*ind).first] ) )
            {
               ind = changed_clique_ubs_inds_vals.erase(ind);
            }
            else
               std::advance(ind, 1);
         }
        reset();
       }

      return { fix_to_zero.end() - fix_to_zero.begin() == cliquelen && cliqueEquation, cliqueEquation && !equationBefore } ;
#endif
   }

   void
   setProbingColumn( int col )
   {
      if( col == -1 )
      {
         probingCol = cliqueind[0];
         for( int i = 0; i < cliquelen; ++i )
         {
            changeUb( cliqueind[i], 0.0 );
         }
      }
      else
      {
         if( col < 0 || col >= static_cast<int>(cliqueind.size()) )
         {
            std::cout<<"\n\nERROR: Bad Index: " << col;
            std::cout<<"\nClique length: " << cliqueind.size();
#ifdef PAPILO_TBB
            std::cout<<"\nTBB is on.";
#else
            std::cout<<"\nTBB is off.";
#endif
            std::cout.flush();
            assert(col >= 0 && col < static_cast<int>(cliqueind.size() ));
         }
         probingCol = cliqueind[col];
         changeLb( probingCol, 1.0 );
         for( int i = 0; i < cliquelen; ++i )
         {
            if( i == col )
               continue;
            changeUb( cliqueind[i], 0.0 );
         }
      }
   }

   void
   resetClique()
   {
    assert( alreadyinitialized == true );
    reset();
    changed_clique_ubs_inds_vals.clear();
    changed_clique_lbs_inds_vals.clear();
    assert(changed_clique_ubs_inds_vals.empty());
    assert(changed_clique_lbs_inds_vals.empty());
    lb_implications.clear();
    ub_implications.clear();
    fix_to_zero.clear();
    cliqueind.clear();
    probingClique = -1;
    cliquelen = -1;
    cliqueEquation = false;
    equationBefore = false;
    fewreductions = false;
    numpropagations = 0;

    alreadyinitialized = false
   }

   void
   activityChanged( ActivityChange actchange, int rowid,
                    RowActivity<REAL>& activity );
   void
   changeLb( int col, REAL newlb );

   void
   changeUb( int col, REAL newub );

   bool
   analyzeImplications();

   void
   propagateDomains();

   bool
   isInfeasible() const
   {
      return infeasible;
   }

   int
   getNumSubstitutions() const
   {
      return static_cast<int>( substitutions.size() );
   }

   Vec<CliqueProbingBoundChg<REAL>>&
   getProbingBoundChanges()
   {
      return boundChanges;
   }


   Vec<CliqueProbingSubstitution<REAL>>&
   getProbingSubstitutions()
   {
      return substitutions;
   }

   const Vec<CliqueProbingBoundChg<REAL>>&
   getProbingBoundChanges() const
   {
      return boundChanges;
   }

   const Vec<CliqueProbingSubstitution<REAL>>&
   getProbingSubstitutions() const
   {
      return substitutions;
   }

   int64_t
   getAmountOfWork() const
   {
      return amountofwork;
   }

   int
   getNumPropagations()
   {
      return numpropagations;
   }

   const Vec<REAL>&
   getProbingLowerBounds() const
   {
      return probing_lower_bounds;
   }

   const Vec<REAL>&
   getProbingUpperBounds() const
   {
      return probing_upper_bounds;
   }

   const Vec<ColFlags>&
   getProbingDomainFlags() const
   {
      return probing_domain_flags;
   }

   void
   clearResults()
   {
      amountofwork = 0;
      boundChanges.clear();
      substitutions.clear();
   }

 private:
   // reference to problem and numerics class
   const Problem<REAL>& problem;
   const Num<REAL>& num;
   REAL minintdomred;
   REAL mincontdomred;

   // datastructures used for probing
   Vec<int> changed_lbs;
   Vec<int> changed_ubs;
   Vec<int> changed_activities;
   Vec<REAL> probing_lower_bounds;
   Vec<REAL> probing_upper_bounds;
   Vec<ColFlags> probing_domain_flags;
   Vec<RowActivity<REAL>> probing_activities;

   std::list<std::pair<int,REAL>> changed_clique_lbs_inds_vals;
   std::list<std::pair<int,REAL>> changed_clique_ubs_inds_vals;

   Vec<int> binary_inds;
   Vec<std::pair<int,int>> lb_implications;
   Vec<std::pair<int,int>> ub_implications;
   Vec<int> fix_to_zero;

   Vec<int> prop_activities;
   Vec<int> next_prop_activities;

   bool infeasible;
   int round;
   int probingCol;
   int probingClique;
   bool probingValue;
   Vec<int> cliqueind;
   int cliquelen;
   bool cliqueEquation;
   bool equationBefore;

   bool alreadyinitialized;

   // datastructures for storing result of probing on one value
   Vec<CliqueProbingBoundChg<REAL>> otherValueImplications;
   bool otherValueInfeasible;
   Vec<bool> infeasibleAssignments;

   // results of probing and statistics
   Vec<CliqueProbingBoundChg<REAL>> boundChanges;
   Vec<CliqueProbingSubstitution<REAL>> substitutions;

   int64_t amountofwork;
   
   const int cliquereductionfactor = 3;
   bool fewreductions;
   int numpropagations = 0;
};

#ifdef PAPILO_USE_EXTERN_TEMPLATES
extern template class CliqueProbingView<double>;
extern template class CliqueProbingView<Quad>;
extern template class CliqueProbingView<Rational>;
#endif

template <typename REAL>
CliqueProbingView<REAL>::CliqueProbingView( const Problem<REAL>& problem_,
                                const Num<REAL>& num_ )
    : problem( problem_ ), num( num_ ),
      probing_lower_bounds( problem_.getLowerBounds() ),
      probing_upper_bounds( problem_.getUpperBounds() ),
      probing_domain_flags( problem_.getColFlags() ),
      probing_activities( problem_.getRowActivities() )
{
   round = -2;
   infeasible = false;
   amountofwork = 0;
   numpropagations = 0;
   probingCol = -1;
   probingClique = -1;
   probingValue = false;
   otherValueInfeasible = false;
   minintdomred = num.getFeasTol() * 1000;
   mincontdomred = 0.3;

   alreadyinitialized = false;
}

template <typename REAL>
void
CliqueProbingView<REAL>::reset()
{
   const Vec<int>& rowsize = problem.getConstraintMatrix().getRowSizes();

   const auto& orig_lbs = problem.getLowerBounds();
   for( int i : changed_lbs )
   {
      if( i < 0 )
      {
         int c = -i - 1;
         assert( !probing_domain_flags[c].test( ColFlag::kLbUseless ) &&
                 problem.getColFlags()[c].test( ColFlag::kLbUseless ) );
         probing_domain_flags[c].set( ColFlag::kLbUseless );
#ifndef NDEBUG
         probing_lower_bounds[c] = orig_lbs[c];
#endif
      }
      else
         probing_lower_bounds[i] = orig_lbs[i];
   }
   changed_lbs.clear();

   const auto& orig_ubs = problem.getUpperBounds();
   for( int i : changed_ubs )
   {
      if( i < 0 )
      {
         int c = -i - 1;
         assert( !probing_domain_flags[c].test( ColFlag::kUbUseless ) &&
                 problem.getColFlags()[c].test( ColFlag::kUbUseless ) );
         probing_domain_flags[c].set( ColFlag::kUbUseless );
#ifndef NDEBUG
         probing_upper_bounds[c] = orig_ubs[c];
#endif
      }
      else
         probing_upper_bounds[i] = orig_ubs[i];
   }
   changed_ubs.clear();

   const auto& orig_activities = problem.getRowActivities();
   for( int i : changed_activities )
   {
      amountofwork += rowsize[i];
      probing_activities[i] = orig_activities[i];
   }
   changed_activities.clear();

   // reset should result in original domains and activities
   assert( std::equal( orig_lbs.begin(), orig_lbs.end(),
                       probing_lower_bounds.begin() ) );
   assert( std::equal( orig_ubs.begin(), orig_ubs.end(),
                       probing_upper_bounds.begin() ) );
   assert( std::equal(
       orig_activities.begin(), orig_activities.end(),
       probing_activities.begin(),
       []( const RowActivity<REAL>& a, const RowActivity<REAL>& b ) {
          return a.ninfmax == b.ninfmax && a.ninfmin == b.ninfmin &&
                 a.min == b.min && a.max == b.max &&
                 a.lastchange == b.lastchange;
       } ) );

   round = -2;
   prop_activities.clear();
   next_prop_activities.clear();
   infeasible = false;
   probingCol = -1;
}

template <typename REAL>
void
CliqueProbingView<REAL>::activityChanged( ActivityChange actchange, int rowid,
                                    RowActivity<REAL>& activity )
{
   const auto& consMatrix = problem.getConstraintMatrix();
   const auto& lhs = consMatrix.getLeftHandSides();
   const auto& rhs = consMatrix.getRightHandSides();
   const auto& rflags = consMatrix.getRowFlags();

   // mark the lastchange fields with round values starting from -2
   // and counting backward By doing this we avoid that we need to
   // alter the lastchange field after copying the original activities
   // since they are always larger or equal to -1
   if( activity.lastchange > -2 )
      changed_activities.push_back(
          rowid ); // activity was changed for the first time

   if( activity.lastchange != round )
      next_prop_activities.push_back( rowid );

   activity.lastchange = round;

   // check if the updated activity is reliable or if it is zero relative to
   // the initial activity
   const RowActivity<REAL>& origactivity = problem.getRowActivities()[rowid];

   bool unreliable;

   if( actchange == ActivityChange::kMin )
      unreliable = ( activity.ninfmin <= 1 && activity.min != 0 &&
                     origactivity.min != 0 &&
                     num.isZero( activity.min / origactivity.min ) );
   else
      unreliable = ( activity.ninfmax <= 1 && activity.max != 0 &&
                     origactivity.max != 0 &&
                     num.isZero( activity.max / origactivity.max ) );

   if( unreliable )
   {
      auto rowvec = problem.getConstraintMatrix().getRowCoefficients( rowid );

      activity = compute_row_activity( rowvec.getValues(), rowvec.getIndices(),
                                       rowvec.getLength(), probing_lower_bounds,
                                       probing_upper_bounds,
                                       probing_domain_flags, round );
   }

   // check for infeasibility
   if( actchange == ActivityChange::kMin && activity.ninfmin == 0 &&
       !rflags[rowid].test( RowFlag::kRhsInf ) &&
       num.isFeasLT( rhs[rowid], activity.min ) &&
       num.isSafeLT( rhs[rowid], activity.min ) )
   {
      Message::debug( this,
                      "[{}:{}] probing on col {} with val {} is infeasible min "
                      "activity is {:.15}, right hand side is {:.15}, and "
                      "original max activity was {:.15}\n",
                      __FILE__, __LINE__, probingCol, probingValue,
                      double( activity.min ), double( rhs[rowid] ),
                      double( problem.getRowActivities()[rowid].min ) );
      infeasible = true;
   }

   if( actchange == ActivityChange::kMax && activity.ninfmax == 0 &&
       !rflags[rowid].test( RowFlag::kLhsInf ) &&
       num.isFeasGT( lhs[rowid], activity.max ) &&
       num.isSafeGT( lhs[rowid], activity.max ) )
   {
      Message::debug( this,
                      "[{}:{}] probing on col {} with val {} is infeasible max "
                      "activity is {:.15}, left hand side is {:.15}, and "
                      "original max activity was {:.15}\n",
                      __FILE__, __LINE__, probingCol, probingValue,
                      double( activity.max ), double( lhs[rowid] ),
                      double( problem.getRowActivities()[rowid].max ) );
      infeasible = true;
   }
}

template <typename REAL>
void
CliqueProbingView<REAL>::changeLb( int col, REAL newlb )
{
   const auto& consMatrix = problem.getConstraintMatrix();
   auto colvec = consMatrix.getColumnCoefficients( col );
   const auto& orig_lbs = problem.getLowerBounds();

   // bound must be tighter than current domains
   bool lbinf = probing_domain_flags[col].test( ColFlag::kLbUseless );
   assert( lbinf || probing_lower_bounds[col] != newlb );

   Message::debug( this, "changing probing lower bound of col {} to {}\n", col,
                   double( newlb ) );

   if( lbinf )
   {
      // bound was not altered yet, store the negative (index + 1) to
      // indicate that the infinity flag was altered
      probing_domain_flags[col].unset( ColFlag::kLbUseless );
      changed_lbs.push_back( -col - 1 );
   }
   else if( probing_lower_bounds[col] == orig_lbs[col] &&
            !problem.getColFlags()[col].test( ColFlag::kLbUseless ) )
      // if bound was not altered yet remember it in the index vector
      changed_lbs.push_back( col );

   // change the bound in the domain vector
   REAL oldlb = probing_lower_bounds[col];
   probing_lower_bounds[col] = newlb;

   // update the probing activities by using the column view
   update_activities_after_boundchange(
       colvec.getValues(), colvec.getIndices(), colvec.getLength(),
       BoundChange::kLower, oldlb, newlb, lbinf, probing_activities,
       [this]( ActivityChange actChange, int rowid,
               RowActivity<REAL>& activity ) {
          activityChanged( actChange, rowid, activity );
       },
       true );
}

template <typename REAL>
void
CliqueProbingView<REAL>::changeUb( int col, REAL newub )
{
   const auto& consMatrix = problem.getConstraintMatrix();
   auto colvec = consMatrix.getColumnCoefficients( col );
   const auto& orig_ubs = problem.getUpperBounds();

   // bound must be tighter than current domains
   bool ubinf = probing_domain_flags[col].test( ColFlag::kUbUseless );
   assert( ubinf || probing_upper_bounds[col] != newub );

   Message::debug( this, "changing probing upper bound of col {} to {}\n", col,
                   double( newub ) );

   if( ubinf )
   {
      // bound was not altered yet, store the negative (index + 1) to
      // indicate that the infinity flag was altered
      probing_domain_flags[col].unset( ColFlag::kUbUseless );
      changed_ubs.push_back( -col - 1 );
   }
   else if( probing_upper_bounds[col] == orig_ubs[col] &&
            !problem.getColFlags()[col].test( ColFlag::kUbUseless ) )
      // if bound was not altered yet remember it in the index vector
      changed_ubs.push_back( col );

   // change the bound in the domain vector
   REAL oldub = probing_upper_bounds[col];
   probing_upper_bounds[col] = newub;

   // update the probing activities by using the column view
   update_activities_after_boundchange(
       colvec.getValues(), colvec.getIndices(), colvec.getLength(),
       BoundChange::kUpper, oldub, newub, ubinf, probing_activities,
       [this]( ActivityChange actChange, int rowid,
               RowActivity<REAL>& activity ) {
          activityChanged( actChange, rowid, activity );
       },
       true );
}


template <typename REAL>
bool
CliqueProbingView<REAL>::analyzeImplications()
{
   for( int ind = 0; ind < static_cast<int>(fix_to_zero.end() - fix_to_zero.begin()); ++ind )
   {
      for(int i = 0; i < cliquelen; ++i )
      {
         if( fix_to_zero[ind] != cliqueind[i] )
            continue;
         /*reset();
         setProbingColumn(i);
         propagateDomains();
         assert( isInfeasible() );*/
      }
      boundChanges.emplace_back(
         CliqueProbingBoundChg<REAL>( true, fix_to_zero[ind], 0.0, -1 ) );
   }
   if( fix_to_zero.end() - fix_to_zero.begin() == cliquelen && cliqueEquation )
   {
      return true;
   }
   if( cliqueEquation )
   {
      /*reset();
      setProbingColumn(-1);
      propagateDomains();
      if( !isInfeasible() && !equationBefore )
      {
         //std::cout<<"\n\nError: Row " << probingClique << " is marked as equation but shoudln't.";
         //std::cout.flush();
      }
      assert( isInfeasible() || equationBefore );*/
   }
   if( fewreductions )
      return false;
   else
   {
      for( typename std::list<std::pair<int,REAL>>::iterator col = changed_clique_lbs_inds_vals.begin();
      col != changed_clique_lbs_inds_vals.end(); std::advance(col,1) )
      {
         /*for(int i = -1 + equationBefore; i < cliquelen; ++i )
         {
            reset();
            setProbingColumn(i);
            propagateDomains();
            if( !isInfeasible() && !num.isGE(probing_upper_bounds[(*col).first], (*col).second) )
            {
               std::cout<<"\nError, Probing ";
               if (i==-1)
                  std::cout<<" all on zero ";
               else
               {   
                  std::cout<<cliqueind[i];
                  std::cout<<" on one and rest on zero ";
               }
               std::cout<<"has weaker lower bounds ";
               std::cout<<probing_lower_bounds[(*col).first];
               std::cout<<" ";
               std::cout<<(*col).second;
               std::cout<<" for variable ";
               std::cout<<(*col).first;
               std::cout.flush();
            }
            assert( isInfeasible() || num.isGE(probing_lower_bounds[(*col).first], (*col).second) );
         }*/
         boundChanges.emplace_back(
            CliqueProbingBoundChg<REAL>( false, (*col).first, (*col).second, cliqueind[0] ) );
      }
      for( typename std::list<std::pair<int,REAL>>::iterator col = changed_clique_ubs_inds_vals.begin();
      col != changed_clique_ubs_inds_vals.end(); std::advance(col,1) )
      {
         /*for(int i = -1 + equationBefore; i < cliquelen; ++i )
         {
            reset();
            setProbingColumn(i);
            propagateDomains();
            if( !isInfeasible() && !num.isLE(probing_upper_bounds[(*col).first], (*col).second) )
            {
               std::cout<<"\nError, Probing ";
               if (i==-1)
                  std::cout<<" all on zero ";
               else
               {   
                  std::cout<<cliqueind[i];
                  std::cout<<" on one and rest on zero ";
               }
               std::cout<<"has weaker upper bounds ";
               std::cout<<probing_upper_bounds[(*col).first];
               std::cout<<" ";
               std::cout<<(*col).second;
               std::cout<<" for variable ";
               std::cout<<(*col).first;
               std::cout.flush();
            }
            assert( isInfeasible() || num.isLE(probing_upper_bounds[(*col).first], (*col).second) );
         }*/
         boundChanges.emplace_back(
            CliqueProbingBoundChg<REAL>( true, (*col).first, (*col).second, cliqueind[0] ) );
      }
      //////std::cout<<"\nBinary inds: " << static_cast<int>(binary_inds.end() - binary_inds.begin());
      for( int ind = 0; ind < static_cast<int>(binary_inds.end() - binary_inds.begin()); ++ind )
      {
         if( lb_implications[ind].first == cliquelen - static_cast<int>(fix_to_zero.size())
            - static_cast<int>(cliqueEquation) && ub_implications[ind].first == 1 
            && ub_implications[ind].second != -1 && binary_inds[ind] != ub_implications[ind].second )
         {
            ////////std::cout<<"\nFOUND SUBSTITUTION";
            substitutions.emplace_back(
               CliqueProbingSubstitution<REAL>( binary_inds[ind], -1.0, ub_implications[ind].second, 1.0 ) );

            /*int i = -1 + equationBefore;
            while( i != cliquelen )
            {
               reset();
               setProbingColumn(i);
               propagateDomains();
               if( i != -1 && cliqueind[i] == ub_implications[ind].second )
                  assert( probing_upper_bounds[binary_inds[ind]] == 0.0 && !isInfeasible() );
               else if( probing_lower_bounds[binary_inds[ind]] != 1.0 && !isInfeasible() )
               {
                  std::cout<<"\nImplicationtest: " << binary_inds[ind] << "\nUbimpsfirst: " << ub_implications[ind].first << " cliquelen: " << cliquelen 
                  << " static_cast<int>(fix_to_zero.size()) " << static_cast<int>(fix_to_zero.size()) <<
                  " static_cast<int>(cliqueEquation) " << static_cast<int>(cliqueEquation) << " lb_implications[ind].first "
                  << lb_implications[ind].first << " lb_implications[ind].second " << lb_implications[ind].second
                  << " ub_implications[ind].second " << ub_implications[ind].second;
                  std::cout.flush();
                  assert( probing_lower_bounds[binary_inds[ind]] == 1.0 || isInfeasible() );
               }
               i +=1;
            }
            reset();*/
      
         }
         else if( ub_implications[ind].first == cliquelen - static_cast<int>(fix_to_zero.size())
            - static_cast<int>(cliqueEquation) && lb_implications[ind].first == 1 
            && lb_implications[ind].second != -1 && binary_inds[ind] != lb_implications[ind].second )
         {
            ////////std::cout<<"\nFOUND SUBSTITUTION";
            substitutions.emplace_back(
               CliqueProbingSubstitution<REAL>( binary_inds[ind], 1.0, lb_implications[ind].second, 0.0 ) );

            /*int i = -1 + equationBefore;
            while( i != cliquelen )
            {
               reset();
               setProbingColumn(i);
               propagateDomains();
               if( i != -1 && cliqueind[i] == lb_implications[ind].second )
                  assert( probing_lower_bounds[binary_inds[ind]] == 1.0 && !isInfeasible() );
               else if( !(probing_upper_bounds[binary_inds[ind]] == 0.0) && !isInfeasible() )
               {
                  std::cout<<"\nImplicationtest: " << binary_inds[ind] << "\nUbimpsfirst: " << ub_implications[ind].first << " cliquelen: " << cliquelen 
                  << " static_cast<int>(fix_to_zero.size()) " << static_cast<int>(fix_to_zero.size()) <<
                  " static_cast<int>(cliqueEquation) " << static_cast<int>(cliqueEquation) << " lb_implications[ind].first "
                  << lb_implications[ind].first << " lb_implications[ind].second " << lb_implications[ind].second
                  << " ub_implications[ind].second " << ub_implications[ind].second;
                  std::cout.flush();
                  assert( probing_upper_bounds[binary_inds[ind]] == 0.0 || isInfeasible() );
               }
               i +=1;
            }
            reset();*/

         }
            //////std::cout<<"\nImplicationtest: " << binary_inds[ind] << "\nUbimpsfirst: " << ub_implications[ind].first << " cliquelen: " << cliquelen 
            //<< " static_cast<int>(fix_to_zero.size()) " << static_cast<int>(fix_to_zero.size()) <<
            //" static_cast<int>(cliqueEquation) " << static_cast<int>(cliqueEquation) << " lb_implications[ind].first "
            //<< lb_implications[ind].first << " lb_implications[ind].second " << lb_implications[ind].second;
      }
      return false;
   }
}

template <typename REAL>
void
CliqueProbingView<REAL>::propagateDomains()
{
   const auto& consMatrix = problem.getConstraintMatrix();
   const auto& lhs = consMatrix.getLeftHandSides();
   const auto& rhs = consMatrix.getRightHandSides();
   const auto& rflags = consMatrix.getRowFlags();

   using std::swap;

   swap( prop_activities, next_prop_activities );
   next_prop_activities.clear();

   while( !prop_activities.empty() )
   {
      int curr_round = round--;

      Message::debug( this,
                      "starting probing propagation round {} on {} rows\n",
                      -curr_round - 2, prop_activities.size() );

      for( int candrow : prop_activities )
      {
         bool propagate = false;

         if( !rflags[candrow].test( RowFlag::kRhsInf ) &&
             probing_activities[candrow].ninfmin <= 1 )
            propagate = true;

         if( !rflags[candrow].test( RowFlag::kLhsInf ) &&
             probing_activities[candrow].ninfmax <= 1 )
            propagate = true;

         if( !propagate )
            continue;

         auto rowvec = consMatrix.getRowCoefficients( candrow );

         propagate_row( num, candrow,
             rowvec.getValues(), rowvec.getIndices(), rowvec.getLength(),
             probing_activities[candrow], lhs[candrow], rhs[candrow],
             rflags[candrow], probing_lower_bounds, probing_upper_bounds,
             probing_domain_flags,
             [this]( BoundChange bndChg, int colid, REAL newbound , int row ) {
                if( num.isHugeVal( newbound ) )
                   return;

                bool isint = probing_domain_flags[colid].test(
                    ColFlag::kIntegral, ColFlag::kImplInt );

                if( bndChg == BoundChange::kLower )
                {
                   if( isint )
                      newbound = num.feasCeil( newbound );

                   if( !probing_domain_flags[colid].test( ColFlag::kUbInf ) &&
                       newbound > probing_upper_bounds[colid] )
                   {
                      if( num.isFeasGT( newbound,
                                        probing_upper_bounds[colid] ) )
                      {
                         Message::debug( this,
                                         "[{}:{}] probing on col {} with "
                                         "val {} is infeasible\n",
                                         __FILE__, __LINE__, probingCol,
                                         probingValue );
                         infeasible = true;
                         return;
                      }

                      newbound = probing_upper_bounds[colid];
                   }

                   REAL delta = newbound - probing_lower_bounds[colid];
                   bool finiteDomain =
                       !probing_domain_flags[colid].test( ColFlag::kUbInf );

                   REAL mindomred = isint ? minintdomred : mincontdomred;

                   // Boundchanges on unbounded variables [lb, inf) or (-inf, ub]
                   // are ignored since their fixing will probably not lead to
                   // additional fixings nor infeasibility
                   // but will trigger a huge sequence of bound changes
                   if( probing_domain_flags[colid].test(
                           ColFlag::kLbUseless ) ||
                       ( finiteDomain && delta > 0 &&
                         ( delta / ( probing_upper_bounds[colid] -
                                     probing_lower_bounds[colid] ) >=
                           mindomred ) ) )
                      changeLb( colid, newbound );
                }
                else
                {
                   assert( bndChg == BoundChange::kUpper );
                   if( isint )
                      newbound = num.feasFloor( newbound );

                   if( !probing_domain_flags[colid].test( ColFlag::kLbInf ) &&
                       newbound < probing_lower_bounds[colid] )
                   {
                      if( num.isFeasLT( newbound,
                                        probing_lower_bounds[colid] ) )
                      {
                         Message::debug( this,
                                         "[{}:{}] probing on col {} with "
                                         "val {} is infeasible\n",
                                         __FILE__, __LINE__, probingCol,
                                         probingValue );
                         infeasible = true;
                         return;
                      }

                      newbound = probing_lower_bounds[colid];
                   }

                   REAL delta = probing_upper_bounds[colid] - newbound;
                   bool finiteDomain =
                       !probing_domain_flags[colid].test( ColFlag::kLbInf );

                   REAL mindomred = isint ? minintdomred : mincontdomred;

                   // Boundchanges on unbounded variables [lb, inf) or (-inf, ub]
                   // are ignored since their fixing will probably not lead to
                   // additional fixings nor infeasibility
                   // but will trigger a huge sequence of bound changes
                   if( probing_domain_flags[colid].test(
                           ColFlag::kUbUseless ) ||
                       ( finiteDomain && delta > 0 &&
                         ( delta / ( probing_upper_bounds[colid] -
                                     probing_lower_bounds[colid] ) >=
                           mindomred ) ) )
                      changeUb( colid, newbound );
                }
             } );
         if( infeasible )
            return;
      }

      swap( prop_activities, next_prop_activities );
      next_prop_activities.clear();
   }
}

} // namespace papilo

#endif
