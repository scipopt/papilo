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

#include "papilo/core/Problem.hpp"
#include "papilo/core/SingleRow.hpp"
#include "papilo/io/Message.hpp"
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

   std::pair<bool,bool>
   probeClique( const int clique, const int*& indices, const int len, const Vec<int>& binary_inds)
   {  
      
      //std::cout<<( "\n\nProbing Clique: ");
      probingClique = clique;
      //std::cout<<( clique);
      //std::cout<<( "\n");
      for( int ind = 0; ind < len; ++ind )
      {
         cliqueind.emplace_back( indices[ind] );
         //std::cout<<(indices[ind]);
         //std::cout<<(" ");
      }
      cliquelen = len;
      assert(len == static_cast<int>(cliqueind.size()));
      bool initbounds = false;
      lb_implications.reserve( static_cast<int>(binary_inds.size()) );
      ub_implications.reserve( static_cast<int>(binary_inds.size()) );
      for( int ind = 0 ; ind != static_cast<int>(binary_inds.size()); ++ind )
      {
        lb_implications.emplace_back( std::pair<int,int> {0,-1} );
        ub_implications.emplace_back( std::pair<int,int> {0,-1} );
      } 

      setProbingColumn(-1);

      cliqueEquation = false;
      propagateDomains();
      if( isInfeasible() )
      {
         //std::cout<<"\nAll 0 infeasible, detected clique equation.\n";
         cliqueEquation = true;
      }
      else
      {
         //std::cout<<("Initializing\n");
         for( int var = 0; var != static_cast<int>(probing_lower_bounds.size()); ++var )
         {
            if( num.isGT(probing_lower_bounds[var], problem.getLowerBounds()[var]) )
            {
               changed_clique_lbs_inds_vals.emplace_back(std::pair<int,REAL>{var, probing_lower_bounds[var]});
               //std::cout<<"The lower bound of ";
               //std::cout<<var;
               //std::cout<<" can be changed after probing everything on zero: ";
               //std::cout<< probing_lower_bounds[var];
               //std::cout<< " ";
               //std::cout<<problem.getLowerBounds()[var];
               //std::cout<<"\n";
            }
            if( num.isLT(probing_upper_bounds[var], problem.getUpperBounds()[var]) )
            {
               changed_clique_ubs_inds_vals.emplace_back(std::pair<int,REAL> {var, probing_upper_bounds[var]});
               //std::cout<<"The upper bound of ";
               //std::cout<<var;
               //std::cout<<" can be changed after probing everything on zero: ";
               //std::cout<< probing_upper_bounds[var];
               //std::cout<< " ";
               //std::cout<<problem.getUpperBounds()[var];
               //std::cout<<"\n";
            }
         }
         initbounds = true;
         //std::cout<<("Initialized\n");
      }
      reset();

      for( int i = 0; i < cliquelen; ++i )
      {
         reset();
         assert( probing_upper_bounds[cliqueind[i]] == 1.0 && probing_lower_bounds[cliqueind[i]] == 0.0 );
         //std::cout<<( "\nSetting probing collumn: ");
         //std::cout<<( cliqueind[i]);
        setProbingColumn(i);
        //std::cout<<( "\nSet probing collumn.");
        //std::cout<<( "\nPropagating.");
        propagateDomains();
      //std::cout<<( "\nPropagated\n");
      //std::cout<<( "Checking for infeasibility\n");
        if( isInfeasible() )
        {
            //std::cout<<"\nInfeasible one assignment, fixing to zero.\n";
            fix_to_zero.emplace_back(probingCol);
            //std::cout<<( "1 is infeasible, added to fix_to_zero and skipped everything else.\n");
            reset();
            continue;
        }
        //std::cout<<( "Changing imps\n");
        for( int ind = 0; ind != static_cast<int>(binary_inds.size()); ++ind )
        {
            assert( ind < static_cast<int>(binary_inds.size()) );
            assert( binary_inds[ind] < static_cast<int>(probing_lower_bounds.size()) );
            if( num.isEq(1.0, probing_lower_bounds[binary_inds[ind]]) )
            {
               assert( ind < static_cast<int>(lb_implications.size()) );
               lb_implications[ind].first += 1;
               lb_implications[ind].second = probingCol;
            }
            assert( ind < static_cast<int>(binary_inds.size()) );
            assert( binary_inds[ind] < static_cast<int>(probing_upper_bounds.size()) );
            if( num.isEq(0.0, probing_upper_bounds[binary_inds[ind]]) )
            {
               assert( ind < static_cast<int>(ub_implications.size()) );
               ub_implications[ind].first += 1;                
               ub_implications[ind].second = probingCol;
            }
        }
        //std::cout<<("Changed imps\n");
        if( initbounds == false )
        {
            //std::cout<<("Initializing\n");
            for( int var = 0; var != static_cast<int>(probing_lower_bounds.size()); ++var )
            {
               if( num.isGT(probing_lower_bounds[var], problem.getLowerBounds()[var]) )
               {
                  changed_clique_lbs_inds_vals.emplace_back(std::pair<int,REAL>{var, probing_lower_bounds[var]});
                  //std::cout<<"The lower bound of ";
                  //std::cout<<var;
                  //std::cout<<" can be changed after probing the current variable on one: ";
                  //std::cout<< probing_lower_bounds[var];
                  //std::cout<< " ";
                  //std::cout<<problem.getLowerBounds()[var];
                  //std::cout<<"\n";
               }
               if( num.isLT(probing_upper_bounds[var], problem.getUpperBounds()[var]) )
               {
                  changed_clique_ubs_inds_vals.emplace_back(std::pair<int,REAL> {var, probing_upper_bounds[var]});
                  //std::cout<<"The upper bound of ";
                  //std::cout<<var;
                  //std::cout<<" can be changed after probing the current variable on one: ";
                  //std::cout<< probing_upper_bounds[var];
                  //std::cout<< " ";
                  //std::cout<<problem.getUpperBounds()[var];
                  //std::cout<<"\n";
               }
            }
            initbounds = true;
            //std::cout<<("Initialized\n");
            reset();
            continue;
        }
        //std::cout<<( "Changing Bounds\n");
        typename std::list<std::pair<int,REAL>>::iterator ind = changed_clique_lbs_inds_vals.begin(); 
         while( ind != changed_clique_lbs_inds_vals.end() )
         {
            if( num.isLT( probing_lower_bounds[(*ind).first], (*ind).second ) )
            {
               //std::cout<<"The probing lower bound of ";
               //std::cout<<(*ind).first;
               //std::cout<<" is weakend after probing the current variable on one: ";
               //std::cout<< probing_lower_bounds[(*ind).first];
               //std::cout<< " ";
               //std::cout<<(*ind).second;
               //std::cout<<"\n";
               (*ind).second = probing_lower_bounds[(*ind).first];
            }
            if( num.isLE(probing_lower_bounds[(*ind).first], problem.getLowerBounds()[(*ind).first] ) )
            {
               //std::cout<<"The probing lower bound of ";
               //std::cout<<(*ind).first;
               //std::cout<<" is deleted after probing the current variable on one: ";
               //std::cout<< probing_lower_bounds[(*ind).first];
               //std::cout<< " ";
               //std::cout<<problem.getLowerBounds()[(*ind).first];
               //std::cout<<"\n";
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
               //std::cout<<"The probing upper bound of ";
               //std::cout<<(*ind).first;
               //std::cout<<" is weakend after probing the current variable on one: ";
               //std::cout<< probing_upper_bounds[(*ind).first];
               //std::cout<< " ";
               //std::cout<<(*ind).second;
               //std::cout<<"\n";
               (*ind).second = probing_upper_bounds[(*ind).first];
            }
            if( num.isGE(probing_upper_bounds[(*ind).first], problem.getUpperBounds()[(*ind).first] ) )
            {  
               //std::cout<<"The probing upper bound of ";
               //std::cout<<(*ind).first;
               //std::cout<<" is deleted after probing the current variable on one: ";
               //std::cout<< probing_upper_bounds[(*ind).first];
               //std::cout<< " ";
               //std::cout<<problem.getUpperBounds()[(*ind).first];
               //std::cout<<"\n";
               ind = changed_clique_ubs_inds_vals.erase(ind);
            }
            else
               std::advance(ind, 1);
         }
         //std::cout<<("Changed bounds.\n");
         //std::cout<<("Resetting probing col:\n");
        reset();
        //std::cout<<( "Reset probing col.\n");
      }
      //if( fix_to_zero.end() - fix_to_zero.begin() == cliquelen && cliqueEquation )
         //std::cout<<"\nInfeasibility due to zero fixings in probe Clique.\n";
      return { fix_to_zero.end() - fix_to_zero.begin() == cliquelen && cliqueEquation, cliqueEquation } ;
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
      reset();
    changed_clique_lbs_inds_vals.clear();
    changed_clique_lbs_inds_vals.clear();
    lb_implications.clear();
    ub_implications.clear();
    fix_to_zero.clear();
    cliqueind.clear();
    probingClique = -1;
    cliquelen = -1;
    cliqueEquation = false;
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

   // datastructures for storing result of probing on one value
   Vec<CliqueProbingBoundChg<REAL>> otherValueImplications;
   bool otherValueInfeasible;
   Vec<bool> infeasibleAssignments;

   // results of probing and statistics
   Vec<CliqueProbingBoundChg<REAL>> boundChanges;
   Vec<CliqueProbingSubstitution<REAL>> substitutions;

   int64_t amountofwork;
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
   probingCol = -1;
   probingClique = -1;
   probingValue = false;
   otherValueInfeasible = false;
   minintdomred = num.getFeasTol() * 1000;
   mincontdomred = 0.3;
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
   //const auto& orig_ubs = problem.getUpperBounds();
   //const auto& orig_lbs = problem.getLowerBounds();
   //const Vec<ColFlags>& orig_domain_flags = problem.getColFlags();
   //std::cout<<( "Analyzing Implications\n");
   if( fix_to_zero.end() - fix_to_zero.begin() == cliquelen && cliqueEquation )
   {
      //std::cout<<"\nInfeasibility due to zero fixings.\n";
      return true;
   }
   for( int ind = 0; ind < static_cast<int>(fix_to_zero.end() - fix_to_zero.begin()); ++ind )
   {
      
      boundChanges.emplace_back(
         CliqueProbingBoundChg<REAL>( true, fix_to_zero[ind], 0.0, -1 ) );
         //std::cout<<"\nFixed to zero: ";
         //std::cout<<fix_to_zero[ind];
   }

   for( typename std::list<std::pair<int,REAL>>::iterator col = changed_clique_lbs_inds_vals.begin(); 
   col != changed_clique_lbs_inds_vals.end(); std::advance(col,1) )
   {
      boundChanges.emplace_back(
         CliqueProbingBoundChg<REAL>( false, (*col).first, (*col).second, cliqueind[0] ) );
         //std::cout<<"\nLower Bound Change: ";
         //std::cout<<(*col).first;
         //std::cout<<" ";
         //std::cout<<(*col).second;
         /*reset();
         setProbingColumn(-1);
         propagateDomains();
         assert( num.isGE( probing_lower_bounds[(*col).first], (*col).second) || isInfeasible() );
         for( int ind = 0; ind < cliquelen; ++ind )
         {
            reset();
            changeUb( cliqueind[ind], 0.0 );
            propagateDomains();
            if( probing_lower_bounds[(*col).first] < (*col).second )
            {
               //std::cout<<"\nProbing ";
               //std::cout<<cliqueind[ind];
               //std::cout<<" to zero has a weaker lower bound than clique probing: ";
               //std::cout<<probing_lower_bounds[(*col).first];
               //std::cout<<" ";
               //std::cout<<(*col).second;
            }
            reset();
            setProbingColumn(ind);
            propagateDomains();
            assert( num.isGE( probing_lower_bounds[(*col).first], (*col).second) || isInfeasible() );
         }*/
   }
   for( typename std::list<std::pair<int,REAL>>::iterator col = changed_clique_ubs_inds_vals.begin(); 
   col != changed_clique_ubs_inds_vals.end(); std::advance(col,1) )
   {
      boundChanges.emplace_back(
         CliqueProbingBoundChg<REAL>( true, (*col).first, (*col).second, cliqueind[0] ) );
         //std::cout<<"\nUpper Bound Change: ";
         //std::cout<<(*col).first;
         //std::cout<<" ";
         //std::cout<<(*col).second;
         /*reset();
         setProbingColumn(-1);
         propagateDomains();
         assert( num.isLE( probing_upper_bounds[(*col).first], (*col).second) || isInfeasible() );
         for( int ind = 0; ind < cliquelen; ++ind )
         {
            reset();
            changeUb( cliqueind[ind], 0.0 );
            propagateDomains();
            if( probing_upper_bounds[(*col).first] > (*col).second )
            {
               //std::cout<<"\nProbing ";
               //std::cout<<cliqueind[ind];
               //std::cout<<" to zero has a weaker upper bound than clique probing: ";
               //std::cout<<probing_upper_bounds[(*col).first];
               //std::cout<<" ";
               //std::cout<<(*col).second;
            }
            reset();
            setProbingColumn(ind);
            propagateDomains();
            assert( (num.isLE( probing_upper_bounds[(*col).first] , (*col).second ) || isInfeasible()) );
         }*/
   }

   for( int ind = 0; ind < static_cast<int>(binary_inds.end() - binary_inds.begin()); ++ind )
   {
      if( lb_implications[ind].first == cliquelen - 1 - static_cast<int>(fix_to_zero.size()) 
          && ub_implications[ind].first == 1 )
      {
         substitutions.emplace_back(
            CliqueProbingSubstitution<REAL>( binary_inds[ind], 1.0, ub_implications[ind].second, 0.0 ) );
         //std::cout<<"\nSubstitution: ";
         //std::cout<<binary_inds[ind];
         //std::cout<<" ";
         //std::cout<<ub_implications[ind].second;
      }
      else if( ub_implications[ind].first == cliquelen - 1 - static_cast<int>(fix_to_zero.size()) 
          && lb_implications[ind].first == 1 )
      {
         substitutions.emplace_back(
            CliqueProbingSubstitution<REAL>( binary_inds[ind], -1.0, lb_implications[ind].second, 1.0 ) );
            //std::cout<<"\nSubstitution: ";
         //std::cout<<binary_inds[ind];
         //std::cout<<" ";
         //std::cout<<lb_implications[ind].second;
      }
   }
   return false;
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
