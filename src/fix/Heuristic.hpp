/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*               This file is part of the program and library                */
/*    PaPILO --- Parallel Presolve for Integer and Linear Optimization       */
/*                                                                           */
/* Copyright (C) 2020-2022 Konrad-Zuse-Zentrum                               */
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

#ifndef FIX_FIX_AND_PROPAGATE_SERVICE_HPP
#define FIX_FIX_AND_PROPAGATE_SERVICE_HPP

#include "fix/FixAndPropagate.hpp"
#include "fix/strategy/FarkasRoundingStrategy.hpp"
#include "fix/strategy/FractionalRoundingStrategy.hpp"
#include "fix/strategy/RandomRoundingStrategy.hpp"
#include "papilo/core/Objective.hpp"
#include "papilo/core/Presolve.hpp"

#include <cassert>
#include <fstream>
#include <string>

using namespace papilo;

template <typename REAL>
class Heuristic
{
   Message msg;
   Num<REAL> num;
   Timer timer;
   Vec<RoundingStrategy<REAL>*> strategies{};
   Vec<Vec<REAL>> int_solutions{};
   Vec<ProbingView<REAL>> views{};
   Vec<REAL> obj_value;
   Vec<bool> infeasible_arr;

 public:
   Heuristic( Message msg_, Num<REAL> num_, Timer& timer_ )
       : msg( msg_ ), num( num_ ), timer( timer_ ), strategies( {} ),
         int_solutions( {} ), views( {} ), obj_value( {} ), infeasible_arr( {} )
   {
   }

   void
   setup( const Problem<REAL>& problem )
   {
#ifdef PAPILO_TBB
      auto s1 = new FarkasRoundingStrategy<REAL>{ 0, num, false };
      auto s2 = new FarkasRoundingStrategy<REAL>{ 0, num, true };
      auto s3 = new FractionalRoundingStrategy<REAL>{ num };
      auto s4 = new RandomRoundingStrategy<REAL>{ 0, num };
      strategies.push_back( s1 );
      strategies.push_back( s2 );
      strategies.push_back( s3 );
      strategies.push_back( s4 );

      Vec<REAL> int_solution{};
      int_solution.resize( problem.getNCols() );

      int_solutions.push_back( { int_solution } );
      int_solutions.push_back( { int_solution } );
      int_solutions.push_back( { int_solution } );
      int_solutions.push_back( { int_solution } );

      views.push_back( { problem, num } );
      views.push_back( { problem, num } );
      views.push_back( { problem, num } );
      views.push_back( { problem, num } );

      infeasible_arr.push_back( true );
      infeasible_arr.push_back( true );
      infeasible_arr.push_back( true );
      infeasible_arr.push_back( true );

      obj_value.push_back( 0 );
      obj_value.push_back( 0 );
      obj_value.push_back( 0 );
      obj_value.push_back( 0 );
#else
      Vec<REAL> int_solution{};
      int_solution.resize( problem.getNCols() );
      auto s1 = new FarkasRoundingStrategy<REAL>{ 0, num, false };
      strategies.push_back( s1 );
      int_solutions.push_back( { int_solution } );
      views.push_back( { problem, num } );
      infeasible_arr.push_back( true );
      obj_value.push_back( 0 );
#endif
   }

   void
   perform_fix_and_propagate( const Vec<REAL>& primal_heur_sol,
                              REAL& best_obj_val,
                              Vec<REAL>& current_best_solution )
   {
      for( auto view : views )
         view.reset();
#ifdef PAPILO_TBB

      tbb::parallel_for(
          tbb::blocked_range<int>( 0, 4 ),
          [&]( const tbb::blocked_range<int>& r )
          {
             for( int i = r.begin(); i != r.end(); ++i )
             {
                FixAndPropagate<REAL> fixAndPropagate{ msg, num, true };
                infeasible_arr[i] = fixAndPropagate.fix_and_propagate(
                    primal_heur_sol, int_solutions[i], *( strategies[i] ),
                    views[i] );
                if( infeasible_arr[i] )
                {
                   obj_value[i] = 0;
                   break;
                }
                StableSum<REAL> sum{};
                for( int j = 0; j < primal_heur_sol.size(); j++ )
                   sum.add( int_solutions[i][j] * views[i].get_obj()[j] );
                obj_value[i] = sum.get();
                msg.info( "Diving {} found obj value {}!\n", i, obj_value[i] );
             }
          } );
#else

      FixAndPropagate<REAL> fixAndPropagate{ msg, num, true };
      infeasible_arr[0] = fixAndPropagate.fix_and_propagate(
          primal_heur_sol, int_solutions[0], *( strategies[0] ), views[0] );
      if( infeasible_arr[0] )
      {
         obj_value[0] = 0;
         return;
      }
      StableSum<REAL> sum{};
      for( int j = 0; j < primal_heur_sol.size(); j++ )
         sum.add( int_solutions[0][j] * views[0].get_obj()[j] );
      obj_value[0] = sum.get();
      msg.info( "Diving {} found obj value {}!\n", 0, obj_value[0] );
#endif
      evaluate( best_obj_val, current_best_solution );
   }

   void
   evaluate( REAL& best_obj_val, Vec<REAL>& current_best_solution )
   {
      bool feasible = std::any_of( infeasible_arr.begin(), infeasible_arr.end(),
                                   []( bool b ) { return !b; } );

      // TODO: copy the best solution;
      if( !feasible )
      {
         msg.info( "Fix and Propagate did not find a feasible solution!\n" );
         return;
      }

      int best_index = -1;
      for( int i = 0; i < obj_value.size(); i++ )
      {
         if( !infeasible_arr[i] &&
             ( num.isLT( obj_value[i], best_obj_val ) ||
               ( current_best_solution.empty() && best_index == -1 ) ) )
         {
            best_index = i;
            best_obj_val = obj_value[i];
         }
      }
      if( best_index == -1 )
      {
         msg.info(
             "Fix and Propagate did not improve the current solution!\n" );
         return;
      }

      if( current_best_solution.empty() )
         msg.info( "Fix and Propagate found an initial solution: {}!\n",
                   best_obj_val );
      else
         msg.info( "Fix and Propagate found a new solution: {}!\n",
                   best_obj_val );

      current_best_solution = int_solutions[best_index];
      assert( best_obj_val == obj_value[best_index] );
   }
};

#endif
