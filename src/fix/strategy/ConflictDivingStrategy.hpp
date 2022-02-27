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

#ifndef FIX_CONFLICT_DIVING_STRATEGY_HPP
#define FIX_CONFLICT_DIVING_STRATEGY_HPP

template <typename REAL>
class ConflictDivingStrategy : public RoundingStrategy<REAL>
{
   RandomGenerator random;
   const Num<REAL> num;
   Vec<int> n_conflict_down_locks;
   Vec<int> n_conflict_up_locks;
   Problem<REAL>& problem;
   int n_cols;
   int n_conflicts;
   Vec<int> n_var_down_locks;
   Vec<int> n_var_up_locks;

 public:
   ConflictDivingStrategy( RandomGenerator random_, Num<REAL> num_,
                           Problem<REAL>& problem_ )
       : random( random_ ), num( num_ ), n_conflict_down_locks( {} ),
         n_conflict_up_locks( {} ), problem( problem_ ),
         n_cols( problem_.getNCols() ), n_conflicts( 0 )
   {
      n_conflict_down_locks.resize( n_cols );
      n_conflict_up_locks.resize( n_cols );

      for( int col = 0; col < n_cols; ++col )
      {
         int n_up_locks = 0;
         int n_down_locks = 0;

         auto colvec = problem.getConstraintMatrix().
                                 getColumnCoefficients( col );

         const REAL* vals = colvec.getValues();
         const int* inds = colvec.getIndices();
         int len = colvec.getLength();
         auto rflags = problem.getRowFlags();

         for( int i = 0; i < len; i++ )
         {
            assert( !rflags[inds[i]].test( RowFlag::kConflictConstraint ) );
            count_locks( vals[i], rflags[inds[i]], n_down_locks,
                         n_up_locks );
         }

         assert( ( n_down_locks > 0 ) || ( n_up_locks > 0 ) );
         n_var_down_locks.push_back( n_down_locks );
         n_var_up_locks.push_back( n_up_locks );
      }
      assert( n_var_down_locks.size() == n_cols );
      assert( n_var_up_locks.size() == n_cols );
   }

   void
   recompute_locks() override
   {
      update_n_conflicts();

      for( int col = 0; col < n_cols; ++col )
      {
         int n_up_locks = 0;
         int n_down_locks = 0;

         auto colvec = problem.getConstraintMatrix().
                                 getColumnCoefficients( col );

         const REAL* vals = colvec.getValues();
         const int* inds = colvec.getIndices();
         int len = colvec.getLength();
         auto rflags = problem.getRowFlags();

         for( int i = 0; i < len; i++ )
            if( rflags[inds[i]].test( RowFlag::kConflictConstraint ) )
               count_locks( vals[i], rflags[inds[i]], n_down_locks,
                            n_up_locks );

         n_conflict_down_locks[col] = n_down_locks;
         n_conflict_up_locks[col] = n_up_locks;
      }

      assert( n_conflict_down_locks.size() == n_cols );
      assert( n_conflict_up_locks.size() == n_cols );
   }

   Fixing<REAL>
   select_rounding_variable( const Vec<REAL>& cont_solution,
                             const ProbingView<REAL>& view ) override
   {
      assert( n_conflict_down_locks.size() == n_cols );
      assert( n_conflict_up_locks.size() == n_cols );

      REAL value = -1;
      int variable = -1;
      REAL score = -1;
      REAL max_score = -1;
      // 0 = round down, 1 = round up
      bool round_up;
      auto obj = view.get_obj();

      for( int i = 0; i < cont_solution.size(); i++ )
      {
         if( num.isIntegral( cont_solution[i] ) ||
             num.isEq( view.getProbingUpperBounds()[i],
                       view.getProbingLowerBounds()[i] ) ||
             !view.is_within_bounds( i, cont_solution[i] ) ||
             !view.is_integer_variable( i ) )
            continue;

         REAL frac = cont_solution[i] - num.epsFloor( cont_solution[i] );
         int n_down_locks = n_conflict_down_locks[i];
         int n_up_locks = n_conflict_up_locks[i];
         // TODO: verify the logic!
         /*
         bool may_round_up = ( ( n_down_locks == 0 ) ||
                               ( n_down_locks > n_up_locks ) ||
                               ( num.isLT( frac, REAL{ 0.5 } ) ) );
         bool may_round_down = ( ( n_up_locks == 0 ) ||
                                 ( n_up_locks > n_down_locks ) ||
                                 ( num.isGT( frac, REAL{ 0.5 } ) ) );
         */
         bool may_round_down = ( n_down_locks == 0 );
         bool may_round_up = ( n_up_locks == 0 );

         score = get_score( i, may_round_down, may_round_up, frac, view,
                            round_up );

         if( num.isGT( score, max_score ) )
         {
            max_score = score;
            variable = i;

            if( round_up )
               value = num.epsCeil( cont_solution[i] );
            else
               value = num.epsFloor( cont_solution[i] );
         }
      }
      return { variable, value };
   }

private:
   void
   update_n_conflicts()
   {
      n_conflicts = 0;
      auto rflags = problem.getRowFlags();

      for( int i = 0; i < problem.getNRows(); i++ )
         if( rflags[i].test( RowFlag::kConflictConstraint ) )
            n_conflicts++;
   }

   REAL
   get_score( const int col, const bool may_round_down, const bool may_round_up,
              const REAL frac, const ProbingView<REAL>& view, bool& round_up )
   {
      REAL score = 0;
      int num_conflict_down_locks = n_conflict_down_locks[col];
      int num_conflict_up_locks = n_conflict_up_locks[col];
      int num_var_down_locks = n_var_down_locks[col];
      int num_var_up_locks = n_var_up_locks[col];
      std::uniform_int_distribution<uint32_t> dist_rounding( 0, 1 );
      std::uniform_real_distribution<> double_dist( 1e-6, 1e-5 );
      REAL threshold = 0.2;
      REAL epsilon = 0.25;
      round_up = 0;

      if( may_round_down || may_round_up )
      {
         if( may_round_down && may_round_up )
         {
            if( num_var_down_locks == num_var_up_locks )
            {
               if( num.isEq( frac, REAL{ 0.5 } ) )
                  round_up = random.get_random_int( dist_rounding );
               else
                  round_up = num.isGT( frac, REAL{ 0.5 } );
            }
            else
               round_up = ( num_var_up_locks > num_var_down_locks );
         }
         else
            round_up = may_round_up;
      }
      else
      {
         assert(!may_round_down);

         if( !num.isEq( num_conflict_down_locks, num_conflict_up_locks ) )
         {
            // TODO: isLT at least gives solutions for mcsched
            round_up = num.isGT( num_conflict_up_locks, num_conflict_down_locks );
         }
         else if( num.isEq( frac, REAL{ 0.5 } ) )
            round_up = random.get_random_int( dist_rounding );
         else
            round_up = num.isGT( frac, REAL{ 0.5 } );
      }

      if( round_up )
         score = num_conflict_up_locks +
                 ( epsilon * num_var_up_locks /
                   ( num_var_up_locks + num_var_down_locks ) ) +
                 random.get_random_double( double_dist );
      else
         score = num_conflict_down_locks +
                 ( epsilon * num_var_down_locks /
                   ( num_var_up_locks + num_var_down_locks ) ) +
                 random.get_random_double( double_dist );

      /* penalize too few conflict locks */
      if( ( ( num_conflict_down_locks + num_conflict_up_locks ) > 0 ) &&
          ( ( num_conflict_down_locks + num_conflict_up_locks ) < ( threshold * n_conflicts ) ) )
         score *= 0.1;

      /* penalize zero conflict locks */
      if( ( num_conflict_down_locks + num_conflict_up_locks ) == 0 )
         score *= 0.01;

      /* penalize too integral values */
      if( num.isLT( frac, REAL{ 0.01 } ) || num.isGT( frac, REAL{ 0.99 } ) )
         score *= 0.01;

      assert( view.is_integer_variable( col ) );
      /* penalize non-binary variables */
      if( !num.isEq( view.getProbingUpperBounds()[col], REAL{ 1.0 } ) ||
          !num.isEq( view.getProbingLowerBounds()[col], REAL{ 0.0 } ) )
         score = -1.0 / score;

      return score;
   }
};

#endif
