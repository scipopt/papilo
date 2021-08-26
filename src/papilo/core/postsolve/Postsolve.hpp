/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*               This file is part of the program and library                */
/*    PaPILO --- Parallel Presolve for Integer and Linear Optimization       */
/*                                                                           */
/* Copyright (C) 2020-2021 Konrad-Zuse-Zentrum                               */
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

#ifndef _PAPILO_CORE_POSTSOLVE_SERVICE_HPP_
#define _PAPILO_CORE_POSTSOLVE_SERVICE_HPP_

#include "PostsolveListener.hpp"
#include "papilo/core/Problem.hpp"
#include "papilo/core/ProblemUpdate.hpp"
#include "papilo/core/postsolve/PostsolveStatus.hpp"
#include "papilo/core/postsolve/PostsolveType.hpp"
#include "papilo/core/postsolve/ReductionType.hpp"
#include "papilo/core/Stored.hpp"
#include "papilo/core/SavedRow.hpp"
#include "papilo/io/MpsWriter.hpp"
#include "papilo/misc/MultiPrecision.hpp"
#include "papilo/misc/Num.hpp"
#include "papilo/misc/StableSum.hpp"
#include "papilo/misc/Vec.hpp"
#include "papilo/misc/compress_vector.hpp"
#include "papilo/misc/dualpostsolve/PrimalDualSolValidation.hpp"
#include "papilo/misc/fmt.hpp"
#include "papilo/misc/tbb.hpp"

#include <fstream>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/tmpdir.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/list.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>

namespace papilo
{

template <typename REAL>
class Postsolve
{
 private:
   Message message{};
   Num<REAL> num{};

   static constexpr int IS_INTEGRAL = static_cast<int>( ColFlag::kIntegral );
   static constexpr int IS_LBINF = static_cast<int>( ColFlag::kLbInf );
   static constexpr int IS_UBINF = static_cast<int>( ColFlag::kUbInf );

 public:
   Postsolve( const Message msg, const Num<REAL> n )
   {
      message = msg;
      num = n;
   };

   PostsolveStatus
   undo( const Solution<REAL>& reducedSolution,
         Solution<REAL>& originalSolution,
         PostsolveListener<REAL> postsolveListener ) const;

 private:
   REAL
   calculate_row_value_for_infinity_column( REAL lhs, REAL rhs, int rowLength,
                                            int column, const int* row_indices,
                                            const REAL* coefficients,
                                            Vec<REAL>& current_solution,
                                            bool is_negative, REAL& coeff_of_column_in_row ) const;

   Problem<REAL>
   calculate_current_problem( const PostsolveListener<REAL>& listener,
                              int current_index ) const;

   void
   copy_from_reduced_to_original(
       const Solution<REAL>& reducedSolution, Solution<REAL>& originalSolution,
       const PostsolveListener<REAL>& postsolveListener) const;

   void
   apply_fix_var_in_original_solution( Solution<REAL>& originalSolution,
                                 const Vec<int>& indices,
                                 const Vec<REAL>& values, int current_index ) const;

   int
   apply_fix_infinity_variable_in_original_solution(
       Solution<REAL>& originalSolution, Vec<int>& indices, Vec<REAL>& values,
       int first, const Problem<REAL>& problem,  Stored<REAL>& stored_bounds  ) const;

   void
   apply_var_bound_change_forced_by_column_in_original_solution(
       Solution<REAL>& originalSolution, const Vec<ReductionType>& types,
       const Vec<int>& start, const Vec<int>& indices, const Vec<REAL>& values,
       int i, int first, Stored<REAL>& stored_bounds ) const;

   void
   apply_parallel_col_to_original_solution( Solution<REAL>& originalSolution,
                                            const Vec<int>& indices,
                                            const Vec<REAL>& values, int first,
                                            int last,
                                            Stored<REAL>& stored ) const;

   void
   apply_row_bound_chang_to_original_solution(
       Solution<REAL>& originalSolution, const Vec<ReductionType>& types,
       const Vec<int>& start, const Vec<int>& indices, const Vec<REAL>& values,
       int i, int first) const;

   void
   apply_substituted_column_to_original_solution(
       Solution<REAL>& originalSolution,  const Vec<int>& indices,
       const Vec<REAL>& values, int first,
       int last ) const;

   VarBasisStatus
   set_basis( int flags, REAL lb, REAL ub, REAL solution,
              bool is_on_bounds ) const;

   bool
   are_the_next_n_types_redundant_rows( const Vec<ReductionType>& types, int i,
                                        int redundant_rows ) const;
   void
   remove_var_from_basis( Solution<REAL>& originalSolution,
                          const Vec<ReductionType>& types,
                          const Vec<int>& start, const Vec<int>& indices,
                          const Vec<REAL>& values, int i,
                          Stored<REAL>& stored_bounds ) const;

   VarBasisStatus
   get_var_basis_status( Stored<REAL>& stored_bounds,
                     int index, REAL val) const;
};

#ifdef PAPILO_USE_EXTERN_TEMPLATES
extern template class Postsolve<double>;
extern template class Postsolve<Quad>;
extern template class Postsolve<Rational>;
#endif

template <typename REAL>
PostsolveStatus
Postsolve<REAL>::undo( const Solution<REAL>& reducedSolution,
                       Solution<REAL>& originalSolution,
                       PostsolveListener<REAL> postsolveListener ) const
{

   PrimalDualSolValidation<REAL> validation{};

   const Vec<REAL>& reducedSol = reducedSolution.primal;

   copy_from_reduced_to_original( reducedSolution, originalSolution,
                                  postsolveListener );


   auto types = postsolveListener.types;
   auto start = postsolveListener.start;
   auto indices = postsolveListener.indices;
   auto values = postsolveListener.values;
   auto origcol_mapping = postsolveListener.origcol_mapping;
   auto origrow_mapping = postsolveListener.origrow_mapping;
   auto problem = postsolveListener.problem;

   bool significant_changes_during_presolve = false;
   // Will be used during dual postsolve for fast access to bound values.
   // TODO: rows bounds are currently not updated during
   Stored<REAL> stored_bounds{ num, (int)postsolveListener.nColsOriginal,
                        (int)postsolveListener.nRowsOriginal,
                        originalSolution.type == SolutionType::kPrimalDual };

   for( int i = postsolveListener.types.size() - 1; i >= 0; --i )
   {
      auto type = types[i];
      int first = start[i];
      int last = start[i + 1];

      switch( type )
      {
      case ReductionType::kColumnDualValue:
         assert( originalSolution.type == SolutionType::kPrimalDual );
         originalSolution.reducedCosts[postsolveListener.indices[first]] =
             postsolveListener.values[postsolveListener.indices[first]];
         break;
      case ReductionType::kRowDualValue:
         assert( originalSolution.type == SolutionType::kPrimalDual );
         originalSolution.dual[postsolveListener.indices[first]] =
             postsolveListener.values[postsolveListener.indices[first]];
         break;
      case ReductionType::kFixedCol:
      {
         apply_fix_var_in_original_solution( originalSolution, indices, values,
                                             first );
         if(SolutionType::kPrimalDual == originalSolution.type)
            stored_bounds.set_bounds_of_variable( indices[first], false, false, values[first],
                                           values[first] );
         break;
      }
      case ReductionType::kFixedInfCol:
      {
         int redundant_rows = apply_fix_infinity_variable_in_original_solution(
             originalSolution, indices, values, first, problem, stored_bounds );
         // skip redundant rows because basis for those is already set and
         // should not be overwritten
         assert(are_the_next_n_types_redundant_rows( types, i, redundant_rows ));
         i -= redundant_rows;
         assert( i >= 0 );
         break;
      }
      case ReductionType::kVarBoundChange:
      {
         assert( originalSolution.type == SolutionType::kPrimalDual );
         apply_var_bound_change_forced_by_column_in_original_solution(
             originalSolution, types, start, indices, values, i, first, stored_bounds );

         bool isLowerBound = indices[first] == 1;
         int col = indices[first + 1];
         REAL old_value = values[first + 2];
         bool isInfinity = indices[first + 2] == 1;
         stored_bounds.set_bound_of_variable(col, isLowerBound, isInfinity, old_value);

         break;
      }
      case ReductionType::kSubstitutedColNoDual:
      {
         int col = indices[first];
         REAL side = values[first];
         REAL colCoef = 0.0;
         StableSum<REAL> sumcols;
         for( int j = first + 1; j < last; ++j )
         {
            if( indices[j] == col )
               colCoef = values[j];
            else
               sumcols.add( originalSolution.primal[indices[j]] * values[j] );
         }
         sumcols.add( -side );

         assert( colCoef != 0.0 );
         originalSolution.primal[col] = ( -sumcols.get() ) / colCoef;
         //TODO: modify stored objective
         break;
      }
      case ReductionType::kSubstitutedCol:
         //TODO: restore information bout the bounds
         apply_substituted_column_to_original_solution(
             originalSolution, indices, values, first, last );
         break;
      case ReductionType::kParallelCol:
         apply_parallel_col_to_original_solution( originalSolution, indices,
                                                  values, first, last, stored_bounds );
         break;
      case ReductionType::kRowBoundChangeForcedByRow:
         assert( originalSolution.type == SolutionType::kPrimalDual );
         apply_row_bound_chang_to_original_solution(
             originalSolution, types, start, indices, values, i, first);
         break;
      case ReductionType::kRowBoundChange:
      {
         assert( originalSolution.type == SolutionType::kPrimalDual );

         // if a row bound change happened because of a substitution skip the
         // verification for this step
         if( i >= 2 and types[i - 1] == ReductionType::kCoefficientChange and
             types[i - 2] == ReductionType::kSubstitutedCol )
            continue;

         bool isLhs = indices[first] == 1;
         bool isInfinity = indices[first + 1];
         int row = (int)values[first];
         REAL new_value = values[first + 1];
         switch( originalSolution.rowBasisStatus[row] )
         {
         case VarBasisStatus::FIXED:
            if( isLhs )
               originalSolution.rowBasisStatus[row] = VarBasisStatus::ON_UPPER;
            //TODO: ZERO missing
            else
               originalSolution.rowBasisStatus[row] = VarBasisStatus::ON_LOWER;
            break;
         case VarBasisStatus::ZERO:
         case VarBasisStatus::ON_LOWER:
            if( isLhs )
               originalSolution.rowBasisStatus[row] = VarBasisStatus::BASIC;
            break;
         case VarBasisStatus::ON_UPPER:
            if( not isLhs )
               originalSolution.rowBasisStatus[row] = VarBasisStatus::BASIC;
            break;
         case VarBasisStatus::BASIC:
            break;
         }

         break;
      }
      case ReductionType::kReducedBoundsCost:
      {
         significant_changes_during_presolve = true;
         assert( originalSolution.type == SolutionType::kPrimalDual );
         // get column bounds
         for( int j = 0; j < postsolveListener.origcol_mapping.size(); j++ )
         {
            int origCol = postsolveListener.origcol_mapping[j];
            int index = first + 2 * j;
            stored_bounds.set_bounds_of_variable(
                origCol, postsolveListener.indices[index] == 1,
                postsolveListener.indices[index + 1] == 1,
                postsolveListener.values[index],
                postsolveListener.values[index + 1] );
         }

         // get row bounds
         int first_row_bounds =
             first + 2 * postsolveListener.origcol_mapping.size();
         for( int k = 0; k < postsolveListener.origrow_mapping.size(); k++ )
         {
            int origRow = postsolveListener.origrow_mapping[k];
            int index = first_row_bounds + 2 * k;
            stored_bounds.set_bounds_of_row( origRow,
                                      postsolveListener.indices[index] == 1,
                                      postsolveListener.indices[index + 1] == 1,
                                      postsolveListener.values[index],
                                      postsolveListener.values[index + 1] );
         }

         break;
      }
      case ReductionType::kCoefficientChange:
      {
         // if a row bound change happened because of a substitution skip the
         // verification for this step
         if( i >= 1 and types[i - 1] == ReductionType::kSubstitutedCol )
            continue;
         break;
      }
      case ReductionType::kRedundantRow:
         assert( originalSolution.type == SolutionType::kPrimalDual );
         if(originalSolution.basisAvailabe)
            originalSolution.rowBasisStatus[indices[first]] = VarBasisStatus::BASIC;
         break;
      case ReductionType::kReasonForRowBoundChangeForcedByRow:
      case ReductionType::kSaveRow:
         assert( originalSolution.type == SolutionType::kPrimalDual );
         break;
      }

#ifndef NDEBUG
      if( reducedSolution.type == SolutionType::kPrimalDual )
      {
         Problem<REAL> problem_at_step_i =
             calculate_current_problem( postsolveListener, i );
         message.info( "Validation of partial ({}) reconstr. sol : ", i );
         validation.verifySolutionAndUpdateSlack( originalSolution, problem_at_step_i );
         assert( not significant_changes_during_presolve or
                 stored_bounds.check_bounds( problem_at_step_i ) );
      }
#endif
   }

   PostsolveStatus status =
       validation.verifySolutionAndUpdateSlack( originalSolution, problem );
   assert( not significant_changes_during_presolve or
           stored_bounds.check_bounds( problem ) );
   if( status == PostsolveStatus::kFailed )
      message.error( "Postsolving solution failed. Please use debug mode to "
                     "obtain more information." );

   return status;
}

template <typename REAL>
bool
Postsolve<REAL>::are_the_next_n_types_redundant_rows(
    const Vec<ReductionType>& types, int i, int redundant_rows ) const
{
   for( int j = 0; j < redundant_rows; j++ )
      if( types[i - j - 1] != ReductionType::kRedundantRow )
         return false;
   return true;
}

template <typename REAL>
void
Postsolve<REAL>::apply_substituted_column_to_original_solution(
    Solution<REAL>& originalSolution,
    const Vec<int>& indices, const Vec<REAL>& values, int first,
    int last ) const
{
   int row = indices[first];
   int row_length = (int)values[first];
   assert( indices[first + 1] == 0 );
   REAL lhs = values[first + 1];
   assert( lhs == values[first + 2] );
   assert( indices[first + 2] == 0 );

   int col = indices[first + 3 + row_length];

   // calculate the primal solution
   REAL colCoef = 0.0;
   StableSum<REAL> sumcols;
   for( int j = first + 3; j < first + 3 + row_length; ++j )
   {
      if( indices[j] == col )
         colCoef = values[j];
      else
         sumcols.add( originalSolution.primal[indices[j]] * values[j] );
   }
   sumcols.add( -lhs );
   assert( colCoef != 0.0 );
   originalSolution.primal[col] = ( -sumcols.get() ) / colCoef;

   // calculate the dual solution
   assert( ( originalSolution.type == SolutionType::kPrimalDual and
             values[first + 3 + row_length] > 0 ) or
           ( originalSolution.type != SolutionType::kPrimalDual and
             values[first + 3 + row_length] == 0 ) );
   if( originalSolution.type == SolutionType::kPrimalDual )
   {

      int col_length = (int)values[first + 3 + row_length];
      assert( indices[first + 4 + row_length] == 0 );
      REAL obj = values[first + 4 + row_length];

      bool ub_infinity = indices[first + 5 + row_length] == 1;
      REAL ub = values[first + 5 + row_length];
      bool lb_infinity = indices[first + 6 + row_length] == 1;
      REAL lb = values[first + 6 + row_length];

      assert( lb_infinity or ub_infinity or num.isGE( ub, lb ) );

      if( ( num.isEq( originalSolution.primal[col], lb ) and
            not lb_infinity ) or
          ( num.isEq( originalSolution.primal[col], ub ) and not ub_infinity ) )
      {
         // adjust the dual solution to the obj coefficient change and calculate
         // the reduced costs
         originalSolution.dual[row] += obj /colCoef;

         StableSum<REAL> sum_dual;
         for( int j = first + 7 + row_length; j < last; ++j )
         {
            sum_dual.add( - originalSolution.dual[indices[j]] * values[j] );
         }
         sum_dual.add( obj );

         originalSolution.reducedCosts[col] = ( sum_dual.get() );


         if( originalSolution.basisAvailabe )
         {
            if(not num.isZero(originalSolution.dual[row]) and
                originalSolution.rowBasisStatus[row] == VarBasisStatus::BASIC)
            {
               originalSolution.rowBasisStatus[row] = VarBasisStatus::FIXED;
               assert( num.isZero( originalSolution.reducedCosts[col] ) );
               originalSolution.varBasisStatus[col] = VarBasisStatus::BASIC;
            }
            else if( num.isEq( originalSolution.primal[col], lb ) and
                not lb_infinity )
            {
               if( num.isZero( lb ) and ub_infinity )
                  originalSolution.varBasisStatus[col] = VarBasisStatus::ZERO;
               else
                  originalSolution.varBasisStatus[col] =
                      VarBasisStatus::ON_LOWER;
            }
            else if( num.isEq( originalSolution.primal[col], ub ) and
                     not ub_infinity )
               originalSolution.varBasisStatus[col] = VarBasisStatus::ON_UPPER;
            else
               assert(false);
         }
      }
      else
      {
         // since variable is not at its bounds the reduced costs are 0 and calculate the dual variable
         originalSolution.reducedCosts[col] = 0;

         REAL rowCoef = 0.0;
         StableSum<REAL> sum_dual;
         for( int j = first + 7 + row_length; j < last; ++j )
         {
            if( indices[j] == row )
               rowCoef = values[j];
            else
               sum_dual.add(- originalSolution.dual[indices[j]] * values[j] );
         }

         assert( rowCoef != 0 );
         sum_dual.add( obj );
         assert( num.isZero(originalSolution.dual[row]) );
         originalSolution.dual[row] = sum_dual.get() / rowCoef;

         if( originalSolution.basisAvailabe )
         {
            if( num.isEq( originalSolution.primal[col], lb ) and
                not lb_infinity )
            {
               if( num.isZero( lb ) )
                  originalSolution.varBasisStatus[col] = VarBasisStatus::ZERO;
               else
                  originalSolution.varBasisStatus[col] =
                      VarBasisStatus::ON_LOWER;
            }
            else if( num.isEq( originalSolution.primal[col], ub ) and
                     not ub_infinity )
               originalSolution.varBasisStatus[col] = VarBasisStatus::ON_UPPER;
            else
            {
               assert( originalSolution.rowBasisStatus[row] ==
                       VarBasisStatus::BASIC );
               // TODO is this correct
               originalSolution.varBasisStatus[col] = VarBasisStatus::BASIC;
               originalSolution.rowBasisStatus[row] = VarBasisStatus::FIXED;
            }
         }
      }
      assert( row_length + col_length + 7 == last - first );

   }
}

template <typename REAL>
void
Postsolve<REAL>::apply_row_bound_chang_to_original_solution(
    Solution<REAL>& originalSolution, const Vec<ReductionType>& types,
    const Vec<int>& start, const Vec<int>& indices, const Vec<REAL>& values,
    int i, int first ) const
{
   bool isLhs = indices[first] == 1;
   int row = (int)values[first];
   REAL new_value = values[first + 1];

   int next_type = i - 1;
   int start_reason = start[next_type];
   assert( types[next_type] ==
           ReductionType::kReasonForRowBoundChangeForcedByRow );
   int remained_row = indices[start_reason];
   int deleted_row = indices[start_reason + 1];
   REAL factor = values[start_reason];
   assert( remained_row == row );
   REAL dual_row_value = originalSolution.dual[remained_row];
   if( ( isLhs and dual_row_value > 0 ) or
       ( not isLhs and dual_row_value < 0 ) )
   {
      originalSolution.dual[deleted_row] = dual_row_value * factor;
      originalSolution.dual[remained_row] = 0;
      if(originalSolution.basisAvailabe)
      {
         assert( originalSolution.rowBasisStatus[deleted_row] ==
                 VarBasisStatus::BASIC );
         originalSolution.rowBasisStatus[deleted_row] =
             originalSolution.rowBasisStatus[remained_row];
         originalSolution.rowBasisStatus[remained_row] = VarBasisStatus::BASIC;
      }
   }
}

template <typename REAL>
void
Postsolve<REAL>::apply_parallel_col_to_original_solution(
    Solution<REAL>& originalSolution, const Vec<int>& indices,
    const Vec<REAL>& values, int first, int last, Stored<REAL>& stored ) const
{
   // calculate values of the parallel cols such that at least one is at its
   // bounds
   assert( last - first == 5 );

   int col1 = indices[first];
   int col1boundFlags = indices[first + 1];
   int col2 = indices[first + 2];
   int col2boundFlags = indices[first + 3];
   const REAL& col1lb = values[first];
   const REAL& col1ub = values[first + 1];
   const REAL& col2lb = values[first + 2];
   const REAL& col2ub = values[first + 3];
   const REAL& col2scale = values[first + 4];
   const REAL& solval = originalSolution.primal[col2];

   REAL col1val;
   REAL col2val;

   if( col1boundFlags & IS_INTEGRAL )
   {
      assert( !( col1boundFlags & IS_LBINF ) );
      assert( !( col1boundFlags & IS_UBINF ) );
      assert( !( col2boundFlags & IS_LBINF ) );
      assert( !( col2boundFlags & IS_UBINF ) );
      assert( col2boundFlags & IS_INTEGRAL );

      col1val = col1lb;

      while( num.isFeasLE( col1val, col1ub ) )
      {
         col2val = solval - col1val * col2scale;

         if( num.isFeasIntegral( col2val ) && num.isFeasGE( col2val, col2lb ) &&
             num.isFeasLE( col2val, col2ub ) )
            break;

         col1val += 1;
      }
   }
   else
   {
      REAL col2valGuess;
      if( !( col2boundFlags & IS_LBINF ) )
         col2valGuess = col2lb;
      else if( !( col2boundFlags & IS_UBINF ) )
         col2valGuess = col2ub;
      else
         col2valGuess = 0;

      col1val = ( solval - col2valGuess ) / col2scale;

      // check if value for column 1 is feasible
      if( !( col1boundFlags & IS_LBINF ) && num.isFeasLT( col1val, col1lb ) )
      {
         // lower bound was violated -> set column 1 to lower bound
         col1val = col1lb;
         // compute new value for column1
         col2val = solval - col2scale * col1val;
      }
      else if( !( col1boundFlags & IS_UBINF ) &&
               num.isFeasGT( col1val, col1ub ) )
      {
         // upper bound was violated -> set column 1 to lower bound
         col1val = col1ub;
         // compute new value for column 2
         col2val = solval - col2scale * col1val;
      }
      else
      {
         // guess was feasible
         col2val = col2valGuess;
      }

      // domains should be valid now,  except for integrality of column
      // 2 which could still be violated
      assert( ( col1boundFlags & IS_LBINF ) ||
              num.isFeasGE( col1val, col1lb ) );
      assert( ( col1boundFlags & IS_UBINF ) ||
              num.isFeasLE( col1val, col1ub ) );
      assert( ( col2boundFlags & IS_LBINF ) ||
              num.isFeasGE( col2val, col2lb ) );
      assert( ( col2boundFlags & IS_UBINF ) ||
              num.isFeasLE( col2val, col2ub ) );

      // maybe integrality is violated now for column 2, then we round
      // further in the direction that we already moved to column 2 to
      if( ( col2boundFlags & IS_INTEGRAL ) && !num.isFeasIntegral( col2val ) )
      {
         // round in the direction that we moved away from when we
         // corrected the bound violation for column 1 otherwise we
         // will violate the bounds of column 1 again
         if( col2val > col2valGuess )
            col2val = ceil( col2val );
         else
            col2val = floor( col2val );

         // recompute value for column 1
         col1val = solval - col1val * col2scale;
      }
   }
   // bounds and integrality should hold now within feasibility
   // tolerance
   assert( ( col1boundFlags & IS_LBINF ) || num.isFeasGE( col1val, col1lb ) );
   assert( ( col1boundFlags & IS_UBINF ) || num.isFeasLE( col1val, col1ub ) );
   assert( !( col1boundFlags & IS_INTEGRAL ) || num.isFeasIntegral( col1val ) );
   assert( ( col2boundFlags & IS_LBINF ) || num.isFeasGE( col2val, col2lb ) );
   assert( ( col2boundFlags & IS_UBINF ) || num.isFeasLE( col2val, col2ub ) );
   assert( !( col2boundFlags & IS_INTEGRAL ) || num.isFeasIntegral( col2val ) );
   assert( num.isFeasEq( solval, col2scale * col1val + col2val ) );

   originalSolution.primal[col1] = col1val;
   originalSolution.primal[col2] = col2val;

   bool col1onBounds =
       ( not( col1boundFlags & IS_UBINF ) and num.isEq( col1val, col1ub ) ) or
       ( not( col1boundFlags & IS_LBINF ) and num.isEq( col1val, col1lb ) );
   bool col2onBounds =
       ( not( col1boundFlags & IS_UBINF ) and num.isEq( col2val, col2ub ) ) or
       ( not( col1boundFlags & IS_LBINF ) and num.isEq( col2val, col2lb ) );
   assert( col1onBounds or col2onBounds );


   if( originalSolution.type == SolutionType::kPrimalDual )
   {
      stored.set_bounds_of_variable(
          col1, ( col1boundFlags & IS_LBINF ) == IS_LBINF,
          ( col1boundFlags & IS_UBINF ) == IS_UBINF, col1lb, col1ub );
      stored.set_bounds_of_variable(
          col2, ( col2boundFlags & IS_LBINF ) == IS_LBINF,
          ( col2boundFlags & IS_UBINF ) == IS_UBINF, col2lb, col2ub );
      if( col1onBounds and col2onBounds )
      {
         if( not num.isZero( originalSolution.reducedCosts[col2] ) )
         {
            assert( num.isZero( originalSolution.reducedCosts[col1] ) );
            originalSolution.reducedCosts[col1] =
                originalSolution.reducedCosts[col2] * col2scale;
         }
         else
         {
            assert( num.isZero( originalSolution.reducedCosts[col2] ) );
            originalSolution.reducedCosts[col2] =
                originalSolution.reducedCosts[col1] / col2scale;
         }
      }

      if(originalSolution.basisAvailabe)
      {

         originalSolution.varBasisStatus[col1] =
             set_basis( col1boundFlags, col1lb, col1ub, col1val, col1onBounds );

         if( col1onBounds and col2onBounds and
             originalSolution.varBasisStatus[col2] == VarBasisStatus::BASIC )
         {
            message.template info("hier\n");
         }
         else
         {
            originalSolution.varBasisStatus[col2] = set_basis(
                col2boundFlags, col2lb, col2ub, col2val, col2onBounds );
         }
         assert(
             originalSolution.varBasisStatus[col1] != VarBasisStatus::BASIC or
             originalSolution.varBasisStatus[col2] != VarBasisStatus::BASIC );
      }
   }
}

template <typename REAL>
VarBasisStatus
Postsolve<REAL>::set_basis( int flags, REAL lb, REAL ub,
                            REAL solution, bool is_on_bounds ) const
{
   if( not is_on_bounds )
      return  VarBasisStatus::BASIC;
   else if( not( flags & IS_UBINF ) and num.isEq( solution, ub ) )
      return VarBasisStatus::ON_UPPER;
   else if( not( flags & IS_LBINF ) and
            ( flags & IS_UBINF ) and num.isEq( solution, lb ) and
            num.isZero( lb ) )
      return VarBasisStatus::ZERO;
   else if( not( flags & IS_LBINF ) and num.isEq( solution, lb ) )
      return VarBasisStatus::ON_LOWER;
   return VarBasisStatus::UNDEFINED;
}

template <typename REAL>
void
Postsolve<REAL>::apply_var_bound_change_forced_by_column_in_original_solution(
    Solution<REAL>& originalSolution, const Vec<ReductionType>& types,
    const Vec<int>& start, const Vec<int>& indices, const Vec<REAL>& values,
    int i, int first, Stored<REAL>& stored_bounds ) const
{

   bool isLowerBound = indices[first] == 1;
   int col = indices[first + 1];
   REAL old_value = values[first + 2];
   REAL new_value = values[first + 1];
   bool was_infinity = indices[first + 2] == 1;

   const REAL reduced_costs = originalSolution.reducedCosts[col];
   bool changes_neg_reduced_costs =
       not isLowerBound and num.isLT( reduced_costs, 0 );
   bool changes_pos_reduced_costs =
       isLowerBound and num.isGT( reduced_costs, 0 );

   int variables_removed_from_basis = 0;

   //calculate the reduced costs if loosens variable at a bound
   if( num.isEq( new_value, originalSolution.primal[col] ) and
       ( changes_neg_reduced_costs or changes_pos_reduced_costs ) )
   {
      assert( not num.isZero( reduced_costs ) );
      SavedRow<REAL> saved_row {num, i, types, start, indices, values, originalSolution.primal};
      int row = saved_row.getRow();
      REAL increasing_value = reduced_costs / saved_row.getCoeffOfCol( col );

      originalSolution.dual[saved_row.getRow()] += increasing_value;

      if( originalSolution.basisAvailabe and
          originalSolution.rowBasisStatus[row] == VarBasisStatus::BASIC and
          not num.isZero( originalSolution.dual[row] ) )
      {
         originalSolution.rowBasisStatus[row] = saved_row.getVBS();
         assert( originalSolution.rowBasisStatus[row] !=
                 VarBasisStatus::BASIC );
         variables_removed_from_basis++;
      }

      for( int j = 0; j < saved_row.getLength(); ++j )
      {
         int col_index = saved_row.getCoeff( j );
         if( col_index == col )
            continue;
         originalSolution.reducedCosts[col_index] -=
             increasing_value * saved_row.getValue( j );

         if( originalSolution.basisAvailabe and
             originalSolution.varBasisStatus[col_index] == VarBasisStatus::BASIC and
             not num.isZero( originalSolution.reducedCosts[col_index] ) )
         {
            originalSolution.varBasisStatus[col_index] = get_var_basis_status(
                stored_bounds, col_index, originalSolution.primal[col_index] );
            assert(originalSolution.varBasisStatus[col_index] != VarBasisStatus::BASIC);
            variables_removed_from_basis++;
         }

         assert( not originalSolution.basisAvailabe or
                 not num.isZero( originalSolution.reducedCosts[col_index] or
                                 originalSolution.varBasisStatus[col_index] !=
                                     VarBasisStatus::BASIC ) );
      }
      if( originalSolution.basisAvailabe and variables_removed_from_basis > 0 )
      {
         originalSolution.varBasisStatus[col] = VarBasisStatus::BASIC;
         variables_removed_from_basis--;
      }
      assert( variables_removed_from_basis == 0 );
      originalSolution.reducedCosts[col] = 0;
   }

   if(originalSolution.basisAvailabe){
      switch(originalSolution.varBasisStatus[col]){
         case VarBasisStatus::FIXED:
            if(isLowerBound)
               originalSolution.varBasisStatus[col] = VarBasisStatus::ON_UPPER;
            else if( was_infinity and num.isZero(new_value) )
               originalSolution.varBasisStatus[col] = VarBasisStatus::ZERO;
            else
               originalSolution.varBasisStatus[col] = VarBasisStatus::ON_LOWER;
            break;
         case VarBasisStatus::ON_LOWER:
         case VarBasisStatus::ZERO:
         {
            if( not isLowerBound )
               break;
            remove_var_from_basis( originalSolution, types, start, indices,
                                   values, i, stored_bounds );
            originalSolution.varBasisStatus[col] = VarBasisStatus::BASIC;
            break;
         }
         case VarBasisStatus::ON_UPPER:
         {
            if( isLowerBound )
               break;
            remove_var_from_basis( originalSolution, types, start, indices,
                                   values, i, stored_bounds );
            originalSolution.varBasisStatus[col] = VarBasisStatus::BASIC;
            break;
         }
         case VarBasisStatus::UNDEFINED:
         case VarBasisStatus::BASIC:
            break;
         }
   }
}

template <typename REAL>
VarBasisStatus
Postsolve<REAL>::get_var_basis_status( Stored<REAL>& stored_bounds,
                                   int index, REAL val ) const
{
   bool isOnUpperBound = stored_bounds.is_on_upper_bound( index, val );
   bool isOnLowerBound = stored_bounds.is_on_lower_bound( index, val );
   if( isOnUpperBound and isOnLowerBound )
      return VarBasisStatus::FIXED;
   else if( isOnUpperBound )
      return VarBasisStatus::ON_UPPER;
   else if( isOnLowerBound and num.isZero( val ) )
      return VarBasisStatus::ZERO;
   else if( isOnLowerBound)
      return VarBasisStatus::ON_LOWER;
   return VarBasisStatus::BASIC;
}


template <typename REAL>
void
Postsolve<REAL>::remove_var_from_basis( Solution<REAL>& originalSolution,
                                        const Vec<ReductionType>& types,
                                        const Vec<int>& start,
                                        const Vec<int>& indices,
                                        const Vec<REAL>& values, int i,
                                        Stored<REAL>& stored_bounds ) const
{
   SavedRow<REAL> saved_row {num, i, types, start, indices, values, originalSolution.primal};

   if( originalSolution.rowBasisStatus[saved_row.getRow()] != VarBasisStatus::BASIC or
       ( not saved_row.is_on_lhs() and not saved_row.is_on_lhs() ) )
   {
      for( int j = 0; j < saved_row.getLength(); ++j )
      {
         int col_index = saved_row.getCoeff(j);
         //check if variable is in Basis and can be removed
         if( originalSolution.varBasisStatus[col_index] ==
             VarBasisStatus::BASIC )
         {
            REAL sol = originalSolution.primal[col_index];
            if( stored_bounds.is_on_lower_bound( col_index, sol ) )
            {
               originalSolution.varBasisStatus[col_index] =
                   VarBasisStatus::ON_LOWER;
               return;
            }
            else if( stored_bounds.is_on_upper_bound( col_index, sol ) )
            {
               originalSolution.varBasisStatus[col_index] =
                   VarBasisStatus::ON_UPPER;
               return;
            }
         }
      }
      assert(false);
   }
   else{
      originalSolution.rowBasisStatus[saved_row.getRow()] = saved_row.getVBS();
      assert( originalSolution.rowBasisStatus[saved_row.getRow()] !=
              VarBasisStatus::BASIC );
   }
}

template <typename REAL>
void
Postsolve<REAL>::apply_fix_var_in_original_solution(
    Solution<REAL>& originalSolution, const Vec<int>& indices,
    const Vec<REAL>& values, int current_index ) const
{
   //fix variable in the primal solution
   int col = indices[current_index];
   originalSolution.primal[col] = values[current_index];

   if( originalSolution.type == SolutionType::kPrimalDual )
   {
      // calculate the reduced costs in the dual z_j = c_j - sum_i a_ij*y_i
      REAL objective_coefficient = values[current_index + 1];

      int col_length = indices[current_index + 1];

      StableSum<REAL> stablesum;
      stablesum.add( objective_coefficient );
      // no need to check for solSetRow because if it is zero then the
      // dual reducedCosts is zero.
      for( int k = 0; k < col_length; ++k )
      {
         int index = current_index + 2 + k;
         int row = indices[index];
         REAL coeff = values[index];
         stablesum.add( -coeff * originalSolution.dual[row] );
      }

      originalSolution.reducedCosts[col] = stablesum.get();
      if( originalSolution.basisAvailabe )
         originalSolution.varBasisStatus[col] = VarBasisStatus::FIXED;
   }
}

template <typename REAL>
int
Postsolve<REAL>::apply_fix_infinity_variable_in_original_solution(
    Solution<REAL>& originalSolution, Vec<int>& indices, Vec<REAL>& values,
    int first, const Problem<REAL>& problem, Stored<REAL>& stored_bounds ) const
{
   //calculate the feasible (minimal) value for the infinity variable
   int col = indices[first];
   REAL bound = values[first + 1];
   int number_rows = indices[first + 1];
   REAL solution = bound;
   int row_counter = 0;
   int current_counter = first + 2;

   bool isNegativeInfinity = values[first] < 0;
   int row_indices[number_rows];
   REAL col_coefficents[number_rows];
   if( isNegativeInfinity )
   {
      while( row_counter < number_rows )
      {
         int length = (int)values[current_counter];
         row_indices[row_counter] = indices[current_counter];

         REAL lhs = values[current_counter + 1];
         REAL rhs = values[current_counter + 2];
         const REAL* coefficients = &values[current_counter + 3];
         const int* col_indices = &indices[current_counter + 3];

         REAL newValue = calculate_row_value_for_infinity_column(
             lhs, rhs, length, col, col_indices, coefficients, originalSolution.primal,
             true, col_coefficents[row_counter]);
         if( num.isLT(newValue, solution) )
         {
            if(originalSolution.basisAvailabe)
            {
               if( num.isGT( col_coefficents[row_counter], 0 ) )
                  originalSolution.rowBasisStatus[row_indices[row_counter]] =
                      VarBasisStatus::ON_UPPER;
               else
                  originalSolution.rowBasisStatus[row_indices[row_counter]] =
                      VarBasisStatus::ON_LOWER;
            }
            solution = newValue;
         }
         else if (originalSolution.basisAvailabe)
            originalSolution.rowBasisStatus[row_indices[row_counter]] = VarBasisStatus::BASIC;

         current_counter += 3 + length;
         row_counter ++;
      }
      if( problem.getColFlags()[col].test( ColFlag::kIntegral ) )
         solution = num.epsFloor( solution );
      originalSolution.primal[col] = solution;
   }
   else
   {
      while( row_counter < number_rows )
      {
         int length = (int)values[current_counter];
         row_indices[row_counter] = indices[current_counter];

         REAL lhs = values[current_counter + 1];
         REAL rhs = values[current_counter + 2];

         const REAL* coefficients = &values[current_counter + 3];
         const int* col_indices = &indices[current_counter + 3];

         REAL newValue = calculate_row_value_for_infinity_column(
             lhs, rhs, length, col, col_indices, coefficients, originalSolution.primal,
             false, col_coefficents[row_counter] );
         if( num.isGT(newValue, solution) )
         {
            if( num.isGT( col_coefficents[row_counter], 0 ) )
               originalSolution.rowBasisStatus[row_indices[row_counter]] =
                   VarBasisStatus::ON_LOWER;
            else
               originalSolution.rowBasisStatus[row_indices[row_counter]] =
                   VarBasisStatus::ON_UPPER;
            solution = newValue;
         }
         else if (originalSolution.basisAvailabe)
            originalSolution.rowBasisStatus[row_indices[row_counter]] = VarBasisStatus::BASIC;

         current_counter += 3 + length;
         row_counter ++;
      }
      if( problem.getColFlags()[col].test( ColFlag::kIntegral ) )
         solution = num.epsCeil( solution );
      originalSolution.primal[col] = solution;
   }

   if( originalSolution.type == SolutionType::kPrimalDual )
   {
      // calculate the reduced costs for the fixed variable with the stored
      // information of calculation of the primal solution
      StableSum<REAL> sum;
      // objective is zero otherwise the problem would be infeasible
      for( int k = 0; k < number_rows; ++k )
         sum.add( -col_coefficents[k] * originalSolution.dual[row_indices[k]] );
      originalSolution.reducedCosts[col] = sum.get();

      if(isNegativeInfinity)
         stored_bounds.set_bounds_of_variable(col, true, false, 0, bound);
      else
         stored_bounds.set_bounds_of_variable(col, false, true, bound, 0);

      if( originalSolution.basisAvailabe )
      {
         if( num.isEq( solution, bound ) )
            if(isNegativeInfinity)
               originalSolution.varBasisStatus[col] = VarBasisStatus::ON_UPPER;
            else
               originalSolution.varBasisStatus[col] = VarBasisStatus::ON_LOWER;
         else
            originalSolution.varBasisStatus[col] = VarBasisStatus::BASIC;
      }
   }
   return number_rows;
}


template <typename REAL>
void
Postsolve<REAL>::copy_from_reduced_to_original(
    const Solution<REAL>& reducedSolution, Solution<REAL>& originalSolution,
    const PostsolveListener<REAL>& postsolveListener) const
{
   if( reducedSolution.type == SolutionType::kPrimalDual )
      originalSolution.type = SolutionType::kPrimalDual;

   originalSolution.primal.clear();
   originalSolution.primal.resize( postsolveListener.nColsOriginal );

   int reduced_columns = (int) reducedSolution.primal.size();
   for( int k = 0; k < reduced_columns; ++k )
      originalSolution.primal[postsolveListener.origcol_mapping[k]] =
          reducedSolution.primal[k];

   if( originalSolution.type == SolutionType::kPrimalDual )
   {
      originalSolution.basisAvailabe = reducedSolution.basisAvailabe;
      int reduced_rows = (int) reducedSolution.dual.size();

      assert( reducedSolution.dual.size() == reduced_rows );
      originalSolution.dual.clear();
      originalSolution.dual.resize( postsolveListener.nRowsOriginal );
      for( int k = 0; k < reduced_rows; k++ )
         originalSolution.dual[postsolveListener.origrow_mapping[k]] =
             reducedSolution.dual[k];

      assert( reducedSolution.reducedCosts.size() == reduced_columns );
      originalSolution.reducedCosts.clear();
      originalSolution.reducedCosts.resize( postsolveListener.nColsOriginal );
      for( int k = 0; k < reduced_columns; k++ )
         originalSolution.reducedCosts[postsolveListener.origcol_mapping[k]] =
             reducedSolution.reducedCosts[k];

      assert( reducedSolution.varBasisStatus.size() == reduced_columns );
      originalSolution.varBasisStatus.clear();
      originalSolution.varBasisStatus.resize( postsolveListener.nColsOriginal, VarBasisStatus::UNDEFINED );
      for( int k = 0; k < reduced_columns; k++ )
         originalSolution.varBasisStatus[postsolveListener.origcol_mapping[k]] =
             reducedSolution.varBasisStatus[k];

      assert( reducedSolution.rowBasisStatus.size() == reduced_rows );

      originalSolution.rowBasisStatus.clear();
      originalSolution.rowBasisStatus.resize( postsolveListener.nRowsOriginal, VarBasisStatus::UNDEFINED );
      for( int k = 0; k < reduced_rows; k++ )
         originalSolution.rowBasisStatus[postsolveListener.origrow_mapping[k]] =
             reducedSolution.rowBasisStatus[k];
   }
}



template <typename REAL>
REAL
Postsolve<REAL>::calculate_row_value_for_infinity_column(
    REAL lhs, REAL rhs, int rowLength, int column, const int* row_indices,
    const REAL* coefficients, Vec<REAL>& current_solution,
    bool is_negative, REAL& coeff_of_column_in_row ) const
{
   StableSum<REAL> stableSum;

   coeff_of_column_in_row = 0;
   for( int l = 0; l < rowLength; l++ )
   {
      int row_index = row_indices[l];
      if( row_index == column )
      {
         coeff_of_column_in_row = coefficients[l];
         continue;
      }

      stableSum.add( -coefficients[l] * current_solution[row_index] );
   }
   if( ( coeff_of_column_in_row > 0 and is_negative ) or
       ( coeff_of_column_in_row < 0 and not is_negative ) )
      stableSum.add( rhs );
   else
      stableSum.add( lhs );
   assert( coeff_of_column_in_row != 0 );
   return ( stableSum.get() / coeff_of_column_in_row );
}

template <typename REAL>
Problem<REAL>
Postsolve<REAL>::calculate_current_problem(
    const PostsolveListener<REAL>& listener, int current_index ) const
{

   auto types = listener.types;
   auto start = listener.start;
   auto indices = listener.indices;
   auto values = listener.values;
   auto origcol_mapping = listener.origcol_mapping;
   auto origrow_mapping = listener.origrow_mapping;
   auto problem = listener.problem;

   Problem<REAL> reduced = Problem<REAL>( problem );
   reduced.recomputeAllActivities();
   Statistics statistics{};
   PresolveOptions presolveOptions{};
   PostsolveListener<REAL> postsolveListener1 =
       PostsolveListener<REAL>( reduced, num );
   ProblemUpdate<REAL> problemUpdate( reduced, postsolveListener1, statistics,
                                      presolveOptions, num, message );

   for( int j = 0; j < current_index; j++ )
   {
      auto type = types[j];
      int first = start[j];
      int last = start[j + 1];
      switch( type )
      {
      case ReductionType::kRedundantRow:
         problemUpdate.markRowRedundant( indices[first] );
         break;
      case ReductionType::kFixedCol:
      case ReductionType::kFixedInfCol:
      {
         // At the moment saves column to the postsolve stack. todo:
         // use notifySavedCol and current index of column on the stack.
         int col = indices[first];
         REAL val = values[first];
         problemUpdate.getProblem().getLowerBounds()[col] = val;
         REAL obj = problemUpdate.getProblem().getObjective().coefficients[col];
         problemUpdate.getProblem().getUpperBounds()[col] = val;
         problemUpdate.getProblem().getColFlags()[col].set( ColFlag::kFixed );
         problemUpdate.addDeletedVar( col );
         problemUpdate.removeFixedCols();
         problemUpdate.clearDeletedCols();
         problemUpdate.getProblem().getObjective().coefficients[col] = obj;
         break;
      }
      case ReductionType::kVarBoundChange:
      {
         bool isLowerBound = indices[first] == 1;
         int col = indices[first + 1];
         REAL old_value = values[first + 2];
         REAL new_value = values[first + 1];
         if( isLowerBound )
         {
            problemUpdate.getProblem().getLowerBounds()[col] = new_value;
            problemUpdate.getProblem().getColFlags()[col].unset(
                ColFlag::kLbInf );
         }
         else
         {
            problemUpdate.getProblem().getUpperBounds()[col] = new_value;
            problemUpdate.getProblem().getColFlags()[col].unset(
                ColFlag::kUbInf );
         }
         break;
      }
      case ReductionType::kRowBoundChange:
      case ReductionType::kRowBoundChangeForcedByRow:
      {
         bool isLhs = indices[first] == 1;
         bool isInfinity = indices[first + 1];
         int row = (int)values[first];
         REAL new_value = values[first + 1];
         if( isLhs )
         {
            if( isInfinity )
            {
               problemUpdate.getProblem()
                   .getConstraintMatrix()
                   .template modifyLeftHandSide<true>( row, num );
            }
            else
            {
               problemUpdate.getProblem()
                   .getConstraintMatrix()
                   .modifyLeftHandSide( row, num, new_value );
            }
         }
         else
         {
            if( isInfinity )
            {
               problemUpdate.getProblem()
                   .getConstraintMatrix()
                   .template modifyRightHandSide<true>( row, num );
            }
            else
            {
               problemUpdate.getProblem()
                   .getConstraintMatrix()
                   .modifyRightHandSide( row, num, new_value );
            }
         }
         break;
      }
      case ReductionType::kCoefficientChange:
         {
//         int row = indices[first];
//         int col = indices[first + 1];
//         REAL value = values[first];
         // TODO:
         //         problemUpdate.matr
         break;
      }
      case ReductionType::kSubstitutedCol:
      case ReductionType::kSubstitutedColNoDual:
      {
         int row = indices[first];
         int row_length = (int)values[first];
         assert( indices[first + 1] == 0 );
         int col = indices[first + 3 + row_length];
         int colsize = problemUpdate.getProblem()
                           .getConstraintMatrix()
                           .getColSizes()[col];

         assert( problemUpdate.getProblem().getRowFlags()[row].test(
             RowFlag::kEquation ) );
         problemUpdate.getProblem().getColFlags()[col].set(
             ColFlag::kSubstituted );
         problemUpdate.getProblem().substituteVarInObj( num, col, row );

         if( colsize > 1 )
         {
//            auto eqRHS = problemUpdate.getProblem()
//                             .getConstraintMatrix()
//                             .getLeftHandSides()[row];

            // TODO: make the changes in the constraint matrix?
            //            problemUpdate.getProblem().getConstraintMatrix().aggregate(
            //                num, col, rowvec, eqRHS,
            //                problemUpdate.getProblem().getVariableDomains(),
            //                intbuffer, realbuffer, tripletbuffer,
            //                changed_activities, problem.getRowActivities(),
            //                singletonRows, singletonColumns, emptyColumns,
            //                stats.nrounds );
         }
         break;
      }
      case ReductionType::kParallelCol:
      {
         int col1 = indices[first];
         int col2 = indices[first + 2];
         const REAL& col2scale = values[first + 4];

         problemUpdate.merge_parallel_columns(
             col1, col2, col2scale, problemUpdate.getConstraintMatrix(),
             problemUpdate.getProblem().getLowerBounds(),
             problemUpdate.getProblem().getUpperBounds(),
             problemUpdate.getProblem().getColFlags() );
      }
      case ReductionType::kReducedBoundsCost:
      case ReductionType::kSaveRow:
      case ReductionType::kReasonForRowBoundChangeForcedByRow:
         break;
      default:
         assert( false );
      }
   }
   MpsWriter<REAL> write{};
   std::string filename = "problem_index_" +
       boost::lexical_cast<std::string>( current_index ) +
       ".mps";

   Vec<int> col_mapping{};
   Vec<int> row_mapping{};
   for( int i = 0; i < problem.getNRows(); i++ )
      row_mapping.push_back( i );
   for( int i = 0; i < problem.getNCols(); i++ )
      col_mapping.push_back( i );

   return reduced;
    }

} // namespace papilo

#endif
