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
                                            bool is_negative ) const;
   void
   verify_current_solution( const Solution<REAL>& originalSolution,
                            PrimalDualSolValidation<REAL>& validation,
                            const PostsolveListener<REAL>& listener,
                            int current_index ) const;
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
   Vec<REAL>& origSol = originalSolution.primal;

   if( reducedSolution.type == SolutionType::kPrimalDual )
   {
      originalSolution.type = SolutionType::kPrimalDual;
   }

   origSol.clear();
   origSol.resize( postsolveListener.nColsOriginal );

   for( int k = 0; k < reducedSol.size(); ++k )
   {
      int origcol = postsolveListener.origcol_mapping[k];
      origSol[origcol] = reducedSol[k];
   }

   if( originalSolution.type == SolutionType::kPrimalDual )
   {
      assert( reducedSolution.reducedCosts.size() ==
              postsolveListener.origcol_mapping.size() );
      originalSolution.reducedCosts.clear();
      originalSolution.reducedCosts.resize( postsolveListener.nColsOriginal );
      for( int k = 0; k < postsolveListener.origcol_mapping.size(); k++ )
      {
         int origcol = postsolveListener.origcol_mapping[k];
         originalSolution.reducedCosts[origcol] =
             reducedSolution.reducedCosts[k];
      }

      assert( reducedSolution.dual.size() ==
              postsolveListener.origrow_mapping.size() );
      originalSolution.dual.clear();
      originalSolution.dual.resize( postsolveListener.nRowsOriginal );
      for( int k = 0; k < postsolveListener.origrow_mapping.size(); k++ )
      {
         int origrow = postsolveListener.origrow_mapping[k];
         originalSolution.dual[origrow] = reducedSolution.dual[k];
      }
   }

   // If problem has been reduced, check solution of reduced problem returned by
   // solver.
   // At the moment not all row and column changes are notified. The check
   // below handles the case when some trivial presolve elimination is applied,
   // but types.size() is still zero.
   if( postsolveListener.origrow_mapping.size() <
           postsolveListener.nRowsOriginal ||
       postsolveListener.origcol_mapping.size() <
           postsolveListener.nColsOriginal )
   {
      // TODO: verify solution
      //      verifySolution(reducedSolution, r)
   }

   // Will be used during dual postsolve for fast access to bound values.
   Vec<REAL> col_cost;
   Vec<REAL> col_lower;
   Vec<REAL> col_upper;
   Vec<REAL> row_lhs;
   Vec<REAL> row_rhs;

   Vec<int> col_infinity_lower;
   Vec<int> col_infinity_upper;
   Vec<int> row_infinity_lhs;
   Vec<int> row_infinity_rhs;

   Vec<int> col_lower_from_row;
   Vec<int> col_upper_from_row;
   Vec<int> row_lower_from_col;
   Vec<int> row_upper_from_col;

   if( originalSolution.type == SolutionType::kPrimalDual )
   {

      col_cost.assign( postsolveListener.nColsOriginal, 0 );
      col_lower.assign( postsolveListener.nColsOriginal, 0 );
      col_upper.assign( postsolveListener.nColsOriginal, 0 );
      row_lhs.assign( postsolveListener.nRowsOriginal, 0 );
      row_rhs.assign( postsolveListener.nRowsOriginal, 0 );
      col_infinity_upper.assign( postsolveListener.nColsOriginal, 1 );
      col_infinity_lower.assign( postsolveListener.nColsOriginal, 1 );
      row_infinity_rhs.assign( postsolveListener.nRowsOriginal, 1 );
      row_infinity_lhs.assign( postsolveListener.nRowsOriginal, 1 );

      col_lower_from_row.assign( postsolveListener.nColsOriginal, -1 );
      col_upper_from_row.assign( postsolveListener.nColsOriginal, -1 );
      row_lower_from_col.assign( postsolveListener.nRowsOriginal, -1 );
      row_upper_from_col.assign( postsolveListener.nRowsOriginal, -1 );
   }

   auto types = postsolveListener.types;
   auto start = postsolveListener.start;
   auto indices = postsolveListener.indices;
   auto values = postsolveListener.values;
   auto origcol_mapping = postsolveListener.origcol_mapping;
   auto origrow_mapping = postsolveListener.origrow_mapping;
   auto problem = postsolveListener.problem;
   for( int i = postsolveListener.types.size() - 1; i >= 0; --i )
   {
      auto type = types[i];
      int first = start[i];
      int last = start[i + 1];

      switch( type )
      {
         // TODO: stack information
      case ReductionType::kReducedBoundsCost:
      {
         assert(originalSolution.type == SolutionType::kPrimalDual);
         // get column bounds
         for( int j = 0; j < postsolveListener.origcol_mapping.size(); j++ )
         {
            int origCol = postsolveListener.origcol_mapping[j];
            int index = first + 2 * j;
            col_lower[origCol] = postsolveListener.values[index];
            col_upper[origCol] = postsolveListener.values[index + 1];
            col_infinity_lower[origCol] = postsolveListener.indices[index];
            col_infinity_upper[origCol] = postsolveListener.indices[index + 1];
         }

         // get row bounds
         int first_row_bounds =
             first + 2 * postsolveListener.origcol_mapping.size();
         for( int k = 0; k < postsolveListener.origrow_mapping.size(); k++ )
         {
            int origRow = postsolveListener.origrow_mapping[k];
            int index = first_row_bounds + 2 * k;
            row_lhs[origRow] = postsolveListener.values[index];
            row_rhs[origRow] = postsolveListener.values[index + 1];
            row_infinity_lhs[origRow] = postsolveListener.indices[index];
            row_infinity_rhs[origRow] = postsolveListener.indices[index + 1];
         }

         // get cost
         int first_cost =
             first_row_bounds + 2 * postsolveListener.origrow_mapping.size();
         for( int j = 0; j < postsolveListener.origcol_mapping.size(); j++ )
         {
            int origcol = postsolveListener.origcol_mapping[j];
            col_cost[origcol] = postsolveListener.values[first_cost + j];
            assert( j == postsolveListener.indices[first_cost + j] );
         }
         break;
      }
      case ReductionType::kColumnDualValue:
         assert(originalSolution.type == SolutionType::kPrimalDual);
         originalSolution.reducedCosts[postsolveListener.indices[first]] =
             postsolveListener.values[postsolveListener.indices[first]];
         break;
      case ReductionType::kRedundantRow:
         assert(originalSolution.type == SolutionType::kPrimalDual);
         break;
      case ReductionType::kDeletedCol:
         assert(originalSolution.type == SolutionType::kPrimalDual);
         break;
      case ReductionType::kRowDualValue:
         assert(originalSolution.type == SolutionType::kPrimalDual);
         originalSolution.dual[postsolveListener.indices[first]] =
             postsolveListener.values[postsolveListener.indices[first]];
         break;
      case ReductionType::kSaveRow:
      {
         assert(originalSolution.type == SolutionType::kPrimalDual);
         // TODO I think this and the SAVE_COL step should just be skipped
         //     we only want to restore redundant rows that have been removed
         //     not all saved rows. Saved row has no logical implications for
         //     postsolve it should just give the deleted row a new index within
         //     the postsolve structure this index would be used by the undo
         //     redundant col I guess.
         int row = indices[first];
         int length = (int)values[first];
         bool lhs_inf = indices[first + 1] == true;
         bool lhs = values[first + 1] == true;
         bool rhs_inf = indices[first + 2] == true;
         bool rhs = values[first + 2] == true;
         REAL coeff[length];
         int col_indices[length];
         for( int j = 0; j < length; ++j )
         {
            col_indices[j] = indices[first + 3 + j];
            coeff[j] = values[first + 3 + j];
         }
         break;
      }
      case ReductionType::kVarBoundChangeForced:
      case ReductionType::kSaveCol:
         assert(originalSolution.type == SolutionType::kPrimalDual);
         break;
      case ReductionType::kFixedCol:
      {
         // At the moment saves column to the postsolve stack. todo:
         // use notifySavedCol and current index of column on the stack.
         int col = indices[first];
         origSol[col] = values[first];

         // todo: checker notify dual value if changed
         if( originalSolution.type == SolutionType::kPrimalDual )
         {
            // get dual reducedCosts z_j = c_j - sum_i a_ij*y_i
            REAL objective_coefficient = values[first + 1];

            // modify changes
            col_cost[col] = objective_coefficient;
            col_infinity_lower[col] = false;
            col_infinity_upper[col] = false;
            col_upper[col] = values[first];
            col_lower[col] = values[first];

            int col_length = indices[first + 1];

            REAL reducedCosts = objective_coefficient;
            // no need to check for solSetRow because if it is zero then the
            // dual reducedCosts is zero.
            for( int k = 0; k < col_length; ++k )
            {
               int index = first + 2 + k;
               int row = indices[index];
               REAL coeff = values[index];
               reducedCosts -= coeff * originalSolution.dual[row];

               if( row_infinity_rhs[row] == false )
                  row_rhs[row] += coeff * origSol[col];
               if( row_infinity_lhs[row] == false )
                  row_lhs[row] += coeff * origSol[col];
            }

            originalSolution.reducedCosts[col] = reducedCosts;
         }
         break;
      }
      case ReductionType::kSingletonRow:
      {
         assert(originalSolution.type == SolutionType::kPrimalDual);

            int row = indices[first];

            row_lhs[row] = values[first + 3];
            row_rhs[row] = values[first + 2];
            row_infinity_lhs[row] = indices[first + 3] == 1;
            row_infinity_rhs[row] = indices[first + 2] == 1;

            //            int col = indices[first + 1];
            //         // code below saves column on stack for the calculation
            //         of dual values.
            //            REAL coeff = values[first + 1];
            //            REAL cost = values[first + 4];
            //            assert( indices[first + 1] == col );
            //            REAL value = cost;
            //            int col_length_minus_one = indices[first + 4];
            // no need to check for solSetRow because if it is zero then the
            // dual value is zero.
            //            for( int k = 0; k < col_length_minus_one; ++k )
            //            {
            //               int index = first + 5 + k;
            //               int col_index = indices[index];
            //               value -= values[index] *
            //               originalSolution.dual[col_index];
            //            }
            //
            //            value = value / coeff;
            //            originalSolution.dual[row] = value;
            //            originalSolution.reducedCosts[col] = 0;
            //            originalSolution.reducedCosts[col] -= coeff * value;


         break;
      }
      case ReductionType::kVarBoundChange:
      {
         assert( originalSolution.type == SolutionType::kPrimalDual );
         bool isLowerBound = indices[first] == 1;
         int col = indices[first + 1];
         REAL old_value = values[first + 2];
         REAL new_value = values[first + 1];
         bool isInfinity = indices[first + 2] == 1;
         if( isLowerBound )
         {
            col_lower[col] = old_value;
            col_infinity_lower[col] = isInfinity;
         }
         else
         {
            col_upper[col] = old_value;
            col_infinity_upper[col] = isInfinity;
         }

         const REAL reduced_costs = originalSolution.reducedCosts[col];
         bool changes_neg_reduced_costs =
             not isLowerBound and num.isLT( reduced_costs, 0 );
         bool changes_pos_reduced_costs =
             isLowerBound and num.isGT( reduced_costs, 0 );

         if( num.isEq( new_value, originalSolution.primal[col] ) and
             ( changes_neg_reduced_costs or changes_pos_reduced_costs ) )
         {
            assert( not num.isZero( reduced_costs ) );
            int next_saved = i - 1;
            assert( types[next_saved] == ReductionType::kSaveRow );
            int first_saved_row = start[next_saved];
            int row = indices[first_saved_row];
            int length = (int)values[first_saved_row];
            for( int j = 0; j < length; ++j )
            {
               int col_index = indices[first_saved_row + 3 + j];
               if( col_index != col )
                  originalSolution.reducedCosts[col_index] -= reduced_costs;
            }

            originalSolution.dual[row] += reduced_costs;
            originalSolution.reducedCosts[col] = 0;
         }
         break;
      }
      case ReductionType::kFixedInfCol:
      {
         // TODO dual case missing
         int column = indices[first];
         REAL bound = values[first + 1];
         REAL solution = bound;

         bool isNegativeInfinity = values[first] < 0;
         if( isNegativeInfinity )
         {
            int current_row_counter = first + 2;
            while( current_row_counter != last )
            {
               int length = (int)values[current_row_counter];

               REAL lhs = values[current_row_counter + 1];
               REAL rhs = values[current_row_counter + 2];
               const REAL* coefficients = &values[current_row_counter + 3];
               const int* row_indices = &indices[current_row_counter + 3];

               REAL newValue = calculate_row_value_for_infinity_column(
                   lhs, rhs, length, column, row_indices, coefficients, origSol,
                   true );
               if( newValue < solution )
               {
                  solution = newValue;
               }
               current_row_counter += 3 + length;
            }
            if( problem.getColFlags()[column].test(
                    papilo::ColFlag::kIntegral ) )
            {
               solution = num.epsFloor( solution );
            }
            origSol[column] = solution;
         }
         else
         {
            int current_row_counter = first + 2;
            while( current_row_counter != last )
            {
               int length = (int)values[current_row_counter];

               REAL lhs = values[current_row_counter + 1];
               REAL rhs = values[current_row_counter + 2];
               const REAL* coefficients = &values[current_row_counter + 3];
               const int* row_indices = &indices[current_row_counter + 3];

               REAL newValue = calculate_row_value_for_infinity_column(
                   lhs, rhs, length, column, row_indices, coefficients, origSol,
                   false );
               if( newValue > solution )
               {
                  solution = newValue;
               }
               current_row_counter += 3 + length;
            }
            if( problem.getColFlags()[column].test(
                    papilo::ColFlag::kIntegral ) )
            {
               solution = num.epsCeil( solution );
            }
            origSol[column] = solution;
         }
         break;
      }
      case ReductionType::kSubstitutedColShort:
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
               sumcols.add( origSol[indices[j]] * values[j] );
         }
         sumcols.add( -side );

         assert( colCoef != 0.0 );
         origSol[col] = ( -sumcols.get() ) / colCoef;
         break;
      }
      case ReductionType::kSubstitutedCol:
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
               sumcols.add( origSol[indices[j]] * values[j] );
         }
         sumcols.add( -lhs );
         assert( colCoef != 0.0 );
         origSol[col] = ( -sumcols.get() ) / colCoef;

         //calculate the dual solution
         assert( (originalSolution.type == SolutionType::kPrimalDual and
                 values[first + 3 + row_length] > 0) or
                 (originalSolution.type != SolutionType::kPrimalDual and
                 values[first + 3 + row_length] == 0 ));
         if( originalSolution.type == SolutionType::kPrimalDual )
         {

            int col_length = (int)values[first + 3 + row_length];
            assert( indices[first + 4 + row_length] == 0 );
            REAL obj = values[first + 4 + row_length];

            REAL rowCoef = 0.0;
            StableSum<REAL> sum_dual;
            for( int j = first + 5 + row_length; j < last; ++j )
            {
               if( indices[j] == row )
                  rowCoef = values[j];
               else
                  sum_dual.add( originalSolution.dual[indices[j]] * values[j] );
            }
            assert(rowCoef != 0);
            sum_dual.add( -obj );
            originalSolution.dual[row] = ( -sum_dual.get() ) / rowCoef;
            assert( row_length + col_length + 5 == last - first );
         }
         break;
      }
      case ReductionType::kParallelCol:
      {
         constexpr int IS_INTEGRAL = static_cast<int>( ColFlag::kIntegral );
         constexpr int IS_LBINF = static_cast<int>( ColFlag::kLbInf );
         constexpr int IS_UBINF = static_cast<int>( ColFlag::kUbInf );

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
         const REAL& solval = origSol[col2];

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

               if( num.isFeasIntegral( col2val ) &&
                   num.isFeasGE( col2val, col2lb ) &&
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
            if( !( col1boundFlags & IS_LBINF ) &&
                num.isFeasLT( col1val, col1lb ) )
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
            if( ( col2boundFlags & IS_INTEGRAL ) &&
                !num.isFeasIntegral( col2val ) )
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
         assert( ( col1boundFlags & IS_LBINF ) ||
                 num.isFeasGE( col1val, col1lb ) );
         assert( ( col1boundFlags & IS_UBINF ) ||
                 num.isFeasLE( col1val, col1ub ) );
         assert( !( col1boundFlags & IS_INTEGRAL ) ||
                 num.isFeasIntegral( col1val ) );
         assert( ( col2boundFlags & IS_LBINF ) ||
                 num.isFeasGE( col2val, col2lb ) );
         assert( ( col2boundFlags & IS_UBINF ) ||
                 num.isFeasLE( col2val, col2ub ) );
         assert( !( col2boundFlags & IS_INTEGRAL ) ||
                 num.isFeasIntegral( col2val ) );
         assert( num.isFeasEq( solval, col2scale * col1val + col2val ) );

         // solSet[col1] = true;
         origSol[col1] = col1val;
         origSol[col2] = col2val;

         break;
      }
      case ReductionType::kRowBoundChange:
      {
         assert( originalSolution.type == SolutionType::kPrimalDual );

         bool isLhs = indices[first] == 1;
         bool isInfinity = indices[first + 1];
         int row = (int)values[first];
         REAL new_value = values[first + 1];
         if( isLhs )
         {
            if( isInfinity )
            {
               row_infinity_lhs[row] = true;
            }
            else
            {
               row_infinity_lhs[row] = false;
               row_lhs[row] = new_value;
            }
         }
         else
         {
            if( isInfinity )
            {
               row_infinity_rhs[row] = true;
            }
            else
            {
               row_infinity_rhs[row] = false;
               row_rhs[row] = new_value;
            }
         }
      }
      }

//#todo only in debug mode
#ifndef NDEBUG
      if( reducedSolution.type == SolutionType::kPrimalDual )
         verify_current_solution( originalSolution, validation,
                                  postsolveListener, i );
#endif
   }

   PostsolveStatus status =
       validation.verifySolution( originalSolution, problem );
   if( status == PostsolveStatus::kFailed )
      message.error( "Postsolving solution failed. Please use debug mode to "
                     "obtain more information." );
   return status;
}
template <typename REAL>
void
Postsolve<REAL>::verify_current_solution(
    const Solution<REAL>& originalSolution,
    PrimalDualSolValidation<REAL>& validation,
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
      case ReductionType::kSingletonRow:
      {
         //         int row = indices[first];
         //         problemUpdate.removeSingletonRow( row );
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
      {
         bool isLhs = indices[first] == 1;
         bool isInfinity = indices[first + 1];
         int row = (int)values[first];
         REAL new_value = values[first + 1];
         if( isLhs )
         {
            if( isInfinity )
            {
               problemUpdate.getProblem().getConstraintMatrix()
                   .template modifyLeftHandSide<true>( row, num );
            }
            else
            {
               problemUpdate.getProblem().getConstraintMatrix().modifyLeftHandSide(
                   row, num, new_value );
            }
         }
         else
         {
            if( isInfinity )
            {
               problemUpdate.getProblem().getConstraintMatrix()
                   .template modifyRightHandSide<true>( row, num );
            }
            else
            {
               problemUpdate.getProblem().getConstraintMatrix().modifyRightHandSide(
                   row, num, new_value );
            }
         }
         break;
      }
      case ReductionType::kSubstitutedCol:
      case ReductionType::kSubstitutedColShort:
      {
         int row = indices[first];
         int row_length = (int)values[first];
         assert( indices[first + 1] == 0 );
         int col = indices[first + 3 + row_length];

         assert(problemUpdate.getProblem().getRowFlags()[row].test(RowFlag::kEquation));
         problemUpdate.getProblem().substituteVarInObj( num, col, row );
         break;
      }
      case ReductionType::kParallelCol:
      case ReductionType::kFixedInfCol:
      case ReductionType::kReducedBoundsCost:
      case ReductionType::kColumnDualValue:
      case ReductionType::kDeletedCol:
      case ReductionType::kRowDualValue:
      case ReductionType::kSaveRow:
      case ReductionType::kSaveCol:
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
   //   write.writeProb( filename, reduced, row_mapping, col_mapping);
   message.info( "Validation of partial ({}) reconstr. sol : ", current_index );
   validation.verifySolution( originalSolution, reduced );
}

template <typename REAL>
REAL
Postsolve<REAL>::calculate_row_value_for_infinity_column(
    REAL lhs, REAL rhs, int rowLength, int column, const int* row_indices,
    const REAL* coefficients, Vec<REAL>& current_solution,
    bool is_negative ) const
{
   StableSum<REAL> stableSum;
   if( ( coefficients[column] > 0 && is_negative ) ||
       ( coefficients[column] < 0 && !is_negative ) )
      stableSum.add( rhs );
   else
      stableSum.add( lhs );
   REAL coeff_of_column_in_row = 0;
   for( int l = 0; l < rowLength; l++ )
   {
      int row_index = row_indices[l];
      if( row_index == column )
      {
         coeff_of_column_in_row = coefficients[l];
         continue;
      }
      // TODO: think about what to do if there are to infinity values ->
      // irrelevant?
      stableSum.add( -coefficients[l] * current_solution[row_index] );
   }
   assert( coeff_of_column_in_row != 0 );
   return ( stableSum.get() / coeff_of_column_in_row );
}

} // namespace papilo

#endif
