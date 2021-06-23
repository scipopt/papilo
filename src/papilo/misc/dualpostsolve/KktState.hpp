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

#include "CheckLevel.hpp"
#include "KktChecker.hpp"
#include "KktInfo.hpp"
#include "kkt_status.hpp"
#include "papilo/core/Solution.hpp"

namespace papilo
{
   /// class to hold all the data needed for the checks. The other class
   /// KktChecker only holds copies of the original and reduced problem data
   template <typename REAL>
   class KktState
   {
    private:
      KktInfo<REAL> info;
      // problem and solution data. rowValues is not a reference because it is
      // calculated internally here. The other three point to the corresponding
      // values in postsolve
      const papilo::Problem<REAL>& problem;

      const papilo::Vec<REAL>& colValues;
      const papilo::Vec<REAL>& colDuals;
      const papilo::Vec<REAL>& rowDuals;
      papilo::Vec<REAL> rowValues;

      const papilo::Vec<uint8_t>& solSetCol;
      const papilo::Vec<uint8_t>& solSetRow;

      // references to problem, for easier checking
      const papilo::SparseStorage<REAL>& matrixRW;
      const papilo::Objective<REAL>& objective;
      const papilo::Vec<REAL>& rowLower;
      const papilo::Vec<REAL>& rowUpper;
      const papilo::Vec<REAL>& colLower;
      const papilo::Vec<REAL>& colUpper;

      // zero tolerance
      // todo set on initialisation and not "static"
      REAL tol = 10e-8;

    public:
      KktState( CheckLevel checker_level, const papilo::Problem<REAL>& prob,
                const papilo::Solution<REAL>& solution,
                const papilo::Vec<uint8_t>& solSetColumns,
                const papilo::Vec<uint8_t>& solSetRows )
          : problem( prob ), colValues( solution.primal ),
            colDuals( solution.reducedCosts ), rowDuals( solution.dual ),
            solSetCol( solSetColumns ), solSetRow( solSetRows ),
            matrixRW( problem.getConstraintMatrix().getConstraintMatrix() ),
            objective( problem.getObjective() ),
            rowLower( problem.getConstraintMatrix().getLeftHandSides() ),
            rowUpper( problem.getConstraintMatrix().getRightHandSides() ),
            colLower( problem.getVariableDomains().lower_bounds ),
            colUpper( problem.getVariableDomains().upper_bounds ),
            level( checker_level )
      {
         const int nRows = problem.getNRows();
         const int nCols = problem.getNCols();

         // We can have a problem with no rows and columns.
         if( nRows == 0 && nCols == 0 )
            return;

         assert( rowLower.size() == nRows );
         assert( rowUpper.size() == nRows );
         assert( colLower.size() == nCols );
         assert( colLower.size() == nCols );

         assert( solSetCol.size() == colValues.size() );
         assert( solSetCol.size() == colDuals.size() || 0 == colDuals.size() );
         assert( solSetRow.size() == rowDuals.size() || 0 == rowDuals.size() );
         assert( solSetRow.size() == rowValues.size() ||
                 0 == rowValues.size() );
      }

      KktState( CheckLevel checker_level, const papilo::Problem<REAL>& prob,
                const papilo::Solution<REAL>& solution,
                const papilo::Vec<uint8_t>& solSetColumns,
                const papilo::Vec<uint8_t>& solSetRows,
                const papilo::Vec<REAL>& rowLower_,
                const papilo::Vec<REAL>& rowUpper_,
                const papilo::Vec<REAL>& colLower_,
                const papilo::Vec<REAL>& colUpper_ )
          : problem( prob ), colValues( solution.primal ),
            colDuals( solution.reducedCosts ), rowDuals( solution.dual ),
            solSetCol( solSetColumns ), solSetRow( solSetRows ),
            matrixRW( problem.getConstraintMatrix().getConstraintMatrix() ),
            objective( problem.getObjective() ), rowLower( rowLower_ ),
            rowUpper( rowUpper_ ), colLower( colLower_ ), colUpper( colUpper_ ),
            level( checker_level )
      {
         const int nRows = problem.getNRows();
         const int nCols = problem.getNCols();

         // We can have a problem with no rows and columns.
         if( nRows == 0 && nCols == 0 )
            return;

         assert( solSetCol.size() == colValues.size() );
         assert( solSetCol.size() == colDuals.size() || 0 == colDuals.size() );
         assert( solSetRow.size() == rowDuals.size() || 0 == rowDuals.size() );
         assert( solSetRow.size() == rowValues.size() ||
                 0 == rowValues.size() );
      }

      int
      getCurrentNRows()
      {
         return std::count_if( solSetRow.begin(), solSetRow.end(),
                               []( uint8_t isset ) { return isset; } );
      }

      int
      getCurrentNCols()
      {
         return std::count_if( solSetCol.begin(), solSetCol.end(),
                               []( uint8_t isset ) { return isset; } );
      }

      kkt_status
      checkPrimalBounds()
      {
         kkt_status status = kkt_status::OK;
         int reduced_index = 0;
         for( unsigned int col = 0; col < solSetCol.size(); col++ )
         {
            if( !solSetCol[col] )
               continue;

            int index =
                ( problem.getNCols() < solSetCol.size() ) ? reduced_index : col;
            if( !isLowerBoundInfinity( problem, index ) &&
                colLower[col] - colValues[col] > tol )
            {
               message.info( "Column {:<3} violates lower column bound.\n",
                             col );
               REAL value = colLower[col] - colValues[col];
               info.primal_col_bounds.values.push_back( value );
               status = kkt_status::Fail_Primal_Bound;
            }

            if( !isUpperBoundInfinity( problem, index ) &&
                colValues[col] - colUpper[col] > tol )
            {
               message.info( "Column {:<3} violates upper column bound.\n",
                             col );
               REAL value = colValues[col] - colUpper[col];
               info.primal_col_bounds.values.push_back( value );
               status = kkt_status::Fail_Primal_Bound;
            }

            reduced_index++;
         }
         return status;
      }

      // called when dual values are being checked too
      kkt_status
      checkLength()
      {
         const int nCols = solSetCol.size();

         if( colLower.size() != nCols || colUpper.size() != nCols ||
             colValues.size() != nCols )
            return kkt_status::Fail_Length;

         const int nRows = solSetRow.size();

         if( rowLower.size() != nRows || rowUpper.size() != nRows ||
             rowValues.size() != nRows )
            return kkt_status::Fail_Length;

         if( level == CheckLevel::After_each_step_and_dual ||
             level == CheckLevel::Primal_and_dual )
            if( colDuals.size() != nCols || rowDuals.size() != nRows )
               return kkt_status::Fail_Length;

         return kkt_status::OK;
      }

      // performs the check with respect to the (postsolved) reduced problem.
      // For the check of primal feasibility of the original problem see
      // CheckPrimalFeasibilityOriginal()
      kkt_status
      checkPrimalFeasibility()
      {
         kkt_status status = checkPrimalBounds();

         int reduced_index = 0;
         for( int row = 0; row < solSetRow.size(); row++ )
         {
            if( !solSetRow[row] )
               continue;

            int index =
                ( problem.getNRows() < solSetRow.size() ) ? reduced_index : row;
            if( !rowLHSInf( problem, index ) &&
                rowLower[row] - rowValues[row] > tol )
            {
               message.info( "Row {:<3} violates row bounds.\n", row );
               REAL value = rowValues[row] - rowUpper[row];
               info.primal_row_bounds.values.push_back( value );
               status = kkt_status::Fail_Primal_Feasibility;
            }
            if( !rowRHSInf( problem, index ) &&
                rowValues[row] - rowUpper[row] > tol )
            {
               message.info( "Row {:<3} violates row bounds.\n", row );
               REAL value = rowValues[row] - rowUpper[row];
               info.primal_row_bounds.values.push_back( value );
               status = kkt_status::Fail_Primal_Feasibility;
            }
            reduced_index++;
         }
         return status;
      }

      kkt_status
      checkDualFeasibility()
      {
         const papilo::Vec<REAL>& lowerBounds =
             problem.getVariableDomains().lower_bounds;
         const papilo::Vec<REAL>& upperBounds =
             problem.getVariableDomains().upper_bounds;

         // check values of z_j are dual feasible
         int reduced_index = 0;
         for( int col = 0; col < solSetCol.size(); col++ )
         {
            if( !solSetCol[col] )
               continue;

            // no lower or upper bound on infinity
            if( isLowerBoundInfinity( problem, reduced_index ) &&
                isUpperBoundInfinity( problem, reduced_index ) )
            {
               if( abs( colDuals[col] ) > tol )
                  return kkt_status::Fail_Dual_Feasibility;
            }
            // non fixed variable at lower bound: x=l and l<u
            else if( abs( colValues[col] - lowerBounds[col] ) < tol &&
                     lowerBounds[col] < upperBounds[col] )
            {
               if( colDuals[col] < 0 && abs( colDuals[col] ) > tol )
                  return kkt_status::Fail_Dual_Feasibility;
            }
            // non fixed variable at upper bound: x=u and l<u
            else if( abs( colValues[col] - upperBounds[col] ) < tol &&
                     lowerBounds[col] < upperBounds[col] )
            {
               if( colDuals[col] > tol )
                  return kkt_status::Fail_Dual_Feasibility;
            }
            reduced_index++;
         }

         std::vector<int> orig_col_index( matrixRW.getNCols(), 0 );
         //      int reduced_col_index = count( solSetCol.begin(), solSetCol.end(), true );
         int reduced_col_index = 0;
         for( int col = 0; col < solSetCol.size(); col++ )
         {
            if( !solSetCol[col] )
               continue;
            orig_col_index[reduced_col_index] = col;
            reduced_col_index++;
         }

         //      TODO
//         assert( matrixRW.getNCols() == reduced_col_index );

         // check values of y_i are dual feasible
         const papilo::Vec<REAL>& lhs =
             problem.getConstraintMatrix().getLeftHandSides();
         const papilo::Vec<REAL>& rhs =
             problem.getConstraintMatrix().getRightHandSides();

         for( int row = 0; row < solSetRow.size(); row++ )
         {
            if( !solSetRow[row] )
               continue;

            // L = Ax = U can be any sign
            if( abs( lhs[row] - rowValues[row] ) < tol &&
                abs( rhs[row] - rowValues[row] ) < tol )
               continue;
            // L = Ax < U
            else if( abs( lhs[row] - rowValues[row] ) < tol &&
                     rowValues[row] < rhs[row] )
            {
               if( rowDuals[row] < -tol )
                  return kkt_status::Fail_Dual_Feasibility;
            }
            // L < Ax = U
            else if( lhs[row] < rowValues[row] &&
                     abs( rowValues[row] - rhs[row] ) < tol )
            {
               if( rowDuals[row] > tol )
                  return kkt_status::Fail_Dual_Feasibility;
            }
            // L < Ax < U
            else if( lhs[row] < ( rowValues[row] + tol ) &&
                     rowValues[row] < ( rhs[row] + tol ) )
            {
               if( abs( rowDuals[row] ) > tol )
                  return kkt_status::Fail_Dual_Feasibility;
            }
         }

         return kkt_status::OK;
      }

      kkt_status
      checkComplementarySlackness()
      {
         int reduced_index = 0;
         for( int col = 0; col < solSetCol.size(); col++ )
         {
            if( !solSetCol[col] )
               continue;

            if( !isLowerBoundInfinity( problem, reduced_index ) )
               if( abs( ( colValues[col] - colLower[col] ) *
                        ( colDuals[col] ) ) > tol &&
                   colValues[col] != colUpper[col] &&
                   abs( colDuals[col] ) > tol )
                  return kkt_status::Fail_Complementary_Slackness;

            if( !isUpperBoundInfinity( problem, reduced_index ) )
               if( abs( ( colUpper[col] - colValues[col] ) *
                        ( colDuals[col] ) ) > tol &&
                   colValues[col] != colLower[col] &&
                   abs( colDuals[col] ) > tol )
                  return kkt_status::Fail_Complementary_Slackness;
         }
         reduced_index++;

         return kkt_status::OK;
      }

      kkt_status
      checkStOfLagrangian()
      {
         // A'y + z = c
         REAL lagrV;

         // getMatrixTranspose() doesn't work because sometimes the matrix of the reduced problem has not been transposed previously. Edit when KKT check is done at each step, rather than just the start and the end of postsolve.
         const papilo::SparseStorage<REAL>& matrixCW =
             problem.getConstraintMatrix().getMatrixTranspose();
         // const SparseStorage<REAL>& matrixCW =
         //    problem.getConstraintMatrix().getConstraintMatrix().getTranspose();

         std::vector<int> orig_row_index( matrixCW.getNCols(), 0 );
         int reduced_row_index = 0;
         for( int row = 0; row < solSetRow.size(); row++ )
         {
            if( !solSetRow[row] )
               continue;
            orig_row_index[reduced_row_index] = row;
            reduced_row_index++;
         }

         // Below is used getNCols() because matrixCW is the transpose.
         assert( matrixCW.getNCols() == reduced_row_index ||
                 matrixCW.getNCols() == solSetRow.size() );

         int reduced_index = 0;
         for( int col = 0; col < solSetCol.size(); col++ )
         {
            if( !solSetCol[col] )
               continue;

            lagrV = 0;

            auto index_range = matrixCW.getRowRanges()[reduced_index];
            for( int k = index_range.start; k < index_range.end; k++ )
            {
               int row = matrixCW.getColumns()[k];
               lagrV += rowDuals[orig_row_index[row]] * matrixCW.getValues()[k];
            }

            lagrV =
                lagrV + colDuals[col] - objective.coefficients[reduced_index];

            if( abs( lagrV ) > tol )
               return kkt_status::Fail_Stationarity_Lagrangian;

            reduced_index++;
         }

         return kkt_status::OK;
      }

      void
      getRowValues()
      {
         const int nRows = problem.getNRows();
         const int nCols = problem.getNCols();

         rowValues.resize( solSetRow.size() );
         REAL rowValue;

         int reduced_row_index = 0;
         auto values = matrixRW.getValues();

         std::vector<int> orig_col_index;
         orig_col_index.reserve( nCols );
         for( int i = 0; i < solSetCol.size(); i++ )
            if( solSetCol[i] )
               orig_col_index.push_back( i );

         bool reduced = false;
         if( nRows < solSetRow.size() || nCols < solSetCol.size() )
         {
            reduced = true;
            assert( orig_col_index.size() == nCols );
         }

         for( int row = 0; row < solSetRow.size(); row++ )
         {
            if( !solSetRow[row] )
               continue;

            rowValue = 0;
            papilo::IndexRange index_range;
            if( reduced )
               index_range = matrixRW.getRowRanges()[reduced_row_index];
            else
               index_range = matrixRW.getRowRanges()[row];

            for( int j = index_range.start; j < index_range.end; j++ )
            {
               int col = matrixRW.getColumns()[j];
               assert( col >= 0 );
               assert( col < (int)nCols );

               // col is index of the column in reduced problem and colValues is
               // expanded.
               if( reduced )
                  rowValue += values[j] * colValues[orig_col_index[col]];
               else
                  rowValue += values[j] * colValues[col];
            }
            rowValues[row] = rowValue;
            reduced_row_index++;
         }
      }

      void
      updateInfo()
      {
         updateKktInfo( getCurrentNCols(), getCurrentNRows(), info );
      }

      void
      reportKktInfo()
      {
         int on = getCurrentNRows();
         int count = solSetRow.size();
         message.info( "KKT check\n" );
         message.info( "         solSetRow: {:>3} / {:>3}\n", on, count );
         on = getCurrentNCols();
         count = solSetCol.size();
         message.info( "         solSetCol: {:>3} / {:>3} \n", on, count );
         message.info( "                    vio / chk    max, sum    \n" );
         message.info( "Primal col bounds:  {:>3} / {:>3}  \n",
                       info.primal_col_bounds.num_violated, info.num_col );
         //               info.primal_col_bounds.max, info.primal_col_bounds.sum );
         message.info( "Primal row bounds:  {:>3} / {:>3}  \n",
                       info.primal_row_bounds.num_violated, info.num_row );

         // todo: dual.
         message.info( "\n" );
      }

      papilo::Message message;
      CheckLevel level;
   };
}
#include "KktCheckHelper.hpp"
#include "KktInfo.hpp"
#include "papilo/core/Problem.hpp"
#include "papilo/core/ProblemBuilder.hpp"
#include "papilo/core/RowFlags.hpp"
#include "papilo/core/Solution.hpp"
#include "papilo/core/SparseStorage.hpp"
#include "papilo/io/Message.hpp"
#include "papilo/misc/VectorUtils.hpp"
#include "papilo/misc/fmt.hpp"
#include <iostream>
