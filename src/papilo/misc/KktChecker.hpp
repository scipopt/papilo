
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

#ifndef _PAPILO_MISC_KKT_CHECKER_HPP_
#define _PAPILO_MISC_KKT_CHECKER_HPP_

#include <iostream>

#include "papilo/core/Problem.hpp"
#include "papilo/core/ProblemBuilder.hpp"
#include "papilo/core/RowFlags.hpp"
#include "papilo/core/Solution.hpp"
#include "papilo/core/SparseStorage.hpp"
#include "papilo/io/Message.hpp"
#include "papilo/misc/KktCheckHelper.hpp"
#include "papilo/misc/VectorUtils.hpp"
#include "papilo/misc/fmt.hpp"

namespace papilo
{

enum kkt_status
{
   OK = 0,
   Fail_Length = 1,
   Fail_Primal_Bound = 2,
   Fail_Primal_Feasibility = 3,
   Fail_Dual_Feasibility = 4,
   Fail_Complementary_Slackness = 5,
   Fail_Stationarity_Lagrangian = 6
};

enum CheckLevel
{
   No_check,
   Check,
   Primal_only,
   Primal_and_dual,
   After_each_step_primal_only,
   After_each_step_and_dual
};

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
   const Problem<REAL>& problem;

   const Vec<REAL>& colValues;
   const Vec<REAL>& colDuals;
   const Vec<REAL>& rowDuals;
   Vec<REAL> rowValues;

   const Vec<uint8_t>& solSetCol;
   const Vec<uint8_t>& solSetRow;

   // references to problem, for easier checking
   const SparseStorage<REAL>& matrixRW;
   const Objective<REAL>& objective;
   const Vec<REAL>& rowLower;
   const Vec<REAL>& rowUpper;
   const Vec<REAL>& colLower;
   const Vec<REAL>& colUpper;

   // zero tolerance todo set on init
   REAL tol = 10e-8;

 public:
   KktState( CheckLevel checker_level, const Problem<REAL>& prob,
             const Solution<REAL>& solution, const Vec<uint8_t>& solSetColumns,
             const Vec<uint8_t>& solSetRows )
       : problem( prob ), colValues( solution.primal ),
         colDuals( solution.col_dual ), rowDuals( solution.row_dual ),
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
      assert( solSetRow.size() == rowValues.size() || 0 == rowValues.size() );
   }

   KktState( CheckLevel checker_level, const Problem<REAL>& prob,
             const Solution<REAL>& solution, const Vec<uint8_t>& solSetColumns,
             const Vec<uint8_t>& solSetRows, const Vec<REAL>& rowLower_,
             const Vec<REAL>& rowUpper_, const Vec<REAL>& colLower_,
             const Vec<REAL>& colUpper_ )
       : problem( prob ), colValues( solution.primal ),
         colDuals( solution.col_dual ), rowDuals( solution.row_dual ),
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
      assert( solSetRow.size() == rowValues.size() || 0 == rowValues.size() );
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
         if( !colLBInf( problem, index ) &&
             colLower[col] - colValues[col] > tol )
         {
            message.info( "Column {:<3} violates lower column bound.\n", col );
            REAL value = colLower[col] - colValues[col];
            info.primal_col_bounds.values.push_back( value );
            status = kkt_status::Fail_Primal_Bound;
         }

         if( !colUBInf( problem, index ) &&
             colValues[col] - colUpper[col] > tol )
         {
            message.info( "Column {:<3} violates upper column bound.\n", col );
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
      const Vec<REAL>& colLower = problem.getVariableDomains().lower_bounds;
      const Vec<REAL>& colUpper = problem.getVariableDomains().upper_bounds;

      // check values of z_j are dual feasible
      int reduced_index = 0;
      for( int col = 0; col < solSetCol.size(); col++ )
      {
         if( !solSetCol[col] )
            continue;

         // no lower or upper bound on column
         if( colLBInf( problem, reduced_index ) &&
             colUBInf( problem, reduced_index ) )
         {
            if( abs( colDuals[col] ) > tol )
               return kkt_status::Fail_Dual_Feasibility;
         }
         // column at lower bound: x=l and l<u
         else if( abs( colValues[col] - colLower[col] ) < tol &&
                  colLower[col] < colUpper[col] )
         {
            if( colDuals[col] < 0 && abs( colDuals[col] ) > tol )
               return kkt_status::Fail_Dual_Feasibility;
         }
         // column at upper bound: x=u and l<u
         else if( abs( colValues[col] - colUpper[col] ) < tol &&
                  colLower[col] < colUpper[col] )
         {
            if( colDuals[col] > tol )
               return kkt_status::Fail_Dual_Feasibility;
         }
         reduced_index++;
      }

      // check values of y_i are dual feasible
      const Vec<REAL>& rowLower =
          problem.getConstraintMatrix().getLeftHandSides();
      const Vec<REAL>& rowUpper =
          problem.getConstraintMatrix().getRightHandSides();

      std::vector<int> orig_col_index( matrixRW.getNCols(), 0 );
      int reduced_col_index = 0;
      for( int col = 0; col < solSetCol.size(); col++ )
      {
         if( !solSetCol[col] )
            continue;
         orig_col_index[reduced_col_index] = col;
         reduced_col_index++;
      }

      assert( matrixRW.getNCols() == reduced_col_index );

      for( int row = 0; row < solSetRow.size(); row++ )
      {
         if( !solSetRow[row] )
            continue;

         // L = Ax = U can be any sign
         if( abs( rowLower[row] - rowValues[row] ) < tol &&
             abs( rowUpper[row] - rowValues[row] ) < tol )
            continue;
         // L = Ax < U
         else if( abs( rowLower[row] - rowValues[row] ) < tol &&
                  rowValues[row] < rowUpper[row] )
         {
            if( rowDuals[row] < -tol )
               return kkt_status::Fail_Dual_Feasibility;
         }
         // L < Ax = U
         else if( rowLower[row] < rowValues[row] &&
                  abs( rowValues[row] - rowUpper[row] ) < tol )
         {
            if( rowDuals[row] > tol )
               return kkt_status::Fail_Dual_Feasibility;
         }
         // L < Ax < U
         else if( rowLower[row] < ( rowValues[row] + tol ) &&
                  rowValues[row] < ( rowUpper[row] + tol ) )
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

         if( !colLBInf( problem, reduced_index ) )
            if( abs( ( colValues[col] - colLower[col] ) * ( colDuals[col] ) ) >
                    tol &&
                colValues[col] != colUpper[col] && abs( colDuals[col] ) > tol )
               return kkt_status::Fail_Complementary_Slackness;

         if( !colUBInf( problem, reduced_index ) )
            if( abs( ( colUpper[col] - colValues[col] ) * ( colDuals[col] ) ) >
                    tol &&
                colValues[col] != colLower[col] && abs( colDuals[col] ) > tol )
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

      // getMatrixTranspose() doesn't work because sometimes the matrix of the
      // reduced problem has not been transposed previously. Edit when KKT
      // check is done at each step, rather than just the start and the end of
      // postsolve.
      const SparseStorage<REAL>& matrixCW =
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

         lagrV = lagrV + colDuals[col] - objective.coefficients[reduced_index];

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
         IndexRange index_range;
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

   Message message;
   CheckLevel level;
};

enum class ProblemType
{
   kOriginal,
   kReduced,
   kPostsolved
};

template <typename REAL, CheckLevel CHECK_LEVEL>
class KktChecker;

template <typename REAL>
class KktChecker<REAL, CheckLevel::No_check>
{
 public:
   CheckLevel level = CheckLevel::No_check;
   using State = int;

   State
   initState( ProblemType type, Solution<REAL>& solution,
              CheckLevel checker_level )
   {
      return 0;
   }

   State
   initState( ProblemType type, Solution<REAL>& solution,
              const Vec<uint8_t>& solSetColumns, const Vec<uint8_t>& solSetRows,
              CheckLevel checker_level )
   {
      return 0;
   }

   State
   initState( ProblemType type, Solution<REAL>& solution,
              const Vec<int>& origcol_map, const Vec<int>& origrow_map,
              CheckLevel checker_level )
   {
      return 0;
   }

   KktChecker() {}

   KktChecker( Problem<REAL> prob ) {}

   void
   checkSolution( State&, bool intermediate = false ) const
   {
   }

   bool
   expandProblem()
   {
   }

   kkt_status
   checkKKT( State& ) const
   {
      return kkt_status::OK;
   }

   void
   setOriginalProblem( const Problem<REAL>& ) const
   {
   }

   void
   setReducedProblem( const Problem<REAL>& ) const
   {
   }

   void
   addRowToProblem( const int row, const int length, const REAL* values,
                    const int* coeffs, const REAL lhs, const REAL rhs,
                    const bool lb_inf, const bool ub_inf )
   {
   }

   void
   undoSubstitutedCol( const int col )
   {
   }

   void
   undoParallelCol( const int col )
   {
   }

   void
   undoFixedCol( const int col, const REAL value )
   {
   }

   void
   undoSingletonRow( const int row )
   {
   }

   void
   undoDeletedCol( const int col )
   {
   }

   void
   undoRedundantRow( const int row )
   {
   }
};

/// class for checking the optimality conditions of the problem during and
/// at end of postsolve. Only contains the original and reduced
/// problem. Any other data required for the check that comes from postsolve is
/// contained in the KktState class
template <typename REAL>
class KktChecker<REAL, CheckLevel::Check>
{
 private:
   Vec<REAL> rowLower_reduced;
   Vec<REAL> rowUpper_reduced;
   Vec<REAL> colLower_reduced;
   Vec<REAL> colUpper_reduced;
   Num<REAL> num;

   // set to false if kkt do not hold
   bool kktHold;

   Problem<REAL> original_problem;
   Problem<REAL> reduced_problem;
   Problem<REAL> problem;

   Vec<uint8_t> solSetRow;
   Vec<uint8_t> solSetCol;

 public:
   CheckLevel level = CheckLevel::Primal_only;
   using State = KktState<REAL>;

   State
   initState( ProblemType type, Solution<REAL>& solution,
              CheckLevel checker_level )
   {
      if( type == ProblemType::kPostsolved )
         return initState( type, solution, solSetCol, solSetRow,
                           checker_level );

      assert( type == ProblemType::kOriginal );
      int ncols = original_problem.getNCols();
      int nrows = original_problem.getNRows();
      solSetCol.resize( ncols );
      solSetCol.assign( ncols, 1 );
      solSetRow.resize( nrows );
      solSetRow.assign( nrows, 1 );

      return initState( type, solution, solSetCol, solSetRow, checker_level );
   }

   State
   initState( ProblemType type, Solution<REAL>& solution,
              const Vec<int>& origcol_map, const Vec<int>& origrow_map,
              CheckLevel checker_level )
   {
      int ncols = original_problem.getNCols();
      int nrows = original_problem.getNRows();
      solSetCol.resize( ncols );
      solSetCol.assign( ncols, 0 );
      solSetRow.resize( nrows );
      solSetRow.assign( nrows, 0 );

      for( int k = 0; k < origcol_map.size(); ++k )
      {
         int origcol = origcol_map[k];
         solSetCol[origcol] = true;
      }

      for( int k = 0; k < origrow_map.size(); ++k )
      {
         int origrow = origrow_map[k];
         solSetRow[origrow] = true;
      }

      return initState( type, solution, solSetCol, solSetRow, checker_level );
   }

   State
   initState( ProblemType type, Solution<REAL>& solution,
              const Vec<uint8_t>& solSetColumns, const Vec<uint8_t>& solSetRows,
              CheckLevel checker_level )
   {
      level = checker_level;

      // todo: make sure constraint matrix transpose is valid since
      //  ...getMatrixTranspose() is used for checking KKT in KktState.

      switch( type )
      {
      case ProblemType::kPostsolved:
      {
         compareMatrixToTranspose( problem.getConstraintMatrix(), num );
         message.info( "\nInitializing check of postsolved solution\n" );
         return State( level, problem, solution, solSetColumns, solSetRows );
      }
      case ProblemType::kOriginal:
      {
         // compares transposes too so no need to call other checks.
         if( level == CheckLevel::After_each_step_and_dual ||
             level == CheckLevel::After_each_step_primal_only )
         {
            bool problems_are_same =
                compareProblems( original_problem, problem, num );
            if( !problems_are_same )
               message.info( "Postsolved problem differs from original.\n" );
         }
         message.info( "Initializing check of original solution\n" );

         if( solSetCol.size() == 0 )
         {
            solSetCol.resize( original_problem.getNCols() );
            solSetCol.assign( original_problem.getNCols(), true );
         }
         assert( solSetCol.size() == original_problem.getNCols() );
         assert( std::all_of( solSetCol.begin(), solSetCol.end(),
                              []( uint8_t isset ) { return isset; } ) );

         if( solSetRow.size() == 0 )
         {
            solSetRow.resize( original_problem.getNRows() );
            solSetRow.assign( original_problem.getNRows(), true );
         }
         assert( solSetRow.size() == original_problem.getNRows() );
         assert( std::all_of( solSetRow.begin(), solSetRow.end(),
                              []( uint8_t isset ) { return isset; } ) );

         // todo: for LP set CheckLevel to prima and dual, once dual postsolve
         // is implemented.
         return State( CheckLevel::Primal_only, original_problem, solution,
                       solSetColumns, solSetRows );
      }
      case ProblemType::kReduced:
      {
         compareMatrixToTranspose( reduced_problem.getConstraintMatrix(), num );

         // if problem type is kReduced expand row / column bound vectors since
         // original_solution is already padded.
         const int nRows = reduced_problem.getNRows();
         const int nCols = reduced_problem.getNCols();

         if( solSetRows.size() > nRows )
         {
            const Vec<REAL> tmp_upper =
                reduced_problem.getConstraintMatrix().getRightHandSides();
            const Vec<REAL> tmp_lower =
                reduced_problem.getConstraintMatrix().getLeftHandSides();

            assert( tmp_upper.size() == nRows );
            assert( tmp_lower.size() == nRows );

            rowUpper_reduced.clear();
            rowLower_reduced.clear();
            rowUpper_reduced.resize( solSetRows.size(), 0 );
            rowLower_reduced.resize( solSetRows.size(), 0 );

            int index = 0;
            for( int k = 0; k < solSetRows.size(); ++k )
            {
               if( solSetRows[k] )
               {
                  rowUpper_reduced[k] = tmp_upper[index];
                  rowLower_reduced[k] = tmp_lower[index];
                  index++;
               }
            }
         }
         else
         {
            assert( solSetRows.size() == nRows );
            rowUpper_reduced =
                reduced_problem.getConstraintMatrix().getRightHandSides();
            rowLower_reduced =
                reduced_problem.getConstraintMatrix().getLeftHandSides();
         }

         if( solSetColumns.size() > nCols )
         {
            Vec<REAL> tmp_upper =
                reduced_problem.getVariableDomains().upper_bounds;
            Vec<REAL> tmp_lower =
                reduced_problem.getVariableDomains().lower_bounds;

            assert( tmp_upper.size() == nCols );
            assert( tmp_lower.size() == nCols );

            colUpper_reduced.clear();
            colLower_reduced.clear();
            colUpper_reduced.resize( solSetColumns.size(), 0 );
            colLower_reduced.resize( solSetColumns.size(), 0 );

            int index = 0;
            for( int k = 0; k < solSetColumns.size(); ++k )
            {
               if( solSetColumns[k] )
               {
                  colUpper_reduced[k] = tmp_upper[index];
                  colLower_reduced[k] = tmp_lower[index];
                  index++;
               }
            }
         }
         else
         {
            assert( solSetColumns.size() == nCols );
            colUpper_reduced =
                reduced_problem.getVariableDomains().upper_bounds;
            colLower_reduced =
                reduced_problem.getVariableDomains().lower_bounds;
         }

         message.info( "\nInitializing check of reduced solution\n" );
         return State( level, reduced_problem, solution, solSetColumns,
                       solSetRows, rowLower_reduced, rowUpper_reduced,
                       colLower_reduced, colUpper_reduced );
      }
      }
   }

   void
   checkSolution( State& state, bool intermediate = false ) const
   {
      state.getRowValues();
      kkt_status status = kkt_status::OK;
      if( state.level == CheckLevel::Primal_only ||
          state.level == CheckLevel::After_each_step_primal_only )
         status = state.checkPrimalFeasibility();
      else
         status = checkKKT( state );

      if( status )
      {
         if( !intermediate )
            message.info( "Check solution: FAIL, status: " );
         else
            message.info( "Check intermediate solution: FAIL, status: " );

         message.info( std::to_string( status ) );
         message.info( "\n" );
      }
      else
      {
         message.info( "Check solution: OK\n" );
      }
      state.updateInfo();
      state.reportKktInfo();
   }

   kkt_status
   checkKKT( State& state ) const
   {
      kkt_status status;
      kkt_status return_status = kkt_status::OK;
      status = state.checkLength();
      if( status != kkt_status::OK )
      {
         message.info( "Solution vector length check failed.\n" );
         return status;
      }

      state.getRowValues();

      status = state.checkPrimalFeasibility();
      if( status != kkt_status::OK )
      {
         return_status = status;
         message.info( "Primal feasibility failed.\n" );
      }

      status = state.checkDualFeasibility();
      if( status != kkt_status::OK )
      {
         return_status = status;
         message.info( "Dual feasibility check failed.\n" );
      }

      status = state.checkComplementarySlackness();
      if( status != kkt_status::OK )
      {
         return_status = status;
         message.info( "Complementary slackness check failed.\n" );
      }

      status = state.checkStOfLagrangian();
      if( status != kkt_status::OK )
      {
         return_status = status;
         message.info( "Stationarity of Lagrangian check failed.\n" );
      }

      return return_status;
   }

   // Expand problem from reduced according to solSet* vectors.
   bool
   expandProblem()
   {
      assert( reduced_problem.getNCols() > -1 );
      assert( original_problem.getNCols() > -1 );
      problem = reduced_problem;

      int ncols = original_problem.getNCols();
      int nrows = original_problem.getNRows();

      ProblemBuilder<REAL> builder;
      builder.setNumCols( ncols );
      builder.setNumRows( nrows );

      std::vector<int> original_index_col;
      original_index_col.reserve( reduced_problem.getNCols() );
      for( int col = 0; col < solSetCol.size(); col++ )
         if( solSetCol[col] )
            original_index_col.push_back( col );
      assert( original_index_col.size() == reduced_problem.getNCols() );

      std::vector<int> original_index_row;
      original_index_row.reserve( reduced_problem.getNRows() );
      for( int row = 0; row < solSetRow.size(); row++ )
         if( solSetRow[row] )
            original_index_row.push_back( row );
      assert( original_index_row.size() == reduced_problem.getNRows() );

      // matrix
      int reduced = 0;
      for( int row = 0; row < solSetRow.size(); row++ )
      {
         const Vec<int>& lengths_o =
             original_problem.getConstraintMatrix().getRowSizes();
         const Vec<int>& lengths_r =
             reduced_problem.getConstraintMatrix().getRowSizes();
         int len_o = lengths_o[row];
         if( solSetRow[row] )
         {
            assert( reduced < lengths_r.size() );
            int len_r = lengths_r[reduced];
            assert( original_index_row[reduced] == row );
            // Saved row at end of presolve.
            const SparseVectorView<REAL>& row_coeff =
                reduced_problem.getConstraintMatrix().getRowCoefficients(
                    reduced );
            assert( row_coeff.getLength() == len_r );
            assert( row_coeff.getLength() > 0 );
            std::vector<int> columns( len_r );
            const int* columns_reduced = row_coeff.getIndices();
            for( int i = 0; i < len_r; i++ )
            {
               assert( len_r = row_coeff.getLength() );
               int col = columns_reduced[i];
               assert( col < original_index_col.size() && col >= 0 );
               columns[i] = original_index_col[col];
            }
            builder.addRowEntries( row, len_r, &columns[0],
                                   row_coeff.getValues() );
            reduced++;
         }
         else
         {
            // todo: // saved row at elimination.
            // empty with size of original.
            int len = len_o;
            std::vector<int> zeros( len, 0 );
            std::vector<REAL> orig_ones( len, 1 );
            const int* cols = &zeros[0];
            const REAL* vals = &orig_ones[0];
            builder.addRowEntries( row, len, cols, vals );
         }
      }
      assert( reduced == reduced_problem.getNRows() );

      // domains and objective
      builder.setObjOffset( reduced_problem.getObjective().offset );
      reduced = 0;
      const Vec<REAL>& reduced_objective =
          reduced_problem.getObjective().coefficients;
      const Vec<REAL>& reduced_lb = reduced_problem.getLowerBounds();
      const Vec<REAL>& reduced_ub = reduced_problem.getUpperBounds();
      const Vec<ColFlags>& reduced_col_flags = reduced_problem.getColFlags();

      for( int col = 0; col < solSetCol.size(); col++ )
      {
         if( solSetCol[col] )
         {
            // add information from reduced problem
            assert( reduced < reduced_col_flags.size() );
            assert( reduced < reduced_lb.size() );
            assert( reduced < reduced_ub.size() );
            assert( reduced < reduced_objective.size() );
            // lb
            if( reduced_col_flags[reduced].test( ColFlag::kLbInf ) )
               builder.setColLbInf( col, true );
            else
               builder.setColLb( col, reduced_lb[reduced] );
            // ub
            if( reduced_col_flags[reduced].test( ColFlag::kUbInf ) )
               builder.setColUbInf( col, true );
            else
               builder.setColUb( col, reduced_ub[reduced] );

            builder.setObj( col, reduced_objective[reduced] );
            reduced++;
         }
         else
         {
            // add infinite bounds and zero cost
            builder.setColLbInf( col, true );
            builder.setColUbInf( col, true );
            builder.setObj( col, 0 );
         }
      }
      assert( reduced == reduced_problem.getNCols() );

      // row bounds
      const Vec<RowFlags>& reduced_flags = reduced_problem.getRowFlags();
      const Vec<REAL>& reduced_lhs =
          reduced_problem.getConstraintMatrix().getLeftHandSides();
      const Vec<REAL>& reduced_rhs =
          reduced_problem.getConstraintMatrix().getRightHandSides();
      reduced = 0;
      for( int row = 0; row < solSetRow.size(); row++ )
      {
         if( solSetRow[row] )
         {
            assert( original_index_row[reduced] == row );
            // lhs
            if( reduced_flags[reduced].test( RowFlag::kLhsInf ) )
               builder.setRowLhsInf( row, true );
            else
               builder.setRowLhs( row, reduced_lhs[reduced] );
            // rhs
            if( reduced_flags[reduced].test( RowFlag::kRhsInf ) )
               builder.setRowRhsInf( row, true );
            else
               builder.setRowRhs( row, reduced_rhs[reduced] );
            reduced++;
         }
         else
         {
            // todo: add bounds of row at elimination?
            builder.setRowLhsInf( row, true );
            builder.setRowRhsInf( row, true );
         }
      }
      assert( reduced == reduced_problem.getNRows() );

      problem = builder.build();

      problem.setVariableNames( original_problem.getVariableNames() );
      problem.setConstraintNames( original_problem.getConstraintNames() );
      problem.setName( original_problem.getName() );

      return true;
   }

   void
   setOriginalProblem( const Problem<REAL>& prob )
   {
      original_problem = prob;
   };
   void
   setReducedProblem( const Problem<REAL>& problem )
   {
      reduced_problem = problem;
   }

   void
   addRowToProblem( const int row, const int length, const REAL* values,
                    const int* coeffs, const REAL lhs, const REAL rhs,
                    const bool lb_inf, const bool ub_inf )
   {
      addRow( problem, row, length, values, coeffs, lhs, rhs, lb_inf, ub_inf );
   }

   void
   undoSubstitutedCol( const int col )
   {
      assert( !solSetCol[col] );
      solSetCol[col] = true;
      // todo: ...
   }

   void
   undoParallelCol( const int col )
   {
      assert( !solSetCol[col] );
      solSetCol[col] = true;
      // todo: ...
   }

   void
   undoFixedCol( const int col, const REAL value )
   {
      assert( !solSetCol[col] );
      solSetCol[col] = true;
      // todo:
      // addCol( problem, col, length, values, coeffs, l, u, lb_inf, ub_inf );
   }

   void
   undoSingletonRow( const int row )
   {
      assert( solSetRow[row] );
      // todo:
      // addRow( problem, row, length, values, coeffs, lhs, rhs, lhs_inf,
      // rhs_inf );
   }

   void
   undoDeletedCol( const int col )
   {
      assert( !solSetCol[col] );
      solSetCol[col] = true;
      // todo:
      // addCol( problem, col, length, values, coeffs, l, u, lb_inf, ub_inf );
   }

   void
   undoRedundantRow( const int row )
   {
      assert( !solSetRow[row] );
      solSetRow[row] = true;
   }

   Message message;
};

template <typename REAL>
bool
rowLHSInf( const Problem<REAL>& problem, const int row )
{
   return problem.getRowFlags()[row].test( RowFlag::kLhsInf );
}

template <typename REAL>
bool
rowRHSInf( const Problem<REAL>& problem, const int row )
{
   return problem.getRowFlags()[row].test( RowFlag::kRhsInf );
}

template <typename REAL>
bool
colLBInf( const Problem<REAL>& problem, const int col )
{
   return problem.getColFlags()[col].test( ColFlag::kLbInf );
}

template <typename REAL>
bool
colUBInf( const Problem<REAL>& problem, const int col )
{
   return problem.getColFlags()[col].test( ColFlag::kUbInf );
}

// for testing
template <typename REAL>
bool
check_solution_feasibility( const Problem<REAL>& prob, Solution<REAL>& sol )
{
   KktChecker<REAL, CheckLevel::Check> chk( prob );

   Vec<uint8_t> colsSet( prob.getNCols(), true );
   Vec<uint8_t> rowsSet( prob.getNRows(), true );

   auto state = chk.initState( true, sol, colsSet, rowsSet, chk.level );
   kkt_status status = state.checkSolution( prob );
   if( status == kkt_status::OK )
      return true;

   return false;
}

} // namespace papilo

#endif