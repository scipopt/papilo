
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
   OK,
   Fail_Length,
   Fail_Primal_Bound,
   Fail_Primal_Feasibility,
   Fail_Dual_Feasibility,
   Fail_Complementary_Slackness,
   Fail_Stationarity_Lagrangian
};

enum CheckLevel
{
   No_check,
   Check,
   Primal_feasibility_only,
   Solver_and_primal_feas,   // expects dual values from solver
   Postsolved_problem_full,  // expects dual values from solver
   After_each_postsolve_step // expects dual values from solver
};

/// class to hold all the data needed for the checks. The other class
/// KktChecker only holds copies of the original and reduced problem data
template <typename REAL>
class KktState
{
 private:
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
         if( !colLBInf( problem, reduced_index ) &&
             colLower[col] - colValues[col] > tol )
         {
            message.info( "Column {:<3} violates lower column bound.\n", col );
            status = kkt_status::Fail_Primal_Bound;
         }

         if( !colUBInf( problem, reduced_index ) &&
             colValues[col] - colUpper[col] > tol )
         {
            message.info( "Column {:<3} violates upper column bound.\n", col );
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
          colValues.size() != nCols || colDuals.size() != nCols )
         return kkt_status::Fail_Length;

      const int nRows = solSetRow.size();

      if( rowLower.size() != nRows || rowUpper.size() != nRows ||
          rowValues.size() != nRows || rowDuals.size() != nRows )
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
      if( status != kkt_status::OK )
         return status;

      int reduced_index = 0;
      for( int row = 0; row < solSetRow.size(); row++ )
      {
         if( !solSetRow[row] )
            continue;

         if( !rowLHSInf( problem, reduced_index ) &&
             rowLower[row] - rowValues[row] > tol )
         {
            message.info( "Row {:<3} violates row bounds.\n", row );
            status = kkt_status::Fail_Primal_Feasibility;
         }
         if( !rowRHSInf( problem, reduced_index ) &&
             rowValues[row] - rowUpper[row] > tol )
         {
            message.info( "Row {:<3} violates row bounds.\n", row );
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
      assert( matrixCW.getNCols() == reduced_row_index );

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
      assert( orig_col_index.size() == nCols );

      for( int i = 0; i < solSetRow.size(); i++ )
      {
         if( !solSetRow[i] )
            continue;

         rowValue = 0;
         auto index_range = matrixRW.getRowRanges()[reduced_row_index];
         for( int j = index_range.start; j < index_range.end; j++ )
         {
            int col = matrixRW.getColumns()[j];
            assert( col >= 0 );
            assert( col < (int)nCols );

            // col is index of the column in reduced problem and colValues is
            // expanded.
            rowValue += values[j] * colValues[orig_col_index[col]];
         }
         rowValues[i] = rowValue;
         reduced_row_index++;
      }
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
   checkSolution( State& ) const
   {
   }

   void
   checkIntermediate( State& ) const
   {
   }

   kkt_status
   checkKKT( State& ) const
   {
      return kkt_status::OK;
   }

   void
   setLevel( State& state, CheckLevel type_of_check ) const
   {
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

 public:
   CheckLevel level = CheckLevel::Primal_feasibility_only;
   using State = KktState<REAL>;

   State
   initState( ProblemType type, Solution<REAL>& solution,
              CheckLevel checker_level )
   {
      int ncols = original_problem.getNCols();
      int nrows = original_problem.getNRows();
      Vec<uint8_t> solSetColumns( ncols, 1 );
      Vec<uint8_t> solSetRows( nrows, 1 );

      return initState( type, solution, solSetColumns, solSetRows,
                        checker_level );
   }

   State
   initState( ProblemType type, Solution<REAL>& solution,
              const Vec<int>& origcol_map, const Vec<int>& origrow_map,
              CheckLevel checker_level )
   {
      int ncols = original_problem.getNCols();
      int nrows = original_problem.getNRows();
      Vec<uint8_t> solSetColumns( ncols, 0 );
      Vec<uint8_t> solSetRows( nrows, 0 );

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

      if( type == ProblemType::kPostsolved )
      {
         compareMatrixToTranspose( problem.getConstraintMatrix(), num );
         message.info( "Initializing check of postsolved solution\n" );
         return State( level, problem, solution, solSetColumns, solSetRows );
      }

      if( type == ProblemType::kOriginal )
      {
         // compares transposes too so no need to call other checks.
         if( level == CheckLevel::After_each_postsolve_step )
         {
            bool problems_are_same =
                compareProblems( original_problem, problem, num );
            assert( problems_are_same );
         }
         message.info( "Initializing check of original solution\n" );
         // todo: assert solSetRows and SolSetColumns are all set.
         return State( level, original_problem, solution, solSetColumns,
                       solSetRows );
      }

      // matrix type is ProblemType::kReduced.
      compareMatrixToTranspose( reduced_problem.getConstraintMatrix(), num );

      // if problem type is REDUCED expand row / column bound vectors since
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
         colUpper_reduced = reduced_problem.getVariableDomains().upper_bounds;
         colLower_reduced = reduced_problem.getVariableDomains().lower_bounds;
      }

      message.info( "Initializing check of reduced solution\n" );
      return State( level, reduced_problem, solution, solSetColumns, solSetRows,
                    rowLower_reduced, rowUpper_reduced, colLower_reduced,
                    colUpper_reduced );
   }

   void
   setLevel( CheckLevel type_of_check )
   {
      level = type_of_check;
   }

   void
   checkSolution( State& state ) const
   {
      state.getRowValues();
      kkt_status status = kkt_status::OK;
      if( state.level == CheckLevel::Primal_feasibility_only )
         status = state.checkPrimalFeasibility();
      else
         status = checkKKT( state );

      if( status )
      {
         message.info( "Check solution: FAIL, status: " );
         message.info( std::to_string( status ) );
         message.info( "\n" );
      }
      else
      {
         message.info( "Check solution: OK\n" );
      }
   }

   void
   checkIntermediate( State& state ) const
   {
      if( state.level != CheckLevel::After_each_postsolve_step )
         return;

      kkt_status status = checkKKT( state );
      if( status )
      {
         message.info( "KKT intermediate FAIL, status: " );
         message.info( std::to_string( status ) );
         message.info( "\n" );
      }
      return;
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

   void
   setOriginalProblem( const Problem<REAL>& prob )
   {
      original_problem = prob;
   };

   void
   setReducedProblem( const Problem<REAL>& problem )
   {
      reduced_problem = problem;

      // todo:
      // assert original_problem is set

      // allocate memory for reduced. for now for each row / col take
      // max of row/col in original and reduced problem.
   }

   void
   addRowToProblem( const int row, const int length, const REAL* values,
                    const int* coeffs, const REAL lhs, const REAL rhs,
                    const bool lb_inf, const bool ub_inf )
   {
      addRow( problem, row, length, values, coeffs, lhs, rhs, lb_inf, ub_inf );
   }

   Message message;

 private:
   // set to false if kkt do not hold
   bool kktHold;

   Problem<REAL> original_problem;
   Problem<REAL> reduced_problem;
   Problem<REAL> problem;

   Vec<uint8_t> solSetRow;
   Vec<uint8_t> solSetCol;
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