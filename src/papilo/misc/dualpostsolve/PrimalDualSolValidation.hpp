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

#ifndef _PAPILO_CORE_PRIMAL_DUAL_SOL_VALIDATION_HPP_
#define _PAPILO_CORE_PRIMAL_DUAL_SOL_VALIDATION_HPP_

#include "papilo/core/Solution.hpp"
#include "papilo/core/postsolve/PostsolveStatus.hpp"

namespace papilo
{
template <typename REAL>
class PrimalDualSolValidation
{

 private:
   const Num<REAL> num;
   Message message{};

   bool
   checkLength( const Solution<REAL>& solution, const Problem<REAL>& problem )
   {
      const int nCols = problem.getNCols();

      bool primal_check = solution.primal.size() != nCols;
      if( solution.type == SolutionType::kPrimalDual )
         return primal_check || solution.reducedCosts.size() != nCols ||
                solution.dual.size() != problem.getNRows();
      return primal_check;
   }

   bool
   checkPrimalBounds( const Vec<REAL>& primalSolution,
                      const Problem<REAL>& problem )
   {
      bool failure = false;

      const Vec<REAL> ub = problem.getUpperBounds();
      const Vec<REAL> lb = problem.getLowerBounds();

      for( unsigned int col = 0; col < problem.getNCols(); col++ )
      {
         if( problem.getColFlags()[col].test( ColFlag::kInactive ) )
            continue;

         if( ( not problem.getColFlags()[col].test( ColFlag::kLbInf ) ) &&
             num.isLT( primalSolution[col], lb[col] ) )
         {
            message.info( "Column {:<3} violates lower column bound.\n", col );
            failure = true;
         }

         if( ( not problem.getColFlags()[col].test( ColFlag::kUbInf ) ) &&
             num.isGT( primalSolution[col], ub[col] ) )
         {
            message.info( "Column {:<3} violates upper column bound.\n", col );
            failure = true;
         }
      }
      return failure;
   }

   bool
   checkPrimalConstraint( const Vec<REAL>& primalSolution,
                          const Problem<REAL>& problem ) const
   {
      const Vec<REAL> rhs = problem.getConstraintMatrix().getRightHandSides();
      const Vec<REAL> lhs = problem.getConstraintMatrix().getLeftHandSides();
      for( int row = 0; row < problem.getNRows(); row++ )
      {
         if( problem.getRowFlags()[row].test( RowFlag::kRedundant ) )
            continue;

         REAL rowValue = 0;
         auto entries = problem.getConstraintMatrix().getRowCoefficients( row );
         for( int j = 0; j < entries.getLength(); j++ )
         {
            int col = entries.getIndices()[j];
            if( problem.getColFlags()[col].test( ColFlag::kFixed ) )
               continue;
            rowValue += entries.getValues()[j] * primalSolution[col];
         }

         if( ( not problem.getRowFlags()[row].test( RowFlag::kLhsInf ) ) &&
             num.isLT( rowValue, lhs[row] ) )
         {
            message.info( "Row {:<3} violates row bounds ({:<3} < {:<3}).\n",
                          row, lhs[row], rowValue );
            return true;
         }
         if( ( not problem.getRowFlags()[row].test( RowFlag::kRhsInf ) ) &&
             num.isGT( rowValue, rhs[row] ) )
         {
            message.info( "Row {:<3} violates row bounds ({:<3} < {:<3}).\n",
                          row, rowValue, rhs[row] );
            return true;
         }
      }
      return false;
   }

   bool
   checkPrimalFeasibility( const Vec<REAL>& primalSolution,
                           const Problem<REAL>& problem )
   {
      return checkPrimalBounds( primalSolution, problem ) and
             checkPrimalConstraint( primalSolution, problem );
   }

   bool
   checkDualFeasibility( const Vec<REAL>& primalSolution,
                         const Vec<REAL>& dualSolution,
                         const Vec<REAL>& reducedCosts,
                         const Problem<REAL>& problem )
   {
      const papilo::Vec<REAL>& lowerBounds = problem.getLowerBounds();
      const papilo::Vec<REAL>& upperBounds = problem.getUpperBounds();

      const Vec<REAL> rhs = problem.getConstraintMatrix().getRightHandSides();
      const Vec<REAL> lhs = problem.getConstraintMatrix().getLeftHandSides();

      for( int variable = 0; variable < problem.getNCols(); variable++ )
      {
         if( problem.getColFlags()[variable].test( ColFlag::kFixed ) )
            continue;

         REAL colValue = 0;

         auto coeff =
             problem.getConstraintMatrix().getColumnCoefficients( variable );
         for( int counter = 0; counter < coeff.getLength(); counter++ )
         {
            REAL value = coeff.getValues()[counter];
            int rowIndex = coeff.getIndices()[counter];
            colValue += dualSolution[rowIndex] * value;
         }

         if( not num.isEq( colValue + reducedCosts[variable],
                           problem.getObjective().coefficients[variable] ) )
            return true;
      }
      return false;

      //      // check values of z_j are dual feasible
      //      for( int col = 0; col < problem.getNCols(); col++ )
      //      {
      //         if( problem.getColFlags()[col].test( ColFlag::kInactive ) )
      //            continue;
      //
      //         // no lower and upper bound on infinity
      //         if( problem.getColFlags()[col].test( ColFlag::kLbInf ) &&
      //             problem.getColFlags()[col].test( ColFlag::kUbInf ) )
      //         {
      //            if( not num.isZero( reducedCosts[col] ) )
      //               return true;
      //         }
      //         // non fixed variable at lower bound: x=l and l<u
      //         else if( num.isEq( primalSolution[col], lowerBounds[col] ) &&
      //                  num.isLT( lowerBounds[col], upperBounds[col] ) )
      //         {
      //            if( num.isLT( reducedCosts[col], 0 ) )
      //               return true;
      //         }
      //         // non fixed variable at upper bound: x=u and l<u
      //         else if( num.isEq( primalSolution[col], upperBounds[col] ) &&
      //                  num.isLT( lowerBounds[col], upperBounds[col] ) )
      //         {
      //            if( num.isGT( reducedCosts[col], 0 ) )
      //               return true;
      //         }
      //      }
      //
      //      // check values of y_i are dual feasible
      //      const papilo::Vec<REAL>& lhs =
      //          problem.getConstraintMatrix().getLeftHandSides();
      //      const papilo::Vec<REAL>& rhs =
      //          problem.getConstraintMatrix().getRightHandSides();
      //
      //      for( int row = 0; row < problem.getNRows(); row++ )
      //      {
      //         if( problem.getRowFlags()[row].test( RowFlag::kRedundant ) )
      //            continue;
      //
      //         REAL rowValue = 0;
      //         auto entries =
      //         problem.getConstraintMatrix().getRowCoefficients( row ); for(
      //         int j = 0; j < entries.getLength(); j++ )
      //         {
      //            int col = entries.getIndices()[j];
      //            if( problem.getColFlags()[col].test( ColFlag::kFixed ) )
      //               continue;
      //            rowValue += entries.getValues()[j] * primalSolution[col];
      //         }
      //
      //         bool isLhsInf = problem.getRowFlags()[row].test(
      //         RowFlag::kLhsInf ); bool isRhsInf =
      //         problem.getRowFlags()[row].test( RowFlag::kRhsInf ); assert(
      //         not( isLhsInf and isRhsInf ) );
      //
      //         // L = Ax = U can be any sign
      //         if( problem.getRowFlags()[row].test( RowFlag::kEquation ) )
      //         {
      //            assert( num.isEq( lhs[row], rowValue ) and
      //                    num.isEq( rhs[row], rowValue ) );
      //            continue;
      //         }
      //         else if( isLhsInf )
      //         {
      //            assert( not isRhsInf );
      //            if( num.isLT( rowValue, rhs[row] ) and  not num.isZero(
      //            dualSolution[row] ) )
      //                  return true;
      //
      //         }
      //         else if( isRhsInf )
      //         {
      //            if( num.isGT( rowValue, lhs[row] ) and  not num.isZero(
      //            dualSolution[row] ) )
      //                  return true;
      //         }
      //         else
      //         {
      //            if( num.isGT( rowValue, lhs[row] ) and
      //                num.isEq( rowValue, rhs[row] ) )
      //            {
      //               if( not num.isZero( dualSolution[row] ) )
      //                  return true;
      //            }
      //            else if( num.isEq( rowValue, lhs[row] ) and
      //                     num.isLT( rowValue, rhs[row] ) )
      //            {
      //               if( not num.isZero( dualSolution[row] ) )
      //                  return true;
      //            }
      //            else
      //            {
      //               assert( num.isGT( rowValue, lhs[row] ) and
      //                       num.isLT( rowValue, rhs[row] ) );
      //               if( num.isZero( dualSolution[row] ) )
      //                  return true;
      //            }
      //         }
      //      }
   }

   bool
   checkObjectiveFunction( const Vec<REAL>& primalSolution,
                           const Vec<REAL>& dualSolution,
                           const Vec<REAL>& reducedCosts,
                           const Problem<REAL>& problem )
   {
      REAL primalObjectiveValue = 0;
      for( int i = 0; i < problem.getNCols(); i++ )
      {
         primalObjectiveValue +=
             primalSolution[i] * problem.getObjective().coefficients[i];
      }
      REAL dualObjectiveValue = 0;
      for( int i = 0; i < problem.getNRows(); i++ )
      {
         REAL dual = dualSolution[i];
         if(dual < 0)
            dualObjectiveValue +=
                dual * problem.getConstraintMatrix().getRightHandSides()[i];
         else
            dualObjectiveValue +=
                dual * problem.getConstraintMatrix().getLeftHandSides()[i];
      }
      for( int i = 0; i < problem.getNCols(); i++ )
      {

         REAL reducedCost = reducedCosts[i];
         if(reducedCost < 0)
            dualObjectiveValue +=
                reducedCost * problem.getUpperBounds()[i];
         else
            dualObjectiveValue +=
                reducedCost * problem.getLowerBounds()[i];

      }
      return not num.isEq( primalObjectiveValue, dualObjectiveValue );
   }

   bool
   checkComplementarySlackness( const Vec<REAL>& primalSolution,
                                const Vec<REAL>& dualSolution,
                                const Vec<REAL>& reducedCosts,
                                const Problem<REAL>& problem )

   {

      const Vec<REAL> lb = problem.getLowerBounds();
      const Vec<REAL> ub = problem.getUpperBounds();

      const Vec<REAL> rhs = problem.getConstraintMatrix().getRightHandSides();
      const Vec<REAL> lhs = problem.getConstraintMatrix().getLeftHandSides();
      for( int row = 0; row < problem.getNRows(); row++ )
      {
         if( problem.getRowFlags()[row].test( RowFlag::kRedundant ) )
            continue;

         REAL rowValue = 0;
         auto entries = problem.getConstraintMatrix().getRowCoefficients( row );
         for( int j = 0; j < entries.getLength(); j++ )
         {
            int col = entries.getIndices()[j];
            if( problem.getColFlags()[col].test( ColFlag::kFixed ) )
               continue;
            rowValue += entries.getValues()[j] * primalSolution[col];
         }

         if( not problem.getRowFlags()[row].test( RowFlag::kLhsInf ) and
             not problem.getRowFlags()[row].test( RowFlag::kRhsInf ) )
         {
            if( num.isGT( lhs[row], rowValue ) and
                num.isLT( rhs[row], rowValue ) and
                not num.isZero( dualSolution[row] ) )
               return true;
         }
         else if( not problem.getRowFlags()[row].test( RowFlag::kLhsInf ) )
         {
            assert( problem.getRowFlags()[row].test( RowFlag::kRhsInf ) );
            if( num.isGT( lhs[row], rowValue ) and
                not num.isZero( dualSolution[row] ) )
               return true;
         }

         if( ( not problem.getRowFlags()[row].test( RowFlag::kLhsInf ) ) &&
             num.isGT( rowValue, lhs[row] ) )
         {
            assert( problem.getRowFlags()[row].test( RowFlag::kRhsInf ) );
            if( num.isLT( rhs[row], rowValue ) and
                not num.isZero( dualSolution[row] ) )
               return true;
         }
      }

      for( int col = 0; col < problem.getNCols(); col++ )
      {
         if( problem.getColFlags()[col].test( ColFlag::kInactive ) )
            continue;

         bool isLbInf = problem.getColFlags()[col].test( ColFlag::kLbInf );
         bool isUbInf = problem.getColFlags()[col].test( ColFlag::kUbInf );
         REAL upperBound = ub[col];
         REAL lowerBound = lb[col];
         REAL reducedCost = reducedCosts[col];
         REAL sol = primalSolution[col];

         // TODO: check this
         if( upperBound == lowerBound and not isLbInf and not isUbInf )
            continue;

         if( not isLbInf and not isUbInf )
         {
            if( num.isGT( sol, lowerBound ) and num.isLT( sol, upperBound ) and
                not num.isZero( reducedCost ) )
               return true;
         }
         else if( not isLbInf )
         {
            assert( isUbInf );
            if( num.isGT( sol, lowerBound ) and not num.isZero( reducedCost ) )
               return true;
         }
         else if( not isUbInf )
         {
            assert( isUbInf );
            if( num.isLT( sol, upperBound ) and not num.isZero( reducedCost ) )
               return true;
         }
      }
      return false;
   }

   //   bool
   //   checkStOfLagrangian( const Vec<REAL>& primalSolution,
   //                        const Vec<REAL>& dualSolution,
   //                        const Vec<REAL>& reducedCosts,
   //                        const Problem<REAL>& problem )
   //
   //   {
   //      // A'y + reduced_costs = c
   //
   //      const papilo::SparseStorage<REAL>& transposed =
   //          problem.getConstraintMatrix().getMatrixTranspose();
   //
   //      std::vector<int> orig_row_index( transposed.getNCols(), 0 );
   //
   //      for( int col = 0; col < transposed.getNCols(); col++ )
   //      {
   //         if( problem.getColFlags()[col].test( ColFlag::kFixed ) )
   //            continue;
   //         REAL lagrV = 0;
   //
   //         auto index_range = transposed.getRowRanges()[col];
   //         for( int k = index_range.start; k < index_range.end; k++ )
   //         {
   //            int row = transposed.getColumns()[k];
   //            lagrV += dualSolution[row] * transposed.getValues()[k];
   //         }
   //
   //         if( not num.isEq( lagrV + reducedCosts[col],
   //                           problem.getObjective().coefficients[col] ) )
   //            return true;
   //      }
   //
   //      return false;
   //   }

 public:
   PostsolveStatus
   verifySolution( const Solution<REAL>& solution,
                   const Problem<REAL>& problem )
   {

      bool failure = checkLength( solution, problem );
      if( failure )
      {
         message.info( "Solution vector length check FAILED.\n" );
         return PostsolveStatus::kFailed;
      }

      failure = checkPrimalFeasibility( solution.primal, problem );
      if( failure )
      {
         message.info( "Primal feasibility check FAILED.\n" );
         return PostsolveStatus::kFailed;
      }

      if( solution.type == SolutionType::kPrimalDual )
      {
         failure = checkDualFeasibility( solution.primal, solution.dual,
                                         solution.reducedCosts, problem );
         if( failure )
         {
            message.info( "Dual feasibility check FAILED.\n" );
            return PostsolveStatus::kFailed;
         }

         failure = checkObjectiveFunction( solution.primal, solution.dual,
                                           solution.reducedCosts, problem );
         if( failure )
         {
            message.info( "Objective function failed.\n" );
            return PostsolveStatus::kFailed;
         }

         failure = checkComplementarySlackness(
             solution.primal, solution.dual, solution.reducedCosts, problem );
         if( failure )
         {
            message.info( "Complementary slack check FAILED.\n" );
            return PostsolveStatus::kFailed;
         }

      }

      message.info( "Solution passed validation\n" );
      return PostsolveStatus::kOk;
   }
};
} // namespace papilo

#endif