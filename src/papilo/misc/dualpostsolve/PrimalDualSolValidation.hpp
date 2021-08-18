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
   REAL maximal_allowed_in_duality_gap = 0;




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
   checkPrimalConstraint( Solution<REAL>& solution,
                          const Problem<REAL>& problem, bool update_slack ) const
   {
      const Vec<REAL> rhs = problem.getConstraintMatrix().getRightHandSides();
      const Vec<REAL> lhs = problem.getConstraintMatrix().getLeftHandSides();
      if( update_slack )
      {
         solution.slack.clear();
         solution.slack.resize( problem.getNRows() );
      }
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
            rowValue += entries.getValues()[j] * solution.primal[col];
         }

         bool lhs_inf = problem.getRowFlags()[row].test( RowFlag::kLhsInf );
         if( ( not lhs_inf ) &&
             num.isLT( rowValue, lhs[row] ) )
         {
            message.info( "Row {:<3} violates row bounds ({:<3} < {:<3}).\n",
                          row, lhs[row], rowValue );
            return true;
         }
         bool rhs_inf = problem.getRowFlags()[row].test( RowFlag::kRhsInf );
         if( ( not rhs_inf ) &&
             num.isGT( rowValue, rhs[row] ) )
         {
            message.info( "Row {:<3} violates row bounds ({:<3} < {:<3}).\n",
                          row, rowValue, rhs[row] );
            return true;
         }
         if( update_slack )
            solution.slack[row] = num.isZero( rowValue ) ? 0 : rowValue;
      }
      return false;
   }

   bool
   checkPrimalFeasibility( Solution<REAL>& solution,
                           const Problem<REAL>& problem, bool update_slack )
   {
      return checkPrimalBounds( solution.primal, problem ) or
             checkPrimalConstraint( solution, problem, update_slack );
   }

   bool
   checkDualFeasibility( const Vec<REAL>& primalSolution,
                         const Vec<REAL>& dualSolution,
                         const Vec<REAL>& reducedCosts,
                         const Vec<VarBasisStatus>& basis,
                         const Problem<REAL>& problem )
   {
      const papilo::Vec<REAL>& lowerBounds = problem.getLowerBounds();
      const papilo::Vec<REAL>& upperBounds = problem.getUpperBounds();

      const Vec<REAL> rhs = problem.getConstraintMatrix().getRightHandSides();
      const Vec<REAL> lhs = problem.getConstraintMatrix().getLeftHandSides();

      for( int variable = 0; variable < problem.getNCols(); variable++ )
      {
         if( problem.getColFlags()[variable].test( ColFlag::kInactive ) )
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
         {
            message.info(
                "Dual row {:<3} violates dual row bounds ({:<3} != {:<3}).\n",
                variable, problem.getObjective().coefficients[variable],
                colValue + reducedCosts[variable],
                problem.getObjective().coefficients[variable] );
            return true;
         }

         bool ub_infinity =
             problem.getColFlags()[variable].test( ColFlag::kUbInf );
         bool lb_infinity =
             problem.getColFlags()[variable].test( ColFlag::kLbInf );
         REAL lb = problem.getLowerBounds()[variable];
         REAL ub = problem.getUpperBounds()[variable];
         REAL sol = primalSolution[variable];

         assert( ub_infinity or lb_infinity or num.isFeasGE( ub, lb ) );
         switch( basis[variable] )
         {
         case VarBasisStatus::FIXED:
            if( ub_infinity or lb_infinity or not num.isEq( lb, ub ) or
                not num.isEq( sol, ub ) )
               return true;
            break;
         case VarBasisStatus::ON_LOWER:
            if( lb_infinity or not num.isEq( sol, lb ) or
                ( ub_infinity and num.isEq( sol, ub ) ) or
                ( num.isZero( lb ) and ub_infinity ) )
               return true;
            break;
         case VarBasisStatus::ON_UPPER:
            if( ub_infinity or not num.isEq( sol, ub ) or
                ( lb_infinity and num.isEq( sol, lb ) ) )
               return true;
            break;
         case VarBasisStatus::ZERO:
            if( lb_infinity or not num.isZero( sol ) or
                not num.isZero( lb ) and not ub_infinity )
               return true;
            break;
         case VarBasisStatus::BASIC:
            if( ( not lb_infinity and num.isEq( sol, lb ) ) or
                ( not ub_infinity and num.isEq( sol, ub ) ) )
               return true;
            break;
         case VarBasisStatus::UNDEFINED:
            return true;
         }
      }
      return false;
   }

   bool
   checkObjectiveFunction( const Vec<REAL>& primalSolution,
                           const Vec<REAL>& dualSolution,
                           const Vec<REAL>& reducedCosts,
                           const Problem<REAL>& problem )
   {
      REAL duality_gap =
          getDualityGap( primalSolution, dualSolution, reducedCosts, problem );
      return not( num.isFeasZero( duality_gap ) or
                  num.isLT( duality_gap, maximal_allowed_in_duality_gap ) );
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
         else if( ( not problem.getRowFlags()[row].test( RowFlag::kLhsInf ) ) &&
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
         if( num.isEq(upperBound, lowerBound) and not isLbInf and not isUbInf )
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
            assert( isLbInf );
            if( num.isLT( sol, upperBound ) and not num.isZero( reducedCost ) )
               return true;
         }
      }
      return false;
   }

 public:
   PostsolveStatus
   verifySolution( Solution<REAL>& solution,
                   const Problem<REAL>& problem, bool update_slack )
   {

      bool failure = checkLength( solution, problem );
      if( failure )
      {
         message.info( "Solution vector length check FAILED.\n" );
         return PostsolveStatus::kFailed;
      }

      failure = checkPrimalFeasibility( solution, problem, update_slack );
      if( failure )
      {
         message.info( "Primal feasibility check FAILED.\n" );
         return PostsolveStatus::kFailed;
      }

      if( solution.type == SolutionType::kPrimalDual )
      {
         if( checkDualFeasibility( solution.primal, solution.dual,
                                   solution.reducedCosts,
                                   solution.varBasisStatus, problem ) )
         {
            message.info( "Dual feasibility check FAILED.\n" );
            failure = true;
         }

         if( checkComplementarySlackness( solution.primal, solution.dual,
                                          solution.reducedCosts, problem ) )
         {
            message.info( "Complementary slack check FAILED.\n" );
            failure = true;
         }

         if( checkObjectiveFunction( solution.primal, solution.dual,
                                     solution.reducedCosts, problem ) )
         {
            message.info( "Objective function failed.\n" );
            //            failure = true;
         }
         if( failure )
            return PostsolveStatus::kFailed;
      }

      message.info( "Solution passed validation\n" );
      return PostsolveStatus::kOk;
   }

   void
   setDualityGap( REAL duality_gap )
   {
      maximal_allowed_in_duality_gap = duality_gap * 10;
   }

   REAL
   getDualityGap( const Vec<REAL>& primalSolution,
                  const Vec<REAL>& dualSolution, const Vec<REAL>& reducedCosts,
                  const Problem<REAL>& problem )
                  {
      StableSum<REAL> primal_objective;
      for( int i = 0; i < problem.getNCols(); i++ )
      {
         primal_objective.add( primalSolution[i] *
         problem.getObjective().coefficients[i] );
      }
      StableSum<REAL> dual_objective;
      for( int i = 0; i < problem.getNRows(); i++ )
      {
         REAL dual = dualSolution[i];
         REAL side;
         if( dual < 0 )
            side = problem.getConstraintMatrix().getRightHandSides()[i];
         else
            side = problem.getConstraintMatrix().getLeftHandSides()[i];
         dual_objective.add( dual * side );
      }
      for( int i = 0; i < problem.getNCols(); i++ )
      {
         REAL reducedCost = reducedCosts[i];
         REAL side;
         if( reducedCost < 0 )
            side = problem.getUpperBounds()[i];
         else
            side = problem.getLowerBounds()[i];
         dual_objective.add( reducedCost * side );
      }
      return primal_objective.get() - dual_objective.get();
                  }

};
} // namespace papilo

#endif