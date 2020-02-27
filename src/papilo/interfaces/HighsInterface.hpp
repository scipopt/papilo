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

#ifndef _PAPILO_INTERFACES_HIGHS_INTERFACE_HPP_
#define _PAPILO_INTERFACES_HIGHS_INTERFACE_HPP_

#include "papilo/misc/String.hpp"
#include "papilo/misc/Vec.hpp"
#include <cassert>
#include <stdexcept>
#include <type_traits>

#include "Highs.h"
#include "papilo/core/Problem.hpp"
#include "papilo/interfaces/SolverInterface.hpp"

namespace papilo
{

template <typename REAL>
class HighsInterface : public SolverInterface<REAL>
{
 private:
   Highs solver;
   HighsOptions opts;
   static constexpr double inf = 1e+25;

 public:
   HighsInterface()
   {
      opts.infinite_bound = inf;
      opts.infinite_cost = inf;
      opts.presolve = "off";
      opts.simplex_scale_strategy = 2;
   }
   void
   setTimeLimit( double tlim ) override
   {
      opts.time_limit = tlim;
   }

   void
   setUp( const Problem<REAL>& problem, const Vec<int>& row_maps,
          const Vec<int>& col_maps ) override
   {
      int ncols = problem.getNCols();
      int nrows = problem.getNRows();
      const Vec<String>& varNames = problem.getVariableNames();
      const Vec<String>& consNames = problem.getConstraintNames();
      const VariableDomains<REAL>& domains = problem.getVariableDomains();
      const Objective<REAL>& obj = problem.getObjective();
      const ConstraintMatrix<REAL>& consMatrix = problem.getConstraintMatrix();
      const auto& lhs_values = consMatrix.getLeftHandSides();
      const auto& rhs_values = consMatrix.getRightHandSides();
      const auto& rflags = problem.getRowFlags();

      opts.simplex_strategy = 3;
      opts.highs_min_threads = 4;
      opts.highs_max_threads = 4;
      opts.simplex_update_limit = 500;
      HighsLp model;

      /* set the objective sense and offset */
      model.sense_ = OBJSENSE_MINIMIZE;
      model.offset_ = double( -obj.offset );

      model.numRow_ = nrows;
      model.numCol_ = ncols;
      model.nnz_ = consMatrix.getNnz();
      model.numInt_ = 0;

      // model.col_names_.resize( ncols );
      model.colCost_.resize( ncols );
      model.colLower_.resize( ncols );
      model.colUpper_.resize( ncols );

      // model.row_names_.resize( numrows );
      model.rowLower_.resize( nrows );
      model.rowUpper_.resize( nrows );

      for( int i = 0; i != nrows; ++i )
      {
         model.rowLower_[i] = rflags[i].test( RowFlag::kLhsInf )
                                  ? -inf
                                  : double( lhs_values[i] );
         model.rowUpper_[i] =
             rflags[i].test( RowFlag::kRhsInf ) ? inf : double( rhs_values[i] );
      }

      model.integrality_.resize( ncols );
      model.Aindex_.resize( model.nnz_ );
      model.Avalue_.resize( model.nnz_ );
      model.Astart_.resize( ncols + 1 );
      model.Astart_[ncols] = model.nnz_;

      int start = 0;

      for( int i = 0; i < ncols; ++i )
      {
         assert( !domains.flags[i].test( ColFlag::kInactive ) );

         model.colLower_[i] = domains.flags[i].test( ColFlag::kLbInf )
                                  ? -inf
                                  : double( domains.lower_bounds[i] );
         model.colUpper_[i] = domains.flags[i].test( ColFlag::kUbInf )
                                  ? inf
                                  : double( domains.upper_bounds[i] );

         model.colCost_[i] = double( obj.coefficients[i] );

         auto colvec = consMatrix.getColumnCoefficients( i );

         int collen = colvec.getLength();
         const int* colrows = colvec.getIndices();
         const REAL* colvals = colvec.getValues();

         model.Astart_[i] = start;

         for( int k = 0; k != collen; ++k )
         {
            model.Avalue_[start + k] = double( colvals[k] );
            model.Aindex_[start + k] = double( colrows[k] );
         }

         start += collen;
      }

      solver.passModel( model );
   }

   void
   setUp( const Problem<REAL>& problem, const Vec<int>& row_maps,
          const Vec<int>& col_maps, const Components& components,
          const ComponentInfo& component ) override
   {
      int ncols = problem.getNCols();
      int nrows = problem.getNRows();
      const Vec<String>& varNames = problem.getVariableNames();
      const Vec<String>& consNames = problem.getConstraintNames();
      const VariableDomains<REAL>& domains = problem.getVariableDomains();
      const Objective<REAL>& obj = problem.getObjective();
      const auto& consMatrix = problem.getConstraintMatrix();
      const auto& lhs_values = consMatrix.getLeftHandSides();
      const auto& rhs_values = consMatrix.getRightHandSides();
      const auto& rflags = problem.getRowFlags();
      const int* rowset = components.getComponentsRows( component.componentid );
      const int* colset = components.getComponentsCols( component.componentid );
      int numrows = components.getComponentsNumRows( component.componentid );
      int numcols = components.getComponentsNumCols( component.componentid );

      HighsLp model;

      /* set the objective sense and offset */
      model.sense_ = OBJSENSE_MINIMIZE;
      model.offset_ = 0;

      model.numRow_ = numrows;
      model.numCol_ = numcols;
      model.nnz_ = component.nnonz;
      model.numInt_ = 0;

      model.colCost_.resize( numcols );
      model.colLower_.resize( numcols );
      model.colUpper_.resize( numcols );

      model.rowLower_.resize( numrows );
      model.rowUpper_.resize( numrows );

      for( int i = 0; i != numrows; ++i )
      {
         int row = rowset[i];

         assert( components.getRowComponentIdx( row ) == i );

         model.rowLower_[i] = rflags[row].test( RowFlag::kLhsInf )
                                  ? -inf
                                  : double( lhs_values[row] );
         model.rowUpper_[i] = rflags[row].test( RowFlag::kRhsInf )
                                  ? inf
                                  : double( rhs_values[row] );
      }

      model.Aindex_.resize( model.nnz_ );
      model.Avalue_.resize( model.nnz_ );
      model.Astart_.resize( numcols + 1 );
      model.Astart_[numcols] = model.nnz_;

      int start = 0;

      for( int i = 0; i != numcols; ++i )
      {
         int col = colset[i];

         assert( components.getColComponentIdx( col ) == i );
         assert( !domains.flags[col].test( ColFlag::kInactive ) );

         model.colLower_[i] = domains.flags[col].test( ColFlag::kLbInf )
                                  ? -inf
                                  : double( domains.lower_bounds[col] );
         model.colUpper_[i] = domains.flags[col].test( ColFlag::kUbInf )
                                  ? inf
                                  : double( domains.upper_bounds[col] );

         model.colCost_[i] = double( obj.coefficients[col] );

         auto colvec = consMatrix.getColumnCoefficients( col );

         int collen = colvec.getLength();
         const int* colrows = colvec.getIndices();
         const REAL* colvals = colvec.getValues();

         model.Astart_[i] = start;

         for( int k = 0; k != collen; ++k )
         {
            model.Avalue_[start + k] = double( colvals[k] );
            model.Aindex_[start + k] =
                components.getRowComponentIdx( colrows[k] );
         }

         start += collen;
      }

      solver.passModel( model );
   }

   void
   solve() override
   {
      solver.passHighsOptions( opts );

      if( solver.run() == HighsStatus::Error )
      {
         this->status = SolverStatus::kError;
         return;
      }

      HighsModelStatus stat = solver.getModelStatus( true );

      switch( stat )
      {
      default:
         this->status = SolverStatus::kError;
         return;
      case HighsModelStatus::PRIMAL_INFEASIBLE:
         this->status = SolverStatus::kInfeasible;
         return;
      case HighsModelStatus::PRIMAL_UNBOUNDED:
         this->status = SolverStatus::kUnbounded;
         return;
      case HighsModelStatus::REACHED_TIME_LIMIT:
      case HighsModelStatus::REACHED_ITERATION_LIMIT:
         this->status = SolverStatus::kInterrupted;
         return;
      case HighsModelStatus::OPTIMAL:
         this->status = SolverStatus::kOptimal;
      }
   }

   void
   setVerbosity( VerbosityLevel verbosity ) override
   {
      switch( verbosity )
      {
      case VerbosityLevel::kQuiet:
         opts.message_level = ML_NONE;
         opts.logfile = nullptr;
         opts.output = nullptr;
         solver.setHighsOutput( nullptr );
         solver.setHighsLogfile( nullptr );
         break;
      case VerbosityLevel::kError:
      case VerbosityLevel::kWarning:
      case VerbosityLevel::kInfo:
      case VerbosityLevel::kExtra:
         opts.message_level = ML_MINIMAL;
      }
   }

   REAL
   getDualBound() override
   {
      if( this->status == SolverStatus::kOptimal )
         return -inf; // todo

      return -inf;
   }

   bool
   getSolution( Solution<REAL>& sol ) override
   {
      const HighsSolution& highsSol = solver.getSolution();
      int numcols = solver.getNumCols();
      int numrows = solver.getNumRows();

      if( highsSol.col_value.size() != numcols ||
          highsSol.col_dual.size() != numcols ||
          highsSol.row_dual.size() != numrows ||
          highsSol.row_value.size() != numrows )
         return false;

      // skip objoffset var
      --numcols;

      // get primal values
      sol.primal.resize( numcols );
      for( int i = 0; i != numcols; ++i )
         sol.primal[i] = highsSol.col_value[i];

      // return if no dual requested
      if( sol.type == SolutionType::kPrimal )
         return true;

      // get reduced costs
      sol.col_dual.resize( numcols );
      for( int i = 0; i != numcols; ++i )
         sol.col_dual[i] = REAL( highsSol.col_dual[i] );

      // get row duals
      sol.row_dual.resize( numrows );
      for( int i = 0; i != numrows; ++i )
         sol.row_dual[i] = REAL( highsSol.row_dual[i] );

      return true;
   }

   bool
   getSolution( const Components& components, int component,
                Solution<REAL>& sol ) override
   {
      if( this->status != SolverStatus::kOptimal )
         return false;

      int numcols = solver.getNumCols();
      int numrows = solver.getNumRows();
      const HighsSolution& highsSol = solver.getSolution();
      if( highsSol.col_value.size() != numcols ||
          highsSol.col_dual.size() != numcols ||
          highsSol.row_dual.size() != numrows ||
          highsSol.row_value.size() != numrows )
         return false;

      assert( components.getComponentsNumCols( component ) == numcols );

      const int* compcols = components.getComponentsCols( component );
      for( int i = 0; i != numcols; ++i )
         sol.primal[compcols[i]] = REAL( highsSol.col_value[i] );

      if( sol.type == SolutionType::kPrimal )
         return true;

      for( int i = 0; i != numcols; ++i )
         sol.col_dual[compcols[i]] = REAL( highsSol.col_dual[i] );

      const int* comprows = components.getComponentsRows( component );

      for( int i = 0; i != numrows; ++i )
         sol.row_dual[comprows[i]] = REAL( highsSol.row_dual[i] );

      return true;
   }

   void
   addParameters( ParameterSet& paramSet ) override
   {
   }

   SolverType
   getType() override
   {
      return SolverType::LP;
   }

   String
   getName() override
   {
      return "HiGHS";
   }

   void
   printDetails() override
   {
   }
};

template <typename REAL>
class HighsFactory : public SolverFactory<REAL>
{
   HighsFactory() {}

 public:
   virtual std::unique_ptr<SolverInterface<REAL>>
   newSolver( VerbosityLevel verbosity ) const override
   {
      auto highs =
          std::unique_ptr<SolverInterface<REAL>>( new HighsInterface<REAL>() );

      highs->setVerbosity( verbosity );

      return std::move( highs );
   }

   static std::unique_ptr<SolverFactory<REAL>>
   create()
   {
      return std::unique_ptr<SolverFactory<REAL>>( new HighsFactory<REAL>() );
   }
};

} // namespace papilo

#endif
