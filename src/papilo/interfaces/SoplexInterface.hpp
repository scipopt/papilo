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

#ifndef _PAPILO_INTERFACES_SOPLEX_INTERFACE_HPP_
#define _PAPILO_INTERFACES_SOPLEX_INTERFACE_HPP_

#include "papilo/misc/String.hpp"
#include "papilo/misc/Vec.hpp"
#include <cassert>
#include <stdexcept>
#include <type_traits>

#include "papilo/core/Problem.hpp"
#include "papilo/interfaces/SolverInterface.hpp"
#include "soplex.h"

namespace papilo
{

template <typename REAL>
class SoplexInterface : public SolverInterface<REAL>
{
 private:
   soplex::SoPlex spx;

 public:
   void
   readSettings( const String& file ) override
   {
      if( !spx.loadSettingsFile( file.c_str() ) )
         this->status = SolverStatus::kError;
   }

   soplex::SoPlex&
   getSoPlex()
   {
      return spx;
   }

   void
   setTimeLimit( double tlim ) override
   {
      using namespace soplex;

      spx.setIntParam( SoPlex::TIMER, SoPlex::TIMER_WALLCLOCK );
      spx.setRealParam( SoPlex::TIMELIMIT, Real( tlim ) );
   }

   void
   setUp( const Problem<REAL>& problem, const Vec<int>& row_maps,
          const Vec<int>& col_maps ) override
   {
      using namespace soplex;

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

      /* set the objective sense and offset */
      spx.setIntParam( SoPlex::OBJSENSE, SoPlex::OBJSENSE_MINIMIZE );

      if( obj.offset != 0 )
         spx.setRealParam( SoPlex::OBJ_OFFSET, Real( obj.offset ) );

      LPRowSet rows( nrows );
      LPColSet cols( ncols );
      DSVector vec( ncols );
      for( int i = 0; i < nrows; ++i )
      {
         Real lhs = rflags[i].test( RowFlag::kLhsInf ) ? -infinity
                                                       : Real( lhs_values[i] );
         Real rhs = rflags[i].test( RowFlag::kRhsInf ) ? infinity
                                                       : Real( rhs_values[i] );

         rows.add( lhs, vec, rhs );
      }

      spx.addRowsReal( rows );

      for( int i = 0; i < ncols; ++i )
      {
         assert( !domains.flags[i].test( ColFlag::kInactive ) );

         Real lb = domains.flags[i].test( ColFlag::kLbInf )
                       ? -infinity
                       : Real( domains.lower_bounds[i] );
         Real ub = domains.flags[i].test( ColFlag::kUbInf )
                       ? infinity
                       : Real( domains.upper_bounds[i] );

         auto colvec = consMatrix.getColumnCoefficients( i );

         int collen = colvec.getLength();
         const int* colrows = colvec.getIndices();
         const REAL* colvals = colvec.getValues();

         vec.clear();

         if( std::is_same<REAL, Real>::value )
         {
            vec.add( collen, colrows, (const Real*)colvals );
         }
         else
         {
            for( int i = 0; i != collen; ++i )
               vec.add( colrows[i], Real( colvals[i] ) );
         }

         cols.add( Real( obj.coefficients[i] ), lb, vec, ub );
      }

      spx.addColsReal( cols );
   }

   void
   setUp( const Problem<REAL>& problem, const Vec<int>& row_maps,
          const Vec<int>& col_maps, const Components& components,
          const ComponentInfo& component ) override
   {
      using namespace soplex;

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

      /* set the objective sense and offset */
      spx.setIntParam( SoPlex::OBJSENSE, SoPlex::OBJSENSE_MINIMIZE );

      LPRowSet rows( numrows );
      LPColSet cols( numcols );
      DSVector vec( numcols );
      for( int i = 0; i != numrows; ++i )
      {
         int row = rowset[i];

         assert( components.getRowComponentIdx( row ) == i );

         Real lhs = rflags[row].test( RowFlag::kLhsInf )
                        ? -infinity
                        : Real( lhs_values[row] );
         Real rhs = rflags[row].test( RowFlag::kRhsInf )
                        ? infinity
                        : Real( rhs_values[row] );

         rows.add( lhs, vec, rhs );
      }

      spx.addRowsReal( rows );

      for( int i = 0; i != numcols; ++i )
      {
         int col = colset[i];

         assert( components.getColComponentIdx( col ) == i );
         assert( !domains.flags[col].test( ColFlag::kInactive ) );

         Real lb = domains.flags[col].test( ColFlag::kLbInf )
                       ? -infinity
                       : Real( domains.lower_bounds[col] );
         Real ub = domains.flags[col].test( ColFlag::kUbInf )
                       ? infinity
                       : Real( domains.upper_bounds[col] );

         auto colvec = consMatrix.getColumnCoefficients( col );

         int collen = colvec.getLength();
         const int* colrows = colvec.getIndices();
         const REAL* colvals = colvec.getValues();

         vec.clear();

         for( int j = 0; j != collen; ++j )
            vec.add( components.getRowComponentIdx( colrows[j] ),
                     Real( colvals[j] ) );

         cols.add( Real( obj.coefficients[col] ), lb, vec, ub );
      }

      spx.addColsReal( cols );
   }

   void
   solve() override
   {
      using namespace soplex;

      assert( this->status != SolverStatus::kError );

      spx.setSettings( spx.settings() );

      SPxSolver::Status stat = spx.optimize();

      switch( stat )
      {
      default:
         this->status = SolverStatus::kError;
         return;
      case SPxSolver::Status::INForUNBD:
         this->status = SolverStatus::kUnbndOrInfeas;
         return;
      case SPxSolver::Status::INFEASIBLE:
         this->status = SolverStatus::kInfeasible;
         return;
      case SPxSolver::Status::UNBOUNDED:
         this->status = SolverStatus::kUnbounded;
         return;
      case SPxSolver::Status::ABORT_CYCLING:
         this->status = SolverStatus::kInterrupted;
         return;
      case SPxSolver::Status::OPTIMAL_UNSCALED_VIOLATIONS:
      case SPxSolver::Status::OPTIMAL:
         this->status = SolverStatus::kOptimal;
      }
   }

   void
   setVerbosity( VerbosityLevel verbosity ) override
   {
      using namespace soplex;

      switch( verbosity )
      {
      case VerbosityLevel::kQuiet:
      case VerbosityLevel::kError:
         spx.setIntParam( SoPlex::VERBOSITY, SoPlex::VERBOSITY_ERROR );
         break;
      case VerbosityLevel::kWarning:
         spx.setIntParam( SoPlex::VERBOSITY, SoPlex::VERBOSITY_WARNING );
         break;
      case VerbosityLevel::kInfo:
         spx.setIntParam( SoPlex::VERBOSITY, SoPlex::VERBOSITY_NORMAL );
         break;
      case VerbosityLevel::kExtra:
         spx.setIntParam( SoPlex::VERBOSITY, SoPlex::VERBOSITY_HIGH );
      }
   }

   REAL
   getDualBound() override
   {
      if( spx.hasPrimal() )
         return spx.objValueReal();
      else
         return -soplex::infinity;
   }

   bool
   getSolution( Solution<REAL>& sol ) override
   {
      Vec<soplex::Real> buffer;

      int numcols = spx.numColsReal();
      buffer.resize( numcols );

      if( !spx.getPrimalReal( buffer.data(), numcols ) )
         return false;

      sol.primal.resize( numcols );
      for( int i = 0; i != numcols; ++i )
         sol.primal[i] = REAL( buffer[i] );

      if( sol.type == SolutionType::kPrimal )
         return true;

      if( !spx.getRedCostReal( buffer.data(), numcols ) )
         return false;

      sol.col_dual.resize( numcols );
      for( int i = 0; i != numcols; ++i )
         sol.col_dual[i] = REAL( buffer[i] );

      int numrows = spx.numRowsReal();

      buffer.resize( numrows );
      if( !spx.getDualReal( buffer.data(), numrows ) )
         return false;

      sol.row_dual.resize( numrows );
      for( int i = 0; i != numrows; ++i )
         sol.row_dual[i] = REAL( buffer[i] );

      return true;
   }

   bool
   getSolution( const Components& components, int component,
                Solution<REAL>& sol ) override
   {
      Vec<soplex::Real> buffer;

      int numcols = spx.numColsReal();
      assert( components.getComponentsNumCols( component ) ==
              spx.numColsReal() );

      buffer.resize( numcols );
      if( !spx.getPrimalReal( buffer.data(), numcols ) )
         return false;

      const int* compcols = components.getComponentsCols( component );
      for( int i = 0; i != numcols; ++i )
         sol.primal[compcols[i]] = REAL( buffer[i] );

      if( sol.type == SolutionType::kPrimal )
         return true;

      if( !spx.getRedCostReal( buffer.data(), numcols ) )
         return false;

      for( int i = 0; i != numcols; ++i )
         sol.col_dual[compcols[i]] = REAL( buffer[i] );

      int numrows = spx.numRowsReal();
      buffer.resize( numrows );
      const int* comprows = components.getComponentsRows( component );

      if( !spx.getDualReal( buffer.data(), numrows ) )
         return false;

      for( int i = 0; i != numrows; ++i )
         sol.row_dual[comprows[i]] = REAL( buffer[i] );

      return true;
   }

   void
   addParameters( ParameterSet& paramSet ) override
   {
      using namespace soplex;

      SoPlex::Settings& settings =
          const_cast<SoPlex::Settings&>( spx.settings() );

      for( int i = 0; i != SoPlex::BOOLPARAM_COUNT; ++i )
      {
         paramSet.addParameter( settings.boolParam.name[i].c_str(),
                                settings.boolParam.description[i].c_str(),
                                settings._boolParamValues[i] );
      }

      for( int i = 0; i != SoPlex::INTPARAM_COUNT; ++i )
      {
         paramSet.addParameter( settings.intParam.name[i].c_str(),
                                settings.intParam.description[i].c_str(),
                                settings._intParamValues[i],
                                settings.intParam.lower[i],
                                settings.intParam.upper[i] );
      }

      for( int i = 0; i != SoPlex::REALPARAM_COUNT; ++i )
      {
         paramSet.addParameter( settings.realParam.name[i].c_str(),
                                settings.realParam.description[i].c_str(),
                                settings._realParamValues[i],
                                settings.realParam.lower[i],
                                settings.realParam.upper[i] );
      }
   }

   SolverType
   getType() override
   {
      return SolverType::LP;
   }

   String
   getName() override
   {
      return "SoPlex";
   }

   void
   printDetails() override
   {
      spx.printStatistics( std::cout );
   }
};

template <typename REAL>
class SoplexFactory : public SolverFactory<REAL>
{
   void ( *soplexsetup )( soplex::SoPlex& soplex, void* usrdata );
   void* soplexsetup_usrdata;

   SoplexFactory( void ( *soplexsetup )( soplex::SoPlex& soplex,
                                         void* usrdata ),
                  void* soplexsetup_usrdata )
       : soplexsetup( soplexsetup ), soplexsetup_usrdata( soplexsetup_usrdata )
   {
   }

 public:
   virtual std::unique_ptr<SolverInterface<REAL>>
   newSolver( VerbosityLevel verbosity ) const override
   {
      auto soplex =
          std::unique_ptr<SolverInterface<REAL>>( new SoplexInterface<REAL>() );

      // set verbosity already before the setup function call
      soplex->setVerbosity( verbosity );

      if( soplexsetup != nullptr )
         soplexsetup(
             static_cast<SoplexInterface<REAL>*>( soplex.get() )->getSoPlex(),
             soplexsetup_usrdata );

      // set verbosity again in case the setup function altered it
      soplex->setVerbosity( verbosity );

      return std::move( soplex );
   }

   static std::unique_ptr<SolverFactory<REAL>>
   create( void ( *soplexsetup )( soplex::SoPlex& soplex,
                                  void* usrdata ) = nullptr,
           void* soplexsetup_usrdata = nullptr )
   {
      return std::unique_ptr<SolverFactory<REAL>>(
          new SoplexFactory<REAL>( soplexsetup, soplexsetup_usrdata ) );
   }
};

} // namespace papilo

#endif
