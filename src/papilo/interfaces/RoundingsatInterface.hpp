/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*               This file is part of the program and library                */
/*    PaPILO --- Parallel Presolve for Integer and Linear Optimization       */
/*                                                                           */
/* Copyright (C) 2020-2025 Zuse Institute Berlin (ZIB)                       */
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

#ifndef _PAPILO_INTERFACES_ROUNDINGSAT_INTERFACE_HPP_
#define _PAPILO_INTERFACES_ROUNDINGSAT_INTERFACE_HPP_

#include "auxiliary.hpp"
#include "globals.hpp"
#include "parsing.hpp"
#include "run.hpp"
#include <csignal>
#include <fstream>
#include <memory>

namespace rs
{
bool asynch_interrupt;
Options options;
Stats stats;
} // namespace rs

namespace papilo
{

template <typename REAL>
class RoundingsatInterface : public SolverInterface<REAL>
{
   //   Num<REAL>& num;
   rs::CeArb objective;
   Vec<int> scaling_row_factor;

   int
   doSetUp( const Problem<REAL>& problem, const Vec<int>& origRowMap,
            const Vec<int>& origColMap )
   {
      Num<REAL> num{};
      objective = rs::run::solver.cePools.takeArb();
      assert( objective->isReset() );
      rs::CeArb input = rs::run::solver.cePools.takeArb();

      const Vec<REAL>& obj = problem.getObjective().coefficients;
      const Vec<REAL>& rhs = problem.getConstraintMatrix().getRightHandSides();
      const Vec<REAL>& lhs = problem.getConstraintMatrix().getLeftHandSides();
      const auto consMatrix = problem.getConstraintMatrix();

      input->reset();
      for( int col = 0; col < problem.getNCols(); ++col )
      {
         assert( num.isIntegral( obj[col] ) );
         if( obj[col] == 0 )
            continue;
         int sat_var_index = get_index_of_variable_name(
             problem.getVariableNames(), origColMap, col );
         if(sat_var_index < 1)
            return -1;
         rs::BigCoef coeff = rs::BigCoef ( obj[col] );
         rs::run::solver.setNbVars( abs( sat_var_index ), true );
         input->addLhs( coeff, sat_var_index );
      }
      input->copyTo( objective );

      for( int row = 0; row < problem.getNRows(); ++row )
      {
         input->reset();
         auto row_coeff =
             problem.getConstraintMatrix().getRowCoefficients( row );
         if( consMatrix.getRowFlags()[row].test( RowFlag::kEquation ) )
         {
            int ret = map_cons_to_lhs( input, row_coeff, row, lhs[row], num, problem, origColMap );
            if( ret < 0 || rs::run::solver.addConstraint( input, rs::Origin::FORMULA )
                    .second == rs::ID_Unsat )
            {
               fmt::print( "An error occurred\n" );
               return -1;
            }
            input->invert();
            const std::pair<rs::ID, rs::ID>& pair =
                rs::run::solver.addConstraint( input, rs::Origin::FORMULA );
            if( pair.second == rs::ID_Unsat )
            {
               fmt::print( "An error occurred\n" );
               return -1;
            }
            continue;
         }

         if( !consMatrix.getRowFlags()[row].test( RowFlag::kLhsInf ) )
         {
            int ret = map_cons_to_lhs( input, row_coeff, row, lhs[row], num, problem, origColMap );
            const std::pair<rs::ID, rs::ID>& pair =
                rs::run::solver.addConstraint( input, rs::Origin::FORMULA );
            if( ret < 0 || pair.second == rs::ID_Unsat )
            {
               fmt::print( "An error occurred\n" );
               return -1;
            }
         }
         if( !consMatrix.getRowFlags()[row].test( RowFlag::kRhsInf ) )
         {
            if(!consMatrix.getRowFlags()[row].test( RowFlag::kLhsInf ))
               input->reset();
            int ret = map_cons_to_rhs( input, row_coeff, row, rhs[row], num, problem, origColMap );
            const std::pair<rs::ID, rs::ID>& pair =
                rs::run::solver.addConstraint( input, rs::Origin::FORMULA );
            if( ret < 0 || pair.second == rs::ID_Unsat )
            {
               fmt::print( "An error occurred\n" );
               return -1;
            }
         }
      }
      //TODO: test symmetries to problem
      for( size_t symmetry_index = 0; symmetry_index < problem.getSymmetries().symmetries.size(); ++symmetry_index )
      {
         assert( problem.test_problem_type(ProblemFlag::kBinary) );
         input->reset();
         auto& symmetry = problem.getSymmetries().symmetries[symmetry_index];
         int col1 = get_index_of_variable_name( problem.getVariableNames(), origColMap,
                                                symmetry.getDominatingCol() );
         int col2 = get_index_of_variable_name( problem.getVariableNames(), origColMap,
                                                symmetry.getDominatedCol() );
         if(col1 < 1 || col2 < 1)
            return -1;
         switch( symmetry.getSymmetryType() )
         {
         case SymmetryType::kXgeY:
            input->addLhs( 1, col1 );
            input->addLhs( -1, col2 );
            input->addRhs( 0 );
            break;
         case SymmetryType::kXplusYge1:
            input->addLhs( 1, col1 );
            input->addLhs( 1, col2 );
            input->addRhs( 1 );
            break;
         default:
            assert( false );
         }

         const std::pair<rs::ID, rs::ID>& pair =
             rs::run::solver.addConstraint( input, rs::Origin::FORMULA );
         if( pair.second == rs::ID_Unsat )
         {
            fmt::print( "An error occurred\n" );
            return -1;
         }
      }
      input->reset();

      return 0;
   }

   int
   get_index_of_variable_name( const Vec<String>& _names,
                               const Vec<int>& origColMap, int col) const
   {
      auto name = _names[origColMap[col]];
      int sat_var_index;
      if (name.empty() || name[0] != 'x')
      {
         sat_var_index = 0;
         fmt::print("RoundingSAT: variablename must start with an 'x'!\n");
      }
      else
         sat_var_index = atoi(name.substr(1).c_str());
      return sat_var_index;
   }

   rs::BigCoef
   to_int( REAL val, REAL scale, Num<REAL>& num )
   {
      if( num.isZero( val ) )
         return 0;
      assert( num.isIntegral( val * scale ) );
      return rs::BigCoef ( val * scale );
   }

   int
   map_cons_to_lhs( rs::CeArb& input, const SparseVectorView<REAL>& row_coeff,
                    int row, REAL lhs, Num<REAL>& num,
                    const Problem<REAL>& problem, const Vec<int>& map )
   {
      bool scale_necessary = !num.isIntegral(lhs);
      for( int j = 0; j < row_coeff.getLength(); ++j )
      {
         if(scale_necessary)
            break;
         scale_necessary = !num.isIntegral(row_coeff.getValues()[j]);
      }
      REAL scale = scale_necessary ? abs(scaling_row_factor[row]) : 1;
      for( int j = 0; j < row_coeff.getLength(); ++j )
      {
         int sat_var_index = get_index_of_variable_name(
             problem.getVariableNames(), map, row_coeff.getIndices()[j]);
         if(sat_var_index < 1)
            return -1;
         assert( num.isIntegral( row_coeff.getValues()[j] * scale )  );
         rs::BigCoef coeff = to_int( row_coeff.getValues()[j], scale, num);
         rs::run::solver.setNbVars( abs( sat_var_index ), true );
         input->addLhs( coeff, sat_var_index );
      }
      input->addRhs( to_int( lhs, scale, num ) );
      return 0;
   }


   int
   map_cons_to_rhs( rs::CeArb& input, const SparseVectorView<REAL>& row_coeff, int row, REAL rhs, Num<REAL>& num,
                    const Problem<REAL>& problem, const Vec<int>& map )
   {
      bool scale_necessary = !num.isIntegral(rhs);
      for( int j = 0; j < row_coeff.getLength(); ++j )
      {
         if(scale_necessary)
            break;
         scale_necessary = !num.isIntegral(row_coeff.getValues()[j]);
      }
      REAL scale = scale_necessary ? abs(scaling_row_factor[row]) : 1;
      for( int j = 0; j < row_coeff.getLength(); ++j )
      {
         int sat_var_index = get_index_of_variable_name(
             problem.getVariableNames(), map, row_coeff.getIndices()[j]);
         if(sat_var_index < 1)
            return -1;
         assert( num.isIntegral( -row_coeff.getValues()[j] * scale )  );
         rs::BigCoef coeff = to_int( -row_coeff.getValues()[j], scale, num);
         rs::run::solver.setNbVars( abs( sat_var_index ), true );
         input->addLhs( coeff, sat_var_index );
      }
      input->addRhs( -to_int( rhs, scale, num ) );
      return 0;
   }

   int
   doSetUp( const Problem<REAL>& problem, const Vec<int>& origRowMap,
            const Vec<int>& origColMap, const Components& components,
            const ComponentInfo& component )
   {
      // TODO: implement
      return -1;
   }

 public:
   RoundingsatInterface()
   {
      // TODO: num = {};
      rs::run::solver.init();
   }

   void
   setUp( const Problem<REAL>& prob, const Vec<int>& row_maps,
          const Vec<int>& col_maps, const Components& components,
          const ComponentInfo& component ) override
   {
      if( doSetUp( prob, row_maps, col_maps, components, component ) != 0 )
         this->status = SolverStatus::kError;
   }

   void
   setNodeLimit( int nodelimit ) override
   {
   }

   void
   setGapLimit( const REAL& gaplim ) override
   {
   }

   void
   setRowScalingFactor( Vec<int> _scaling_row_factor )
   {
      assert(scaling_row_factor.size() <= _scaling_row_factor.size());
      scaling_row_factor = _scaling_row_factor;
   }

   void
   setSoftTimeLimit( double tlim ) override
   {
   }

   void
   setTimeLimit( double tlim ) override
   {
      rs::options.time_limit.parse( std::to_string(tlim));
   }

   void
   setVerbosity( VerbosityLevel verbosity ) override
   {
   }

   void
   setUp( const Problem<REAL>& prob, const Vec<int>& row_maps,
          const Vec<int>& col_maps ) override
   {
      if( doSetUp( prob, row_maps, col_maps ) != 0 )
         this->status = SolverStatus::kError;
   }

   void
   solve() override
   {
      if(this->status == SolverStatus::kError)
         return ;

      rs::run::solver.initLP( objective );
      rs::run::run( objective );

      bool satisfiable = rs::run::solver.foundSolution();
      bool time_expired = rs::stats.getTime() > rs::options.time_limit.get();
      if( time_expired && !satisfiable )
         // postsolving not possible -> so return error to avoid it
         this->status = SolverStatus::kError;
      else if( time_expired )
         this->status = SolverStatus::kInterrupted;
      else
         this->status =
             satisfiable ? SolverStatus::kOptimal : SolverStatus::kInfeasible;
   }

   REAL
   getDualBound() override
   {
      return 0;
   }

   bool
   getSolution( Solution<REAL>& sol_buffer, PostsolveStorage<REAL>& postsolve ) override
   {

      std::vector<rs::Lit>& sol = rs::run::solver.lastSol;
      Vec<REAL> primal{};

      int ncols = postsolve.origcol_mapping.size();
      auto names = postsolve.getOriginalProblem().getVariableNames();
      auto origcol_mapping = postsolve.origcol_mapping;
      primal.resize( ncols );

      for( int red_col = 0; red_col < ncols ; ++red_col )
      {
         int orig_col = get_index_of_variable_name( names, origcol_mapping, red_col );
         assert( abs( sol[orig_col] ) == orig_col );
         primal[red_col] = sol[orig_col] > 0 ? 1 : 0;
      }
      sol_buffer = Solution<REAL>( primal );
      return true;
   }

   bool
   getSolution( const Components& components, int component,
                Solution<REAL>& solbuffer ) override
   {

      return false;
   }

   SolverType
   getType() override
   {
      return SolverType::PseudoBoolean;
   }

   String
   getName() override
   {
      return "RoundingSat";
   }

   void
   printDetails() override
   {
   }

   bool
   is_dual_solution_available() override
   {
      return false;
   }

   void
   addParameters( ParameterSet& paramSet ) override
   {
   }

   ~RoundingsatInterface() = default;
};

template <typename REAL>
class RoundingsatFactory : public SolverFactory<REAL>
{
   RoundingsatFactory() = default;

 public:
   std::unique_ptr<SolverInterface<REAL>>
   newSolver( VerbosityLevel verbosity ) const override
   {
      auto satsolver = std::unique_ptr<SolverInterface<REAL>>(
          new RoundingsatInterface<REAL>() );

      auto res = std::move( satsolver );
      return res;
   }

   virtual void
   add_parameters( ParameterSet& parameter ) const
   {
   }

   static std::unique_ptr<SolverFactory<REAL>>
   create()
   {
      return std::unique_ptr<SolverFactory<REAL>>(
          new RoundingsatFactory<REAL>() );
   }
};

} // namespace papilo

#endif
