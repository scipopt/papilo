/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*               This file is part of the program and library                */
/*    PaPILO --- Parallel Presolve for Integer and Linear Optimization       */
/*                                                                           */
/* Copyright (C) 2020-2023 Konrad-Zuse-Zentrum                               */
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
      assert( objective->isReset() );
      rs::CeArb input = rs::run::solver.cePools.takeArb();

      const Vec<REAL>& obj = problem.getObjective().coefficients;
      const Vec<REAL>& rhs = problem.getConstraintMatrix().getRightHandSides();
      const Vec<REAL>& lhs = problem.getConstraintMatrix().getLeftHandSides();
      const auto consMatrix = problem.getConstraintMatrix();

      input->reset();
      for( int col = 0; col < problem.getNCols(); ++col )
      {
         input->reset();
         assert( num.isIntegral( obj[col] ) );
         if( obj[col] == 0 )
            continue;
         int sat_var_index = col + 1;
         int coeff = to_int( obj[col] );
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
            map_cons_to_lhs( input, row_coeff, row );
            input->addRhs( to_int( lhs[row] * scaling_row_factor[row] ) );
            if( rs::run::solver.addConstraint( input, rs::Origin::FORMULA )
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
            map_cons_to_lhs( input, row_coeff, row );
            input->addRhs( to_int( lhs[row] * scaling_row_factor[row] ) );
            const std::pair<rs::ID, rs::ID>& pair =
                rs::run::solver.addConstraint( input, rs::Origin::FORMULA );
            if( pair.second == rs::ID_Unsat )
            {
               fmt::print( "An error occurred\n" );
               return -1;
            }
         }
         if( !consMatrix.getRowFlags()[row].test( RowFlag::kRhsInf ) )
         {
            map_cons_to_lhs( input, row_coeff, row );
            input->addRhs( to_int( rhs[row] * scaling_row_factor[row] ) );
            input->invert();
            const std::pair<rs::ID, rs::ID>& pair =
                rs::run::solver.addConstraint( input, rs::Origin::FORMULA );
            if( pair.second == rs::ID_Unsat )
            {
               fmt::print( "An error occurred\n" );
               return -1;
            }
         }
      }
      return 0;
   }

   int
   to_int( REAL val )
   {
      //TODO: num
      Num<REAL> num;
      assert( num.isIntegral(val) );
      return num.round_to_int(val);
   }

   void
   map_cons_to_lhs( rs::CeArb& input, const SparseVectorView<REAL>& row_coeff, int row )
   {
      Num<REAL> num{};
      for( int j = 0; j < row_coeff.getLength(); ++j )
      {
         int sat_var_index = row_coeff.getIndices()[j] + 1;
         assert( num.isIntegral( row_coeff.getValues()[j] * scaling_row_factor[row] )  );
         rs::bigint coeff = to_int( row_coeff.getValues()[j] * scaling_row_factor[row]);
         rs::run::solver.setNbVars( abs( sat_var_index ), true );
         input->addLhs( coeff, sat_var_index );
      }
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
      objective = rs::run::solver.cePools.takeArb();
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

      // TODO how to detect (infeasibility or) time limit reached
      bool satisfiable = rs::run::solver.foundSolution();
      this->status =
          satisfiable ? SolverStatus::kOptimal : SolverStatus::kError;
   }

   REAL
   getDualBound() override
   {
      return 0;
   }

   bool
   getSolution( Solution<REAL>& solbuffer ) override
   {

      std::vector<rs::Lit>& sol = rs::run::solver.lastSol;
      Vec<REAL> primal{};
      primal.resize( sol.size() - 1 );
      for( int v = 1; v < sol.size(); ++v )
      {
         assert( abs( sol[v] ) == v );
         primal[v - 1] = sol[v] > 0 ? 1 : 0;
      }
      solbuffer = Solution<REAL>( primal );
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

      return std::move( satsolver );
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
