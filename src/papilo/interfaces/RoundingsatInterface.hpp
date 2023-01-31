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

#include "Solver.hpp"
#include "typedefs.hpp"
#include "run.hpp"
#include <csignal>
#include <fstream>
#include <memory>

namespace papilo {

   template<typename REAL>
   class RoundingsatInterface : public SolverInterface<REAL> {

      int
      doSetUp(const Problem<REAL> &problem, const Vec<int> &origRowMap,
              const Vec<int> &origColMap) {

         rs::run::solver.init( );
         rs::CeArb objective = rs::run::solver.cePools.takeArb( );
         assert(objective->isReset( ));
         rs::CeArb input = rs::run::solver.cePools.takeArb( );

         const Vec<String> &varNames = problem.getVariableNames( );
         const VariableDomains<REAL> &domains = problem.getVariableDomains( );
         const Vec<REAL> &obj = problem.getObjective( ).coefficients;
         const Vec<REAL> &rhs = problem.getConstraintMatrix( ).getRightHandSides( );
         const Vec<REAL> &lhs = problem.getConstraintMatrix( ).getLeftHandSides( );
         const auto consMatrix = problem.getConstraintMatrix( );

         for( int row = 0; row < problem.getNRows( ); ++row )
         {
            auto row_coeff = problem.getConstraintMatrix( ).getRowCoefficients(row);
            if( !consMatrix.getRowFlags( )[ row ].test(RowFlag::kLhsInf))
            {
               for( int j = 0; j < row_coeff.getLength( ); ++j )
               {
                  String var_name = varNames[ origColMap[ row_coeff.getIndices( )[ j ]]];
                  assert(row_coeff.getIndices( )[ j ] != 0);
                  assert(domains.isBinary(row_coeff.getIndices( )[ j ]));
                  if( var_name.empty( ) || var_name[ 0 ] != 'x' )
                  {
                     fmt::print("Invalid literal token: {}", var_name);
                     return -1;
                  }

                  rs::Lit var_index = atoi(var_name.substr(1).c_str( ));
                  if( var_index < 1 )
                  {
                     fmt::print("Variable {} less than 1\n ", var_name);
                     return -1;
                  }
                  //TODO:
                  // assert(num.isInteger(row_coeff.getIndices()[j]));
                  int coeff = ( int ) row_coeff.getValues( )[ j ];
                  if( coeff < 0 )
                     var_index = -var_index;
                  rs::run::solver.setNbVars(std::abs(var_index), true);
                  input->addLhs(coeff, var_index);
               }
               input->addRhs(lhs[ row ]);
               if( rs::run::solver.addConstraint(input, rs::Origin::FORMULA).second == rs::ID_Unsat )
               {
                  fmt::print("An error occurred\n");
                  return -1;
               }
            }
            if( !consMatrix.getRowFlags( )[ row ].test(RowFlag::kRhsInf))
            {
               for( int j = 0; j < row_coeff.getLength( ); ++j )
               {
                  String var_name = varNames[ origColMap[ row_coeff.getIndices( )[ j ]]];
                  assert(row_coeff.getIndices( )[ j ] != 0);
                  assert(domains.isBinary(row_coeff.getIndices( )[ j ]));
                  if( var_name.empty( ) || var_name[ 0 ] != 'x' )
                  {
                     fmt::print("Invalid literal token: {}", var_name);
                     return -1;
                  }
                  rs::Lit var_index = atoi(var_name.substr(1).c_str( ));
                  if( var_index < 1 )
                  {
                     fmt::print("Variable {} less than 1\n ", var_name);
                     return -1;
                  }
                  //TODO:
                  // assert(num.isInteger(row_coeff.getIndices()[j]));
                  int coeff = ( int ) row_coeff.getValues( )[ j ];
                  if( coeff < 0 )
                     var_index = -var_index;
                  rs::run::solver.setNbVars(std::abs(var_index), true);
                  input->addLhs(coeff, var_index);
               }
               input->addRhs(rhs[ row ]);
               input->invert( );
               if( rs::run::solver.addConstraint(input, rs::Origin::FORMULA).second == rs::ID_Unsat )
               {
                  fmt::print("An error occurred\n");
                  return -1;
               }
            }

         }
      }

      int
      doSetUp(const Problem<REAL> &problem, const Vec<int> &origRowMap,
              const Vec<int> &origColMap, const Components &components,
              const ComponentInfo &component) {
      }

   public:
      RoundingsatInterface( ) { }

      void
      setUp(const Problem<REAL> &prob, const Vec<int> &row_maps,
            const Vec<int> &col_maps, const Components &components,
            const ComponentInfo &component) override {
         if( doSetUp(prob, row_maps, col_maps, components, component) != 0 )
            this->status = SolverStatus::kError;
      }

      void
      setNodeLimit(int
                   num) override {
      }

      void
      setGapLimit(
            const REAL &gaplim) override {
      }

      void
      setSoftTimeLimit(double
                       tlim) override {
      }

      void
      setTimeLimit(double
                   tlim) override {
      }

      void
      setVerbosity(VerbosityLevel
                   verbosity) override {
      }

      void
      setUp(const Problem<REAL> &prob, const Vec<int> &row_maps,
            const Vec<int> &col_maps) override {
         if( doSetUp(prob, row_maps, col_maps) != 0 )
            this->status = SolverStatus::kError;
      }

      void
      solve( )
      override {
      }

      REAL
      getDualBound( )
      override {
         return 0;
      }

      bool
      getSolution(Solution<REAL> &solbuffer)
      override {
         return false;
      }

      bool
      getSolution(const Components &components, int component,
                  Solution<REAL> &solbuffer) override {

         return false;
      }

      SolverType
      getType( )
      override {
         return SolverType::PseudoBoolean;
      }

      String
      getName( )
      override {
         return "RoundingSat";
      }

      void
      printDetails( )
      override {
      }

      bool
      is_dual_solution_available( )
      override {
         return false;
      }

      void
      addParameters(ParameterSet &paramSet)
      override {
      }

      ~RoundingsatInterface( ) =
      default;

   };

   template<typename REAL>
   class RoundingsatFactory : public SolverFactory<REAL> {
      RoundingsatFactory( ) = default;

   public:
      std::unique_ptr<SolverInterface<REAL>>
      newSolver(VerbosityLevel verbosity) const override {
         auto satsolver = std::unique_ptr<SolverInterface<REAL>>(
               new RoundingsatInterface<REAL>( ));

         return std::move(satsolver);
      }

      virtual void
      add_parameters(ParameterSet &parameter) const {
      }

      static std::unique_ptr<SolverFactory<REAL>>
      create( ) {
         return std::unique_ptr<SolverFactory<REAL>>(
               new RoundingsatFactory<REAL>( ));
      }
   };

} // namespace papilo

#endif
