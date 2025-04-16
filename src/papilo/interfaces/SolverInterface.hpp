/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*               This file is part of the program and library                */
/*    PaPILO --- Parallel Presolve for Integer and Linear Optimization       */
/*                                                                           */
/* Copyright (C) 2020-2025 Zuse Institute Berlin (ZIB)                       */
/*                                                                           */
/* Licensed under the Apache License, Version 2.0 (the "License");           */
/* you may not use this file except in compliance with the License.          */
/* You may obtain a copy of the License at                                   */
/*                                                                           */
/*     http://www.apache.org/licenses/LICENSE-2.0                            */
/*                                                                           */
/* Unless required by applicable law or agreed to in writing, software       */
/* distributed under the License is distributed on an "AS IS" BASIS,         */
/* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  */
/* See the License for the specific language governing permissions and       */
/* limitations under the License.                                            */
/*                                                                           */
/* You should have received a copy of the Apache-2.0 license                 */
/* along with PaPILO; see the file LICENSE. If not visit scipopt.org.        */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef _PAPILO_INTERFACES_SOLVER_INTERFACE_HPP_
#define _PAPILO_INTERFACES_SOLVER_INTERFACE_HPP_

#include "papilo/core/Components.hpp"
#include "papilo/io/Message.hpp"
#include "papilo/misc/ParameterSet.hpp"
#include "papilo/misc/String.hpp"
#include "papilo/misc/Vec.hpp"
#include <memory>
#include <string>

namespace papilo
{

enum class SolverType : int
{
   LP,
   MIP,
   PseudoBoolean
};

enum class SolverStatus : int
{
   kInit,
   kOptimal,
   kInfeasible,
   kUnbounded,
   kUnbndOrInfeas,
   kInterrupted,
   kError
};

template <typename REAL>
class SolverInterface
{
 protected:
   SolverStatus status;

 public:
   SolverInterface() : status( SolverStatus::kInit ) {}

   virtual void
   setUp( const Problem<REAL>& prob, const Vec<int>& row_maps,
          const Vec<int>& col_maps ) = 0;

   virtual void
   setUp( const Problem<REAL>& prob, const Vec<int>& row_maps,
          const Vec<int>& col_maps, const Components& components,
          const ComponentInfo& component ) = 0;

   virtual void
   solve() = 0;

   virtual SolverType
   getType() = 0;

   virtual String
   getName() = 0;

   virtual void
   printDetails()
   {
   }

   virtual void
   readSettings( const String& file )
   {
   }

   SolverStatus
   getStatus()
   {
      return status;
   }

   virtual void
   setNodeLimit( int num )
   {
   }

   virtual void
   setGapLimit( const REAL& gaplim )
   {
   }

   virtual void
   setSoftTimeLimit( double tlim )
   {
   }

   virtual void
   setTimeLimit( double tlim ) = 0;

   virtual void
   setVerbosity( VerbosityLevel verbosity ) = 0;

   virtual void
   setRowScalingFactor( Vec<int> scaling_factor ) {}

   virtual bool
   getSolution( Solution<REAL>& sol, PostsolveStorage<REAL>& postsolve ) = 0;

   virtual bool
   getSolution( const Components& components, int component,
                Solution<REAL>& sol ) = 0;

   virtual REAL
   getDualBound() = 0;

   virtual bool
   is_dual_solution_available() = 0;

   virtual void
   addParameters( ParameterSet& paramSet )
   {
   }

   virtual ~SolverInterface() {}
};

template <typename REAL>
class SolverFactory
{
 public:
   virtual std::unique_ptr<SolverInterface<REAL>>
   newSolver( VerbosityLevel verbosity = VerbosityLevel::kQuiet ) const = 0;

   virtual void
   add_parameters( ParameterSet& parameter ) const = 0;

   virtual ~SolverFactory() {}
};

} // namespace papilo

#endif
