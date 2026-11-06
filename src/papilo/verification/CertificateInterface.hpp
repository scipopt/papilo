/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*               This file is part of the program and library                */
/*    PaPILO --- Parallel Presolve for Integer and Linear Optimization       */
/*                                                                           */
/* Copyright (C) 2020-2026 Zuse Institute Berlin (ZIB)                       */
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

#ifndef _PAPILO_VERI_CERTIFICATE_INTERFACE_HPP_
#define _PAPILO_VERI_CERTIFICATE_INTERFACE_HPP_

#include "papilo/core/Problem.hpp"
#include "papilo/misc/Num.hpp"
#include "papilo/misc/Vec.hpp"
#include "papilo/misc/compress_vector.hpp"
#include "papilo/misc/fmt.hpp"
#include "papilo/verification/ArgumentType.hpp"

namespace papilo
{

/// type to store necessary data for post solve
template <typename REAL>
class CertificateInterface
{

 public:
   CertificateInterface() = default;

   virtual void
   start_transaction() = 0;

   virtual void
   end_transaction( const Problem<REAL>& problem,
                    const Vec<int>& var_mapping, const Vec<int>& dirty_row_states ) = 0;

   virtual void
   print_header() = 0;

   virtual const Vec<int>&
   getRowScalingFactor() = 0;

   virtual void
   flush() = 0;

   virtual void
   change_upper_bound( REAL val, int col, const Problem<REAL>& problem,
                       const Vec<int>& var_mapping, MatrixBuffer<REAL>& matrix_buffer,
                       ArgumentType argument = ArgumentType::kPrimal ) = 0;

   virtual void
   change_lower_bound( REAL val, int col, const Problem<REAL>& problem,
                       const Vec<int>& var_mapping, MatrixBuffer<REAL>& matrix_buffer,
                       ArgumentType argument = ArgumentType::kPrimal ) = 0;

   virtual void
   dominating_columns( int dominating_column, int dominated_column,
                       const Vec<String>& names,
                       const Vec<int>& var_mapping) = 0;

   virtual void
   add_probing_reasoning( bool is_upper, int causing_col, int col,
                       const Vec<String>& names,
                       const Vec<int>& var_mapping) = 0;

   virtual void
   change_rhs( int row, REAL val, const SparseVectorView<REAL>& data,
               const Vec<String>& names, const Vec<int>& var_mapping,
               ArgumentType argument = ArgumentType::kPrimal ) = 0;

   virtual void
   change_lhs( int row, REAL val, const SparseVectorView<REAL>& data,
               const Vec<String>& names, const Vec<int>& var_mapping,
               ArgumentType argument = ArgumentType::kPrimal ) = 0;

   virtual void
   store_gcd( int row, REAL gcd ) = 0;

   virtual void
   store_parallel_row( int row ) = 0;

   virtual void
   store_implied_bound( int row, REAL lowerbound ) = 0;

   virtual void
   change_rhs_parallel_row( int row, REAL val, int parallel_row,
                            const Problem<REAL>& problem,
                            const Vec<int>& var_mapping ) = 0;

   virtual void
   change_lhs_parallel_row( int row, REAL val, int parallel_row,
                            const Problem<REAL>& problem ) = 0;

   virtual void
   change_lhs_inf( int row ) = 0;

   virtual void
   change_rhs_inf( int row ) = 0;

   virtual void
   mark_row_redundant( int row, const Problem<REAL>& currentProblem, ArgumentType argument = ArgumentType::kPrimal  ) = 0;

   virtual void
   log_forcing_row ( int row ) = 0;

   virtual void
   change_matrix_entry( int row, int col, REAL new_val,
                        const SparseVectorView<REAL>& data, RowFlags& rflags,
                        REAL lhs, REAL rhs, const Vec<String>& names,
                        const Vec<int>& var_mapping, bool is_next_reduction_matrix_entry,
                        ArgumentType argument ) = 0;

   virtual void
   substitute( int col, int row, REAL old_obj_coeff, const Problem<REAL>& currentProblem, const Vec<int>& var_mapping, ArgumentType argument  ) = 0;

   virtual void
   substitute_col_singleton_implied( int col, int row, REAL old_obj_coeff, const Problem<REAL>& currentProblem, const Vec<int>& var_mapping  ) = 0;

   virtual void
   substitute( int col, const SparseVectorView<REAL>& equality, REAL offset, REAL old_obj_coeff,
               const Problem<REAL>& currentProblem, const Vec<String>& names,
               const Vec<int>& var_mapping ) = 0;

   virtual void
   sparsify( int eqrow, int candrow, REAL scale,
             const Problem<REAL>& currentProblem ) = 0;

   virtual void
   log_solution( const Solution<REAL>& orig_solution,
                 const Vec<String>& names, REAL origobj ) = 0;

   virtual void
   symmetries(
       const SymmetryStorage& symmetries, const Vec<String>& names,
       const Vec<int>& var_mapping ) = 0;

   virtual void
   compress( const Vec<int>& rowmapping, const Vec<int>& colmapping,
             bool full = false ) = 0;

   virtual void
   setInfeasibleCause(int col) = 0;

   virtual void
   infeasible( ) { };

   virtual void
   end_proof( ) { };

   virtual void
   infeasible( const Vec<int>& colmapping, const Vec<String>& names ){ };

   virtual ~CertificateInterface() = default;
};

} // namespace papilo

#endif
