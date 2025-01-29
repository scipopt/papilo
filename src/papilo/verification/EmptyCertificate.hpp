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

#ifndef _PAPILO_VERI_EMPTY_CERTIFICATE_HPP_
#define _PAPILO_VERI_EMPTY_CERTIFICATE_HPP_

#include "papilo/verification/CertificateInterface.hpp"

namespace papilo
{

/// type to store necessary data for post solve
template <typename REAL>
class EmptyCertificate : public CertificateInterface<REAL>
{

 private:
   Vec<int> vec{};

 public:
   EmptyCertificate() = default;

   void
   print_header(){};

   void
   start_transaction(){};

   void
   end_transaction( const Problem<REAL>& problem,
                    const Vec<int>& var_mapping, const Vec<int>& dirty_row_states ){};

   void
   flush(){};

   const Vec<int>&
   getRowScalingFactor()
   {
      return vec;
   }

   void
   change_upper_bound( REAL val, int col, const Problem<REAL>& problem,
                       const Vec<int>& var_mapping, MatrixBuffer<REAL>& matrix_buffer,
                       ArgumentType argument = ArgumentType::kPrimal )
   {
   }

   void
   change_lower_bound( REAL val, int col, const Problem<REAL>& problem,
                       const Vec<int>& var_mapping, MatrixBuffer<REAL>& matrix_buffer,
                       ArgumentType argument = ArgumentType::kPrimal )
   {
   }

   void
   dominating_columns( int dominating_column, int dominated_column,
                       const Vec<String>& names, const Vec<int>& var_mapping )
   {
   }


   void
   add_probing_reasoning( bool is_upper, int causing_col, int col,
                          const Vec<String>& names,
                          const Vec<int>& var_mapping) {}
   void
   change_rhs( int row, REAL val, const SparseVectorView<REAL>& data,
               const Vec<String>& names, const Vec<int>& var_mapping,
               ArgumentType argument = ArgumentType::kPrimal )
   {
   }

   void
   change_lhs( int row, REAL val, const SparseVectorView<REAL>& data,
               const Vec<String>& names, const Vec<int>& var_mapping,
               ArgumentType argument = ArgumentType::kPrimal )
   {
   }

   void
   change_rhs_parallel_row( int row, REAL val, int parallel_row,
                            const Problem<REAL>& problem,
                            const Vec<int>& var_mapping )
   {
   }

   void
   change_lhs_parallel_row( int row, REAL val, int parallel_row,
                            const Problem<REAL>& problem )
   {
   }

   void
   store_gcd( int row, REAL gcd ) {};

   void
   store_parallel_row( int row ) {};

   void
   store_implied_bound( int row, REAL lowerbound ) {};

   void
   change_lhs_inf( int row ){}

   void
   change_rhs_inf( int row ){}

   void
   log_forcing_row ( int row ) {}

   void
   mark_row_redundant( int row, const Problem<REAL>& currentProblem, ArgumentType argument = ArgumentType::kPrimal  )
   {
   }

   void
   change_matrix_entry( int row, int col, REAL new_val,
                        const SparseVectorView<REAL>& data, RowFlags& rflags,
                        REAL lhs, REAL rhs, const Vec<String>& names,
                        const Vec<int>& var_mapping, bool is_next_reduction_matrix_entry, ArgumentType argument ){};

   void
   substitute( int col, int row, REAL old_obj_coeff, const Problem<REAL>& currentProblem, const Vec<int>& var_mapping, ArgumentType argument  ){};

   void
   substitute_col_singleton_implied( int col, int row, REAL old_obj_coeff, const Problem<REAL>& currentProblem, const Vec<int>& var_mapping  ) {};

   void
   substitute( int col, const SparseVectorView<REAL>& equality, REAL offset, REAL old_obj_coeff,
               const Problem<REAL>& currentProblem, const Vec<String>& names,
               const Vec<int>& var_mapping )   {
   }

   void
   sparsify( int eqrow, int candrow, REAL scale,
             const Problem<REAL>& currentProblem )
   {
   }

   void
   symmetries(
       const SymmetryStorage& symmetries, const Vec<String>& names,
       const Vec<int>& var_mapping ) {};

   void
   log_solution( const Solution<REAL>& orig_solution,
                 const Vec<String>& names, REAL origobj ){};

   void
   setInfeasibleCause(int col){};

   void
   compress( const Vec<int>& rowmapping, const Vec<int>& colmapping,
             bool full = false )
   {
   }

};

} // namespace papilo

#endif
