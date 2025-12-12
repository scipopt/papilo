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

#ifndef _PAPILO_PRESOLVERS_GF2_HPP_
#define _PAPILO_PRESOLVERS_GF2_HPP_

#include "papilo/core/PresolveMethod.hpp"
#include "papilo/core/Problem.hpp"
#include "papilo/core/ProblemUpdate.hpp"

namespace papilo
{

#if GF2_PRESOLVE_DEBUG
#define NOT_GF2( reason, ... )                                                 \
   do                                                                          \
   {                                                                           \
      printf( "NO : Cons %d is not gf2: " reason "\n", cstr_idx,               \
              ##__VA_ARGS__ );                                                 \
      goto not_valid;                                                          \
   } while( 0 )
#else
#define NOT_GF2( reason, ... )                                                 \
   do                                                                          \
   {                                                                           \
      goto not_valid;                                                          \
   } while( 0 )
#endif

/// presolver to fix continuous variables whose bounds are very close
template <typename REAL>
class GF2 final : public PresolveMethod<REAL>
{
 public:
   GF2() : PresolveMethod<REAL>()
   {
      this->setName( "gf2" );
      this->setTiming( PresolverTiming::kMedium );
      this->setType( PresolverType::kIntegralCols );
   }

   bool
   check_variables_and_coeff_in_constraint(
       const Vec<ColFlags>& col_flags, const Vec<REAL>& lower_bounds,
       const Vec<REAL>& upper_bounds,
       std::unordered_map<size_t, size_t> gf2_bin_vars,
       std::unordered_map<size_t, size_t> gf2_key_vars,
       REAL integrality_tolerance, int& key_var_idx, REAL& key_var_coeff,
       std::vector<std::pair<size_t, REAL>> constraint_bin_vars,
       const int* row_indices, const REAL* row_values, int row_length );
   PresolveStatus
   execute( const Problem<REAL>& problem,
            const ProblemUpdate<REAL>& problemUpdate, const Num<REAL>& num,
            Reductions<REAL>& reductions, const Timer& timer,
            int& reason_of_infeasibility ) override;

   struct gf2_constraint_t
   {
      size_t cstr_idx;
      Vec<std::pair<size_t, REAL>> bin_vars;
      std::pair<size_t, REAL> key_var;
      size_t rhs; // 0 or 1

      gf2_constraint_t() = default;
      gf2_constraint_t( size_t cstr_idx,
                        std::vector<std::pair<size_t, REAL>> bin_vars,
                        std::pair<size_t, REAL> key_var, size_t rhs )
          : cstr_idx( cstr_idx ), bin_vars( std::move( bin_vars ) ),
            key_var( key_var ), rhs( rhs )
      {
      }
      gf2_constraint_t( const gf2_constraint_t& other ) = default;
      gf2_constraint_t( gf2_constraint_t&& other ) noexcept = default;
      gf2_constraint_t&
      operator=( const gf2_constraint_t& other ) = default;
      gf2_constraint_t&
      operator=( gf2_constraint_t&& other ) noexcept = default;
   };

};

#ifdef PAPILO_USE_EXTERN_TEMPLATES
extern template class GF2<double>;
extern template class GF2<Quad>;
extern template class GF2<Rational>;
#endif

template <typename i_t>
static inline i_t
positive_modulo( i_t i, i_t n )
{
   return ( i % n + n ) % n;
}

// this is kind-of a stopgap implementation (as in practice MIPLIB2017 only
// contains a couple of GF2 problems and they're small) but cuDSS could be used
// for this since A is likely to be sparse and low-bandwidth (i think?) unlikely
// to occur in real-world problems however. doubt it'd be worth the effort
// trashes A and b, return true if solved
static bool
gf2_solve( std::vector<std::vector<int>>& A, std::vector<int>& b,
           std::vector<int>& x )
{
   int i, j, k;
   const int N = A.size();
   for( i = 0; i < N; i++ )
   {
      // Find pivot
      int pivot = -1;
      for( j = i; j < N; j++ )
      {
         if( A[j][i] )
         {
            pivot = j;
            break;
         }
      }
      if( pivot == -1 )
         return false; // No solution

      // Swap current row with pivot row if needed
      if( pivot != i )
      {
         for( k = 0; k < N; k++ )
         {
            const int temp = A[i][k];
            A[i][k] = A[pivot][k];
            A[pivot][k] = temp;
         }
         const int temp = b[i];
         b[i] = b[pivot];
         b[pivot] = temp;
      }

      // Eliminate downwards
      for( j = i + 1; j < N; j++ )
      {
         if( A[j][i] )
         {
            for( k = i; k < N; k++ )
               A[j][k] ^= A[i][k];
            b[j] ^= b[i];
         }
      }
   }

   // Back-substitution
   for( i = N - 1; i >= 0; i-- )
   {
      x[i] = b[i];
      for( j = i + 1; j < N; j++ )
         x[i] ^= ( A[i][j] & x[j] );
      if( !A[i][i] && x[i] )
         return false; // No solution
   }
   return true; // Success
}

template <typename REAL>
bool
GF2<REAL>::check_variables_and_coeff_in_constraint(
    const Vec<ColFlags>& col_flags, const Vec<REAL>& lower_bounds,
    const Vec<REAL>& upper_bounds,
    std::unordered_map<size_t, size_t> gf2_bin_vars,
    std::unordered_map<size_t, size_t> gf2_key_vars, REAL integrality_tolerance,
    int& key_var_idx, REAL& key_var_coeff,
    std::vector<std::pair<size_t, REAL>> constraint_bin_vars,
    const int* row_indices, const REAL* row_values, int row_length )
{
   // TODO:
   Num<REAL> num{};
   for( int j = 0; j < row_length; ++j )
   {
      if( !num.isIntegral( row_values[j] ) )
         return false;

      int var_idx = row_indices[j];

      REAL coeff = num.round( row_values[j] );

      // Check if variable is integer
      if( !col_flags[var_idx].test( ColFlag::kIntegral ) )
         return false;

      bool is_binary = col_flags[var_idx].test( ColFlag::kLbInf ) ? false
                       : col_flags[var_idx].test( ColFlag::kUbInf )
                           ? false
                           : ( lower_bounds[var_idx] == 0.0 &&
                               upper_bounds[var_idx] == 1.0 );

      bool abs2 = coeff == 2 || coeff == -2;
      bool abs1 = coeff == 1 || coeff == -1;
      // Check coefficient constraints
      if( is_binary && !abs1 && !abs2 )
         return false;
      if( !is_binary && !abs2 )
         return false;

      // Key variable (coefficient of 2)
      if( abs2 )
      {
         if( key_var_idx != -1 )
            return false;

         key_var_idx = var_idx;
         key_var_coeff = coeff;
         gf2_key_vars.insert( { var_idx, gf2_key_vars.size() } );
      }
      else
      {
         // Binary variable
         constraint_bin_vars.push_back( { var_idx, coeff } );
         gf2_bin_vars.insert( { var_idx, gf2_bin_vars.size() } );
      }
   }

   if( key_var_idx == -1 )
      return false;
   return true;
}

template <typename REAL>
PresolveStatus
GF2<REAL>::execute( const Problem<REAL>& problem,
                    const ProblemUpdate<REAL>& problemUpdate,
                    const Num<REAL>& num, Reductions<REAL>& reductions,
                    const Timer& timer, int& reason_of_infeasibility )
{
   const auto& constraint_matrix = problem.getConstraintMatrix();
   const auto& lhs_values = constraint_matrix.getLeftHandSides();
   const auto& rhs_values = constraint_matrix.getRightHandSides();
   const auto& row_flags = constraint_matrix.getRowFlags();
   const auto& domains = problem.getVariableDomains();
   const auto& col_flags = domains.flags;
   const auto& lower_bounds = domains.lower_bounds;
   const auto& upper_bounds = domains.upper_bounds;

   const int num_rows = constraint_matrix.getNRows();

   std::unordered_map<size_t, size_t> gf2_bin_vars;
   std::unordered_map<size_t, size_t> gf2_key_vars;
   std::vector<gf2_constraint_t> gf2_constraints;

   const REAL integrality_tolerance = num.getFeasTol();

   for( int cstr_idx = 0; cstr_idx < num_rows; ++cstr_idx )
   {
      int key_var_idx = -1;
      REAL key_var_coeff = 0.0;

      std::vector<std::pair<size_t, REAL>> constraint_bin_vars;

      // Check constraint coefficients
      auto row_coeff = constraint_matrix.getRowCoefficients( cstr_idx );
      const int* row_indices = row_coeff.getIndices();
      const REAL* row_values = row_coeff.getValues();
      const int row_length = row_coeff.getLength();
      REAL rhs = num.round( lhs_values[cstr_idx] );

      // Check if this is an equality constraint
      if( !problem.getRowFlags()[cstr_idx].test( RowFlag::kEquation ) )
         continue;

      // Only accept 0, 1, -1 as rhs
      if( rhs != 0.0 && rhs != 1.0 && rhs != -1.0 )
         continue;

      if( !check_variables_and_coeff_in_constraint(
              col_flags, lower_bounds, upper_bounds, gf2_bin_vars, gf2_key_vars,
              integrality_tolerance, key_var_idx, key_var_coeff,
              constraint_bin_vars, row_indices, row_values, row_length ) )
         continue;

      gf2_constraints.emplace_back(
          static_cast<size_t>( cstr_idx ), std::move( constraint_bin_vars ),
          std::pair<size_t, REAL>{ key_var_idx, key_var_coeff },
          positive_modulo( static_cast<int>( rhs ), 2 ) );
   }

   // If no GF2 constraints found, return unchanged
   if( gf2_constraints.empty() )
      return PresolveStatus::kUnchanged;

   // Skip if that would cause computational explosion (O(n^3) with simple
   // gaussian elimination)
   if( gf2_constraints.size() > 1000 )
      return PresolveStatus::kUnchanged;

   // Validate structure
   if( gf2_key_vars.size() != gf2_constraints.size() ||
       gf2_bin_vars.size() != gf2_constraints.size() )
      return PresolveStatus::kUnchanged;

   // Create inverse mappings
   std::unordered_map<size_t, size_t> gf2_bin_vars_invmap;
   for( const auto& [var_idx, gf2_idx] : gf2_bin_vars )
      gf2_bin_vars_invmap.insert( { gf2_idx, var_idx } );

   // Build binary matrix
   // Could be a flat vector but. oh well. in practice N is small
   std::vector<std::vector<int>> A(
       gf2_constraints.size(), std::vector<int>( gf2_constraints.size(), 0 ) );
   std::vector<int> b( gf2_constraints.size() );
   for( const auto& cons : gf2_constraints )
   {
      for( auto [bin_var, _] : cons.bin_vars )
      {
         A[cons.cstr_idx][gf2_bin_vars[bin_var]] = 1;
      }
      b[cons.cstr_idx] = cons.rhs;
   }

   std::vector<int> solution( gf2_constraints.size() );
   bool feasible = gf2_solve( A, b, solution );
   if( !feasible )
   {
      return PresolveStatus::kInfeasible;
   }

   std::unordered_map<size_t, REAL> fixings;
   // Fix binary variables
   for( size_t sol_idx = 0; sol_idx < gf2_constraints.size(); ++sol_idx )
      fixings[gf2_bin_vars_invmap[sol_idx]] = solution[sol_idx];

   // Compute fixings for key variables by solving for the constraint
   for( const auto& cons : gf2_constraints )
   {
      auto [key_var_idx, key_var_coeff] = cons.key_var;
      REAL constraint_rhs = lhs_values[cons.cstr_idx]; // equality constraint
      REAL lhs = -constraint_rhs;
      for( auto [bin_var, coeff] : cons.bin_vars )
      {
         lhs += fixings[bin_var] * coeff;
      }
      fixings[key_var_idx] = num.round( -lhs / key_var_coeff );
   }

   auto status = PresolveStatus::kUnchanged;
   TransactionGuard<REAL> rg{ reductions };
   for( const auto& [var_idx, fixing] : fixings )
   {
      if( num.isZero( fixing ) )
         reductions.fixCol( var_idx, 0 );
      else
         reductions.fixCol( var_idx, fixing );
      status = PresolveStatus::kReduced;
   }

   return status;
}

} // namespace papilo

#endif
