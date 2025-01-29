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

#ifndef PAPILO_VERI_VERI_PB_HPP_
#define PAPILO_VERI_VERI_PB_HPP_

//turn on to enable some comments and extended checks
//#define VERIPB_DEBUG

#define VERIPB_VERSION 2

#include "papilo/Config.hpp"
#include "papilo/core/Problem.hpp"
#include "papilo/misc/Hash.hpp"
#include "papilo/misc/Num.hpp"
#include "papilo/misc/Vec.hpp"
#include "papilo/misc/compress_vector.hpp"
#include "papilo/misc/fmt.hpp"
#include "papilo/verification/ArgumentType.hpp"
#include "papilo/verification/CertificateInterface.hpp"

namespace papilo
{

static const char* const NEGATED = "~";
static const char* const COMMENT = "* ";
static const char* const CONCLUSION = "conclusion ";
static const char* const OUTPUT = "output ";
static const char* const NONE = "NONE";
#if VERIPB_VERSION >= 2
static const char* const DELETE_CONS = "delc ";
#else
static const char* const DELETE_CONS = "del id ";
#endif
static const char* const OBJECTIVE_DIFF = "obju diff ";
static const char* const POL = "pol ";
static const char* const RUP = "rup ";
static const char* const RED = "red ";
static const char* const EQUAL_CHECK = "e ";
static const char* const SATURATION = "s";
static const char* const WEAKENING = "w";
static const int UNKNOWN = -1;
static const char* const MOVE_LAST_CONS_TO_CORE = "core id -1\n";

template <typename REAL>
class VeriPb : public CertificateInterface<REAL>
{
 public:

   Num<REAL> num;
   std::ofstream proof_out;

   int propagation_option;
   int status = 0; // 1 = solved, -1 = infeasible -2 = finished;

   Objective<REAL> stored_objective;


   /// mapping constraint from PaPILO to constraint ids from PaPILO to VeriPb
   Vec<int> rhs_row_mapping;
   Vec<int> lhs_row_mapping;

   /// information required for optimization problems
   bool is_optimization_problem = false;
#if VERIPB_VERSION == 1
   bool verification_possible = true;
#endif
   HashMap<int, Vec<int>> substitutions;
   int cause = -1;
   int row_forcing_propagation = -1;

   ////
   int stored_dominating_col = UNKNOWN;
   int stored_dominated_col = UNKNOWN;
   Vec<int> weakened_columns{};
   std::pair<int, int> row_with_gcd = {UNKNOWN, UNKNOWN};
   int row_implying_lb = UNKNOWN;
   int row_implying_ub = UNKNOWN;
   int parallel_remaining_row = UNKNOWN;

   /// PaPILO does not care about the integrality of the coefficient
   /// therefore store scale factors to ensure the integrality
   Vec<int> scale_factor;


   Vec<int> fixed_variable;


   /// this holds the id of the next generated constraint of VeriPB
   int next_constraint_id = 0;

   /// skip reductions
   int skip_deleting_rhs_constraint_id = UNKNOWN;
   int skip_deleting_lhs_constraint_id = UNKNOWN;
   int skip_changing_rhs = UNKNOWN;
   int skip_changing_lhs = UNKNOWN;
   bool saturation_already_called = false;



#ifdef VERIPB_DEBUG
   int validate_row = UNKNOWN;
#endif

   /// this stores the changed entries for the current row (otherwise refer to
   /// the matrix_buffer)
   HashMap<int, int> changed_entries_during_current_tsxs{};

   const Vec<int>&
   getRowScalingFactor() override
   {
      return scale_factor;
   }

   VeriPb( const Problem<REAL>& _problem, const Num<REAL>& _num, PresolveOptions options ) :
         num( _num ), propagation_option(options.veripb_propagation_option), substitutions( {})
   {
      rhs_row_mapping.reserve( _problem.getNRows() );
      lhs_row_mapping.reserve( _problem.getNRows() );
      scale_factor.reserve( _problem.getNRows() );

      fixed_variable.reserve(_problem.getNCols());
      std::fill_n(std::back_inserter(fixed_variable), _problem.getNCols(), 0);


      stored_objective = Objective<REAL>(_problem.getObjective());
      auto coefficients = stored_objective.coefficients;

      for( int i = 0; i < _problem.getNRows(); ++i )
      {
         scale_factor.push_back( 1 );
         if( !_problem.getRowFlags()[i].test( RowFlag::kLhsInf ) )
         {
            next_constraint_id++;
            lhs_row_mapping.push_back( next_constraint_id );
         }
         else
            lhs_row_mapping.push_back( UNKNOWN );
         if( !_problem.getRowFlags()[i].test( RowFlag::kRhsInf ) )
         {
            next_constraint_id++;
            rhs_row_mapping.push_back( next_constraint_id );
         }
         else
            rhs_row_mapping.push_back( UNKNOWN );
      }
      for( int i = 0; i < _problem.getNCols(); ++i )
      {
         if( coefficients[i] != 0 )
         {
            is_optimization_problem = true;
#if VERIPB_VERSION == 1
            verification_possible = true;
            fmt::print("VeriPB 1.0 does not support objective update rule. Hence verification not entirely supported for optimization problems.\n");
#endif
            break;
         }
      }

      assert( rhs_row_mapping.size() == lhs_row_mapping.size() );
      assert( static_cast<int>(rhs_row_mapping.size()) == _problem.getNRows() );

      auto problem_name = _problem.getName();
      int length = problem_name.length();
      int ending = 4;
#ifdef PAPILO_USE_BOOST_IOSTREAMS_WITH_ZLIB
      if( problem_name.substr( length - 3 ) == ".gz" )
         ending = 7;
#endif
#ifdef PAPILO_USE_BOOST_IOSTREAMS_WITH_BZIP2
      if( problem_name.substr( length - 4 ) == ".bz2" )
         ending = 8;
#endif
      proof_out =
          std::ofstream( problem_name.substr( 0, length - ending ) + ".pbp" );
   }

   void
   print_header() override
   {
#if VERIPB_VERSION == 1
      proof_out << "pseudo-Boolean proof version 1.0\n";
#else
      proof_out << "pseudo-Boolean proof version 2.0\n";
#endif
      proof_out << COMMENT << "Log files generated by PaPILO "
                << PAPILO_VERSION_MAJOR << "." << PAPILO_VERSION_MINOR << "."
                << PAPILO_VERSION_PATCH;
#ifdef PAPILO_GITHASH_AVAILABLE
      proof_out << " [GitHash: " << PAPILO_GITHASH << " ]";
#endif
      proof_out << "\n";
      proof_out << "f " << next_constraint_id << "\n";
#if VERIPB_VERSION == 1
      if(!verification_possible)
         proof_out << "* Verification currently not possible for optimization problems!\n";
#endif
      proof_out << std::fixed;
   };

   void
   start_transaction() override
   {
      skip_changing_lhs = UNKNOWN;
      skip_changing_rhs = UNKNOWN;
      skip_deleting_lhs_constraint_id = UNKNOWN;
      skip_deleting_rhs_constraint_id = UNKNOWN;
      stored_dominating_col = UNKNOWN;
      stored_dominated_col = UNKNOWN;
      changed_entries_during_current_tsxs.clear();
      saturation_already_called = false;
      row_with_gcd = {UNKNOWN, UNKNOWN};
      row_implying_lb = UNKNOWN;
      weakened_columns.clear();
      row_implying_ub = UNKNOWN;
      parallel_remaining_row = UNKNOWN;
#ifdef VERIPB_DEBUG
      validate_row = UNKNOWN;
#endif
   };

   void
   end_transaction( const Problem<REAL>& problem,
                    const Vec<int>& var_mapping, const Vec<int>& dirty_row_states ) override
   {
      if( row_with_gcd.first != UNKNOWN )
      {
         int row = row_with_gcd.first;
         assert( ( lhs_row_mapping[row] == UNKNOWN &&
                   rhs_row_mapping[row] != UNKNOWN ) ||
                 ( rhs_row_mapping[row] == UNKNOWN &&
                   lhs_row_mapping[row] != UNKNOWN ) );
         if( lhs_row_mapping[row] != UNKNOWN )
            change_lhs( row, row_with_gcd.second,
                        problem.getConstraintMatrix().getRowCoefficients( row ),
                        problem.getVariableNames(), var_mapping,
                        ArgumentType::kWeakening );
         else
            change_rhs( row, row_with_gcd.second,
                        problem.getConstraintMatrix().getRowCoefficients( row ),
                        problem.getVariableNames(), var_mapping,
                        ArgumentType::kWeakening );
      }
      assert(weakened_columns.empty());
#ifdef VERIPB_DEBUG
      if( validate_row != UNKNOWN )
      {
         verify_changed_row( validate_row, problem, var_mapping, dirty_row_states );
         validate_row = UNKNOWN;
      }
      proof_out << "*end trans " << next_constraint_id << "\n";
#endif
   };

   void
   flush() override
   {
      proof_out.flush();
   };

   void
   change_upper_bound( REAL val, int col, const Problem<REAL>& problem,
                       const Vec<int>& var_mapping, MatrixBuffer<REAL>& matrix_buffer,
                       ArgumentType argument = ArgumentType::kPrimal ) override
   {
#if VERIPB_VERSION == 1
      if( !verification_possible )
         return;
#endif
      next_constraint_id++;
      assert( val == 0 );
      const Vec<String>& names = problem.getVariableNames();
      int orig_col = var_mapping[col];
      switch( argument )
      {
      case ArgumentType::kPropagation:
         if(propagation_option == 1)
         {
            assert( row_forcing_propagation != -1 );
            propagate_row( row_forcing_propagation, col, val, false, problem,
                           var_mapping );
            break;
         }
         proof_out << RUP << "1 " << NEGATED << names[orig_col] << " >= 1 ;\n";
         break;
      case ArgumentType::kPrimal:
         if( stored_dominated_col == orig_col)
         {
            assert(stored_dominating_col != UNKNOWN);
            proof_out << RED << "1 " << NEGATED << names[orig_col] << " >= 1 ; "
                      << names[orig_col] << " -> 0 " << names[stored_dominating_col] << " -> 0";
#if VERIPB_VERSION == 1
            add_substitutions_fix_to_witness( names, stored_dominated_col, val == 1 );
            add_substitutions_fix_to_witness( names, stored_dominating_col, val == 1 );
#endif
            proof_out << "\n";
            break;
         }

         proof_out << RUP << "1 " << NEGATED << names[orig_col] << " >= 1 ;\n";
         break;
      case ArgumentType::kAggregation:
      case ArgumentType::kDual:
      case ArgumentType::kSymmetry:
         proof_out << RED << "1 " << NEGATED << names[orig_col] << " >= 1 ; "
                   << names[orig_col] << " -> 0";
#if VERIPB_VERSION == 1
         add_substitutions_fix_to_witness( names, var_mapping[col], val == 1 );
#endif
         proof_out << "\n";
         break;
      default:
         assert( false );
         return;
      }
#if VERIPB_VERSION == 2
      proof_out << MOVE_LAST_CONS_TO_CORE;
#endif
      substitutions.erase(var_mapping[col]);
      int cons_id_fixing = next_constraint_id;
      auto col_coeff =
          problem.getConstraintMatrix().getColumnCoefficients( col );
      for( int row_index = 0; row_index < col_coeff.getLength(); row_index++ )
      {
         int row = col_coeff.getIndices()[row_index];
         if( problem.getRowFlags()[row].test( RowFlag::kRedundant ) )
         {
            assert( rhs_row_mapping[row] == UNKNOWN );
            assert( lhs_row_mapping[row] == UNKNOWN );
            continue;
         }
         assert( num.isIntegral( col_coeff.getValues()[row_index] *
                                 scale_factor[row] ) );
         REAL unscaled_row_value;
         // check if the matrix coefficients were updated since the last call
         const MatrixEntry<REAL>* entry = matrix_buffer.template findEntry<false>( row, col);
         if( entry != nullptr )
            unscaled_row_value = entry->val;
         else
            unscaled_row_value = col_coeff.getValues()[row_index];
         assert(num.isIntegral(unscaled_row_value * scale_factor[row]));
         int row_value =
             cast_to_long( unscaled_row_value * scale_factor[row] );
         if( !problem.getRowFlags()[row].test( RowFlag::kLhsInf ) )
         {
            next_constraint_id++;
            assert( lhs_row_mapping[row] != UNKNOWN );
            if( row_value < 0 )
               proof_out << POL << lhs_row_mapping[row] << " "
                         << names[orig_col] << " " << abs( row_value )
                         << " * +\n";
            else
               proof_out << POL << lhs_row_mapping[row] << " " << cons_id_fixing
                         << " " << abs( row_value ) << " * +\n";
#if VERIPB_VERSION >= 2
            proof_out << MOVE_LAST_CONS_TO_CORE;
#endif
            proof_out << DELETE_CONS << lhs_row_mapping[row];
            lhs_row_mapping[row] = next_constraint_id;
#if VERIPB_VERSION >= 2
            if( problem.getConstraintMatrix().getRowCoefficients(row).getLength() >1 )
            {
               proof_out << " ; ; begin \n\t";
               if( row_value < 0 )
                  proof_out << POL << lhs_row_mapping[row] << " "
                            << cons_id_fixing << " " << abs( row_value )
                            << " * +\n";
               else
                  proof_out << POL << lhs_row_mapping[row] << " "
                            << names[orig_col] << " " << abs( row_value )
                            << " * +\n";
               proof_out << "end";
               next_constraint_id += 2;
            }
#endif
            proof_out << "\n";
         }
         if( !problem.getRowFlags()[row].test( RowFlag::kRhsInf ) )
         {
            next_constraint_id++;
            assert( rhs_row_mapping[row] != UNKNOWN );
            if( row_value > 0 )
               proof_out << POL << rhs_row_mapping[row] << " "
                         << names[orig_col] << " " << abs( row_value )
                         << " * +\n";
            else
               proof_out << POL << rhs_row_mapping[row] << " " << cons_id_fixing
                         << " " << abs( row_value ) << " * +\n";
#if VERIPB_VERSION >= 2
            proof_out << MOVE_LAST_CONS_TO_CORE;
#endif
            proof_out << DELETE_CONS << rhs_row_mapping[row];
            rhs_row_mapping[row] = next_constraint_id;
#if VERIPB_VERSION >= 2
            if( problem.getConstraintMatrix().getRowCoefficients(row).getLength() >1 )
            {
               proof_out << " ; ; begin \n\t";
               if( row_value > 0 )
                  proof_out << POL << rhs_row_mapping[row] << " "
                            << cons_id_fixing << " " << abs( row_value )
                            << " * +\n";
               else
                  proof_out << POL << rhs_row_mapping[row] << " "
                            << names[orig_col] << " " << abs( row_value )
                            << " * +\n";
               proof_out << "end";
               next_constraint_id += 2;
            }
#endif
            proof_out << "\n";
         }
      }

#if VERIPB_VERSION == 2
      const auto obj_coeff = cast_to_long(stored_objective.coefficients[col]);
      assert( obj_coeff == problem.getObjective().coefficients[col] );
      if( obj_coeff != 0)
         proof_out << OBJECTIVE_DIFF << (-obj_coeff ) << " " << names[orig_col] << ";\n";
      stored_objective.coefficients[col] = 0;
      fixed_variable[col] = -1;
#endif
   }

   void
   change_lower_bound( REAL val, int col, const Problem<REAL>& problem,
                       const Vec<int>& var_mapping, MatrixBuffer<REAL>& matrix_buffer,
                       ArgumentType argument = ArgumentType::kPrimal ) override
   {
#if VERIPB_VERSION == 1
      if( !verification_possible )
         return;
#endif
      next_constraint_id++;
      assert( val == 1 );
      const Vec<String>& names = problem.getVariableNames();
      int orig_col = var_mapping[col];
      switch( argument )
      {
      case ArgumentType::kPropagation:
         if( propagation_option == 1)
         {
            assert( row_forcing_propagation != -1 );
            propagate_row( row_forcing_propagation, col, val, true, problem,
                           var_mapping );
            break;
         }
         else
         {
            proof_out << RUP << "1 " << names[orig_col]
                      << " >= " << cast_to_long( val ) << " ;\n";
            break;
         }
      case ArgumentType::kPrimal:
         if( stored_dominating_col == orig_col)
         {
            assert(stored_dominated_col != UNKNOWN);
            proof_out << RED << "1 " << names[orig_col]
                      << " >= " << cast_to_long( val ) << " ; "
                      << names[orig_col] << " -> " << cast_to_long( val ) << " "
                      << names[stored_dominated_col] << " -> 1";
#if VERIPB_VERSION == 1
            add_substitutions_fix_to_witness( names, stored_dominated_col, val == 1 );
            add_substitutions_fix_to_witness( names, stored_dominating_col, val == 1 );
#endif
            proof_out << "\n";
            break;
         }
         proof_out << RUP << "1 " << names[orig_col]
                   << " >= " << cast_to_long( val ) << " ;\n";
         break;
      case ArgumentType::kAggregation:
      case ArgumentType::kDual:
      case ArgumentType::kSymmetry:
         proof_out << RED << "1 " << names[orig_col]
                   << " >= " << cast_to_long( val ) << " ; "
                   << names[orig_col] << " -> " << cast_to_long( val );
#if VERIPB_VERSION == 1
         add_substitutions_fix_to_witness( names, var_mapping[col], val == 1 );
#endif
         proof_out << "\n";
         break;
      default:
         assert( false );
         return;
      }
#if VERIPB_VERSION == 2
      proof_out << MOVE_LAST_CONS_TO_CORE;
#endif
      substitutions.erase(var_mapping[col]);
      int cons_id_fixing = next_constraint_id;
      auto col_coeff =
          problem.getConstraintMatrix().getColumnCoefficients( col );
      for( int row_index = 0; row_index < col_coeff.getLength(); row_index++ )
      {
         int row = col_coeff.getIndices()[row_index];
         if( problem.getRowFlags()[row].test( RowFlag::kRedundant ) )
         {
            assert( rhs_row_mapping[row] == UNKNOWN );
            assert( lhs_row_mapping[row] == UNKNOWN );
            continue;
         }
         assert( num.isIntegral( col_coeff.getValues()[row_index] *
                                 scale_factor[row] ) );
         int row_value = cast_to_long( col_coeff.getValues()[row_index] *
                                           scale_factor[row] );
         if( !problem.getRowFlags()[row].test( RowFlag::kLhsInf ) )
         {
            next_constraint_id++;
            assert( lhs_row_mapping[row] != UNKNOWN );
            if( row_value > 0 )
               proof_out << POL << lhs_row_mapping[row] << " " << NEGATED
                         << names[orig_col] << " " << abs( row_value )
                         << " * +\n";
            else
               proof_out << POL << lhs_row_mapping[row] << " " << cons_id_fixing
                         << " " << abs( row_value ) << " * +\n";
#if VERIPB_VERSION>= 2
            proof_out << MOVE_LAST_CONS_TO_CORE;
#endif
            proof_out << DELETE_CONS << lhs_row_mapping[row];
            lhs_row_mapping[row] = next_constraint_id;
#if VERIPB_VERSION >= 2
            if( problem.getConstraintMatrix().getRowCoefficients(row).getLength() >1 )
            {
               proof_out << " ; ; begin \n\t";
               if( row_value > 0 )
                  proof_out << POL << lhs_row_mapping[row] << " "
                            << cons_id_fixing << " " << abs( row_value )
                            << " * +\n";
               else
                  proof_out << POL << lhs_row_mapping[row] << " " << NEGATED
                            << names[orig_col] << " " << abs( row_value )
                            << " * +\n";
               proof_out << "end";
               next_constraint_id += 2;
            }
#endif
            proof_out << "\n";
         }
         if( !problem.getRowFlags()[row].test( RowFlag::kRhsInf ) )
         {
            next_constraint_id++;
            assert( rhs_row_mapping[row] != UNKNOWN );
            if( row_value < 0 )
               proof_out << POL << rhs_row_mapping[row] << " " << NEGATED
                         << names[orig_col] << " " << abs( row_value )
                         << " * +\n";
            else
               proof_out << POL << rhs_row_mapping[row] << " " << cons_id_fixing
                         << " " << abs( row_value ) << " * +\n";
#if VERIPB_VERSION>= 2
            proof_out << MOVE_LAST_CONS_TO_CORE;
#endif
            proof_out << DELETE_CONS << rhs_row_mapping[row];
            rhs_row_mapping[row] = next_constraint_id;
#if VERIPB_VERSION >= 2
            if( problem.getConstraintMatrix().getRowCoefficients(row).getLength() >1 )
            {
               proof_out << " ; ; begin \n\t";
               if( row_value < 0 )
                  proof_out << POL << rhs_row_mapping[row] << " "
                            << cons_id_fixing << " " << abs( row_value )
                            << " * +\n";
               else
                  proof_out << POL << rhs_row_mapping[row] << " " << NEGATED
                            << names[orig_col] << " " << abs( row_value )
                            << " * +\n";
               proof_out << "end";
               next_constraint_id += 2;
            }
#endif
            proof_out << "\n";
         }
      }
#if VERIPB_VERSION >= 2
      const auto obj_coeff = cast_to_long(stored_objective.coefficients[col]);
      assert( obj_coeff == problem.getObjective().coefficients[col] );
      if( obj_coeff != 0)
      {
         proof_out << OBJECTIVE_DIFF << (-obj_coeff ) << " " << names[orig_col] << " " << cast_to_long(obj_coeff * val)<< " ;\n";
         stored_objective.offset += obj_coeff * val;
      }
      stored_objective.coefficients[col] = 0;
      fixed_variable[col] = 1;
#endif
   }

   void
   dominating_columns( int dominating_column, int dominated_column,
                       const Vec<String>& names, const Vec<int>& var_mapping) override
   {
#if VERIPB_VERSION == 1
      if( !verification_possible )
         return;
#endif
      next_constraint_id++;
      stored_dominating_col = var_mapping[dominating_column];
      stored_dominated_col = var_mapping[dominated_column];
      const auto name_dominating = names[stored_dominating_col];
      const auto name_dominated = names[stored_dominated_col];
      proof_out << RED << "1 " << name_dominating << " +1 " << NEGATED
                << name_dominated << " >= 1 ; " << name_dominating << " -> "
                << name_dominated << " " << name_dominated << " -> "
                << name_dominating;
#if VERIPB_VERSION == 1
      add_substitutions_to_witness(names, stored_dominating_col, stored_dominated_col);
#endif
      proof_out << "\n";
   }

   void
   add_probing_reasoning( bool is_upper, int causing_col, int col,
                          const Vec<String>& names,
                          const Vec<int>& var_mapping) override
   {
#if VERIPB_VERSION == 1
      if( !verification_possible )
         return;
#endif
      const auto& name_causing_col = names[var_mapping[causing_col]];
      const auto& name_col = names[var_mapping[col]];
      next_constraint_id++;
      proof_out << RUP << "1 " << name_causing_col << " +1 ";
      if(is_upper)
         proof_out << NEGATED;
      proof_out << name_col << " >= 1;\n";
      next_constraint_id++;
      proof_out << RUP << "1 " << NEGATED << name_causing_col << " +1 ";
      if(is_upper)
         proof_out << NEGATED;
      proof_out << name_col << " >= 1;\n";
   }

   void
   store_gcd( int row, REAL gcd ) override
   {
#if VERIPB_VERSION == 1
      if( !verification_possible )
         return;
#endif
      assert( num.isIntegral( gcd ) );
      row_with_gcd = {row, cast_to_long( gcd )};
   };

   void
   store_parallel_row( int row ) override
   {
      parallel_remaining_row = row;
   };

   void
   store_implied_bound( int row, REAL lowerbound ) override
   {
      assert(lowerbound == 1 || lowerbound == 0 );
      if( lowerbound == 1 )
         row_implying_lb = row;
      else
         row_implying_ub = row;
   };


   void
   change_rhs( int row, REAL val, const SparseVectorView<REAL>& data,
               const Vec<String>& names, const Vec<int>& var_mapping,
               ArgumentType argument = ArgumentType::kPrimal ) override
   {
#if VERIPB_VERSION == 1
      if( !verification_possible )
         return;
#endif
      if( skip_changing_rhs == row )
      {
         skip_changing_rhs = UNKNOWN;
         return;
      }
      next_constraint_id++;
      switch( argument )
      {
      case ArgumentType::kPrimal:
      case ArgumentType::kPropagation:
      case ArgumentType::kAggregation:
      case ArgumentType::kDual:
      case ArgumentType::kSymmetry:
      case ArgumentType::kSaturation:
      {
         assert( num.isIntegral( val * scale_factor[row] ) );
         proof_out << RUP;
         int offset = 0;
         for( int i = 0; i < data.getLength(); i++ )
         {
            int unscaled_coeff = cast_to_long( data.getValues()[i] );
            assert( unscaled_coeff != 0 );
            auto found = changed_entries_during_current_tsxs.find( data.getIndices()[i] );
            if( found != changed_entries_during_current_tsxs.end() )
            {
               unscaled_coeff = found->second;
               if( unscaled_coeff == 0 )
                  continue;
            }
            if( i != 0 )
               proof_out << " +";
            int coeff = unscaled_coeff * scale_factor[row];
            proof_out << abs( coeff ) << " ";
            if( coeff > 0 )
            {
               offset += coeff;
               proof_out << NEGATED;
            }
            proof_out << names[var_mapping[data.getIndices()[i]]];
         }
         proof_out << " >=  "
                   << abs( offset ) -
                          cast_to_long( val ) * scale_factor[row]
                   << ";\n";
         break;
      }
      case ArgumentType::kWeakening:
      {
         assert( rhs_row_mapping[row] != UNKNOWN );
         assert( row_with_gcd.first == row );
         int gcd = row_with_gcd.second;
         proof_out << POL << rhs_row_mapping[row] << " " << gcd << " d " << gcd << " *\n";
#ifdef VERIPB_DEBUG
         assert( validate_row == row || validate_row == UNKNOWN );
         validate_row = row;
#endif
         row_with_gcd = {UNKNOWN, UNKNOWN};
         break;
      }
      default:
         assert( false );
      }
#if VERIPB_VERSION >= 2
      proof_out << MOVE_LAST_CONS_TO_CORE;
#endif
      proof_out << DELETE_CONS << rhs_row_mapping[row] << "\n";
      rhs_row_mapping[row] = next_constraint_id;
   }

   void
   change_lhs( int row, REAL val, const SparseVectorView<REAL>& data,
               const Vec<String>& names, const Vec<int>& var_mapping, ArgumentType argument = ArgumentType::kPrimal ) override
   {
#if VERIPB_VERSION == 1
      if( !verification_possible )
         return;
#endif
      if( skip_changing_lhs == row )
      {
         skip_changing_lhs = UNKNOWN;
         return;
      }
      next_constraint_id++;
      switch( argument )
      {
      case ArgumentType::kPrimal:
      case ArgumentType::kPropagation:
      case ArgumentType::kAggregation:
      case ArgumentType::kDual:
      case ArgumentType::kSymmetry:
      case ArgumentType::kSaturation:
      {
         assert( num.isIntegral( val * scale_factor[row] ) );
         proof_out << RUP;
         int offset = 0;
         for( int i = 0; i < data.getLength(); i++ )
         {
            int unscaled_coeff = cast_to_long( data.getValues()[i] );
            assert( unscaled_coeff != 0 );
            auto found = changed_entries_during_current_tsxs.find( data.getIndices()[i] );
            if( found != changed_entries_during_current_tsxs.end() )
            {
               unscaled_coeff = found->second;
               if( unscaled_coeff == 0 )
                  continue;
            }
            if( i != 0 )
               proof_out << " +";
            int coeff = unscaled_coeff * scale_factor[row];
            proof_out << abs( coeff ) << " ";
            if( coeff < 0 )
            {
               proof_out << NEGATED;
               offset += coeff;
            }
            proof_out << names[var_mapping[data.getIndices()[i]]];
         }

         proof_out << " >=  "
                   << cast_to_long( val ) * scale_factor[row] + abs( offset )
                   << ";\n";
         break;
      }
      case ArgumentType::kWeakening:
      {
         assert( lhs_row_mapping[row] != UNKNOWN );
         assert( row_with_gcd.first == row );
         int gcd = row_with_gcd.second;
         proof_out << POL << lhs_row_mapping[row] << " " << gcd << " d " << gcd << " *\n";
#ifdef VERIPB_DEBUG
         assert( validate_row == row || validate_row == UNKNOWN );
         validate_row = row;
#endif
         row_with_gcd = {UNKNOWN, UNKNOWN};
         break;
      }
      default:
         assert( false );
      }
#if VERIPB_VERSION >= 2
      proof_out << MOVE_LAST_CONS_TO_CORE;
#endif
      proof_out << DELETE_CONS << lhs_row_mapping[row] << "\n";
      lhs_row_mapping[row] = next_constraint_id;
   }

   void
   change_rhs_parallel_row( int row, REAL val, int parallel_row,
                            const Problem<REAL>& problem,
                            const Vec<int>& var_mapping ) override
   {
#if VERIPB_VERSION == 1
      if( !verification_possible )
         return;
#endif
      REAL factor_row = problem.getConstraintMatrix().getRowCoefficients( row ).getValues()[0] * scale_factor[row];
      REAL factor_parallel = problem.getConstraintMatrix().getRowCoefficients( parallel_row )
                                 .getValues()[0] *
                             scale_factor[parallel_row];
      REAL factor = factor_row / factor_parallel;
      assert( abs( factor ) >= 1 );
#ifdef VERIPB_DEBUG
      proof_out << COMMENT;
      if( factor < 0 )
         proof_out << lhs_row_mapping[parallel_row];
      else
         proof_out << rhs_row_mapping[parallel_row];
      proof_out << " is parallel to " << rhs_row_mapping[row] << "/"
                << lhs_row_mapping[row] << " are parallel.\n";
#endif
      if( abs( factor ) == 1 )
      {
         // rhs_row_mapping[row] can be KNOWN for example if Singleton relaxed a
         // constraint.
         if( rhs_row_mapping[row] != UNKNOWN )
         {
            proof_out << DELETE_CONS << rhs_row_mapping[row];
            if( factor == 1 )
               rhs_row_mapping[row] = rhs_row_mapping[parallel_row];
            else
               rhs_row_mapping[row] = lhs_row_mapping[parallel_row];
#if VERIPB_VERSION >= 2
            int used_row = rhs_row_mapping[parallel_row];
            if( factor < 0)
               used_row = lhs_row_mapping[parallel_row];
            proof_out << " ; ; begin\n\t" << POL << used_row << " -1 + \nend -1";
            next_constraint_id+=2;
#endif
            proof_out << "\n";
         }
         else
         {
            if( factor == 1 )
               rhs_row_mapping[row] = rhs_row_mapping[parallel_row];
            else
               rhs_row_mapping[row] = lhs_row_mapping[parallel_row];
         }
         if(factor < 0)
            skip_deleting_rhs_constraint_id = -rhs_row_mapping[row];
         else
            skip_deleting_rhs_constraint_id = rhs_row_mapping[row];
      }
      else
      {
         if( factor > 0 )
         {
            bool is_not_integral = false;
            if( !num.isIntegral( factor ) )
            {
               is_not_integral = true;
               factor = factor_row;
            }
            assert( rhs_row_mapping[parallel_row] != UNKNOWN );
            next_constraint_id++;
            proof_out << POL << rhs_row_mapping[parallel_row] << " " <<  (int) factor << " *\n";
#if VERIPB_VERSION >= 2
            proof_out << MOVE_LAST_CONS_TO_CORE;
#endif
            if( rhs_row_mapping[row] != UNKNOWN )
            {
               proof_out << DELETE_CONS << rhs_row_mapping[row] << "\n";
               rhs_row_mapping[row] = next_constraint_id;
#if VERIPB_VERSION >= 2
               int used_row = rhs_row_mapping[parallel_row];
               proof_out << " ; ; begin\n\t" << POL << used_row << " " << factor <<" * -1 + \nend -1";
               next_constraint_id += 2;
#endif
               proof_out << "\n";
            }
            else
               rhs_row_mapping[row] = next_constraint_id;
            //TODO: handle also lhs

            // scale also lhs
            if( lhs_row_mapping[row] != UNKNOWN && is_not_integral )
            {
               next_constraint_id++;
               proof_out << POL << lhs_row_mapping[row] << " " << (int) factor_parallel << " *\n";
#if VERIPB_VERSION >= 2
               proof_out << MOVE_LAST_CONS_TO_CORE;
#endif
               proof_out << DELETE_CONS << lhs_row_mapping[row] ;
               lhs_row_mapping[row] = next_constraint_id;
#if VERIPB_VERSION >= 2
               int used_row = lhs_row_mapping[parallel_row];
               proof_out << " ; ; begin\n\t" << POL << used_row << " " << cast_to_long(factor) <<" * -1 + \nend -1";
               next_constraint_id += 2;
#endif
               proof_out << "\n";
               scale_factor[row] *= cast_to_long( abs( factor_parallel ) );
            }
         }
         else
         {
            bool is_not_integral = false;
            if( !num.isIntegral( factor ) )
            {
               is_not_integral = true;
               factor = factor_row;
            }
            next_constraint_id++;
            assert( lhs_row_mapping[parallel_row] != UNKNOWN );
            proof_out << POL << lhs_row_mapping[parallel_row] << " " << (int) abs( factor ) << " *\n";
#if VERIPB_VERSION >= 2
            proof_out << MOVE_LAST_CONS_TO_CORE;
#endif
            if( rhs_row_mapping[row] != UNKNOWN )
            {
               proof_out << DELETE_CONS << rhs_row_mapping[row];
               rhs_row_mapping[row] = next_constraint_id;
#if VERIPB_VERSION >= 2
               int used_row = lhs_row_mapping[parallel_row];
               proof_out << " ; ; begin\n\t" << POL << used_row << " " << (int) abs( factor ) <<" * -1 + \nend -1";
               next_constraint_id += 2;
#endif
               proof_out << "\n";
            }
            else
            {
               rhs_row_mapping[row] = next_constraint_id;
            }
            // scale also lhs
            if( lhs_row_mapping[row] != UNKNOWN && is_not_integral )
            {
               next_constraint_id++;
               proof_out << POL << lhs_row_mapping[row] << " " << cast_to_long( abs( factor_parallel ) ) << " *\n";
#if VERIPB_VERSION >= 2
               proof_out << MOVE_LAST_CONS_TO_CORE;
#endif
               proof_out << DELETE_CONS << lhs_row_mapping[row];
               lhs_row_mapping[row] = next_constraint_id;
#if VERIPB_VERSION >= 2
               int used_row;
               if( factor > 0)
                  used_row = rhs_row_mapping[row];
               else
                  used_row = lhs_row_mapping[row];
               proof_out << " ; ; begin\n\t" << POL << used_row  <<" -1 " << cast_to_long (abs( factor_parallel )) << " * + \nend -1";
               next_constraint_id += 2;
#endif
               proof_out << "\n";
               scale_factor[row] *= cast_to_long( abs( factor_parallel ) );
            }
         }
      }
   }

   void
   change_lhs_parallel_row( int row, REAL val, int parallel_row,
                            const Problem<REAL>& problem ) override
   {
#if VERIPB_VERSION == 1
      if( !verification_possible )
         return;
#endif
      REAL factor_row = problem.getConstraintMatrix()
                            .getRowCoefficients( row )
                            .getValues()[0] *
                        scale_factor[row];
      REAL factor_parallel = problem.getConstraintMatrix()
                                 .getRowCoefficients( parallel_row )
                                 .getValues()[0] *
                             scale_factor[parallel_row];
      REAL factor = factor_row / factor_parallel;
      assert( abs( factor ) >= 1 );
#ifdef VERIPB_DEBUG
      proof_out << COMMENT;
      if( factor > 0 )
         proof_out << lhs_row_mapping[parallel_row];
      else
         proof_out << rhs_row_mapping[parallel_row];
      proof_out << " is parallel to " << rhs_row_mapping[row] << "/"
                << lhs_row_mapping[row] << " are parallel.\n";
#endif
      // shift the constraint ids
      if( abs( factor ) == 1 )
      {
         // lhs_row_mapping[row] can be KNOWN for example if Singleton relaxed a
         // constraint.
         if( lhs_row_mapping[row] != UNKNOWN )
         {
            proof_out << DELETE_CONS << lhs_row_mapping[row];
            if( factor == 1 )
               lhs_row_mapping[row] = lhs_row_mapping[parallel_row];
            else
               lhs_row_mapping[row] = rhs_row_mapping[parallel_row];
#if VERIPB_VERSION >= 2
            int used_row = lhs_row_mapping[parallel_row];
            if( factor < 0)
               used_row = rhs_row_mapping[parallel_row];
            proof_out << " ; ; begin\n" << POL << used_row << " -1 + \nend -1";
            next_constraint_id+=2;
#endif
            proof_out << "\n";
         }
         else
         {
            assert( ( lhs_row_mapping[row] != UNKNOWN && factor == 1 ) ||
                    ( factor == -1 && rhs_row_mapping[row] != UNKNOWN ) );
            assert( factor == 1 || factor == -1 );
            if( factor == 1 )
               lhs_row_mapping[row] = lhs_row_mapping[parallel_row];
            else
               lhs_row_mapping[row] = rhs_row_mapping[parallel_row];
         }
         if(factor > 0 )
            skip_deleting_lhs_constraint_id = lhs_row_mapping[row];
         else
            skip_deleting_rhs_constraint_id = -rhs_row_mapping[row];

      }
      else
      {
         if( factor > 0 )
         {
            bool is_not_integral = false;
            if( !num.isIntegral( factor ) )
            {
               is_not_integral = true;
               factor = factor_row;
            }
            next_constraint_id++;
            proof_out << POL << lhs_row_mapping[parallel_row] << " " << (int) factor << " *\n";
            proof_out << MOVE_LAST_CONS_TO_CORE;
            if( lhs_row_mapping[row] != UNKNOWN )
            {
               proof_out << DELETE_CONS << lhs_row_mapping[row];
               lhs_row_mapping[row] = next_constraint_id;
#if VERIPB_VERSION >= 2
               int used_row = lhs_row_mapping[parallel_row];
               if( factor < 0)
                  used_row = rhs_row_mapping[parallel_row];
               proof_out << " ; ; begin\n" << POL << used_row << " " << cast_to_long( factor ) <<" * -1 + \nend -1";
               next_constraint_id += 2;
#endif
               proof_out << "\n";
            }
            else
               lhs_row_mapping[row] = next_constraint_id;
            // scale also rhs
            if( rhs_row_mapping[row] != UNKNOWN && is_not_integral )
            {
               next_constraint_id++;
               proof_out << POL << rhs_row_mapping[row] << " "
                         << (int) factor_parallel << " *\n";
#if VERIPB_VERSION >= 2
               proof_out << MOVE_LAST_CONS_TO_CORE;
#endif
               proof_out << DELETE_CONS << rhs_row_mapping[row];
               rhs_row_mapping[row] = next_constraint_id;
#if VERIPB_VERSION >= 2
               int used_row = rhs_row_mapping[parallel_row];
               if( factor < 0)
                  used_row = lhs_row_mapping[parallel_row];
               proof_out << " ; ; begin\n" << POL << used_row << " " << cast_to_long( factor ) <<" * -1 + \nend -1";
               next_constraint_id += 2;
#endif
               proof_out << "\n";
               scale_factor[row] *= cast_to_long( abs( factor_parallel ) );
            }
         }
         else
         {
            bool is_not_integral = false;
            if( !num.isIntegral( factor ) )
            {
               is_not_integral = true;
               factor = factor_row;
            }
            next_constraint_id++;
            assert( rhs_row_mapping[parallel_row] != UNKNOWN );
            proof_out << POL << rhs_row_mapping[parallel_row] << " " << (int) abs( factor ) << " *\n";
#if VERIPB_VERSION >= 2
            proof_out << MOVE_LAST_CONS_TO_CORE;
#endif
            if( lhs_row_mapping[row] != UNKNOWN )
            {
               proof_out << DELETE_CONS << lhs_row_mapping[row];
               lhs_row_mapping[row] = next_constraint_id;
#if VERIPB_VERSION >= 2
               int used_row = lhs_row_mapping[parallel_row];
               proof_out << " ; ; begin\n\t" << POL << used_row << " " << (int) abs( factor ) <<" * -1 + \nend -1";
               next_constraint_id += 2;
#endif
               proof_out << "\n";
            }
            else
               lhs_row_mapping[row] = next_constraint_id;
            // scale also rhs
            if( rhs_row_mapping[row] != UNKNOWN && is_not_integral )
            {
               next_constraint_id++;
               proof_out << POL << rhs_row_mapping[row] << " " << (int) abs( factor_parallel ) << " *\n";
               proof_out << MOVE_LAST_CONS_TO_CORE;
               proof_out << DELETE_CONS << rhs_row_mapping[row];
               rhs_row_mapping[row] = next_constraint_id;
#if VERIPB_VERSION >= 2
               int used_row = lhs_row_mapping[parallel_row];
               proof_out << " ; ; begin\n\t" << POL << used_row << " " << (int) abs( factor ) <<" * -1 + \nend -1";
               next_constraint_id += 2;
#endif
               proof_out << "\n";
               scale_factor[row] *= cast_to_long( abs( factor_parallel ) );
            }
         }
      }
   }

   void
   change_lhs_inf( int row ) override
   {
#if VERIPB_VERSION == 1
      if( !verification_possible )
         return;
#endif
      proof_out << DELETE_CONS << lhs_row_mapping[row] << "\n";
      lhs_row_mapping[row] = UNKNOWN;
   }

   void
   change_rhs_inf( int row ) override
   {
#if VERIPB_VERSION == 1
      if( !verification_possible )
         return;
#endif
      proof_out << DELETE_CONS << rhs_row_mapping[row] << "\n";
      rhs_row_mapping[row] = UNKNOWN;
   }

   void
   change_matrix_entry( int row, int col, REAL new_val,
                        const SparseVectorView<REAL>& data, RowFlags& rflags,
                        REAL lhs, REAL rhs, const Vec<String>& names,
                        const Vec<int>& var_mapping, bool is_next_reduction_matrix_entry, ArgumentType argument ) override
   {
#if VERIPB_VERSION == 1
      if( !verification_possible )
         return;
#endif
      // remove singleton variable from equation
      changed_entries_during_current_tsxs.emplace( col, cast_to_long( new_val ) );
      if( argument == ArgumentType::kAggregation )
      {
         assert( lhs == rhs );
         assert( !rflags.test( RowFlag::kLhsInf ) );
         assert( !rflags.test( RowFlag::kRhsInf ) );
         assert( new_val == 0 );
         assert( skip_changing_lhs == UNKNOWN );
         assert( skip_changing_rhs == UNKNOWN );
         skip_changing_lhs = UNKNOWN;
         skip_changing_rhs = UNKNOWN;
         int old_value = 0;
         for( int i = 0; i < data.getLength(); i++ )
            if( data.getIndices()[i] == col )
            {
               assert(
                   num.isIntegral( data.getValues()[i] * scale_factor[row] ) );
               old_value =
                   cast_to_long( data.getValues()[i] * scale_factor[row] );
            }
         assert( old_value != 0 );

         const auto& name = names[var_mapping[col]];
         assert( num.isIntegral( new_val ) );
         int diff = old_value - cast_to_long( new_val );
         if( !rflags.test( RowFlag::kLhsInf ) )
         {
            next_constraint_id++;
            assert( lhs_row_mapping[row] != UNKNOWN );
            if( old_value > 0 )
               proof_out << POL << lhs_row_mapping[row] << " " << NEGATED
                         << name << " " << abs( diff ) << " * +\n";
            else
               proof_out << POL << lhs_row_mapping[row] << " " << name << " "
                         << abs( diff ) << " * +\n";
            skip_changing_lhs = row;
#if VERIPB_VERSION>= 2
            proof_out << MOVE_LAST_CONS_TO_CORE;
#endif
            proof_out << DELETE_CONS << lhs_row_mapping[row];
            lhs_row_mapping[row] = next_constraint_id;
#if VERIPB_VERSION >= 2
            if( old_value > 0 )
               proof_out << " ; " << name << " -> 1" ;
            else
               proof_out << " ; " << name << " -> 0" ;
#endif
            proof_out << "\n";
         }
         if( !rflags.test( RowFlag::kRhsInf ) )
         {
            next_constraint_id++;
            assert( rhs_row_mapping[row] != UNKNOWN );
            skip_changing_rhs = row;
            if( old_value < 0 )
            {
               proof_out << POL << rhs_row_mapping[row] << " " << NEGATED
                         << name << " " << abs( diff ) << " * +\n";
            }
            else
            {
               proof_out << POL << rhs_row_mapping[row] << " " << name << " "
                         << abs( diff ) << " * +\n";
            }
#if VERIPB_VERSION>= 2
            proof_out << MOVE_LAST_CONS_TO_CORE;
#endif
            proof_out << DELETE_CONS << rhs_row_mapping[row];
            rhs_row_mapping[row] = next_constraint_id;
#if VERIPB_VERSION >= 2
            if( old_value < 0 )
               proof_out << " ; " << name << " -> 1" ;
            else
               proof_out << " ; " << name << " -> 0" ;
#endif
            proof_out << "\n";
         }

         return;
      }
      else if( argument == ArgumentType::kSaturation )
      {
         if( saturation_already_called )
            return;
         next_constraint_id++;
         assert( !rflags.test( RowFlag::kEquation ) );
         assert( rflags.test( RowFlag::kRhsInf ) ||
                 rflags.test( RowFlag::kLhsInf ) );
         proof_out << POL;
         if( !rflags.test( RowFlag::kRhsInf ) )
         {
            assert( rhs_row_mapping[row] != UNKNOWN );
            proof_out << rhs_row_mapping[row] << " ";
            skip_changing_rhs = next_constraint_id;
         }
         else
         {
            assert( lhs_row_mapping[row] != UNKNOWN );
            proof_out << lhs_row_mapping[row] << " ";
            skip_changing_lhs = next_constraint_id;
         }
         proof_out << SATURATION << "\n";
#if VERIPB_VERSION >= 2
         proof_out << MOVE_LAST_CONS_TO_CORE;
#endif
         if( !rflags.test( RowFlag::kRhsInf ) )
         {
            proof_out << DELETE_CONS << rhs_row_mapping[row];
            rhs_row_mapping[row] = next_constraint_id;
#if VERIPB_VERSION >= 2
            proof_out << " ; ; begin \n\t" << POL << rhs_row_mapping[row] << " -1 +\nend -1";
            next_constraint_id += 2;
#endif
            proof_out << "\n";
         }
         else
         {
            proof_out << DELETE_CONS << lhs_row_mapping[row];
            lhs_row_mapping[row] = next_constraint_id;
#if VERIPB_VERSION >= 2
            proof_out << " ; ; begin \n\t" << POL << lhs_row_mapping[row] << " -1 +\nend -1";
            next_constraint_id += 2;
#endif
            proof_out << "\n";
         }
         skip_changing_lhs = row;
         skip_changing_rhs = row;
         saturation_already_called = true;
#ifdef VERIPB_DEBUG
         assert( validate_row == row || validate_row == UNKNOWN );
         validate_row = row;
#endif
      }
      else if( argument == ArgumentType::kWeakening )
      {

         assert( new_val == 0 );
         weakened_columns.push_back(col);
         if( is_next_reduction_matrix_entry )
         {
            return;
         }
         assert( rflags.test( RowFlag::kRhsInf ) ||
                 rflags.test( RowFlag::kLhsInf ) );
         next_constraint_id++;
         proof_out << POL;
         if( rhs_row_mapping[row] != UNKNOWN )
         {
            assert( !rflags.test( RowFlag::kRhsInf ) );
            proof_out << rhs_row_mapping[row] << " ";
         }
         else
         {
            assert( !rflags.test( RowFlag::kLhsInf ) );
            proof_out << lhs_row_mapping[row] << " ";
         }
         for( int c : weakened_columns)
            proof_out << names[var_mapping[c]] << " " << WEAKENING << " ";
         proof_out << "\n";
         weakened_columns.clear();
#if VERIPB_VERSION >= 2
         proof_out << MOVE_LAST_CONS_TO_CORE;
#endif
         if( rhs_row_mapping[row] != UNKNOWN )
         {
            proof_out << DELETE_CONS << rhs_row_mapping[row];
            rhs_row_mapping[row] = next_constraint_id;
#if VERIPB_VERSION >= 2
            proof_out << " ; ; begin\n\t" << POL << next_constraint_id << " " << row_with_gcd.second << " d " << row_with_gcd.second << " * -1 + \nend -1" ;
            next_constraint_id+=2;
#endif
            proof_out << "\n";
         }
         else
         {
            proof_out << DELETE_CONS << lhs_row_mapping[row];
            lhs_row_mapping[row] = next_constraint_id;
#if VERIPB_VERSION >= 2
            proof_out << " ; ; begin\n\t" << POL << next_constraint_id << " " << row_with_gcd.second << " d " << row_with_gcd.second << " * -1 + \nend -1" ;
            next_constraint_id+=2;
#endif
            proof_out << "\n";

         }
#ifdef VERIPB_DEBUG
         assert( validate_row == row || validate_row == UNKNOWN );
         validate_row = row;
#endif
      }
      else
      {
         assert( false );
      }
   }

   void
   sparsify( int eqrow, int candrow, REAL scale,
             const Problem<REAL>& currentProblem ) override
   {
#if VERIPB_VERSION == 1
      if( !verification_possible )
         return;
#endif
      const ConstraintMatrix<REAL>& matrix =
          currentProblem.getConstraintMatrix();
      int scale_eqrow = scale_factor[eqrow];
      int scale_candrow = scale_factor[candrow];
      assert( scale != 0 );
      REAL scale_updated = scale * scale_candrow / scale_eqrow;
      if( num.isIntegral( scale_updated ) )
      {
         int int_scale_updated = cast_to_long( scale_updated );
         if( !matrix.getRowFlags()[candrow].test( RowFlag::kRhsInf ) )
         {
            next_constraint_id++;
            assert( rhs_row_mapping[candrow] != UNKNOWN );
            assert( rhs_row_mapping[eqrow] != UNKNOWN );
            if( int_scale_updated > 0 )
               proof_out << POL << rhs_row_mapping[eqrow] << " "
                         << abs( int_scale_updated ) << " * "
                         << rhs_row_mapping[candrow] << " +\n";
            else
               proof_out << POL << lhs_row_mapping[eqrow] << " "
                         << abs( int_scale_updated ) << " * "
                         << rhs_row_mapping[candrow] << " +\n";
#if VERIPB_VERSION >= 2
            proof_out << MOVE_LAST_CONS_TO_CORE;
#endif
            proof_out << DELETE_CONS << rhs_row_mapping[candrow];
            rhs_row_mapping[candrow] = next_constraint_id;
#if VERIPB_VERSION >= 2
            proof_out << " ; ; begin \n\t";
            if( int_scale_updated > 0 )
               proof_out << POL << lhs_row_mapping[eqrow] << " "
                         << abs( int_scale_updated ) << " * "
                         << rhs_row_mapping[candrow] << " +\n";
            else
               proof_out << POL << rhs_row_mapping[eqrow] << " "
                         << abs( int_scale_updated ) << " * "
                         << rhs_row_mapping[candrow] << " +\n";
            proof_out << "end";
            next_constraint_id += 2;
#endif
            proof_out << "\n";
         }
         if( !matrix.getRowFlags()[candrow].test( RowFlag::kLhsInf ) )
         {
            next_constraint_id++;
            assert( lhs_row_mapping[candrow] != UNKNOWN );
            assert( lhs_row_mapping[eqrow] != UNKNOWN );
            if( int_scale_updated > 0 )
               proof_out << POL << lhs_row_mapping[eqrow] << " "
                         << abs( int_scale_updated ) << " * "
                         << lhs_row_mapping[candrow] << " +\n";
            else
               proof_out << POL << rhs_row_mapping[eqrow] << " "
                         << abs( int_scale_updated ) << " * "
                         << lhs_row_mapping[candrow] << " +\n";
#if VERIPB_VERSION >= 2
            proof_out << MOVE_LAST_CONS_TO_CORE;
#endif
            proof_out << DELETE_CONS << lhs_row_mapping[candrow];
            lhs_row_mapping[candrow] = next_constraint_id;
#if VERIPB_VERSION >= 2
            proof_out << " ; ; begin \n\t";
            if( int_scale_updated > 0 )
               proof_out << POL << rhs_row_mapping[eqrow] << " "
                         << abs( int_scale_updated ) << " * "
                         << lhs_row_mapping[candrow] << " +\n";
            else
               proof_out << POL << lhs_row_mapping[eqrow] << " "
                         << abs( int_scale_updated ) << " * "
                         << lhs_row_mapping[candrow] << " +\n";
            proof_out << "end";
            next_constraint_id += 2;
#endif
            proof_out << "\n";
         }
      }
      else if( num.isIntegral( 1.0 / scale_updated ) )
      {
         int int_scale_updated = cast_to_long( 1.0 / scale_updated );

         if( !matrix.getRowFlags()[candrow].test( RowFlag::kRhsInf ) )
         {
            next_constraint_id++;
            assert( rhs_row_mapping[candrow] != UNKNOWN );
            assert( rhs_row_mapping[eqrow] != UNKNOWN );
            if( int_scale_updated > 0 )
               proof_out << POL << rhs_row_mapping[candrow] << " "
                         << abs( int_scale_updated ) << " * "
                         << rhs_row_mapping[eqrow] << " +\n";
            else
               proof_out << POL << rhs_row_mapping[candrow] << " "
                         << abs( int_scale_updated ) << " * "
                         << lhs_row_mapping[eqrow] << " +\n";
#if VERIPB_VERSION>= 2
            proof_out << MOVE_LAST_CONS_TO_CORE;
#endif
            proof_out << DELETE_CONS << rhs_row_mapping[candrow] << "";
            rhs_row_mapping[candrow] = next_constraint_id;
#if VERIPB_VERSION >= 2
            proof_out << " ; ; begin \n\t";
            if( int_scale_updated > 0 )
               proof_out << POL << lhs_row_mapping[eqrow] << " "
                         << next_constraint_id  << " + "
                         << abs( int_scale_updated ) << " d \n";
            else
               proof_out << POL << rhs_row_mapping[eqrow] << " "
                         << next_constraint_id  << " + "
                         << abs( int_scale_updated ) << " d \n";
            proof_out << "end";
            next_constraint_id += 2;
#endif
            proof_out << "\n";
         }
         if( !matrix.getRowFlags()[candrow].test( RowFlag::kLhsInf ) )
         {
            next_constraint_id++;
            assert( lhs_row_mapping[candrow] != UNKNOWN );
            assert( lhs_row_mapping[eqrow] != UNKNOWN );
            if( int_scale_updated > 0 )
               proof_out << POL << lhs_row_mapping[candrow] << " "
                         << abs( int_scale_updated ) << " * "
                         << lhs_row_mapping[eqrow] << " +\n";
            else
               proof_out << POL << lhs_row_mapping[candrow] << " "
                         << abs( int_scale_updated ) << " * "
                         << rhs_row_mapping[eqrow] << " +\n";
#if VERIPB_VERSION>= 2
            proof_out << MOVE_LAST_CONS_TO_CORE;
#endif
            proof_out << DELETE_CONS << lhs_row_mapping[candrow];
            lhs_row_mapping[candrow] = next_constraint_id;
#if VERIPB_VERSION >= 2
            proof_out << " ; ; begin \n\t";
            if( int_scale_updated > 0 )
               proof_out << POL << rhs_row_mapping[eqrow] << " "
                         << next_constraint_id  << " + "
                         << abs( int_scale_updated ) << " d \n";
            else
               proof_out << POL << lhs_row_mapping[eqrow] << " "
                         << next_constraint_id  << " + "
                         << abs( int_scale_updated ) << " d \n";

            proof_out << "end";
            next_constraint_id += 2;
#endif
            proof_out << "\n";
         }
         scale_factor[candrow] *= abs(int_scale_updated);
      }
      else
      {
         auto frac = sparsify_convert_scale_to_frac( eqrow, candrow, scale, matrix );
         assert( frac.second / frac.first == -scale );
         int frac_eqrow =
             abs( cast_to_long( frac.second * scale_candrow ) );
         int frac_candrow = abs( cast_to_long( frac.first * scale_eqrow ) );

         if( !matrix.getRowFlags()[candrow].test( RowFlag::kRhsInf ) )
         {
            next_constraint_id++;
            assert( rhs_row_mapping[candrow] != UNKNOWN );
            assert( rhs_row_mapping[eqrow] != UNKNOWN );
            if( scale > 0 )
               proof_out << POL << rhs_row_mapping[candrow] << " "
                         << frac_candrow << " * " << rhs_row_mapping[eqrow]
                         << " " << frac_eqrow << " * +\n";
            else
               proof_out << POL << rhs_row_mapping[candrow] << " "
                         << frac_candrow << " * " << lhs_row_mapping[eqrow]
                         << " " << frac_eqrow << " * +\n";
#if VERIPB_VERSION>= 2
            proof_out << MOVE_LAST_CONS_TO_CORE;
#endif
            proof_out << DELETE_CONS << rhs_row_mapping[candrow];
            rhs_row_mapping[candrow] = next_constraint_id;
#if VERIPB_VERSION >= 2
            proof_out << " ; ; begin \n\t";
            if( scale > 0 )
               proof_out << POL << rhs_row_mapping[candrow] << " " << lhs_row_mapping[eqrow] << " " << frac_eqrow << " * + "
                         << frac_candrow << " d "  <<"\n";
            else
               proof_out << POL << rhs_row_mapping[candrow] << " " << rhs_row_mapping[eqrow] << " " << frac_eqrow << " * + "
                         << frac_candrow << " d "  <<"\n";
            proof_out << "end";
            next_constraint_id += 2;
#endif
            proof_out << "\n";
         }
         if( !matrix.getRowFlags()[candrow].test( RowFlag::kLhsInf ) )
         {
            next_constraint_id++;
            assert( lhs_row_mapping[candrow] != UNKNOWN );
            assert( lhs_row_mapping[eqrow] != UNKNOWN );
            if( scale > 0 )
               proof_out << POL << lhs_row_mapping[candrow] << " "
                         << frac_candrow << " * " << lhs_row_mapping[eqrow]
                         << " " << frac_eqrow << " * +\n";
            else
               proof_out << POL << lhs_row_mapping[candrow] << " "
                         << frac_candrow << " * " << rhs_row_mapping[eqrow]
                         << " " << frac_eqrow << " * +\n";
#if VERIPB_VERSION>= 2
            proof_out << MOVE_LAST_CONS_TO_CORE;
#endif
            proof_out << DELETE_CONS << lhs_row_mapping[candrow];
            lhs_row_mapping[candrow] = next_constraint_id;
#if VERIPB_VERSION >= 2
            proof_out << " ; ; begin \n\t";
            if( scale > 0 )
               proof_out << POL << lhs_row_mapping[candrow] << " " << rhs_row_mapping[eqrow] << " " << frac_eqrow << " * + "
                         << frac_candrow << " d "  <<"\n";
            else
               proof_out << POL << lhs_row_mapping[candrow] << " " << lhs_row_mapping[eqrow] << " " << frac_eqrow << " * + "
                         << frac_candrow << " d "  <<"\n";
            proof_out << "end";
            next_constraint_id += 2;
#endif
            proof_out << "\n";
         }
         scale_factor[candrow] *= abs(frac_candrow);
      }
   }

   void
   substitute( int col, const SparseVectorView<REAL>& equality, REAL offset, REAL old_obj_coeff,
               const Problem<REAL>& currentProblem, const Vec<String>& names,
               const Vec<int>& var_mapping ) override {
#if VERIPB_VERSION == 1
      if( !verification_possible )
         return;
#endif
      assert( num.isIntegral( offset ) );
      const REAL* values = equality.getValues();
      const int* indices = equality.getIndices();
      assert( equality.getLength() == 2 );
      assert( num.isIntegral( values[0] ) && num.isIntegral( values[1] ) );
      assert( cast_to_long( values[0] ) != 0 );
      assert( cast_to_long( values[1] ) != 0 );
      REAL substitute_factor = indices[0] == col ? values[0] : values[1];

      next_constraint_id++;
      int orig_index_0 = var_mapping[indices[0]];
      int orig_index_1 = var_mapping[indices[1]];

      int first_constraint_id = next_constraint_id;
#ifdef VERIPB_DEBUG
      proof_out << COMMENT << "postsolve stack : row id " << next_constraint_id << "\n";
#endif
      proof_out << RUP;

      int lhs = cast_to_long( offset );
      proof_out << abs( cast_to_long( values[0] ) ) << " ";
      if( values[0] < 0 )
      {
         proof_out << NEGATED;
         lhs += abs( cast_to_long( values[0] ) );
      }
      proof_out << names[orig_index_0] << " +"
                << abs( cast_to_long( values[1] ) ) << " ";
      if( values[1] < 0 )
      {
         proof_out << NEGATED;
         lhs += abs( cast_to_long( values[1] ) );
      }
      proof_out << names[orig_index_1] << " >= " << lhs << ";\n";
      int lhs_id = next_constraint_id;
#if VERIPB_VERSION >= 2
      proof_out << MOVE_LAST_CONS_TO_CORE;
#endif

      next_constraint_id++;
      int second_constraint_id = next_constraint_id;

#ifdef VERIPB_DEBUG
      proof_out << COMMENT << "postsolve stack : row id " << next_constraint_id << "\n";
#endif
      proof_out << RUP;
      int rhs = -cast_to_long( offset );
      proof_out << abs( cast_to_long( values[0] ) ) << " ";
      if( values[0] > 0 )
      {
         proof_out << NEGATED;
         rhs += abs( cast_to_long( values[0] ) );
      }
      proof_out << names[orig_index_0] << " +"
                << abs( cast_to_long( values[1] ) ) << " ";
      if( values[1] > 0 )
      {
         proof_out << NEGATED;
         rhs += abs( cast_to_long( values[1] ) );
      }
      proof_out << names[orig_index_1] << " >= " << rhs << ";\n";

#if VERIPB_VERSION >= 2
      proof_out << MOVE_LAST_CONS_TO_CORE;
#endif

      substitute( col, substitute_factor, lhs_id, next_constraint_id,
                  currentProblem );
#if VERIPB_VERSION == 1
      if(!is_optimization_problem)
      {
#endif
#if VERIPB_VERSION >= 2
         apply_substitution_to_objective(col, equality, offset);
         if(old_obj_coeff != 0)
         {
            proof_out << OBJECTIVE_DIFF ;
            for( int i=0 ; i < 2; i++)
               if(indices[i] == col)
                  proof_out << cast_to_long( -old_obj_coeff ) << " " << names[var_mapping[indices[i]]] << " ";
               else
                  proof_out << cast_to_long(-old_obj_coeff * values[0]/values[1]) << " " << names[var_mapping[indices[i]]] << " ";
            proof_out << cast_to_long(offset * old_obj_coeff * values[0]/values[1]) << ";";
            if ( abs(old_obj_coeff) != 1 )
            {
               proof_out << " ; begin\n\tproofgoal #1\n\t\t" << POL;
               if(old_obj_coeff/ substitute_factor < 0 )
                  proof_out << first_constraint_id << " " << cast_to_long( abs( old_obj_coeff ) ) << " * " << " -1 " << cast_to_long( abs( substitute_factor ) ) << " * +";
               else
                  proof_out << second_constraint_id << " " << cast_to_long( abs( old_obj_coeff ) ) << " * " << " -1 " << cast_to_long( abs( substitute_factor ) ) << " * +";
               proof_out << "\t\nend -1\n\tproofgoal #2\n\t\t" << POL;
               if(old_obj_coeff/ substitute_factor > 0 )
                  proof_out << first_constraint_id << " " << cast_to_long( abs( old_obj_coeff ) ) << " * " << " -1 " << cast_to_long( abs( substitute_factor ) ) << " * +";
               else
                  proof_out << second_constraint_id << " " << cast_to_long( abs( old_obj_coeff ) ) << " * " << " -1 " << cast_to_long( abs( substitute_factor ) ) << " * +";
               proof_out << "\t\nend -1\nend";
               next_constraint_id += 4;
            }
            proof_out << "\n";
         }
#endif

         proof_out << DELETE_CONS << first_constraint_id << " ; ";
#if VERIPB_VERSION >= 2
         int index = indices[0] == col ? 0 : 1;
         const auto& name = names[var_mapping[indices[index]]];
         int value = (values[index] > 0) ? 1 : 0;
         proof_out << name << " -> " << value;
#endif
         proof_out << "\n";

         proof_out << DELETE_CONS << second_constraint_id << " ; ";
#if VERIPB_VERSION >= 2
         value = (values[index] > 0) ? 0 : 1;
         proof_out << name << " -> " << value;
#endif
         proof_out << "\n";
#if VERIPB_VERSION == 1
      }
      else
         store_substitution( values[1], values[0], orig_index_1, orig_index_0 );
#endif
   }

   void
   substitute_col_singleton_implied( int col, int row, REAL old_obj_coeff, const Problem<REAL>& currentProblem, const Vec<int>& var_mapping  ) override
   {
#if VERIPB_VERSION == 1
      if( !verification_possible || (matrix.getColumnCoefficients( col ).getLength() == 1 && !is_optimization_problem))
         return;
#endif
      const ConstraintMatrix<REAL>& matrix =
          currentProblem.getConstraintMatrix();
      auto col_vec = matrix.getColumnCoefficients( col );
      auto row_data = matrix.getRowCoefficients( row );
#if VERIPB_VERSION == 1
      if( currentProblem.getConstraintMatrix().getRowSizes()[row] > 2 && is_optimization_problem)
      {
         fmt::print("Verification currently not possible for multi-aggregations for optimization problem!\n");
         proof_out << "* Verification currently not possible for multi-aggregations for optimization problem!\n";
         verification_possible = false;
         return;
      }
#endif
      REAL substitute_factor = get_coeff_in_col_vec(row, col_vec);

      const String name = currentProblem.getVariableNames()[var_mapping[col]];
      assert(col_vec.getLength() == 1);
      apply_substitution_to_objective(col, row_data, matrix.getLeftHandSides()[row]);
      if( old_obj_coeff != 0)
      {
         proof_out << OBJECTIVE_DIFF << cast_to_long( -old_obj_coeff ) << " " << name << " ";
         REAL factor = old_obj_coeff / substitute_factor;
         REAL offset = currentProblem.getConstraintMatrix().getRightHandSides()[row];
         for(int i=0; i<row_data.getLength(); i++)
         {
            int index = row_data.getIndices()[i];
            if( index == col || fixed_variable[index]  == -1 )
               continue;
            if(fixed_variable[index]  == 1)
            {
               offset -= row_data.getValues()[i];
               continue;
            }
            proof_out  << cast_to_long( -factor * row_data.getValues()[i] ) << " " << currentProblem.getVariableNames()[var_mapping[index]] << " ";
         }
         proof_out << cast_to_long(offset * factor) << ";";

         if( abs(old_obj_coeff) != 1)
         {
            proof_out << " ; begin\n\tproofgoal #1\n\t\t" << POL;
            if(old_obj_coeff/ substitute_factor < 0 )
               proof_out << lhs_row_mapping[row] << " " << cast_to_long( abs( old_obj_coeff ) ) << " * " << " -1 " << cast_to_long( abs( substitute_factor ) ) << " * +";
            else
               proof_out << rhs_row_mapping[row] << " " << cast_to_long( abs( old_obj_coeff ) ) << " * " << " -1 " << cast_to_long( abs( substitute_factor ) ) << " * +";
            proof_out << "\nend -1\n\tproofgoal #2\n\t\t" << POL;
            if(old_obj_coeff/ substitute_factor > 0 )
               proof_out << lhs_row_mapping[row] << " " << cast_to_long( abs( old_obj_coeff ) ) << " * " << " -1 " << cast_to_long( abs( substitute_factor ) ) << " * +";
            else
               proof_out << rhs_row_mapping[row] << " " << cast_to_long( abs( old_obj_coeff ) ) << " * " << " -1 " << cast_to_long( abs( substitute_factor ) ) << " * +";
            proof_out << "\nend -1\nend";
            next_constraint_id += 4;
         }
         proof_out << "\n";
      }
      proof_out << DELETE_CONS << lhs_row_mapping[row];
#if VERIPB_VERSION >= 2
      if( substitute_factor > 0)
         proof_out << " ; " << name << " -> 1";
      else
         proof_out << " ; " << name << " -> 0";
#endif
      proof_out << "\n";
      proof_out << DELETE_CONS << rhs_row_mapping[row];
#if VERIPB_VERSION >= 2
      if( substitute_factor < 0)
         proof_out << " ; " << name << " -> 1";
      else
         proof_out << " ; " << name << " -> 0";
#endif
      proof_out << "\n";
      skip_deleting_rhs_constraint_id = rhs_row_mapping[row];
      skip_deleting_lhs_constraint_id = lhs_row_mapping[row];
   };

   void
   substitute( int col, int substituted_row, REAL old_obj_coeff, const Problem<REAL>& currentProblem, const Vec<int>& var_mapping, ArgumentType argument  )  override {
#if VERIPB_VERSION == 1
      if( !verification_possible || (matrix.getColumnCoefficients( col ).getLength() == 1 && !is_optimization_problem))
         return;
#endif
      const ConstraintMatrix<REAL>& matrix =
          currentProblem.getConstraintMatrix();
      auto col_vec = matrix.getColumnCoefficients( col );
      auto row_data = matrix.getRowCoefficients( substituted_row );

#if VERIPB_VERSION == 1
      if( currentProblem.getConstraintMatrix().getRowSizes()[substituted_row] > 2 && is_optimization_problem)
      {
         fmt::print("Verification currently not possible for multi-aggregations for optimization problem!\n");
         proof_out << "* Verification currently not possible for multi-aggregations for optimization problem!\n";
         verification_possible = false;
         return;
      }
#endif
      REAL substitute_factor = get_coeff_in_col_vec( substituted_row, col_vec );
      assert(substitute_factor != 0);
      const String name = currentProblem.getVariableNames()[var_mapping[col]];
#if VERIPB_VERSION >= 2
      int implied_constraint_lhs = -1;
      int implied_constraint_rhs = -1;
      if( row_implying_lb != UNKNOWN || row_implying_ub != UNKNOWN)
      {
         if( substitute_factor > 0)
         {
            proof_out << POL << lhs_row_mapping[substituted_row] << " " << NEGATED << name << " " << abs( cast_to_long( substitute_factor )) << " * +\n";
            proof_out << MOVE_LAST_CONS_TO_CORE;
            proof_out << POL << rhs_row_mapping[substituted_row] << " " << name << " " << abs( cast_to_long( substitute_factor )) << " * +\n";
            proof_out << MOVE_LAST_CONS_TO_CORE;
            implied_constraint_lhs = next_constraint_id + 1;
            implied_constraint_rhs = next_constraint_id + 2;
         }
         else
         {
            proof_out << POL << rhs_row_mapping[substituted_row] << " " << NEGATED << name << " " << abs( cast_to_long( substitute_factor )) << " * +\n";
            proof_out << MOVE_LAST_CONS_TO_CORE;
            proof_out << POL << lhs_row_mapping[substituted_row] << " " << name << " " << abs( cast_to_long( substitute_factor )) << " * +\n";
            proof_out << MOVE_LAST_CONS_TO_CORE;
            implied_constraint_rhs = next_constraint_id + 1;
            implied_constraint_lhs = next_constraint_id + 2;
         }
         next_constraint_id++;
         next_constraint_id++;
      }
#endif
      if( col_vec.getLength() != 1)
      {
         substitute( col, substitute_factor, lhs_row_mapping[substituted_row],
                     rhs_row_mapping[substituted_row], currentProblem,
                     substituted_row );
      }
      else if (argument != ArgumentType::kAggregation)
      {

         skip_deleting_lhs_constraint_id = lhs_row_mapping[substituted_row];
         skip_deleting_rhs_constraint_id = rhs_row_mapping[substituted_row];
      }
      assert( !matrix.getRowFlags()[substituted_row].test( RowFlag::kRhsInf ) );
      assert( !matrix.getRowFlags()[substituted_row].test( RowFlag::kLhsInf ) );
#if VERIPB_VERSION == 1
      if(!is_optimization_problem)
      {
#endif
#if VERIPB_VERSION >= 2
         assert(matrix.getLeftHandSides()[substituted_row] == matrix.getRightHandSides()[substituted_row]);
         apply_substitution_to_objective(col, row_data, matrix.getLeftHandSides()[substituted_row]);
         if( old_obj_coeff != 0)
         {
            proof_out << OBJECTIVE_DIFF << cast_to_long( -old_obj_coeff ) << " " << name << " ";
            REAL factor = old_obj_coeff / substitute_factor;
            REAL offset = currentProblem.getConstraintMatrix().getRightHandSides()[substituted_row];
            for( int i = 0; i < row_data.getLength(); i++ )
            {
               int index = row_data.getIndices()[i];
               if( index == col || fixed_variable[index]  == -1  )
                  continue;
               if(fixed_variable[index]  == 1)
               {
                  offset -= row_data.getValues()[i];
                  continue;
               }
               proof_out  << cast_to_long(( -factor ) * row_data.getValues()[i] ) << " " << currentProblem.getVariableNames()[var_mapping[index]] << " ";
            }
            proof_out << cast_to_long( offset * factor) << ";";
            if( abs(old_obj_coeff) != 1)
            {
               proof_out << " ; begin\n\tproofgoal #1\n\t\t" << POL;
               if(old_obj_coeff/ substitute_factor < 0 )
                  proof_out << lhs_row_mapping[substituted_row] << " " << cast_to_long( abs( old_obj_coeff ) ) << " * " << " -1 " << cast_to_long( abs( substitute_factor ) ) << " * +";
               else
                  proof_out << rhs_row_mapping[substituted_row] << " " << cast_to_long( abs( old_obj_coeff ) ) << " * " << " -1 " << cast_to_long( abs( substitute_factor ) ) << " * +";
               proof_out << "\nend -1\n\tproofgoal #2\n\t\t" << POL;
               if(old_obj_coeff/ substitute_factor > 0 )
                  proof_out << lhs_row_mapping[substituted_row] << " " << cast_to_long( abs( old_obj_coeff ) ) << " * " << " -1 " << cast_to_long( abs( substitute_factor ) ) << " * +";
               else
                  proof_out << rhs_row_mapping[substituted_row] << " " << cast_to_long( abs( old_obj_coeff ) ) << " * " << " -1 " << cast_to_long( abs( substitute_factor ) ) << " * +";
               proof_out << "\nend -1\nend";
               next_constraint_id += 4;
            }
            proof_out << "\n";

         }
         if (argument == ArgumentType::kAggregation)
            return;
#endif
#ifdef VERIPB_DEBUG
         proof_out << COMMENT << "postsolve stack : row id "
                   << lhs_row_mapping[substituted_row] << "\n";
#endif
         proof_out << DELETE_CONS << rhs_row_mapping[substituted_row];
#if VERIPB_VERSION >= 2
         if( substitute_factor > 0)
            proof_out << " ; " << name << " -> 0";
         else
            proof_out << " ; " << name << " -> 1";
         if( row_implying_lb != UNKNOWN || row_implying_ub != UNKNOWN)
         {
            proof_out << " ; begin\n\t" << POL << implied_constraint_lhs << " " ;
                if( substitute_factor < 0)
               proof_out << NEGATED ;
            proof_out << name << " " << abs( cast_to_long( substitute_factor )) << " * +\nend";
            next_constraint_id += 2;
         }
#endif
         proof_out << "\n";

#ifdef VERIPB_DEBUG
         proof_out << COMMENT << "postsolve stack : row id "
                << rhs_row_mapping[substituted_row] << "\n";
#endif
         proof_out << DELETE_CONS << lhs_row_mapping[substituted_row];
#if VERIPB_VERSION >= 2
         if( substitute_factor < 0)
            proof_out << " ; " <<  name << " -> 0";
         else
            proof_out << " ; " << name << " -> 1";

         if( row_implying_lb != UNKNOWN || row_implying_ub != UNKNOWN)
         {
            proof_out << " ; begin\n\t" << POL << implied_constraint_rhs << " " ;
            if( substitute_factor > 0)
               proof_out << NEGATED ;
            proof_out << name << " " << abs( cast_to_long( substitute_factor )) << " * +\nend";
            next_constraint_id += 2;
         }
#endif
         proof_out << "\n";

         if( row_implying_lb != UNKNOWN || row_implying_ub != UNKNOWN)
         {
            /**
             * lb implied externally + coeff > 0 -> proof rhs
             * ub implied externally + coeff > 0 -> proof lhs
             * lb implied externally + coeff < 0 -> proof lhs
             * ub implied externally + coeff < 0 -> proof rhs
             */
            if( row_implying_lb != UNKNOWN && substitute_factor > 0)
            {
               REAL coeff_lb = get_coeff_in_col_vec( row_implying_lb, col_vec );
               int contradiction = coeff_lb > 0 ? lhs_row_mapping[row_implying_lb]: rhs_row_mapping[row_implying_lb];
               assert(contradiction != UNKNOWN);
               proof_out << DELETE_CONS << implied_constraint_rhs << " ; ; begin\n\t" << POL << contradiction
                         << " -1 +\nend\n";
               next_constraint_id += 2;
            }
            else if ( row_implying_ub != UNKNOWN && substitute_factor < 0)
            {
               REAL coeff_ub = get_coeff_in_col_vec( row_implying_ub, col_vec );
               int contradiction = coeff_ub > 0 ? rhs_row_mapping[row_implying_ub]: lhs_row_mapping[row_implying_ub];
               assert(contradiction != UNKNOWN);
               proof_out << DELETE_CONS << implied_constraint_rhs << " ; ; begin\n\t" << POL << contradiction << " -1 +\nend\n";
               next_constraint_id += 2;
            }
            else
               proof_out << DELETE_CONS << implied_constraint_rhs << "\n";

            if( row_implying_ub != UNKNOWN && substitute_factor > 0)
            {

               assert( (get_coeff_in_col_vec( row_implying_ub, col_vec ) >0  && rhs_row_mapping[row_implying_ub] != UNKNOWN)
                       || (lhs_row_mapping[row_implying_ub] != UNKNOWN));
               proof_out << DELETE_CONS << implied_constraint_lhs << " ; ; begin\n\t" << POL << rhs_row_mapping[row_implying_ub] << " -1 +\nend\n";
               next_constraint_id += 2;
            }
            else if ( row_implying_lb != UNKNOWN && substitute_factor < 0)
            {
               assert( (get_coeff_in_col_vec( row_implying_lb, col_vec ) >0  && lhs_row_mapping[row_implying_ub] != UNKNOWN)
                       || (rhs_row_mapping[row_implying_ub] != UNKNOWN));
               proof_out << DELETE_CONS << implied_constraint_lhs << " ; ; begin\n\t" << POL << rhs_row_mapping[row_implying_lb] << " -1 +\nend\n";
               next_constraint_id += 2;
            }
            else
               proof_out << DELETE_CONS << implied_constraint_lhs << "\n";
         }

#if VERIPB_VERSION == 1
      }
      else
      {
         auto row = currentProblem.getConstraintMatrix().getRowCoefficients(substituted_row);
         assert(row.getLength() == 2);
         int col2 = row.getIndices()[0] == col ? row.getIndices()[1] : row.getIndices()[0];
         store_substitution(row.getValues()[0], row.getValues()[1],
                             var_mapping[col], var_mapping[col2]);
         if( col_vec.getLength() == 1 )
         {
            skip_deleting_lhs_constraint_id = lhs_row_mapping[substituted_row];
            skip_deleting_rhs_constraint_id = rhs_row_mapping[substituted_row];
         }
      }
#endif
   }

   void
   mark_row_redundant( int row, const Problem<REAL>& currentProblem, ArgumentType argument = ArgumentType::kPrimal ) override
   {
#if VERIPB_VERSION == 1
      if( !verification_possible )
         return;
#endif
      if(status == -2)
         return;
      assert( lhs_row_mapping[row] != UNKNOWN ||
              rhs_row_mapping[row] != UNKNOWN );
      if( lhs_row_mapping[row] != UNKNOWN )
      {
         if( lhs_row_mapping[row] == skip_deleting_lhs_constraint_id )
            skip_deleting_lhs_constraint_id = UNKNOWN;
         else if( lhs_row_mapping[row] == -skip_deleting_rhs_constraint_id )
            skip_deleting_rhs_constraint_id = UNKNOWN;
         else
         {
            proof_out << DELETE_CONS << lhs_row_mapping[row];
            if(argument == ArgumentType::kRedundant)
            {
               assert(parallel_remaining_row != UNKNOWN);
               int coeff_remaining =
                   cast_to_long(
                       currentProblem.getConstraintMatrix()
                           .getRowCoefficients( parallel_remaining_row )
                           .getValues()[0] ) * scale_factor[parallel_remaining_row];
               int coeff =
                   cast_to_long( currentProblem.getConstraintMatrix()
                                         .getRowCoefficients( row )
                                         .getValues()[0] ) * scale_factor[row];
               if(abs(coeff/coeff_remaining) != 1)
               {
                  int use_row_to_proof =lhs_row_mapping[parallel_remaining_row];
                  if( coeff * 1.0 / coeff_remaining < 0 )
                     use_row_to_proof = rhs_row_mapping[parallel_remaining_row];
                  proof_out << " ; ; begin\n\t" << POL << use_row_to_proof
                            << " " << abs( coeff ) << " * -1 "
                            << abs( coeff_remaining ) << " * +\nend -1";
                  next_constraint_id += 2;
               }
            }
            proof_out << "\n";
         }
         lhs_row_mapping[row] = UNKNOWN;
      }
      if( rhs_row_mapping[row] != UNKNOWN )
      {
         if( rhs_row_mapping[row] == -skip_deleting_lhs_constraint_id )
            skip_deleting_lhs_constraint_id = UNKNOWN;
         else if( rhs_row_mapping[row] == skip_deleting_rhs_constraint_id )
            skip_deleting_rhs_constraint_id = UNKNOWN;
         else
         {
            proof_out << DELETE_CONS << rhs_row_mapping[row];
            if( argument == ArgumentType::kRedundant )
            {
               assert( parallel_remaining_row != UNKNOWN );
               int coeff_remaining =
                   cast_to_long(
                       currentProblem.getConstraintMatrix()
                           .getRowCoefficients( parallel_remaining_row )
                           .getValues()[0] )* scale_factor[parallel_remaining_row];
               int coeff =
                   cast_to_long( currentProblem.getConstraintMatrix()
                                         .getRowCoefficients( row )
                                         .getValues()[0] )
                        * scale_factor[row];
               if(abs(coeff/coeff_remaining) != 1)
               {
                  int use_row_to_proof =rhs_row_mapping[parallel_remaining_row];
                  if( coeff *1.0/ coeff_remaining < 0 )
                     use_row_to_proof = lhs_row_mapping[parallel_remaining_row];
                  proof_out << " ; ; begin\n\t" << POL << use_row_to_proof
                            << " " << abs( coeff ) << " * -1 "
                            << abs( coeff_remaining ) << " * +\nend -1";
                  next_constraint_id+=2;
               }
            }
            proof_out << "\n";
         }
         rhs_row_mapping[row] = UNKNOWN;
      }
   }

   void
   log_solution( const Solution<REAL>& orig_solution, const Vec<String>& names, REAL origobj ) override
   {
#if VERIPB_VERSION == 1
      if( !verification_possible )
         return;
      proof_out << "o";
#else
      if( is_optimization_problem )
         proof_out << "o";
      else
         proof_out << "sol";
#endif
      next_constraint_id++;
      for( unsigned int i = 0; i < orig_solution.primal.size(); i++ )
      {
         assert( orig_solution.primal[i] == 0 || orig_solution.primal[i] == 1 );
         proof_out << " ";
         if( orig_solution.primal[i] == 0 )
            proof_out << NEGATED;
         proof_out << names[i];
      }
      next_constraint_id++;
      proof_out << "\n";
#if VERIPB_VERSION == 1
      proof_out << "u >= 1 ;\n";
      proof_out << "c " << next_constraint_id << "\n";
#else
      status = 1;
      end_proof((int)origobj);
#endif
   };

   void
   setInfeasibleCause(int col) override
   {
      cause = col;
   }

   void
   infeasible( ) override
   {
      if( status == -2)
         return;
#if VERIPB_VERSION == 1
      if( !verification_possible )
         return;
#endif
      next_constraint_id++;
      proof_out << "u >= 1 ;\n";
#if VERIPB_VERSION == 1
      proof_out << "c " << next_constraint_id << "\n";
#else
      status = -1;
      end_proof();
#endif
   };

   void
   infeasible( const Vec<int>& colmapping, const Vec<String>& names ) override
   {
#if VERIPB_VERSION == 1
      if( !verification_possible )
         return;
#endif
      if( status == -2)
         return;
      if(cause != -1)
      {
         next_constraint_id++;
         proof_out << RUP << "1 " << names[colmapping[cause]] << " >= 1 ;\n";
      }
      next_constraint_id++;
      proof_out << "u >= 1 ;\n";
#if VERIPB_VERSION == 1
      proof_out << "c " << next_constraint_id << "\n";
#else
      status = -1;
      end_proof();
#endif
   };

   void
   end_proof( ) override
   {
      end_proof( 0 ) ;
   }



   void
   symmetries( const SymmetryStorage& symmetries, const Vec<String>& names,
               const Vec<int>& var_mapping ) override
   {
#if VERIPB_VERSION == 1
      if( !verification_possible )
         return;
#endif
      if( symmetries.symmetries.empty() )
         return;
#ifdef VERIPB_DEBUG
      proof_out << COMMENT << "symmetries: \n";
#endif
      for(Symmetry symmetry: symmetries.symmetries)
      {
         const auto& name_dominating = names[var_mapping[symmetry.getDominatingCol()]];
         const auto& name_dominated = names[var_mapping[symmetry.getDominatedCol()]];
         switch( symmetry.getSymmetryType() )
         {
         case SymmetryType::kXgeY:
            proof_out << RED << "1 " << name_dominating << " +1 " << NEGATED
                      << name_dominated << " >= 1 ; " << name_dominating << " -> "
                      << name_dominated << " " << name_dominated << " -> "
                      << name_dominating << "\n";
            break;
         case SymmetryType::kXplusYge1:
            proof_out << RED << "1 " << name_dominating << " +1 "
                      << name_dominated << " >= 1 ; " << name_dominating << " -> ~"
                      << name_dominated << " " << name_dominated << " -> ~"
                      << name_dominating << "\n";
            break;
         default:
            assert(false);
         }

      }
   };

   void
   compress( const Vec<int>& rowmapping, const Vec<int>& colmapping,
             bool full = false ) override
   {
      flush();
#ifdef PAPILO_TBB
      tbb::parallel_invoke(
          [this, &rowmapping, full]()
          {
             compress_vector( rowmapping, lhs_row_mapping );
             if( full )
                lhs_row_mapping.shrink_to_fit();
          },
          [this, &rowmapping, full]()
          {
             compress_vector( rowmapping, scale_factor );
             if( full )
                scale_factor.shrink_to_fit();
          },
          [this, &colmapping, full]()
          {
             compress_vector( colmapping, fixed_variable );
             if( full )
                fixed_variable.shrink_to_fit();
          },
          [this, &colmapping, full]()
          {
             REAL count = 0;
             for(auto v: stored_objective.coefficients)
                count += v;
             compress_vector( colmapping, stored_objective.coefficients );
             REAL count2 = 0;
             for(auto v: stored_objective.coefficients)
                count2 += v;
             if( full )
                stored_objective.coefficients.shrink_to_fit();
          },
          [this, &rowmapping, full]()
          {
             compress_vector( rowmapping, rhs_row_mapping );
             if( full )
                rhs_row_mapping.shrink_to_fit();
          } );
#else
      compress_vector( rowmapping, lhs_row_mapping );
      compress_vector( rowmapping, rhs_row_mapping );
      compress_vector( rowmapping, scale_factor );
      compress_vector( colmapping, stored_objective.coefficients );
      if( full )
      {
         rhs_row_mapping.shrink_to_fit();
         lhs_row_mapping.shrink_to_fit();
         scale_factor.shrink_to_fit();
         stored_objective.coefficients.shrink_to_fit();
      }
#endif
   }

   void
   log_forcing_row ( int row ) override {
      row_forcing_propagation = row;
   }


 private:

   REAL
   get_coeff_in_col_vec( int substituted_row,
                       const SparseVectorView<REAL>& col_vec )
   {
      for( int i = 0; i < col_vec.getLength(); i++ )
         if( col_vec.getIndices()[i] == substituted_row )
            return col_vec.getValues()[i] * scale_factor[substituted_row];
      assert(false);
      return 0;
   };


   void
   end_proof( int obj )
   {
      if( status == -2)
         return;
#if VERIPB_VERSION >= 2
      proof_out << OUTPUT << NONE << " \n";
      proof_out << CONCLUSION;
      if(is_optimization_problem)
      {
         if( status > 0 )
            proof_out << "BOUNDS " << obj << " " << obj;
         else if( status < 0 )
            proof_out << " BOUNDS INF INF";
         else
            proof_out << NONE;
      }
      else
      {
         if( status > 0 )
            proof_out << "SAT";
         else if( status < 0 )
            proof_out << "UNSAT";
         else
            proof_out << NONE;
      }
      proof_out << "\n";
      proof_out << "end pseudo-Boolean proof\n";
      status = -2;
#endif
   };

#if VERIPB_VERSION >= 2
   
   void
   apply_substitution_to_objective(int sub_col, const SparseVectorView<REAL>& equality, REAL rhs)
   {
      if(stored_objective.coefficients[sub_col] == 0 )
         return;
      REAL factor = 0;
      const REAL* values = equality.getValues();
      const int* indices = equality.getIndices();
      for(int i=0; i<equality.getLength(); i++)
         if( indices[i] == sub_col)
         {
            factor = stored_objective.coefficients[sub_col]/ values[i];
            break;
         }
      for(int i=0; i<equality.getLength(); i++)
      {
         if( indices[i] == sub_col || fixed_variable[indices[i]]  == -1)
            continue;
         if( fixed_variable[indices[i]]  == 1 )
         {
            stored_objective.offset -= factor * values[i];
            continue;
         }
         stored_objective.coefficients[indices[i]] -= factor * values[i];
      }
      stored_objective.offset += rhs * factor ;
      stored_objective.coefficients[sub_col] = 0;
   }



#endif

   long
   cast_to_long( const REAL& x )
   {
      return (long) floor( REAL( x + REAL( 0.5 ) ) );
   }

#if VERIPB_VERSION == 1
   void
   store_substitution( const REAL value_0, const REAL value_1, int orig_index_0,
                       int orig_index_1 )
   {
      if (  substitutions.find(orig_index_1) == substitutions.end() )
      {
         Vec<int> v {};
         substitutions.insert_or_assign(orig_index_1, v);
      }
      int sign = -(int) value_1 / (int) value_0;
      assert( abs(sign) == 1);
      int index = orig_index_0;
      if (sign < 0 )
         index = -orig_index_0 -1;

      auto it1 = substitutions.find(orig_index_1);
      assert(it1 != substitutions.end());
      it1->second.push_back(index);

      auto it0 = substitutions.find(orig_index_0);
      if (  it0 == substitutions.end() )
         return;
      for( const auto& subs_vars : it0->second )
      {
         if(sign > 0)
            it1->second.push_back(subs_vars);
         else
         {
            int v = convert_substitution_to_col( subs_vars );
            if( subs_vars >= 0 )
                v = -subs_vars - 1 ;
            it1->second.push_back( v );
         }
      }
   }
#endif

   void
   substitute( int col, REAL substitute_factor, int lhs_id, int rhs_id,
               const Problem<REAL>& currentProblem, int skip_row_id = UNKNOWN )
   {
      const ConstraintMatrix<REAL>& matrix =
          currentProblem.getConstraintMatrix();
      auto col_vec = matrix.getColumnCoefficients( col );

      for( int i = 0; i < col_vec.getLength(); i++ )
      {
         int row = col_vec.getIndices()[i];
         if( row == skip_row_id )
            continue;
         if( matrix.getRowFlags()[row].test( RowFlag::kRedundant ) )
         {
            assert( rhs_row_mapping[row] == UNKNOWN );
            assert( lhs_row_mapping[row] == UNKNOWN );
            continue;
         }
         REAL factor = col_vec.getValues()[i] * abs(scale_factor[row]);
         if( num.isIntegral( factor / substitute_factor ) )
         {
            int val = cast_to_long( factor / substitute_factor );
            if( !matrix.getRowFlags()[row].test( RowFlag::kRhsInf ) )
            {
               next_constraint_id++;
               assert( rhs_row_mapping[row] != UNKNOWN );
               if( substitute_factor * factor > 0 )
                  proof_out << POL << lhs_id << " " << val << " * "
                            << rhs_row_mapping[row] << " +\n";
               else
                  proof_out << POL << rhs_id << " " << abs( val ) << " * "
                            << rhs_row_mapping[row] << " +\n";
#if VERIPB_VERSION >= 2
               proof_out << MOVE_LAST_CONS_TO_CORE;
#endif
               proof_out << DELETE_CONS << rhs_row_mapping[row];
               rhs_row_mapping[row] = next_constraint_id;
#if VERIPB_VERSION >= 2
               proof_out << " ; ; begin \n\t";
               if( substitute_factor * factor > 0 )
                  proof_out << POL << " " << rhs_row_mapping[row] << " " << rhs_id <<  " " << cast_to_long( abs( val ) ) << " * + \n";
               else
                  proof_out << POL << " " << rhs_row_mapping[row] << " " << lhs_id << " " << cast_to_long( abs( val ) ) << " * + \n";

               proof_out << "end";
               next_constraint_id += 2;
#endif
               proof_out << "\n";
            }
            if( !matrix.getRowFlags()[row].test( RowFlag::kLhsInf ) )
            {
               next_constraint_id++;
               assert( lhs_row_mapping[row] != UNKNOWN );
               if( substitute_factor * factor > 0 )
                  proof_out << POL << rhs_id << " " << val << " * "
                            << lhs_row_mapping[row] << " +\n";
               else
                  proof_out << POL << lhs_id << " " << abs( val ) << " * "
                            << lhs_row_mapping[row] << " +\n";
#if VERIPB_VERSION >= 2
               proof_out << MOVE_LAST_CONS_TO_CORE;
#endif
               proof_out << DELETE_CONS << lhs_row_mapping[row];
               lhs_row_mapping[row] = next_constraint_id;
#if VERIPB_VERSION >= 2
               proof_out << " ; ; begin \n\t";
               if( substitute_factor * factor > 0 )
                  proof_out << POL << lhs_row_mapping[row] << " " << lhs_id << " " << cast_to_long( abs( val ) ) << " * + \n";
               else
                  proof_out << POL << lhs_row_mapping[row] << " " << rhs_id << " " << cast_to_long( abs( val ) ) << " * + \n";
               proof_out << "end";
               next_constraint_id += 2;
#endif
               proof_out << "\n";
            }
         }
         else if( num.isIntegral( substitute_factor / factor ) )
         {
            scale_factor[row] *=
                cast_to_long( abs( substitute_factor / factor ) );
            int val = abs( cast_to_long( substitute_factor / factor ) );
            assert( val > 0 );
            if( !matrix.getRowFlags()[row].test( RowFlag::kRhsInf ) )
            {
               next_constraint_id++;
               assert( rhs_row_mapping[row] != UNKNOWN );
               if( substitute_factor * factor > 0 )
                  proof_out << POL << rhs_row_mapping[row] << " " << val << " * " << lhs_id << " +\n";
               else
                  proof_out << POL << rhs_row_mapping[row] << " " << val << " * " << rhs_id << " +\n";
#if VERIPB_VERSION >= 2
               proof_out << MOVE_LAST_CONS_TO_CORE;
#endif
               proof_out << DELETE_CONS << rhs_row_mapping[row];
               rhs_row_mapping[row] = next_constraint_id;
#if VERIPB_VERSION >= 2
               proof_out << " ; ; begin \n\t";
               if( substitute_factor * factor > 0 )
                  proof_out << POL << rhs_row_mapping[row] << " " << rhs_id << " + "  << cast_to_long( abs( val ) ) << " d\n";
               else
                  proof_out << POL << rhs_row_mapping[row] << " " << lhs_id << " + " << cast_to_long( abs( val ) ) << " d\n";
               proof_out << "end";
               next_constraint_id += 2;
#endif
               proof_out << "\n";

            }
            if( !matrix.getRowFlags()[row].test( RowFlag::kLhsInf ) )
            {
               next_constraint_id++;
               assert( lhs_row_mapping[row] != UNKNOWN );
               if( substitute_factor * factor > 0 )
                  proof_out << POL << lhs_row_mapping[row] << " " << val
                            << " * " << rhs_id << " +\n";
               else
                  proof_out << POL << lhs_row_mapping[row] << " " << val
                            << " * " << lhs_id << " +\n";
#if VERIPB_VERSION >= 2
               proof_out << MOVE_LAST_CONS_TO_CORE;
#endif
               proof_out << DELETE_CONS << lhs_row_mapping[row];
               lhs_row_mapping[row] = next_constraint_id;
#if VERIPB_VERSION >= 2
               proof_out << " ; ; begin \n\t";
               if( substitute_factor * factor < 0 )
                  proof_out << POL << lhs_row_mapping[row] << " " << rhs_id << " + " << cast_to_long( abs( val ) ) << " d\n";
               else
                  proof_out << POL << lhs_row_mapping[row] << " " << lhs_id << " + " << cast_to_long( abs( val ) ) << " d\n";
               proof_out << "end";
               next_constraint_id += 2;
#endif
               proof_out << "\n";
            }
         }
         else
         {
            assert( num.isIntegral( substitute_factor ) );
            assert( num.isIntegral( factor ) );
            scale_factor[row] *= cast_to_long( abs( substitute_factor ) );
            int val = abs( cast_to_long( factor ) );
            int val2 = abs( cast_to_long( substitute_factor ) );

            if( !matrix.getRowFlags()[row].test( RowFlag::kRhsInf ) )
            {
               next_constraint_id++;
               assert( rhs_row_mapping[row] != UNKNOWN );
               if( substitute_factor * factor > 0 )
                  proof_out << POL << lhs_id << " " << val << " * "
                            << rhs_row_mapping[row] << " " << val2 << " * +\n";
               else
                  proof_out << POL << rhs_id << " " << val << " * "
                            << rhs_row_mapping[row] << " " << val2 << " * +\n";
#if VERIPB_VERSION >= 2
               proof_out << MOVE_LAST_CONS_TO_CORE;
#endif
               proof_out << DELETE_CONS << rhs_row_mapping[row];
               rhs_row_mapping[row] = next_constraint_id;
#if VERIPB_VERSION >= 2
               proof_out << " ; ; begin \n\t";
               if( substitute_factor * factor > 0 )
                  proof_out << POL << rhs_id << " " << cast_to_long( abs( val ) ) << " * " << rhs_row_mapping[row] << " + " << cast_to_long( abs( val2 ) ) << " d\n";
               else
                  proof_out << POL << lhs_id << " " << cast_to_long( abs( val ) ) << " * " << rhs_row_mapping[row] << " + " << cast_to_long( abs( val2 ) ) << " d\n";
               proof_out << "end";
               next_constraint_id += 2;
#endif
               proof_out << "\n";
            }
            if( !matrix.getRowFlags()[row].test( RowFlag::kLhsInf ) )
            {
               next_constraint_id++;
               assert( lhs_row_mapping[row] != UNKNOWN );
               if( substitute_factor * factor > 0 )
                  proof_out << POL << rhs_id << " " << val << " * "
                            << lhs_row_mapping[row] << " " << val2 << " * +\n";
               else
                  proof_out << POL << lhs_id << " " << val << " * "
                            << lhs_row_mapping[row] << " " << val2 << " * +\n";
#if VERIPB_VERSION >= 2
               proof_out << MOVE_LAST_CONS_TO_CORE;
#endif
               proof_out << DELETE_CONS << lhs_row_mapping[row];
               lhs_row_mapping[row] = next_constraint_id;
#if VERIPB_VERSION >= 2
               proof_out << " ; ; begin \n\t";
               if( substitute_factor * factor > 0 )
                  proof_out << POL << lhs_id << " " << cast_to_long( abs( val ) ) << " * " << lhs_row_mapping[row] << " + " << cast_to_long( abs( val2 ) ) << " d\n";
               else
                  proof_out << POL << rhs_id << " " << cast_to_long( abs( val ) ) << " * " << lhs_row_mapping[row] << " + " << cast_to_long( abs( val2 ) ) << " d\n";
               proof_out << "end";
               next_constraint_id += 2;
#endif
               proof_out << "\n";
            }
         }
      }

   };

   std::pair<REAL, REAL>
   sparsify_convert_scale_to_frac( int eqrow, int candrow, REAL scale,
                                   const ConstraintMatrix<REAL>& matrix ) const
   {
      auto data_eq_row = matrix.getRowCoefficients( eqrow );
      auto data_cand_row = matrix.getRowCoefficients( candrow );
      int counter_eq_row = 0;
      int counter_cand_row = 0;
      while( true )
      {
         assert( counter_eq_row < data_eq_row.getLength() );
         int col_index_eq = data_eq_row.getIndices()[counter_eq_row];
         if( counter_cand_row >= data_cand_row.getLength() )
         {
            REAL val_eq = data_eq_row.getValues()[counter_eq_row];
            return { val_eq, val_eq * -scale };
         }
         int col_index_cand = data_cand_row.getIndices()[counter_cand_row];
         if( col_index_eq == col_index_cand)
         {
            counter_eq_row++;
            counter_cand_row++;
         }
         else if( col_index_eq > col_index_cand )
            counter_cand_row++;
         else
         {
            REAL val_eq = data_eq_row.getValues()[counter_eq_row];
            return { val_eq, val_eq * -scale };
         }
      }
      assert(false);
      return {-1, -1};
   }

#ifdef VERIPB_DEBUG
   void
   verify_changed_row( int row, const Problem<REAL>& problem, const Vec<int>& var_mapping, const Vec<int>& dirty_row_states )
   {
      auto constraintMatrix = problem.getConstraintMatrix();
      auto data = constraintMatrix.getRowCoefficients( validate_row );
      const RowFlags& rflags = constraintMatrix.getRowFlags()[validate_row];
      const REAL lhs = constraintMatrix.getLeftHandSides()[validate_row];
      const REAL rhs = constraintMatrix.getRightHandSides()[validate_row];
      const Vec<String>& names = problem.getVariableNames();
      assert( rhs_row_mapping[row] != UNKNOWN ||
                lhs_row_mapping[row] != UNKNOWN );
      if( lhs_row_mapping[row] != UNKNOWN )
      {
         if (std::find(dirty_row_states.begin(), dirty_row_states.end(), row) != dirty_row_states.end())
            proof_out << COMMENT << " dirty row state ";
         proof_out << EQUAL_CHECK << lhs_row_mapping[row] << " ";
         int offset = 0;
         for( int i = 0; i < data.getLength(); i++ )
         {

            REAL value = data.getValues()[i];
            int index = data.getIndices()[i];
            auto found = changed_entries_during_current_tsxs.find( index );
            if( found != changed_entries_during_current_tsxs.end() )
            {
               value = found->second;
               if( value == 0 )
                  continue;
            }
            if(problem.getUpperBounds()[index] == problem.getLowerBounds()[index])
            {
               assert(num.isIntegral(problem.getUpperBounds()[index] * value * scale_factor[row]));
               offset -= num.round_to_int(problem.getUpperBounds()[index] * value * scale_factor[row]);
               continue;
            }
            if( i != 0 )
               proof_out << "+";
            proof_out << num.round_to_int( abs( value ) * scale_factor[row] )
                      << " ";
            if( value < 0 )
            {
               offset += num.round_to_int( value );
               proof_out << NEGATED;
            }
            assert( var_mapping.size() > index );
            proof_out << names[var_mapping[index]] << " ";
         }
         proof_out << ">= "
                   << num.round_to_int( lhs ) * scale_factor[row] +
                          abs( offset )
                   << ";\n";
      }
      if( rhs_row_mapping[row] != UNKNOWN )
      {
         if (std::find(dirty_row_states.begin(), dirty_row_states.end(), row) != dirty_row_states.end())
            proof_out << COMMENT << " dirty row state ";
         proof_out << EQUAL_CHECK << rhs_row_mapping[row] << " ";
         int offset = 0;
         for( int i = 0; i < data.getLength(); i++ )
         {
            REAL value = data.getValues()[i];
            int index = data.getIndices()[i];

            auto found = changed_entries_during_current_tsxs.find( index );
            if( found != changed_entries_during_current_tsxs.end() )
            {
               value = found->second;
               if( value == 0 )
                  continue;
            }
            if(problem.getUpperBounds()[index] == problem.getLowerBounds()[index])
            {
               assert(num.isIntegral(problem.getUpperBounds()[index] * value * scale_factor[row]));
               offset -= num.round_to_int(problem.getUpperBounds()[index] * value * scale_factor[row]);
               continue;
            }
            if( i != 0 )
               proof_out << "+";
            proof_out << num.round_to_int( abs( value ) * scale_factor[row] )
                      << " ";
            if( value > 0 )
            {
               offset += num.round_to_int( value );
               proof_out << NEGATED;
            }
            else
            assert( var_mapping.size() > index );
            proof_out << names[var_mapping[index]] << " ";
         }
         proof_out << ">= "
                   << abs( offset ) -
                          num.round_to_int( rhs ) * scale_factor[row]
                   << ";\n";
      }
   }
#endif

#if VERIPB_VERSION == 1
   void
   add_substitutions_fix_to_witness(const Vec<String>& names, int orig_col_1, bool var)
   {
      if(!is_optimization_problem )
         return;
      auto it = substitutions.find(orig_col_1);
      if(it != substitutions.end())
      {
         for( const auto& substituted_vars : it->second )
         {
            proof_out << " " << names[convert_substitution_to_col(substituted_vars)] << " -> ";
            if( substituted_vars > 0 )
            proof_out << ( var ? "1" : "0" );
            else
            proof_out << ( var ? "0" : "1" );
         }
      }
   }

   void
   add_substitutions_to_witness(const Vec<String>& names, int orig_col_1, int orig_col_2)
   {
      if(!is_optimization_problem )
         return;
      auto it1 = substitutions.find(orig_col_1);
      if(it1 != substitutions.end())
      {
         for( const auto& substituted_vars : it1->second )
         {
            proof_out << " " << names[convert_substitution_to_col(substituted_vars)] << " -> ";
            if( substituted_vars < 0 )
               proof_out << NEGATED;
            proof_out << names[orig_col_2];
         }
      }
      auto it2 = substitutions.find(orig_col_2);
      if(it2 != substitutions.end())
      {
         for( const auto& substituted_vars : it2->second )
         {
            proof_out << " " << names[convert_substitution_to_col(substituted_vars)] << " -> ";
            if( substituted_vars < 0 )
               proof_out << NEGATED;
            proof_out << names[orig_col_1];
         }
      }
   }
#endif

   int
   convert_substitution_to_col( const int cols ) const
   {
      if(cols >= 0)
         return cols;
      return abs(cols + 1);
   }

   void
   propagate_row( int row, int col, REAL val, bool is_lb,
                  const Problem<REAL>& problem, const Vec<int>& var_mapping )
   {
      proof_out << POL << " ";
      const Vec<String>& names = problem.getVariableNames();
      const SparseVectorView<REAL>& row_data = problem.getConstraintMatrix().getRowCoefficients( row );
      const REAL* values = row_data.getValues();
      const int* indices = row_data.getIndices();
      auto& col_flags = problem.getVariableDomains().flags;
      bool is_lhs = false;
      if(lhs_row_mapping[row] != UNKNOWN && rhs_row_mapping[row] != UNKNOWN)
      {    
         REAL col_coef = 0;
         for(int i = 0; i < row_data.getLength(); i++)
            if( indices[i] == col )
            {
               col_coef = values[i];
               break;
            }
         assert( col_coef != 0);
         is_lhs = (is_lb && col_coef > 0) || (!is_lb && col_coef < 0);
      }
      else if(lhs_row_mapping[row] != UNKNOWN)
      {
         assert(rhs_row_mapping[row] == UNKNOWN);
         is_lhs = true;
      }
      else
      {
         assert(lhs_row_mapping[row] == UNKNOWN);
         is_lhs = false;
      }
      if(is_lhs)
         proof_out << lhs_row_mapping[row];
      else
         proof_out << rhs_row_mapping[row];
      proof_out <<  " " ;
      REAL col_coef = 0;
      for(int i = 0; i < row_data.getLength(); i++)
      {
         int c = indices[i];
         if(c == col)
         {
            col_coef = values[i];
            continue;
         }
         if( col_flags[c].test(ColFlag::kFixed, ColFlag::kSubstituted) )
            continue ;
         assert(num.isIntegral( values[i]));
         if(!( (values[i] < 0 && is_lhs) || (values[i] > 0 && !is_lhs)) )
            proof_out << NEGATED;
         proof_out << names[var_mapping[c]] << " " << cast_to_long( abs( values[i] ) ) << " * + ";
      }
      assert(col_coef != 0);
      proof_out << cast_to_long( abs( col_coef ) ) << " d\n";
      assert(
          (col_coef > 0 && is_lb && lhs_row_mapping[row]!= UNKNOWN) ||
          (col_coef < 0 && is_lb && rhs_row_mapping[row]!= UNKNOWN) ||
          (col_coef > 0 && !is_lb && rhs_row_mapping[row]!= UNKNOWN) ||
          (col_coef < 0 && !is_lb && lhs_row_mapping[row]!= UNKNOWN)
                  );
#ifdef VERIPB_DEBUG
      int orig_col = var_mapping[col];
      if(is_lb)
         proof_out << EQUAL_CHECK << next_constraint_id <<  " +1 " << names[orig_col] << " >= 1 ;\n";
      else 
         proof_out << EQUAL_CHECK << next_constraint_id << " +1 " << NEGATED << names[orig_col] << " >= 1 ;\n";

#endif
   }
};

} // namespace papilo

#endif
