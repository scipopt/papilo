/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*               This file is part of the program and library                */
/*    PaPILO --- Parallel Presolve for Integer and Linear Optimization       */
/*                                                                           */
/* Copyright (C) 2020-2022 Konrad-Zuse-Zentrum                               */
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

#ifndef _PAPILO_VERI_VERI_PB_HPP_
#define _PAPILO_VERI_VERI_PB_HPP_

#define VERIPB_DEBUG

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
static const char* const DELETE_CONS = "del id ";
static const char* const RUP = "rup ";
static const char* const RED = "red ";
static const char* const POL = "pol ";
static const char* const EQUAL_CHECK = "e ";
static const char* const SATURATION = "s";
static const char* const WEAKENING = "w";
static const int UNKNOWN = -1;

template <typename REAL>
class VeriPb : public CertificateInterface<REAL>
{
 public:
   unsigned int nRowsOriginal{};
   std::ofstream proof_out;

   /// mapping constraint from PaPILO to constraint ids from PaPILO to VeriPb
   Vec<int> rhs_row_mapping;
   Vec<int> lhs_row_mapping;

   /// PaPILO does not care about the integrality of the coefficient
   /// therefore store scale factors to ensure the integrality
   Vec<int> scale_factor;

   const Vec<int>&
   getRowScalingFactor() const
   {
      return scale_factor;
   }

   /// this holds the id of the next generated constraint of VeriPB
   int next_constraint_id = 0;

   //TODO: set Num
   Num<REAL> num;

   int skip_deleting_rhs_constraint_id = UNKNOWN;
   int skip_deleting_lhs_constraint_id = UNKNOWN;
   int skip_changing_rhs = UNKNOWN;
   int skip_changing_lhs = UNKNOWN;

   // TODO: maybe this is not sufficient and matrix buffer needs to be handed over
   HashMap<int, int> changed_entries {};



   VeriPb() = default;

   VeriPb( const Problem<REAL>& _problem, const Num<REAL>& _num )
       : num( _num ), nRowsOriginal( _problem.getNRows() )
   {
      rhs_row_mapping.reserve( nRowsOriginal );
      lhs_row_mapping.reserve( nRowsOriginal );
      scale_factor.reserve( nRowsOriginal );

      for( unsigned int i = 0; i < nRowsOriginal; ++i )
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
      assert( rhs_row_mapping.size() == lhs_row_mapping.size() );
      assert( rhs_row_mapping.size() == nRowsOriginal );

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
      proof_out = std::ofstream(
             problem_name.substr( 0, length - ending ) + ".pbp" );
   }

   void
   print_header()
   {
      proof_out << "pseudo-Boolean proof version 1.0\n";
      proof_out << COMMENT << "Log files generated by PaPILO\n";
      proof_out << COMMENT << "Be aware that this is currently an experimental feature\n";
      proof_out << "f " << next_constraint_id << "\n";
      proof_out << std::fixed;
   };

   void
   start_transaction() {
       skip_changing_lhs = UNKNOWN;
       skip_changing_rhs = UNKNOWN;
       skip_deleting_lhs_constraint_id = UNKNOWN;
       skip_deleting_rhs_constraint_id = UNKNOWN;
       changed_entries.clear();
   };

   void
   end_transaction() {};

   void
   flush()
   {
      proof_out.flush();
   };

   void
   change_upper_bound( REAL val, int col, const Problem<REAL>& problem, const Vec<int>& var_mapping,
                       ArgumentType argument = ArgumentType::kPrimal)
   {
      next_constraint_id++;
      assert( val == 0 );
      const Vec<String>& names = problem.getVariableNames();
      int orig_col = var_mapping[col];
      switch( argument )
      {
      case ArgumentType::kPrimal:
         proof_out << RUP << "1 " << NEGATED << names[orig_col] << " >= 1 ;\n";
         break;
      case ArgumentType::kAggregation:
      case ArgumentType::kDual:
      case ArgumentType::kSymmetry:
         proof_out << RED << "1 " << NEGATED << names[orig_col] << " >= 1 ; " << names[orig_col] << " -> 0\n";
         break;
      default:
         assert( false );
         return;
      }
      int cons_id_fixing = next_constraint_id;

      auto col_coeff = problem.getConstraintMatrix().getColumnCoefficients( col );
      for( int row_index = 0; row_index < col_coeff.getLength(); row_index++ )
      {
         int row = col_coeff.getIndices()[row_index];
         if( problem.getRowFlags()[row].test( RowFlag::kRedundant ) )
         {
            assert( rhs_row_mapping[row] == UNKNOWN );
            assert( lhs_row_mapping[row] == UNKNOWN );
            continue;
         }
         assert(num.isIntegral(col_coeff.getValues()[row_index] * scale_factor[row]));
         int row_value = num.round_to_int(col_coeff.getValues()[row_index] * scale_factor[row]);
         if( !problem.getRowFlags()[row].test( RowFlag::kLhsInf ) )
         {
            next_constraint_id++;
            assert(lhs_row_mapping[row] != UNKNOWN);
            int cons = cons_id_fixing;
            if( row_value < 0)
               proof_out << POL << lhs_row_mapping[row] << " " << names[orig_col] << " " << abs(row_value) << " * + \n";
            else
               proof_out << POL << lhs_row_mapping[row] << " " << cons_id_fixing << " " << abs(row_value) << " * + \n";
            proof_out << DELETE_CONS << lhs_row_mapping[row] << "\n";
            lhs_row_mapping[row] = next_constraint_id;
         }
         if( !problem.getRowFlags()[row].test( RowFlag::kRhsInf ) )
         {
            next_constraint_id++;
            assert(rhs_row_mapping[row] != UNKNOWN);
            int cons = cons_id_fixing;
            if( row_value > 0)
               proof_out << POL << rhs_row_mapping[row] << " " << names[orig_col] << " " << abs(row_value) << " * + \n";
            else
               proof_out << POL << rhs_row_mapping[row] << " " << cons_id_fixing << " " << abs(row_value) << " * + \n";
            proof_out << DELETE_CONS << rhs_row_mapping[row] << "\n";
            rhs_row_mapping[row] = next_constraint_id;
         }
      }
   }

   void
   change_lower_bound( REAL val, int col, const Problem<REAL>& problem,
                       const Vec<int>& var_mapping,
                       ArgumentType argument = ArgumentType::kPrimal )
   {
      next_constraint_id++;
      assert( val == 1 );
      const Vec<String>& names = problem.getVariableNames();
      int orig_col = var_mapping[col];
      switch( argument )
      {
      case ArgumentType::kPrimal:
         proof_out << RUP << "1 " << names[orig_col]
                   << " >= " << num.round_to_int( val ) << " ;\n";
         break;
      case ArgumentType::kAggregation:
      case ArgumentType::kDual:
      case ArgumentType::kSymmetry:
         proof_out << RED << "1 " << names[orig_col]
                   << " >= " << num.round_to_int( val ) << " ; "
                   << names[orig_col] << " -> " << num.round_to_int( val )
                   << "\n";
         break;
      default:
         assert( false );
         return;
      }
      int cons_id_fixing = next_constraint_id;
      auto col_coeff = problem.getConstraintMatrix().getColumnCoefficients( col );
      for( int row_index = 0; row_index < col_coeff.getLength(); row_index++ )
      {
         int row = col_coeff.getIndices()[row_index];
         assert(num.isIntegral(col_coeff.getValues()[row_index] * scale_factor[row]));
         int row_value = num.round_to_int(col_coeff.getValues()[row_index] * scale_factor[row]);
         if( !problem.getRowFlags()[row].test( RowFlag::kLhsInf ) )
         {
            next_constraint_id++;
            assert(lhs_row_mapping[row] != UNKNOWN);
            int cons = cons_id_fixing;
            if( row_value > 0)
               proof_out << POL << lhs_row_mapping[row] << " " << NEGATED << names[orig_col] << " " << abs(row_value) << " * + \n";
            else
               proof_out << POL << lhs_row_mapping[row] << " " << cons_id_fixing  << " " << abs(row_value) << " * + \n";
            proof_out << DELETE_CONS << lhs_row_mapping[row] << "\n";
            lhs_row_mapping[row] = next_constraint_id;
         }
         if( !problem.getRowFlags()[row].test( RowFlag::kRhsInf ) )
         {
            next_constraint_id++;
            assert(rhs_row_mapping[row] != UNKNOWN);
            if( row_value < 0)
               proof_out << POL << rhs_row_mapping[row] << " " << NEGATED  << names[orig_col] << " " << abs(row_value) << " * + \n";
            else
               proof_out << POL << rhs_row_mapping[row] << " " << cons_id_fixing << " " << abs(row_value) << " * + \n";
            proof_out << DELETE_CONS << rhs_row_mapping[row] << "\n";
            rhs_row_mapping[row] = next_constraint_id;
         }
      }
   }

   void
   dominating_columns( int dominating_column, int dominated_column,
                       const Vec<String>& names, const Vec<int>& var_mapping )
   {
      next_constraint_id ++;
      auto name_dominating = names[var_mapping[dominating_column]];
      auto name_dominated = names[var_mapping[dominated_column]];
      proof_out << RED << "1 " << name_dominating << " +1 " << NEGATED << name_dominated
                << " >= 1 ; " << name_dominating << " -> " << name_dominated
                << " " << name_dominated << " -> " << name_dominating << "\n";
   }

   void
   change_rhs( int row, REAL val, const SparseVectorView<REAL>& data,
               const Vec<String>& names, const Vec<int>& var_mapping )
   {
      if(skip_changing_rhs == row)
      {
         skip_changing_rhs = UNKNOWN;
         return ;
      }
      assert( num.isIntegral( val * scale_factor[row] ) );
      next_constraint_id++;
      proof_out << RUP;
      int offset = 0;
      for( int i = 0; i < data.getLength(); i++ )
      {
         int unscaled_coeff = num.round_to_int( data.getValues()[i] );
         assert( unscaled_coeff != 0 );
         auto found = changed_entries.find( data.getIndices()[i] );
         if( found != changed_entries.end() )
         {
            unscaled_coeff = found->second;
            if(unscaled_coeff == 0)
               continue ;
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
      proof_out << " >=  " << abs(offset) - num.round_to_int( val ) * scale_factor[row]
                << ";\n";
      rhs_row_mapping[row] = next_constraint_id;
   }

   void
   change_lhs( int row, REAL val, const SparseVectorView<REAL>& data,
               const Vec<String>& names, const Vec<int>& var_mapping )
   {
      if(skip_changing_lhs == row)
      {
         skip_changing_lhs = UNKNOWN;
         return ;
      }
      assert( num.isIntegral( val * scale_factor[row] ) );
      next_constraint_id++;
      proof_out << RUP;
      int offset = 0;
      for( int i = 0; i < data.getLength(); i++ )
      {
         int unscaled_coeff = num.round_to_int( data.getValues()[i] );
         assert( unscaled_coeff != 0 );
         auto found = changed_entries.find( data.getIndices()[i] );
         if( found != changed_entries.end() )
         {
            unscaled_coeff = found->second;
            if(unscaled_coeff == 0)
               continue ;
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

      proof_out << " >=  " << num.round_to_int( val ) * scale_factor[row] + abs(offset)
                << ";\n";

      lhs_row_mapping[row] = next_constraint_id;
   }

   void
   change_rhs_parallel_row( int row, REAL val, int parallel_row,  const Problem<REAL>& problem, const Vec<int>& var_mapping)
   {
      REAL factor_row = problem.getConstraintMatrix()
                            .getRowCoefficients( row )
                            .getValues()[0] * scale_factor[row];
      REAL factor_parallel = problem.getConstraintMatrix()
                                 .getRowCoefficients( parallel_row )
                                 .getValues()[0] * scale_factor[parallel_row];
      REAL factor = factor_row / factor_parallel;
      assert( abs( factor ) >= 1 );
      proof_out << COMMENT ;
      if( factor < 0 )
         proof_out << lhs_row_mapping[parallel_row] ;
      else
         proof_out << rhs_row_mapping[parallel_row] ;
      proof_out << " is parallel to " << rhs_row_mapping[row] << "/"
                << lhs_row_mapping[row] << " are parallel.\n";

      if( abs(factor) == 1 )
      {
         // rhs_row_mapping[row] can be KNOWN for example if Singleton relaxed a constraint.
         if(rhs_row_mapping[row] != UNKNOWN)
            proof_out << DELETE_CONS << rhs_row_mapping[row] << "\n";
         if(factor == 1)
            rhs_row_mapping[row] = rhs_row_mapping[parallel_row];
         else
            rhs_row_mapping[row] = lhs_row_mapping[parallel_row];
         skip_deleting_lhs_constraint_id = rhs_row_mapping[row];
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
            assert(rhs_row_mapping[parallel_row] != UNKNOWN);
            next_constraint_id++;
            proof_out << POL << rhs_row_mapping[parallel_row] << " " << factor
                      << " *\n";
            if( rhs_row_mapping[row] != UNKNOWN )
               proof_out << DELETE_CONS << rhs_row_mapping[row] << "\n";
            rhs_row_mapping[row] = next_constraint_id;
            //scale also lhs
            if( lhs_row_mapping[row] != UNKNOWN && is_not_integral )
            {
               next_constraint_id++;
               proof_out << POL << lhs_row_mapping[row] << " " << factor_parallel
                         << " *\n";
               proof_out << DELETE_CONS << lhs_row_mapping[row] << "\n";
               lhs_row_mapping[row] = next_constraint_id;
               scale_factor[row] *= num.round_to_int(factor_parallel);
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
            assert(lhs_row_mapping[parallel_row] != UNKNOWN);
            proof_out << POL << lhs_row_mapping[parallel_row] << " " << abs(factor)
                      << " *\n";
            if( rhs_row_mapping[row] != UNKNOWN )
               proof_out << DELETE_CONS << rhs_row_mapping[row] << "\n";
            rhs_row_mapping[row] = next_constraint_id;
            //scale also lhs
            if( lhs_row_mapping[row] != UNKNOWN && is_not_integral )
            {
               next_constraint_id++;
               proof_out << POL << lhs_row_mapping[row] << " " << abs(factor_parallel)
                         << " *\n";
               proof_out << DELETE_CONS << lhs_row_mapping[row] << "\n";
               lhs_row_mapping[row] = next_constraint_id;
               scale_factor[row] *= num.round_to_int(abs(factor_parallel));
            }
         }
      }

   }

   void
   change_lhs_parallel_row( int row, REAL val, int parallel_row,  const Problem<REAL>& problem)
   {
      REAL factor_row = problem.getConstraintMatrix()
                                .getRowCoefficients( row )
                                .getValues()[0] * scale_factor[row];
      REAL factor_parallel = problem.getConstraintMatrix()
                                 .getRowCoefficients( parallel_row )
                                 .getValues()[0] * scale_factor[parallel_row];
      REAL factor = factor_row / factor_parallel;
      assert( abs( factor ) >= 1 );
      proof_out << COMMENT ;
      if( factor > 0 )
         proof_out << lhs_row_mapping[parallel_row] ;
      else
         proof_out << rhs_row_mapping[parallel_row] ;
      proof_out << " is parallel to " << rhs_row_mapping[row] << "/"
                                   << lhs_row_mapping[row] << " are parallel.\n";

      // shift the constraint ids
      if( abs(factor) == 1 )
      {
         // lhs_row_mapping[row] can be KNOWN for example if Singleton relaxed a constraint.
         if(lhs_row_mapping[row] != UNKNOWN)
            proof_out << DELETE_CONS << lhs_row_mapping[row] << "\n";
         assert( ( lhs_row_mapping[row] != UNKNOWN && factor == 1 ) ||
                 ( factor == -1 && rhs_row_mapping[row] != UNKNOWN ) );
         assert( factor == 1 || factor == -1 );
         if(factor == 1)
            lhs_row_mapping[row] = lhs_row_mapping[parallel_row];
         else
            lhs_row_mapping[row] = rhs_row_mapping[parallel_row];
         skip_deleting_lhs_constraint_id = lhs_row_mapping[row];
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
            proof_out << POL << lhs_row_mapping[parallel_row] << " " << factor
                      << " *\n";
            if( lhs_row_mapping[row] != UNKNOWN )
               proof_out << DELETE_CONS << lhs_row_mapping[row] << "\n";
            lhs_row_mapping[row] = next_constraint_id;
            //scale also rhs
            if( rhs_row_mapping[row] != UNKNOWN && is_not_integral )
            {
               next_constraint_id++;
               proof_out << POL << rhs_row_mapping[row] << " " << factor_parallel
                         << " *\n";
               proof_out << DELETE_CONS << rhs_row_mapping[row] << "\n";
               rhs_row_mapping[row] = next_constraint_id;
               scale_factor[row] *= num.round_to_int(factor_parallel);
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
            assert(rhs_row_mapping[parallel_row] != UNKNOWN);
            proof_out << POL << rhs_row_mapping[parallel_row] << " " << abs(factor)
                      << " *\n";
            if( lhs_row_mapping[row] != UNKNOWN )
               proof_out << DELETE_CONS << lhs_row_mapping[row] << "\n";
            lhs_row_mapping[row] = next_constraint_id;
            //scale also rhs
            if( rhs_row_mapping[row] != UNKNOWN && is_not_integral )
            {
               next_constraint_id++;
               proof_out << POL << rhs_row_mapping[row] << " " << abs(factor_parallel)
                         << " *\n";
               proof_out << DELETE_CONS << rhs_row_mapping[row] << "\n";
               rhs_row_mapping[row] = next_constraint_id;
               scale_factor[row] *= num.round_to_int(abs(factor_parallel));
            }
         }
      }
}

   void
   change_lhs_inf( int row )
   {
      proof_out << DELETE_CONS << lhs_row_mapping[row] << "\n";
      lhs_row_mapping[row] = UNKNOWN;
   }

   void
   change_rhs_inf( int row )
   {
      proof_out << DELETE_CONS << rhs_row_mapping[row] << "\n";
      rhs_row_mapping[row] = UNKNOWN;
   }

   void
   change_matrix_entry( int row, int col, REAL new_val,
               const SparseVectorView<REAL>& data, RowFlags& rflags, REAL lhs,
               REAL rhs, const Vec<String>& names, const Vec<int>& var_mapping, ArgumentType argument )
   {
      // remove singleton variable from equation
      changed_entries.emplace( col, num.round_to_int( new_val ) );
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
                   num.round_to_int( data.getValues()[i] * scale_factor[row] );
            }
         assert( old_value != 0 );

         auto name = names[var_mapping[col]];
         assert( num.isIntegral( new_val ) );
         int diff = old_value - num.round_to_int( new_val );
         if( !rflags.test( RowFlag::kLhsInf ) )
         {
            next_constraint_id++;
            assert( lhs_row_mapping[row] != UNKNOWN );
            if( old_value > 0 )
               proof_out << POL << lhs_row_mapping[row] << " " << NEGATED  << name << " "
                         << abs( diff ) << " * + \n";
            else
               proof_out << POL << lhs_row_mapping[row] << " " << name << " "
                         << abs( diff ) << " * + \n";
            skip_changing_lhs = row;
            proof_out << DELETE_CONS << lhs_row_mapping[row] << "\n";
            lhs_row_mapping[row] = next_constraint_id;
         }
         if( !rflags.test( RowFlag::kRhsInf ) )
         {
            next_constraint_id++;
            assert( rhs_row_mapping[row] != UNKNOWN );
            skip_changing_rhs = row;
            if( old_value < 0 )
            {
               proof_out << POL << rhs_row_mapping[row] << " " << NEGATED  << name << " "
                         << abs( diff ) << " * + \n";
            }
            else
            {
               proof_out << POL << rhs_row_mapping[row] << " " << name << " "
                         << abs( diff ) << " * + \n";
            }
            proof_out << DELETE_CONS << rhs_row_mapping[row] << "\n";
            rhs_row_mapping[row] = next_constraint_id;
         }

         return;
      }
      else if( argument == ArgumentType::kSaturation )
      {
         next_constraint_id++;
         assert(!rflags.test( RowFlag::kEquation ));
         assert(rflags.test( RowFlag::kRhsInf ) || rflags.test( RowFlag::kLhsInf ));
         proof_out << POL;
         if(!rflags.test( RowFlag::kRhsInf ))
         {
            assert(rhs_row_mapping[row] != UNKNOWN);
            proof_out << rhs_row_mapping[row] << " ";
            skip_changing_rhs = next_constraint_id;
         }
         else
         {
            assert(lhs_row_mapping[row] != UNKNOWN);
            proof_out << lhs_row_mapping[row] << " ";
            skip_changing_lhs = next_constraint_id;
         }
         proof_out << SATURATION << "\n";
         if( !rflags.test( RowFlag::kRhsInf ) )
         {
            proof_out << DELETE_CONS << rhs_row_mapping[row] << "\n";
            rhs_row_mapping[row] = next_constraint_id;
         }
         else
         {
            proof_out << DELETE_CONS << lhs_row_mapping[row] << "\n";
            lhs_row_mapping[row] = next_constraint_id;
         }
#ifdef VERIPB_DEBUG
         verify_changed_row( row, col, new_val, data, rflags, lhs, rhs,
                             names, var_mapping );
#endif
      }
      else if( argument == ArgumentType::kWeakening )
      {
         assert( new_val == 0);
         assert( rflags.test(RowFlag::kRhsInf) || rflags.test(RowFlag::kLhsInf));
         next_constraint_id++;
         proof_out << POL;
         if(rhs_row_mapping[row] != UNKNOWN)
         {
            assert(!rflags.test(RowFlag::kRhsInf));
            proof_out << rhs_row_mapping[row] << " ";
         }
         else
         {
            assert( !rflags.test( RowFlag::kLhsInf ) );
            proof_out << lhs_row_mapping[row] << " ";
         }
         proof_out << names[var_mapping[col]] << " " << WEAKENING << "\n";
         if(rhs_row_mapping[row] != UNKNOWN)
         {
            proof_out << DELETE_CONS << rhs_row_mapping[row] << "\n";
            rhs_row_mapping[row] = next_constraint_id;         }
         else
         {
            proof_out << DELETE_CONS << lhs_row_mapping[row] << "\n";
            lhs_row_mapping[row] = next_constraint_id;
         }
#ifdef VERIPB_DEBUG
         verify_changed_row( row, col, new_val, data, rflags, lhs, rhs,
                             names, var_mapping );
#endif
      }
      else
      {
         assert( false );
      }

//      assert( num.isIntegral( new_val * scale_factor[row] ) );
//      if( !rflags.test( RowFlag::kLhsInf ) )
//      {
//         next_constraint_id++;
//         proof_out << RUP;
//         int offset = 0;
//         for( int i = 0; i < data.getLength(); i++ )
//         {
//
//            int index = data.getIndices()[i];
//            if( index == col )
//            {
//               if( new_val == 0 )
//                  continue;
//               if( i != 0 )
//                  proof_out << " +";
//               proof_out << abs( num.round(new_val) * scale_factor[row] ) << " ";
//
//               if( new_val < 0 )
//               {
//                  proof_out << NEGATED;
//                  offset -= num.round_to_int(new_val);
//               }
//            }
//            else
//            {
//               if( i != 0 )
//                  proof_out << " +";
//               int val = num.round_to_int( data.getValues()[i] * scale_factor[row] );
//               proof_out << abs( val ) << " ";
//               if( val < 0 )
//               {
//                  proof_out << NEGATED;
//                  offset -= val;
//               }
//            }
//            proof_out << names[var_mapping[index]];
//         }
//         proof_out << " >= " << offset + num.round_to_int( lhs * scale_factor[row] )
//                   << ";\n";
//
//         proof_out << DELETE_CONS << lhs_row_mapping[row] << "\n";
//
//         lhs_row_mapping[row] = next_constraint_id;
//      }
//      if( !rflags.test( RowFlag::kRhsInf ) )
//      {
//         next_constraint_id++;
//         int offset = 0;
//         proof_out<< RUP ;
//         for( int i = 0; i < data.getLength(); i++ )
//         {
//            int index = data.getIndices()[i];
//            if( index == col )
//            {
//               if( new_val == 0 )
//                  continue;
//               if( i != 0 )
//                  proof_out << " +";
//               proof_out << abs( num.round_to_int(new_val) * scale_factor[row] ) << " ";
//
//               if( new_val > 0 )
//               {
//                  offset += num.round_to_int(new_val);
//                  proof_out << NEGATED;
//               }
//            }
//            else
//            {
//               if( i != 0 )
//                  proof_out << " +";
//               int val = num.round_to_int( data.getValues()[i] * scale_factor[row] );
//               proof_out << abs( val ) << " ";
//               if( val > 0 )
//               {
//                  offset += val;
//                  proof_out << NEGATED;
//               }
//            }
//            proof_out << names[var_mapping[index]];
//         }
//
//         proof_out << " >= " << offset - num.round_to_int( rhs * scale_factor[row] )
//                   << ";\n";
//         proof_out << DELETE_CONS << rhs_row_mapping[row] << "\n";
//         rhs_row_mapping[row] = next_constraint_id;
//      }
   }

   void
   sparsify( int eqrow, int candrow, REAL scale,
             const Problem<REAL>& currentProblem )
   {

      const ConstraintMatrix<REAL>& matrix =
          currentProblem.getConstraintMatrix();
      int scale_eqrow = scale_factor[eqrow];
      int scale_candrow = scale_factor[candrow];
      assert( scale != 0 );
      REAL scale_updated = scale * scale_candrow / scale_eqrow;
      if( num.isIntegral( scale_updated ) )
      {
         int int_scale_updated = num.round_to_int(scale_updated);
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
            proof_out << DELETE_CONS << rhs_row_mapping[candrow] << "\n";
            rhs_row_mapping[candrow] = next_constraint_id;
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
            proof_out << DELETE_CONS << lhs_row_mapping[candrow] << "\n";

            lhs_row_mapping[candrow] = next_constraint_id;
         }
      }
      else if( num.isIntegral( 1.0 / scale_updated ) )
      {
         int int_scale_updated = num.round_to_int( 1.0 / scale_updated );

         if( !matrix.getRowFlags()[candrow].test( RowFlag::kRhsInf ) )
         {
            next_constraint_id++;
            assert( rhs_row_mapping[candrow] != UNKNOWN );
            assert( rhs_row_mapping[eqrow] != UNKNOWN );
            if( int_scale_updated > 0 )
               proof_out << POL << rhs_row_mapping[candrow] << " "
                         << abs( int_scale_updated ) << " * "
                         << rhs_row_mapping[candrow] << " +\n";
            else
               proof_out << POL << rhs_row_mapping[candrow] << " "
                         << abs( int_scale_updated ) << " * "
                         << lhs_row_mapping[candrow] << " +\n";
            proof_out << DELETE_CONS << rhs_row_mapping[candrow] << "\n";
            rhs_row_mapping[candrow] = next_constraint_id;
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

            proof_out << DELETE_CONS << lhs_row_mapping[candrow] << "\n";

            lhs_row_mapping[candrow] = next_constraint_id;
         }
         scale_factor[candrow] *= int_scale_updated;
      }
      else
      {
         auto frac =
             sparsify_convert_scale_to_frac( eqrow, candrow, matrix );
         assert( frac.second / frac.first == -scale );
         int frac_eqrow = abs( num.round_to_int( frac.second * scale_candrow ) );
         int frac_candrow = abs( num.round_to_int( frac.first * scale_eqrow ) );

         if( !matrix.getRowFlags()[candrow].test( RowFlag::kRhsInf ) )
         {
            next_constraint_id++;
            assert( rhs_row_mapping[candrow] != UNKNOWN );
            assert( rhs_row_mapping[eqrow] != UNKNOWN );
            if( scale > 0 )
               proof_out << POL << rhs_row_mapping[candrow] << " "
                         << frac_candrow << " * " << rhs_row_mapping[eqrow]
                         << " * " << frac_eqrow << " +\n";
            else
               proof_out << POL << rhs_row_mapping[candrow] << " "
                         << frac_candrow << " * " << lhs_row_mapping[eqrow]
                         << " * " << frac_eqrow << " +\n";
            proof_out << DELETE_CONS << rhs_row_mapping[candrow] << "\n";
            rhs_row_mapping[candrow] = next_constraint_id;
         }
         if( !matrix.getRowFlags()[candrow].test( RowFlag::kLhsInf ) )
         {
            next_constraint_id++;
            assert( lhs_row_mapping[candrow] != UNKNOWN );
            assert( lhs_row_mapping[eqrow] != UNKNOWN );
            if( scale > 0 )
               proof_out << POL << lhs_row_mapping[candrow] << " "
                         << frac_candrow << " * " << lhs_row_mapping[eqrow]
                         << " * " << frac_eqrow << " +\n";
            else
               proof_out << POL << lhs_row_mapping[candrow] << " "
                         << frac_candrow << " * " << rhs_row_mapping[eqrow]
                         << " * " << frac_eqrow << " +\n";
            proof_out << DELETE_CONS << lhs_row_mapping[candrow] << "\n";

            lhs_row_mapping[candrow] = next_constraint_id;
         }
         scale_factor[candrow] *= frac_candrow;
      }
   }

   void
   substitute( int col, const SparseVectorView<REAL>& equality, REAL offset,
               const Problem<REAL>& currentProblem, const Vec<String>& names,
               const Vec<int>& var_mapping )
   {
      assert( num.isIntegral( offset ) );
      const REAL* values = equality.getValues();
      const int* indices = equality.getIndices();
      assert( equality.getLength() == 2 );
      assert( num.isIntegral( values[0] ) && num.isIntegral( values[1] ) );
      assert(num.round_to_int(values[0]) != 0);
      assert(num.round_to_int(values[1]) != 0);
      REAL substitute_factor = indices[0] == col ? values[0] : values[1];

      next_constraint_id++;
      int first_constraint_id = next_constraint_id;
      proof_out << COMMENT << "postsolve stack : row id " << next_constraint_id << "\n";
      proof_out << RUP;

      int lhs = num.round_to_int( offset);
      proof_out << abs( num.round_to_int( values[0]) ) << " ";
      if(values[0] < 0)
      {
         proof_out << NEGATED;
         lhs += abs(num.round_to_int( values[0]));
      }
      proof_out << names[var_mapping[indices[0]]] << " +" << abs( num.round_to_int( values[1] )) << " ";
      if(values[1] < 0)
      {
         proof_out << NEGATED;
         lhs += abs(num.round_to_int( values[1]));
      }
      proof_out << names[var_mapping[indices[1]]] << " >= " << lhs << ";\n";
      int lhs_id = next_constraint_id;

      next_constraint_id++;
      int second_constraint_id = next_constraint_id;

      proof_out << COMMENT << "postsolve stack : row id " << next_constraint_id << "\n";

      proof_out << RUP;
      int rhs = - num.round_to_int( offset);
      proof_out << abs( num.round_to_int(values[0]) ) << " ";
      if(values[0] > 0)
      {
         proof_out << NEGATED;
         rhs += abs(num.round_to_int( values[0]));
      }
      proof_out << names[var_mapping[indices[0]]] << " +" << abs( num.round_to_int( values[1] )) << " ";
      if(values[1] > 0)
      {
         proof_out << NEGATED;
         rhs += abs(num.round_to_int(values[1]));
      }
      proof_out << names[var_mapping[indices[1]]] << " >= " << rhs << ";\n";

      substitute( col, substitute_factor, lhs_id, next_constraint_id,
                  currentProblem );
      proof_out << DELETE_CONS << first_constraint_id << "\n";
      proof_out << DELETE_CONS << second_constraint_id << "\n";

   }

   void
   substitute( int col, int substituted_row,
               const Problem<REAL>& currentProblem )
   {
      const ConstraintMatrix<REAL>& matrix =
          currentProblem.getConstraintMatrix();
      auto col_vec = matrix.getColumnCoefficients( col );
      auto row_data = matrix.getRowCoefficients( substituted_row );

      // TODO: if singleton column is substiuted in objective -> continue there is nothing to do in this case
      if(col_vec.getLength() == 1)
         return ;

      REAL substitute_factor = 0;
      for( int i = 0; i < col_vec.getLength(); i++ )
      {
         if( col_vec.getIndices()[i] == substituted_row )
         {
            substitute_factor =
                col_vec.getValues()[i] * scale_factor[substituted_row];
            break;
         }
      }
      substitute( col, substitute_factor, lhs_row_mapping[substituted_row],
                  rhs_row_mapping[substituted_row], currentProblem,
                  substituted_row );
      assert( !matrix.getRowFlags()[substituted_row].test( RowFlag::kRhsInf ) );
      assert( !matrix.getRowFlags()[substituted_row].test( RowFlag::kLhsInf ) );
      proof_out << COMMENT << "postsolve stack : row id "
                << rhs_row_mapping[substituted_row] << "\n";
      proof_out << COMMENT << "postsolve stack : row id "
                << lhs_row_mapping[substituted_row] << "\n";
      proof_out << DELETE_CONS << rhs_row_mapping[substituted_row] << "\n";
      proof_out << DELETE_CONS << lhs_row_mapping[substituted_row] << "\n";
   };

   void
   mark_row_redundant( int row )
   {
      assert( lhs_row_mapping[row] != UNKNOWN || rhs_row_mapping[row] != UNKNOWN );
      if( lhs_row_mapping[row] != UNKNOWN )
      {
         if(lhs_row_mapping[row] == skip_deleting_lhs_constraint_id)
            skip_deleting_lhs_constraint_id = UNKNOWN;
         else
            proof_out << DELETE_CONS << lhs_row_mapping[row] << "\n";
         lhs_row_mapping[row] = UNKNOWN;

      }
      if( rhs_row_mapping[row] != UNKNOWN )
      {
         if(rhs_row_mapping[row] == skip_deleting_rhs_constraint_id)
            skip_deleting_rhs_constraint_id = UNKNOWN;
         else
            proof_out << DELETE_CONS << rhs_row_mapping[row] << "\n";
         rhs_row_mapping[row] = UNKNOWN;

      }
   }

   void
   log_solution( const Solution<REAL>& orig_solution, const Vec<String>& names )
   {
      proof_out << "o";
      next_constraint_id++;
      for( int i = 0; i < orig_solution.primal.size(); i++ )
      {
         assert( orig_solution.primal[i] == 0 || orig_solution.primal[i] == 1 );
         proof_out << " ";
         if( orig_solution.primal[i] == 0 )
            proof_out << NEGATED;
         proof_out << names[i];
      }
      next_constraint_id++;
      proof_out << "\n";
      proof_out << "u >= 1 ;\n";
      proof_out << "c " << next_constraint_id << "\n";
   };

   void
   compress( const Vec<int>& rowmapping, const Vec<int>& colmapping,
             bool full = false )
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
      if( full )
      {
         rhs_row_mapping.shrink_to_fit();
         lhs_row_mapping.shrink_to_fit();
         scale_factor.shrink_to_fit();
      }
#endif
   }

 private:
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
         REAL factor = col_vec.getValues()[i] * scale_factor[row];
         auto data = matrix.getRowCoefficients( row );
         if( num.isIntegral( factor / substitute_factor ) )
         {
            int val = num.round_to_int( factor / substitute_factor );
            if( !matrix.getRowFlags()[row].test( RowFlag::kRhsInf ) )
            {
               next_constraint_id++;
               assert( rhs_row_mapping[row] != UNKNOWN );
               if( substitute_factor * factor > 0 )
                  proof_out << POL << lhs_id << " " << val << " * "
                            << rhs_row_mapping[row] << " +\n";
               else
                  proof_out << POL << rhs_id << " " << abs(val) << " * "
                            << rhs_row_mapping[row] << " +\n";

               proof_out << DELETE_CONS << rhs_row_mapping[row] << "\n";
               rhs_row_mapping[row] = next_constraint_id;
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
               proof_out << DELETE_CONS << lhs_row_mapping[row] << "\n";
               lhs_row_mapping[row] = next_constraint_id;
            }
         }
         else if( num.isIntegral( substitute_factor / factor ) )
         {
            scale_factor[row] *= num.round_to_int( substitute_factor / factor );
            int val = abs( num.round_to_int( substitute_factor / factor ) );
            assert(val > 0);
            if( !matrix.getRowFlags()[row].test( RowFlag::kRhsInf ) )
            {
               next_constraint_id++;
               assert( rhs_row_mapping[row] != UNKNOWN );
               if( substitute_factor * factor > 0 )
                  proof_out << POL << rhs_row_mapping[row] << " " << val
                            << " * " << lhs_id << " +\n";
               else
                  proof_out << POL << rhs_row_mapping[row] << " " << val
                            << " * " << rhs_id << " +\n";

               proof_out << DELETE_CONS << rhs_row_mapping[row] << "\n";
               rhs_row_mapping[row] = next_constraint_id;
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
               proof_out << DELETE_CONS << lhs_row_mapping[row] << "\n";
               lhs_row_mapping[row] = next_constraint_id;
            }
         }
         else
         {
            assert( num.isIntegral( substitute_factor ) );
            assert( num.isIntegral( factor ) );
            scale_factor[row] *= num.round_to_int(substitute_factor);
            int val = abs( num.round_to_int(factor) );
            int val2 = abs( num.round_to_int(substitute_factor) );

            if( !matrix.getRowFlags()[row].test( RowFlag::kRhsInf ) )
            {
               next_constraint_id++;
               assert( rhs_row_mapping[row] != UNKNOWN );
               if( substitute_factor * factor > 0 )
                  proof_out << POL << lhs_id << " " << val << " * "
                            << rhs_row_mapping[row] << " " << val2 << " +\n";
               else
                  proof_out << POL << rhs_id << " " << val << " * "
                            << rhs_row_mapping[row] << " " << val2 << " +\n";
               proof_out << DELETE_CONS << rhs_row_mapping[row] << "\n";
               rhs_row_mapping[row] = next_constraint_id;
            }
            if( !matrix.getRowFlags()[row].test( RowFlag::kLhsInf ) )
            {
               next_constraint_id++;
               assert( lhs_row_mapping[row] != UNKNOWN );
               if( substitute_factor * factor > 0 )
                  proof_out << POL << rhs_id << " " << val << " * "
                         << lhs_row_mapping[row] << " " << val2 << " +\n";
               else
                  proof_out << POL << lhs_id << " " << val << " * "
                            << lhs_row_mapping[row] << " " << val2 << " +\n";

               proof_out << DELETE_CONS << lhs_row_mapping[row] << "\n";
               lhs_row_mapping[row] = next_constraint_id;
            }
         }
      }
   };

   std::pair<REAL, REAL>
   sparsify_convert_scale_to_frac( int eqrow, int candrow,
                                   const ConstraintMatrix<REAL>& matrix ) const
   {
      auto data_eq_row = matrix.getRowCoefficients( eqrow );
      auto data_cand_row = matrix.getRowCoefficients( candrow );
      int index_eq_row = 0;
      int index_cand_row = 0;
      while( true )
      {
         assert( index_eq_row < data_eq_row.getLength() );
         assert( index_cand_row < data_cand_row.getLength() );
         if( data_eq_row.getIndices()[index_eq_row] ==
             data_cand_row.getIndices()[index_cand_row] )
         {
            index_eq_row++;
            index_cand_row++;
            continue;
         }
         if( data_eq_row.getIndices()[index_eq_row] <
             data_cand_row.getIndices()[index_cand_row] )
            break;
         if( data_eq_row.getIndices()[index_cand_row] <
             data_cand_row.getIndices()[index_eq_row] )
            index_eq_row++;
      }
      return { data_eq_row.getValues()[index_eq_row],
               data_cand_row.getValues()[index_cand_row] };
   }

   void
   verify_changed_row( int row, int col, REAL new_val,
               const SparseVectorView<REAL>& data, RowFlags& rflags, REAL lhs,
               REAL rhs, const Vec<String>& names, const Vec<int>& var_mapping )
   {
      if( lhs_row_mapping[row] != UNKNOWN )
      {
         proof_out << EQUAL_CHECK << lhs_row_mapping[row] << " ";
         int offset = 0;
         for( int i = 0; i < data.getLength(); i++ )
         {

            REAL value = data.getValues()[i];
            auto found = changed_entries.find(data.getIndices()[i]);
            if( found != changed_entries.end() )
            {
               value = found->second;
               if(value == 0)
                  continue ;
            }
            if( i != 0 )
               proof_out << "+";
            proof_out << num.round_to_int(abs( value ) * scale_factor[row]) << " ";
            if( value < 0 )
            {
               offset += num.round_to_int( value );
               proof_out << NEGATED;
            }
            assert( var_mapping.size() > data.getIndices()[i] );
            proof_out << names[var_mapping[data.getIndices()[i]]] << " ";
         }
         proof_out << " >= " << num.round_to_int( lhs ) * scale_factor[row] + abs(offset) << ";\n";
      }
      if( rhs_row_mapping[row] != UNKNOWN )
      {
         proof_out << EQUAL_CHECK << rhs_row_mapping[row] << " ";
         int offset = 0;
         for( int i = 0; i < data.getLength(); i++ )
         {
            REAL value = data.getValues()[i];
            auto found = changed_entries.find(data.getIndices()[i]);
            {
               value = found->second;
               if(value == 0)
                  continue ;
            }
            if( i != 0 )
               proof_out << "+";
            proof_out << num.round_to_int(abs( value ) * scale_factor[row]) << " ";
            if( value < 0 )
            {
               offset += num.round_to_int( value );
               proof_out << NEGATED;
            }
            proof_out << names[var_mapping[data.getIndices()[i]]] << " ";
         }
         proof_out << " >= " << abs(offset) - num.round_to_int( rhs ) * scale_factor[row]
                   << ";\n";
      }
   }
};

} // namespace papilo

#endif
