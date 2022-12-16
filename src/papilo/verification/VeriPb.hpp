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

#include "papilo/core/Problem.hpp"
#include "papilo/misc/Num.hpp"
#include "papilo/misc/Vec.hpp"
#include "papilo/misc/compress_vector.hpp"
#include "papilo/misc/fmt.hpp"
#include "papilo/verification/ArgumentType.hpp"
#include "papilo/verification/CertificateInterface.hpp"



namespace papilo
{

/// type to store necessary data for post solve
template <typename REAL>
class VeriPb : public CertificateInterface<REAL>
{
 public:
   unsigned int nRowsOriginal;


   /// mapping constraint ids from PaPILO to VeriPb
   /// since VeriPb only supports >= each equation is mapped to 2 constraints
   Vec<int> rhs_row_mapping;
   Vec<int> lhs_row_mapping;

   /// it may be necessary to scale constraints to secure positive integer cons
   Vec<int> scale_factor;

   /// this holds the id of the constraints of VeriPB
   int next_constraint_id = 0;

   Num<REAL> num;
   Message msg;

   VeriPb() = default;

   VeriPb( const Problem<REAL>& _problem, const Num<REAL>& _num,
           const Message& _msg )
       : num( _num ), msg( _msg ), nRowsOriginal(_problem.getNRows())
   {
      rhs_row_mapping.reserve( nRowsOriginal );
      lhs_row_mapping.reserve( nRowsOriginal );
      scale_factor.reserve( nRowsOriginal );

      for( unsigned int i = 0; i < nRowsOriginal; ++i )
      {
         scale_factor.push_back( 1 );
         if(!_problem.getRowFlags()[i].test(RowFlag::kLhsInf))
         {
            next_constraint_id++;
            lhs_row_mapping.push_back( next_constraint_id );
         }
         else
            lhs_row_mapping.push_back( -1 );
         if(!_problem.getRowFlags()[i].test(RowFlag::kRhsInf))
         {
            next_constraint_id++;
            rhs_row_mapping.push_back( next_constraint_id );
         }
         else
            rhs_row_mapping.push_back( -1 );
      }
      assert( rhs_row_mapping.size() == lhs_row_mapping.size() );
      assert( rhs_row_mapping.size() == nRowsOriginal );
   }

   void
   print_header()
   {
      msg.info( "pseudo-Boolean proof version 1.1\n" );
      msg.info( "f {}\n", next_constraint_id );
   };

   void
   change_upper_bound( REAL val, const String& name,
                       ArgumentType argument = ArgumentType::kPrimal )
   {
      next_constraint_id++;
      // VeriPb can only handle >= constraint and they must start with variables
      // -> invert variable
      assert( val == 0 );
      switch( argument )
      {
      case ArgumentType::kPrimal:
         msg.info( "rup 1 ~{} >= 1 ;\n", name, (int)val );
         break;
      case ArgumentType::kDual:
         msg.info( "red 1 ~{} >= 1 ; {} -> 0\n", name, (int)val, name );
         break;
      case ArgumentType::kSymmetry:
         assert( false );
         break;
      default:
         assert( false );
      }

   }


   void
   change_lower_bound(  REAL val, const String& name, ArgumentType argument = ArgumentType::kPrimal)
   {
      next_constraint_id++;
      assert( val == 1 );
      switch( argument )
      {
      case ArgumentType::kPrimal:
         msg.info( "rup 1 {} >= {} ;\n", name, (int) val );
         break;
      case ArgumentType::kDual:
         msg.info( "red 1 {} >= {} ; {} -> {}\n", name, (int) val, name,
                   (int) val );
         break;
      case ArgumentType::kSymmetry:
         assert( false );
         break;
      default:
         assert( false );
      }
   }



   void
   change_rhs( int row, REAL val, const SparseVectorView<REAL>& data,
               const Vec<String>& names, const Vec<int>& var_mapping )
   {
      next_constraint_id++;
      fmt::print( "rup" );
      for( int i = 0; i < data.getLength(); i++ )
      {
         fmt::print( " ~{} {}", names[var_mapping[data.getIndices()[i]]],
                     (int) ( ( -1 ) * data.getValues()[i] ) );
         if( i != data.getLength() - 1 )
            fmt::print( " +" );
      }

      fmt::print( " >= {};\n", (int)( val + 1 ) );
      rhs_row_mapping[row] = next_constraint_id;
   }

   void
   change_lhs( int row, REAL val, const SparseVectorView<REAL>& data,
               const Vec<String>& names, const Vec<int>& var_mapping )
   {
      next_constraint_id++;
      fmt::print( "rup" );
      for( int i = 0; i < data.getLength(); i++ )
      {
         fmt::print( " {} {}", names[var_mapping[data.getIndices()[i]]],
                     (int) data.getValues()[i] );
         if( i != data.getLength() - 1 )
            fmt::print( " +" );
      }
      fmt::print( " >= {};\n", (int)( val ) );
      rhs_row_mapping[row] = next_constraint_id;
   }

   void
   update_row( int row, int col, REAL new_val,  const SparseVectorView<REAL>& data,
               RowFlags& rflags, REAL lhs, REAL rhs,
               const Vec<String>& names, const Vec<int>& var_mapping )
   {
      if( !rflags.test( RowFlag::kLhsInf ) )
      {
         next_constraint_id++;
         fmt::print( "rup " );
         for( int i = 0; i < data.getLength(); i++ )
         {
            if(data.getIndices()[i] == col)
            {
               if(new_val == 0)
                  continue ;
               if( i != 0 )
                  fmt::print( " +" );
               fmt::print( "{} {}", (int)( new_val ), names[col] );
            }
            else
            {
               if( i != 0 )
                  fmt::print( " +" );
               fmt::print( "{} {}",
                           (int)( data.getValues()[i] ), names[var_mapping[data.getIndices()[i]]] );
            }
         }

         fmt::print( " >= {} ;\n", (int)lhs );
         fmt::print( "del id {}\n", lhs_row_mapping[row] );
         lhs_row_mapping[row] = next_constraint_id;
      }
      if( !rflags.test( RowFlag::kRhsInf ) )
      {
         next_constraint_id++;
         fmt::print( "rup" );
         for( int i = 0; i < data.getLength(); i++ )
         {
            if( data.getIndices()[i] == col )
            {
               if( new_val == 0 )
                  continue;
               if( i != 0 )
                  fmt::print( " +" );
               fmt::print( "{} ~{}", (int)( new_val ),  names[col] );
            }
            else
            {
               if( i != 0 )
                  fmt::print( " +" );
               fmt::print( "{} ~{}", (int)( ( -1 ) * data.getValues()[i] ), names[var_mapping[data.getIndices()[i]]] );
            }
         }

         fmt::print( " >= {} ;\n", (int)( rhs + 1 ) );
         fmt::print( "del id {}\n", rhs_row_mapping[row] );
         rhs_row_mapping[row] = next_constraint_id;
      }
   }

   void
   substitute( int col, int row, const Problem<REAL>& currentProblem, const Vec<String>& names, const Vec<int>& var_mapping )
   {
      const ConstraintMatrix<REAL>& matrix =
          currentProblem.getConstraintMatrix();
      auto col_vec = matrix.getColumnCoefficients( col );
      auto row_data = matrix.getRowCoefficients(row);

      REAL factor = 0;
      for( int i = 0; i < col_vec.getLength(); i++ )
      {
         if(col_vec.getIndices()[i] == row)
         {
            factor = col_vec.getValues()[i];
            break;
         }
      }
      assert( factor != 0);
      for( int i = 0; i < col_vec.getLength(); i++ )
      {
         if(col_vec.getIndices()[i] == row)
            continue ;
         REAL div = col_vec.getValues()[i];
         assert( div != 0 && factor != 0);
         int current_row = col_vec.getIndices()[i];
         auto data = matrix.getRowCoefficients(current_row);
         if( num.isIntegral( div * scale_factor[i] / factor / scale_factor[row] ) )
         {
            assert( div * scale_factor[i] / factor / scale_factor[row]  > 0 );
            if( !matrix.getRowFlags()[current_row].test( RowFlag::kRhsInf ) )
            {
               next_constraint_id++;
               msg.info( "pol {} {} * {} +\n",
                         (int)lhs_row_mapping[row], (int)( div / factor ),
                         (int)rhs_row_mapping[current_row] );
               msg.info( "del id {}\n", (int)rhs_row_mapping[current_row] );
               rhs_row_mapping[current_row] = next_constraint_id;
            }
            if( !matrix.getRowFlags()[current_row].test( RowFlag::kLhsInf ) )
            {
               next_constraint_id++;
               msg.info( "pol {} {} * {} +\n",
                         (int)rhs_row_mapping[row], (int)( div / factor ),
                         (int)lhs_row_mapping[current_row] );
               msg.info( "del id {}\n", (int)lhs_row_mapping[current_row] );
               lhs_row_mapping[current_row] = next_constraint_id;
            }
         }
         else if( num.isIntegral( factor * scale_factor[row] / div / scale_factor[i] ) )
         {
            assert( factor / div > 0 );
            if( !matrix.getRowFlags()[current_row].test( RowFlag::kRhsInf ) )
            {
               next_constraint_id++;
               msg.info( "pol {} {} * {} +\n",
                         (int)rhs_row_mapping[current_row],
                         (int)( factor / div ), (int)lhs_row_mapping[row] );
               msg.info( "del id {}\n", (int)rhs_row_mapping[current_row] );
               rhs_row_mapping[current_row] = next_constraint_id;
            }
            if( !matrix.getRowFlags()[current_row].test( RowFlag::kLhsInf ) )
            {
               next_constraint_id++;
               msg.info( "pol {} {} * {} +\n",
                         (int)lhs_row_mapping[current_row],
                         (int)( factor / div ), (int)rhs_row_mapping[row] );
               msg.info( "del id {}\n", (int)lhs_row_mapping[current_row] );
               lhs_row_mapping[current_row] = next_constraint_id;
            }
         }
         else
         {
            assert( num.isIntegral( factor ) );
            scale_factor[row] =  (int) factor;
            if( !matrix.getRowFlags()[current_row].test( RowFlag::kRhsInf ) )
            {
               next_constraint_id++;
               msg.info( "pol {} {} * {} {} * +\n",
                         (int)lhs_row_mapping[row], (int)( div ),
                         (int)rhs_row_mapping[current_row], int (factor) );
               msg.info( "del id {}\n", (int)rhs_row_mapping[current_row] );
               rhs_row_mapping[current_row] = next_constraint_id;
            }
            if( !matrix.getRowFlags()[current_row].test( RowFlag::kLhsInf ) )
            {
               next_constraint_id++;
               msg.info( "pol {} {} * {} {} * +\n",
                         (int)rhs_row_mapping[row], (int) ( div ),
                         (int)lhs_row_mapping[current_row], (int) (factor) );
               msg.info( "del id {}\n", (int)lhs_row_mapping[current_row] );
               lhs_row_mapping[current_row] = next_constraint_id;
            }
         }
      }
      assert( !matrix.getRowFlags()[row].test( RowFlag::kRhsInf ) );
      assert( !matrix.getRowFlags()[row].test( RowFlag::kLhsInf ) );
      msg.info( "del id {}\n", (int)rhs_row_mapping[row] );
      msg.info( "del id {}\n", (int)lhs_row_mapping[row] );
   };

   void
   mark_row_redundant( int row )
   {
      assert( lhs_row_mapping[row] != -1 || rhs_row_mapping[row] != -1 );
      if( lhs_row_mapping[row] != -1 )
      {
         msg.info( "del id {}\n", (int) lhs_row_mapping[row] );
         lhs_row_mapping[row] = -1;
      }
      if( rhs_row_mapping[row] != -1 )
      {
         msg.info( "del id {}\n", (int) rhs_row_mapping[row] );
         rhs_row_mapping[row] = -1;
      }
   }

   void
   compress( const Vec<int>& rowmapping, const Vec<int>& colmapping,
             bool full = false )
   {
#ifdef PAPILO_TBB
      tbb::parallel_invoke(
          [this, &rowmapping, full]() {
             // update information about rows that is stored by index
             compress_vector( rowmapping, lhs_row_mapping );
             if( full )
                lhs_row_mapping.shrink_to_fit();
          },
          [this, &rowmapping, full]() {
             // update information about rows that is stored by index
             compress_vector( rowmapping, scale_factor );
             if( full )
                scale_factor.shrink_to_fit();
          },
          [this, &rowmapping, full]() {
             // update information about rows that is stored by index
             compress_vector( rowmapping, rhs_row_mapping );
             if( full )
                rhs_row_mapping.shrink_to_fit();
          } );
#else
      compress_vector( rowmapping, lhs_row_mapping );
      compress_vector( rowmapping, rhs_row_mapping );
      if( full )
      {
         rhs_row_mapping.shrink_to_fit();
         lhs_row_mapping.shrink_to_fit();
      }
#endif
   }


};



} // namespace papilo

#endif
