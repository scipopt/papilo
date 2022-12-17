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

   // TODO: test change_rhs and change_lhs
   void
   change_rhs( int row, REAL val, const SparseVectorView<REAL>& data,
               const Vec<String>& names, const Vec<int>& var_mapping )
   {
      assert( num.isIntegral( val * scale_factor[row] ) );
      next_constraint_id++;
      fmt::print( "rup" );
      int offset = 0;
      for( int i = 0; i < data.getLength(); i++ )
      {
         int coeff = (int) data.getValues()[i] * scale_factor[row];
         const String& varname = names[var_mapping[data.getIndices()[i]]];
         assert( coeff != 0 );
         fmt::print( "{}", abs(coeff) );
         if( coeff < 0 )
            offset -= coeff;
         else
            fmt::print( " ~" );
         fmt::print( "{}",  varname );
         if( i != data.getLength() - 1 )
            fmt::print( " +" );
      }

      fmt::print( " >= {};\n", (int)( val * scale_factor[row] ) + offset );
      rhs_row_mapping[row] = next_constraint_id;
   }

   // TODO: test change_rhs and change_lhs
   void
   change_lhs( int row, REAL val, const SparseVectorView<REAL>& data,
               const Vec<String>& names, const Vec<int>& var_mapping )
   {
      assert( num.isIntegral( val * scale_factor[row] ) );
      next_constraint_id++;
      fmt::print( "rup" );
      int offset = 0;
      for( int i = 0; i < data.getLength(); i++ )
      {
         int coeff = (int) data.getValues()[i] * scale_factor[row];
         const String& varname = names[var_mapping[data.getIndices()[i]]];
         assert( coeff != 0 );
         fmt::print( "{} ", abs(coeff) );
         if( coeff < 0 )
         {
            fmt::print( "~" );
            offset -= coeff;
         }
         fmt::print( "{}",  varname );
         if( i != data.getLength() - 1 )
            fmt::print( " +" );
      }
      fmt::print( " >= {};\n", (int)( val * scale_factor[row]) + offset );
      rhs_row_mapping[row] = next_constraint_id;
   }

   //TODO test for negative values in the coeff
   void
   update_row( int row, int col, REAL new_val,  const SparseVectorView<REAL>& data,
               RowFlags& rflags, REAL lhs, REAL rhs,
               const Vec<String>& names, const Vec<int>& var_mapping )
   {
      assert( num.isIntegral( new_val * scale_factor[row] ) );
      if( !rflags.test( RowFlag::kLhsInf ) )
      {
         next_constraint_id++;
         fmt::print( "rup " );
         int offset = 0;
         for( int i = 0; i < data.getLength(); i++ )
         {

            if(data.getIndices()[i] == col)
            {
               if(new_val == 0)
                  continue ;
               if( i != 0 )
                  fmt::print( " +" );
               fmt::print( "{} ", abs((int) new_val * scale_factor[row]) );
               if( new_val < 0 )
               {
                  fmt::print( "~" );
                  offset -= (int) new_val * scale_factor[row];
               }
            }
            else
            {
               if( i != 0 )
                  fmt::print( " +" );
               int val = (int)( data.getValues()[i] * scale_factor[row] );
               fmt::print( "{} ", abs(val) );
               if(val < 0 )
               {
                  fmt::print( "~" );
                  offset -= val;
               }
            }
            fmt::print( "{}",  names[var_mapping[col]] );

         }
         fmt::print( " >= {} ;\n", (int)lhs *  scale_factor[row] + offset );
         fmt::print( "del id {}\n", lhs_row_mapping[row] );
         lhs_row_mapping[row] = next_constraint_id;
      }
      if( !rflags.test( RowFlag::kRhsInf ) )
      {
         next_constraint_id++;
         int offset = 0;
         fmt::print( "rup" );
         for( int i = 0; i < data.getLength(); i++ )
         {
            if( data.getIndices()[i] == col )
            {
               if( new_val == 0 )
                  continue;
               fmt::print( "{} ", abs((int) new_val * scale_factor[row]) );
               if( new_val < 0 )
                  offset += (int) new_val * scale_factor[row];
               else
                  fmt::print( "~" );
            }
            else
            {
               if( i != 0 )
                  fmt::print( " +" );
               int val = (int)( data.getValues()[i] * scale_factor[row] );
               fmt::print( "{} ", abs( val ) );
               if( val < 0 )
                  offset -= val;
               else
                  fmt::print( "~" );
            }
         }
         fmt::print( "{}",  names[var_mapping[col]] );
         fmt::print( " >= {} ;\n", (int)( rhs * scale_factor[row] ) + offset);
         fmt::print( "del id {}\n", rhs_row_mapping[row] );
         rhs_row_mapping[row] = next_constraint_id;
      }
   }


   void
   sparsify( int eqrow, int candrow, REAL scale, const Problem<REAL>& currentProblem )
   {
      const ConstraintMatrix<REAL>& matrix =
          currentProblem.getConstraintMatrix();
      //TODO: test
      //TODO: handle scale_factor
      assert( scale != 0 );
      if( ( num.isIntegral( scale ) || num.isIntegral( 1.0 / scale ) ) && scale > 0 )
         {
         int fac =  num.isIntegral( scale ) ? (int) scale : (int)(1.0 / scale) ;
         assert( fac > 0 );

         if( !matrix.getRowFlags()[candrow].test( RowFlag::kRhsInf ) )
         {
            next_constraint_id++;
            msg.info( "pol {} {} * {} +\n", (int)lhs_row_mapping[eqrow],
                      (int)( fac ), (int)rhs_row_mapping[candrow] );
            msg.info( "del id {}\n", (int)rhs_row_mapping[candrow] );
            rhs_row_mapping[candrow] = next_constraint_id;
         }
         if( !matrix.getRowFlags()[candrow].test( RowFlag::kLhsInf ) )
         {
            next_constraint_id++;
            msg.info( "pol {} {} * {} +\n", (int)rhs_row_mapping[eqrow] ,
                      (int)( fac ), (int)lhs_row_mapping[candrow] );
            msg.info( "del id {}\n", (int)lhs_row_mapping[candrow] );
            lhs_row_mapping[candrow] = next_constraint_id;
         }

      }
      else if( (num.isIntegral( 1.0 / scale ) || num.isIntegral( scale )) && scale < 0 )
      {
         int fac = num.isIntegral( scale ) ? (int) -scale : (int)(-1.0 / scale);
         assert( fac > 0 );

         if( !matrix.getRowFlags()[candrow].test( RowFlag::kRhsInf ) )
         {
            next_constraint_id++;
            msg.info( "pol {} {} * {} +\n", (int)rhs_row_mapping[candrow],
                      (int)( fac ), (int)lhs_row_mapping[eqrow] );
            msg.info( "del id {}\n", (int)rhs_row_mapping[candrow] );
            rhs_row_mapping[candrow] = next_constraint_id;
         }
         if( !matrix.getRowFlags()[candrow].test( RowFlag::kLhsInf ) )
         {
            next_constraint_id++;
            msg.info( "pol {} {} * {} +\n", (int)lhs_row_mapping[candrow] ,
                      (int)( fac ), (int)rhs_row_mapping[eqrow] );
            msg.info( "del id {}\n", (int)lhs_row_mapping[candrow] );
            lhs_row_mapping[candrow] = next_constraint_id;
         }
         scale_factor[candrow] *= fac;
      }
      else
      {
         // TODO:
      }
   }

   // TODO: test it
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
      REAL substitute_factor = indices[0] == col ? values[0]: values[1];

      next_constraint_id++;
      msg.info( "* postsolve stack : row id {}\n", next_constraint_id );
      msg.info( "rup {} {} + {} {} >= {};\n", (int)( values[0] ),
                var_mapping[indices[0]], (int)( values[1] ),
                var_mapping[indices[1]], (int)( offset ) );

      int lhs_id = next_constraint_id;
      next_constraint_id++;
      msg.info( "* postsolve stack : row id {}\n", next_constraint_id );
      msg.info( "rup {} ~{} + {} ~{} >= {};\n", (int)( values[0] ),
                var_mapping[indices[0]], (int)( values[1] ),
                var_mapping[indices[1]], (int)( offset ) );

      substitute( col, substitute_factor, lhs_id, next_constraint_id,
                  currentProblem );
   }

   void
   substitute( int col, int substituted_row, const Problem<REAL>& currentProblem )
   {
      const ConstraintMatrix<REAL>& matrix =
          currentProblem.getConstraintMatrix();
      auto col_vec = matrix.getColumnCoefficients( col );
      auto row_data = matrix.getRowCoefficients( substituted_row );

      REAL substitute_factor = 0;
      for( int i = 0; i < col_vec.getLength(); i++ )
      {
         if(col_vec.getIndices()[i] == substituted_row )
         {
            substitute_factor = col_vec.getValues()[i] * scale_factor[substituted_row];
            break;
         }
      }
      substitute( col, substitute_factor, lhs_row_mapping[substituted_row],
                  rhs_row_mapping[substituted_row], currentProblem,
                  substituted_row );
      assert( !matrix.getRowFlags()[substituted_row].test( RowFlag::kRhsInf ) );
      assert( !matrix.getRowFlags()[substituted_row].test( RowFlag::kLhsInf ) );
      msg.info( "* postsolve stack : row id {}\n", (int)rhs_row_mapping[substituted_row] );
      msg.info( "* postsolve stack : row id {}\n", (int)lhs_row_mapping[substituted_row] );
//      msg.info( "del id {}\n", (int) rhs_row_mapping[row] );
//      msg.info( "del id {}\n", (int) lhs_row_mapping[row] );
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
      compress_vector( rowmapping, scale_factor );
      if( full )
      {
         rhs_row_mapping.shrink_to_fit();
         lhs_row_mapping.shrink_to_fit();
         scale_factor.shrink_to_fit();
      }
#endif
   }

 private :

     void
     substitute( int col, REAL substitute_factor, int lhs_id, int rhs_id, const Problem<REAL>& currentProblem, int skip_row_id = -1 )
     {
      const ConstraintMatrix<REAL>& matrix =
          currentProblem.getConstraintMatrix();
      auto col_vec = matrix.getColumnCoefficients( col );

      for( int i = 0; i < col_vec.getLength(); i++ )
      {
         int row = col_vec.getIndices()[i];
         if( row == skip_row_id )
            continue ;
         REAL factor = col_vec.getValues()[i] * scale_factor[row];
         auto data = matrix.getRowCoefficients( row );
         if( num.isIntegral( factor / substitute_factor ) )
         {
            int val = (int)( factor / substitute_factor );
            if( !matrix.getRowFlags()[row].test( RowFlag::kRhsInf ) )
            {
               next_constraint_id++;
               if( substitute_factor * factor > 0)
                  msg.info( "pol {} {} * {} +\n", lhs_id, val,
                            rhs_row_mapping[row] );
               else
                  msg.info( "pol {} {} * {} +\n", lhs_id, val,
                            lhs_row_mapping[row] );
               msg.info( "del id {}\n", rhs_row_mapping[row] );
               rhs_row_mapping[row] = next_constraint_id;
            }
            if( !matrix.getRowFlags()[row].test( RowFlag::kLhsInf ) )
            {
               next_constraint_id++;
               if( substitute_factor * factor > 0 )
                  msg.info( "pol {} {} * {} +\n", rhs_id, val,
                            lhs_row_mapping[row] );
               else
                  msg.info( "pol {} {} * {} +\n", rhs_id, val,
                            rhs_row_mapping[row] );
               msg.info( "del id {}\n", lhs_row_mapping[row] );
               lhs_row_mapping[row] = next_constraint_id;
            }
         }
         else if( num.isIntegral( substitute_factor / factor ) )
         {
            scale_factor[row] *= (int) ( substitute_factor / factor);
            int val = abs( (int)( substitute_factor / factor ));
            if( !matrix.getRowFlags()[row].test( RowFlag::kRhsInf ) )
            {
               next_constraint_id++;

               if( substitute_factor * factor > 0 )
                  msg.info( "pol {} {} * {} +\n", rhs_row_mapping[row],
                            val,
                            lhs_id );
               else
                  msg.info( "pol {} {} * {} +\n", rhs_row_mapping[row],
                            val, rhs_id );

               msg.info( "del id {}\n", rhs_row_mapping[row] );
               rhs_row_mapping[row] = next_constraint_id;
            }
            if( !matrix.getRowFlags()[row].test( RowFlag::kLhsInf ) )
            {
               next_constraint_id++;
               if( substitute_factor * factor > 0 )
                  msg.info( "pol {} {} * {} +\n", lhs_row_mapping[row], val,
                            rhs_id );
               else
                  msg.info( "pol {} {} * {} +\n", lhs_row_mapping[row], val,
                            lhs_id );
               msg.info( "del id {}\n", lhs_row_mapping[row] );
               lhs_row_mapping[row] = next_constraint_id;
            }
         }
         else
         {
            assert( num.isIntegral( substitute_factor ) );
            assert( num.isIntegral( factor ) );
            scale_factor[row] *= (int)substitute_factor;
            int val = abs( (int) factor );
            int val2 = abs( (int)substitute_factor );

            if( !matrix.getRowFlags()[row].test( RowFlag::kRhsInf ) )
            {
               next_constraint_id++;
               if( substitute_factor * factor > 0 )
               {
                  msg.info( "pol {} {} * {} {} * +\n",
                            (int)lhs_id, val,
                            (int)rhs_row_mapping[row], val2 );
               }
               else
                  msg.info( "pol {} {} * {} {} * +\n",
                            (int)rhs_id, val,
                            (int)rhs_row_mapping[row], val2 );
               msg.info( "del id {}\n", (int)rhs_row_mapping[row] );
               rhs_row_mapping[row] = next_constraint_id;
            }
            if( !matrix.getRowFlags()[row].test( RowFlag::kLhsInf ) )
            {
               next_constraint_id++;
               if( substitute_factor * factor > 0 )
                  msg.info( "pol {} {} * {} {} * +\n",
                            (int)rhs_id, val,
                            (int)lhs_row_mapping[row], val2 );
               else
                  msg.info( "pol {} {} * {} {} * +\n",
                            (int)lhs_id, val,
                            (int)lhs_row_mapping[row], val2 );

               msg.info( "del id {}\n", (int)lhs_row_mapping[row] );
               lhs_row_mapping[row] = next_constraint_id;
            }
         }
      }

     };


};



} // namespace papilo

#endif
