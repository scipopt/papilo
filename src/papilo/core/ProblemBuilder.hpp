/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*               This file is part of the program and library                */
/*    PaPILO --- Parallel Presolve for Integer and Linear Optimization       */
/*                                                                           */
/* Copyright (C) 2020  Konrad-Zuse-Zentrum                                   */
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

#ifndef _PAPILO_CORE_PROBLEM_BUILDER_HPP_
#define _PAPILO_CORE_PROBLEM_BUILDER_HPP_

namespace papilo
{

#include "papilo/core/MatrixBuffer.hpp"
#include "papilo/core/Problem.hpp"
#include "papilo/misc/String.hpp"
#include "papilo/misc/Vec.hpp"

template <typename REAL>
class ProblemBuilder
{
 public:
   /// Sets the number of columns to the given value. The information of columns
   /// that already exist is kept, new columns are continuous and fixed to zero.
   /// If the number of columns is reduced, then the information for the removed
   /// columns is lost.
   void
   setNumCols( int ncols )
   {
      // allocate column information
      obj.coefficients.resize( ncols );
      domains.lower_bounds.resize( ncols );
      domains.upper_bounds.resize( ncols );
      domains.flags.resize( ncols );
      colnames.resize( ncols );
   }

   /// Sets the number of rows to the given value. The information of rows that
   /// already exist is kept, new rows are qualities with right hand side zero.
   /// If the number of rows is reduced, then the information for the removed
   /// rows is lost.
   void
   setNumRows( int nrows )
   {
      // allocate row information
      lhs.resize( nrows );
      rhs.resize( nrows );
      rflags.resize( nrows );
      rownames.resize( nrows );
   }

   /// Returns the current number of rows
   int
   getNumRows() const
   {
      return static_cast<int>( rflags.size() );
   }

   /// Returns the current number of cols
   int
   getNumCols() const
   {
      return static_cast<int>( domains.flags.size() );
   }

   /// reserve storage for the given number of non-zeros
   void
   reserve( int nnz, int nrows, int ncols )
   {
      matrix_buffer.reserve( nnz );

      // reserve space for row information
      lhs.reserve( nrows );
      rhs.reserve( nrows );
      rflags.reserve( nrows );
      rownames.reserve( nrows );

      // reserve space for column information
      obj.coefficients.reserve( ncols );
      domains.lower_bounds.reserve( ncols );
      domains.upper_bounds.reserve( ncols );
      domains.flags.reserve( ncols );
      colnames.reserve( ncols );
   }

   /// change the objective coefficient of a column
   void
   setObj( int col, REAL val )
   {
      obj.coefficients[col] = std::move( val );
   }

   /// change the objectives constant offset
   void
   setObjOffset( REAL val )
   {
      obj.offset = std::move( val );
   }

   void
   setColLbInf( int col, bool isInfinite )
   {
      if( isInfinite )
         domains.flags[col].set( ColFlag::kLbInf );
      else
         domains.flags[col].unset( ColFlag::kLbInf );
   }

   void
   setColUbInf( int col, bool isInfinite )
   {
      if( isInfinite )
         domains.flags[col].set( ColFlag::kUbInf );
      else
         domains.flags[col].unset( ColFlag::kUbInf );
   }

   void
   setColLb( int col, REAL lb )
   {
      domains.lower_bounds[col] = std::move( lb );
   }

   void
   setColUb( int col, REAL ub )
   {
      domains.upper_bounds[col] = std::move( ub );
   }

   void
   setColIntegral( int col, bool isIntegral )
   {
      if( isIntegral )
         domains.flags[col].set( ColFlag::kIntegral );
      else
         domains.flags[col].unset( ColFlag::kIntegral );
   }

   void
   setRowLhsInf( int row, bool isInfinite )
   {
      if( isInfinite )
         rflags[row].set( RowFlag::kLhsInf );
      else
         rflags[row].unset( RowFlag::kLhsInf );
   }

   void
   setRowRhsInf( int row, bool isInfinite )
   {
      if( isInfinite )
         rflags[row].set( RowFlag::kRhsInf );
      else
         rflags[row].unset( RowFlag::kRhsInf );
   }

   void
   setRowLhs( int row, REAL lhsval )
   {
      lhs[row] = std::move( lhsval );
   }

   void
   setRowRhs( int row, REAL rhsval )
   {
      rhs[row] = std::move( rhsval );
   }

   /// add a nonzero entry for the given row and column
   void
   addEntry( int row, int col, const REAL& val )
   {
      assert( val != 0 );
      matrix_buffer.addEntry( row, col, val );
   }

   /// add the nonzero entries for the given row
   template <typename R>
   void
   addRowEntries( int row, int len, const int* cols, const R* vals )
   {
      for( int i = 0; i != len; ++i )
      {
         assert( vals[i] != 0 );
         matrix_buffer.addEntry( row, cols[i], REAL{ vals[i] } );
      }
   }

   template <typename Str>
   void
   setRowName( int row, Str&& name )
   {
      rownames[row] = String( name );
   }

   template <typename Str>
   void
   setColName( int col, Str&& name )
   {
      colnames[col] = String( name );
   }

   template <typename Str>
   void
   setProblemName( Str&& name )
   {
      probname = String( name );
   }

   /// add the nonzero entries for the given column
   template <typename R>
   void
   addColEntries( int col, int len, const int* rows, const R* vals )
   {
      for( int i = 0; i != len; ++i )
      {
         assert( vals[i] != 0 );
         matrix_buffer.addEntry( rows[i], col, REAL{ vals[i] } );
      }
   }

   Problem<REAL>
   build()
   {
      Problem<REAL> problem;

      int nrows = lhs.size();
      int ncols = obj.coefficients.size();

      problem.setName( std::move( probname ) );

      problem.setConstraintMatrix( ConstraintMatrix<REAL>{
          matrix_buffer.buildCSR( nrows, ncols ),
          matrix_buffer.buildCSC( nrows, ncols ), std::move( lhs ),
          std::move( rhs ), std::move( rflags ) } );

      matrix_buffer.clear();

      problem.setObjective( std::move( obj ) );
      problem.setVariableDomains( std::move( domains ) );
      problem.setVariableNames( std::move( colnames ) );
      problem.setConstraintNames( std::move( rownames ) );

      return problem;
   }

 private:
   MatrixBuffer<REAL> matrix_buffer;
   Objective<REAL> obj;
   VariableDomains<REAL> domains;
   Vec<REAL> lhs;
   Vec<REAL> rhs;
   Vec<RowFlags> rflags;
   Vec<String> rownames;
   Vec<String> colnames;
   String probname;
};

} // namespace papilo

#endif