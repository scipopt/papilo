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

#ifndef _IO_MPS_WRITER_
#define _IO_MPS_WRITER_

#include "papilo/core/Problem.hpp"
#include "papilo/misc/Vec.hpp"
#include "papilo/misc/fmt.hpp"
#include <boost/algorithm/string/predicate.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <cmath>
#include <fstream>
#include <iostream>

namespace papilo
{

/// Writer to write problem structures into an mps file
template <typename REAL>
struct MpsWriter
{
   static void
   writeProb( const std::string& filename, const Problem<REAL>& prob,
              const Vec<int>& row_mapping, const Vec<int>& col_mapping )
   {
      const ConstraintMatrix<REAL>& consmatrix = prob.getConstraintMatrix();
      const Vec<std::string>& consnames = prob.getConstraintNames();
      const Vec<std::string>& varnames = prob.getVariableNames();
      const Vec<REAL>& lhs = consmatrix.getLeftHandSides();
      const Vec<REAL>& rhs = consmatrix.getRightHandSides();
      const Objective<REAL>& obj = prob.getObjective();
      const Vec<ColFlags>& col_flags = prob.getColFlags();
      const Vec<RowFlags>& row_flags = prob.getRowFlags();

      std::ofstream file( filename, std::ofstream::out );
      boost::iostreams::filtering_ostream out;

      if( boost::algorithm::ends_with( filename, ".gz" ) )
         out.push( boost::iostreams::gzip_compressor() );
      else if( boost::algorithm::ends_with( filename, ".bz2" ) )
         out.push( boost::iostreams::bzip2_compressor() );

      out.push( file );

      fmt::print( out, "*ROWS:         {}\n", consmatrix.getNRows() );
      fmt::print( out, "*COLUMNS:      {}\n", consmatrix.getNCols() );
      fmt::print( out, "*INTEGER:      {}\n", prob.getNumIntegralCols() );
      fmt::print( out, "*NONZERO:      {}\n*\n*\n", consmatrix.getNnz() );

      fmt::print( out, "NAME          {}\n", prob.getName() );
      fmt::print( out, "ROWS\n" );
      fmt::print( out, " N  OBJ\n" );
      bool hasRangedRow = false;
      for( int i = 0; i < consmatrix.getNRows(); ++i )
      {
         // discard redundant rows when writing problem
         assert( !consmatrix.isRowRedundant( i ) );

         char type;

         if( row_flags[i].test( RowFlag::LHS_INF ) &&
             row_flags[i].test( RowFlag::RHS_INF ) )
            type = 'N';
         else if( row_flags[i].test( RowFlag::RHS_INF ) )
            type = 'G';
         else if( row_flags[i].test( RowFlag::LHS_INF ) )
            type = 'L';
         else
         {
            if( !row_flags[i].test( RowFlag::EQUALITY ) )
               hasRangedRow = true;
            type = 'E';
         }

         fmt::print( out, " {}  {}\n", type, consnames[row_mapping[i]] );
      }

      fmt::print( out, "COLUMNS\n" );

      int hasintegral = prob.getNumIntegralCols() != 0;

      for( int integral = 0; integral <= hasintegral; ++integral )
      {
         if( integral )
            fmt::print( out,
                        "    MARK0000  'MARKER'                 'INTORG'\n" );

         for( int i = 0; i < consmatrix.getNCols(); ++i )
         {
            if( col_flags[i].test( ColFlag::INACTIVE ) )
               continue;

            if( ( !col_flags[i].test( ColFlag::INTEGRAL ) && integral ) ||
                ( col_flags[i].test( ColFlag::INTEGRAL ) && !integral ) )
               continue;

            assert(
                !col_flags[i].test( ColFlag::FIXED, ColFlag::SUBSTITUTED ) );

            if( obj.coefficients[i] != 0.0 )
            {
               fmt::print( out, "    {: <9} OBJ       {:.15}\n",
                           varnames[col_mapping[i]],
                           double( obj.coefficients[i] ) );
            }

            SparseVectorView<REAL> column =
                consmatrix.getColumnCoefficients( i );

            const int* rowinds = column.getIndices();
            const REAL* colvals = column.getValues();
            int len = column.getLength();

            for( int j = 0; j < len; ++j )
            {
               int r = rowinds[j];

               // discard redundant rows when writing problem
               if( consmatrix.isRowRedundant( r ) )
                  continue;

               // normal row
               fmt::print( out, "    {: <9} {: <9} {:.15}\n",
                           varnames[col_mapping[i]], consnames[row_mapping[r]],
                           double( colvals[j] ) );
            }
         }

         if( integral )
            fmt::print( out,
                        "    MARK0000  'MARKER'                 'INTEND'\n" );
      }

      const Vec<REAL>& lower_bounds = prob.getLowerBounds();
      const Vec<REAL>& upper_bounds = prob.getUpperBounds();

      fmt::print( out, "RHS\n" );

      if( obj.offset != 0 )
      {
         if( obj.offset != REAL{0.0} )
            fmt::print( out, "    B         {: <9} {:.15}\n", "OBJ",
                        double( -obj.offset ) );
      }

      for( int i = 0; i < consmatrix.getNRows(); ++i )
      {
         // discard redundant rows when writing problem
         if( consmatrix.isRowRedundant( i ) )
            continue;

         if( row_flags[i].test( RowFlag::LHS_INF ) &&
             row_flags[i].test( RowFlag::RHS_INF ) )
            continue;

         if( row_flags[i].test( RowFlag::LHS_INF ) )
         {
            if( rhs[i] != REAL{0.0} )
               fmt::print( out, "    B         {: <9} {:.15}\n",
                           consnames[row_mapping[i]], double( rhs[i] ) );
         }
         else
         {
            if( lhs[i] != REAL{0.0} )
               fmt::print( out, "    B         {: <9} {:.15}\n",
                           consnames[row_mapping[i]], double( lhs[i] ) );
         }
      }

      if( hasRangedRow )
      {
         fmt::print( out, "RANGES\n" );
         for( int i = 0; i < consmatrix.getNRows(); ++i )
         {
            if( row_flags[i].test( RowFlag::LHS_INF, RowFlag::RHS_INF,
                                   RowFlag::EQUALITY, RowFlag::REDUNDANT ) )
               continue;

            double rangeval = double( rhs[i] - lhs[i] );

            if( rangeval != 0 )
            {
               fmt::print( out, "    B         {: <9} {:.15}\n",
                           consnames[row_mapping[i]], rangeval );
            }
         }
      }

      fmt::print( out, "BOUNDS\n" );

      for( int i = 0; i < consmatrix.getNCols(); ++i )
      {
         if( col_flags[i].test( ColFlag::INACTIVE ) )
            continue;

         if( !col_flags[i].test( ColFlag::LB_INF ) &&
             !col_flags[i].test( ColFlag::UB_INF ) &&
             lower_bounds[i] == upper_bounds[i] )
         {
            fmt::print( out, " FX BND       {: <9} {:.15}\n",
                        varnames[col_mapping[i]], double( lower_bounds[i] ) );
         }
         else
         {
            if( col_flags[i].test( ColFlag::LB_INF ) || lower_bounds[i] != 0.0 )
            {
               if( col_flags[i].test( ColFlag::LB_INF ) )
                  fmt::print( out, " MI BND       {}\n",
                              varnames[col_mapping[i]] );
               else
                  fmt::print( out, " LO BND       {: <9} {:.15}\n",
                              varnames[col_mapping[i]],
                              double( lower_bounds[i] ) );
            }

            if( !col_flags[i].test( ColFlag::UB_INF ) )
               fmt::print( out, " UP BND       {: <9} {:.15}\n",
                           varnames[col_mapping[i]],
                           double( upper_bounds[i] ) );
            else
               fmt::print( out, " PL BND       {: <9}\n",
                           varnames[col_mapping[i]] );
         }
      }
      fmt::print( out, "ENDATA\n" );
   }
};

} // namespace papilo

#endif
