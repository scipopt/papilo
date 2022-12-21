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

#ifndef _PAPILO_IO_OPB_WRITER_
#define _PAPILO_IO_OPB_WRITER_

#include "papilo/Config.hpp"
#include "papilo/core/Problem.hpp"
#include "papilo/misc/Vec.hpp"
#include "papilo/misc/fmt.hpp"
#include <boost/algorithm/string/predicate.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <cmath>
#include <fstream>
#include <iostream>

#ifdef PAPILO_USE_BOOST_IOSTREAMS_WITH_BZIP2
#include <boost/iostreams/filter/bzip2.hpp>
#endif
#ifdef PAPILO_USE_BOOST_IOSTREAMS_WITH_ZLIB
#include <boost/iostreams/filter/gzip.hpp>
#endif

namespace papilo
{

/// Writer to write problem structures into an opb file
template <typename REAL>
struct OpbWriter
{
   static bool
   writeProb( const std::string& filename, const Problem<REAL>& prob,
              const Vec<int>& col_mapping )
   {
      const ConstraintMatrix<REAL>& matrix = prob.getConstraintMatrix();
      const Vec<std::string>& consnames = prob.getConstraintNames();
      const Vec<std::string>& varnames = prob.getVariableNames();
      const Vec<REAL>& lhs = matrix.getLeftHandSides();
      const Vec<REAL>& rhs = matrix.getRightHandSides();
      const Objective<REAL>& obj = prob.getObjective();
      const Vec<ColFlags>& col_flags = prob.getColFlags();
      const Vec<RowFlags>& row_flags = prob.getRowFlags();

      std::ofstream file( filename, std::ofstream::out );
      boost::iostreams::filtering_ostream out;

#ifdef PAPILO_USE_BOOST_IOSTREAMS_WITH_ZLIB
      if( boost::algorithm::ends_with( filename, ".gz" ) )
         out.push( boost::iostreams::gzip_compressor() );
#endif

#ifdef PAPILO_USE_BOOST_IOSTREAMS_WITH_BZIP2
      if( boost::algorithm::ends_with( filename, ".bz2" ) )
         out.push( boost::iostreams::bzip2_compressor() );
#endif

      if( prob.getNumContinuousCols() > 0 )
      {
         fmt::print("Problem contains continuous variables. Opb is not the write format.");
         return false;
      }
      for( int i = 0; i < prob.getNumIntegralCols(); ++i )
      {
         if( !prob.getVariableDomains().isBinary( i ) )
         {
            fmt::print("Problem contains non binary integer variables. Opb is not the write format.");
            return false;
         }
      }
      for( int i = 0; i < prob.getNRows(); ++i )
      {
         auto vector = matrix.getRowCoefficients(i);
         for( int j = 0; j < vector.getLength(); j++ )
            //TODO: check values
            //if( !num.isIntegral( vector.getValues()[j] ) )
            {
            //   fmt::print( "Matrix contains fractional values. Opb is not the write format." );
            //   return false;
            }
      }
      out.push( file );

      fmt::print( out, "*NAME          {}\n", prob.getName() );
      fmt::print( out, "*ROWS:         {}\n", matrix.getNRows() );
      fmt::print( out, "*VARIABLES:    {}\n", prob.getNumIntegralCols() );
      fmt::print( out, "*NONZERO:      {}\n*\n*\n", matrix.getNnz() );

      bool obj_has_nonzeros = false;
      if( obj.offset != 0 )
      {
         for( auto coef : obj.coefficients )
            if( coef != 0 )
            {
               obj_has_nonzeros = true;
               break;
            }
      }
      else
         obj_has_nonzeros = true;

      if( obj_has_nonzeros )
      {
         fmt::print(out, "min: ");
         for( int i = 0; i < prob.getNCols(); i++ )
         {
            REAL coef = obj.coefficients[i];
            if( coef == 0 )
               continue;
            fmt::print( out, "{}{} {} ", coef > 0 ? "+" : "-", abs( (int)coef ),
                        varnames[col_mapping[i]] );
         }
         int obj_offset = (int) prob.getObjective().offset;
         if( obj_offset != 0 )
            fmt::print( out, "{}{} ", obj_offset > 0 ? "+" : "-", abs( obj_offset ) );
         fmt::print( ";\n" );
      }

      //TODO: scaling and more
      for( int row = 0; row < matrix.getNRows(); ++row )
      {
         assert( !matrix.isRowRedundant( row ) );
         char type;
         if( row_flags[row].test( RowFlag::kEquation ) &&
             row_flags[row].test( RowFlag::kRhsInf ) )
         {
            assert( !row_flags[row].test( RowFlag::kLhsInf ) );
            auto vector = matrix.getRowCoefficients( row );
            for( int j = 0; j < vector.getLength(); j++ )
            {
               REAL val = vector.getValues()[j];
               assert( !val == 0 );
               fmt::print( out, "{}{} {} ", val > 0 ? "+" : "-",
                           abs( (int)val ),
                           varnames[col_mapping[vector.getIndices()[j]]] );
            }
            if(row_flags[row].test( RowFlag::kEquation ))
               fmt::print(out, "= ");
            else
               fmt::print(out, ">= ");
            fmt::print(out, " {} ;\n", lhs[row]);
         }
         else{
            assert(row_flags[row].test( RowFlag::kLhsInf ));
            assert(!row_flags[row].test( RowFlag::kRhsInf ));
            auto vector = matrix.getRowCoefficients( row );
            for( int j = 0; j < vector.getLength(); j++ )
            {
               REAL val = vector.getValues()[j];
               assert( !val == 0 );
               fmt::print( out, "{}{} {} ", val < 0 ? "+" : "-",
                           abs( (int)val ),
                           varnames[col_mapping[vector.getIndices()[j]]] );
            }
            if(row_flags[row].test( RowFlag::kEquation ))
               fmt::print(out, "= ");
            else
               fmt::print(out, ">= ");
            fmt::print(out, " {} ;\n", rhs[row]);
         }

      }
      return true;
   }
};

} // namespace papilo

#endif
