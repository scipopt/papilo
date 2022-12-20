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

#ifndef _PAPILO_IO_MPS_WRITER_
#define _PAPILO_IO_MPS_WRITER_

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
struct MpsWriter
{
   static bool
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

#ifdef PAPILO_USE_BOOST_IOSTREAMS_WITH_ZLIB
      if( boost::algorithm::ends_with( filename, ".gz" ) )
         out.push( boost::iostreams::gzip_compressor() );
#endif

#ifdef PAPILO_USE_BOOST_IOSTREAMS_WITH_BZIP2
      if( boost::algorithm::ends_with( filename, ".bz2" ) )
         out.push( boost::iostreams::bzip2_compressor() );
#endif

      if( prob.getNumContinuousCols() > 0 )
         return false;
      for( int i = 0; i < prob.getNumIntegralCols(); ++i )
      {
         if( !prob.getVariableDomains().isBinary( i ) )
            return false;
      }

      out.push( file );

      fmt::print( out, "*ROWS:         {}\n", consmatrix.getNRows() );
      fmt::print( out, "*INTEGER:      {}\n", prob.getNumIntegralCols() );
      fmt::print( out, "*NONZERO:      {}\n*\n*\n", consmatrix.getNnz() );

      fmt::print( out, "min: ", consmatrix.getNnz() );

      bool obj_has_nonzeros = false;
      if( prob.getObjective().offset != 0 )
      {
         for( auto coef : prob.getObjective().coefficients )
            if( coef != 0 )
            {
               obj_has_nonzeros = true;
               break;
            }
      }
      else
         obj_has_nonzeros = true;

      if( objective_has_nonzeros( prob.getObjective() ) )
      {
         int offset = 0;
         for( int i = 0; i < prob.getNCols(); i++ )
         {
            REAL coef = prob.getObjective().coefficients[i];
            if( coef == 0 )
               continue;
            if( coef < 0 )
               offset += abs( (int)coef );
            fmt::print( out, "+{} {} ", coef > 0 ? "" : "~", abs( (int)coef ),
                        prob.getConstraintNames()[col_mapping[i]] );
         }
         assert( prob.getObjective().offset >= 0 );
         int obj_offset = (int) prob.getObjective().offset + offset;
         if( obj_offset != 0 )
            fmt::print( out, "+{} ", abs( obj_offset ) );
         fmt::print( "\n;" );
      };

      //TODO: scaling and more
      for( int i = 0; i < consmatrix.getNRows(); ++i )
      {
         assert( !consmatrix.isRowRedundant( i ) );
         char type;
         if( row_flags[i].test( RowFlag::kEquation ) &&
             row_flags[i].test( RowFlag::kRhsInf ) )
         {
            assert( row_flags[i].test( RowFlag::kLhsInf ) );
         }

         if( row_flags[i].test( RowFlag::kLhsInf ) &&
             row_flags[i].test( RowFlag::kRhsInf ) )
            type = 'N';
         else if( row_flags[i].test( RowFlag::kRhsInf ) )
            type = 'G';
         else if( row_flags[i].test( RowFlag::kLhsInf ) )
            type = 'L';
         else
         {
            type = 'E';
         }

         fmt::print( out, " {}  {}\n", type, consnames[row_mapping[i]] );
      }
      return true;
   }
};

} // namespace papilo

#endif
