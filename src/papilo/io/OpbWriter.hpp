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
              const Vec<int>& col_mapping, const Vec<int>& row_scaling,
              const Num<REAL>& num)
   {
      const ConstraintMatrix<REAL>& matrix = prob.getConstraintMatrix();
      const Vec<std::string>& varnames = prob.getVariableNames();
      const Vec<REAL>& lhs = matrix.getLeftHandSides();
      const Vec<REAL>& rhs = matrix.getRightHandSides();
      const Objective<REAL>& obj = prob.getObjective();
      const Vec<RowFlags>& row_flags = prob.getRowFlags();
      const Vec<Symmetry>& sym = prob.getSymmetries().symmetries;

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
         fmt::print( "Problem contains continuous variables. Opb is not the "
                     "write format." );
         return false;
      }
      for( int i = 0; i < prob.getNumIntegralCols(); ++i )
      {
         if( !prob.getVariableDomains().isBinary( i ) )
         {
            fmt::print( "Problem contains non binary integer variables. Opb is "
                        "not the write format." );
            return false;
         }
      }
      for( int i = 0; i < prob.getNRows(); ++i )
      {
         auto vector = matrix.getRowCoefficients( i );
         for( int j = 0; j < vector.getLength(); j++ )
            if( !num.isIntegral( vector.getValues()[j] * row_scaling[i] ) )
            {
               fmt::print( "Matrix contains fractional values. Opb is not the "
                           "write format." );
               return false;
            }
         if( !num.isIntegral( lhs[i] * row_scaling[i] ) ||
             !num.isIntegral( rhs[i] * row_scaling[i] ) )
         {
            fmt::print( "Lhs/Rhs contains fractional values. Opb is not the "
                        "correct format.\n" );
            return false;
         }
      }
      out.push( file );

      fmt::print( out, "* #variable= {} #constraint= {}\n",
                  getVars( prob, col_mapping ), getRows( prob ) );
      fmt::print( out, "* Objective Offset {}\n", boost::multiprecision::cpp_int( abs(prob.getObjective().offset) ).str() );

      bool obj_has_nonzeros = false;
      if( obj.offset == 0 )
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
         fmt::print( out, "min: " );
         for( int i = 0; i < prob.getNCols(); i++ )
         {
            REAL coef = obj.coefficients[i];
            if( coef == 0 )
               continue;
            fmt::print( out, "{}{} {} ", coef > 0 ? "+" : "-", boost::multiprecision::cpp_int( abs(coef) ).str(),
                        varnames[col_mapping[i]] );
         }
         fmt::print( out, ";\n" );
      }

      for( int row = 0; row < matrix.getNRows(); ++row )
      {
         assert( !matrix.isRowRedundant( row ) );
         auto vector = matrix.getRowCoefficients( row );

         bool scale_necessary = false;
         for( int j = 0; j < vector.getLength(); ++j )
         {
            if(scale_necessary)
               break;
            scale_necessary = !num.isIntegral(vector.getValues()[j]);
         }
         if( row_flags[row].test( RowFlag::kEquation ) ||
             !row_flags[row].test( RowFlag::kLhsInf ) )
         {
            scale_necessary = scale_necessary || !num.isIntegral(lhs[row]);
            REAL scale = scale_necessary ? abs(row_scaling[row]) : 1;

            assert( !row_flags[row].test( RowFlag::kLhsInf ) );
            for( int j = 0; j < vector.getLength(); j++ )
            {
               REAL val = vector.getValues()[j] * scale;
               assert( val != 0 );
               assert( num.isIntegral(val ) );
               fmt::print( out, "{}{} {} ", val > 0 ? "+" : "-",
                           abs( cast_to_long( val ) ),
                           varnames[col_mapping[vector.getIndices()[j]]] );
            }
            assert(num.isIntegral( lhs[row] * scale ));
            if( row_flags[row].test( RowFlag::kEquation ) )
               fmt::print( out, "= " );
            else
               fmt::print( out, ">= " );
            fmt::print( out, " {} ;\n", cast_to_long( lhs[row] * scale ) );
         }
         if(!row_flags[row].test( RowFlag::kEquation ) &&
             !row_flags[row].test( RowFlag::kRhsInf ))
         {
            assert( !row_flags[row].test( RowFlag::kRhsInf ) );

            scale_necessary = scale_necessary || !num.isIntegral(rhs[row]);
            REAL scale = scale_necessary ? abs(row_scaling[row]) : 1;
            for( int j = 0; j < vector.getLength(); j++ )
            {
               REAL val = vector.getValues()[j] * scale;
               assert( val != 0 );
               assert( num.isIntegral(val ) );
               fmt::print( out, "{}{} {} ", val < 0 ? "+" : "-",
                           boost::multiprecision::cpp_int( abs( val ) ).str(),
                           varnames[col_mapping[vector.getIndices()[j]]] );
            }
            assert(num.isIntegral( rhs[row] * scale ));
            fmt::print( out, ">= {} ;\n",
                        -cast_to_long( rhs[row] * scale ) );
         }
      }
      for( auto sym_c : sym )
      {
         switch( sym_c.getSymmetryType() )
         {
         case SymmetryType::kXgeY:
            fmt::print( out, "+1 {} +1 ~{} >= 1\n",
                        varnames[col_mapping[sym_c.getDominatingCol()]],
                        varnames[col_mapping[sym_c.getDominatedCol()]] );
            break;
         case SymmetryType::kXplusYge1:
            fmt::print( out, "+1 {} +1 {} >= 1\n",
                        varnames[col_mapping[sym_c.getDominatingCol()]],
                        varnames[col_mapping[sym_c.getDominatedCol()]] );
            break;
         default:
            assert( false );
         }

      }

      return true;
   }

   static int
   getVars( const Problem<REAL>& prob,
                 const Vec<int>& col_mapping)
   {
      int var = -1;
      for( int i = 0; i < prob.getNCols(); i++ )
      {
         auto var_name = prob.getVariableNames()[col_mapping[i]];
         assert( !var_name.empty() );
         assert( var_name[0] == 'x' );
         var_name = var_name.substr( 1 );
         int v = atoi( var_name.c_str() );
         if( v >= var )
            var = v;
      }
      assert(var > 0);
      return var;
   }

   static int
   getRows( const Problem<REAL>& prob)
   {
      int rows = 0;
      for(int i =0; i< prob.getNRows(); i++)
      {
         if(prob.getRowFlags()[i].test(RowFlag::kEquation))
         {
            rows++;
            continue ;
         }
         if(!prob.getRowFlags()[i].test(RowFlag::kLhsInf))
            rows++;
         if(!prob.getRowFlags()[i].test(RowFlag::kRhsInf))
            rows++;
      }
      return rows;

   }

   static long
   cast_to_long( const REAL& x )
   {
      return (long) ( REAL( x + REAL( 0.5 ) ) );
   }
};

} // namespace papilo

#endif
