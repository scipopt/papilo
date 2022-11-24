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

#ifndef _PAPILO_IO_OPB_PARSER_HPP_
#define _PAPILO_IO_OPB_PARSER_HPP_

#include "papilo/Config.hpp"
#include "papilo/core/ConstraintMatrix.hpp"
#include "papilo/core/Objective.hpp"
#include "papilo/core/Problem.hpp"
#include "papilo/core/VariableDomains.hpp"
#include "papilo/external/pdqsort/pdqsort.h"
#include "papilo/io/BoundType.hpp"
#include "papilo/io/ParseKey.hpp"
#include "papilo/misc/Flags.hpp"
#include "papilo/misc/Hash.hpp"
#include "papilo/misc/Num.hpp"
#include <algorithm>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/optional.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/utility/string_ref.hpp>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <memory>
#include <regex>
#include <tuple>
#include <utility>

#ifdef PAPILO_USE_BOOST_IOSTREAMS_WITH_BZIP2
#include <boost/iostreams/filter/bzip2.hpp>
#endif
#ifdef PAPILO_USE_BOOST_IOSTREAMS_WITH_ZLIB
#include <boost/iostreams/filter/gzip.hpp>
#endif

/* This file reader parses the @a opb format and is also used by the @a wbo
 * reader for the @a wbo format. For a detailed description of this format see
 *
 * - http://www.cril.univ-artois.fr/PB07/solver_req.html
 * - http://www.cril.univ-artois.fr/PB10/format.pdf
 */

namespace papilo
{

/// Parser for mps files in fixed and free format
template <typename REAL>
class OpbParser
{
   static_assert(
       num_traits<typename RealParseType<REAL>::type>::is_floating_point,
       "the parse type must be a floating point type" );

 public:
   static boost::optional<Problem<REAL>>
   loadProblem( const std::string& filename )
   {
      OpbParser<REAL> parser;

      Problem<REAL> problem;

      if( !parser.parseFile( filename ) )
         return boost::none;

      assert( parser.nnz >= 0 );

      Vec<REAL> obj_vec( size_t( parser.nCols ), REAL{ 0.0 } );

      for( auto i : parser.coeffobj )
         obj_vec[i.first] = i.second;

      problem.setObjective( std::move( obj_vec ), parser.objoffset );
      problem.setConstraintMatrix(
          SparseStorage<REAL>{ std::move( parser.entries ), parser.nRows,
                               parser.nCols },
          std::move( parser.rowlhs ), std::move( parser.rowrhs ),
          std::move( parser.row_flags ) );
      problem.setVariableDomains( std::move( parser.lb4cols ),
                                  std::move( parser.ub4cols ),
                                  std::move( parser.col_flags ) );
      problem.setVariableNames( std::move( parser.colnames ) );
      problem.setName( std::move( filename ) );
      problem.setConstraintNames( std::move( parser.rownames ) );

      problem.setInputTolerance(
          REAL{ pow( typename RealParseType<REAL>::type{ 10 },
                     -std::numeric_limits<
                         typename RealParseType<REAL>::type>::digits10 ) } );
      return problem;
   }

 private:
   OpbParser() = default;

   /// load LP from MPS file as transposed triplet matrix
   bool
   parseFile( const std::string& filename );

   bool
   parse( boost::iostreams::filtering_istream& file );

   void
   printErrorMessage( ParseKey keyword )
   {
      switch( keyword )
      {
      case ParseKey::kRows:
         std::cerr << "read error in section ROWS " << std::endl;
         break;
      case ParseKey::kCols:
         std::cerr << "read error in section COLUMNS " << std::endl;
         break;
      case ParseKey::kRhs:
         std::cerr << "read error in section RHS " << std::endl;
         break;
      case ParseKey::kBounds:
         std::cerr << "read error in section BOUNDS " << std::endl;
         break;
      case ParseKey::kRanges:
         std::cerr << "read error in section RANGES " << std::endl;
         break;
      default:
         std::cerr << "undefined read error " << std::endl;
         break;
      }
   };

   /*
    * data for mps problem
    */

   Vec<Triplet<REAL>> entries;
   Vec<std::pair<int, REAL>> coeffobj;
   Vec<REAL> rowlhs;
   Vec<REAL> rowrhs;
   Vec<std::string> rownames;
   Vec<std::string> colnames;

   HashMap<std::string, int> rowname2idx;
   HashMap<std::string, int> colname2idx;
   Vec<REAL> lb4cols;
   Vec<REAL> ub4cols;
   Vec<BoundType> row_type;
   Vec<RowFlags> row_flags;
   Vec<ColFlags> col_flags;
   REAL objoffset = 0;

   int nCols = 0;
   int nRows = 0;
   int nnz = -1;

   ParseKey
   parseRows( const std::string& line );

   ParseKey
   parseObjective( const std::string& line );

   Vec<String>
   split( const char* str );

   bool
   isNumeric( const String& token );

   void
   add_binary_variable( const String& name );

   REAL
   to_int( const String& token ) const;
};

template <typename REAL>
bool
OpbParser<REAL>::parseFile( const std::string& filename )
{
   std::ifstream file( filename, std::ifstream::in );
   boost::iostreams::filtering_istream in;

   if( !file )
      return false;

#ifdef PAPILO_USE_BOOST_IOSTREAMS_WITH_ZLIB
   if( boost::algorithm::ends_with( filename, ".gz" ) )
      in.push( boost::iostreams::gzip_decompressor() );
#endif

#ifdef PAPILO_USE_BOOST_IOSTREAMS_WITH_BZIP2
   if( boost::algorithm::ends_with( filename, ".bz2" ) )
      in.push( boost::iostreams::bzip2_decompressor() );
#endif

   in.push( file );

   return parse( in );
}

template <typename REAL>
bool
OpbParser<REAL>::parse( boost::iostreams::filtering_istream& file )
{
   nnz = 0;
   std::string strline;

   while( getline( file, strline ) )
   {
      if( strline.rfind( '*', 0 ) == 0 || strline.empty() )
         continue;
      else if( strline.rfind( "min:", 0 ) == 0 )
      {
         auto type = parseObjective( strline );
         if( type == ParseKey::kFail )
            return false;
      }
      else
      {
         auto type = parseRows( strline );
         if( type == ParseKey::kFail )
            return false;
      }
   }

   assert( row_type.size() == unsigned( nRows ) );
   assert( nCols == colname2idx.size() );
   assert( nRows == rowname2idx.size() );

   return true;
}

template <typename REAL>
ParseKey
OpbParser<REAL>::parseRows( const std::string& line )
{
   auto tokens = split( line.c_str() );
   bool init_last_value = false;
   bool operator_occurred = false;
   REAL last_value;
   REAL offset = REAL{ 0 };
   rownames.push_back( std::to_string( nRows ) );
   rowname2idx.insert( { std::to_string( nRows ), nRows } );
   for( String token : tokens )
   {
      if( token == "<=" )
      {
         row_type.push_back( BoundType::kLE );
         RowFlags flags{};
         flags.unset( RowFlag::kRhsInf );
         flags.set( RowFlag::kLhsInf );
         row_flags.push_back( flags );
         operator_occurred = true;
      }
      else if( token == ">=" )
      {
         row_type.push_back( BoundType::kGE );
         RowFlags flags{};
         flags.set( RowFlag::kRhsInf );
         flags.unset( RowFlag::kLhsInf );
         row_flags.push_back( flags );
         operator_occurred = true;
      }
      else if( token == "=" )
      {
         row_type.push_back( BoundType::kEq );
         RowFlags flags{};
         flags.unset( RowFlag::kRhsInf );
         flags.unset( RowFlag::kLhsInf );
         flags.set( RowFlag::kEquation );
         row_flags.push_back( flags );
         operator_occurred = true;
      }
      else if( isNumeric( token ) )
      {
         init_last_value = true;
         if( operator_occurred )
         {
            if( row_type[row_type.size() - 1] == BoundType::kEq )
            {
               REAL val = to_int( token ) - offset;
               rowrhs.push_back( val );
               rowlhs.push_back( val );
            }
            else if( row_type[row_type.size() - 1] == BoundType::kGE )
            {
               rowlhs.push_back( to_int( token ) - offset );
               rowrhs.push_back( REAL{ 0 } );
            }
            else
            {
               assert( row_type[row_type.size() - 1] == BoundType::kLE );
               rowlhs.push_back( REAL{ 0 } );
               rowrhs.push_back( to_int( token ) - offset );
            }
         }
         else
            last_value = to_int( token );
      }
      else if( token == ";" )
      {
         assert( operator_occurred );
         assert( init_last_value );
         break;
      }
      else
      {
         assert( !operator_occurred );
         std::string& name = token;
         bool is_negated = false;
         if( token[0] == '~' )
         {
            is_negated = true;
            offset += last_value;
            name = token.substr( 1, token.size() );
         }
         auto iterator = colname2idx.find( name );
         int col = iterator->second;
         if( iterator == colname2idx.end() )
         {
            col = nCols;
            add_binary_variable( name );
            coeffobj.push_back( { nCols, REAL{ 0 } } );
         }
         if( !init_last_value )
         {
            fmt::print(
                "PaPILO does not support non-linear pseudo-boolean equations" );
            return ParseKey::kFail;
         }
         entries.push_back(
             { nRows, col, is_negated ? -last_value : last_value } );
         nnz++;
         init_last_value = false;
      }
   }
   nRows++;
   assert( operator_occurred );
   assert( init_last_value );
   assert( rowlhs.size() == rowrhs.size() );
   assert( rowlhs.size() == row_flags.size() );
   assert( rowlhs.size() == row_type.size() );
   assert( rowlhs.size() == rownames.size() );
   assert( rowlhs.size() == rowname2idx.size() );
   return ParseKey::kNone;
}

template <typename REAL>
ParseKey
OpbParser<REAL>::parseObjective( const std::string& line )
{
   auto tokens = split( line.c_str() );
   bool is_last_token_numeric = false;
   int counter = 0;
   REAL val;
   for( String token : tokens )
   {
      // TODO: token ~
      if( token == ( "min:" ) )
      {
         assert( counter == 0 );
         counter++;
      }
      else if( token == ";" )
         break;
      else if( isNumeric( token ) )
      {
         assert( !is_last_token_numeric );
         is_last_token_numeric = true;
         counter++;
         val = to_int( token );
      }
      else
      {
         assert( colname2idx.empty() ||
                 colname2idx.find( token ) == colname2idx.end() );

         bool is_negated = false;
         if( token[0] == '~' )
         {
            coeffobj.push_back( { nCols, -val } );
            objoffset += val;
            add_binary_variable( token.substr( 1, token.size() ) );
         }
         else
         {
            coeffobj.push_back( { nCols, val } );
            add_binary_variable( token );
         }
         counter++;
         if( !is_last_token_numeric )
         {
            fmt::print(
                "PaPILO does not support non-linear pseudo-boolean equations" );
            return ParseKey::kFail;
         }
         assert( coeffobj.size() == colnames.size() );
         is_last_token_numeric = false;
      }
   }
   assert( !is_last_token_numeric );
   return ParseKey::kNone;
}

template <typename REAL>
REAL
OpbParser<REAL>::to_int( const String& token ) const
{
   int val = std::stoi( token );
   assert( std::stof( token ) == val );
   return REAL{ (double)val };
}

template <typename REAL>
void
OpbParser<REAL>::add_binary_variable( const String& name )
{
   colnames.push_back( name );
   colname2idx.emplace( name, nCols );
   lb4cols.push_back( REAL{ 0 } );
   ub4cols.push_back( REAL{ 1 } );
   ColFlags flags{};
   flags.set( ColFlag::kIntegral );
   flags.unset( ColFlag::kUbInf );
   flags.unset( ColFlag::kLbInf );
   col_flags.push_back( flags );
   nCols++;
}

template <typename REAL>
Vec<String>
OpbParser<REAL>::split( const char* str )
{
   Vec<String> tokens;
   char c1 = ' ';
   char c2 = '\t';

   do
   {
      const char* begin = str;

      while( *str != c1 && *str != c2 && *str )
         str++;

      tokens.emplace_back( begin, str );

      while( ( *str == c1 || *str == c2 ) && *str )
         str++;

   } while( 0 != *str );

   return tokens;
}

template <typename REAL>
bool
OpbParser<REAL>::isNumeric( const String& token )
{
   return std::regex_match( token,
                            std::regex( "(\\+|-)?[0-9]*(\\.?([0-9]+))$" ) );
}

} // namespace papilo

#endif /* _PARSING_MPS_PARSER_HPP_ */
