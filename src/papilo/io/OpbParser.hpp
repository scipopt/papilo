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

 public:
   static boost::optional<Problem<REAL>>
   loadProblem( const std::string& filename )
   {
      OpbParser<REAL> parser;

      Problem<REAL> problem;

      if( !parser.parseFile( filename ) )
         return boost::none;

      assert( parser.nnz >= 0 );

      assert(static_cast<int>(parser.coeffobj.size()) == parser.nCols);

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

      problem.set_problem_type( ProblemFlag::kMixedInteger );
      problem.set_problem_type( ProblemFlag::kInteger );
      problem.set_problem_type( ProblemFlag::kBinary );

      return problem;
   }

 private:
   OpbParser() = default;

   /// load LP from MPS file as transposed triplet matrix
   bool
   parseFile( const std::string& filename );

   bool
   parse( boost::iostreams::filtering_istream& file );

   /*
    * data for opb problem
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
   parseRows( std::string& line );

   ParseKey
   parseObjective( std::string& line );

   void
   add_binary_variable( const String& name );
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

      if( strline[0] == '*' || strline.empty() )
         continue;
      for (char& c : strline)
         if (c == ';') c = ' ';
      if( strline.substr(0, 4) == "min:" )
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
   assert( nCols == static_cast<int>(colname2idx.size()) );
   assert( nRows == static_cast<int>(rowname2idx.size()) );

   return true;
}

template <typename REAL>
ParseKey
OpbParser<REAL>::parseRows( std::string& line )
{
   rownames.push_back( std::to_string( nRows ) );
   rowname2idx.insert( { std::to_string( nRows ), nRows } );

   unsigned long pos = line.find(">=");
   std::string line_rhs;
   REAL offset = 0;
   if( pos == std::string::npos )
   {
      pos = line.find( '=' );
      row_type.push_back( BoundType::kEq );
      RowFlags flags{};
      flags.unset( RowFlag::kRhsInf );
      flags.unset( RowFlag::kLhsInf );
      flags.set( RowFlag::kEquation );
      row_flags.push_back( flags );
      line_rhs = line.substr( pos + 1 );

   }
   else
   {
      row_type.push_back( BoundType::kGE );
      RowFlags flags{};
      flags.set( RowFlag::kRhsInf );
      flags.unset( RowFlag::kLhsInf );
      row_flags.push_back( flags );
      line_rhs = line.substr( pos + 2 );
   }
   boost::trim(line_rhs);
   line = line.substr(0, pos);
   assert( pos != std::string::npos );
   assert( line.find( "<=" ) == std::string::npos );

   std::istringstream is(line);
   std::vector<std::string> tokens;
   std::string tmp;
   while (is >> tmp)
      tokens.push_back(tmp);

   if (tokens.size() % 2 != 0)
   {
      fmt::print(
          "PaPILO does not support non-linear pseudo-boolean equations\n" );
      return ParseKey::kFail;
   }
   for (int i = 0; i < (long long)tokens.size(); i += 2)
      if (find(tokens[i].begin(), tokens[i].end(), 'x') != tokens[i].end())
      {
         fmt::print(
             "PaPILO does not support non-linear pseudo-boolean equations\n" );
         return ParseKey::kFail;
      }

   std::pair<bool, REAL> result;
   for( int counter = 0; counter < (long long)tokens.size(); counter += 2 )
   {
      std::string s_coef = tokens[counter];
      std::string var = tokens[counter + 1];
      result = parse_number<REAL>( s_coef );
      if( result.first )
      {
         fmt::print("Could not parse coefficient {}\n", s_coef);
         return ParseKey::kFail;
      }
      REAL coef = result.second;
      bool negated = false;
      if( !var.empty() && var[0] == '~' )
      {
         negated = true;
         var = var.substr( 1 );
         offset += coef;
      }
      if( var.empty() || var[0] != 'x' )
      {
         fmt::print( "Variable must start with 'x'\n" );
         return ParseKey::kFail;
      }

      auto iterator = colname2idx.find( var );
      int col;
      if( iterator == colname2idx.end() )
      {
         col = nCols;
         add_binary_variable( var );
         coeffobj.push_back( { col, REAL{ 0 } } );
      }
      else
         col = iterator->second;
      entries.push_back( { nRows, col, negated ? -coef : coef } );
      nnz++;
   }
   result = parse_number<REAL>( line_rhs );
   if( result.first )
   {
      fmt::print("Could not parse side {}\n", line_rhs);
      return ParseKey::kFail;
   }
   REAL rhs = result.second;

   if( row_type[row_type.size() - 1] == BoundType::kEq )
   {
      REAL val = rhs - offset;
      rowrhs.push_back( val );
      rowlhs.push_back( val );
   }
   else if( row_type[row_type.size() - 1] == BoundType::kGE )
   {
      rowlhs.push_back( rhs - offset );
      rowrhs.push_back( REAL{ 0 } );
   }
   nRows++;
   assert( rowlhs.size() == rowrhs.size() );
   assert( rowlhs.size() == row_flags.size() );
   assert( rowlhs.size() == row_type.size() );
   assert( rowlhs.size() == rownames.size() );
   assert( rowlhs.size() == rowname2idx.size() );
   return ParseKey::kNone;
}

template <typename REAL>
ParseKey
OpbParser<REAL>::parseObjective( std::string& line )
{
   assert(line.substr(0, 4) == "min:");

   line = line.substr( 4 );

   std::istringstream is(line);
   std::vector<std::string> tokens;
   std::string tmp;
   while (is >> tmp)
      tokens.push_back(tmp);

   if (tokens.size() % 2 != 0)
   {
      fmt::print(
          "PaPILO does not support non-linear pseudo-boolean equations\n" );
      return ParseKey::kFail;
   }
   for (int i = 0; i < (long long)tokens.size(); i += 2)
      if (find(tokens[i].begin(), tokens[i].end(), 'x') != tokens[i].end())
      {
         fmt::print(
             "PaPILO does not support non-linear pseudo-boolean equations\n" );
         return ParseKey::kFail;
      }

   REAL offset = 0;
   std::pair<bool, REAL> result;
   for( int counter = 0; counter < (long long)tokens.size(); counter += 2 )
   {
      std::string s_coef = tokens[counter];
      std::string var = tokens[counter + 1];
      result = parse_number<REAL>( s_coef );
      if( result.first )
      {
         fmt::print("Could not parse objective {}\n", s_coef);
         return ParseKey::kFail;
      }
      REAL coef = result.second;
      bool negated = false;
      if( !var.empty() && var[0] == '~' )
      {
         negated = true;
         var = var.substr( 1 );
         offset += coef;
      }
      if( var.empty() || var[0] != 'x' )
      {
         fmt::print( "Variable must start with 'x'\n" );
         return ParseKey::kFail;
      }

      if( negated )
         objoffset += coef;

      coeffobj.push_back( { nCols, negated ? -coef: coef } );
      assert(colname2idx.find( var ) == colname2idx.end());
      add_binary_variable( var);
   }
   return ParseKey::kNone;
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

} // namespace papilo

#endif /* _PARSING_MPS_PARSER_HPP_ */
