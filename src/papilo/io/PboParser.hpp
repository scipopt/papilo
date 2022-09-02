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

#ifndef _PAPILO_IO_PBO_PARSER_HPP_
#define _PAPILO_IO_PBO_PARSER_HPP_

#include "papilo/Config.hpp"
#include "papilo/core/ConstraintMatrix.hpp"
#include "papilo/core/Objective.hpp"
#include "papilo/core/Problem.hpp"
#include "papilo/core/VariableDomains.hpp"
#include "papilo/misc/Flags.hpp"
#include "papilo/misc/Hash.hpp"
#include "papilo/misc/Num.hpp"
#include "papilo/external/pdqsort/pdqsort.h"
#include <algorithm>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/optional.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/utility/string_ref.hpp>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <memory>
#include <tuple>
#include <utility>

#ifdef PAPILO_USE_BOOST_IOSTREAMS_WITH_BZIP2
#include <boost/iostreams/filter/bzip2.hpp>
#endif
#ifdef PAPILO_USE_BOOST_IOSTREAMS_WITH_ZLIB
#include <boost/iostreams/filter/gzip.hpp>
#endif

namespace papilo
{

template <typename REAL, bool isfptype = num_traits<REAL>::is_floating_point>
struct RealParseType
{
   using type = double;
};

template <typename REAL>
struct RealParseType<REAL, true>
{
   using type = REAL;
};

/// Parser for pbo files derived from MpsParser.hpp
template <typename REAL>
class PboParser
{
   static_assert(
       num_traits<typename RealParseType<REAL>::type>::is_floating_point,
       "the parse type must be a floating point type" ); 
       // TODO replace this so fractional values are required instead.

 public:
   static boost::optional<Problem<REAL>>
   loadProblem( const std::string& filename )
   {
      PboParser<REAL> parser;

      Problem<REAL> problem;

      if( !parser.parseFile( filename ) )
         return boost::none;

      assert( parser.nnz >= 0 );

      Vec<REAL> obj_vec( size_t( parser.nCols ), REAL{ 0.0 } );

      for( auto i : parser.coeffobj )
         obj_vec[i.first] = i.second;

      problem.setObjective( std::move( obj_vec ), parser.objoffset );
      problem.setConstraintMatrix(
          SparseStorage<REAL>{ std::move( parser.entries ), parser.nCols,
                               parser.nRows, true },
          std::move( parser.rowlhs ), std::move( parser.rowrhs ),
          std::move( parser.row_flags ), true );
      problem.setVariableDomains( std::move( Vec<REAL> vect(n, REAL(0)) ),
                                  std::move( Vec<REAL> vect(n, REAL(1)) ),
                                  std::move( Vec<ColFlags> vect(n, kIntegral) ) );
      problem.setVariableNames( std::move( parser.colnames ) );
      problem.setName( std::move( filename ) );
      // We do not have ConstraintNames in PBO
      //problem.setConstraintNames( std::move( parser.rownames ) );

      /* Not sure what to do with InputTolerance
         problem.setInputTolerance(
            REAL{ pow( typename RealParseType<REAL>::type{ 10 },
                       -std::numeric_limits<
                           typename RealParseType<REAL>::type>::digits10 ) } );
       */
      return problem;
   }

 private:
   PboParser() {}

   /// load LP from PBO file as transposed triplet matrix
   bool
   parseFile( const std::string& filename );

   bool
   parse( boost::iostreams::filtering_istream& file );

   /* Try to comply with http://www.cril.univ-artois.fr/PB16/format.pdf 
    * but not rely on competition specific rules
    * data for pbo problem
    */

   Vec<Triplet<REAL>> entries;
   Vec<std::pair<int, REAL>> coeffobj;
   Vec<REAL> rowlhs;
   Vec<REAL> rowrhs;
   Vec<std::string> colnames;

   HashMap<std::string, int> colname2idx;
   Vec<RowFlags> row_flags;
   REAL objoffset = REAL(0);

   int nCols = 0;
   int nRows = 0;
   int nnz = -1;

   Vec<std::pair<int, REAL>> parseRow(std::string& trimmedstrline);

};

template <typename REAL>
bool
PboParser<REAL>::parseFile( const std::string& filename )
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
std::pair<Vec<std::pair<int, REAL>>,REAL> parseRow(std::string& trimmedstrline)
{
   const std::string whitespace = " ";
   auto beginSpace = trimmedstrline.find_first_of(whitespace);
   while (beginSpace != std::string::npos)
   {
        const auto endSpace = trimmedstrline.find_first_not_of(whitespace, beginSpace);
        const auto range = endSpace - beginSpace;

        trimmedstrline.replace(beginSpace, range, fill);

        const auto newStart = beginSpace + fill.length();
        beginSpace = trimmedstrline.find_first_of(whitespace, newStart);
   } // having a nicer string to work with makes me comfortable i get it right. Loop maybe O(n^2)
   Vec<std::pair<int, REAL>> result;

   REAL rhsoff = REAL(0);

   std::stringstream row(beginSpace);

   int degree = 0;
   int variable_index;
   REAL weight;

   // I am unfamiliar with error handling conventions in this code base.

   while(getline(row, token, ' '))
   // You can use space as line break and the getline the result.
   {
      if (token == "+")
      {
         if (degree != 2) variable_index = -1;
         degree = 0;
         result.push_back(std::make_pair(variable_index, weight));
         continue;
      } 
      else if ((token == ">=") || (token == "="))
      {
         if ((degree != 2) || (std::string::npos != getline(row, token, ' '))) variable_index = -1;
         degree = 0;
         result.push_back(std::make_pair(variable_index, weight));
         getline(row, token, ' ');
         std::istringstream(token) >> weight;
         rhsoff += weight;

         break;
      }
      if (degree == 0)
      {
         std::istringstream(token) >> weight;
      } 
      else if (degree == 1)
      {
         if(token.starts_with('~'))
         {
            // lhs <= a*~x = a*(1-x) = a*1 - a*x <=> lhs -a <= -a*x 
            weight = -weight; 
            // weight is initialized here since degree == 0 branch must run before
            rhsoff += weight;
            token.erase(0,1)
         }
         if(colname2idx.count(token) == 0)
         {
            colname2idx.insert(std::pair<std::string,int>(token,nCols++));
            colnames.push_back( colname );
            assert(colnames[nCols-1] != token);

         } 
         variable_index = colname2idx[token];
      }

      degree++;
   }
   if(degree != 2) {
      result[result.begin].first = -1;
      return std::make_pair(uniqresult, rhsoff)
   }
   bool older_var (std::pair<int, REAL> i, std::pair<int, REAL> j) { return (i.first<j.first); }

   std::sort(result.begin,result.end, older_var);
   Vec<std::pair<int, REAL>> uniqresult;

   int last_index = -2;
   
   for (const auto& pair : row)
   {
      if (last_index == pair.first)
      {
         uniqresult[uniqresult.end] += pair.second;
      } else uniqresult.push_back(pair);
   } 

   return std::make_pair(uniqresult, rhsoff)

}

template <typename REAL>
bool
PboParser<REAL>::parse( boost::iostreams::filtering_istream& file )
{
   nnz = 0;
   bool has_objective = false;
   std::pair<Vec<std::pair<int, REAL>>,REAL> unpack_helper;
   // parsing loop
   std::string line;


   while(std::getline(file, line)){
      if (line[0] == '*' || line.empty()) continue;
      if (line[0] == 'm' && line[1] == 'i' && line[2] == 'n' && line[3] == ':')
      {
         const auto strBegin = line.find_first_not_of(" ", 4); 
         const auto strEnd = line.find_last_not_of(" ;");
         // being a bit liberal in what is accepted
         const auto strRange = strEnd - strBegin + 1;

         line = line.substr(strBegin, strRange);
         //[coeffobj, objoffset]
         unpack_helper = parseRow(line);
         coeffobj = unpack_helper.first;
         objoffset = - unpack_helper.second; // not sure about the sign exactly
         for (const auto& pair : coeffobj)
         {  
            if (pair.first == -1) {
               std::cerr << "Objective contains non-linear and currently unsupported constraint or is malformed" << std::endl;
               return false;
            } 
         }

         break; // objective may only be first non comment line
      }
      break;
   }
   do {
      if (line[0] == '*' || line.empty()) continue;

      Vec<std::pair<int, REAL> row;
      int rhs; 
      const auto strBegin = line.find_first_not_of(" "); 
      const auto strEnd = line.find_last_not_of(" ;");
      // being a bit liberal in what is accepted
      const auto strRange = strEnd - strBegin + 1;
      line = line.substr(strBegin, strRange);     

      //[row, lhs] 
      unpack_helper = parseRow(line);
      row = unpack_helper.first;
      lhs = unpack_helper.second;      
      
      for (const auto& pair : row)
      {  
         // This implementation assumes there are no constraints like 
         // a1*x1 + a2*~x1 +a3*x1 < bi as (x1, i, a1) and (x1, i, -a2) and (x1, i, a3) 
         // would all be added to entries which might be bad
      if (pair.first == -1) {
         std::cerr << "The " << nRows <<" constraint contains non-linear and currently unsupported constraint or is malformed" << std::endl;
         return false;
      } 
         entries.push_back(
            std::make_tuple( pair.first, nRows, pair.second ) );
         nnz++;
      }
      // a1 x1 + a2 x2 = b;
      if (line.find("=") != std::string::npos) 
      {
         rowlhs.push_back( lhs );
         rowrhs.push_back( lhs );
         row_flags.emplace_back( RowFlag::kEquation );

      }
      // a1 x1 + a2 x2 >= b;
      // interpret as lhs = b <= a1*x1 + a2*x2 
      else if (line.find(">=") != std::string::npos) 
      {
         rowlhs.push_back( lhs ); 
         // Not sure what to put here as Floating point in does not exist 
         // for rational types i think.
         rowrhs.push_back( REAL{ 0.0 } );
         row_flags.emplace_back( RowFlag::kRhsInf );
      }
      else 
      {
         return false;
      }
      nRows++;
   } while(std::getline(file,line))

   // those asserts might be off by one or something
   assert(nRows == rowname2idx.size());
   assert(nRows == rowlhs.size());
   assert(nRows == rowrhs.size());
   assert(nRows == row_flags.size());

   return true;
}

} // namespace papilo

#endif /* _PARSING_PBO_PARSER_HPP_ */
