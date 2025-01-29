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

#ifndef _PAPILO_IO_SOL_PARSER_HPP_
#define _PAPILO_IO_SOL_PARSER_HPP_

#include "papilo/Config.hpp"
#include "papilo/misc/Hash.hpp"
#include "papilo/misc/Num.hpp"
#include "papilo/misc/Vec.hpp"
#include <algorithm>
#include <boost/algorithm/string.hpp>
#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <sstream>

namespace papilo
{

template <typename REAL>
struct SolParser
{

   static bool
   read( const std::string& filename, const Vec<int>& origcol_mapping,
         const Vec<String>& colnames, Vec<REAL>& solution_vector )
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

      HashMap<String, int> nameToCol;

      for( size_t i = 0; i != origcol_mapping.size(); ++i )
      {
         int origcol = origcol_mapping[i];
         nameToCol.emplace( colnames[origcol], i );
      }

      solution_vector.resize( origcol_mapping.size(), REAL{ 0 } );
      String strline;

      skip_header( colnames, in, strline );

      std::pair<bool, REAL> result;
      do
      {
         auto tokens = split( strline.c_str() );
         assert( !tokens.empty() );

         auto it = nameToCol.find( tokens[0] );
         if( it != nameToCol.end() )
         {
            assert( tokens.size() > 1 );
            result = parse_number<REAL>( tokens[1] );
            if( result.first )
            {
               fmt::print("Could not parse solution {}\n", tokens[1]);
               return false;
            }
            solution_vector[it->second] = result.second;
         }
         else if(strline.empty()){}
         else
         {
            fmt::print( stderr,
                        "WARNING: skipping unknown column {} in solution\n",
                        tokens[0] );
         }
      } while( getline( in, strline ) );

      return true;
   }

   static bool
   read_basis( const std::string& filename, const PostsolveStorage<REAL>& ps,
         Vec<VarBasisStatus>& var_basis, Vec<VarBasisStatus>& row_basis )
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

      HashMap<String, int> nameToCol;

      auto var_names = ps.getOriginalProblem().getVariableNames();
      for( size_t i = 0; i != ps.origcol_mapping.size(); ++i )
      {
         int origcol = ps.origcol_mapping[i];
         nameToCol.emplace( var_names[origcol], i );
      }

      HashMap<String, int> nameToRow;

      for( size_t i = 0; i != ps.origrow_mapping.size(); ++i )
      {
         int origcol = ps.origrow_mapping[i];
         nameToRow.emplace( ps.getOriginalProblem().getConstraintNames()[origcol], i );
      }

      var_basis.resize( ps.origcol_mapping.size(),  VarBasisStatus::ON_LOWER);
      for( int i = 0; i < ps.problem.getNCols(); i++ )
      {
         if(ps.problem.getColFlags()[i].test(ColFlag::kUbInf) &&
             ps.problem.getColFlags()[i].test(ColFlag::kLbInf) )
            var_basis[ps.origcol_mapping[i]] = VarBasisStatus::ZERO;
      }
      row_basis.resize( ps.origrow_mapping.size(),  VarBasisStatus::BASIC);
      String strline;

      skip_header( var_names, in, strline );

      do
      {
         auto tokens = split( strline.c_str() );
         assert( !tokens.empty() );
         if(strline.rfind("ENDATA") == 0)
            break;
         auto it = nameToCol.find( tokens[2] );
         if( it != nameToCol.end() )
         {
            assert( tokens.size() > 1 );
            if( tokens[1] == "UL" )
            {
               assert( tokens.size() == 3 );
               var_basis[it->second] = VarBasisStatus::ON_UPPER;
            }
            else if( tokens[1] == "LL" )
            {
               assert( tokens.size() == 3 );
               var_basis[it->second] = VarBasisStatus::ON_LOWER;
            }
            else if( tokens[1] == "XL" )
            {
               var_basis[it->second] = VarBasisStatus::BASIC;
               assert( tokens.size() == 4 );
               auto it_row = nameToRow.find( tokens[3] );
               if( it_row != nameToCol.end() )
                  row_basis[it_row->second] = VarBasisStatus::ON_LOWER;
               else
               {
                  fmt::print( stderr,
                              "WARNING: skipping unknown row {} in solution\n",
                              tokens[3] );
               }

            }
            else if( tokens[1] == "XU" )
            {
               var_basis[it->second] = VarBasisStatus::BASIC;
               assert( tokens.size() == 4 );
               auto it_row = nameToRow.find( tokens[3] );
               if( it_row != nameToCol.end() )
                  row_basis[it_row->second] = VarBasisStatus::ON_UPPER;
               else
                  fmt::print( stderr,
                              "WARNING: skipping unknown row {} in solution\n",
                              tokens[3] );
            }
            else
            {
               fmt::print(
                   stderr,
                   "WARNING: skipping unknown basistype {} in solution\n",
                   tokens[1] );
            }
         }
         else
         {
            fmt::print( stderr,
                        "WARNING: skipping unknown column {} in solution\n",
                        tokens[2] );
         }
      } while( getline( in, strline ) );

      return true;
   }


 private:

   static void
   skip_header( const Vec<String>& colnames,
                boost::iostreams::filtering_istream& filteringIstream,
                String& strline )
   {
      while(getline( filteringIstream, strline ))
      {
         for(const auto & colname : colnames)
          {
            if( strline.rfind( colname ) != ULLONG_MAX)
               return;
         }
      }
   }

   Vec<String> static split( const char* str )
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
};

} // namespace papilo

#endif
