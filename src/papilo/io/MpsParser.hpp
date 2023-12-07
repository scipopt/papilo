/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*               This file is part of the program and library                */
/*    PaPILO --- Parallel Presolve for Integer and Linear Optimization       */
/*                                                                           */
/* Copyright (C) 2020-2023 Konrad-Zuse-Zentrum                               */
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

#ifndef _PAPILO_IO_MPS_PARSER_HPP_
#define _PAPILO_IO_MPS_PARSER_HPP_

#include "BoundType.hpp"
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

/// Parser for mps files in fixed and free format
template <typename REAL>
class MpsParser
{
   static_assert(
       num_traits<typename RealParseType<REAL>::type>::is_floating_point,
       "the parse type must be a floating point type" );

 public:
   static boost::optional<Problem<REAL>>
   loadProblem( const std::string& filename )
   {
      MpsParser<REAL> parser;

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
      problem.setVariableDomains( std::move( parser.lb4cols ),
                                  std::move( parser.ub4cols ),
                                  std::move( parser.col_flags ) );
      problem.setVariableNames( std::move( parser.colnames ) );
      problem.setName( std::move( filename ) );
      problem.setConstraintNames( std::move( parser.rownames ) );

      problem.set_problem_type( ProblemFlag::kMixedInteger );
      if(problem.getNumIntegralCols() == 0 )
         problem.set_problem_type( ProblemFlag::kLinear );
      if(problem.getNumContinuousCols() == 0 )
      {
         //TODO check if PseudoBoolean
         problem.set_problem_type( ProblemFlag::kInteger );
      }

      problem.setInputTolerance(
          REAL{ pow( typename RealParseType<REAL>::type{ 10 },
                     -std::numeric_limits<
                         typename RealParseType<REAL>::type>::digits10 ) } );
      return problem;
   }

 private:
   MpsParser() {}

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

   /// checks first word of strline and wraps it by it_begin and it_end
   ParseKey
   checkFirstWord( std::string& strline, std::string::iterator& it,
                   boost::string_ref& word_ref ) const;

   ParseKey
   parseDefault( boost::iostreams::filtering_istream& file ) const;

   ParseKey
   parseRows( boost::iostreams::filtering_istream& file,
              Vec<BoundType>& rowtype );

   ParseKey
   parseCols( boost::iostreams::filtering_istream& file,
              const Vec<BoundType>& rowtype );

   ParseKey
   parseRhs( boost::iostreams::filtering_istream& file );

   ParseKey
   parseRanges( boost::iostreams::filtering_istream& file );

   ParseKey
   parseBounds( boost::iostreams::filtering_istream& file );
};

template <typename REAL>
ParseKey
MpsParser<REAL>::checkFirstWord( std::string& strline,
                                 std::string::iterator& it,
                                 boost::string_ref& word_ref ) const
{
   using namespace boost::spirit;

   it = strline.begin() + strline.find_first_not_of( " " );
   std::string::iterator it_start = it;

   // TODO: Daniel
   qi::parse( it, strline.end(), qi::lexeme[+qi::graph] );

   const std::size_t length = std::distance( it_start, it );

   boost::string_ref word( &( *it_start ), length );

   word_ref = word;

   if( word.front() == 'R' ) // todo
   {
      if( word == "ROWS" )
         return ParseKey::kRows;
      else if( word == "RHS" )
         return ParseKey::kRhs;
      else if( word == "RANGES" )
         return ParseKey::kRanges;
      else
         return ParseKey::kNone;
   }
   else if( word == "COLUMNS" )
      return ParseKey::kCols;
   else if( word == "BOUNDS" )
      return ParseKey::kBounds;
   else if( word == "ENDATA" )
      return ParseKey::kEnd;
   else
      return ParseKey::kNone;
}

template <typename REAL>
ParseKey
MpsParser<REAL>::parseDefault( boost::iostreams::filtering_istream& file ) const
{
   std::string strline;
   getline( file, strline );

   std::string::iterator it;
   boost::string_ref word_ref;
   return checkFirstWord( strline, it, word_ref );
}

template <typename REAL>
ParseKey
MpsParser<REAL>::parseRows( boost::iostreams::filtering_istream& file,
                            Vec<BoundType>& rowtype )
{
   using namespace boost::spirit;

   std::string strline;
   size_t nrows = 0;
   bool hasobj = false;

   while( getline( file, strline ) )
   {
      bool isobj = false;
      std::string::iterator it;
      boost::string_ref word_ref;
      ParseKey key = checkFirstWord( strline, it, word_ref );

      // start of new section?
      if( key != ParseKey::kNone )
      {
         nRows = int( nrows );
         if( !hasobj )
         {
            std::cout << "WARNING: no objective row found" << std::endl;
            rowname2idx.emplace( "artificial_empty_objective", -1 );
         }

         return key;
      }

      if( word_ref.front() == 'G' )
      {
         rowlhs.push_back( REAL{ 0.0 } );
         rowrhs.push_back( REAL{ 0.0 } );
         row_flags.emplace_back( RowFlag::kRhsInf );
         rowtype.push_back( BoundType::kGE );
      }
      else if( word_ref.front() == 'E' )
      {
         rowlhs.push_back( REAL{ 0.0 } );
         rowrhs.push_back( REAL{ 0.0 } );
         row_flags.emplace_back( RowFlag::kEquation );
         rowtype.push_back( BoundType::kEq );
      }
      else if( word_ref.front() == 'L' )
      {
         rowlhs.push_back( REAL{ 0.0 } );
         rowrhs.push_back( REAL{ 0.0 } );
         row_flags.emplace_back( RowFlag::kLhsInf );
         rowtype.push_back( BoundType::kLE );
      }
      // todo properly treat multiple free rows
      else if( word_ref.front() == 'N' )
      {
         if( hasobj )
         {
            rowlhs.push_back( REAL{ 0.0 } );
            rowrhs.push_back( REAL{ 0.0 } );
            RowFlags rowf;
            rowf.set( RowFlag::kLhsInf, RowFlag::kRhsInf );
            row_flags.emplace_back( rowf );
            rowtype.push_back( BoundType::kLE );
         }
         else
         {
            isobj = true;
            hasobj = true;
         }
      }
      else if( word_ref.empty() ) // empty line
         continue;
      else
         return ParseKey::kFail;

      std::string rowname = ""; // todo use ref

      // get row name
      qi::phrase_parse( it, strline.end(), qi::lexeme[+qi::graph], ascii::space,
                        rowname ); // todo use ref

      // todo whitespace in name possible?
      auto ret = rowname2idx.emplace( rowname, isobj ? ( -1 ) : ( nrows++ ) );

      if( !isobj )
         rownames.push_back( rowname );

      if( !ret.second )
      {
         std::cerr << "duplicate row " << rowname << std::endl;
         return ParseKey::kFail;
      }
   }

   return ParseKey::kFail;
}

template <typename REAL>
ParseKey
MpsParser<REAL>::parseCols( boost::iostreams::filtering_istream& file,
                            const Vec<BoundType>& rowtype )
{
   using namespace boost::spirit;

   std::string colname = "";
   std::string strline;
   int rowidx;
   int ncols = 0;
   int colstart = 0;
   bool integral_cols = false;

   auto parsename = [&rowidx, this]( std::string name ) {
      auto mit = rowname2idx.find( name );

      assert( mit != rowname2idx.end() );
      rowidx = mit->second;

      if( rowidx >= 0 )
         this->nnz++;
      else
         assert( -1 == rowidx );
   };

   auto addtuple = [&rowidx, &ncols,
                    this]( typename RealParseType<REAL>::type coeff ) {
      if( rowidx >= 0 )
         entries.push_back(
             std::make_tuple( ncols - 1, rowidx, REAL{ coeff } ) );
      else
         coeffobj.push_back( std::make_pair( ncols - 1, REAL{ coeff } ) );
   };

   while( getline( file, strline ) )
   {
      std::string::iterator it;
      boost::string_ref word_ref;
      ParseKey key = checkFirstWord( strline, it, word_ref );

      // start of new section?
      if( key != ParseKey::kNone )
      {
         if( ncols > 1 )
            pdqsort( entries.begin() + colstart, entries.end(),
                     []( Triplet<REAL> a, Triplet<REAL> b ) {
                        return std::get<1>( b ) > std::get<1>( a );
                     } );

         return key;
      }

      // check for integrality marker
      std::string marker = ""; // todo use ref
      std::string::iterator it2 = it;

      qi::phrase_parse( it2, strline.end(), qi::lexeme[+qi::graph],
                        ascii::space, marker );

      if( marker == "'MARKER'" )
      {
         marker = "";
         qi::phrase_parse( it2, strline.end(), qi::lexeme[+qi::graph],
                           ascii::space, marker );

         if( ( integral_cols && marker != "'INTEND'" ) ||
             ( !integral_cols && marker != "'INTORG'" ) )
         {
            std::cerr << "integrality marker error " << std::endl;
            return ParseKey::kFail;
         }
         integral_cols = !integral_cols;

         continue;
      }

      // new column?
      if( !( word_ref == colname ) )
      {
         if( word_ref.empty() ) // empty line
            continue;

         colname = word_ref.to_string();
         auto ret = colname2idx.emplace( colname, ncols++ );
         colnames.push_back( colname );

         if( !ret.second )
         {
            std::cerr << "duplicate column " << std::endl;
            return ParseKey::kFail;
         }

         assert( lb4cols.size() == col_flags.size() );

         col_flags.emplace_back( integral_cols ? ColFlag::kIntegral
                                               : ColFlag::kNone );

         // initialize with default bounds
         if( integral_cols )
         {
            lb4cols.push_back( REAL{ 0.0 } );
            ub4cols.push_back( REAL{ 1.0 } );
         }
         else
         {
            lb4cols.push_back( REAL{ 0.0 } );
            ub4cols.push_back( REAL{ 0.0 } );
            col_flags.back().set( ColFlag::kUbInf );
         }

         assert( col_flags.size() == lb4cols.size() );

         if( ncols > 1 )
            pdqsort( entries.begin() + colstart, entries.end(),
                     []( Triplet<REAL> a, Triplet<REAL> b ) {
                        return std::get<1>( b ) > std::get<1>( a );
                     } );

         colstart = entries.size();
      }

      assert( ncols > 0 );

      if( !qi::phrase_parse(
              it, strline.end(),
              +( qi::lexeme[qi::as_string[+qi::graph][( parsename )]] >>
                 qi::real_parser<typename RealParseType<REAL>::type>()[(
                     addtuple )] ),
              ascii::space ) )
         return ParseKey::kFail;
   }

   return ParseKey::kFail;
}

template <typename REAL>
ParseKey
MpsParser<REAL>::parseRanges( boost::iostreams::filtering_istream& file )
{
   using namespace boost::spirit;
   std::string strline;
   assert( rowrhs.size() == rowlhs.size() );

   while( getline( file, strline ) )
   {
      std::string::iterator it;
      boost::string_ref word_ref;
      ParseKey key = checkFirstWord( strline, it, word_ref );

      // start of new section?
      if( key != ParseKey::kNone && key != ParseKey::kRanges )
         return key;

      if( word_ref.empty() )
         continue;

      int rowidx;

      auto parsename = [&rowidx, this]( std::string name ) {
         auto mit = rowname2idx.find( name );

         assert( mit != rowname2idx.end() );
         rowidx = mit->second;

         assert( rowidx >= 0 && rowidx < nRows );
      };

      auto addrange = [&rowidx,
                       this]( typename RealParseType<REAL>::type val ) {
         assert( size_t( rowidx ) < rowrhs.size() );

         if( row_type[rowidx] == BoundType::kGE )
         {
            row_flags[rowidx].unset( RowFlag::kRhsInf );
            rowrhs[rowidx] = rowlhs[rowidx] + REAL(abs( val ));
         }
         else if( row_type[rowidx] == BoundType::kLE )
         {
            row_flags[rowidx].unset( RowFlag::kLhsInf );
            rowlhs[rowidx] = rowrhs[rowidx] - REAL(abs( val ));
         }
         else
         {
            assert( row_type[rowidx] == BoundType::kEq );
            assert( rowrhs[rowidx] == rowlhs[rowidx] );
            assert( row_flags[rowidx].test(RowFlag::kEquation) );

            if( val > REAL{ 0.0 } )
            {
               row_flags[rowidx].unset( RowFlag::kEquation );
               rowrhs[rowidx] = rowrhs[rowidx] + REAL( val );
            }
            else if( val < REAL{ 0.0 } )
            {
               rowlhs[rowidx] = rowlhs[rowidx] + REAL( val );
               row_flags[rowidx].unset( RowFlag::kEquation );
            }
         }
      };

      // compulsory part
      if( !qi::phrase_parse(
              it, strline.end(),
              +( qi::lexeme[qi::as_string[+qi::graph][( parsename )]] >>
                 qi::real_parser<typename RealParseType<REAL>::type>()[(
                     addrange )] ),
              ascii::space ) )
         return ParseKey::kFail;

      // optional part todo don't replicate code
      qi::phrase_parse(
          it, strline.end(),
          +( qi::lexeme[qi::as_string[+qi::graph][( parsename )]] >>
             qi::real_parser<typename RealParseType<REAL>::type>()[(
                 addrange )] ),
          ascii::space );
   }

   return ParseKey::kFail;
}

template <typename REAL>
ParseKey
MpsParser<REAL>::parseRhs( boost::iostreams::filtering_istream& file )
{
   using namespace boost::spirit;
   std::string strline;

   while( getline( file, strline ) )
   {
      std::string::iterator it;
      boost::string_ref word_ref;
      ParseKey key = checkFirstWord( strline, it, word_ref );

      // start of new section?
      if( key != ParseKey::kNone && key != ParseKey::kRhs )
         return key;

      if( word_ref.empty() )
         continue;

      int rowidx;

      auto parsename = [&rowidx, this]( std::string name ) {
         auto mit = rowname2idx.find( name );

         assert( mit != rowname2idx.end() );
         rowidx = mit->second;

         assert( rowidx >= -1 );
         assert( rowidx < nRows );
      };

      auto addrhs = [&rowidx, this]( typename RealParseType<REAL>::type val ) {
         if( rowidx == -1 )
         {
            objoffset = -REAL{ val };
            return;
         }
         if( row_type[rowidx] == BoundType::kEq ||
             row_type[rowidx] == BoundType::kLE )
         {
            assert( size_t( rowidx ) < rowrhs.size() );
            rowrhs[rowidx] = REAL{ val };
            row_flags[rowidx].unset( RowFlag::kRhsInf );
         }

         if( row_type[rowidx] == BoundType::kEq ||
             row_type[rowidx] == BoundType::kGE )
         {
            assert( size_t( rowidx ) < rowlhs.size() );
            rowlhs[rowidx] = REAL{ val };
            row_flags[rowidx].unset( RowFlag::kLhsInf );
         }
      };

      // Documentation Link to qi:
      // https://www.boost.org/doc/libs/1_66_0/libs/spirit/doc/html/spirit/qi/tutorials/warming_up.html
      // +: Parse a one or more times
      // lexeme[a]: Disable skip parsing for a, does pre-skipping
      // as_string: Force atomic assignment for string attributes
      // graph: Matches a character based on the equivalent of std::isgraph in the current character set
      if( !qi::phrase_parse(
              it, strline.end(),
              +( qi::lexeme[qi::as_string[+qi::graph][( parsename )]] >>
                 qi::real_parser<typename RealParseType<REAL>::type>()[(
                     addrhs )] ),
              ascii::space ) )
         return ParseKey::kFail;
   }

   return ParseKey::kFail;
}

template <typename REAL>
ParseKey
MpsParser<REAL>::parseBounds( boost::iostreams::filtering_istream& file )
{
   using namespace boost::spirit;
   std::string strline;

   Vec<bool> ub_is_default( lb4cols.size(), true );
   Vec<bool> lb_is_default( lb4cols.size(), true );

   while( getline( file, strline ) )
   {
      std::string::iterator it;
      boost::string_ref word_ref;
      ParseKey key = checkFirstWord( strline, it, word_ref );

      // start of new section?
      if( key != ParseKey::kNone )
         return key;

      if( word_ref.empty() )
         continue;

      bool islb = false;
      bool isub = false;
      bool isintegral = false;
      bool isdefaultbound = false;

      if( word_ref == "UP" ) // lower bound
         isub = true;
      else if( word_ref == "LO" ) // upper bound
         islb = true;
      else if( word_ref == "FX" ) // fixed
      {
         islb = true;
         isub = true;
      }
      else if( word_ref == "MI" ) // infinite lower bound
      {
         islb = true;
         isdefaultbound = true;
      }
      else if( word_ref == "PL" ) // infinite upper bound (redundant)
      {
         isub = true;
         isdefaultbound = true;
      }
      else if( word_ref == "BV" ) // binary
      {
         isintegral = true;
         isdefaultbound = true;
         islb = true;
         isub = true;
      }
      else if( word_ref == "LI" ) // integer lower bound
      {
         islb = true;
         isintegral = true;
      }
      else if( word_ref == "UI" ) // integer upper bound
      {
         isub = true;
         isintegral = true;
      }
      else if( word_ref == "FR" ) // free variable
      {
         islb = true;
         isub = true;
         isdefaultbound = true;
      }
      else
      {
         if( word_ref == "INDICATORS" )
            std::cerr << "PaPILO does not support INDICATORS in the MPS file!!"<< std::endl;
         else
            std::cerr << "unknown bound type " << word_ref << std::endl;
         return ParseKey::kFail;
      }

      // parse over next word
      qi::phrase_parse( it, strline.end(), qi::lexeme[+qi::graph],
                        ascii::space );

      int colidx;

      auto parsename = [&colidx, this]( std::string name ) {
         auto mit = colname2idx.find( name );
         assert( mit != colname2idx.end() );
         colidx = mit->second;
         assert( colidx >= 0 );
      };

      if( isdefaultbound )
      {
         if( !qi::phrase_parse(
                 it, strline.end(),
                 ( qi::lexeme[qi::as_string[+qi::graph][( parsename )]] ),
                 ascii::space ) )
            return ParseKey::kFail;

         if( isintegral ) // binary
         {
            if( islb )
               lb4cols[colidx] = REAL{ 0.0 };
            if( isub )
            {
               col_flags[colidx].unset( ColFlag::kUbInf );
               ub4cols[colidx] = REAL{ 1.0 };
            }
            col_flags[colidx].set( ColFlag::kIntegral );
         }
         else
         {
            if( islb )
               col_flags[colidx].set( ColFlag::kLbInf );
            if( isub )
               col_flags[colidx].set( ColFlag::kUbInf );
         }
         continue;
      }

      if( !qi::phrase_parse(
              it, strline.end(),
              +( qi::lexeme[qi::as_string[+qi::graph][( parsename )]] >>
                 qi::real_parser<typename RealParseType<REAL>::type>()[(
                     [&ub_is_default, &lb_is_default, &colidx, &islb, &isub,
                      &isintegral,
                      this]( typename RealParseType<REAL>::type val ) {
                        if( islb )
                        {
                           lb4cols[colidx] = REAL{ val };
                           lb_is_default[colidx] = false;
                           col_flags[colidx].unset( ColFlag::kLbInf );
                        }
                        if( isub )
                        {
                           ub4cols[colidx] = REAL{ val };
                           ub_is_default[colidx] = false;
                           col_flags[colidx].unset( ColFlag::kUbInf );
                        }

                        if( isintegral )
                           col_flags[colidx].set( ColFlag::kIntegral );

                        if( col_flags[colidx].test( ColFlag::kIntegral ) )
                        {
                           col_flags[colidx].set( ColFlag::kIntegral );
                           if( !islb && lb_is_default[colidx] )
                              lb4cols[colidx] = REAL{ 0.0 };
                           if( !isub && ub_is_default[colidx] )
                              col_flags[colidx].set( ColFlag::kUbInf );
                        }
                     } )] ),
              ascii::space ) )
         return ParseKey::kFail;
   }

   return ParseKey::kFail;
}

template <typename REAL>
bool
MpsParser<REAL>::parseFile( const std::string& filename )
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
MpsParser<REAL>::parse( boost::iostreams::filtering_istream& file )
{
   nnz = 0;
   ParseKey keyword = ParseKey::kNone;
   ParseKey keyword_old = ParseKey::kNone;

   // parsing loop
   while( keyword != ParseKey::kFail && keyword != ParseKey::kEnd &&
          !file.eof() && file.good() )
   {
      keyword_old = keyword;
      switch( keyword )
      {
      case ParseKey::kRows:
         keyword = parseRows( file, row_type );
         break;
      case ParseKey::kCols:
         keyword = parseCols( file, row_type );
         break;
      case ParseKey::kRhs:
         keyword = parseRhs( file );
         break;
      case ParseKey::kRanges:
         keyword = parseRanges( file );
         break;
      case ParseKey::kBounds:
         keyword = parseBounds( file );
         break;
      case ParseKey::kFail:
         break;
      default:
         keyword = parseDefault( file );
         break;
      }
   }

   if( keyword == ParseKey::kFail || keyword != ParseKey::kEnd )
   {
      printErrorMessage( keyword_old );
      return false;
   }

   assert( row_type.size() == unsigned( nRows ) );

   nCols = colname2idx.size();
   nRows = rowname2idx.size() - 1; // subtract obj row

   return true;
}

} // namespace papilo

#endif /* _PARSING_MPS_PARSER_HPP_ */
