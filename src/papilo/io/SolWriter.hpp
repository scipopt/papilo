/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*               This file is part of the program and library                */
/*    PaPILO --- Parallel Presolve for Integer and Linear Optimization       */
/*                                                                           */
/* Copyright (C) 2020-2021 Konrad-Zuse-Zentrum                               */
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

#ifndef _PAPILO_IO_SOL_WRITER_HPP_
#define _PAPILO_IO_SOL_WRITER_HPP_

#include "papilo/Config.hpp"
#include "papilo/misc/Vec.hpp"
#include "papilo/misc/fmt.hpp"
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

/// Writer to write problem structures into an mps file
template <typename REAL>
struct SolWriter
{
   static void
   writePrimalSol( const std::string& filename, const Vec<REAL>& sol,
                   const Vec<REAL>& objective, const REAL& solobj,
                   const Vec<std::string>& colnames )
   {
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

      out.push( file );

      fmt::print( out, "{: <50} {: <18.15}\n", "=obj=", double( solobj ) );

      for( int i = 0; i != sol.size(); ++i )
      {
         if( sol[i] != 0.0 )
         {
            fmt::print( out, "{: <50} {: <18.15}   obj({:.15})\n", colnames[i],
                        double( sol[i] ), double( objective[i] ) );
         }
      }
   }

   static void
   writeDualSol( const std::string& filename, const Vec<REAL>& sol,
                 const Vec<REAL>& rhs, const Vec<REAL>& lhs,
                 const REAL& obj_value, const Vec<std::string>& row_names )
   {
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

      out.push( file );

      fmt::print( out, "{: <50} {: <18.15}\n", "=obj=", double( obj_value ) );

      for( int i = 0; i < sol.size(); ++i )
      {
         if( sol[i] != 0.0 )
         {
            REAL objective = lhs[i];
            if( sol[i] < 0 )
               objective = rhs[i];
            fmt::print( out, "{: <50} {: <18.15}   obj({:.15})\n", row_names[i],
                        double( sol[i] ), double( objective ) );
         }
      }
   }

   static void
   writeReducedCostsSol( const std::string& filename, const Vec<REAL>& sol,
                         const Vec<REAL>& ub, const Vec<REAL>& lb,
                         const REAL& solobj, const Vec<std::string>& col_names )
   {
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

      out.push( file );

      fmt::print( out, "{: <50} {: <18.15}\n", "=obj=", double( solobj ) );

      for( int i = 0; i < sol.size(); ++i )
      {
         if( sol[i] != 0.0 )
         {
            REAL objective = lb[i];
            if( sol[i] < 0 )
               objective = ub[i];
            fmt::print( out, "{: <50} {: <18.15}   obj({:.15})\n", col_names[i],
                        double( sol[i] ), double( objective ) );
         }
      }
   }
};

} // namespace papilo

#endif