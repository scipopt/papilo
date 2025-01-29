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

#ifndef _PAPILO_IO_PARSER_HPP_
#define _PAPILO_IO_PARSER_HPP_

#include "papilo/Config.hpp"
#include "papilo/io/MpsParser.hpp"
#include "papilo/io/OpbParser.hpp"
#ifdef PAPILO_USE_BOOST_IOSTREAMS_WITH_BZIP2
#include <boost/iostreams/filter/bzip2.hpp>
#endif
#ifdef PAPILO_USE_BOOST_IOSTREAMS_WITH_ZLIB
#include <boost/iostreams/filter/gzip.hpp>
#endif

namespace papilo
{

template <typename REAL>
class Parser
{

 public:
   static boost::optional<Problem<REAL>>
   loadProblem( const std::string& filename )
   {
      if( filename.find(".mps") != std::string::npos)
         return MpsParser<REAL>::loadProblem( filename );
      else if( filename.find(".opb") != std::string::npos)
         return OpbParser<REAL>::loadProblem( filename );
      else
         return boost::none;
   }
};
}; // namespace papilo

#endif /* _PARSING_MPS_PARSER_HPP_ */
