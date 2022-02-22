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

#include "papilo/io/MpsParser.hpp"
#include "papilo/io/MpsWriter.hpp"
#include "papilo/core/Presolve.hpp"
#include "papilo/core/postsolve/Postsolve.hpp"
#include "papilo/core/postsolve/PostsolveStorage.hpp"
#include "boost/algorithm/string/trim.hpp"
#include <cassert>
#include <string>

using namespace papilo;

void*
presolve( const char* filename, const char* reduced_filename,
          const char* setting_filename, int* status, int* result )
{

   std::string filename_as_string( filename );
   std::string setting_as_string( setting_filename );
   boost::optional<Problem<double>> prob;
   {
      prob = MpsParser<double>::loadProblem( filename_as_string );
   }
   if( !prob )
   {
      fmt::print( "error loading problem {}\n", filename );
      *result = 1;
      return nullptr;
   }

   Presolve<double> presolve{};
   presolve.addDefaultPresolvers();
   ParameterSet paramSet = presolve.getParameters();

   if( !setting_as_string.empty() )
   {
      std::ifstream input( setting_filename );
      if( input )
      {
         String theoptionstr;
         String thevaluestr;
         for( String line; getline( input, line ); )
         {
            std::size_t pos = line.find_first_of( '#' );
            if( pos != String::npos )
               line = line.substr( 0, pos );

            pos = line.find_first_of( '=' );

            if( pos == String::npos )
               continue;

            theoptionstr = line.substr( 0, pos - 1 );
            thevaluestr = line.substr( pos + 1 );

            boost::algorithm::trim( theoptionstr );
            boost::algorithm::trim( thevaluestr );

            try
            {
               paramSet.parseParameter( theoptionstr.c_str(),
                                        thevaluestr.c_str() );
               fmt::print( "set {} = {}\n", theoptionstr, thevaluestr );
            }
            catch( const std::exception& e )
            {
               fmt::print( "parameter '{}' could not be set: {}\n", line,
                           e.what() );
            }
         }
      }
      else
      {
         fmt::print( "could not read parameter file '{}'\n",
                     setting_filename );
      }
   }

   *result = 0;
   auto problem = new Problem<double>( prob.get() );
   auto presolve_result = presolve.apply( *problem, false );
   switch( presolve_result.status )
   {

   case PresolveStatus::kUnchanged:
      *status = 1;
      return nullptr;
   case PresolveStatus::kReduced:
      *status = 0;
      break;
   case PresolveStatus::kUnbndOrInfeas:
      *status = 2;
      return nullptr;
   case PresolveStatus::kUnbounded:
      *status = 3;
      return nullptr;
   case PresolveStatus::kInfeasible:
      *status = 4;
      return nullptr;
   }
   std::string reduced_filename_as_string( reduced_filename );
   MpsWriter<double>::writeProb( reduced_filename_as_string, *problem,
                                 presolve_result.postsolve.origrow_mapping,
                                 presolve_result.postsolve.origcol_mapping );
   return new PostsolveStorage<double>( presolve_result.postsolve );
}

void
delete_postsolve_storage( void* postsolve_storage_ptr )
{
   auto postsolve_storage =
       (PostsolveStorage<double>*)( postsolve_storage_ptr );
   delete postsolve_storage;
}

void
postsolve( void* postsolve_storage_ptr, double* red_solution, int red_cols, double* org_solution,
           int org_cols )
{
   auto postsolve_storage =
       (PostsolveStorage<double>*)( postsolve_storage_ptr );
   assert( postsolve_storage->origcol_mapping.size() == red_cols );
   assert( postsolve_storage->getOriginalProblem().getNCols() == org_cols );
   Vec<double> sol_vec( red_solution, red_solution + red_cols );
   Vec<double> res_vec( org_solution, org_solution + org_cols );
   Solution<double> sol( sol_vec );
   Solution<double> res( res_vec );

   Postsolve<double> postsolve{ {}, postsolve_storage->num };
   PostsolveStatus status = postsolve.undo( sol, res, *postsolve_storage );
   assert( status == PostsolveStatus::kOk );
   org_solution = &res.primal[0];
}
