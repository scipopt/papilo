/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*               This file is part of the program and library                */
/*    PaPILO --- Parallel Presolve for Integer and Linear Optimization       */
/*                                                                           */
/* Copyright (C) 2020  Konrad-Zuse-Zentrum                                   */
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

#include "papilo/core/ConstraintMatrix.hpp"
#include "papilo/core/Objective.hpp"
#include "papilo/core/Problem.hpp"
#include "papilo/core/VariableDomains.hpp"
#include "papilo/io/MpsParser.hpp"
#include "papilo/misc/fmt.hpp"

using namespace papilo;

static bool
check_duplicates( const Problem<double>& prob1, const Problem<double>& prob2 )
{
   // hier prob1 und prob2 vergleichen
   int ncols = prob1.getNCols();

   if( ncols != prob2.getNCols() )
   {
      // kein duplikat: unterschiedlich viele variablen
      return false;
   }

   const VariableDomains<double>& vd1 = prob1.getVariableDomains();
   const VariableDomains<double>& vd2 = prob2.getVariableDomains();

   for( int i = 0; i < ncols; ++i )
   {
      if( vd1.flags[i].test( ColFlag::kIntegral ) !=
          vd2.flags[i].test( ColFlag::kIntegral ) )
      {
         // kein duplikat: eine variable ist ganzzahlig die andere nicht
         return false;
      }

      if( vd1.flags[i].test( ColFlag::kUbInf ) !=
          vd2.flags[i].test( ColFlag::kUbInf ) )
      {
         // kein duplikat: ein upper bound ist +infinity, der andere nicht
         return false;
      }

      if( vd1.flags[i].test( ColFlag::kLbInf ) !=
          vd2.flags[i].test( ColFlag::kLbInf ) )
      {
         // kein duplikat: ein lower bound ist -infinity, der andere nicht
         return false;
      }

      if( !vd1.flags[i].test( ColFlag::kLbInf ) &&
          vd1.lower_bounds[i] != vd2.lower_bounds[i] )
      {
         assert( !vd2.flags[i].test( ColFlag::kLbInf ) );
         // kein duplikat: lower bounds sind endlich aber unterschiedlich
         return false;
      }

      if( !vd1.flags[i].test( ColFlag::kUbInf ) &&
          vd1.upper_bounds[i] != vd2.upper_bounds[i] )
      {
         assert( !vd2.flags[i].test( ColFlag::kUbInf ) );
         // kein duplikat: upper bounds sind endlich aber unterschiedlich
         return false;
      }
   }

   return false;
}

int
main( int argc, char* argv[] )
{

   assert( argc == 3 );

   Problem<double> prob1 = MpsParser<double>::loadProblem( argv[1] );
   Problem<double> prob2 = MpsParser<double>::loadProblem( argv[2] );

   fmt::print( "duplicates: {}\n", check_duplicates( prob1, prob2 ) );

   return 0;
}
