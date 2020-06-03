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

// Returns True if cols in given permutation are same
static bool
check_cols( const VariableDomains<double> vd1, const VariableDomains<double> vd2, Vec<int> perm1, const Vec<int> perm2, int ncols )
{

   for( int i = 0; i < ncols; ++i )
   {
      int i1 = perm1[i];
      int i2 = perm2[i];

      if( vd1.flags[i1].test( ColFlag::kIntegral ) !=
          vd2.flags[i2].test( ColFlag::kIntegral ) )
      {
         // kein duplikat: eine variable ist ganzzahlig die andere nicht
         return false;
      }

      if( vd1.flags[i1].test( ColFlag::kUbInf ) !=
          vd2.flags[i2].test( ColFlag::kUbInf ) )
      {
         // kein duplikat: ein upper bound ist +infinity, der andere nicht
         return false;
      }

      if( vd1.flags[i1].test( ColFlag::kLbInf ) !=
          vd2.flags[i2].test( ColFlag::kLbInf ) )
      {
         // kein duplikat: ein lower bound ist -infinity, der andere nicht
         return false;
      }

      if( !vd1.flags[i1].test( ColFlag::kLbInf ) &&
          vd1.lower_bounds[i1] != vd2.lower_bounds[i2] )
      {
         assert( !vd2.flags[i2].test( ColFlag::kLbInf ) );
         // kein duplikat: lower bounds sind endlich aber unterschiedlich
         return false;
      }

      if( !vd1.flags[i1].test( ColFlag::kUbInf ) &&
          vd1.upper_bounds[i1] != vd2.upper_bounds[i2] )
      {
         assert( !vd2.flags[i2].test( ColFlag::kUbInf ) );
         // kein duplikat: upper bounds sind endlich aber unterschiedlich
         return false;
      }
   }
   return true;
}

// Returns True if rows are in
static bool
check_rows( const ConstraintMatrix<double> cm1, const ConstraintMatrix<double> cm2, Vec<int> permrow1, const Vec<int> permrow2, int nrows )
{

   for( int i = 0; i < nrows; ++i)
   {
      int i1 = permrow1[i];
      int i2 = permrow2[i];

   }
}

static bool
check_duplicates( const Problem<double>& prob1, const Problem<double>& prob2 )
{
   // Check for columns
   int ncols = prob1.getNCols();

   if( ncols != prob2.getNCols() )
   {
      // kein duplikat: unterschiedlich viele variablen
      fmt::print("not same number of variables\n");
      return false;
   }

   const VariableDomains<double>& vd1 = prob1.getVariableDomains();
   const VariableDomains<double>& vd2 = prob2.getVariableDomains();

   Vec<int> nocolperm(ncols);
   // std::iota(noperm.begin(), noperm.end(), 0);
   std::generate(nocolperm.begin(), nocolperm.end(), [] {
      static int i = 0;
      return i++;
   });

   check_cols(vd1, vd2, nocolperm, nocolperm, ncols);

   // Check for rows
   // First assume for being same you need to have same rows (even though not true)
   int nrows = prob1.getNRows();

   if( nrows = prob2.getNRows() )
   {
      fmt::print("not same number of rows");
      return false;
   }

   Vec<int> norowperm(nrows);
   std::generate(norowperm.begin(), norowperm.end(), [] {
      static int i = 0;
      return i++;
   });

   const ConstraintMatrix<double> cm1 = prob1.getConstraintMatrix();
   const ConstraintMatrix<double> cm2 = prob2.getConstraintMatrix();

   check_rows(cm1, cm2, norowperm, norowperm, nrows);

   // All checks passed
   return true;
}

int
main( int argc, char* argv[] )
{

   assert( argc == 3 );

   Problem<double> prob1 = MpsParser<double>::loadProblem( argv[1] );
   Problem<double> prob2 = MpsParser<double>::loadProblem( argv[2] );

   fmt::print( "duplicates: \n");
   fmt::print( "{}", check_duplicates( prob1, prob2 ) );

   return 0;
}
