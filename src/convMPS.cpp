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
#include "papilo/misc/Hash.hpp"
#include "papilo/misc/Vec.hpp"
#include "papilo/misc/fmt.hpp"
#include "papilo/misc/tbb.hpp"
#include "pdqsort/pdqsort.h"
#include "tbb/concurrent_unordered_set.h"
#include <algorithm>

using namespace papilo;

static void
convMPS( const Problem<double>& prob )
{
   /// Print all the relevant stuff for problem creation as cpp commands

   // Get all relevant data
   int nCols = prob.getNCols();
   int nRows = prob.getNRows();
   const Objective<double>& obj = prob.getObjective();
   const ConstraintMatrix<double>& cm = prob.getConstraintMatrix();
   Vec<double> rowlhs = cm.getLeftHandSides();
   Vec<double> rowrhs = cm.getRightHandSides();
   Vec<RowFlags> row_flags = cm.getRowFlags();
   const int nnz = cm.getNnz();
   const VariableDomains<double> vd = prob.getVariableDomains();

   // Data structures
   fmt::print( "   // enum declaration, only needed once\n" );
   fmt::print( "   enum class boundtype{{ kLE, kEq, kGE }};\n" );

   // Variables
   fmt::print( "   // Variable declaration\n" );
   fmt::print( "   int nCols = {}; int nRows = {};\n", nCols, nRows );
   fmt::print( "   Vec<std::string> rownames;\n" );
   fmt::print( "   Vec<std::string> colnames;\n\n" );
   fmt::print( "   HashMap<std::string, int> rowname2idx;\n" );
   fmt::print( "   HashMap<std::string, int> colname2idx;\n" );
   fmt::print( "   Vec<double> lb4cols;\n" );
   fmt::print( "   Vec<double> ub4cols;\n" );
   fmt::print( "   Vec<boundtype> row_type;\n" );
   fmt::print( "   Vec<RowFlags> row_flags({});\n", nRows );
   fmt::print( "   Vec<ColFlags> col_flags;\n" );
   fmt::print( "   Problem<double> problem;\n" );

   // Objective
   fmt::print( "   // Objective\n" );
   fmt::print( "   Vec<double> coeffobj{{ " );
   for( double coeff: obj.coefficients )
      fmt::print ("{},", coeff );
   fmt::print( "}};\n" );
   fmt::print( "   problem.setObjective( coeffobj, {} );\n\n", obj.offset );

   // Constraint Matrix
   fmt::print( "   // Constraint Matrix\n" );
   fmt::print( "   Vec<std::tuple<int, int, double>> entries{{" );
   // iterate through every row and add columns
   for( int r = 0; r < nRows; ++r )
   {
      const SparseVectorView<double>& rc = cm.getRowCoefficients( r );
      const int len = rc.getLength();
      const int* indices = rc.getIndices();
      const double* vals = rc.getValues();
      for(int i = 0; i < len; ++i)
         fmt::print( "{{{},{},{}}},", r, *(indices + i), *(vals + i) );
   }
   fmt::print( "}};\n" );
   // iterate through every row and set lhs and rhs
   fmt::print( "   Vec<double> rowlhs{{" );
   for( int r = 0; r < nRows; ++r )
      fmt::print( "{},", rowlhs[r] );
   fmt::print( "   }};\n" );
   fmt::print( "   Vec<double> rowrhs{{" );
   for( int r = 0; r < nRows; ++r )
      fmt::print( "{},", rowrhs[r] );
   fmt::print( "   }};\n" );
   // Rowflags
   fmt::print( "   " );
   for( int r = 0; r < nRows; ++r )
   {
      fmt::print( "row_flags[{}].set(", r );

      if( row_flags[r].test( RowFlag::NONE ) )
         fmt::print( "RowFlag::NONE" );
      else
      {
         bool somethingset = false;
         if( row_flags[r].test( RowFlag::kLhsInf ) )
         {
            fmt::print( "RowFlag::kLhsInf" );
            somethingset = true;
         }
         if( row_flags[r].test( RowFlag::kRhsInf ) )
         {
            somethingset ? fmt::print( ",RowFlag::kRhsInf" ) : fmt::print( "RowFlag::kRhsInf" );
            somethingset = true;
         }
         if( row_flags[r].test( RowFlag::kEquation ) )
         {
            somethingset ? fmt::print( ",RowFlag::kEquation" ) : fmt::print( "RowFlag::kEquation" );
            somethingset = true;
         }
         if( row_flags[r].test( RowFlag::kIntegral ) )
         {
            somethingset ? fmt::print( ",RowFlag::kIntegral" ) : fmt::print( "RowFlag::kIntegral" );
            somethingset = true;
         }
         assert( somethingset == true );
      }

      fmt::print( ");" ); // no newline so the file does not get clustered with rowflags...
   }
   fmt::print( "\n" );
   fmt::print( "   problem.setConstraintMatrix( SparseStorage<double>{{ entries, nCols, nRows, false }} , rowlhs, rowrhs, row_flags, false );" );

   // Version 2: The problem builder
   fmt::print( "\n\n\n   ///PROBLEM BUILDER CODE\n" );
   fmt::print( "   //int nCols = {}; int nRows = {};\n", nCols, nRows );
   fmt::print( "   ProblemBuilder<double> pb;\n" );

   // Set all needed things
   fmt::print( "   pb.reserve( {},{},{} );\n", nnz, nRows, nCols );
   // Obj
   fmt::print( "   //Vec<double> coeffobj{{ " );
   for( double coeff: obj.coefficients )
      fmt::print ("{},", coeff );
   fmt::print( "}};\n" );
   fmt::print( "   pb.setObjAll( coeffobj );\n" );
   fmt::print( "   pb.setObjOffset( {} );\n", obj.offset );
   // Columns
   fmt::print( "   Vec<double> lbs = {{" );
   for( int c = 0; c < nCols; ++c )
      fmt::print( "{},", vd.lower_bounds[c] );
   fmt::print( "}};\n" );
   fmt::print( "   Vec<bool> lbInf = {{" );
   for( int c = 0; c < nCols; ++c )
      fmt::print( "{},", vd.flags[c].test( ColFlag::kLbInf ) );
   fmt::print( "}};\n" );
   fmt::print( "   Vec<double> ubs = {{" );
   for( int c = 0; c < nCols; ++c )
      fmt::print( "{},", vd.upper_bounds[c] );
   fmt::print( "}};\n" );
   fmt::print( "   Vec<bool> ubInf = {{" );
   for( int c = 0; c < nCols; ++c )
      fmt::print( "{},", vd.flags[c].test( ColFlag::kUbInf ) );
   fmt::print( "}};\n" );
   fmt::print( "   Vec<bool> isIntegral = {{" );
   for( int c = 0; c < nCols; ++c )
      fmt::print( "{},", vd.flags[c].test( ColFlag::kIntegral ) );
   fmt::print( "}};\n" );
   fmt::print( "   pb.setColLbAll( lbs );\n" );
   fmt::print( "   pb.setColLbInfAll( lbInf );\n" );
   fmt::print( "   pb.setColUbAll( ubs );\n" );
   fmt::print( "   pb.setColUbInfAll( ubInf );\n" );
   fmt::print( "   pb.setColIntegralAll( isIntegral );\n" );
   // Rows


}

int
main( int argc, char* argv[] )
{
   if( argc != 2 )
   {
      fmt::print( "usage:\n" );
      fmt::print( "./convMPS instance1.mps         - create array of cpp code to load instance.mps to papilo\n" );
      return 1;
   }
   assert( argc == 2 );

   Problem<double> prob = MpsParser<double>::loadProblem( argv[1] );

   convMPS( prob );

   return 0;
}
