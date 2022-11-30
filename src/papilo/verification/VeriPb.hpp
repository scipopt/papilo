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

#ifndef _PAPILO_CORE_VERI_PB_HPP_
#define _PAPILO_CORE_VERI_PB_HPP_

#include "papilo/core/Problem.hpp"
#include "papilo/core/postsolve/PostsolveType.hpp"
#include "papilo/core/postsolve/ReductionType.hpp"
#include "papilo/misc/MultiPrecision.hpp"
#include "papilo/misc/Num.hpp"
#include "papilo/misc/Vec.hpp"
#include "papilo/misc/compress_vector.hpp"
#include "papilo/misc/fmt.hpp"



namespace papilo
{


/// type to store necessary data for post solve
template <typename REAL>
class VeriPb
{
 public:
   unsigned int nRowsOriginal;


   /// mapping of reduced problems row indices to row indices in the original
   /// problem
   Vec<int> rhs_row_mapping;

   Vec<int> lhs_row_mapping;

   // number of cons in VeriPb
   int veri_pb_rows = 0;

   Num<REAL> num;
   Message msg;

   VeriPb() = default;

   VeriPb( const Problem<REAL>& _problem, const Num<REAL>& _num,
           const Message& _msg )
       : num( _num ), msg( _msg )
   {
      nRowsOriginal = _problem.getNRows();

      rhs_row_mapping.reserve( nRowsOriginal );
      lhs_row_mapping.reserve( nRowsOriginal );

      for( unsigned int i = 0; i < nRowsOriginal; ++i )
      {
         if(!_problem.getRowFlags()[i].test(RowFlag::kLhsInf))
         {
            veri_pb_rows++;
            lhs_row_mapping.push_back( veri_pb_rows );
         }
         if(!_problem.getRowFlags()[i].test(RowFlag::kLhsInf))
         {
            veri_pb_rows++;
            rhs_row_mapping.push_back( veri_pb_rows );
         }
      }
   }

   void
   print_header()
   {
      msg.info( "pseudo-Boolean proof version 1.1\n" );
      msg.info( "f {}\n", veri_pb_rows );
   };

   //TODO: information if dual reduction
   void
   change_upper_bound( REAL val, const String& name )
   {
      // VeriPb can only handle >= constraint and they must start with variables
      // -> invert variable
      veri_pb_rows++;
      assert( val == 0 );
      msg.info("rup 1 ~{} >= 1 ;\n", name, (int) val);
   }

   //TODO: information if dual reduction
   void
   change_lower_bound(  REAL val, const String& name )
   {
      veri_pb_rows++;
      assert( val == 1 );
      msg.info("rup 1 {} >= {} ;\n", name, (int) val);
   }

   //TODO: information if dual reduction
   void
   fix_var( REAL val,  const String& name )
   {
      veri_pb_rows++;
      msg.info("rup 1 {} = {} ;\n", name, (int) val);
   }

   //TODO: compress the mappings
   void
   compress( const Vec<int>& rowmapping, const Vec<int>& colmapping,
             bool full = false )
   {
#ifdef PAPILO_TBB
      tbb::parallel_invoke(
          [this, &rowmapping, full]() {
             // update information about rows that is stored by index
             compress_vector( rowmapping, lhs_row_mapping );
             if( full )
                lhs_row_mapping.shrink_to_fit();
          },
          [this, &rowmapping, full]() {
             // update information about rows that is stored by index
             compress_vector( rowmapping, rhs_row_mapping );
             if( full )
                rhs_row_mapping.shrink_to_fit();
          } );
#else
      compress_vector( rowmapping, lhs_row_mapping );
      compress_vector( rowmapping, rhs_row_mapping );
      if( full )
      {
         rhs_row_mapping.shrink_to_fit();
         lhs_row_mapping.shrink_to_fit();
      }
#endif
   }


 private:


};



} // namespace papilo

#endif
