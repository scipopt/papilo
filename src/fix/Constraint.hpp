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

#ifndef _PAPILO_CORE_CONSTRAINT_HPP_
#define _PAPILO_CORE_CONSTRAINT_HPP_

#include "papilo/core/Problem.hpp"
#include "papilo/misc/Num.hpp"
#include "papilo/misc/Hash.hpp"

namespace papilo
{

template <typename REAL>
class Constraint
{

 public:
   const SparseVectorView<REAL>&
   get_data() const
   {
      return data;
   }
   const RowFlags&
   get_row_flag() const
   {
      return row_flag;
   }
   REAL
   get_lhs() const
   {
      return lhs;
   }
   REAL
   get_rhs() const
   {
      return rhs;
   }

   unsigned int
   get_hash(){

      Hasher<unsigned int> hasher( data.getLength() );
      if( data.getLength() > 1 )
      {
         REAL scale = REAL( 2.0 / ( 1.0 + sqrt( 5.0 ) ) ) / data.getValues()[0];
         for( int j = 1; j < data.getLength(); ++j )
            hasher.addValue( Num<REAL>::hashCode( data.getValues()[j] * scale ) );
      }

      return hasher.getHash();
   }

 private:
   SparseVectorView<REAL> data;
   RowFlags row_flag;
   REAL lhs;
   REAL rhs;

 public:
   Constraint( SparseVectorView<REAL> data_, RowFlags row_flag_, REAL lhs_,
               REAL rhs_ )
       : data( data_ ), row_flag( row_flag_ ), lhs( lhs_ ), rhs( rhs_ )
   {
      assert( !row_flag.test( RowFlag::kEquation ) || lhs == rhs );
   }
};

} // namespace papilo
#endif
