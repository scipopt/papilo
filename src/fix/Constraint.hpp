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

namespace papilo
{

template <typename REAL>
class Constraint
{

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
  SparseVectorView<REAL>
  get_data()
  {
    return data;
  }

  RowFlags
  get_row_flags()
  {
    return row_flag;
  }
  REAL
  get_lhs()
  {
    return lhs;
  }
  REAL
  get_rhs()
  {
    return rhs;
  }
};

} // namespace papilo
#endif
