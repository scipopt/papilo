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

#ifndef _PAPILO_VERI_CERTIFICATE_INTERFACE_HPP_
#define _PAPILO_VERI_CERTIFICATE_INTERFACE_HPP_

#include "papilo/core/Problem.hpp"
#include "papilo/misc/Num.hpp"
#include "papilo/misc/Vec.hpp"
#include "papilo/misc/compress_vector.hpp"
#include "papilo/misc/fmt.hpp"
#include "papilo/verification/ArgumentType.hpp"

namespace papilo
{

/// type to store necessary data for post solve
template <typename REAL>
class CertificateInterface
{

 public:
   CertificateInterface() = default;

   virtual void
   print_header() = 0;

   virtual void
   change_upper_bound( REAL val, const String& name,
                       ArgumentType argument = ArgumentType::kPrimal ) = 0;

   virtual void
   change_lower_bound( REAL val, const String& name,
                       ArgumentType argument = ArgumentType::kPrimal ) = 0;

   virtual void
   change_rhs( int row, REAL val, const SparseVectorView<REAL>& data,
               const Vec<String>& names, const Vec<int>& var_mapping ) = 0;

   virtual void
   change_lhs( int row, REAL val, const SparseVectorView<REAL>& data,
               const Vec<String>& names, const Vec<int>& var_mapping ) = 0;

   virtual void
   mark_row_redundant( int row ) = 0;

   virtual void
   update_row( int row, int col, REAL new_val,  const SparseVectorView<REAL>& data,
               RowFlags& rflags, REAL lhs, REAL rhs,
               const Vec<String>& names, const Vec<int>& var_mapping ) = 0;

   virtual void
   substitute( int col, int row,
               const Problem<REAL>& currentProblem, const Vec<String>& names, const Vec<int>& var_mapping ) = 0;

   virtual void
   compress( const Vec<int>& rowmapping, const Vec<int>& colmapping,
             bool full = false ) = 0;

   virtual ~CertificateInterface() = default;

};

} // namespace papilo

#endif
