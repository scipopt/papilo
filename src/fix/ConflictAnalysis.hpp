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

#include "papilo/io/Message.hpp"
#include "papilo/misc/Num.hpp"
#include "papilo/core/RowFlags.hpp"
#include <cassert>
#include <cmath>
#include <fstream>

#include <cassert>
#include <fstream>

namespace papilo
{

template <typename REAL>
class ConflictAnalysis
{
   Message msg;
   Num<REAL> num;

 public:
   ConflictAnalysis( Message _msg, Num<REAL> _num ) : msg( _msg ), num( _num )
   {
   }

   bool
   perform_conflict_analysis( Vec<int> length, Vec<int*> indices,
                              Vec<REAL*> values, Vec<RowFlags> flags,
                              Vec<REAL> lhs, Vec<REAL> rhs )
   {
      // TODO: to be implemented
      //  should return a list of constraint to be added to the builder
      return true;
   }

   bool
   perform_conflict_analysis( )
   {
      msg.info("function call is dummy and waited to be implemented above");
      return true;
   }
};

} // namespace papilo
