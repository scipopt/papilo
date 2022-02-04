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

#ifndef SRC_FIX_ANALYZECONFLICT_HPP
#define SRC_FIX_ANALYZECONFLICT_HPP

#include "papilo/core/ProbingView.hpp"

#include <cassert>
#include <fstream>
#include <string>

using namespace papilo;

template <typename REAL>
class AnalyzeConflict
{
   Message msg;
   ProbingView<REAL> probing_view;
   
 public:
   // Constructor
   AnalyzeConflict(Message msg_, ProbingView<REAL> view_)
   : msg( msg_ ), probing_view( view_ )
   {
   }

   void
   analyze_conflict()
   {
        msg.info( "Start conflict analysis \n");
   }

   void
   resolve_bound_change()
   {

   }

   void
   add_conflict_constraint()
   {

   }
};

#endif /* SRC_FIX_ANALYZECONFLICT_HPP */
