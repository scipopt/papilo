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

#ifndef _CORE_SOLUTION_HPP_
#define _CORE_SOLUTION_HPP_

namespace papilo
{

enum class SolutionType
{
   PRIMAL_ONLY,
   PRIMAL_AND_DUAL
};

template <typename REAL>
class Solution
{
 public:
   SolutionType type;
   Vec<REAL> primal;
   Vec<REAL> col_dual;
   Vec<REAL> row_dual;

   // Default type primal only.
   Solution() : type( SolutionType::PRIMAL_ONLY ) {}

   Solution( SolutionType type_ ) : type( type_ ) {}

   Solution( SolutionType type_, Vec<REAL> values )
       : type( type_ ), primal( std::move( values ) )
   {
   }

   Solution( Vec<REAL> values )
       : type( SolutionType::PRIMAL_ONLY ), primal( std::move( values ) )
   {
   }

   Solution( Vec<REAL> primal_values, Vec<REAL> dual_col_values,
             Vec<REAL> dual_row_values )
       : type( SolutionType::PRIMAL_AND_DUAL ),
         primal( std::move( primal_values ) ),
         col_dual( std::move( dual_col_values ) ),
         row_dual( std::move( dual_row_values ) )
   {
   }
};

/*
template <typename REAL>
class PrimalSolution : public Solution<REAL>
{
 public:
   PrimalSolution() : Solution<REAL>( SolutionType::PRIMAL_ONLY ) {}
   PrimalSolution( std::vector<REAL> values )
       : Solution<REAL>( SolutionType::PRIMAL_ONLY, values )
   {
   }
};

template <typename REAL>
class PrimalDualSolution : public Solution<REAL>
{
 public:
   std::vector<REAL> col_dual;
   std::vector<REAL> row_dual;

   PrimalDualSolution() : Solution<REAL>( SolutionType::PRIMAL_AND_DUAL ) {}
   PrimalDualSolution( std::vector<REAL> primal_values )
       : Solution<REAL>( SolutionType::PRIMAL_ONLY, primal_values )
   {
   }
};
*/

} // namespace papilo

#endif