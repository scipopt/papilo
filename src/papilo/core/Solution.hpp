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

#ifndef _PAPILO_CORE_SOLUTION_HPP_
#define _PAPILO_CORE_SOLUTION_HPP_

namespace papilo
{

enum class SolutionType
{
   kPrimal,
   kPrimalDual
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
   Solution() : type( SolutionType::kPrimal ) {}

   Solution( SolutionType type_ ) : type( type_ ) {}

   Solution( SolutionType type_, Vec<REAL> values )
       : type( type_ ), primal( std::move( values ) )
   {
   }

   Solution( Vec<REAL> values )
       : type( SolutionType::kPrimal ), primal( std::move( values ) )
   {
   }

   Solution( Vec<REAL> primal_values, Vec<REAL> dual_col_values,
             Vec<REAL> dual_row_values )
       : type( SolutionType::kPrimalDual ),
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
   PrimalSolution() : Solution<REAL>( SolutionType::kPrimal ) {}
   PrimalSolution( std::vector<REAL> values )
       : Solution<REAL>( SolutionType::kPrimal, values )
   {
   }
};

template <typename REAL>
class PrimalDualSolution : public Solution<REAL>
{
 public:
   std::vector<REAL> col_dual;
   std::vector<REAL> row_dual;

   PrimalDualSolution() : Solution<REAL>( SolutionType::kPrimalDual ) {}
   PrimalDualSolution( std::vector<REAL> primal_values )
       : Solution<REAL>( SolutionType::kPrimal, primal_values )
   {
   }
};
*/

} // namespace papilo

#endif