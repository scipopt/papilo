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

#ifndef _PAPILO_MISC_NUMERICALSTATISTICS_HPP_
#define _PAPILO_MISC_NUMERICALSTATISTICS_HPP_

#include "papilo/core/ConstraintMatrix.hpp"
#include "papilo/core/Problem.hpp"
#include "papilo/misc/fmt.hpp"
#include <cmath>

namespace papilo
{

template <typename REAL>
struct Num_stats
{
   REAL matrixMin;
   REAL matrixMax;
   REAL objMin;
   REAL objMax;
   REAL boundsMin;
   REAL boundsMax;
   bool boundsMaxInf;
   REAL rhsMin;
   REAL rhsMax;
   REAL lhsMin;
   REAL lhsMax;
   REAL dynamism;
   REAL rowDynamism;
   REAL colDynamism;
};

template <typename REAL>
class NumericalStatistics
{
public:
   NumericalStatistics(Problem<REAL>& p)
      : stats(Num_stats<REAL>())
      , prob(p)
   {
      // Set all values in Num_stats

      const ConstraintMatrix<REAL>& cm = prob.getConstraintMatrix();
      const VariableDomains<REAL>& vd = prob.getVariableDomains();

      int nrows = cm.getNRows();
      int ncols = cm.getNCols();

      // matrixMin, matrixMax, dynamism, bounds

      stats.matrixMin = 0.0;
      stats.matrixMax = 0.0;
      stats.rowDynamism = 0.0;


      // Row dynamism, matrixMin/Max
      for( int r = 0; r < nrows; ++r )
      {
         const SparseVectorView<REAL>& row = cm.getRowCoefficients(r);
         std::pair<REAL,REAL> minmax = row.getMinMaxAbsValue();

         stats.matrixMax = std::max( minmax.second, stats.matrixMax);
         if( r == 0 ) stats.matrixMin = stats.matrixMax;
         else stats.matrixMin = std::min( minmax.first, stats.matrixMin);

         REAL dyn = minmax.second / minmax.first;
         stats.rowDynamism = std::max( dyn , stats.rowDynamism );
      }


      stats.colDynamism = 0.0;
      stats.boundsMaxInf = false;
      stats.boundsMax = 0.0;
      stats.boundsMin = 0.0;
      bool boundsMinSet = false;

      // Column dynamism, Variable Bounds
      for( int c = 0; c < ncols; ++c )
      {
         // Column dynamism
         const SparseVectorView<REAL>& col = cm.getColumnCoefficients(c);
         std::pair<REAL,REAL> minmax = col.getMinMaxAbsValue();

         REAL dyn = minmax.second / minmax.first;
         stats.colDynamism = std::max( dyn, stats.colDynamism );

         // Bounds

         // Handle case where first variables are unbounded
         if( !boundsMinSet )
         {
            boundsMinSet = true;
            if( !vd.flags[c].test( ColFlag::kLbInf ) && !vd.flags[c].test( ColFlag::kUbInf ) )
               stats.boundsMin = std::min( abs( vd.lower_bounds[c] ), abs( vd.upper_bounds[c] ) );
            else if( !vd.flags[c].test( ColFlag::kLbInf ) )
               stats.boundsMin = abs( vd.lower_bounds[c] );
            else if( !vd.flags[c].test( ColFlag::kUbInf ) )
               stats.boundsMin = abs( vd.upper_bounds[c] );
            else
               boundsMinSet = false;
         }
         else
         {
            if( !vd.flags[c].test( ColFlag::kLbInf ) && !vd.flags[c].test( ColFlag::kUbInf ) )
               stats.boundsMin = std::min( stats.boundsMin,
                                           REAL( std::min( abs( vd.lower_bounds[c] ),
                                                           abs( vd.upper_bounds[c] )
                                                           ) )
                                           );
            else if( !vd.flags[c].test( ColFlag::kLbInf ) )
               stats.boundsMin = std::min( stats.boundsMin,
                                           REAL( abs( vd.lower_bounds[c] ) )
                                           );
            else if( !vd.flags[c].test( ColFlag::kUbInf ) )
               stats.boundsMin = std::min( stats.boundsMin,
                                           REAL( abs( vd.upper_bounds[c] ) )
                                           );
         }

         if( !stats.boundsMaxInf )
         {
            if( vd.flags[c].test( ColFlag::kLbInf ) || vd.flags[c].test( ColFlag::kUbInf ) )
               stats.boundsMaxInf = true;
            else
               stats.boundsMax = std::max( stats.boundsMax,
                                           REAL( std::max( abs( vd.lower_bounds[c]),
                                                           abs( vd.upper_bounds[c] )
                                                           ) )
                                           );
         }

      }

      stats.dynamism = stats.matrixMax/stats.matrixMin;

      // Objective
      const Objective<REAL>& obj = prob.getObjective();

      stats.objMax = 0.0;
      stats.objMin = 0.0;
      bool objMinSet = false;

      for( int i = 0; i < obj.coefficients.size(); ++i )
      {
         if( obj.coefficients[i] != 0 )
         {
            stats.objMax = std::max( stats.objMax,
                                     REAL( abs( obj.coefficients[i] ) )
                                     );
            if( !objMinSet )
            {
               stats.objMin = abs( obj.coefficients[i] );
               objMinSet = true;
            }
            else
               stats.objMin = std::min( stats.objMin,
                                        REAL ( abs( obj.coefficients[i] ) )
                                        );
         }
      }

      fmt::print("Matrix[{},{}], bounds: [{},{}], obj: [{},{}], dyn: {} dynCol: {}, dynRow: {}",
                 double(stats.matrixMax),
                 double(stats.matrixMin),
                 double(stats.boundsMin),
                 double(stats.boundsMax),
                 double(stats.objMin),
                 double(stats.objMax),
                 double(stats.dynamism),
                 double(stats.colDynamism),
                 double(stats.rowDynamism)
                 );

   }


private:
   Num_stats<REAL> stats;
   Problem<REAL>& prob;
};

} // namespace papilo

#endif
