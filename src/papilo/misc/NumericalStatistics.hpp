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

      REAL minabsval = 0.0;
      REAL maxabsval = 0.0;
      REAL maxRowDyn = 0.0;
      REAL maxColDyn = 0.0;

      stats.boundsMaxInf = false;
      stats.boundsMax = 0.0;
      stats.boundsMin = 0.0;


      for( int r = 0; r < nrows; ++r )
      {
         const SparseVectorView<REAL>& row = cm.getRowCoefficients(r);
         std::pair<REAL,REAL> minmax = row.getMinMaxAbsValue();

         maxabsval = std::max( minmax.second, maxabsval);
         if( r == 0 ) minabsval = maxabsval;
         else minabsval = std::min( minmax.first, minabsval);

         REAL dyn = minmax.second / minmax.first;
         maxRowDyn = std::max( dyn , maxRowDyn );
      }

      for( int c = 0; c < ncols; ++c )
      {
         // Column dynamism
         const SparseVectorView<REAL>& col = cm.getColumnCoefficients(c);
         std::pair<REAL,REAL> minmax = col.getMinMaxAbsValue();

         REAL dyn = minmax.second / minmax.first;
         maxColDyn = std::max( dyn, maxColDyn );

         // Bounds
         if( c == 0 )
            stats.boundsMin = std::min( abs( vd.lower_bounds[c] ), abs( vd.upper_bounds[c] ) );
         else
            stats.boundsMin = std::min( stats.boundsMin,
                                        REAL( std::min( abs( vd.lower_bounds[c] ), abs( vd.upper_bounds[c] ) ) )
                                        );

         if( !stats.boundsMaxInf )
         {
            if( vd.flags[c].test( ColFlag::kLbInf ) || vd.flags[c].test( ColFlag::kUbInf ) )
               stats.boundsMaxInf = true;
            else
               stats.boundsMax = std::max( stats.boundsMax,
                                           REAL( std::min( abs( vd.lower_bounds[c]), abs( vd.upper_bounds[c] ) ) )
                                           );
         }

      }

      stats.matrixMin = minabsval;
      stats.matrixMax = maxabsval;
      stats.dynamism = maxabsval/minabsval;
      stats.rowDynamism = maxRowDyn;
      stats.colDynamism = maxColDyn;

      // Objective
      const Objective<REAL>& obj = prob.getObjective();

      REAL maxAbsObj = 0.0;
      REAL minAbsObj = 0.0; // == 0 beachten
      for( int i = 0; i < obj.coefficients.size(); ++i )
      {
         maxAbsObj = std::max( REAL( abs( obj.coefficients[i] ) ), maxAbsObj );
         if( i == 0 ) minAbsObj = abs( obj.coefficients[i] );
         else minAbsObj = std::min( REAL ( abs( obj.coefficients[i] ) ), minAbsObj );
      }
      stats.objMax = maxAbsObj;
      stats.objMin = minAbsObj;

      fmt::print("max {}, min {}, dyn {}, bounds: [{},{}], obj: [{},{}]",
                 double(maxabsval),
                 double(minabsval),
                 double(stats.dynamism),
                 double(stats.boundsMin),
                 double(stats.boundsMax),
                 double(stats.objMin),
                 double(stats.objMax)
                 );

   }


private:
   Num_stats<REAL> stats;
   Problem<REAL>& prob;
};

} // namespace papilo

#endif
