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
   REAL rhsMin;
   REAL rhsMax;
   REAL lhsMin;
   REAL lhsMax;
   REAL dynamism;
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

      ConstraintMatrix<REAL>& cm = prob.getConstraintMatrix();

      int nrows = cm.getNRows();
      int ncols = cm.getNCols();

      REAL minabsval;
      REAL maxabsval = 0.0;
      REAL maxdyn = 0.0;
      if( nrows < ncols )
      {
         for( int r = 0; r < nrows; ++r)
         {
            const SparseVectorView<REAL>& row = cm.getRowCoefficients(r);
            // technically you loop 3 times over row -> not efficient
            maxabsval = std::max( row.getMaxAbsValue(), maxabsval );
            maxdyn = std::max( row.getDynamism(), maxdyn);
            if( r == 0 ) minabsval = maxabsval;
            else minabsval = std::min( row.getMinAbsValue(), minabsval);
         }
      }
      else
      {
         for( int c = 0; c < ncols; ++c)
         {
            const SparseVectorView<REAL>& col = cm.getColumnCoefficients(c);
            maxabsval = std::max( col.getMaxAbsValue(), maxabsval );
            maxdyn = std::max( col.getDynamism(), maxdyn);
            if( c == 0 ) minabsval = maxabsval;
            else minabsval = std::min( col.getMinAbsValue(), minabsval);
         }
      }
      stats.matrixMin = minabsval;
      stats.matrixMax = maxabsval;
      stats.dynamism = maxdyn;
      // todo: dynamism for col/row missing sometimes. Gotta fix that
      fmt::print("max {}, min {}, dyn {}", double(maxabsval), double(minabsval), double(maxdyn));

   }


private:
   Num_stats<REAL> stats;
   Problem<REAL>& prob;
};

} // namespace papilo

#endif
