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

#include "papilo/core/Objective.hpp"
#include "papilo/core/Presolve.hpp"
#include "papilo/misc/OptionsParser.hpp"
#include "papilo/misc/VersionLogger.hpp"

#include <cassert>
#include <fstream>

using namespace papilo;

template <typename REAL>
class VectorMultiplication
{

 public:
   VectorMultiplication() = default;

   Vec<REAL>
   multiplication( const ConstraintMatrix<REAL>& matrix,
                   const Vec<REAL>& scalar, const Vec<REAL>& substract )
   {
      assert( matrix.getNRows() == substract.size() );
      assert( matrix.getNCols() == scalar.size() );
      Vec<REAL> result( substract );
#ifdef PAPILO_TBB
      tbb::parallel_for( tbb::blocked_range<int>( 0, matrix.getNRows() ),
                         [&]( const tbb::blocked_range<int>& r )
                         {
                            for( int i = r.begin(); i < r.end(); ++i )
#else
      for( int i = 0; i < matrix.getNRows(); ++i )
#endif
                            {
                               auto coeff = matrix.getRowCoefficients( i );
                               StableSum<REAL> aux( -substract[i] );
                               for( int j = 0; j < coeff.getLength(); j++ )
                                  aux.add( coeff.getValues()[j] *
                                           scalar[coeff.getIndices()[j]] );
                               result[i] = aux.get();
                            }
#ifdef PAPILO_TBB
                         } );
#endif
      return result;
   }
};
