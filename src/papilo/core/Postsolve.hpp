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

#ifndef _PAPILO_CORE_POSTSOLVE_HPP_
#define _PAPILO_CORE_POSTSOLVE_HPP_

#include "papilo/core/Problem.hpp"
#include "papilo/misc/KktChecker.hpp"
#include "papilo/misc/MultiPrecision.hpp"
#include "papilo/misc/Num.hpp"
#include "papilo/misc/StableSum.hpp"
#include "papilo/misc/Vec.hpp"
#include "papilo/misc/compress_vector.hpp"
#include "papilo/misc/fmt.hpp"
#include "papilo/misc/tbb.hpp"
#include <fstream>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/tmpdir.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/list.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>

namespace papilo
{

/// possible types of post solving
enum class PostsolveType : int
{
   kPrimal = 0,
   kFull = 1,
};

enum class ReductionType : int
{
   kFixedCol = 0,
   kSubstitutedCol = 1,
   kParallelCol = 2,
   kSaveRow = 3,
   kSaveCol = 4
};

enum class PostsolveStatus : int
{
   kOk,
   kFail
};

// forward declarations
template <typename REAL>
class SparseVectorView;

struct IndexRange;

/// type to store necessary data for post solve
template <typename REAL>
class Postsolve
{
 public:
   unsigned int nColsOriginal;
   unsigned int nRowsOriginal;

   /// mapping of reduced problems column indices to column indices in the
   /// original problem
   Vec<int> origcol_mapping;

   /// mapping of reduced problems row indices to row indices in the original
   /// problem
   Vec<int> origrow_mapping;

   // set to full for development of postsolve,
   // later will not be default value
   // PostsolveType postsolveType = PostsolveType::kFull;
   PostsolveType postsolveType = PostsolveType::kPrimal;

   Vec<ReductionType> types;
   Vec<int> indices;
   Vec<REAL> values;
   Vec<int> start;

   Problem<REAL> problem;

   Num<REAL> num;

// #define CHECK_KKT
#ifndef CHECK_KKT
   using Kkt = KktChecker<REAL, CheckLevel::No_check>;
#else
   using Kkt = KktChecker<REAL, CheckLevel::Check>;
#endif

   mutable Kkt checker;

   Postsolve() {}

   Postsolve( const Problem<REAL>& problem, const Num<REAL>& num )
       : problem( problem ), num( num )
   {
      int nrows = problem.getNRows();
      int ncols = problem.getNCols();

      origrow_mapping.reserve( nrows );
      origrow_mapping.reserve( ncols );

      for( int i = 0; i < nrows; ++i )
         origrow_mapping.push_back( i );

      for( int i = 0; i < ncols; ++i )
         origcol_mapping.push_back( i );

      nColsOriginal = ncols;
      nRowsOriginal = nrows;

      start.push_back( 0 );

      // release excess storage in original problem copy
      this->problem.compress( true );
   }

   int
   notifySavedRow( int row, const SparseVectorView<REAL>& coefficients,
                   REAL lhs, REAL rhs, const RowFlags& flags );

   void
   notifyModifiedRow( int row );

   void
   notifyFixedCol( int col, REAL val );

   void
   notifySubstitution( int col, SparseVectorView<REAL> equalityLHS,
                       REAL equalityRHS );

   /// col1 = col2scale * col2 and are merged into a new column y = col2 +
   /// col2scale * col1 which takes over the index of col2
   void
   notifyParallelCols( int col1, bool col1integral, bool col1lbinf,
                       const REAL& col1lb, bool col1ubinf, const REAL& col1ub,
                       int col2, bool col2integral, bool col2lbinf,
                       const REAL& col2lb, bool col2ubinf, const REAL& col2ub,
                       const REAL& col2scale );

   void
   compress( const Vec<int>& rowmapping, const Vec<int>& colmapping,
             bool full = false )
   {
      tbb::parallel_invoke(
          [this, &colmapping, full]() {
             compress_vector( colmapping, origcol_mapping );
             if( full )
                origcol_mapping.shrink_to_fit();
          },
          [this, &rowmapping, full]() {
             // update information about rows that is stored by index
             compress_vector( rowmapping, origrow_mapping );
             if( full )
                origrow_mapping.shrink_to_fit();
          } );
   }

   template <typename Archive>
   void
   serialize( Archive& ar, const unsigned int version )
   {
      ar& nColsOriginal;
      ar& nRowsOriginal;
      ar& origcol_mapping;
      ar& origrow_mapping;
      ar& postsolveType;
      ar& types;
      ar& indices;
      ar& values;
      ar& start;

      ar& problem;

      ar& num;
   }

   PostsolveStatus
   undo( const Solution<REAL>& reducedSolution,
         Solution<REAL>& originalSolution ) const;

   const Problem<REAL>&
   getOriginalProblem() const
   {
      return problem;
   }

   const Num<REAL>&
   getNum() const
   {
      return num;
   }

   Kkt&
   getChecker()
   {
      return checker;
   }

 private:
   void
   finishNotify()
   {
      assert( types.size() == start.size() );
      assert( values.size() == indices.size() );
      start.push_back( values.size() );
   }

   Vec<int> row_stack_index;
};

#ifdef PAPILO_USE_EXTERN_TEMPLATES
extern template class Postsolve<double>;
extern template class Postsolve<Quad>;
extern template class Postsolve<Rational>;
#endif

template <typename REAL>
void
Postsolve<REAL>::notifyFixedCol( int col, REAL val )
{
   types.push_back( ReductionType::kFixedCol );
   indices.push_back( origcol_mapping[col] );
   values.push_back( val );

   finishNotify();
}

// template <typename REAL>
// void
// Postsolve<REAL>::notifySavedRow( int row,
//                                  const SparseVectorView<REAL>& coefficients,
//                                  REAL lhs, REAL rhs, const RowFlags& flags )
// {
//    const REAL* coefs = coefficients.getValues();
//    const int* columns = coefficients.getIndices();
//    const int length = coefficients.getLength();

//    types.push_back( ReductionType::kSaveRow );
//    indices.push_back( origrow_mapping[row] );
//    values.push_back( (double)length );

//    // LB
//    if( flags.test( RowFlag::kLhsInf ) )
//       indices.push_back( 1 );
//    else
//       indices.push_back( 0 );
//    values.push_back( lhs );

//    // UB
//    if( flags.test( RowFlag::kRhsInf ) )
//       indices.push_back( 1 );
//    else
//       indices.push_back( 0 );
//    values.push_back( rhs );

//    for( int i = 0; i < length; ++i )
//    {
//       indices.push_back( columns[i] );
//       values.push_back( coefs[i] );
//    }

//    finishNotify();
// }

template <typename REAL>
void
Postsolve<REAL>::notifyModifiedRow( int row )
{
   int origrow = origrow_mapping[row];
   row_stack_index[origrow] = -1;
}

template <typename REAL>
int
Postsolve<REAL>::notifySavedRow( int row,
                                 const SparseVectorView<REAL>& coefficients,
                                 REAL lhs, REAL rhs, const RowFlags& flags )
{
   // initialize arrays if necessary
   if( row_stack_index.size() == 0 )
   {
      int nrows = problem.getNRows();
      row_stack_index.resize( nrows, -1 );
   }

   // check if row is valid on the postsolve stack
   if( row_stack_index[row] >= 0 )
      return row_stack_index[row];

   const REAL* coefs = coefficients.getValues();
   const int* columns = coefficients.getIndices();
   const int length = coefficients.getLength();

   types.push_back( ReductionType::kSaveRow );

   int stack_index = indices.size();
   indices.push_back( origrow_mapping[row] );
   values.push_back( (double)length );

   // LB
   if( flags.test( RowFlag::kLhsInf ) )
      indices.push_back( 1 );
   else
      indices.push_back( 0 );
   values.push_back( lhs );

   // UB
   if( flags.test( RowFlag::kRhsInf ) )
      indices.push_back( 1 );
   else
      indices.push_back( 0 );
   values.push_back( rhs );

   for( int i = 0; i < length; ++i )
   {
      indices.push_back( columns[i] );
      values.push_back( coefs[i] );
   }

   finishNotify();

   row_stack_index[row] = stack_index;

   return stack_index;
}

template <typename REAL>
void
Postsolve<REAL>::notifySubstitution( int col,
                                     SparseVectorView<REAL> equalityLHS,
                                     REAL equalityRHS )
{
   const REAL* coefs = equalityLHS.getValues();
   const int* columns = equalityLHS.getIndices();
   const int length = equalityLHS.getLength();
   assert( length > 1 );

   types.push_back( ReductionType::kSubstitutedCol );
   values.push_back( equalityRHS );
   // values.insert( values.end(), coefs, coefs + length );
   indices.push_back( origcol_mapping[col] );
   for( int i = 0; i < length; ++i )
   {
      indices.push_back( origcol_mapping[columns[i]] );
      values.push_back( coefs[i] );
   }

   finishNotify();
}

/// col1 = col2scale * col2 and are merged into a new column y = col2 +
/// col2scale * col1 which takes over the index of col2
template <typename REAL>
void
Postsolve<REAL>::notifyParallelCols( int col1, bool col1integral,
                                     bool col1lbinf, const REAL& col1lb,
                                     bool col1ubinf, const REAL& col1ub,
                                     int col2, bool col2integral,
                                     bool col2lbinf, const REAL& col2lb,
                                     bool col2ubinf, const REAL& col2ub,
                                     const REAL& col2scale )
{
   // encode the finiteness of the bounds in one integer and store it as
   // value for column 1
   int col1BoundFlags = 0;
   int col2BoundFlags = 0;

   if( col1integral )
      col1BoundFlags |= static_cast<int>( ColFlag::kIntegral );
   if( col1lbinf )
      col1BoundFlags |= static_cast<int>( ColFlag::kLbInf );
   if( col1ubinf )
      col1BoundFlags |= static_cast<int>( ColFlag::kUbInf );
   if( col2integral )
      col2BoundFlags |= static_cast<int>( ColFlag::kIntegral );
   if( col2lbinf )
      col2BoundFlags |= static_cast<int>( ColFlag::kLbInf );
   if( col2ubinf )
      col2BoundFlags |= static_cast<int>( ColFlag::kUbInf );

   // add all information
   indices.push_back( origcol_mapping[col1] );
   indices.push_back( col1BoundFlags );
   indices.push_back( origcol_mapping[col2] );
   indices.push_back( col2BoundFlags );
   indices.push_back( -1 ); // last index slot is not used
   values.push_back( col1lb );
   values.push_back( col1ub );
   values.push_back( col2lb );
   values.push_back( col2ub );
   values.push_back( col2scale );

   // add the range and the type of the reduction
   types.push_back( ReductionType::kParallelCol );

   finishNotify();
}

template <typename REAL>
PostsolveStatus
Postsolve<REAL>::undo( const Solution<REAL>& reducedSolution,
                       Solution<REAL>& originalSolution ) const
{
   const Vec<REAL>& reducedSol = reducedSolution.primal;
   Vec<REAL>& origSol = originalSolution.primal;

   if( reducedSolution.type == SolutionType::kPrimalDual )
   {
      originalSolution.type = SolutionType::kPrimalDual;
   }

   origSol.clear();
   origSol.resize( nColsOriginal );

   for( int k = 0; k < reducedSol.size(); ++k )
   {
      int origcol = origcol_mapping[k];
      origSol[origcol] = reducedSol[k];
   }

   if( originalSolution.type == SolutionType::kPrimalDual )
   {
      assert( reducedSolution.col_dual.size() == origcol_mapping.size() );
      originalSolution.col_dual.clear();
      originalSolution.col_dual.resize( nColsOriginal );
      for( int k = 0; k < origcol_mapping.size(); k++ )
      {
         int origcol = origcol_mapping[k];
         originalSolution.col_dual[origcol] = reducedSolution.col_dual[k];
      }

      assert( reducedSolution.row_dual.size() == origrow_mapping.size() );
      originalSolution.row_dual.clear();
      originalSolution.row_dual.resize( nRowsOriginal );
      for( int k = 0; k < origrow_mapping.size(); k++ )
      {
         int origrow = origrow_mapping[k];
         originalSolution.row_dual[origrow] = reducedSolution.row_dual[k];
      }
   }

   // If problem has been reduced, check solution of reduced problem returned by
   // solver.
   // At the moment not all row and column changes are notified. The check
   // below handles the case when some trivial presolve elimination is applied,
   // but types.size() is still zero.
   if( ( reducedSolution.row_dual.size() != 0 &&
         nRowsOriginal != reducedSolution.row_dual.size() ) ||
       nColsOriginal != reducedSolution.primal.size() )
   {
      // originalSoluiton is already the reduced solution padded with zeros
      auto kktState =
          checker.initState( ProblemType::kReduced, originalSolution,
                             origcol_mapping, origrow_mapping, checker.level );

      if( originalSolution.type == SolutionType::kPrimalDual )
         checker.level = CheckLevel::Solver_and_primal_feas;
      //   checker.setLevel( CheckLevel::After_each_postsolve_step);

      checker.checkSolution( kktState );
   }

   for( int i = types.size() - 1; i >= 0; --i )
   {
      auto type = types[i];
      int first = start[i];
      int last = start[i + 1];

      switch( type )
      {
      case ReductionType::kSaveRow:
      {
         int row = indices[first];
         int length = (int)values[first];
         bool lb_inf = false;
         bool ub_inf = false;
         if( indices[1] )
            lb_inf = true;
         if( indices[2] )
            ub_inf = true;

         checker.addRowToProblem( row, length, &values[3], &indices[3],
                                  values[1], values[2], lb_inf, ub_inf );
         break;
      }
      case ReductionType::kFixedCol:
      {
         int col = indices[first];
         // todo: move to checker
         // assert( !solSet[col] );
         // solSet[col] = true;
         origSol[col] = values[first];
         break;
      }
      case ReductionType::kSubstitutedCol:
      {
         int col = indices[first];
         // assert( !solSet[col] );
         REAL side = values[first];
         REAL colCoef = 0.0;
         StableSum<REAL> sumcols;
         for( int j = first + 1; j < last; ++j )
         {
            if( indices[j] == col )
               colCoef = values[j];
            else
            {
               // assert( solSet[indices[j]] );
               sumcols.add( origSol[indices[j]] * values[j] );
            }
         }
         sumcols.add( -side );

         assert( colCoef != 0.0 );
         origSol[col] = ( -sumcols.get() ) / colCoef;
         // solSet[col] = true;
         break;
      }
      case ReductionType::kParallelCol:
      {
         constexpr int IS_INTEGRAL = static_cast<int>( ColFlag::kIntegral );
         constexpr int IS_LBINF = static_cast<int>( ColFlag::kLbInf );
         constexpr int IS_UBINF = static_cast<int>( ColFlag::kUbInf );

         assert( last - first == 5 );

         int col1 = indices[first];
         int col1boundFlags = indices[first + 1];
         int col2 = indices[first + 2];
         int col2boundFlags = indices[first + 3];
         const REAL& col1lb = values[first];
         const REAL& col1ub = values[first + 1];
         const REAL& col2lb = values[first + 2];
         const REAL& col2ub = values[first + 3];
         const REAL& col2scale = values[first + 4];
         const REAL& solval = origSol[col2];
         // assert( !solSet[col1] );
         // assert( solSet[col2] );

         // fmt::print( "uncrushing solval {} for parallel cols with scale
         // {}: col1 "
         //            "([{},{}], {}) col2 ([{},{}], {})\n",
         //            solval,
         //            col2scale,
         //            col1boundFlags & IS_LBINF
         //                ? -std::numeric_limits<double>::infinity()
         //                : double( col1lb ),
         //            col1boundFlags& IS_UBINF
         //                ? std::numeric_limits<double>::infinity()
         //                : double( col1ub ),
         //            col1boundFlags& IS_INTEGRAL ? "int." : "cont.",
         //            col2boundFlags& IS_LBINF
         //                ? -std::numeric_limits<double>::infinity()
         //                : double( col2lb ),
         //            col2boundFlags& IS_UBINF
         //                ? std::numeric_limits<double>::infinity()
         //                : double( col2ub ),
         //            col2boundFlags& IS_INTEGRAL ? "int." : "cont." );

         REAL col1val;
         REAL col2val;

         if( col1boundFlags & IS_INTEGRAL )
         {
            assert( !( col1boundFlags & IS_LBINF ) );
            assert( !( col1boundFlags & IS_UBINF ) );
            assert( !( col2boundFlags & IS_LBINF ) );
            assert( !( col2boundFlags & IS_UBINF ) );
            assert( col2boundFlags & IS_INTEGRAL );

            col1val = col1lb;

            while( num.isFeasLE( col1val, col1ub ) )
            {
               col2val = solval - col1val * col2scale;

               if( num.isFeasIntegral( col2val ) &&
                   num.isFeasGE( col2val, col2lb ) &&
                   num.isFeasLE( col2val, col2ub ) )
                  break;

               col1val += 1;
            }
         }
         else
         {
            REAL col2valGuess;
            if( !( col2boundFlags & IS_LBINF ) )
               col2valGuess = col2lb;
            else if( !( col2boundFlags & IS_UBINF ) )
               col2valGuess = col2ub;
            else
               col2valGuess = 0;

            col1val = ( solval - col2valGuess ) / col2scale;

            // check if value for column 1 is feasible
            if( !( col1boundFlags & IS_LBINF ) &&
                num.isFeasLT( col1val, col1lb ) )
            {
               // lower bound was violated -> set column 1 to lower bound
               col1val = col1lb;
               // compute new value for column1
               col2val = solval - col2scale * col1val;
            }
            else if( !( col1boundFlags & IS_UBINF ) &&
                     num.isFeasGT( col1val, col1ub ) )
            {
               // upper bound was violated -> set column 1 to lower bound
               col1val = col1ub;
               // compute new value for column 2
               col2val = solval - col2scale * col1val;
            }
            else
            {
               // guess was feasible
               col2val = col2valGuess;
            }

            // domains should be valid now,  except for integrality of column
            // 2 which could still be violated
            assert( ( col1boundFlags & IS_LBINF ) ||
                    num.isFeasGE( col1val, col1lb ) );
            assert( ( col1boundFlags & IS_UBINF ) ||
                    num.isFeasLE( col1val, col1ub ) );
            assert( ( col2boundFlags & IS_LBINF ) ||
                    num.isFeasGE( col2val, col2lb ) );
            assert( ( col2boundFlags & IS_UBINF ) ||
                    num.isFeasLE( col2val, col2ub ) );

            // maybe integrality is violated now for column 2, then we round
            // further in the direction that we already moved to column 2 to
            if( ( col2boundFlags & IS_INTEGRAL ) &&
                !num.isFeasIntegral( col2val ) )
            {
               // round in the direction that we moved away from when we
               // corrected the bound violation for column 1 otherwise we
               // will violate the bounds of column 1 again
               if( col2val > col2valGuess )
                  col2val = ceil( col2val );
               else
                  col2val = floor( col2val );

               // recompute value for column 1
               col1val = solval - col1val * col2scale;
            }
         }
         // bounds and integrality should hold now within feasibility
         // tolerance
         assert( ( col1boundFlags & IS_LBINF ) ||
                 num.isFeasGE( col1val, col1lb ) );
         assert( ( col1boundFlags & IS_UBINF ) ||
                 num.isFeasLE( col1val, col1ub ) );
         assert( !( col1boundFlags & IS_INTEGRAL ) ||
                 num.isFeasIntegral( col1val ) );
         assert( ( col2boundFlags & IS_LBINF ) ||
                 num.isFeasGE( col2val, col2lb ) );
         assert( ( col2boundFlags & IS_UBINF ) ||
                 num.isFeasLE( col2val, col2ub ) );
         assert( !( col2boundFlags & IS_INTEGRAL ) ||
                 num.isFeasIntegral( col2val ) );
         assert( num.isFeasEq( solval, col2scale * col1val + col2val ) );

         // solSet[col1] = true;
         origSol[col1] = col1val;
         origSol[col2] = col2val;

         break;
      }
      }

      // intermediate kkt check
      if( checker.level == CheckLevel::After_each_postsolve_step )
      {
         auto kktStatePostsolvedProblem = checker.initState(
             ProblemType::kPostsolved, originalSolution, checker.level );
         checker.checkIntermediate( kktStatePostsolvedProblem );
      }
   }

   // todo: move to checker
   // assert( std::all_of( solSet.begin(), solSet.end(),
   //                      []( uint8_t isset ) { return isset; } ) );

   auto kktStateOriginalProblem = checker.initState(
       ProblemType::kOriginal, originalSolution, checker.level );
   checker.checkSolution( kktStateOriginalProblem );

   return PostsolveStatus::kOk;
}

} // namespace papilo

#endif
