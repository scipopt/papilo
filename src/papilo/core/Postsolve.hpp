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
   kFixedCol,
   kSubstitutedCol,
   kParallelCol,
   kSingletonRow,
   kRedundantRow,
   kDeletedCol,
   kBoundChange,
   kColumnDualValue,
   kRowDualValue,
   kReducedBoundsCost,
   kSaveRow,
   kSaveCol
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

   Postsolve( const int nrows, const int ncols )
   {
      origrow_mapping.reserve( nrows );
      origrow_mapping.reserve( ncols );

      for( int i = 0; i < nrows; ++i )
         origrow_mapping.push_back( i );

      for( int i = 0; i < ncols; ++i )
         origcol_mapping.push_back( i );

      nColsOriginal = ncols;
      nRowsOriginal = nrows;

      start.push_back( 0 );
   }

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
   notifyRedundantRow( const int row );

   void
   notifyDeletedCol( const int col );

   void
   notifyBoundChange( const bool is_row, const bool is_lower, const int col,
                      const int row, const REAL old_bound,
                      const REAL new_bound );

   void
   notifyReducedBoundsAndCost( const Vec<REAL>& col_lb, const Vec<REAL>& col_ub,
                               const Vec<REAL>& row_lb, const Vec<REAL>& row_ub,
                               const Vec<REAL>& cost,
                               const Vec<RowFlags>& row_flags,
                               const Vec<ColFlags>& col_flags );

   // todo: modify with colvec and col cost so if dual postsolve
   // col values are added so we can get dual value
   void
   notifyFixedCol( const int col, const REAL val,
                   const SparseVectorView<REAL>& colvec,
                   const Vec<REAL>& cost );

   void
   notifySingletonRow( const int row, const int col, const REAL coeff,
                       const Vec<REAL>& cost,
                       const SparseVectorView<REAL>& colvec );

   void
   notifyDualValue( bool is_column_dual, int index, REAL value );

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

   // TODO add mechanism for saving columns as well
   Vec<int> row_stack_index;
};

#ifdef PAPILO_USE_EXTERN_TEMPLATES
extern template class Postsolve<double>;
extern template class Postsolve<Quad>;
extern template class Postsolve<Rational>;
#endif

template <typename REAL>
void
Postsolve<REAL>::notifyRedundantRow( const int row )
{
   // TODO2: actually this should not require to have the row stored on
   // postsolve at all but you might need the row for the checker. Therefore
   // the row should only be stored in the checker. To make this easier
   // the checker should just store a reference to the problem
   // that is passed in during construction of the postsolve, which wil always
   // correspond to the current reduced problem. Then we do not need to pass
   // extra information to the checker and it can store the additional
   // information by itself.

   // Apart from that the value should be able to store a
   // dual value so that this function can also be used when components are
   // solved individually. For that case the postsolve will just restore exactly
   // that dual value and not care about reduced costs as they would be fixed
   // separately to the values returned by the solver for one component.

   types.push_back( ReductionType::kRedundantRow );
   indices.push_back(
       origrow_mapping[row] ); // TODO: this was an error! this must map to
                               // original space, check at other places too! I
                               // added it for here
   values.push_back( 0 );

   finishNotify();
}

template <typename REAL>
void
Postsolve<REAL>::notifyDeletedCol( const int col )
{
   // TODO I think we do not need notifyDeletedCol. A column is deleted when it
   // is fixed or substituted or a parallel column. But all those have their own
   // postsolve notify function that must handle all necessary information.
   types.push_back( ReductionType::kDeletedCol );
   indices.push_back( col );
   values.push_back( 0 );

   finishNotify();
}

template <typename REAL>
void
Postsolve<REAL>::notifyBoundChange( const bool is_row, const bool is_lower,
                                    const int row, const int col,
                                    const REAL old_bound, const REAL new_bound )
{
   // TODO, this is not needed due to the bound relaxing strategy I'll
   // add for constraint propagation, instead there should only be a function
   // notifyForcingRow. This is called for the case where a row forces a column
   // upper bound to its lower bound and the column is fixed as a result, or the
   // other way around.
   types.push_back( ReductionType::kBoundChange );
   if( is_row && is_lower )
      indices.push_back( 0 );
   else if( is_row )
      indices.push_back( 1 );
   else if( is_lower )
      indices.push_back( 2 );
   else
      indices.push_back( 3 );
   values.push_back( 0 );

   indices.push_back( row );
   indices.push_back( col );
   values.push_back( old_bound );
   values.push_back( new_bound );

   finishNotify();
}

template <typename REAL>
void
Postsolve<REAL>::notifyReducedBoundsAndCost(
    const Vec<REAL>& col_lb, const Vec<REAL>& col_ub, const Vec<REAL>& row_lb,
    const Vec<REAL>& row_ub, const Vec<REAL>& cost,
    const Vec<RowFlags>& row_flags, const Vec<ColFlags>& col_flags )
{
   // TODO for what is this notification required? Can you add comments?
   // the postsolve stores the original problem. The notify functions are not
   // for the checker, the checker must get around without notifies and is only
   // informed about changes from within postsolve notify functions. As
   // mentioned in the above comment, the checker can store a reference that
   // always contains the current reduced problem

   types.push_back( ReductionType::kReducedBoundsCost );

   // would be better to only pass finite values, not all
   // col bounds
   for( int col = 0; col < col_lb.size(); col++ )
   {
      int flag_lb = 0;
      int flag_ub = 0;
      if( col_flags[col].test( ColFlag::LB_INF ) )
         flag_lb |= static_cast<int>( ColFlag::LB_INF );
      if( col_flags[col].test( ColFlag::UB_INF ) )
         flag_ub |= static_cast<int>( ColFlag::UB_INF );
      indices.push_back( flag_lb );
      values.push_back( col_lb[col] );
      indices.push_back( flag_ub );
      values.push_back( col_ub[col] );
   }

   // row bounds
   for( int row = 0; row < row_lb.size(); row++ )
   {
      int flag_lb = 0;
      int flag_ub = 0;
      if( row_flags[row].test( RowFlag::LHS_INF ) )
         flag_lb |= static_cast<int>( RowFlag::LHS_INF );
      if( row_flags[row].test( RowFlag::RHS_INF ) )
         flag_ub |= static_cast<int>( RowFlag::RHS_INF );
      indices.push_back( flag_lb );
      values.push_back( row_lb[row] );
      indices.push_back( flag_ub );
      values.push_back( row_ub[row] );
   }

   // col cost
   for( int col = 0; col < cost.size(); col++ )
   {
      indices.push_back( col );
      values.push_back( cost[col] );
   }

   finishNotify();
}

template <typename REAL>
void
Postsolve<REAL>::notifyFixedCol( int col, const REAL val,
                                 const SparseVectorView<REAL>& colvec,
                                 const Vec<REAL>& cost )
{
   types.push_back( ReductionType::kFixedCol );
   indices.push_back( origcol_mapping[col] );
   values.push_back( val );

   if( postsolveType == PostsolveType::FULL )
   {
      // TODO this should probably use the saveCol mechanism if the column
      // values are needed
      const int length = colvec.getLength();
      indices.push_back( length );
      values.push_back( cost[origcol_mapping[col]] );

      const REAL* vals = colvec.getValues();
      const int* inds = colvec.getIndices();

      for( int j = 0; j < length; j++ )
      {
         indices.push_back( inds[j] );
         values.push_back( vals[j] );
      }
   }

   finishNotify();
}

template <typename REAL>
void
Postsolve<REAL>::notifySingletonRow( const int row, const int col,
                                     const REAL coeff, const Vec<REAL>& cost,
                                     const SparseVectorView<REAL>& colvec )
{
   // TODO: in general it would be good for all the notify functions, now that
   // there are so many to have comments that say which information needs to be
   // stored for which postsolve type
   types.push_back( ReductionType::kSingletonRow );
   indices.push_back( origrow_mapping[row] );
   values.push_back( 0 );
   indices.push_back( origcol_mapping[col] );
   values.push_back( coeff );

   const int length = colvec.getLength();
   indices.push_back( length - 1 );
   values.push_back( cost[origcol_mapping[col]] );

   const REAL* vals = colvec.getValues();
   const int* inds = colvec.getIndices();

   for( int j = 0; j < length; j++ )
   {
      if( inds[j] != row )
      {
         indices.push_back( inds[j] );
         values.push_back( vals[j] );
      }
   }

   finishNotify();
}

template <typename REAL>
void
Postsolve<REAL>::notifyDualValue( bool is_column_dual, int index, REAL value )
{
   // TODO, for which reduction is this notify function for?
   // Pushing zero so I don't modity finishNotify()'s assert (for the moment)
   if( is_column_dual )
      types.push_back( ReductionType::kColumnDualValue );
   else
      types.push_back( ReductionType::kRowDualValue );

   indices.push_back( index );
   values.push_back( value );
   finishNotify();
}

// TODO remove dead code if not needed anymore
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
   // TODO: depending on the postsolve type I guess we need to save also the
   //       column, but for this branch lets focus on a working dual postsolve
   //       only for the trivial presolve and lets add tests for that. Simple
   //       tests could basically just read the MIP  instances in test/instances
   //       folder, then discard integrality information, and apply a presolve
   //       procedure that only uses trivial presolve. I think when this is done
   //       I will add a flag to the presolvers for which postsolve type they
   //       are compatible and default all to only primal. Then we can work on
   //       adding more and more presolvers in later merge requests.
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
   if( origrow_mapping.size() < nRowsOriginal ||
       origcol_mapping.size() < nColsOriginal )
   {
      CheckLevel level = CheckLevel::Primal_only;
      if( originalSolution.type == SolutionType::PRIMAL_AND_DUAL )
         CheckLevel level = CheckLevel::Primal_and_dual;

      auto kktState =
          checker.initState( ProblemType::kReduced, originalSolution,
                             origcol_mapping, origrow_mapping, level );

      checker.checkSolution( kktState );
      checker.level = CheckLevel::After_each_step_primal_only;

      checker.expandProblem();
      auto kktState_expand = checker.initState(
          ProblemType::POSTSOLVED, originalSolution, checker.level );
      checker.checkSolution( kktState_expand );
   }

   // Will be used during dual postsolve for fast access to bound values.
   Vec<REAL> col_cost;
   Vec<REAL> col_lower;
   Vec<REAL> col_upper;
   Vec<REAL> row_lower;
   Vec<REAL> row_upper;

   Vec<int> col_bound_lower;
   Vec<int> col_bound_upper;
   Vec<int> row_bound_lower;
   Vec<int> row_bound_upper;

   Vec<int> col_lower_from_row;
   Vec<int> col_upper_from_row;
   Vec<int> row_lower_from_col;
   Vec<int> row_upper_from_col;

   if( originalSolution.type == SolutionType::PRIMAL_AND_DUAL )
   {

      col_cost.assign( nColsOriginal, 0 );
      col_lower.assign( nColsOriginal, 0 );
      col_upper.assign( nColsOriginal, 0 );
      row_lower.assign( nRowsOriginal, 0 );
      row_upper.assign( nRowsOriginal, 0 );
      col_bound_upper.assign( nColsOriginal, 0 );
      col_bound_lower.assign( nColsOriginal, 0 );
      row_bound_upper.assign( nRowsOriginal, 0 );
      row_bound_lower.assign( nRowsOriginal, 0 );

      // // expand vectors.
      // Vec<REAL> tmp_col_cost = col_cost;
      // Vec<REAL> tmp_col_lower = col_lower;
      // Vec<REAL> tmp_col_upper = col_upper;

      // col_cost.clear();
      // col_lower.clear();
      // col_upper.clear();
      // col_cost.resize( nColsOriginal );
      // col_lower.resize( nColsOriginal );
      // col_upper.resize( nColsOriginal );

      // for( int k = 0; k < origcol_mapping.size(); k++ )
      // {
      //    int origcol = origcol_mapping[k];
      //    col_cost[origcol] = tmp_col_cost[k];
      //    col_lower[origcol] = tmp_col_lower[k];
      //    col_upper[origcol] = tmp_col_upper[k];
      // }

      // Vec<RowActivity<REAL>> tmp_activity = activity;
      // Vec<RowFlags> tmp_flags = flags;

      // // what are flags and activities initialized to?
      // Vec<RowActivity<REAL>> expanded_activity;
      // Vec<RowFlags> expanded_flags;

      // int current = 0;
      // for( int k = 0; k < origrow_mapping.size(); k++ )
      // {
      //    int origrow = origrow_mapping[k];
      //    while( origrow > current )
      //    {
      //       RowActivity<REAL> activity;
      //       expanded_activity.push_back( activity );
      //       RowFlags flags;
      //       // todo: set to infinite sce rows have been removed so are
      //       // redundant.
      //       expanded_flags.push_back( flags );
      //       current++;
      //    }
      //    flags[origrow] = tmp_flags[k];
      //    activity[origrow] = tmp_activity[k];
      //    current++;
      // }
      // flags = expanded_flags;
      // activity = expanded_activity;

      col_lower_from_row.assign( nColsOriginal, -1 );
      col_upper_from_row.assign( nColsOriginal, -1 );
      row_lower_from_col.assign( nRowsOriginal, -1 );
      row_upper_from_col.assign( nRowsOriginal, -1 );
   }

   for( int i = types.size() - 1; i >= 0; --i )
   {
      auto type = types[i];
      int first = start[i];
      int last = start[i + 1];

      switch( type )
      {
      case ReductionType::kReducedBoundsCost:
      {
         // get column bounds
         for( int j = 0; j < origcol_mapping.size(); j++ )
         {
            int origcol = origcol_mapping[j];
            int index = first + 2 * j;
            col_lower[origcol] = values[index];
            col_upper[origcol] = values[index + 1];
            col_bound_lower[origcol] = indices[index];
            col_bound_upper[origcol] = indices[index + 1];
         }

         // get row bounds
         int first_row_bounds = first + 2 * origcol_mapping.size();
         for( int i = 0; i < origrow_mapping.size(); i++ )
         {
            int origrow = origrow_mapping[i];
            int index = first_row_bounds + 2 * i;
            row_lower[origrow] = values[index];
            row_upper[origrow] = values[index + 1];
            row_bound_lower[origrow] = indices[index];
            row_bound_upper[origrow] = indices[index + 1];
         }

         // get cost
         int first_cost = first_row_bounds + 2 * origrow_mapping.size();
         for( int j = 0; j < origcol_mapping.size(); j++ )
         {
            int origcol = origcol_mapping[j];
            col_cost[origcol] = values[first_cost + j];
            assert( j == indices[first_cost + j] );
         }
         break;
      }
      case ReductionType::kColumnDualValue:
         originalSolution.col_dual[indices[first]] = values[indices[first]];
         break;
      case ReductionType::kRedundantRow:
         checker.undoRedundantRow( indices[first] );
         break;
      case ReductionType::kDeletedCol:
         checker.undoDeletedCol( indices[first] );
         break;
      case ReductionType::kRowDualValue:
         originalSolution.row_dual[indices[first]] = values[indices[first]];
         break;
      case ReductionType::kSaveRow:
      {
         // TODO I think this and the SAVE_COL step should just be skipped
         //     we only want to restore redundant rows that have been removed
         //     not all saved rows. Saved row has no logical implications for
         //     postsolve it should just give the deleted row a new index within
         //     the postsolve structure this index would be used by the undo
         //     redundant col I guess.
         int row = indices[first];
         int length = (int)values[first];
         bool lb_inf = false;
         bool ub_inf = false;
         if( indices[first + 1] )
            lb_inf = true;
         if( indices[first + 2] )
            ub_inf = true;

         checker.addRowToProblem( row, length, &values[first + 3],
                                  &indices[first + 3], values[first + 1],
                                  values[first + 2], lb_inf, ub_inf );
         break;
      }
      case ReductionType::SAVE_COL:
      {
         break;
      }
      case ReductionType::kFixedCol:
      {
         // At the moment saves column to the postsolve stack. todo:
         // use notifySavedCol and current index of column on the stack.
         int col = indices[first];
         origSol[col] = values[first];
         checker.undoFixedCol( col, values[first] );
         // todo: checker notify dual value if changed
         if( originalSolution.type == SolutionType::kPrimalDual )
         {
            // get dual value z_j = c_j - sum_i a_ij*y_i
            REAL value = 0;

            REAL objective_coefficient = values[first + 1];
            int col_length = indices[first + 1];

            // no need to check for solSetRow because if it iz zero then the
            // row_dual value is zero.
            for( int i = 0; i < col_length; ++i )
            {
               int index = first + 1 + i;
               value = value - values[index] *
                                   originalSolution.row_dual[indices[index]];
            }
            value = value + objective_coefficient;
            originalSolution.col_dual[col] = value;
         }
         break;
      }
      case ReductionType::kSingletonRow:
      {
         int row = indices[first];
         int col = indices[first + 1];
         if( originalSolution.type == SolutionType::kPrimalDual )
         // code below saves column on stack for the calculation of dual values.
         // todo: use column on stack
         {
            REAL coeff = values[first + 1];
            REAL cost = values[first + 2];
            assert( indices[first + 1] == col );

            REAL value = 0;
            int col_length_minus_one = indices[first + 2];

            // no need to check for solSetRow because if it iz zero then the
            // row_dual value is zero.
            for( int i = 0; i < col_length_minus_one; ++i )
            {
               int index = first + 3 + i;
               value = value + values[index] *
                                   originalSolution.row_dual[indices[index]];
            }

            value = cost - value;
            value = value / coeff;

            originalSolution.row_dual[row] = value;
            originalSolution.col_dual[col] = 0;
         }
         checker.undoSingletonRow( row );
         break;
      }
      case ReductionType::kBoundChange:
      {
         // -1 column cost
         // 0 primal row lower
         // 1 primal row upper
         // 2 primal col lower
         // 3 primal col upper
         // .. dual strong
         // .. dual weak
         int type = indices[first];
         int row = indices[first + 1];
         int col = indices[first + 2];
         REAL old_value = values[first + 1];
         REAL new_value = values[first + 2];

         // todo: also set flags for rows and columns which have a bound now
         switch( type )
         {
         case -1:
            col_cost[col] = old_value;
            break;
         case 0:
            row_lower[row] = old_value;
            row_lower_from_col[row] = col;
            break;
         case 1:
            row_upper[row] = old_value;
            row_upper_from_col[row] = col;
            break;
         case 2:
            col_lower[col] = old_value;
            col_lower_from_row[col] = row;
            break;
         case 3:
            col_upper[col] = old_value;
            col_upper_from_row[col] = row;
            break;
         }

         break;
      }
         // todo: modify, unused right now
         // case ReductionType::kRedundantRow:
         //{
         // if( originalSolution.type == SolutionType::PRIMAL_AND_DUAL )
         // {
         //    int row = indices[first];
         //    int col_lo = indices[first + 1];
         //    int col_up = indices[first + 2];
         //    // todo: it happens sometimes like deteq8, debug.
         //    assert( ( col_lo == -1 ) != ( col_up == -1 ) );

         //    REAL coeff = 0;
         //    int col = 0;
         //    if( col_lo != -1 )
         //    {
         //       col = col_lo;
         //       coeff = values[first + 1];
         //    }
         //    else
         //    {
         //       col = col_up;
         //       coeff = values[first + 2];
         //    }

         //    originalSolution.row_dual[row] =
         //        originalSolution.row_dual[row] +
         //        originalSolution.col_dual[col] / coeff;
         //    originalSolution.col_dual[col] = 0;

         //    solSetRow[row] = true;
         // }
         //         break;
         //     }
      case ReductionType::kSubstitutedCol:
      {
         int col = indices[first];
         checker.undoSubstitutedCol( col );

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

         checker.undoParallelCol( col1 );
         // solSet[col1] = true;
         origSol[col1] = col1val;
         origSol[col2] = col2val;

         break;
      }
      }

      // intermediate kkt check
      if( checker.level == CheckLevel::After_each_step_primal_only ||
          checker.level == CheckLevel::After_each_step_and_dual )
      {
         auto kktStatePostsolvedProblem = checker.initState(

             ProblemType::kPostsolved, originalSolution, checker.level );
         checker.checkSolution( kktStatePostsolvedProblem, true );
      }
   }

   auto kktStateOriginalProblem = checker.initState(
       ProblemType::kOriginal, originalSolution, checker.level );
   checker.checkSolution( kktStateOriginalProblem );

   return PostsolveStatus::kOk;
}

} // namespace papilo

#endif
