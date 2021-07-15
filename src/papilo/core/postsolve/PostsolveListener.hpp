/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*               This file is part of the program and library                */
/*    PaPILO --- Parallel Presolve for Integer and Linear Optimization       */
/*                                                                           */
/* Copyright (C) 2020-2021 Konrad-Zuse-Zentrum                               */
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

#ifndef _PAPILO_CORE_POSTSOLVE_LISTENER_HPP_
#define _PAPILO_CORE_POSTSOLVE_LISTENER_HPP_

#include "papilo/core/Problem.hpp"
#include "papilo/misc/MultiPrecision.hpp"
#include "papilo/misc/Num.hpp"
#include "papilo/misc/StableSum.hpp"
#include "papilo/misc/Vec.hpp"
#include "papilo/misc/compress_vector.hpp"
#include "papilo/misc/fmt.hpp"
#include "papilo/misc/tbb.hpp"
#include <fstream>
#include "papilo/misc/dualpostsolve/PrimalDualSolValidation.hpp"
#include "papilo/core/postsolve/PostsolveType.hpp"
#include "papilo/core/postsolve/ReductionType.hpp"

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/tmpdir.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/list.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>

namespace papilo
{

template <typename REAL>
class SparseVectorView;

struct IndexRange;

/// type to store necessary data for post solve
template <typename REAL>
class PostsolveListener
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

   PostsolveListener() = default;

   PostsolveListener( const int nrows, const int ncols )
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

   PostsolveListener( const Problem<REAL>& problem, const Num<REAL>& num )
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

   void
   notifyRedundantRow( const int row );

   void
   notifyDeletedCol( const int col );

   void
   notifyVarBoundChange( const bool isLowerBound,
                         const int row, const int col, const REAL oldBound,
                         bool isInfinity, const REAL newBound )
   ;

   void
   notifyRowBoundChange( const bool isLhs,
                         const int row, const int col,
                         const REAL oldBound, const REAL newBound );


   void
   notifyReducedBoundsAndCost( const Vec<REAL>& col_lb, const Vec<REAL>& col_ub,
                               const Vec<REAL>& row_lhs, const Vec<REAL>& row_rhs,
                               const Vec<REAL>& coefficients,
                               const Vec<RowFlags>& row_flags,
                               const Vec<ColFlags>& col_flags );

   // todo: modify with colvec and col cost so if dual postsolve
   // col values are added so we can get dual value
   void
   notifyFixedCol( const int col, const REAL val,
                   const SparseVectorView<REAL>& colvec,
                   const Vec<REAL>& cost );

   void
   notifySingletonRow( const int row, const int col,
                       const REAL coeff, const Vec<REAL>& cost,
                       const SparseVectorView<REAL>& colvec,
                       const REAL lhs, bool isLhsInfinity,
                       const REAL rhs, bool isRhsInfinity )
   ;

   void
   notifyDualValue( bool is_column_dual, int index, REAL value );

   void
   notifyFixedInfCol( int col, REAL val, REAL bound,
                      const Problem<REAL>& currentProblem );

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
   void
   push_back_row( int row, const Problem<REAL>& currentProblem );

};

#ifdef PAPILO_USE_EXTERN_TEMPLATES
extern template class PostsolveListener<double>;
extern template class PostsolveListener<Quad>;
extern template class PostsolveListener<Rational>;
#endif

template <typename REAL>
void
PostsolveListener<REAL>::notifyRedundantRow( const int row )
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
PostsolveListener<REAL>::notifyDeletedCol( const int col )
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
PostsolveListener<REAL>::notifyVarBoundChange( const bool isLowerBound,
                                       const int row, const int col, const REAL oldBound,
                                       bool isInfinity, const REAL newBound )
{
   // TODO, this is not needed due to the bound relaxing strategy I'll
   // add for constraint propagation, instead there should only be a function
   // notifyForcingRow. This is called for the case where a row forces a column
   // upper bound to its lower bound and the column is fixed as a result, or the
   // other way around.
   types.push_back( ReductionType::kVarBoundChange );
   if( isLowerBound )
      indices.push_back( 1 );
   else
      indices.push_back( 0 );
   values.push_back( 0 );

   indices.push_back( col );
   values.push_back( newBound );

   indices.push_back( isInfinity );
   values.push_back( oldBound );

   finishNotify();
}

template <typename REAL>
void
PostsolveListener<REAL>::notifyRowBoundChange( const bool isLhs,
                                       const int row, const int col,
                                       const REAL oldBound, const REAL newBound )
{
   // TODO, this is not needed due to the bound relaxing strategy I'll
   // add for constraint propagation, instead there should only be a function
   // notifyForcingRow. This is called for the case where a row forces a column
   // upper bound to its lower bound and the column is fixed as a result, or the
   // other way around.
   types.push_back( ReductionType::kRowBoundChange );
   if( isLhs )
      indices.push_back( 1 );
   else
      indices.push_back( 0 );
   values.push_back( 0 );
   indices.push_back( 0 );
   indices.push_back( col );
   values.push_back( oldBound );
   values.push_back( newBound );

   finishNotify();
}

template <typename REAL>
void
PostsolveListener<REAL>::notifyReducedBoundsAndCost(
    const Vec<REAL>& col_lb, const Vec<REAL>& col_ub, const Vec<REAL>& row_lhs,
    const Vec<REAL>& row_rhs, const Vec<REAL>& coefficients,
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
      if( col_flags[col].test( ColFlag::kLbInf ) )
         flag_lb |= static_cast<int>( ColFlag::kLbInf );
      if( col_flags[col].test( ColFlag::kUbInf ) )
         flag_ub |= static_cast<int>( ColFlag::kUbInf );
      indices.push_back( flag_lb );
      values.push_back( col_lb[col] );
      indices.push_back( flag_ub );
      values.push_back( col_ub[col] );
   }

   // row bounds
   for( int row = 0; row < row_lhs.size(); row++ )
   {
      int flag_lb = 0;
      int flag_ub = 0;
      if( row_flags[row].test( RowFlag::kLhsInf ) )
         flag_lb |= static_cast<int>( RowFlag::kLhsInf );
      if( row_flags[row].test( RowFlag::kRhsInf ) )
         flag_ub |= static_cast<int>( RowFlag::kRhsInf );
      indices.push_back( flag_lb );
      values.push_back( row_lhs[row] );
      indices.push_back( flag_ub );
      values.push_back( row_rhs[row] );
   }

   // col coefficients
   for( int col = 0; col < coefficients.size(); col++ )
   {
      indices.push_back( col );
      values.push_back( coefficients[col] );
   }

   finishNotify();
}

template <typename REAL>
void
PostsolveListener<REAL>::notifyFixedCol( int col, const REAL val,
                                 const SparseVectorView<REAL>& colvec,
                                 const Vec<REAL>& cost )
{
   types.push_back( ReductionType::kFixedCol );
   indices.push_back( origcol_mapping[col] );
   values.push_back( val );

   if( postsolveType == PostsolveType::kFull )
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

/**
 * If a singleton row is deleted aka converted to an lower or upper bound,
 * this function saves the information to recalculate the dual solution.
 * This function should only be called if the singletonRow implies an tighter
 * bound. (TO BE CHECKED)
 * In this case the c^T - y^T *A needs to be zero because the original bound is not going to be met.
 * Therefore save the current column vector to recalculate.
 * (c^T-(y^T*A\{col}))/a_col.
 * @tparam REAL
 * @param row row index of singleton row
 * @param col column index of singleton row
 * @param coeff a_col
 * @param cost obj_col = c^T
 * @param colvec A\{col}
 */
template <typename REAL>
void
PostsolveListener<REAL>::notifySingletonRow( const int row, const int col,
                                     const REAL coeff, const Vec<REAL>& cost,
                                     const SparseVectorView<REAL>& colvec,
                                     const REAL lhs, bool isLhsInfinity,
                                     const REAL rhs, bool isRhsInfinity )
{
   types.push_back( ReductionType::kSingletonRow );
   indices.push_back( origrow_mapping[row] );
   values.push_back( 0 );
   indices.push_back( origcol_mapping[col] );
   values.push_back( coeff );

   indices.push_back( isLhsInfinity );
   values.push_back( lhs );
   indices.push_back( isRhsInfinity );
   values.push_back( rhs );

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
PostsolveListener<REAL>::notifyDualValue( bool is_column_dual, int index, REAL value )
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
// PostsolveListener<REAL>::notifySavedRow( int row,
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
PostsolveListener<REAL>::push_back_row( int row, const Problem<REAL>& currentProblem )
{
   const auto& coefficients =
       currentProblem.getConstraintMatrix().getRowCoefficients( row );
   REAL lhs = currentProblem.getConstraintMatrix().getLeftHandSides()[row];
   REAL rhs = currentProblem.getConstraintMatrix().getRightHandSides()[row];
   const auto& flags = currentProblem.getConstraintMatrix().getRowFlags()[row];

   const REAL* coefs = coefficients.getValues();
   const int* columns = coefficients.getIndices();
   const int length = coefficients.getLength();

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
      indices.push_back( origcol_mapping[columns[i]] );
      values.push_back( coefs[i] );
   }
}

template <typename REAL>
void
PostsolveListener<REAL>::notifyFixedInfCol( int col, REAL val, REAL bound,
                                    const Problem<REAL>& currentProblem )
{
   types.push_back( ReductionType::kFixedInfCol );
   indices.push_back( origcol_mapping[col] );
   values.push_back( val );
   indices.push_back( 0 );
   values.push_back( bound );

   const auto& coefficients =
       currentProblem.getConstraintMatrix().getColumnCoefficients( col );
   const int* row_indices = coefficients.getIndices();

   for( int i = 0; i < coefficients.getLength(); i++ )
      push_back_row( row_indices[i], currentProblem );

   finishNotify();
}

template <typename REAL>
void
PostsolveListener<REAL>::notifySubstitution( int col,
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
PostsolveListener<REAL>::notifyParallelCols( int col1, bool col1integral,
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




} // namespace papilo

#endif
