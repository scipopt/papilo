/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*               This file is part of the program and library                */
/*    PaPILO --- Parallel Presolve for Integer and Linear Optimization       */
/*                                                                           */
/* Copyright (C) 2020-2024 Zuse Institute Berlin (ZIB)                       */
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

#ifndef _PAPILO_PRESOLVERS_CliqueMerging_HPP_
#define _PAPILO_PRESOLVERS_CliqueMerging_HPP_

#include "papilo/core/PresolveMethod.hpp"
#include "papilo/core/Problem.hpp"
#include "papilo/core/ProblemUpdate.hpp"
#include <vector>
#include <algorithm>
#ifdef PAPILO_TBB
#include "papilo/misc/tbb.hpp"
#endif

namespace papilo
{

template <typename REAL>
class CliqueMerging : public PresolveMethod<REAL>
{
 public:
   CliqueMerging() : PresolveMethod<REAL>()
   {
      this->setName( "cliquemerging" );
      this->setTiming( PresolverTiming::kExhaustive );
      this->setType( PresolverType::kIntegralCols );
   }

   PresolveStatus
   execute( const Problem<REAL>& problem,
            const ProblemUpdate<REAL>& problemUpdate, const Num<REAL>& num,
            Reductions<REAL>& reductions, const Timer& timer,
            int& reason_of_infeasibility ) override;

   Vec<int>
   greedyClique( const ConstraintMatrix<REAL>& matrix,
               const int cliqueRow );

   int
   getNeighbour( const ConstraintMatrix<REAL>& matrix,
               int col, int neighbournumber );

   int
   getNeighbourhoodSize( const ConstraintMatrix<REAL>& matrix, 
                        int col);

   bool
   isNeighbour( const ConstraintMatrix<REAL>& matrix, 
               int u, int v );

   bool
   isCovered( const ConstraintMatrix<REAL>& matrix, 
            const int row, const Vec<int> newClique, 
            const int clique);

};

#ifdef PAPILO_USE_EXTERN_TEMPLATES
extern template class CliqueMerging<double>;
extern template class CliqueMerging<Quad>;
extern template class CliqueMerging<Rational>;
#endif

template <typename REAL>
int
CliqueMerging<REAL>::getNeighbour( const ConstraintMatrix<REAL>& matrix,
                                    int col, int neighbournumber )
{
   const auto colvec = matrix.getColumnCoefficients( col );
   const int* colinds = colvec.getIndices();
   const std::vector<RowFlags> rowFlags = matrix.getRowFlags();
   int neighbour = -1;
   int neighboursskipped = 0;
   for( int i = 0; i < colvec.getLength(); ++i )
   {
      if( !rowFlags[colinds[i]].test(RowFlag::kClique) )
         continue;
      auto rowvec = matrix.getRowCoefficients( colinds[i] );
      if( neighbournumber - neighboursskipped >= rowvec.getLength() - 1 )
      {
         neighboursskipped += rowvec.getLength() - 1;
         continue;
      }
      else
      {
         const auto rowinds = rowvec.getIndices();
         if( rowinds[neighbournumber - neighboursskipped] >= col )
            neighbour = rowinds[neighbournumber - neighboursskipped + 1];
         else
            neighbour = rowinds[neighbournumber - neighboursskipped];
         break;
      }
   }
   assert( neighbour >= 0);
   return neighbour ;
}

template <typename REAL>
int
CliqueMerging<REAL>::getNeighbourhoodSize( const ConstraintMatrix<REAL>& matrix, 
                                    int col)
{
   const auto colvec = matrix.getColumnCoefficients( col );
   const int* inds = colvec.getIndices();
   const std::vector<RowFlags> rowFlags = matrix.getRowFlags();
   int neighbournumber = 0;
   for( int i = 0; i < colvec.getLength(); ++i )
   {
      if( rowFlags[inds[i]].test( RowFlag::kClique ))
      {
         const auto rowvec = matrix.getRowCoefficients(inds[i]);
         neighbournumber += rowvec.getLength() - 1;
      }
   }
   return neighbournumber ;
}

template <typename REAL>
bool
CliqueMerging<REAL>::isNeighbour( const ConstraintMatrix<REAL>& matrix, 
                                    int u, int v )
{
   const auto colvecu = matrix.getColumnCoefficients( u );
   const auto colvecv = matrix.getColumnCoefficients( v );
   const int* indu = colvecu.getIndices();
   const int* indv = colvecv.getIndices();
   const int lv = colvecv.getLength();
   std::vector<RowFlags> rowFlags = matrix.getRowFlags();
   int j = 0;
   for( int i = 0; i < colvecu.getLength(); ++i )
   {
      while( indv[j] < indu[i] && j < lv )
      {
         j += 1;
      }
      if( j == lv)
         break;
      if( rowFlags[indu[i]].test( RowFlag::kClique ) && indv[j] == indu[i] )
         return( true );
   }
   return( false );
}
            
template <typename REAL>
Vec<int>
CliqueMerging<REAL>::greedyClique( const ConstraintMatrix<REAL>& matrix,
                                const int cliqueRow )
{
   Vec<int> newClique;
   const auto rowvec = matrix.getRowCoefficients( cliqueRow );
   const auto indices = rowvec.getIndices();
   const int neighlength = getNeighbourhoodSize( matrix, indices[0] );
   for( int i = 0; i < neighlength; ++i )
   {
      int potNeighbour = getNeighbour( matrix, indices[0], i );
      bool isIn = true;
      bool alreadyIn = false;
      for( int j = 0; j < rowvec.getLength(); ++j )
      {
         if( indices[j] == potNeighbour )
         {
            alreadyIn = true;
            break;
         }
      }
      for( int j = 0; j < newClique.end() - newClique.begin(); ++j )
      {
         if( newClique[j] == potNeighbour )
         {
            alreadyIn = true;
            break;
         }
      }
      if( alreadyIn )
         continue;
      for( int j = 0; j < rowvec.getLength(); ++j )
      {
         if( !isNeighbour( matrix, indices[j], potNeighbour) )
         {
            isIn = false;
            break;
         }
      }
      if( isIn )
      {
         for( int j = 0; j < newClique.end() - newClique.begin(); ++j )
         {
            if( !isNeighbour( matrix, newClique[j], potNeighbour) )
            {
               isIn = false;
               break;
            }
         }
      }
      if( isIn )
      {
         newClique.push_back(potNeighbour);
      }
   }
   std::sort(newClique.begin(), newClique.end());
   return( newClique );
}

template <typename REAL>
bool
CliqueMerging<REAL>::isCovered( const ConstraintMatrix<REAL>& matrix, 
         const int row, const Vec<int> newClique, 
         const int clique)
{
   const auto rowvec1 = matrix.getRowCoefficients( clique );
   const auto indices1 = rowvec1.getIndices();
   const auto rowvec2 = matrix.getRowCoefficients( row );
   const auto indices2 = rowvec2.getIndices();
   int i1 = 0;
   int i3 = 0;
   const int l1 = rowvec1.getLength();
   const int l3 = newClique.end() - newClique.begin();
   const int l2 = rowvec2.getLength();
   for( int i2 = 0; i2 < l2; ++i2 )
   {
      while( i1 < l1 )
      {
         if( indices1[i1] >= indices2[i2] )
            break;
         i1 += 1;
      }
      while( i3 < l3 )
      {
         if( newClique[i3] >= indices2[i2] )
            break;
         i3 += 1;
      }
      if( ( i1 == l1 && i3 == l3 ) )
         return( false );
      if( i1 < l1 && i3 < l3 )
      {
         if( newClique[i3] > indices2[i2] && indices1[i1] > indices2[i2] )
            return( false );
      }
      else if ( i1 < l1 )
      {
         if( indices1[i1] > indices2[i2] )
            return( false );
      }
      else
      {
         if( newClique[i3] > indices2[i2] )
            return( false );
      }
      
   }
   return( true );
}

template <typename REAL>
PresolveStatus
CliqueMerging<REAL>::execute( const Problem<REAL>& problem,
                              const ProblemUpdate<REAL>& problemUpdate,
                              const Num<REAL>& num,
                              Reductions<REAL>& reductions, const Timer& timer,
                              int& reason_of_infeasibility )
{

   const auto& matrix = problem.getConstraintMatrix();

   std::vector<RowFlags> rowFlags = matrix.getRowFlags();

   const int nrows = matrix.getNRows();

   PresolveStatus result = PresolveStatus::kUnchanged;

   Vec<int> Cliques;

   for( int row = 0; row < nrows; ++row )
   {
      std::pair<bool, bool> cliqueCheck =
          problem.is_clique_or_sos1( matrix, row, num );
      if( cliqueCheck.first & !cliqueCheck.second )
      {
         assert( matrix.getRowCoefficients(row).getLength() != 0);
         Cliques.push_back( row );
      }
   }

   Vec<bool> completedCliques(Cliques.end() - Cliques.begin(), false);
   for( int clique = 0; clique < Cliques.end() - Cliques.begin(); ++clique )
   {
      if( completedCliques[clique])
         continue;
      Vec<int> newClique = greedyClique( matrix, Cliques[clique] );
      Vec<int> coveredCliques;
      for( int cl = 0; cl < Cliques.end()- Cliques.begin(); ++cl )
      {
         if( cl != clique && isCovered( matrix, Cliques[cl], newClique, Cliques[clique] ) )
            coveredCliques.push_back( cl );
      }
      if( coveredCliques.end() - coveredCliques.begin() > 0 )
      {
         auto rowvec = matrix.getRowCoefficients( Cliques[clique] );
         auto indices = rowvec.getIndices();
         auto val = rowvec.getValues()[0];
         Vec<int> newVertices;
         for( int cl = 0; cl < coveredCliques.end() - coveredCliques.begin(); ++cl )
         {
            auto rowvector = matrix.getRowCoefficients( Cliques[coveredCliques[cl]] );
            auto rowindex = rowvector.getIndices();
            for( int vertex = 0; vertex < rowvector.getLength(); ++vertex )
            {
               bool vertexIn = false;
               for( int node = 0; node < rowvec.getLength(); ++node )
               {
                  if( indices[node] == rowindex[vertex] )
                  {
                     vertexIn = true;
                     break;
                  }
               }
               for( int node = 0; node < newVertices.end() - newVertices.begin(); ++ node )
               {
                  if( vertexIn || newVertices[node] == rowindex[vertex] )
                  {
                     vertexIn = true;
                     break;
                  }
               }
               if( !vertexIn )
                  newVertices.push_back(rowindex[vertex]);
            }
         }
         result = PresolveStatus::kReduced;
         TransactionGuard<REAL> tg{ reductions };
         reductions.lockRow(Cliques[clique]);
         for( int i = 0; i < coveredCliques.end() - coveredCliques.begin(); ++i )
         {
            reductions.lockRow(Cliques[coveredCliques[i]]);
            completedCliques[coveredCliques[i]] = true;
         }
        for( int i = 0; i < Cliques.end() - Cliques.begin(); ++i)
        {
         reductions.lockRow(Cliques[i]);
         auto r = matrix.getRowCoefficients( Cliques[i] );
         auto indi = r.getIndices();
         for( int j = 0; j < r.getLength(); ++j )
         {
            reductions.lockCol(indi[j]);
            reductions.lockColBounds(indi[j]);
         }
        }
         for( int i = 0; i < newVertices.end() - newVertices.begin(); ++i )
         {
            reductions.lockColBounds( newVertices[i] );
            reductions.lockCol( newVertices[i] );
         }
         for( int i = 0; i < rowvec.getLength(); ++i )
         {
            reductions.lockColBounds( indices[i] );
            reductions.lockCol( indices[i] );
         }
         for( int i = 0; i < coveredCliques.end() - coveredCliques.begin(); ++i )
         {
            reductions.markRowRedundant(Cliques[coveredCliques[i]]);
         }
         auto lb = problem.getLowerBounds();
         auto ub = problem.getUpperBounds();
         for( int i = 0; i < newVertices.end() - newVertices.begin(); ++i )
         {
            assert( num.isFeasGT(val * ( ub[indices[0]] - lb[indices[0]] )
                                      * ( ub[newVertices[i]] - lb[newVertices[i]] ),0.0) || num.isFeasLT( val * ( ub[indices[0]] - lb[indices[0]] )
                                      * ( ub[newVertices[i]] - lb[newVertices[i]] ), 0.0) );
            assert( !num.isEq(0.0,val * ( ub[indices[0]] - lb[indices[0]] )
                                      * ( ub[newVertices[i]] - lb[newVertices[i]] )));
            auto rv = matrix.getRowCoefficients( Cliques[i] );
            auto vl = rv.getValues();
            auto inc = rv.getIndices();
            for( int x = 0; x < rv.getLength(); ++x )
            {
               assert( inc[x] != newVertices[i] );
               assert( vl[x] == 1.0 || vl[x] == -1.0 );
            }
            reductions.changeMatrixEntry( Cliques[clique], newVertices[i], val * ( ub[indices[0]] - lb[indices[0]] )
                                      * ( ub[newVertices[i]] - lb[newVertices[i]] ) );
         }
         return result ;
      }
   }
   return result ;
}

} // namespace papilo

#endif